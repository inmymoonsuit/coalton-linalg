"""Independently verify sparse LDL output artifacts without importing the implementation.

All numerical checks are derived from Ax=b and PAP^T=LDL^T.  The small
fixture oracle uses exact Fraction arithmetic.  No sparse-library code is used.
"""
from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
import shutil
import sys


class VerificationError(Exception):
    """An artifact does not satisfy its declared format or numerical contract."""


def require(condition, message):
    if not condition:
        raise VerificationError(message)


def finite_float(token):
    value = float(token.replace("D", "e").replace("d", "e"))
    require(math.isfinite(value), "non-finite numeric artifact")
    return value


def read_vector(path, n, integers=False):
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    require(len(lines) == n, f"{path.name}: expected {n} lines, got {len(lines)}")
    require(all(len(line.split()) == 1 for line in lines), f"{path.name}: not a vector")
    return [int(x) if integers else finite_float(x) for x in lines]


def read_triplets(path, strict=False):
    lines = [line.split() for line in path.read_text().splitlines() if line.strip()]
    require(bool(lines) and len(lines[0]) == 2, f"{path.name}: malformed header")
    n, nnz = map(int, lines[0])
    require(n >= 1 and nnz >= 0, f"{path.name}: invalid dimensions")
    require(len(lines) - 1 == nnz, f"{path.name}: entry count mismatch")
    entries = {}
    previous = (-1, -1)
    for line in lines[1:]:
        require(len(line) == 3, f"{path.name}: malformed triplet")
        row, col = map(int, line[:2])
        value = finite_float(line[2])
        require(0 <= col <= row < n, f"{path.name}: entry outside lower triangle")
        require(not strict or row > col, f"{path.name}: expected strict lower triangle")
        require((col, row) > previous, f"{path.name}: unsorted or duplicate entry")
        previous = (col, row)
        entries[row, col] = value
    if not strict:
        require(all((i, i) in entries for i in range(n)), f"{path.name}: missing diagonal")
    return n, entries


def metadata(path):
    if not path.exists():
        return {}
    return dict(line.split(": ", 1) for line in path.read_text().splitlines() if ": " in line)


def exact_symbolic(n, pattern):
    # Derived independently from elimination-graph adjacency updates.
    neighbours = [set() for _ in range(n)]
    for row, col in pattern:
        if row != col:
            neighbours[row].add(col)
            neighbours[col].add(row)
    live = set(range(n))
    order = []
    frontiers = []
    while live:
        vertex = min(live, key=lambda i: (len(neighbours[i]), i))
        adjacent = sorted(neighbours[vertex])
        order.append(vertex)
        frontiers.append(adjacent)
        for i in adjacent:
            neighbours[i].remove(vertex)
            neighbours[i].update(j for j in adjacent if j != i)
        live.remove(vertex)
    inverse = {old: new for new, old in enumerate(order)}
    structure = {(inverse[row], col) for col, adjacent in enumerate(frontiers) for row in adjacent}
    return order, structure


def exact_factor(n, entries, shift=0):
    """Small independent dense rational oracle; returns order,L,D or refusal."""
    require(n <= 32, "exact refusal oracle is limited to 32 rows")
    order, _ = exact_symbolic(n, entries)
    matrix = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for (r, c), value in entries.items():
        matrix[r][c] = matrix[c][r] = Fraction(value)
    for i in range(n):
        matrix[i][i] -= Fraction(shift)
    lower = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    diagonal = []
    for j in range(n):
        pivot = matrix[order[j]][order[j]] - sum(
            lower[j][k] ** 2 * diagonal[k] for k in range(j)
        )
        if pivot <= 0:
            return order, lower, diagonal, order[j]
        diagonal.append(pivot)
        for i in range(j + 1, n):
            lower[i][j] = (
                matrix[order[i]][order[j]]
                - sum(lower[i][k] * lower[j][k] * diagonal[k] for k in range(j))
            ) / pivot
    return order, lower, diagonal, None


def check_factor(n, entries, output, prefix, shift=0.0):
    order = read_vector(output / f"{prefix}.perm.out.txt", n, integers=True)
    require(sorted(order) == list(range(n)), f"{prefix}: permutation is not a bijection")
    size, lower = read_triplets(output / f"{prefix}.L.out.txt", strict=True)
    require(size == n, f"{prefix}: factor size differs from A")
    if n <= 32:
        expected_order, expected_structure = exact_symbolic(n, entries)
        require(order == expected_order, f"{prefix}: wrong minimum-degree order or tie break")
        require(set(lower) == expected_structure, f"{prefix}: symbolic factor structure differs, including structural zeros")
    diagonal = read_vector(output / f"{prefix}.D.out.txt", n)
    require(all(value > 0 for value in diagonal), f"{prefix}: D is not strictly positive")
    inverse = [0] * n
    for new, old in enumerate(order):
        inverse[old] = new

    # Keep every produced entry, including cancellations outside A's structure.
    difference = {}
    for (old_r, old_c), value in entries.items():
        r, c = sorted((inverse[old_r], inverse[old_c]), reverse=True)
        difference[r, c] = -value
    for j in range(n):
        difference[j, j] += shift
    columns = [[(j, 1.0)] for j in range(n)]
    for (row, col), value in lower.items():
        columns[col].append((row, value))
    for j, column in enumerate(columns):
        for a, (row, value) in enumerate(column):
            scaled = value * diagonal[j]
            for col, other in column[:a + 1]:
                key = (row, col)
                difference[key] = difference.get(key, 0.0) + scaled * other
    require(all(math.isfinite(value) for value in difference.values()),
            f"{prefix}: non-finite reconstruction intermediate")
    error = max(map(abs, difference.values()), default=0.0)
    scale = max(map(abs, entries.values()), default=0.0)
    require(math.isfinite(error), f"{prefix}: reconstruction overflowed")
    require(error <= 1e-12 * scale,
            f"{prefix}: reconstruction {error:.17g} > {1e-12 * scale:.17g}")
    return {"nnzL": len(lower) + n, "rebuild": error,
            "rebuild_scaled": error / scale if scale else error,
            "order": order}


def refusal_check(n, entries, metadata_values, output, prefix, shift=0.0):
    refusal_path = output / (f"{prefix}.out.txt" if prefix.endswith(".shift-refused") else f"{prefix}.refused.out.txt")
    refusal = read_vector(refusal_path, 1, integers=True)[0]
    require(0 <= refusal < n, f"{prefix}: refusal is not a caller column")
    if n <= 32:
        _, _, _, expected = exact_factor(n, entries, shift)
        require(expected is not None, f"{prefix}: refused a positive definite matrix")
        require(refusal == expected,
                f"{prefix}: refused caller column {refusal}, exact oracle expects {expected}")
    elif "smallest eigenvalue" in metadata_values:
        require(float(metadata_values["smallest eigenvalue"]) <= shift,
                f"{prefix}: refused below the independently recorded smallest eigenvalue")
    else:
        raise VerificationError(f"{prefix}: cannot independently validate a large refusal")
    return refusal


def verify_fixture(fixture, output, name, require_shifts=False, enforce_fill=False):
    n, entries = read_triplets(fixture / f"{name}.A.txt")
    b = read_vector(fixture / f"{name}.b.txt", n)
    meta = metadata(fixture / f"{name}.meta.txt")
    files = [output / f"{name}.{suffix}.out.txt" for suffix in ("x", "perm", "L", "D")]
    refused = output / f"{name}.refused.out.txt"
    if refused.exists():
        require(not any(path.exists() for path in files), f"{name}: stale factor artifacts beside refusal")
        value = refusal_check(n, entries, meta, output, name)
        return {"name": name, "n": n, "refused": value, "passed": True}
    require(meta.get("expected") != "refused", f"{name}: expected a refusal")
    result = check_factor(n, entries, output, name)
    x = read_vector(files[0], n)
    products = [[] for _ in range(n)]
    row_abs = [[] for _ in range(n)]
    for (r, c), value in entries.items():
        products[r].append(value * x[c])
        row_abs[r].append(abs(value))
        if r != c:
            products[c].append(value * x[r])
            row_abs[c].append(abs(value))
    residual = max(abs(math.fsum(products[i]) - b[i]) for i in range(n))
    anorm = max(math.fsum(row) for row in row_abs)
    xnorm = max(map(abs, x))
    bnorm = max(map(abs, b))
    scale = anorm * xnorm + bnorm
    require(math.isfinite(residual) and math.isfinite(scale), f"{name}: residual check overflowed")
    require(residual <= 1e-12 * scale, f"{name}: residual {residual:.17g} > {1e-12 * scale:.17g}")
    result.update(name=name, n=n, nnzA=len(entries), residual=residual,
                  backward_error=residual / scale if scale else residual)
    reference = fixture / f"{name}.x.txt"
    if reference.exists():
        target = read_vector(reference, n)
        result["reference_x_relative"] = max(abs(a - b) for a, b in zip(x, target)) / max(
            max(map(abs, target)), sys.float_info.min)
    if "reference fill" in meta:
        target_fill = int(meta["reference fill"].split()[0])
        result["fill_ratio"] = result["nnzL"] / target_fill
        result["fill_target_met"] = result["fill_ratio"] <= 1.25
        if enforce_fill:
            require(result["fill_target_met"], f"{name}: fill ratio {result['fill_ratio']:.6g} > 1.25")
    kept_key = "shift that must still factorise with every pivot positive"
    refused_key = "shift that must be refused"
    if kept_key in meta:
        prefix = f"{name}.shift-kept"
        if (output / f"{prefix}.D.out.txt").exists():
            kept = check_factor(n, entries, output, prefix, float(meta[kept_key]))
            result["shift_kept"] = True
            result["shift_rebuild"] = kept["rebuild"]
        else:
            require(not require_shifts, f"{name}: missing independent kept-shift artifacts")
            result["shift_kept"] = "not independently checked"
    if refused_key in meta:
        prefix = f"{name}.shift-refused"
        if (output / f"{prefix}.out.txt").exists():
            result["shift_refused"] = refusal_check(n, entries, meta, output, prefix, float(meta[refused_key]))
            result["shift_refusal_column_check"] = "exact rational" if n <= 32 else "range only; matrix indefiniteness from fixture eigenvalue"
        else:
            require(not require_shifts, f"{name}: missing independent refused-shift artifact")
            result["shift_refused"] = "not independently checked"
    result.pop("order")
    result["passed"] = True
    return result


def write_vector(path, values):
    path.write_text("".join(f"{float(x):.17g}\n" for x in values))


def write_fixture(folder, name, matrix, x, structural=(), shifts=None):
    n = len(matrix)
    pattern = {(i, j) for i in range(n) for j in range(i + 1) if matrix[i][j] != 0 or i == j}
    pattern.update(structural)
    entries = {(i, j): Fraction(matrix[i][j]) for i, j in pattern}
    ordered = sorted(entries, key=lambda key: (key[1], key[0]))
    text = f"{n} {len(entries)}\n" + "".join(
        f"{i} {j} {float(entries[i, j]):.17g}\n" for i, j in ordered)
    (folder / f"{name}.A.txt").write_text(text)
    b = [sum(Fraction(matrix[i][j]) * x[j] for j in range(n)) for i in range(n)]
    write_vector(folder / f"{name}.b.txt", b)
    order, lower, diagonal, refusal = exact_factor(n, entries)
    meta = [f"name: {name}", f"expected: {'success' if refusal is None else 'refused'}"]
    if refusal is None:
        write_vector(folder / f"{name}.x.txt", x)
    else:
        meta.append(f"expected refusal: {refusal}")
    if shifts:
        meta.extend((f"shift that must still factorise with every pivot positive: {shifts[0]:.17g}",
                     f"shift that must be refused: {shifts[1]:.17g}"))
    (folder / f"{name}.meta.txt").write_text("\n".join(meta) + "\n")
    return {"name": name, "n": n, "expected": "success" if refusal is None else "refused",
            "refused_column": refusal, "order": order}


def generate(folder):
    folder.mkdir(parents=True, exist_ok=True)
    cases = []
    def add(name, matrix, x=None, **kwargs):
        if x is None:
            x = [Fraction((i % 3) - 1, 2) for i in range(len(matrix))]
        cases.append(write_fixture(folder, name, matrix, x, **kwargs))
    add("scalar", [[4]], [Fraction(3, 2)], shifts=(1.0, 5.0))
    add("scalarzero", [[0]], [Fraction(0)])
    add("scalarnegative", [[-2]], [Fraction(0)])
    add("diagonal", [[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 4, 0], [0, 0, 0, 8]],
        shifts=(0.5, 1.5))
    add("structuralzero", [[4, 0, 0], [0, 3, 1], [0, 1, 2]], structural=((1, 0), (2, 0)))
    path = [[4, -1, 0, 0, 0], [-1, 4, -1, 0, 0], [0, -1, 4, -1, 0],
            [0, 0, -1, 4, -1], [0, 0, 0, -1, 4]]
    add("path", path, shifts=(1.0, 5.0))
    relabel = [2, 4, 0, 3, 1]
    add("path_relabelled", [[path[i][j] for j in relabel] for i in relabel],
        [Fraction((i % 3) - 1, 2) for i in relabel])
    star = [[5 if i == j == 0 else 2 if i == j else 1 if min(i, j) == 0 else 0
             for j in range(5)] for i in range(5)]
    add("star", star)
    add("star_refusal", [[2 if i == j == 0 else 2 if i == j else 1 if min(i, j) == 0 else 0
                         for j in range(5)] for i in range(5)])
    add("singularpair", [[1, 1], [1, 1]])
    add("indefinitepair", [[1, 2], [2, 1]])
    add("zero_leading", [[0, 1], [1, 2]])
    add("near_singular", [[1, 1], [1, 1 + Fraction(1, 2**40)]], [Fraction(1), Fraction(-1)])
    add("tiny_scale", [[Fraction(4, 2**400), Fraction(1, 2**400)],
                       [Fraction(1, 2**400), Fraction(2, 2**400)]])
    add("large_scale", [[4 * 2**400, 2**400], [2**400, 2 * 2**400]])
    add("disconnected", [[4, 1, 0, 0], [1, 2, 0, 0], [0, 0, 3, -1], [0, 0, -1, 2]])
    cycle = [[4 if i == j else -1 if (i - j) % 7 in (1, 6) else 0
              for j in range(7)] for i in range(7)]
    add("cycle7", cycle)
    reindex = [4, 0, 6, 2, 5, 1, 3]
    add("cycle7_relabelled", [[cycle[i][j] for j in reindex] for i in reindex],
        [Fraction((i % 3) - 1, 2) for i in reindex])
    grid = [[5 if i == j else -1 if abs(i // 3 - j // 3) + abs(i % 3 - j % 3) == 1 else 0
             for j in range(9)] for i in range(9)]
    add("grid3x3", grid)
    seed_lower = [[Fraction(1) if i == j else Fraction(((i * 5 + j * 3) % 7) - 3, 8)
                   if i > j else Fraction(0) for j in range(8)] for i in range(8)]
    dense = [[sum(seed_lower[i][k] * (k + 1) * seed_lower[j][k] for k in range(8))
              for j in range(8)] for i in range(8)]
    add("dense_dyadic8", dense)
    (folder / "manifest.json").write_text(json.dumps(cases, indent=2) + "\n")
    return cases


def write_oracle_outputs(fixtures, output, name, shift=0.0, suffix=""):
    n, entries = read_triplets(fixtures / f"{name}.A.txt")
    order, lower, diagonal, refusal = exact_factor(n, entries, shift)
    prefix = name + suffix
    output.mkdir(parents=True, exist_ok=True)
    if refusal is not None:
        (output / (f"{prefix}.out.txt" if suffix == ".shift-refused" else f"{prefix}.refused.out.txt")).write_text(f"{refusal}\n")
        return
    (output / f"{prefix}.perm.out.txt").write_text("".join(f"{i}\n" for i in order))
    _, structure = exact_symbolic(n, entries)
    nonzero = [(i, j, lower[i][j]) for i, j in sorted(structure, key=lambda p: (p[1], p[0]))]
    (output / f"{prefix}.L.out.txt").write_text(
        f"{n} {len(nonzero)}\n" + "".join(f"{i} {j} {float(v):.17g}\n" for i, j, v in nonzero))
    write_vector(output / f"{prefix}.D.out.txt", diagonal)
    if not suffix:
        shutil.copyfile(fixtures / f"{name}.x.txt", output / f"{prefix}.x.out.txt")


def selftest(folder):
    cases = generate(folder / "fixtures")
    fixtures = folder / "fixtures"
    output = folder / "oracle"
    for case in cases:
        name = case["name"]
        write_oracle_outputs(fixtures, output, name)
        meta = metadata(fixtures / f"{name}.meta.txt")
        for key, suffix in (
            ("shift that must still factorise with every pivot positive", ".shift-kept"),
            ("shift that must be refused", ".shift-refused"),
        ):
            if key in meta:
                write_oracle_outputs(fixtures, output, name, float(meta[key]), suffix)
        verify_fixture(fixtures, output, name, require_shifts=True)

    # Plant corrupt artifacts in separate retained directories and require red checks.
    mutations = [
        ("wrong_solution", "path", "x", lambda text: "999\n" + "\n".join(text.splitlines()[1:]) + "\n"),
        ("wrong_diagonal", "path", "D", lambda text: "999\n" + "\n".join(text.splitlines()[1:]) + "\n"),
        ("duplicate_permutation", "path", "perm", lambda text: text.splitlines()[1] + "\n" + "\n".join(text.splitlines()[1:]) + "\n"),
        ("nonfinite", "path", "D", lambda text: "nan\n" + "\n".join(text.splitlines()[1:]) + "\n"),
        ("wrong_refusal", "star_refusal", "refused", lambda text: "1\n"),
        ("false_refusal", "scalar", "refused", lambda text: "0\n"),
        ("missing_structural_zero", "structuralzero", "L",
         lambda text: text.splitlines()[0].split()[0] + " " + str(int(text.splitlines()[0].split()[1]) - 1)
         + "\n" + "\n".join(text.splitlines()[2:]) + "\n"),
        ("negative_pivot", "path", "D", lambda text: "-1\n" + "\n".join(text.splitlines()[1:]) + "\n"),
        ("unshifted_kept_factor", "scalar", "shift-kept.D", lambda text: "4\n"),
        ("overflowed_reconstruction", "path", "L",
         lambda text: text.splitlines()[0] + "\n" + " ".join(text.splitlines()[1].split()[:2])
         + " 1e308\n" + "\n".join(text.splitlines()[2:]) + "\n"),
    ]
    evidence = []
    for label, name, suffix, mutate in mutations:
        destination = folder / label
        destination.mkdir(parents=True, exist_ok=True)
        for path in output.glob(f"{name}.*.out.txt"):
            shutil.copyfile(path, destination / path.name)
        target = destination / f"{name}.{suffix}.out.txt"
        target.write_text(mutate(target.read_text() if target.exists() else ""))
        if label == "false_refusal":
            for field in ("x", "perm", "L", "D"):
                (destination / f"{name}.{field}.out.txt").unlink()
        try:
            verify_fixture(fixtures, destination, name, require_shifts=True)
        except VerificationError as exc:
            evidence.append({"control": label, "result": "rejected", "reason": str(exc)})
        else:
            raise AssertionError(f"planted defect escaped: {label}")
    (folder / "evidence.json").write_text(json.dumps(
        {"oracle_cases_passed": len(cases), "planted_controls": evidence}, indent=2) + "\n")
    print(json.dumps({"oracle_cases_passed": len(cases), "planted_controls_rejected": len(evidence)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    gen = sub.add_parser("generate")
    gen.add_argument("folder", type=Path)
    own = sub.add_parser("selftest")
    own.add_argument("folder", type=Path)
    check = sub.add_parser("verify")
    check.add_argument("fixtures", type=Path)
    check.add_argument("output", type=Path)
    check.add_argument("--name", action="append", help="Fixture basename; repeat to select a subset.")
    check.add_argument("--require-shifts", action="store_true")
    check.add_argument("--enforce-fill", action="store_true")
    check.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.command == "generate":
        print(json.dumps(generate(args.folder)))
    elif args.command == "selftest":
        selftest(args.folder)
    else:
        names = args.name or [path.name[:-6] for path in sorted(args.fixtures.glob("*.A.txt"))]
        require(bool(names), "no fixtures selected")
        results = []
        for name in names:
            try:
                result = verify_fixture(args.fixtures, args.output, name,
                                        args.require_shifts, args.enforce_fill)
            except (VerificationError, ValueError, OSError, OverflowError) as exc:
                result = {"name": name, "passed": False, "error": str(exc)}
            results.append(result)
            print(json.dumps(result, sort_keys=True))
        if args.json:
            args.json.write_text(json.dumps(results, indent=2) + "\n")
        return 0 if all(result["passed"] for result in results) else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
