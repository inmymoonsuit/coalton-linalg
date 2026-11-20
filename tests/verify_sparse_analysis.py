"""Generate and independently verify exact sparse-analysis structure.

The library is never imported or executed here. Small/general cases use the
existing set-graph oracle. Large analytic families have directly derived orders
and fill edges; selftest cross-checks those formulas against the set oracle.
Timing and allocation come from a separate Lisp process; neither measures peak
memory. See tools/sparse-analysis-tests.lisp for the integer-only file protocol.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from verify import exact_symbolic

DEFAULT_SEED = 0x3A2026
ANALYTIC = {"chain", "isolated", "star", "chain-blocks", "clique-blocks",
            "bipartite-blocks", "k2r", "hub-cycles"}
OUTPUT_SUFFIXES = (".perm.out.txt", ".columns.out.txt", ".rows.out.txt")


class AnalysisMismatch(ValueError):
    """An input, exported structure, or measurement violates its contract."""


def require(condition, message):
    if not condition:
        raise AnalysisMismatch(message)


def edge(a, b):
    require(a != b, "Graph edges must be off-diagonal.")
    return max(a, b), min(a, b)


def definition(name, family, parameters, suite="small"):
    return {"name": name, "family": family, "parameters": parameters, "suite": suite}


def definitions(seed=DEFAULT_SEED):
    """A fixed roster; labels and seeded choices are part of the generated input."""
    result = []
    for n in range(5):
        for mask in range(1 << (n * (n - 1) // 2)):
            result.append(definition(f"enum-n{n}-m{mask:03x}", "enumerated", {"n": n, "mask": mask}))
    for index in range(64):
        result.append(definition(f"seeded-{index:02d}", "random",
                                 {"n": 5 + index % 28, "seed": seed + index,
                                  "sixteenths": (0, 1, 4, 8, 15, 16)[index % 6]}))
    result.extend([
        definition("chain-33", "chain", {"n": 33}),
        definition("star-33", "star", {"n": 33}),
        definition("cycle-17", "cycle", {"n": 17}),
        definition("wheel-17", "wheel", {"n": 17}),
        definition("grid-5x7", "grid", {"height": 5, "width": 7}),
        definition("postorder-tie-trap", "explicit", {"n": 7,
                   "edges": [[0, 5], [0, 6], [5, 1], [5, 3], [6, 2], [6, 4]]}),
        definition("bipartite-4x7", "bipartite", {"left": 4, "right": 7}),
        definition("hub-cycles-3", "hub-cycles", {"cycles": 3}),
    ])
    result.extend([
        definition("chain-4097", "chain", {"n": 4097}, "large"),
        definition("chain-8192", "chain", {"n": 8192}, "large"),
        definition("chain-16384", "chain", {"n": 16384}, "large"),
        definition("chain-65536", "chain", {"n": 65536}, "large"),
        definition("isolated-16384", "isolated", {"n": 16384}, "large"),
        definition("star-16384", "star", {"n": 16384}, "large"),
        definition("disconnected-chains-8192", "chain-blocks", {"block": 64, "blocks": 128}, "large"),
        definition("disconnected-cliques-8192", "clique-blocks", {"block": 16, "blocks": 512}, "large"),
        definition("bipartite-fill-4480", "bipartite-blocks", {"side": 32, "blocks": 70}, "large"),
        definition("grid-65x65", "grid", {"height": 65, "width": 65}, "large"),
        definition("wheel-4097", "wheel", {"n": 4097}, "large"),
        definition("k2r-8192", "k2r", {"right": 8192}, "large"),
        definition("hub-cycles-2048", "hub-cycles", {"cycles": 2048}, "large"),
    ])
    return result


def graph(case):
    """Construct only the declared input graph; numerical values are absent."""
    family, p = case["family"], case["parameters"]
    edges = set()
    if family in {"chain", "isolated", "star", "cycle", "wheel", "enumerated", "random", "explicit"}:
        n = p["n"]
    elif family in {"chain-blocks", "clique-blocks"}:
        n = p["block"] * p["blocks"]
    elif family == "bipartite-blocks":
        n = 2 * p["side"] * p["blocks"]
    elif family == "grid":
        n = p["height"] * p["width"]
    elif family == "bipartite":
        n = p["left"] + p["right"]
    elif family == "k2r":
        n = 2 + p["right"]
    elif family == "hub-cycles":
        n = 1 + 3 * p["cycles"]
    else:
        raise AnalysisMismatch(f"Unknown graph family: {family}")
    require(isinstance(n, int) and 0 <= n <= 1000000, "Invalid graph dimension.")

    if family == "chain":
        edges.update((i, i - 1) for i in range(1, n))
    elif family == "star":
        edges.update((i, 0) for i in range(1, n))
    elif family == "cycle":
        require(n >= 3, "A cycle needs at least three vertices.")
        edges.update((i, i - 1) for i in range(1, n))
        edges.add((n - 1, 0))
    elif family == "wheel":
        require(n >= 4, "A wheel needs at least four vertices.")
        edges.update((i, 0) for i in range(1, n))
        edges.update((i, i - 1) for i in range(2, n))
        edges.add((n - 1, 1))
    elif family == "enumerated":
        candidates = [(i, j) for i in range(n) for j in range(i)]
        require(0 <= p["mask"] < 1 << len(candidates), "Invalid graph mask.")
        edges.update(e for bit, e in enumerate(candidates) if p["mask"] & (1 << bit))
    elif family == "random":
        rng = random.Random(p["seed"])
        require(0 <= p["sixteenths"] <= 16, "Invalid edge probability.")
        edges.update((i, j) for i in range(n) for j in range(i)
                     if rng.randrange(16) < p["sixteenths"])
    elif family == "explicit":
        edges.update(edge(a, b) for a, b in p["edges"])
    elif family in {"chain-blocks", "clique-blocks"}:
        block = p["block"]
        require(block >= 1 and p["blocks"] >= 0, "Invalid disconnected block dimensions.")
        for base in range(0, n, block):
            if family == "chain-blocks":
                edges.update((i, i - 1) for i in range(base + 1, base + block))
            else:
                edges.update((i, j) for i in range(base, base + block) for j in range(base, i))
    elif family == "bipartite-blocks":
        side = p["side"]
        require(side >= 1 and p["blocks"] >= 0, "Invalid bipartite block dimensions.")
        for base in range(0, n, 2 * side):
            edges.update((base + side + j, base + i) for i in range(side) for j in range(side))
    elif family == "grid":
        height, width = p["height"], p["width"]
        require(height >= 1 and width >= 1, "Invalid grid dimensions.")
        for row in range(height):
            for column in range(width):
                i = row * width + column
                if row:
                    edges.add((i, i - width))
                if column:
                    edges.add((i, i - 1))
    elif family in {"bipartite", "k2r"}:
        left = p["left"] if family == "bipartite" else 2
        right = p["right"]
        require(left >= 0 and right >= 0, "Invalid bipartite dimensions.")
        edges.update((left + j, i) for i in range(left) for j in range(right))
    elif family == "hub-cycles":
        require(p["cycles"] >= 0, "Invalid cycle count.")
        for k in range(p["cycles"]):
            a, b, c = 3 * k + 1, 3 * k + 2, 3 * k + 3
            edges.update((edge(0, a), edge(a, b), edge(b, c), edge(c, 0)))
    require(all(0 <= column < row < n for row, column in edges), "Graph edge outside its dimension.")
    return n, edges


def csc(n, pairs, include_diagonal=False):
    columns = [[] for _ in range(n)]
    if include_diagonal:
        for i in range(n):
            columns[i].append(i)
    for row, column in pairs:
        columns[column].append(row)
    pointers, rows = [0], []
    for column in columns:
        rows.extend(sorted(column))
        pointers.append(len(rows))
    return pointers, rows


def pattern_bytes(n, edges):
    pointers, rows = csc(n, edges, include_diagonal=True)
    text = f"{n} {len(rows)}\n"
    text += " ".join(map(str, pointers)) + "\n"
    text += " ".join(map(str, rows)) + "\n"
    return text.encode("ascii")


def analytic_order_and_fill(case, n, edges):
    """Closed forms for selected families, independently cross-checked in selftest."""
    family, p = case["family"], case["parameters"]
    order = list(range(n))
    filled = set(edges)
    if family == "star" and n > 1:
        order = list(range(1, n - 1)) + [0, n - 1]
    elif family == "bipartite-blocks":
        side = p["side"]
        for base in range(0, n, 2 * side):
            filled.update((i, j) for i in range(base + side, base + 2 * side)
                          for j in range(base + side, i))
    elif family == "k2r":
        right = p["right"]
        if right >= 3:
            order = list(range(2, n - 1)) + [0, 1, n - 1]
            filled.add((1, 0))
        elif right == 2:
            filled.add((3, 2))
    elif family == "hub-cycles" and p["cycles"]:
        cycles = p["cycles"]
        order = list(range(1, 3 * (cycles - 1) + 1)) + [0, 3 * cycles - 2, 3 * cycles - 1, 3 * cycles]
        for k in range(cycles - 1):
            filled.add((3 * k + 2, 0))
        filled.add((3 * cycles, 3 * cycles - 2))
    inverse = {old: new for new, old in enumerate(order)}
    require(sorted(order) == list(range(n)), "Analytic order is not a permutation.")
    structure = {edge(inverse[row], inverse[column]) for row, column in filled}
    return order, structure


def expected_structure(case, n, edges):
    if case["family"] in ANALYTIC:
        order, structure = analytic_order_and_fill(case, n, edges)
    else:
        order, structure = exact_symbolic(n, edges)
    pointers, rows = csc(n, structure)
    return order, pointers, rows


def read_pattern(path):
    try:
        tokens = list(map(int, path.read_text(encoding="ascii").split()))
    except (ValueError, UnicodeError) as error:
        raise AnalysisMismatch(f"Non-integer pattern file: {path.name}") from error
    require(len(tokens) >= 3, f"Missing pattern header: {path.name}")
    n, nnz = tokens[:2]
    require(0 <= n <= 1000000 and n <= nnz <= 10000000, "Pattern dimensions outside reader limits.")
    require(len(tokens) == 2 + n + 1 + nnz, "Pattern token count differs from header.")
    pointers, rows = tokens[2:3 + n], tokens[3 + n:]
    require(pointers[0] == 0 and pointers[-1] == nnz, "Invalid input CSC endpoints.")
    edges = set()
    for column in range(n):
        begin, end = pointers[column:column + 2]
        require(0 <= begin < end <= nnz, "Input CSC column lacks its diagonal or has invalid pointers.")
        entries = rows[begin:end]
        require(entries[0] == column and all(column <= r < n for r in entries),
                "Input CSC diagonal or row bound is invalid.")
        require(all(a < b for a, b in zip(entries, entries[1:])), "Input CSC rows are not strictly ordered.")
        edges.update((row, column) for row in entries[1:])
    return n, edges


def read_indices(path):
    try:
        return [int(token) for token in path.read_text(encoding="ascii").split()]
    except (ValueError, UnicodeError) as error:
        raise AnalysisMismatch(f"Non-integer analysis output: {path.name}") from error


def derived_tree_checks(n, pointers, rows):
    """Check forest/ancestor identities derived from L, not an unexported production tree."""
    parent = [rows[pointers[j]] if pointers[j] < pointers[j + 1] else n for j in range(n)]
    children = [[] for _ in range(n)]
    roots = []
    for child, ancestor in enumerate(parent):
        if ancestor == n:
            roots.append(child)
        else:
            require(child < ancestor < n, "Derived parent is not a later elimination vertex.")
            children[ancestor].append(child)
    postorder = []
    stack = [(root, False) for root in reversed(roots)]
    while stack:
        vertex, visited = stack.pop()
        if visited:
            postorder.append(vertex)
        else:
            stack.append((vertex, True))
            stack.extend((child, False) for child in reversed(children[vertex]))
    require(sorted(postorder) == list(range(n)), "Derived postorder does not cover the forest.")
    position = [0] * n
    for index, vertex in enumerate(postorder):
        position[vertex] = index
    subtree_size = [1] * n
    for child, ancestor in enumerate(parent):
        if ancestor != n:
            subtree_size[ancestor] += subtree_size[child]
            require(position[child] < position[ancestor], "Derived postorder violates parent precedence.")
    for column in range(n):
        for row in rows[pointers[column]:pointers[column + 1]]:
            start = position[row] - subtree_size[row] + 1
            require(start <= position[column] <= position[row], "A factor row is not an ancestor of its column.")
    depth = [0] * n
    for vertex in range(n - 1, -1, -1):
        if parent[vertex] != n:
            depth[vertex] = depth[parent[vertex]] + 1
    encoded = ",".join(map(str, postorder)).encode("ascii")
    return {"derived_tree_roots": len(roots), "derived_tree_max_depth": max(depth, default=0),
            "derived_postorder_sha256": hashlib.sha256(encoded).hexdigest()}


def validate_structure(n, actual, expected):
    permutation, pointers, rows = actual
    require(len(permutation) == n and sorted(permutation) == list(range(n)),
            "Permutation is not a bijection.")
    require(len(pointers) == n + 1 and pointers[0] == 0 and pointers[-1] == len(rows),
            "CSC pointers have invalid length or endpoints.")
    for column in range(n):
        begin, end = pointers[column:column + 2]
        require(0 <= begin <= end <= len(rows), "CSC pointers are not monotone/in bounds.")
        entries = rows[begin:end]
        require(all(column < row < n for row in entries), "Factor rows are not strictly below the diagonal.")
        require(all(a < b for a, b in zip(entries, entries[1:])), "Factor rows are not strictly increasing.")
    require(permutation == expected[0], "Permutation differs from exact minimum-degree tie order.")
    require(pointers == expected[1] and rows == expected[2], "Factor pattern differs from independent exact fill.")
    return derived_tree_checks(n, pointers, rows)


def read_metrics(path, n, input_nnz, strict_nnz):
    result = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        require(": " in line, "Malformed analysis measurement line.")
        key, value = line.split(": ", 1)
        require(key not in result, "Duplicate analysis measurement key.")
        result[key] = value
    require(result.get("format") == "1" and result.get("status") == "success", "Analysis did not report success.")
    for key, expected in (("size", n), ("input_nnz", input_nnz), ("strict_l_nnz", strict_nnz)):
        require(int(result[key]) == expected, f"Analysis metadata {key} differs from exported structure.")
    ticks = int(result["analysis_ticks"])
    units = int(result["internal_time_units_per_second"])
    allocated = int(result["analysis_bytes_consed"])
    seconds = float(result["analysis_seconds"])
    require(ticks >= 0 and units > 0 and allocated >= 0 and math.isfinite(seconds),
            "Invalid analysis measurement.")
    require(abs(seconds - ticks / units) <= 1e-8, "Analysis seconds do not match the recorded clock ticks.")
    require("not peak memory" in result["measurement_scope"], "Allocation scope is missing.")
    return {"analysis_seconds": seconds, "analysis_ticks": ticks,
            "internal_time_units_per_second": units, "analysis_bytes_consed": allocated}


def selected(cases, names):
    if not names:
        return cases
    requested = set(names)
    found = [case for case in cases if case["name"] in requested]
    require({case["name"] for case in found} == requested, "Requested case is absent from this suite.")
    return found


def generate(folder, suite, seed, names):
    cases = [case for case in definitions(seed) if suite == "all" or case["suite"] == suite]
    cases = selected(cases, names)
    require(bool(cases), "Empty analysis suite.")
    folder.mkdir(parents=True, exist_ok=True)
    generated = []
    for case in cases:
        n, edges = graph(case)
        contents = pattern_bytes(n, edges)
        (folder / (case["name"] + ".pattern.txt")).write_bytes(contents)
        generated.append({**case, "size": n, "input_nnz": n + len(edges),
                          "input_sha256": hashlib.sha256(contents).hexdigest(),
                          "oracle": "analytic-order-and-fill" if case["family"] in ANALYTIC else "exact-set-graph"})
    manifest = {"format": 1, "generator_revision": 1, "seed": seed, "cases": generated}
    (folder / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"SPARSE-ANALYSIS-GENERATED cases={len(cases)} above4096={sum(c['size'] > 4096 for c in generated)}")


def verify(folder, output, names=None):
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("format") == 1 and manifest.get("generator_revision") == 1, "Unsupported analysis manifest.")
    all_cases = manifest["cases"]
    roster = [case["name"] for case in all_cases]
    require(len(roster) == len(set(roster)) and all(re.fullmatch(r"[a-z0-9_-]+", name) for name in roster),
            "Duplicate or unsafe analysis case name.")
    require({p.name for p in folder.glob("*.pattern.txt")} == {name + ".pattern.txt" for name in roster},
            "Input pattern roster differs from manifest; use a separate directory for a different suite.")
    cases = selected(all_cases, names)
    require(bool(cases), "Empty analysis verification scope.")
    expected_outputs = {case["name"] + suffix for case in cases for suffix in OUTPUT_SUFFIXES}
    require({p.name for p in output.glob("*.out.txt")} == expected_outputs,
            "Analysis output roster differs from the requested case roster.")
    expected_meta = {case["name"] + ".analysis.meta.txt" for case in cases}
    require({p.name for p in output.glob("*.analysis.meta.txt")} == expected_meta,
            "Analysis metadata roster differs from the requested case roster.")
    records = []
    for case in cases:
        name = case["name"]
        path = folder / (name + ".pattern.txt")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == case["input_sha256"], f"Input bytes changed: {name}")
        n, edges = read_pattern(path)
        declared_n, declared_edges = graph(case)
        require(n == declared_n == case["size"] and edges == declared_edges, f"Input graph differs from declared family: {name}")
        require(n + len(edges) == case["input_nnz"], f"Input count differs: {name}")
        expected = expected_structure(case, n, edges)
        actual = tuple(read_indices(output / (name + suffix)) for suffix in OUTPUT_SUFFIXES)
        tree = validate_structure(n, actual, expected)
        metrics = read_metrics(output / (name + ".analysis.meta.txt"), n, n + len(edges), len(actual[2]))
        records.append({"name": name, "family": case["family"], "size": n,
                        "input_nnz": n + len(edges), "strict_l_nnz": len(actual[2]),
                        "fill_edges": len(actual[2]) - len(edges), "oracle": case["oracle"], **tree, **metrics})
    return {"format": 1, "cases": len(records), "above4096": sum(row["size"] > 4096 for row in records),
            "measurement_scope": "Lisp analyze time and program-wide allocation; no peak-memory measurement",
            "tree_scope": "Derived from exported L; no private production tree or postorder is inspected",
            "results": records}


def selftest():
    crosschecks = 0
    families = []
    for n in range(13):
        for family in ("chain", "isolated", "star"):
            families.append(definition("formula", family, {"n": n}))
    for block in (1, 2, 3, 7):
        for blocks in (0, 1, 2, 3):
            for family in ("chain-blocks", "clique-blocks"):
                families.append(definition("formula", family, {"block": block, "blocks": blocks}))
    for side in (1, 2, 3, 5):
        for blocks in (0, 1, 2, 3):
            families.append(definition("formula", "bipartite-blocks", {"side": side, "blocks": blocks}))
    for count in range(9):
        families.append(definition("formula", "k2r", {"right": count}))
        families.append(definition("formula", "hub-cycles", {"cycles": count}))
    for case in families:
        n, edges = graph(case)
        analytic = analytic_order_and_fill(case, n, edges)
        reference = exact_symbolic(n, edges)
        require(analytic == reference, f"Analytic family differs from independent set oracle: {case}")
        pointers, rows = csc(n, analytic[1])
        derived_tree_checks(n, pointers, rows)
        crosschecks += 1
    for case in definitions():
        if case["suite"] == "small":
            n, edges = graph(case)
            expected = expected_structure(case, n, edges)
            validate_structure(n, expected, expected)

    controls = []

    def rejected(name, n, actual, expected, message):
        try:
            validate_structure(n, actual, expected)
        except AnalysisMismatch as error:
            require(message in str(error), f"{name} failed for an unintended reason: {error}")
            controls.append({"name": name, "rejected": str(error)})
        else:
            raise AssertionError(f"Planted defect passed: {name}")

    case = definition("tie", "isolated", {"n": 3})
    n, edges = graph(case)
    good = expected_structure(case, n, edges)
    rejected("wrong-minimum-degree-tie", n, ([1, 0, 2], good[1], good[2]), good, "Permutation differs")

    case = definition("fill", "cycle", {"n": 4})
    n, edges = graph(case)
    good = expected_structure(case, n, edges)
    order, structure = exact_symbolic(n, edges)
    require((3, 1) in structure and (3, 1) not in edges, "Cycle control has no planted fill edge.")
    pointers, rows = csc(n, structure - {(3, 1)})
    rejected("missing-fill-edge", n, (order, pointers, rows), good, "Factor pattern differs")

    case = definition("extra", "chain", {"n": 4})
    n, edges = graph(case)
    good = expected_structure(case, n, edges)
    pointers, rows = csc(n, edges | {(2, 0)})
    rejected("extra-fill-edge", n, (good[0], pointers, rows), good, "Factor pattern differs")
    broken_pointers = good[1][:]
    broken_pointers[0] = 1
    rejected("malformed-column-pointers", n, (good[0], broken_pointers, good[2]), good, "CSC pointers")
    return {"format": 1, "analytic_crosschecks": crosschecks, "planted_controls": controls,
            "small_cases": sum(case["suite"] == "small" for case in definitions())}


def write_report(path, report):
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    generate_parser = commands.add_parser("generate")
    generate_parser.add_argument("folder", type=Path)
    generate_parser.add_argument("--suite", choices=("small", "large", "all"), default="all")
    generate_parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    generate_parser.add_argument("--case", action="append")
    verify_parser = commands.add_parser("verify")
    verify_parser.add_argument("input", type=Path)
    verify_parser.add_argument("output", type=Path)
    verify_parser.add_argument("--case", action="append")
    verify_parser.add_argument("--json", type=Path)
    selftest_parser = commands.add_parser("selftest")
    selftest_parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "generate":
            generate(args.folder, args.suite, args.seed, args.case)
        elif args.command == "verify":
            report = verify(args.input, args.output, args.case)
            write_report(args.json, report)
            print(f"SPARSE-ANALYSIS-VERIFIED cases={report['cases']} above4096={report['above4096']}")
        else:
            report = selftest()
            write_report(args.json, report)
            print(f"SPARSE-ANALYSIS-CONTROLS-PASS analytic={report['analytic_crosschecks']} controls={len(report['planted_controls'])} small={report['small_cases']}")
        return 0
    except (AnalysisMismatch, AssertionError, OSError, ValueError, KeyError, TypeError) as error:
        print(f"SPARSE-ANALYSIS-FAIL: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())