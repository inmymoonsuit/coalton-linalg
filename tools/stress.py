"""Seeded independent factorization checks; numerical oracle only, no runtime dependency.

Generate small exact/dyadic cases, verify emitted factors with verify.py, and check
known solutions plus permutation/power-of-two relations. Floating-point boundary
observations belong to tests/property-tests.lisp, not this strict exact-oracle lane.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import verify as oracle


SCHEMA = 1
DEFAULT_SEEDS = 12
FAMILIES = ("dominant", "disconnected", "known_ldl", "cancellation",
            "near_singular", "zero_diagonal", "indefinite")
EXPONENTS = (-300, 300)


class Stream:
    """A specified 32-bit recurrence makes fixture selection version-independent."""

    def __init__(self, seed):
        self.state = seed & 0xFFFFFFFF

    def below(self, count):
        self.state = (1664525 * self.state + 1013904223) & 0xFFFFFFFF
        return self.state % count


def write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def digest(path):
    # Text normalization allows LF and CRLF checkouts without changing fixtures.
    return hashlib.sha256(path.read_text(encoding="utf-8").replace("\r\n", "\n")
                          .encode("utf-8")).hexdigest()


def dyadic_power(exponent):
    return Fraction(2 ** exponent) if exponent >= 0 else Fraction(1, 2 ** -exponent)


def dense_product(lower, diagonal):
    n = len(diagonal)
    return [[sum(lower[i][k] * diagonal[k] * lower[j][k] for k in range(n))
             for j in range(n)] for i in range(n)]


def make_family(family, seed):
    rng = Stream(seed + 17 * (FAMILIES.index(family) + 1))
    n = 3 + seed % 10
    structural = set()
    details = {}
    if family in ("dominant", "disconnected"):
        a = [[Fraction(0) for _ in range(n)] for _ in range(n)]
        split = max(1, n // 2)
        for i in range(n):
            for j in range(i):
                same_block = (i < split) == (j < split)
                if (family == "dominant" or same_block) and rng.below(3) == 0:
                    value = Fraction(rng.below(7) - 3, 4)
                    a[i][j] = a[j][i] = value
                    structural.add((i, j))
        for i in range(n):
            a[i][i] = sum(abs(x) for x in a[i]) + 2 + Fraction(rng.below(4), 4)
        details["strict_diagonal_dominance_margin_at_least"] = 2
    elif family == "known_ldl":
        n = 3 + seed % 8
        lower = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(i):
                if rng.below(3) != 0:
                    lower[i][j] = Fraction(rng.below(5) - 2, 32)
        diagonal = [Fraction(4 + rng.below(5), 4) for _ in range(n)]
        a = dense_product(lower, diagonal)
        details["construction"] = "exact dyadic unit-lower times positive D times transpose"
    elif family == "cancellation":
        n = 3 + seed % 4
        lower = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
        lower[1][0] = Fraction(1 if seed % 2 else -1, 2)
        lower[2][0] = Fraction(1, 2)
        lower[2][1] = -lower[2][0] * lower[1][0]
        for i in range(3, n):
            lower[i][i - 1] = Fraction(1 if i % 2 else -1, 8)
        a = dense_product(lower, [Fraction(1)] * n)
        assert a[2][1] == 0
        structural.add((2, 1))
        details["exact_cancellation_entry"] = [2, 1]
    elif family == "near_singular":
        n = 2
        exponent = 20 + seed % 12
        delta = dyadic_power(-exponent)
        a = [[Fraction(1), Fraction(1)], [Fraction(1), 1 + delta]]
        details["positive_margin_exponent"] = -exponent
        details["boundary_class"] = "representable positive margin, away from input rounding tie"
    else:
        a = [[Fraction(int(i == j) * (1 + i % 3)) for j in range(n)] for i in range(n)]
        if family == "zero_diagonal":
            a[seed % n][seed % n] = Fraction(0)
        elif family == "indefinite":
            a[0][0] = a[1][1] = Fraction(1)
            a[0][1] = a[1][0] = Fraction(2)
        else:
            raise ValueError(f"unknown family {family}")
    x = [Fraction((1 if (i + seed) % 2 else -1) * (1 + (i * 3 + seed) % 4), 4)
         for i in range(n)]
    # A conservative test budget, not an advertised condition-number estimate.
    tolerance = (128 * sys.float_info.epsilon * 2 ** (20 + seed % 12)
                 if family == "near_singular" else 2e-10)
    return a, x, structural, details, tolerance


def permutation(n, seed):
    indices = list(range(n))
    stream = Stream(seed ^ 0xA534D7)
    for end in range(n - 1, 0, -1):
        at = stream.below(end + 1)
        indices[end], indices[at] = indices[at], indices[end]
    if indices == list(range(n)) and n > 1:
        indices = indices[1:] + indices[:1]
    return indices


def transformed(matrix, solution, structural, order, exponent):
    scale = dyadic_power(exponent)
    inverse = {old: new for new, old in enumerate(order)}
    converted = {tuple(sorted((inverse[i], inverse[j]), reverse=True)) for i, j in structural}
    return ([[scale * matrix[i][j] for j in order] for i in order],
            [solution[i] for i in order], converted)


def generate(folder, seeds=DEFAULT_SEEDS):
    oracle.require(1 <= seeds <= 1000, "seed count must be between 1 and 1000")
    folder.mkdir(parents=True, exist_ok=True)
    cases = []
    for seed in range(seeds):
        for family in FAMILIES:
            matrix, x, structural, details, tolerance = make_family(family, seed)
            base = f"s{seed:04d}_{family}_base"
            n = len(matrix)
            variants = [("base", list(range(n)), 0),
                        ("permuted", permutation(n, seed), 0),
                        ("small", list(range(n)), EXPONENTS[0]),
                        ("large", list(range(n)), EXPONENTS[1])]
            for variant, order, exponent in variants:
                name = f"s{seed:04d}_{family}_{variant}"
                a, target, zeros = transformed(matrix, x, structural, order, exponent)
                shifts = None
                if family not in ("zero_diagonal", "indefinite"):
                    scale = dyadic_power(exponent)
                    shifts = (-float(scale), float((max(matrix[i][i] for i in range(n)) + 1) * scale))
                item = oracle.write_fixture(folder, name, a, target, zeros, shifts)
                item.update(seed=seed, family=family, variant=variant, base=base,
                            permutation=order, scale_exponent=exponent,
                            forward_tolerance=tolerance, construction=details)
                paths = sorted(folder.glob(f"{name}.*.txt"))
                item["input_sha256"] = {path.name: digest(path) for path in paths}
                cases.append(item)
    expected = {case["name"] for case in cases}
    found = {path.name[:-6] for path in folder.glob("*.A.txt")}
    oracle.require(found == expected, "generation folder contains unrelated/stale matrices")
    manifest = {"schema": SCHEMA, "seed_count": seeds, "case_count": len(cases),
                "generator": "specified 32-bit linear recurrence; seeds start at zero",
                "classification": "strict exact/dyadic input lane; no ambiguous boundary signs",
                "families": dict(Counter(item["family"] for item in cases)),
                "expected": dict(Counter(item["expected"] for item in cases)),
                "cases": cases}
    write_json(folder / "stress-manifest.json", manifest)
    return manifest


def load_manifest(folder):
    manifest = json.loads((folder / "stress-manifest.json").read_text(encoding="utf-8"))
    oracle.require(manifest.get("schema") == SCHEMA, "unsupported stress manifest")
    cases = manifest["cases"]
    oracle.require(len(cases) == manifest["case_count"] and bool(cases), "invalid case count")
    names = [case["name"] for case in cases]
    oracle.require(all(re.fullmatch(r"[a-zA-Z0-9_]+", name) for name in names),
                   "invalid case name")
    oracle.require(len(set(names)) == len(names), "duplicate case name")
    found = {path.name[:-6] for path in folder.glob("*.A.txt")}
    oracle.require(found == set(names), "manifest does not cover every input matrix")
    for case in cases:
        for filename, expected in case["input_sha256"].items():
            oracle.require(Path(filename).name == filename, "invalid manifest input path")
            oracle.require(digest(folder / filename) == expected,
                           f"{case['name']}: input digest mismatch for {filename}")
    return manifest


def close_vector(actual, expected, tolerance, label):
    oracle.require(len(actual) == len(expected), f"{label}: length mismatch")
    oracle.require(all(math.isfinite(x) for x in actual), f"{label}: nonfinite values")
    scale = max((abs(x) for x in expected), default=0.0)
    error = max((abs(a - b) for a, b in zip(actual, expected)), default=0.0)
    oracle.require(error <= tolerance * max(scale, sys.float_info.min),
                   f"{label}: relative solution error {error / max(scale, sys.float_info.min):.8g}")
    return error / max(scale, sys.float_info.min)


def close_diagonal(actual, expected, tolerance, label):
    """Compare each positive pivot relative to itself, without a larger-entry floor."""
    oracle.require(len(actual) == len(expected), f"{label}: length mismatch")
    worst = 0.0
    for index, (a, e) in enumerate(zip(actual, expected)):
        oracle.require(math.isfinite(a) and math.isfinite(e) and a > 0 and e > 0,
                       f"{label}: pivot {index} must be finite and positive")
        relative = abs(a - e) / e
        oracle.require(math.isfinite(relative) and relative <= tolerance,
                       f"{label}: pivot {index} relative scale error {relative:.8g}")
        worst = max(worst, relative)
    return worst


def verify_all(folder, output):
    manifest = load_manifest(folder)
    cases = manifest["cases"]
    results = {}
    worst_forward = worst_backward = 0.0
    for case in cases:
        name = case["name"]
        result = oracle.verify_fixture(folder, output, name, require_shifts=True)
        oracle.require(("refused" in result) == (case["expected"] == "refused"),
                       f"{name}: outcome differs from the exact input classification")
        if case["expected"] == "success":
            actual = oracle.read_vector(output / f"{name}.x.out.txt", case["n"])
            expected = oracle.read_vector(folder / f"{name}.x.txt", case["n"])
            relative = close_vector(actual, expected, case["forward_tolerance"], name)
            result["forward_error"] = relative
            worst_forward = max(worst_forward, relative)
            worst_backward = max(worst_backward, result["backward_error"])
        results[name] = result

    relations = 0
    for case in cases:
        if case["variant"] == "base":
            continue
        name, base, n = case["name"], case["base"], case["n"]
        base_n, base_entries = oracle.read_triplets(folder / f"{base}.A.txt")
        actual_n, actual_entries = oracle.read_triplets(folder / f"{name}.A.txt")
        order = case["permutation"]
        oracle.require(base_n == actual_n == n and sorted(order) == list(range(n)),
                       f"{name}: invalid transformation dimensions")
        inverse = {old: new for new, old in enumerate(order)}
        scale = math.ldexp(1.0, case["scale_exponent"])
        expected_entries = {tuple(sorted((inverse[i], inverse[j]), reverse=True)): value * scale
                            for (i, j), value in base_entries.items()}
        oracle.require(actual_entries == expected_entries, f"{name}: matrix transformation differs")
        base_b = oracle.read_vector(folder / f"{base}.b.txt", n)
        actual_b = oracle.read_vector(folder / f"{name}.b.txt", n)
        oracle.require(actual_b == [base_b[i] * scale for i in order],
                       f"{name}: RHS transformation differs")
        oracle.require(("refused" in results[name]) == ("refused" in results[base]),
                       f"{name}: success/refusal is not transformation-invariant")
        if case["expected"] == "success":
            base_x = oracle.read_vector(output / f"{base}.x.out.txt", n)
            actual_x = oracle.read_vector(output / f"{name}.x.out.txt", n)
            close_vector(actual_x, [base_x[i] for i in order],
                         2 * case["forward_tolerance"], f"{name}: solution transformation")
            if case["scale_exponent"]:
                base_d = oracle.read_vector(output / f"{base}.D.out.txt", n)
                actual_d = oracle.read_vector(output / f"{name}.D.out.txt", n)
                close_diagonal([value / scale for value in actual_d], base_d, 2e-12,
                             f"{name}: diagonal scale")
                _, base_l = oracle.read_triplets(output / f"{base}.L.out.txt", strict=True)
                _, actual_l = oracle.read_triplets(output / f"{name}.L.out.txt", strict=True)
                oracle.require(set(base_l) == set(actual_l), f"{name}: scaled L pattern differs")
                close_vector([actual_l[key] for key in sorted(base_l)],
                             [base_l[key] for key in sorted(base_l)], 2e-12,
                             f"{name}: scaled L values")
        relations += 1
    return {"schema": SCHEMA, "passed": True, "case_count": len(cases),
            "seed_count": manifest["seed_count"], "families": manifest["families"],
            "outcomes": dict(Counter("refused" if "refused" in r else "success"
                                     for r in results.values())),
            "metamorphic_relations": relations,
            "max_forward_error": worst_forward, "max_backward_error": worst_backward}


def oracle_outputs(fixtures, output, cases):
    for case in cases:
        name = case["name"]
        oracle.write_oracle_outputs(fixtures, output, name)
        meta = oracle.metadata(fixtures / f"{name}.meta.txt")
        for key, suffix in (
            ("shift that must still factorise with every pivot positive", ".shift-kept"),
            ("shift that must be refused", ".shift-refused"),
        ):
            if key in meta:
                oracle.write_oracle_outputs(fixtures, output, name, float(meta[key]), suffix)


def selftest(folder):
    fixtures, output = folder / "fixtures", folder / "oracle"
    manifest = generate(fixtures, seeds=1)
    oracle_outputs(fixtures, output, manifest["cases"])
    baseline = verify_all(fixtures, output)
    evidence = []
    # Exercise the complete checker; preserve each corrupted artifact separately.
    controls = [
        ("permuted_solution", output / "s0000_dominant_permuted.x.out.txt", "999\n"),
        ("scaled_solution", output / "s0000_dominant_large.x.out.txt", "-999\n"),
        ("scaled_diagonal", output / "s0000_dominant_small.D.out.txt", "1\n"),
        ("input_transformation", fixtures / "s0000_dominant_permuted.b.txt", "999\n"),
    ]
    for label, target, first in controls:
        original = target.read_text(encoding="utf-8")
        corrupted = first + "\n".join(original.splitlines()[1:]) + "\n"
        destination = folder / "controls" / label
        destination.mkdir(parents=True, exist_ok=True)
        (destination / target.name).write_text(corrupted, encoding="utf-8")
        try:
            target.write_text(corrupted, encoding="utf-8")
            try:
                verify_all(fixtures, output)
            except (oracle.VerificationError, ValueError, OSError, OverflowError) as exc:
                evidence.append({"control": label, "rejected": True, "reason": str(exc)})
            else:
                raise AssertionError(f"planted defect escaped: {label}")
        finally:
            target.write_text(original, encoding="utf-8")
    # Missing outputs must not reduce the scope silently.
    target = output / "s0000_known_ldl_base.x.out.txt"
    original = target.read_bytes()
    try:
        target.unlink()
        try:
            verify_all(fixtures, output)
        except (oracle.VerificationError, ValueError, OSError, OverflowError) as exc:
            evidence.append({"control": "missing_output", "rejected": True, "reason": str(exc)})
        else:
            raise AssertionError("missing artifact escaped")
    finally:
        target.write_bytes(original)
    expected_pivots = [1.0, math.ldexp(1.0, -100)]
    corrupted_pivots = [1.0, math.ldexp(1.0, -99)]
    # The former whole-vector check must accept this corruption, exposing its blind spot.
    global_error = close_vector(corrupted_pivots, expected_pivots, 2e-12,
                                "tiny-pivot control: former global comparison")
    try:
        close_diagonal(corrupted_pivots, expected_pivots, 2e-12,
                       "tiny-pivot control: componentwise comparison")
    except oracle.VerificationError as exc:
        evidence.append({"control": "tiny_pivot_scaling", "rejected": True,
                         "former_global_comparison_accepted": True,
                         "former_global_relative_error": global_error,
                         "expected_pivots": expected_pivots, "actual_pivots": corrupted_pivots,
                         "reason": str(exc)})
    else:
        raise AssertionError("a factor-of-two tiny-pivot error escaped")
    repeated = generate(folder / "regenerated", seeds=1)
    original_names = sorted(path.name for path in fixtures.iterdir() if path.is_file())
    repeated_names = sorted(path.name for path in (folder / "regenerated").iterdir() if path.is_file())
    oracle.require(original_names == repeated_names, "generator file roster is not reproducible")
    oracle.require(all(digest(fixtures / name) == digest(folder / "regenerated" / name)
                       for name in original_names), "generator contents are not reproducible")
    oracle.require(manifest == repeated, "generator manifest is not reproducible")
    verify_all(fixtures, output)
    result = {"baseline": baseline, "planted_controls": evidence,
              "reproducible_files": len(original_names)}
    write_json(folder / "evidence.json", result)
    return {"selftest_passed": True, "oracle_cases": baseline["case_count"],
            "planted_controls_rejected": len(evidence), "reproducible_files": len(original_names)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    gen = commands.add_parser("generate")
    gen.add_argument("folder", type=Path)
    gen.add_argument("--seeds", type=int, default=DEFAULT_SEEDS)
    check = commands.add_parser("verify")
    check.add_argument("fixtures", type=Path)
    check.add_argument("output", type=Path)
    check.add_argument("--json", type=Path)
    own = commands.add_parser("selftest")
    own.add_argument("folder", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "generate":
            result = generate(args.folder, args.seeds)
            result = {key: result[key] for key in ("case_count", "seed_count", "families", "expected")}
        elif args.command == "verify":
            result = verify_all(args.fixtures, args.output)
            if args.json:
                write_json(args.json, result)
        else:
            result = selftest(args.folder)
        print(json.dumps(result, sort_keys=True))
        return 0
    except (oracle.VerificationError, ValueError, OSError, OverflowError, AssertionError) as exc:
        print(json.dumps({"passed": False, "error": str(exc)}, sort_keys=True))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
