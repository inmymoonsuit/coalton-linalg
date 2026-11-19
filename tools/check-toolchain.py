"""Check an installed offline toolchain against reviewed source-content pins.

This tool never installs, downloads, compiles, or modifies dependencies. A source
check proves file contents; --metadata additionally checks the ASDF systems that
the separate Lisp metadata probe actually loaded.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tempfile


class Mismatch(ValueError):
    """An installed input or loaded system differs from the lock."""


EXCLUDED_DIRECTORIES = {".git", "__pycache__", ".cache"}
EXCLUDED_SUFFIXES = {".fasl", ".pyc", ".log"}
EXCLUDED_FILES = {"system-index.txt"}


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def manifest(root: Path, selected: list[str] | None = None) -> dict:
    """Hash all source inputs, including their relative names and sizes."""
    if not root.is_dir():
        raise Mismatch(f"Missing installed root: {root}")
    paths = [root / name for name in selected] if selected is not None else root.rglob("*")
    result = {}
    for path in paths:
        relative = path.relative_to(root)
        if any(part in EXCLUDED_DIRECTORIES for part in relative.parts):
            continue
        if selected is None and (path.name in EXCLUDED_FILES or path.suffix.lower() in EXCLUDED_SUFFIXES):
            continue
        if not path.is_file():
            if selected is not None:
                raise Mismatch(f"Missing pinned file: {path}")
            continue
        if not path.resolve().is_relative_to(root.resolve()):
            raise Mismatch(f"Pinned input escapes its root: {path}")
        result[relative.as_posix()] = {"sha256": digest(path), "bytes": path.stat().st_size}
    return dict(sorted(result.items()))


def resolve_roots(lock: dict, supplied: dict[str, Path]) -> dict[str, Path]:
    result = {}
    for name, component in lock["components"].items():
        base = supplied[component["base"]].resolve()
        root = (base / component["path"]).resolve()
        if not root.is_relative_to(base):
            raise Mismatch(f"Component escapes its supplied root: {name}")
        result[name] = root
    return result


def read_metadata(path: Path) -> dict:
    result = {"systems": {}, "environment": {}}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        fields = line.split("\t")
        if fields[0] in {"SBCL", "ASDF"} and len(fields) == 2:
            if fields[0].lower() in result:
                raise Mismatch(f"Duplicate metadata field: {fields[0]}")
            result[fields[0].lower()] = fields[1]
        elif fields[0] == "ENV" and len(fields) == 3:
            result["environment"][fields[1]] = fields[2]
        elif fields[0] == "SYSTEM" and len(fields) == 4:
            if fields[1] in result["systems"]:
                raise Mismatch(f"Duplicate loaded system: {fields[1]}")
            result["systems"][fields[1]] = {"version": fields[2], "source": fields[3]}
        else:
            raise Mismatch(f"Unrecognized metadata row: {line!r}")
    if "sbcl" not in result or "asdf" not in result or not result["systems"]:
        raise Mismatch("Metadata is incomplete.")
    return result


def normalize_metadata(metadata: dict, roots: dict[str, Path], library_root: Path) -> dict:
    """Replace machine paths with pinned component names; refuse unknown sources."""
    normalized = {key: metadata[key] for key in ("sbcl", "asdf", "environment")}
    normalized["systems"] = {}
    for name, system in sorted(metadata["systems"].items()):
        raw_source = system["source"]
        if name == "coalton-linalg" or name.startswith("coalton-linalg/"):
            source = Path(raw_source).resolve()
            if not source.is_relative_to(library_root.resolve()):
                raise Mismatch(f"Project system resolved outside this checkout: {name}: {source}")
            continue
        if raw_source == "<built-in>":
            if name not in {"asdf", "asdf-package-system", "uiop"} and not name.startswith("sb-"):
                raise Mismatch(f"Unexpected built-in system: {name}")
            normalized["systems"][name] = {"version": system["version"], "builtin": True}
            continue
        source = Path(raw_source).resolve()
        matches = [(len(root.parts), component, root) for component, root in roots.items()
                   if source.is_relative_to(root)]
        if not matches:
            raise Mismatch(f"Loaded system has an unpinned source: {name}: {source}")
        _, component, root = max(matches)
        relative = source.relative_to(root).as_posix()
        if not source.is_file():
            raise Mismatch(f"Loaded system source is absent: {name}: {source}")
        normalized["systems"][name] = {
            "version": system["version"], "component": component, "source": relative
        }
    return normalized


def verify(lock: dict, roots: dict[str, Path], metadata: dict | None = None,
           library_root: Path | None = None) -> tuple[int, int]:
    if lock.get("format") != 1:
        raise Mismatch("Unsupported lock format.")
    file_count = 0
    for name, component in lock["components"].items():
        actual = manifest(roots[name], component.get("selected"))
        expected = component["files"]
        if actual != expected:
            missing = sorted(expected.keys() - actual.keys())
            extra = sorted(actual.keys() - expected.keys())
            changed = sorted(key for key in actual.keys() & expected.keys()
                             if actual[key] != expected[key])
            raise Mismatch(f"{name}: missing={missing[:5]}, extra={extra[:5]}, changed={changed[:5]}")
        file_count += len(actual)
    system_count = 0
    if metadata is not None:
        if not lock.get("runtime"):
            raise Mismatch("The lock has no recorded ASDF runtime resolution.")
        normalized = normalize_metadata(metadata, roots, library_root)
        if normalized != lock["runtime"]:
            raise Mismatch("Loaded ASDF systems, versions, sources, or compiler environment differ from the lock.")
        for name, system in normalized["systems"].items():
            if not system.get("builtin"):
                pinned_files = lock["components"][system["component"]]["files"]
                if system["source"] not in pinned_files:
                    raise Mismatch(f"Loaded source was excluded from the content pin: {name}")
        system_count = len(normalized["systems"])
    return file_count, system_count


def selftest(folder: Path) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="toolchain-", dir=folder) as scratch:
        root = Path(scratch)
        source = root / "dependency.asd"
        source.write_text("original\n", encoding="utf-8")
        component = {"files": manifest(root)}
        lock = {"format": 1, "components": {"dependency": component}, "runtime": {
            "sbcl": "test", "asdf": "test", "environment": {"COALTON_ENV": "release"},
            "systems": {"dependency": {"version": "1", "component": "dependency",
                                       "source": "dependency.asd"}}}}
        roots = {"dependency": root}
        metadata = {"sbcl": "test", "asdf": "test", "environment": {"COALTON_ENV": "release"},
                    "systems": {"dependency": {"version": "1", "source": str(source)}}}
        verify(lock, roots, metadata, root)
        rejected = 0

        def expect_mismatch() -> None:
            nonlocal rejected
            try:
                verify(lock, roots, metadata, root)
            except Mismatch:
                rejected += 1
            else:
                raise AssertionError("Planted toolchain change was accepted.")

        source.write_text("modified\n", encoding="utf-8")
        expect_mismatch()
        source.unlink()
        expect_mismatch()
        source.write_text("original\n", encoding="utf-8")
        extra = root / "extra.lisp"
        extra.write_text("extra", encoding="utf-8")
        expect_mismatch()
        extra.unlink()
        metadata["systems"]["dependency"]["version"] = "2"
        expect_mismatch()
        metadata["systems"]["dependency"]["version"] = "1"
        metadata["systems"]["dependency"]["source"] = str(root.parent / "unpinned.asd")
        expect_mismatch()
        metadata["systems"]["dependency"]["source"] = str(source)
        metadata["environment"]["COALTON_ENV"] = "development"
        expect_mismatch()
        metadata["environment"]["COALTON_ENV"] = "release"
        verify(lock, roots, metadata, root)
        print(f"TOOLCHAIN-CONTROLS-PASS {rejected} planted mismatches rejected")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, default=Path(__file__).resolve().parents[1] / "toolchain.lock.json")
    parser.add_argument("--sbcl-root", type=Path)
    parser.add_argument("--coalton-root", type=Path)
    parser.add_argument("--quicklisp-root", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--selftest", type=Path, metavar="SCRATCH_FOLDER")
    args = parser.parse_args()
    try:
        if args.selftest is not None:
            selftest(args.selftest)
            return 0
        supplied = {name: getattr(args, name + "_root") for name in ("sbcl", "coalton", "quicklisp")}
        if any(value is None for value in supplied.values()):
            parser.error("--sbcl-root, --coalton-root and --quicklisp-root are required")
        lock = json.loads(args.lock.read_text(encoding="utf-8"))
        roots = resolve_roots(lock, supplied)
        metadata = read_metadata(args.metadata) if args.metadata else None
        files, systems = verify(lock, roots, metadata, Path(__file__).resolve().parents[1])
        print(f"TOOLCHAIN-CONTENT-PASS {len(roots)} components, {files} files")
        if metadata is None:
            print("ASDF resolution not checked; provide fresh --metadata from the Lisp probe.")
        else:
            print(f"TOOLCHAIN-RESOLUTION-PASS {systems} loaded dependency systems")
        return 0
    except (Mismatch, OSError, ValueError, KeyError) as error:
        print(f"TOOLCHAIN-FAIL {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
