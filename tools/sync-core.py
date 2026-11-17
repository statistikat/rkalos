#!/usr/bin/env python3
"""Export a committed, allowlisted Kalos snapshot, or verify one without Git.

Lock schema 2 records the core's three-section manifest (tools/source-files.json:
core, fixtures, evidence) as one snapshot entry per section. The historical
validation/fixtures/ tree and the newer validation/evidence/ tree both map to
the SAME inst/tinytest/fixtures/ layout (each section's files keep their
post-prefix paths), so the R test paths never moved; the entries keep their
section provenance (and source prefix) distinct in the lock. Lock schema 1
(single fixtures section) stays readable by --verify; schema 2 is rejected by
older tools rather than silently misverified.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import subprocess
import tarfile


def sha(data):
    return hashlib.sha256(data).hexdigest()


def safe_path(value):
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or str(path) != value:
        raise ValueError(f"Unsafe manifest path: {value}")
    return value


def files_under(root):
    result = {}
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Unexpected snapshot symlink: {path}")
        if path.is_file():
            result[path.relative_to(root).as_posix()] = sha(path.read_bytes())
    return result


def merge_sectioned(snapshots):
    """Group section entries by destination directory, checking for conflicts."""
    grouped = {}
    for section in snapshots:
        directory = section["directory"]
        files = grouped.setdefault(directory, {})
        for path, checksum in section["files"].items():
            if path in files and files[path] != checksum:
                raise ValueError(
                    f"Snapshot sections disagree on {directory}/{path}")
            files[path] = checksum
    return grouped


def verify(destination):
    lock_path = destination / "core.lock.json"
    if not lock_path.exists():
        lock_path = destination / "inst/core.lock.json"
    lock = json.loads(lock_path.read_text())
    if lock["schema_version"] not in (1, 2):
        raise ValueError("Unsupported source lock version")
    for section in lock["snapshots"]:
        for path in section["files"]:
            safe_path(path)
    for directory, expected in merge_sectioned(lock["snapshots"]).items():
        root = destination / safe_path(directory)
        actual = files_under(root)
        if actual != expected:
            missing = sorted(expected.keys() - actual.keys())
            extra = sorted(actual.keys() - expected.keys())
            changed = sorted(p for p in actual.keys() & expected.keys() if actual[p] != expected[p])
            raise ValueError(f"Snapshot mismatch in {root}: missing={missing}, unexpected={extra}, changed={changed}")
    if lock["schema_version"] == 1:
        detail = ", ".join(
            f"{len(s['files'])} files in {s['directory']}" for s in lock["snapshots"])
    else:
        detail = ", ".join(
            f"{len(s['files'])} {s['section']} files in {s['directory']}"
            for s in lock["snapshots"])
    print(f"Verified {lock['commit']} (lock schema {lock['schema_version']}): {detail}")


def export(core, revision, kind, destination):
    def git(*args):
        return subprocess.check_output(["git", "-C", str(core), *args])
    commit = git("rev-parse", "--verify", revision + "^{commit}").decode().strip()
    inventory = json.loads(git("show", commit + ":tools/source-files.json"))
    fixture_sections = [("fixtures", "validation/fixtures/")]
    if "evidence" in inventory:
        fixture_sections.append(("evidence", "validation/evidence/"))
    selected = inventory["core"] + [name for section, _ in fixture_sections
        for name in inventory.get(section, [])] if kind == "r" else list(inventory["core"])
    for name in selected:
        safe_path(name)
    if len(set(selected)) != len(selected):
        raise ValueError("Duplicate allowlist entries")
    archive = git("archive", "--format=tar", commit, "--", *selected)
    with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
        contents = {m.name: tar.extractfile(m).read() for m in tar.getmembers() if m.isfile()}
    if set(contents) != set(selected):
        raise ValueError("The committed export does not match its allowlist")
    core_files = {name: contents[name] for name in inventory["core"]}
    core_files["include/kalos/revision.hpp"] = (
        f'// Generated from core.lock.json; do not edit.\n#define KALOS_CORE_REVISION "{commit}"\n'
    ).encode()
    # Section entries: (directory, section, source_prefix, files). The fixture
    # and evidence sections share one destination tree with identical layouts.
    snapshots = [("src/kalos" if kind == "r" else "vendor/kalos", "core", "", core_files)]
    if kind == "r":
        seen = {}
        for section, prefix in fixture_sections:
            files = {name.removeprefix(prefix): contents[name] for name in inventory[section]}
            for path, data in files.items():
                if path in seen and seen[path] != data:
                    raise ValueError(f"fixtures and evidence disagree on {path}")
                seen[path] = data
            snapshots.append(("inst/tinytest/fixtures", section, prefix, files))
    previous_path = destination / "core.lock.json"
    previous = json.loads(previous_path.read_text()) if previous_path.exists() else None
    if previous:
        verify(destination)  # Refuse to discard local or unexpected snapshot work.
    lock = {"schema_version": 2, "repository": "https://github.com/statistikat/kalos",
            "commit": commit, "version": inventory["version"],
            "source_manifest_sha256": sha(git("show", commit + ":tools/source-files.json")),
            "generated_files": ["include/kalos/revision.hpp"], "snapshots": []}
    previous_files = merge_sectioned(previous["snapshots"]) if previous else {}
    written = {}
    for directory, section, prefix, files in snapshots:
        root = destination / directory
        if root.exists() and not previous:
            # The extracted historical fixture directory already contains exact originals.
            actual = files_under(root)
            if any(path not in files or checksum != sha(files[path]) for path, checksum in actual.items()):
                raise ValueError(f"Refusing to overwrite unrecognized contents in {root}")
        if previous:
            old = previous_files.get(directory, {})
            keep = set().union(*[set(f) for d2, _, _, f in snapshots if d2 == directory])
            for removed in old.keys() - keep:
                (root / safe_path(removed)).unlink()
        for name, data in files.items():
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        written[directory] = written.get(directory, 0) + 1
        entry = {"directory": directory, "section": section, "files":
            {name: sha(data) for name, data in sorted(files.items())}}
        if prefix:
            entry["source_prefix"] = prefix
        lock["snapshots"].append(entry)
    previous_path.write_text(json.dumps(lock, indent=2) + "\n")
    if kind == "r":
        (destination / "inst/core.lock.json").write_bytes(previous_path.read_bytes())
    verify(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--core", type=Path)
    parser.add_argument("--revision", default="HEAD")
    parser.add_argument("--kind", choices=["r", "python"])
    parser.add_argument("--destination", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    if args.verify:
        verify(args.destination)
    else:
        if not args.core or not args.kind:
            parser.error("Export requires --core and --kind")
        export(args.core, args.revision, args.kind, args.destination)


if __name__ == "__main__":
    main()
