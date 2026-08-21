#!/usr/bin/env python3
"""Fail the build when a pinned package resolved to something else.

pip exits 0 after installing a version that satisfies no pin we care about, so
the pins in base.txt are only a request until something checks the result. The
CUDA line matters here: a torch built for a newer driver than the host leaves
torch.cuda unavailable at runtime while the build and nvidia-smi both look fine.

torch.cuda itself is not checked, because the builder has no GPU.
"""

import re
import sys
from importlib import metadata
from pathlib import Path

REQUIREMENTS = Path(__file__).with_name("base.txt")
PIN = re.compile(r"^(?P<name>[A-Za-z0-9][A-Za-z0-9._-]*)==(?P<version>[^\s;#]+)")


def pinned(path):
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "-")):
            continue
        match = PIN.match(line)
        if match:
            yield match.group("name"), match.group("version")


def main():
    pins = dict(pinned(REQUIREMENTS))
    if not pins:
        sys.exit(f"No pins parsed from {REQUIREMENTS}; the check would pass vacuously.")

    failures = []
    for name, expected in sorted(pins.items()):
        try:
            found = metadata.version(name)
        except metadata.PackageNotFoundError:
            failures.append(f"  {name}: pinned {expected}, not installed")
            continue
        # A local version such as 2.11.0+cu128 must match exactly, including the
        # +cu suffix, which is the whole point of pinning it.
        if found != expected:
            failures.append(f"  {name}: pinned {expected}, installed {found}")

    print(f"verify_base: checked {len(pins)} pinned package(s)")
    if failures:
        sys.exit("verify_base: pinned versions did not resolve as requested:\n"
                 + "\n".join(failures))
    print("verify_base: all pinned versions match")


if __name__ == "__main__":
    main()
