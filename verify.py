#!/usr/bin/env python3
"""Verify the file hashes and run both circuit-distance witness checkers."""

import hashlib
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
FILES = {
    "gen_bb_memory.py",
    "verify_stim.py",
    "verify_transcription.py",
    "witness_Ztype_16faults.json",
    "witness_Xtype_16faults.json",
    "witness_A23swap_Ztype_16faults.json",
    "witness_A23swap_Xtype_16faults.json",
}
FINGERPRINT = "8002fb86c839b990f501ceb0237149c7a94ecc77b55b9f76da4cd052d8689769"


def check_files():
    seen = set()
    for line in (ROOT / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        expected, name = line.split()
        if name not in FILES or name in seen:
            raise ValueError(f"Unexpected or duplicate checksum entry: {name}")
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"Checksum mismatch: {name}")
        seen.add(name)
    if seen != FILES:
        raise ValueError("SHA256SUMS does not list all seven verification files")
    print("File integrity: all seven SHA-256 hashes match.", flush=True)


def run_checker(name):
    print(f"\nRunning {name}", flush=True)
    # Isolated mode ignores PYTHONOPTIMIZE so the original assertions stay active.
    result = subprocess.run(
        [sys.executable, "-I", "-u", str(ROOT / name)],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
    )
    print(result.stdout, end="", flush=True)
    # The transcription checker reports failure in its output, not its exit code.
    if result.returncode != 0 or result.stdout.rstrip().splitlines()[-1:] != ["ALL CHECKS PASS"]:
        raise ValueError(f"{name} did not pass")
    if name == "verify_stim.py" and f"sha256 {FINGERPRINT} " not in result.stdout:
        raise ValueError("The two-gross detector-signature fingerprint does not match")


def main():
    try:
        check_files()
        run_checker("verify_stim.py")
        run_checker("verify_transcription.py")
    except (OSError, ValueError) as error:
        print(f"\nVERIFICATION FAILED: {error}", file=sys.stderr)
        return 1
    print("\nVERIFIED: all four 16-fault witnesses pass (circuit distance <= 16).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
