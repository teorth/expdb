#!/usr/bin/env python3
"""Check (or fix) trailing whitespace and final newlines on tracked text files.

Enforces the on-save formatting rules declared in .vscode/settings.json:
  - files.trimTrailingWhitespace : no spaces/tabs before a line terminator
  - files.insertFinalNewline     : file ends with a newline
  - files.trimFinalNewlines      : no blank lines at end of file

Line endings (CRLF vs LF) are preserved per file: VS Code only trims
horizontal whitespace on save, it does not rewrite an existing file's EOL
style. Run with --fix to rewrite offending files, otherwise the script only
reports them and exits non-zero if any are found (used by CI).
"""

import re
import subprocess
import sys
from pathlib import Path

# Text file types that the editor settings apply to.
INCLUDE_EXT = {
    ".py", ".tex", ".lean", ".md", ".yml", ".yaml", ".json", ".sty",
    ".sh", ".scss", ".css", ".cfg", ".bib", ".toml", ".html", ".bat",
}

# Generated / vendored files that are maintained by tooling, not by hand.
EXCLUDE = {
    "lake-manifest.json",
    "home_page/Gemfile.lock",
}

# Trailing spaces/tabs before a line terminator (CRLF or LF), per line.
_TRAILING = re.compile(rb"[ \t]+(\r?\n)")


def normalize(data: bytes) -> bytes:
    if b"\x00" in data:  # binary safety net
        return data
    eol = b"\r\n" if b"\r\n" in data else b"\n"
    data = _TRAILING.sub(rb"\1", data)          # trim trailing whitespace
    body = data.rstrip(b" \t\r\n")               # drop trailing blank lines/ws
    return body + eol if body else b""           # exactly one final newline


def repo_files():
    out = subprocess.run(
        ["git", "ls-files"], capture_output=True, text=True, check=True
    ).stdout.splitlines()
    for rel in out:
        if rel in EXCLUDE:
            continue
        if Path(rel).suffix.lower() in INCLUDE_EXT:
            yield rel


def main() -> int:
    fix = "--fix" in sys.argv[1:]
    root = Path(
        subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    )
    offenders = []
    for rel in repo_files():
        path = root / rel
        data = path.read_bytes()
        fixed = normalize(data)
        if fixed != data:
            offenders.append(rel)
            if fix:
                path.write_bytes(fixed)

    if not offenders:
        return 0
    verb = "Fixed" if fix else "Would fix"
    print(f"{verb} {len(offenders)} file(s):")
    for rel in offenders:
        print(f"  {rel}")
    if not fix:
        print("\nRun `python scripts/check_whitespace.py --fix` to normalize.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
