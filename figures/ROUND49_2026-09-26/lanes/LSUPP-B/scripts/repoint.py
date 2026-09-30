"""Copy a builder script from an earlier round into this lane and apply exact-string edits, each recorded (old -> new) in the sheet's
scripts/R49_EDITS.md. Every old string must occur exactly once (the R40 convention). Never touches the source file."""
import hashlib, os, shutil, sys


def apply_edits(src, dst, edits, notes_md):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    txt = open(src, encoding="utf-8").read()
    sha = hashlib.sha256(txt.encode()).hexdigest()
    lines = [f"\n## {os.path.basename(dst)}  (from {src}, sha256 {sha[:16]})"]
    for old, new in edits:
        n = txt.count(old)
        assert n == 1, f"{dst}: expected exactly one occurrence of {old[:80]!r}, found {n}"
        txt = txt.replace(old, new)
        lines.append(f"- `{old.strip()[:110]}` -> `{new.strip()[:110]}`")
    open(dst, "w", encoding="utf-8").write(txt)
    with open(notes_md, "a", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return dst
