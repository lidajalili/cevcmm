#!/usr/bin/env python3
"""
fix_citations_v2.py

Replaces the citation line-by-line (not as one big block) so it's robust to
small whitespace / wrap differences between files. The Applications line
grows from 1 line to 2 lines with an embedded '\\n#\\' ' so the added
"Statistical Association." tail sits on its own roxygen-formatted line.

Run once from the package root:
    python3 fix_citations_v2.py

Then in RStudio:
    devtools::document()
    devtools::build_readme()
    devtools::build_vignettes()
    pkgdown::build_site()
    devtools::build_manual()
"""

import os

SKIP_DIRS = {
    'docs', 'doc', 'man',
    '.Rproj.user', '.git', '.github',
    'renv', '.Rcheck',
}
SKIP_FILES = {'fix_citations.py', 'fix_citations_v2.py'}
PROCESS_EXTS = ('.R', '.Rmd', '.md', '.yml', '.yaml')
EXTRA_FILES = {'DESCRIPTION'}

# ----------------------------------------------------------------------------
# Roxygen (.R files): 3 independent single-line replacements.
# Line C's replacement contains "\n#' " so it grows to 2 lines cleanly.
# ----------------------------------------------------------------------------

R_LINE_REPLACEMENTS = [
    # Author line
    (
        "Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient",
        "Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient"
    ),
    # Title line (drop the hyphen, "Effects" -> "Effect")
    (
        "Varying Coefficient Mixed-Effects Models: Methodology, Theory, and",
        "Varying Coefficient Mixed Effect Models: Methodology, Theory, and"
    ),
    # Venue line -> arXiv reference (1 line becomes 2)
    (
        "Applications. Journal of the American Statistical Association (under review).",
        "Applications. arXiv:2511.12732; under review at Journal of the American\n#' Statistical Association."
    ),
]

# ----------------------------------------------------------------------------
# Non-R text files (Rmd, md, yml, DESCRIPTION): also cover the fully-inline
# long-form citation that may appear in a single line (markdown italicized).
# ----------------------------------------------------------------------------

NONR_REPLACEMENTS = [
    # Long inline form with markdown italics (Rmd, md)
    (
        "Lin, L.-H. and Jalili, L. (2026). *Scalable and Communication-Efficient Varying Coefficient Mixed-Effects Models: Methodology, Theory, and Applications.* Journal of the American Statistical Association (under review).",
        "Jalili, L. and Lin, L.-H. (2025). *Scalable and Communication-Efficient Varying Coefficient Mixed Effect Models: Methodology, Theory, and Applications.* arXiv:2511.12732; under review at *Journal of the American Statistical Association*."
    ),
    # Long inline form no italics (plain md, yml)
    (
        "Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient Varying Coefficient Mixed-Effects Models: Methodology, Theory, and Applications. Journal of the American Statistical Association (under review).",
        "Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient Varying Coefficient Mixed Effect Models: Methodology, Theory, and Applications. arXiv:2511.12732; under review at Journal of the American Statistical Association."
    ),
]

# ----------------------------------------------------------------------------
# Inline mentions of the paper (all file types). Ordered most-specific-first.
# ----------------------------------------------------------------------------

INLINE_REPLACEMENTS = [
    # Multi-line inline citations wrapped across two roxygen lines
    ("Lin and Jalili,\n#' 2026)", "Jalili and Lin,\n#' 2025)"),
    ("Lin & Jalili,\n#' 2026)",   "Jalili and Lin,\n#' 2025)"),

    # Single-line, no closing paren (matches "(2026)", "(2026, Sec. 3)", etc.)
    ("Lin and Jalili (2026", "Jalili and Lin (2025"),
    ("Lin & Jalili (2026",   "Jalili and Lin (2025"),

    # Single-line, comma-year form
    ("Lin and Jalili, 2026", "Jalili and Lin, 2025"),
    ("Lin & Jalili, 2026",   "Jalili and Lin, 2025"),
]


def process_file(path):
    """Return True if the file was changed, False otherwise."""
    try:
        with open(path, 'r', encoding='utf-8') as f:
            content = f.read()
    except (UnicodeDecodeError, PermissionError):
        return False
    original = content

    is_r = path.endswith('.R')

    if is_r:
        # Line-by-line replacements for the 3-line roxygen bibliographic block.
        # Third line grows to 2 lines (added "\n#' Statistical Association.")
        for old, new in R_LINE_REPLACEMENTS:
            content = content.replace(old, new)
    else:
        # Long-form single-line inline citations (README, NEWS, Rmd, yml)
        for old, new in NONR_REPLACEMENTS:
            content = content.replace(old, new)

    # Inline mentions apply to every file type
    for old, new in INLINE_REPLACEMENTS:
        content = content.replace(old, new)

    if content != original:
        with open(path, 'w', encoding='utf-8') as f:
            f.write(content)
        return True
    return False


def should_process(fname):
    if fname in SKIP_FILES:
        return False
    if fname in EXTRA_FILES:
        return True
    return fname.endswith(PROCESS_EXTS)


def main():
    changed = []
    scanned = 0
    for root, dirs, files in os.walk('.'):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for fname in files:
            if not should_process(fname):
                continue
            scanned += 1
            path = os.path.join(root, fname)
            if process_file(path):
                changed.append(path)
                print(f"  fixed: {path}")

    print()
    print(f"Scanned {scanned} files.")
    print(f"Changed {len(changed)} files.")
    if changed:
        print()
        print("Next: in R:")
        print("    devtools::document()")
        print("    devtools::build_readme()")
        print("    devtools::build_vignettes()")
        print("    pkgdown::build_site()")
        print("    devtools::build_manual()")
        print()
        print("Then verify no leftovers:")
        print("    grep -rn 'Lin, L.-H. and Jalili' . \\")
        print("         --exclude-dir=.Rproj.user --exclude-dir=.git \\")
        print("         --exclude=fix_citations.py --exclude=fix_citations_v2.py")


if __name__ == '__main__':
    main()
