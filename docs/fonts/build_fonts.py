"""Build the site's web fonts from the Nerd Fonts release.

Every text on the site is set in a Nerd Font: Noto Sans for the body, Noto
Serif for headings and RobotoMono for code. Each patched file is about 3 MB,
almost all of it icons, so each face is cut into woff2 subsets by Unicode
range, as Google Fonts serves them, and a browser fetches only the subsets a
page uses: ``latin`` (Google's Latin set), ``latin-ext`` (European accented
letters), ``greek``, ``symbols`` (super- and subscripts, letterlike forms,
arrows, math operators) and ``ext`` (Cyrillic, Vietnamese, IPA and the rest).
The icons are the same in every patched font,
so they are served once, from Symbols Nerd Font, and only fetched by a page
that uses one.

usage: pip install fonttools brotli; python docs/fonts/build_fonts.py

Writes the woff2 files and the font licenses next to this script, and the
@font-face rules to docs/stylesheets/fonts.css.
"""

import urllib.request
from pathlib import Path

from fontTools import subset
from fontTools.ttLib import TTFont

TAG = "v3.5.1"
RAW = f"https://raw.githubusercontent.com/ryanoasis/nerd-fonts/{TAG}/patched-fonts/"
HERE = Path(__file__).resolve().parent
CSS = HERE.parent / "stylesheets" / "fonts.css"
CACHE = HERE / ".cache"

# (CSS family, weight, style, file under patched-fonts/)
FACES = [
    ("Noto Sans Nerd Font Propo", 400, "normal", "Noto/Sans/NotoSansNerdFontPropo-Regular.ttf"),
    ("Noto Sans Nerd Font Propo", 400, "italic", "Noto/Sans/NotoSansNerdFontPropo-Italic.ttf"),
    ("Noto Sans Nerd Font Propo", 600, "normal", "Noto/Sans/NotoSansNerdFontPropo-SemiBold.ttf"),
    ("Noto Sans Nerd Font Propo", 700, "normal", "Noto/Sans/NotoSansNerdFontPropo-Bold.ttf"),
    ("Noto Sans Nerd Font Propo", 700, "italic", "Noto/Sans/NotoSansNerdFontPropo-BoldItalic.ttf"),
    ("Noto Serif Nerd Font Propo", 600, "normal", "Noto/Serif/NotoSerifNerdFontPropo-SemiBold.ttf"),
    ("Noto Serif Nerd Font Propo", 700, "normal", "Noto/Serif/NotoSerifNerdFontPropo-Bold.ttf"),
    ("RobotoMono Nerd Font", 400, "normal", "RobotoMono/RobotoMonoNerdFont-Regular.ttf"),
    ("RobotoMono Nerd Font", 400, "italic", "RobotoMono/RobotoMonoNerdFont-Italic.ttf"),
    ("RobotoMono Nerd Font", 700, "normal", "RobotoMono/RobotoMonoNerdFont-Bold.ttf"),
]
SYMBOLS = ("Symbols Nerd Font", "NerdFontsSymbolsOnly/SymbolsNerdFont-Regular.ttf")
LICENSES = {
    "OFL-NotoSans.txt": "Noto/Sans/OFL.txt",
    "OFL-NotoSerif.txt": "Noto/Serif/OFL.txt",
    "LICENSE-RobotoMono.txt": "RobotoMono/LICENSE.txt",
    "LICENSE-SymbolsNerdFont.txt": "NerdFontsSymbolsOnly/LICENSE",
}

# Nerd Fonts puts its icons in the Private Use Areas.
ICONS = [(0xE000, 0xF8FF), (0xF0000, 0x10FFFF)]
# Each subset takes its ranges minus those of the subsets before it; ``ext``
# takes whatever is left outside the icons.
SUBSETS = [
    ("latin", [
        (0x0000, 0x00FF), (0x0131, 0x0131), (0x0152, 0x0153), (0x02BB, 0x02BC),
        (0x02C6, 0x02C6), (0x02DA, 0x02DA), (0x02DC, 0x02DC), (0x0304, 0x0304),
        (0x0308, 0x0308), (0x0329, 0x0329), (0x2000, 0x206F), (0x20AC, 0x20AC),
        (0x2122, 0x2122), (0x2190, 0x2193), (0x2212, 0x2212), (0x2215, 0x2215),
        (0xFEFF, 0xFEFF), (0xFFFD, 0xFFFD),
    ]),
    ("latin-ext", [(0x0100, 0x024F)]),
    ("greek", [(0x0370, 0x03FF), (0x1F00, 0x1FFF)]),
    ("symbols", [(0x2070, 0x209F), (0x2100, 0x23FF), (0x25A0, 0x25FF)]),
]


def subtract(spans, taken):
    """The parts of ``spans`` outside every range in ``taken``."""
    out = []
    for lo, hi in spans:
        for t_lo, t_hi in sorted(taken):
            if t_hi < lo or t_lo > hi:
                continue
            if t_lo > lo:
                out.append((lo, t_lo - 1))
            lo = t_hi + 1
            if lo > hi:
                break
        if lo <= hi:
            out.append((lo, hi))
    return out


def partition():
    taken, parts = list(ICONS), []
    for name, spans in SUBSETS:
        part = subtract(spans, taken)
        parts.append((name, part))
        taken += part
    parts.append(("ext", subtract([(0x0000, 0x10FFFF)], taken)))
    return parts


def unicode_range(spans):
    return ", ".join(f"U+{lo:04X}" if lo == hi else f"U+{lo:04X}-{hi:04X}" for lo, hi in spans)


def fetch(path):
    target = CACHE / path
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(RAW + path) as response:
            target.write_bytes(response.read())
    return target


def cut(source, spans, out):
    """Write the glyphs of ``source`` mapped from ``spans`` to ``out`` as woff2.

    Returns the file size, or None when ``source`` maps nothing in ``spans``.
    """
    font = TTFont(source)
    wanted = [cp for cp in font.getBestCmap() if any(lo <= cp <= hi for lo, hi in spans)]
    if not wanted:
        return None
    options = subset.Options()
    options.flavor = "woff2"
    options.layout_features = ["*"]
    options.name_IDs = ["*"]
    options.name_languages = ["*"]
    options.drop_tables += ["PfEd"]  # FontForge's private notes table
    subsetter = subset.Subsetter(options)
    subsetter.populate(unicodes=wanted)
    subsetter.subset(font)
    subset.save_font(font, out, options)
    return out.stat().st_size


def face(family, weight, style, url, spans):
    return (
        "@font-face {\n"
        f'  font-family: "{family}";\n'
        f"  font-style: {style};\n"
        f"  font-weight: {weight};\n"
        "  font-display: swap;\n"
        f'  src: url("../fonts/{url}") format("woff2");\n'
        f"  unicode-range: {unicode_range(spans)};\n"
        "}\n"
    )


def main():
    rules, parts = [], partition()
    for family, weight, style, path in FACES:
        source = fetch(path)
        for name, spans in parts:
            out = HERE / f"{Path(path).stem}.{name}.woff2"
            size = cut(source, spans, out)
            if size is None:
                continue
            print(f"{out.name:52} {size / 1024:7.1f} KiB")
            rules.append(face(family, weight, style, out.name, spans))
    family, path = SYMBOLS
    out = HERE / f"{Path(path).stem}.woff2"
    size = cut(fetch(path), ICONS, out)
    print(f"{out.name:52} {size / 1024:7.1f} KiB")
    rules.append(face(family, "100 900", "normal", out.name, ICONS))

    for name, path in LICENSES.items():
        (HERE / name).write_bytes(fetch(path).read_bytes())

    header = (
        f"/* Generated by docs/fonts/build_fonts.py from Nerd Fonts {TAG}; do not edit.\n"
        " * Noto Sans and Noto Serif: SIL OFL 1.1. RobotoMono: Apache 2.0.\n"
        " * Symbols Nerd Font: MIT. License texts are in docs/fonts/. */\n\n"
    )
    CSS.write_text(header + "\n".join(rules))


if __name__ == "__main__":
    main()
