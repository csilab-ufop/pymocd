"""MkDocs hooks for the pymocd site.

* Copy the generated ``pymocd.pyi`` stub where mkdocstrings reads it, rewriting
  pyo3's ``__new__`` constructors as ``__init__``.
* Copy the logo and the favicon's otter face from ``res/`` and derive their
  dark-mode variants.
* Fail the build when a public name of the stub has no ``::: pymocd.<name>``
  directive under ``docs/api/``, or a directive names something the stub lacks.
* Render the API summary on ``api/index.md`` from the stub.
* Check every parameter table marked ``<!-- pymocd:defaults ... -->`` against
  the defaults in the stub, so a page cannot drift from the code.
* Fail the build on a ```python block that markdown-exec would not execute,
  so every Python example on the site is run against the current package.
* Point the edit button of generated gallery pages at the example's source.
* List the ``<!-- TODO(docs): ... -->`` markers left in the pages.
"""

import ast
import logging
import math
import re
import shutil
from fractions import Fraction
from pathlib import Path

from mkdocs.exceptions import PluginError

ROOT = Path(__file__).resolve().parent.parent
STUB_DIR = ROOT / ".docs-stub"
REPO = "https://github.com/csilab-ufop/pymocd"

log = logging.getLogger("mkdocs.hooks.pymocd")

# res/ colour -> dark-mode replacement: the navy and blue strokes of the logo
# vanish on the slate background, so they are lifted to lighter tints.
DARK_COLOURS = {
    "#0F3D75": "#8DB3E8",
    "#0B4DA2": "#5C9CF0",
    "#C7D3E6": "#5A6B85",
}

API_PAGES = ("api/detectors.md", "api/fronts.md", "api/metrics.md", "api/generators.md")
DIRECTIVE = re.compile(r"^:::\s+pymocd\.(\w+)\s*$", re.M)
DEFAULTS_MARKER = re.compile(r"<!--\s*pymocd:defaults\s+([\w. ]+?)\s*-->")
TODO_MARKER = re.compile(r"<!--\s*TODO\(docs\):\s*(.*?)\s*-->", re.S)
FENCE = re.compile(r"^\s*(`{3,}|~{3,})(.*)$")
EXECUTED = re.compile(r'\bexec="(?:on|true|1)"')
ABBREVIATIONS = ("al.", "e.g.", "i.e.", "Eq.", "Eqs.", "Sec.", "vs.", "cf.", "Fig.")

_api = {}
_documented = {}
_todos = []


def on_startup(command, dirty):
    stub = ROOT / "pymocd.pyi"
    if not stub.exists():
        raise FileNotFoundError("pymocd.pyi not found — run `make stubs` first")
    STUB_DIR.mkdir(exist_ok=True)
    text = re.sub(
        r"def __new__\(cls(.*)\) -> [\w.]+:( \.\.\.)?$",
        r"def __init__(self\1) -> None:\2",
        stub.read_text(),
        flags=re.M,
    )
    (STUB_DIR / "pymocd.py").write_text(text)

    assets = ROOT / "docs" / "assets"
    assets.mkdir(exist_ok=True)
    # The header shows the favicon's otter face. Its colour must follow the
    # site's palette toggle, not the OS scheme, so the favicon's dark-mode
    # media query is dropped and a dark variant is derived like the logo's.
    face = re.sub(
        r"^\s*@media \(prefers-color-scheme: dark\).*\n",
        "",
        (ROOT / "res" / "favicon.svg").read_text(),
        flags=re.M,
    )
    for name, svg in (("logo", (ROOT / "res" / "logo.svg").read_text()), ("face", face)):
        (assets / f"{name}.svg").write_text(svg)
        for light, dark in DARK_COLOURS.items():
            svg = re.sub(re.escape(light), dark, svg, flags=re.I)
        (assets / f"{name}-dark.svg").write_text(svg)
    for name in ("logo.png", "favicon.svg", "favicon-32.png", "apple-touch-icon.png"):
        shutil.copy(ROOT / "res" / name, assets / name)


def on_config(config):
    _api.clear()
    _documented.clear()
    _todos.clear()
    tree = ast.parse((STUB_DIR / "pymocd.py").read_text())
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and not node.name.startswith("_"):
            _api[node.name] = node

    for rel in API_PAGES:
        for name in DIRECTIVE.findall((ROOT / "docs" / rel).read_text()):
            _documented.setdefault(name, rel)
    missing = sorted(set(_api) - set(_documented))
    unknown = sorted(set(_documented) - set(_api))
    if missing or unknown:
        raise PluginError(
            "API reference out of sync with pymocd.pyi — "
            f"undocumented: {missing or 'none'}; not in the stub: {unknown or 'none'}"
        )
    return config


def on_page_markdown(markdown, page, config, files):
    src = page.file.src_uri
    for todo in TODO_MARKER.findall(markdown):
        _todos.append((src, " ".join(todo.split())))
    if src == "api/index.md":
        markdown = markdown.replace("<!-- pymocd:api-summary -->", _api_summary())
    for match in DEFAULTS_MARKER.finditer(markdown):
        _check_defaults(src, match.group(1).split(), markdown[match.end():])
    if not src.startswith("examples/"):
        _check_executed(src, markdown)
    return markdown


def on_page_context(context, page, config, nav):
    src = page.file.src_uri
    if src.startswith("examples/"):
        rel = src.removeprefix("examples/")
        source = rel.replace("index.md", "README.md") if rel.endswith("index.md") else rel.replace(".md", ".py")
        page.edit_url = f"{REPO}/edit/master/examples/{source}"
    return context


def on_post_build(config):
    site = Path(config["site_dir"])
    html = "".join(p.read_text() for p in (site / "api").rglob("*.html"))
    absent = sorted(n for n in _api if f'id="pymocd.{n}"' not in html)
    if absent:
        raise PluginError(f"API reference pages render no anchor for: {absent}")
    if _todos:
        log.info("%d TODO(docs) marker(s) left in the pages:", len(_todos))
        for src, text in _todos:
            log.info("  %s: %s", src, text)


def _check_executed(src, markdown):
    """Raise for a top-level ```python fence without exec="on" (gallery pages
    are executed by mkdocs-gallery and skipped; nested fences are examples)."""
    open_fence, unexecuted = None, []
    for number, line in enumerate(markdown.splitlines(), 1):
        match = FENCE.match(line)
        if not match:
            continue
        marker, info = match.groups()
        if open_fence is None:
            open_fence = marker
            lang = info.strip().split(" ")[0].lower() if info.strip() else ""
            if lang in ("python", "py") and not EXECUTED.search(info):
                unexecuted.append(number)
        elif marker[0] == open_fence[0] and len(marker) >= len(open_fence) and not info.strip():
            open_fence = None
    if unexecuted:
        raise PluginError(
            f'{src}: ```python block(s) at line(s) {unexecuted} lack exec="on"; every Python '
            "example must run at build time (use ```text for code that must not run)"
        )


def _summary(doc):
    """First sentence of a docstring, as one line of Markdown."""
    para = " ".join(doc.strip().split("\n\n")[0].split()).replace("``", "`")
    words, out = para.split(" "), []
    for word in words:
        out.append(word)
        if word.endswith(".") and not word.endswith(ABBREVIATIONS):
            break
    return " ".join(out)


def _api_summary():
    by_page = {}
    for name, rel in _documented.items():
        by_page.setdefault(rel, []).append(name)
    blocks = []
    for rel in API_PAGES:
        text = (ROOT / "docs" / rel).read_text()
        title = re.search(r"^# (.+)$", text, re.M).group(1).strip()
        order = DIRECTIVE.findall(text)
        rows = ["| Name | Kind | Summary |", "|---|---|---|"]
        for name in sorted(by_page.get(rel, []), key=order.index):
            node = _api[name]
            kind = "class" if isinstance(node, ast.ClassDef) else "function"
            doc = _summary(ast.get_docstring(node) or "")
            rows.append(f"| [`{name}`][pymocd.{name}] | {kind} | {doc} |")
        page = rel.removeprefix("api/")
        blocks.append(f"## [{title}]({page})\n\n" + "\n".join(rows))
    return "\n\n".join(blocks)


def _signature(name):
    node = _api.get(name)
    if node is None:
        return None
    if isinstance(node, ast.ClassDef):
        node = next(n for n in node.body if isinstance(n, ast.FunctionDef) and n.name == "__init__")
    args = [a.arg for a in node.args.args if a.arg not in ("self", "cls", "graph")]
    defaults = node.args.defaults
    named = node.args.args[len(node.args.args) - len(defaults):]
    return args, {a.arg: ast.literal_eval(d) for a, d in zip(named, defaults)}


def _parse_default(text):
    text = text.strip()
    frac = re.fullmatch(r"(\d+)\s*/\s*(\d+)", text)
    if frac:
        return float(Fraction(int(frac.group(1)), int(frac.group(2))))
    return ast.literal_eval(text)


def _check_defaults(src, targets, rest):
    sigs = {}
    for target in targets:
        sig = _signature(target.removeprefix("pymocd."))
        if sig is None:
            raise PluginError(f"{src}: defaults table names {target}, which is not in pymocd.pyi")
        sigs[target] = sig

    rows, started = {}, False
    for line in rest.lstrip("\n").splitlines():
        if not line.startswith("|"):
            if started:
                break
            continue
        started = True
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        name = re.fullmatch(r"`(\w+)`", cells[0])
        default = re.fullmatch(r"`([^`]+)`", cells[1]) if len(cells) > 1 else None
        if name and default:
            rows[name.group(1)] = default.group(1)

    errors = []
    for target, (params, defaults) in sigs.items():
        for param in params:
            if param not in rows:
                errors.append(f"{target}: parameter `{param}` missing from the table")
                continue
            want = defaults.get(param)
            try:
                got = _parse_default(rows[param])
            except (ValueError, SyntaxError):
                errors.append(f"{target}: `{param}` default `{rows[param]}` is not a literal")
                continue
            same = (
                math.isclose(got, want, rel_tol=1e-12)
                if isinstance(got, (int, float)) and isinstance(want, (int, float))
                else got == want
            )
            if not same:
                errors.append(f"{target}: `{param}` documented as `{rows[param]}`, stub says `{want!r}`")
    extra = set(rows) - {p for params, _ in sigs.values() for p in params}
    errors += [f"table lists `{p}`, which no target takes" for p in sorted(extra)]
    if errors:
        raise PluginError(f"{src}: parameter table disagrees with pymocd.pyi:\n  " + "\n  ".join(errors))
