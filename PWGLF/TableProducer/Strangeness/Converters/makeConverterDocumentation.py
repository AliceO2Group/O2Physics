#!/usr/bin/env python3
"""
converter_tables.py -- summarise O2Physics converter tasks.

Scans a directory of O2Physics workflow sources (by default the PWGLF
strangeness Converters directory), works out which AOD tables each task
subscribes to and which it produces, and writes the result as a Markdown
table (README.md) or as JSON/CSV.

Typical use:

    ./converter_tables.py                       # print Markdown to stdout
    ./converter_tables.py -o README.md          # write the README in place
    ./converter_tables.py --check               # CI: fail if README is stale
    ./converter_tables.py --format json

The parser is regex based: it does not need a compiler, O2 headers or a
build environment, and it deliberately errs on the side of reporting what
it cannot parse (see the --strict flag) rather than silently dropping it.

Recognised constructs
---------------------
  outputs   Produces<aod::T>, Spawns<aod::T>, Builds<aod::T>
  inputs    arguments of every void process*(...) member function,
            including soa::Join<...>, soa::Filtered<...>, ::iterator,
            and file-scope `using X = ...;` aliases
  extras    PROCESS_SWITCH(...) descriptions and default on/off state,
            the descriptive comment above the struct,
            resolution of unversioned table names (StraEvSels) to the
            concrete version they currently alias (StraEvSels_006)

Written for the ALICE O2Physics repository; stdlib only, Python >= 3.8.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import re
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# --------------------------------------------------------------------------
# Source cleaning
# --------------------------------------------------------------------------

_STRING_RE = re.compile(r'"(?:\\.|[^"\\])*"')


def strip_comments(text: str) -> str:
    """Remove // and /* */ comments, preserving line structure.

    Line count is preserved so that reported line numbers stay meaningful.
    String literals are protected so that "http://..." is not eaten.
    """
    out = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c == '"':
            m = _STRING_RE.match(text, i)
            if m:
                out.append(m.group(0))
                i = m.end()
                continue
            out.append(c)
            i += 1
        elif c == "'" and i + 2 < n:
            j = text.find("'", i + 1)
            while j != -1 and text[j - 1] == "\\":
                j = text.find("'", j + 1)
            if j != -1:
                out.append(text[i:j + 1])
                i = j + 1
                continue
            out.append(c)
            i += 1
        elif text.startswith("//", i):
            j = text.find("\n", i)
            i = n if j == -1 else j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            chunk = text[i:] if j == -1 else text[i:j + 2]
            out.append("\n" * chunk.count("\n"))
            i = n if j == -1 else j + 2
        else:
            out.append(c)
            i += 1
    return "".join(out)


# --------------------------------------------------------------------------
# Small C++-ish helpers
# --------------------------------------------------------------------------

TABLE_RE = re.compile(r"\b(?:o2\s*::\s*)?aod\s*::\s*(\w+)")


def match_paren(text: str, open_pos: int) -> int:
    """Index of the ')' matching the '(' at open_pos (-1 if unbalanced)."""
    depth = 0
    for i in range(open_pos, len(text)):
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                return i
    return -1


def split_top_level(text: str, sep: str = ",") -> List[str]:
    """Split on `sep`, ignoring separators nested in <>, (), [] or {}."""
    parts, depth, cur = [], 0, []
    for ch in text:
        if ch in "<([{":
            depth += 1
        elif ch in ">)]}":
            depth -= 1
        if ch == sep and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return [p.strip() for p in parts if p.strip()]


# --------------------------------------------------------------------------
# Data model: `using X = X_00N;` aliases and versioned declarations
# --------------------------------------------------------------------------

ALIAS_RE = re.compile(r"^\s*using\s+(\w+)\s*=\s*([^;]+);", re.M)
VERSIONED_RE = re.compile(
    r"DECLARE_SOA_TABLE_VERSIONED\s*\(\s*(\w+)\s*,\s*\"(\w+)\"\s*,\s*\"([^\"]+)\"\s*,\s*(\d+)"
)
PLAIN_TABLE_RE = re.compile(
    r"DECLARE_SOA_TABLE(?:_FULL|_STAGED)?\s*\(\s*(\w+)\s*,\s*\"(\w+)\"\s*,\s*\"([^\"]+)\""
)


@dataclass
class DataModel:
    """Maps unversioned table names to the version they currently alias."""

    alias: Dict[str, str] = field(default_factory=dict)   # StraEvSels -> StraEvSels_006
    declared: Dict[str, str] = field(default_factory=dict)  # name -> "AOD/STRAEVSELS/6"

    def resolve(self, name: str) -> Optional[str]:
        """Return Name_00N if `name` is an alias for a versioned table."""
        target = self.alias.get(name)
        if not target or target == name:
            return None
        # only report genuine version aliases (StraEvSels -> StraEvSels_006),
        # not renames such as V0Cores -> V0CoresBase
        return target if re.fullmatch(re.escape(name) + r"_\d+", target) else None

    @classmethod
    def from_paths(cls, paths: Iterable[Path]) -> "DataModel":
        dm = cls()
        for path in paths:
            if not path.exists():
                continue
            headers = [path] if path.is_file() else sorted(path.rglob("*.h"))
            for header in headers:
                try:
                    src = strip_comments(header.read_text(errors="replace"))
                except OSError:
                    continue
                for name, origin, desc, ver in VERSIONED_RE.findall(src):
                    dm.declared[name] = f"{origin}/{desc}/v{ver}"
                for name, origin, desc in PLAIN_TABLE_RE.findall(src):
                    dm.declared.setdefault(name, f"{origin}/{desc}")
                for lhs, rhs in ALIAS_RE.findall(src):
                    rhs = rhs.strip()
                    if re.fullmatch(r"(?:o2::)?(?:aod::)?\w+", rhs):
                        dm.alias[lhs] = rhs.split("::")[-1]
        # follow alias chains (A = B, B = B_002)
        for name in list(dm.alias):
            seen, cur = {name}, dm.alias[name]
            while cur in dm.alias and cur not in seen:
                seen.add(cur)
                cur = dm.alias[cur]
            dm.alias[name] = cur
        return dm


# --------------------------------------------------------------------------
# Parsed representation of one workflow file
# --------------------------------------------------------------------------

@dataclass
class Subscription:
    """One argument of a process function: one or more joined tables."""

    tables: List[str]
    joined: bool = False
    filtered: bool = False
    iterator: bool = False
    raw: str = ""


@dataclass
class ProcessFn:
    name: str
    subscriptions: List[Subscription] = field(default_factory=list)
    description: str = ""
    default_on: Optional[bool] = None  # None: no PROCESS_SWITCH (always on)


@dataclass
class Task:
    file: str
    struct: str
    comment: str = ""
    outputs: List[str] = field(default_factory=list)
    processes: List[ProcessFn] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def inputs(self) -> List[str]:
        seen, flat = set(), []
        for p in self.processes:
            for s in p.subscriptions:
                for t in s.tables:
                    if t not in seen:
                        seen.add(t)
                        flat.append(t)
        return flat


STRUCT_RE = re.compile(r"^\s*struct\s+(\w+)\s*(?::[^{]*)?\{", re.M)
PRODUCES_RE = re.compile(r"\b(Produces|Spawns|Builds)\s*<\s*([^>]+?)\s*>")
PROCESS_RE = re.compile(r"\b(?:void|auto)\s+(process\w*)\s*\(")
SWITCH_RE = re.compile(
    r"PROCESS_SWITCH\s*\(\s*(\w+)\s*,\s*(\w+)\s*,\s*\"([^\"]*)\"\s*,\s*(true|false)\s*\)"
)
LOCAL_ALIAS_RE = re.compile(r"^\s*using\s+(\w+)\s*=\s*([^;]+);", re.M)


def struct_body(src: str, start: int) -> str:
    """Return the body of the struct whose opening brace follows `start`."""
    brace = src.find("{", start)
    if brace == -1:
        return ""
    depth = 0
    for i in range(brace, len(src)):
        if src[i] == "{":
            depth += 1
        elif src[i] == "}":
            depth -= 1
            if depth == 0:
                return src[brace + 1:i]
    return src[brace + 1:]


def leading_comment(raw: str, struct_pos: int) -> str:
    """The // comment block immediately above a struct, as one line."""
    lines = raw[:struct_pos].splitlines()
    picked: List[str] = []
    for line in reversed(lines):
        s = line.strip()
        if not s:
            if picked:
                break
            continue
        if s.startswith("//"):
            body = s.lstrip("/").strip()
            if body:
                picked.append(body)
            continue
        break
    return " ".join(reversed(picked)).strip()


def parse_subscription(arg: str, aliases: Dict[str, str]) -> Optional[Subscription]:
    """Turn one process-function argument into a Subscription."""
    expr = arg
    # expand file-scope using aliases, repeatedly (aliases may nest)
    for _ in range(8):
        expanded = re.sub(
            r"\b(\w+)\b",
            lambda m: aliases.get(m.group(1), m.group(1)),
            expr,
        )
        if expanded == expr:
            break
        expr = expanded

    tables = []
    for name in TABLE_RE.findall(expr):
        if name not in tables:
            tables.append(name)
    if not tables:
        return None
    return Subscription(
        tables=tables,
        joined="Join" in expr and len(tables) > 1,
        filtered="Filtered" in expr,
        iterator="::iterator" in expr.replace(" ", ""),
        raw=" ".join(arg.split()),
    )


def parse_file(path: Path) -> List[Task]:
    raw = path.read_text(errors="replace")
    src = strip_comments(raw)

    aliases: Dict[str, str] = {}
    for lhs, rhs in LOCAL_ALIAS_RE.findall(src):
        aliases[lhs] = rhs.strip()

    switches = {(s[0], s[1]): (s[2], s[3] == "true") for s in SWITCH_RE.findall(src)}

    tasks: List[Task] = []
    for m in STRUCT_RE.finditer(src):
        struct = m.group(1)
        body = struct_body(src, m.start())
        # the descriptive comment lives in the uncleaned source
        raw_pos = raw.find(f"struct {struct}")
        task = Task(
            file=path.name,
            struct=struct,
            comment=leading_comment(raw, raw_pos) if raw_pos != -1 else "",
        )

        for kind, tmpl in PRODUCES_RE.findall(body):
            name = tmpl.split("::")[-1].strip()
            label = name if kind == "Produces" else f"{name} ({kind.lower()})"
            if label not in task.outputs:
                task.outputs.append(label)

        for pm in PROCESS_RE.finditer(body):
            open_paren = body.index("(", pm.end() - 1)
            close = match_paren(body, open_paren)
            if close == -1:
                task.warnings.append(f"unbalanced parentheses in {pm.group(1)}()")
                continue
            fn = ProcessFn(name=pm.group(1))
            for arg in split_top_level(body[open_paren + 1:close]):
                sub = parse_subscription(arg, aliases)
                if sub:
                    fn.subscriptions.append(sub)
                elif arg.strip():
                    task.warnings.append(
                        f"{fn.name}(): unrecognised argument '{' '.join(arg.split())}'"
                    )
            desc, default_on = switches.get((struct, fn.name), (None, None))
            fn.description = desc or ""
            fn.default_on = default_on
            task.processes.append(fn)

        if not task.outputs:
            task.warnings.append("no Produces<> found")
        if not task.processes:
            task.warnings.append("no process() function found")
        tasks.append(task)
    return tasks


# --------------------------------------------------------------------------
# Consistency checks on the sources themselves
# --------------------------------------------------------------------------

VERSION_SUFFIX_RE = re.compile(r"^(.*?)_(\d+)$")
COMMENT_RANGE_RE = re.compile(
    r"(?:version\s+|from\s+|v)?(\d{1,3})\s*(?:to|->|-->|=>)\s*v?(\d{1,3})", re.I
)


def split_version(name: str) -> Tuple[str, Optional[int]]:
    """'StraEvSels_006' -> ('StraEvSels', 6); 'V0Extras' -> ('V0Extras', None)."""
    m = VERSION_SUFFIX_RE.match(name)
    return (m.group(1), int(m.group(2))) if m else (name, None)


def family(task: "Task") -> str:
    """Table family a converter belongs to, taken from its first output."""
    if task.outputs:
        return split_version(task.outputs[0].split(" ")[0])[0]
    return "other"


def lint_task(task: "Task", dm: DataModel) -> List[str]:
    """Flag descriptive comments that disagree with the parsed code.

    Several converters in the repository carry copy-pasted comments (e.g.
    'Converts V0 version 001 to 002' on a task that converts StraRawCents
    003 to 004). These are harmless to the build but actively mislead the
    reader, which is exactly who this summary is for.
    """
    issues: List[str] = []
    if not task.comment:
        return issues

    out_names = [o.split(" ")[0] for o in task.outputs]
    in_names = [t for p in task.processes for s in p.subscriptions for t in s.tables]
    out_ver = next((v for _, v in map(split_version, out_names) if v is not None), None)
    in_ver = next((v for _, v in map(split_version, in_names) if v is not None), None)

    m = COMMENT_RANGE_RE.search(task.comment)
    if m and in_ver is not None and out_ver is not None:
        said = (int(m.group(1)), int(m.group(2)))
        if said != (in_ver, out_ver):
            issues.append(
                f"comment says {said[0]} -> {said[1]} but the task converts "
                f"{in_ver:03d} -> {out_ver:03d}"
            )

    # a table family named in the comment that the task never touches
    mentioned = set(re.findall(r"\b([A-Z][A-Za-z0-9]{3,})\b", task.comment))
    touched = {split_version(n)[0] for n in in_names + out_names}
    touched |= {n for n in in_names + out_names}
    stray = {w for w in mentioned
             if (w in dm.declared or w in dm.alias)
             and split_version(w)[0] not in touched and w not in touched}
    if stray:
        issues.append("comment mentions " + ", ".join(sorted(stray))
                      + " which the task does not use")
    return issues


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------

def fmt_subscription(sub: Subscription, dm: DataModel, resolve: bool) -> str:
    names = []
    for t in sub.tables:
        target = dm.resolve(t) if resolve else None
        names.append(f"`{t}`" + (f" *(={target})*" if target else ""))
    text = " + ".join(names)
    if sub.filtered:
        text += " *(filtered)*"
    return text


def fmt_inputs(task: Task, dm: DataModel, resolve: bool) -> str:
    chunks = []
    multi = len(task.processes) > 1
    for p in task.processes:
        if not p.subscriptions:
            body = "*(none)*"
        else:
            body = ", ".join(fmt_subscription(s, dm, resolve) for s in p.subscriptions)
        if multi:
            flag = "" if p.default_on is None else ("" if p.default_on else " *(off by default)*")
            chunks.append(f"**{p.name}**: {body}{flag}")
        else:
            chunks.append(body)
    return "<br>".join(chunks) if chunks else "*(none)*"


def fmt_outputs(task: Task, dm: DataModel, resolve: bool) -> str:
    if not task.outputs:
        return "*(none)*"
    out = []
    for o in task.outputs:
        base = o.split(" ")[0]
        target = dm.resolve(base) if resolve else None
        out.append(f"`{o}`" + (f" *(={target})*" if target else ""))
    return "<br>".join(out)


HEADER = """<!-- ############################################################ -->
<!-- This file is generated by converter_tables.py -- do not edit. -->
<!-- Regenerate with:  ./converter_tables.py -o README.md          -->
<!-- ############################################################ -->

# Strangeness table converters

Converter tasks migrate old versions of the strangeness derived-data tables to
newer ones, so that analyses written against the current data model can still
run over older derived data. Each task reads one or more legacy tables and
produces the corresponding newer version, filling any column that did not exist
in the old table with a dummy value.

Add the workflow you need to your job with `o2-analysis-lf-<name>`; only the
converters whose *input* tables are actually present in your derived data need
to be enabled.

The table below is generated automatically from the sources.
"""

FOOTER = """
Notes:

* Tables joined with `+` are subscribed together in a single argument
  (`soa::Join<...>`); a comma separates independent arguments.
* `Name *(=Name_00N)*` means the task subscribes to (or produces) the
  unversioned alias, which currently resolves to version `N` in the data model.
* Tasks with more than one `process` function list each one separately;
  a switch marked *off by default* must be enabled explicitly.
"""


def _table(buf: io.StringIO, tasks: List[Task], dm: DataModel,
           resolve: bool, show_comments: bool) -> None:
    cols = ["Converter", "Reads", "Produces"]
    if show_comments:
        cols.append("Description")
    buf.write("| " + " | ".join(cols) + " |\n")
    buf.write("|" + "|".join(["---"] * len(cols)) + "|\n")
    for t in sorted(tasks, key=lambda x: x.file):
        row = [
            f"`{t.file}`",
            fmt_inputs(t, dm, resolve),
            fmt_outputs(t, dm, resolve),
        ]
        if show_comments:
            row.append(t.comment or "")
        buf.write("| " + " | ".join(c.replace("|", "\\|") for c in row) + " |\n")


def render_markdown(tasks: List[Task], dm: DataModel, resolve: bool,
                    show_comments: bool, with_header: bool,
                    group: bool = False) -> str:
    buf = io.StringIO()
    if with_header:
        buf.write(HEADER)
        buf.write("\n")
    if not group:
        _table(buf, tasks, dm, resolve, show_comments)
    else:
        groups: Dict[str, List[Task]] = {}
        for t in tasks:
            groups.setdefault(family(t), []).append(t)
        for fam in sorted(groups):
            buf.write(f"\n## {fam}\n\n")
            _table(buf, groups[fam], dm, resolve, show_comments)
    if with_header:
        buf.write(FOOTER)
    return buf.getvalue()


def render_json(tasks: List[Task], dm: DataModel) -> str:
    payload = []
    for t in sorted(tasks, key=lambda x: x.file):
        d = asdict(t)
        d["inputs"] = t.inputs
        d["resolved"] = {n: dm.resolve(n) for n in t.inputs + t.outputs if dm.resolve(n)}
        payload.append(d)
    return json.dumps(payload, indent=2)


def render_csv(tasks: List[Task]) -> str:
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["file", "struct", "process", "reads", "produces", "description"])
    for t in sorted(tasks, key=lambda x: x.file):
        for p in t.processes or [ProcessFn(name="process")]:
            reads = "; ".join(" + ".join(s.tables) for s in p.subscriptions)
            w.writerow([t.file, t.struct, p.name, reads,
                        "; ".join(t.outputs), t.comment])
    return buf.getvalue()


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

# DEFAULT_DIR = "PWGLF/TableProducer/Strangeness/Converters"
DEFAULT_DIR = "."


def find_repo_root(start: Path) -> Optional[Path]:
    for parent in [start.resolve()] + list(start.resolve().parents):
        if (parent / ".git").exists() or (parent / "PWGLF").is_dir():
            return parent
    return None


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Summarise the input/output tables of O2Physics converter tasks."
    )
    ap.add_argument("directory", nargs="?", default=DEFAULT_DIR,
                    help=f"directory of .cxx workflows (default: {DEFAULT_DIR})")
    ap.add_argument("-o", "--output", help="write to this file instead of stdout")
    ap.add_argument("--format", choices=["markdown", "json", "csv"], default="markdown")
    ap.add_argument("--datamodel", action="append", default=[],
                    help="extra header file or directory to scan for table "
                         "declarations and aliases (repeatable)")
    ap.add_argument("--no-resolve", action="store_true",
                    help="do not annotate unversioned names with their current version")
    ap.add_argument("--no-comments", action="store_true",
                    help="omit the Description column")
    ap.add_argument("--group", action="store_true",
                    help="split the table into one section per table family")
    ap.add_argument("--lint", action="store_true",
                    help="also report source comments that contradict the code")
    ap.add_argument("--bare", action="store_true",
                    help="emit only the Markdown table, without the surrounding text")
    ap.add_argument("--check", action="store_true",
                    help="exit 2 if --output exists and differs (for CI)")
    ap.add_argument("--strict", action="store_true",
                    help="exit 1 if any file could not be fully parsed")
    args = ap.parse_args(argv)

    directory = Path(args.directory)
    if not directory.is_dir():
        print(f"error: {directory} is not a directory", file=sys.stderr)
        return 1

    sources = sorted(p for p in directory.glob("*.cxx"))
    if not sources:
        print(f"error: no .cxx files in {directory}", file=sys.stderr)
        return 1

    tasks: List[Task] = []
    for src in sources:
        tasks.extend(parse_file(src))

    dm_paths = [Path(p) for p in args.datamodel]
    if not dm_paths:
        root = find_repo_root(directory)
        if root:
            dm_paths = [root / "PWGLF" / "DataModel", root / "Common" / "DataModel"]
    dm = DataModel.from_paths(dm_paths)

    if args.format == "json":
        text = render_json(tasks, dm)
    elif args.format == "csv":
        text = render_csv(tasks)
    else:
        text = render_markdown(tasks, dm, not args.no_resolve,
                               not args.no_comments, not args.bare,
                               group=args.group)

    problems = [(t.file, w) for t in tasks for w in t.warnings]
    if args.lint:
        for t in tasks:
            problems.extend((t.file, f"[lint] {i}") for i in lint_task(t, dm))
    for f, w in problems:
        print(f"warning: {f}: {w}", file=sys.stderr)

    if args.output:
        out = Path(args.output)
        if args.check:
            existing = out.read_text() if out.exists() else None
            if existing != text:
                print(f"{out} is out of date; regenerate with "
                      f"`{Path(sys.argv[0]).name} -o {out}`", file=sys.stderr)
                return 2
            print(f"{out} is up to date ({len(tasks)} tasks).", file=sys.stderr)
            return 0
        out.write_text(text)
        print(f"wrote {out} ({len(tasks)} tasks, "
              f"{len(sources)} files).", file=sys.stderr)
    else:
        sys.stdout.write(text)

    return 1 if (args.strict and problems) else 0


if __name__ == "__main__":
    raise SystemExit(main())
