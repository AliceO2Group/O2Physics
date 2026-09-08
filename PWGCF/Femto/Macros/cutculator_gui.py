#!/usr/bin/env python3

# Copyright 2019-2025 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

"""!
@brief  CutCulator GUI (Compute bitmask for selecting particles in the Femto Framework)
@author Anton Riedel <anton.riedel@cern.ch>, Technical University of Munich
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import argparse
import itertools
import json
import os
from typing import Any, Dict, List

try:
    import ROOT  # pylint: disable=import-error

    ROOT.gROOT.SetBatch(True)
except ImportError:
    ROOT = None

VALUE_DELIM = "___"
SECTION_DELIM = ":::"

# Above this number of enumerated variations the export asks for confirmation
# before building the list (the product over cuts grows multiplicatively).
VARIATION_WARN_LIMIT = 5000

# ── Colours ──────────────────────────────────────────────────────────────────
BG = "#1e1e2e"
BG_CARD = "#2a2a3e"
BG_HOVER = "#3a3a50"
ACCENT = "#89b4fa"  # blue  – minimal
ACCENT_OPT = "#a6e3a1"  # green – optional
ACCENT_REJ = "#f38ba8"  # red   – rejection / neutral
ACCENT_ALWAYS = "#cba6f7"  # purple – loosest/always-true (no bit set)
ACCENT_FILTER = "#f9e2af"  # yellow – filter values (read-only)
FG = "#cdd6f4"
FG_DIM = "#6c7086"
BORDER = "#45475a"
SEL_BG = "#313244"

FONT_TITLE = ("Inter", 15, "bold")
FONT_HEAD = ("Inter", 11, "bold")
FONT_BODY = ("Inter", 10)
FONT_MONO = ("JetBrains Mono", 11, "bold")
FONT_SMALL = ("Inter", 9)

# Leading column width shared by every row type (checkbox glyph, loosest-bin
# spacer) so all rows' text starts at exactly the same x position.
CHECK_COL_WIDTH = 3


# ── Helpers ───────────────────────────────────────────────────────────────────
def parse_bin_label(label):
    result = {}
    for sec in label.split(SECTION_DELIM):
        if VALUE_DELIM in sec:
            k, v = sec.split(VALUE_DELIM, 1)
            result[k.strip()] = v.strip()
    return result


def format_value_with_comment(b):
    val = b.get("Value", "")
    comment = b.get("Comment", "")
    if comment and comment.upper() != "X":
        return f"{val}  ({comment})"
    return val


def bit_position_int(b):
    """BitPosition as an int, or -1 for the always-applied (X) bins."""
    pos = b.get("BitPosition", "X")
    return int(pos) if pos.upper() != "X" else -1


def bin_type(b):
    """Return 'minimal', 'optional', 'neutral', or 'skip'."""
    is_min = b.get("MinimalCut", "0") == "1"
    is_opt = b.get("OptionalCut", "0") == "1"
    pos = b.get("BitPosition", "X")
    if pos.upper() == "X":
        return "skip"
    if is_min and not is_opt:
        return "minimal"
    if is_opt:
        return "optional"
    return "neutral"


def variation_kind(b):
    """Axis kind used for grouping systematic variations.

    Identical to bin_type() except that an always-applied minimal floor
    (BitPosition=="X", which bin_type() calls 'skip') is folded onto the
    'minimal' axis of its selection: it is the loosest rung of the same
    threshold ladder, so it is a legitimate variation point alongside the
    stricter, bit-carrying rungs — it just contributes no bit.
    """
    kind = bin_type(b)
    if kind == "skip" and b.get("MinimalCut", "0") == "1" and b.get("OptionalCut", "0") == "0":
        return "minimal"
    return kind


def has_only_skipped_minimal_bins(group):
    """
    Returns True when a group has at least one minimal bin but ALL of them have
    BitPosition=="X" — meaning the cut is the loosest (always-true) selection
    and no bit needs to be set. The group should still be displayed so the
    user knows the cut exists.
    """
    minimal = [b for b in group if b.get("MinimalCut", "0") == "1" and b.get("OptionalCut", "0") == "0"]
    return bool(minimal) and all(b.get("BitPosition", "X").upper() == "X" for b in minimal)


def get_skipped_minimal_bin_items(group):
    """
    Returns (index_in_group, bin) for the minimal bins whose BitPosition=="X" —
    i.e. the loosest threshold(s), always in effect unless a stricter minimal
    bit is set. They carry no bit, but they are still a selectable variation
    point (see variation_kind), so the index is needed to key self._vars.
    """
    return [
        (i, b)
        for i, b in enumerate(group)
        if b.get("MinimalCut", "0") == "1"
        and b.get("OptionalCut", "0") == "0"
        and b.get("BitPosition", "X").upper() == "X"
    ]


def get_skipped_minimal_bins(group):
    """The bins from get_skipped_minimal_bin_items(), without their indices."""
    return [b for _i, b in get_skipped_minimal_bin_items(group)]


def is_filter_group(bins):
    """
    Determine whether the parsed bins belong to a filter histogram (fixed,
    non-selectable pre-filter values) rather than a selection/bitmask histogram.
    Filter bins carry a "FilterName" key; selection bins carry "SelectionName".
    """
    for b in bins:
        if "FilterName" in b:
            return True
        if "SelectionName" in b:
            return False
    return False


def load_bins_from_hist(hist):
    nbins = hist.GetNbinsX()
    bins = []
    for i in range(1, nbins - 2 + 1):
        label = hist.GetXaxis().GetBinLabel(i)
        if not label:
            continue
        b = parse_bin_label(label)
        b["_bin_index"] = i
        bins.append(b)

    is_filter = is_filter_group(bins)
    if is_filter:
        # filter histograms are a flat list of fixed name/value pairs — no grouping needed
        return bins, True

    # group by SelectionName preserving order
    groups = {}
    for b in bins:
        name = b.get("SelectionName", f"unknown_{b['_bin_index']}")
        groups.setdefault(name, []).append(b)
    return groups, False


def resolve_file_path(path):
    """Absolute, ~-expanded path for a local file.

    Remote specifications (root://, alien://, http://, …) are returned
    untouched — os.path.abspath() would prepend the cwd and corrupt them.
    """
    if "://" in path:
        return path
    return os.path.abspath(os.path.expanduser(path))


def unique_key(mapping, base):
    """Return `base`, or `base #2`, `base #3`, … if it is already taken.

    Needed because two different axes (or an axis and an always-applied floor)
    can carry the same SelectionName; JSON objects cannot hold duplicate keys,
    so the collision has to be resolved rather than silently overwritten.
    """
    if base not in mapping:
        return base
    n = 2
    while f"{base} #{n}" in mapping:
        n += 1
    return f"{base} #{n}"


# ── Main Application ──────────────────────────────────────────────────────────
class CutCulatorApp(tk.Tk):
    def __init__(self, rootfile=None, tdir="femto-producer"):
        super().__init__()
        self.title("CutCulator")
        self.configure(bg=BG)
        self.geometry("1200x820")
        self.minsize(900, 600)
        self.resizable(True, True)

        self._rootfile_path = rootfile
        self._tdir_path = tdir
        self._root_file = None
        self._hist = None
        self._groups: Dict[str, List[Dict[str, Any]]] = {}  # SelectionName → list[bin_dict]
        self._filter_bins: List[Dict[str, Any]] = []  # list[bin_dict], used when the histogram is a filter histogram
        self._is_filter_hist = False
        self._vars = {}  # (SelectionName, idx) → BooleanVar
        self._check_labels = {}  # (SelectionName, idx) → Label (custom checkbox glyph)

        self._build_ui()

        if rootfile:
            self._load_file(rootfile)

    # ── UI skeleton ──────────────────────────────────────────────────────────
    def _build_ui(self):
        # ── top bar ──
        top = tk.Frame(self, bg=BG, pady=10, padx=18)
        top.pack(fill="x")

        tk.Label(top, text="CutCulator", font=FONT_TITLE, bg=BG, fg=ACCENT).pack(side="left")

        btn_frame = tk.Frame(top, bg=BG)
        btn_frame.pack(side="right")
        self._btn_open = self._make_button(btn_frame, "Open ROOT file", self._open_file_dialog, ACCENT)
        self._btn_open.pack(side="left", padx=4)

        # ── file + histogram selector bar ──
        bar = tk.Frame(self, bg=BG_CARD, pady=8, padx=18)
        bar.pack(fill="x")

        tk.Label(bar, text="File:", font=FONT_BODY, bg=BG_CARD, fg=FG_DIM).pack(side="left")
        self._lbl_file = tk.Label(bar, text="(none)", font=FONT_BODY, bg=BG_CARD, fg=FG, wraplength=400, anchor="w")
        self._lbl_file.pack(side="left", padx=6)

        tk.Label(bar, text="Histogram:", font=FONT_BODY, bg=BG_CARD, fg=FG_DIM).pack(side="left", padx=(20, 0))
        self._hist_var = tk.StringVar()
        self._hist_combo = ttk.Combobox(bar, textvariable=self._hist_var, state="disabled", width=30, font=FONT_BODY)
        self._hist_combo.pack(side="left", padx=6)
        self._hist_combo.bind("<<ComboboxSelected>>", self._on_hist_selected)

        self._style_combobox()

        # ── legend ──
        self._legend = tk.Frame(self, bg=BG, pady=4, padx=18)
        self._legend.pack(fill="x")
        self._build_legend_selection()

        # ── bottom result bar (packed to bottom BEFORE the main area so the
        #    main area fills remaining space correctly) ──
        self._bottom = tk.Frame(self, bg=BG_CARD, pady=10, padx=18)
        self._bottom.pack(fill="x", side="bottom")

        tk.Label(self._bottom, text="Bitmask:", font=FONT_HEAD, bg=BG_CARD, fg=FG).pack(side="left")

        self._lbl_dec = tk.Label(self._bottom, text="—", font=FONT_MONO, bg=BG_CARD, fg=ACCENT, width=14, anchor="w")
        self._lbl_dec.pack(side="left", padx=8)

        self._lbl_hex = tk.Label(
            self._bottom, text="—", font=FONT_MONO, bg=BG_CARD, fg=ACCENT_OPT, width=14, anchor="w"
        )
        self._lbl_hex.pack(side="left", padx=8)

        self._lbl_bin = tk.Label(self._bottom, text="—", font=FONT_MONO, bg=BG_CARD, fg=FG_DIM, anchor="w")
        self._lbl_bin.pack(side="left", padx=8)

        # live count of enumerated systematic variations for the current checkboxes
        self._lbl_var_count = tk.Label(self._bottom, text="", font=FONT_SMALL, bg=BG_CARD, fg=ACCENT_ALWAYS, anchor="w")
        self._lbl_var_count.pack(side="left", padx=12)

        copy_frame = tk.Frame(self._bottom, bg=BG_CARD)
        copy_frame.pack(side="right")
        self._make_button(copy_frame, "Copy Dec", lambda: self._copy(self._lbl_dec["text"]), ACCENT).pack(
            side="left", padx=3
        )
        self._make_button(copy_frame, "Copy Hex", lambda: self._copy(self._lbl_hex["text"]), ACCENT_OPT).pack(
            side="left", padx=3
        )
        self._make_button(copy_frame, "Export variations", self._export_variations, ACCENT_ALWAYS).pack(
            side="left", padx=3
        )

        # ── main paned area: selection cards (left) + summary sidebar (right) ──
        self._paned = tk.PanedWindow(self, orient="horizontal", bg=BORDER, sashwidth=4, sashrelief="flat")
        self._paned.pack(fill="both", expand=True, padx=12, pady=4)

        # ── scrollable selection area (left) ──
        outer = tk.Frame(self._paned, bg=BG)
        self._paned.add(outer, stretch="always", minsize=360)

        self._canvas = tk.Canvas(outer, bg=BG, highlightthickness=0, bd=0)
        vsb = ttk.Scrollbar(outer, orient="vertical", command=self._canvas.yview)
        self._canvas.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        self._canvas.pack(side="left", fill="both", expand=True)

        self._inner = tk.Frame(self._canvas, bg=BG)
        self._canvas_win = self._canvas.create_window((0, 0), window=self._inner, anchor="nw")
        self._inner.bind("<Configure>", self._on_inner_configure)
        self._canvas.bind("<Configure>", self._on_canvas_configure)
        self._canvas.bind_all("<MouseWheel>", self._on_mousewheel)
        self._canvas.bind_all("<Button-4>", self._on_mousewheel)
        self._canvas.bind_all("<Button-5>", self._on_mousewheel)

        # ── summary sidebar (right) ──
        self._summary_outer = tk.Frame(self._paned, bg=BG_CARD)
        self._paned.add(self._summary_outer, stretch="never", minsize=280, width=320)

        summary_hdr = tk.Frame(self._summary_outer, bg=BG_CARD, pady=8, padx=10)
        summary_hdr.pack(fill="x")
        self._lbl_summary_title = tk.Label(summary_hdr, text="Selected cuts", font=FONT_HEAD, bg=BG_CARD, fg=FG)
        self._lbl_summary_title.pack(anchor="w")
        self._lbl_summary_empty = tk.Label(summary_hdr, text="(none)", font=FONT_SMALL, bg=BG_CARD, fg=FG_DIM)
        self._lbl_summary_empty.pack(anchor="w", pady=(2, 0))

        tk.Frame(self._summary_outer, bg=BORDER, height=1).pack(fill="x")

        summary_scroll_outer = tk.Frame(self._summary_outer, bg=BG_CARD)
        summary_scroll_outer.pack(fill="both", expand=True)

        self._summary_canvas = tk.Canvas(summary_scroll_outer, bg=BG_CARD, highlightthickness=0, bd=0)
        summary_vsb = ttk.Scrollbar(summary_scroll_outer, orient="vertical", command=self._summary_canvas.yview)
        self._summary_canvas.configure(yscrollcommand=summary_vsb.set)
        summary_vsb.pack(side="right", fill="y")
        self._summary_canvas.pack(side="left", fill="both", expand=True)

        self._summary_inner = tk.Frame(self._summary_canvas, bg=BG_CARD)
        self._summary_canvas_win = self._summary_canvas.create_window((0, 0), window=self._summary_inner, anchor="nw")
        self._summary_inner.bind("<Configure>", self._on_summary_inner_configure)
        self._summary_canvas.bind("<Configure>", self._on_summary_canvas_configure)

    def _build_legend_selection(self):
        for w in self._legend.winfo_children():
            w.destroy()
        for color, label in [
            (ACCENT, "Minimal"),
            (ACCENT_OPT, "Optional"),
            (ACCENT_REJ, "Neutral"),
            (ACCENT_ALWAYS, "Loosest (no bit set)"),
        ]:
            dot = tk.Label(self._legend, text="●", font=FONT_BODY, bg=BG, fg=color)
            dot.pack(side="left")
            tk.Label(self._legend, text=label, font=FONT_SMALL, bg=BG, fg=FG_DIM).pack(side="left", padx=(2, 12))

    def _build_legend_filter(self):
        for w in self._legend.winfo_children():
            w.destroy()
        dot = tk.Label(self._legend, text="●", font=FONT_BODY, bg=BG, fg=ACCENT_FILTER)
        dot.pack(side="left")
        tk.Label(
            self._legend,
            text="Filter value (fixed, already applied — not selectable)",
            font=FONT_SMALL,
            bg=BG,
            fg=FG_DIM,
        ).pack(side="left", padx=(2, 12))

    # ── Combobox styling ──────────────────────────────────────────────────────
    def _style_combobox(self):
        s = ttk.Style(self)
        s.theme_use("clam")
        s.configure(
            "TCombobox",
            fieldbackground=BG_CARD,
            background=BG_CARD,
            foreground=FG,
            selectbackground=SEL_BG,
            selectforeground=FG,
            bordercolor=BORDER,
            arrowcolor=ACCENT,
        )
        s.map("TCombobox", fieldbackground=[("readonly", BG_CARD)])
        s.configure("TScrollbar", troughcolor=BG, background=BORDER)

    # ── Canvas resize helpers ─────────────────────────────────────────────────
    def _on_inner_configure(self, _e=None):
        self._canvas.configure(scrollregion=self._canvas.bbox("all"))

    def _on_canvas_configure(self, e):
        self._canvas.itemconfig(self._canvas_win, width=e.width)

    def _on_summary_inner_configure(self, _e=None):
        self._summary_canvas.configure(scrollregion=self._summary_canvas.bbox("all"))

    def _on_summary_canvas_configure(self, e):
        self._summary_canvas.itemconfig(self._summary_canvas_win, width=e.width)

    def _on_mousewheel(self, e):
        if e.num == 4:
            self._canvas.yview_scroll(-1, "units")
        elif e.num == 5:
            self._canvas.yview_scroll(1, "units")
        else:
            self._canvas.yview_scroll(int(-1 * (e.delta / 120)), "units")

    # ── File / histogram loading ──────────────────────────────────────────────
    def _open_file_dialog(self):
        path = filedialog.askopenfilename(
            title="Open ROOT file", filetypes=[("ROOT files", "*.root"), ("All files", "*.*")]
        )
        if path:
            self._load_file(path)

    def _load_file(self, path):
        if ROOT is None:
            messagebox.showerror("Missing dependency", "PyROOT is not available. Please install ROOT.")
            return
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            messagebox.showerror("Error", f"Cannot open ROOT file:\n{path}")
            return
        self._root_file = f
        # Keep the attribute in sync with what is actually loaded — it was
        # previously only ever set from the CLI argument, so a file opened via
        # the dialog left it pointing at the old (or None) path. Resolved to an
        # absolute path so a relative CLI argument (./AnalysisResults.root) still
        # gives a reproducible provenance record in the exported JSON.
        self._rootfile_path = resolve_file_path(path)
        self._lbl_file.config(text=self._rootfile_path)

        d = f.Get(self._tdir_path)
        if not d:
            messagebox.showerror("Error", f"Directory '{self._tdir_path}' not found in file.")
            return

        histograms = [k.GetName() for k in d.GetListOfKeys() if k.ReadObj().InheritsFrom("TH1")]
        if not histograms:
            messagebox.showwarning("Warning", "No TH1 histograms found in directory.")
            return

        self._hist_combo.config(values=histograms, state="readonly")
        self._hist_combo.current(0)
        self._on_hist_selected()

    def _on_hist_selected(self, _e=None):
        if self._root_file is None:
            return
        d = self._root_file.Get(self._tdir_path)
        if not d:
            return
        hname = self._hist_var.get()
        hist = d.Get(hname)
        if not hist:
            return
        self._hist = hist
        self._vars = {}
        self._check_labels = {}

        # Always start from a clean summary panel — forget()-ing the pane only
        # hides it, it does not destroy previously built rows.
        self._clear_summary()
        self._lbl_var_count.config(text="")

        parsed, is_filter = load_bins_from_hist(hist)
        self._is_filter_hist = is_filter

        if is_filter:
            # load_bins_from_hist returns a flat list of bin dicts when is_filter
            # is True. This check narrows `parsed` from the dict|list union the
            # type checker infers from the function's two return shapes (same
            # purpose an assert would serve, but unlike assert this is never
            # stripped out under python -O, so it stays a real runtime guard).
            if not isinstance(parsed, list):
                raise TypeError("load_bins_from_hist() must return a list of bin dicts when is_filter is True")
            self._filter_bins = parsed
            self._groups = {}
            self._build_legend_filter()
            self._build_filter_view()
            self._set_bitmask_bar_visible(False)
            self._set_summary_visible(False)
        else:
            # load_bins_from_hist returns a dict of SelectionName -> bins when
            # is_filter is False; same narrowing purpose as above.
            if not isinstance(parsed, dict):
                raise TypeError("load_bins_from_hist() must return a dict of bins when is_filter is False")
            self._groups = parsed
            self._filter_bins = []
            self._build_legend_selection()
            self._build_selections()
            self._set_bitmask_bar_visible(True)
            self._set_summary_visible(True)
            self._update_bitmask()

    def _clear_summary(self):
        """Destroy all summary rows and reset to the empty state. Called on every
        histogram switch so stale rows from a previous selection histogram never
        linger while a filter histogram (which never calls _rebuild_summary) is
        being shown."""
        for w in self._summary_inner.winfo_children():
            w.destroy()
        self._lbl_summary_empty.config(text="(none)")
        self._on_summary_inner_configure()

    def _set_bitmask_bar_visible(self, visible):
        # self._bottom is a plain pack()ed frame (never added to the PanedWindow),
        # so toggling it with pack()/pack_forget() is fine here.
        if visible:
            self._bottom.pack(fill="x", side="bottom")
        else:
            self._bottom.pack_forget()

    def _set_summary_visible(self, visible):
        # self._summary_outer IS managed by the PanedWindow, so it must be
        # shown/hidden via the PanedWindow's own add()/forget() API — calling
        # pack() on it directly conflicts with the PanedWindow geometry manager
        # and produces broken layout (this was the bug: the summary pane would
        # swallow the whole window and the selection cards would disappear).
        panes = self._paned.panes()
        summary_path = str(self._summary_outer)
        if visible:
            if summary_path not in panes:
                self._paned.add(self._summary_outer, stretch="never", minsize=280, width=320)
        else:
            if summary_path in panes:
                self._paned.forget(self._summary_outer)

    # ── Filter (read-only) view ───────────────────────────────────────────────
    def _build_filter_view(self):
        for w in self._inner.winfo_children():
            w.destroy()

        card = tk.Frame(self._inner, bg=BG_CARD, bd=0, highlightthickness=1, highlightbackground=BORDER)
        card.pack(fill="x", padx=10, pady=5, ipadx=8, ipady=6)

        hdr = tk.Frame(card, bg=BG_CARD)
        hdr.pack(fill="x", padx=6, pady=(4, 2))
        tk.Label(hdr, text="Filter values", font=FONT_HEAD, bg=BG_CARD, fg=FG).pack(side="left")
        tk.Label(
            hdr,
            text="fixed values already applied when the producer ran — nothing to select",
            font=FONT_SMALL,
            bg=BG_CARD,
            fg=FG_DIM,
        ).pack(side="left", padx=10)

        tk.Frame(card, bg=BORDER, height=1).pack(fill="x", padx=6, pady=2)

        rows_frame = tk.Frame(card, bg=BG_CARD)
        rows_frame.pack(fill="x", padx=6, pady=4)

        for b in self._filter_bins:
            name = b.get("FilterName", "unknown")
            value = b.get("FilterValue", "")

            row = tk.Frame(rows_frame, bg=BG_CARD)
            row.pack(fill="x", pady=1)

            tk.Label(row, text="●", font=FONT_BODY, bg=BG_CARD, fg=ACCENT_FILTER).pack(side="left", padx=(0, 4))
            tk.Label(row, text=name, font=FONT_BODY, bg=BG_CARD, fg=FG, anchor="w").pack(side="left")
            tk.Label(row, text=value, font=FONT_BODY, bg=BG_CARD, fg=FG_DIM, anchor="e").pack(side="right", padx=4)

        self._on_inner_configure()

    # ── Selection cards ───────────────────────────────────────────────────────
    def _build_selections(self):
        for w in self._inner.winfo_children():
            w.destroy()

        for sel_name, group in self._groups.items():
            self._build_group_card(sel_name, group)

        self._on_inner_configure()

    def _build_group_card(self, sel_name, group):
        # categorise into selectable bins
        minimal = [(i, b) for i, b in enumerate(group) if bin_type(b) == "minimal"]
        optional = [(i, b) for i, b in enumerate(group) if bin_type(b) == "optional"]
        neutral = [(i, b) for i, b in enumerate(group) if bin_type(b) == "neutral"]

        # detect loosest/always-true: minimal bins exist but all have BitPosition==X
        loosest_only = has_only_skipped_minimal_bins(group)

        # skip entirely if there is truly nothing to show
        if not (minimal or optional or neutral or loosest_only):
            return

        card = tk.Frame(self._inner, bg=BG_CARD, bd=0, highlightthickness=1, highlightbackground=BORDER)
        card.pack(fill="x", padx=10, pady=5, ipadx=8, ipady=6)

        # header row
        hdr = tk.Frame(card, bg=BG_CARD)
        hdr.pack(fill="x", padx=6, pady=(4, 2))
        tk.Label(hdr, text=sel_name, font=FONT_HEAD, bg=BG_CARD, fg=FG).pack(side="left")

        if loosest_only:
            # All minimal bins are BitPosition==X → loosest cut, no bit set
            tk.Label(
                hdr,
                text="loosest selection — no bit set",
                font=FONT_SMALL,
                bg=BG_CARD,
                fg=ACCENT_ALWAYS,
            ).pack(side="left", padx=10)
        else:
            # show the loosest (most permissive) minimal threshold as a hint
            skipped_minimal = get_skipped_minimal_bins(group)
            if skipped_minimal:
                loosest_val = format_value_with_comment(skipped_minimal[0])
                tk.Label(
                    hdr,
                    text=f"minimal cut  →  always applied floor: {loosest_val}",
                    font=FONT_SMALL,
                    bg=BG_CARD,
                    fg=FG_DIM,
                ).pack(side="left", padx=10)
            elif minimal:
                loosest_val = format_value_with_comment(minimal[0][1])
                tk.Label(
                    hdr,
                    text=f"minimal cut  →  loosest selection: {loosest_val}",
                    font=FONT_SMALL,
                    bg=BG_CARD,
                    fg=FG_DIM,
                ).pack(side="left", padx=10)
            elif optional:
                tk.Label(hdr, text="optional", font=FONT_SMALL, bg=BG_CARD, fg=ACCENT_OPT).pack(side="left", padx=10)
            elif neutral:
                tk.Label(hdr, text="neutral", font=FONT_SMALL, bg=BG_CARD, fg=ACCENT_REJ).pack(side="left", padx=10)

        # separator
        tk.Frame(card, bg=BORDER, height=1).pack(fill="x", padx=6, pady=2)

        # bins
        bins_frame = tk.Frame(card, bg=BG_CARD)
        bins_frame.pack(fill="x", padx=6, pady=4)

        if loosest_only:
            # The whole cut is a single always-true rung — there is no stricter
            # alternative to vary against, so nothing to check: informational only.
            for b in group:
                if b.get("MinimalCut", "0") == "1" and b.get("OptionalCut", "0") == "0":
                    self._build_loosest_row(bins_frame, None, None, b)
        else:
            # The loosest bin is always in effect when no stricter minimal bit is
            # set — which makes it the loosest rung of the ladder and therefore a
            # variation point in its own right. It gets a checkbox like any other
            # rung; it just contributes no bit when picked.
            for i, b in get_skipped_minimal_bin_items(group):
                self._build_loosest_row(bins_frame, sel_name, i, b)
            for i, b in minimal:
                self._build_bin_row(bins_frame, sel_name, i, b, "minimal")
            for i, b in optional:
                self._build_bin_row(bins_frame, sel_name, i, b, "optional")
            for i, b in neutral:
                self._build_bin_row(bins_frame, sel_name, i, b, "neutral")

    def _build_loosest_row(self, parent, sel_name, idx, b):
        """Row for a loosest/always-true cut — never sets a bit.

        With sel_name/idx given the row is checkable, so the loosest rung can be
        included as a systematic variation; with both None it is purely
        informational (used when the cut has no stricter rung to vary against).

        Either way the leading column is the same widget type/width (a Label with
        width=CHECK_COL_WIDTH) as _build_bin_row's checkbox glyph, so all row
        types line up pixel-for-pixel. A native tk.Checkbutton's indicator box has
        no queryable pixel width, so matching it with a plain spacer Label never
        aligns reliably — using an identical Label everywhere is the only way to
        guarantee it.
        """
        label_text = format_value_with_comment(b)
        checkable = sel_name is not None and idx is not None

        row = tk.Frame(parent, bg=BG_CARD, cursor="hand2" if checkable else "")
        row.pack(fill="x", pady=1)

        tk.Label(row, text="●", font=FONT_BODY, bg=BG_CARD, fg=ACCENT_ALWAYS).pack(side="left", padx=(0, 4))

        check_lbl = tk.Label(
            row,
            text="[ ]" if checkable else " ",
            font=FONT_BODY,
            bg=BG_CARD,
            fg=FG_DIM,
            width=CHECK_COL_WIDTH,
            anchor="w",
        )
        check_lbl.pack(side="left")

        text_lbl = tk.Label(row, text=label_text, font=FONT_BODY, bg=BG_CARD, fg=FG_DIM, anchor="w")
        text_lbl.pack(side="left", fill="x", expand=True)

        tk.Label(row, text="no bit", font=FONT_SMALL, bg=BG_CARD, fg=ACCENT_ALWAYS, width=8).pack(side="right", padx=4)

        if not checkable:
            return

        var = tk.BooleanVar(value=False)
        self._vars[(sel_name, idx)] = var
        self._check_labels[(sel_name, idx)] = check_lbl

        def toggle(_e=None):
            var.set(not var.get())
            check_lbl.config(text="[x]" if var.get() else "[ ]", fg=ACCENT_ALWAYS if var.get() else FG_DIM)
            text_lbl.config(fg=FG if var.get() else FG_DIM)
            self._update_bitmask()

        for w in (row, check_lbl, text_lbl):
            w.bind("<Button-1>", toggle)

    def _build_bin_row(self, parent, sel_name, idx, b, kind):
        color = {"minimal": ACCENT, "optional": ACCENT_OPT, "neutral": ACCENT_REJ}[kind]
        label_text = format_value_with_comment(b)

        var = tk.BooleanVar(value=False)
        self._vars[(sel_name, idx)] = var

        row = tk.Frame(parent, bg=BG_CARD, cursor="hand2")
        row.pack(fill="x", pady=1)

        tk.Label(row, text="●", font=FONT_BODY, bg=BG_CARD, fg=color).pack(side="left", padx=(0, 4))

        # Custom checkbox glyph, built from the SAME Label(width=CHECK_COL_WIDTH)
        # pattern as the loosest-row spacer above — this is what actually fixes
        # the alignment: both row types share one widget recipe for their
        # leading column instead of trying to match a native Checkbutton's
        # unmeasurable indicator box.
        check_lbl = tk.Label(row, text="[ ]", font=FONT_BODY, bg=BG_CARD, fg=FG_DIM, width=CHECK_COL_WIDTH, anchor="w")
        check_lbl.pack(side="left")
        self._check_labels[(sel_name, idx)] = check_lbl

        text_lbl = tk.Label(row, text=label_text, font=FONT_BODY, bg=BG_CARD, fg=FG, anchor="w")
        text_lbl.pack(side="left", fill="x", expand=True)

        # bit-position badge
        pos = b.get("BitPosition", "X")
        if pos.upper() != "X":
            tk.Label(row, text=f"bit {pos}", font=FONT_SMALL, bg=BG_CARD, fg=FG_DIM, width=8).pack(side="right", padx=4)

        def toggle(_e=None):
            var.set(not var.get())
            check_lbl.config(text="[x]" if var.get() else "[ ]", fg=color if var.get() else FG_DIM)
            self._update_bitmask()

        for w in (row, check_lbl, text_lbl):
            w.bind("<Button-1>", toggle)

    # ── Bitmask computation + summary update ──────────────────────────────────
    def _update_bitmask(self):
        if self._is_filter_hist:
            return

        bitmask = 0
        checked = []  # (sel_name, kind, pos_int, label_text, pos_str, color, is_loosest)

        for (sel_name, idx), var in self._vars.items():
            if not var.get():
                continue
            b = self._groups[sel_name][idx]
            pos = b.get("BitPosition", "X")
            # a checked loosest rung sets no bit; variation_kind still files it
            # under its selection's minimal axis
            is_loosest = pos.upper() == "X"
            kind = variation_kind(b)
            color = (
                ACCENT_ALWAYS
                if is_loosest
                else {"minimal": ACCENT, "optional": ACCENT_OPT, "neutral": ACCENT_REJ}.get(kind, FG_DIM)
            )

            if not is_loosest:
                bitmask |= 1 << int(pos)

            label_text = format_value_with_comment(b)
            checked.append((sel_name, kind, bit_position_int(b), label_text, pos, color, is_loosest))

        # Within the same (sel_name, kind) group, multiple checked bits represent a
        # threshold ladder for one observable, not independent alternatives — only
        # the strictest one (highest BitPosition) is the meaningful cut, so collapse
        # to that single entry for display. The loosest rung sorts to -1 and so
        # loses to any bit-carrying rung of the same ladder, which is correct: the
        # bitmask bar shows one working point, not the scan.
        tightest_by_group = {}
        for sel_name, kind, pos_int, label_text, pos, color, is_loosest in checked:
            key = (sel_name, kind)
            if key not in tightest_by_group or pos_int > tightest_by_group[key][0]:
                tightest_by_group[key] = (pos_int, label_text, pos, color, is_loosest)

        selected_entries = [
            (sel_name, label_text, pos, color, is_loosest)
            for (sel_name, kind), (pos_int, label_text, pos, color, is_loosest) in tightest_by_group.items()
        ]

        # Selection names that already have a checked "minimal" entry — for those,
        # the checked threshold is strictly stricter than the always-applied floor
        # (that's what makes it a selectable, higher-bit-position option in the
        # first place), so the floor is logically subsumed and showing it alongside
        # the checked value is redundant. Only surface the floor when nothing
        # stricter has been chosen for that selection.
        sel_names_with_checked_minimal = {
            sel_name for (sel_name, kind) in tightest_by_group.keys() if kind == "minimal"
        }

        # Loosest/always-applied bins: shown only when no stricter minimal cut for
        # the same selection is currently checked — otherwise this entry is fully
        # dominated by the checked one and would just duplicate the same cut.
        for sel_name, group in self._groups.items():
            if sel_name in sel_names_with_checked_minimal:
                continue
            for b in get_skipped_minimal_bins(group):
                label_text = format_value_with_comment(b)
                selected_entries.append((sel_name, label_text, "X", ACCENT_ALWAYS, True))

        self._lbl_dec.config(text=str(bitmask))
        self._lbl_hex.config(text=hex(bitmask))
        self._lbl_bin.config(text=bin(bitmask))

        self._rebuild_summary(selected_entries)
        self._update_variation_count()

    def _rebuild_summary(self, entries):
        """Rebuild the selected-cuts summary sidebar. Narrow column → each value
        gets its own row, stacked under the selection name, rather than packed
        side-by-side in a single line."""
        self._clear_summary()

        if not entries:
            return

        self._lbl_summary_empty.config(text="")

        # group by selection name, preserving is_loosest flag
        by_name = {}
        for sel_name, label_text, pos, color, is_loosest in entries:
            by_name.setdefault(sel_name, []).append((label_text, pos, color, is_loosest))

        for sel_name, items in by_name.items():
            block = tk.Frame(self._summary_inner, bg=BG_CARD)
            block.pack(fill="x", padx=10, pady=(4, 2))

            tk.Label(
                block,
                text=sel_name,
                font=FONT_BODY,
                bg=BG_CARD,
                fg=FG,
                anchor="w",
                wraplength=280,
                justify="left",
            ).pack(fill="x")

            for label_text, pos, color, is_loosest in items:
                entry_row = tk.Frame(block, bg=BG_CARD)
                entry_row.pack(fill="x", padx=(14, 0), pady=1)

                tk.Label(entry_row, text="●", font=FONT_SMALL, bg=BG_CARD, fg=color).pack(side="left")
                tk.Label(
                    entry_row,
                    text=label_text,
                    font=FONT_SMALL,
                    bg=BG_CARD,
                    fg=FG_DIM if is_loosest else FG,
                    anchor="w",
                    wraplength=220,
                    justify="left",
                ).pack(side="left", padx=(4, 0), fill="x", expand=True)

                if is_loosest:
                    tk.Label(
                        entry_row,
                        text="no bit",
                        font=FONT_SMALL,
                        bg=BG_CARD,
                        fg=ACCENT_ALWAYS,
                    ).pack(side="right")
                elif pos.upper() != "X":
                    tk.Label(
                        entry_row,
                        text=f"bit {pos}",
                        font=FONT_SMALL,
                        bg=BG_CARD,
                        fg=FG_DIM,
                    ).pack(side="right")

            tk.Frame(self._summary_inner, bg=BORDER, height=1).pack(fill="x", padx=10, pady=(4, 0))

        self._on_summary_inner_configure()

    # ── Systematic variations ─────────────────────────────────────────────────
    def _variation_axes(self):
        """Group the currently checked bins into one variation axis per
        (SelectionName, kind).

        Semantics: within one axis the checked bins are treated as MUTUALLY
        EXCLUSIVE alternatives — each one is a separate systematic variation of
        the same observable, not an additional cut applied on top. This holds for
        minimal and optional cuts alike, and includes a checked loosest rung
        (BitPosition=="X"), which is just the variation that sets no bit.

        Returns an ordered dict {(sel_name, kind): [bin_dict, ...]}, each list
        sorted by ascending BitPosition (loosest → strictest). Insertion order
        follows the card order on screen, so the enumeration is reproducible.
        """
        axes: Dict[Any, List[Dict[str, Any]]] = {}
        for (sel_name, idx), var in self._vars.items():
            if not var.get():
                continue
            b = self._groups[sel_name][idx]
            axes.setdefault((sel_name, variation_kind(b)), []).append(b)

        for bins in axes.values():
            bins.sort(key=bit_position_int)
        return axes

    @staticmethod
    def _axis_key_labels(axes):
        """JSON key for each axis: the plain SelectionName, disambiguated with the
        cut kind only when the same SelectionName carries more than one axis."""
        counts: Dict[str, int] = {}
        for sel_name, _kind in axes:
            counts[sel_name] = counts.get(sel_name, 0) + 1
        return {
            (sel_name, kind): (sel_name if counts[sel_name] == 1 else f"{sel_name} [{kind}]") for sel_name, kind in axes
        }

    def _always_applied_bins(self, axes):
        """The BitPosition==X minimal floors that are in effect for every variation.

        Skipped for any SelectionName that has a checked minimal axis: either the
        picked rung sits at a higher bit position and is therefore strictly
        tighter (floor subsumed), or the floor itself is the picked rung and is
        already carried by the axis. Emitting it here too would duplicate the key.
        """
        sel_with_minimal_axis = {sel_name for sel_name, kind in axes if kind == "minimal"}
        out = []
        for sel_name, group in self._groups.items():
            if sel_name in sel_with_minimal_axis:
                continue
            for b in get_skipped_minimal_bins(group):
                out.append((sel_name, b))
        return out

    def _variation_count(self, axes=None):
        axes = self._variation_axes() if axes is None else axes
        if not axes:
            return 0
        n = 1
        for bins in axes.values():
            n *= len(bins)
        return n

    def _update_variation_count(self):
        n = self._variation_count()
        if n <= 1:
            # 0 = nothing checked, 1 = a single fixed working point (no variation)
            self._lbl_var_count.config(text="")
        else:
            self._lbl_var_count.config(text=f"→ {n} variations")

    def _compute_variations(self):
        """Cartesian product over all variation axes → one dict per combination."""
        axes = self._variation_axes()
        if not axes:
            return []

        labels = self._axis_key_labels(axes)
        axis_items = list(axes.items())
        floors = self._always_applied_bins(axes)

        variations = []
        for combo in itertools.product(*[bins for _key, bins in axis_items]):
            bitmask = 0
            entry: Dict[str, Any] = {}
            for (key, _bins), b in zip(axis_items, combo):
                pos = b.get("BitPosition", "X")
                if pos.upper() != "X":
                    bitmask |= 1 << int(pos)
                entry[unique_key(entry, labels[key])] = b.get("Value", "")

            # floors carry no bit but define the cut actually in effect, so they
            # belong in the record of each variation
            for sel_name, b in floors:
                base = sel_name if sel_name not in entry else f"{sel_name} [floor]"
                entry[unique_key(entry, base)] = b.get("Value", "")

            variations.append(
                {
                    "index": len(variations),
                    "bitmask": bitmask,
                    "bitmask_hex": hex(bitmask),
                    **entry,
                }
            )
        return variations

    def _export_variations(self):
        if self._is_filter_hist:
            messagebox.showinfo(
                "Not applicable",
                "This is a filter histogram — its values are fixed and carry no bits, so there is nothing to vary.",
            )
            return

        axes = self._variation_axes()
        if not axes:
            messagebox.showwarning(
                "Nothing selected",
                "Check the options you want to scan first.\n\n"
                "Options checked within the same cut are enumerated as alternative "
                "variations; the export is the product across all cuts.",
            )
            return

        n = self._variation_count(axes)
        if n > VARIATION_WARN_LIMIT and not messagebox.askyesno(
            "Many variations", f"This will enumerate {n} variations. Continue?"
        ):
            return

        variations = self._compute_variations()

        path = filedialog.asksaveasfilename(
            title="Save systematic variations",
            defaultextension=".json",
            initialfile="cut_variations.json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return

        payload = {
            "file": self._rootfile_path,
            "directory": self._tdir_path,
            "histogram": self._hist_var.get(),
            "n_variations": len(variations),
            "variations": variations,
        }

        try:
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2)
                fh.write("\n")
        except OSError as exc:
            messagebox.showerror("Error", f"Could not write file:\n{exc}")
            return

        messagebox.showinfo("Exported", f"Wrote {len(variations)} variations to\n{path}")

    # ── Utilities ─────────────────────────────────────────────────────────────
    def _copy(self, text):
        self.clipboard_clear()
        self.clipboard_append(text)
        self.update()

    @staticmethod
    def _make_button(parent, text, cmd, color):
        return tk.Button(
            parent,
            text=text,
            command=cmd,
            font=FONT_BODY,
            bg=BG_CARD,
            fg=color,
            activebackground=BG_HOVER,
            activeforeground=color,
            relief="flat",
            bd=0,
            padx=10,
            pady=4,
            highlightthickness=1,
            highlightbackground=color,
            cursor="hand2",
        )


# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CutCulator GUI")
    parser.add_argument("rootfile", nargs="?", default=None, help="Path to ROOT file (optional, can be opened via GUI)")
    parser.add_argument("--dir", default="femto-producer", help="TDirectory path in ROOT file")
    args = parser.parse_args()

    app = CutCulatorApp(rootfile=args.rootfile, tdir=args.dir)
    app.mainloop()
