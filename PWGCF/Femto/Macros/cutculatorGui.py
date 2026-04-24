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
import sys

try:
    import ROOT

    ROOT.gROOT.SetBatch(True)
except ImportError:
    ROOT = None

VALUE_DELIM = "___"
SECTION_DELIM = ":::"

# ── Colours ──────────────────────────────────────────────────────────────────
BG = "#1e1e2e"
BG_CARD = "#2a2a3e"
BG_HOVER = "#3a3a50"
ACCENT = "#89b4fa"  # blue  – minimal
ACCENT_OPT = "#a6e3a1"  # green – optional
ACCENT_REJ = "#f38ba8"  # red   – rejection / neutral
FG = "#cdd6f4"
FG_DIM = "#6c7086"
BORDER = "#45475a"
SEL_BG = "#313244"

FONT_TITLE = ("Inter", 15, "bold")
FONT_HEAD = ("Inter", 11, "bold")
FONT_BODY = ("Inter", 10)
FONT_MONO = ("JetBrains Mono", 11, "bold")
FONT_SMALL = ("Inter", 9)


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
    # group by SelectionName preserving order
    groups = {}
    for b in bins:
        name = b.get("SelectionName", f"unknown_{b['_bin_index']}")
        groups.setdefault(name, []).append(b)
    return groups


# ── Main Application ──────────────────────────────────────────────────────────
class CutCulatorApp(tk.Tk):
    def __init__(self, rootfile=None, tdir="femto-producer"):
        super().__init__()
        self.title("CutCulator")
        self.configure(bg=BG)
        self.geometry("920x720")
        self.minsize(780, 560)
        self.resizable(True, True)

        self._rootfile_path = rootfile
        self._tdir_path = tdir
        self._root_file = None
        self._hist = None
        self._groups = {}  # SelectionName → list[bin_dict]
        self._vars = {}  # (SelectionName, idx) → BooleanVar

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
        legend = tk.Frame(self, bg=BG, pady=4, padx=18)
        legend.pack(fill="x")
        for color, label in [(ACCENT, "Minimal"), (ACCENT_OPT, "Optional"), (ACCENT_REJ, "Neutral")]:
            dot = tk.Label(legend, text="●", font=FONT_BODY, bg=BG, fg=color)
            dot.pack(side="left")
            tk.Label(legend, text=label, font=FONT_SMALL, bg=BG, fg=FG_DIM).pack(side="left", padx=(2, 12))

        # ── scrollable selection area ──
        outer = tk.Frame(self, bg=BG)
        outer.pack(fill="both", expand=True, padx=12, pady=4)

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

        # ── bottom result bar ──
        bottom = tk.Frame(self, bg=BG_CARD, pady=10, padx=18)
        bottom.pack(fill="x", side="bottom")

        tk.Label(bottom, text="Bitmask:", font=FONT_HEAD, bg=BG_CARD, fg=FG).pack(side="left")

        self._lbl_dec = tk.Label(bottom, text="—", font=FONT_MONO, bg=BG_CARD, fg=ACCENT, width=14, anchor="w")
        self._lbl_dec.pack(side="left", padx=8)

        self._lbl_hex = tk.Label(bottom, text="—", font=FONT_MONO, bg=BG_CARD, fg=ACCENT_OPT, width=14, anchor="w")
        self._lbl_hex.pack(side="left", padx=8)

        self._lbl_bin = tk.Label(bottom, text="—", font=FONT_MONO, bg=BG_CARD, fg=FG_DIM, anchor="w")
        self._lbl_bin.pack(side="left", padx=8)

        copy_frame = tk.Frame(bottom, bg=BG_CARD)
        copy_frame.pack(side="right")
        self._make_button(copy_frame, "Copy Dec", lambda: self._copy(self._lbl_dec["text"]), ACCENT).pack(
            side="left", padx=3
        )
        self._make_button(copy_frame, "Copy Hex", lambda: self._copy(self._lbl_hex["text"]), ACCENT_OPT).pack(
            side="left", padx=3
        )

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
        self._lbl_file.config(text=path)

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
        self._groups = load_bins_from_hist(hist)
        self._vars = {}
        self._build_selections()
        self._update_bitmask()

    # ── Selection cards ───────────────────────────────────────────────────────
    def _build_selections(self):
        for w in self._inner.winfo_children():
            w.destroy()

        for sel_name, group in self._groups.items():
            self._build_group_card(sel_name, group)

        self._on_inner_configure()

    def _build_group_card(self, sel_name, group):
        # categorise
        minimal = [(i, b) for i, b in enumerate(group) if bin_type(b) == "minimal"]
        optional = [(i, b) for i, b in enumerate(group) if bin_type(b) == "optional"]
        neutral = [(i, b) for i, b in enumerate(group) if bin_type(b) == "neutral"]

        if not (minimal or optional or neutral):
            return  # nothing to show (all "skip")

        card = tk.Frame(self._inner, bg=BG_CARD, bd=0, highlightthickness=1, highlightbackground=BORDER)
        card.pack(fill="x", padx=10, pady=5, ipadx=8, ipady=6)

        # header row
        hdr = tk.Frame(card, bg=BG_CARD)
        hdr.pack(fill="x", padx=6, pady=(4, 2))
        tk.Label(hdr, text=sel_name, font=FONT_HEAD, bg=BG_CARD, fg=FG).pack(side="left")

        # show the loosest (most permissive) minimal threshold as a hint.
        # the truly loosest threshold has mSkipMostPermissiveBit=true so its
        # BitPosition is "X" — it lands in the "skip" category. check there first,
        # then fall back to the loosest bit-carrying minimal bin.
        skipped_minimal = [
            b for b in group if b.get("MinimalCut", "0") == "1" and b.get("BitPosition", "X").upper() == "X"
        ]
        if skipped_minimal:
            loosest_val = format_value_with_comment(skipped_minimal[0])
            tk.Label(
                hdr, text=f"minimal cut  →  loosest selection: {loosest_val}", font=FONT_SMALL, bg=BG_CARD, fg=FG_DIM
            ).pack(side="left", padx=10)
        elif minimal:
            loosest_val = format_value_with_comment(minimal[0][1])
            tk.Label(
                hdr, text=f"minimal cut  →  loosest selection: {loosest_val}", font=FONT_SMALL, bg=BG_CARD, fg=FG_DIM
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

        for i, b in minimal:
            self._build_bin_row(bins_frame, sel_name, i, b, "minimal")
        for i, b in optional:
            self._build_bin_row(bins_frame, sel_name, i, b, "optional")
        for i, b in neutral:
            self._build_bin_row(bins_frame, sel_name, i, b, "neutral")

    def _build_bin_row(self, parent, sel_name, idx, b, kind):
        color = {"minimal": ACCENT, "optional": ACCENT_OPT, "neutral": ACCENT_REJ}[kind]
        label_text = format_value_with_comment(b)

        var = tk.BooleanVar(value=False)

        self._vars[(sel_name, idx)] = var

        row = tk.Frame(parent, bg=BG_CARD)
        row.pack(fill="x", pady=1)

        # coloured dot
        tk.Label(row, text="●", font=FONT_BODY, bg=BG_CARD, fg=color).pack(side="left", padx=(0, 4))

        # checkbox styled as a toggle button
        cb = tk.Checkbutton(
            row,
            text=label_text,
            variable=var,
            font=FONT_BODY,
            bg=BG_CARD,
            fg=FG,
            activebackground=BG_HOVER,
            activeforeground=FG,
            selectcolor=SEL_BG,
            relief="flat",
            bd=0,
            highlightthickness=0,
            cursor="hand2",
            command=self._update_bitmask,
        )
        cb.pack(side="left", fill="x", expand=True)

        # bit-position badge
        pos = b.get("BitPosition", "X")
        if pos.upper() != "X":
            tk.Label(row, text=f"bit {pos}", font=FONT_SMALL, bg=BG_CARD, fg=FG_DIM, width=8).pack(side="right", padx=4)

    # ── Bitmask computation ───────────────────────────────────────────────────
    def _update_bitmask(self):
        bitmask = 0
        for (sel_name, idx), var in self._vars.items():
            if not var.get():
                continue
            b = self._groups[sel_name][idx]
            pos = b.get("BitPosition", "X")
            if pos.upper() == "X":
                continue
            bitmask |= 1 << int(pos)

        self._lbl_dec.config(text=str(bitmask))
        self._lbl_hex.config(text=hex(bitmask))
        self._lbl_bin.config(text=bin(bitmask))

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
