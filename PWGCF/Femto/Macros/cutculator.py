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
@brief  CutCulator (Compute bitmask for selecting particles in the Femto Framework)
@author Anton Riedel <anton.riedel@cern.ch>, Technical University of Munich
"""

import ROOT
import argparse

VALUE_DELIM = "___"
SECTION_DELIM = ":::"


def parse_bin_label(label):
    """Parse a bin label into a dictionary."""
    result = {}
    sections = label.split(SECTION_DELIM)
    for sec in sections:
        if VALUE_DELIM not in sec:
            continue
        key, value = sec.split(VALUE_DELIM, 1)
        result[key] = value
    return result


def format_value_with_comment(b):
    """Return Value plus optional (comment=...) suffix."""
    val = b.get("Value", "")
    comment = b.get("Comment", "")
    if comment and comment.upper() != "X":
        return f"{val} (comment={comment})"
    return val


def ask_user_selection(group):
    """
    Prompt user to select bin(s) for this selection group.
    - If minimal selections contain exactly 1 entry → auto-select it.
    - Optional selections remain user-selectable.
    """
    selection_name = group[0].get("SelectionName", "unknown")

    # Separate minimal and optional bins
    minimal_bins = [b for b in group if b.get("MinimalCut", "0") == "1" and b.get("OptionalCut", "0") == "0"]
    optional_bins = [b for b in group if b.get("OptionalCut", "0") == "1"]

    selected_bins = []

    # ----- Minimal selection -----
    if minimal_bins:
        if len(minimal_bins) == 1:
            only = minimal_bins[0]
            print(
                f"\nSelection: {selection_name} — only one minimal option → auto-selecting: "
                f"{format_value_with_comment(only)}"
            )
            selected_bins.append(only)
        else:
            print(f"\nSelection: {selection_name}")
            for idx, b in enumerate(minimal_bins):
                print(f"  [{idx}] {format_value_with_comment(b)}")
            while True:
                sel_input = input("Enter index for minimal cut (0 = loosest minimal): ")
                if sel_input.strip() == "":
                    sel_input = "0"
                try:
                    sel_idx = int(sel_input)
                    if 0 <= sel_idx < len(minimal_bins):
                        choice = minimal_bins[sel_idx]
                        selected_bins.append(choice)
                        print(f"Selected: {format_value_with_comment(choice)}")
                        break
                except ValueError:
                    pass
                print("Invalid input. Please enter a valid index.")

    # ----- Optional selection -----
    if optional_bins:
        print(f"\nSelection: {selection_name} (optional selection, 0 to skip)")
        for idx, b in enumerate(optional_bins, start=1):
            print(f"  [{idx}] {format_value_with_comment(b)}")

        while True:
            sel_input = input("Enter indices separated by space (0 to skip): ")
            if not sel_input.strip() or sel_input.strip() == "0":
                print("Selected: (skipped)")
                break

            try:
                indices = [int(x) for x in sel_input.split()]
                if all(0 <= i <= len(optional_bins) for i in indices):
                    chosen = []
                    for i in indices:
                        if i != 0:
                            b = optional_bins[i - 1]
                            selected_bins.append(b)
                            chosen.append(format_value_with_comment(b))

                    print("Selected: " + ", ".join(chosen))
                    break
            except ValueError:
                pass

            print("Invalid input. Please enter valid indices separated by space.")

    return selected_bins


def main(rootfile_path, tdir_path="femto-producer"):
    print(f"Opening ROOT file: {rootfile_path}")
    f = ROOT.TFile.Open(rootfile_path)
    if not f:
        print("Cannot open ROOT file")
        return

    print(f"Accessing directory: {tdir_path}")
    d = f.Get(tdir_path)
    if not d:
        print(f"Cannot access directory {tdir_path}")
        return

    histograms = [k.GetName() for k in d.GetListOfKeys() if k.ReadObj().InheritsFrom("TH1")]
    if not histograms:
        print("No histograms found")
        return

    print("\nHistograms found in directory:")
    for i, hname in enumerate(histograms):
        print(f"  [{i}] {hname}")
    hidx = int(input("\nSelect histogram index: "))
    hname = histograms[hidx]
    hist = d.Get(hname)
    nbins = hist.GetNbinsX()
    print(f"\nUsing histogram: {hname}")
    print(f"Histogram contains {nbins} bins.\n")

    # parse all bins, ignoring the last 2 special bins
    bins = []
    for i in range(1, nbins - 2 + 1):
        label = hist.GetXaxis().GetBinLabel(i)
        if not label:
            continue
        bdict = parse_bin_label(label)
        bdict["_bin_index"] = i
        bins.append(bdict)

    # group by SelectionName
    groups = {}
    for b in bins:
        sel_name = b.get("SelectionName", f"unknown_{b['_bin_index']}")
        groups.setdefault(sel_name, []).append(b)

    selected_bins = []

    for group in groups.values():
        res = ask_user_selection(group)
        if res:
            selected_bins.extend(res)

    # compute bitmask from selected bins
    bitmask = 0
    for b in selected_bins:
        pos = b.get("BitPosition", "")
        if pos.upper() == "X":
            continue
        bitmask |= 1 << int(pos)

    print("\n=======================================")
    print("Summary of your selections:")
    print("=======================================\n")

    summary = {}
    for b in selected_bins:
        sel = b.get("SelectionName", "unknown")
        summary.setdefault(sel, []).append(format_value_with_comment(b))

    for sel, values in summary.items():
        print(f"  {sel}: {', '.join(values)}")

    print("\nFinal selected bitmask:")
    print(f"  Decimal: {bitmask}")
    print(f"  Binary:  {bin(bitmask)}")
    print(f"  Hex:     {hex(bitmask)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", help="Path to ROOT file")
    parser.add_argument("--dir", default="femto-producer", help="TDirectory path in ROOT file")
    args = parser.parse_args()
    main(args.rootfile, args.dir)
