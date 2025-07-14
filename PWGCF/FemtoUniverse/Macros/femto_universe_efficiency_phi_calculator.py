#!/usr/bin/env python3

# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

"""
A tool to calculate efficiency and upload it to CCDB.
Author: Dawid Karpiński (dawid.karpinski@cern.ch)
Modified by: Zuzanna Chochulska (zchochul@cern.ch)
"""

import argparse
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import ROOT  # pylint: disable=import-error


class CustomHelpFormatter(argparse.HelpFormatter):
    "Add default value to help format"

    def _get_help_string(self, action):
        help_str = action.help
        if help_str is not None and action.default not in [argparse.SUPPRESS, None]:
            help_str += f" (default: {action.default})"
        return help_str


parser = argparse.ArgumentParser(
    description="A tool to calculate efficiency and upload it to CCDB",
    formatter_class=CustomHelpFormatter,
)
parser.add_argument(
    "--alien-path",
    type=Path,
    help="path to train run's directory in Alien with analysis results "
    "[example: /alice/cern.ch/user/a/alihyperloop/outputs/0033/332611/70301]",
)
parser.add_argument(
    "--mc-reco",
    type=str,
    nargs="+",
    help="paths to MC Reco histograms, separated by space [example: task/mcreco_one/hPt task/mcreco_two/hPt]",
    required=True,
)
parser.add_argument(
    "--mc-truth",
    type=str,
    nargs="+",
    help="paths to MC Truth histograms, separated by space [example: task/mctruth_one/hPt task/mctruth_one/hPt]",
    required=True,
)
parser.add_argument(
    "--ccdb-path",
    type=str,
    help="location in CCDB to where objects will be uploaded",
    required=True,
)
parser.add_argument(
    "--ccdb-url",
    type=str,
    help="URL to CCDB",
    default="http://ccdb-test.cern.ch:8080",
)
parser.add_argument(
    "--ccdb-labels",
    type=str,
    nargs="+",
    help="custom labels to add to objects' metadata in CCDB [example: label1 label2]",
    default=[],
)
parser.add_argument(
    "--ccdb-lifetime",
    type=int,
    help="how long should objects in CCDB remain valid (milliseconds)",
    default=365 * 24 * 60 * 60 * 1000,  # one year
)
args = parser.parse_args()

if len(args.mc_reco) != len(args.mc_truth):
    print("[!] Provided number of histograms with MC Reco must match MC Truth", file=sys.stderr)
    sys.exit(1)

if len(args.ccdb_labels) > 0 and len(args.ccdb_labels) != len(args.mc_reco):
    print("[!] You must provide labels for all particles", file=sys.stderr)
    sys.exit(1)

if len(args.ccdb_labels) == 0:
    # if flag is not provided, fill with empty strings to match size
    args.ccdb_labels = [""] * len(args.mc_reco)

ANALYSIS_RESULTS = "AnalysisResults.root"
results_path = args.alien_path / ANALYSIS_RESULTS
job_id = results_path.parent.name
train_number = results_path.parent.parent.name

tmp_dir = Path(tempfile.gettempdir())

res_dest = tmp_dir / f"{train_number}-{job_id}-{ANALYSIS_RESULTS}"
eff_dest = tmp_dir / f"{train_number}-{job_id}-Efficiency.root"

# get file from alien
if not res_dest.is_file():
    print(f"[↓] Downloading analysis results from Alien to '{res_dest}' ...", file=sys.stderr)
    ROOT.TGrid.Connect("alien://")
    try:
        subprocess.run(
            ["alien_cp", results_path, "file://" + str(res_dest)],
            capture_output=True,
            check=True,
        )
        print("[-] Download complete!", file=sys.stderr)
    except subprocess.CalledProcessError as error:
        print(f"[!] Error while downloading results file: {error.stderr}", file=sys.stderr)
        sys.exit(1)
else:
    print(
        f"[-] Skipping download from Alien, since '{res_dest}' is already present",
        file=sys.stderr,
    )

print()

histos_to_upload = []

# get reco & truth histos
with (
    ROOT.TFile.Open(res_dest.as_uri()) as res_file,
    ROOT.TFile.Open(eff_dest.as_uri(), "recreate") as eff_file,
):

    for idx, (mc_reco, mc_truth) in enumerate(zip(args.mc_reco, args.mc_truth)):
        hist_reco = res_file.Get(mc_reco)
        if not hist_reco:
            print(f"[!] Cannot find MC Reco histogram in '{mc_reco}', aborting", file=sys.stderr)
            sys.exit(1)

        hist_truth = res_file.Get(mc_truth)
        if not hist_truth:
            print(f"[!] Cannot find MC Truth histogram in '{mc_truth}', aborting", file=sys.stderr)
            sys.exit(1)

        hist_reco.Rebin(4)
        hist_truth.Rebin(4)

        num_bins = hist_reco.GetNbinsX()
        x_max = hist_reco.GetXaxis().GetBinLowEdge(num_bins)
        x_max += hist_reco.GetXaxis().GetBinWidth(num_bins)

        hist_name = f"Efficiency_part{idx + 1}"

        # calculate efficiency
        eff = ROOT.TH1F(hist_name, "", num_bins, 0, x_max)
        for bin_idx in range(1, num_bins + 1):  # Bins start at 1 in ROOT
            denom = hist_truth.GetBinContent(bin_idx)
            if idx == 0:
                denom *= 0.489
            eff.SetBinContent(bin_idx, hist_reco.GetBinContent(bin_idx) / denom if denom > 0 else 0)

        # save efficiency object to file
        eff_file.WriteObject(eff, hist_name)
        histos_to_upload.append(hist_name)

if len(histos_to_upload) == 0:
    print("[-] Exiting, since there is nothing to upload", file=sys.stderr)
    sys.exit(1)

# upload objects to ccdb
try:
    for idx, key in enumerate(histos_to_upload):
        timestamp_start = int(time.time() * 1000)
        timestamp_end = timestamp_start + args.ccdb_lifetime

        print(f"[↑] Uploading {key} to {args.ccdb_url} ... ", file=sys.stderr, end="")
        upload_cmd = [
            "o2-ccdb-upload",
            "--file",
            eff_dest.as_uri(),
            "--host",
            args.ccdb_url,
            "--key",
            key,
            "--path",
            args.ccdb_path,
            "--starttimestamp",
            str(timestamp_start),
            "--endtimestamp",
            str(timestamp_end),
            "--meta",
            f"trainNumber={train_number};label={args.ccdb_labels[idx]}",
        ]
        result = subprocess.run(upload_cmd, capture_output=True, check=True)

        if "html" in result.stdout.decode("utf-8"):
            print(
                f"\n[!] Something went wrong with upload request: {result.stdout.decode('utf-8')}",
                file=sys.stderr,
            )
            sys.exit(1)

        if result.stderr:
            print(f"\n[!] Error while uploading: {result.stderr.decode('utf-8')}", file=sys.stderr)
            sys.exit(1)

        print("complete!", file=sys.stderr)

except subprocess.CalledProcessError as error:
    print(f"\n[!] Error while uploading: {error.stderr.decode('utf-8')}", file=sys.stderr)
    sys.exit(1)

print()
print("[✓] Success!", file=sys.stderr)
