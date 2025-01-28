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
"""

import argparse
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import ROOT  # pylint: disable=import-error

parser = argparse.ArgumentParser(description="A tool to calculate efficiency and upload it to CCDB")
parser.add_argument(
    "run_dir",
    type=Path,
    help="Path to run directory with analysis results",
)
parser.add_argument(
    "--task",
    "-t",
    type=str,
    help="Task name to query histograms from",
    required=True,
)
parser.add_argument(
    "--ccdb-path",
    "-P",
    type=str,
    help="Path where to save in CCDB",
    required=True,
)
parser.add_argument(
    "--ccdb-url",
    "-U",
    type=str,
    help="Path where to save in CCDB",
    default="http://ccdb-test.cern.ch:8080",
)
parser.add_argument(
    "--ccdb-labels",
    "-L",
    type=str,
    nargs="+",
    help="Labels for objects' metadata in CCDB",
    default=[""] * 2,
)
parser.add_argument(
    "--ccdb-lifetime",
    "-T",
    type=int,
    help="How long should the objects' validity in CCDB last (default: 1 year)",
    default=365 * 24 * 60 * 60 * 1000,  # one year
)
args = parser.parse_args()

ANALYSIS_RESULTS = "AnalysisResults.root"
results_path = args.run_dir / ANALYSIS_RESULTS
job_id = results_path.parent.name
train_number = results_path.parent.parent.name

tmp_dir = Path(tempfile.gettempdir())

res_dest = tmp_dir / f"{train_number}-{job_id}-{ANALYSIS_RESULTS}"
eff_dest = tmp_dir / f"{train_number}-{job_id}-Efficiency.root"

# get file from alien
if not res_dest.is_file():
    print(f"[↓] Downloading analysis results from Alien to {res_dest} ...", file=sys.stderr)
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
        sys.exit(0)
else:
    print(
        f"[-] Skipping download from Alien, since {res_dest} is already present",
        file=sys.stderr,
    )

particles = {1: "one", 2: "two"}
histos_to_upload = []

# get reco & truth histos
with (
    ROOT.TFile.Open(res_dest.as_uri()) as res_file,
    ROOT.TFile.Open(eff_dest.as_uri(), "recreate") as eff_file,
):
    for idx, part_num in particles.items():
        reco = res_file.Get(f"{args.task}/Tracks_{part_num}_MC/hPt")
        truth = res_file.Get(f"{args.task}/MCTruthTracks_{part_num}/hPt")
        if not reco and not truth:
            print(
                f"[-] No MC Reco nor MC Truth histogram found for particle {part_num}",
                file=sys.stderr,
            )
            continue

        num_bins = reco.GetNbinsX()
        x_max = reco.GetXaxis().GetBinLowEdge(num_bins) + reco.GetXaxis().GetBinWidth(num_bins)

        hist_name = f"Efficiency_part{idx}"

        # calculate efficiency
        eff = ROOT.TH1F(hist_name, "", num_bins, 0, x_max)
        for bin_idx in range(len(eff)):
            denom = truth.GetBinContent(bin_idx)
            eff.SetBinContent(bin_idx, reco.GetBinContent(bin_idx) / denom if denom > 0 else 0)

        # save efficiency object to file
        eff_file.WriteObject(eff, hist_name)
        histos_to_upload.append(hist_name)

if len(histos_to_upload) == 0:
    print("[-] Exiting, since there is nothing to upload", file=sys.stderr)
    sys.exit(0)

# upload objects to ccdb
try:
    for idx, key in enumerate(histos_to_upload):
        timestamp_start = int(time.time() * 1000)
        timestamp_end = timestamp_start + args.ccdb_lifetime

        print(f"[↑] Uploading {key} to CCDB ... ", file=sys.stderr, end="")
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
        if result.stderr:
            print(f"[!] Error while uploading: {result.stderr.decode('utf-8')}", file=sys.stderr)
            sys.exit(0)

        print("complete!", file=sys.stderr)

except subprocess.CalledProcessError as error:
    print(f"[!] Error while uploading: {error.stderr.decode('utf-8')}", file=sys.stderr)
    sys.exit(0)

print("[✓] Success!", file=sys.stderr)
