#!/usr/bin/env python3

import subprocess
import tempfile
import time
from pathlib import Path

import ROOT

# import argparse
# parser = argparse.ArgumentParser(description="a tool to calculate efficiency and upload it to CCDB")
# args = parser.parse_args()

analysis_results = "AnalysisResults.root"
alien_path = Path("/alice/cern.ch/user/a/alihyperloop/outputs/0033/332611/70301") / analysis_results
job_dir = alien_path.parent.name

tmp_dir = Path(tempfile.gettempdir())

dest_in = tmp_dir / f"{job_dir}-{analysis_results}"
dest_out = tmp_dir / f"{job_dir}-Efficiency.root"

# get file from alien
if not dest_in.is_file():
    ROOT.TGrid.Connect("alien://")
    try:
        subprocess.run(
            ["alien_cp", alien_path, "file://" + str(dest_in)],
            capture_output=True,
            check=True,
        )
    except subprocess.CalledProcessError as error:
        print(f"[!] Error while downloading results file: {error.stderr}")

# get reco & truth histos
with ROOT.TFile.Open(dest_in.as_uri()) as infile:
    truth = infile["femto-universe-pair-task-track-track-extended/MCTruthTracks_one/hPt"]
    reco = infile["femto-universe-pair-task-track-track-extended/Tracks_one_MC/hPt"]

    # calculate efficiency
    eff = ROOT.TH1F("Efficiency/part1", "Efficiency origin/generated part2; p_{T} (GeV/c); Efficiency", 100, 0, 4)
    for bin_idx in range(len(eff)):
        denom = truth.GetBinContent(bin_idx)
        eff.SetBinContent(bin_idx, reco.GetBinContent(bin_idx) / denom if denom > 0 else 0)

    # save efficiency histo
    with ROOT.TFile.Open(dest_out.as_uri(), "recreate") as outfile:
        outfile.WriteObject(eff, "Efficiency_part1")

# upload object to ccdb
ccdb_path = "Users/d/dkarpins/Efficiency"
ccdb_url = "http://ccdb-test.cern.ch:8080"

timestamp_start = int(time.time() * 1000)
timestamp_end = 365 * 24 * 60 * 60 * 60 * 1000

try:
    print("[^] Uploading to CCDB...")
    upload_cmd = [
        "o2-ccdb-upload",
        "--file",
        dest_out.as_uri(),
        "--host",
        ccdb_url,
        "--key",
        "Efficiency_part1",
        "--path",
        ccdb_path,
        "--starttimestamp",
        str(timestamp_start),
        "--endtimestamp",
        str(timestamp_end),
        "--meta",
        "label=proton1",
    ]
    result = subprocess.run(upload_cmd, capture_output=True, check=True)
    print(result.stdout.decode("utf-8"))
    if result.stderr:
        print(f"[!] Error while uploading: {result.stderr.decode('utf-8')}")

except subprocess.CalledProcessError as error:
    print(f"[!] Error while uploading: {error.stderr.decode('utf-8')}")
