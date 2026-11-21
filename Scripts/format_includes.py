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

"""!
@brief  Fix include style
@author Vít Kučera <vit.kucera@cern.ch>, Inha University
@date   2026-10-08

NB: Run before sorting.

Usage: format_includes.py FILE [FILE ...]
"""

import re
import sys

LOCAL = ('"', '"')
EXTERNAL = ("<", ">")

INCLUDE = re.compile(r"(\s*#include\s+)(\S+)")


def fix_line(line: str) -> str:
    m = INCLUDE.match(line)
    if not m:
        return line
    pre, tok = m.groups()
    rest = line[m.end() :]
    h = tok[1:-1]

    if re.match(r"(PWG[A-Z]{2}|Common|ALICE3|DPG|EventFiltering|PID|Tools|Tutorials)/.*\.h", h):
        d = LOCAL  # O2Physics
    elif re.match(
        r"(Algorithm|CCDB|Common[A-Z]|DataFormats|DCAFitter|Detectors|EMCAL|FDD|Field|Framework|FT0|FV0|GlobalTracking|GPU|ITS|MathUtils|MCH|MFT|MID|PHOS|ReconstructionDataFormats|SimulationDataFormat|TOF|TPC|ZDC).*/.*\.h",
        h,
    ):
        d = EXTERNAL  # O2
    elif re.match(r"(T[A-Z]|Math/|Roo[A-Z])[A-Za-z0-9/]+\.h", h):
        d = EXTERNAL  # ROOT
    elif re.match(r"KF[A-Z][A-Za-z0-9]+\.h", h):
        d = EXTERNAL  # KFParticle
    elif re.match(r"(fastjet/|onnxruntime)", h):
        d = EXTERNAL  # FastJet, ONNX
    elif re.match(r".*DataModel/", h):
        d = LOCAL  # incomplete path to DataModel
    elif re.match(r"([A-Za-z0-9_]+/)+[A-Za-z0-9_]+\.h", h):
        d = EXTERNAL  # other third-party
    elif re.match(r'".*\.', tok):
        return line  # other local-looking file
    elif re.match(r"[a-z_]+\.h", h):
        d = EXTERNAL  # C system
    elif re.match(r"[a-z_/]+$", h):
        d = EXTERNAL  # C++ system (whole string)
    else:
        return line

    return pre + d[0] + h + d[1] + rest


def process(path: str):
    # newline='' keeps the original line endings untouched.
    with open(path, newline="") as f:
        lines = f.readlines()
    new_lines = [fix_line(line) for line in lines]
    if new_lines != lines:
        with open(path, "w", newline="") as f:
            f.writelines(new_lines)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for p in sys.argv[1:]:
        process(p)
