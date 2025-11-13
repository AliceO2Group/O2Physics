#!/bin/bash
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

OPTION="$1"
EXECUTABLE="o2-analysis-pid-feature-extractor"

o2-analysis-timestamp "${OPTION}" | \
  o2-analysis-event-selection "${OPTION}" | \
  o2-analysis-tracks-extra-v002-converter "${OPTION}" | \
  o2-analysis-track-propagation "${OPTION}" | \
  o2-analysis-pid-tpc-base "${OPTION}" | \
  o2-analysis-pid-tpc "${OPTION}" | \
  o2-analysis-pid-tof-base "${OPTION}" | \
  o2-analysis-pid-tof "${OPTION}" | \
  o2-analysis-pid-tof-beta "${OPTION}" | \
  o2-analysis-multiplicity-table "${OPTION}" | \
  o2-analysis-mccollision-converter "${OPTION}" | \
  "${EXECUTABLE}" "${OPTION}"