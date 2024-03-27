// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "FastJetUtilities.h"

void fastjetutilities::setFastJetUserInfo(std::vector<fastjet::PseudoJet>& constituents, int index, int status)
{
  fastjet_user_info* user_info = new fastjet_user_info(status, index); // FIXME: can setting this as a pointer be avoided?
  constituents.back().set_user_info(user_info);
  if (index != -99999999) { // FIXME: needed for constituent subtraction as user_info is not propagated, but need to be quite careful to make sure indices dont overlap between tracks, clusters and HF candidates. Current solution might not be optimal
    int i = index;
    if (status == static_cast<int>(JetConstituentStatus::track)) {
      i = i + 1;
    }
    if (status == static_cast<int>(JetConstituentStatus::cluster)) {
      i = -1 * (i + 1);
    }
    if (status == static_cast<int>(JetConstituentStatus::candidateHF)) {
      i = 0;
    }
    // TODO : add V0s
    constituents.back().set_user_index(i); // FIXME: needed for constituent subtraction, but need to be quite careful to make sure indices dont overlap between tracks, clusters and HF candidates. Current solution might not be optimal
  }
}
