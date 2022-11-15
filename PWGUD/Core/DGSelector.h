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

#ifndef O2_ANALYSIS_DGSELECTOR_H_
#define O2_ANALYSIS_DGSELECTOR_H_

#include "Framework/Logger.h"
#include "TDatabasePDG.h"
#include "PWGUD/Core/DGCutparHolder.h"

// -----------------------------------------------------------------------------
// add here Selectors for different types of diffractive events
// Selector for Double Gap events
class DGSelector {
 public:
  // constructor/destructor
  DGSelector();
  ~DGSelector();
  
  void Print()
  {
    LOGF(info, "In DGSelector"); 
  }
  
  // Function to check if collisions passes filter
  template <typename CC, typename BCs, typename TCs, typename FWs>
  int IsSelected(DGCutparHolder diffCuts, CC const& collision, BCs& bcRange, TCs& tracks, FWs& fwdtracks);

  //template <typename BCs, typename TCs, typename FWs>
  //int IsSelected(DGCutparHolder diffCuts, BCs& bcRange, TCs& tracks, FWs& fwdtracks);

 private:
  TDatabasePDG* fPDG;
};

#endif // O2_ANALYSIS_DGSELECTOR_H_
