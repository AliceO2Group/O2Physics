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

#ifndef PWGUD_CORE_UPCHELPERS_H_
#define PWGUD_CORE_UPCHELPERS_H_

#include "UPCCutparHolder.h"

#include "PWGUD/DataModel/UDTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"

#include "CommonConstants/LHCConstants.h"
#include "Framework/AnalysisDataModel.h"

#include "TLorentzVector.h"

using BCsWithBcSels = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;

using ForwardTracks = o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov>;

using BarrelTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA,
                                   o2::aod::pidTPCFullEl, o2::aod::pidTPCFullMu, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr,
                                   o2::aod::TOFSignal, o2::aod::pidTOFbeta,
                                   o2::aod::pidTOFFullEl, o2::aod::pidTOFFullMu, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr>;

// namespace with helpers for UPC track skimming and candidate production
namespace upchelpers
{

// forward tracks selection counters
enum FwdSels {
  kFwdSelAll = 0,
  kFwdSelPt,
  kFwdSelEta,
  kFwdSelRabs,
  kFwdSelpDCA,
  kFwdSelChi2,
  kNFwdSels
};

// central barrel tracks selection counters
enum BarrelSels {
  kBarrelSelAll = 0,
  kBarrelSelHasTOF,
  kBarrelSelPt,
  kBarrelSelEta,
  kBarrelSelITSNCls,
  kBarrelSelITSChi2,
  kBarrelSelTPCNCls,
  kBarrelSelTPCChi2,
  kBarrelSelDCAXY,
  kBarrelSelDCAZ,
  kAmbiguous,
  kNBarrelSels
};

// FIT info holder
struct FITInfo {
  float ampFT0A = -1;
  float ampFT0C = -1;
  float timeFT0A = -999.;
  float timeFT0C = -999.;
  uint8_t triggerMaskFT0 = 0;
  float ampFDDA = -1;
  float ampFDDC = -1;
  float timeFDDA = -999.;
  float timeFDDC = -999.;
  uint8_t triggerMaskFDD = 0;
  float ampFV0A = -1;
  float timeFV0A = -999.;
  uint8_t triggerMaskFV0A = 0;
  int32_t BBFT0Apf = 0;
  int32_t BGFT0Apf = 0;
  int32_t BBFT0Cpf = 0;
  int32_t BGFT0Cpf = 0;
  int32_t BBFV0Apf = 0;
  int32_t BGFV0Apf = 0;
  int32_t BBFDDApf = 0;
  int32_t BGFDDApf = 0;
  int32_t BBFDDCpf = 0;
  int32_t BGFDDCpf = 0;
  int32_t distClosestBcTOR = 999;
  int32_t distClosestBcTSC = 999;
  int32_t distClosestBcTVX = 999;
  int32_t distClosestBcV0A = 999;
  int32_t distClosestBcT0A = 999;
};

template <typename T, typename TSelectorsArray>
void applyFwdCuts(UPCCutparHolder& upcCuts, const T& track, TSelectorsArray& fwdSelectors)
{
  fwdSelectors[kFwdSelPt] = track.pt() > upcCuts.getFwdPtLow() && track.pt() < upcCuts.getFwdPtHigh();                                                        // check pt
  fwdSelectors[kFwdSelEta] = track.eta() > upcCuts.getFwdEtaLow() && track.eta() < upcCuts.getFwdEtaHigh();                                                   // check pseudorapidity
  fwdSelectors[kFwdSelRabs] = track.rAtAbsorberEnd() > upcCuts.getMuonRAtAbsorberEndLow() && track.rAtAbsorberEnd() < upcCuts.getMuonRAtAbsorberEndHigh();    // check muon R
  fwdSelectors[kFwdSelpDCA] = track.rAtAbsorberEnd() < 26.5 ? track.pDca() < upcCuts.getMuonPDcaHighFirst() : track.pDca() < upcCuts.getMuonPDcaHighSecond(); // check pDCA
  fwdSelectors[kFwdSelChi2] = track.chi2() > upcCuts.getFwdChi2Low() && track.chi2() < upcCuts.getFwdChi2High();                                              // check chi2
}

template <typename T, typename TSelectorsArray>
void applyBarrelCuts(UPCCutparHolder& upcCuts, const T& track, TSelectorsArray& barrelSelectors)
{
  barrelSelectors[kAmbiguous] = true;
  if (upcCuts.getAmbigSwitch())
    barrelSelectors[kAmbiguous] = track.isPVContributor();
  barrelSelectors[kBarrelSelHasTOF] = true;
  if (upcCuts.getRequireTOF())
    barrelSelectors[kBarrelSelHasTOF] = track.hasTOF();                                                           // require TOF match if needed
  barrelSelectors[kBarrelSelPt] = track.pt() > upcCuts.getBarPtLow() && track.pt() < upcCuts.getBarPtHigh();      // check Pt cuts
  barrelSelectors[kBarrelSelEta] = track.eta() > upcCuts.getBarEtaLow() && track.eta() < upcCuts.getBarEtaHigh(); // check pseudorapidity cuts

  // check ITS cuts
  barrelSelectors[kBarrelSelITSNCls] = track.itsNCls() >= static_cast<uint8_t>(upcCuts.getITSNClusLow()) && track.itsNCls() <= static_cast<uint8_t>(upcCuts.getITSNClusHigh());
  barrelSelectors[kBarrelSelITSChi2] = track.itsChi2NCl() > upcCuts.getITSChi2Low() && track.itsChi2NCl() < upcCuts.getITSChi2High();

  // check TPC cuts
  barrelSelectors[kBarrelSelTPCNCls] = track.tpcNClsFound() > static_cast<int16_t>(upcCuts.getTPCNClsLow()) && track.tpcNClsFound() < static_cast<int16_t>(upcCuts.getTPCNClsHigh());
  barrelSelectors[kBarrelSelTPCChi2] = track.tpcChi2NCl() > upcCuts.getTPCChi2Low() && track.tpcChi2NCl() < upcCuts.getTPCChi2High();

  // check DCA
  barrelSelectors[kBarrelSelDCAZ] = track.dcaZ() > upcCuts.getDcaZLow() && track.dcaZ() < upcCuts.getDcaZHigh();

  if (upcCuts.getCheckMaxDcaXY()) {
    float dca = track.dcaXY();
    float maxDCA = 0.0105f + 0.0350f / pow(track.pt(), 1.1f);
    barrelSelectors[kBarrelSelDCAXY] = dca < maxDCA;
  } else {
    barrelSelectors[kBarrelSelDCAXY] = true;
  }
}

} // namespace upchelpers

#endif // PWGUD_CORE_UPCHELPERS_H_
