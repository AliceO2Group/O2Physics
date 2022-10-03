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

///
/// \file   qaPIDTOF.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Header file for QA tasks of the TOF PID quantities
///

#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Task to produce the TOF QA plots
struct tofPidQa {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr std::string_view hexpected[Np] = {"expected/El", "expected/Mu", "expected/Pi",
                                                     "expected/Ka", "expected/Pr", "expected/De",
                                                     "expected/Tr", "expected/He", "expected/Al"};
  static constexpr std::string_view hdelta[Np] = {"delta/El", "delta/Mu", "delta/Pi",
                                                  "delta/Ka", "delta/Pr", "delta/De",
                                                  "delta/Tr", "delta/He", "delta/Al"};
  static constexpr std::string_view hdelta_pt_pos[Np] = {"delta/pt/pos/El", "delta/pt/pos/Mu", "delta/pt/pos/Pi",
                                                         "delta/pt/pos/Ka", "delta/pt/pos/Pr", "delta/pt/pos/De",
                                                         "delta/pt/pos/Tr", "delta/pt/pos/He", "delta/pt/pos/Al"};
  static constexpr std::string_view hdelta_pt_neg[Np] = {"delta/pt/neg/El", "delta/pt/neg/Mu", "delta/pt/neg/Pi",
                                                         "delta/pt/neg/Ka", "delta/pt/neg/Pr", "delta/pt/neg/De",
                                                         "delta/pt/neg/Tr", "delta/pt/neg/He", "delta/pt/neg/Al"};
  // Ev. Time fill
  static constexpr std::string_view hdelta_evtime_fill[Np] = {"delta/evtime/fill/El", "delta/evtime/fill/Mu", "delta/evtime/fill/Pi",
                                                              "delta/evtime/fill/Ka", "delta/evtime/fill/Pr", "delta/evtime/fill/De",
                                                              "delta/evtime/fill/Tr", "delta/evtime/fill/He", "delta/evtime/fill/Al"};
  static constexpr std::string_view hdelta_pt_pos_evtime_fill[Np] = {"delta/pt/pos/evtime/fill/El", "delta/pt/pos/evtime/fill/Mu", "delta/pt/pos/evtime/fill/Pi",
                                                                     "delta/pt/pos/evtime/fill/Ka", "delta/pt/pos/evtime/fill/Pr", "delta/pt/pos/evtime/fill/De",
                                                                     "delta/pt/pos/evtime/fill/Tr", "delta/pt/pos/evtime/fill/He", "delta/pt/pos/evtime/fill/Al"};
  static constexpr std::string_view hdelta_pt_neg_evtime_fill[Np] = {"delta/pt/neg/evtime/fill/El", "delta/pt/neg/evtime/fill/Mu", "delta/pt/neg/evtime/fill/Pi",
                                                                     "delta/pt/neg/evtime/fill/Ka", "delta/pt/neg/evtime/fill/Pr", "delta/pt/neg/evtime/fill/De",
                                                                     "delta/pt/neg/evtime/fill/Tr", "delta/pt/neg/evtime/fill/He", "delta/pt/neg/evtime/fill/Al"};
  // Ev. Time TOF
  static constexpr std::string_view hdelta_evtime_tof[Np] = {"delta/evtime/tof/El", "delta/evtime/tof/Mu", "delta/evtime/tof/Pi",
                                                             "delta/evtime/tof/Ka", "delta/evtime/tof/Pr", "delta/evtime/tof/De",
                                                             "delta/evtime/tof/Tr", "delta/evtime/tof/He", "delta/evtime/tof/Al"};
  static constexpr std::string_view hdelta_pt_pos_evtime_tof[Np] = {"delta/pt/pos/evtime/tof/El", "delta/pt/pos/evtime/tof/Mu", "delta/pt/pos/evtime/tof/Pi",
                                                                    "delta/pt/pos/evtime/tof/Ka", "delta/pt/pos/evtime/tof/Pr", "delta/pt/pos/evtime/tof/De",
                                                                    "delta/pt/pos/evtime/tof/Tr", "delta/pt/pos/evtime/tof/He", "delta/pt/pos/evtime/tof/Al"};
  static constexpr std::string_view hdelta_pt_neg_evtime_tof[Np] = {"delta/pt/neg/evtime/tof/El", "delta/pt/neg/evtime/tof/Mu", "delta/pt/neg/evtime/tof/Pi",
                                                                    "delta/pt/neg/evtime/tof/Ka", "delta/pt/neg/evtime/tof/Pr", "delta/pt/neg/evtime/tof/De",
                                                                    "delta/pt/neg/evtime/tof/Tr", "delta/pt/neg/evtime/tof/He", "delta/pt/neg/evtime/tof/Al"};
  // Ev. Time FT0
  static constexpr std::string_view hdelta_evtime_ft0[Np] = {"delta/evtime/ft0/El", "delta/evtime/ft0/Mu", "delta/evtime/ft0/Pi",
                                                             "delta/evtime/ft0/Ka", "delta/evtime/ft0/Pr", "delta/evtime/ft0/De",
                                                             "delta/evtime/ft0/Tr", "delta/evtime/ft0/He", "delta/evtime/ft0/Al"};
  static constexpr std::string_view hdelta_pt_pos_evtime_ft0[Np] = {"delta/pt/pos/evtime/ft0/El", "delta/pt/pos/evtime/ft0/Mu", "delta/pt/pos/evtime/ft0/Pi",
                                                                    "delta/pt/pos/evtime/ft0/Ka", "delta/pt/pos/evtime/ft0/Pr", "delta/pt/pos/evtime/ft0/De",
                                                                    "delta/pt/pos/evtime/ft0/Tr", "delta/pt/pos/evtime/ft0/He", "delta/pt/pos/evtime/ft0/Al"};
  static constexpr std::string_view hdelta_pt_neg_evtime_ft0[Np] = {"delta/pt/neg/evtime/ft0/El", "delta/pt/neg/evtime/ft0/Mu", "delta/pt/neg/evtime/ft0/Pi",
                                                                    "delta/pt/neg/evtime/ft0/Ka", "delta/pt/neg/evtime/ft0/Pr", "delta/pt/neg/evtime/ft0/De",
                                                                    "delta/pt/neg/evtime/ft0/Tr", "delta/pt/neg/evtime/ft0/He", "delta/pt/neg/evtime/ft0/Al"};
  // Ev. Time TOF+FT0
  static constexpr std::string_view hdelta_evtime_tofft0[Np] = {"delta/evtime/tofft0/El", "delta/evtime/tofft0/Mu", "delta/evtime/tofft0/Pi",
                                                                "delta/evtime/tofft0/Ka", "delta/evtime/tofft0/Pr", "delta/evtime/tofft0/De",
                                                                "delta/evtime/tofft0/Tr", "delta/evtime/tofft0/He", "delta/evtime/tofft0/Al"};
  static constexpr std::string_view hdelta_pt_pos_evtime_tofft0[Np] = {"delta/pt/pos/evtime/tofft0/El", "delta/pt/pos/evtime/tofft0/Mu", "delta/pt/pos/evtime/tofft0/Pi",
                                                                       "delta/pt/pos/evtime/tofft0/Ka", "delta/pt/pos/evtime/tofft0/Pr", "delta/pt/pos/evtime/tofft0/De",
                                                                       "delta/pt/pos/evtime/tofft0/Tr", "delta/pt/pos/evtime/tofft0/He", "delta/pt/pos/evtime/tofft0/Al"};
  static constexpr std::string_view hdelta_pt_neg_evtime_tofft0[Np] = {"delta/pt/neg/evtime/tofft0/El", "delta/pt/neg/evtime/tofft0/Mu", "delta/pt/neg/evtime/tofft0/Pi",
                                                                       "delta/pt/neg/evtime/tofft0/Ka", "delta/pt/neg/evtime/tofft0/Pr", "delta/pt/neg/evtime/tofft0/De",
                                                                       "delta/pt/neg/evtime/tofft0/Tr", "delta/pt/neg/evtime/tofft0/He", "delta/pt/neg/evtime/tofft0/Al"};
  static constexpr std::string_view hexpsigma[Np] = {"expsigma/El", "expsigma/Mu", "expsigma/Pi",
                                                     "expsigma/Ka", "expsigma/Pr", "expsigma/De",
                                                     "expsigma/Tr", "expsigma/He", "expsigma/Al"};
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi",
                                                   "nsigma/Ka", "nsigma/Pr", "nsigma/De",
                                                   "nsigma/Tr", "nsigma/He", "nsigma/Al"};
  static constexpr std::string_view hnsigma_pt[Np] = {"nsigma/pt/El", "nsigma/pt/Mu", "nsigma/pt/Pi",
                                                      "nsigma/pt/Ka", "nsigma/pt/Pr", "nsigma/pt/De",
                                                      "nsigma/pt/Tr", "nsigma/pt/He", "nsigma/pt/Al"};
  static constexpr std::string_view hnsigma_pt_pos[Np] = {"nsigma/pt/pos/El", "nsigma/pt/pos/Mu", "nsigma/pt/pos/Pi",
                                                          "nsigma/pt/pos/Ka", "nsigma/pt/pos/Pr", "nsigma/pt/pos/De",
                                                          "nsigma/pt/pos/Tr", "nsigma/pt/pos/He", "nsigma/pt/pos/Al"};
  static constexpr std::string_view hnsigma_pt_neg[Np] = {"nsigma/pt/neg/El", "nsigma/pt/neg/Mu", "nsigma/pt/neg/Pi",
                                                          "nsigma/pt/neg/Ka", "nsigma/pt/neg/Pr", "nsigma/pt/neg/De",
                                                          "nsigma/pt/neg/Tr", "nsigma/pt/neg/He", "nsigma/pt/neg/Al"};
  // Ev. Time fill
  static constexpr std::string_view hnsigma_evtime_fill[Np] = {"nsigma/evtime/fill/El", "nsigma/evtime/fill/Mu", "nsigma/evtime/fill/Pi",
                                                               "nsigma/evtime/fill/Ka", "nsigma/evtime/fill/Pr", "nsigma/evtime/fill/De",
                                                               "nsigma/evtime/fill/Tr", "nsigma/evtime/fill/He", "nsigma/evtime/fill/Al"};
  static constexpr std::string_view hnsigma_pt_evtime_fill[Np] = {"nsigma/pt/evtime/fill/El", "nsigma/pt/evtime/fill/Mu", "nsigma/pt/evtime/fill/Pi",
                                                                  "nsigma/pt/evtime/fill/Ka", "nsigma/pt/evtime/fill/Pr", "nsigma/pt/evtime/fill/De",
                                                                  "nsigma/pt/evtime/fill/Tr", "nsigma/pt/evtime/fill/He", "nsigma/pt/evtime/fill/Al"};
  static constexpr std::string_view hnsigma_pt_pos_evtime_fill[Np] = {"nsigma/pt/pos/evtime/fill/El", "nsigma/pt/pos/evtime/fill/Mu", "nsigma/pt/pos/evtime/fill/Pi",
                                                                      "nsigma/pt/pos/evtime/fill/Ka", "nsigma/pt/pos/evtime/fill/Pr", "nsigma/pt/pos/evtime/fill/De",
                                                                      "nsigma/pt/pos/evtime/fill/Tr", "nsigma/pt/pos/evtime/fill/He", "nsigma/pt/pos/evtime/fill/Al"};
  static constexpr std::string_view hnsigma_pt_neg_evtime_fill[Np] = {"nsigma/pt/neg/evtime/fill/El", "nsigma/pt/neg/evtime/fill/Mu", "nsigma/pt/neg/evtime/fill/Pi",
                                                                      "nsigma/pt/neg/evtime/fill/Ka", "nsigma/pt/neg/evtime/fill/Pr", "nsigma/pt/neg/evtime/fill/De",
                                                                      "nsigma/pt/neg/evtime/fill/Tr", "nsigma/pt/neg/evtime/fill/He", "nsigma/pt/neg/evtime/fill/Al"};
  // Ev. Time TOF
  static constexpr std::string_view hnsigma_evtime_tof[Np] = {"nsigma/evtime/tof/El", "nsigma/evtime/tof/Mu", "nsigma/evtime/tof/Pi",
                                                              "nsigma/evtime/tof/Ka", "nsigma/evtime/tof/Pr", "nsigma/evtime/tof/De",
                                                              "nsigma/evtime/tof/Tr", "nsigma/evtime/tof/He", "nsigma/evtime/tof/Al"};
  static constexpr std::string_view hnsigma_pt_evtime_tof[Np] = {"nsigma/pt/evtime/tof/El", "nsigma/pt/evtime/tof/Mu", "nsigma/pt/evtime/tof/Pi",
                                                                 "nsigma/pt/evtime/tof/Ka", "nsigma/pt/evtime/tof/Pr", "nsigma/pt/evtime/tof/De",
                                                                 "nsigma/pt/evtime/tof/Tr", "nsigma/pt/evtime/tof/He", "nsigma/pt/evtime/tof/Al"};
  static constexpr std::string_view hnsigma_pt_pos_evtime_tof[Np] = {"nsigma/pt/pos/evtime/tof/El", "nsigma/pt/pos/evtime/tof/Mu", "nsigma/pt/pos/evtime/tof/Pi",
                                                                     "nsigma/pt/pos/evtime/tof/Ka", "nsigma/pt/pos/evtime/tof/Pr", "nsigma/pt/pos/evtime/tof/De",
                                                                     "nsigma/pt/pos/evtime/tof/Tr", "nsigma/pt/pos/evtime/tof/He", "nsigma/pt/pos/evtime/tof/Al"};
  static constexpr std::string_view hnsigma_pt_neg_evtime_tof[Np] = {"nsigma/pt/neg/evtime/tof/El", "nsigma/pt/neg/evtime/tof/Mu", "nsigma/pt/neg/evtime/tof/Pi",
                                                                     "nsigma/pt/neg/evtime/tof/Ka", "nsigma/pt/neg/evtime/tof/Pr", "nsigma/pt/neg/evtime/tof/De",
                                                                     "nsigma/pt/neg/evtime/tof/Tr", "nsigma/pt/neg/evtime/tof/He", "nsigma/pt/neg/evtime/tof/Al"};
  // Ev. Time FT0
  static constexpr std::string_view hnsigma_evtime_ft0[Np] = {"nsigma/evtime/ft0/El", "nsigma/evtime/ft0/Mu", "nsigma/evtime/ft0/Pi",
                                                              "nsigma/evtime/ft0/Ka", "nsigma/evtime/ft0/Pr", "nsigma/evtime/ft0/De",
                                                              "nsigma/evtime/ft0/Tr", "nsigma/evtime/ft0/He", "nsigma/evtime/ft0/Al"};
  static constexpr std::string_view hnsigma_pt_evtime_ft0[Np] = {"nsigma/pt/evtime/ft0/El", "nsigma/pt/evtime/ft0/Mu", "nsigma/pt/evtime/ft0/Pi",
                                                                 "nsigma/pt/evtime/ft0/Ka", "nsigma/pt/evtime/ft0/Pr", "nsigma/pt/evtime/ft0/De",
                                                                 "nsigma/pt/evtime/ft0/Tr", "nsigma/pt/evtime/ft0/He", "nsigma/pt/evtime/ft0/Al"};
  static constexpr std::string_view hnsigma_pt_pos_evtime_ft0[Np] = {"nsigma/pt/pos/evtime/ft0/El", "nsigma/pt/pos/evtime/ft0/Mu", "nsigma/pt/pos/evtime/ft0/Pi",
                                                                     "nsigma/pt/pos/evtime/ft0/Ka", "nsigma/pt/pos/evtime/ft0/Pr", "nsigma/pt/pos/evtime/ft0/De",
                                                                     "nsigma/pt/pos/evtime/ft0/Tr", "nsigma/pt/pos/evtime/ft0/He", "nsigma/pt/pos/evtime/ft0/Al"};
  static constexpr std::string_view hnsigma_pt_neg_evtime_ft0[Np] = {"nsigma/pt/neg/evtime/ft0/El", "nsigma/pt/neg/evtime/ft0/Mu", "nsigma/pt/neg/evtime/ft0/Pi",
                                                                     "nsigma/pt/neg/evtime/ft0/Ka", "nsigma/pt/neg/evtime/ft0/Pr", "nsigma/pt/neg/evtime/ft0/De",
                                                                     "nsigma/pt/neg/evtime/ft0/Tr", "nsigma/pt/neg/evtime/ft0/He", "nsigma/pt/neg/evtime/ft0/Al"};
  // Ev. Time TOF+FT0
  static constexpr std::string_view hnsigma_evtime_tofft0[Np] = {"nsigma/evtime/tofft0/El", "nsigma/evtime/tofft0/Mu", "nsigma/evtime/tofft0/Pi",
                                                                 "nsigma/evtime/tofft0/Ka", "nsigma/evtime/tofft0/Pr", "nsigma/evtime/tofft0/De",
                                                                 "nsigma/evtime/tofft0/Tr", "nsigma/evtime/tofft0/He", "nsigma/evtime/tofft0/Al"};
  static constexpr std::string_view hnsigma_pt_evtime_tofft0[Np] = {"nsigma/pt/evtime/tofft0/El", "nsigma/pt/evtime/tofft0/Mu", "nsigma/pt/evtime/tofft0/Pi",
                                                                    "nsigma/pt/evtime/tofft0/Ka", "nsigma/pt/evtime/tofft0/Pr", "nsigma/pt/evtime/tofft0/De",
                                                                    "nsigma/pt/evtime/tofft0/Tr", "nsigma/pt/evtime/tofft0/He", "nsigma/pt/evtime/tofft0/Al"};
  static constexpr std::string_view hnsigma_pt_pos_evtime_tofft0[Np] = {"nsigma/pt/pos/evtime/tofft0/El", "nsigma/pt/pos/evtime/tofft0/Mu", "nsigma/pt/pos/evtime/tofft0/Pi",
                                                                        "nsigma/pt/pos/evtime/tofft0/Ka", "nsigma/pt/pos/evtime/tofft0/Pr", "nsigma/pt/pos/evtime/tofft0/De",
                                                                        "nsigma/pt/pos/evtime/tofft0/Tr", "nsigma/pt/pos/evtime/tofft0/He", "nsigma/pt/pos/evtime/tofft0/Al"};
  static constexpr std::string_view hnsigma_pt_neg_evtime_tofft0[Np] = {"nsigma/pt/neg/evtime/tofft0/El", "nsigma/pt/neg/evtime/tofft0/Mu", "nsigma/pt/neg/evtime/tofft0/Pi",
                                                                        "nsigma/pt/neg/evtime/tofft0/Ka", "nsigma/pt/neg/evtime/tofft0/Pr", "nsigma/pt/neg/evtime/tofft0/De",
                                                                        "nsigma/pt/neg/evtime/tofft0/Tr", "nsigma/pt/neg/evtime/tofft0/He", "nsigma/pt/neg/evtime/tofft0/Al"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> logAxis{"logAxis", 0, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 400, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.1f, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 5.f, "Maximum momentum in range"};
  Configurable<int> nBinsDelta{"nBinsDelta", 200, "Number of bins for the Delta"};
  Configurable<float> minDelta{"minDelta", -1000.f, "Minimum Delta in range"};
  Configurable<float> maxDelta{"maxDelta", 1000.f, "Maximum Delta in range"};
  Configurable<int> nBinsExpSigma{"nBinsExpSigma", 200, "Number of bins for the ExpSigma"};
  Configurable<float> minExpSigma{"minExpSigma", 0.f, "Minimum ExpSigma in range"};
  Configurable<float> maxExpSigma{"maxExpSigma", 200.f, "Maximum ExpSigma in range"};
  Configurable<int> nBinsNSigma{"nBinsNSigma", 200, "Number of bins for the NSigma"};
  Configurable<float> minNSigma{"minNSigma", -10.f, "Minimum NSigma in range"};
  Configurable<float> maxNSigma{"maxNSigma", 10.f, "Maximum NSigma in range"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply event selection cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<bool> applyTrackCut{"applyTrackCut", false, "Flag to apply standard track cuts"};
  Configurable<bool> applyRapidityCut{"applyRapidityCut", false, "Flag to apply rapidity cut"};
  Configurable<bool> enableEvTimeSplitting{"enableEvTimeSplitting", false, "Flag to enable histograms splitting depending on the Event Time used"};

  template <o2::track::PID::ID id>
  void initPerParticle(const AxisSpec& pAxis, const AxisSpec& ptAxis)
  {
    static_assert(id >= 0 && id <= PID::Alpha && "Particle index outside limits");
    bool enableFullHistos = false;
    int enabledProcesses = 0;
    switch (id) { // Skipping disabled particles
#define particleCase(particleId)                                 \
  case PID::particleId:                                          \
    if (!doprocess##particleId && !doprocessFull##particleId) {  \
      return;                                                    \
    }                                                            \
    if (doprocess##particleId) {                                 \
      enabledProcesses++;                                        \
    }                                                            \
    if (doprocessFull##particleId) {                             \
      enableFullHistos = true;                                   \
      enabledProcesses++;                                        \
    }                                                            \
    LOGF(info, "Enabled TOF QA for %s %s", #particleId, pT[id]); \
    break;

      particleCase(Electron);
      particleCase(Muon);
      particleCase(Pion);
      particleCase(Kaon);
      particleCase(Proton);
      particleCase(Deuteron);
      particleCase(Triton);
      particleCase(Helium3);
      particleCase(Alpha);
#undef particleCase
    }
    if (enabledProcesses != 1) {
      LOG(fatal) << "Cannot enable more than one process function per particle, check and retry!";
    }

    // NSigma
    const char* axisTitle = Form("N_{#sigma}^{TOF}(%s)", pT[id]);
    const AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, axisTitle};
    histos.add(hnsigma[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
    histos.add(hnsigma_pt[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigma_pt_pos[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigma_pt_neg[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    if (enableEvTimeSplitting) {
      histos.add(hnsigma_evtime_fill[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
      histos.add(hnsigma_pt_evtime_fill[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
      histos.add(hnsigma_pt_pos_evtime_fill[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
      histos.add(hnsigma_pt_neg_evtime_fill[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});

      histos.add(hnsigma_evtime_tof[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
      histos.add(hnsigma_pt_evtime_tof[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
      histos.add(hnsigma_pt_pos_evtime_tof[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
      histos.add(hnsigma_pt_neg_evtime_tof[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});

      histos.add(hnsigma_evtime_ft0[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
      histos.add(hnsigma_pt_evtime_ft0[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
      histos.add(hnsigma_pt_pos_evtime_ft0[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
      histos.add(hnsigma_pt_neg_evtime_ft0[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});

      histos.add(hnsigma_evtime_tofft0[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
      histos.add(hnsigma_pt_evtime_tofft0[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
      histos.add(hnsigma_pt_pos_evtime_tofft0[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
      histos.add(hnsigma_pt_neg_evtime_tofft0[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    }

    if (!enableFullHistos) { // Enabling only NSigma for tiny tables
      return;
    }

    // Exp signal
    const AxisSpec expAxis{1000, 0, 2e6, Form("t_{exp}(%s) (ps)", pT[id])};
    histos.add(hexpected[id].data(), "", kTH2F, {pAxis, expAxis});

    // Signal - Expected signal
    const AxisSpec deltaAxis{nBinsDelta, minDelta, maxDelta, Form("t-t_{ev}-t_{exp}(%s) (ps)", pT[id])};
    axisTitle = Form("#Delta^{TOF}(%s)", pT[id]);
    histos.add(hdelta[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});
    histos.add(hdelta_pt_pos[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});
    histos.add(hdelta_pt_neg[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});

    if (enableEvTimeSplitting) {
      histos.add(hdelta_evtime_fill[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});
      histos.add(hdelta_pt_pos_evtime_fill[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});
      histos.add(hdelta_pt_neg_evtime_fill[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});

      histos.add(hdelta_evtime_tof[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});
      histos.add(hdelta_pt_pos_evtime_tof[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});
      histos.add(hdelta_pt_neg_evtime_tof[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});

      histos.add(hdelta_evtime_ft0[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});
      histos.add(hdelta_pt_pos_evtime_ft0[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});
      histos.add(hdelta_pt_neg_evtime_ft0[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});

      histos.add(hdelta_evtime_tofft0[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});
      histos.add(hdelta_pt_pos_evtime_tofft0[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});
      histos.add(hdelta_pt_neg_evtime_tofft0[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});
    }

    // Exp Sigma
    const AxisSpec expSigmaAxis{nBinsExpSigma, minExpSigma, maxExpSigma, Form("Exp_{#sigma}^{TOF}(%s) (ps)", pT[id])};
    histos.add(hexpsigma[id].data(), "", kTH2F, {pAxis, expSigmaAxis});
  }

  void init(o2::framework::InitContext&)
  {
    const AxisSpec multAxis{100, 0, 100, "TOF multiplicity"};
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec etaAxis{100, -1, 1, "#it{#eta}"};
    const AxisSpec phiAxis{100, 0, TMath::TwoPi(), "#it{#phi}"};
    const AxisSpec colTimeAxis{100, -2000, 2000, "Collision time (ps)"};
    const AxisSpec colTimeResoAxis{100, 0, 1000, "#sigma_{Collision time} (ps)"};
    const AxisSpec lAxis{100, 0, 500, "Track length (cm)"};
    const AxisSpec ptResoAxis{100, 0, 0.1, "#sigma_{#it{p}_{T}}"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    AxisSpec pExpAxis{nBinsP, minP, maxP, "#it{p}_{Exp. TOF} (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogarithmic();
      pAxis.makeLogarithmic();
      pExpAxis.makeLogarithmic();
    }
    const AxisSpec tofAxis{10000, 0, 2e6, "TOF Signal (ps)"};

    // Event properties
    auto h = histos.add<TH1>("event/evsel", "", kTH1F, {{10, 0.5, 10.5, "Ev. Sel."}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Passed ev. sel.");
    h->GetXaxis()->SetBinLabel(3, "Passed mult.");
    h->GetXaxis()->SetBinLabel(4, "Passed vtx Z");

    h = histos.add<TH1>("event/trackselection", "", kTH1F, {{10, 0.5, 10.5, "Selection passed"}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "isGlobalTrack");
    h->GetXaxis()->SetBinLabel(3, "hasITS");
    h->GetXaxis()->SetBinLabel(4, "hasTPC");
    h->GetXaxis()->SetBinLabel(5, "hasTOF");

    histos.add("event/vertexz", "", kTH1F, {vtxZAxis});
    h = histos.add<TH1>("event/particlehypo", "", kTH1F, {{10, 0, 10, "PID in tracking"}});
    for (int i = 0; i < 9; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, PID::getName(i));
    }
    histos.add("event/trackmultiplicity", "", kTH1F, {multAxis});
    histos.add("event/tofmultiplicity", "", kTH1F, {multAxis});

    histos.add("event/evtime/colltime", "collisionTime()", kTH1F, {colTimeAxis});
    histos.add("event/evtime/colltimereso", "collisionTimeRes()", kTH2F, {multAxis, colTimeResoAxis});
    histos.add("event/evtime/undef", "Undefined event time", kTH1F, {colTimeAxis});
    histos.add("event/evtime/undefreso", "Undefined event time reso.", kTH2F, {multAxis, colTimeResoAxis});
    histos.add("event/evtime/avail", "Available event time", kTH1F, {colTimeAxis});
    histos.add("event/evtime/availreso", "Available event time reso.", kTH2F, {multAxis, colTimeResoAxis});
    histos.add("event/evtime/ft0tof", "FT0+TOF event time", kTH1F, {colTimeAxis});
    histos.add("event/evtime/ft0tofreso", "FT0+TOF event time reso.", kTH2F, {multAxis, colTimeResoAxis});
    histos.add("event/evtime/tof", "TOF event time", kTH1F, {colTimeAxis});
    histos.add("event/evtime/tofreso", "TOF event time reso.", kTH2F, {multAxis, colTimeResoAxis});
    histos.add("event/evtime/ft0", "FT0 event time", kTH1F, {colTimeAxis});
    histos.add("event/evtime/ft0reso", "FT0 event time reso.", kTH2F, {multAxis, colTimeResoAxis});

    histos.add("event/tofsignal", "TOF signal", kTH2F, {pAxis, tofAxis});
    histos.add("event/tofsignalunassigned", "TOF signal (unassigned tracks)", kTH2F, {pAxis, tofAxis});
    histos.add("event/pexp", "", kTH2F, {pAxis, pExpAxis});
    histos.add("event/eta", "", kTH1F, {etaAxis});
    histos.add("event/phi", "", kTH1F, {phiAxis});
    histos.add("event/etaphi", "", kTH2F, {etaAxis, phiAxis});
    histos.add("event/length", "", kTH1F, {lAxis});
    histos.add("event/pt", "", kTH1F, {ptAxis});
    histos.add("event/p", "", kTH1F, {pAxis});
    // histos.add("event/ptreso", "", kTH2F, {pAxis, ptResoAxis});

    static_for<0, 8>([&](auto i) {
      initPerParticle<i>(pAxis, ptAxis);
    });
  }

  template <bool fillHistograms, typename CollisionType, typename TrackType>
  bool isEventSelected(const CollisionType& collision, const TrackType& tracks)
  {

    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 1);
    }
    if (applyEvSel == 1) {
      if (!collision.sel7()) {
        return false;
      }
    } else if (applyEvSel == 2) {
      if (!collision.sel8()) {
        return false;
      }
    }

    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 2);
    }

    // Computing Multiplicity first
    int ntracks = 0;
    int tofmult = 0;
    float evtime = 0.f;
    float evtimereso = 0.f;
    int evtimeflag = 0;

    if constexpr (fillHistograms) {
      for (auto t : tracks) {
        if (applyTrackCut && !t.isGlobalTrack()) {
          continue;
        }
        ntracks += 1;
        if (!t.hasTOF()) { // Skipping tracks without TOF
          continue;
        }
        tofmult++;
        evtime = t.tofEvTime();
        evtimereso = t.tofEvTimeErr();
        evtimeflag = 0;
        if (t.isEvTimeDefined()) {
          evtimeflag = 1;
        }
        if (t.isEvTimeTOF() && t.isEvTimeT0AC()) {
          evtimeflag = 2;
        } else if (t.isEvTimeTOF()) {
          evtimeflag = 3;
        } else if (t.isEvTimeT0AC()) {
          evtimeflag = 4;
        }
      }
      histos.fill(HIST("event/evsel"), 3);
    }
    if (abs(collision.posZ()) > 10.f) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 4);
      histos.fill(HIST("event/vertexz"), collision.posZ());
      histos.fill(HIST("event/trackmultiplicity"), ntracks);
      histos.fill(HIST("event/tofmultiplicity"), tofmult);

      histos.fill(HIST("event/evtime/colltime"), collision.collisionTime() * 1000.f);
      histos.fill(HIST("event/evtime/colltimereso"), tofmult, collision.collisionTimeRes() * 1000.f);

      switch (evtimeflag) {
        case 0:
          histos.fill(HIST("event/evtime/undef"), evtime);
          histos.fill(HIST("event/evtime/undefreso"), tofmult, evtimereso);
          break;
        case 1:
          histos.fill(HIST("event/evtime/avail"), evtime);
          histos.fill(HIST("event/evtime/availreso"), tofmult, evtimereso);
          break;
        case 2:
          histos.fill(HIST("event/evtime/ft0tof"), evtime);
          histos.fill(HIST("event/evtime/ft0tofreso"), tofmult, evtimereso);
          break;
        case 3:
          histos.fill(HIST("event/evtime/tof"), evtime);
          histos.fill(HIST("event/evtime/tofreso"), tofmult, evtimereso);
          break;
        case 4:
          histos.fill(HIST("event/evtime/tof"), evtime);
          histos.fill(HIST("event/evtime/tofreso"), tofmult, evtimereso);
          break;
        default:
          LOG(fatal) << "Unrecognized Event time flag";
          break;
      }
    }
    return true;
  }

  template <bool fillHistograms, typename CollisionType, typename TrackType>
  bool isTrackSelected(const CollisionType& collision, const TrackType& track)
  {
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 1.f);
    }
    if (!track.isGlobalTrack()) { // Skipping non global tracks
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 2.f);
    }
    if (!track.hasITS()) { // Skipping tracks without ITS
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 3.f);
    }
    if (!track.hasTPC()) { // Skipping tracks without TPC
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 4.f);
    }
    if (!track.hasTOF()) { // Skipping tracks without TOF
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 5.f);
      histos.fill(HIST("event/particlehypo"), track.pidForTracking());
      if (track.has_collision()) {
        histos.fill(HIST("event/tofsignal"), track.p(), track.tofSignal());
      } else {
        histos.fill(HIST("event/tofsignalunassigned"), track.p(), track.tofSignal());
      }
      histos.fill(HIST("event/pexp"), track.p(), track.tofExpMom());
      histos.fill(HIST("event/eta"), track.eta());
      histos.fill(HIST("event/phi"), track.phi());
      histos.fill(HIST("event/etaphi"), track.eta(), track.phi());
      histos.fill(HIST("event/length"), track.length());
      histos.fill(HIST("event/pt"), track.pt());
      histos.fill(HIST("event/p"), track.p());
      // histos.fill(HIST("event/ptreso"), track.p(), track.sigma1Pt() * track.pt() * track.pt());
    }
    return true;
  }

  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  void process(CollisionCandidate const& collision,
               soa::Join<aod::Tracks, aod::TracksExtra,
                         aod::TOFSignal, aod::TOFEvTime, aod::pidEvTimeFlags,
                         aod::TrackSelection> const& tracks)
  {
    isEventSelected<true>(collision, tracks);
    for (auto t : tracks) {
      isTrackSelected<true>(collision, t);
    }
  }

  template <o2::track::PID::ID id, bool fillFullHistograms,
            typename TrackType>
  void processSingleParticle(CollisionCandidate const& collision,
                             TrackType const& tracks)
  {
    if (!isEventSelected<false>(collision, tracks)) {
      return;
    }

    for (auto t : tracks) {
      if (!isTrackSelected<false>(collision, t)) {
        continue;
      }

      if (applyRapidityCut) {
        if (abs(t.rapidity(PID::getMass(id))) > 0.5) {
          continue;
        }
      }

      const auto nsigma = o2::aod::pidutils::tofNSigma<id>(t);
      histos.fill(HIST(hnsigma[id]), t.p(), nsigma);
      histos.fill(HIST(hnsigma_pt[id]), t.pt(), nsigma);
      if (t.sign() > 0) {
        histos.fill(HIST(hnsigma_pt_pos[id]), t.pt(), nsigma);
      } else {
        histos.fill(HIST(hnsigma_pt_neg[id]), t.pt(), nsigma);
      }
      // Filling info split per ev. time
      if (enableEvTimeSplitting) {
        if (t.isEvTimeTOF() && t.isEvTimeT0AC()) { // TOF + FT0 Ev. Time
          histos.fill(HIST(hnsigma_evtime_tofft0[id]), t.p(), nsigma);
          histos.fill(HIST(hnsigma_pt_evtime_tofft0[id]), t.pt(), nsigma);
          if (t.sign() > 0) {
            histos.fill(HIST(hnsigma_pt_pos_evtime_tofft0[id]), t.pt(), nsigma);
          } else {
            histos.fill(HIST(hnsigma_pt_neg_evtime_tofft0[id]), t.pt(), nsigma);
          }
        } else if (t.isEvTimeT0AC()) { // FT0 Ev. Time
          histos.fill(HIST(hnsigma_evtime_ft0[id]), t.p(), nsigma);
          histos.fill(HIST(hnsigma_pt_evtime_ft0[id]), t.pt(), nsigma);
          if (t.sign() > 0) {
            histos.fill(HIST(hnsigma_pt_pos_evtime_ft0[id]), t.pt(), nsigma);
          } else {
            histos.fill(HIST(hnsigma_pt_neg_evtime_ft0[id]), t.pt(), nsigma);
          }
        } else if (t.isEvTimeTOF()) { // TOF Ev. Time
          histos.fill(HIST(hnsigma_evtime_tof[id]), t.p(), nsigma);
          histos.fill(HIST(hnsigma_pt_evtime_tof[id]), t.pt(), nsigma);
          if (t.sign() > 0) {
            histos.fill(HIST(hnsigma_pt_pos_evtime_tof[id]), t.pt(), nsigma);
          } else {
            histos.fill(HIST(hnsigma_pt_neg_evtime_tof[id]), t.pt(), nsigma);
          }
        } else { // No Ev. Time -> Fill Ev. Time
          histos.fill(HIST(hnsigma_evtime_fill[id]), t.p(), nsigma);
          histos.fill(HIST(hnsigma_pt_evtime_fill[id]), t.pt(), nsigma);
          if (t.sign() > 0) {
            histos.fill(HIST(hnsigma_pt_pos_evtime_fill[id]), t.pt(), nsigma);
          } else {
            histos.fill(HIST(hnsigma_pt_neg_evtime_fill[id]), t.pt(), nsigma);
          }
        }
      }
      if constexpr (fillFullHistograms) {
        const float tof = t.tofSignal() - t.tofEvTime();
        const auto diff = o2::aod::pidutils::tofExpSignalDiff<id>(t);
        histos.fill(HIST(hexpected[id]), t.p(), tof - diff);
        histos.fill(HIST(hdelta[id]), t.p(), diff);
        if (t.sign() > 0) {
          histos.fill(HIST(hdelta_pt_pos[id]), t.p(), diff);
        } else {
          histos.fill(HIST(hdelta_pt_neg[id]), t.p(), diff);
        }
        // Filling info split per ev. time
        if (enableEvTimeSplitting) {
          if (t.isEvTimeTOF() && t.isEvTimeT0AC()) { // TOF + FT0 Ev. Time
            histos.fill(HIST(hdelta_evtime_tofft0[id]), t.p(), diff);
            if (t.sign() > 0) {
              histos.fill(HIST(hdelta_pt_pos_evtime_tofft0[id]), t.p(), diff);
            } else {
              histos.fill(HIST(hdelta_pt_neg_evtime_tofft0[id]), t.p(), diff);
            }
          } else if (t.isEvTimeT0AC()) { // FT0 Ev. Time
            histos.fill(HIST(hdelta_evtime_ft0[id]), t.p(), diff);
            if (t.sign() > 0) {
              histos.fill(HIST(hdelta_pt_pos_evtime_ft0[id]), t.p(), diff);
            } else {
              histos.fill(HIST(hdelta_pt_neg_evtime_ft0[id]), t.p(), diff);
            }
          } else if (t.isEvTimeTOF()) { // TOF Ev. Time
            histos.fill(HIST(hdelta_evtime_tof[id]), t.p(), diff);
            if (t.sign() > 0) {
              histos.fill(HIST(hdelta_pt_pos_evtime_tof[id]), t.p(), diff);
            } else {
              histos.fill(HIST(hdelta_pt_neg_evtime_tof[id]), t.p(), diff);
            }
          } else { // No Ev. Time -> Fill Ev. Time
            histos.fill(HIST(hdelta_evtime_fill[id]), t.p(), diff);
            if (t.sign() > 0) {
              histos.fill(HIST(hdelta_pt_pos_evtime_fill[id]), t.p(), diff);
            } else {
              histos.fill(HIST(hdelta_pt_neg_evtime_fill[id]), t.p(), diff);
            }
          }
        }
        histos.fill(HIST(hexpsigma[id]), t.p(), o2::aod::pidutils::tofExpSigma<id>(t));
      }
    }
  }

  // QA of nsigma only tables
#define makeProcessFunction(inputPid, particleId)                                         \
  void process##particleId(CollisionCandidate const& collision,                           \
                           soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,  \
                                     aod::pidEvTimeFlags, aod::TOFSignal, aod::TOFEvTime, \
                                     inputPid> const& tracks)                             \
  {                                                                                       \
    processSingleParticle<PID::particleId, false>(collision, tracks);                     \
  }                                                                                       \
  PROCESS_SWITCH(tofPidQa, process##particleId, Form("Process for the %s hypothesis for TOF NSigma QA", #particleId), false);

  makeProcessFunction(aod::pidTOFEl, Electron);
  makeProcessFunction(aod::pidTOFMu, Muon);
  makeProcessFunction(aod::pidTOFPi, Pion);
  makeProcessFunction(aod::pidTOFKa, Kaon);
  makeProcessFunction(aod::pidTOFPr, Proton);
  makeProcessFunction(aod::pidTOFDe, Deuteron);
  makeProcessFunction(aod::pidTOFTr, Triton);
  makeProcessFunction(aod::pidTOFHe, Helium3);
  makeProcessFunction(aod::pidTOFAl, Alpha);
#undef makeProcessFunction

// QA of full tables
#define makeProcessFunction(inputPid, particleId)                                             \
  void processFull##particleId(CollisionCandidate const& collision,                           \
                               soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,  \
                                         aod::pidEvTimeFlags, aod::TOFSignal, aod::TOFEvTime, \
                                         inputPid> const& tracks)                             \
  {                                                                                           \
    processSingleParticle<PID::particleId, true>(collision, tracks);                          \
  }                                                                                           \
  PROCESS_SWITCH(tofPidQa, processFull##particleId, Form("Process for the %s hypothesis for full TOF PID QA", #particleId), false);

  makeProcessFunction(aod::pidTOFFullEl, Electron);
  makeProcessFunction(aod::pidTOFFullMu, Muon);
  makeProcessFunction(aod::pidTOFFullPi, Pion);
  makeProcessFunction(aod::pidTOFFullKa, Kaon);
  makeProcessFunction(aod::pidTOFFullPr, Proton);
  makeProcessFunction(aod::pidTOFFullDe, Deuteron);
  makeProcessFunction(aod::pidTOFFullTr, Triton);
  makeProcessFunction(aod::pidTOFFullHe, Helium3);
  makeProcessFunction(aod::pidTOFFullAl, Alpha);
#undef makeProcessFunction
};

/// Task to produce the TOF QA plots for Beta
struct tofPidBetaQa {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr std::string_view hexpected[Np] = {"expected/El", "expected/Mu", "expected/Pi",
                                                     "expected/Ka", "expected/Pr", "expected/De",
                                                     "expected/Tr", "expected/He", "expected/Al"};
  static constexpr std::string_view hdelta[Np] = {"delta/El", "delta/Mu", "delta/Pi",
                                                  "delta/Ka", "delta/Pr", "delta/De",
                                                  "delta/Tr", "delta/He", "delta/Al"};
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi",
                                                   "nsigma/Ka", "nsigma/Pr", "nsigma/De",
                                                   "nsigma/Tr", "nsigma/He", "nsigma/Al"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> logAxis{"logAxis", 0, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 400, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.1f, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 5.f, "Maximum momentum in range"};
  Configurable<int> nBinsBeta{"nBinsBeta", 4000, "Number of bins for the beta"};
  Configurable<float> minBeta{"minBeta", 0, "Minimum beta in range"};
  Configurable<float> maxBeta{"maxBeta", 2.f, "Maximum beta in range"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply event selection cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<bool> applyTrackCut{"applyTrackCut", false, "Flag to apply standard track cuts"};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec tofAxis{10000, 0, 2e6, "TOF Signal"};
    const AxisSpec betaAxis{nBinsBeta, minBeta, maxBeta, "TOF #beta"};
    const AxisSpec massAxis{1000, 0, 3, "TOF mass (GeV/#it{c}^{2})"};
    const AxisSpec etaAxis{100, -2, 2, "#it{#eta}"};
    const AxisSpec colTimeAxis{100, -2000, 2000, "Collision time (ps)"};
    const AxisSpec lAxis{100, 0, 500, "Track length (cm)"};
    const AxisSpec ptResoAxis{100, 0, 0.1, "#sigma_{#it{p}_{T}}"};
    const AxisSpec pAxisPosNeg{2 * nBinsP, -maxP, maxP, "#it{p}/z (GeV/#it{c})"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogarithmic();
      pAxis.makeLogarithmic();
    }

    // Event properties
    histos.add("event/tofsignal", "", HistType::kTH2F, {pAxis, tofAxis});
    histos.add("event/tofmass", "TOF mass", HistType::kTH1F, {massAxis});
    histos.add("event/tofmassEvTimeTOF", "TOF mass Ev. Time TOF", HistType::kTH1F, {massAxis});
    histos.add("event/tofmassEvTimeTOFOnly", "TOF mass Ev. Time TOF Only", HistType::kTH1F, {massAxis});
    histos.add("event/tofmassEvTimeT0AC", "TOF mass Ev. Time T0AC", HistType::kTH1F, {massAxis});
    histos.add("event/tofmassEvTimeT0ACOnly", "TOF mass Ev. Time T0AC Only", HistType::kTH1F, {massAxis});
    histos.add("event/tofbeta", "", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("event/tofbetaEvTimeTOF", "Ev. Time TOF", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("event/tofbetaEvTimeTOFOnly", "Ev. Time TOF Only", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("event/tofbetaEvTimeT0AC", "Ev. Time T0AC", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("event/tofbetaEvTimeT0ACOnly", "Ev. Time T0AC Only", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("event/signedtofbeta", "", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("event/signedtofbetaEvTimeTOF", "Ev. Time TOF", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("event/signedtofbetaEvTimeTOFOnly", "Ev. Time TOF Only", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("event/signedtofbetaEvTimeT0AC", "Ev. Time T0AC", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("event/signedtofbetaEvTimeT0ACOnly", "Ev. Time T0AC Only", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("event/eta", "", HistType::kTH1F, {etaAxis});
    histos.add("event/length", "", HistType::kTH1F, {lAxis});
    histos.add("event/pt", "", HistType::kTH1F, {ptAxis});
    histos.add("event/p", "", HistType::kTH1F, {pAxis});
    auto h = histos.add<TH1>("event/evsel", "", kTH1F, {{10, 0.5, 10.5, "Ev. Sel."}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Passed ev. sel.");
    h->GetXaxis()->SetBinLabel(3, "Passed mult.");
    h->GetXaxis()->SetBinLabel(4, "Passed vtx Z");

    h = histos.add<TH1>("event/trackselection", "", kTH1F, {{10, 0.5, 10.5, "Selection passed"}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "hasTOF");
    h->GetXaxis()->SetBinLabel(3, "isGlobalTrack");
  }

  template <uint8_t i, typename T>
  void fillParticleHistos(const T& t, const float tof, const float exp_diff, const float nsigma)
  {
    histos.fill(HIST(hexpected[i]), t.p(), tof - exp_diff);
    histos.fill(HIST(hdelta[i]), t.p(), exp_diff);
    histos.fill(HIST(hnsigma[i]), t.p(), nsigma);
  }

  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  void process(CollisionCandidate const& collision,
               soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTOFbeta, aod::pidTOFmass, aod::TrackSelection, aod::TOFSignal, aod::pidEvTimeFlags> const& tracks)
  {

    histos.fill(HIST("event/evsel"), 1);
    if (applyEvSel == 1) {
      if (!collision.sel7()) {
        return;
      }
    } else if (applyEvSel == 2) {
      if (!collision.sel8()) {
        return;
      }
    }

    histos.fill(HIST("event/evsel"), 2);

    // Computing Multiplicity first
    float ntracks = 0;
    for (auto t : tracks) {
      if (applyTrackCut && !t.isGlobalTrack()) {
        continue;
      }
      ntracks += 1;
    }
    histos.fill(HIST("event/evsel"), 3);
    if (abs(collision.posZ()) > 10.f) {
      return;
    }

    histos.fill(HIST("event/evsel"), 4);

    for (auto const& track : tracks) {
      histos.fill(HIST("event/trackselection"), 1.f);
      if (!track.hasTOF()) { // Skipping tracks without TOF
        continue;
      }
      histos.fill(HIST("event/trackselection"), 2.f);
      if (!track.isGlobalTrack()) {
        continue;
      }
      histos.fill(HIST("event/trackselection"), 3.f);
      if (track.isEvTimeTOF()) {
        histos.fill(HIST("event/tofmassEvTimeTOF"), track.mass());
        histos.fill(HIST("event/tofbetaEvTimeTOF"), track.p(), track.beta());
        histos.fill(HIST("event/signedtofbetaEvTimeTOF"), track.p() * track.sign(), track.beta());
      }
      if (track.isEvTimeTOF() && !track.isEvTimeT0AC()) {
        histos.fill(HIST("event/tofmassEvTimeTOFOnly"), track.mass());
        histos.fill(HIST("event/tofbetaEvTimeTOFOnly"), track.p(), track.beta());
        histos.fill(HIST("event/signedtofbetaEvTimeTOFOnly"), track.p() * track.sign(), track.beta());
      }
      if (track.isEvTimeT0AC()) {
        histos.fill(HIST("event/tofmassEvTimeT0AC"), track.mass());
        histos.fill(HIST("event/tofbetaEvTimeT0AC"), track.p(), track.beta());
        histos.fill(HIST("event/signedtofbetaEvTimeT0AC"), track.p() * track.sign(), track.beta());
      }
      if (track.isEvTimeT0AC() && !track.isEvTimeTOF()) {
        histos.fill(HIST("event/tofmassEvTimeT0ACOnly"), track.mass());
        histos.fill(HIST("event/tofbetaEvTimeT0ACOnly"), track.p(), track.beta());
        histos.fill(HIST("event/signedtofbetaEvTimeT0ACOnly"), track.p() * track.sign(), track.beta());
      }
      histos.fill(HIST("event/tofmass"), track.mass());
      histos.fill(HIST("event/tofbeta"), track.p(), track.beta());
      histos.fill(HIST("event/signedtofbeta"), track.p() * track.sign(), track.beta());
      histos.fill(HIST("event/length"), track.length());
      histos.fill(HIST("event/eta"), track.eta());
      histos.fill(HIST("event/tofsignal"), track.p(), track.tofSignal());
      histos.fill(HIST("event/pt"), track.pt());
      histos.fill(HIST("event/p"), track.p());
    }
  }
};
