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
/// \author Alberto Caliva (alberto.caliva@cern.ch)
/// \since June 27, 2023

#include <TMath.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TDatabasePDG.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/DataTypes.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "ReconstructionDataFormats/DCA.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct AntimatterAbsorptionHMPID {

  // Histograms
  HistogramRegistry registryQC{"registryQC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryDA{"registryDA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Track Selection Parameters
  Configurable<float> pmin{"pmin", 0.1, "pmin"};
  Configurable<float> pmax{"pmax", 3.0, "pmax"};
  Configurable<float> etaMin{"etaMin", -0.8, "etaMin"};
  Configurable<float> etaMax{"etaMax", +0.8, "etaMax"};
  Configurable<float> phiMin{"phiMin", 0.0, "phiMin"};
  Configurable<float> phiMax{"phiMax", 2.0 * TMath::Pi(), "phiMax"};
  Configurable<float> nsigmaTPCMin{"nsigmaTPCMin", -3.0, "nsigmaTPCMin"};
  Configurable<float> nsigmaTPCMax{"nsigmaTPCMax", +3.0, "nsigmaTPCMax"};
  Configurable<float> nsigmaTOFMin{"nsigmaTOFMin", -3.0, "nsigmaTOFMin"};
  Configurable<float> nsigmaTOFMax{"nsigmaTOFMax", +3.5, "nsigmaTOFMax"};
  Configurable<float> minReqClusterITS{"minReqClusterITS", 4.0, "min number of clusters required in ITS"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 50.0f, "minTPCnClsFound"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.5f, "maxDCAxy"};
  Configurable<float> maxDCAz{"maxDCAz", 0.5f, "maxDCAz"};
  Configurable<bool> use_hmpid_mom{"use_hmpid_mom", true, "use hmpid momentum"};

  void init(o2::framework::InitContext&)
  {
    // Event Counter
    registryQC.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{10, 0, 10, "counter"}});

    // HMPID Maps
    registryDA.add("hmpidXYpos", "hmpidXYpos", HistType::kTH2F, {{300, 0.0, 300.0, "x_{MIP}"}, {300, 0.0, 300.0, "y_{MIP}"}});
    registryDA.add("hmpidXYneg", "hmpidXYneg", HistType::kTH2F, {{300, 0.0, 300.0, "x_{MIP}"}, {300, 0.0, 300.0, "y_{MIP}"}});
    registryDA.add("mom_corr_pos", "mom_corr_pos", HistType::kTH2F, {{150, 0.0, 3.0, "p_{vtx}"}, {150, 0.0, 3.0, "p_{hmpid}"}});
    registryDA.add("mom_corr_neg", "mom_corr_neg", HistType::kTH2F, {{150, 0.0, 3.0, "p_{vtx}"}, {150, 0.0, 3.0, "p_{hmpid}"}});

    registryDA.add("hmpidEtaPhiMomPos", "hmpidEtaPhiMomPos", HistType::kTH3F, {{29, 0.1, 3.0, "p (GeV/c)"}, {180, -0.9, 0.9, "#eta"}, {200, 0.0, TMath::Pi(), "#phi"}});
    registryDA.add("hmpidEtaPhiMomNeg", "hmpidEtaPhiMomNeg", HistType::kTH3F, {{29, 0.1, 3.0, "p (GeV/c)"}, {180, -0.9, 0.9, "#eta"}, {200, 0.0, TMath::Pi(), "#phi"}});

    // Pion Pos
    registryDA.add("incomingPi_Pos_8cm", "incomingPi_Pos_8cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("incomingPi_Pos_4cm", "incomingPi_Pos_4cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("survivingPi_Pos_8cm", "survivingPi_Pos_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("survivingPi_Pos_4cm", "survivingPi_Pos_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("Pi_Pos_Q_8cm", "Pi_Pos_Q_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("Pi_Pos_Q_4cm", "Pi_Pos_Q_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("Pi_Pos_ClsSize_8cm", "Pi_Pos_ClsSize_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("Pi_Pos_ClsSize_4cm", "Pi_Pos_ClsSize_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("Pi_Pos_momentum", "Pi_Pos_momentum", HistType::kTH2F, {{100, 0.0, 3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {100, 0.0, 3.0, "#it{p}_{mhpid} (GeV/#it{c})"}});

    // Pion Neg
    registryDA.add("incomingPi_Neg_8cm", "incomingPi_Neg_8cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("incomingPi_Neg_4cm", "incomingPi_Neg_4cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("survivingPi_Neg_8cm", "survivingPi_Neg_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("survivingPi_Neg_4cm", "survivingPi_Neg_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("Pi_Neg_Q_8cm", "Pi_Neg_Q_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("Pi_Neg_Q_4cm", "Pi_Neg_Q_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("Pi_Neg_ClsSize_8cm", "Pi_Neg_ClsSize_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("Pi_Neg_ClsSize_4cm", "Pi_Neg_ClsSize_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("Pi_Neg_momentum", "Pi_Neg_momentum", HistType::kTH2F, {{100, 0.0, 3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {100, 0.0, 3.0, "#it{p}_{mhpid} (GeV/#it{c})"}});

    // Kaon Pos
    registryDA.add("incomingKa_Pos_8cm", "incomingKa_Pos_8cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("incomingKa_Pos_4cm", "incomingKa_Pos_4cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("survivingKa_Pos_8cm", "survivingKa_Pos_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("survivingKa_Pos_4cm", "survivingKa_Pos_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("Ka_Pos_Q_8cm", "Ka_Pos_Q_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("Ka_Pos_Q_4cm", "Ka_Pos_Q_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("Ka_Pos_ClsSize_8cm", "Ka_Pos_ClsSize_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("Ka_Pos_ClsSize_4cm", "Ka_Pos_ClsSize_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("Ka_Pos_momentum", "Ka_Pos_momentum", HistType::kTH2F, {{100, 0.0, 3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {100, 0.0, 3.0, "#it{p}_{mhpid} (GeV/#it{c})"}});

    // Kaon Neg
    registryDA.add("incomingKa_Neg_8cm", "incomingKa_Neg_8cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("incomingKa_Neg_4cm", "incomingKa_Neg_4cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("survivingKa_Neg_8cm", "survivingKa_Neg_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("survivingKa_Neg_4cm", "survivingKa_Neg_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("Ka_Neg_Q_8cm", "Ka_Neg_Q_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("Ka_Neg_Q_4cm", "Ka_Neg_Q_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("Ka_Neg_ClsSize_8cm", "Ka_Neg_ClsSize_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("Ka_Neg_ClsSize_4cm", "Ka_Neg_ClsSize_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("Ka_Neg_momentum", "Ka_Neg_momentum", HistType::kTH2F, {{100, 0.0, 3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {100, 0.0, 3.0, "#it{p}_{mhpid} (GeV/#it{c})"}});

    // Proton Pos
    registryDA.add("incomingPr_Pos_8cm", "incomingPr_Pos_8cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("incomingPr_Pos_4cm", "incomingPr_Pos_4cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("survivingPr_Pos_8cm", "survivingPr_Pos_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("survivingPr_Pos_4cm", "survivingPr_Pos_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("Pr_Pos_Q_8cm", "Pr_Pos_Q_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("Pr_Pos_Q_4cm", "Pr_Pos_Q_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("Pr_Pos_ClsSize_8cm", "Pr_Pos_ClsSize_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("Pr_Pos_ClsSize_4cm", "Pr_Pos_ClsSize_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("Pr_Pos_momentum", "Pr_Pos_momentum", HistType::kTH2F, {{100, 0.0, 3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {100, 0.0, 3.0, "#it{p}_{mhpid} (GeV/#it{c})"}});

    // Proton Neg
    registryDA.add("incomingPr_Neg_8cm", "incomingPr_Neg_8cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("incomingPr_Neg_4cm", "incomingPr_Neg_4cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("survivingPr_Neg_8cm", "survivingPr_Neg_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("survivingPr_Neg_4cm", "survivingPr_Neg_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("Pr_Neg_Q_8cm", "Pr_Neg_Q_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("Pr_Neg_Q_4cm", "Pr_Neg_Q_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("Pr_Neg_ClsSize_8cm", "Pr_Neg_ClsSize_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("Pr_Neg_ClsSize_4cm", "Pr_Neg_ClsSize_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("Pr_Neg_momentum", "Pr_Neg_momentum", HistType::kTH2F, {{100, 0.0, 3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {100, 0.0, 3.0, "#it{p}_{mhpid} (GeV/#it{c})"}});

    // Deuteron Pos
    registryDA.add("incomingDe_Pos_8cm", "incomingDe_Pos_8cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("incomingDe_Pos_4cm", "incomingDe_Pos_4cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("survivingDe_Pos_8cm", "survivingDe_Pos_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("survivingDe_Pos_4cm", "survivingDe_Pos_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("De_Pos_Q_8cm", "De_Pos_Q_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("De_Pos_Q_4cm", "De_Pos_Q_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("De_Pos_ClsSize_8cm", "De_Pos_ClsSize_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("De_Pos_ClsSize_4cm", "De_Pos_ClsSize_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("De_Pos_momentum", "De_Pos_momentum", HistType::kTH2F, {{100, 0.0, 3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {100, 0.0, 3.0, "#it{p}_{mhpid} (GeV/#it{c})"}});

    // Deuteron Neg
    registryDA.add("incomingDe_Neg_8cm", "incomingDe_Neg_8cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("incomingDe_Neg_4cm", "incomingDe_Neg_4cm", HistType::kTH1F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("survivingDe_Neg_8cm", "survivingDe_Neg_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("survivingDe_Neg_4cm", "survivingDe_Neg_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("De_Neg_Q_8cm", "De_Neg_Q_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("De_Neg_Q_4cm", "De_Neg_Q_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("De_Neg_ClsSize_8cm", "De_Neg_ClsSize_8cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("De_Neg_ClsSize_4cm", "De_Neg_ClsSize_4cm", HistType::kTH2F, {{290, 0.1, 3.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Cls size"}});
    registryDA.add("De_Neg_momentum", "De_Neg_momentum", HistType::kTH2F, {{100, 0.0, 3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {100, 0.0, 3.0, "#it{p}_{mhpid} (GeV/#it{c})"}});
  }

  // Single-Track Selection
  template <typename trackType>
  bool passedTrackSelection(const trackType& track)
  {
    if (!track.hasITS())
      return false;
    if (!track.hasTPC())
      return false;
    if (!track.hasTOF())
      return false;
    if (!track.passedITSRefit())
      return false;
    if (!track.passedTPCRefit())
      return false;
    if (track.itsNCls() < minReqClusterITS)
      return false;
    if (track.tpcNClsFound() < minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.itsChi2NCl() > maxChi2ITS)
      return false;
    if (TMath::Abs(track.dcaXY()) > maxDCAxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > maxDCAz)
      return false;
    /*
      if (track.eta() < etaMin)
      return false;
    if (track.eta() > etaMax)
      return false;
    if (track.phi() < phiMin)
      return false;
    if (track.phi() > phiMax)
      return false;*/

    return true;
  }

  // Particle Identification (Pions)
  template <typename pionCandidate>
  bool passedPionSelection(const pionCandidate& track)
  {
    if (track.tpcNSigmaPi() < nsigmaTPCMin)
      return false;
    if (track.tpcNSigmaPi() > nsigmaTPCMax)
      return false;
    if (track.tofNSigmaPi() < nsigmaTOFMin)
      return false;
    if (track.tofNSigmaPi() > nsigmaTOFMax)
      return false;

    return true;
  }

  // Particle Identification (Kaons)
  template <typename kaonCandidate>
  bool passedKaonSelection(const kaonCandidate& track)
  {
    if (track.tpcNSigmaKa() < nsigmaTPCMin)
      return false;
    if (track.tpcNSigmaKa() > nsigmaTPCMax)
      return false;
    if (track.tofNSigmaKa() < nsigmaTOFMin)
      return false;
    if (track.tofNSigmaKa() > nsigmaTOFMax)
      return false;

    return true;
  }

  // Particle Identification (Protons)
  template <typename protonCandidate>
  bool passedProtonSelection(const protonCandidate& track)
  {
    if (track.tpcNSigmaPr() < nsigmaTPCMin)
      return false;
    if (track.tpcNSigmaPr() > nsigmaTPCMax)
      return false;
    if (track.tofNSigmaPr() < nsigmaTOFMin)
      return false;
    if (track.tofNSigmaPr() > nsigmaTOFMax)
      return false;

    return true;
  }

  // Particle Identification (Deuterons)
  template <typename deuteronCandidate>
  bool passedDeuteronSelection(const deuteronCandidate& track)
  {
    if (track.tpcNSigmaDe() < nsigmaTPCMin)
      return false;
    if (track.tpcNSigmaDe() > nsigmaTPCMax)
      return false;
    if (track.tofNSigmaDe() < nsigmaTOFMin)
      return false;
    if (track.tofNSigmaDe() > nsigmaTOFMax)
      return false;

    return true;
  }

  // Full Tracks
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe>;

  // Process Data
  void processData(o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator const& event, o2::aod::HMPIDs const& hmpids)
  {
    // Event Selection
    registryQC.fill(HIST("number_of_events_data"), 0.5);

    if (!event.sel8())
      return;

    // Event Counter
    registryQC.fill(HIST("number_of_events_data"), 1.5);

    if (abs(event.posZ()) > 10.0)
      return;

    // Event Counter
    registryQC.fill(HIST("number_of_events_data"), 2.5);

    for (const auto& hmpid : hmpids) {

      // Get Track
      const auto& track = hmpid.track_as<FullTracks>();

      // Track Momentum
      float momentum = track.p();
      if (use_hmpid_mom)
        momentum = hmpid.hmpidMom();

      // Track Selection
      if (!passedTrackSelection(track))
        continue;

      if (track.sign() > 0) {
        registryDA.fill(HIST("hmpidXYpos"), hmpid.hmpidXMip(), hmpid.hmpidYMip());
        registryDA.fill(HIST("hmpidEtaPhiMomPos"), track.p(), track.eta(), track.phi());
        registryDA.fill(HIST("mom_corr_pos"), track.p(), hmpid.hmpidMom());
      }
      if (track.sign() < 0) {
        registryDA.fill(HIST("hmpidXYneg"), hmpid.hmpidXMip(), hmpid.hmpidYMip());
        registryDA.fill(HIST("hmpidEtaPhiMomNeg"), track.p(), track.eta(), track.phi());
        registryDA.fill(HIST("mom_corr_neg"), track.p(), hmpid.hmpidMom());
      }

      // Absorber
      bool hmpidAbs8cm = true;
      bool hmpidAbs4cm = true;

      // Distance between extrapolated and matched point
      float dx = hmpid.hmpidXTrack() - hmpid.hmpidXMip();
      float dy = hmpid.hmpidYTrack() - hmpid.hmpidYMip();
      float dr = sqrt(dx * dx + dy * dy);

      // Fill Histograms for Positive Pions
      if (passedPionSelection(track) && track.sign() > 0) {

        if (hmpidAbs8cm) {
          registryDA.fill(HIST("incomingPi_Pos_8cm"), momentum);
          registryDA.fill(HIST("survivingPi_Pos_8cm"), momentum, dr);
          registryDA.fill(HIST("Pi_Pos_Q_8cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("Pi_Pos_ClsSize_8cm"), momentum, hmpid.hmpidClusSize());
        }
        if (hmpidAbs4cm) {
          registryDA.fill(HIST("incomingPi_Pos_4cm"), momentum);
          registryDA.fill(HIST("survivingPi_Pos_4cm"), momentum, dr);
          registryDA.fill(HIST("Pi_Pos_Q_4cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("Pi_Pos_ClsSize_4cm"), momentum, hmpid.hmpidClusSize());
        }
      }

      // Fill Histograms for Negative Pions
      if (passedPionSelection(track) && track.sign() < 0) {

        if (hmpidAbs8cm) {
          registryDA.fill(HIST("incomingPi_Neg_8cm"), momentum);
          registryDA.fill(HIST("survivingPi_Neg_8cm"), momentum, dr);
          registryDA.fill(HIST("Pi_Neg_Q_8cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("Pi_Neg_ClsSize_8cm"), momentum, hmpid.hmpidClusSize());
        }
        if (hmpidAbs4cm) {
          registryDA.fill(HIST("incomingPi_Neg_4cm"), momentum);
          registryDA.fill(HIST("survivingPi_Neg_4cm"), momentum, dr);
          registryDA.fill(HIST("Pi_Neg_Q_4cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("Pi_Neg_ClsSize_4cm"), momentum, hmpid.hmpidClusSize());
        }
      }

      // Fill Histograms for Positive Kaons
      if (passedKaonSelection(track) && track.sign() > 0) {

        if (hmpidAbs8cm) {
          registryDA.fill(HIST("incomingKa_Pos_8cm"), momentum);
          registryDA.fill(HIST("survivingKa_Pos_8cm"), momentum, dr);
          registryDA.fill(HIST("Ka_Pos_Q_8cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("Ka_Pos_ClsSize_8cm"), momentum, hmpid.hmpidClusSize());
        }
        if (hmpidAbs4cm) {
          registryDA.fill(HIST("incomingKa_Pos_4cm"), momentum);
          registryDA.fill(HIST("survivingKa_Pos_4cm"), momentum, dr);
          registryDA.fill(HIST("Ka_Pos_Q_4cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("Ka_Pos_ClsSize_4cm"), momentum, hmpid.hmpidClusSize());
        }
      }

      // Fill Histograms for Negative Kaons
      if (passedKaonSelection(track) && track.sign() < 0) {

        if (hmpidAbs8cm) {
          registryDA.fill(HIST("incomingKa_Neg_8cm"), momentum);
          registryDA.fill(HIST("survivingKa_Neg_8cm"), momentum, dr);
          registryDA.fill(HIST("Ka_Neg_Q_8cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("Ka_Neg_ClsSize_8cm"), momentum, hmpid.hmpidClusSize());
        }
        if (hmpidAbs4cm) {
          registryDA.fill(HIST("incomingKa_Neg_4cm"), momentum);
          registryDA.fill(HIST("survivingKa_Neg_4cm"), momentum, dr);
          registryDA.fill(HIST("Ka_Neg_Q_4cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("Ka_Neg_ClsSize_4cm"), momentum, hmpid.hmpidClusSize());
        }
      }

      // Fill Histograms for Positive Protons
      if (passedProtonSelection(track) && track.sign() > 0) {

        if (hmpidAbs8cm) {
          registryDA.fill(HIST("incomingPr_Pos_8cm"), momentum);
          registryDA.fill(HIST("survivingPr_Pos_8cm"), momentum, dr);
          registryDA.fill(HIST("Pr_Pos_Q_8cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("Pr_Pos_ClsSize_8cm"), momentum, hmpid.hmpidClusSize());
        }
        if (hmpidAbs4cm) {
          registryDA.fill(HIST("incomingPr_Pos_4cm"), momentum);
          registryDA.fill(HIST("survivingPr_Pos_4cm"), momentum, dr);
          registryDA.fill(HIST("Pr_Pos_Q_4cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("Pr_Pos_ClsSize_4cm"), momentum, hmpid.hmpidClusSize());
        }
      }

      // Fill Histograms for Negative Protons
      if (passedProtonSelection(track) && track.sign() < 0) {

        if (hmpidAbs8cm) {
          registryDA.fill(HIST("incomingPr_Neg_8cm"), momentum);
          registryDA.fill(HIST("survivingPr_Neg_8cm"), momentum, dr);
          registryDA.fill(HIST("Pr_Neg_Q_8cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("Pr_Neg_ClsSize_8cm"), momentum, hmpid.hmpidClusSize());
        }
        if (hmpidAbs4cm) {
          registryDA.fill(HIST("incomingPr_Neg_4cm"), momentum);
          registryDA.fill(HIST("survivingPr_Neg_4cm"), momentum, dr);
          registryDA.fill(HIST("Pr_Neg_Q_4cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("Pr_Neg_ClsSize_4cm"), momentum, hmpid.hmpidClusSize());
        }
      }

      // Fill Histograms for Positive Deuterons
      if (passedDeuteronSelection(track) && track.sign() > 0) {

        if (hmpidAbs8cm) {
          registryDA.fill(HIST("incomingDe_Pos_8cm"), momentum);
          registryDA.fill(HIST("survivingDe_Pos_8cm"), momentum, dr);
          registryDA.fill(HIST("De_Pos_Q_8cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("De_Pos_ClsSize_8cm"), momentum, hmpid.hmpidClusSize());
        }
        if (hmpidAbs4cm) {
          registryDA.fill(HIST("incomingDe_Pos_4cm"), momentum);
          registryDA.fill(HIST("survivingDe_Pos_4cm"), momentum, dr);
          registryDA.fill(HIST("De_Pos_Q_4cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("De_Pos_ClsSize_4cm"), momentum, hmpid.hmpidClusSize());
        }
      }

      // Fill Histograms for Negative Deuterons
      if (passedDeuteronSelection(track) && track.sign() < 0) {

        if (hmpidAbs8cm) {
          registryDA.fill(HIST("incomingDe_Neg_8cm"), momentum);
          registryDA.fill(HIST("survivingDe_Neg_8cm"), momentum, dr);
          registryDA.fill(HIST("De_Neg_Q_8cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("De_Neg_ClsSize_8cm"), momentum, hmpid.hmpidClusSize());
        }
        if (hmpidAbs4cm) {
          registryDA.fill(HIST("incomingDe_Neg_4cm"), momentum);
          registryDA.fill(HIST("survivingDe_Neg_4cm"), momentum, dr);
          registryDA.fill(HIST("De_Neg_Q_4cm"), momentum, hmpid.hmpidQMip());
          registryDA.fill(HIST("De_Neg_ClsSize_4cm"), momentum, hmpid.hmpidClusSize());
        }
      }
    }
  }
  PROCESS_SWITCH(AntimatterAbsorptionHMPID, processData, "process data", true);
};

//*************************************************************************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<AntimatterAbsorptionHMPID>(cfgc)};
}
