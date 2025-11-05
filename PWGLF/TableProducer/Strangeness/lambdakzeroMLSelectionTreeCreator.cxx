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
//
//  *+-+*+-+*+-+*+-+*+-+*+-+*
//  V0 ML Tree Creator task
//  *+-+*+-+*+-+*+-+*+-+*+-+*
//
//  This task loops over a set of V0 indices and
//  creates a TTree for ML training and testing.
//
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    gianni.shigeru.setoue.liveraro@cern.ch
//    romain.schotter@cern.ch
//    david.dobrigkeit.chinellato@cern.ch
//

#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h> // C system
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>   // C++ system
#include <cmath>   // C++ system
#include <cstdlib> // C++ system

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using std::cout;
using std::endl;

using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

struct lambdakzeroMLSelectionTreeCreator {
  Produces<aod::V0MLCandidates> v0MLCandidates;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Base selection criteria

  // Lambda standard criteria::
  Configurable<float> Lambdav0cospa{"Lambdav0cospa", 0.90, "min V0 CosPA"};
  Configurable<float> Lambdadcav0dau{"Lambdadcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
  Configurable<float> Lambdadcanegtopv{"Lambdadcanegtopv", .05, "min DCA Neg To PV (cm)"};
  Configurable<float> Lambdadcapostopv{"Lambdadcapostopv", .05, "min DCA Pos To PV (cm)"};
  Configurable<float> Lambdav0radius{"Lambdav0radius", 1.5, "minimum V0 radius (cm)"};
  Configurable<float> LambdaWindow{"LambdaWindow", 0.01, "Mass window around expected (in GeV/c2)"};

  // Anti-Lambda standard criteria::
  Configurable<float> AntiLambdav0cospa{"AntiLambdav0cospa", 0.90, "min V0 CosPA"};
  Configurable<float> AntiLambdadcav0dau{"AntiLambdadcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
  Configurable<float> AntiLambdadcanegtopv{"AntiLambdadcanegtopv", .05, "min DCA Neg To PV (cm)"};
  Configurable<float> AntiLambdadcapostopv{"AntiLambdadcapostopv", .05, "min DCA Pos To PV (cm)"};
  Configurable<float> AntiLambdav0radius{"AntiLambdav0radius", 1.5, "minimum V0 radius (cm)"};
  Configurable<float> AntiLambdaWindow{"AntiLambdaWindow", 0.01, "Mass window around expected (in GeV/c2)"};

  // Photon standard criteria:
  Configurable<float> PhotonMinRadius{"PhotonMinRadius", 1.0, "minimum photon conversion radius (cm)"};
  Configurable<float> PhotonMaxRadius{"PhotonMaxRadius", 200, "maximum photon conversion radius (cm)"};
  Configurable<float> PhotonMinPt{"PhotonMinPt", 0.001, "minimum photon pT (GeV/c)"};
  Configurable<float> PhotonMaxMass{"PhotonMaxMass", 0.2, "Maximum photon mass (GeV/c^{2})"};
  Configurable<float> PhotonMaxqt{"PhotonMaxqt", 0.1, "Maximum photon qt value (AP plot) (GeV/c)"};
  Configurable<float> PhotonMaxalpha{"PhotonMaxalpha", 1.0, "Max photon alpha absolute value (AP plot)"};
  Configurable<float> PhotonWindow{"PhotonWindow", 0.01, "Mass window around expected (in GeV/c2)"};

  // KZeroShort standard criteria:
  // TODO: update criteria
  Configurable<float> K0Shortv0cospa{"K0Shortv0cospa", 0.90, "min V0 CosPA"};
  Configurable<float> K0Shortdcav0dau{"K0Shortdcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
  Configurable<float> K0Shortdcanegtopv{"K0Shortdcanegtopv", .05, "min DCA Neg To PV (cm)"};
  Configurable<float> K0Shortdcapostopv{"K0Shortdcapostopv", .05, "min DCA Pos To PV (cm)"};
  Configurable<float> K0Shortv0radius{"K0Shortv0radius", 1.5, "minimum V0 radius (cm)"};
  Configurable<float> K0ShortWindow{"K0ShortWindow", 0.01, "Mass window around expected (in GeV/c2)"};

  Configurable<int> RejectHypothesis{"RejectHypothesis", -1, "SelHypothesis to reject"};
  Configurable<bool> SelMCCand{"SelMCCand", true, "Select candidate based on PDG Code"};
  Configurable<bool> SelMCCandMother{"SelMCCandMother", false, "Select candidate based on PDGCodeMother"};
  Configurable<int> SelPDGCode{"SelPDGCode", 22, "PDGCode to select candidate"};
  Configurable<int> SelPDGCodeMother{"SelPDGCodeMother", 3212, "PDGCodeMother to select candidate"};

  // Axis:
  ConfigurableAxis centralityAxis{"centralityAxis", {100, 0.0f, 100.0f}, ""};
  ConfigurableAxis candSelectionAxis{"candSelectionAxis", {16, 0.0f, 16.0f}, ""};

  void init(InitContext const&)
  {
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {centralityAxis});
    histos.add("hCandSelection", "hCandSelection", kTH1F, {candSelectionAxis});
  }

  // Helper struct to pass v0 information
  struct {
    int posITSCls;
    int negITSCls;
    uint32_t posITSClSize;
    uint32_t negITSClSize;
    uint8_t posTPCRows;
    uint8_t negTPCRows;
    float posTPCSigmaPi;
    float negTPCSigmaPi;
    float posTPCSigmaPr;
    float negTPCSigmaPr;
    float posTPCSigmaEl;
    float negTPCSigmaEl;
    float TOFSigmaLaPr;
    float TOFSigmaLaPi;
    float TOFSigmaALaPi;
    float TOFSigmaALaPr;
    float TOFSigmaK0PiPlus;
    float TOFSigmaK0PiMinus;
    float LambdaMass;
    float AntiLambdaMass;
    float GammaMass;
    float KZeroShortMass;
    float pT;
    float qt;
    float alpha;
    float posEta;
    float negEta;
    float v0Eta;
    float Z;
    float v0radius;
    float PA;
    float dcapostopv;
    float dcanegtopv;
    float dcaV0daughters;
    float dcav0topv;
    float PsiPair;
    uint8_t v0type;
    float centrality;
    uint8_t SelHypothesis;
    bool isLambda;
    bool isAntiLambda;
    bool isGamma;
    bool isKZeroShort;
    int PDGCodeMother;
  } Candidate;

  // Process candidate and store properties in object
  template <typename TV0Object, typename TCollision>
  void processCandidate(TCollision const& coll, TV0Object const& cand)
  {
    bool lConsistentWithLambda = true;
    bool lConsistentWithAntiLambda = true;
    bool lConsistentWithGamma = true;
    bool lConsistentWithK0Short = true;

    auto posTrackExtra = cand.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = cand.template negTrackExtra_as<dauTracks>();

    // Track quality
    Candidate.posITSCls = posTrackExtra.itsNCls();
    Candidate.negITSCls = negTrackExtra.itsNCls();
    Candidate.posITSClSize = posTrackExtra.itsClusterSizes();
    Candidate.negITSClSize = negTrackExtra.itsClusterSizes();
    Candidate.posTPCRows = posTrackExtra.tpcCrossedRows();
    Candidate.negTPCRows = negTrackExtra.tpcCrossedRows();

    // TPC PID
    Candidate.posTPCSigmaPi = posTrackExtra.tpcNSigmaPi();
    Candidate.negTPCSigmaPi = negTrackExtra.tpcNSigmaPi();
    Candidate.posTPCSigmaPr = posTrackExtra.tpcNSigmaPr();
    Candidate.negTPCSigmaPr = negTrackExtra.tpcNSigmaPr();
    Candidate.posTPCSigmaEl = posTrackExtra.tpcNSigmaEl();
    Candidate.negTPCSigmaEl = negTrackExtra.tpcNSigmaEl();

    // TOF PID:
    Candidate.TOFSigmaLaPr = cand.tofNSigmaLaPr();
    Candidate.TOFSigmaLaPi = cand.tofNSigmaLaPi();
    Candidate.TOFSigmaALaPi = cand.tofNSigmaALaPi();
    Candidate.TOFSigmaALaPr = cand.tofNSigmaALaPr();
    Candidate.TOFSigmaK0PiPlus = cand.tofNSigmaK0PiPlus();
    Candidate.TOFSigmaK0PiMinus = cand.tofNSigmaK0PiMinus();

    // General
    Candidate.LambdaMass = cand.mLambda();
    Candidate.AntiLambdaMass = cand.mAntiLambda();
    Candidate.GammaMass = cand.mGamma();
    Candidate.KZeroShortMass = cand.mK0Short();
    Candidate.pT = cand.pt();
    Candidate.qt = cand.qtarm();
    Candidate.alpha = cand.alpha();
    Candidate.posEta = cand.positiveeta();
    Candidate.negEta = cand.negativeeta();
    Candidate.v0Eta = cand.eta();
    Candidate.Z = cand.z();

    // Topological
    Candidate.v0radius = cand.v0radius();
    Candidate.PA = TMath::ACos(cand.v0cosPA());
    Candidate.dcapostopv = cand.dcapostopv();
    Candidate.dcanegtopv = cand.dcanegtopv();
    Candidate.dcaV0daughters = cand.dcaV0daughters();
    Candidate.dcav0topv = cand.dcav0topv();
    Candidate.PsiPair = cand.psipair();

    // Debug/Aditional
    Candidate.v0type = cand.v0Type();
    Candidate.centrality = coll.centFT0C();

    // Applying selections and saving hypothesis
    if ((std::abs(cand.mLambda() - 1.115683) > LambdaWindow) || (cand.v0radius() < Lambdav0radius) || (cand.v0cosPA() < Lambdav0cospa) || (TMath::Abs(cand.dcapostopv()) < Lambdadcapostopv) || (TMath::Abs(cand.dcanegtopv()) < Lambdadcanegtopv) || (cand.dcaV0daughters() > Lambdadcav0dau))
      lConsistentWithLambda = false;
    if ((std::abs(cand.mAntiLambda() - 1.115683) > AntiLambdaWindow) || (cand.v0radius() < AntiLambdav0radius) || (cand.v0cosPA() < AntiLambdav0cospa) || (TMath::Abs(cand.dcapostopv()) < AntiLambdadcapostopv) || (TMath::Abs(cand.dcanegtopv()) < AntiLambdadcanegtopv) || (cand.dcaV0daughters() > AntiLambdadcav0dau))
      lConsistentWithAntiLambda = false;
    if ((std::abs(cand.mGamma()) > PhotonWindow) || (cand.v0radius() < PhotonMinRadius) || (cand.v0radius() > PhotonMaxRadius) || (cand.pt() < PhotonMinPt) || (cand.qtarm() > PhotonMaxqt) || (TMath::Abs(cand.alpha()) > PhotonMaxalpha))
      lConsistentWithGamma = false;
    if ((std::abs(cand.mK0Short() - 0.497) > K0ShortWindow) || (cand.v0radius() < K0Shortv0radius) || (cand.v0cosPA() < K0Shortv0cospa) || (TMath::Abs(cand.dcapostopv()) < K0Shortdcapostopv) || (TMath::Abs(cand.dcanegtopv()) < K0Shortdcanegtopv) || (cand.dcaV0daughters() > K0Shortdcav0dau))
      lConsistentWithK0Short = false;

    Candidate.SelHypothesis = lConsistentWithLambda << 0 | lConsistentWithAntiLambda << 1 | lConsistentWithGamma << 2 | lConsistentWithK0Short << 3;
    // 1: Consistent with Lambda only, 2: Consistent with Anti-Lambda only
    // 3: Consistent with Lambda and Anti-Lambda, 4: Consistent with Gamma only
    // 5: Consistent with Lambda and Gamma, 6: Consistent with Anti-Lambda and Gamma
    // 7: Consistent with Lambda, Anti-Lambda, and Gamma, 8: Consistent with K0Short only
    // 9: Consistent with Lambda and K0Short, 10: Consistent with Anti-Lambda and K0Short
    // 11: Consistent with Lambda, Anti-Lambda, and K0Short, 12: Consistent with Gamma and K0Short
    // 13: Consistent with Lambda, Gamma, and K0Short, 14: Consistent with Anti-Lambda, Gamma, and K0Short
    // 15: Consistent with Lambda, Anti-Lambda, Gamma, and K0Short

    histos.fill(HIST("hCandSelection"), Candidate.SelHypothesis);

    if (Candidate.SelHypothesis == RejectHypothesis)
      return;

    // MC flags
    Candidate.isLambda = false;
    Candidate.isAntiLambda = false;
    Candidate.isGamma = false;
    Candidate.isKZeroShort = false;
    Candidate.PDGCodeMother = -1;

    if constexpr (requires { cand.pdgCode(); }) {
      if (SelMCCand && (cand.pdgCode() != SelPDGCode))
        return;

      Candidate.isLambda = (cand.pdgCode() == 3122);
      Candidate.isAntiLambda = (cand.pdgCode() == -3122);
      Candidate.isGamma = (cand.pdgCode() == 22);
      Candidate.isKZeroShort = (cand.pdgCode() == 310);
    }

    if constexpr (requires { cand.pdgCodeMother(); }) {
      if (SelMCCandMother && (cand.pdgCodeMother() != SelPDGCodeMother))
        return;

      Candidate.PDGCodeMother = cand.pdgCodeMother();
    }

    // Filling TTree for ML analysis
    v0MLCandidates(Candidate.posITSCls, Candidate.negITSCls, Candidate.posITSClSize, Candidate.negITSClSize, Candidate.posTPCRows, Candidate.negTPCRows,
                   Candidate.posTPCSigmaPi, Candidate.negTPCSigmaPi, Candidate.posTPCSigmaPr, Candidate.negTPCSigmaPr,
                   Candidate.posTPCSigmaEl, Candidate.negTPCSigmaEl, Candidate.TOFSigmaLaPr, Candidate.TOFSigmaLaPi,
                   Candidate.TOFSigmaALaPi, Candidate.TOFSigmaALaPr, Candidate.TOFSigmaK0PiPlus, Candidate.TOFSigmaK0PiMinus,
                   Candidate.LambdaMass, Candidate.AntiLambdaMass, Candidate.GammaMass, Candidate.KZeroShortMass, Candidate.pT,
                   Candidate.qt, Candidate.alpha, Candidate.posEta, Candidate.negEta, Candidate.v0Eta, Candidate.Z,
                   Candidate.v0radius, Candidate.PA, Candidate.dcapostopv, Candidate.dcanegtopv, Candidate.dcaV0daughters, Candidate.dcav0topv, Candidate.PsiPair,
                   Candidate.v0type, Candidate.centrality, Candidate.SelHypothesis, Candidate.isLambda, Candidate.isAntiLambda, Candidate.isGamma, Candidate.isKZeroShort, Candidate.PDGCodeMother);
  }

  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents>::iterator const& coll, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFNSigmas> const& v0s, dauTracks const&)
  {
    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    for (auto& cand : v0s) { // looping over lambdas
      processCandidate(coll, cand);
    }
  }

  void processSimData(soa::Join<aod::StraCollisions, aod::StraCents>::iterator const& coll, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0MCDatas, aod::V0TOFNSigmas> const& v0s, dauTracks const&)
  {
    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    for (auto& cand : v0s) { // looping over lambdas
      processCandidate(coll, cand);
    }
  }

  PROCESS_SWITCH(lambdakzeroMLSelectionTreeCreator, processRealData, "Produce Run 3 v0 tables (real data)", false);
  PROCESS_SWITCH(lambdakzeroMLSelectionTreeCreator, processSimData, "Produce Run 3 v0 tables (simulation)", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdakzeroMLSelectionTreeCreator>(cfgc)};
}
