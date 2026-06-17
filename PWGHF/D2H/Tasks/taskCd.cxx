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

/// \file taskCd.cxx
/// \brief cd± → d± K∓ π±  analysis task
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg Universiity

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>
#include <TPDGCode.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <string>
#include <string_view>
#include <vector> // std::vector

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::hf_evsel;

namespace o2::aod
{
namespace full
{
// Candidate kinematics
DECLARE_SOA_COLUMN(MassCd, massCd, float);                          //! Invariant mass of cd candidate (GeV/c^2)
DECLARE_SOA_COLUMN(MassLc, massLc, float);                          //! Invariant mass of lc candidate (GeV/c^2)
DECLARE_SOA_COLUMN(Pt, pt, float);                                  //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                                //! eta of candidate (GeV/c)
DECLARE_SOA_COLUMN(Phi, phi, float);                                //! phi of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                    //! rapidity of generated particle
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);                      //! Transverse momentum of prong 0 (GeV/c)
DECLARE_SOA_COLUMN(PxProng0, pxProng0, float);                      //! Px of prong 0 (GeV/c)
DECLARE_SOA_COLUMN(PyProng0, pyProng0, float);                      //! Py of prong 0 (GeV/c)
DECLARE_SOA_COLUMN(PzProng0, pzProng0, float);                      //! Pz of prong 0 (GeV/c)
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);                      //! Transverse momentum of prong 1 (GeV/c)
DECLARE_SOA_COLUMN(PxProng1, pxProng1, float);                      //! Px of prong 1 (GeV/c)
DECLARE_SOA_COLUMN(PyProng1, pyProng1, float);                      //! Py of prong 1 (GeV/c)
DECLARE_SOA_COLUMN(PzProng1, pzProng1, float);                      //! Pz of prong 1 (GeV/c)
DECLARE_SOA_COLUMN(PtProng2, ptProng2, float);                      //! Transverse momentum of prong 2 (GeV/c)
DECLARE_SOA_COLUMN(PxProng2, pxProng2, float);                      //! Px of prong 2 (GeV/c)
DECLARE_SOA_COLUMN(PyProng2, pyProng2, float);                      //! Py of prong 2 (GeV/c)
DECLARE_SOA_COLUMN(PzProng2, pzProng2, float);                      //! Pz of prong 2 (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameter0, impactParameter0, float);      //! Impact parameter (DCA to PV) of prong 0 (cm)
DECLARE_SOA_COLUMN(ImpactParameter1, impactParameter1, float);      //! Impact parameter (DCA to PV) of prong 1 (cm)
DECLARE_SOA_COLUMN(ImpactParameter2, impactParameter2, float);      //! Impact parameter (DCA to PV) of prong 2 (cm)
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                //! Decay length (3D) of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);            //! Decay length in transverse plane (cm)
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                //! Cosine of pointing angle (3D)
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                            //! Cosine of pointing angle in XY plane
DECLARE_SOA_COLUMN(Chi2PCA, chi2PCA, float);                        //! chi2PCA
DECLARE_SOA_COLUMN(NSigmaTpcDe, nSigmaTpcDe, float);                //! TPC nσ for deuteron hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcPr, nSigmaTpcPr, float);                //! TPC nσ for proton hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcKa, nSigmaTpcKa, float);                //! TPC nσ for kaon hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcPi, nSigmaTpcPi, float);                //! TPC nσ for pion hypothesis
DECLARE_SOA_COLUMN(NSigmaItsDe, nSigmaItsDe, float);                //! ITS nσ for deuteron hypothesis
DECLARE_SOA_COLUMN(NSigmaTofDe, nSigmaTofDe, float);                //! TOF nσ for deuteron hypothesis
DECLARE_SOA_COLUMN(NSigmaTofKa, nSigmaTofKa, float);                //! TOF nσ for kaon hypothesis
DECLARE_SOA_COLUMN(NSigmaTofPi, nSigmaTofPi, float);                //! TOF nσ for pion hypothesis
DECLARE_SOA_COLUMN(NItsClusters, nItsClusters, float);              //! Number of ITS clusters used in the track fit
DECLARE_SOA_COLUMN(NItsNClusterSize, nItsNClusterSize, float);      //! Number of ITS clusters size used in the track fit
DECLARE_SOA_COLUMN(NTpcClusters, nTpcClusters, float);              //! Number of TPC clusters used in the track fit
DECLARE_SOA_COLUMN(NTpcSignalsDe, nTpcSignalsDe, float);            //! Number of TPC signas for deuteron
DECLARE_SOA_COLUMN(NTpcSignalsPi, nTpcSignalsPi, float);            //! Number of TPC signas for pion
DECLARE_SOA_COLUMN(NTpcSignalsKa, nTpcSignalsKa, float);            //! Number of TPC signas for kaon
DECLARE_SOA_COLUMN(NItsSignalsDe, nItsSignalsDe, float);            //! Number of ITS signas
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);     //! Candidates falg
DECLARE_SOA_COLUMN(CandidateSign, candidateSign, int8_t);           //! Candidates sign
DECLARE_SOA_COLUMN(FlagMc, flagMc, int8_t);                         //! MC matching flag
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               //! MC origin for reconstructed candidates
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t); //! Resonant MC decay channel for reconstructed candidates
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               //! MC origin for generated particles
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t); //! Resonant MC decay channel for generated candidates
DECLARE_SOA_COLUMN(CtGen, ctGen, float);                            //! Generated ct computed wrt to c-deuteron production vertex, which can be either PV (prompt) or B-hadron decay vertex (non-prompt)
DECLARE_SOA_COLUMN(Cent, cent, float);                              //! Centrality
DECLARE_SOA_COLUMN(VtxZ, vtxZ, float);                              //! Vertex Z
DECLARE_SOA_COLUMN(GIndexCol, gIndexCol, int);                      //! Global index for the collision
DECLARE_SOA_COLUMN(McCollisionId, mcCollisionId, int);              //! Global index for the MC collision
DECLARE_SOA_COLUMN(TimeStamp, timeStamp, int64_t);                  //! Timestamp for the collision
} // namespace full

// Lite table
DECLARE_SOA_TABLE(HfCandCdLite, "AOD", "HFCANDCDLITE",
                  full::MassCd,
                  full::MassLc,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::PtProng0,
                  full::PtProng1,
                  full::PtProng2,
                  full::ImpactParameter0,
                  full::ImpactParameter1,
                  full::ImpactParameter2,
                  full::DecayLength,
                  full::Cpa,
                  full::Chi2PCA,
                  full::NSigmaTpcDe,
                  full::NSigmaTpcPr,
                  full::NSigmaItsDe,
                  full::NSigmaTofDe,
                  full::CandidateSelFlag,
                  full::CandidateSign,
                  full::FlagMc,
                  full::OriginMcRec,
                  full::FlagMcDecayChanRec,
                  full::CtGen,
                  full::Cent);

// full table for local Rotation & Event Mixing
DECLARE_SOA_TABLE(HfCandCdFull, "AOD", "HFCANDCDFULL",
                  full::PxProng0,
                  full::PyProng0,
                  full::PzProng0,
                  full::PxProng1,
                  full::PyProng1,
                  full::PzProng1,
                  full::PxProng2,
                  full::PyProng2,
                  full::PzProng2,
                  full::ImpactParameter0,
                  full::ImpactParameter1,
                  full::ImpactParameter2,
                  full::DecayLength,
                  full::Cpa,
                  full::Chi2PCA,
                  full::NSigmaTpcDe,
                  full::NSigmaTpcPr,
                  full::NSigmaItsDe,
                  full::NSigmaTofDe,
                  full::NSigmaTpcPi,
                  full::NSigmaTofPi,
                  full::NSigmaTpcKa,
                  full::NSigmaTofKa,
                  full::CandidateSelFlag,
                  full::CandidateSign,
                  full::FlagMc,
                  full::OriginMcRec,
                  full::FlagMcDecayChanRec,
                  full::CtGen,
                  full::Cent,
                  full::VtxZ,
                  full::GIndexCol,
                  full::TimeStamp);

DECLARE_SOA_TABLE(HfCandCdGen, "AOD", "HFCANDCDGEN",
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::FlagMc,
                  full::OriginMcGen,
                  full::FlagMcDecayChanGen,
                  full::CtGen,
                  full::Cent,
                  full::VtxZ,
                  full::McCollisionId);
} // namespace o2::aod

struct HfTaskCd {

  Produces<o2::aod::HfCandCdLite> rowCandCdLite;
  Produces<o2::aod::HfCandCdFull> rowCandCdFull;
  Produces<o2::aod::HfCandCdGen> rowCandCdGen;

  Configurable<int> selectionFlagCd{"selectionFlagCd", 1, "Selection Flag for Cd"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_cd_to_de_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<bool> fillTHn{"fillTHn", false, "fill THn"};
  Configurable<bool> fillCandLiteTree{"fillCandLiteTree", false, "Flag to fill candiates lite tree"};
  Configurable<bool> fillCandFullTree{"fillCandFullTree", false, "Flag to fill candiates full tree"};
  Configurable<bool> cfgUseTofPidForDeuteron{"cfgUseTofPidForDeuteron", false, "Use TOF PID for deuteron candidates"};
  Configurable<bool> cfgCutOnDeuteronDcaOrdering{"cfgCutOnDeuteronDcaOrdering", false, "Require deuteron DCA to be smaller than kaon and pion DCAs"};
  Configurable<float> cfgMinDeuteronDcaPreselection{"cfgMinDeuteronDcaPreselection", 0.004, "Minimum deuteron DCA for preselection (cm)"};
  Configurable<float> cfgMaxDeuteronTofPidPreselection{"cfgMaxDeuteronTofPidPreselection", 5, "Maximum |nSigma TOF| for deuteron preselection"};

  SliceCache cache;

  using CollisionsWEvSel = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsMc = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  using CollisionsWithEvSelFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsMcWithEvSelFT0C = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithEvSelFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using CollisionsMcWithEvSelFT0M = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Ms>;

  using CdCandidates = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelCd, aod::HfCand3ProngWPidPiKaDe>>;
  using CdCandidatesMc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelCd, aod::HfCand3ProngWPidPiKaDe, aod::HfCand3ProngMcRec>>;
  using McParticles3ProngMatched = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;
  using HFTracks = soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullDe>;
  using HFTracksMc = soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullDe, aod::McTrackLabels>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_cd::isSelCdToDeKPi >= selectionFlagCd || aod::hf_sel_candidate_cd::isSelCdToPiKDe >= selectionFlagCd;
  Preslice<aod::HfCand3Prong> candCdPerCollision = aod::hf_cand::collisionId;
  PresliceUnsorted<aod::McCollisionLabels> colPerMcCollision = aod::mcparticle::mcCollisionId;

  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {72, 0, 36}, ""};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {400, 2.4, 4.4}, ""};
  ConfigurableAxis thnConfigAxisPtProng{"thnConfigAxisPtProng", {100, 0, 20}, ""};
  ConfigurableAxis thnConfigAxisChi2PCA{"thnConfigAxisChi2PCA", {100, 0, 20}, ""};
  ConfigurableAxis thnConfigAxisDecLength{"thnConfigAxisDecLength", {10, 0, 0.05}, ""};
  ConfigurableAxis thnConfigAxisCPA{"thnConfigAxisCPA", {20, 0.8, 1}, ""};
  ConfigurableAxis thnConfigAxisCentrality{"thnConfigAxisCentrality", {100, 0, 100}, ""};
  ConfigurableAxis thnAxisRapidity{"thnAxisRapidity", {20, -1, 1}, "Cand. rapidity bins"};
  ConfigurableAxis thnConfigAxisGenPtB{"thnConfigAxisGenPtB", {1000, 0, 100}, "Gen Pt B"};
  ConfigurableAxis thnConfigAxisNumPvContr{"thnConfigAxisNumPvContr", {200, -0.5, 199.5}, "Number of PV contributors"};
  ConfigurableAxis thnConfigAxisCt{"thnConfigAxisCt", {500, 0., 5000.}, ""};

  constexpr static std::string_view SignalFolders[] = {"signal", "prompt", "nonprompt"};
  constexpr static std::string_view SignalSuffixes[] = {"", "Prompt", "NonPrompt"};
  const float cmToMum = 1.e4; 

  enum SignalClasses : int {
    Signal = 0,
    Prompt,
    NonPrompt
  };

  HistogramRegistry registry{
    "registry",
    {/// mass candidate
     {"Data/hMass", "3-prong candidates;inv. mass (de K #pi) (GeV/#it{c}^{2})", {HistType::kTH1D, {{400, 2.4, 4.4}}}},
     /// pT
     {"Data/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}}},
     {"Data/hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}}},
     {"Data/hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}}},
     {"Data/hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}}},
     /// DCAxy to prim. vertex prongs
     {"Data/hd0Prong0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1D, {{600, -0.4, 0.4}}}},
     {"Data/hd0Prong1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1D, {{600, -0.4, 0.4}}}},
     {"Data/hd0Prong2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1D, {{600, -0.4, 0.4}}}},
     /// decay length candidate
     {"Data/hDecLength", "3-prong candidates;decay length (cm);entries", {HistType::kTH1D, {{400, 0., 1.}}}},
     /// decay length xy candidate
     {"Data/hDecLengthxy", "3-prong candidates;decay length xy (cm);entries", {HistType::kTH1D, {{400, 0., 1.}}}},
     /// cosine of pointing angle
     {"Data/hCPA", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1D, {{110, -1.1, 1.1}}}},
     /// cosine of pointing angle xy
     {"Data/hCPAxy", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1D, {{110, -1.1, 1.1}}}},
     /// Chi 2 PCA to sec. vertex
     {"Data/hDca2", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1D, {{400, 0., 20.}}}},
     /// eta
     {"Data/hEta", "3-prong candidates;#it{#eta};entries", {HistType::kTH1D, {{100, -2., 2.}}}},
     /// phi
     {"Data/hPhi", "3-prong candidates;#it{#Phi};entries", {HistType::kTH1D, {{100, 0., 6.3}}}}}};

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    std::array<bool, 6> doprocess{doprocessDataStd, doprocessDataStdWithFT0C, doprocessDataStdWithFT0M, doprocessMcStd, doprocessMcStdWithFT0C, doprocessMcStdWithFT0M};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "no or more than one process function enabled! Please check your configuration!");
    }
    const bool isData = doprocessDataStd || doprocessDataStdWithFT0C || doprocessDataStdWithFT0M;

    auto addHistogramsRec = [&](const std::string& histoName, const std::string& xAxisTitle, const std::string& yAxisTitle, const HistogramConfigSpec& configSpec) {
      if (!isData) {
        registry.add(("MC/reconstructed/signal/" + histoName + "RecSig").c_str(), ("3-prong candidates (matched);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
        registry.add(("MC/reconstructed/prompt/" + histoName + "RecSigPrompt").c_str(), ("3-prong candidates (matched, prompt);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
        registry.add(("MC/reconstructed/nonprompt/" + histoName + "RecSigNonPrompt").c_str(), ("3-prong candidates (matched, non-prompt);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
      }
    };

    auto addHistogramsGen = [&](const std::string& histoName, const std::string& xAxisTitle, const std::string& yAxisTitle, const HistogramConfigSpec& configSpec) {
      if (!isData) {
        registry.add(("MC/generated/signal/" + histoName + "Gen").c_str(), ("MC particles (matched);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
        registry.add(("MC/generated/prompt/" + histoName + "GenPrompt").c_str(), ("MC particles (matched, prompt);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
        registry.add(("MC/generated/nonprompt/" + histoName + "GenNonPrompt").c_str(), ("MC particles (matched, non-prompt);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
      }
    };

    addHistogramsRec("hMass", "inv. mass (de K #pi) (GeV/#it{c}^{2})", "", {HistType::kTH1D, {{400, 2.4, 4.4}}});
    addHistogramsRec("hPt", "#it{p}_{T}^{rec.} (GeV/#it{c})", "entries", {HistType::kTH1D, {{360, 0., 36.}}});
    addHistogramsGen("hPt", "#it{p}_{T}^{gen.} (GeV/#it{c})", "entries", {HistType::kTH1D, {{360, 0., 36.}}});
    if (!isData) {
      registry.add("MC/generated/signal/hPtGenSig", "3-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}});
    }
    addHistogramsRec("hPtProng0", "prong 0 #it{p}_{T} (GeV/#it{c})", "entries", {HistType::kTH1D, {{360, 0., 36.}}});
    addHistogramsRec("hPtProng1", "prong 1 #it{p}_{T} (GeV/#it{c})", "entries", {HistType::kTH1D, {{360, 0., 36.}}});
    addHistogramsRec("hPtProng2", "prong 2 #it{p}_{T} (GeV/#it{c})", "entries", {HistType::kTH1D, {{360, 0., 36.}}});
    addHistogramsRec("hd0Prong0", "prong 0 DCAxy to prim. vertex (cm)", "entries", {HistType::kTH1D, {{600, -0.4, 0.4}}});
    addHistogramsRec("hd0Prong1", "prong 1 DCAxy to prim. vertex (cm)", "entries", {HistType::kTH1D, {{600, -0.4, 0.4}}});
    addHistogramsRec("hd0Prong2", "prong 2 DCAxy to prim. vertex (cm)", "entries", {HistType::kTH1D, {{600, -0.4, 0.4}}});
    addHistogramsRec("hDecLength", "decay length (cm)", "entries", {HistType::kTH1D, {{400, 0., 1.}}});
    addHistogramsRec("hDecLengthxy", "decay length xy (cm)", "entries", {HistType::kTH1D, {{400, 0., 1.}}});
    addHistogramsRec("hCPA", "cosine of pointing angle", "entries", {HistType::kTH1D, {{110, -1.1, 1.1}}});
    addHistogramsRec("hCPAxy", "cosine of pointing angle xy", "entries", {HistType::kTH1D, {{110, -1.1, 1.1}}});
    addHistogramsRec("hDca2", "prong Chi2PCA to sec. vertex (cm)", "entries", {HistType::kTH1D, {{400, 0., 20.}}});
    addHistogramsRec("hEta", "#it{#eta}", "entries", {HistType::kTH1D, {{100, -2., 2.}}});
    addHistogramsGen("hEta", "#it{#eta}", "entries", {HistType::kTH1D, {{100, -2., 2.}}});
    addHistogramsGen("hY", "#it{y}", "entries", {HistType::kTH1D, {{100, -2., 2.}}});
    addHistogramsRec("hPhi", "#it{#Phi}", "entries", {HistType::kTH1D, {{100, 0., 6.3}}});
    addHistogramsGen("hPhi", "#it{#Phi}", "entries", {HistType::kTH1D, {{100, 0., 6.3}}});

    /// mass candidate
    if (isData) {
      registry.add("Data/hMassVsPtVsNPvContributors", "3-prong candidates;inv. mass (de K #pi) (GeV/#it{c}^{2}); p_{T}; Number of PV contributors", {HistType::kTH3F, {{400, 2.4, 4.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}, {500, 0., 5000.}}});
    }
    addHistogramsRec("hMassVsPt", "inv. mass (de K #pi) (GeV/#it{c}^{2})", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 2.4, 4.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hMassVsPt", "3-prong candidates;inv. mass (de K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{400, 2.4, 4.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// DCAxy to prim. vertex prongs
    registry.add("Data/hd0VsPtProng0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// decay length candidate
    registry.add("Data/hDecLengthVsPt", "3-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// decay length xy candidate
    registry.add("Data/hDecLengthxyVsPt", "3-prong candidates;decay length xy(cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// cosine of pointing angle
    registry.add("Data/hCPAVsPt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// cosine of pointing angle xy
    registry.add("Data/hCPAxyVsPt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// Chi 2 PCA to sec. vertex
    registry.add("Data/hDca2VsPt", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{400, 0., 20.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// eta
    registry.add("Data/hEtaVsPt", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    addHistogramsGen("hEtaVsPt", "#it{#eta}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    addHistogramsGen("hYVsPt", "#it{y}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// phi
    registry.add("Data/hPhiVsPt", "3-prong candidates;candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    addHistogramsGen("hPhiVsPt", "#it{#Phi}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 6.3}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// selection status
    registry.add("hSelectionStatus", "3-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// impact parameter error
    registry.add("Data/hImpParErrProng0", "3-prong candidates;prong 0 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErrProng1", "3-prong candidates;prong 1 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErrProng2", "3-prong candidates;prong 2 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hNsigmaTPCDeVsP", "deuteron;#it{p} (GeV/#it{c}); n#sigma^{TPC}_{d}", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    registry.add("Data/hNsigmaTPCPrVsP", "proton;#it{p} (GeV/#it{c}); n#sigma^{TPC}_{p}", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    registry.add("Data/hNsigmaTOFDeVsP", "deuteron;#it{p} (GeV/#it{c}); n#sigma^{TOF}_{d}", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    registry.add("Data/hNsigmaITSDeVsP", "deuteron;#it{p} (GeV/#it{c}); n#sigma^{ITS}_{d}", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    registry.add("Data/hTPCSignalDeVsP", "deuteron;#it{p} (GeV/#it{c}); TPC signals", {HistType::kTH2F, {{200, -10.f, 10.f}, {2000, 0, 2000}}});
    registry.add("Data/hTPCSignalPiVsP", "Pion;#it{p} (GeV/#it{c}); TPC signals", {HistType::kTH2F, {{200, -10.f, 10.f}, {2000, 0, 2000}}});
    registry.add("Data/hTPCSignalKaVsP", "Kaon;#it{p} (GeV/#it{c}); TPC signals", {HistType::kTH2F, {{200, -10.f, 10.f}, {2000, 0, 2000}}});
    registry.add("Data/hITSSignalDeVsP", "deuteron;#it{p} (GeV/#it{c}); ITS signals", {HistType::kTH2F, {{200, -10.f, 10.f}, {20, 0, 20}}});
    registry.add("Data/hNsigmaTPCPiVsP", "Pion;#it{p} (GeV/#it{c});n#sigma^{TPC}_{pi};", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    registry.add("Data/hNsigmaTOFPiVsP", "Pion;#it{p} (GeV/#it{c});n#sigma^{TOF}_{pi};", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    registry.add("Data/hNsigmaTPCKaVsP", "Kaon;#it{p} (GeV/#it{c}); n#sigma^{TPC}_{Kaon}", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    registry.add("Data/hNsigmaTOFKaVsP", "Kaon;#it{p} (GeV/#it{c}); n#sigma^{TOF}_{Kaon}", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    addHistogramsRec("hd0VsPtProng0", "prong 0 DCAxy to prim. vertex (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    addHistogramsRec("hd0VsPtProng1", "prong 1 DCAxy to prim. vertex (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    addHistogramsRec("hd0VsPtProng2", "prong 2 DCAxy to prim. vertex (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    addHistogramsRec("hDecLengthVsPt", "decay length (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    addHistogramsRec("hDecLengthxyVsPt", "decay length xy (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    addHistogramsRec("hCPAVsPt", "cosine of pointing angle", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    addHistogramsRec("hCPAxyVsPt", "cosine of pointing angle xy", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    addHistogramsRec("hDca2VsPt", "prong Chi2PCA to sec. vertex (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 20.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    addHistogramsRec("hEtaVsPt", "candidate #it{#eta}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    addHistogramsRec("hPhiVsPt", "candidate #it{#Phi}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 6.3}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    if (fillTHn) {
      const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (de K #pi) (GeV/#it{c}^{2})"};
      const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T}(C_{d}^{+}) (GeV/#it{c})"};
      const AxisSpec thnAxisPtProng0{thnConfigAxisPtProng, "#it{p}_{T}(prong0) (GeV/#it{c})"};
      const AxisSpec thnAxisPtProng1{thnConfigAxisPtProng, "#it{p}_{T}(prong1) (GeV/#it{c})"};
      const AxisSpec thnAxisPtProng2{thnConfigAxisPtProng, "#it{p}_{T}(prong2) (GeV/#it{c})"};
      const AxisSpec thnAxisChi2PCA{thnConfigAxisChi2PCA, "Chi2PCA to sec. vertex (cm)"};
      const AxisSpec thnAxisDecLength{thnConfigAxisDecLength, "decay length (cm)"};
      const AxisSpec thnAxisCPA{thnConfigAxisCPA, "cosine of pointing angle"};
      const AxisSpec thnAxisCentrality{thnConfigAxisCentrality, "centrality (FT0C)"};
      const AxisSpec thnAxisY{thnAxisRapidity, "rapidity"};
      const AxisSpec thnAxisPtB{thnConfigAxisGenPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
      const AxisSpec thnAxisCt{thnConfigAxisCt, "#it{ct} (#mum)"};
      const AxisSpec thnAxisTracklets{thnConfigAxisNumPvContr, "Number of PV contributors"};

      std::vector axesStd{thnAxisMass, thnAxisPt, thnAxisPtProng0, thnAxisPtProng1, thnAxisPtProng2, thnAxisChi2PCA, thnAxisDecLength, thnAxisCPA, thnAxisCentrality};
      std::vector axesGen{thnAxisPt, thnAxisCentrality, thnAxisY, thnAxisTracklets, thnAxisCt, thnAxisPtB};
      registry.add("hnCdVars", isData ? "THn for Reconstructed Cd candidates for data" : "THn for Reconstructed Cd candidates for MC", HistType::kTHnSparseF, axesStd);
      if (!isData) {
        registry.add("hnCdVarsGen", "THn for Generated Cd", HistType::kTHnSparseF, axesGen);
      }
    }
  }

  // taken from: https://github.com/AliceO2Group/O2Physics/blob/master/EventFiltering/PWGCF/CFFilterAll.cxx
  template <typename T>
  float itsSignal(T const& track)
  {
    uint32_t clsizeflag = track.itsClusterSizes();
    auto clSizeLayer0 = (clsizeflag >> (0 * 4)) & 0xf;
    auto clSizeLayer1 = (clsizeflag >> (1 * 4)) & 0xf;
    auto clSizeLayer2 = (clsizeflag >> (2 * 4)) & 0xf;
    auto clSizeLayer3 = (clsizeflag >> (3 * 4)) & 0xf;
    auto clSizeLayer4 = (clsizeflag >> (4 * 4)) & 0xf;
    auto clSizeLayer5 = (clsizeflag >> (5 * 4)) & 0xf;
    auto clSizeLayer6 = (clsizeflag >> (6 * 4)) & 0xf;
    int numLayers = 7;
    int sumClusterSizes = clSizeLayer1 + clSizeLayer2 + clSizeLayer3 + clSizeLayer4 + clSizeLayer5 + clSizeLayer6 + clSizeLayer0;
    float cosLamnda = 1. / std::cosh(track.eta());
    return (static_cast<float>(sumClusterSizes) / numLayers) * cosLamnda;
  };

  /// Helper function for filling MC reconstructed histograms for prompt, nonprompt and common signal
  template <int SignalType, typename CandidateType>
  void fillHistogramsRecSig(CandidateType const& candidate)
  {
    const auto& mcParticleProng0 = candidate.template prong0_as<HFTracksMc>().template mcParticle_as<McParticles3ProngMatched>();
    const auto pdgCodeProng0 = std::abs(mcParticleProng0.pdgCode());
    if ((candidate.isSelCdToDeKPi() >= selectionFlagCd) && pdgCodeProng0 == o2::constants::physics::Pdg::kDeuteron) {
      registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hMassRecSig") + HIST(SignalSuffixes[SignalType]), HfHelper::invMassCdToDeKPi(candidate));
      registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hMassVsPtRecSig") + HIST(SignalSuffixes[SignalType]), HfHelper::invMassCdToDeKPi(candidate), candidate.pt());
    }
    if ((candidate.isSelCdToPiKDe() >= selectionFlagCd) && pdgCodeProng0 == kPiPlus) {
      registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hMassRecSig") + HIST(SignalSuffixes[SignalType]), HfHelper::invMassCdToPiKDe(candidate));
      registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hMassVsPtRecSig") + HIST(SignalSuffixes[SignalType]), HfHelper::invMassCdToPiKDe(candidate), candidate.pt());
    }
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.pt());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hPtProng0RecSig") + HIST(SignalSuffixes[SignalType]), candidate.ptProng0());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hPtProng1RecSig") + HIST(SignalSuffixes[SignalType]), candidate.ptProng1());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hPtProng2RecSig") + HIST(SignalSuffixes[SignalType]), candidate.ptProng2());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hd0Prong0RecSig") + HIST(SignalSuffixes[SignalType]), candidate.impactParameter0());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hd0Prong1RecSig") + HIST(SignalSuffixes[SignalType]), candidate.impactParameter1());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hd0Prong2RecSig") + HIST(SignalSuffixes[SignalType]), candidate.impactParameter2());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hd0VsPtProng0RecSig") + HIST(SignalSuffixes[SignalType]), candidate.impactParameter0(), candidate.pt());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hd0VsPtProng1RecSig") + HIST(SignalSuffixes[SignalType]), candidate.impactParameter1(), candidate.pt());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hd0VsPtProng2RecSig") + HIST(SignalSuffixes[SignalType]), candidate.impactParameter2(), candidate.pt());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hDecLengthRecSig") + HIST(SignalSuffixes[SignalType]), candidate.decayLength());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hDecLengthVsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.decayLength(), candidate.pt());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hDecLengthxyRecSig") + HIST(SignalSuffixes[SignalType]), candidate.decayLengthXY());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hDecLengthxyVsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.decayLengthXY(), candidate.pt());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hCPARecSig") + HIST(SignalSuffixes[SignalType]), candidate.cpa());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hCPAVsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.cpa(), candidate.pt());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hCPAxyRecSig") + HIST(SignalSuffixes[SignalType]), candidate.cpaXY());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hCPAxyVsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.cpaXY(), candidate.pt());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hDca2RecSig") + HIST(SignalSuffixes[SignalType]), candidate.chi2PCA());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hDca2VsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.chi2PCA(), candidate.pt());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hEtaRecSig") + HIST(SignalSuffixes[SignalType]), candidate.eta());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hEtaVsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.eta(), candidate.pt());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hPhiRecSig") + HIST(SignalSuffixes[SignalType]), candidate.phi());
    registry.fill(HIST("MC/reconstructed/") + HIST(SignalFolders[SignalType]) + HIST("/hPhiVsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.phi(), candidate.pt());
  }

  template <typename CollType, typename CandCdMcRec, typename CandCdMcGen, typename TrackWithItsType, typename BcType>
  void fillHistosMcRec(CollType const& collision, CandCdMcRec const& candidates, CandCdMcGen const& mcParticles, TrackWithItsType const& tracksWithItsPid, BcType const& /*bcs*/)
  {
    const auto thisCollId = collision.globalIndex();
    const auto& groupedCdCandidates = candidates.sliceBy(candCdPerCollision, thisCollId);
    const auto bc = collision.template bc_as<BcType>();
    const int64_t timeStamp = bc.timestamp();

    for (const auto& candidate : groupedCdCandidates) {
      if (candidate.flagMcMatchRec() == 0) { // we skip combinatorial background
        continue;
      }
      if (!TESTBIT(candidate.hfflag(), aod::hf_cand_3prong::DecayType::CdToDeKPi)) {
        continue;
      }
      const auto yCd = RecoDecay::y(candidate.pVector(), o2::constants::physics::MassCDeuteron);
      if (yCandRecoMax >= 0. && std::abs(yCd) > yCandRecoMax) {
        continue;
      }

      float ctGen{-1.f}, ptGen{-1.f};
      int pdgCodeProng0{0};
      if (candidate.flagMcMatchRec() == hf_decay::hf_cand_3prong::DecayChannelMain::CDeuteronToDeKPi) {
        const auto& mcParticleProng0 = candidate.template prong0_as<HFTracksMc>().template mcParticle_as<CandCdMcGen>();
        pdgCodeProng0 = std::abs(mcParticleProng0.pdgCode());
        const auto indexMother = RecoDecay::getMother(mcParticles, mcParticleProng0, o2::constants::physics::Pdg::kCDeuteron, true);
        const auto particleMother = mcParticles.rawIteratorAt(indexMother);
        ctGen = RecoDecay::ct(std::array{particleMother.px(), particleMother.py(), particleMother.pz()}, RecoDecay::distance(std::array{particleMother.vx(), particleMother.vy(), particleMother.vz()}, std::array{mcParticleProng0.vx(), mcParticleProng0.vy(), mcParticleProng0.vz()}), o2::constants::physics::MassCDeuteron) * cmToMum;
        ptGen = particleMother.pt();
      }

      if (fillCandLiteTree || fillCandFullTree) {
        float invMassCd = 0.f;
        float invMassLc = 0.f;
        int candFlag = -999;
        int candSign = -999;

        float nSigmaTpcDe = 0.f, nSigmaTpcKa = 0.f, nSigmaTpcPi = 0.f, nSigmaTpcPr = 0.f;
        float nSigmaItsDe = 0.f;
        float nSigmaTofDe = 0.f, nSigmaTofKa = 0.f, nSigmaTofPi = 0.f;

        float dcaDeuteron = 0.f, dcaKaon = 0.f, dcaPion = 0.f;

        const bool selDeKPi = (candidate.isSelCdToDeKPi() >= selectionFlagCd);
        const bool selPiKDe = (candidate.isSelCdToPiKDe() >= selectionFlagCd);

        auto prong1 = candidate.template prong1_as<HFTracksMc>();

        auto prong0Its = tracksWithItsPid.iteratorAt(candidate.prong0Id() - tracksWithItsPid.offset());
        auto prong2Its = tracksWithItsPid.iteratorAt(candidate.prong2Id() - tracksWithItsPid.offset());

        candSign = static_cast<int8_t>(-prong1.sign());
        nSigmaTpcKa = candidate.nSigTpcKa1();
        nSigmaTofKa = candidate.nSigTofKa1();

        if (selDeKPi) {
          invMassCd = HfHelper::invMassCdToDeKPi(candidate);
          invMassLc = HfHelper::invMassLcToPKPi(candidate);
          candFlag = 1;
          nSigmaTpcDe = candidate.nSigTpcDe0();
          nSigmaTpcPr = candidate.nSigTpcPr0();
          nSigmaTofDe = candidate.nSigTofDe0();
          nSigmaTpcPi = candidate.nSigTpcPi2();
          nSigmaTofPi = candidate.nSigTofPi2();
          nSigmaItsDe = prong0Its.itsNSigmaDe();
          dcaDeuteron = candidate.impactParameter0();
          dcaKaon = candidate.impactParameter1();
          dcaPion = candidate.impactParameter2();
        } else if (selPiKDe) {
          invMassCd = HfHelper::invMassCdToPiKDe(candidate);
          invMassLc = HfHelper::invMassLcToPiKP(candidate);
          candFlag = -1;
          nSigmaTpcDe = candidate.nSigTpcDe2();
          nSigmaTpcPr = candidate.nSigTpcPr2();
          nSigmaTofDe = candidate.nSigTofDe2();
          nSigmaTpcPi = candidate.nSigTpcPi0();
          nSigmaTofPi = candidate.nSigTofPi0();
          nSigmaItsDe = prong2Its.itsNSigmaDe();
          dcaDeuteron = candidate.impactParameter2();
          dcaKaon = candidate.impactParameter1();
          dcaPion = candidate.impactParameter0();
        }

        if (cfgUseTofPidForDeuteron && std::abs(nSigmaTofDe) > cfgMaxDeuteronTofPidPreselection) {
          continue;
        }
        if (std::abs(dcaDeuteron) < cfgMinDeuteronDcaPreselection) {
          continue;
        }
        if (cfgCutOnDeuteronDcaOrdering && (std::abs(dcaDeuteron) > std::abs(dcaKaon) || std::abs(dcaDeuteron) > std::abs(dcaPion))) {
          continue;
        }

        if (fillCandLiteTree) {
          rowCandCdLite(
            invMassCd,
            invMassLc,
            candidate.pt(),
            candidate.eta(),
            candidate.phi(),
            candidate.ptProng0(),
            candidate.ptProng1(),
            candidate.ptProng2(),
            candidate.impactParameter0(),
            candidate.impactParameter1(),
            candidate.impactParameter2(),
            candidate.decayLength(),
            candidate.cpa(),
            candidate.chi2PCA(),
            nSigmaTpcDe,
            nSigmaTpcPr,
            nSigmaItsDe,
            nSigmaTofDe,
            candFlag,
            candSign,
            candidate.flagMcMatchRec(),
            candidate.originMcRec(),
            candidate.flagMcDecayChanRec(),
            ctGen,
            o2::hf_centrality::getCentralityColl(collision));
        }

        if (fillCandFullTree) {
          rowCandCdFull(
            candidate.pxProng0(),
            candidate.pyProng0(),
            candidate.pzProng0(),
            candidate.pxProng1(),
            candidate.pyProng1(),
            candidate.pzProng1(),
            candidate.pxProng2(),
            candidate.pyProng2(),
            candidate.pzProng2(),
            candidate.impactParameter0(),
            candidate.impactParameter1(),
            candidate.impactParameter2(),
            candidate.decayLength(),
            candidate.cpa(),
            candidate.chi2PCA(),
            nSigmaTpcDe,
            nSigmaTpcPr,
            nSigmaItsDe,
            nSigmaTofDe,
            nSigmaTpcPi,
            nSigmaTofPi,
            nSigmaTpcKa,
            nSigmaTofKa,
            candFlag,
            candSign,
            candidate.flagMcMatchRec(),
            candidate.originMcRec(),
            candidate.flagMcDecayChanRec(),
            ctGen,
            o2::hf_centrality::getCentralityColl(collision),
            collision.posZ(),
            collision.globalIndex(),
            timeStamp);
        }
      }

      if (std::abs(candidate.flagMcMatchRec()) != hf_decay::hf_cand_3prong::DecayChannelMain::CDeuteronToDeKPi) {
        continue;
      }

      registry.fill(HIST("MC/generated/signal/hPtGenSig"), ptGen);

      fillHistogramsRecSig<Signal>(candidate);
      if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
        fillHistogramsRecSig<Prompt>(candidate);
      } else if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
        fillHistogramsRecSig<NonPrompt>(candidate);
      }

      if (fillTHn) {
        const float cent = o2::hf_centrality::getCentralityColl(collision);
        auto fillTHnRecSig = [&](bool isDeKPi) {
          const auto massCd = isDeKPi ? HfHelper::invMassCdToDeKPi(candidate) : HfHelper::invMassCdToPiKDe(candidate);
          std::vector<double> valuesToFill{massCd, candidate.pt(), candidate.ptProng0(), candidate.ptProng1(), candidate.ptProng2(), candidate.chi2PCA(), candidate.decayLength(), candidate.cpa(), cent};
          registry.get<THnSparse>(HIST("hnCdVars"))->Fill(valuesToFill.data());
        };
        if ((candidate.isSelCdToDeKPi() >= selectionFlagCd) && pdgCodeProng0 == o2::constants::physics::Pdg::kDeuteron) {
          fillTHnRecSig(true);
        }
        if ((candidate.isSelCdToPiKDe() >= selectionFlagCd) && pdgCodeProng0 == kPiPlus) {
          fillTHnRecSig(false);
        }
      }
    }
  }

  template <int SignalType, typename ParticleType>
  void fillHistogramsGen(ParticleType const& particle, float yGen)
  {
    registry.fill(HIST("MC/generated/") + HIST(SignalFolders[SignalType]) + HIST("/hPtGen") + HIST(SignalSuffixes[SignalType]), particle.pt());
    registry.fill(HIST("MC/generated/") + HIST(SignalFolders[SignalType]) + HIST("/hEtaGen") + HIST(SignalSuffixes[SignalType]), particle.eta());
    registry.fill(HIST("MC/generated/") + HIST(SignalFolders[SignalType]) + HIST("/hYGen") + HIST(SignalSuffixes[SignalType]), yGen);
    registry.fill(HIST("MC/generated/") + HIST(SignalFolders[SignalType]) + HIST("/hPhiGen") + HIST(SignalSuffixes[SignalType]), particle.phi());
    registry.fill(HIST("MC/generated/") + HIST(SignalFolders[SignalType]) + HIST("/hEtaVsPtGen") + HIST(SignalSuffixes[SignalType]), particle.eta(), particle.pt());
    registry.fill(HIST("MC/generated/") + HIST(SignalFolders[SignalType]) + HIST("/hYVsPtGen") + HIST(SignalSuffixes[SignalType]), yGen, particle.pt());
    registry.fill(HIST("MC/generated/") + HIST(SignalFolders[SignalType]) + HIST("/hPhiVsPtGen") + HIST(SignalSuffixes[SignalType]), particle.phi(), particle.pt());
  }

  template <typename CandCdMcGen, typename Coll>
  void fillHistosMcGen(CandCdMcGen const& mcParticles, Coll const& recoCollisions)
  {
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) != hf_decay::hf_cand_3prong::DecayChannelMain::CDeuteronToDeKPi) {
        continue;
      }
      const auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassCDeuteron);
      if (yCandGenMax >= 0. && std::abs(yGen) > yCandGenMax) {
        continue;
      }
      const auto& recoCollsPerMcColl = recoCollisions.sliceBy(colPerMcCollision, particle.mcCollision().globalIndex());
      unsigned int numPvContributors = 0;
      float vtxZ = -999.f;
      for (const auto& recCol : recoCollsPerMcColl) {
        numPvContributors = recCol.numContrib() > numPvContributors ? recCol.numContrib() : numPvContributors;
        vtxZ = recCol.posZ();
      }
      const float cent = o2::hf_centrality::getCentralityGenColl(recoCollsPerMcColl);
      const float ptGenB = particle.originMcGen() == RecoDecay::OriginType::Prompt ? -1.f : mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt();
      const auto firstDau = particle.template daughters_as<CandCdMcGen>().begin();
      const float ctGen = RecoDecay::ct(std::array{particle.px(), particle.py(), particle.pz()}, RecoDecay::distance(std::array{particle.vx(), particle.vy(), particle.vz()}, std::array{firstDau.vx(), firstDau.vy(), firstDau.vz()}), o2::constants::physics::MassCDeuteron) * cmToMum;

      fillHistogramsGen<Signal>(particle, yGen);
      if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
        fillHistogramsGen<Prompt>(particle, yGen);
      } else if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
        fillHistogramsGen<NonPrompt>(particle, yGen);
      }

      if (fillTHn) {
        std::vector<double> valuesToFill{particle.pt(), cent, yGen, static_cast<double>(numPvContributors), ctGen, ptGenB};
        registry.get<THnSparse>(HIST("hnCdVarsGen"))->Fill(valuesToFill.data());
      }

      rowCandCdGen(
        particle.pt(),
        particle.eta(),
        particle.phi(),
        yGen,
        particle.flagMcMatchGen(),
        particle.originMcGen(),
        particle.flagMcDecayChanGen(),
        ctGen,
        cent,
        vtxZ,
        particle.mcCollision().globalIndex());
    }
  }

  /// Fill histograms for real data
  template <typename CollType, typename CandType, typename TrackType, typename TrackWithItsType, typename BcType>
  void fillHistosData(CollType const& collision, CandType const& candidates, TrackType const& /*tracks*/, TrackWithItsType const& tracksWithItsPid, BcType const& /*bcs*/)
  {
    auto thisCollId = collision.globalIndex();
    auto groupedCdCandidates = candidates.sliceBy(candCdPerCollision, thisCollId);
    auto numPvContributors = collision.numContrib();
    auto bc = collision.template bc_as<BcType>();
    int64_t timeStamp = bc.timestamp();

    for (const auto& candidate : groupedCdCandidates) {
      if (!TESTBIT(candidate.hfflag(), aod::hf_cand_3prong::DecayType::CdToDeKPi)) {
        continue;
      }

      const auto pt = candidate.pt();
      const auto eta = candidate.eta();
      const auto phi = candidate.phi();
      const auto ptProng0 = candidate.ptProng0();
      const auto ptProng1 = candidate.ptProng1();
      const auto ptProng2 = candidate.ptProng2();
      const auto decayLength = candidate.decayLength();
      const auto decayLengthXY = candidate.decayLengthXY();
      const auto chi2PCA = candidate.chi2PCA();
      const auto cpa = candidate.cpa();
      const auto cpaXY = candidate.cpaXY();
      float invMassCd = 0.f;
      float invMassLc = 0.f;
      if (candidate.isSelCdToDeKPi() >= selectionFlagCd) {
        invMassCd = HfHelper::invMassCdToDeKPi(candidate);
        invMassLc = HfHelper::invMassLcToPKPi(candidate);
      }
      if (candidate.isSelCdToPiKDe() >= selectionFlagCd) {
        invMassCd = HfHelper::invMassCdToPiKDe(candidate);
        invMassLc = HfHelper::invMassLcToPiKP(candidate);
      }

      if (candidate.isSelCdToDeKPi() >= selectionFlagCd) {
        registry.fill(HIST("Data/hMass"), HfHelper::invMassCdToDeKPi(candidate));
        registry.fill(HIST("Data/hMassVsPtVsNPvContributors"), HfHelper::invMassCdToDeKPi(candidate), pt, numPvContributors);
        registry.fill(HIST("Data/hMassVsPt"), HfHelper::invMassCdToDeKPi(candidate), pt);
      }
      if (candidate.isSelCdToPiKDe() >= selectionFlagCd) {
        registry.fill(HIST("Data/hMass"), HfHelper::invMassCdToPiKDe(candidate));
        registry.fill(HIST("Data/hMassVsPtVsNPvContributors"), HfHelper::invMassCdToPiKDe(candidate), pt, numPvContributors);
        registry.fill(HIST("Data/hMassVsPt"), HfHelper::invMassCdToPiKDe(candidate), pt);
      }
      registry.fill(HIST("Data/hPt"), pt);
      registry.fill(HIST("Data/hPtProng0"), ptProng0);
      registry.fill(HIST("Data/hPtProng1"), ptProng1);
      registry.fill(HIST("Data/hPtProng2"), ptProng2);
      registry.fill(HIST("Data/hd0Prong0"), candidate.impactParameter0());
      registry.fill(HIST("Data/hd0Prong1"), candidate.impactParameter1());
      registry.fill(HIST("Data/hd0Prong2"), candidate.impactParameter2());
      registry.fill(HIST("Data/hd0VsPtProng0"), candidate.impactParameter0(), pt);
      registry.fill(HIST("Data/hd0VsPtProng1"), candidate.impactParameter1(), pt);
      registry.fill(HIST("Data/hd0VsPtProng2"), candidate.impactParameter2(), pt);
      registry.fill(HIST("Data/hDecLength"), decayLength);
      registry.fill(HIST("Data/hDecLengthVsPt"), decayLength, pt);
      registry.fill(HIST("Data/hDecLengthxy"), decayLengthXY);
      registry.fill(HIST("Data/hDecLengthxyVsPt"), decayLengthXY, pt);
      registry.fill(HIST("Data/hCPA"), cpa);
      registry.fill(HIST("Data/hCPAVsPt"), cpa, pt);
      registry.fill(HIST("Data/hCPAxy"), cpaXY);
      registry.fill(HIST("Data/hCPAxyVsPt"), cpaXY, pt);
      registry.fill(HIST("Data/hDca2"), chi2PCA);
      registry.fill(HIST("Data/hDca2VsPt"), chi2PCA, pt);
      registry.fill(HIST("Data/hEta"), eta);
      registry.fill(HIST("Data/hEtaVsPt"), eta, pt);
      registry.fill(HIST("Data/hPhi"), phi);
      registry.fill(HIST("Data/hPhiVsPt"), phi, pt);
      registry.fill(HIST("hSelectionStatus"), candidate.isSelCdToDeKPi(), pt);
      registry.fill(HIST("hSelectionStatus"), candidate.isSelCdToPiKDe(), pt);
      registry.fill(HIST("Data/hImpParErrProng0"), candidate.errorImpactParameter0(), pt);
      registry.fill(HIST("Data/hImpParErrProng1"), candidate.errorImpactParameter1(), pt);
      registry.fill(HIST("Data/hImpParErrProng2"), candidate.errorImpactParameter2(), pt);

      float const cent = o2::hf_centrality::getCentralityColl(collision);

      if (fillTHn) {
        double massCd(-1);
        if (candidate.isSelCdToDeKPi() >= selectionFlagCd) {
          massCd = HfHelper::invMassCdToDeKPi(candidate);
          std::vector<double> valuesToFill{massCd, pt, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, cent};
          registry.get<THnSparse>(HIST("hnCdVars"))->Fill(valuesToFill.data());
        }
        if (candidate.isSelCdToPiKDe() >= selectionFlagCd) {
          massCd = HfHelper::invMassCdToPiKDe(candidate);
          std::vector<double> valuesToFill{massCd, pt, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, cent};
          registry.get<THnSparse>(HIST("hnCdVars"))->Fill(valuesToFill.data());
        }
      }

      if (fillCandLiteTree || fillCandFullTree) {

        int candFlag = -999;
        int candSign = -999;

        float nSigmaTpcDe = 0.f, nSigmaTpcKa = 0.f, nSigmaTpcPi = 0.f, nSigmaTpcPr = 0.f;
        float nSigmaItsDe = 0.f;
        float nSigmaTofDe = 0.f, nSigmaTofKa = 0.f, nSigmaTofPi = 0.f;

        float dcaDeuteron = 0.f, dcaKaon = 0.f, dcaPion = 0.f;
        // int itsNClusterSizeDe = 0;

        float tpcSignalsDe = 0.f;
        float tpcSignalsPi = 0.f;
        float tpcSignalsKa = 0.f;

        float itsSignalsDe = 0.f;

        float pSignedDe = -999.f;
        float pSignedPi = -999.f;

        nSigmaTpcKa = candidate.nSigTpcKa1();
        nSigmaTofKa = candidate.nSigTofKa1();

        const bool selDeKPi = (candidate.isSelCdToDeKPi() >= selectionFlagCd);
        const bool selPiKDe = (candidate.isSelCdToPiKDe() >= selectionFlagCd);

        auto prong0 = candidate.template prong0_as<TrackType>();
        auto prong1 = candidate.template prong1_as<TrackType>();
        auto prong2 = candidate.template prong2_as<TrackType>();

        auto prong0Its = tracksWithItsPid.iteratorAt(candidate.prong0Id() - tracksWithItsPid.offset());
        auto prong2Its = tracksWithItsPid.iteratorAt(candidate.prong2Id() - tracksWithItsPid.offset());

        candSign = static_cast<int8_t>(-prong1.sign());

        tpcSignalsKa = prong1.tpcSignal();

        if (selDeKPi) {
          candFlag = 1;
          pSignedDe = prong0.tpcInnerParam() * prong0.sign();
          pSignedPi = prong2.tpcInnerParam() * prong2.sign();
          nSigmaTpcDe = candidate.nSigTpcDe0();
          nSigmaTpcPr = candidate.nSigTpcPr0();
          nSigmaTofDe = candidate.nSigTofDe0();
          nSigmaTpcPi = candidate.nSigTpcPi2();
          nSigmaTofPi = candidate.nSigTofPi2();
          nSigmaItsDe = prong0Its.itsNSigmaDe();
          // itsNClusterSizeDe = prong0.itsClusterSizes();
          tpcSignalsDe = prong0.tpcSignal();
          tpcSignalsPi = prong2.tpcSignal();
          itsSignalsDe = itsSignal(prong0);

          dcaDeuteron = candidate.impactParameter0();
          dcaKaon = candidate.impactParameter1();
          dcaPion = candidate.impactParameter2();
        } else if (selPiKDe) {
          candFlag = -1;
          pSignedDe = prong2.tpcInnerParam() * prong2.sign();
          pSignedPi = prong0.tpcInnerParam() * prong0.sign();
          nSigmaTpcDe = candidate.nSigTpcDe2();
          nSigmaTpcPr = candidate.nSigTpcPr2();
          nSigmaTofDe = candidate.nSigTofDe2();
          nSigmaTpcPi = candidate.nSigTpcPi0();
          nSigmaTofPi = candidate.nSigTofPi0();
          nSigmaItsDe = prong2Its.itsNSigmaDe();
          // itsNClusterSizeDe = prong2.itsClusterSizes();
          tpcSignalsDe = prong2.tpcSignal();
          tpcSignalsPi = prong0.tpcSignal();
          itsSignalsDe = itsSignal(prong2);

          dcaDeuteron = candidate.impactParameter2();
          dcaKaon = candidate.impactParameter1();
          dcaPion = candidate.impactParameter0();
        }

        //  PID QA
        registry.fill(HIST("Data/hNsigmaTPCDeVsP"), pSignedDe, nSigmaTpcDe);
        registry.fill(HIST("Data/hNsigmaTPCPrVsP"), pSignedDe, nSigmaTpcPr);
        registry.fill(HIST("Data/hNsigmaTOFDeVsP"), pSignedDe, nSigmaTofDe);
        registry.fill(HIST("Data/hNsigmaITSDeVsP"), pSignedDe, nSigmaItsDe);
        registry.fill(HIST("Data/hTPCSignalDeVsP"), pSignedDe, tpcSignalsDe);
        registry.fill(HIST("Data/hTPCSignalPiVsP"), pSignedPi, tpcSignalsPi);
        registry.fill(HIST("Data/hTPCSignalKaVsP"), prong1.tpcInnerParam() * prong1.sign(), tpcSignalsKa);
        registry.fill(HIST("Data/hITSSignalDeVsP"), pSignedDe, itsSignalsDe);
        registry.fill(HIST("Data/hNsigmaTPCPiVsP"), pSignedPi, nSigmaTpcPi);
        registry.fill(HIST("Data/hNsigmaTOFPiVsP"), pSignedPi, nSigmaTofPi);
        registry.fill(HIST("Data/hNsigmaTPCKaVsP"), prong1.tpcInnerParam() * prong1.sign(), nSigmaTpcKa);
        registry.fill(HIST("Data/hNsigmaTOFKaVsP"), prong1.tpcInnerParam() * prong1.sign(), nSigmaTofKa);

        if (cfgUseTofPidForDeuteron && std::abs(nSigmaTofDe) > cfgMaxDeuteronTofPidPreselection) {
          continue;
        }
        if (std::abs(dcaDeuteron) < cfgMinDeuteronDcaPreselection) {
          continue;
        }
        if (cfgCutOnDeuteronDcaOrdering && (std::abs(dcaDeuteron) > std::abs(dcaKaon) || std::abs(dcaDeuteron) > std::abs(dcaPion))) {
          continue;
        }
        if (fillCandLiteTree) {

          rowCandCdLite(
            invMassCd,
            invMassLc,
            pt,
            eta,
            phi,
            ptProng0,
            ptProng1,
            ptProng2,
            candidate.impactParameter0(),
            candidate.impactParameter1(),
            candidate.impactParameter2(),
            decayLength,
            cpa,
            chi2PCA,
            nSigmaTpcDe,
            nSigmaTpcPr,
            nSigmaItsDe,
            nSigmaTofDe,
            candFlag,
            candSign,
            0,
            0,
            -1,
            -1.f,
            cent);
        }

        if (fillCandFullTree) {

          rowCandCdFull(
            candidate.pxProng0(),
            candidate.pyProng0(),
            candidate.pzProng0(),
            candidate.pxProng1(),
            candidate.pyProng1(),
            candidate.pzProng1(),
            candidate.pxProng2(),
            candidate.pyProng2(),
            candidate.pzProng2(),
            candidate.impactParameter0(),
            candidate.impactParameter1(),
            candidate.impactParameter2(),
            decayLength,
            cpa,
            chi2PCA,
            nSigmaTpcDe,
            nSigmaTpcPr,
            nSigmaItsDe,
            nSigmaTofDe,
            nSigmaTpcPi,
            nSigmaTofPi,
            nSigmaTpcKa,
            nSigmaTofKa,
            candFlag,
            candSign,
            0,
            0,
            -1,
            -1.f,
            cent,
            collision.posZ(),
            collision.globalIndex(),
            timeStamp);
        }
      }
    }
  }
  /// Run the analysis on real data
  template <typename CollType, typename CandType, typename TrackType, typename TrackWithItsType, typename BcType>
  void runAnalysisPerCollisionData(CollType const& collisions,
                                   CandType const& candidates,
                                   TrackType const& tracks,
                                   TrackWithItsType const& tracksWithItsPid,
                                   BcType const& bcs)
  {

    for (const auto& collision : collisions) {
      fillHistosData(collision, candidates, tracks, tracksWithItsPid, bcs);
    }
  }

  template <typename CollType, typename CandType, typename CandCdMcGen, typename TrackWithItsType, typename BcType>
  void runAnalysisPerCollisionMc(CollType const& collisions,
                                 CandType const& candidates,
                                 CandCdMcGen const& mcParticles,
                                 TrackWithItsType const& tracksWithItsPid,
                                 BcType const& bcs)
  {
    for (const auto& collision : collisions) {
      fillHistosMcRec(collision, candidates, mcParticles, tracksWithItsPid, bcs);
    }
    fillHistosMcGen(mcParticles, collisions);
  }

  void processDataStd(CollisionsWEvSel const& collisions,
                      CdCandidates const& selectedCdCandidates,
                      HFTracks const& tracks,
                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // inlcude ITS PID information
    auto tracksWithItsPid = soa::Attach<HFTracks, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);
    runAnalysisPerCollisionData(collisions, selectedCdCandidates, tracks, tracksWithItsPid, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfTaskCd, processDataStd, "Process Data with the standard method", true);

  void processDataStdWithFT0C(CollisionsWithEvSelFT0C const& collisions,
                              CdCandidates const& selectedCdCandidates,
                              HFTracks const& tracks,
                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // inlcude ITS PID information
    auto tracksWithItsPid = soa::Attach<HFTracks, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);
    runAnalysisPerCollisionData(collisions, selectedCdCandidates, tracks, tracksWithItsPid, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfTaskCd, processDataStdWithFT0C, "Process real data with the standard method and with FT0C centrality", false);

  void processDataStdWithFT0M(CollisionsWithEvSelFT0M const& collisions,
                              CdCandidates const& selectedCdCandidates,
                              HFTracks const& tracks,
                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // inlcude ITS PID information
    auto tracksWithItsPid = soa::Attach<HFTracks, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);
    runAnalysisPerCollisionData(collisions, selectedCdCandidates, tracks, tracksWithItsPid, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfTaskCd, processDataStdWithFT0M, "Process real data with the standard method and with FT0M centrality", false);

  void processMcStd(CollisionsMc const& collisions,
                    CdCandidatesMc const& selectedCdCandidatesMc,
                    McParticles3ProngMatched const& mcParticles,
                    aod::McCollisions const&,
                    HFTracksMc const& tracks,
                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    auto tracksWithItsPid = soa::Attach<HFTracksMc, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);
    runAnalysisPerCollisionMc(collisions, selectedCdCandidatesMc, mcParticles, tracksWithItsPid, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfTaskCd, processMcStd, "Process MC with the standard method", false);

  void processMcStdWithFT0C(CollisionsMcWithEvSelFT0C const& collisions,
                            CdCandidatesMc const& selectedCdCandidatesMc,
                            McParticles3ProngMatched const& mcParticles,
                            aod::McCollisions const&,
                            HFTracksMc const& tracks,
                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    auto tracksWithItsPid = soa::Attach<HFTracksMc, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);
    runAnalysisPerCollisionMc(collisions, selectedCdCandidatesMc, mcParticles, tracksWithItsPid, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfTaskCd, processMcStdWithFT0C, "Process MC with the standard method with FT0C centrality", false);

  void processMcStdWithFT0M(CollisionsMcWithEvSelFT0M const& collisions,
                            CdCandidatesMc const& selectedCdCandidatesMc,
                            McParticles3ProngMatched const& mcParticles,
                            aod::McCollisions const&,
                            HFTracksMc const& tracks,
                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    auto tracksWithItsPid = soa::Attach<HFTracksMc, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);
    runAnalysisPerCollisionMc(collisions, selectedCdCandidatesMc, mcParticles, tracksWithItsPid, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfTaskCd, processMcStdWithFT0M, "Process MC with the standard method with FT0M centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCd>(cfgc)};
}
