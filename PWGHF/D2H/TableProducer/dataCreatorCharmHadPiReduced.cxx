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

/// \file dataCreatorCharmHadPiReduced.cxx
/// \brief Creation of CharmHadron-Pi pairs for Beauty hadron analyses
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Fabio Catalano <fabio.catalano@cern.ch>, CERN
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg University

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannelsLegacy.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/Track.h>

#include <TH1.h>
#include <TPDGCode.h>
#include <TString.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_trkcandsel;

enum Event : uint8_t {
  Processed = 0,
  NoCharmHadPiSelected,
  CharmHadPiSelected,
  kNEvent
};

enum DecayChannel : uint8_t {
  B0ToDminusPi = 0,
  BplusToD0barPi,
  BsToDsminusPi,
  LbToLcplusPi,
  B0ToDstarPi
};

enum WrongCollisionType : uint8_t {
  None = 0,
  WrongAssociation,
  SplitCollision,
};

/// Creation of CharmHad-Pi pairs for Beauty hadrons
struct HfDataCreatorCharmHadPiReduced {
  // Produces AOD tables to store track information
  // collision related tables
  struct : ProducesGroup {
    Produces<aod::HfRedCollisions> hfReducedCollision;
    Produces<aod::HfRedCollCents> hfReducedCollCentrality;
    Produces<aod::HfRedQvectors> hfReducedQvector;
    Produces<aod::HfRedCollExtras> hfReducedCollExtra;
    Produces<aod::HfOrigColCounts> hfCollisionCounter;
    // Pi bachelor related tables
    Produces<aod::HfRedTrackBases> hfTrackPion;
    Produces<aod::HfRedTracksCov> hfTrackCovPion;
    Produces<aod::HfRedTracksPid> hfTrackPidPion;
    Produces<aod::HfRedTracksMom> hfTrackMomPion;
    // charm hadron related tables
    Produces<aod::HfRed2Prongs> hfCand2Prong;
    Produces<aod::HfRed2ProngsCov> hfCand2ProngCov;
    Produces<aod::HfRed2ProngsMl> hfCand2ProngMl;
    Produces<aod::HfRed3Prongs> hfCand3Prong;
    Produces<aod::HfRed3ProngsCov> hfCand3ProngCov;
    Produces<aod::HfRed3ProngsMl> hfCand3ProngMl;
    Produces<aod::HfRedMomDDaugs> hfMomDMesDaugs;
    // D* related tables
    Produces<aod::HfRedSoftPiBases> hfTrackSoftPion;
    Produces<aod::HfRedSoftPiCov> hfTrackCovSoftPion;
    Produces<aod::HfRedSoftPiPid> hfTrackPidSoftPion;
    // PID tables for charm-hadron candidate daughter tracks
    Produces<aod::HfRedPidDau0s> hfCandPidProng0;
    Produces<aod::HfRedPidDau1s> hfCandPidProng1;
    Produces<aod::HfRedPidDau2s> hfCandPidProng2;

    // B-hadron config and MC related tables
    Produces<aod::HfCandB0Configs> rowCandidateConfigB0;
    Produces<aod::HfMcRecRedDpPis> rowHfDPiMcRecReduced;
    Produces<aod::HfMcCheckDpPis> rowHfDPiMcCheckReduced;
    Produces<aod::HfMcRecRedDStarPis> rowHfDStarPiMcRecReduced;
    Produces<aod::HfMcGenRedB0s> rowHfB0McGenReduced;

    Produces<aod::HfCandBpConfigs> rowCandidateConfigBplus;
    Produces<aod::HfMcRecRedD0Pis> rowHfD0PiMcRecReduced;
    Produces<aod::HfMcCheckD0Pis> rowHfD0PiMcCheckReduced;
    Produces<aod::HfMcGenRedBps> rowHfBpMcGenReduced;

    Produces<aod::HfCandBsConfigs> rowCandidateConfigBs;
    Produces<aod::HfMcRecRedDsPis> rowHfDsPiMcRecReduced;
    Produces<aod::HfMcCheckDsPis> rowHfDsPiMcCheckReduced;
    Produces<aod::HfMcGenRedBss> rowHfBsMcGenReduced;

    Produces<aod::HfCandLbConfigs> rowCandidateConfigLb;
    Produces<aod::HfMcRecRedLcPis> rowHfLcPiMcRecReduced;
    Produces<aod::HfMcCheckLcPis> rowHfLcPiMcCheckReduced;
    Produces<aod::HfMcGenRedLbs> rowHfLbMcGenReduced;
  } tables;

  // generic configurables
  struct : o2::framework::ConfigurableGroup {
    // event selection
    Configurable<bool> skipRejectedCollisions{"skipRejectedCollisions", true, "skips collisions rejected by the event selection, instead of flagging only"};
    // magnetic field setting from CCDB
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
    // pair selection
    Configurable<double> invMassWindowCharmHadPi{"invMassWindowCharmHadPi", 0.3, "invariant-mass window for CharmHad-Pi pair preselections (GeV/c2)"};
    // MC extra
    Configurable<bool> checkDecayTypeMc{"checkDecayTypeMc", false, "flag to enable MC checks on decay type"};
  } configs;
  // vertexing
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
    Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
    Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
    Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
    Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any B0 is smaller than this"};
    Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  } vertexConfigurations;
  // selection
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> usePionIsGlobalTrackWoDCA{"usePionIsGlobalTrackWoDCA", true, "check isGlobalTrackWoDCA status for pions, for Run3 studies"};
    Configurable<double> ptPionMin{"ptPionMin", 0.5, "minimum pion pT threshold (GeV/c)"};
    Configurable<double> etaPionMax{"etaPionMax", 0.8, "maximum pion absolute eta threshold"};
    Configurable<std::vector<double>> binsPtPion{"binsPtPion", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for pion DCA XY pT-dependent cut"};
    Configurable<LabeledArray<double>> cutsTrackPionDCA{"cutsTrackPionDCA", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for pions"};
  } trackPionConfigurations;
  // HF flags
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for D+"};
    Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
    Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
    Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
    Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
    Configurable<bool> selectionFlagDstar{"selectionFlagDstar", true, "Selection Flag for D* decay to D0 Pi"};
  } hfflagConfigurations;

  o2::hf_evsel::HfEventSelection hfEvSel;
  o2::hf_evsel::HfEventSelectionMc hfEvSelMc;

  // CCDB service
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber{};

  // O2DatabasePDG service
  Service<o2::framework::O2DatabasePDG> pdg;

  double massC{0.};
  double massB{0.};
  double invMass2ChHadPiMin{0.};
  double invMass2ChHadPiMax{0.};
  double bz{0.};
  static constexpr std::size_t NDaughtersDs{2u};
  static constexpr std::size_t NDaughtersDstar{2u};
  bool isHfCandBhadConfigFilled = false;

  // Fitter to redo D-vertex to get extrapolated daughter tracks (2/3-prong vertex filter)
  o2::vertexing::DCAFitterN<3> df3;
  o2::vertexing::DCAFitterN<2> df2;

  using TracksPid = soa::Join<aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>; // TODO: revert to pion only once the Nsigma variables for the charm-hadron candidate daughters are in the candidate table for 3 prongs too
  using TracksPidWithSel = soa::Join<aod::TracksWCovDcaExtra, TracksPid, aod::TrackSelection>;
  using TracksPidWithSelAndMc = soa::Join<TracksPidWithSel, aod::McTrackLabels>;

  using CandsDplusFiltered = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKa, aod::HfSelDplusToPiKPi>>;
  using CandsDplusFilteredWithMl = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKa, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandsDsFiltered = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKa, aod::HfSelDsToKKPi>>;
  using CandsDsFilteredWithMl = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKa, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandsD0Filtered = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0>>;
  using CandsD0FilteredWithMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0, aod::HfMlD0>>;
  using CandsLcFiltered = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelLc>>;
  using CandsLcFilteredWithMl = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelLc, aod::HfMlLcToPKPi>>;
  using CandsDstarFiltered = soa::Filtered<soa::Join<aod::HfCandDstarsWPid, aod::HfD0FromDstar, aod::HfSelDstarToD0Pi>>;
  using CandsDstarFilteredWithMl = soa::Filtered<soa::Join<aod::HfCandDstarsWPid, aod::HfD0FromDstar, aod::HfSelDstarToD0Pi, aod::HfMlDstarToD0Pi>>;

  using CollisionsWCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs>;
  using CollisionsWCentAndMcLabels = soa::Join<CollisionsWCent, aod::McCollisionLabels>;
  using CollisionsWCentAndQvectors = soa::Join<CollisionsWCent, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorTPCposs, aod::QvectorTPCnegs, aod::QvectorTPCalls>;
  using BCsInfo = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;

  Filter filterSelectDplusCandidates = (aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= hfflagConfigurations.selectionFlagDplus);
  Filter filterSelectDsCandidates = (aod::hf_sel_candidate_ds::isSelDsToKKPi >= hfflagConfigurations.selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= hfflagConfigurations.selectionFlagDs);
  Filter filterSelectDzeroCandidates = (aod::hf_sel_candidate_d0::isSelD0 >= hfflagConfigurations.selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= hfflagConfigurations.selectionFlagD0bar);
  Filter filterSelectLcCandidates = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= hfflagConfigurations.selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= hfflagConfigurations.selectionFlagLc);
  Filter filterSelectDstarCandidates = (aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == hfflagConfigurations.selectionFlagDstar);

  struct : PresliceGroup {
    Preslice<CandsDplusFiltered> candsDplusPerCollision = aod::track_association::collisionId;
    Preslice<CandsDplusFilteredWithMl> candsDplusPerCollisionWithMl = aod::track_association::collisionId;
    Preslice<CandsDsFiltered> candsDsPerCollision = aod::track_association::collisionId;
    Preslice<CandsDsFilteredWithMl> candsDsPerCollisionWithMl = aod::track_association::collisionId;
    Preslice<CandsD0Filtered> candsD0PerCollision = aod::track_association::collisionId;
    Preslice<CandsD0FilteredWithMl> candsD0PerCollisionWithMl = aod::track_association::collisionId;
    Preslice<CandsLcFiltered> candsLcPerCollision = aod::track_association::collisionId;
    Preslice<CandsLcFilteredWithMl> candsLcPerCollisionWithMl = aod::track_association::collisionId;
    Preslice<CandsDstarFiltered> candsDstarPerCollision = aod::track_association::collisionId;
    Preslice<CandsDstarFilteredWithMl> candsDstarPerCollisionWithMl = aod::track_association::collisionId;
    Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
    Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;
    PresliceUnsorted<CollisionsWCentAndMcLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  } preslices;

  std::shared_ptr<TH1> hCandidatesD0, hCandidatesDPlus, hCandidatesDs, hCandidatesLc, hCandidatesD0FromDstar;
  HistogramRegistry registry{"registry"};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  std::array<int, 2> arrPDGResonantDsPhiPi = {kPhi, kPiPlus};      // Ds± → Phi π±
  std::array<int, 2> arrPDGResonantDKstarK = {kK0Star892, kKPlus}; // Ds± → K*(892)0bar K± and D± → K*(892)0bar K±

  void init(InitContext& initContext)
  {
    std::array<int, 28> doProcess = {doprocessDplusPiData, doprocessDplusPiDataWithMl, doprocessDplusPiDataWithQvec, doprocessDplusPiDataWithMlAndQvec, doprocessDplusPiMc, doprocessDplusPiMcWithMl,
                                     doprocessDsPiData, doprocessDsPiDataWithMl, doprocessDsPiDataWithQvec, doprocessDsPiDataWithMlAndQvec, doprocessDsPiMc, doprocessDsPiMcWithMl,
                                     doprocessD0PiData, doprocessD0PiDataWithMl, doprocessD0PiDataWithQvec, doprocessD0PiDataWithMlAndQvec, doprocessD0PiMc, doprocessD0PiMcWithMl,
                                     doprocessLcPiData, doprocessLcPiDataWithMl, doprocessLcPiMc, doprocessLcPiMcWithMl,
                                     doprocessDstarPiData, doprocessDstarPiDataWithMl, doprocessDstarPiDataWithQvec, doprocessDstarPiDataWithMlAndQvec, doprocessDstarPiMc, doprocessDstarPiMcWithMl};
    if (std::accumulate(doProcess.begin(), doProcess.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function can be enabled at a time, please fix your configuration!");
    }

    // invariant-mass window cut
    if (doprocessDplusPiData || doprocessDplusPiDataWithMl || doprocessDplusPiDataWithQvec || doprocessDplusPiDataWithMlAndQvec || doprocessDplusPiMc || doprocessDplusPiMcWithMl) {
      massC = MassDMinus;
      massB = MassB0;
    } else if (doprocessDsPiData || doprocessDsPiDataWithMl || doprocessDsPiDataWithQvec || doprocessDsPiDataWithMlAndQvec || doprocessDsPiMc || doprocessDsPiMcWithMl) {
      massC = MassDS;
      massB = MassBS;
    } else if (doprocessD0PiData || doprocessD0PiDataWithMl || doprocessD0PiDataWithQvec || doprocessD0PiDataWithMlAndQvec || doprocessD0PiMc || doprocessD0PiMcWithMl) {
      massC = MassD0;
      massB = MassBPlus;
    } else if (doprocessLcPiData || doprocessLcPiDataWithMl || doprocessLcPiMc || doprocessLcPiMcWithMl) {
      massC = MassLambdaCPlus;
      massB = MassLambdaB0;
    } else if (doprocessDstarPiData || doprocessDstarPiDataWithMl || doprocessDstarPiDataWithQvec || doprocessDstarPiDataWithMlAndQvec || doprocessDstarPiMc || doprocessDstarPiMcWithMl) {
      massC = MassDStar;
      massB = MassB0;
    }
    invMass2ChHadPiMin = (massB - configs.invMassWindowCharmHadPi) * (massB - configs.invMassWindowCharmHadPi);
    invMass2ChHadPiMax = (massB + configs.invMassWindowCharmHadPi) * (massB + configs.invMassWindowCharmHadPi);

    // Initialize fitter
    if (doprocessDplusPiData || doprocessDplusPiDataWithMl || doprocessDplusPiDataWithQvec || doprocessDplusPiDataWithMlAndQvec || doprocessDplusPiMc || doprocessDplusPiMcWithMl ||
        doprocessDsPiData || doprocessDsPiDataWithMl || doprocessDsPiDataWithQvec || doprocessDsPiDataWithMlAndQvec || doprocessDsPiMc || doprocessDsPiMcWithMl ||
        doprocessLcPiData || doprocessLcPiDataWithMl || doprocessLcPiMc || doprocessLcPiMcWithMl) {
      df3.setPropagateToPCA(vertexConfigurations.propagateToPCA);
      df3.setMaxR(vertexConfigurations.maxR);
      df3.setMaxDZIni(vertexConfigurations.maxDZIni);
      df3.setMinParamChange(vertexConfigurations.minParamChange);
      df3.setMinRelChi2Change(vertexConfigurations.minRelChi2Change);
      df3.setUseAbsDCA(vertexConfigurations.useAbsDCA);
      df3.setWeightedFinalPCA(vertexConfigurations.useWeightedFinalPCA);
      df3.setMatCorrType(noMatCorr);
    } else if (doprocessD0PiData || doprocessD0PiDataWithMl || doprocessD0PiDataWithQvec || doprocessD0PiDataWithMlAndQvec || doprocessD0PiMc || doprocessD0PiMcWithMl ||
               doprocessDstarPiData || doprocessDstarPiDataWithMl || doprocessDstarPiDataWithQvec || doprocessDstarPiDataWithMlAndQvec || doprocessDstarPiMc || doprocessDstarPiMcWithMl) {
      df2.setPropagateToPCA(vertexConfigurations.propagateToPCA);
      df2.setMaxR(vertexConfigurations.maxR);
      df2.setMaxDZIni(vertexConfigurations.maxDZIni);
      df2.setMinParamChange(vertexConfigurations.minParamChange);
      df2.setMinRelChi2Change(vertexConfigurations.minRelChi2Change);
      df2.setUseAbsDCA(vertexConfigurations.useAbsDCA);
      df2.setWeightedFinalPCA(vertexConfigurations.useWeightedFinalPCA);
      df2.setMatCorrType(noMatCorr);
    }

    // Configure CCDB access
    ccdb->setURL(configs.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    runNumber = 0;

    // histograms
    constexpr int NumBinsEvents = kNEvent;
    std::string labels[NumBinsEvents];
    labels[Event::Processed] = "processed";
    labels[Event::NoCharmHadPiSelected] = "without CharmHad-Pi pairs";
    labels[Event::CharmHadPiSelected] = "with CharmHad-Pi pairs";
    static const AxisSpec axisEvents = {NumBinsEvents, 0.5, NumBinsEvents + 0.5, ""};
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {axisEvents});
    for (int iBin = 0; iBin < NumBinsEvents; iBin++) {
      registry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    std::string charmHadTitle;
    std::string histMassTitle;
    if (doprocessDplusPiData || doprocessDplusPiDataWithMl || doprocessDplusPiDataWithQvec || doprocessDplusPiDataWithMlAndQvec || doprocessDplusPiMc || doprocessDplusPiMcWithMl) {
      charmHadTitle = "D^{#plus}";
      histMassTitle = "Dplus";
      registry.add("hMassDplus", "D^{#plus} candidates; #it{M}(K#pi#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
    } else if (doprocessDsPiData || doprocessDsPiDataWithMl || doprocessDsPiDataWithQvec || doprocessDsPiDataWithMlAndQvec || doprocessDsPiMc || doprocessDsPiMcWithMl) {
      charmHadTitle = "D_{s}^{#plus}";
      histMassTitle = "Ds";
      registry.add("hMassDsToKKPi", "D_{s}^{#plus} to KKpi candidates; #it{M}(KK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      registry.add("hMassDsToPiKK", "D_{s}^{#plus} to piKK candidates; #it{M}(KK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
    } else if (doprocessD0PiData || doprocessD0PiDataWithMl || doprocessD0PiDataWithQvec || doprocessD0PiDataWithMlAndQvec || doprocessD0PiMc || doprocessD0PiMcWithMl) {
      charmHadTitle = "D^{0}";
      histMassTitle = "D0";
      registry.add("hMassD0", "D^{0} candidates; #it{M}(K#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      registry.add("hMassD0bar", "#overline{D}^{0} candidates; #it{M}(K#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
    } else if (doprocessLcPiData || doprocessLcPiDataWithMl || doprocessLcPiMc || doprocessLcPiMcWithMl) {
      charmHadTitle = "#Lambda_{c}^{+}";
      histMassTitle = "Lc";
      registry.add("hMassLcToPKPi", "#Lambda_{c}^{+} to KKpi candidates; #it{M}(pK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      registry.add("hMassLcToPiKP", "#Lambda_{c}^{+} to piKK candidates; #it{M}(#piKp) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
    } else if (doprocessDstarPiData || doprocessDstarPiDataWithMl || doprocessDstarPiDataWithQvec || doprocessDstarPiDataWithMlAndQvec || doprocessDstarPiMc || doprocessDstarPiMcWithMl) {
      charmHadTitle = "D^{*}";
      histMassTitle = "Dstar";
      registry.add("hMassDstarToD0Pi", "D^{*} candidates; #it{M}(K#pi#pi) - #it{M}(K#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 1.}}});
    }

    registry.add(Form("hPt%s", histMassTitle.data()), Form("%s candidates candidates;%s candidate #it{p}_{T} (GeV/#it{c});entries", charmHadTitle.data(), charmHadTitle.data()), {HistType::kTH1D, {{100, 0., 10.}}});
    registry.add("hPtPion", "#pi^{#plus} candidates;#pi^{#plus} candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{100, 0., 10.}}});
    registry.add(Form("hCpa%s", histMassTitle.data()), Form("%s candidates;%s cosine of pointing angle;entries", charmHadTitle.data(), charmHadTitle.data()), {HistType::kTH1D, {{110, -1.1, 1.1}}});

    /// candidate monitoring
    hCandidatesD0 = registry.add<TH1>("hCandidatesD0", "D0 candidate counter", {HistType::kTH1D, {axisCands}});
    hCandidatesDPlus = registry.add<TH1>("hCandidatesDPlus", "Dplus candidate counter", {HistType::kTH1D, {axisCands}});
    hCandidatesDs = registry.add<TH1>("hCandidatesDs", "Ds candidate counter", {HistType::kTH1D, {axisCands}});
    hCandidatesLc = registry.add<TH1>("hCandidatesLc", "Lc candidate counter", {HistType::kTH1D, {axisCands}});
    hCandidatesD0FromDstar = registry.add<TH1>("hCandidatesD0FromDstar", "D0 from D* candidate counter", {HistType::kTH1D, {axisCands}});

    setLabelHistoCands(hCandidatesD0);
    setLabelHistoCands(hCandidatesDPlus);
    setLabelHistoCands(hCandidatesDs);
    setLabelHistoCands(hCandidatesLc);
    setLabelHistoCands(hCandidatesD0FromDstar);

    // init HF event selection helper
    hfEvSel.init(registry, &zorroSummary);
    if (doprocessDplusPiMc || doprocessDplusPiMcWithMl ||
        doprocessDsPiMc || doprocessDsPiMcWithMl ||
        doprocessD0PiMc || doprocessD0PiMcWithMl ||
        doprocessLcPiMc || doprocessLcPiMcWithMl ||
        doprocessDstarPiMc || doprocessDstarPiMcWithMl) {
      const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
      for (const DeviceSpec& device : workflows.devices) {
        if (device.name == "hf-data-creator-charm-had-pi-reduced") {
          // init HF event selection helper
          hfEvSelMc.init(device, registry);
          break;
        }
      }
    }
  }

  /// Pion selection (D Pi <-- B0)
  /// \param trackPion is a track with the pion hypothesis
  /// \param trackParCovPion is the track parametrisation of the pion
  /// \param dcaPion is the 2-D array with track DCAs of the pion
  /// \param charmDautracks charm-hadron daughter tracks
  /// \return true if trackPion passes all cuts
  template <typename T1, typename T2, typename T3>
  bool isPionSelected(const T1& trackPion, const T2& trackParCovPion, const T3& dcaPion, const std::vector<T1>& charmDautracks)
  {
    // check isGlobalTrackWoDCA status for pions if wanted
    if (trackPionConfigurations.usePionIsGlobalTrackWoDCA && !trackPion.isGlobalTrackWoDCA()) {
      return false;
    }
    // minimum pT and eta selection
    if (trackParCovPion.getPt() < trackPionConfigurations.ptPionMin || std::abs(trackParCovPion.getEta()) > trackPionConfigurations.etaPionMax || !isSelectedTrackDCA(trackParCovPion, dcaPion, trackPionConfigurations.binsPtPion, trackPionConfigurations.cutsTrackPionDCA)) {
      return false;
    }
    // reject pions that are charm-hadron daughters
    for (const auto& track : charmDautracks) {
      if (trackPion.globalIndex() == track.globalIndex()) {
        return false;
      }
    }

    return true;
  }

  /// Calculates the index of the collision with the maximum number of contributions.
  ///\param collisions are the collisions to search through.
  ///\return The index of the collision with the maximum number of contributions.
  template <typename CColl>
  int64_t getIndexCollisionMaxNumContrib(const CColl& collisions)
  {
    unsigned maxNumContrib = 0;
    int64_t indexCollisionMaxNumContrib = -1;
    for (const auto& collision : collisions) {
      if (collision.numContrib() > maxNumContrib) {
        maxNumContrib = collision.numContrib();
        indexCollisionMaxNumContrib = collision.globalIndex();
      }
    }
    return indexCollisionMaxNumContrib;
  }

  /// Checks if the B meson is associated with a different collision than the one it was generated in
  /// \param particleMother is the mother particle
  /// \param collision is the reconstructed collision
  /// \param indexCollisionMaxNumContrib is the index of the collision associated with a given MC collision with the largest number of contributors.
  /// \param flagWrongCollision is the flag indicating if whether the associated collision is incorrect.
  template <typename PParticle, typename CColl>
  void checkWrongCollision(const PParticle& particleMother,
                           const CColl& collision,
                           const int64_t& indexCollisionMaxNumContrib,
                           int8_t& flagWrongCollision)
  {

    if (particleMother.mcCollision().globalIndex() != collision.mcCollisionId()) {
      flagWrongCollision = WrongCollisionType::WrongAssociation;
    } else {
      if (collision.globalIndex() != indexCollisionMaxNumContrib) {
        flagWrongCollision = WrongCollisionType::SplitCollision;
      }
    }
  }

  /// Function for filling MC reco information in the tables
  /// \param particlesMc is the table with MC particles
  /// \param vecDaughtersB is the vector with all daughter tracks (bachelor pion in last position)
  /// \param indexHfCandCharm is the index of the charm-hadron candidate
  /// \param selectedTracksPion is the map with the indices of selected bachelor pion tracks
  template <uint8_t DecChannel, typename CColl, typename PParticles, typename TTrack>
  void fillMcRecoInfo(const CColl& collision,
                      const PParticles& particlesMc,
                      const std::vector<TTrack>& vecDaughtersB,
                      int& indexHfCandCharm,
                      std::map<int64_t, int64_t> selectedTracksPion,
                      const int64_t indexCollisionMaxNumContrib)
  {

    // we check the MC matching to be stored
    int8_t sign{0};
    int8_t signD{0};
    int8_t flag{0};
    int8_t flagWrongCollision{WrongCollisionType::None};
    int8_t debug{0};
    int pdgCodeBeautyMother{-1};
    int pdgCodeCharmMother{-1};
    int pdgCodeProng0{0};
    int pdgCodeProng1{0};
    int pdgCodeProng2{0};
    int pdgCodeProng3{0};
    float motherPt{-1.f};

    if constexpr (DecChannel == DecayChannel::B0ToDminusPi) {
      // B0 → D- π+ → (π- K+ π-) π+
      auto indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kB0, std::array{-kPiPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 3);
      if (indexRec > -1) {
        // D- → π- K+ π-
        // Printf("Checking D- → π- K+ π-");
        indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, Pdg::kDMinus, std::array{-kPiPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
        if (indexRec > -1) {
          flag = sign * BIT(hf_cand_b0::DecayTypeMc::B0ToDplusPiToPiKPiPi);
        } else {
          debug = 1;
          LOGF(debug, "B0 decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }

        auto indexMother = RecoDecay::getMother(particlesMc, vecDaughtersB.back().template mcParticle_as<PParticles>(), Pdg::kB0, true);
        if (indexMother >= 0) {
          auto particleMother = particlesMc.rawIteratorAt(indexMother);
          motherPt = particleMother.pt();
          checkWrongCollision(particleMother, collision, indexCollisionMaxNumContrib, flagWrongCollision);
        }
      }

      // additional checks for correlated backgrounds
      if (configs.checkDecayTypeMc) {
        // B0 → Ds- π+ → (K- K+ π-) π+
        if (!flag) {
          indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kB0, std::array{-kKPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 3);
          if (indexRec > -1) {
            // Ds- → K- K+ π-
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, -Pdg::kDS, std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
            if (indexRec > -1) {
              flag = sign * BIT(hf_cand_b0::DecayTypeMc::B0ToDsPiToKKPiPi);
            }
          }
        }
        // Bs → Ds- π+ → (K- K+ π-) π+
        if (!flag) {
          indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kBS, std::array{-kKPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 3);
          if (indexRec > -1) {
            // Ds- → K- K+ π-
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, -Pdg::kDS, std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
            if (indexRec > -1) {
              flag = sign * BIT(hf_cand_b0::DecayTypeMc::BsToDsPiToKKPiPi);
            }
          }
        }
        // B0 → D- K+ → (π- K+ π-) K+
        if (!flag) {
          indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kB0, std::array{-kPiPlus, +kKPlus, -kPiPlus, +kKPlus}, true, &sign, 3);
          if (indexRec > -1) {
            // D- → π- K+ π-
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, Pdg::kDMinus, std::array{-kPiPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
            if (indexRec > -1) {
              flag = sign * BIT(hf_cand_b0::DecayTypeMc::B0ToDplusKToPiKPiK);
            }
          }
        }
        // Partly reconstructed decays, i.e. the 4 prongs have a common b-hadron ancestor
        // convention: final state particles are prong0,1,2,3
        if (!flag) {
          auto particleProng0 = vecDaughtersB[0].mcParticle();
          auto particleProng1 = vecDaughtersB[1].mcParticle();
          auto particleProng2 = vecDaughtersB[2].mcParticle();
          auto particleProng3 = vecDaughtersB[3].mcParticle();
          // b-hadron hypothesis
          std::array<int, 3> const bHadronMotherHypos = {Pdg::kB0, Pdg::kBS, Pdg::kLambdaB0};
          // c-hadron hypothesis
          std::array<int, 4> const cHadronMotherHypos = {Pdg::kDPlus, Pdg::kDS, Pdg::kDStar, Pdg::kLambdaCPlus};

          for (const auto& bHadronMotherHypo : bHadronMotherHypos) {
            int const index0Mother = RecoDecay::getMother(particlesMc, particleProng0, bHadronMotherHypo, true);
            int const index1Mother = RecoDecay::getMother(particlesMc, particleProng1, bHadronMotherHypo, true);
            int const index2Mother = RecoDecay::getMother(particlesMc, particleProng2, bHadronMotherHypo, true);
            int const index3Mother = RecoDecay::getMother(particlesMc, particleProng3, bHadronMotherHypo, true);

            // look for common b-hadron ancestor
            if (index0Mother > -1 && index1Mother > -1 && index2Mother > -1 && index3Mother > -1) {
              if (index0Mother == index1Mother && index1Mother == index2Mother && index2Mother == index3Mother) {
                flag = BIT(hf_cand_b0::DecayTypeMc::PartlyRecoDecay);
                pdgCodeBeautyMother = particlesMc.rawIteratorAt(index0Mother).pdgCode();
                pdgCodeCharmMother = 0;
                pdgCodeProng0 = particleProng0.pdgCode();
                pdgCodeProng1 = particleProng1.pdgCode();
                pdgCodeProng2 = particleProng2.pdgCode();
                pdgCodeProng3 = particleProng3.pdgCode();
                // look for common c-hadron mother among prongs 0, 1 and 2
                for (const auto& cHadronMotherHypo : cHadronMotherHypos) {
                  int8_t depthMax = 2;
                  if (cHadronMotherHypo == Pdg::kDStar) { // to include D* -> D π0/γ and D* -> D0 π
                    depthMax += 1;
                  }
                  int const index0CharmMother = RecoDecay::getMother(particlesMc, particleProng0, cHadronMotherHypo, true, &sign, depthMax);
                  int const index1CharmMother = RecoDecay::getMother(particlesMc, particleProng1, cHadronMotherHypo, true, &sign, depthMax);
                  int const index2CharmMother = RecoDecay::getMother(particlesMc, particleProng2, cHadronMotherHypo, true, &sign, depthMax);
                  if (index0CharmMother > -1 && index1CharmMother > -1 && index2CharmMother > -1) {
                    if (index0CharmMother == index1CharmMother && index1CharmMother == index2CharmMother) {
                      // pdgCodeCharmMother =
                      //   Pdg::kDPlus (if D+ is the mother and does not come from D*+)
                      //   Pdg::kDPlus + Pdg::kDStar (if D+ is the mother and D*+ -> D+ π0/γ)
                      //   Pdg::kDStar (if D*+ is the mother and D*+ -> D0 π+)
                      //   Pdg::kDS (if Ds is the mother)
                      //   Pdg::kLambdaCPlus (if Λc+ is the mother)
                      pdgCodeCharmMother += std::abs(particlesMc.rawIteratorAt(index0CharmMother).pdgCode());
                    }
                  }
                }
                break;
              }
            }
          }
        }
        tables.rowHfDPiMcCheckReduced(pdgCodeBeautyMother, pdgCodeCharmMother, pdgCodeProng0, pdgCodeProng1, pdgCodeProng2, pdgCodeProng3);
      }
      tables.rowHfDPiMcRecReduced(indexHfCandCharm, selectedTracksPion[vecDaughtersB.back().globalIndex()], flag, flagWrongCollision, debug, motherPt);
    } else if constexpr (DecChannel == DecayChannel::BsToDsminusPi) {
      // Bs → Ds- π+ → (K- K+ π-) π+
      auto indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kBS, std::array{-kKPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 3);
      if (indexRec > -1) {
        // Ds- → K- K+ π-
        indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, -Pdg::kDS, std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
        if (indexRec > -1) {
          std::vector<int> arrDaughDsIndex;
          std::array<int, 2> arrPDGDaughDs{};
          RecoDecay::getDaughters(particlesMc.rawIteratorAt(indexRec), &arrDaughDsIndex, std::array{0}, 1);
          if (arrDaughDsIndex.size() == NDaughtersDs) {
            for (auto iProng = 0u; iProng < arrDaughDsIndex.size(); ++iProng) {
              auto daughI = particlesMc.rawIteratorAt(arrDaughDsIndex[iProng]);
              arrPDGDaughDs[iProng] = std::abs(daughI.pdgCode());
            }
            // Ds- → Phi π- → K- K+ π- and Ds- → K0* K- → K- K+ π-
            if ((arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[0] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[1]) || (arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[1] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[0])) {
              flag = sign * BIT(hf_cand_bs::DecayTypeMc::BsToDsPiToPhiPiPiToKKPiPi);
            } else if ((arrPDGDaughDs[0] == arrPDGResonantDKstarK[0] && arrPDGDaughDs[1] == arrPDGResonantDKstarK[1]) || (arrPDGDaughDs[0] == arrPDGResonantDKstarK[1] && arrPDGDaughDs[1] == arrPDGResonantDKstarK[0])) {
              flag = sign * BIT(hf_cand_bs::DecayTypeMc::BsToDsPiToK0starKPiToKKPiPi);
            }
          }
        } else {
          debug = 1;
          LOGF(debug, "Bs decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }

        auto indexMother = RecoDecay::getMother(particlesMc, vecDaughtersB.back().template mcParticle_as<PParticles>(), Pdg::kBS, true);
        if (indexMother >= 0) {
          auto particleMother = particlesMc.rawIteratorAt(indexMother);
          motherPt = particleMother.pt();
          checkWrongCollision(particleMother, collision, indexCollisionMaxNumContrib, flagWrongCollision);
        }
      }

      // additional checks for correlated backgrounds
      if (configs.checkDecayTypeMc) {
        // B0 → Ds- π+ → (K- K+ π-) π+
        if (!flag) {
          indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kB0, std::array{-kKPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 3);
          if (indexRec > -1) {
            // Ds- → K- K+ π-
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, -Pdg::kDS, std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
            if (indexRec > -1) {
              std::vector<int> arrDaughDsIndex;
              std::array<int, 2> arrPDGDaughDs{};
              RecoDecay::getDaughters(particlesMc.rawIteratorAt(indexRec), &arrDaughDsIndex, std::array{0}, 1);
              if (arrDaughDsIndex.size() == NDaughtersDs) {
                for (auto iProng = 0u; iProng < arrDaughDsIndex.size(); ++iProng) {
                  auto daughI = particlesMc.rawIteratorAt(arrDaughDsIndex[iProng]);
                  arrPDGDaughDs[iProng] = std::abs(daughI.pdgCode());
                }
                // Ds- → Phi π- → K- K+ π- and Ds- → K0* K- → K- K+ π-
                if ((arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[0] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[1]) || (arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[1] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[0])) {
                  flag = sign * BIT(hf_cand_bs::DecayTypeMc::B0ToDsPiToPhiPiPiToKKPiPi);
                } else if ((arrPDGDaughDs[0] == arrPDGResonantDKstarK[0] && arrPDGDaughDs[1] == arrPDGResonantDKstarK[1]) || (arrPDGDaughDs[0] == arrPDGResonantDKstarK[1] && arrPDGDaughDs[1] == arrPDGResonantDKstarK[0])) {
                  flag = sign * BIT(hf_cand_bs::DecayTypeMc::B0ToDsPiToK0starKPiToKKPiPi);
                }
              }
            }
          }
        }
        // Bs → Ds- K+ → (K- K+ π-) K+
        if (!flag) {
          indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kBS, std::array{-kKPlus, +kKPlus, -kPiPlus, +kKPlus}, true, &sign, 3);
          if (indexRec > -1) {
            // Ds- → K- K+ π-
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, -Pdg::kDS, std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
            if (indexRec > -1) {
              std::vector<int> arrDaughDsIndex;
              std::array<int, 2> arrPDGDaughDs{};
              RecoDecay::getDaughters(particlesMc.rawIteratorAt(indexRec), &arrDaughDsIndex, std::array{0}, 1);
              if (arrDaughDsIndex.size() == NDaughtersDs) {
                for (auto iProng = 0u; iProng < arrDaughDsIndex.size(); ++iProng) {
                  auto daughI = particlesMc.rawIteratorAt(arrDaughDsIndex[iProng]);
                  arrPDGDaughDs[iProng] = std::abs(daughI.pdgCode());
                }
                // Ds- → Phi π- → K- K+ π- and Ds- → K0* K- → K- K+ π-
                if ((arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[0] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[1]) || (arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[1] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[0])) {
                  flag = sign * BIT(hf_cand_bs::DecayTypeMc::BsToDsKToPhiPiKToKKPiK);
                } else if ((arrPDGDaughDs[0] == arrPDGResonantDKstarK[0] && arrPDGDaughDs[1] == arrPDGResonantDKstarK[1]) || (arrPDGDaughDs[0] == arrPDGResonantDKstarK[1] && arrPDGDaughDs[1] == arrPDGResonantDKstarK[0])) {
                  flag = sign * BIT(hf_cand_bs::DecayTypeMc::BsToDsKToK0starKKToKKPiK);
                }
              }
            } else {
              debug = 1;
              LOGF(debug, "Bs decays in the expected final state but the condition on the intermediate state is not fulfilled");
            }

            auto indexMother = RecoDecay::getMother(particlesMc, vecDaughtersB.back().template mcParticle_as<PParticles>(), Pdg::kBS, true);
            if (indexMother >= 0) {
              auto particleMother = particlesMc.rawIteratorAt(indexMother);
              motherPt = particleMother.pt();
              checkWrongCollision(particleMother, collision, indexCollisionMaxNumContrib, flagWrongCollision);
            }
          }
        }
        // Partly reconstructed decays, i.e. the 4 prongs have a common b-hadron ancestor
        // convention: final state particles are prong0,1,2,3
        if (!flag) {
          auto particleProng0 = vecDaughtersB[0].mcParticle();
          auto particleProng1 = vecDaughtersB[1].mcParticle();
          auto particleProng2 = vecDaughtersB[2].mcParticle();
          auto particleProng3 = vecDaughtersB[3].mcParticle();
          // b-hadron hypothesis
          std::array<int, 3> const bHadronMotherHypos = {Pdg::kB0, Pdg::kBS, Pdg::kLambdaB0};
          // c-hadron hypothesis
          std::array<int, 5> const cHadronMotherHypos = {Pdg::kDPlus, Pdg::kDS, Pdg::kDStar, Pdg::kDSStar, Pdg::kLambdaCPlus};

          for (const auto& bHadronMotherHypo : bHadronMotherHypos) {
            int const index0Mother = RecoDecay::getMother(particlesMc, particleProng0, bHadronMotherHypo, true);
            int const index1Mother = RecoDecay::getMother(particlesMc, particleProng1, bHadronMotherHypo, true);
            int const index2Mother = RecoDecay::getMother(particlesMc, particleProng2, bHadronMotherHypo, true);
            int const index3Mother = RecoDecay::getMother(particlesMc, particleProng3, bHadronMotherHypo, true);

            // look for common b-hadron ancestor
            if (index0Mother > -1 && index1Mother > -1 && index2Mother > -1 && index3Mother > -1) {
              if (index0Mother == index1Mother && index1Mother == index2Mother && index2Mother == index3Mother) {
                flag = BIT(hf_cand_bs::DecayTypeMc::PartlyRecoDecay);
                pdgCodeBeautyMother = particlesMc.rawIteratorAt(index0Mother).pdgCode();
                pdgCodeCharmMother = 0;
                pdgCodeProng0 = particleProng0.pdgCode();
                pdgCodeProng1 = particleProng1.pdgCode();
                pdgCodeProng2 = particleProng2.pdgCode();
                pdgCodeProng3 = particleProng3.pdgCode();
                // look for common c-hadron mother among prongs 0, 1 and 2
                for (const auto& cHadronMotherHypo : cHadronMotherHypos) {
                  int8_t depthMax = 2;
                  if (cHadronMotherHypo == Pdg::kDStar || cHadronMotherHypo == Pdg::kDSStar) { // to include D* -> D π0/γ, D* -> D0 π, and Ds* -> Ds π0/γ
                    depthMax += 1;
                  }
                  int const index0CharmMother = RecoDecay::getMother(particlesMc, particleProng0, cHadronMotherHypo, true, &sign, depthMax);
                  int const index1CharmMother = RecoDecay::getMother(particlesMc, particleProng1, cHadronMotherHypo, true, &sign, depthMax);
                  int const index2CharmMother = RecoDecay::getMother(particlesMc, particleProng2, cHadronMotherHypo, true, &sign, depthMax);
                  if (index0CharmMother > -1 && index1CharmMother > -1 && index2CharmMother > -1) {
                    if (index0CharmMother == index1CharmMother && index1CharmMother == index2CharmMother) {
                      // pdgCodeCharmMother =
                      //   Pdg::kDPlus (if D+ is the mother and does not come from D*+)
                      //   Pdg::kDPlus + Pdg::kDStar (if D+ is the mother and D*+ -> D+ π0/γ)
                      //   Pdg::kDStar (if D*+ is the mother and D*+ -> D0 π+)
                      //   Pdg::kDS (if Ds is the mother and does not come from Ds*)
                      //   Pdg::kDS + Pdg::kDSStar (if Ds is the mother and Ds* -> Ds π0/γ)
                      //   Pdg::kLambdaCPlus (if Λc+ is the mother)
                      pdgCodeCharmMother += std::abs(particlesMc.rawIteratorAt(index0CharmMother).pdgCode());
                    }
                  }
                }
                break;
              }
            }
          }
        }
        tables.rowHfDsPiMcCheckReduced(pdgCodeBeautyMother, pdgCodeCharmMother, pdgCodeProng0, pdgCodeProng1, pdgCodeProng2, pdgCodeProng3);
      }
      tables.rowHfDsPiMcRecReduced(indexHfCandCharm, selectedTracksPion[vecDaughtersB.back().globalIndex()], flag, flagWrongCollision, debug, motherPt);
    } else if constexpr (DecChannel == DecayChannel::BplusToD0barPi) {
      // B+ → D0(bar) π+ → (K+ π-) π+
      auto indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, Pdg::kBPlus, std::array{+kPiPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
      if (indexRec > -1) {
        // D0(bar) → K+ π-;
        indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1]}, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign, 1);
        if (indexRec > -1) {
          flag = sign * BIT(hf_cand_bplus::DecayTypeMc::BplusToD0PiToKPiPi);
        } else {
          debug = 1;
          LOGF(debug, "B+ decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }

        auto indexMother = RecoDecay::getMother(particlesMc, vecDaughtersB.back().template mcParticle_as<PParticles>(), Pdg::kBPlus, true);
        if (indexMother >= 0) {
          auto particleMother = particlesMc.rawIteratorAt(indexMother);
          motherPt = particleMother.pt();
          checkWrongCollision(particleMother, collision, indexCollisionMaxNumContrib, flagWrongCollision);
        }
      }
      // additional checks for correlated backgrounds
      if (configs.checkDecayTypeMc) {
        if (!flag) {
          // B+ → D0(bar) K+ → (K+ π-) K+
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, Pdg::kBPlus, std::array{+kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
          if (indexRec > -1) {
            // D0(bar) → K+ π-;
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1]}, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign, 1);
            if (indexRec > -1) {
              flag = sign * BIT(hf_cand_bplus::DecayTypeMc::BplusToD0KToKPiK);
            } else {
              debug = 1;
              LOGF(debug, "B+ decays in the expected final state but the condition on the intermediate state is not fulfilled");
            }

            auto indexMother = RecoDecay::getMother(particlesMc, vecDaughtersB.back().template mcParticle_as<PParticles>(), Pdg::kBPlus, true);
            if (indexMother >= 0) {
              auto particleMother = particlesMc.rawIteratorAt(indexMother);
              motherPt = particleMother.pt();
              checkWrongCollision(particleMother, collision, indexCollisionMaxNumContrib, flagWrongCollision);
            }
          }
        }
        // Partly reconstructed decays, i.e. the 3 prongs have a common b-hadron ancestor
        // convention: final state particles are prong0,1,2
        if (!flag) {
          auto particleProng0 = vecDaughtersB[0].mcParticle();
          auto particleProng1 = vecDaughtersB[1].mcParticle();
          auto particleProng2 = vecDaughtersB[2].mcParticle();
          // b-hadron hypothesis
          std::array<int, 4> const bHadronMotherHypos = {Pdg::kBPlus, Pdg::kB0, Pdg::kBS, Pdg::kLambdaB0};
          // c-hadron hypothesis
          std::array<int, 7> const cHadronMotherHypos = {Pdg::kD0, Pdg::kDPlus, Pdg::kDS, Pdg::kDStar, 423, Pdg::kDSStar, Pdg::kLambdaCPlus};

          for (const auto& bHadronMotherHypo : bHadronMotherHypos) {
            int const index0Mother = RecoDecay::getMother(particlesMc, particleProng0, bHadronMotherHypo, true);
            int const index1Mother = RecoDecay::getMother(particlesMc, particleProng1, bHadronMotherHypo, true);
            int const index2Mother = RecoDecay::getMother(particlesMc, particleProng2, bHadronMotherHypo, true);

            // look for common b-hadron ancestor
            if (index0Mother > -1 && index1Mother > -1 && index2Mother > -1) {
              if (index0Mother == index1Mother && index1Mother == index2Mother) {
                flag = BIT(hf_cand_bplus::DecayTypeMc::PartlyRecoDecay);
                pdgCodeBeautyMother = particlesMc.rawIteratorAt(index0Mother).pdgCode();
                pdgCodeCharmMother = 0;
                pdgCodeProng0 = particleProng0.pdgCode();
                pdgCodeProng1 = particleProng1.pdgCode();
                pdgCodeProng2 = particleProng2.pdgCode();
                // look for common c-hadron mother among prongs 0, 1 and 2
                for (const auto& cHadronMotherHypo : cHadronMotherHypos) {
                  int8_t depthMax = 2;
                  if (cHadronMotherHypo == Pdg::kDStar || cHadronMotherHypo == Pdg::kDStar0 || cHadronMotherHypo == Pdg::kDSStar) { // to include D* -> D π0/γ, D* -> D0 π, and Ds* -> Ds π0/γ
                    depthMax += 1;
                  }
                  int const index0CharmMother = RecoDecay::getMother(particlesMc, particleProng0, cHadronMotherHypo, true, &sign, depthMax);
                  int const index1CharmMother = RecoDecay::getMother(particlesMc, particleProng1, cHadronMotherHypo, true, &sign, depthMax);
                  if (index0CharmMother > -1 && index1CharmMother > -1) {
                    if (index0CharmMother == index1CharmMother) {
                      // pdgCodeCharmMother =
                      //   Pdg::kDPlus (if D+ is the mother and does not come from D*+)
                      //   Pdg::kDPlus + Pdg::kDStar (if D+ is the mother and D*+ -> D+ π0/γ)
                      //   Pdg::kDStar (if D*+ is the mother and D*+ -> D0 π+)
                      //   Pdg::kDS (if Ds is the mother and does not come from Ds*)
                      //   Pdg::kDS + Pdg::kDSStar (if Ds is the mother and Ds* -> Ds π0/γ)
                      //   Pdg::kLambdaCPlus (if Λc+ is the mother)
                      pdgCodeCharmMother += std::abs(particlesMc.rawIteratorAt(index0CharmMother).pdgCode());
                    }
                  }
                }
                break;
              }
            }
          }
        }
        tables.rowHfD0PiMcCheckReduced(pdgCodeBeautyMother, pdgCodeCharmMother, pdgCodeProng0, pdgCodeProng1, pdgCodeProng2);
      }
      tables.rowHfD0PiMcRecReduced(indexHfCandCharm, selectedTracksPion[vecDaughtersB.back().globalIndex()], flag, flagWrongCollision, debug, motherPt);
    } else if constexpr (DecChannel == DecayChannel::LbToLcplusPi) {
      // Lb → Lc+ π- → (p K- π+) π-
      auto indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kLambdaB0, std::array{+kProton, -kKPlus, +kPiPlus, -kPiPlus}, true, &sign, 3);
      if (indexRec > -1) {
        // Lc+ → p K- π+
        // Printf("Checking Lc+ → p K- π+");
        indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec > -1) {
          flag = sign * BIT(hf_cand_lb::DecayTypeMc::LbToLcPiToPKPiPi);
        } else {
          debug = 1;
          LOGF(debug, "Lb decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }

        auto indexMother = RecoDecay::getMother(particlesMc, vecDaughtersB.back().template mcParticle_as<PParticles>(), Pdg::kLambdaB0, true);
        if (indexMother >= 0) {
          auto particleMother = particlesMc.rawIteratorAt(indexMother);
          motherPt = particleMother.pt();
          checkWrongCollision(particleMother, collision, indexCollisionMaxNumContrib, flagWrongCollision);
        }
      }

      // additional checks for correlated backgrounds
      if (configs.checkDecayTypeMc) {
        // Lb → Lc+ K- →  (p K- π+) K-
        if (!flag) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kLambdaB0, std::array{+kProton, -kKPlus, +kPiPlus, -kKPlus}, true, &sign, 3);
          if (indexRec > -1) {
            //  Lc+ → p K- π+
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
            if (indexRec > -1) {
              flag = sign * BIT(hf_cand_lb::DecayTypeMc::LbToLcKToPKPiK);
            }
          }
        }
        // B0 → D- π+ → (π- K+ π-) π+
        if (!flag) {
          indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kB0, std::array{-kPiPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 3);
          if (indexRec > -1) {
            // D- → (π- K+ π-
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, Pdg::kDMinus, std::array{-kPiPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
            if (indexRec > -1) {
              flag = sign * BIT(hf_cand_lb::DecayTypeMc::B0ToDplusPiToPiKPiPi);
            }
          }
        }

        // Partly reconstructed decays, i.e. the 4 prongs have a common b-hadron ancestor
        // convention: final state particles are prong0,1,2,3
        if (!flag) {
          auto particleProng0 = vecDaughtersB[0].mcParticle();
          auto particleProng1 = vecDaughtersB[1].mcParticle();
          auto particleProng2 = vecDaughtersB[2].mcParticle();
          auto particleProng3 = vecDaughtersB[3].mcParticle();
          // b-hadron hypothesis
          std::array<int, 3> const bHadronMotherHypos = {Pdg::kB0, Pdg::kBS, Pdg::kLambdaB0};
          // c-hadron hypothesis
          std::array<int, 4> const cHadronMotherHypos = {Pdg::kDPlus, Pdg::kDS, Pdg::kDStar, Pdg::kLambdaCPlus};

          for (const auto& bHadronMotherHypo : bHadronMotherHypos) {
            int const index0Mother = RecoDecay::getMother(particlesMc, particleProng0, bHadronMotherHypo, true);
            int const index1Mother = RecoDecay::getMother(particlesMc, particleProng1, bHadronMotherHypo, true);
            int const index2Mother = RecoDecay::getMother(particlesMc, particleProng2, bHadronMotherHypo, true);
            int const index3Mother = RecoDecay::getMother(particlesMc, particleProng3, bHadronMotherHypo, true);

            // look for common b-hadron ancestor
            if (index0Mother > -1 && index1Mother > -1 && index2Mother > -1 && index3Mother > -1) {
              if (index0Mother == index1Mother && index1Mother == index2Mother && index2Mother == index3Mother) {
                flag = BIT(hf_cand_lb::DecayTypeMc::PartlyRecoDecay);
                pdgCodeBeautyMother = particlesMc.rawIteratorAt(index0Mother).pdgCode();
                pdgCodeCharmMother = 0;
                pdgCodeProng0 = particleProng0.pdgCode();
                pdgCodeProng1 = particleProng1.pdgCode();
                pdgCodeProng2 = particleProng2.pdgCode();
                pdgCodeProng3 = particleProng3.pdgCode();
                // look for common c-hadron mother among prongs 0, 1 and 2
                for (const auto& cHadronMotherHypo : cHadronMotherHypos) {
                  int8_t depthMax = 2;
                  if (cHadronMotherHypo == Pdg::kDStar) { // to include D* -> D π0/γ and D* -> D0 π
                    depthMax += 1;
                  }
                  int const index0CharmMother = RecoDecay::getMother(particlesMc, particleProng0, cHadronMotherHypo, true, &sign, depthMax);
                  int const index1CharmMother = RecoDecay::getMother(particlesMc, particleProng1, cHadronMotherHypo, true, &sign, depthMax);
                  int const index2CharmMother = RecoDecay::getMother(particlesMc, particleProng2, cHadronMotherHypo, true, &sign, depthMax);
                  if (index0CharmMother > -1 && index1CharmMother > -1 && index2CharmMother > -1) {
                    if (index0CharmMother == index1CharmMother && index1CharmMother == index2CharmMother) {
                      // pdgCodeCharmMother =
                      //   Pdg::kDPlus (if D+ is the mother and does not come from D*+)
                      //   Pdg::kDPlus + Pdg::kDStar (if D+ is the mother and D*+ -> D+ π0/γ)
                      //   Pdg::kDStar (if D*+ is the mother and D*+ -> D0 π+)
                      //   Pdg::kDS (if Ds is the mother)
                      //   Pdg::kLambdaCPlus (if Λc+ is the mother)
                      pdgCodeCharmMother += std::abs(particlesMc.rawIteratorAt(index0CharmMother).pdgCode());
                    }
                  }
                }
                break; // Early exit: found a valid decay chain with common b-hadron mother
              }
            }
          }
        }
        tables.rowHfLcPiMcCheckReduced(pdgCodeBeautyMother, pdgCodeCharmMother, pdgCodeProng0, pdgCodeProng1, pdgCodeProng2, pdgCodeProng3);
      }
      tables.rowHfLcPiMcRecReduced(indexHfCandCharm, selectedTracksPion[vecDaughtersB.back().globalIndex()], flag, flagWrongCollision, debug, motherPt);
    } else if constexpr (DecChannel == DecayChannel::B0ToDstarPi) {
      // B0 → D*+ π- → (D0 π+) π- → (K- π+ π+) π-
      auto indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kB0, std::array{+kKPlus, -kPiPlus, -kPiPlus, +kPiPlus}, true, &sign, 4);
      if (indexRec > -1) {
        // D*+ → (D0 π+) → K- π+ π+
        indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, +Pdg::kDStar, std::array{-kKPlus, +kPiPlus, +kPiPlus}, true, &signD, 3);
        if (indexRec > -1) {
          std::vector<int> arrDaughDstarIndex;
          RecoDecay::getDaughters(particlesMc.rawIteratorAt(indexRec), &arrDaughDstarIndex, std::array{0}, 1);
          if (arrDaughDstarIndex.size() == NDaughtersDstar) {
            bool matchD0{false};
            for (const int iProng : arrDaughDstarIndex) { // o2-linter: disable=const-ref-in-for-loop (it is preferable to copy values in case of simple variables)
              auto daughI = particlesMc.rawIteratorAt(iProng);
              if (std::abs(daughI.pdgCode()) == Pdg::kD0) {
                matchD0 = RecoDecay::isMatchedMCGen(particlesMc, daughI, +Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &signD, 2);
              }
            }
            if (matchD0) {
              flag = sign * BIT(hf_cand_b0::DecayTypeMc::B0ToDstarPiToD0PiPiToKPiPiPi);
            } else {
              debug = 1;
              LOGF(debug, "B0 decays in the expected final state but the condition on D* intermediate state is not fulfilled");
            }
          }
        } else {
          debug = 1;
          LOGF(debug, "B0 decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }

        auto indexMother = RecoDecay::getMother(particlesMc, vecDaughtersB.back().template mcParticle_as<PParticles>(), Pdg::kB0, true);
        if (indexMother >= 0) {
          auto particleMother = particlesMc.rawIteratorAt(indexMother);
          motherPt = particleMother.pt();
          checkWrongCollision(particleMother, collision, indexCollisionMaxNumContrib, flagWrongCollision);
        }
      }
      tables.rowHfDStarPiMcRecReduced(indexHfCandCharm, selectedTracksPion[vecDaughtersB.back().globalIndex()], flag, flagWrongCollision, debug, motherPt);
    }
  }

  template <bool DoMc, bool WithMl, uint8_t DecChannel, bool WithQvec, typename PParticles, typename TTracks, typename CCharmCands, typename Coll, typename BBCs>
  void runDataCreation(Coll const& collision,
                       CCharmCands const& candsC,
                       aod::TrackAssoc const& trackIndices,
                       TTracks const&,
                       PParticles const& particlesMc,
                       uint64_t const& indexCollisionMaxNumContrib,
                       BBCs const&,
                       int& zvtxColl,
                       int& sel8Coll,
                       int& zvtxAndSel8Coll,
                       int& zvtxAndSel8CollAndSoftTrig,
                       int& allSelColl)
  {
    registry.fill(HIST("hEvents"), 1 + Event::Processed);

    const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, BBCs>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
    if (configs.skipRejectedCollisions && hfRejMap != 0) {
      return;
    }

    // helpers for ReducedTables filling
    int const indexHfReducedCollision = tables.hfReducedCollision.lastIndex() + 1;
    // std::map where the key is the track.globalIndex() and
    // the value is the track index in the table of the selected pions
    std::map<int64_t, int64_t> selectedTracksPion;
    bool fillHfReducedCollision = false;

    auto primaryVertex = getPrimaryVertex(collision);

    // Set the magnetic field from ccdb.
    // The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
    // but this is not true when running on Run2 data/MC already converted into AO2Ds.
    auto bc = collision.template bc_as<BBCs>();
    if (runNumber != bc.runNumber()) {
      LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
      auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(configs.ccdbPathGrpMag, bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      runNumber = bc.runNumber();
    }
    df2.setBz(bz);
    df3.setBz(bz);

    auto thisCollId = collision.globalIndex();
    for (const auto& candC : candsC) {
      int indexHfCandCharm{-1};
      float invMassC0{-1.f}, invMassC1{-1.f};
      if constexpr (DecChannel == DecayChannel::B0ToDminusPi) {
        indexHfCandCharm = tables.hfCand3Prong.lastIndex() + 1;
        invMassC0 = HfHelper::invMassDplusToPiKPi(candC);
        registry.fill(HIST("hMassDplus"), invMassC0);
        registry.fill(HIST("hPtDplus"), candC.pt());
        registry.fill(HIST("hCpaDplus"), candC.cpa());
      } else if constexpr (DecChannel == DecayChannel::BsToDsminusPi) {
        indexHfCandCharm = tables.hfCand3Prong.lastIndex() + 1;
        if (candC.isSelDsToKKPi() >= hfflagConfigurations.selectionFlagDs) {
          invMassC0 = HfHelper::invMassDsToKKPi(candC);
          registry.fill(HIST("hMassDsToKKPi"), invMassC0);
        }
        if (candC.isSelDsToPiKK() >= hfflagConfigurations.selectionFlagDs) {
          invMassC1 = HfHelper::invMassDsToPiKK(candC);
          registry.fill(HIST("hMassDsToPiKK"), invMassC1);
        }
        registry.fill(HIST("hPtDs"), candC.pt());
        registry.fill(HIST("hCpaDs"), candC.cpa());
      } else if constexpr (DecChannel == DecayChannel::BplusToD0barPi) {
        indexHfCandCharm = tables.hfCand2Prong.lastIndex() + 1;
        if (candC.isSelD0() >= hfflagConfigurations.selectionFlagD0) {
          invMassC0 = HfHelper::invMassD0ToPiK(candC);
          registry.fill(HIST("hMassD0"), invMassC0);
        }
        if (candC.isSelD0bar() >= hfflagConfigurations.selectionFlagD0bar) {
          invMassC1 = HfHelper::invMassD0barToKPi(candC);
          registry.fill(HIST("hMassD0bar"), invMassC1);
        }
        registry.fill(HIST("hPtD0"), candC.pt());
        registry.fill(HIST("hCpaD0"), candC.cpa());
      } else if constexpr (DecChannel == DecayChannel::LbToLcplusPi) {
        indexHfCandCharm = tables.hfCand3Prong.lastIndex() + 1;
        if (candC.isSelLcToPKPi() >= hfflagConfigurations.selectionFlagLc) {
          invMassC0 = HfHelper::invMassLcToPKPi(candC);
          registry.fill(HIST("hMassLcToPKPi"), invMassC0);
        }
        if (candC.isSelLcToPiKP() >= hfflagConfigurations.selectionFlagLc) {
          invMassC1 = HfHelper::invMassLcToPiKP(candC);
          registry.fill(HIST("hMassLcToPiKP"), invMassC1);
        }
        registry.fill(HIST("hPtLc"), candC.pt());
        registry.fill(HIST("hCpaLc"), candC.cpa());
      } else if constexpr (DecChannel == DecayChannel::B0ToDstarPi) {
        indexHfCandCharm = tables.hfCand2Prong.lastIndex() + 1;
        if (candC.signSoftPi() > 0) {
          invMassC0 = candC.invMassDstar() - candC.invMassD0();
        } else {
          invMassC0 = candC.invMassAntiDstar() - candC.invMassD0Bar();
        }
        registry.fill(HIST("hMassDstarToD0Pi"), invMassC0);
        registry.fill(HIST("hPtDstar"), candC.pt());
        registry.fill(HIST("hCpaDstar"), candC.cpaD0());
      }
      bool fillHfCandCharm = false;

      std::vector<typename TTracks::iterator> charmHadDauTracks{candC.template prong0_as<TTracks>(), candC.template prong1_as<TTracks>()};
      o2::track::TrackParCov trackParCov0 = getTrackParCov(charmHadDauTracks[0]);
      o2::track::TrackParCov trackParCov1 = getTrackParCov(charmHadDauTracks[1]);
      o2::track::TrackParCov trackParCov2{};

      std::array<float, 3> pVec0 = charmHadDauTracks[0].pVector();
      std::array<float, 3> pVec1 = charmHadDauTracks[1].pVector();
      std::array<float, 3> pVec2{};

      auto dca0 = o2::dataformats::DCA(charmHadDauTracks[0].dcaXY(), charmHadDauTracks[0].dcaZ(), charmHadDauTracks[0].cYY(), charmHadDauTracks[0].cZY(), charmHadDauTracks[0].cZZ());
      auto dca1 = o2::dataformats::DCA(charmHadDauTracks[1].dcaXY(), charmHadDauTracks[1].dcaZ(), charmHadDauTracks[1].cYY(), charmHadDauTracks[1].cZY(), charmHadDauTracks[1].cZZ());

      // repropagate tracks to this collision if needed
      if (charmHadDauTracks[0].collisionId() != thisCollId) {
        trackParCov0.propagateToDCA(primaryVertex, bz, &dca0);
      }

      if (charmHadDauTracks[1].collisionId() != thisCollId) {
        trackParCov1.propagateToDCA(primaryVertex, bz, &dca1);
      }

      // third track, if it's a 3-prong
      if constexpr (DecChannel == DecayChannel::B0ToDminusPi || DecChannel == DecayChannel::BsToDsminusPi || DecChannel == DecayChannel::LbToLcplusPi || DecChannel == DecayChannel::B0ToDstarPi) {
        if constexpr (DecChannel == DecayChannel::B0ToDstarPi) {
          charmHadDauTracks.push_back(candC.template prongPi_as<TTracks>()); // Soft pion from D* decay
        } else {
          charmHadDauTracks.push_back(candC.template prong2_as<TTracks>());
        }
        trackParCov2 = getTrackParCov(charmHadDauTracks[2]);
        pVec2 = charmHadDauTracks[2].pVector();
        auto dca2 = o2::dataformats::DCA(charmHadDauTracks[2].dcaXY(), charmHadDauTracks[2].dcaZ(), charmHadDauTracks[2].cYY(), charmHadDauTracks[2].cZY(), charmHadDauTracks[2].cZZ());

        if (charmHadDauTracks[2].collisionId() != thisCollId) {
          trackParCov2.propagateToDCA(primaryVertex, bz, &dca2);
        }
      }

      // ---------------------------------
      // reconstruct charm candidate secondary vertex
      o2::track::TrackParCov trackParCovCharmHad{};
      std::array<float, 3> pVecCharm{};
      if constexpr (DecChannel == DecayChannel::B0ToDminusPi || DecChannel == DecayChannel::BsToDsminusPi || DecChannel == DecayChannel::LbToLcplusPi) { // D∓ → π∓ K± π∓ and Ds∓ → K∓ K± π∓ and Lc∓ → p∓ K± π∓

        if constexpr (DecChannel == DecayChannel::B0ToDminusPi) {
          hCandidatesDPlus->Fill(SVFitting::BeforeFit);
        } else if constexpr (DecChannel == DecayChannel::BsToDsminusPi) {
          hCandidatesDs->Fill(SVFitting::BeforeFit);
        } else if constexpr (DecChannel == DecayChannel::LbToLcplusPi) {
          hCandidatesLc->Fill(SVFitting::BeforeFit);
        }

        try {
          if (df3.process(trackParCov0, trackParCov1, trackParCov2) == 0) {
            continue;
          }
        } catch (const std::runtime_error& error) {
          LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
          if constexpr (DecChannel == DecayChannel::B0ToDminusPi) {
            hCandidatesDPlus->Fill(SVFitting::Fail);
          } else if constexpr (DecChannel == DecayChannel::BsToDsminusPi) {
            hCandidatesDs->Fill(SVFitting::Fail);
          } else if constexpr (DecChannel == DecayChannel::LbToLcplusPi) {
            hCandidatesLc->Fill(SVFitting::Fail);
          }
          continue;
        }
        if constexpr (DecChannel == DecayChannel::B0ToDminusPi) {
          hCandidatesDPlus->Fill(SVFitting::FitOk);
        } else if constexpr (DecChannel == DecayChannel::BsToDsminusPi) {
          hCandidatesDs->Fill(SVFitting::FitOk);
        } else if constexpr (DecChannel == DecayChannel::LbToLcplusPi) {
          hCandidatesLc->Fill(SVFitting::FitOk);
        }

        auto secondaryVertexCharm = df3.getPCACandidate();
        trackParCov0.propagateTo(secondaryVertexCharm[0], bz);
        trackParCov1.propagateTo(secondaryVertexCharm[0], bz);
        trackParCov2.propagateTo(secondaryVertexCharm[0], bz);
        df3.getTrack(0).getPxPyPzGlo(pVec0);
        df3.getTrack(1).getPxPyPzGlo(pVec1);
        df3.getTrack(2).getPxPyPzGlo(pVec2);
        pVecCharm = RecoDecay::pVec(pVec0, pVec1, pVec2);
        trackParCovCharmHad = df3.createParentTrackParCov();
        trackParCovCharmHad.setAbsCharge(charmHadDauTracks[1].sign());   // to be sure
      } else if constexpr (DecChannel == DecayChannel::BplusToD0barPi) { // D0(bar) → K± π∓

        hCandidatesD0->Fill(SVFitting::BeforeFit);
        try {
          if (df2.process(trackParCov0, trackParCov1) == 0) {
            continue;
          }
        } catch (const std::runtime_error& error) {
          LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
          hCandidatesD0->Fill(SVFitting::Fail);
          continue;
        }
        hCandidatesD0->Fill(SVFitting::FitOk);

        auto secondaryVertexCharm = df2.getPCACandidate();
        trackParCov0.propagateTo(secondaryVertexCharm[0], bz);
        trackParCov1.propagateTo(secondaryVertexCharm[0], bz);
        df2.getTrack(0).getPxPyPzGlo(pVec0);
        df2.getTrack(1).getPxPyPzGlo(pVec1);
        pVecCharm = RecoDecay::pVec(pVec0, pVec1);
        trackParCovCharmHad = df2.createParentTrackParCov();
        trackParCovCharmHad.setAbsCharge(0); // to be sure
      } else if constexpr (DecChannel == DecayChannel::B0ToDstarPi) {

        hCandidatesD0FromDstar->Fill(SVFitting::BeforeFit);
        try {
          // D0 vertex
          if (df2.process(trackParCov0, trackParCov1) == 0) {
            continue;
          }
        } catch (const std::runtime_error& error) {
          LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
          hCandidatesD0FromDstar->Fill(SVFitting::Fail);
          continue;
        }
        hCandidatesD0FromDstar->Fill(SVFitting::FitOk);
        auto secondaryVertexCharm = df2.getPCACandidate();
        trackParCov0.propagateTo(secondaryVertexCharm[0], bz);
        trackParCov1.propagateTo(secondaryVertexCharm[0], bz);
        df2.getTrack(0).getPxPyPzGlo(pVec0);
        df2.getTrack(1).getPxPyPzGlo(pVec1);
        pVecCharm = RecoDecay::pVec(pVec0, pVec1, pVec2);
        trackParCovCharmHad = df2.createParentTrackParCov();
        trackParCovCharmHad.setAbsCharge(0); // to be sure
      }

      float ptDauMin = 1.e6, etaDauMin = 999.f, chi2TpcDauMax = -1.f;
      int nItsClsDauMin = 8, nTpcCrossRowsDauMin = 200;
      for (const auto& charmHadTrack : charmHadDauTracks) {
        if (charmHadTrack.pt() < ptDauMin) {
          ptDauMin = charmHadTrack.pt();
        }
        if (std::abs(charmHadTrack.eta()) < etaDauMin) {
          etaDauMin = std::abs(charmHadTrack.eta());
        }
        if (charmHadTrack.itsNCls() < nItsClsDauMin) {
          nItsClsDauMin = charmHadTrack.itsNCls();
        }
        if (charmHadTrack.tpcNClsCrossedRows() < nTpcCrossRowsDauMin) {
          nTpcCrossRowsDauMin = charmHadTrack.tpcNClsCrossedRows();
        }
        if (charmHadTrack.tpcChi2NCl() > chi2TpcDauMax) {
          chi2TpcDauMax = charmHadTrack.tpcChi2NCl();
        }
      }

      for (const auto& trackId : trackIndices) {
        auto trackPion = trackId.template track_as<TTracks>();

        // apply selections on pion tracks
        auto trackParCovPion = getTrackParCov(trackPion);
        std::array<float, 2> dcaPion{trackPion.dcaXY(), trackPion.dcaZ()};
        std::array<float, 3> pVecPion = trackPion.pVector();
        if (trackPion.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovPion, 2.f, noMatCorr, &dcaPion);
          getPxPyPz(trackParCovPion, pVecPion);
        }

        // reject pi D with same sign as D
        if constexpr (DecChannel == DecayChannel::B0ToDminusPi || DecChannel == DecayChannel::BsToDsminusPi || DecChannel == DecayChannel::LbToLcplusPi) { // D∓ → π∓ K± π∓ and Ds∓ → K∓ K± π∓ and Lc∓ → p∓ K± π∓
          if (trackPion.sign() * charmHadDauTracks[0].sign() > 0) {
            continue;
          }
        } else if constexpr (DecChannel == DecayChannel::BplusToD0barPi) { // D0(bar) → K± π∓
          if (!((candC.isSelD0() >= hfflagConfigurations.selectionFlagD0 && trackPion.sign() < 0) || (candC.isSelD0bar() >= hfflagConfigurations.selectionFlagD0bar && trackPion.sign() > 0))) {
            continue;
          }
        } else if constexpr (DecChannel == DecayChannel::B0ToDstarPi) { // D*+ → D0 π+
          if (trackPion.sign() * charmHadDauTracks.back().sign() > 0) {
            continue;
          }
        }

        // apply selections on pion tracks
        if (!isPionSelected(trackPion, trackParCovPion, dcaPion, charmHadDauTracks)) {
          continue;
        }

        registry.fill(HIST("hPtPion"), trackParCovPion.getPt());
        // compute invariant mass square and apply selection
        auto invMass2DPi = RecoDecay::m2(std::array{pVecCharm, pVecPion}, std::array{massC, MassPiPlus});
        if ((invMass2DPi < invMass2ChHadPiMin) || (invMass2DPi > invMass2ChHadPiMax)) {
          continue;
        }

        // fill Pion tracks table
        // if information on track already stored, go to next track
        if (!selectedTracksPion.count(trackPion.globalIndex())) {
          tables.hfTrackPion(trackPion.globalIndex(), indexHfReducedCollision,
                             trackParCovPion.getX(), trackParCovPion.getAlpha(),
                             trackParCovPion.getY(), trackParCovPion.getZ(), trackParCovPion.getSnp(),
                             trackParCovPion.getTgl(), trackParCovPion.getQ2Pt(),
                             trackPion.itsNCls(), trackPion.tpcNClsCrossedRows(), trackPion.tpcChi2NCl());
          tables.hfTrackCovPion(trackParCovPion.getSigmaY2(), trackParCovPion.getSigmaZY(), trackParCovPion.getSigmaZ2(),
                                trackParCovPion.getSigmaSnpY(), trackParCovPion.getSigmaSnpZ(),
                                trackParCovPion.getSigmaSnp2(), trackParCovPion.getSigmaTglY(), trackParCovPion.getSigmaTglZ(),
                                trackParCovPion.getSigmaTglSnp(), trackParCovPion.getSigmaTgl2(),
                                trackParCovPion.getSigma1PtY(), trackParCovPion.getSigma1PtZ(), trackParCovPion.getSigma1PtSnp(),
                                trackParCovPion.getSigma1PtTgl(), trackParCovPion.getSigma1Pt2());
          tables.hfTrackPidPion(trackPion.hasTPC(), trackPion.hasTOF(),
                                trackPion.tpcNSigmaPi(), trackPion.tofNSigmaPi());
          tables.hfTrackMomPion(pVecPion[0], pVecPion[1], pVecPion[2], trackPion.sign());
          // add trackPion.globalIndex() to a list
          // to keep memory of the pions filled in the table and avoid refilling them if they are paired to another D candidate
          // and keep track of their index in tables.hfTrackPion for McRec purposes
          selectedTracksPion[trackPion.globalIndex()] = tables.hfTrackPion.lastIndex();
        }

        if constexpr (DoMc) {
          std::vector<typename TTracks::iterator> beautyHadDauTracks{};
          beautyHadDauTracks.reserve(charmHadDauTracks.size());
          for (const auto& track : charmHadDauTracks) {
            beautyHadDauTracks.push_back(track);
          }
          beautyHadDauTracks.push_back(trackPion);
          fillMcRecoInfo<DecChannel>(collision, particlesMc, beautyHadDauTracks, indexHfCandCharm, selectedTracksPion, indexCollisionMaxNumContrib);
        }
        fillHfCandCharm = true;
      } // pion loop
      if (fillHfCandCharm) { // fill candCplus table only once per D candidate
        constexpr std::size_t NSizeMLScore{3u};
        if constexpr (DecChannel == DecayChannel::B0ToDminusPi || DecChannel == DecayChannel::BsToDsminusPi || DecChannel == DecayChannel::LbToLcplusPi) { // D∓ → π∓ K± π∓ and Ds∓ → K∓ K± π∓ and Lc∓ → p∓ K± π∓
          tables.hfCand3Prong(charmHadDauTracks[0].globalIndex(), charmHadDauTracks[1].globalIndex(), charmHadDauTracks[2].globalIndex(),
                              indexHfReducedCollision,
                              trackParCovCharmHad.getX(), trackParCovCharmHad.getAlpha(),
                              trackParCovCharmHad.getY(), trackParCovCharmHad.getZ(), trackParCovCharmHad.getSnp(),
                              trackParCovCharmHad.getTgl(), trackParCovCharmHad.getQ2Pt(),
                              candC.xSecondaryVertex(), candC.ySecondaryVertex(), candC.zSecondaryVertex(), invMassC0, invMassC1,
                              ptDauMin, etaDauMin, nItsClsDauMin, nTpcCrossRowsDauMin, chi2TpcDauMax);
          tables.hfCand3ProngCov(trackParCovCharmHad.getSigmaY2(), trackParCovCharmHad.getSigmaZY(), trackParCovCharmHad.getSigmaZ2(),
                                 trackParCovCharmHad.getSigmaSnpY(), trackParCovCharmHad.getSigmaSnpZ(),
                                 trackParCovCharmHad.getSigmaSnp2(), trackParCovCharmHad.getSigmaTglY(), trackParCovCharmHad.getSigmaTglZ(),
                                 trackParCovCharmHad.getSigmaTglSnp(), trackParCovCharmHad.getSigmaTgl2(),
                                 trackParCovCharmHad.getSigma1PtY(), trackParCovCharmHad.getSigma1PtZ(), trackParCovCharmHad.getSigma1PtSnp(),
                                 trackParCovCharmHad.getSigma1PtTgl(), trackParCovCharmHad.getSigma1Pt2());
          float nSigmaTpcPr0{-999.f}, nSigmaTpcPr1{-999.f}, nSigmaTpcPr2{-999.f};
          float nSigmaTofPr0{-999.f}, nSigmaTofPr1{-999.f}, nSigmaTofPr2{-999.f};
          if constexpr (DecChannel == DecayChannel::LbToLcplusPi) {
            /// assign non-dummy values only for Lb->LcPi analysis
            nSigmaTpcPr0 = candC.nSigTpcPr0();
            nSigmaTpcPr1 = candC.nSigTpcPr1();
            nSigmaTpcPr2 = candC.nSigTpcPr2();
            nSigmaTofPr0 = candC.nSigTofPr0();
            nSigmaTofPr1 = candC.nSigTofPr1();
            nSigmaTofPr2 = candC.nSigTofPr2();
          }
          tables.hfCandPidProng0(candC.nSigTpcPi0(), candC.nSigTofPi0(), candC.nSigTpcKa0(), candC.nSigTofKa0(), nSigmaTpcPr0, nSigmaTofPr0, charmHadDauTracks[0].hasTOF(), charmHadDauTracks[0].hasTPC());
          tables.hfCandPidProng1(candC.nSigTpcPi1(), candC.nSigTofPi1(), candC.nSigTpcKa1(), candC.nSigTofKa1(), nSigmaTpcPr1, nSigmaTofPr1, charmHadDauTracks[1].hasTOF(), charmHadDauTracks[1].hasTPC());
          tables.hfCandPidProng2(candC.nSigTpcPi2(), candC.nSigTofPi2(), candC.nSigTpcKa2(), candC.nSigTofKa2(), nSigmaTpcPr2, nSigmaTofPr2, charmHadDauTracks[2].hasTOF(), charmHadDauTracks[2].hasTPC());
          if constexpr (WithMl) {
            std::array<float, 6> mlScores = {-1.f, -1.f, -1.f, -1.f, -1.f, -1.f};
            if constexpr (DecChannel == DecayChannel::B0ToDminusPi) {
              tables.hfCand3ProngMl(candC.mlProbDplusToPiKPi()[0], candC.mlProbDplusToPiKPi()[1], candC.mlProbDplusToPiKPi()[2], -1., -1., -1.);
            } else if constexpr (DecChannel == DecayChannel::BsToDsminusPi) {
              if (candC.mlProbDsToKKPi().size() == NSizeMLScore) {
                std::copy(candC.mlProbDsToKKPi().begin(), candC.mlProbDsToKKPi().end(), mlScores.begin());
              }
              if (candC.mlProbDsToPiKK().size() == NSizeMLScore) {
                std::copy(candC.mlProbDsToPiKK().begin(), candC.mlProbDsToPiKK().end(), mlScores.begin() + 3);
              }
              tables.hfCand3ProngMl(mlScores[0], mlScores[1], mlScores[2], mlScores[3], mlScores[4], mlScores[5]);
            } else if constexpr (DecChannel == DecayChannel::LbToLcplusPi) {
              if (candC.mlProbLcToPKPi().size() == NSizeMLScore) {
                std::copy(candC.mlProbLcToPKPi().begin(), candC.mlProbLcToPKPi().end(), mlScores.begin());
              }
              if (candC.mlProbLcToPiKP().size() == NSizeMLScore) {
                std::copy(candC.mlProbLcToPiKP().begin(), candC.mlProbLcToPiKP().end(), mlScores.begin() + 3);
              }
              tables.hfCand3ProngMl(mlScores[0], mlScores[1], mlScores[2], mlScores[3], mlScores[4], mlScores[5]);
            }
          }
        } else if constexpr (DecChannel == DecayChannel::BplusToD0barPi) { // D0(bar) → K± π∓
          tables.hfCand2Prong(charmHadDauTracks[0].globalIndex(), charmHadDauTracks[1].globalIndex(),
                              indexHfReducedCollision,
                              trackParCovCharmHad.getX(), trackParCovCharmHad.getAlpha(),
                              trackParCovCharmHad.getY(), trackParCovCharmHad.getZ(), trackParCovCharmHad.getSnp(),
                              trackParCovCharmHad.getTgl(), trackParCovCharmHad.getQ2Pt(),
                              candC.xSecondaryVertex(), candC.ySecondaryVertex(), candC.zSecondaryVertex(), invMassC0, invMassC1,
                              ptDauMin, etaDauMin, nItsClsDauMin, nTpcCrossRowsDauMin, chi2TpcDauMax);
          tables.hfCand2ProngCov(trackParCovCharmHad.getSigmaY2(), trackParCovCharmHad.getSigmaZY(), trackParCovCharmHad.getSigmaZ2(),
                                 trackParCovCharmHad.getSigmaSnpY(), trackParCovCharmHad.getSigmaSnpZ(),
                                 trackParCovCharmHad.getSigmaSnp2(), trackParCovCharmHad.getSigmaTglY(), trackParCovCharmHad.getSigmaTglZ(),
                                 trackParCovCharmHad.getSigmaTglSnp(), trackParCovCharmHad.getSigmaTgl2(),
                                 trackParCovCharmHad.getSigma1PtY(), trackParCovCharmHad.getSigma1PtZ(), trackParCovCharmHad.getSigma1PtSnp(),
                                 trackParCovCharmHad.getSigma1PtTgl(), trackParCovCharmHad.getSigma1Pt2());
          tables.hfCandPidProng0(candC.nSigTpcPi0(), candC.nSigTofPi0(), candC.nSigTpcKa0(), candC.nSigTofKa0(), 0., 0., charmHadDauTracks[0].hasTOF(), charmHadDauTracks[0].hasTPC());
          tables.hfCandPidProng1(candC.nSigTpcPi1(), candC.nSigTofPi1(), candC.nSigTpcKa1(), candC.nSigTofKa1(), 0., 0., charmHadDauTracks[1].hasTOF(), charmHadDauTracks[1].hasTPC());
          if constexpr (WithMl) {
            std::array<float, 6> mlScores = {-1.f, -1.f, -1.f, -1.f, -1.f, -1.f};
            if (candC.mlProbD0().size() == NSizeMLScore) {
              std::copy(candC.mlProbD0().begin(), candC.mlProbD0().end(), mlScores.begin());
            }
            if (candC.mlProbD0bar().size() == NSizeMLScore) {
              std::copy(candC.mlProbD0bar().begin(), candC.mlProbD0bar().end(), mlScores.begin() + 3);
            }
            tables.hfCand2ProngMl(mlScores[0], mlScores[1], mlScores[2], mlScores[3], mlScores[4], mlScores[5]);
          }
        } else if constexpr (DecChannel == DecayChannel::B0ToDstarPi) {
          tables.hfCand2Prong(charmHadDauTracks[0].globalIndex(), charmHadDauTracks[1].globalIndex(),
                              indexHfReducedCollision,
                              trackParCovCharmHad.getX(), trackParCovCharmHad.getAlpha(),
                              trackParCovCharmHad.getY(), trackParCovCharmHad.getZ(), trackParCovCharmHad.getSnp(),
                              trackParCovCharmHad.getTgl(), trackParCovCharmHad.getQ2Pt(),
                              candC.xSecondaryVertexD0(), candC.ySecondaryVertexD0(), candC.zSecondaryVertexD0(), invMassC0, invMassC1,
                              ptDauMin, etaDauMin, nItsClsDauMin, nTpcCrossRowsDauMin, chi2TpcDauMax);
          tables.hfCand2ProngCov(trackParCovCharmHad.getSigmaY2(), trackParCovCharmHad.getSigmaZY(), trackParCovCharmHad.getSigmaZ2(),
                                 trackParCovCharmHad.getSigmaSnpY(), trackParCovCharmHad.getSigmaSnpZ(),
                                 trackParCovCharmHad.getSigmaSnp2(), trackParCovCharmHad.getSigmaTglY(), trackParCovCharmHad.getSigmaTglZ(),
                                 trackParCovCharmHad.getSigmaTglSnp(), trackParCovCharmHad.getSigmaTgl2(),
                                 trackParCovCharmHad.getSigma1PtY(), trackParCovCharmHad.getSigma1PtZ(), trackParCovCharmHad.getSigma1PtSnp(),
                                 trackParCovCharmHad.getSigma1PtTgl(), trackParCovCharmHad.getSigma1Pt2());
          float nSigmaTpcPr0{-999.f}, nSigmaTpcPr1{-999.f};
          float nSigmaTofPr0{-999.f}, nSigmaTofPr1{-999.f};
          tables.hfCandPidProng0(candC.nSigTpcPi0(), candC.nSigTofPi0(), candC.nSigTpcKa0(), candC.nSigTofKa0(), nSigmaTpcPr0, nSigmaTofPr0, charmHadDauTracks[0].hasTOF(), charmHadDauTracks[0].hasTPC());
          tables.hfCandPidProng1(candC.nSigTpcPi1(), candC.nSigTofPi1(), candC.nSigTpcKa1(), candC.nSigTofKa1(), nSigmaTpcPr1, nSigmaTofPr1, charmHadDauTracks[1].hasTOF(), charmHadDauTracks[1].hasTPC());

          // Soft pion tables
          auto trackSoftPion = charmHadDauTracks.back();
          auto trackParCovSoftPion = getTrackParCov(trackSoftPion);
          std::array<float, 2> dcaSoftPion{trackSoftPion.dcaXY(), trackSoftPion.dcaZ()};
          std::array<float, 3> pVecSoftPion = trackSoftPion.pVector();
          if (trackSoftPion.collisionId() != thisCollId) {
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovSoftPion, 2.f, noMatCorr, &dcaSoftPion);
            getPxPyPz(trackParCovSoftPion, pVecSoftPion);
          }
          tables.hfTrackSoftPion(trackSoftPion.globalIndex(), indexHfReducedCollision,
                                 trackParCovSoftPion.getX(), trackParCovSoftPion.getAlpha(),
                                 trackParCovSoftPion.getY(), trackParCovSoftPion.getZ(), trackParCovSoftPion.getSnp(),
                                 trackParCovSoftPion.getTgl(), trackParCovSoftPion.getQ2Pt(),
                                 trackSoftPion.itsNCls(), trackSoftPion.tpcNClsCrossedRows(), trackSoftPion.tpcChi2NCl());
          tables.hfTrackCovSoftPion(trackParCovSoftPion.getSigmaY2(), trackParCovSoftPion.getSigmaZY(), trackParCovSoftPion.getSigmaZ2(),
                                    trackParCovSoftPion.getSigmaSnpY(), trackParCovSoftPion.getSigmaSnpZ(),
                                    trackParCovSoftPion.getSigmaSnp2(), trackParCovSoftPion.getSigmaTglY(), trackParCovSoftPion.getSigmaTglZ(),
                                    trackParCovSoftPion.getSigmaTglSnp(), trackParCovSoftPion.getSigmaTgl2(),
                                    trackParCovSoftPion.getSigma1PtY(), trackParCovSoftPion.getSigma1PtZ(), trackParCovSoftPion.getSigma1PtSnp(),
                                    trackParCovSoftPion.getSigma1PtTgl(), trackParCovSoftPion.getSigma1Pt2());
          tables.hfTrackPidSoftPion(candC.nSigTpcPi2(), candC.nSigTofPi2(), charmHadDauTracks[2].hasTOF(), charmHadDauTracks[2].hasTPC());
          if constexpr (WithMl) {
            std::array<float, 6> mlScores = {-1.f, -1.f, -1.f, -1.f, -1.f, -1.f};
            if (candC.mlProbDstarToD0Pi().size() == NSizeMLScore) {
              std::copy(candC.mlProbDstarToD0Pi().begin(), candC.mlProbDstarToD0Pi().end(), mlScores.begin());
            }
            tables.hfCand3ProngMl(mlScores[0], mlScores[1], mlScores[2], mlScores[3], mlScores[4], mlScores[5]);
          }

          tables.hfMomDMesDaugs(pVec0[0], pVec0[1], pVec0[2],
                                pVec1[0], pVec1[1], pVec1[2],
                                pVecSoftPion[0], pVecSoftPion[1], pVecSoftPion[2]);
        }
        fillHfReducedCollision = true;
      }
    } // candsC loop

    if (!fillHfReducedCollision) {
      registry.fill(HIST("hEvents"), 1 + Event::NoCharmHadPiSelected);
      return;
    }
    registry.fill(HIST("hEvents"), 1 + Event::CharmHadPiSelected);

    // fill collision table if it contains a DPi pair a minima
    tables.hfReducedCollision(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), hfRejMap, bz);
    tables.hfReducedCollExtra(collision.covXX(), collision.covXY(), collision.covYY(),
                              collision.covXZ(), collision.covYZ(), collision.covZZ());
    tables.hfReducedCollCentrality(collision.centFT0C(), collision.centFT0M(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange());
    if constexpr (WithQvec) {
      tables.hfReducedQvector(collision.qvecFT0CRe(), collision.qvecFT0CIm(), collision.sumAmplFT0C(),
                              collision.qvecFT0ARe(), collision.qvecFT0AIm(), collision.sumAmplFT0A(),
                              collision.qvecFT0MRe(), collision.qvecFT0MIm(), collision.sumAmplFT0M(),
                              collision.qvecTPCposRe(), collision.qvecTPCposIm(), collision.nTrkTPCpos(),
                              collision.qvecTPCnegRe(), collision.qvecTPCnegIm(), collision.nTrkTPCneg(),
                              collision.qvecTPCallRe(), collision.qvecTPCallIm(), collision.nTrkTPCall());
    }
  }

  template <uint8_t DecayChannel>
  void runMcGen(aod::McCollision const& mcCollision,
                aod::McParticles const& particlesMc,
                CollisionsWCentAndMcLabels const& collisions,
                BCsInfo const&)
  {
    // Check event selection
    float centDummy{-1.f}, centFT0C{-1.f}, centFT0M{-1.f};
    const auto collSlice = collisions.sliceBy(preslices.colPerMcCollision, mcCollision.globalIndex());
    const auto hfRejMap = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, o2::hf_centrality::CentralityEstimator::None>(mcCollision, collSlice, centDummy);
    if (configs.skipRejectedCollisions && hfRejMap != 0) {
      return;
    }

    // get centrality
    using TMult = uint16_t; // type of numContrib
    TMult multiplicity{};
    for (const auto& collision : collSlice) {
      const TMult collMult = collision.numContrib();
      if (collMult > multiplicity) {
        centFT0C = collision.centFT0C();
        centFT0M = collision.centFT0M();
        multiplicity = collMult;
      }
    }

    const auto mcParticlesPerMcColl = particlesMc.sliceBy(preslices.mcParticlesPerMcCollision, mcCollision.globalIndex());

    // Match generated particles.
    for (const auto& particle : mcParticlesPerMcColl) {
      int8_t sign{0};
      int8_t flag{0};
      if constexpr (DecayChannel == DecayChannel::B0ToDminusPi) {
        // B0 → D- π+
        if (RecoDecay::isMatchedMCGen<true>(particlesMc, particle, Pdg::kB0, std::array{-static_cast<int>(Pdg::kDPlus), +kPiPlus}, true)) {
          // Match D- -> π- K+ π-
          auto candCMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
          // Printf("Checking D- -> π- K+ π-");
          if (RecoDecay::isMatchedMCGen(particlesMc, candCMC, -static_cast<int>(Pdg::kDPlus), std::array{-kPiPlus, +kKPlus, -kPiPlus}, true, &sign, 2)) {
            flag = sign * BIT(hf_cand_b0::DecayType::B0ToDPi);
          }
        }

        // save information for B0 task
        if (!TESTBIT(std::abs(flag), hf_cand_b0::DecayType::B0ToDPi)) {
          continue;
        }

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), massB);
        auto etaParticle = particle.eta();

        std::array<float, 2> ptProngs{};
        std::array<float, 2> yProngs{};
        std::array<float, 2> etaProngs{};
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }
        tables.rowHfB0McGenReduced(flag, -1 /*channel*/, ptParticle, yParticle, etaParticle,
                                   ptProngs[0], yProngs[0], etaProngs[0],
                                   ptProngs[1], yProngs[1], etaProngs[1], hfRejMap, centFT0C, centFT0M);
      } else if constexpr (DecayChannel == DecayChannel::BsToDsminusPi) {
        // Bs → Ds- π+
        if (RecoDecay::isMatchedMCGen<true>(particlesMc, particle, Pdg::kBS, std::array{-static_cast<int>(Pdg::kDS), +kPiPlus}, true)) {
          // Match Ds- -> π- K+ π-
          auto candCMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(particlesMc, candCMC, -static_cast<int>(Pdg::kDS), std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2)) {
            std::vector<int> arrDaughDsIndex;
            std::array<int, 2> arrPDGDaughDs{};
            RecoDecay::getDaughters(candCMC, &arrDaughDsIndex, std::array{0}, 1);
            if (arrDaughDsIndex.size() == NDaughtersDs) {
              for (auto jProng = 0u; jProng < arrDaughDsIndex.size(); ++jProng) {
                auto daughJ = particlesMc.rawIteratorAt(arrDaughDsIndex[jProng]);
                arrPDGDaughDs[jProng] = std::abs(daughJ.pdgCode());
              }
              // Ds- → Phi π- → K- K+ π- and Ds- → K0* K- → K- K+ π-
              if ((arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[0] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[1]) || (arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[1] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[0])) {
                flag = sign * BIT(hf_cand_bs::DecayTypeMc::BsToDsPiToPhiPiPiToKKPiPi);
              } else if ((arrPDGDaughDs[0] == arrPDGResonantDKstarK[0] && arrPDGDaughDs[1] == arrPDGResonantDKstarK[1]) || (arrPDGDaughDs[0] == arrPDGResonantDKstarK[1] && arrPDGDaughDs[1] == arrPDGResonantDKstarK[0])) {
                flag = sign * BIT(hf_cand_bs::DecayTypeMc::BsToDsPiToK0starKPiToKKPiPi);
              }
            }
          }
        }

        // additional checks for correlated backgrounds
        if (configs.checkDecayTypeMc) {
          // B0 → Ds- π+
          if (!flag) {
            if (RecoDecay::isMatchedMCGen<true>(particlesMc, particle, Pdg::kB0, std::array{-static_cast<int>(Pdg::kDS), +kPiPlus}, true)) {
              // Match Ds- -> π- K+ π-
              auto candCMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
              if (RecoDecay::isMatchedMCGen(particlesMc, candCMC, -static_cast<int>(Pdg::kDS), std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2)) {
                std::vector<int> arrDaughDsIndex;
                std::array<int, 2> arrPDGDaughDs{};
                RecoDecay::getDaughters(candCMC, &arrDaughDsIndex, std::array{0}, 1);
                if (arrDaughDsIndex.size() == NDaughtersDs) {
                  for (auto jProng = 0u; jProng < arrDaughDsIndex.size(); ++jProng) {
                    auto daughJ = particlesMc.rawIteratorAt(arrDaughDsIndex[jProng]);
                    arrPDGDaughDs[jProng] = std::abs(daughJ.pdgCode());
                  }
                  // Ds- → Phi π- → K- K+ π- and Ds- → K0* K- → K- K+ π-
                  if ((arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[0] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[1]) || (arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[1] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[0])) {
                    flag = sign * BIT(hf_cand_bs::DecayTypeMc::B0ToDsPiToPhiPiPiToKKPiPi);
                  } else if ((arrPDGDaughDs[0] == arrPDGResonantDKstarK[0] && arrPDGDaughDs[1] == arrPDGResonantDKstarK[1]) || (arrPDGDaughDs[0] == arrPDGResonantDKstarK[1] && arrPDGDaughDs[1] == arrPDGResonantDKstarK[0])) {
                    flag = sign * BIT(hf_cand_bs::DecayTypeMc::B0ToDsPiToK0starKPiToKKPiPi);
                  }
                }
              }
            }
          }
        }

        // save information for Bs task
        if (!TESTBIT(std::abs(flag), hf_cand_bs::DecayTypeMc::BsToDsPiToPhiPiPiToKKPiPi) && !TESTBIT(std::abs(flag), hf_cand_bs::DecayTypeMc::BsToDsPiToK0starKPiToKKPiPi) &&
            !TESTBIT(std::abs(flag), hf_cand_bs::DecayTypeMc::B0ToDsPiToPhiPiPiToKKPiPi) && !TESTBIT(std::abs(flag), hf_cand_bs::DecayTypeMc::B0ToDsPiToK0starKPiToKKPiPi)) {
          continue;
        }

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), massB);
        auto etaParticle = particle.eta();

        std::array<float, 2> ptProngs{};
        std::array<float, 2> yProngs{};
        std::array<float, 2> etaProngs{};
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }
        tables.rowHfBsMcGenReduced(flag, -1 /*channel*/, ptParticle, yParticle, etaParticle,
                                   ptProngs[0], yProngs[0], etaProngs[0],
                                   ptProngs[1], yProngs[1], etaProngs[1], hfRejMap, centFT0C, centFT0M);
      } else if constexpr (DecayChannel == DecayChannel::BplusToD0barPi) {
        // B+ → D0bar π+
        if (RecoDecay::isMatchedMCGen(particlesMc, particle, Pdg::kBPlus, std::array{-static_cast<int>(Pdg::kD0), +kPiPlus}, true)) {
          // Match D0bar -> π- K+
          auto candD0MC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
          // Printf("Checking D0bar -> π- K+");
          if (RecoDecay::isMatchedMCGen(particlesMc, candD0MC, static_cast<int>(Pdg::kD0), std::array{+kPiPlus, -kKPlus}, true, &sign)) {
            flag = sign * BIT(hf_cand_bplus::DecayType::BplusToD0Pi);
          }
        }

        // save information for B+ task
        if (!TESTBIT(std::abs(flag), hf_cand_bplus::DecayType::BplusToD0Pi)) {
          continue;
        }

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), massB);
        auto etaParticle = particle.eta();

        std::array<float, 2> ptProngs{};
        std::array<float, 2> yProngs{};
        std::array<float, 2> etaProngs{};
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }
        tables.rowHfBpMcGenReduced(flag, -1 /*channel*/, ptParticle, yParticle, etaParticle,
                                   ptProngs[0], yProngs[0], etaProngs[0],
                                   ptProngs[1], yProngs[1], etaProngs[1], hfRejMap, centFT0C, centFT0M);
      } else if constexpr (DecayChannel == DecayChannel::LbToLcplusPi) {
        // Lb → Lc+ π-
        if (RecoDecay::isMatchedMCGen<true>(particlesMc, particle, Pdg::kLambdaB0, std::array{static_cast<int>(Pdg::kLambdaCPlus), -kPiPlus}, true)) {
          // Match Lc+ → p K- π+
          auto candCMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
          // Printf("Checking Lc+ → p K- π+");
          if (RecoDecay::isMatchedMCGen(particlesMc, candCMC, static_cast<int>(Pdg::kLambdaCPlus), std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2)) {
            flag = sign * BIT(hf_cand_lb::DecayType::LbToLcPi);
          }
        }

        // save information for Lc task
        if (!TESTBIT(std::abs(flag), hf_cand_lb::DecayType::LbToLcPi)) {
          continue;
        }

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), massB);
        auto etaParticle = particle.eta();

        std::array<float, 2> ptProngs{};
        std::array<float, 2> yProngs{};
        std::array<float, 2> etaProngs{};
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }
        tables.rowHfLbMcGenReduced(flag, ptParticle, yParticle, etaParticle,
                                   ptProngs[0], yProngs[0], etaProngs[0],
                                   ptProngs[1], yProngs[1], etaProngs[1], hfRejMap, centFT0C, centFT0M);
      } else if constexpr (DecayChannel == DecayChannel::B0ToDstarPi) {
        // B0 → D* π+
        if (RecoDecay::isMatchedMCGen<true>(particlesMc, particle, Pdg::kB0, std::array{-static_cast<int>(Pdg::kDStar), +kPiPlus}, true)) {
          // Match D- -> π- K+ π-
          auto candCMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
          // Printf("Checking D- -> π- K+ π-");
          if (RecoDecay::isMatchedMCGen(particlesMc, candCMC, +static_cast<int>(Pdg::kDStar), std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 3)) {
            flag = sign * BIT(hf_cand_b0::DecayType::B0ToDstarPi);
          }
        }

        // save information for B0 task
        if (!TESTBIT(std::abs(flag), hf_cand_b0::DecayType::B0ToDstarPi)) {
          continue;
        }

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), massB);
        auto etaParticle = particle.eta();

        std::array<float, 2> ptProngs{};
        std::array<float, 2> yProngs{};
        std::array<float, 2> etaProngs{};
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }
        tables.rowHfB0McGenReduced(flag, -1 /*channel*/, ptParticle, yParticle, etaParticle,
                                   ptProngs[0], yProngs[0], etaProngs[0],
                                   ptProngs[1], yProngs[1], etaProngs[1], hfRejMap, centFT0C, centFT0M);
      }
    } // gen
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // PROCESS FUNCTIONS FOR DATA

  void processDplusPiData(CollisionsWCent const& collisions,
                          CandsDplusFiltered const& candsC,
                          aod::TrackAssoc const& trackIndices,
                          TracksPidWithSel const& tracks,
                          aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigB0(hfflagConfigurations.selectionFlagDplus.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDplusPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::B0ToDminusPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDplusPiData, "Process DplusPi without MC info and without ML info", true);

  void processDplusPiDataWithMl(CollisionsWCent const& collisions,
                                CandsDplusFilteredWithMl const& candsC,
                                aod::TrackAssoc const& trackIndices,
                                TracksPidWithSel const& tracks,
                                aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigB0(hfflagConfigurations.selectionFlagDplus.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDplusPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::B0ToDminusPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDplusPiDataWithMl, "Process DplusPi without MC info and with ML info", false);

  void processDplusPiDataWithQvec(CollisionsWCentAndQvectors const& collisions,
                                  CandsDplusFiltered const& candsC,
                                  aod::TrackAssoc const& trackIndices,
                                  TracksPidWithSel const& tracks,
                                  aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigB0(hfflagConfigurations.selectionFlagDplus.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDplusPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::B0ToDminusPi, true>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDplusPiDataWithQvec, "Process DplusPi without MC info, without ML info and with Q-vectors", false);

  void processDplusPiDataWithMlAndQvec(CollisionsWCentAndQvectors const& collisions,
                                       CandsDplusFilteredWithMl const& candsC,
                                       aod::TrackAssoc const& trackIndices,
                                       TracksPidWithSel const& tracks,
                                       aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigB0(hfflagConfigurations.selectionFlagDplus.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDplusPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::B0ToDminusPi, true>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDplusPiDataWithMlAndQvec, "Process DplusPi without MC info, with ML info and with Q-vectors", false);

  void processDstarPiData(CollisionsWCent const& collisions,
                          CandsDstarFiltered const& candsC,
                          aod::TrackAssoc const& trackIndices,
                          TracksPidWithSel const& tracks,
                          aod::BCsWithTimestamps const& bcs,
                          aod::Hf2Prongs const&)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigB0(hfflagConfigurations.selectionFlagDstar.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDstarPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::B0ToDstarPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDstarPiData, "Process DstarPi without MC info and without ML info", false);

  void processDstarPiDataWithMl(CollisionsWCent const& collisions,
                                CandsDstarFilteredWithMl const& candsC,
                                aod::TrackAssoc const& trackIndices,
                                TracksPidWithSel const& tracks,
                                aod::BCsWithTimestamps const& bcs,
                                aod::Hf2Prongs const&)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigB0(hfflagConfigurations.selectionFlagDstar.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDstarPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::B0ToDstarPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDstarPiDataWithMl, "Process DstarPi without MC info and with ML info", false);

  void processDstarPiDataWithQvec(CollisionsWCentAndQvectors const& collisions,
                                  CandsDstarFiltered const& candsC,
                                  aod::TrackAssoc const& trackIndices,
                                  TracksPidWithSel const& tracks,
                                  aod::BCsWithTimestamps const& bcs,
                                  aod::Hf2Prongs const&)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigB0(hfflagConfigurations.selectionFlagDstar.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDstarPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::B0ToDstarPi, true>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDstarPiDataWithQvec, "Process DstarPi without MC info, without ML info and with Q-vectors", false);

  void processDstarPiDataWithMlAndQvec(CollisionsWCentAndQvectors const& collisions,
                                       CandsDstarFilteredWithMl const& candsC,
                                       aod::TrackAssoc const& trackIndices,
                                       TracksPidWithSel const& tracks,
                                       aod::BCsWithTimestamps const& bcs,
                                       aod::Hf2Prongs const&)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigB0(hfflagConfigurations.selectionFlagDstar.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDstarPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::B0ToDstarPi, true>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDstarPiDataWithMlAndQvec, "Process DstarPi without MC info, with ML info and with Q-vectors", false);

  void processDsPiData(CollisionsWCent const& collisions,
                       CandsDsFiltered const& candsC,
                       aod::TrackAssoc const& trackIndices,
                       TracksPidWithSel const& tracks,
                       aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for Bs workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigBs(hfflagConfigurations.selectionFlagDs.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDsPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::BsToDsminusPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDsPiData, "Process DsPi without MC info and without ML info", false);

  void processDsPiDataWithMl(CollisionsWCent const& collisions,
                             CandsDsFilteredWithMl const& candsC,
                             aod::TrackAssoc const& trackIndices,
                             TracksPidWithSel const& tracks,
                             aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for Bs workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigBs(hfflagConfigurations.selectionFlagDs.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDsPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::BsToDsminusPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDsPiDataWithMl, "Process DsPi without MC info and with ML info", false);

  void processDsPiDataWithQvec(CollisionsWCentAndQvectors const& collisions,
                               CandsDsFiltered const& candsC,
                               aod::TrackAssoc const& trackIndices,
                               TracksPidWithSel const& tracks,
                               aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for Bs workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigBs(hfflagConfigurations.selectionFlagDs.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDsPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::BsToDsminusPi, true>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDsPiDataWithQvec, "Process DsPi without MC info, without ML info and with Q-vectors", false);

  void processDsPiDataWithMlAndQvec(CollisionsWCentAndQvectors const& collisions,
                                    CandsDsFilteredWithMl const& candsC,
                                    aod::TrackAssoc const& trackIndices,
                                    TracksPidWithSel const& tracks,
                                    aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for Bs workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigBs(hfflagConfigurations.selectionFlagDs.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDsPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::BsToDsminusPi, true>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDsPiDataWithMlAndQvec, "Process DsPi without MC info, with ML info and Q-vectors", false);

  void processD0PiData(CollisionsWCent const& collisions,
                       CandsD0Filtered const& candsC,
                       aod::TrackAssoc const& trackIndices,
                       TracksPidWithSel const& tracks,
                       aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigBplus(hfflagConfigurations.selectionFlagD0.value, hfflagConfigurations.selectionFlagD0bar.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsD0PerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::BplusToD0barPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processD0PiData, "Process D0Pi without MC info and without ML info", false);

  void processD0PiDataWithMl(CollisionsWCent const& collisions,
                             CandsD0FilteredWithMl const& candsC,
                             aod::TrackAssoc const& trackIndices,
                             TracksPidWithSel const& tracks,
                             aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigBplus(hfflagConfigurations.selectionFlagD0.value, hfflagConfigurations.selectionFlagD0bar.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsD0PerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::BplusToD0barPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processD0PiDataWithMl, "Process D0Pi without MC info and with ML info", false);

  void processD0PiDataWithQvec(CollisionsWCentAndQvectors const& collisions,
                               CandsD0Filtered const& candsC,
                               aod::TrackAssoc const& trackIndices,
                               TracksPidWithSel const& tracks,
                               aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigBplus(hfflagConfigurations.selectionFlagD0.value, hfflagConfigurations.selectionFlagD0bar.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsD0PerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::BplusToD0barPi, true>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processD0PiDataWithQvec, "Process D0Pi without MC info, without ML info, and with Q-vectors", false);

  void processD0PiDataWithMlAndQvec(CollisionsWCentAndQvectors const& collisions,
                                    CandsD0FilteredWithMl const& candsC,
                                    aod::TrackAssoc const& trackIndices,
                                    TracksPidWithSel const& tracks,
                                    aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigBplus(hfflagConfigurations.selectionFlagD0.value, hfflagConfigurations.selectionFlagD0bar.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsD0PerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::BplusToD0barPi, true>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processD0PiDataWithMlAndQvec, "Process D0Pi without MC info, with ML info, and with Q-vectors", false);

  void processLcPiData(CollisionsWCent const& collisions,
                       CandsLcFiltered const& candsC,
                       aod::TrackAssoc const& trackIndices,
                       TracksPidWithSel const& tracks,
                       aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for Lb workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigLb(hfflagConfigurations.selectionFlagLc.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsLcPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::LbToLcplusPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processLcPiData, "Process LcPi without MC info and without ML info", false);

  void processLcPiDataWithMl(CollisionsWCent const& collisions,
                             CandsLcFilteredWithMl const& candsC,
                             aod::TrackAssoc const& trackIndices,
                             TracksPidWithSel const& tracks,
                             aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for Lb workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigLb(hfflagConfigurations.selectionFlagLc.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsLcPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::LbToLcplusPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processLcPiDataWithMl, "Process LcPi without MC info and with ML info", false);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // PROCESS FUNCTIONS FOR MC

  void processDplusPiMc(CollisionsWCentAndMcLabels const& collisions,
                        CandsDplusFiltered const& candsC,
                        aod::TrackAssoc const& trackIndices,
                        TracksPidWithSelAndMc const& tracks,
                        aod::McParticles const& particlesMc,
                        BCsInfo const& bcs,
                        McCollisions const& mcCollisions)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigB0(hfflagConfigurations.selectionFlagDplus.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDplusPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(preslices.colPerMcCollision, collision.mcCollisionId());
      int64_t const indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, false, DecayChannel::B0ToDminusPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    for (const auto& mcCollision : mcCollisions) {
      runMcGen<DecayChannel::B0ToDminusPi>(mcCollision, particlesMc, collisions, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDplusPiMc, "Process DplusPi with MC info and without ML info", false);

  void processDplusPiMcWithMl(CollisionsWCentAndMcLabels const& collisions,
                              CandsDplusFilteredWithMl const& candsC,
                              aod::TrackAssoc const& trackIndices,
                              TracksPidWithSelAndMc const& tracks,
                              aod::McParticles const& particlesMc,
                              BCsInfo const& bcs,
                              McCollisions const& mcCollisions)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigB0(hfflagConfigurations.selectionFlagDplus.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDplusPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(preslices.colPerMcCollision, collision.mcCollisionId());
      int64_t const indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, true, DecayChannel::B0ToDminusPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    for (const auto& mcCollision : mcCollisions) {
      runMcGen<DecayChannel::B0ToDminusPi>(mcCollision, particlesMc, collisions, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDplusPiMcWithMl, "Process DplusPi with MC info and with ML info", false);

  void processDstarPiMc(CollisionsWCentAndMcLabels const& collisions,
                        CandsDstarFiltered const& candsC,
                        aod::TrackAssoc const& trackIndices,
                        TracksPidWithSelAndMc const& tracks,
                        aod::McParticles const& particlesMc,
                        BCsInfo const& bcs,
                        McCollisions const& mcCollisions,
                        aod::Hf2Prongs const&)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigB0(hfflagConfigurations.selectionFlagDstar.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDstarPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(preslices.colPerMcCollision, collision.mcCollisionId());
      int64_t const indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, false, DecayChannel::B0ToDstarPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    for (const auto& mcCollision : mcCollisions) {
      runMcGen<DecayChannel::B0ToDstarPi>(mcCollision, particlesMc, collisions, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDstarPiMc, "Process DstarPi with MC info and without ML info", false);

  void processDstarPiMcWithMl(CollisionsWCentAndMcLabels const& collisions,
                              CandsDstarFilteredWithMl const& candsC,
                              aod::TrackAssoc const& trackIndices,
                              TracksPidWithSelAndMc const& tracks,
                              aod::McParticles const& particlesMc,
                              BCsInfo const& bcs,
                              McCollisions const& mcCollisions,
                              aod::Hf2Prongs const&)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigB0(hfflagConfigurations.selectionFlagDstar.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDstarPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(preslices.colPerMcCollision, collision.mcCollisionId());
      int64_t const indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, true, DecayChannel::B0ToDstarPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    for (const auto& mcCollision : mcCollisions) {
      runMcGen<DecayChannel::B0ToDstarPi>(mcCollision, particlesMc, collisions, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDstarPiMcWithMl, "Process DstarPi with MC info and with ML info", false);

  void processDsPiMc(CollisionsWCentAndMcLabels const& collisions,
                     CandsDsFiltered const& candsC,
                     aod::TrackAssoc const& trackIndices,
                     TracksPidWithSelAndMc const& tracks,
                     aod::McParticles const& particlesMc,
                     BCsInfo const& bcs,
                     McCollisions const& mcCollisions)
  {
    // store configurables needed for Bs workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigBs(hfflagConfigurations.selectionFlagDs.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDsPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(preslices.colPerMcCollision, collision.mcCollisionId());
      int64_t const indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, false, DecayChannel::BsToDsminusPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    for (const auto& mcCollision : mcCollisions) {
      runMcGen<DecayChannel::BsToDsminusPi>(mcCollision, particlesMc, collisions, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDsPiMc, "Process DsPi with MC info and without ML info", false);

  void processDsPiMcWithMl(CollisionsWCentAndMcLabels const& collisions,
                           CandsDsFilteredWithMl const& candsC,
                           aod::TrackAssoc const& trackIndices,
                           TracksPidWithSelAndMc const& tracks,
                           aod::McParticles const& particlesMc,
                           BCsInfo const& bcs,
                           McCollisions const& mcCollisions)
  {
    // store configurables needed for Bs workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigBs(hfflagConfigurations.selectionFlagDs.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsDsPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(preslices.colPerMcCollision, collision.mcCollisionId());
      int64_t const indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, true, DecayChannel::BsToDsminusPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    for (const auto& mcCollision : mcCollisions) {
      runMcGen<DecayChannel::BsToDsminusPi>(mcCollision, particlesMc, collisions, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDsPiMcWithMl, "Process DsPi with MC info and with ML info", false);

  void processD0PiMc(CollisionsWCentAndMcLabels const& collisions,
                     CandsD0Filtered const& candsC,
                     aod::TrackAssoc const& trackIndices,
                     TracksPidWithSelAndMc const& tracks,
                     aod::McParticles const& particlesMc,
                     BCsInfo const& bcs,
                     McCollisions const& mcCollisions)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigBplus(hfflagConfigurations.selectionFlagD0.value, hfflagConfigurations.selectionFlagD0bar.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsD0PerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(preslices.colPerMcCollision, collision.mcCollisionId());
      int64_t const indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, false, DecayChannel::BplusToD0barPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    for (const auto& mcCollision : mcCollisions) {
      runMcGen<DecayChannel::BplusToD0barPi>(mcCollision, particlesMc, collisions, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processD0PiMc, "Process D0Pi with MC info and without ML info", false);

  void processD0PiMcWithMl(CollisionsWCentAndMcLabels const& collisions,
                           CandsD0FilteredWithMl const& candsC,
                           aod::TrackAssoc const& trackIndices,
                           TracksPidWithSelAndMc const& tracks,
                           aod::McParticles const& particlesMc,
                           BCsInfo const& bcs,
                           McCollisions const& mcCollisions)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigBplus(hfflagConfigurations.selectionFlagD0.value, hfflagConfigurations.selectionFlagD0bar.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsD0PerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(preslices.colPerMcCollision, collision.mcCollisionId());
      int64_t const indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, true, DecayChannel::BplusToD0barPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    for (const auto& mcCollision : mcCollisions) {
      runMcGen<DecayChannel::BplusToD0barPi>(mcCollision, particlesMc, collisions, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processD0PiMcWithMl, "Process D0Pi with MC info and with ML info", false);

  void processLcPiMc(CollisionsWCentAndMcLabels const& collisions,
                     CandsLcFiltered const& candsC,
                     aod::TrackAssoc const& trackIndices,
                     TracksPidWithSelAndMc const& tracks,
                     aod::McParticles const& particlesMc,
                     BCsInfo const& bcs,
                     McCollisions const& mcCollisions)
  {
    // store configurables needed for Lb workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigLb(hfflagConfigurations.selectionFlagDplus.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsLcPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(preslices.colPerMcCollision, collision.mcCollisionId());
      int64_t const indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, false, DecayChannel::LbToLcplusPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    for (const auto& mcCollision : mcCollisions) {
      runMcGen<DecayChannel::LbToLcplusPi>(mcCollision, particlesMc, collisions, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processLcPiMc, "Process LcPi with MC info and without ML info", false);

  void processLcPiMcWithMl(CollisionsWCentAndMcLabels const& collisions,
                           CandsLcFilteredWithMl const& candsC,
                           aod::TrackAssoc const& trackIndices,
                           TracksPidWithSelAndMc const& tracks,
                           aod::McParticles const& particlesMc,
                           BCsInfo const& bcs,
                           McCollisions const& mcCollisions)
  {
    // store configurables needed for Lb workflow
    if (!isHfCandBhadConfigFilled) {
      tables.rowCandidateConfigLb(hfflagConfigurations.selectionFlagLc.value, configs.invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(preslices.candsLcPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(preslices.trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(preslices.colPerMcCollision, collision.mcCollisionId());
      int64_t const indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, true, DecayChannel::LbToLcplusPi, false>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    }
    // handle normalization by the right number of collisions
    tables.hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    for (const auto& mcCollision : mcCollisions) {
      runMcGen<DecayChannel::LbToLcplusPi>(mcCollision, particlesMc, collisions, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processLcPiMcWithMl, "Process LcPi with MC info and with ML info", false);

}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDataCreatorCharmHadPiReduced>(cfgc)};
}
