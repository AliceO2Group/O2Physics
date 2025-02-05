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

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"

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
  BsToDsminusPi
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
  Produces<aod::HfRedCollisions> hfReducedCollision;
  Produces<aod::HfRedCollExtras> hfReducedCollExtra;
  Produces<aod::HfOrigColCounts> hfCollisionCounter;
  // Pi bachelor related tables
  Produces<aod::HfRedTrackBases> hfTrackPion;
  Produces<aod::HfRedTracksCov> hfTrackCovPion;
  Produces<aod::HfRedTracksPid> hfTrackPidPion;
  // charm hadron related tables
  Produces<aod::HfRed2Prongs> hfCand2Prong;
  Produces<aod::HfRed2ProngsCov> hfCand2ProngCov;
  Produces<aod::HfRed2ProngsMl> hfCand2ProngMl;
  Produces<aod::HfRed3Prongs> hfCand3Prong;
  Produces<aod::HfRed3ProngsCov> hfCand3ProngCov;
  Produces<aod::HfRed3ProngsMl> hfCand3ProngMl;
  // PID tables for charm-hadron candidate daughter tracks
  Produces<aod::HfRedPidDau0s> hfCandPidProng0;
  Produces<aod::HfRedPidDau1s> hfCandPidProng1;
  Produces<aod::HfRedPidDau2s> hfCandPidProng2;

  // B-hadron config and MC related tables
  Produces<aod::HfCandB0Configs> rowCandidateConfigB0;
  Produces<aod::HfMcRecRedDpPis> rowHfDPiMcRecReduced;
  Produces<aod::HfMcCheckDpPis> rowHfDPiMcCheckReduced;
  Produces<aod::HfMcGenRedB0s> rowHfB0McGenReduced;

  Produces<aod::HfCandBpConfigs> rowCandidateConfigBplus;
  Produces<aod::HfMcRecRedD0Pis> rowHfD0PiMcRecReduced;
  Produces<aod::HfMcCheckD0Pis> rowHfD0PiMcCheckReduced;
  Produces<aod::HfMcGenRedBps> rowHfBpMcGenReduced;

  Produces<aod::HfCandBsConfigs> rowCandidateConfigBs;
  Produces<aod::HfMcRecRedDsPis> rowHfDsPiMcRecReduced;
  Produces<aod::HfMcCheckDsPis> rowHfDsPiMcCheckReduced;
  Produces<aod::HfMcGenRedBss> rowHfBsMcGenReduced;

  // vertexing
  // Configurable<double> bz{"bz", 5., "magnetic field"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any B0 is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // selection
  Configurable<bool> usePionIsGlobalTrackWoDCA{"usePionIsGlobalTrackWoDCA", true, "check isGlobalTrackWoDCA status for pions, for Run3 studies"};
  Configurable<double> ptPionMin{"ptPionMin", 0.5, "minimum pion pT threshold (GeV/c)"};
  Configurable<std::vector<double>> binsPtPion{"binsPtPion", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for pion DCA XY pT-dependent cut"};
  Configurable<LabeledArray<double>> cutsTrackPionDCA{"cutsTrackPionDCA", {hf_cuts_single_track::cutsTrack[0], hf_cuts_single_track::nBinsPtTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for pions"};
  Configurable<double> invMassWindowCharmHadPi{"invMassWindowCharmHadPi", 0.3, "invariant-mass window for CharmHad-Pi pair preselections (GeV/c2)"};
  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for D+"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};

  // magnetic field setting from CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  // MC extra
  Configurable<bool> checkDecayTypeMc{"checkDecayTypeMc", false, "flag to enable MC checks on decay type"};

  HfHelper hfHelper;
  o2::hf_evsel::HfEventSelection hfEvSel;

  // CCDB service
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber;

  // O2DatabasePDG service
  Service<o2::framework::O2DatabasePDG> pdg;

  double massPi{0.};
  double massC{0.};
  double massB{0.};
  double invMass2ChHadPiMin{0.};
  double invMass2ChHadPiMax{0.};
  double bz{0.};

  bool isHfCandBhadConfigFilled = false;

  // Fitter to redo D-vertex to get extrapolated daughter tracks (2/3-prong vertex filter)
  o2::vertexing::DCAFitterN<3> df3;
  o2::vertexing::DCAFitterN<2> df2;

  using TracksPid = soa::Join<aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa>; // TODO: revert to pion only once the Nsigma variables for the charm-hadron candidate daughters are in the candidate table for 3 prongs too
  using TracksPidWithSel = soa::Join<aod::TracksWCovDcaExtra, TracksPid, aod::TrackSelection>;
  using TracksPidWithSelAndMc = soa::Join<TracksPidWithSel, aod::McTrackLabels>;

  using CandsDplusFiltered = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using CandsDplusFilteredWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandsDsFiltered = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandsDsFilteredWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandsD0Filtered = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0>>;
  using CandsD0FilteredWithMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0, aod::HfMlD0>>;

  using CollisionsWMcLabels = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;

  Filter filterSelectDplusCandidates = (aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus);
  Filter filterSelectDsCandidates = (aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs);
  Filter filterSelectDzeroCandidates = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar);

  Preslice<CandsDplusFiltered> candsDplusPerCollision = aod::track_association::collisionId;
  Preslice<CandsDplusFilteredWithMl> candsDplusPerCollisionWithMl = aod::track_association::collisionId;
  Preslice<CandsDsFiltered> candsDsPerCollision = aod::track_association::collisionId;
  Preslice<CandsDsFilteredWithMl> candsDsPerCollisionWithMl = aod::track_association::collisionId;
  Preslice<CandsD0Filtered> candsD0PerCollision = aod::track_association::collisionId;
  Preslice<CandsD0FilteredWithMl> candsD0PerCollisionWithMl = aod::track_association::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  PresliceUnsorted<CollisionsWMcLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  std::shared_ptr<TH1> hCandidatesD0, hCandidatesDPlus, hCandidatesDs;
  HistogramRegistry registry{"registry"};

  std::array<int, 2> arrPDGResonantDsPhiPi = {kPhi, kPiPlus};      // Ds± → Phi π±
  std::array<int, 2> arrPDGResonantDKstarK = {kK0Star892, kKPlus}; // Ds± → K*(892)0bar K± and D± → K*(892)0bar K±

  void init(InitContext const&)
  {
    std::array<int, 12> doProcess = {doprocessDplusPiData, doprocessDplusPiDataWithMl, doprocessDplusPiMc, doprocessDplusPiMcWithMl,
                                     doprocessDsPiData, doprocessDsPiDataWithMl, doprocessDsPiMc, doprocessDsPiMcWithMl,
                                     doprocessD0PiData, doprocessD0PiDataWithMl, doprocessD0PiMc, doprocessD0PiMcWithMl};
    if (std::accumulate(doProcess.begin(), doProcess.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function can be enabled at a time, please fix your configuration!");
    }

    // invariant-mass window cut
    massPi = MassPiPlus;
    if (doprocessDplusPiData || doprocessDplusPiDataWithMl || doprocessDplusPiMc || doprocessDplusPiMcWithMl) {
      massC = MassDMinus;
      massB = MassB0;
    } else if (doprocessDsPiData || doprocessDsPiDataWithMl || doprocessDsPiMc || doprocessDsPiMcWithMl) {
      massC = MassDS;
      massB = MassBS;
    } else if (doprocessD0PiData || doprocessD0PiDataWithMl || doprocessD0PiMc || doprocessD0PiMcWithMl) {
      massC = MassD0;
      massB = MassBPlus;
    }
    invMass2ChHadPiMin = (massB - invMassWindowCharmHadPi) * (massB - invMassWindowCharmHadPi);
    invMass2ChHadPiMax = (massB + invMassWindowCharmHadPi) * (massB + invMassWindowCharmHadPi);

    // Initialize fitter
    if (doprocessDplusPiData || doprocessDplusPiDataWithMl || doprocessDplusPiMc || doprocessDplusPiMcWithMl ||
        doprocessDsPiData || doprocessDsPiDataWithMl || doprocessDsPiMc || doprocessDsPiMcWithMl) {
      df3.setPropagateToPCA(propagateToPCA);
      df3.setMaxR(maxR);
      df3.setMaxDZIni(maxDZIni);
      df3.setMinParamChange(minParamChange);
      df3.setMinRelChi2Change(minRelChi2Change);
      df3.setUseAbsDCA(useAbsDCA);
      df3.setWeightedFinalPCA(useWeightedFinalPCA);
      df3.setMatCorrType(noMatCorr);
    } else if (doprocessD0PiData || doprocessD0PiDataWithMl || doprocessD0PiMc || doprocessD0PiMcWithMl) {
      df2.setPropagateToPCA(propagateToPCA);
      df2.setMaxR(maxR);
      df2.setMaxDZIni(maxDZIni);
      df2.setMinParamChange(minParamChange);
      df2.setMinRelChi2Change(minRelChi2Change);
      df2.setUseAbsDCA(useAbsDCA);
      df2.setWeightedFinalPCA(useWeightedFinalPCA);
      df2.setMatCorrType(noMatCorr);
    }

    // Configure CCDB access
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    runNumber = 0;

    // histograms
    constexpr int kNBinsEvents = kNEvent;
    std::string labels[kNBinsEvents];
    labels[Event::Processed] = "processed";
    labels[Event::NoCharmHadPiSelected] = "without CharmHad-Pi pairs";
    labels[Event::CharmHadPiSelected] = "with CharmHad-Pi pairs";
    static const AxisSpec axisEvents = {kNBinsEvents, 0.5, kNBinsEvents + 0.5, ""};
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {axisEvents});
    for (int iBin = 0; iBin < kNBinsEvents; iBin++) {
      registry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    std::string charmHadTitle = "";
    std::string histMassTitle = "";
    if (doprocessDplusPiData || doprocessDplusPiDataWithMl || doprocessDplusPiMc || doprocessDplusPiMcWithMl) {
      charmHadTitle = "D^{#plus}";
      histMassTitle = "Dplus";
      registry.add("hMassDplus", "D^{#plus} candidates; #it{M}(K#pi#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
    } else if (doprocessDsPiData || doprocessDsPiDataWithMl || doprocessDsPiMc || doprocessDsPiMcWithMl) {
      charmHadTitle = "D_{s}^{#plus}";
      histMassTitle = "Ds";
      registry.add("hMassDsToKKPi", "D_{s}^{#plus} to KKpi candidates; #it{M}(KK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
      registry.add("hMassDsToPiKK", "D_{s}^{#plus} to piKK candidates; #it{M}(KK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
    } else if (doprocessD0PiData || doprocessD0PiDataWithMl || doprocessD0PiMc || doprocessD0PiMcWithMl) {
      charmHadTitle = "D^{0}";
      histMassTitle = "D0";
      registry.add("hMassD0", "D^{0} candidates; #it{M}(K#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
      registry.add("hMassD0bar", "#overline{D}^{0} candidates; #it{M}(K#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
    }

    registry.add(Form("hPt%s", histMassTitle.data()), Form("%s candidates candidates;%s candidate #it{p}_{T} (GeV/#it{c});entries", charmHadTitle.data(), charmHadTitle.data()), {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hPtPion", "#pi^{#plus} candidates;#pi^{#plus} candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add(Form("hCpa%s", histMassTitle.data()), Form("%s candidates;%s cosine of pointing angle;entries", charmHadTitle.data(), charmHadTitle.data()), {HistType::kTH1F, {{110, -1.1, 1.1}}});

    /// candidate monitoring
    hCandidatesD0 = registry.add<TH1>("hCandidatesD0", "D0 candidate counter", {HistType::kTH1D, {axisCands}});
    hCandidatesDPlus = registry.add<TH1>("hCandidatesDPlus", "Dplus candidate counter", {HistType::kTH1D, {axisCands}});
    hCandidatesDs = registry.add<TH1>("hCandidatesDs", "Ds candidate counter", {HistType::kTH1D, {axisCands}});
    setLabelHistoCands(hCandidatesD0);
    setLabelHistoCands(hCandidatesDPlus);
    setLabelHistoCands(hCandidatesDs);
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
    if (usePionIsGlobalTrackWoDCA && !trackPion.isGlobalTrackWoDCA()) {
      return false;
    }
    // minimum pT selection
    if (trackParCovPion.getPt() < ptPionMin || !isSelectedTrackDCA(trackParCovPion, dcaPion)) {
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

  /// Single-track cuts for pions on dcaXY
  /// \param trackPar is the track parametrisation
  /// \param dca is the 2-D array with track DCAs
  /// \return true if track passes all cuts
  template <typename T1, typename T2>
  bool isSelectedTrackDCA(const T1& trackPar, const T2& dca)
  {
    auto pTBinTrack = findBin(binsPtPion, trackPar.getPt());
    if (pTBinTrack == -1) {
      return false;
    }

    if (std::abs(dca[0]) < cutsTrackPionDCA->get(pTBinTrack, "min_dcaxytoprimary")) {
      return false; // minimum DCAxy
    }
    if (std::abs(dca[0]) > cutsTrackPionDCA->get(pTBinTrack, "max_dcaxytoprimary")) {
      return false; // maximum DCAxy
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
  template <uint8_t decChannel, typename CColl, typename PParticles, typename TTrack>
  void fillMcRecoInfo(const CColl& collision,
                      const PParticles& particlesMc,
                      const std::vector<TTrack>& vecDaughtersB,
                      int& indexHfCandCharm,
                      std::map<int64_t, int64_t> selectedTracksPion,
                      const int64_t indexCollisionMaxNumContrib)
  {

    // we check the MC matching to be stored
    int8_t sign{0};
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

    if constexpr (decChannel == DecayChannel::B0ToDminusPi) {
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
      if (checkDecayTypeMc) {
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
        // Partly reconstructed decays, i.e. the 4 prongs have a common b-hadron ancestor
        // convention: final state particles are prong0,1,2,3
        if (!flag) {
          auto particleProng0 = vecDaughtersB[0].mcParticle();
          auto particleProng1 = vecDaughtersB[1].mcParticle();
          auto particleProng2 = vecDaughtersB[2].mcParticle();
          auto particleProng3 = vecDaughtersB[3].mcParticle();
          // b-hadron hypothesis
          std::array<int, 3> bHadronMotherHypos = {Pdg::kB0, Pdg::kBS, Pdg::kLambdaB0};
          // c-hadron hypothesis
          std::array<int, 4> cHadronMotherHypos = {Pdg::kDPlus, Pdg::kDS, Pdg::kDStar, Pdg::kLambdaCPlus};

          for (const auto& bHadronMotherHypo : bHadronMotherHypos) {
            int index0Mother = RecoDecay::getMother(particlesMc, particleProng0, bHadronMotherHypo, true);
            int index1Mother = RecoDecay::getMother(particlesMc, particleProng1, bHadronMotherHypo, true);
            int index2Mother = RecoDecay::getMother(particlesMc, particleProng2, bHadronMotherHypo, true);
            int index3Mother = RecoDecay::getMother(particlesMc, particleProng3, bHadronMotherHypo, true);

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
                  int index0CharmMother = RecoDecay::getMother(particlesMc, particleProng0, cHadronMotherHypo, true, &sign, depthMax);
                  int index1CharmMother = RecoDecay::getMother(particlesMc, particleProng1, cHadronMotherHypo, true, &sign, depthMax);
                  int index2CharmMother = RecoDecay::getMother(particlesMc, particleProng2, cHadronMotherHypo, true, &sign, depthMax);
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
        rowHfDPiMcCheckReduced(pdgCodeBeautyMother, pdgCodeCharmMother, pdgCodeProng0, pdgCodeProng1, pdgCodeProng2, pdgCodeProng3);
      }
      rowHfDPiMcRecReduced(indexHfCandCharm, selectedTracksPion[vecDaughtersB.back().globalIndex()], flag, flagWrongCollision, debug, motherPt);
    } else if constexpr (decChannel == DecayChannel::BsToDsminusPi) {
      // Bs → Ds- π+ → (K- K+ π-) π+
      auto indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kBS, std::array{-kKPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 3);
      if (indexRec > -1) {
        // Ds- → K- K+ π-
        indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, -Pdg::kDS, std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
        if (indexRec > -1) {
          std::vector<int> arrDaughDsIndex;
          std::array<int, 2> arrPDGDaughDs;
          RecoDecay::getDaughters(particlesMc.rawIteratorAt(indexRec), &arrDaughDsIndex, std::array{0}, 1);
          if (arrDaughDsIndex.size() == 2) {
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
      if (checkDecayTypeMc) {
        // B0 → Ds- π+ → (K- K+ π-) π+
        if (!flag) {
          indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kB0, std::array{-kKPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 3);
          if (indexRec > -1) {
            // Ds- → K- K+ π-
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, -Pdg::kDS, std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
            if (indexRec > -1) {
              std::vector<int> arrDaughDsIndex;
              std::array<int, 2> arrPDGDaughDs;
              RecoDecay::getDaughters(particlesMc.rawIteratorAt(indexRec), &arrDaughDsIndex, std::array{0}, 1);
              if (arrDaughDsIndex.size() == 2) {
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
        // Partly reconstructed decays, i.e. the 4 prongs have a common b-hadron ancestor
        // convention: final state particles are prong0,1,2,3
        if (!flag) {
          auto particleProng0 = vecDaughtersB[0].mcParticle();
          auto particleProng1 = vecDaughtersB[1].mcParticle();
          auto particleProng2 = vecDaughtersB[2].mcParticle();
          auto particleProng3 = vecDaughtersB[3].mcParticle();
          // b-hadron hypothesis
          std::array<int, 3> bHadronMotherHypos = {Pdg::kB0, Pdg::kBS, Pdg::kLambdaB0};
          // c-hadron hypothesis
          std::array<int, 5> cHadronMotherHypos = {Pdg::kDPlus, Pdg::kDS, Pdg::kDStar, Pdg::kDSStar, Pdg::kLambdaCPlus};

          for (const auto& bHadronMotherHypo : bHadronMotherHypos) {
            int index0Mother = RecoDecay::getMother(particlesMc, particleProng0, bHadronMotherHypo, true);
            int index1Mother = RecoDecay::getMother(particlesMc, particleProng1, bHadronMotherHypo, true);
            int index2Mother = RecoDecay::getMother(particlesMc, particleProng2, bHadronMotherHypo, true);
            int index3Mother = RecoDecay::getMother(particlesMc, particleProng3, bHadronMotherHypo, true);

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
                  int index0CharmMother = RecoDecay::getMother(particlesMc, particleProng0, cHadronMotherHypo, true, &sign, depthMax);
                  int index1CharmMother = RecoDecay::getMother(particlesMc, particleProng1, cHadronMotherHypo, true, &sign, depthMax);
                  int index2CharmMother = RecoDecay::getMother(particlesMc, particleProng2, cHadronMotherHypo, true, &sign, depthMax);
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
        rowHfDsPiMcCheckReduced(pdgCodeBeautyMother, pdgCodeCharmMother, pdgCodeProng0, pdgCodeProng1, pdgCodeProng2, pdgCodeProng3);
      }
      rowHfDsPiMcRecReduced(indexHfCandCharm, selectedTracksPion[vecDaughtersB.back().globalIndex()], flag, flagWrongCollision, debug, motherPt);
    } else if constexpr (decChannel == DecayChannel::BplusToD0barPi) {
      // B+ → D0(bar) π+ → (K+ π-) π+
      auto indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, Pdg::kBPlus, std::array{+kPiPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
      if (indexRec > -1) {
        // D0(bar) → K+ π-;
        indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1]}, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign, 1);
        if (indexRec > -1) {
          flag = sign * BIT(hf_cand_bplus::DecayType::BplusToD0Pi);
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
      if (checkDecayTypeMc) {
        // Partly reconstructed decays, i.e. the 3 prongs have a common b-hadron ancestor
        // convention: final state particles are prong0,1,2
        if (!flag) {
          auto particleProng0 = vecDaughtersB[0].mcParticle();
          auto particleProng1 = vecDaughtersB[1].mcParticle();
          auto particleProng2 = vecDaughtersB[2].mcParticle();
          // b-hadron hypothesis
          std::array<int, 3> bHadronMotherHypos = {Pdg::kBPlus, Pdg::kB0, Pdg::kBS};

          for (const auto& bHadronMotherHypo : bHadronMotherHypos) {
            int index0Mother = RecoDecay::getMother(particlesMc, particleProng0, bHadronMotherHypo, true);
            int index1Mother = RecoDecay::getMother(particlesMc, particleProng1, bHadronMotherHypo, true);
            int index2Mother = RecoDecay::getMother(particlesMc, particleProng2, bHadronMotherHypo, true);

            // look for common b-hadron ancestor
            if (index0Mother > -1 && index1Mother > -1 && index2Mother > -1) {
              if (index0Mother == index1Mother && index1Mother == index2Mother) {
                flag = BIT(hf_cand_bplus::DecayTypeMc::PartlyRecoDecay);
                pdgCodeBeautyMother = particlesMc.rawIteratorAt(index0Mother).pdgCode();
                pdgCodeCharmMother = 0;
                pdgCodeProng0 = particleProng0.pdgCode();
                pdgCodeProng1 = particleProng1.pdgCode();
                pdgCodeProng2 = particleProng2.pdgCode();
                break;
              }
            }
          }
        }
        rowHfD0PiMcCheckReduced(pdgCodeBeautyMother, pdgCodeProng0, pdgCodeProng1, pdgCodeProng2);
      }
      rowHfD0PiMcRecReduced(indexHfCandCharm, selectedTracksPion[vecDaughtersB.back().globalIndex()], flag, flagWrongCollision, debug, motherPt);
    }
  }

  template <bool doMc, bool withMl, uint8_t decChannel, typename PParticles, typename TTracks, typename CCharmCands, typename Coll>
  void runDataCreation(Coll const& collision,
                       CCharmCands const& candsC,
                       aod::TrackAssoc const& trackIndices,
                       TTracks const&,
                       PParticles const& particlesMc,
                       uint64_t const& indexCollisionMaxNumContrib,
                       aod::BCsWithTimestamps const&)
  {
    // helpers for ReducedTables filling
    int indexHfReducedCollision = hfReducedCollision.lastIndex() + 1;
    // std::map where the key is the track.globalIndex() and
    // the value is the track index in the table of the selected pions
    std::map<int64_t, int64_t> selectedTracksPion;
    bool fillHfReducedCollision = false;

    auto primaryVertex = getPrimaryVertex(collision);

    // Set the magnetic field from ccdb.
    // The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
    // but this is not true when running on Run2 data/MC already converted into AO2Ds.
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrpMag, bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      runNumber = bc.runNumber();
    }
    df3.setBz(bz);

    auto thisCollId = collision.globalIndex();
    for (const auto& candC : candsC) {
      int indexHfCandCharm{-1};
      float invMassC0{-1.f}, invMassC1{-1.f};
      if constexpr (decChannel == DecayChannel::B0ToDminusPi) {
        indexHfCandCharm = hfCand3Prong.lastIndex() + 1;
        invMassC0 = hfHelper.invMassDplusToPiKPi(candC);
        registry.fill(HIST("hMassDplus"), invMassC0);
        registry.fill(HIST("hPtDplus"), candC.pt());
        registry.fill(HIST("hCpaDplus"), candC.cpa());
      } else if constexpr (decChannel == DecayChannel::BsToDsminusPi) {
        indexHfCandCharm = hfCand3Prong.lastIndex() + 1;
        if (candC.isSelDsToKKPi() >= selectionFlagDs) {
          invMassC0 = hfHelper.invMassDsToKKPi(candC);
          registry.fill(HIST("hMassDsToKKPi"), invMassC0);
        }
        if (candC.isSelDsToPiKK() >= selectionFlagDs) {
          invMassC1 = hfHelper.invMassDsToPiKK(candC);
          registry.fill(HIST("hMassDsToPiKK"), invMassC1);
        }
        registry.fill(HIST("hPtDs"), candC.pt());
        registry.fill(HIST("hCpaDs"), candC.cpa());
      } else if constexpr (decChannel == DecayChannel::BplusToD0barPi) {
        indexHfCandCharm = hfCand2Prong.lastIndex() + 1;
        if (candC.isSelD0() >= selectionFlagD0) {
          invMassC0 = hfHelper.invMassD0ToPiK(candC);
          registry.fill(HIST("hMassD0"), invMassC0);
        }
        if (candC.isSelD0bar() >= selectionFlagD0bar) {
          invMassC1 = hfHelper.invMassD0barToKPi(candC);
          registry.fill(HIST("hMassD0bar"), invMassC1);
        }
        registry.fill(HIST("hPtD0"), candC.pt());
        registry.fill(HIST("hCpaD0"), candC.cpa());
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
      if constexpr (decChannel == DecayChannel::B0ToDminusPi || decChannel == DecayChannel::BsToDsminusPi) {
        charmHadDauTracks.push_back(candC.template prong2_as<TTracks>());
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
      if constexpr (decChannel == DecayChannel::B0ToDminusPi || decChannel == DecayChannel::BsToDsminusPi) { // D∓ → π∓ K± π∓ and Ds∓ → K∓ K± π∓

        if constexpr (decChannel == DecayChannel::B0ToDminusPi) {
          hCandidatesDPlus->Fill(SVFitting::BeforeFit);
        } else {
          hCandidatesDs->Fill(SVFitting::BeforeFit);
        }

        try {
          if (df3.process(trackParCov0, trackParCov1, trackParCov2) == 0) {
            continue;
          }
        } catch (const std::runtime_error& error) {
          LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
          if constexpr (decChannel == DecayChannel::B0ToDminusPi) {
            hCandidatesDPlus->Fill(SVFitting::Fail);
          } else {
            hCandidatesDs->Fill(SVFitting::Fail);
          }
          continue;
        }
        if constexpr (decChannel == DecayChannel::B0ToDminusPi) {
          hCandidatesDPlus->Fill(SVFitting::FitOk);
        } else {
          hCandidatesDs->Fill(SVFitting::FitOk);
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
      } else if constexpr (decChannel == DecayChannel::BplusToD0barPi) { // D0(bar) → K± π∓

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
        o2::gpu::gpustd::array<float, 2> dcaPion{trackPion.dcaXY(), trackPion.dcaZ()};
        std::array<float, 3> pVecPion = trackPion.pVector();
        if (trackPion.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovPion, 2.f, noMatCorr, &dcaPion);
          getPxPyPz(trackParCovPion, pVecPion);
        }

        // reject pi D with same sign as D
        if constexpr (decChannel == DecayChannel::B0ToDminusPi || decChannel == DecayChannel::BsToDsminusPi) { // D∓ → π∓ K± π∓ and Ds∓ → K∓ K± π∓
          if (trackPion.sign() * charmHadDauTracks[0].sign() > 0) {
            continue;
          }
        } else if constexpr (decChannel == DecayChannel::BplusToD0barPi) { // D0(bar) → K± π∓
          if (!((candC.isSelD0() >= selectionFlagD0 && trackPion.sign() < 0) || (candC.isSelD0bar() >= selectionFlagD0bar && trackPion.sign() > 0))) {
            continue;
          }
        }

        // apply selections on pion tracks
        if (!isPionSelected(trackPion, trackParCovPion, dcaPion, charmHadDauTracks)) {
          continue;
        }

        registry.fill(HIST("hPtPion"), trackParCovPion.getPt());
        // compute invariant mass square and apply selection
        auto invMass2DPi = RecoDecay::m2(std::array{pVecCharm, pVecPion}, std::array{massC, massPi});
        if ((invMass2DPi < invMass2ChHadPiMin) || (invMass2DPi > invMass2ChHadPiMax)) {
          continue;
        }

        // fill Pion tracks table
        // if information on track already stored, go to next track
        if (!selectedTracksPion.count(trackPion.globalIndex())) {
          hfTrackPion(trackPion.globalIndex(), indexHfReducedCollision,
                      trackParCovPion.getX(), trackParCovPion.getAlpha(),
                      trackParCovPion.getY(), trackParCovPion.getZ(), trackParCovPion.getSnp(),
                      trackParCovPion.getTgl(), trackParCovPion.getQ2Pt(),
                      trackPion.itsNCls(), trackPion.tpcNClsCrossedRows(), trackPion.tpcChi2NCl());
          hfTrackCovPion(trackParCovPion.getSigmaY2(), trackParCovPion.getSigmaZY(), trackParCovPion.getSigmaZ2(),
                         trackParCovPion.getSigmaSnpY(), trackParCovPion.getSigmaSnpZ(),
                         trackParCovPion.getSigmaSnp2(), trackParCovPion.getSigmaTglY(), trackParCovPion.getSigmaTglZ(),
                         trackParCovPion.getSigmaTglSnp(), trackParCovPion.getSigmaTgl2(),
                         trackParCovPion.getSigma1PtY(), trackParCovPion.getSigma1PtZ(), trackParCovPion.getSigma1PtSnp(),
                         trackParCovPion.getSigma1PtTgl(), trackParCovPion.getSigma1Pt2());
          hfTrackPidPion(trackPion.hasTPC(), trackPion.hasTOF(),
                         trackPion.tpcNSigmaPi(), trackPion.tofNSigmaPi());
          // add trackPion.globalIndex() to a list
          // to keep memory of the pions filled in the table and avoid refilling them if they are paired to another D candidate
          // and keep track of their index in hfTrackPion for McRec purposes
          selectedTracksPion[trackPion.globalIndex()] = hfTrackPion.lastIndex();
        }

        if constexpr (doMc) {
          std::vector<typename TTracks::iterator> beautyHadDauTracks{};
          for (const auto& track : charmHadDauTracks) {
            beautyHadDauTracks.push_back(track);
          }
          beautyHadDauTracks.push_back(trackPion);
          fillMcRecoInfo<decChannel>(collision, particlesMc, beautyHadDauTracks, indexHfCandCharm, selectedTracksPion, indexCollisionMaxNumContrib);
        }
        fillHfCandCharm = true;
      }                                                           // pion loop
      if (fillHfCandCharm) {                                      // fill candCplus table only once per D candidate
        if constexpr (decChannel == DecayChannel::B0ToDminusPi || decChannel == DecayChannel::BsToDsminusPi) { // D∓ → π∓ K± π∓ and Ds∓ → K∓ K± π∓
          hfCand3Prong(charmHadDauTracks[0].globalIndex(), charmHadDauTracks[1].globalIndex(), charmHadDauTracks[2].globalIndex(),
                       indexHfReducedCollision,
                       trackParCovCharmHad.getX(), trackParCovCharmHad.getAlpha(),
                       trackParCovCharmHad.getY(), trackParCovCharmHad.getZ(), trackParCovCharmHad.getSnp(),
                       trackParCovCharmHad.getTgl(), trackParCovCharmHad.getQ2Pt(),
                       candC.xSecondaryVertex(), candC.ySecondaryVertex(), candC.zSecondaryVertex(), invMassC0, invMassC1,
                       ptDauMin, etaDauMin, nItsClsDauMin, nTpcCrossRowsDauMin, chi2TpcDauMax);
          hfCand3ProngCov(trackParCovCharmHad.getSigmaY2(), trackParCovCharmHad.getSigmaZY(), trackParCovCharmHad.getSigmaZ2(),
                          trackParCovCharmHad.getSigmaSnpY(), trackParCovCharmHad.getSigmaSnpZ(),
                          trackParCovCharmHad.getSigmaSnp2(), trackParCovCharmHad.getSigmaTglY(), trackParCovCharmHad.getSigmaTglZ(),
                          trackParCovCharmHad.getSigmaTglSnp(), trackParCovCharmHad.getSigmaTgl2(),
                          trackParCovCharmHad.getSigma1PtY(), trackParCovCharmHad.getSigma1PtZ(), trackParCovCharmHad.getSigma1PtSnp(),
                          trackParCovCharmHad.getSigma1PtTgl(), trackParCovCharmHad.getSigma1Pt2());
          hfCandPidProng0(charmHadDauTracks[0].tpcNSigmaPi(), charmHadDauTracks[0].tofNSigmaPi(), charmHadDauTracks[0].tpcNSigmaKa(), charmHadDauTracks[0].tofNSigmaKa(), charmHadDauTracks[0].hasTOF(), charmHadDauTracks[0].hasTPC());
          hfCandPidProng1(charmHadDauTracks[1].tpcNSigmaPi(), charmHadDauTracks[1].tofNSigmaPi(), charmHadDauTracks[1].tpcNSigmaKa(), charmHadDauTracks[1].tofNSigmaKa(), charmHadDauTracks[1].hasTOF(), charmHadDauTracks[1].hasTPC());
          hfCandPidProng2(charmHadDauTracks[2].tpcNSigmaPi(), charmHadDauTracks[2].tofNSigmaPi(), charmHadDauTracks[2].tpcNSigmaKa(), charmHadDauTracks[2].tofNSigmaKa(), charmHadDauTracks[2].hasTOF(), charmHadDauTracks[2].hasTPC());
          if constexpr (withMl) {
            if constexpr (decChannel == DecayChannel::B0ToDminusPi) {
              hfCand3ProngMl(candC.mlProbDplusToPiKPi()[0], candC.mlProbDplusToPiKPi()[1], candC.mlProbDplusToPiKPi()[2], -1., -1., -1.);
            } else {
              std::array<float, 6> mlScores = {-1.f, -1.f, -1.f, -1.f, -1.f, -1.f};
              if (candC.mlProbDsToKKPi().size() == 3) {
                std::copy(candC.mlProbDsToKKPi().begin(), candC.mlProbDsToKKPi().end(), mlScores.begin());
              }
              if (candC.mlProbDsToPiKK().size() == 3) {
                std::copy(candC.mlProbDsToPiKK().begin(), candC.mlProbDsToPiKK().end(), mlScores.begin() + 3);
              }
              hfCand3ProngMl(mlScores[0], mlScores[1], mlScores[2], mlScores[3], mlScores[4], mlScores[5]);
            }
          }
        } else if constexpr (decChannel == DecayChannel::BplusToD0barPi) { // D0(bar) → K± π∓
          hfCand2Prong(charmHadDauTracks[0].globalIndex(), charmHadDauTracks[1].globalIndex(),
                       indexHfReducedCollision,
                       trackParCovCharmHad.getX(), trackParCovCharmHad.getAlpha(),
                       trackParCovCharmHad.getY(), trackParCovCharmHad.getZ(), trackParCovCharmHad.getSnp(),
                       trackParCovCharmHad.getTgl(), trackParCovCharmHad.getQ2Pt(),
                       candC.xSecondaryVertex(), candC.ySecondaryVertex(), candC.zSecondaryVertex(), invMassC0, invMassC1,
                       ptDauMin, etaDauMin, nItsClsDauMin, nTpcCrossRowsDauMin, chi2TpcDauMax);
          hfCand2ProngCov(trackParCovCharmHad.getSigmaY2(), trackParCovCharmHad.getSigmaZY(), trackParCovCharmHad.getSigmaZ2(),
                          trackParCovCharmHad.getSigmaSnpY(), trackParCovCharmHad.getSigmaSnpZ(),
                          trackParCovCharmHad.getSigmaSnp2(), trackParCovCharmHad.getSigmaTglY(), trackParCovCharmHad.getSigmaTglZ(),
                          trackParCovCharmHad.getSigmaTglSnp(), trackParCovCharmHad.getSigmaTgl2(),
                          trackParCovCharmHad.getSigma1PtY(), trackParCovCharmHad.getSigma1PtZ(), trackParCovCharmHad.getSigma1PtSnp(),
                          trackParCovCharmHad.getSigma1PtTgl(), trackParCovCharmHad.getSigma1Pt2());
          hfCandPidProng0(candC.nSigTpcPi0(), candC.nSigTofPi0(), candC.nSigTpcKa0(), candC.nSigTofKa0(), charmHadDauTracks[0].hasTOF(), charmHadDauTracks[0].hasTPC());
          hfCandPidProng1(candC.nSigTpcPi1(), candC.nSigTofPi1(), candC.nSigTpcKa1(), candC.nSigTofKa1(), charmHadDauTracks[1].hasTOF(), charmHadDauTracks[1].hasTPC());
          if constexpr (withMl) {
            std::array<float, 6> mlScores = {-1.f, -1.f, -1.f, -1.f, -1.f, -1.f};
            if (candC.mlProbD0().size() == 3) {
              std::copy(candC.mlProbD0().begin(), candC.mlProbD0().end(), mlScores.begin());
            }
            if (candC.mlProbD0bar().size() == 3) {
              std::copy(candC.mlProbD0bar().begin(), candC.mlProbD0bar().end(), mlScores.begin() + 3);
            }
            hfCand2ProngMl(mlScores[0], mlScores[1], mlScores[2], mlScores[3], mlScores[4], mlScores[5]);
          }
        }
        fillHfReducedCollision = true;
      }
    } // candsC loop

    registry.fill(HIST("hEvents"), 1 + Event::Processed);
    if (!fillHfReducedCollision) {
      registry.fill(HIST("hEvents"), 1 + Event::NoCharmHadPiSelected);
      return;
    }
    registry.fill(HIST("hEvents"), 1 + Event::CharmHadPiSelected);
    float centrality = -1.f;
    uint16_t hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
    // fill collision table if it contains a DPi pair a minima
    hfReducedCollision(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), hfRejMap, bz);
    hfReducedCollExtra(collision.covXX(), collision.covXY(), collision.covYY(),
                       collision.covXZ(), collision.covYZ(), collision.covZZ());
  }

  template <uint8_t decayChannel>
  void runMcGen(aod::McParticles const& particlesMc)
  {
    // Match generated particles.
    for (const auto& particle : particlesMc) {
      int8_t sign{0};
      int8_t flag{0};
      if constexpr (decayChannel == DecayChannel::B0ToDminusPi) {
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

        std::array<float, 2> ptProngs;
        std::array<float, 2> yProngs;
        std::array<float, 2> etaProngs;
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }
        rowHfB0McGenReduced(flag, ptParticle, yParticle, etaParticle,
                            ptProngs[0], yProngs[0], etaProngs[0],
                            ptProngs[1], yProngs[1], etaProngs[1]);
      } else if constexpr (decayChannel == DecayChannel::BsToDsminusPi) {
        // Bs → Ds- π+
        if (RecoDecay::isMatchedMCGen<true>(particlesMc, particle, Pdg::kBS, std::array{-static_cast<int>(Pdg::kDS), +kPiPlus}, true)) {
          // Match Ds- -> π- K+ π-
          auto candCMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(particlesMc, candCMC, -static_cast<int>(Pdg::kDS), std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2)) {
            std::vector<int> arrDaughDsIndex;
            std::array<int, 2> arrPDGDaughDs;
            RecoDecay::getDaughters(candCMC, &arrDaughDsIndex, std::array{0}, 1);
            if (arrDaughDsIndex.size() == 2) {
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
        if (checkDecayTypeMc) {
          // B0 → Ds- π+
          if (!flag) {
            if (RecoDecay::isMatchedMCGen<true>(particlesMc, particle, Pdg::kB0, std::array{-static_cast<int>(Pdg::kDS), +kPiPlus}, true)) {
              // Match Ds- -> π- K+ π-
              auto candCMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
              if (RecoDecay::isMatchedMCGen(particlesMc, candCMC, -static_cast<int>(Pdg::kDS), std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2)) {
                std::vector<int> arrDaughDsIndex;
                std::array<int, 2> arrPDGDaughDs;
                RecoDecay::getDaughters(candCMC, &arrDaughDsIndex, std::array{0}, 1);
                if (arrDaughDsIndex.size() == 2) {
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

        std::array<float, 2> ptProngs;
        std::array<float, 2> yProngs;
        std::array<float, 2> etaProngs;
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }
        rowHfBsMcGenReduced(flag, ptParticle, yParticle, etaParticle,
                            ptProngs[0], yProngs[0], etaProngs[0],
                            ptProngs[1], yProngs[1], etaProngs[1]);
      } else if constexpr (decayChannel == DecayChannel::BplusToD0barPi) {
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

        std::array<float, 2> ptProngs;
        std::array<float, 2> yProngs;
        std::array<float, 2> etaProngs;
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }
        rowHfBpMcGenReduced(flag, ptParticle, yParticle, etaParticle,
                            ptProngs[0], yProngs[0], etaProngs[0],
                            ptProngs[1], yProngs[1], etaProngs[1]);
      }
    } // gen
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // PROCESS FUNCTIONS FOR DATA

  void processDplusPiData(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                          CandsDplusFiltered const& candsC,
                          aod::TrackAssoc const& trackIndices,
                          TracksPidWithSel const& tracks,
                          aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigB0(selectionFlagDplus.value, invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(candsDplusPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::B0ToDminusPi>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDplusPiData, "Process DplusPi without MC info and without ML info", true);

  void processDplusPiDataWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                CandsDplusFilteredWithMl const& candsC,
                                aod::TrackAssoc const& trackIndices,
                                TracksPidWithSel const& tracks,
                                aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigB0(selectionFlagDplus.value, invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(candsDplusPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::B0ToDminusPi>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDplusPiDataWithMl, "Process DplusPi without MC info and with ML info", false);

  void processDsPiData(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                       CandsDsFiltered const& candsC,
                       aod::TrackAssoc const& trackIndices,
                       TracksPidWithSel const& tracks,
                       aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for Bs workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBs(selectionFlagDs.value, invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(candsDplusPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::BsToDsminusPi>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDsPiData, "Process DsPi without MC info and without ML info", true);

  void processDsPiDataWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                             CandsDsFilteredWithMl const& candsC,
                             aod::TrackAssoc const& trackIndices,
                             TracksPidWithSel const& tracks,
                             aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for Bs workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBs(selectionFlagDs.value, invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(candsDplusPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::BsToDsminusPi>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDsPiDataWithMl, "Process DsPi without MC info and with ML info", false);

  void processD0PiData(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                       CandsD0Filtered const& candsC,
                       aod::TrackAssoc const& trackIndices,
                       TracksPidWithSel const& tracks,
                       aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBplus(selectionFlagD0.value, selectionFlagD0bar.value, invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(candsD0PerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::BplusToD0barPi>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processD0PiData, "Process D0Pi without MC info and without ML info", false);

  void processD0PiDataWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                             CandsD0FilteredWithMl const& candsC,
                             aod::TrackAssoc const& trackIndices,
                             TracksPidWithSel const& tracks,
                             aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBplus(selectionFlagD0.value, selectionFlagD0bar.value, invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(candsD0PerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::BplusToD0barPi>(collision, candsCThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processD0PiDataWithMl, "Process D0Pi without MC info and with ML info", false);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // PROCESS FUNCTIONS FOR MC

  void processDplusPiMc(CollisionsWMcLabels const& collisions,
                        CandsDplusFiltered const& candsC,
                        aod::TrackAssoc const& trackIndices,
                        TracksPidWithSelAndMc const& tracks,
                        aod::McParticles const& particlesMc,
                        aod::BCsWithTimestamps const& bcs,
                        McCollisions const&)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigB0(selectionFlagDplus.value, invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(candsDplusPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(colPerMcCollision, collision.mcCollisionId());
      int64_t indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, false, DecayChannel::B0ToDminusPi>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    runMcGen<DecayChannel::B0ToDminusPi>(particlesMc);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDplusPiMc, "Process DplusPi with MC info and without ML info", false);

  void processDplusPiMcWithMl(CollisionsWMcLabels const& collisions,
                              CandsDplusFilteredWithMl const& candsC,
                              aod::TrackAssoc const& trackIndices,
                              TracksPidWithSelAndMc const& tracks,
                              aod::McParticles const& particlesMc,
                              aod::BCsWithTimestamps const& bcs,
                              McCollisions const&)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigB0(selectionFlagDplus.value, invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(candsDplusPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(colPerMcCollision, collision.mcCollisionId());
      int64_t indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, true, DecayChannel::B0ToDminusPi>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    runMcGen<DecayChannel::B0ToDminusPi>(particlesMc);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDplusPiMcWithMl, "Process DplusPi with MC info and with ML info", false);

  void processDsPiMc(CollisionsWMcLabels const& collisions,
                     CandsDsFiltered const& candsC,
                     aod::TrackAssoc const& trackIndices,
                     TracksPidWithSelAndMc const& tracks,
                     aod::McParticles const& particlesMc,
                     aod::BCsWithTimestamps const& bcs,
                     McCollisions const&)
  {
    // store configurables needed for Bs workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBs(selectionFlagDs.value, invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(candsDplusPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(colPerMcCollision, collision.mcCollisionId());
      int64_t indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, false, DecayChannel::BsToDsminusPi>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    runMcGen<DecayChannel::BsToDsminusPi>(particlesMc);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDsPiMc, "Process DsPi with MC info and without ML info", false);

  void processDsPiMcWithMl(CollisionsWMcLabels const& collisions,
                           CandsDsFilteredWithMl const& candsC,
                           aod::TrackAssoc const& trackIndices,
                           TracksPidWithSelAndMc const& tracks,
                           aod::McParticles const& particlesMc,
                           aod::BCsWithTimestamps const& bcs,
                           McCollisions const&)
  {
    // store configurables needed for Bs workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBs(selectionFlagDs.value, invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(candsDplusPerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(colPerMcCollision, collision.mcCollisionId());
      int64_t indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, true, DecayChannel::BsToDsminusPi>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    runMcGen<DecayChannel::BsToDsminusPi>(particlesMc);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processDsPiMcWithMl, "Process DsPi with MC info and with ML info", false);

  void processD0PiMc(CollisionsWMcLabels const& collisions,
                     CandsD0Filtered const& candsC,
                     aod::TrackAssoc const& trackIndices,
                     TracksPidWithSelAndMc const& tracks,
                     aod::McParticles const& particlesMc,
                     aod::BCsWithTimestamps const& bcs,
                     McCollisions const&)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBplus(selectionFlagD0.value, selectionFlagD0bar.value, invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(candsD0PerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(colPerMcCollision, collision.mcCollisionId());
      int64_t indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, false, DecayChannel::BplusToD0barPi>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    runMcGen<DecayChannel::BplusToD0barPi>(particlesMc);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processD0PiMc, "Process D0Pi with MC info and without ML info", false);

  void processD0PiMcWithMl(CollisionsWMcLabels const& collisions,
                           CandsD0FilteredWithMl const& candsC,
                           aod::TrackAssoc const& trackIndices,
                           TracksPidWithSelAndMc const& tracks,
                           aod::McParticles const& particlesMc,
                           aod::BCsWithTimestamps const& bcs,
                           McCollisions const&)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBplus(selectionFlagD0.value, selectionFlagD0bar.value, invMassWindowCharmHadPi.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsCThisColl = candsC.sliceBy(candsD0PerCollisionWithMl, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(colPerMcCollision, collision.mcCollisionId());
      int64_t indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, true, DecayChannel::BplusToD0barPi>(collision, candsCThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    runMcGen<DecayChannel::BplusToD0barPi>(particlesMc);
  }
  PROCESS_SWITCH(HfDataCreatorCharmHadPiReduced, processD0PiMcWithMl, "Process D0Pi with MC info and with ML info", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDataCreatorCharmHadPiReduced>(cfgc)};
}
