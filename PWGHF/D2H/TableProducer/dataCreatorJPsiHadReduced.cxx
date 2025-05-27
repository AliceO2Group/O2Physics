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

/// \file dataCreatorJPsiHadReduced.cxx
/// \brief Creation of J/Psi-LF hadron pairs for Beauty hadron analyses
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università degli Studi and INFN Torino
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

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
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Qvectors.h"

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

enum Event : uint8_t { // TODO: check if needed
  Processed = 0,
  NoCharmHadPiSelected,
  CharmHadPiSelected,
  kNEvent
};

enum DecayChannel : uint8_t {
  B0ToJPsiK0Star = 0,
  BplusToJPsiK,
  BsToJPsiPhi
};

enum WrongCollisionType : uint8_t {
  None = 0,
  WrongAssociation,
  SplitCollision,
};

/// Creation of JPsi-Had pairs for Beauty hadrons
struct HfDataCreatorJPsiHadReduced {
  // Produces AOD tables to store track information
  // collision related tables
  Produces<aod::HfRedCollisions> hfReducedCollision;
  Produces<aod::HfRedCollCents> hfReducedCollCentrality;
  Produces<aod::HfRedQvectors> hfReducedQvector;
  Produces<aod::HfRedCollExtras> hfReducedCollExtra;
  Produces<aod::HfOrigColCounts> hfCollisionCounter;
  // J/Psi related tables
  Produces<aod::HfRedJPsis> hfJPsi;
  Produces<aod::HfRedJPsiCov> hfRedJPsiCov;
  // Ka bachelor related tables
  Produces<aod::HfRedBachProng0Bases> hfTrackLfDau0;
  Produces<aod::HfRedBachProng0Cov> hfTrackCovLfDau0;
  Produces<aod::HfRedBachProng1Bases> hfTrackLfDau1;
  Produces<aod::HfRedBachProng1Cov> hfTrackCovLfDau1;
  // MC related tables
  Produces<aod::HfMcRecRedJPsiKs> rowHfJPsiKMcRecReduced;
  Produces<aod::HfMcRecRedJPsiPhis> rowHfJPsiPhiMcRecReduced;
  Produces<aod::HfMcGenRedBps> rowHfBpMcGenReduced;
  Produces<aod::HfMcGenRedBss> rowHfBsMcGenReduced;

  Produces<aod::HfCfgBpToJPsi> rowCandidateConfigBplus;
  Produces<aod::HfCfgBsToJPsi> rowCandidateConfigBs;

  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any B0 is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};

  Configurable<bool> runJPsiToee{"runJPsiToee", false, "Run analysis for J/Psi to ee (debug)"};
  Configurable<double> ptJpsiMin{"ptJpsiMin", 0., "Lower bound of J/Psi pT"};
  Configurable<double> ptJpsiMax{"ptJpsiMax", 50., "Upper bound of J/Psi pT"};
  Configurable<bool> useTrackIsGlobalTrackWoDCA{"useTrackIsGlobalTrackWoDCA", true, "check isGlobalTrackWoDCA status for the bachelor tracks"};
  Configurable<double> ptTrackMin{"ptTrackMin", 0.5, "minimum bachelor track pT threshold (GeV/c)"};
  Configurable<double> absEtaTrackMax{"absEtaTrackMax", 0.8, "maximum bachelor track absolute eta threshold"};
  Configurable<std::vector<double>> binsPtTrack{"binsPtTrack", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for bachelor track DCA XY pT-dependent cut"};
  Configurable<LabeledArray<double>> cutsTrackDCA{"cutsTrackDCA", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for bachelor track"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_jpsi_to_mu_mu::vecBinsPt}, "J/Psi pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_jpsi_to_mu_mu::Cuts[0], hf_cuts_jpsi_to_mu_mu::NBinsPt, hf_cuts_jpsi_to_mu_mu::NCutVars, hf_cuts_jpsi_to_mu_mu::labelsPt, hf_cuts_jpsi_to_mu_mu::labelsCutVar}, "J/Psi candidate selection per pT bin"};

  Configurable<double> invMassWindowJPsiHad{"invMassWindowJPsiHad", 0.3, "invariant-mass window for JPsi-Had pair preselections (GeV/c2)"};
  Configurable<double> deltaMPhiMax{"deltaMPhiMax", 0.02, "invariant-mass window for phi preselections (GeV/c2) (only for Bs->J/PsiPhi)"};
  // magnetic field setting from CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  HfHelper hfHelper;

  // CCDB service
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // O2DatabasePDG service
  Service<o2::framework::O2DatabasePDG> pdg;

  using TracksPid = soa::Join<aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;
  using TracksPidWithSel = soa::Join<aod::TracksWCovDcaExtra, TracksPid, aod::TrackSelection>;
  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;
  using TracksPidWithSelAndMc = soa::Join<TracksPidWithSel, aod::McTrackLabels>;
  using CollisionsWCMcLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

  Preslice<aod::HfCand2ProngWPid> candsJPsiPerCollision = aod::track_association::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  PresliceUnsorted<CollisionsWCMcLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber;
  double bz{0.};
  double invMass2JPsiHadMin, invMass2JPsiHadMax;
  bool isHfCandBhadConfigFilled = false;

  o2::hf_evsel::HfEventSelection hfEvSel;
  o2::vertexing::DCAFitterN<2> df2;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    std::array<int, 4> doProcess = {doprocessJPsiKData, doprocessJPsiKMc, doprocessJPsiPhiData, doprocessJPsiPhiMc}; // TODO: check completeness
    if (std::accumulate(doProcess.begin(), doProcess.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function can be enabled at a time, please fix your configuration!");
    }

    // Set up the histogram registry
    constexpr int kNBinsSelections = 2 + aod::SelectionStep::RecoTopol;
    std::string labels[kNBinsSelections];
    labels[0] = "No selection";
    labels[1 + aod::SelectionStep::RecoSkims] = "Skims selection";
    labels[1 + aod::SelectionStep::RecoTopol] = "Skims & Topological selections";
    static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections + 0.5, ""};
    registry.add("hSelectionsJPsi", "J/Psi selection;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSelections, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
      registry.get<TH2>(HIST("hSelectionsJPsi"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    constexpr int kNBinsEvents = kNEvent;
    std::string labelsEvents[kNBinsEvents];
    labelsEvents[Event::Processed] = "processed";
    labelsEvents[Event::NoCharmHadPiSelected] = "without CharmHad-Pi pairs";
    labelsEvents[Event::CharmHadPiSelected] = "with CharmHad-Pi pairs";
    static const AxisSpec axisEvents = {kNBinsEvents, 0.5, kNBinsEvents + 0.5, ""};
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {axisEvents});
    for (int iBin = 0; iBin < kNBinsEvents; iBin++) {
      registry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(iBin + 1, labelsEvents[iBin].data());
    }

    registry.add("hMassJPsi", "J/Psi mass;#it{M}_{#mu#mu} (GeV/#it{c}^{2});Counts", {HistType::kTH1F, {{600, 2.8, 3.4, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtJPsi", "J/Psi #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", {HistType::kTH1F, {{(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCpaJPsi", "J/Psi cos#theta_{p};J/Psi cos#theta_{p};Counts", {HistType::kTH1F, {{200, -1., 1, "J/Psi cos#theta_{p}"}}});
    std::shared_ptr<TH1> hFitCandidatesJPsi = registry.add<TH1>("hFitCandidatesJPsi", "JPsi candidate counter", {HistType::kTH1D, {axisCands}});
    setLabelHistoCands(hFitCandidatesJPsi);
    if (doprocessJPsiKData || doprocessJPsiKMc) {
      registry.add("hPtKaon", "Kaon #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", {HistType::kTH1F, {{100, 0., 10.}}});
      registry.add("hMassJPsiKaon", "J/Psi Kaon mass;#it{M}_{J/#PsiK} (GeV/#it{c}^{2});Counts", {HistType::kTH1F, {{800, 4.9, 5.7}}});
    } else if (doprocessJPsiPhiData || doprocessJPsiPhiMc) {
      registry.add("hPtPhi", "Phi #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", {HistType::kTH1F, {{100, 0., 10.}}});
      registry.add("hMassPhi", "Phi mass;#it{M}_{KK} (GeV/#it{c}^{2});Counts", {HistType::kTH1F, {{200, 0.9, 1.2}}});
      registry.add("hMassJPsiPhi", "J/Psi Phi mass;#it{M}_{J/#Psi#phi} (GeV/#it{c}^{2});Counts", {HistType::kTH1F, {{800, 4.9, 5.7}}});
      std::shared_ptr<TH1> hFitCandidatesPhi = registry.add<TH1>("hFitCandidatesPhi", "Phi candidate counter", {HistType::kTH1D, {axisCands}});
      setLabelHistoCands(hFitCandidatesPhi);
    }

    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);
    df2.setMatCorrType(noMatCorr);

    // Configure CCDB access
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    runNumber = 0;

    if (doprocessJPsiKData || doprocessJPsiKMc) {
      invMass2JPsiHadMin = (MassBPlus - invMassWindowJPsiHad) * (MassBPlus - invMassWindowJPsiHad); // TODO: change mass bplus for other species
      invMass2JPsiHadMax = (MassBPlus + invMassWindowJPsiHad) * (MassBPlus + invMassWindowJPsiHad); // TODO: change mass bplus for other species
    } else if (doprocessJPsiPhiData || doprocessJPsiPhiMc) {
      invMass2JPsiHadMin = (MassBS - invMassWindowJPsiHad) * (MassBS - invMassWindowJPsiHad); // TODO: change mass bplus for other species
      invMass2JPsiHadMax = (MassBS + invMassWindowJPsiHad) * (MassBS + invMassWindowJPsiHad); // TODO: change mass bplus for other species
    }
  }

  /// Topological cuts
  /// \param candidate is candidate
  /// \param trackPos is the positive track
  /// \param trackNeg is the negative track
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2>
  bool selectionTopol(const T1& candidate, const T2& trackPos, const T2& trackNeg)
  {
    auto candpT = candidate.pt();
    auto pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < ptJpsiMin || candpT >= ptJpsiMax) {
      return false;
    }

    // cut on μ+ μ− (e+e−) invariant mass
    if (runJPsiToee && std::abs(hfHelper.invMassJpsiToEE(candidate) - o2::constants::physics::MassJPsi) > cuts->get(pTBin, "m")) {
      return false;
    } else if (!runJPsiToee && std::abs(hfHelper.invMassJpsiToMuMu(candidate) - o2::constants::physics::MassJPsi) > cuts->get(pTBin, "m")) {
      return false;
    }

    // cut on daughter pT (same cut used for both channels)
    if (trackNeg.pt() < cuts->get(pTBin, "pT mu") || trackPos.pt() < cuts->get(pTBin, "pT mu")) {
      return false;
    }

    // decay length
    if (candidate.decayLength() < cuts->get(pTBin, "decay length")) {
      return false;
    }

    // decay length in XY plane
    if (candidate.decayLengthXY() < cuts->get(pTBin, "decay length xy")) {
      return false;
    }

    // cosine of pointing angle
    if (candidate.cpa() < cuts->get(pTBin, "cpa")) {
      return false;
    }

    // cosine of pointing angle XY
    if (candidate.cpaXY() < cuts->get(pTBin, "cpa xy")) {
      return false;
    }

    // product of daughter impact parameters
    if (candidate.impactParameterProduct() > cuts->get(pTBin, "d0xd0")) {
      return false;
    }

    return true;
  }

  /// Single-track cuts for kaons on dcaXY
  /// \param trackPar is the track parametrisation
  /// \param dca is the 2-D array with track DCAs
  /// \return true if track passes all cuts
  template <typename T1, typename T2>
  bool isSelectedTrackDCA(const T1& trackPar, const T2& dca)
  {
    auto pTBinTrack = findBin(binsPtTrack, trackPar.getPt());
    if (pTBinTrack == -1) {
      return false;
    }

    if (std::abs(dca[0]) < cutsTrackDCA->get(pTBinTrack, "min_dcaxytoprimary")) {
      return false; // minimum DCAxy
    }
    if (std::abs(dca[0]) > cutsTrackDCA->get(pTBinTrack, "max_dcaxytoprimary")) {
      return false; // maximum DCAxy
    }
    return true;
  }

  /// Kaon selection (J/Psi K+ <-- B+)
  /// \param track is the considered track
  /// \param trackParCov is the track parametrisation
  /// \param dca is the 2-D array with track DCAs
  /// \param jPsiDautracks J/Psi daughter tracks
  /// \return true if track passes all cuts
  template <typename T1, typename T2, typename T3>
  bool isTrackSelected(const T1& track, const T2& trackParCov, const T3& dca, const std::vector<T1>& jPsiDautracks)
  {
    // check isGlobalTrackWoDCA status for kaons if wanted
    if (useTrackIsGlobalTrackWoDCA && !track.isGlobalTrackWoDCA()) {
      return false;
    }
    // minimum pT, eta, and DCA selection
    if (trackParCov.getPt() < ptTrackMin || std::abs(trackParCov.getEta()) > absEtaTrackMax || !isSelectedTrackDCA(trackParCov, dca)) {
      return false;
    }
    // reject kaons that are J/Psi daughters
    for (const auto& trackJPsi : jPsiDautracks) {
      if (track.globalIndex() == trackJPsi.globalIndex()) {
        return false;
      }
    }

    return true;
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
  /// \param vecDaughtersB is the vector with all daughter tracks (JPsi daughters in first position)
  /// \param indexHfCandJPsi is the index of the JPsi candidate
  /// \param selectedTracksBach is the map with the indices of selected bachelor pion tracks
  template <uint8_t decChannel, typename CColl, typename PParticles, typename TTrack>
  void fillMcRecoInfo(const CColl& collision,
                      const PParticles& particlesMc,
                      const std::vector<TTrack>& vecDaughtersB,
                      int& indexHfCandJPsi,
                      std::array<std::map<int64_t, int64_t>, 2> selectedTracksBach,
                      const int64_t indexCollisionMaxNumContrib)
  {

    // we check the MC matching to be stored
    int8_t sign{0};
    int8_t flag{0};
    int8_t flagWrongCollision{WrongCollisionType::None};
    int8_t debug{0};
    float motherPt{-1.f};

    if constexpr (decChannel == DecayChannel::BplusToJPsiK) {
      // B+ → J/Psi K+ → (µ+µ-) K+
      int indexRec = -1;
      if (!runJPsiToee) {
        indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, Pdg::kBPlus, std::array{-kMuonMinus, +kMuonMinus, +kKPlus}, true, &sign, 3);
      } else {
        indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, Pdg::kBPlus, std::array{-kElectron, +kElectron, +kKPlus}, true, &sign, 3);
      }
      if (indexRec > -1) {
        // J/Psi → µ+µ-
        if (!runJPsiToee) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1]}, Pdg::kJPsi, std::array{-kMuonMinus, +kMuonMinus}, true, &sign, 1);
        } else {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1]}, Pdg::kJPsi, std::array{-kElectron, +kElectron}, true, &sign, 1);
        }
        if (indexRec > -1) {
          flag = sign * BIT(static_cast<uint8_t>(hf_cand_bplus::DecayTypeBToJPsiMc::BplusToJPsiKToMuMuK));
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
      rowHfJPsiKMcRecReduced(indexHfCandJPsi, selectedTracksBach[0][vecDaughtersB.back().globalIndex()], flag, flagWrongCollision, debug, motherPt);
    } else if constexpr (decChannel == DecayChannel::BsToJPsiPhi) {
      // Bs → J/Psi phi → (µ+µ-) (K+K-)
      int indexRec = -1;
      if (!runJPsiToee) {
        indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kBS, std::array{-kMuonMinus, +kMuonMinus, +kKPlus, -kKPlus}, true, &sign, 4);
      } else {
        indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kBS, std::array{-kElectron, +kElectron, +kKPlus, -kKPlus}, true, &sign, 4);
      }
      if (indexRec > -1) {
        // J/Psi → µ+µ-
        if (!runJPsiToee) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1]}, Pdg::kJPsi, std::array{-kMuonMinus, +kMuonMinus}, true, &sign, 1);
        } else {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1]}, Pdg::kJPsi, std::array{-kElectron, +kElectron}, true, &sign, 1);
        }
        if (indexRec > -1) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kPhi, std::array{-kKPlus, +kKPlus}, true, &sign, 1);
          if (indexRec > -1) {
            flag = sign * BIT(static_cast<uint8_t>(hf_cand_bs::DecayTypeBToJPsiMc::BsToJPsiPhiToMuMuKK));
          } else {
            debug = 1;
            LOGF(debug, "Bs decays in the expected final state but the condition on the phi intermediate state is not fulfilled");
          }
        } else {
          debug = 1;
          LOGF(debug, "Bs decays in the expected final state but the condition on the J/Psi intermediate state is not fulfilled");
        }

        auto indexMother = RecoDecay::getMother(particlesMc, vecDaughtersB.back().template mcParticle_as<PParticles>(), Pdg::kBS, true);
        if (indexMother >= 0) {
          auto particleMother = particlesMc.rawIteratorAt(indexMother);
          motherPt = particleMother.pt();
          checkWrongCollision(particleMother, collision, indexCollisionMaxNumContrib, flagWrongCollision);
        }
      }
      rowHfJPsiPhiMcRecReduced(indexHfCandJPsi, selectedTracksBach[0][vecDaughtersB.back().globalIndex()], selectedTracksBach[1][vecDaughtersB.back().globalIndex()], flag, flagWrongCollision, debug, motherPt);
    }
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

  template <uint8_t decChannel>
  void runMcGen(aod::McParticles const& particlesMc)
  {
    // Match generated particles.
    for (const auto& particle : particlesMc) {
      int8_t sign{0};
      int8_t flag{0};
      if constexpr (decChannel == DecayChannel::BplusToJPsiK) {
        // B+ → J/Psi K+ → (µ+µ-) K+
        if (RecoDecay::isMatchedMCGen<false>(particlesMc, particle, Pdg::kBPlus, std::array{static_cast<int>(Pdg::kJPsi), +kKPlus}, true)) {
          // Match J/Psi -> µ+µ-
          auto candJPsiMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
          // Printf("Checking J/Psi -> µ+µ-");
          if (!runJPsiToee) {
            if (RecoDecay::isMatchedMCGen(particlesMc, candJPsiMC, static_cast<int>(Pdg::kJPsi), std::array{-kMuonMinus, +kMuonMinus}, true, &sign, 2)) {
              flag = sign * BIT(hf_cand_bplus::DecayType::BplusToJPsiK);
            }
          } else { // debug
            // Printf("Checking J/Psi -> e+e-");
            if (RecoDecay::isMatchedMCGen(particlesMc, candJPsiMC, static_cast<int>(Pdg::kJPsi), std::array{-kElectron, +kElectron}, true, &sign, 2)) {
              flag = sign * BIT(hf_cand_bplus::DecayType::BplusToJPsiK);
            }
          }
        }

        // save information for B+ task
        if (!TESTBIT(std::abs(flag), hf_cand_bplus::DecayType::BplusToJPsiK)) {
          continue;
        }

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), MassBPlus);
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
      } else if constexpr (decChannel == DecayChannel::BsToJPsiPhi) {
        // Bs → J/Psi phi → (µ+µ-) (K+K-)
        if (RecoDecay::isMatchedMCGen<true>(particlesMc, particle, Pdg::kBS, std::array{static_cast<int>(Pdg::kJPsi), static_cast<int>(Pdg::kPhi)}, true)) {
          // Match J/Psi -> µ+µ- and phi -> K+K-
          auto candJPsiMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
          auto candPhiMC = particlesMc.rawIteratorAt(particle.daughtersIds().back());
          // Printf("Checking J/Psi -> µ+µ- and phi -> K+K-");
          if (!runJPsiToee) {
            if (RecoDecay::isMatchedMCGen(particlesMc, candJPsiMC, static_cast<int>(Pdg::kJPsi), std::array{-kMuonMinus, +kMuonMinus}, true, &sign, 2) &&
                RecoDecay::isMatchedMCGen(particlesMc, candPhiMC, static_cast<int>(Pdg::kPhi), std::array{-kKPlus, +kKPlus}, true, &sign, 2)) {
              flag = sign * BIT(hf_cand_bs::DecayType::BsToJPsiPhi);
            }
          } else { // debug
            // Printf("Checking J/Psi -> e+e- and phi -> K+K-");
            if (RecoDecay::isMatchedMCGen(particlesMc, candJPsiMC, static_cast<int>(Pdg::kJPsi), std::array{-kElectron, +kElectron}, true, &sign, 2) &&
                RecoDecay::isMatchedMCGen(particlesMc, candPhiMC, static_cast<int>(Pdg::kPhi), std::array{-kKPlus, +kKPlus}, true, &sign, 2)) {
              flag = sign * BIT(hf_cand_bs::DecayType::BsToJPsiPhi);
            }
          }
        }

        // save information for Bs task
        if (!TESTBIT(std::abs(flag), hf_cand_bs::DecayType::BsToJPsiPhi)) {
          continue;
        }

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), MassBPlus);
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
      }
    } // gen
  }

  // JPsi candidate selection
  template <bool doMc, uint8_t decChannel, typename Coll, typename JPsiCands, typename TTracks, typename PParticles>
  void runDataCreation(Coll const& collision,
                       JPsiCands const& candsJPsi,
                       aod::TrackAssoc const& trackIndices,
                       TTracks const&,
                       PParticles const& particlesMc,
                       uint64_t const& indexCollisionMaxNumContrib,
                       aod::BCsWithTimestamps const&)
  {

    // helpers for ReducedTables filling
    int indexHfReducedCollision = hfReducedCollision.lastIndex() + 1;
    // std::map where the key is the track.globalIndex() and
    // the value is the track index in the table of the selected tracks
    std::map<int64_t, int64_t> selectedTracksBach;
    std::map<int64_t, int64_t> selectedTracksBach2; // for the second daughter (for B0 and Bs)

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
    df2.setBz(bz);

    auto thisCollId = collision.globalIndex();
    // looping over 2-prong candidates
    for (const auto& candidate : candsJPsi) {

      // Apply the selections on the J/Psi candidates
      registry.fill(HIST("hSelectionsJPsi"), 1, candidate.pt());

      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::JpsiToMuMu)) {
        continue;
      }
      registry.fill(HIST("hSelectionsJPsi"), 2 + aod::SelectionStep::RecoSkims, candidate.pt());

      auto trackPos = candidate.template prong0_as<TTracks>(); // positive daughter
      auto trackNeg = candidate.template prong1_as<TTracks>(); // negative daughter

      auto trackPosParCov = getTrackParCov(trackPos);
      auto trackNegParCov = getTrackParCov(trackNeg);

      // topological selection
      if (!selectionTopol(candidate, trackPos, trackNeg)) {
        continue;
      }
      registry.fill(HIST("hSelectionsJPsi"), 2 + aod::SelectionStep::RecoTopol, candidate.pt());

      int indexHfCandJPsi = hfJPsi.lastIndex() + 1;
      float invMassJPsi{0.f};
      if (runJPsiToee) {
        invMassJPsi = hfHelper.invMassJpsiToEE(candidate);
      } else {
        invMassJPsi = hfHelper.invMassJpsiToMuMu(candidate);
      }
      registry.fill(HIST("hMassJPsi"), invMassJPsi);
      registry.fill(HIST("hPtJPsi"), candidate.pt());
      registry.fill(HIST("hCpaJPsi"), candidate.cpa());

      bool fillHfCandJPsi = false;

      std::vector<typename TTracks::iterator> jPsiDauTracks{trackPos, trackNeg};

      std::array<float, 3> pVec0 = jPsiDauTracks[0].pVector();
      std::array<float, 3> pVec1 = jPsiDauTracks[1].pVector();

      auto dca0 = o2::dataformats::DCA(jPsiDauTracks[0].dcaXY(), jPsiDauTracks[0].dcaZ(), jPsiDauTracks[0].cYY(), jPsiDauTracks[0].cZY(), jPsiDauTracks[0].cZZ());
      auto dca1 = o2::dataformats::DCA(jPsiDauTracks[1].dcaXY(), jPsiDauTracks[1].dcaZ(), jPsiDauTracks[1].cYY(), jPsiDauTracks[1].cZY(), jPsiDauTracks[1].cZZ());

      // repropagate tracks to this collision if needed
      if (jPsiDauTracks[0].collisionId() != thisCollId) {
        trackPosParCov.propagateToDCA(primaryVertex, bz, &dca0);
      }

      if (jPsiDauTracks[1].collisionId() != thisCollId) {
        trackNegParCov.propagateToDCA(primaryVertex, bz, &dca1);
      }

      // ---------------------------------
      // reconstruct J/Psi candidate secondary vertex
      o2::track::TrackParCov trackParCovJPsi{};
      std::array<float, 3> pVecJPsi{};
      registry.fill(HIST("hFitCandidatesJPsi"), SVFitting::BeforeFit);
      try {
        if (df2.process(trackPosParCov, trackNegParCov) == 0) {
          continue;
        }
      } catch (const std::runtime_error& error) {
        LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
        registry.fill(HIST("hFitCandidatesJPsi"), SVFitting::Fail);
        continue;
      }
      registry.fill(HIST("hFitCandidatesJPsi"), SVFitting::FitOk);

      auto secondaryVertexJPsi = df2.getPCACandidate();
      trackPosParCov.propagateTo(secondaryVertexJPsi[0], bz);
      trackNegParCov.propagateTo(secondaryVertexJPsi[0], bz);
      df2.getTrack(0).getPxPyPzGlo(pVec0);
      df2.getTrack(1).getPxPyPzGlo(pVec1);
      pVecJPsi = RecoDecay::pVec(pVec0, pVec1);
      trackParCovJPsi = df2.createParentTrackParCov();
      trackParCovJPsi.setAbsCharge(0); // to be sure

      // TODO: add single track information (min eta, min ITS/TPC clusters, etc.)
      double invMass2JPsiHad{0.};
      for (const auto& trackId : trackIndices) {
        auto trackBach = trackId.template track_as<TTracks>();

        // apply selections on bachelor tracks
        auto trackParCovBach = getTrackParCov(trackBach);
        o2::gpu::gpustd::array<float, 2> dcaBach{trackBach.dcaXY(), trackBach.dcaZ()};
        std::array<float, 3> pVecBach = trackBach.pVector();
        if (trackBach.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovBach, 2.f, noMatCorr, &dcaBach);
          getPxPyPz(trackParCovBach, pVecBach);
        }

        // apply selections on bachelor tracks
        if (!isTrackSelected(trackBach, trackParCovBach, dcaBach, jPsiDauTracks)) {
          continue;
        }

        if constexpr (decChannel == DecayChannel::BplusToJPsiK) {
          registry.fill(HIST("hPtKaon"), trackParCovBach.getPt());
          // compute invariant mass square and apply selection
          invMass2JPsiHad = RecoDecay::m2(std::array{pVecJPsi, pVecBach}, std::array{MassJPsi, MassKPlus});
          if ((invMass2JPsiHad < invMass2JPsiHadMin) || (invMass2JPsiHad > invMass2JPsiHadMax)) {
            continue;
          }
          registry.fill(HIST("hMassJPsiKaon"), std::sqrt(invMass2JPsiHad));

          // fill Kaon tracks table
          // if information on track already stored, go to next track
          if (!selectedTracksBach.count(trackBach.globalIndex())) {
            hfTrackLfDau0(trackBach.globalIndex(), indexHfReducedCollision,
                          trackParCovBach.getX(), trackParCovBach.getAlpha(),
                          trackParCovBach.getY(), trackParCovBach.getZ(), trackParCovBach.getSnp(),
                          trackParCovBach.getTgl(), trackParCovBach.getQ2Pt(),
                          trackBach.itsNCls(), trackBach.tpcNClsCrossedRows(), trackBach.tpcChi2NCl(),
                          trackBach.hasTPC(), trackBach.hasTOF(),
                          trackBach.tpcNSigmaPi(), trackBach.tofNSigmaPi(),
                          trackBach.tpcNSigmaKa(), trackBach.tofNSigmaKa(),
                          trackBach.tpcNSigmaPr(), trackBach.tofNSigmaPr());
            hfTrackCovLfDau0(trackParCovBach.getSigmaY2(), trackParCovBach.getSigmaZY(), trackParCovBach.getSigmaZ2(),
                             trackParCovBach.getSigmaSnpY(), trackParCovBach.getSigmaSnpZ(),
                             trackParCovBach.getSigmaSnp2(), trackParCovBach.getSigmaTglY(), trackParCovBach.getSigmaTglZ(),
                             trackParCovBach.getSigmaTglSnp(), trackParCovBach.getSigmaTgl2(),
                             trackParCovBach.getSigma1PtY(), trackParCovBach.getSigma1PtZ(), trackParCovBach.getSigma1PtSnp(),
                             trackParCovBach.getSigma1PtTgl(), trackParCovBach.getSigma1Pt2());
            // add trackBach.globalIndex() to a list
            // to keep memory of the pions filled in the table and avoid refilling them if they are paired to another JPsi candidate
            // and keep track of their index in hfTrackLfDau0 for McRec purposes
            selectedTracksBach[trackBach.globalIndex()] = hfTrackLfDau0.lastIndex();
          }

          if constexpr (doMc) {
            std::vector<typename TTracks::iterator> beautyHadDauTracks{};
            for (const auto& track : jPsiDauTracks) {
              beautyHadDauTracks.push_back(track);
            }
            beautyHadDauTracks.push_back(trackBach);
            fillMcRecoInfo<DecayChannel::BplusToJPsiK>(collision, particlesMc, beautyHadDauTracks, indexHfCandJPsi, std::array<std::map<int64_t, int64_t>, 2>{selectedTracksBach}, indexCollisionMaxNumContrib);
          }
          fillHfCandJPsi = true;
        } else if constexpr (decChannel == DecayChannel::BsToJPsiPhi) {
          for (auto trackBachId2 = trackId + 1; trackBachId2 != trackIndices.end(); ++trackBachId2) {
            auto trackBach2 = trackBachId2.template track_as<TTracks>();
            auto trackBach2ParCov = getTrackParCov(trackBach2);

            o2::gpu::gpustd::array<float, 2> dcaBach2{trackBach2.dcaXY(), trackBach2.dcaZ()};
            std::array<float, 3> pVecBach2 = trackBach2.pVector();
            if (trackBach2.collisionId() != thisCollId) {
              o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackBach2ParCov, 2.f, noMatCorr, &dcaBach2);
              getPxPyPz(trackBach2ParCov, pVecBach2);
            }

            // apply selections on bachelor tracks
            if (!isTrackSelected(trackBach2, trackBach2ParCov, dcaBach2, jPsiDauTracks)) {
              continue;
            }

            auto dcaBach0 = o2::dataformats::DCA(trackBach.dcaXY(), trackBach.dcaZ(), trackBach.cYY(), trackBach.cZY(), trackBach.cZZ());
            auto dcaBach1 = o2::dataformats::DCA(trackBach2.dcaXY(), trackBach2.dcaZ(), trackBach2.cYY(), trackBach2.cZY(), trackBach2.cZZ());

            // repropagate tracks to this collision if needed
            if (trackBach2.collisionId() != thisCollId) {
              trackBach2ParCov.propagateToDCA(primaryVertex, bz, &dcaBach1);
            }

            // ---------------------------------
            // reconstruct phi candidate secondary vertex
            o2::track::TrackParCov trackParCovPhi{};
            std::array<float, 3> pVecPhi{};
            registry.fill(HIST("hFitCandidatesPhi"), SVFitting::BeforeFit);
            try {
              if (df2.process(trackParCovBach, trackBach2ParCov) == 0) {
                continue;
              }
            } catch (const std::runtime_error& error) {
              LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
              registry.fill(HIST("hFitCandidatesPhi"), SVFitting::Fail);
              continue;
            }
            registry.fill(HIST("hFitCandidatesPhi"), SVFitting::FitOk);

            auto secondaryVertexPhi = df2.getPCACandidate();
            trackParCovBach.propagateTo(secondaryVertexPhi[0], bz);
            trackBach2ParCov.propagateTo(secondaryVertexPhi[0], bz);
            df2.getTrack(0).getPxPyPzGlo(pVecBach);
            df2.getTrack(1).getPxPyPzGlo(pVecBach2);
            pVecPhi = RecoDecay::pVec(pVecBach, pVecBach2);
            trackParCovPhi = df2.createParentTrackParCov();
            trackParCovPhi.setAbsCharge(0); // to be sure
            auto invMassPhi = RecoDecay::m(std::array{pVecBach, pVecBach2}, std::array{MassKPlus, MassKPlus});

            if (std::abs(invMassPhi - MassPhi) > deltaMPhiMax) {
              continue;
            }

            registry.fill(HIST("hPtPhi"), RecoDecay::pt(pVecBach, pVecBach2));
            registry.fill(HIST("hMassPhi"), RecoDecay::m(std::array{pVecBach, pVecBach2}, std::array{MassKPlus, MassKPlus}));
            invMass2JPsiHad = RecoDecay::m2(std::array{pVecJPsi, pVecPhi}, std::array{MassJPsi, MassPhi});
            if ((invMass2JPsiHad < invMass2JPsiHadMin) || (invMass2JPsiHad > invMass2JPsiHadMax)) {
              continue;
            }
            registry.fill(HIST("hMassJPsiPhi"), std::sqrt(invMass2JPsiHad));

            // fill daughter tracks table
            // if information on track already stored, go to next track
            if (!selectedTracksBach.count(trackBach.globalIndex())) {
              hfTrackLfDau0(trackBach.globalIndex(), indexHfReducedCollision,
                            trackParCovBach.getX(), trackParCovBach.getAlpha(),
                            trackParCovBach.getY(), trackParCovBach.getZ(), trackParCovBach.getSnp(),
                            trackParCovBach.getTgl(), trackParCovBach.getQ2Pt(),
                            trackBach.itsNCls(), trackBach.tpcNClsCrossedRows(), trackBach.tpcChi2NCl(),
                            trackBach.hasTPC(), trackBach.hasTOF(),
                            trackBach.tpcNSigmaPi(), trackBach.tofNSigmaPi(),
                            trackBach.tpcNSigmaKa(), trackBach.tofNSigmaKa(),
                            trackBach.tpcNSigmaPr(), trackBach.tofNSigmaPr());
              hfTrackCovLfDau0(trackParCovBach.getSigmaY2(), trackParCovBach.getSigmaZY(), trackParCovBach.getSigmaZ2(),
                               trackParCovBach.getSigmaSnpY(), trackParCovBach.getSigmaSnpZ(),
                               trackParCovBach.getSigmaSnp2(), trackParCovBach.getSigmaTglY(), trackParCovBach.getSigmaTglZ(),
                               trackParCovBach.getSigmaTglSnp(), trackParCovBach.getSigmaTgl2(),
                               trackParCovBach.getSigma1PtY(), trackParCovBach.getSigma1PtZ(), trackParCovBach.getSigma1PtSnp(),
                               trackParCovBach.getSigma1PtTgl(), trackParCovBach.getSigma1Pt2());
              // add trackBach.globalIndex() to a list
              // to keep memory of the pions filled in the table and avoid refilling them if they are paired to another JPsi candidate
              // and keep track of their index in hfTrackLfDau0 for McRec purposes
              selectedTracksBach[trackBach.globalIndex()] = hfTrackLfDau0.lastIndex();
            }

            // fill daughter tracks table
            // if information on track already stored, go to next track
            if (!selectedTracksBach2.count(trackBach2.globalIndex())) {
              hfTrackLfDau1(trackBach2.globalIndex(), indexHfReducedCollision,
                            trackBach2ParCov.getX(), trackBach2ParCov.getAlpha(),
                            trackBach2ParCov.getY(), trackBach2ParCov.getZ(), trackBach2ParCov.getSnp(),
                            trackBach2ParCov.getTgl(), trackBach2ParCov.getQ2Pt(),
                            trackBach2.itsNCls(), trackBach2.tpcNClsCrossedRows(), trackBach2.tpcChi2NCl(),
                            trackBach2.hasTPC(), trackBach2.hasTOF(),
                            trackBach2.tpcNSigmaPi(), trackBach2.tofNSigmaPi(),
                            trackBach2.tpcNSigmaKa(), trackBach2.tofNSigmaKa(),
                            trackBach2.tpcNSigmaPr(), trackBach2.tofNSigmaPr());
              hfTrackCovLfDau1(trackBach2ParCov.getSigmaY2(), trackBach2ParCov.getSigmaZY(), trackBach2ParCov.getSigmaZ2(),
                               trackBach2ParCov.getSigmaSnpY(), trackBach2ParCov.getSigmaSnpZ(),
                               trackBach2ParCov.getSigmaSnp2(), trackBach2ParCov.getSigmaTglY(), trackBach2ParCov.getSigmaTglZ(),
                               trackBach2ParCov.getSigmaTglSnp(), trackBach2ParCov.getSigmaTgl2(),
                               trackBach2ParCov.getSigma1PtY(), trackBach2ParCov.getSigma1PtZ(), trackBach2ParCov.getSigma1PtSnp(),
                               trackBach2ParCov.getSigma1PtTgl(), trackBach2ParCov.getSigma1Pt2());
              // add trackBach2.globalIndex() to a list
              // to keep memory of the pions filled in the table and avoid refilling them if they are paired to another JPsi candidate
              // and keep track of their index in hfTrackLfDau1 for McRec purposes
              selectedTracksBach2[trackBach2.globalIndex()] = hfTrackLfDau1.lastIndex();
            }

            if constexpr (doMc) {
              std::vector<typename TTracks::iterator> beautyHadDauTracks{};
              for (const auto& track : jPsiDauTracks) {
                beautyHadDauTracks.push_back(track);
              }
              beautyHadDauTracks.push_back(trackBach);
              fillMcRecoInfo<DecayChannel::BplusToJPsiK>(collision, particlesMc, beautyHadDauTracks, indexHfCandJPsi, std::array<std::map<int64_t, int64_t>, 2>{selectedTracksBach, selectedTracksBach2}, indexCollisionMaxNumContrib);
            }
            fillHfCandJPsi = true;
          }
        }
      } // kaon loop
      if (fillHfCandJPsi) { // fill JPsi table only once per JPsi candidate
        double invMassJpsi{0.};
        if (runJPsiToee) {
          invMassJpsi = hfHelper.invMassJpsiToEE(candidate);
        } else {
          invMassJpsi = hfHelper.invMassJpsiToMuMu(candidate);
        }
        hfJPsi(trackPos.globalIndex(), trackNeg.globalIndex(),
               indexHfReducedCollision,
               candidate.xSecondaryVertex(), candidate.ySecondaryVertex(), candidate.zSecondaryVertex(),
               invMassJpsi,
               trackPosParCov.getX(), trackNegParCov.getX(),
               trackPosParCov.getY(), trackNegParCov.getY(),
               trackPosParCov.getZ(), trackNegParCov.getZ(),
               trackPosParCov.getAlpha(), trackNegParCov.getAlpha(),
               trackPosParCov.getSnp(), trackNegParCov.getSnp(),
               trackPosParCov.getTgl(), trackNegParCov.getTgl(),
               trackPosParCov.getQ2Pt(), trackNegParCov.getQ2Pt()); // Q/pT
        hfRedJPsiCov(trackPosParCov.getSigmaY2(), trackNegParCov.getSigmaY2(),
                     trackPosParCov.getSigmaZY(), trackNegParCov.getSigmaZY(),
                     trackPosParCov.getSigmaZ2(), trackNegParCov.getSigmaZ2(),
                     trackPosParCov.getSigmaSnpY(), trackNegParCov.getSigmaSnpY(),
                     trackPosParCov.getSigmaSnpZ(), trackNegParCov.getSigmaSnpZ(),
                     trackPosParCov.getSigmaSnp2(), trackNegParCov.getSigmaSnp2(),
                     trackPosParCov.getSigmaTglY(), trackNegParCov.getSigmaTglY(),
                     trackPosParCov.getSigmaTglZ(), trackNegParCov.getSigmaTglZ(),
                     trackPosParCov.getSigmaTglSnp(), trackNegParCov.getSigmaTglSnp(),
                     trackPosParCov.getSigmaTgl2(), trackNegParCov.getSigmaTgl2(),
                     trackPosParCov.getSigma1PtY(), trackNegParCov.getSigma1PtY(),
                     trackPosParCov.getSigma1PtZ(), trackNegParCov.getSigma1PtZ(),
                     trackPosParCov.getSigma1PtSnp(), trackNegParCov.getSigma1PtSnp(),
                     trackPosParCov.getSigma1PtTgl(), trackNegParCov.getSigma1PtTgl(),
                     trackPosParCov.getSigma1Pt2(), trackNegParCov.getSigma1Pt2());
        fillHfReducedCollision = true;
      }
    } // candsJPsi loop

    registry.fill(HIST("hEvents"), 1 + Event::Processed);
    if (!fillHfReducedCollision) {
      registry.fill(HIST("hEvents"), 1 + Event::NoCharmHadPiSelected);
      return;
    }
    registry.fill(HIST("hEvents"), 1 + Event::CharmHadPiSelected);
    float centrality = -1.f;
    uint16_t hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
    // fill collision table if it contains a J/Psi K pair at minimum
    hfReducedCollision(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), hfRejMap, bz);
    hfReducedCollExtra(collision.covXX(), collision.covXY(), collision.covYY(),
                       collision.covXZ(), collision.covYZ(), collision.covZZ());
    // hfReducedCollCentrality(collision.centFT0C(), collision.centFT0M(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange()); // TODO: add
    // if constexpr (withQvec) {
    //   hfReducedQvector(collision.qvecFT0CRe(), collision.qvecFT0CIm(), collision.sumAmplFT0C(),
    //                    collision.qvecFT0ARe(), collision.qvecFT0AIm(), collision.sumAmplFT0A(),
    //                    collision.qvecFT0MRe(), collision.qvecFT0MIm(), collision.sumAmplFT0M(),
    //                    collision.qvecTPCposRe(), collision.qvecTPCposIm(), collision.nTrkTPCpos(),
    //                    collision.qvecTPCnegRe(), collision.qvecTPCnegIm(), collision.nTrkTPCneg(),
    //                    collision.qvecTPCallRe(), collision.qvecTPCallIm(), collision.nTrkTPCall());
    // }
  }

  void processJPsiKData(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                        aod::HfCand2ProngWPid const& candsJPsi,
                        aod::TrackAssoc const& trackIndices,
                        TracksPidWithSel const& tracks,
                        aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBplus(invMassWindowJPsiHad.value);
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
      auto candsJPsiThisColl = candsJPsi.sliceBy(candsJPsiPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, DecayChannel::BplusToJPsiK>(collision, candsJPsiThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorJPsiHadReduced, processJPsiKData, "Process J/Psi K without MC info", true);

  void processJPsiPhiData(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                          aod::HfCand2ProngWPid const& candsJPsi,
                          aod::TrackAssoc const& trackIndices,
                          TracksPidWithSel const& tracks,
                          aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBs(invMassWindowJPsiHad.value);
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
      auto candsJPsiThisColl = candsJPsi.sliceBy(candsJPsiPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, DecayChannel::BsToJPsiPhi>(collision, candsJPsiThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorJPsiHadReduced, processJPsiPhiData, "Process J/Psi phi without MC info", false);

  void processJPsiKMc(CollisionsWCMcLabels const& collisions,
                      aod::HfCand2ProngWPid const& candsJPsi,
                      aod::TrackAssoc const& trackIndices,
                      TracksPidWithSelAndMc const& tracks,
                      aod::McParticles const& particlesMc,
                      aod::BCsWithTimestamps const& bcs,
                      McCollisions const&)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBplus(invMassWindowJPsiHad.value);
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
      auto candsJPsiThisColl = candsJPsi.sliceBy(candsJPsiPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(colPerMcCollision, collision.mcCollisionId());
      int64_t indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, DecayChannel::BplusToJPsiK>(collision, candsJPsiThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    runMcGen<DecayChannel::BplusToJPsiK>(particlesMc);
  }
  PROCESS_SWITCH(HfDataCreatorJPsiHadReduced, processJPsiKMc, "Process J/Psi K with MC info", false);

  void processJPsiPhiMc(CollisionsWCMcLabels const& collisions,
                        aod::HfCand2ProngWPid const& candsJPsi,
                        aod::TrackAssoc const& trackIndices,
                        TracksPidWithSelAndMc const& tracks,
                        aod::McParticles const& particlesMc,
                        aod::BCsWithTimestamps const& bcs,
                        McCollisions const&)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBs(invMassWindowJPsiHad.value);
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
      auto candsJPsiThisColl = candsJPsi.sliceBy(candsJPsiPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(colPerMcCollision, collision.mcCollisionId());
      int64_t indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, DecayChannel::BsToJPsiPhi>(collision, candsJPsiThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    runMcGen<DecayChannel::BsToJPsiPhi>(particlesMc);
  }
  PROCESS_SWITCH(HfDataCreatorJPsiHadReduced, processJPsiPhiMc, "Process J/Psi phi with MC info", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDataCreatorJPsiHadReduced>(cfgc)};
}
