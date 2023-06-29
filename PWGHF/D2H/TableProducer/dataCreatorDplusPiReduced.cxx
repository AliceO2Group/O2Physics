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

/// \file dataCreatorDplusPiReduced.cxx
/// \brief Creation of Dplus-Pi pairs
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

// event types
enum Event : uint8_t {
  Processed = 0,
  NoDPiSelected,
  DPiSelected,
  kNEvent
};

/// Creation of Dplus-Pi pairs
struct HfDataCreatorDplusPiReduced {
  // Produces AOD tables to store track information
  Produces<aod::HfReducedCollisions> hfReducedCollision;
  Produces<aod::HfOriginalCollisionsCounter> hfCollisionCounter;
  Produces<aod::HfTracksReduced> hfTrackPion;
  Produces<aod::HfTracksPidReduced> hfTrackPidPion;
  Produces<aod::HfCand3ProngReduced> hfCand3Prong;
  Produces<aod::HfCandB0Config> rowCandidateConfig;

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
  Configurable<double> invMassWindowB0{"invMassWindowB0", 0.3, "invariant-mass window for B0 candidates"};
  Configurable<int> selectionFlagD{"selectionFlagD", 1, "Selection Flag for D"};

  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber;

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massD = RecoDecay::getMassPDG(pdg::Code::kDMinus);
  double massB0 = RecoDecay::getMassPDG(pdg::Code::kB0);
  double massDPi{0.};
  double invMassD{0.};
  double bz{0.};

  bool isHfCandB0ConfigFilled = false;

  // Fitter to redo D-vertex to get extrapolated daughter tracks (3-prong vertex filter)
  o2::vertexing::DCAFitterN<3> df3;

  using TracksPIDWithSel = soa::Join<aod::BigTracksPIDExtended, aod::TrackSelection>;
  using CandsDFiltered = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagD);

  Preslice<CandsDFiltered> candsDPerCollision = aod::track_association::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    // histograms
    constexpr int kNBinsEvents = kNEvent;
    std::string labels[kNBinsEvents];
    labels[Event::Processed] = "processed";
    labels[Event::NoDPiSelected] = "without DPi pairs";
    labels[Event::DPiSelected] = "with DPi pairs";
    static const AxisSpec axisEvents = {kNBinsEvents, 0.5, kNBinsEvents + 0.5, ""};
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {axisEvents});
    for (int iBin = 0; iBin < kNBinsEvents; iBin++) {
      registry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    registry.add("hMassDToPiKPi", "D^{#minus} candidates;inv. mass (p^{#minus} K^{#plus} #pi^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
    registry.add("hPtD", "D^{#minus} candidates;D^{#minus} candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hPtPion", "#pi^{#plus} candidates;#pi^{#plus} candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hCPAD", "D^{#minus} candidates;D^{#minus} cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});

    // Initialize fitter
    df3.setPropagateToPCA(propagateToPCA);
    df3.setMaxR(maxR);
    df3.setMaxDZIni(maxDZIni);
    df3.setMinParamChange(minParamChange);
    df3.setMinRelChi2Change(minRelChi2Change);
    df3.setUseAbsDCA(useAbsDCA);
    df3.setWeightedFinalPCA(useWeightedFinalPCA);

    // Configure CCDB access
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }

  /// Pion selection (D Pi <-- B0)
  /// \param trackPion is a track with the pion hypothesis
  /// \param track0 is prong0 of selected D candidate
  /// \param track1 is prong1 of selected D candidate
  /// \param track2 is prong2 of selected D candidate
  /// \return true if trackPion passes all cuts
  template <typename T1, typename T2>
  bool isPionSelected(const T1& trackPion, const T2& track0, const T2& track1, const T2& track2)
  {
    // check isGlobalTrackWoDCA status for pions if wanted
    if (usePionIsGlobalTrackWoDCA && !trackPion.isGlobalTrackWoDCA()) {
      return false;
    }
    // minimum pT selection
    if (trackPion.pt() < ptPionMin || !isSelectedTrackDCA(trackPion)) {
      return false;
    }
    // reject pions that are D daughters
    if (trackPion.globalIndex() == track0.globalIndex() || trackPion.globalIndex() == track1.globalIndex() || trackPion.globalIndex() == track2.globalIndex()) {
      return false;
    }
    // reject pi D with same sign as D
    if (trackPion.sign() * track0.sign() > 0) {
      return false;
    }
    return true;
  }

  /// Single-track cuts for pions on dcaXY
  /// \param track is a track
  /// \return true if track passes all cuts
  template <typename T>
  bool isSelectedTrackDCA(const T& track)
  {
    auto pTBinTrack = findBin(binsPtPion, track.pt());
    if (pTBinTrack == -1) {
      return false;
    }

    if (std::abs(track.dcaXY()) < cutsTrackPionDCA->get(pTBinTrack, "min_dcaxytoprimary")) {
      return false; // minimum DCAxy
    }
    if (std::abs(track.dcaXY()) > cutsTrackPionDCA->get(pTBinTrack, "max_dcaxytoprimary")) {
      return false; // maximum DCAxy
    }
    return true;
  }

  void process(aod::Collisions const& collisions,
               CandsDFiltered const& candsD,
               aod::TrackAssoc const& trackIndices,
               TracksPIDWithSel const&,
               aod::BCsWithTimestamps const&)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandB0ConfigFilled) {
      rowCandidateConfig(selectionFlagD.value);
      isHfCandB0ConfigFilled = true;
    }

    static int aodNumber = 0;
    aodNumber++;

    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize());

    static int nCol = 0;
    for (const auto& collision : collisions) {
      nCol++;

      // helpers for ReducedTables filling
      int hfReducedCollisionIndex = hfReducedCollision.lastIndex() + 1;
      std::vector<int64_t> selectedTracksPion;
      selectedTracksPion.reserve(trackIndices.size());
      bool fillHfReducedCollision = false;

      auto primaryVertex = getPrimaryVertex(collision);
      if (nCol % 10000 == 0) {
        LOG(debug) << nCol << " collisions parsed";
      }

      // Set the magnetic field from ccdb.
      // The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      // but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      df3.setBz(bz);

      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      for (const auto& candD : candsDThisColl) {
        bool fillHfCand3Prong = false;
        float invMassD;

        registry.fill(HIST("hMassDToPiKPi"), invMassDplusToPiKPi(candD));
        registry.fill(HIST("hPtD"), candD.pt());
        registry.fill(HIST("hCPAD"), candD.cpa());

        // track0 <-> pi, track1 <-> K, track2 <-> pi
        auto track0 = candD.prong0_as<TracksPIDWithSel>();
        auto track1 = candD.prong1_as<TracksPIDWithSel>();
        auto track2 = candD.prong2_as<TracksPIDWithSel>();
        auto trackParCov0 = getTrackParCov(track0);
        auto trackParCov1 = getTrackParCov(track1);
        auto trackParCov2 = getTrackParCov(track2);

        std::array<float, 3> pVec0 = {track0.px(), track0.py(), track0.pz()};
        std::array<float, 3> pVec1 = {track1.px(), track1.py(), track1.pz()};
        std::array<float, 3> pVec2 = {track2.px(), track2.py(), track2.pz()};

        auto dca0 = o2::dataformats::DCA(track0.dcaXY(), track0.dcaZ(), track0.cYY(), track0.cZY(), track0.cZZ());
        auto dca1 = o2::dataformats::DCA(track1.dcaXY(), track1.dcaZ(), track1.cYY(), track1.cZY(), track1.cZZ());
        auto dca2 = o2::dataformats::DCA(track2.dcaXY(), track2.dcaZ(), track2.cYY(), track2.cZY(), track2.cZZ());

        // repropagate tracks to this collision if needed
        if (track0.collisionId() != thisCollId) {
          trackParCov0.propagateToDCA(primaryVertex, bz, &dca0);
        }

        if (track1.collisionId() != thisCollId) {
          trackParCov1.propagateToDCA(primaryVertex, bz, &dca1);
        }

        if (track2.collisionId() != thisCollId) {
          trackParCov2.propagateToDCA(primaryVertex, bz, &dca2);
        }

        // ---------------------------------
        // reconstruct 3-prong secondary vertex (D±)
        if (df3.process(trackParCov0, trackParCov1, trackParCov2) == 0) {
          continue;
        }

        const auto& secondaryVertexD = df3.getPCACandidate();
        // propagate the 3 prongs to the secondary vertex
        trackParCov0.propagateTo(secondaryVertexD[0], bz);
        trackParCov1.propagateTo(secondaryVertexD[0], bz);
        trackParCov2.propagateTo(secondaryVertexD[0], bz);

        // update pVec of tracks
        df3.getTrack(0).getPxPyPzGlo(pVec0);
        df3.getTrack(1).getPxPyPzGlo(pVec1);
        df3.getTrack(2).getPxPyPzGlo(pVec2);

        // D∓ → π∓ K± π∓
        std::array<float, 3> pVecPiK = RecoDecay::pVec(pVec0, pVec1);
        std::array<float, 3> pVecD = RecoDecay::pVec(pVec0, pVec1, pVec2);
        auto trackParCovPiK = o2::dataformats::V0(df3.getPCACandidatePos(), pVecPiK, df3.calcPCACovMatrixFlat(),
                                                  trackParCov0, trackParCov1, {0, 0}, {0, 0});
        auto trackParCovD = o2::dataformats::V0(df3.getPCACandidatePos(), pVecD, df3.calcPCACovMatrixFlat(),
                                                trackParCovPiK, trackParCov2, {0, 0}, {0, 0});

        auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
        for (const auto& trackId : trackIdsThisCollision) {
          auto trackPion = trackId.track_as<TracksPIDWithSel>();

          // apply selections on pion tracks
          if (!isPionSelected(trackPion, track0, track1, track2) || !isSelectedTrackDCA(trackPion)) {
            continue;
          }
          registry.fill(HIST("hPtPion"), trackPion.pt());
          std::array<float, 3> pVecPion = {trackPion.px(), trackPion.py(), trackPion.pz()};
          // compute invariant mass and apply selection
          massDPi = RecoDecay::m(std::array{pVecD, pVecPion}, std::array{massD, massPi});
          if (std::abs(massDPi - massB0) > invMassWindowB0) {
            continue;
          }

          invMassD = hf_cand_3prong_reduced::invMassDplusToPiKPi(pVec0, pVec1, pVec2);

          // fill Pion tracks table
          // if information on track already stored, go to next track
          if (!std::count(selectedTracksPion.begin(), selectedTracksPion.end(), trackPion.globalIndex())) {
            hfTrackPion(trackPion.globalIndex(), hfReducedCollisionIndex,
                        trackPion.x(), trackPion.alpha(),
                        trackPion.y(), trackPion.z(), trackPion.snp(),
                        trackPion.tgl(), trackPion.signed1Pt(),
                        trackPion.cYY(), trackPion.cZY(), trackPion.cZZ(),
                        trackPion.cSnpY(), trackPion.cSnpZ(),
                        trackPion.cSnpSnp(), trackPion.cTglY(), trackPion.cTglZ(),
                        trackPion.cTglSnp(), trackPion.cTglTgl(),
                        trackPion.c1PtY(), trackPion.c1PtZ(), trackPion.c1PtSnp(),
                        trackPion.c1PtTgl(), trackPion.c1Pt21Pt2(),
                        trackPion.px(), trackPion.py(), trackPion.pz());
            hfTrackPidPion(hfReducedCollisionIndex,
                           trackPion.pt(),
                           trackPion.hasTPC(), trackPion.hasTOF(),
                           trackPion.tpcNSigmaEl(), trackPion.tpcNSigmaMu(), trackPion.tpcNSigmaPi(), trackPion.tpcNSigmaKa(), trackPion.tpcNSigmaPr(),
                           trackPion.tofNSigmaEl(), trackPion.tofNSigmaMu(), trackPion.tofNSigmaPi(), trackPion.tofNSigmaKa(), trackPion.tofNSigmaPr());
            // add trackPion.globalIndex() to a list
            // to keep memory of the pions filled in the table and avoid refilling them if they are paired to another D candidate
            selectedTracksPion.emplace_back(trackPion.globalIndex());
          }
          fillHfCand3Prong = true;
        }                       // pion loop
        if (fillHfCand3Prong) { // fill candDplus table only once per D candidate
          hfCand3Prong(track0.globalIndex(), track1.globalIndex(), track2.globalIndex(),
                       hfReducedCollisionIndex,
                       trackParCovD.getX(), trackParCovD.getAlpha(),
                       trackParCovD.getY(), trackParCovD.getZ(), trackParCovD.getSnp(),
                       trackParCovD.getTgl(), trackParCovD.getQ2Pt(),
                       trackParCovD.getSigmaY2(), trackParCovD.getSigmaZY(), trackParCovD.getSigmaZ2(),
                       trackParCovD.getSigmaSnpY(), trackParCovD.getSigmaSnpZ(),
                       trackParCovD.getSigmaSnp2(), trackParCovD.getSigmaTglY(), trackParCovD.getSigmaTglZ(),
                       trackParCovD.getSigmaTglSnp(), trackParCovD.getSigmaTgl2(),
                       trackParCovD.getSigma1PtY(), trackParCovD.getSigma1PtZ(), trackParCovD.getSigma1PtSnp(),
                       trackParCovD.getSigma1PtTgl(), trackParCovD.getSigma1Pt2(),
                       pVecD[0], pVecD[1], pVecD[2],
                       candD.cpa(),
                       candD.decayLength(),
                       invMassD);
          fillHfReducedCollision = true;
        }
      } // candsD loop
      registry.fill(HIST("hEvents"), 1 + Event::Processed);
      if (!fillHfReducedCollision) {
        registry.fill(HIST("hEvents"), 1 + Event::NoDPiSelected);
        continue;
      }
      registry.fill(HIST("hEvents"), 1 + Event::DPiSelected);
      // fill collision table if it contains a DPi pair a minima
      hfReducedCollision(collision.posX(), collision.posY(), collision.posZ(),
                         collision.covXX(), collision.covXY(), collision.covYY(),
                         collision.covXZ(), collision.covYZ(), collision.covZZ(),
                         bz);
    } // collision
  }   // process
};    // struct

/// Performs MC matching.
struct HfDataCreatorDplusPiReducedMc {
  Produces<aod::HfDPiMcRecReduced> rowHfDPiMcRecReduced;
  Produces<aod::HfB0McGenReduced> rowHfB0McGenReduced;

  void init(InitContext const&) {}

  void processMc(aod::HfCand3ProngReduced const& candsD,
                 aod::HfTracksReduced const& tracksPion,
                 aod::BigTracksMC const&,
                 aod::McParticles const& particlesMc)
  {
    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t debug = 0;

    for (const auto& candD : candsD) {
      auto arrayDaughtersD = array{candD.prong0_as<aod::BigTracksMC>(),
                                   candD.prong1_as<aod::BigTracksMC>(),
                                   candD.prong2_as<aod::BigTracksMC>()};

      for (const auto& trackPion : tracksPion) {
        if (trackPion.hfReducedCollisionId() != candD.hfReducedCollisionId()) {
          continue;
        }
        // const auto& trackId = trackPion.globalIndex();
        auto arrayDaughtersB0 = array{candD.prong0_as<aod::BigTracksMC>(),
                                      candD.prong1_as<aod::BigTracksMC>(),
                                      candD.prong2_as<aod::BigTracksMC>(),
                                      trackPion.track_as<aod::BigTracksMC>()};
        // B0 → D- π+ → (π- K+ π-) π+
        // Printf("Checking B0 → D- π+");
        indexRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughtersB0, pdg::Code::kB0, array{-kPiPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec > -1) {
          // D- → π- K+ π-
          // Printf("Checking D- → π- K+ π-");
          indexRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughtersD, pdg::Code::kDMinus, array{-kPiPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
          if (indexRec > -1) {
            flag = sign * BIT(hf_cand_b0::DecayType::B0ToDPi);
          } else {
            debug = 1;
            LOGF(info, "WARNING: B0 decays in the expected final state but the condition on the intermediate state is not fulfilled");
          }
        }
        auto indexMother = RecoDecay::getMother(particlesMc, trackPion.track_as<aod::BigTracksMC>().mcParticle_as<aod::McParticles>(), pdg::Code::kB0, true);
        auto particleMother = particlesMc.rawIteratorAt(indexMother);

        rowHfDPiMcRecReduced(candD.globalIndex(), trackPion.globalIndex(), flag, origin, debug, particleMother.pt());
      }
    } // rec

    // Match generated particles.
    for (auto const& particle : particlesMc) {
      // Printf("New gen. candidate");
      flag = 0;
      origin = 0;
      // B0 → D- π+
      if (RecoDecay::isMatchedMCGen(particlesMc, particle, pdg::Code::kB0, array{-static_cast<int>(pdg::Code::kDPlus), +kPiPlus}, true)) {
        // Match D- -> π- K+ π-
        auto candDMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
        // Printf("Checking D- -> π- K+ π-");
        if (RecoDecay::isMatchedMCGen(particlesMc, candDMC, -static_cast<int>(pdg::Code::kDPlus), array{-kPiPlus, +kKPlus, -kPiPlus}, true, &sign)) {
          flag = sign * BIT(hf_cand_b0::DecayType::B0ToDPi);
        }
      }

      // save information for B0 task
      if (!TESTBIT(std::abs(flag), hf_cand_b0::DecayType::B0ToDPi)) {
        continue;
      }

      auto ptParticle = particle.pt();
      auto yParticle = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(pdg::Code::kB0));
      auto etaParticle = particle.eta();

      std::array<float, 2> ptProngs;
      std::array<float, 2> yProngs;
      std::array<float, 2> etaProngs;
      int counter = 0;
      for (auto const& daught : particle.daughters_as<aod::McParticles>()) {
        ptProngs[counter] = daught.pt();
        etaProngs[counter] = daught.eta();
        yProngs[counter] = RecoDecay::y(array{daught.px(), daught.py(), daught.pz()}, RecoDecay::getMassPDG(daught.pdgCode()));
        counter++;
      }
      rowHfB0McGenReduced(flag, origin,
                          ptParticle, yParticle, etaParticle,
                          ptProngs[0], yProngs[0], etaProngs[0],
                          ptProngs[1], yProngs[1], etaProngs[1]);
    } // gen

  } // processMc
  PROCESS_SWITCH(HfDataCreatorDplusPiReducedMc, processMc, "Process MC", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfDataCreatorDplusPiReduced>(cfgc),
    adaptAnalysisTask<HfDataCreatorDplusPiReducedMc>(cfgc)};
}
