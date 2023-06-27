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

/// \file dataCreatorD0PiReduced.cxx
/// \brief Creation of D0-Pi pairs
///
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociation.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

// event types
enum Event : uint8_t {
  Processed = 0,
  NoD0PiSelected,
  D0PiSelected,
  kNEvent
};

/// Creation of D0-Pi pairs
struct HfDataCreatorD0PiReduced {
  // Produces AOD tables to store track information
  Produces<aod::HfReducedCollisions> hfReducedCollision;
  Produces<aod::HfOriginalCollisionsCounter> hfCollisionCounter;
  Produces<aod::HfTracksReduced> hfTrackPion;
  Produces<aod::HfTracksPidReduced> hfTrackPidPion;
  Produces<aod::HfCand2ProngReduced> hfCand2Prong;
  Produces<aod::HfCandBpConfig> rowCandidateConfig;

  // vertexing
  // Configurable<double> bz{"bz", 5., "magnetic field"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any B+ is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // selection
  Configurable<bool> usePionIsGlobalTrackWoDCA{"usePionIsGlobalTrackWoDCA", true, "check isGlobalTrackWoDCA status for pions, for Run3 studies"};
  Configurable<double> ptPionMin{"ptPionMin", 0.5, "minimum pion pT threshold (GeV/c)"};
  Configurable<std::vector<double>> binsPtPion{"binsPtPion", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for pion DCA XY pT-dependent cut"};
  Configurable<LabeledArray<double>> cutsTrackPionDCA{"cutsTrackPionDCA", {hf_cuts_single_track::cutsTrack[0], hf_cuts_single_track::nBinsPtTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for pions"};
  Configurable<double> invMassWindowBplus{"invMassWindowBplus", 0.3, "invariant-mass window for B+ candidates"};
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGeo{"ccdbPathGeo", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber;

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massD0 = RecoDecay::getMassPDG(pdg::Code::kD0);
  double massBplus = RecoDecay::getMassPDG(pdg::Code::kBPlus);
  double massD0Pi{0.};
  double invMassD0{0.};
  double bz{0.};

  bool isHfCandBplusConfigFilled = false;

  // Fitter to redo D0-vertex to get extrapolated daughter tracks (2-prong vertex filter)
  o2::vertexing::DCAFitterN<2> df2;

  using TracksPIDWithSel = soa::Join<aod::BigTracksPIDExtended, aod::TrackSelection>;
  using CandsDFiltered = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar);

  Preslice<CandsDFiltered> candsDPerCollision = aod::track_association::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    // histograms
    constexpr int kNBinsEvents = kNEvent;
    std::string labels[kNBinsEvents];
    labels[Event::Processed] = "processed";
    labels[Event::NoD0PiSelected] = "without D0Pi pairs";
    labels[Event::D0PiSelected] = "with D0Pi pairs";
    static const AxisSpec axisEvents = {kNBinsEvents, 0.5, kNBinsEvents + 0.5, ""};
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {axisEvents});
    for (int iBin = 0; iBin < kNBinsEvents; iBin++) {
      registry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    registry.add("hMassD0ToKPi", "D^{0}} candidates;inv. mass (K^{#minus} #pi^{#plus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
    registry.add("hPtD0", "D^{0} candidates;D^{0} candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hPtPion", "#pi^{#plus} candidates;#pi^{#plus} candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hCPAD0", "D^{0} candidates;D^{0} cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});

    // Initialize fitter
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);

    // Configure CCDB access
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(ccdbPathGeo);
    }
    runNumber = 0;
  }

  /// Pion selection (D0 Pi <-- B+)
  /// \param trackPion is a track with the pion hypothesis
  /// \param track0 is prong0 of selected D0 candidate
  /// \param track1 is prong1 of selected D0 candidate
  /// \param candD0 is the D0 candidate
  /// \return true if trackPion passes all cuts
  template <typename T1, typename T2, typename T3>
  bool isPionSelected(const T1& trackPion, const T2& track0, const T2& track1, const T3& candD0)
  {
    // check isGlobalTrackWoDCA status for pions if wanted
    if (usePionIsGlobalTrackWoDCA && !trackPion.isGlobalTrackWoDCA()) {
      return false;
    }
    // minimum pT selection
    if (trackPion.pt() < ptPionMin || !isSelectedTrackDCA(trackPion)) {
      return false;
    }
    // reject pion not compatible with D0/D0bar hypothesis
    if (!((candD0.isSelD0() >= selectionFlagD0 && trackPion.sign() < 0) || (candD0.isSelD0bar() >= selectionFlagD0bar && trackPion.sign() > 0))) {
      // Printf("D0: %d, D0bar%d, sign: %d", candD0.isSelD0(), candD0.isSelD0bar(), track.sign());
      return false;
    }
    // reject pions that are D daughters
    if (trackPion.globalIndex() == track0.globalIndex() || trackPion.globalIndex() == track1.globalIndex()) {
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
               CandsDFiltered const& candsD0,
               aod::TrackAssoc const& trackIndices,
               TracksPIDWithSel const&,
               aod::BCsWithTimestamps const&)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBplusConfigFilled) {
      rowCandidateConfig(selectionFlagD0.value, selectionFlagD0bar.value);
      isHfCandBplusConfigFilled = true;
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
      df2.setBz(bz);

      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD0.sliceBy(candsDPerCollision, thisCollId);
      for (const auto& candD0 : candsDThisColl) {
        bool fillHfCand2Prong = false;
        float invMassD0;

        if (candD0.isSelD0() >= selectionFlagD0) {
          invMassD0 = invMassD0ToPiK(candD0);
          registry.fill(HIST("hMassD0ToKPi"), invMassD0);
        }
        if (candD0.isSelD0bar() >= selectionFlagD0bar) {
          invMassD0 = invMassD0barToKPi(candD0);
          registry.fill(HIST("hMassD0ToKPi"), invMassD0);
        }
        registry.fill(HIST("hPtD0"), candD0.pt());
        registry.fill(HIST("hCPAD0"), candD0.cpa());

        auto track0 = candD0.prong0_as<TracksPIDWithSel>();
        auto track1 = candD0.prong1_as<TracksPIDWithSel>();
        auto trackParCov0 = getTrackParCov(track0);
        auto trackParCov1 = getTrackParCov(track1);

        std::array<float, 3> pVec0 = {track0.px(), track0.py(), track0.pz()};
        std::array<float, 3> pVec1 = {track1.px(), track1.py(), track1.pz()};

        auto dca0 = o2::dataformats::DCA(track0.dcaXY(), track0.dcaZ(), track0.cYY(), track0.cZY(), track0.cZZ());
        auto dca1 = o2::dataformats::DCA(track1.dcaXY(), track1.dcaZ(), track1.cYY(), track1.cZY(), track1.cZZ());

        // repropagate tracks to this collision if needed
        if (track0.collisionId() != thisCollId) {
          trackParCov0.propagateToDCA(primaryVertex, bz, &dca0);
        }

        if (track1.collisionId() != thisCollId) {
          trackParCov1.propagateToDCA(primaryVertex, bz, &dca1);
        }

        // ---------------------------------
        // reconstruct 2-prong secondary vertex (D0)
        if (df2.process(trackParCov0, trackParCov1) == 0) {
          continue;
        }

        const auto& secondaryVertexD0 = df2.getPCACandidate();
        // propagate the 3 prongs to the secondary vertex
        trackParCov0.propagateTo(secondaryVertexD0[0], bz);
        trackParCov1.propagateTo(secondaryVertexD0[0], bz);

        // update pVec of tracks
        df2.getTrack(0).getPxPyPzGlo(pVec0);
        df2.getTrack(1).getPxPyPzGlo(pVec1);

        // D0(bar) → π∓ K±
        std::array<float, 3> pVecD0 = RecoDecay::pVec(pVec0, pVec1);
        auto trackParCovD0 = o2::dataformats::V0(df2.getPCACandidatePos(), pVecD0, df2.calcPCACovMatrixFlat(),
                                                 trackParCov0, trackParCov1, {0, 0}, {0, 0});

        auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
        for (const auto& trackId : trackIdsThisCollision) {
          auto trackPion = trackId.track_as<TracksPIDWithSel>();

          // apply selections on pion tracks
          if (!isPionSelected(trackPion, track0, track1, candD0) || !isSelectedTrackDCA(trackPion)) {
            continue;
          }
          registry.fill(HIST("hPtPion"), trackPion.pt());
          std::array<float, 3> pVecPion = {trackPion.px(), trackPion.py(), trackPion.pz()};
          // compute invariant mass and apply selection
          massD0Pi = RecoDecay::m(std::array{pVecD0, pVecPion}, std::array{massD0, massPi});
          if (std::abs(massD0Pi - massBplus) > invMassWindowBplus) {
            continue;
          }

          /* if (candD0.isSelD0() >= selectionFlagD0) {
            invMassD0 = hf_cand_2prong_reduced::invMassD0barToKPi(pVec0, pVec1);
          }
          if (candD0.isSelD0bar() >= selectionFlagD0bar) {
            invMassD0 = hf_cand_2prong_reduced::invMassD0ToPiK(pVec0, pVec1);
          } */
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
          fillHfCand2Prong = true;
        }                       // pion loop
        if (fillHfCand2Prong) { // fill candD0 table only once per D0 candidate
          hfCand2Prong(track0.globalIndex(), track1.globalIndex(),
                       hfReducedCollisionIndex,
                       trackParCovD0.getX(), trackParCovD0.getAlpha(),
                       trackParCovD0.getY(), trackParCovD0.getZ(), trackParCovD0.getSnp(),
                       trackParCovD0.getTgl(), trackParCovD0.getQ2Pt(),
                       trackParCovD0.getSigmaY2(), trackParCovD0.getSigmaZY(), trackParCovD0.getSigmaZ2(),
                       trackParCovD0.getSigmaSnpY(), trackParCovD0.getSigmaSnpZ(),
                       trackParCovD0.getSigmaSnp2(), trackParCovD0.getSigmaTglY(), trackParCovD0.getSigmaTglZ(),
                       trackParCovD0.getSigmaTglSnp(), trackParCovD0.getSigmaTgl2(),
                       trackParCovD0.getSigma1PtY(), trackParCovD0.getSigma1PtZ(), trackParCovD0.getSigma1PtSnp(),
                       trackParCovD0.getSigma1PtTgl(), trackParCovD0.getSigma1Pt2(),
                       pVecD0[0], pVecD0[1], pVecD0[2],
                       candD0.cpa(),
                       candD0.decayLength(),
                       invMassD0);
          fillHfReducedCollision = true;
        }
      } // candsD loop
      registry.fill(HIST("hEvents"), 1 + Event::Processed);
      if (!fillHfReducedCollision) {
        registry.fill(HIST("hEvents"), 1 + Event::NoD0PiSelected);
        continue;
      }
      registry.fill(HIST("hEvents"), 1 + Event::D0PiSelected);
      // fill collision table if it contains a D0Pi pair a minima
      hfReducedCollision(collision.posX(), collision.posY(), collision.posZ(),
                         collision.covXX(), collision.covXY(), collision.covYY(),
                         collision.covXZ(), collision.covYZ(), collision.covZZ(),
                         bz);
    } // collision
  }   // process
};    // struct

/// Performs MC matching.
struct HfDataCreatorD0PiReducedMc {
  Produces<aod::HfD0PiMcRecReduced> rowHfD0PiMcRecReduced;
  Produces<aod::HfBpMcGenReduced> rowHfBPMcGenReduced;

  void init(InitContext const&) {}

  void processMc(aod::HfCand2ProngReduced const& candsD0,
                 aod::HfTracksReduced const& tracksPion,
                 aod::BigTracksMC const&,
                 aod::McParticles const& particlesMc)
  {
    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t debug = 0;

    for (const auto& candD0 : candsD0) {
      auto arrayDaughtersD0 = array{candD0.prong0_as<aod::BigTracksMC>(),
                                    candD0.prong1_as<aod::BigTracksMC>()};

      for (const auto& trackPion : tracksPion) {
        if (trackPion.hfReducedCollisionId() != candD0.hfReducedCollisionId()) {
          continue;
        }
        // const auto& trackId = trackPion.globalIndex();
        auto arrayDaughtersBplus = array{candD0.prong0_as<aod::BigTracksMC>(),
                                         candD0.prong1_as<aod::BigTracksMC>(),
                                         trackPion.track_as<aod::BigTracksMC>()};
        // B+ → D0(bar) π+ → (K+ π-) π+
        // Printf("Checking B+ → D0bar π+");
        indexRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughtersBplus, pdg::Code::kBPlus, array{+kPiPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
        if (indexRec > -1) {
          // D0bar → K+ π-
          // Printf("Checking D0bar → K+ π-");
          indexRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughtersD0, -pdg::Code::kD0, array{-kPiPlus, +kKPlus}, true, &sign, 1);
          if (indexRec > -1) {
            flag = sign * BIT(hf_cand_bplus::DecayType::BplusToD0Pi);
          } else {
            debug = 1;
            LOGF(info, "WARNING: B+ decays in the expected final state but the condition on the intermediate state is not fulfilled");
          }
        }
        auto indexMother = RecoDecay::getMother(particlesMc, trackPion.track_as<aod::BigTracksMC>().mcParticle_as<aod::McParticles>(), pdg::Code::kBPlus, true);
        auto particleMother = particlesMc.rawIteratorAt(indexMother);

        rowHfD0PiMcRecReduced(candD0.globalIndex(), trackPion.globalIndex(), flag, origin, particleMother.pt());
      }
    } // rec

    // Match generated particles.
    for (auto const& particle : particlesMc) {
      // Printf("New gen. candidate");
      flag = 0;
      origin = 0;
      // B+ → D0bar π+
      if (RecoDecay::isMatchedMCGen(particlesMc, particle, pdg::Code::kBPlus, array{-static_cast<int>(pdg::Code::kD0), +kPiPlus}, true)) {
        // Match D0bar -> π- K+
        auto candD0MC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
        // Printf("Checking D0bar -> π- K+");
        if (RecoDecay::isMatchedMCGen(particlesMc, candD0MC, -static_cast<int>(pdg::Code::kD0), array{-kPiPlus, +kKPlus}, true, &sign)) {
          flag = sign * BIT(hf_cand_bplus::DecayType::BplusToD0Pi);
        }
      }

      // save information for B+ task
      if (!TESTBIT(std::abs(flag), hf_cand_bplus::DecayType::BplusToD0Pi)) {
        continue;
      }

      auto ptParticle = particle.pt();
      auto yParticle = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(pdg::Code::kBPlus));
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
      rowHfBPMcGenReduced(flag, origin,
                          ptParticle, yParticle, etaParticle,
                          ptProngs[0], yProngs[0], etaProngs[0],
                          ptProngs[1], yProngs[1], etaProngs[1]);
    } // gen

  } // processMc
  PROCESS_SWITCH(HfDataCreatorD0PiReducedMc, processMc, "Process MC", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfDataCreatorD0PiReduced>(cfgc),
    adaptAnalysisTask<HfDataCreatorD0PiReducedMc>(cfgc)};
}
