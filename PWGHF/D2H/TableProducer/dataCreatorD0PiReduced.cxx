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

#include <map>

#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::analysis;
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
  Produces<aod::HfRedCollisions> hfReducedCollision;
  Produces<aod::HfOrigColCounts> hfCollisionCounter;
  Produces<aod::HfRedTrackBases> hfTrackPion;
  Produces<aod::HfRedTracksCov> hfTrackCovPion;
  Produces<aod::HfRedTracksPid> hfTrackPidPion;
  Produces<aod::HfRed2Prongs> hfCand2Prong;
  Produces<aod::HfRed2ProngsCov> hfCand2ProngCov;
  Produces<aod::HfCandBpConfigs> rowCandidateConfig;
  Produces<aod::HfMcRecRedD0Pis> rowHfD0PiMcRecReduced;
  Produces<aod::HfMcGenRedBps> rowHfBpMcGenReduced;

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
  Configurable<double> invMassWindowD0Pi{"invMassWindowD0Pi", 0.3, "invariant-mass window for D0Pi pair preselections (GeV/c2)"};
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGeo{"ccdbPathGeo", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  HfHelper hfHelper;

  // CCDB service
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber;

  // O2DatabasePDG service
  Service<o2::framework::O2DatabasePDG> pdg;

  double massPi{0.};
  double massD0{0.};
  double massBplus{0.};
  double invMass2D0PiMin{0.};
  double invMass2D0PiMax{0.};
  double bz{0.};

  bool isHfCandBplusConfigFilled = false;

  // Fitter to redo D0-vertex to get extrapolated daughter tracks (2-prong vertex filter)
  o2::vertexing::DCAFitterN<2> df2;

  using TracksPidAll = soa::Join<aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using TracksPIDWithSel = soa::Join<aod::TracksWCovDcaExtra, TracksPidAll, aod::TrackSelection>;
  using TracksPIDWithSelAndMc = soa::Join<TracksPIDWithSel, aod::McTrackLabels>;
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
    runNumber = 0;

    // invariant-mass window cut
    massPi = o2::analysis::pdg::MassPiPlus;
    massD0 = o2::analysis::pdg::MassD0;
    massBplus = o2::analysis::pdg::MassBPlus;
    invMass2D0PiMin = (massBplus - invMassWindowD0Pi) * (massBplus - invMassWindowD0Pi);
    invMass2D0PiMax = (massBplus + invMassWindowD0Pi) * (massBplus + invMassWindowD0Pi);
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
      // LOGF(debug, "D0: %d, D0bar%d, sign: %d", candD0.isSelD0(), candD0.isSelD0bar(), track.sign());
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

  template <bool doMc = false, typename P, typename T>
  void runDataCreation(aod::Collisions const& collisions,
                       CandsDFiltered const& candsD0,
                       aod::TrackAssoc const& trackIndices,
                       T const& tracks,
                       P const& particlesMc,
                       aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBplusConfigFilled) {
      rowCandidateConfig(selectionFlagD0.value, selectionFlagD0bar.value, invMassWindowD0Pi.value);
      isHfCandBplusConfigFilled = true;
    }

    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize());

    for (const auto& collision : collisions) {

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
        int indexHfCand2Prong = hfCand2Prong.lastIndex() + 1;
        bool fillHfCand2Prong = false;
        float invMassD0{-1.f}, invMassD0bar{-1.f};

        if (candD0.isSelD0() >= selectionFlagD0) {
          invMassD0 = hfHelper.invMassD0ToPiK(candD0);
          registry.fill(HIST("hMassD0ToKPi"), invMassD0);
        }
        if (candD0.isSelD0bar() >= selectionFlagD0bar) {
          invMassD0bar = hfHelper.invMassD0barToKPi(candD0);
          registry.fill(HIST("hMassD0ToKPi"), invMassD0bar);
        }
        registry.fill(HIST("hPtD0"), candD0.pt());
        registry.fill(HIST("hCPAD0"), candD0.cpa());

        auto track0 = candD0.template prong0_as<T>();
        auto track1 = candD0.template prong1_as<T>();
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
        // propagate the 2 prongs to the secondary vertex
        trackParCov0.propagateTo(secondaryVertexD0[0], bz);
        trackParCov1.propagateTo(secondaryVertexD0[0], bz);

        // update pVec of tracks
        df2.getTrack(0).getPxPyPzGlo(pVec0);
        df2.getTrack(1).getPxPyPzGlo(pVec1);

        // D0(bar) → π∓ K±
        std::array<float, 3> pVecD0 = RecoDecay::pVec(pVec0, pVec1);
        auto trackParCovD0 = o2::dataformats::V0(df2.getPCACandidatePos(), pVecD0, df2.calcPCACovMatrixFlat(), trackParCov0, trackParCov1);

        auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
        for (const auto& trackId : trackIdsThisCollision) {
          auto trackPion = trackId.template track_as<T>();

          // apply selections on pion tracks
          if (!isPionSelected(trackPion, track0, track1, candD0) || !isSelectedTrackDCA(trackPion)) {
            continue;
          }
          registry.fill(HIST("hPtPion"), trackPion.pt());
          std::array<float, 3> pVecPion = {trackPion.px(), trackPion.py(), trackPion.pz()};
          // compute invariant mass square and apply selection
          auto invMass2D0Pi = RecoDecay::m2(std::array{pVecD0, pVecPion}, std::array{massD0, massPi});
          if ((invMass2D0Pi < invMass2D0PiMin) || (invMass2D0Pi > invMass2D0PiMax)) {
            continue;
          }

          // fill Pion tracks table
          // if information on track already stored, go to next track
          if (!selectedTracksPion.count(trackPion.globalIndex())) {
            hfTrackPion(trackPion.globalIndex(), indexHfReducedCollision,
                        trackPion.x(), trackPion.alpha(),
                        trackPion.y(), trackPion.z(), trackPion.snp(),
                        trackPion.tgl(), trackPion.signed1Pt());
            hfTrackCovPion(trackPion.cYY(), trackPion.cZY(), trackPion.cZZ(),
                           trackPion.cSnpY(), trackPion.cSnpZ(),
                           trackPion.cSnpSnp(), trackPion.cTglY(), trackPion.cTglZ(),
                           trackPion.cTglSnp(), trackPion.cTglTgl(),
                           trackPion.c1PtY(), trackPion.c1PtZ(), trackPion.c1PtSnp(),
                           trackPion.c1PtTgl(), trackPion.c1Pt21Pt2());
            hfTrackPidPion(trackPion.hasTPC(), trackPion.hasTOF(),
                           trackPion.tpcNSigmaPi(), trackPion.tofNSigmaPi());
            // add trackPion.globalIndex() to a list
            // to keep memory of the pions filled in the table and avoid refilling them if they are paired to another D candidate
            // and keep track of their index in hfTrackPion for McRec purposes
            selectedTracksPion[trackPion.globalIndex()] = hfTrackPion.lastIndex();
          }

          if constexpr (doMc) {
            // we check the MC matching to be stored
            auto arrayDaughtersD0 = std::array{track0, track1};
            auto arrayDaughtersBplus = std::array{track0, track1, trackPion};
            int8_t sign{0};
            int8_t flag{0};
            // B+ → D0(bar) π+ → (K+ π-) π+
            // Printf("Checking B+ → D0bar π+");
            auto indexRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughtersBplus, pdg::Code::kBPlus, std::array{+kPiPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
            if (indexRec > -1) {
              // D0bar → K+ π-
              // Printf("Checking D0bar → K+ π-");
              indexRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughtersD0, pdg::Code::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign, 1);
              if (indexRec > -1) {
                flag = sign * BIT(hf_cand_bplus::DecayType::BplusToD0Pi);
              } else {
                LOGF(info, "WARNING: B+ decays in the expected final state but the condition on the intermediate state is not fulfilled");
              }
            }
            auto indexMother = RecoDecay::getMother(particlesMc, trackPion.template mcParticle_as<P>(), pdg::Code::kBPlus, true);
            auto particleMother = particlesMc.rawIteratorAt(indexMother);

            rowHfD0PiMcRecReduced(indexHfCand2Prong, selectedTracksPion[trackPion.globalIndex()], flag, particleMother.pt());
          }
          fillHfCand2Prong = true;
        }                       // pion loop
        if (fillHfCand2Prong) { // fill candD0 table only once per D0 candidate
          hfCand2Prong(track0.globalIndex(), track1.globalIndex(),
                       indexHfReducedCollision,
                       trackParCovD0.getX(), trackParCovD0.getAlpha(),
                       trackParCovD0.getY(), trackParCovD0.getZ(), trackParCovD0.getSnp(),
                       trackParCovD0.getTgl(), trackParCovD0.getQ2Pt(),
                       candD0.xSecondaryVertex(), candD0.ySecondaryVertex(), candD0.zSecondaryVertex(), invMassD0, invMassD0bar);
          hfCand2ProngCov(trackParCovD0.getSigmaY2(), trackParCovD0.getSigmaZY(), trackParCovD0.getSigmaZ2(),
                          trackParCovD0.getSigmaSnpY(), trackParCovD0.getSigmaSnpZ(),
                          trackParCovD0.getSigmaSnp2(), trackParCovD0.getSigmaTglY(), trackParCovD0.getSigmaTglZ(),
                          trackParCovD0.getSigmaTglSnp(), trackParCovD0.getSigmaTgl2(),
                          trackParCovD0.getSigma1PtY(), trackParCovD0.getSigma1PtZ(), trackParCovD0.getSigma1PtSnp(),
                          trackParCovD0.getSigma1PtTgl(), trackParCovD0.getSigma1Pt2());
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

    if constexpr (doMc) {
      // Match generated particles.
      for (const auto& particle : particlesMc) {
        int8_t sign{0};
        int8_t flag{0};
        // B+ → D0bar π+
        if (RecoDecay::isMatchedMCGen(particlesMc, particle, pdg::Code::kBPlus, std::array{static_cast<int>(pdg::Code::kD0), +kPiPlus}, true)) {
          // Match D0bar -> π- K+
          auto candD0MC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
          // Printf("Checking D0bar -> π- K+");
          if (RecoDecay::isMatchedMCGen(particlesMc, candD0MC, static_cast<int>(pdg::Code::kD0), std::array{+kPiPlus, -kKPlus}, true, &sign)) {
            flag = sign * BIT(hf_cand_bplus::DecayType::BplusToD0Pi);
          }
        }

        // save information for B+ task
        if (!TESTBIT(std::abs(flag), hf_cand_bplus::DecayType::BplusToD0Pi)) {
          continue;
        }

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, massBplus);
        auto etaParticle = particle.eta();

        std::array<float, 2> ptProngs;
        std::array<float, 2> yProngs;
        std::array<float, 2> etaProngs;
        int counter = 0;
        for (const auto& daught : particle.template daughters_as<P>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(std::array{daught.px(), daught.py(), daught.pz()}, pdg->Mass(daught.pdgCode()));
          counter++;
        }
        rowHfBpMcGenReduced(flag, ptParticle, yParticle, etaParticle,
                            ptProngs[0], yProngs[0], etaProngs[0],
                            ptProngs[1], yProngs[1], etaProngs[1]);
      } // gen
    }
  }

  void processData(aod::Collisions const& collisions,
                   CandsDFiltered const& candsD0,
                   aod::TrackAssoc const& trackIndices,
                   TracksPIDWithSel const& tracks,
                   aod::BCsWithTimestamps const& bcs)
  {
    runDataCreation<false>(collisions, candsD0, trackIndices, tracks, tracks, bcs);
  }
  PROCESS_SWITCH(HfDataCreatorD0PiReduced, processData, "Process without MC info", true);

  void processMc(aod::Collisions const& collisions,
                 CandsDFiltered const& candsD0,
                 aod::TrackAssoc const& trackIndices,
                 TracksPIDWithSelAndMc const& tracks,
                 aod::McParticles const& particlesMc,
                 aod::BCsWithTimestamps const& bcs)
  {
    runDataCreation<true>(collisions, candsD0, trackIndices, tracks, particlesMc, bcs);
  }
  PROCESS_SWITCH(HfDataCreatorD0PiReduced, processMc, "Process with MC info", false);

}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDataCreatorD0PiReduced>(cfgc)};
}
