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

/// \file candidateCreatorB0.cxx
/// \brief Reconstruction of B0 candidates
/// \note Adapted from candidateCreatorXicc.cxx
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include <vector>
#include <string>
#include <memory>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"
#include "PWGHF/Utils/utilsMcGen.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_trkcandsel;

/// Reconstruction of B0 candidates
struct HfCandidateCreatorB0 {
  Produces<aod::HfCandB0Base> rowCandidateBase;     // table defined in CandidateReconstructionTables.h
  Produces<aod::HfCandB0Prongs> rowCandidateProngs; // table defined in CandidateReconstructionTables.h

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
  Configurable<LabeledArray<double>> cutsTrackPionDCA{"cutsTrackPionDCA", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for pions"};
  Configurable<double> invMassWindowB0{"invMassWindowB0", 0.3, "invariant-mass window for B0 candidates"};
  Configurable<int> selectionFlagD{"selectionFlagD", 1, "Selection Flag for D"};
  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  HfHelper hfHelper;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber;

  double massPi{0.};
  double massD{0.};
  double massB0{0.};
  double invMass2DPiMin{0.};
  double invMass2DPiMax{0.};
  double bz{0.};

  // Fitter for B vertex (2-prong vertex filter)
  o2::vertexing::DCAFitterN<2> df2;
  // Fitter to redo D-vertex to get extrapolated daughter tracks (3-prong vertex filter)
  o2::vertexing::DCAFitterN<3> df3;

  using TracksWithSel = soa::Join<aod::TracksWCovDca, aod::TrackSelection>;
  using CandsDFiltered = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagD);

  Preslice<CandsDFiltered> candsDPerCollision = aod::track_association::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  std::shared_ptr<TH1> hCandidatesD, hCandidatesB;
  HistogramRegistry registry{"registry"};

  OutputObj<TH1F> hMassDToPiKPi{TH1F("hMassDToPiKPi", "D^{#minus} candidates;inv. mass (p^{#minus} K^{#plus} #pi^{#minus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtD{TH1F("hPtD", "D^{#minus} candidates;D^{#minus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hPtPion{TH1F("hPtPion", "#pi^{#plus} candidates;#pi^{#plus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hCPAD{TH1F("hCPAD", "D^{#minus} candidates;D^{#minus} cosine of pointing angle;entries", 110, -1.1, 1.1)};
  OutputObj<TH1F> hMassB0ToDPi{TH1F("hMassB0ToDPi", "2-prong candidates;inv. mass (B^{0} #rightarrow D^{#minus}#pi^{#plus} #rightarrow #pi^{#minus}K^{#plus}#pi^{#minus}#pi^{#plus}) (GeV/#it{c}^{2});entries", 500, 3., 8.)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};

  void init(InitContext const&)
  {
    // invariant-mass window cut
    massPi = MassPiPlus;
    massD = MassDMinus;
    massB0 = MassB0;
    invMass2DPiMin = (massB0 - invMassWindowB0) * (massB0 - invMassWindowB0);
    invMass2DPiMax = (massB0 + invMassWindowB0) * (massB0 + invMassWindowB0);

    // Initialise fitter for B vertex (2-prong vertex filter)
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);

    // Initial fitter to redo D-vertex to get extrapolated daughter tracks (3-prong vertex filter)
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

    /// candidate monitoring
    hCandidatesD = registry.add<TH1>("hCandidatesD", "D candidate counter", {HistType::kTH1D, {axisCands}});
    hCandidatesB = registry.add<TH1>("hCandidatesB", "B candidate counter", {HistType::kTH1D, {axisCands}});
    setLabelHistoCands(hCandidatesD);
    setLabelHistoCands(hCandidatesB);
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
               TracksWithSel const&,
               aod::BCsWithTimestamps const&)
  {

    static int ncol = 0;

    for (const auto& collision : collisions) {
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();

      if (ncol % 10000 == 0) {
        LOG(debug) << ncol << " collisions parsed";
      }
      ncol++;

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      df2.setBz(bz);
      df3.setBz(bz);

      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);

      for (const auto& candD : candsDThisColl) { // start loop over filtered D candidates indices as associated to this collision in candidateCreator3Prong.cxx
        hMassDToPiKPi->Fill(hfHelper.invMassDplusToPiKPi(candD), candD.pt());
        hPtD->Fill(candD.pt());
        hCPAD->Fill(candD.cpa());

        // track0 <-> pi, track1 <-> K, track2 <-> pi
        auto track0 = candD.prong0_as<TracksWithSel>();
        auto track1 = candD.prong1_as<TracksWithSel>();
        auto track2 = candD.prong2_as<TracksWithSel>();
        auto trackParCov0 = getTrackParCov(track0);
        auto trackParCov1 = getTrackParCov(track1);
        auto trackParCov2 = getTrackParCov(track2);

        std::array<float, 3> pVec0 = track0.pVector();
        std::array<float, 3> pVec1 = track1.pVector();
        std::array<float, 3> pVec2 = track2.pVector();

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
        hCandidatesD->Fill(SVFitting::BeforeFit);
        try {
          if (df3.process(trackParCov0, trackParCov1, trackParCov2) == 0) {
            continue;
          }
        } catch (const std::runtime_error& error) {
          LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN for D cannot work, skipping the candidate.";
          hCandidatesD->Fill(SVFitting::Fail);
          continue;
        }
        hCandidatesD->Fill(SVFitting::FitOk);

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
        auto trackParCovPiK = o2::dataformats::V0(df3.getPCACandidatePos(), pVecPiK, df3.calcPCACovMatrixFlat(), trackParCov0, trackParCov1);
        auto trackParCovD = o2::dataformats::V0(df3.getPCACandidatePos(), pVecD, df3.calcPCACovMatrixFlat(), trackParCovPiK, trackParCov2);

        int indexTrack0 = track0.globalIndex();
        int indexTrack1 = track1.globalIndex();
        int indexTrack2 = track2.globalIndex();

        auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);

        for (const auto& trackId : trackIdsThisCollision) { // start loop over track indices associated to this collision
          auto trackPion = trackId.track_as<TracksWithSel>();

          // check isGlobalTrackWoDCA status for pions if wanted
          if (usePionIsGlobalTrackWoDCA && !trackPion.isGlobalTrackWoDCA()) {
            continue;
          }

          // minimum pT selection
          if (trackPion.pt() < ptPionMin || !isSelectedTrackDCA(trackPion)) {
            continue;
          }
          // reject pions that are D daughters
          if (trackPion.globalIndex() == indexTrack0 || trackPion.globalIndex() == indexTrack1 || trackPion.globalIndex() == indexTrack2) {
            continue;
          }
          // reject pi and D with same sign
          if (trackPion.sign() * track0.sign() > 0) {
            continue;
          }

          hPtPion->Fill(trackPion.pt());
          std::array<float, 3> pVecPion = trackPion.pVector();
          auto trackParCovPi = getTrackParCov(trackPion);

          // ---------------------------------
          // reconstruct the 2-prong B0 vertex
          hCandidatesB->Fill(SVFitting::BeforeFit);
          try {
            if (df2.process(trackParCovD, trackParCovPi) == 0) {
              continue;
            }
          } catch (const std::runtime_error& error) {
            LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN for B cannot work, skipping the candidate.";
            hCandidatesB->Fill(SVFitting::Fail);
            continue;
          }
          hCandidatesB->Fill(SVFitting::FitOk);

          // calculate relevant properties
          const auto& secondaryVertexB0 = df2.getPCACandidate();
          auto chi2PCA = df2.getChi2AtPCACandidate();
          auto covMatrixPCA = df2.calcPCACovMatrixFlat();
          hCovSVXX->Fill(covMatrixPCA[0]);
          hCovPVXX->Fill(covMatrixPV[0]);

          // propagate D and Pi to the B0 vertex
          df2.propagateTracksToVertex();
          // track.getPxPyPzGlo(pVec) modifies pVec of track
          df2.getTrack(0).getPxPyPzGlo(pVecD);    // momentum of D at the B0 vertex
          df2.getTrack(1).getPxPyPzGlo(pVecPion); // momentum of Pi at the B0 vertex

          // calculate invariant mass square and apply selection
          auto invMass2DPi = RecoDecay::m2(std::array{pVecD, pVecPion}, std::array{massD, massPi});
          if ((invMass2DPi < invMass2DPiMin) || (invMass2DPi > invMass2DPiMax)) {
            continue;
          }
          hMassB0ToDPi->Fill(std::sqrt(invMass2DPi));

          // compute impact parameters of D and Pi
          o2::dataformats::DCA dcaD;
          o2::dataformats::DCA dcaPion;
          trackParCovD.propagateToDCA(primaryVertex, bz, &dcaD);
          trackParCovPi.propagateToDCA(primaryVertex, bz, &dcaPion);

          // get uncertainty of the decay length
          double phi, theta;
          // getPointDirection modifies phi and theta
          getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexB0, phi, theta);
          auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
          auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

          // fill the candidate table for the B0 here:
          rowCandidateBase(thisCollId,
                           collision.posX(), collision.posY(), collision.posZ(),
                           secondaryVertexB0[0], secondaryVertexB0[1], secondaryVertexB0[2],
                           errorDecayLength, errorDecayLengthXY,
                           chi2PCA,
                           pVecD[0], pVecD[1], pVecD[2],
                           pVecPion[0], pVecPion[1], pVecPion[2],
                           dcaD.getY(), dcaPion.getY(),
                           std::sqrt(dcaD.getSigmaY2()), std::sqrt(dcaPion.getSigmaY2()));

          rowCandidateProngs(candD.globalIndex(), trackPion.globalIndex());
        } // pi loop
      } // D loop
    } // collision loop
  } // process
}; // struct

/// Extends the base table with expression columns and performs MC matching.
struct HfCandidateCreatorB0Expressions {
  Spawns<aod::HfCandB0Ext> rowCandidateB0;
  Produces<aod::HfCandB0McRec> rowMcMatchRec; // table defined in CandidateReconstructionTables.h
  Produces<aod::HfCandB0McGen> rowMcMatchGen; // table defined in CandidateReconstructionTables.h

  void init(InitContext const&) {}

  void processMc(aod::HfCand3Prong const&,
                 aod::TracksWMc const&,
                 aod::McParticles const& mcParticles,
                 aod::HfCandB0Prongs const& candsB0)
  {

    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t debug = 0;

    // Match reconstructed candidates.
    for (const auto& candidate : candsB0) {
      flag = 0;
      origin = 0;
      debug = 0;
      auto candD = candidate.prong0();
      auto arrayDaughtersB0 = std::array{candD.prong0_as<aod::TracksWMc>(),
                                         candD.prong1_as<aod::TracksWMc>(),
                                         candD.prong2_as<aod::TracksWMc>(),
                                         candidate.prong1_as<aod::TracksWMc>()};
      auto arrayDaughtersD = std::array{candD.prong0_as<aod::TracksWMc>(),
                                        candD.prong1_as<aod::TracksWMc>(),
                                        candD.prong2_as<aod::TracksWMc>()};

      // B0 → D- π+ → (π- K+ π-) π+
      indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersB0, Pdg::kB0, std::array{-kPiPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 3);
      if (indexRec > -1) {
        // D- → π- K+ π-
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersD, Pdg::kDMinus, std::array{-kPiPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
        if (indexRec > -1) {
          flag = sign * BIT(hf_cand_b0::DecayTypeMc::B0ToDplusPiToPiKPiPi);
        } else {
          debug = 1;
          LOGF(debug, "WARNING: B0 in decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }
      }

      // B0 → Ds- π+ → (K- K+ π-) π+
      if (!flag) {
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersB0, Pdg::kB0, std::array{-kKPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 3);
        if (indexRec > -1) {
          // Ds- → K- K+ π-
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersD, -Pdg::kDS, std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
          if (indexRec > -1) {
            flag = sign * BIT(hf_cand_b0::DecayTypeMc::B0ToDsPiToKKPiPi);
          }
        }
      }

      // Partly reconstructed decays, i.e. the 4 prongs have a common b-hadron ancestor
      // convention: final state particles are prong0,1,2,3
      if (!flag) {
        auto particleProng0 = arrayDaughtersB0[0].mcParticle();
        auto particleProng1 = arrayDaughtersB0[1].mcParticle();
        auto particleProng2 = arrayDaughtersB0[2].mcParticle();
        auto particleProng3 = arrayDaughtersB0[3].mcParticle();
        // b-hadron hypothesis
        std::array<int, 3> bHadronMotherHypos = {Pdg::kB0, Pdg::kBS, Pdg::kLambdaB0};

        for (const auto& bHadronMotherHypo : bHadronMotherHypos) {
          int index0Mother = RecoDecay::getMother(mcParticles, particleProng0, bHadronMotherHypo, true);
          int index1Mother = RecoDecay::getMother(mcParticles, particleProng1, bHadronMotherHypo, true);
          int index2Mother = RecoDecay::getMother(mcParticles, particleProng2, bHadronMotherHypo, true);
          int index3Mother = RecoDecay::getMother(mcParticles, particleProng3, bHadronMotherHypo, true);

          // look for common b-hadron ancestor
          if (index0Mother > -1 && index1Mother > -1 && index2Mother > -1 && index3Mother > -1) {
            if (index0Mother == index1Mother && index1Mother == index2Mother && index2Mother == index3Mother) {
              flag = BIT(hf_cand_b0::DecayTypeMc::PartlyRecoDecay);
              break;
            }
          }
        }
      }

      rowMcMatchRec(flag, origin, debug);
    } // rec

    hf_mc_gen::fillMcMatchGenB0(mcParticles, rowMcMatchGen); // gen
  } // processMc
  PROCESS_SWITCH(HfCandidateCreatorB0Expressions, processMc, "Process MC", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfCandidateCreatorB0>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorB0Expressions>(cfgc)};

  return workflow;
}
