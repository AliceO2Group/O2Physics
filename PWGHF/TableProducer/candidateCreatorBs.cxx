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

/// \file candidateCreatorBs.cxx
/// \brief Reconstruction of Bs candidates
/// \note Adapted from candidateCreatorB0.cxx
///
/// \author Phil Stahlhut <phil.lennart.stahlhut@cern.ch>

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/V0.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
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
using namespace o2::hf_decay::hf_cand_beauty;

/// Reconstruction of Bs candidates
struct HfCandidateCreatorBs {
  Produces<aod::HfCandBsBase> rowCandidateBase;     // table defined in CandidateReconstructionTables.h
  Produces<aod::HfCandBsProngs> rowCandidateProngs; // table defined in CandidateReconstructionTables.h

  // vertexing
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCADs{"useAbsDCADs", false, "Minimise abs. distance rather than chi2 for 3-prong Ds vertex"};
  Configurable<bool> useAbsDCABs{"useAbsDCABs", true, "Minimise abs. distance rather than chi2 for 2-prong Bs vertex"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 4., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any Bs is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // selection
  Configurable<bool> usePionIsGlobalTrackWoDCA{"usePionIsGlobalTrackWoDCA", true, "check isGlobalTrackWoDCA status for pions, for Run3 studies"};
  Configurable<double> ptPionMin{"ptPionMin", 0.5, "minimum pion pT threshold (GeV/c)"};
  Configurable<std::vector<double>> binsPtPion{"binsPtPion", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for pion DCA XY pT-dependent cut"};
  Configurable<LabeledArray<double>> cutsTrackPionDCA{"cutsTrackPionDCA", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for pions"};
  Configurable<double> invMassWindowBs{"invMassWindowBs", 0.3, "invariant-mass window for Bs candidates"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 1, "Selection Flag for Ds"};
  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  o2::vertexing::DCAFitterN<2> df2; // 2-prong vertex fitter
  o2::vertexing::DCAFitterN<3> df3; // 3-prong vertex fitter
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  int runNumber{0};
  double massDsPi{0.};
  double bz{0.};

  using TracksWithSel = soa::Join<aod::TracksWCovDca, aod::TrackSelection>;
  using CandsDsFiltered = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs);

  Preslice<CandsDsFiltered> candsDsPerCollision = aod::track_association::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  std::shared_ptr<TH1> hCandidatesD, hCandidatesB;
  HistogramRegistry registry{"registry"};

  OutputObj<TH1F> hMassDsToKKPi{TH1F("hMassDsToKKPi", "D_{s} candidates;inv. mass (K K #pi) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtDs{TH1F("hPtDs", "D_{s} candidates;D_{s} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hPtPion{TH1F("hPtPion", "#pi candidates;#pi candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hCPADs{TH1F("hCPADs", "D_{s} candidates;D_{s} cosine of pointing angle;entries", 110, -1.1, 1.1)};
  OutputObj<TH1F> hMassBsToDsPi{TH1F("hMassBsToDsPi", "2-prong candidates;inv. mass (B_{s} #rightarrow D_{s}#pi #rightarrow KK#pi#pi) (GeV/#it{c}^{2});entries", 500, 3., 8.)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};

  void init(InitContext const&)
  {
    // Initialise fitter for Bs vertex (2-prong vertex fitter)
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCABs);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);

    // Initialise fitter to redo Ds-vertex to get extrapolated daughter tracks (3-prong vertex fitter)
    df3.setPropagateToPCA(propagateToPCA);
    df3.setMaxR(maxR);
    df3.setMaxDZIni(maxDZIni);
    df3.setMinParamChange(minParamChange);
    df3.setMinRelChi2Change(minRelChi2Change);
    df3.setUseAbsDCA(useAbsDCADs);
    df3.setWeightedFinalPCA(useWeightedFinalPCA);

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
               CandsDsFiltered const& candsDs,
               aod::TrackAssoc const& trackIndices,
               TracksWithSel const&,
               aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();

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
      auto candsDsThisColl = candsDs.sliceBy(candsDsPerCollision, thisCollId);

      for (const auto& candDs : candsDsThisColl) { // start loop over filtered Ds candidates indices as associated to this collision in candidateCreator3Prong.cxx

        // track0 <-> K, track1 <-> K, track2 <-> pi
        auto track0 = candDs.prong0_as<TracksWithSel>();
        auto track1 = candDs.prong1_as<TracksWithSel>();
        auto track2 = candDs.prong2_as<TracksWithSel>();
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
        // reconstruct 3-prong secondary vertex (Ds±)
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

        const auto& secondaryVertexDs = df3.getPCACandidate();
        // propagate the 3 prongs to the secondary vertex
        trackParCov0.propagateTo(secondaryVertexDs[0], bz);
        trackParCov1.propagateTo(secondaryVertexDs[0], bz);
        trackParCov2.propagateTo(secondaryVertexDs[0], bz);

        // update pVec of tracks
        df3.getTrack(0).getPxPyPzGlo(pVec0);
        df3.getTrack(1).getPxPyPzGlo(pVec1);
        df3.getTrack(2).getPxPyPzGlo(pVec2);

        // Ds∓ → K∓ K± π∓
        std::array<float, 3> const pVecKK = RecoDecay::pVec(pVec0, pVec1);
        std::array<float, 3> pVecDs = RecoDecay::pVec(pVec0, pVec1, pVec2);
        auto trackParCovKK = o2::dataformats::V0(df3.getPCACandidatePos(), pVecKK, df3.calcPCACovMatrixFlat(),
                                                 trackParCov0, trackParCov1);
        auto trackParCovDs = o2::dataformats::V0(df3.getPCACandidatePos(), pVecDs, df3.calcPCACovMatrixFlat(),
                                                 trackParCovKK, trackParCov2);

        int const indexTrack0 = track0.globalIndex();
        int const indexTrack1 = track1.globalIndex();
        int const indexTrack2 = track2.globalIndex();

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
          // reject pions that are Ds daughters
          if (trackPion.globalIndex() == indexTrack0 || trackPion.globalIndex() == indexTrack1 || trackPion.globalIndex() == indexTrack2) {
            continue;
          }
          // reject pi and Ds with same sign
          if (trackPion.sign() * track0.sign() > 0) {
            continue;
          }

          std::array<float, 3> pVecPion = trackPion.pVector();
          auto trackParCovPi = getTrackParCov(trackPion);

          // ---------------------------------
          // reconstruct the 2-prong Bs vertex
          hCandidatesB->Fill(SVFitting::BeforeFit);
          try {
            if (df2.process(trackParCovDs, trackParCovPi) == 0) {
              continue;
            }
          } catch (const std::runtime_error& error) {
            LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN for B cannot work, skipping the candidate.";
            hCandidatesB->Fill(SVFitting::Fail);
            continue;
          }
          hCandidatesB->Fill(SVFitting::FitOk);

          // calculate relevant properties
          const auto& secondaryVertexBs = df2.getPCACandidate();
          auto chi2PCA = df2.getChi2AtPCACandidate();
          auto covMatrixPCA = df2.calcPCACovMatrixFlat();

          // propagate Ds and Pi to the Bs vertex
          df2.propagateTracksToVertex();
          // track.getPxPyPzGlo(pVec) modifies pVec of track
          df2.getTrack(0).getPxPyPzGlo(pVecDs);   // momentum of Ds at the Bs vertex
          df2.getTrack(1).getPxPyPzGlo(pVecPion); // momentum of Pi at the Bs vertex

          // calculate invariant mass and apply selection
          massDsPi = RecoDecay::m(std::array{pVecDs, pVecPion}, std::array{MassDSBar, MassPiPlus});
          if (std::abs(massDsPi - MassBS) > invMassWindowBs) {
            continue;
          }

          // compute impact parameters of Ds and Pi
          o2::dataformats::DCA dcaDs;
          o2::dataformats::DCA dcaPion;
          trackParCovDs.propagateToDCA(primaryVertex, bz, &dcaDs);
          trackParCovPi.propagateToDCA(primaryVertex, bz, &dcaPion);

          // get uncertainty of the decay length
          double phi{}, theta{};
          // getPointDirection modifies phi and theta
          getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexBs, phi, theta);
          auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
          auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

          // fill output histograms for Bs candidates
          hMassDsToKKPi->Fill(HfHelper::invMassDsToKKPi(candDs), candDs.pt());
          hCovSVXX->Fill(covMatrixPCA[0]);
          hCovPVXX->Fill(covMatrixPV[0]);
          hMassBsToDsPi->Fill(massDsPi);
          hPtDs->Fill(candDs.pt());
          hCPADs->Fill(candDs.cpa());
          hPtPion->Fill(trackPion.pt());

          // fill the candidate table for the Bs here:
          rowCandidateBase(thisCollId,
                           collision.posX(), collision.posY(), collision.posZ(),
                           secondaryVertexBs[0], secondaryVertexBs[1], secondaryVertexBs[2],
                           errorDecayLength, errorDecayLengthXY,
                           chi2PCA,
                           pVecDs[0], pVecDs[1], pVecDs[2],
                           pVecPion[0], pVecPion[1], pVecPion[2],
                           dcaDs.getY(), dcaPion.getY(),
                           std::sqrt(dcaDs.getSigmaY2()), std::sqrt(dcaPion.getSigmaY2()));

          rowCandidateProngs(candDs.globalIndex(), trackPion.globalIndex());
        } // pi loop
      } // Ds loop
    } // collision loop
  } // process
}; // struct

/// Extends the base table with expression columns and performs MC matching.
struct HfCandidateCreatorBsExpressions {
  Spawns<aod::HfCandBsExt> rowCandidateBs;
  Produces<aod::HfCandBsMcRec> rowMcMatchRec; // table defined in CandidateReconstructionTables.h
  Produces<aod::HfCandBsMcGen> rowMcMatchGen; // table defined in CandidateReconstructionTables.h

  void init(InitContext const&) {}

  void processMc(aod::HfCand3Prong const&,
                 aod::TracksWMc const&,
                 aod::McParticles const& mcParticles,
                 aod::HfCandBsProngs const& candsBs)
  {
    int indexRec = -1;
    int8_t sign = 0;
    int8_t flagChannelMain = 0;
    int8_t flagChannelReso = 0;
    std::vector<int> arrDaughDsIndex{};
    std::array<int, 2> arrPDGDaughDs{};
    std::array<int, 2> arrPDGResonantDsPhiPi = {Pdg::kPhi, kPiPlus}; // Ds± → Phi π±

    // Match reconstructed candidates.
    for (const auto& candidate : candsBs) {
      flagChannelMain = 0;
      flagChannelReso = 0;
      arrDaughDsIndex.clear();
      auto candDs = candidate.prong0();
      auto arrayDaughtersBs = std::array{candDs.prong0_as<aod::TracksWMc>(),
                                         candDs.prong1_as<aod::TracksWMc>(),
                                         candDs.prong2_as<aod::TracksWMc>(),
                                         candidate.prong1_as<aod::TracksWMc>()};
      auto arrayDaughtersDs = std::array{candDs.prong0_as<aod::TracksWMc>(),
                                         candDs.prong1_as<aod::TracksWMc>(),
                                         candDs.prong2_as<aod::TracksWMc>()};

      // Checking Bs0(bar) → Ds∓ π± → (K- K+ π∓) π±
      indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersBs, Pdg::kBS, std::array{-kKPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 3);
      if (indexRec > -1) {
        // Checking Ds∓ → K- K+ π∓
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersDs, Pdg::kDSBar, std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2);
        if (indexRec > -1) {
          RecoDecay::getDaughters(mcParticles.rawIteratorAt(indexRec), &arrDaughDsIndex, std::array{0}, 1);
          if (arrDaughDsIndex.size() == 2) {
            for (auto iProng = 0u; iProng < arrDaughDsIndex.size(); ++iProng) {
              auto daughI = mcParticles.rawIteratorAt(arrDaughDsIndex[iProng]);
              arrPDGDaughDs[iProng] = std::abs(daughI.pdgCode());
            }
            if ((arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[0] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[1]) || (arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[1] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[0])) {
              flagChannelMain = sign * DecayChannelMain::BsToDsPi;
            }
          }
        }
      }

      if (flagChannelMain == 0) {
        // Checking B0(bar) → Ds± π∓ → (K- K+ π±) π∓
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersBs, Pdg::kB0, std::array{-kKPlus, +kKPlus, +kPiPlus, -kPiPlus}, true, &sign, 3);
        if (indexRec > -1) {
          // Checking Ds± → K- K+ π±
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersDs, Pdg::kDS, std::array{-kKPlus, +kKPlus, +kPiPlus}, true, &sign, 2);
          if (indexRec > -1) {
            RecoDecay::getDaughters(mcParticles.rawIteratorAt(indexRec), &arrDaughDsIndex, std::array{0}, 1);
            if (arrDaughDsIndex.size() == 2) {
              for (auto iProng = 0u; iProng < arrDaughDsIndex.size(); ++iProng) {
                auto daughI = mcParticles.rawIteratorAt(arrDaughDsIndex[iProng]);
                arrPDGDaughDs[iProng] = std::abs(daughI.pdgCode());
              }
              if ((arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[0] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[1]) || (arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[1] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[0])) {
                flagChannelMain = sign * DecayChannelMain::B0ToDsPi;
              }
            }
          }
        }
      }

      // Partly reconstructed decays, i.e. the 4 prongs have a common b-hadron ancestor
      // convention: final state particles are prong0,1,2,3
      if (flagChannelMain == 0) {
        auto particleProng0 = arrayDaughtersBs[0].mcParticle();
        auto particleProng1 = arrayDaughtersBs[1].mcParticle();
        auto particleProng2 = arrayDaughtersBs[2].mcParticle();
        auto particleProng3 = arrayDaughtersBs[3].mcParticle();
        // b-hadron hypothesis
        std::array<int, 4> const bHadronMotherHypos = {Pdg::kB0, Pdg::kBPlus, Pdg::kBS, Pdg::kLambdaB0};

        for (const auto& bHadronMotherHypo : bHadronMotherHypos) {
          int const index0Mother = RecoDecay::getMother(mcParticles, particleProng0, bHadronMotherHypo, true);
          int const index1Mother = RecoDecay::getMother(mcParticles, particleProng1, bHadronMotherHypo, true);
          int const index2Mother = RecoDecay::getMother(mcParticles, particleProng2, bHadronMotherHypo, true);
          int const index3Mother = RecoDecay::getMother(mcParticles, particleProng3, bHadronMotherHypo, true);

          // look for common b-hadron ancestor
          if (index0Mother > -1 && index0Mother == index1Mother && index1Mother == index2Mother && index2Mother == index3Mother) {
            flagChannelMain = hf_cand_bs::DecayTypeMc::PartlyRecoDecay; // FIXME
            break;
          }
        }
      }

      rowMcMatchRec(flagChannelMain, flagChannelReso);
    } // rec

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      flagChannelMain = 0;
      flagChannelReso = 0;
      arrDaughDsIndex.clear();

      // Checking Bs0(bar) → Ds∓ π± → (K- K+ π∓) π±
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kBS, std::array{+Pdg::kDSBar, +kPiPlus}, true)) {
        // Checking Ds∓ → K- K+ π∓
        auto candDsMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
        if (RecoDecay::isMatchedMCGen(mcParticles, candDsMC, Pdg::kDSBar, std::array{-kKPlus, +kKPlus, -kPiPlus}, true, &sign, 2)) {
          RecoDecay::getDaughters(candDsMC, &arrDaughDsIndex, std::array{0}, 1);
          if (arrDaughDsIndex.size() == 2) {
            for (auto jProng = 0u; jProng < arrDaughDsIndex.size(); ++jProng) {
              auto daughJ = mcParticles.rawIteratorAt(arrDaughDsIndex[jProng]);
              arrPDGDaughDs[jProng] = std::abs(daughJ.pdgCode());
            }
            if ((arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[0] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[1]) || (arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[1] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[0])) {
              flagChannelMain = sign * DecayChannelMain::BsToDsPi;
            }
          }
        }
      }

      if (flagChannelMain == 0) {
        // Checking B0(bar) → Ds± π∓ → (K- K+ π±) π∓
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kB0, std::array{+Pdg::kDS, -kPiPlus}, true)) {
          // Checking Ds± → K- K+ π±
          auto candDsMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(mcParticles, candDsMC, Pdg::kDS, std::array{-kKPlus, +kKPlus, +kPiPlus}, true, &sign, 2)) {
            RecoDecay::getDaughters(candDsMC, &arrDaughDsIndex, std::array{0}, 1);
            if (arrDaughDsIndex.size() == 2) {
              for (auto jProng = 0u; jProng < arrDaughDsIndex.size(); ++jProng) {
                auto daughJ = mcParticles.rawIteratorAt(arrDaughDsIndex[jProng]);
                arrPDGDaughDs[jProng] = std::abs(daughJ.pdgCode());
              }
              if ((arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[0] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[1]) || (arrPDGDaughDs[0] == arrPDGResonantDsPhiPi[1] && arrPDGDaughDs[1] == arrPDGResonantDsPhiPi[0])) {
                flagChannelMain = sign * DecayChannelMain::B0ToDsPi;
              }
            }
          }
        }
      }

      rowMcMatchGen(flagChannelMain, flagChannelReso);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfCandidateCreatorBsExpressions, processMc, "Process MC", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfCandidateCreatorBs>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorBsExpressions>(cfgc)};

  return workflow;
}
