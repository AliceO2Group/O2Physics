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

/// \file candidateCreatorBplus.cxx
/// \brief Reconstruction of B± → D0bar(D0) π± candidates
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN
/// \author Deepa Thomas <deepa.thomas@cern.ch>, UT Austin
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari & INFN, Sezione di Bari

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsMcGen.h"
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

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::hf_trkcandsel;
using namespace o2::aod;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_decay::hf_cand_beauty;

/// Reconstruction of B± → D0bar(D0) π± → (K± π∓) π±
struct HfCandidateCreatorBplus {
  Produces<aod::HfCandBplusBase> rowCandidateBase;     // table defined in CandidateReconstructionTables.h
  Produces<aod::HfCandBplusProngs> rowCandidateProngs; // table defined in CandidateReconstructionTables.h

  // vertexing parameters
  // Configurable<double> bz{"bz", 5., "magnetic field"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 5., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 999, "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
  // selection
  Configurable<bool> usePionIsGlobalTrackWoDCA{"usePionIsGlobalTrackWoDCA", true, "check isGlobalTrackWoDCA status for pions, for Run3 studies"};
  Configurable<double> ptPionMin{"ptPionMin", 0.2, "minimum pion pT threshold (GeV/c)"};
  Configurable<std::vector<double>> binsPtPion{"binsPtPion", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for pion DCA XY pT-dependent cut"};
  Configurable<LabeledArray<double>> cutsTrackPionDCA{"cutsTrackPionDCA", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for pions"};
  Configurable<double> invMassWindowBplus{"invMassWindowBplus", 0.3, "invariant-mass window for B^{+} candidates"};
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<double> etaTrackMax{"etaTrackMax", -1, "max. bach track. pseudorapidity"};
  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber{};

  double invMass2D0PiMin{0.};
  double invMass2D0PiMax{0.};
  double bz{0.};

  // Fitter for B vertex
  o2::vertexing::DCAFitterN<2> dfB;
  // Fitter to redo D-vertex to get extrapolated daughter tracks
  o2::vertexing::DCAFitterN<2> df;

  using TracksWithSel = soa::Join<aod::TracksWCovDca, aod::TrackSelection>;
  using CandsDFiltered = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar);
  Preslice<CandsDFiltered> candsDPerCollision = aod::track_association::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  std::shared_ptr<TH1> hCandidatesD, hCandidatesB;
  HistogramRegistry registry{"registry"};

  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH1F> hRapidityD0{TH1F("hRapidityD0", "D0 candidates;#it{y};entries", 100, -2, 2)};
  OutputObj<TH1F> hEtaPi{TH1F("hEtaPi", "Pion track;#it{#eta};entries", 400, -2., 2.)};
  OutputObj<TH1F> hMassBplusToD0Pi{TH1F("hMassBplusToD0Pi", "2-prong candidates;inv. mass (B^{+} #rightarrow #bar{D^{0}}#pi^{+}) (GeV/#it{c}^{2});entries", 500, 3., 8.)};

  void init(InitContext const&)
  {
    // invariant-mass window cut
    invMass2D0PiMin = (MassBPlus - invMassWindowBplus) * (MassBPlus - invMassWindowBplus);
    invMass2D0PiMax = (MassBPlus + invMassWindowBplus) * (MassBPlus + invMassWindowBplus);

    // Initialise fitter for B vertex
    dfB.setPropagateToPCA(propagateToPCA);
    dfB.setMaxR(maxR);
    dfB.setMinParamChange(minParamChange);
    dfB.setMinRelChi2Change(minRelChi2Change);
    dfB.setUseAbsDCA(true);
    dfB.setWeightedFinalPCA(useWeightedFinalPCA);

    // Initial fitter to redo D-vertex to get extrapolated daughter tracks
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

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
  bool isSelectedTrack(const T& track)
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

    for (const auto& collision : collisions) {
      auto primaryVertex = getPrimaryVertex(collision);
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
      df.setBz(bz);
      dfB.setBz(bz);

      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);

      // loop over pairs of track indices
      for (const auto& candD0 : candsDThisColl) {

        if (!TESTBIT(candD0.hfflag(), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          continue;
        }
        if (yCandMax >= 0. && std::abs(HfHelper::yD0(candD0)) > yCandMax) {
          continue;
        }

        hRapidityD0->Fill(HfHelper::yD0(candD0));

        // track0 <-> pi, track1 <-> K
        auto prong0 = candD0.prong0_as<TracksWithSel>();
        auto prong1 = candD0.prong1_as<TracksWithSel>();
        auto trackParCovProng0 = getTrackParCov(prong0);
        auto trackParCovProng1 = getTrackParCov(prong1);

        auto dca0 = o2::dataformats::DCA(prong0.dcaXY(), prong0.dcaZ(), prong0.cYY(), prong0.cZY(), prong0.cZZ());
        auto dca1 = o2::dataformats::DCA(prong1.dcaXY(), prong1.dcaZ(), prong1.cYY(), prong1.cZY(), prong1.cZZ());

        // repropagate tracks to this collision if needed
        if (prong0.collisionId() != thisCollId) {
          trackParCovProng0.propagateToDCA(primaryVertex, bz, &dca0);
        }

        if (prong1.collisionId() != thisCollId) {
          trackParCovProng1.propagateToDCA(primaryVertex, bz, &dca1);
        }

        // reconstruct D0 secondary vertex
        hCandidatesD->Fill(SVFitting::BeforeFit);
        try {
          if (df.process(trackParCovProng0, trackParCovProng1) == 0) {
            continue;
          }
        } catch (const std::runtime_error& error) {
          LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN for D cannot work, skipping the candidate.";
          hCandidatesD->Fill(SVFitting::Fail);
          continue;
        }
        hCandidatesD->Fill(SVFitting::FitOk);

        const auto& vertexD0 = df.getPCACandidatePos();
        trackParCovProng0.propagateTo(vertexD0[0], bz);
        trackParCovProng1.propagateTo(vertexD0[0], bz);
        // Get pVec of tracks
        std::array<float, 3> pVec0 = {0};
        std::array<float, 3> pVec1 = {0};
        df.getTrack(0).getPxPyPzGlo(pVec0);
        df.getTrack(1).getPxPyPzGlo(pVec1);
        // Get D0 momentum
        std::array<float, 3> const pVecD = RecoDecay::pVec(pVec0, pVec1);

        // build a D0 neutral track
        auto trackD0 = o2::dataformats::V0(vertexD0, pVecD, df.calcPCACovMatrixFlat(), trackParCovProng0, trackParCovProng1);

        int const indexTrack0 = prong0.globalIndex();
        int const indexTrack1 = prong1.globalIndex();

        auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);

        // loop over tracks pi
        for (const auto& trackId : trackIdsThisCollision) { // start loop over track indices associated to this collision
          auto trackPion = trackId.track_as<TracksWithSel>();

          // check isGlobalTrackWoDCA status for pions if wanted
          if (usePionIsGlobalTrackWoDCA && !trackPion.isGlobalTrackWoDCA()) {
            continue;
          }

          // Select D0pi- and D0(bar)pi+ pairs only
          if ((candD0.isSelD0() < selectionFlagD0 || trackPion.sign() >= 0) && (candD0.isSelD0bar() < selectionFlagD0bar || trackPion.sign() <= 0)) {
            // LOGF(debug, "D0: %d, D0bar%d, sign: %d", candD0.isSelD0(), candD0.isSelD0bar(), track.sign());
            continue;
          }

          // minimum pT selection
          if (trackPion.pt() < ptPionMin || !isSelectedTrack(trackPion)) {
            continue;
          }

          if (indexTrack0 == trackPion.globalIndex() || indexTrack1 == trackPion.globalIndex()) {
            continue; // different id between D0 daughters and bachelor track
          }

          if (etaTrackMax >= 0. && std::abs(trackPion.eta()) > etaTrackMax) {
            continue;
          }

          hEtaPi->Fill(trackPion.eta());

          auto trackParCovPi = getTrackParCov(trackPion);
          std::array<float, 3> pVecD0 = {0., 0., 0.};
          std::array<float, 3> pVecBach = {0., 0., 0.};

          // find the DCA between the D0 and the bachelor track, for B+
          hCandidatesB->Fill(SVFitting::BeforeFit);
          try {
            if (dfB.process(trackD0, trackParCovPi) == 0) {
              continue;
            }
          } catch (const std::runtime_error& error) {
            LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN for B cannot work, skipping the candidate.";
            hCandidatesB->Fill(SVFitting::Fail);
            continue;
          }
          hCandidatesB->Fill(SVFitting::FitOk);

          // get D and Pi tracks (propagated to the B+ vertex if propagateToPCA==true)
          trackD0.getPxPyPzGlo(pVecD0);         // momentum of D0 at the B+ vertex
          trackParCovPi.getPxPyPzGlo(pVecBach); // momentum of pi+ at the B+ vertex

          const auto& secVertexBplus = dfB.getPCACandidate();
          auto chi2PCA = dfB.getChi2AtPCACandidate();
          auto covMatrixPCA = dfB.calcPCACovMatrixFlat();
          hCovSVXX->Fill(covMatrixPCA[0]); // FIXME: Calculation of errorDecayLength(XY) gives wrong values without this line.

          // get track impact parameters
          // This modifies track momenta!
          auto covMatrixPV = primaryVertex.getCov();
          hCovPVXX->Fill(covMatrixPV[0]);
          o2::dataformats::DCA impactParameter0;
          o2::dataformats::DCA impactParameter1;
          trackD0.propagateToDCA(primaryVertex, bz, &impactParameter0);
          trackParCovPi.propagateToDCA(primaryVertex, bz, &impactParameter1);

          // get uncertainty of the decay length
          double phi, theta;
          getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secVertexBplus, phi, theta);
          auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
          auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

          // compute invariant mass square and apply selection
          auto invMass2D0Pi = RecoDecay::m2(std::array{pVecD0, pVecBach}, std::array{MassD0, MassPiPlus});
          if ((invMass2D0Pi < invMass2D0PiMin) || (invMass2D0Pi > invMass2D0PiMax)) {
            continue;
          }
          hMassBplusToD0Pi->Fill(std::sqrt(invMass2D0Pi));

          // fill candidate table rows
          rowCandidateBase(collision.globalIndex(),
                           collision.posX(), collision.posY(), collision.posZ(),
                           secVertexBplus[0], secVertexBplus[1], secVertexBplus[2],
                           errorDecayLength, errorDecayLengthXY,
                           chi2PCA,
                           pVecD0[0], pVecD0[1], pVecD0[2],
                           pVecBach[0], pVecBach[1], pVecBach[2],
                           impactParameter0.getY(), impactParameter1.getY(),
                           std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()));

          rowCandidateProngs(candD0.globalIndex(), trackPion.globalIndex()); // index D0 and bachelor
        } // track loop
      } // D0 cand loop
    } // collision
  } // process
}; // struct

/// Extends the base table with expression columns and performs MC matching
struct HfCandidateCreatorBplusExpressions {
  Spawns<aod::HfCandBplusExt> rowCandidateBPlus;
  Produces<aod::HfCandBplusMcRec> rowMcMatchRec;
  Produces<aod::HfCandBplusMcGen> rowMcMatchGen;

  void init(InitContext const&) {}

  void processMc(aod::HfCand2Prong const&,
                 aod::TracksWMc const&,
                 aod::McParticles const& mcParticles,
                 aod::HfCandBplusProngs const& candsBplus)
  {

    int indexRec = -1, indexRecD0 = -1;
    int8_t signB = 0, signD0 = 0;
    int8_t flagChannelMain = 0;
    int8_t flagChannelReso = 0;
    int8_t origin = 0;

    // Match reconstructed candidates.
    // Spawned table can be used directly
    for (const auto& candidate : candsBplus) {

      flagChannelMain = 0;
      flagChannelReso = 0;
      origin = 0;
      auto candDaughterD0 = candidate.prong0();
      auto arrayDaughtersD0 = std::array{candDaughterD0.prong0_as<aod::TracksWMc>(), candDaughterD0.prong1_as<aod::TracksWMc>()};
      auto arrayDaughters = std::array{candidate.prong1_as<aod::TracksWMc>(), candDaughterD0.prong0_as<aod::TracksWMc>(), candDaughterD0.prong1_as<aod::TracksWMc>()};

      // B± → D0bar(D0) π± → (K± π∓) π±
      indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kBPlus, std::array{+kPiPlus, +kKPlus, -kPiPlus}, true, &signB, 2);
      indexRecD0 = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersD0, -Pdg::kD0, std::array{+kKPlus, -kPiPlus}, true, &signD0, 1);

      if (indexRecD0 > -1 && indexRec > -1) {
        flagChannelMain = signB * DecayChannelMain::BplusToD0Pi;
      }
      rowMcMatchRec(flagChannelMain, flagChannelReso, origin);
    }
    hf_mc_gen::fillMcMatchGenBplus(mcParticles, rowMcMatchGen); // gen
  } // process
  PROCESS_SWITCH(HfCandidateCreatorBplusExpressions, processMc, "Process MC", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorBplus>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorBplusExpressions>(cfgc)};
}
