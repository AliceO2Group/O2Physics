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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "DCAFitter/DCAFitterN.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "ReconstructionDataFormats/V0.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::aod::hf_cand_bplus;

/// Reconstruction of B± → D0bar(D0) π± → (K± π∓) π±
struct HfCandidateCreatorBplus {
  Produces<aod::HfCandBplusBase> rowCandidateBase;

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
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<double> etaTrackMax{"etaTrackMax", -1, "max. bach track. pseudorapidity"};
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
  double massD0 = RecoDecay::getMassPDG(pdg::Code::kDMinus);
  double massD0Pi = 0.;
  double bz = 0.;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar);

  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH1F> hNEvents{TH1F("hNEvents", "Number of events;Nevents;entries", 1, 0., 1)};
  OutputObj<TH1F> hRapidityD0{TH1F("hRapidityD0", "D0 candidates;#it{y};entries", 100, -2, 2)};
  OutputObj<TH1F> hEtaPi{TH1F("hEtaPi", "Pion track;#it{#eta};entries", 400, 2, 2)};
  OutputObj<TH1F> hMassBplusToD0Pi{TH1F("hMassBplusToD0Pi", "2-prong candidates;inv. mass (B^{+} #rightarrow #bar{D^{0}}#pi^{+}) (GeV/#it{c}^{2});entries", 500, 3., 8.)};

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(ccdbPathGeo);
    }
    runNumber = 0;
  }

  void process(aod::Collision const& collisions,
               soa::Filtered<soa::Join<aod::HfCand2Prong,
                                       aod::HfSelD0>> const& candidates,
               aod::BigTracks const& tracks,
               aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    hNEvents->Fill(0);

    // Initialise fitter for B vertex
    o2::vertexing::DCAFitterN<2> bfitter;
    // bfitter.setBz(bz);
    bfitter.setPropagateToPCA(propagateToPCA);
    bfitter.setMaxR(maxR);
    bfitter.setMinParamChange(minParamChange);
    bfitter.setMinRelChi2Change(minRelChi2Change);
    bfitter.setUseAbsDCA(useAbsDCA);
    bfitter.setWeightedFinalPCA(useWeightedFinalPCA);

    // Initial fitter to redo D-vertex to get extrapolated daughter tracks
    o2::vertexing::DCAFitterN<2> df;
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

    // loop over pairs of track indices
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(yD0(candidate)) > yCandMax) {
        continue;
      }

      hRapidityD0->Fill(yD0(candidate));

      const std::array<float, 3> vertexD0 = {candidate.xSecondaryVertex(), candidate.ySecondaryVertex(), candidate.zSecondaryVertex()};
      const std::array<float, 3> momentumD0 = {candidate.px(), candidate.py(), candidate.pz()};

      auto prong0 = candidate.prong0_as<aod::BigTracks>();
      auto prong1 = candidate.prong1_as<aod::BigTracks>();
      auto prong0TrackParCov = getTrackParCov(prong0);
      auto prong1TrackParCov = getTrackParCov(prong1);
      auto collision = prong0.collision();

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = prong0.collision().bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      df.setBz(bz);

      // LOGF(info, "All track: %d (prong0); %d (prong1)", candidate.prong0().globalIndex(), candidate.prong1().globalIndex());
      // LOGF(info, "All track pT: %f (prong0); %f (prong1)", prong0.pt(), prong1.pt());

      // reconstruct D0 secondary vertex
      if (df.process(prong0TrackParCov, prong1TrackParCov) == 0) {
        continue;
      }

      prong0TrackParCov.propagateTo(candidate.xSecondaryVertex(), bz);
      prong1TrackParCov.propagateTo(candidate.xSecondaryVertex(), bz);
      // std::cout << df.getTrack(0).getX() << " "<< "secVx=" << candidate.xSecondaryVertex() << std::endl;

      const std::array<float, 6> pCovMatrixD0 = df.calcPCACovMatrixFlat();
      // build a D0 neutral track
      auto trackD0 = o2::dataformats::V0(vertexD0, momentumD0, pCovMatrixD0, prong0TrackParCov, prong1TrackParCov, {0, 0}, {0, 0});

      // loop over tracks for pi selection
      // auto count = 0;
      for (auto& track : tracks) {
        // if(count % 100 == 0){
        // LOGF(info, "Col: %d (cand); %d (track)", candidate.collisionId(), track.collisionId());
        //   count++;
        //  }
        bfitter.setBz(bz);

        if (etaTrackMax >= 0. && std::abs(track.eta()) > etaTrackMax) {
          continue;
        }

        hEtaPi->Fill(track.eta());

        if (candidate.prong0Id() == track.globalIndex() || candidate.prong1Id() == track.globalIndex()) {
          continue; // daughter track id and bachelor track id not the same
        }

        // Select D0pi- and D0(bar)pi+ pairs only
        if (!((candidate.isSelD0() >= selectionFlagD0 && track.sign() < 0) || (candidate.isSelD0bar() >= selectionFlagD0bar && track.sign() > 0))) {
          // Printf("D0: %d, D0bar%d, sign: %d", candidate.isSelD0(), candidate.isSelD0bar(), track.sign());
          continue;
        }

        auto trackBach = getTrackParCov(track);
        std::array<float, 3> pVecD0 = {0., 0., 0.};
        std::array<float, 3> pVecBach = {0., 0., 0.};
        std::array<float, 3> pVecBCand = {0., 0., 0.};

        // find the DCA between the D0 and the bachelor track, for B+
        if (bfitter.process(trackD0, trackBach) == 0) {
          continue;
        }

        bfitter.propagateTracksToVertex();          // propagate the bachelor and D0 to the B+ vertex
        bfitter.getTrack(0).getPxPyPzGlo(pVecD0);   // momentum of D0 at the B+ vertex
        bfitter.getTrack(1).getPxPyPzGlo(pVecBach); // momentum of pi+ at the B+ vertex
        const auto& BSecVertex = bfitter.getPCACandidate();
        auto chi2PCA = bfitter.getChi2AtPCACandidate();
        auto covMatrixPCA = bfitter.calcPCACovMatrixFlat();
        hCovSVXX->Fill(covMatrixPCA[0]); // FIXME: Calculation of errorDecayLength(XY) gives wrong values without this line.

        pVecBCand = RecoDecay::pVec(pVecD0, pVecBach);

        // get track impact parameters
        // This modifies track momenta!
        auto primaryVertex = getPrimaryVertex(collision);
        auto covMatrixPV = primaryVertex.getCov();
        hCovPVXX->Fill(covMatrixPV[0]);
        o2::dataformats::DCA impactParameter0;
        o2::dataformats::DCA impactParameter1;

        bfitter.getTrack(0).propagateToDCA(primaryVertex, bz, &impactParameter0);
        bfitter.getTrack(1).propagateToDCA(primaryVertex, bz, &impactParameter1);

        // get uncertainty of the decay length
        double phi, theta;
        getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, BSecVertex, phi, theta);
        auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        int hfFlag = 1 << hf_cand_bplus::DecayType::BplusToD0Pi;

        // calculate invariant mass and fill the Invariant Mass control plot
        auto arrayMomenta = array{pVecD0, pVecBach};
        massD0Pi = RecoDecay::m(std::move(arrayMomenta), array{massD0, massPi});
        hMassBplusToD0Pi->Fill(massD0Pi);
        // B+ candidate invariant mass selction window
        if (massD0Pi < 4.5 || massD0Pi > 6.) {
          continue;
        }

        // fill candidate table rows
        rowCandidateBase(collision.globalIndex(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         BSecVertex[0], BSecVertex[1], BSecVertex[2],
                         errorDecayLength, errorDecayLengthXY,
                         chi2PCA,
                         pVecD0[0], pVecD0[1], pVecD0[2],
                         pVecBach[0], pVecBach[1], pVecBach[2],
                         impactParameter0.getY(), impactParameter1.getY(),
                         std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
                         candidate.globalIndex(), track.globalIndex(), // index D0 and bachelor
                         hfFlag);
      } // track loop
    }   // D0 cand loop
  }     // process
};      // struct

/// Extends the base table with expression columns.
struct HfCandidateCreatorBplusExpressions {
  Spawns<aod::HfCandBplusExt> rowCandidateBPlus;

  void init(InitContext const&) {}
};

/// Performs MC matching.
struct HfCandidateCreatorBplusMc {
  Produces<aod::HfCandBplusMcRec> rowMcMatchRec;
  Produces<aod::HfCandBplusMcGen> rowMcMatchGen;

  void processMc(aod::HfCandBplus const& candidates,
                 aod::HfCand2Prong const&,
                 aod::BigTracksMC const& tracks,
                 aod::McParticles const& particlesMC)
  {
    int indexRec = -1, indexRecD0 = -1;
    int8_t signB = 0, signD0 = 0;
    int8_t flag = 0;
    int kD0pdg = pdg::Code::kD0;

    // Match reconstructed candidates.
    for (auto& candidate : candidates) {
      // Printf("New rec. candidate");

      flag = 0;
      auto candDaughterD0 = candidate.prong0_as<aod::HfCand2Prong>();
      auto arrayDaughtersD0 = array{candDaughterD0.prong0_as<aod::BigTracksMC>(), candDaughterD0.prong1_as<aod::BigTracksMC>()};
      auto arrayDaughters = array{candidate.prong1_as<aod::BigTracksMC>(), candDaughterD0.prong0_as<aod::BigTracksMC>(), candDaughterD0.prong1_as<aod::BigTracksMC>()};

      // B± → D0bar(D0) π± → (K± π∓) π±
      // Printf("Checking B± → D0(bar) π±");
      indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kBPlus, array{+kPiPlus, +kKPlus, -kPiPlus}, true, &signB, 2);
      indexRecD0 = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersD0, pdg::Code::kD0, array{-kKPlus, +kPiPlus}, true, &signD0, 1);

      if (indexRecD0 > -1 && indexRec > -1) {
        flag = signB * (1 << hf_cand_bplus::DecayType::BplusToD0Pi);
      }
      rowMcMatchRec(flag);
    }

    // Match generated particles.
    for (auto& particle : particlesMC) {
      // Printf("New gen. candidate");
      flag = 0;
      signB = 0;
      signD0 = 0;
      int indexGenD0 = -1;

      // B± → D0bar(D0) π± → (K± π∓) π±
      // Printf("Checking B± → D0(bar) π±");
      std::vector<int> arrayDaughterB;
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kBPlus, array{-kD0pdg, +kPiPlus}, true, &signB, 1, &arrayDaughterB)) {
        // D0(bar) → π± K∓
        // Printf("Checking D0(bar) → π± K∓");
        for (auto iD : arrayDaughterB) {
          auto candDaughterMC = particlesMC.rawIteratorAt(iD);
          if (std::abs(candDaughterMC.pdgCode()) == kD0pdg) {
            indexGenD0 = RecoDecay::isMatchedMCGen(particlesMC, candDaughterMC, pdg::Code::kD0, array{-kKPlus, +kPiPlus}, true, &signD0, 1);
          }
        }
        if (indexGenD0 > -1) {
          flag = signB * (1 << hf_cand_bplus::DecayType::BplusToD0Pi);
        }
      }
      rowMcMatchGen(flag);
    } // B candidate
  }   // process
  PROCESS_SWITCH(HfCandidateCreatorBplusMc, processMc, "Process MC", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfCandidateCreatorBplus>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorBplusExpressions>(cfgc)};
  workflow.push_back(adaptAnalysisTask<HfCandidateCreatorBplusMc>(cfgc));
  return workflow;
}
