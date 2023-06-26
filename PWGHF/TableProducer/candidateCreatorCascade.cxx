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

/// \file candidateCreatorCascade.cxx
/// \brief Reconstruction of heavy-flavour cascade decay candidates
///
/// \author Chiara Zampolli, <Chiara.Zampolli@cern.ch>, CERN
///         Paul Buehler, <paul.buehler@oeaw.ac.at>, Vienna

#include "Common/Core/trackUtilities.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsDebugLcToK0sP.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_2prong;

// #define MY_DEBUG

#ifdef MY_DEBUG
using MyBigTracks = aod::BigTracksMC;
#define MY_DEBUG_MSG(condition, cmd) \
  if (condition) {                   \
    cmd;                             \
  }
#else
using MyBigTracks = aod::BigTracks;
#define MY_DEBUG_MSG(condition, cmd)
#endif

/// Reconstruction of heavy-flavour cascade decay candidates
struct HfCandidateCreatorCascade {
  Produces<aod::HfCandCascBase> rowCandidateBase;

  // vertexing
  // Configurable<double> bz{"bz", 5., "magnetic field"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", true, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "fill validation histograms"};
  Configurable<bool> silenceV0DataWarning{"silenceV0DataWarning", false, "do not print a warning for not found V0s and silently skip them"};

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

  // for debugging
#ifdef MY_DEBUG
  Configurable<std::vector<int>> indexK0Spos{"indexK0Spos", {729, 2866, 4754, 5457, 6891, 7824, 9243, 9810}, "indices of K0S positive daughters, for debug"};
  Configurable<std::vector<int>> indexK0Sneg{"indexK0Sneg", {730, 2867, 4755, 5458, 6892, 7825, 9244, 9811}, "indices of K0S negative daughters, for debug"};
  Configurable<std::vector<int>> indexProton{"indexProton", {717, 2810, 4393, 5442, 6769, 7793, 9002, 9789}, "indices of protons, for debug"};
#endif

  double massP = RecoDecay::getMassPDG(kProton);
  double massK0s = RecoDecay::getMassPDG(kK0Short);
  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massLc = RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus);
  double mass2K0sP{0.};
  double bz = 0.;

  OutputObj<TH1F> hMass2{TH1F("hMass2", "2-prong candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }

  void process(aod::Collisions const&,
               aod::HfCascades const& rowsTrackIndexCasc,
               MyBigTracks const&,
               aod::V0sLinked const&,
               aod::V0Datas const&,
               aod::BCsWithTimestamps const&
#ifdef MY_DEBUG
               ,
               aod::McParticles& mcParticles
#endif
  )
  {
    // 2-prong vertex fitter
    o2::vertexing::DCAFitterN<2> df;
    // df.setBz(bz);
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

    // loop over pairs of track indeces
    for (const auto& casc : rowsTrackIndexCasc) {

      const auto& bach = casc.prong0_as<MyBigTracks>();
      LOGF(debug, "V0 %d in HF cascade %d.", casc.v0Id(), casc.globalIndex());
      if (!casc.has_v0()) {
        LOGF(error, "V0 not there for HF cascade %d. Skipping candidate.", casc.globalIndex());
        continue;
      }
      if (!casc.v0_as<aod::V0sLinked>().has_v0Data()) {
        if (!silenceV0DataWarning) {
          LOGF(warning, "V0Data not there for V0 %d in HF cascade %d. Skipping candidate.", casc.v0Id(), casc.globalIndex());
        }
        continue;
      }
      LOGF(debug, "V0Data ID: %d", casc.v0_as<aod::V0sLinked>().v0DataId());
      const auto& v0 = casc.v0_as<aod::V0sLinked>().v0Data();
      const auto& trackV0DaughPos = v0.posTrack_as<MyBigTracks>();
      const auto& trackV0DaughNeg = v0.negTrack_as<MyBigTracks>();

      auto collision = casc.collision();

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        // df.setBz(bz); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
      }
      df.setBz(bz);

#ifdef MY_DEBUG
      auto indexBach = bach.mcParticleId();
      auto indexV0DaughPos = trackV0DaughPos.mcParticleId();
      auto indexV0DaughNeg = trackV0DaughNeg.mcParticleId();
      bool isLc = isLcK0SpFunc(indexBach, indexV0DaughPos, indexV0DaughNeg, indexProton, indexK0Spos, indexK0Sneg);
#endif

      MY_DEBUG_MSG(isLc, LOG(info) << "Processing the Lc with proton " << indexBach << " trackV0DaughPos " << indexV0DaughPos << " trackV0DaughNeg " << indexV0DaughNeg);

      auto trackParCovBach = getTrackParCov(bach);
      auto trackParCovV0DaughPos = getTrackParCov(trackV0DaughPos); // check that MyBigTracks does not need TracksDCA!
      auto trackParCovV0DaughNeg = getTrackParCov(trackV0DaughNeg); // check that MyBigTracks does not need TracksDCA!
      trackParCovV0DaughPos.propagateTo(v0.posX(), bz);             // propagate the track to the X closest to the V0 vertex
      trackParCovV0DaughNeg.propagateTo(v0.negX(), bz);             // propagate the track to the X closest to the V0 vertex
      const std::array<float, 3> vertexV0 = {v0.x(), v0.y(), v0.z()};
      const std::array<float, 3> momentumV0 = {v0.px(), v0.py(), v0.pz()};
      // we build the neutral track to then build the cascade
      auto trackV0 = o2::dataformats::V0(vertexV0, momentumV0, {0, 0, 0, 0, 0, 0}, trackParCovV0DaughPos, trackParCovV0DaughNeg, {0, 0}, {0, 0}); // build the V0 track (indices for v0 daughters set to 0 for now)

      // reconstruct the cascade secondary vertex
      if (df.process(trackV0, trackParCovBach) == 0) {
        MY_DEBUG_MSG(isLc, LOG(info) << "Vertexing failed for Lc candidate");
        //  if (isLc) {
        // LOG(info) << "Vertexing failed for Lc with proton " << indexBach << " trackV0DaughPos " << indexV0DaughPos << " trackV0DaughNeg " << indexV0DaughNeg;
        //}
        continue;
      } else {
        // LOG(info) << "Vertexing succeeded for Lc candidate";
      }

      const auto& secondaryVertex = df.getPCACandidate();
      auto chi2PCA = df.getChi2AtPCACandidate();
      auto covMatrixPCA = df.calcPCACovMatrixFlat();
      hCovSVXX->Fill(covMatrixPCA[0]); // FIXME: Calculation of errorDecayLength(XY) gives wrong values without this line.
      // do I have to call "df.propagateTracksToVertex();"?
      auto trackParVarV0 = df.getTrack(0);
      auto trackParVarBach = df.getTrack(1);

      // get track momenta
      array<float, 3> pVecV0;
      array<float, 3> pVecBach;
      trackParVarV0.getPxPyPzGlo(pVecV0);
      trackParVarBach.getPxPyPzGlo(pVecBach);

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      hCovPVXX->Fill(covMatrixPV[0]);
      o2::dataformats::DCA impactParameterV0;
      o2::dataformats::DCA impactParameterBach;
      trackParVarV0.propagateToDCA(primaryVertex, bz, &impactParameterV0); // we do this wrt the primary vtx
      trackParVarBach.propagateToDCA(primaryVertex, bz, &impactParameterBach);

      // get uncertainty of the decay length
      double phi, theta;
      getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

      // fill candidate table rows
      MY_DEBUG_MSG(isLc, LOG(info) << "IT IS A Lc! Filling for Lc with proton " << indexBach << " trackV0DaughPos " << indexV0DaughPos << " trackV0DaughNeg " << indexV0DaughNeg);
      rowCandidateBase(collision.globalIndex(),
                       collision.posX(), collision.posY(), collision.posZ(),
                       secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                       errorDecayLength, errorDecayLengthXY,
                       chi2PCA,
                       pVecBach[0], pVecBach[1], pVecBach[2],
                       pVecV0[0], pVecV0[1], pVecV0[2],
                       impactParameterBach.getY(), impactParameterV0.getY(),
                       std::sqrt(impactParameterBach.getSigmaY2()), std::sqrt(impactParameterV0.getSigmaY2()),
                       casc.prong0Id(), casc.v0Id(),
                       v0.x(), v0.y(), v0.z(),
                       // v0.posTrack(), v0.negTrack(), // why this was not fine?
                       trackV0DaughPos.globalIndex(), trackV0DaughNeg.globalIndex(),
                       v0.pxpos(), v0.pypos(), v0.pzpos(),
                       v0.pxneg(), v0.pyneg(), v0.pzneg(),
                       v0.dcaV0daughters(),
                       v0.dcapostopv(),
                       v0.dcanegtopv());

      // fill histograms
      if (fillHistograms) {
        // calculate invariant masses
        mass2K0sP = RecoDecay::m(array{pVecBach, pVecV0}, array{massP, massK0s});
        hMass2->Fill(mass2K0sP);
      }
    }
  }
};

/// Performs MC matching.
struct HfCandidateCreatorCascadeMc {
  Spawns<aod::HfCandCascExt> rowCandidateCasc;
  Produces<aod::HfCandCascadeMcRec> rowMcMatchRec;
  Produces<aod::HfCandCascadeMcGen> rowMcMatchGen;

#ifdef MY_DEBUG
  Configurable<std::vector<int>> indexK0Spos{"indexK0Spos", {729, 2866, 4754, 5457, 6891, 7824, 9243, 9810}, "indices of K0S positive daughters, for debug"};
  Configurable<std::vector<int>> indexK0Sneg{"indexK0Sneg", {730, 2867, 4755, 5458, 6892, 7825, 9244, 9811}, "indices of K0S negative daughters, for debug"};
  Configurable<std::vector<int>> indexProton{"indexProton", {717, 2810, 4393, 5442, 6769, 7793, 9002, 9789}, "indices of protons, for debug"};
#endif

  void processMc(aod::BigTracksMC const& tracks,
                 aod::McParticles const& particlesMC)
  {
    int8_t sign = 0;
    int8_t origin = 0;
    int indexRec = -1;
    std::vector<int> arrDaughLcIndex;
    std::array<int, 3> arrDaughLcPDG;
    std::array<int, 3> arrDaughLcPDGRef = {2212, 211, -211};

    // Match reconstructed candidates.
    rowCandidateCasc->bindExternalIndices(&tracks);
    for (const auto& candidate : *rowCandidateCasc) {

      origin = 0;

      const auto& bach = candidate.prong0_as<aod::BigTracksMC>();
      const auto& trackV0DaughPos = candidate.posTrack_as<aod::BigTracksMC>();
      const auto& trackV0DaughNeg = candidate.negTrack_as<aod::BigTracksMC>();

      auto arrayDaughtersV0 = array{trackV0DaughPos, trackV0DaughNeg};
      auto arrayDaughtersLc = array{bach, trackV0DaughPos, trackV0DaughNeg};

      // First we check the K0s
      LOG(debug) << "\n";
      LOG(debug) << "Checking MC for candidate!";
      LOG(debug) << "Looking for K0s";
#ifdef MY_DEBUG
      auto indexV0DaughPos = trackV0DaughPos.mcParticleId();
      auto indexV0DaughNeg = trackV0DaughNeg.mcParticleId();
      auto indexBach = bach.mcParticleId();
      bool isLc = isLcK0SpFunc(indexBach, indexV0DaughPos, indexV0DaughNeg, indexProton, indexK0Spos, indexK0Sneg);
      bool isK0SfromLc = isK0SfromLcFunc(indexV0DaughPos, indexV0DaughNeg, indexK0Spos, indexK0Sneg);
#endif
      MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "correct K0S in the Lc daughters: posTrack --> " << indexV0DaughPos << ", negTrack --> " << indexV0DaughNeg);

      // if (isLc) {
      RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersV0, kK0Short, array{+kPiPlus, -kPiPlus}, true, &sign, 1); // does it matter the "acceptAntiParticle" in the K0s case? In principle, there is no anti-K0s

      if (sign != 0) { // we have already positively checked the K0s
        // then we check the Lc
        MY_DEBUG_MSG(sign, LOG(info) << "K0S was correct! now we check the Lc");
        MY_DEBUG_MSG(sign, LOG(info) << "index proton = " << indexBach);
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersLc, pdg::Code::kLambdaCPlus, array{+kProton, +kPiPlus, -kPiPlus}, true, &sign, 3); // 3-levels Lc --> p + K0 --> p + K0s --> p + pi+ pi-
        MY_DEBUG_MSG(sign, LOG(info) << "Lc found with sign " << sign; printf("\n"));
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (sign != 0) {
        auto particle = particlesMC.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
      }

      rowMcMatchRec(sign, origin);
    }
    //}

    // Match generated particles.
    for (const auto& particle : particlesMC) {
      origin = 0;
      // checking if I have a Lc --> K0S + p
      RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kLambdaCPlus, array{+kProton, +kK0Short}, true, &sign, 2);
      if (sign != 0) {
        MY_DEBUG_MSG(sign, LOG(info) << "Lc in K0S p");
        arrDaughLcIndex.clear();
        // checking that the final daughters (decay depth = 3) are p, pi+, pi-
        RecoDecay::getDaughters(particle, &arrDaughLcIndex, arrDaughLcPDGRef, 3); // best would be to check the K0S daughters
        if (arrDaughLcIndex.size() == 3) {
          for (std::size_t iProng = 0; iProng < arrDaughLcIndex.size(); ++iProng) {
            auto daughI = particlesMC.rawIteratorAt(arrDaughLcIndex[iProng]);
            arrDaughLcPDG[iProng] = daughI.pdgCode();
          }
          if (!(arrDaughLcPDG[0] == arrDaughLcPDGRef[0] && arrDaughLcPDG[1] == arrDaughLcPDGRef[1] && arrDaughLcPDG[2] == arrDaughLcPDGRef[2])) { // this should be the condition, first bach, then v0
            sign = 0;
          } else {
            LOG(debug) << "Lc --> K0S+p found in MC table";
          }
          MY_DEBUG_MSG(sign == 0, LOG(info) << "Pity, the three final daughters are not p, pi+, pi-, but " << arrDaughLcPDG[0] << ", " << arrDaughLcPDG[1] << ", " << arrDaughLcPDG[2]);
        }
      }
      // Check whether the particle is non-prompt (from a b quark).
      if (sign != 0) {
        origin = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
      }
      rowMcMatchGen(sign, origin);
    }
  }

  PROCESS_SWITCH(HfCandidateCreatorCascadeMc, processMc, "Process MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorCascade>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorCascadeMc>(cfgc)};
}
