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

/// \file candidateCreator3Prong.cxx
/// \brief Reconstruction of heavy-flavour 3-prong decay candidates
/// \note Extended from candidateCreator2Prong.cxx
///
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "Common/Core/trackUtilities.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_3prong;

/// Reconstruction of heavy-flavour 3-prong decay candidates
struct HfCandidateCreator3Prong {
  Produces<aod::HfCand3ProngBase> rowCandidateBase;

  // vertexing
  Configurable<bool> doPvRefit{"doPvRefit", false, "do PV refit excluding the candidate daughters, if contributors"};
  // Configurable<double> bz{"bz", 5., "magnetic field"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "do validation plots"};
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

  float toMicrometers = 10000.; // from cm to µm

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massK = RecoDecay::getMassPDG(kKPlus);
  double massPiKPi{0.};
  double bz = 0.;

  OutputObj<TH1F> hMass3{TH1F("hMass3", "3-prong candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", 500, 1.6, 2.1)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "3-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "3-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH1F> hCovPVYY{TH1F("hCovPVYY", "3-prong candidates;YY element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVYY{TH1F("hCovSVYY", "3-prong candidates;YY element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH1F> hCovPVXZ{TH1F("hCovPVXZ", "3-prong candidates;XZ element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, -1.e-4, 1.e-4)};
  OutputObj<TH1F> hCovSVXZ{TH1F("hCovSVXZ", "3-prong candidates;XZ element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, -1.e-4, 0.2)};
  OutputObj<TH1F> hCovPVZZ{TH1F("hCovPVZZ", "3-prong candidates;ZZ element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVZZ{TH1F("hCovSVZZ", "3-prong candidates;ZZ element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH2F> hDcaXYProngs{TH2F("hDcaXYProngs", "DCAxy of 3-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{xy}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};
  OutputObj<TH2F> hDcaZProngs{TH2F("hDcaZProngs", "DCAz of 3-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{z}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }

  void process(aod::Collisions const& collisions,
               soa::Join<aod::Hf3Prongs, aod::HfPvRefit3Prong> const& rowsTrackIndexProng3,
               aod::BigTracks const& tracks,
               aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // 3-prong vertex fitter
    o2::vertexing::DCAFitterN<3> df;
    // df.setBz(bz);
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

    // loop over triplets of track indices
    for (const auto& rowTrackIndexProng3 : rowsTrackIndexProng3) {
      auto track0 = rowTrackIndexProng3.prong0_as<aod::BigTracks>();
      auto track1 = rowTrackIndexProng3.prong1_as<aod::BigTracks>();
      auto track2 = rowTrackIndexProng3.prong2_as<aod::BigTracks>();
      auto trackParVar0 = getTrackParCov(track0);
      auto trackParVar1 = getTrackParCov(track1);
      auto trackParVar2 = getTrackParCov(track2);
      auto collision = rowTrackIndexProng3.collision();

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
        // df.setBz(bz); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
        // df.print();
      }
      df.setBz(bz);

      // reconstruct the 3-prong secondary vertex
      if (df.process(trackParVar0, trackParVar1, trackParVar2) == 0) {
        continue;
      }
      const auto& secondaryVertex = df.getPCACandidate();
      auto chi2PCA = df.getChi2AtPCACandidate();
      auto covMatrixPCA = df.calcPCACovMatrixFlat();
      hCovSVXX->Fill(covMatrixPCA[0]); // FIXME: Calculation of errorDecayLength(XY) gives wrong values without this line.
      hCovSVYY->Fill(covMatrixPCA[2]);
      hCovSVXZ->Fill(covMatrixPCA[3]);
      hCovSVZZ->Fill(covMatrixPCA[5]);
      trackParVar0 = df.getTrack(0);
      trackParVar1 = df.getTrack(1);
      trackParVar2 = df.getTrack(2);

      // get track momenta
      array<float, 3> pvec0;
      array<float, 3> pvec1;
      array<float, 3> pvec2;
      trackParVar0.getPxPyPzGlo(pvec0);
      trackParVar1.getPxPyPzGlo(pvec1);
      trackParVar2.getPxPyPzGlo(pvec2);

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      if (doPvRefit) {
        /// use PV refit
        /// Using it in the rowCandidateBase all dynamic columns shall take it into account
        // coordinates
        primaryVertex.setX(rowTrackIndexProng3.pvRefitX());
        primaryVertex.setY(rowTrackIndexProng3.pvRefitY());
        primaryVertex.setZ(rowTrackIndexProng3.pvRefitZ());
        // covariance matrix
        primaryVertex.setSigmaX2(rowTrackIndexProng3.pvRefitSigmaX2());
        primaryVertex.setSigmaXY(rowTrackIndexProng3.pvRefitSigmaXY());
        primaryVertex.setSigmaY2(rowTrackIndexProng3.pvRefitSigmaY2());
        primaryVertex.setSigmaXZ(rowTrackIndexProng3.pvRefitSigmaXZ());
        primaryVertex.setSigmaYZ(rowTrackIndexProng3.pvRefitSigmaYZ());
        primaryVertex.setSigmaZ2(rowTrackIndexProng3.pvRefitSigmaZ2());
        covMatrixPV = primaryVertex.getCov();
      }
      hCovPVXX->Fill(covMatrixPV[0]);
      hCovPVYY->Fill(covMatrixPV[2]);
      hCovPVXZ->Fill(covMatrixPV[3]);
      hCovPVZZ->Fill(covMatrixPV[5]);
      o2::dataformats::DCA impactParameter0;
      o2::dataformats::DCA impactParameter1;
      o2::dataformats::DCA impactParameter2;
      trackParVar0.propagateToDCA(primaryVertex, bz, &impactParameter0);
      trackParVar1.propagateToDCA(primaryVertex, bz, &impactParameter1);
      trackParVar2.propagateToDCA(primaryVertex, bz, &impactParameter2);
      hDcaXYProngs->Fill(track0.pt(), impactParameter0.getY() * toMicrometers);
      hDcaXYProngs->Fill(track1.pt(), impactParameter1.getY() * toMicrometers);
      hDcaXYProngs->Fill(track2.pt(), impactParameter2.getY() * toMicrometers);
      hDcaZProngs->Fill(track0.pt(), impactParameter0.getZ() * toMicrometers);
      hDcaZProngs->Fill(track1.pt(), impactParameter1.getZ() * toMicrometers);
      hDcaZProngs->Fill(track2.pt(), impactParameter2.getZ() * toMicrometers);

      // get uncertainty of the decay length
      double phi, theta;
      getPointDirection(array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

      // fill candidate table rows
      rowCandidateBase(collision.globalIndex(),
                       primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                       secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                       errorDecayLength, errorDecayLengthXY,
                       chi2PCA,
                       pvec0[0], pvec0[1], pvec0[2],
                       pvec1[0], pvec1[1], pvec1[2],
                       pvec2[0], pvec2[1], pvec2[2],
                       impactParameter0.getY(), impactParameter1.getY(), impactParameter2.getY(),
                       std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()), std::sqrt(impactParameter2.getSigmaY2()),
                       rowTrackIndexProng3.prong0Id(), rowTrackIndexProng3.prong1Id(), rowTrackIndexProng3.prong2Id(),
                       rowTrackIndexProng3.hfflag());

      // fill histograms
      if (fillHistograms) {
        // calculate invariant mass
        auto arrayMomenta = array{pvec0, pvec1, pvec2};
        massPiKPi = RecoDecay::m(std::move(arrayMomenta), array{massPi, massK, massPi});
        hMass3->Fill(massPiKPi);
      }
    }
  }
};

/// Extends the base table with expression columns.
struct HfCandidateCreator3ProngExpressions {
  Spawns<aod::HfCand3ProngExt> rowCandidateProng3;
  Produces<aod::HfCand3ProngMcRec> rowMcMatchRec;
  Produces<aod::HfCand3ProngMcGen> rowMcMatchGen;

  void init(InitContext const&) {}

  /// Performs MC matching.
  void processMc(aod::BigTracksMC const& tracks,
                 aod::McParticles const& particlesMC)
  {
    rowCandidateProng3->bindExternalIndices(&tracks);

    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t swapping = 0;
    int8_t channel = 0;
    std::vector<int> arrDaughIndex;
    std::array<int, 2> arrPDGDaugh;
    std::array<int, 2> arrPDGResonant1 = {kProton, 313};  // Λc± → p± K*
    std::array<int, 2> arrPDGResonant2 = {2224, kKPlus};  // Λc± → Δ(1232)±± K∓
    std::array<int, 2> arrPDGResonant3 = {3124, kPiPlus}; // Λc± → Λ(1520) π±
    std::array<int, 2> arrPDGResonantDsPhiPi = {333, kPiPlus}; // Ds± → Phi π±
    std::array<int, 2> arrPDGResonantDsKstarK = {313, kKPlus}; // Ds± → K*(892)0bar K±

    // Match reconstructed candidates.
    // Spawned table can be used directly
    for (auto& candidate : *rowCandidateProng3) {
      // Printf("New rec. candidate");
      flag = 0;
      origin = 0;
      swapping = 0;
      channel = 0;
      arrDaughIndex.clear();
      auto arrayDaughters = array{candidate.prong0_as<aod::BigTracksMC>(), candidate.prong1_as<aod::BigTracksMC>(), candidate.prong2_as<aod::BigTracksMC>()};

      // D± → π± K∓ π±
      // Printf("Checking D± → π± K∓ π±");
      indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kDPlus, array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2);
      if (indexRec > -1) {
        flag = sign * (1 << DecayType::DplusToPiKPi);
      }

      // Ds± → K± K∓ π±
      if (flag == 0) {
        // Printf("Checking Ds± → K± K∓ π±");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kDS, array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec > -1) {
          flag = sign * (1 << DecayType::DsToKKPi);
          if (arrayDaughters[0].has_mcParticle()) {
            swapping = int8_t(std::abs(arrayDaughters[0].mcParticle().pdgCode()) == kPiPlus);
          }
          RecoDecay::getDaughters(particlesMC.rawIteratorAt(indexRec), &arrDaughIndex, array{0}, 1);
          if (arrDaughIndex.size() == 2) {
            for (auto iProng = 0u; iProng < arrDaughIndex.size(); ++iProng) {
              auto daughI = particlesMC.rawIteratorAt(arrDaughIndex[iProng]);
              arrPDGDaugh[iProng] = std::abs(daughI.pdgCode());
            }
            if ((arrPDGDaugh[0] == arrPDGResonantDsPhiPi[0] && arrPDGDaugh[1] == arrPDGResonantDsPhiPi[1]) || (arrPDGDaugh[0] == arrPDGResonantDsPhiPi[1] && arrPDGDaugh[1] == arrPDGResonantDsPhiPi[0])) {
              channel = DecayChannelDs::PhiPi;
            } else if ((arrPDGDaugh[0] == arrPDGResonantDsKstarK[0] && arrPDGDaugh[1] == arrPDGResonantDsKstarK[1]) || (arrPDGDaugh[0] == arrPDGResonantDsKstarK[1] && arrPDGDaugh[1] == arrPDGResonantDsKstarK[0])) {
              channel = DecayChannelDs::K0starK;
            }
          }
        }
      }

      // Λc± → p± K∓ π±
      if (flag == 0) {
        // Printf("Checking Λc± → p± K∓ π±");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec > -1) {
          flag = sign * (1 << DecayType::LcToPKPi);

          // Printf("Flagging the different Λc± → p± K∓ π± decay channels");
          if (arrayDaughters[0].has_mcParticle()) {
            swapping = int8_t(std::abs(arrayDaughters[0].mcParticle().pdgCode()) == kPiPlus);
          }
          RecoDecay::getDaughters(particlesMC.rawIteratorAt(indexRec), &arrDaughIndex, array{0}, 1);
          if (arrDaughIndex.size() == 2) {
            for (auto iProng = 0u; iProng < arrDaughIndex.size(); ++iProng) {
              auto daughI = particlesMC.rawIteratorAt(arrDaughIndex[iProng]);
              arrPDGDaugh[iProng] = std::abs(daughI.pdgCode());
            }
            if ((arrPDGDaugh[0] == arrPDGResonant1[0] && arrPDGDaugh[1] == arrPDGResonant1[1]) || (arrPDGDaugh[0] == arrPDGResonant1[1] && arrPDGDaugh[1] == arrPDGResonant1[0])) {
              channel = 1;
            } else if ((arrPDGDaugh[0] == arrPDGResonant2[0] && arrPDGDaugh[1] == arrPDGResonant2[1]) || (arrPDGDaugh[0] == arrPDGResonant2[1] && arrPDGDaugh[1] == arrPDGResonant2[0])) {
              channel = 2;
            } else if ((arrPDGDaugh[0] == arrPDGResonant3[0] && arrPDGDaugh[1] == arrPDGResonant3[1]) || (arrPDGDaugh[0] == arrPDGResonant3[1] && arrPDGDaugh[1] == arrPDGResonant3[0])) {
              channel = 3;
            }
          }
        }
      }

      // Ξc± → p± K∓ π±
      if (flag == 0) {
        // Printf("Checking Ξc± → p± K∓ π±");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kXiCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec > -1) {
          flag = sign * (1 << DecayType::XicToPKPi);
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flag != 0) {
        auto particle = particlesMC.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
      }

      rowMcMatchRec(flag, origin, swapping, channel);
    }

    // Match generated particles.
    for (auto& particle : particlesMC) {
      // Printf("New gen. candidate");
      flag = 0;
      origin = 0;
      channel = 0;
      arrDaughIndex.clear();

      // D± → π± K∓ π±
      // Printf("Checking D± → π± K∓ π±");
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kDPlus, array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2)) {
        flag = sign * (1 << DecayType::DplusToPiKPi);
      }

      // Ds± → K± K∓ π±
      if (flag == 0) {
        // Printf("Checking Ds± → K± K∓ π±");
        if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kDS, array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          flag = sign * (1 << DecayType::DsToKKPi);
          RecoDecay::getDaughters(particle, &arrDaughIndex, array{0}, 1);
          if (arrDaughIndex.size() == 2) {
            for (auto jProng = 0u; jProng < arrDaughIndex.size(); ++jProng) {
              auto daughJ = particlesMC.rawIteratorAt(arrDaughIndex[jProng]);
              arrPDGDaugh[jProng] = std::abs(daughJ.pdgCode());
            }
            if ((arrPDGDaugh[0] == arrPDGResonantDsPhiPi[0] && arrPDGDaugh[1] == arrPDGResonantDsPhiPi[1]) || (arrPDGDaugh[0] == arrPDGResonantDsPhiPi[1] && arrPDGDaugh[1] == arrPDGResonantDsPhiPi[0])) {
              channel = DecayChannelDs::PhiPi;
            } else if ((arrPDGDaugh[0] == arrPDGResonantDsKstarK[0] && arrPDGDaugh[1] == arrPDGResonantDsKstarK[1]) || (arrPDGDaugh[0] == arrPDGResonantDsKstarK[1] && arrPDGDaugh[1] == arrPDGResonantDsKstarK[0])) {
              channel = DecayChannelDs::K0starK;
            }
          }
        }
      }

      // Λc± → p± K∓ π±
      if (flag == 0) {
        // Printf("Checking Λc± → p± K∓ π±");
        if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          flag = sign * (1 << DecayType::LcToPKPi);

          // Printf("Flagging the different Λc± → p± K∓ π± decay channels");
          RecoDecay::getDaughters(particle, &arrDaughIndex, array{0}, 1);
          if (arrDaughIndex.size() == 2) {
            for (auto jProng = 0u; jProng < arrDaughIndex.size(); ++jProng) {
              auto daughJ = particlesMC.rawIteratorAt(arrDaughIndex[jProng]);
              arrPDGDaugh[jProng] = std::abs(daughJ.pdgCode());
            }
            if ((arrPDGDaugh[0] == arrPDGResonant1[0] && arrPDGDaugh[1] == arrPDGResonant1[1]) || (arrPDGDaugh[0] == arrPDGResonant1[1] && arrPDGDaugh[1] == arrPDGResonant1[0])) {
              channel = 1;
            } else if ((arrPDGDaugh[0] == arrPDGResonant2[0] && arrPDGDaugh[1] == arrPDGResonant2[1]) || (arrPDGDaugh[0] == arrPDGResonant2[1] && arrPDGDaugh[1] == arrPDGResonant2[0])) {
              channel = 2;
            } else if ((arrPDGDaugh[0] == arrPDGResonant3[0] && arrPDGDaugh[1] == arrPDGResonant3[1]) || (arrPDGDaugh[0] == arrPDGResonant3[1] && arrPDGDaugh[1] == arrPDGResonant3[0])) {
              channel = 3;
            }
          }
        }
      }

      // Ξc± → p± K∓ π±
      if (flag == 0) {
        // Printf("Checking Ξc± → p± K∓ π±");
        if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kXiCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          flag = sign * (1 << DecayType::XicToPKPi);
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flag != 0) {
        origin = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
      }

      rowMcMatchGen(flag, origin, channel);
    }
  }

  PROCESS_SWITCH(HfCandidateCreator3ProngExpressions, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreator3Prong>(cfgc, TaskName{"hf-candidate-creator-3prong"}),
    adaptAnalysisTask<HfCandidateCreator3ProngExpressions>(cfgc, TaskName{"hf-candidate-creator-3prong-expressions"})};
}
