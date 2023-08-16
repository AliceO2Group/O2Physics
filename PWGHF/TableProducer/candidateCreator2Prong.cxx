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

/// \file candidateCreator2Prong.cxx
/// \brief Reconstruction of heavy-flavour 2-prong decay candidates
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Pengzhong Lu <pengzhong.lu@cern.ch>, GSI Darmstadt, USTC

#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Common/Core/trackUtilities.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

/// includes KFParticle
#include "Tools/KFparticle/KFUtilities.h"
#include <KFPTrack.h>
#include <KFParticleBase.h>
#include <KFPVertex.h>
#include <KFParticle.h>
#include <KFVertex.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_2prong;

/// Reconstruction of heavy-flavour 2-prong decay candidates
struct HfCandidateCreator2Prong {
  Produces<aod::HfCand2ProngBase> rowCandidateBase;
  Produces<aod::HfCand2ProngKF> rowCandidateKF;

  // vertexing
  Configurable<bool> constrainKfToPv{"constrainKfToPv", true, "constraint KFParticle to PV"};
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
  double massPiK{0.};
  double massKPi{0.};
  double bz = 0.;

  OutputObj<TH1F> hMass2{TH1F("hMass2", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hKFMassFailed{TH1F("hKFMassFailed", ";inv. mass with topol fit failed;entries", 500, 0., 5.)}; // 0 for D0 failed, 1 for D0bar failed
  OutputObj<TH1F> hMass2Topol{TH1F("hMass2Topol", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", 600, -1.0, 5.)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH1F> hCovPVYY{TH1F("hCovPVYY", "2-prong candidates;YY element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVYY{TH1F("hCovSVYY", "2-prong candidates;YY element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH1F> hCovPVXZ{TH1F("hCovPVXZ", "2-prong candidates;XZ element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, -1.e-4, 1.e-4)};
  OutputObj<TH1F> hCovSVXZ{TH1F("hCovSVXZ", "2-prong candidates;XZ element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, -1.e-4, 0.2)};
  OutputObj<TH1F> hCovPVZZ{TH1F("hCovPVZZ", "2-prong candidates;ZZ element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVZZ{TH1F("hCovSVZZ", "2-prong candidates;ZZ element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH2F> hDcaXYProngs{TH2F("hDcaXYProngs", "DCAxy of 2-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{xy}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};
  OutputObj<TH2F> hDcaZProngs{TH2F("hDcaZProngs", "DCAz of 2-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{z}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};
  OutputObj<TH1F> hUseKForDCAFitter{TH1F("hUseKForDCAFitter", ";Use KF or DCAFitterN;entries", 2, -0.5, 1.5)}; // 0 for KF, 1 for DCAFitterN

  void init(InitContext const&)
  {
    std::array<bool, 2> doprocessDF{doprocessPvRefit, doprocessNoPvRefit};
    std::array<bool, 2> doprocessKF{doprocessPvRefitWithKF, doprocessNoPvRefitWithKF};
    if ((std::accumulate(doprocessDF.begin(), doprocessDF.end(), 0) + std::accumulate(doprocessKF.begin(), doprocessKF.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }

  template <bool doPvRefit = false, typename Cand, typename TTracks>
  void runCreator2Prong(aod::Collisions const& collisions,
                        Cand const& rowsTrackIndexProng2,
                        TTracks const& tracks,
                        aod::BCsWithTimestamps const& bcWithTimeStamps)
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
    hUseKForDCAFitter->Fill(1.0);

    // loop over pairs of track indices
    for (const auto& rowTrackIndexProng2 : rowsTrackIndexProng2) {
      auto track0 = rowTrackIndexProng2.template prong0_as<TTracks>();
      auto track1 = rowTrackIndexProng2.template prong1_as<TTracks>();
      auto trackParVarPos1 = getTrackParCov(track0);
      auto trackParVarNeg1 = getTrackParCov(track1);
      auto collision = rowTrackIndexProng2.collision();

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
        // df.setBz(bz); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
        // df.print();
      }

      df.setBz(bz);

      // reconstruct the 2-prong secondary vertex
      if (df.process(trackParVarPos1, trackParVarNeg1) == 0) {
        continue;
      }
      const auto& secondaryVertex = df.getPCACandidate();
      auto chi2PCA = df.getChi2AtPCACandidate();
      auto covMatrixPCA = df.calcPCACovMatrixFlat();
      hCovSVXX->Fill(covMatrixPCA[0]); // FIXME: Calculation of errorDecayLength(XY) gives wrong values without this line.
      hCovSVYY->Fill(covMatrixPCA[2]);
      hCovSVXZ->Fill(covMatrixPCA[3]);
      hCovSVZZ->Fill(covMatrixPCA[5]);
      auto trackParVar0 = df.getTrack(0);
      auto trackParVar1 = df.getTrack(1);

      // get track momenta
      array<float, 3> pvec0;
      array<float, 3> pvec1;
      trackParVar0.getPxPyPzGlo(pvec0);
      trackParVar1.getPxPyPzGlo(pvec1);

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      if constexpr (doPvRefit) {
        /// use PV refit
        /// Using it in the rowCandidateBase all dynamic columns shall take it into account
        // coordinates
        primaryVertex.setX(rowTrackIndexProng2.pvRefitX());
        primaryVertex.setY(rowTrackIndexProng2.pvRefitY());
        primaryVertex.setZ(rowTrackIndexProng2.pvRefitZ());
        // covariance matrix
        primaryVertex.setSigmaX2(rowTrackIndexProng2.pvRefitSigmaX2());
        primaryVertex.setSigmaXY(rowTrackIndexProng2.pvRefitSigmaXY());
        primaryVertex.setSigmaY2(rowTrackIndexProng2.pvRefitSigmaY2());
        primaryVertex.setSigmaXZ(rowTrackIndexProng2.pvRefitSigmaXZ());
        primaryVertex.setSigmaYZ(rowTrackIndexProng2.pvRefitSigmaYZ());
        primaryVertex.setSigmaZ2(rowTrackIndexProng2.pvRefitSigmaZ2());
        covMatrixPV = primaryVertex.getCov();
      }
      hCovPVXX->Fill(covMatrixPV[0]);
      hCovPVYY->Fill(covMatrixPV[2]);
      hCovPVXZ->Fill(covMatrixPV[3]);
      hCovPVZZ->Fill(covMatrixPV[5]);
      o2::dataformats::DCA impactParameter0;
      o2::dataformats::DCA impactParameter1;
      trackParVar0.propagateToDCA(primaryVertex, bz, &impactParameter0);
      trackParVar1.propagateToDCA(primaryVertex, bz, &impactParameter1);
      hDcaXYProngs->Fill(track0.pt(), impactParameter0.getY() * toMicrometers);
      hDcaXYProngs->Fill(track1.pt(), impactParameter1.getY() * toMicrometers);
      hDcaZProngs->Fill(track0.pt(), impactParameter0.getZ() * toMicrometers);
      hDcaZProngs->Fill(track1.pt(), impactParameter1.getZ() * toMicrometers);

      // get uncertainty of the decay length
      double phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
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
                       impactParameter0.getY(), impactParameter1.getY(),
                       std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
                       rowTrackIndexProng2.prong0Id(), rowTrackIndexProng2.prong1Id(),
                       rowTrackIndexProng2.hfflag());

      // fill histograms
      if (fillHistograms) {
        // calculate invariant masses
        auto arrayMomenta = array{pvec0, pvec1};
        massPiK = RecoDecay::m(arrayMomenta, array{massPi, massK});
        massKPi = RecoDecay::m(arrayMomenta, array{massK, massPi});
        hMass2->Fill(massPiK);
        hMass2->Fill(massKPi);
      }
    }
  }

  template <bool doPvRefit = false, typename Cand, typename TTracks>
  void runCreator2ProngWithKF(aod::Collisions const& collisions,
                              Cand const& rowsTrackIndexProng2,
                              TTracks const& tracks,
                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    hUseKForDCAFitter->Fill(0.0);

    for (const auto& rowTrackIndexProng2 : rowsTrackIndexProng2) {
      auto track0 = rowTrackIndexProng2.template prong0_as<TTracks>();
      auto track1 = rowTrackIndexProng2.template prong1_as<TTracks>();
      auto trackParVarPos1 = getTrackParCov(track0);
      auto trackParVarNeg1 = getTrackParCov(track1);
      auto collision = rowTrackIndexProng2.collision();

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
        // df.setBz(bz); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
        // df.print();
      }
      float covMatrixPV[6];

      KFParticle::SetField(bz);
      KFPVertex kfpVertex = createKFPVertexFromCollision(collision);

      if constexpr (doPvRefit) {
        /// use PV refit
        /// Using it in the rowCandidateBase all dynamic columns shall take it into account
        // coordinates
        kfpVertex.SetXYZ(rowTrackIndexProng2.pvRefitX(), rowTrackIndexProng2.pvRefitY(), rowTrackIndexProng2.pvRefitZ());
        // covariance matrix
        kfpVertex.SetCovarianceMatrix(rowTrackIndexProng2.pvRefitSigmaX2(), rowTrackIndexProng2.pvRefitSigmaXY(), rowTrackIndexProng2.pvRefitSigmaY2(), rowTrackIndexProng2.pvRefitSigmaXZ(), rowTrackIndexProng2.pvRefitSigmaYZ(), rowTrackIndexProng2.pvRefitSigmaZ2());
      }
      kfpVertex.GetCovarianceMatrix(covMatrixPV);
      KFParticle KFPV(kfpVertex);
      hCovPVXX->Fill(covMatrixPV[0]);
      hCovPVYY->Fill(covMatrixPV[2]);
      hCovPVXZ->Fill(covMatrixPV[3]);
      hCovPVZZ->Fill(covMatrixPV[5]);

      KFPTrack kfpTrack0;
      KFPTrack kfpTrack1;

      kfpTrack0 = createKFPTrackFromTrack(track0);
      kfpTrack1 = createKFPTrackFromTrack(track1);

      KFParticle kfPosPion(kfpTrack0, kPiPlus);
      KFParticle kfNegPion(kfpTrack1, kPiPlus);
      KFParticle kfPosKaon(kfpTrack0, kKPlus);
      KFParticle kfNegKaon(kfpTrack1, kKPlus);

      float impactParameter0XY = 0., errImpactParameter0XY = 0., impactParameter1XY = 0., errImpactParameter1XY = 0.;
      if (!kfPosPion.GetDistanceFromVertexXY(KFPV, impactParameter0XY, errImpactParameter0XY)) {
        hDcaXYProngs->Fill(track0.pt(), impactParameter0XY * toMicrometers);
        hDcaZProngs->Fill(track0.pt(), std::sqrt(kfPosPion.GetDistanceFromVertex(KFPV) * kfPosPion.GetDistanceFromVertex(KFPV) - impactParameter0XY * impactParameter0XY) * toMicrometers);
      } else {
        hDcaXYProngs->Fill(track0.pt(), -999.);
        hDcaZProngs->Fill(track0.pt(), -999.);
      }
      if (!kfNegPion.GetDistanceFromVertexXY(KFPV, impactParameter1XY, errImpactParameter1XY)) {
        hDcaXYProngs->Fill(track1.pt(), impactParameter1XY * toMicrometers);
        hDcaZProngs->Fill(track1.pt(), std::sqrt(kfNegPion.GetDistanceFromVertex(KFPV) * kfNegPion.GetDistanceFromVertex(KFPV) - impactParameter1XY * impactParameter1XY) * toMicrometers);
      } else {
        hDcaXYProngs->Fill(track1.pt(), -999.);
        hDcaZProngs->Fill(track1.pt(), -999.);
      }

      KFParticle kfCandD0;
      const KFParticle* kfDaughtersD0[2] = {&kfPosPion, &kfNegKaon};
      kfCandD0.SetConstructMethod(2);
      kfCandD0.Construct(kfDaughtersD0, 2);
      KFParticle kfCandD0bar;
      const KFParticle* kfDaughtersD0bar[2] = {&kfNegPion, &kfPosKaon};
      kfCandD0bar.SetConstructMethod(2);
      kfCandD0bar.Construct(kfDaughtersD0bar, 2);

      auto massD0 = kfCandD0.GetMass();
      auto massD0bar = kfCandD0bar.GetMass();

      hCovSVXX->Fill(kfCandD0.Covariance(0, 0));
      hCovSVYY->Fill(kfCandD0.Covariance(1, 1));
      hCovSVXZ->Fill(kfCandD0.Covariance(2, 0));
      hCovSVZZ->Fill(kfCandD0.Covariance(2, 2));
      auto covMatrixSV = kfCandD0.CovarianceMatrix();

      double phi, theta;
      getPointDirection(std::array{KFPV.GetX(), KFPV.GetY(), KFPV.GetZ()}, array{kfCandD0.GetX(), kfCandD0.GetY(), kfCandD0.GetZ()}, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSV, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixSV, phi, 0.));

      float topolChi2PerNdfD0 = -999.;
      KFParticle kfCandD0Topol2PV;
      if (constrainKfToPv) {
        kfCandD0Topol2PV = kfCandD0;
        kfCandD0Topol2PV.SetProductionVertex(KFPV);
        topolChi2PerNdfD0 = kfCandD0Topol2PV.GetChi2() / kfCandD0Topol2PV.GetNDF();
      }

      // fill candidate table rows
      rowCandidateBase(collision.globalIndex(),
                       KFPV.GetX(), KFPV.GetY(), KFPV.GetZ(),
                       kfCandD0.GetX(), kfCandD0.GetY(), kfCandD0.GetZ(),
                       errorDecayLength, errorDecayLengthXY,   // TODO: much different from the DCAFitterN one
                       kfCandD0.GetChi2() / kfCandD0.GetNDF(), // TODO: to make sure it should be chi2 only or chi2/ndf, much different from the DCAFitterN one
                       kfPosPion.GetPx(), kfPosPion.GetPy(), kfPosPion.GetPz(),
                       kfNegKaon.GetPx(), kfNegKaon.GetPy(), kfNegKaon.GetPz(),
                       impactParameter0XY, impactParameter1XY,
                       errImpactParameter0XY, errImpactParameter1XY,
                       rowTrackIndexProng2.prong0Id(), rowTrackIndexProng2.prong1Id(),
                       rowTrackIndexProng2.hfflag());
      rowCandidateKF(topolChi2PerNdfD0,
                     massD0, massD0bar);

      // fill histograms
      if (fillHistograms) {
        hMass2->Fill(massD0);
        hMass2->Fill(massD0bar);
        if (constrainKfToPv) {
          float massD0Topol = 0., errMassD0Topol = 0.;
          if (kfCandD0Topol2PV.GetMass(massD0Topol, errMassD0Topol))
            hKFMassFailed->Fill(massD0); // Set up D0 topological constraints to PV failed
          hMass2Topol->Fill(massD0Topol);
        }
      }
    }
  }

  void processPvRefit(aod::Collisions const& collisions,
                      soa::Join<aod::Hf2Prongs, aod::HfPvRefit2Prong> const& rowsTrackIndexProng2,
                      aod::TracksWCov const& tracks,
                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2Prong<true>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefit, "Run candidate creator with PV refit", false);

  void processNoPvRefit(aod::Collisions const& collisions,
                        aod::Hf2Prongs const& rowsTrackIndexProng2,
                        aod::TracksWCov const& tracks,
                        aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2Prong(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefit, "Run candidate creator without PV refit", true);

  void processPvRefitWithKF(aod::Collisions const& collisions,
                            soa::Join<aod::Hf2Prongs, aod::HfPvRefit2Prong> const& rowsTrackIndexProng2,
                            soa::Join<aod::TracksWCov, aod::TracksExtra> const& tracks,
                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKF<true>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithKF, "Run candidate creator with PV refit", false);

  void processNoPvRefitWithKF(aod::Collisions const& collisions,
                              aod::Hf2Prongs const& rowsTrackIndexProng2,
                              soa::Join<aod::TracksWCov, aod::TracksExtra> const& tracks,
                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKF(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefitWithKF, "Run candidate creator without PV refit", false);
};

/// Extends the base table with expression columns.
struct HfCandidateCreator2ProngExpressions {
  Spawns<aod::HfCand2ProngExt> rowCandidateProng2;
  Produces<aod::HfCand2ProngMcRec> rowMcMatchRec;
  Produces<aod::HfCand2ProngMcGen> rowMcMatchGen;

  void init(InitContext const&) {}

  /// Performs MC matching.
  void processMc(aod::TracksWMc const& tracks,
                 aod::McParticles const& particlesMC)
  {
    rowCandidateProng2->bindExternalIndices(&tracks);

    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;

    // Match reconstructed candidates.
    // Spawned table can be used directly
    for (auto& candidate : *rowCandidateProng2) {
      // Printf("New rec. candidate");
      flag = 0;
      origin = 0;
      auto arrayDaughters = array{candidate.prong0_as<aod::TracksWMc>(), candidate.prong1_as<aod::TracksWMc>()};

      // D0(bar) → π± K∓
      // Printf("Checking D0(bar) → π± K∓");
      indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kD0, array{+kPiPlus, -kKPlus}, true, &sign);
      if (indexRec > -1) {
        flag = sign * (1 << DecayType::D0ToPiK);
      }

      // J/ψ → e+ e−
      if (flag == 0) {
        // Printf("Checking J/ψ → e+ e−");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kJPsi, array{+kElectron, -kElectron}, true);
        if (indexRec > -1) {
          flag = 1 << DecayType::JpsiToEE;
        }
      }

      // J/ψ → μ+ μ−
      if (flag == 0) {
        // Printf("Checking J/ψ → μ+ μ−");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kJPsi, array{+kMuonPlus, -kMuonPlus}, true);
        if (indexRec > -1) {
          flag = 1 << DecayType::JpsiToMuMu;
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flag != 0) {
        auto particle = particlesMC.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
      }

      rowMcMatchRec(flag, origin);
    }

    // Match generated particles.
    for (auto& particle : particlesMC) {
      // Printf("New gen. candidate");
      flag = 0;
      origin = 0;

      // D0(bar) → π± K∓
      // Printf("Checking D0(bar) → π± K∓");
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kD0, array{+kPiPlus, -kKPlus}, true, &sign)) {
        flag = sign * (1 << DecayType::D0ToPiK);
      }

      // J/ψ → e+ e−
      if (flag == 0) {
        // Printf("Checking J/ψ → e+ e−");
        if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kJPsi, array{+kElectron, -kElectron}, true)) {
          flag = 1 << DecayType::JpsiToEE;
        }
      }

      // J/ψ → μ+ μ−
      if (flag == 0) {
        // Printf("Checking J/ψ → μ+ μ−");
        if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kJPsi, array{+kMuonPlus, -kMuonPlus}, true)) {
          flag = 1 << DecayType::JpsiToMuMu;
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flag != 0) {
        origin = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
      }

      rowMcMatchGen(flag, origin);
    }
  }

  PROCESS_SWITCH(HfCandidateCreator2ProngExpressions, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreator2Prong>(cfgc, TaskName{"hf-candidate-creator-2prong"}),
    adaptAnalysisTask<HfCandidateCreator2ProngExpressions>(cfgc, TaskName{"hf-candidate-creator-2prong-expressions"})};
}
