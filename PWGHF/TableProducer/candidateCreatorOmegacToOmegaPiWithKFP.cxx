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
//
/// \file candidateCreatorOmegacToOmegaPiWithKFP.cxx
/// \brief Reconstruction of Omegac0  -> omega pi candidates using KFParticle
/// \author Yunfan Liu <yunfan.liu@cern.ch>

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/V0.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

/// includes KFParticle
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"


#ifndef HomogeneousField
#define HomogeneousField
#endif

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::aod::cascdata;
using namespace o2::aod::v0data;
using namespace o2::aod::hf_track_index;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Reconstruction of omegac0 candidates
struct HfCandidateCreatorOmegacToOmegaPiWithKfp {
  Produces<aod::HfCandOmegac> rowCandidate;
  Produces<aod::HfOmegacKf> kfCandidateData;

  Configurable<bool> doDCAfitter{"doDCAfitter", false, "do DCAfitte before KFParicle and using DCAfitter updated track"};
  Configurable<bool> doPvRefit{"doPvRefit", false, "set to true if you do PV refit in trackIndexSkimCreator.cxx"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", true, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> maxDXYIni{"maxDXYIni", 4., "reject (if>0) PCA candidate if tracks DXY exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<double> maxChi2{"maxChi2", 100., "discard vertices with chi2/Nprongs > this (or sum{DCAi^2}/Nprongs for abs. distance minimization)"};


  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  // V0 cuts
  Configurable<float> lambdaMassWindow{"lambdaMassWindow", 0.0075, "Distance from Lambda mass"};

  // cascade cuts
  Configurable<double> massToleranceCascade{"massToleranceCascade", 0.01, "Invariant mass tolerance for cascade"};

  // for KF particle operation
  Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "KF Construct Method"};
  Configurable<bool> kfUseV0MassConstraint{"kfUseV0MassConstraint", false, "KF: use Lambda mass constraint"};
  Configurable<bool> kfUseCascadeMassConstraint{"kfUseCascadeMassConstraint", false, "KF: use Cascade mass constraint - WARNING: not adequate for inv mass analysis of Xi"};
  Configurable<bool> kfDoDCAFitterPreMinimV0{"kfDoDCAFitterPreMinimV0", false, "KF: do DCAFitter pre-optimization before KF fit to include material corrections for V0"};
  Configurable<bool> kfDoDCAFitterPreMinimCasc{"kfDoDCAFitterPreMinimCasc", false, "KF: do DCAFitter pre-optimization before KF fit to include material corrections for Omega"};
  Configurable<bool> kfDoDCAFitterPreMinimCharmBaryon{"kfDoDCAFitterPreMinimCharmBaryon", true, "KF: do DCAFitter pre-optimization before KF fit to include material corrections for OmegaC0"};
  // hist
  HistogramRegistry registry{"registry",{}};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::vertexing::DCAFitterN<2> df;

  float bz;
  int mRunNumber;

  using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::HfSelCollision>>;
  using MyTracks = soa::Join<aod::TracksWCovDcaExtra, aod::HfPvRefitTrack, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;
  using FilteredHfTrackAssocSel = soa::Filtered<soa::Join<aod::TrackAssoc, aod::HfSelTrack>>;
  using MyKfCascTable = soa::Join<KFCascDatas, aod::KFCascCovs>;
  using MyV0Table = soa::Join<aod::V0Datas, aod::V0Covs>;
  using MySkimIdx = HfCascLf2Prongs;

  Filter filterSelectCollisions = (aod::hf_sel_collision::whyRejectColl == static_cast<uint16_t>(0)); // filter to use only HF selected collisions
  Filter filterSelectTrackIds = (aod::hf_sel_track::isSelProng > static_cast<uint32_t>(0));

  
  Preslice<FilteredHfTrackAssocSel> trackIndicesPerCollision = aod::track_association::collisionId; // aod::hf_track_association::collisionId
  Preslice<MyKfCascTable> KFCascadesPerCollision = aod::cascdata::collisionId;
  Preslice<MySkimIdx> candidatesPerCollision = hf_track_index::collisionId;

  // Helper struct to pass  information
  struct {
    float nSigmaTPCPiFromOmegac;
    float nSigmaTOFPiFromOmegac;
    float nSigmaTPCKaFromCasc;
    float nSigmaTOFKaFromCasc;
    float nSigmaTPCPiFromV0;
    float nSigmaTPCPrFromV0;
    float dcaCascDau;
    float chi2GeoV0;
    float ldlV0;
    float chi2TopoV0ToPv;
    float chi2GeoCasc;
    float ldlCasc;
    float chi2TopoCascToPv;
    float decayLenXYLambda;
    float decayLenXYCasc;
    float cosPaV0ToCasc;     // PA
    float cosPaV0ToPv;       // PA
    float cosPaCascToOmegac; // PA
    float cosPaCascToPv;     // PA
    float massV0;
    float massCasc;
    float ptPiFromOmegac;
    float ptOmegac;
    float rapOmegac;
    float massOmegac;
    float cosThetaStarPiFromOmegac;
    float chi2TopoPiFromOmegacToPv;
    float kfDcaXYPiFromOmegac;
    float chi2TopoV0ToCasc;
    float chi2TopoCascToOmegac;
    float decayLenXYOmegac;
    float chi2GeoOmegac;
    float dcaV0Dau;
    float kfDcaCascDau;
    float kfDcaOmegacDau;
    float kfDcaXYCascToPv;
    float chi2TopoOmegacToPv;
    float cosPaOmegacToPv; // PA
    float ldlOmegac;
    float ctOmegac;
    float chi2MassV0;
    float chi2MassCasc;
    float etaOmegac;
  } kfOmegac0Candidate;

  OutputObj<TH1F> hCollision{TH1F("hCollision", "Counter Collison;entries", 500, -2, 2)};
  OutputObj<TH1F> hPionNumber{TH1F("hPionNumber", "Pion Number;entries", 500, 0, 5)};
  OutputObj<TH1F> hInvMassCharmBaryon{TH1F("hInvMassOmegac", "Omegac candidate invariant mass;inv mass;entries", 500, 2.2, 3.1)};
  OutputObj<TH1F> hKfInvMassCharmBaryon{TH1F("hKfInvMassOmegac", "KFParticle Omegac candidate invariant mass;inv mass;entries", 500, 2.2, 3.1)};
  OutputObj<TH1F> hFitterStatus{TH1F("hFitterStatus", "Charm DCAFitter status;status;entries", 3, -0.5, 2.5)};                     // 0 --> vertex(es) found, 1 --> exception found, 2 --> no vertex found (but no exception)
  OutputObj<TH1F> hCandidateCounter{TH1F("hCandidateCounter", "Candidate counter wrt derived data;status;entries", 4, -0.5, 3.5)}; // 0 --> candidates in derived data table, 1 --> candidates passing testbit selection, 2 --> candidates passing fitter step 3 --> candidates filled in new table

  void init(InitContext const&)
  {
    registry.add("hKFParticleV0Chi2", "hKFParticleV0Chi2", kTH1F, {{1000, -0.10f, 10.0f}});
    registry.add("hKFParticleCascChi2", "hKFParticleCascChi2", kTH1F, {{1000, -0.1f, 10.0f}});
    registry.add("hKFParticleOmegaC0Chi2", "hKFParticleOmegaC0Chi2", kTH1F, {{1000, -0.1f, 10.0f}});

    registry.add("hKFParticleV0TopoChi2", "hKFParticleV0TopoChi2", kTH1F, {{1000, -0.10f, 10.0f}});
    registry.add("hKFParticleCascTopoChi2", "hKFParticleCascTopoChi2", kTH1F, {{1000, -0.1f, 10.0f}});
    registry.add("hKFParticleCascBachTopoChi2", "hKFParticleCascBachTopoChi2", kTH1F, {{1000, -0.1f, 10.0f}});

    registry.add("hKfLambda_ldl", "hKfLambda_ldl", kTH1F, {{1000, 0.0f, 100.0f}});
    registry.add("hKfOmega_ldl", "hKfOmega_ldl", kTH1F, {{1000, 0.0f, 100.0f}});
    registry.add("hKfOmegaC0_ldl", "hKfOmegaC0_ldl", kTH1F, {{1000, 0.0f, 100.0f}});
    registry.add("hDcaXYCascadeToPVKf", "hDcaXYCascadeToPVKf", kTH1F, {{1000, 0.0f, 2.0f}});

    registry.add("hInvMassOmegaMinus", "hInvMassOmegaMinus", kTH1F, {{1000, 1.6f, 2.0f}});

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    mRunNumber = 0;
    bz = 0;
	
	// 2-prong vertex fitter to build the omegac vertex
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMaxDXYIni(maxDXYIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setMaxChi2(maxChi2);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

  };
  
  // TrackParCov to KF converter
  // FIXME: could be an utility somewhere else
  // from Carolina Reetz (thank you!)
  template <typename T>
  KFParticle createKFParticleFromTrackParCov(const o2::track::TrackParametrizationWithError<T>& trackparCov, int charge, float mass)
  {
    std::array<T, 3> xyz, pxpypz;
    float xyzpxpypz[6];
    trackparCov.getPxPyPzGlo(pxpypz);
    trackparCov.getXYZGlo(xyz);
    for (int i{0}; i < 3; ++i) {
      xyzpxpypz[i] = xyz[i];
      xyzpxpypz[i + 3] = pxpypz[i];
    }

    std::array<float, 21> cv;
    try {
      trackparCov.getCovXYZPxPyPzGlo(cv);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to get cov matrix from TrackParCov" << e.what();
    }

    KFParticle kfPart;
    float Mini, SigmaMini, M, SigmaM;
    kfPart.GetMass(Mini, SigmaMini);
    LOG(debug) << "Daughter KFParticle mass before creation: " << Mini << " +- " << SigmaMini;

    try {
      kfPart.Create(xyzpxpypz, cv.data(), charge, mass);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to create KFParticle from daughter TrackParCov" << e.what();
    }

    kfPart.GetMass(M, SigmaM);
    LOG(debug) << "Daughter KFParticle mass after creation: " << M << " +- " << SigmaM;
    return kfPart;
  }

  // KF to TrackParCov converter
  // FIXME: could be an utility somewhere else
  // from Carolina Reetz (thank you!)
  o2::track::TrackParCov getTrackParCovFromKFP(const KFParticle& kfParticle, const o2::track::PID pid, const int sign)
  {
    o2::gpu::gpustd::array<float, 3> xyz, pxpypz;
    o2::gpu::gpustd::array<float, 21> cv;

    // get parameters from kfParticle
    xyz[0] = kfParticle.GetX();
    xyz[1] = kfParticle.GetY();
    xyz[2] = kfParticle.GetZ();
    pxpypz[0] = kfParticle.GetPx();
    pxpypz[1] = kfParticle.GetPy();
    pxpypz[2] = kfParticle.GetPz();

    // set covariance matrix elements (lower triangle)
    for (int i = 0; i < 21; i++) {
      cv[i] = kfParticle.GetCovariance(i);
    }

    // create TrackParCov track
    o2::track::TrackParCov track = o2::track::TrackParCov(xyz, pxpypz, cv, sign, true, pid);
    return track;
  }

  void processKfIdxCombinatorics(SelectedCollisions const& collisions,
                                 aod::BCsWithTimestamps const& bcWithTimeStamps,
                                 MyTracks const& tracks,
                                 FilteredHfTrackAssocSel const& trackIndices,
                                 MyKfCascTable const& kfcascades,
                                 MyV0Table const&,
                                 aod::V0sLinked const&)
  {

    for (const auto& collision : collisions) {

      // set the magnetic field from CCDB
      auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
	  initCCDB(bc, mRunNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
      auto magneticField = o2::base::Propagator::Instance()->getNominalBz(); // z component
      df.setBz(magneticField);
      KFParticle::SetField(magneticField);
	  
      // loop over cascades reconstructed by cascadebuilder.cxx
      auto thisCollId = collision.globalIndex();
      auto groupedCascades = kfcascades.sliceBy(KFCascadesPerCollision, thisCollId);
      hCollision->Fill(1);

      for (const auto& casc : groupedCascades) {

        //----------------accessing particles in the decay chain-------------
        // Get cascade daughter track
        // cascade daughter - V0
        int v0index = casc.cascadeId();
        // cascade daughter - charged particle
        auto trackOmegaDauCharged = casc.bachelor_as<MyTracks>(); // kaon <- omega track from MyTracks table //
        // V0 positive daughter
        auto trackV0Dau0 = casc.posTrack_as<MyTracks>(); // p <- V0 track (positive track) from MyTracks table
        // V0 negative daughter
        auto trackV0Dau1 = casc.negTrack_as<MyTracks>(); // pion <- V0 track (negative track) from MyTracks table

        auto bachCharge = trackOmegaDauCharged.signed1Pt() > 0 ? +1 : -1;

        //-------------------------- V0 info---------------------------
        // pseudorapidity
        double pseudorapV0Dau0 = trackV0Dau0.eta();
        double pseudorapV0Dau1 = trackV0Dau1.eta();
        //
        //// pion & p TrackParCov
        auto trackParCovV0Dau0 = getTrackParCov(trackV0Dau0);
        auto trackParCovV0Dau1 = getTrackParCov(trackV0Dau1);
        // kaon <- casc TrackParCov
        auto omegaDauChargedTrackParCov = getTrackParCov(trackOmegaDauCharged);
		//
        float massBachelorKaon = o2::constants::physics::MassKaonCharged;
        float massNegTrack,massPosTrack;
        if (bachCharge < 0) {
          massPosTrack = o2::constants::physics::MassProton;
          massNegTrack = o2::constants::physics::MassPionCharged;
        } else {
          massPosTrack = o2::constants::physics::MassPionCharged;
          massNegTrack = o2::constants::physics::MassProton;
        }
		
        //__________________________________________
        //*>~<* step 1 : V0 with dca fitter, uses material corrections implicitly
        // This is optional - move close to minima and therefore take material
        if (doDCAfitter) {
            int nCand = 0;
            try {
              nCand = df.process(trackParCovV0Dau0, trackParCovV0Dau1);
            } catch (...) {
              LOG(error) << "Exception caught in V0 DCA fitter process call!";
              continue;
            }
            if (nCand == 0) {
              continue;
            }
            // re-acquire from DCA fitter

            trackParCovV0Dau0 = df.getTrack(0);
            trackParCovV0Dau1 = df.getTrack(1);
        }
        //__________________________________________
        //*>~<* step 2 : construct V0 with KF
        // Create KF particle

        KFParticle kfpPos = createKFParticleFromTrackParCov(trackParCovV0Dau0, trackParCovV0Dau0.getCharge(), massPosTrack);
        KFParticle kfpNeg = createKFParticleFromTrackParCov(trackParCovV0Dau1, trackParCovV0Dau1.getCharge(), massNegTrack);
		
        const KFParticle* V0Daughters[2] = {&kfpPos, &kfpNeg};
        // construct V0
        KFParticle KFV0;
        KFV0.SetConstructMethod(kfConstructMethod);
        try {
          KFV0.Construct(V0Daughters, 2);
        } catch (std::runtime_error& e) {
          LOG(debug) << "Failed to construct cascade V0 from daughter tracks: " << e.what();
          continue;
        }
        registry.fill(HIST("hKFParticleV0Chi2"), KFV0.GetChi2());
        // mass window cut on lambda before mass constraint
        float massLam, sigLam;
        KFV0.GetMass(massLam, sigLam);
        if (TMath::Abs(massLam - MassLambda0) > lambdaMassWindow)
          continue;
		
        if (kfUseV0MassConstraint) {
          KFV0.SetNonlinearMassConstraint(o2::constants::physics::MassLambda);
        }
        KFParticle KFV0_m = KFV0;
        KFV0_m.SetNonlinearMassConstraint(o2::constants::physics::MassLambda);
        // check
        if ((KFV0.GetNDF() <= 0 || KFV0.GetChi2() <= 0))
          continue;
        if (sigLam <= 0)
          continue;

        // V0 constructed, now recovering TrackParCov for dca fitter minimization (with material correction)
        o2::track::TrackParCov v0TrackParCov = getTrackParCovFromKFP(KFV0, o2::track::PID::Lambda, 0);
        v0TrackParCov.setAbsCharge(0);
        KFV0.TransportToDecayVertex();

        if (doDCAfitter) {
          //__________________________________________
          //*>~<* step 3 : Cascade with dca fitter
          // OmegaMinus with DCAFitter
            int nCandCascade = 0;
            try {
              nCandCascade = df.process(v0TrackParCov, omegaDauChargedTrackParCov);
            } catch (...) {
              LOG(error) << "Exception caught in DCA fitter process call!";
              continue;
            }
            if (nCandCascade == 0)
              continue;

            v0TrackParCov = df.getTrack(0);
            omegaDauChargedTrackParCov = df.getTrack(1);
        }
		
        //__________________________________________
        //*>~<* step 4 : reconstruc cascade(Omega) with KF
        // create KF particle for V0 and bachelor
        KFParticle kfpV0 = KFV0;
        if (doDCAfitter) {
          kfpV0 = createKFParticleFromTrackParCov(v0TrackParCov, 0, massLam);
        }
        KFParticle kfpBachKaon = createKFParticleFromTrackParCov(omegaDauChargedTrackParCov, bachCharge, massBachelorKaon);
        const KFParticle* OmegaDaugthers[2] = {&kfpBachKaon, &kfpV0};

        //// construct cascade
        KFParticle KFOmega;
        KFOmega.SetConstructMethod(kfConstructMethod);
        try {
          KFOmega.Construct(OmegaDaugthers, 2);
        } catch (std::runtime_error& e) {
          LOG(debug) << "Failed to construct omega from V0 and bachelor track: " << e.what();
          continue;
        }

        float massCasc, sigCasc;
        KFOmega.GetMass(massCasc, sigCasc);

        if (kfUseCascadeMassConstraint) {
          // set mass constraint if requested
          // WARNING: this is only adequate for decay chains, i.e. XiC -> Xi or OmegaC -> Omega
          KFOmega.SetNonlinearMassConstraint(o2::constants::physics::MassOmegaMinus);
        }

        KFParticle KFOmega_m = KFOmega;
        KFOmega_m.SetNonlinearMassConstraint(o2::constants::physics::MassOmegaMinus);

        // get DCA of updated daughters at vertex
        KFOmega.TransportToDecayVertex();

        registry.fill(HIST("hInvMassOmegaMinus"), massCasc);
        registry.fill(HIST("hKFParticleCascChi2"), KFOmega.GetChi2());

        omegaDauChargedTrackParCov = getTrackParCovFromKFP(kfpBachKaon, o2::track::PID::Kaon, bachCharge);
        auto trackCasc = getTrackParCovFromKFP(KFOmega, o2::track::PID::OmegaMinus, bachCharge);
		
        omegaDauChargedTrackParCov.setAbsCharge(1);
        omegaDauChargedTrackParCov.setPID(o2::track::PID::Kaon);
        trackCasc.setAbsCharge(1);
        trackCasc.setPID(o2::track::PID::OmegaMinus);

        // info from KFParticle
        std::array<float, 3> pVecV0 = {KFV0.GetPx(), KFV0.GetPy(), KFV0.GetPz()}; // pVec stands for vector containing the 3-momentum components
        std::array<float, 3> vertexV0 = {KFV0.GetX(), KFV0.GetY(), KFV0.GetZ()};
        std::array<float, 3> pVecV0Dau0 = {kfpPos.GetPx(), kfpPos.GetPy(), kfpPos.GetPz()};
        std::array<float, 3> pVecV0Dau1 = {kfpNeg.GetPx(), kfpNeg.GetPy(), kfpNeg.GetPz()};

        //-----------------------------reconstruct cascade track-----------------------------
        // pseudorapidity
        double pseudorapKaFromCas = trackOmegaDauCharged.eta();

        // info from KFParticle 
        std::array<float, 3> vertexCasc = {KFOmega.GetX(), KFOmega.GetY(), KFOmega.GetZ()};
        std::array<float, 3> pVecCasc = {KFOmega.GetPx(), KFOmega.GetPy(), KFOmega.GetPz()};
		
        std::array<float, 3> pVecPionFromCasc = {casc.pxbach(), casc.pybach(), casc.pzbach()};

        //-------------------combining cascade and pion tracks--------------------------
        auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
        for (const auto& trackIndexPion : groupedTrackIndices) {

          // use bachelor selections from HfTrackIndexSkimCreatorTagSelTracks --> bit =2 is CandidateType::CandV0bachelor
          if (!TESTBIT(trackIndexPion.isSelProng(), 2)) {
            continue;
          }
          // OmegaC0 daughter track
          auto trackPion = trackIndexPion.track_as<MyTracks>(); //

          // ask for opposite sign daughters (charm baryon daughters)
          if (trackPion.sign() * trackOmegaDauCharged.sign() >= 0) {
            continue;
          }

          // pseudorapidity
          double pseudorapPiFromCharmBaryon = trackPion.eta();

          // charm bachelor pion track to be processed with DCAFitter
          o2::track::TrackParCov trackParVarPi = getTrackParCov(trackPion);
		  //__________________________________________
          //*>~<* step 5 : OmegaC0 with dca fitter
          // OmegaC0 with DCAFitter
          if (doDCAfitter) {
              // reconstruct charm baryon with DCAFitter
              int nVtxFromFitterCharmBaryon = 0;
              try {
                nVtxFromFitterCharmBaryon = df.process(trackCasc, trackParVarPi);
              } catch (...) {
                LOG(error) << "Exception caught in charm DCA fitter process call!";
                hFitterStatus->Fill(1);
                continue;
              }
              if (nVtxFromFitterCharmBaryon == 0) {
                continue;
              }
              trackCasc = df.getTrack(0);
              trackParVarPi = df.getTrack(1);
          }
          hFitterStatus->Fill(0);

          //__________________________________________
          //*>~<* step 5 : reconstruc Omegac0 with KF
          // create KF particle for V0 and bachelor
          KFParticle kfpCasc = KFOmega;
          if (doDCAfitter) {
            KFParticle kfpCasc = createKFParticleFromTrackParCov(trackCasc, trackCasc.getCharge(), massCasc);
          }
          KFParticle kfpBachPion = createKFParticleFromTrackParCov(trackParVarPi, trackParVarPi.getCharge(), o2::constants::physics::MassPionCharged);
          const KFParticle* OmegaC0Daugthers[2] = {&kfpBachPion, &kfpCasc};

          // construct OmegaC0
          KFParticle KFOmegaC0;
          KFOmegaC0.SetConstructMethod(kfConstructMethod);
          try {
            KFOmegaC0.Construct(OmegaC0Daugthers, 2);
          } catch (std::runtime_error& e) {
            LOG(debug) << "Failed to construct OmegaC0 from V0 and bachelor track: " << e.what();
            continue;
          }

          float massOmegaC0, sigOmegaC0;
          KFOmegaC0.GetMass(massOmegaC0, sigOmegaC0);
          registry.fill(HIST("hKFParticleOmegaC0Chi2"), KFOmegaC0.GetChi2());

          // computing info
          std::array<float, 3> pVecCascAsD;
          std::array<float, 3> pVecPionFromCharmBaryon;
		 
		  pVecPionFromCharmBaryon[0] =  kfpBachPion.GetPx();
		  pVecPionFromCharmBaryon[1] =  kfpBachPion.GetPy();
		  pVecPionFromCharmBaryon[2] =  kfpBachPion.GetPz();
		  pVecCascAsD[0] =  kfpCasc.GetPx();
		  pVecCascAsD[1] =  kfpCasc.GetPy();
		  pVecCascAsD[2] =  kfpCasc.GetPz();
		  if (doDCAfitter){
		    df.propagateTracksToVertex();
            df.getTrack(0).getPxPyPzGlo(pVecCascAsD);
            df.getTrack(1).getPxPyPzGlo(pVecPionFromCharmBaryon);
		  }

          float dcaxyV0Dau0 = trackV0Dau0.dcaXY();
          float dcaxyV0Dau1 = trackV0Dau1.dcaXY();
          float dcaxyKaFromCasc = trackOmegaDauCharged.dcaXY();

          // DCAz (computed with propagateToDCABxByBz method)
          float dcazV0Dau0 = trackV0Dau0.dcaZ();
          float dcazV0Dau1 = trackV0Dau1.dcaZ();
          float dcazKaFromCasc = trackOmegaDauCharged.dcaZ();

          // primary vertex of the collision
          auto primaryVertex = getPrimaryVertex(collision); // get the associated covariance matrix with auto covMatrixPV = primaryVertex.getCov();
          std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};

          // KFPVertex
          KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
          if (doPvRefit && ((trackPion.pvRefitSigmaX2() != 1e10f) || (trackPion.pvRefitSigmaY2() != 1e10f) || (trackPion.pvRefitSigmaZ2() != 1e10f))) { // if I asked for PV refit in trackIndexSkimCreator.cxx
                                                                                                                                                        // KF PVrefit
            kfpVertex.SetXYZ(trackPion.pvRefitX(), trackPion.pvRefitY(), trackPion.pvRefitZ());
            kfpVertex.SetCovarianceMatrix(trackPion.pvRefitSigmaX2(), trackPion.pvRefitSigmaXY(), trackPion.pvRefitSigmaY2(), trackPion.pvRefitSigmaXZ(), trackPion.pvRefitSigmaYZ(), trackPion.pvRefitSigmaZ2());
            
			// o2::dataformats::VertexBase Pvtx;
            primaryVertex.setX(trackPion.pvRefitX());
            primaryVertex.setY(trackPion.pvRefitY());
            primaryVertex.setZ(trackPion.pvRefitZ());
            primaryVertex.setCov(trackPion.pvRefitSigmaX2(), trackPion.pvRefitSigmaXY(), trackPion.pvRefitSigmaY2(), trackPion.pvRefitSigmaXZ(), trackPion.pvRefitSigmaYZ(), trackPion.pvRefitSigmaZ2());

            o2::dataformats::DCA impactParameterV0Dau0;
            o2::dataformats::DCA impactParameterV0Dau1;
            o2::dataformats::DCA impactParameterKaFromCasc;
            o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0Dau0, 2.f, matCorr, &impactParameterV0Dau0);
            o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0Dau1, 2.f, matCorr, &impactParameterV0Dau1);
            o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, omegaDauChargedTrackParCov, 2.f, matCorr, &impactParameterKaFromCasc);
            dcaxyV0Dau0 = impactParameterV0Dau0.getY();
            dcaxyV0Dau1 = impactParameterV0Dau1.getY();
            dcaxyKaFromCasc = impactParameterKaFromCasc.getY();
            dcazV0Dau0 = impactParameterV0Dau0.getZ();
            dcazV0Dau1 = impactParameterV0Dau1.getZ();
            dcazKaFromCasc = impactParameterKaFromCasc.getZ();
          }
          KFParticle KFPV(kfpVertex);
		  
		  pvCoord = {KFPV.GetX(), KFPV.GetY(), KFPV.GetZ()};
		  std::array<float, 3> pVecCharmBaryon = {pVecCascAsD[0] + pVecPionFromCharmBaryon[0], pVecCascAsD[1] + pVecPionFromCharmBaryon[1], pVecCascAsD[2] + pVecPionFromCharmBaryon[2]};
          std::array<float, 3> coordVtxCharmBaryon = {KFOmegaC0.GetX(), KFOmegaC0.GetY(), KFOmegaC0.GetZ()};
          auto covVtxCharmBaryon = KFOmegaC0.CovarianceMatrix();
		  float covMatrixPV[6];
          kfpVertex.GetCovarianceMatrix(covMatrixPV);
		  
          // impact parameters
          o2::dataformats::DCA impactParameterCasc;
          o2::dataformats::DCA impactParameterPiFromCharmBaryon;
          o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackCasc, 2.f, matCorr, &impactParameterCasc);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarPi, 2.f, matCorr, &impactParameterPiFromCharmBaryon);
          float impactParPiFromCharmBaryonXY = impactParameterPiFromCharmBaryon.getY();
          float impactParPiFromCharmBaryonZ = impactParameterPiFromCharmBaryon.getZ();

          // invariant mass under the hypothesis of particles ID corresponding to the decay chain
          double mLambda = casc.mLambda(); // from LF table, V0 mass under lambda hypothesis
          double mCasc = casc.mOmega();
          const std::array<double, 2> arrMassCharmBaryon = {MassOmegaMinus, MassPiPlus};
          double mCharmBaryon = RecoDecay::m(std::array{pVecCascAsD, pVecPionFromCharmBaryon}, arrMassCharmBaryon);

          // computing cosPA
          double cpaV0 = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
          double cpaCharmBaryon = RecoDecay::cpa(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
          double cpaCasc = RecoDecay::cpa(coordVtxCharmBaryon, vertexCasc, pVecCasc);
          double cpaxyV0 = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);
          double cpaxyCharmBaryon = RecoDecay::cpaXY(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
          double cpaxyCasc = RecoDecay::cpaXY(coordVtxCharmBaryon, vertexCasc, pVecCasc);

          // computing decay length and ctau
          double decLenCharmBaryon = RecoDecay::distance(pvCoord, coordVtxCharmBaryon);
          double decLenCascade = RecoDecay::distance(coordVtxCharmBaryon, vertexCasc);
          double decLenV0 = RecoDecay::distance(vertexCasc, vertexV0);

          double phiCharmBaryon, thetaCharmBaryon;
          getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, coordVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon);
          auto errorDecayLengthCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, thetaCharmBaryon) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon));
          auto errorDecayLengthXYCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, 0.) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, 0.));

          double ctOmegac = RecoDecay::ct(pVecCharmBaryon, decLenCharmBaryon, MassOmegaC0);
          double ctCascade = RecoDecay::ct(pVecCasc, decLenCascade, MassOmegaMinus);
          double ctV0 = RecoDecay::ct(pVecV0, decLenV0, MassLambda0);

          // computing eta
          double pseudorapCharmBaryon = RecoDecay::eta(pVecCharmBaryon);
          double pseudorapCascade = RecoDecay::eta(pVecCasc);
          double pseudorapV0 = RecoDecay::eta(pVecV0);

          // DCA between daughters
          float dcaCascDau = casc.dcacascdaughters();
          float dcaV0Dau = casc.dcaV0daughters();
          float dcaCharmBaryonDau = std::sqrt(df.getChi2AtPCACandidate());

          // set hfFlag
          int hfFlag = 1 << aod::hf_cand_omegac::DecayType::DecayToOmegaPi;

          //// KFParticle table information

          kfOmegac0Candidate.dcaCascDau = casc.dcacascdaughters();
          kfOmegac0Candidate.dcaV0Dau = casc.dcaV0daughters();
          // PID information
          auto trackPrFromV0 = trackV0Dau0;
          auto trackPiFromV0 = trackV0Dau1;
          auto trackKaFromCasc = trackOmegaDauCharged;
          auto trackPiFromOmegac = trackPion;

          if (bachCharge > 0) {
            trackPrFromV0 = trackV0Dau1;
            trackPiFromV0 = trackV0Dau0;
          }

          kfOmegac0Candidate.nSigmaTPCPiFromV0 = trackPiFromV0.tpcNSigmaPi();
          kfOmegac0Candidate.nSigmaTPCPrFromV0 = trackPrFromV0.tpcNSigmaPr();
          kfOmegac0Candidate.nSigmaTOFKaFromCasc = trackKaFromCasc.tofNSigmaKa();
          kfOmegac0Candidate.nSigmaTPCKaFromCasc = trackKaFromCasc.tpcNSigmaKa();
          kfOmegac0Candidate.nSigmaTOFPiFromOmegac = trackPiFromOmegac.tofNSigmaPi();
          kfOmegac0Candidate.nSigmaTPCPiFromOmegac = trackPiFromOmegac.tpcNSigmaPi();

          // KFParticle vriables
          //  KFparticle to decay vertext

          KFParticle kfpNegToV0 = kfpNeg;
          KFParticle kfpPosToV0 = kfpPos;
          kfpNegToV0.SetProductionVertex(KFV0);
          kfpPosToV0.SetProductionVertex(KFV0);

          KFParticle kfpBachKaonToOmega = kfpBachKaon;
          KFParticle kfpV0ToCasc = kfpV0;
          kfpBachKaonToOmega.SetProductionVertex(KFOmega);
          kfpV0ToCasc.SetProductionVertex(KFOmega);

          KFParticle kfpCascToOmegaC = kfpCasc;
          KFParticle kfpBachPionToOmegaC = kfpBachPion;
          kfpBachPionToOmegaC.SetProductionVertex(KFOmegaC0);
          kfpCascToOmegaC.SetProductionVertex(KFOmegaC0);

          // KFParticle to PV
          KFParticle kfpV0ToPv = kfpV0;
          KFParticle kfpCascToPv = kfpCasc;
          KFParticle kfpOmegacToPv = KFOmegaC0;
          KFParticle kfpPiFromOmegacToPv = kfpBachPion;

          kfpV0ToPv.SetProductionVertex(KFPV);
          kfpCascToPv.SetProductionVertex(KFPV);
          kfpOmegacToPv.SetProductionVertex(KFPV);
          kfpPiFromOmegacToPv.SetProductionVertex(KFPV);

          // KF geochi2
          kfOmegac0Candidate.chi2GeoV0 = KFV0.GetChi2();
          auto v0NDF = KFV0.GetNDF();
          auto v0Chi2OverNdf = kfOmegac0Candidate.chi2GeoV0 / v0NDF;

          kfOmegac0Candidate.chi2GeoCasc = KFOmega.GetChi2();
          auto cascNDF = KFOmega.GetNDF();
          auto cascChi2OverNdf = kfOmegac0Candidate.chi2GeoCasc / cascNDF;

          kfOmegac0Candidate.chi2GeoOmegac = KFOmegaC0.GetChi2();
          auto charmbaryonNDF = KFOmegaC0.GetNDF();
          auto charmbaryonChi2OverNdf = kfOmegac0Candidate.chi2GeoOmegac / charmbaryonNDF;

          kfOmegac0Candidate.chi2MassV0 = KFV0_m.GetChi2();
          auto v0NDF_m = KFV0_m.GetNDF();
          auto v0Chi2OverNdf_m = kfOmegac0Candidate.chi2MassV0 / v0NDF_m;

          kfOmegac0Candidate.chi2MassCasc = KFOmega_m.GetChi2();
          auto cascNDF_m = KFOmega_m.GetNDF();
          auto cascChi2OverNdf_m = kfOmegac0Candidate.chi2MassCasc / cascNDF_m;

          // KF topo Chi2
          kfOmegac0Candidate.chi2TopoV0ToPv = kfpV0ToPv.GetChi2();
          kfOmegac0Candidate.chi2TopoCascToPv = kfpCascToPv.GetChi2();
          kfOmegac0Candidate.chi2TopoPiFromOmegacToPv = kfpPiFromOmegacToPv.GetChi2();
          kfOmegac0Candidate.chi2TopoOmegacToPv = kfpOmegacToPv.GetChi2();

          auto cascBachTopoChi2 = kfpBachKaonToOmega.GetChi2();
          kfOmegac0Candidate.chi2TopoV0ToCasc = kfpV0ToCasc.GetChi2();
          kfOmegac0Candidate.chi2TopoCascToOmegac = kfpCascToOmegaC.GetChi2();

          // KF ldl
          kfOmegac0Candidate.ldlV0 = ldlFromKF(KFV0, KFPV);
          kfOmegac0Candidate.ldlCasc = ldlFromKF(KFOmega, KFPV);
          kfOmegac0Candidate.ldlOmegac = ldlFromKF(KFOmegaC0, KFPV);

          // KF dca
          kfOmegac0Candidate.kfDcaXYPiFromOmegac = kfpBachPion.GetDistanceFromVertexXY(KFPV);
          kfOmegac0Candidate.kfDcaCascDau = kfpBachKaon.GetDistanceFromParticle(kfpV0);
          kfOmegac0Candidate.kfDcaOmegacDau = kfpBachPion.GetDistanceFromParticle(kfpCasc);
          kfOmegac0Candidate.kfDcaXYCascToPv = kfpCasc.GetDistanceFromVertexXY(KFPV);

          // KF decay length
          float DecayLxy_Lam, err_DecayLxy_Lam;
          kfpV0ToCasc.GetDecayLengthXY(DecayLxy_Lam, err_DecayLxy_Lam);
          kfOmegac0Candidate.decayLenXYLambda = DecayLxy_Lam;

          float DecayLxy_Casc, err_DecayLxy_Casc;
          kfpCascToPv.GetDecayLengthXY(DecayLxy_Casc, err_DecayLxy_Casc);
          kfOmegac0Candidate.decayLenXYCasc = DecayLxy_Casc;

          float DecayLxy_Omegac0, err_DecayLxy_Omegac0;
          kfpOmegacToPv.GetDecayLengthXY(DecayLxy_Omegac0, err_DecayLxy_Omegac0);
          kfOmegac0Candidate.decayLenXYOmegac = DecayLxy_Omegac0;

          // KF cosPA
          kfOmegac0Candidate.cosPaV0ToPv = cpaFromKF(kfpV0, KFPV);
          kfOmegac0Candidate.cosPaCascToPv = cpaFromKF(kfpCasc, KFPV);
          kfOmegac0Candidate.cosPaOmegacToPv = cpaFromKF(KFOmegaC0, KFPV);

          kfOmegac0Candidate.cosPaV0ToCasc = cpaFromKF(kfpV0, kfpCasc);
          kfOmegac0Candidate.cosPaCascToOmegac = cpaFromKF(kfpCasc, KFOmegaC0);

          // KF mass
          kfOmegac0Candidate.massV0 = massLam;
          kfOmegac0Candidate.massCasc = massCasc;
          kfOmegac0Candidate.massOmegac = massOmegaC0;

          // KF pT
          kfOmegac0Candidate.ptPiFromOmegac = trackPion.pt();
          kfOmegac0Candidate.ptOmegac = kfpOmegacToPv.GetPt();

          // KF rapidity
          kfOmegac0Candidate.rapOmegac = kfpOmegacToPv.GetRapidity();

          // KF cosThetaStar
          kfOmegac0Candidate.cosThetaStarPiFromOmegac = cosThetaStarFromKF(0, 4332, 211, 3312, kfpBachPionToOmegaC, kfpCascToOmegaC);

          // KF ct
          float ct_Omegac0 = 0., err_ct_Omegac0 = 0.;
          kfpOmegacToPv.GetLifeTime(ct_Omegac0, err_ct_Omegac0);
          kfOmegac0Candidate.ctOmegac = ct_Omegac0;

          // KF eta
          kfOmegac0Candidate.etaOmegac = kfpOmegacToPv.GetEta();

          // fill KF hist
          registry.fill(HIST("hKFParticleCascBachTopoChi2"), cascBachTopoChi2);
          registry.fill(HIST("hKFParticleV0TopoChi2"), kfOmegac0Candidate.chi2TopoV0ToCasc);
          registry.fill(HIST("hKFParticleCascTopoChi2"), kfOmegac0Candidate.chi2TopoCascToOmegac);

          registry.fill(HIST("hKfLambda_ldl"), kfOmegac0Candidate.ldlV0);
          registry.fill(HIST("hKfOmega_ldl"), kfOmegac0Candidate.ldlCasc);
          registry.fill(HIST("hKfOmegaC0_ldl"), kfOmegac0Candidate.ldlOmegac);
          registry.fill(HIST("hDcaXYCascadeToPVKf"), kfOmegac0Candidate.kfDcaXYCascToPv);

          // fill test histograms
          hInvMassCharmBaryon->Fill(mCharmBaryon);
          hKfInvMassCharmBaryon->Fill(massOmegaC0);
          // fill the table
          rowCandidate(collision.globalIndex(),
                       pvCoord[0], pvCoord[1], pvCoord[2],
                       vertexCasc[0], vertexCasc[1], vertexCasc[2],
                       vertexV0[0], vertexV0[1], vertexV0[2],
                       trackOmegaDauCharged.sign(),
                       covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
                       pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
                       pVecCasc[0], pVecCasc[1], pVecCasc[2],
                       pVecPionFromCharmBaryon[0], pVecPionFromCharmBaryon[1], pVecPionFromCharmBaryon[2],
                       pVecV0[0], pVecV0[1], pVecV0[2],
                       pVecPionFromCasc[0], pVecPionFromCasc[1], pVecPionFromCasc[2],
                       pVecV0Dau0[0], pVecV0Dau0[1], pVecV0Dau0[2],
                       pVecV0Dau1[0], pVecV0Dau1[1], pVecV0Dau1[2],
                       impactParameterCasc.getY(), impactParPiFromCharmBaryonXY,
                       impactParameterCasc.getZ(), impactParPiFromCharmBaryonZ,
                       std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterPiFromCharmBaryon.getSigmaY2()),
                       v0index, casc.posTrackId(), casc.negTrackId(),
                       casc.globalIndex(), trackPion.globalIndex(), trackOmegaDauCharged.globalIndex(),
                       mLambda, mCasc, mCharmBaryon,
                       cpaV0, cpaCharmBaryon, cpaCasc, cpaxyV0, cpaxyCharmBaryon, cpaxyCasc,
                       ctOmegac, ctCascade, ctV0,
                       pseudorapV0Dau0, pseudorapV0Dau1, pseudorapKaFromCas, pseudorapPiFromCharmBaryon,
                       pseudorapCharmBaryon, pseudorapCascade, pseudorapV0,
                       dcaxyV0Dau0, dcaxyV0Dau1, dcaxyKaFromCasc,
                       dcazV0Dau0, dcazV0Dau1, dcazKaFromCasc,
                       dcaCascDau, dcaV0Dau, dcaCharmBaryonDau,
                       decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon,
                       hfFlag);

          kfCandidateData(kfOmegac0Candidate.nSigmaTPCPiFromOmegac, kfOmegac0Candidate.nSigmaTOFPiFromOmegac,
                          kfOmegac0Candidate.nSigmaTOFKaFromCasc, kfOmegac0Candidate.nSigmaTOFKaFromCasc,
                          kfOmegac0Candidate.nSigmaTPCPiFromV0, kfOmegac0Candidate.nSigmaTPCPrFromV0,
                          kfOmegac0Candidate.kfDcaXYPiFromOmegac, kfOmegac0Candidate.kfDcaCascDau, kfOmegac0Candidate.kfDcaOmegacDau, kfOmegac0Candidate.kfDcaXYCascToPv,
                          kfOmegac0Candidate.chi2GeoV0, kfOmegac0Candidate.chi2GeoCasc, kfOmegac0Candidate.chi2GeoOmegac, kfOmegac0Candidate.chi2MassV0, kfOmegac0Candidate.chi2MassCasc,
                          kfOmegac0Candidate.ldlV0, kfOmegac0Candidate.ldlCasc, kfOmegac0Candidate.ldlOmegac,
                          kfOmegac0Candidate.chi2TopoV0ToPv, kfOmegac0Candidate.chi2TopoCascToPv, kfOmegac0Candidate.chi2TopoPiFromOmegacToPv, kfOmegac0Candidate.chi2TopoOmegacToPv,
                          kfOmegac0Candidate.chi2TopoV0ToCasc, kfOmegac0Candidate.chi2TopoCascToOmegac,
                          kfOmegac0Candidate.decayLenXYLambda, kfOmegac0Candidate.decayLenXYCasc, kfOmegac0Candidate.decayLenXYOmegac,
                          kfOmegac0Candidate.cosPaV0ToCasc, kfOmegac0Candidate.cosPaV0ToPv, kfOmegac0Candidate.cosPaCascToOmegac, kfOmegac0Candidate.cosPaCascToPv, kfOmegac0Candidate.cosPaOmegacToPv,
                          kfOmegac0Candidate.massV0, kfOmegac0Candidate.massCasc, kfOmegac0Candidate.massOmegac,
                          kfOmegac0Candidate.rapOmegac, kfOmegac0Candidate.ptPiFromOmegac, kfOmegac0Candidate.ptOmegac,
                          kfOmegac0Candidate.cosThetaStarPiFromOmegac,
                          kfOmegac0Candidate.ctOmegac,
                          kfOmegac0Candidate.etaOmegac,
                          v0NDF, cascNDF, charmbaryonNDF, v0NDF_m, cascNDF_m,
                          v0Chi2OverNdf, cascChi2OverNdf, charmbaryonChi2OverNdf, v0Chi2OverNdf_m, cascChi2OverNdf_m);
        } // loop over pions
      } // loop over cascades
    }   // close loop collisions

  } // end of process
  PROCESS_SWITCH(HfCandidateCreatorOmegacToOmegaPiWithKfp, processKfIdxCombinatorics, "Do index combinatorics with KFParticle", true);
}; // end of struct

/// Performs MC matching.
struct HfCandidateCreatorOmegacToOmegaPiWithKfpMc {
  Produces<aod::HfOmegacMCRec> rowMCMatchRec;
  Produces<aod::HfOmegacMCGen> rowMCMatchGen;

  Configurable<bool> matchOmegacMc{"matchOmegacMc", false, "Do MC matching for Omegac0"};

  void init(InitContext const&) {}

  void processDoNoMc(aod::Collisions::iterator const& collision)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(HfCandidateCreatorOmegacToOmegaPiWithKfpMc, processDoNoMc, "Do not run MC process function", true);

  void processMc(aod::HfCandOmegac const& candidates,
                 aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles,
                 aod::McCollisionLabels const&)
  {
    float ptCharmBaryonGen = -999.;
    float etaCharmBaryonGen = -999.;
    int indexRec = -1;
    int indexRecCharmBaryon = -1;
    int8_t sign = -9;
    int8_t flag = 0;
    int8_t origin = 0; // to be used for prompt/non prompt
    int8_t debug = 0;
    int8_t debugGenCharmBar = 0;
    int8_t debugGenOmega = 0;
    int8_t debugGenLambda = 0;

    //  Match reconstructed candidates.
    for (const auto& candidate : candidates) {
      flag = 0;
      origin = RecoDecay::OriginType::None;
      debug = 0;
      bool collisionMatched = false;
      auto arrayDaughters = std::array{candidate.piFromCharmBaryon_as<aod::TracksWMc>(), // pi <- charm baryon
                                       candidate.bachelor_as<aod::TracksWMc>(),          // Kaon <- cascade
                                       candidate.posTrack_as<aod::TracksWMc>(),          // p <- lambda
                                       candidate.negTrack_as<aod::TracksWMc>()};         // pi <- lambda
      auto arrayDaughtersCasc = std::array{candidate.bachelor_as<aod::TracksWMc>(),
                                           candidate.posTrack_as<aod::TracksWMc>(),
                                           candidate.negTrack_as<aod::TracksWMc>()};
      auto arrayDaughtersV0 = std::array{candidate.posTrack_as<aod::TracksWMc>(),
                                         candidate.negTrack_as<aod::TracksWMc>()};

      // Omegac matching
      if (matchOmegacMc) {
        // Omegac → pi pi pi p
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, int(kOmegaC0), std::array{int(kPiPlus), int(kKMinus), int(kProton), int(kPiMinus)}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Omega- → kaon pi p
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersCasc, int(kOmegaMinus), std::array{int(kKMinus), int(kProton), int(kPiMinus)}, true, &sign, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersV0, int(kLambda0), std::array{int(kProton), int(kPiMinus)}, true, &sign, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << aod::hf_cand_omegac::DecayType::OmegaczeroToOmegaPi);
              collisionMatched = candidate.collision_as<aod::McCollisionLabels>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
            }
          }
        }

        // Check whether the charm baryon is non-prompt (from a b quark).
        if (flag != 0) {
          auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, true);
        }
      }
      if (debug == 2 || debug == 3) {
        LOGF(info, "WARNING: Charm baryon decays in the expected final state but the condition on the intermediate states are not fulfilled");
      }
      rowMCMatchRec(flag, debug, origin, collisionMatched);

    } // close loop over candidates

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      ptCharmBaryonGen = -999.;
      etaCharmBaryonGen = -999.;
      flag = 0;
      sign = -9;
      debugGenCharmBar = 0;
      debugGenOmega = 0;
      debugGenLambda = 0;
      origin = RecoDecay::OriginType::None;
      if (matchOmegacMc) {
        //  Omegac → Omega pi
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, int(kOmegaC0), std::array{int(kOmegaMinus), int(kPiPlus)}, true, &sign)) {
          debugGenCharmBar = 1;
          ptCharmBaryonGen = particle.pt();
          etaCharmBaryonGen = particle.eta();
          // Omega -> Lambda kaon
          auto cascMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(mcParticles, cascMC, int(kOmegaMinus), std::array{int(kLambda0), int(kKMinus)}, true)) {
            debugGenOmega = 1;
            // Lambda -> p pi
            auto v0MC = mcParticles.rawIteratorAt(cascMC.daughtersIds().front());
            if (RecoDecay::isMatchedMCGen(mcParticles, v0MC, int(kLambda0), std::array{int(kProton), int(kPiMinus)}, true)) {
              debugGenLambda = 1;
              flag = sign * (1 << aod::hf_cand_omegac::DecayType::OmegaczeroToOmegaPi);
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark)
        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, true);
        }
      }

      rowMCMatchGen(flag, debugGenCharmBar, debugGenOmega, debugGenLambda, ptCharmBaryonGen, etaCharmBaryonGen, origin);
    }
  }                                                                                          // close process
  PROCESS_SWITCH(HfCandidateCreatorOmegacToOmegaPiWithKfpMc, processMc, "Process MC", true); // false
};                                                                                           // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorOmegacToOmegaPiWithKfp>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorOmegacToOmegaPiWithKfpMc>(cfgc)};
}
