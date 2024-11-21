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

/// \brief Builder task for 3-body decay reconstruction (p + pion + bachelor)
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>
/// \author Carolina Reetz <c.reetz@cern.ch> (KFParticle specific part)

#include <cmath>
#include <array>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "PWGLF/DataModel/pidTOFGeneric.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/PID/PIDTOF.h"
#include "TableHelper.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsTPC/BetheBlochAleph.h"

#ifndef HomogeneousField
#define HomogeneousField
#endif

// includes KFParticle
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using FullTracksExtPIDIU = soa::Join<FullTracksExtIU, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe>;

using ColwithEvTimes = o2::soa::Join<aod::Collisions, aod::EvSels, aod::EvTimeTOFFT0>;
using FullCols = o2::soa::Join<ColwithEvTimes, aod::CentFT0Cs>;
using TrackExtIUwithEvTimes = soa::Join<FullTracksExtIU, aod::EvTimeTOFFT0ForTrack>;
using TrackExtPIDIUwithEvTimes = soa::Join<FullTracksExtPIDIU, aod::EvTimeTOFFT0ForTrack>;

using MCLabeledTracksIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

struct vtxCandidate {
  int track0Id;
  int track1Id;
  int track2Id;
  int collisionId;
  int decay3bodyId;
  float vtxPos[3];
  float track0P[3];
  float track1P[3];
  float track2P[3];
  float dcadaughters;
  float daudcaxytopv[3]; // 0 - proton, 1 - pion, 2 - bachelor
  float daudcatopv[3];   // 0 - proton, 1 - pion, 2 - bachelor
  float bachelortofNsigma;
};

struct decay3bodyBuilder {

  Produces<aod::StoredVtx3BodyDatas> vtx3bodydata;
  Produces<aod::KFVtx3BodyDatas> kfvtx3bodydata;
  Produces<aod::KFVtx3BodyDatasLite> kfvtx3bodydatalite;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  std::vector<vtxCandidate> vtxCandidates;

  // Configurables
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};

  enum hyp3body { kH3L = 0,
                  kH4L,
                  kHe4L,
                  kHe5L,
                  kNHyp3body };

  enum vtxstep { kVtxAll = 0,
                 kVtxTPCNcls,
                 kVtxPIDCut,
                 kVtxhasSV,
                 kVtxDcaDau,
                 kVtxCosPA,
                 kNVtxSteps };

  enum kfvtxstep { kKfVtxAll = 0,
                   kKfVtxCollIds,
                   kKfVtxCharge,
                   kKfVtxEta,
                   kKfVtxTPCNcls,
                   kKfVtxTPCRows,
                   kKfVtxTPCPID,
                   kKfVtxDCAxyPV,
                   kKfVtxDCAzPV,
                   kKfVtxV0MassConst,
                   kKfVtxhasSV,
                   kKfVtxDcaDau,
                   kKfVtxDcaDauVtx,
                   kKfVtxDauPt,
                   kKfVtxRap,
                   kKfVtxPt,
                   kKfVtxMass,
                   kKfVtxCosPA,
                   kKfVtxCosPAXY,
                   kKfVtxChi2geo,
                   kKfVtxTopoConstr,
                   kKfVtxChi2topo,
                   kKfNVtxSteps };

  HistogramRegistry registry{"registry", {}};

  // hypothesis
  Configurable<int> motherhyp{"motherhyp", 0, "hypothesis of the 3body decayed particle"}; // corresponds to hyp3body
  int bachelorcharge = 1;                                                                  // to be updated in Init base on the hypothesis
  // o2::aod::pidtofgeneric::TofPidNewCollision<ColwithEvTimes::iterator, TrackExtPIDIUwithEvTimes::iterator> bachelorTOFPID; // to be updated in Init base on the hypothesis
  o2::aod::pidtofgeneric::TofPidNewCollision<TrackExtPIDIUwithEvTimes::iterator> bachelorTOFPID; // to be updated in Init base on the hypothesis

  // Selection criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<int> mintpcNCls{"mintpcNCls", 70, "min tpc Nclusters"};
  Configurable<float> minCosPA3body{"minCosPA3body", 0.9, "minCosPA3body"};
  Configurable<float> dcavtxdau{"dcavtxdau", 1.0, "DCA Vtx Daughters"};
  Configurable<bool> enablePidCut{"enablePidCut", 0, "enable function checkPIDH3L"};
  Configurable<float> TofPidNsigmaMin{"TofPidNsigmaMin", -5, "TofPidNsigmaMin"};
  Configurable<float> TofPidNsigmaMax{"TofPidNsigmaMax", 5, "TofPidNsigmaMax"};
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
  Configurable<float> minBachPUseTOF{"minBachPUseTOF", 1, "minBachP Enable TOF PID"};

  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  // CCDB TOF PID paras
  Configurable<int64_t> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};
  Configurable<std::string> paramFileName{"paramFileName", "", "Path to the parametrization object. If empty the parametrization is not taken from file"};
  Configurable<std::string> parametrizationPath{"parametrizationPath", "TOF/Calib/Params", "Path of the TOF parametrization on the CCDB or in the file, if the paramFileName is not empty"};
  Configurable<std::string> passName{"passName", "", "Name of the pass inside of the CCDB parameter collection. If empty, the automatically deceted from metadata (to be implemented!!!)"};
  Configurable<std::string> timeShiftCCDBPath{"timeShiftCCDBPath", "", "Path of the TOF time shift vs eta. If empty none is taken"};
  Configurable<bool> loadResponseFromCCDB{"loadResponseFromCCDB", false, "Flag to load the response from the CCDB"};
  Configurable<bool> fatalOnPassNotAvailable{"fatalOnPassNotAvailable", true, "Flag to throw a fatal if the pass is not available in the retrieved CCDB object"};
  // for KFParticle reconstruction
  struct : ConfigurableGroup {
    Configurable<bool> fillCandidateLiteTable{"kfparticleConfigurations.fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
    Configurable<bool> doSel8selection{"kfparticleConfigurations.doSel8selection", true, "flag for sel8 event selection"};
    Configurable<bool> doPosZselection{"kfparticleConfigurations.doPosZselection", true, "flag for posZ event selection"};
    Configurable<bool> doDCAFitterPreMinimum{"kfparticleConfigurations.doDCAFitterPreMinimum", false, "do DCAFitter pre-optimization before KF fit to include material corrections for decay3body vertex"};
    Configurable<bool> doTrackQA{"kfparticleConfigurations.doTrackQA", false, "Flag to fill QA histograms for daughter tracks."};
    Configurable<bool> doVertexQA{"kfparticleConfigurations.doVertexQA", false, "Flag to fill QA histograms for KFParticle PV."};
    Configurable<bool> useLambdaMassConstraint{"kfparticleConfigurations.useLambdaMassConstraint", false, "Apply Lambda mass constraint on proton-pion vertex"};
    Configurable<float> maxEta{"kfparticleConfigurations.maxEta", 0.9, "Maximum eta for daughter tracks"};
    Configurable<bool> useTPCforPion{"kfparticleConfigurations.useTPCforPion", true, "Flag to ask for TPC info for pion track (PID, nClusters), false: pion track can be ITS only"};
    Configurable<float> mintpcNClsProton{"kfparticleConfigurations.mintpcNClsProton", 70, "Minimum number of TPC clusters for proton track"};
    Configurable<float> mintpcNClsPion{"kfparticleConfigurations.mintpcNClsPion", 70, "Minimum number of TPC clusters for pion track"};
    Configurable<float> mintpcNClsBach{"kfparticleConfigurations.mintpcNClsBach", 70, "Minimum number of TPC clusters for bachelor track"};
    Configurable<float> mintpcCrossedRows{"kfparticleConfigurations.mintpcCrossedRows", 70, "Minimum number of TPC crossed rows for proton and deuteron track"};
    Configurable<float> mintpcCrossedRowsPion{"kfparticleConfigurations.mintpcCrossedRowsPion", 70, "Minimum number of TPC crossed rows for pion track"};
    Configurable<float> minPtProton{"kfparticleConfigurations.minPtProton", 0.1, "Minimum pT of proton track"};
    Configurable<float> maxPtProton{"kfparticleConfigurations.maxPtProton", 10, "Maximum pT of proton track"};
    Configurable<float> minPtPion{"kfparticleConfigurations.minPtPion", 0.1, "Minimum pT of pion track"};
    Configurable<float> maxPtPion{"kfparticleConfigurations.maxPtPion", 10, "Maximum pT of pion track"};
    Configurable<float> minPtDeuteron{"kfparticleConfigurations.minPtDeuteron", 0.1, "Minimum pT of deuteron track"};
    Configurable<float> maxPtDeuteron{"kfparticleConfigurations.maxPtDeuteron", 10, "Maximum pT of deuteron track"};
    Configurable<float> mindcaXYPionPV{"kfparticleConfigurations.mindcaXYPionPV", 0.1, "Minimum DCA XY of the pion daughter track to the PV"};
    Configurable<float> mindcaXYProtonPV{"kfparticleConfigurations.mindcaXYProtonPV", 0.1, "Minimum DCA XY of the proton daughter track to the PV"};
    Configurable<float> mindcaZPionPV{"kfparticleConfigurations.mindcaZPionPV", 0.1, "Minimum DCA Z of the pion daughter track to the PV"};
    Configurable<float> mindcaZProtonPV{"kfparticleConfigurations.mindcaZProtonPV", 0.1, "Minimum DCA Z of the proton daughter track to the PV"};
    Configurable<float> maxtpcnSigma{"kfparticleConfigurations.maxtpcnSigma", 5., "Maximum nSigma TPC for daughter tracks"};
    Configurable<float> maxDcaProDeu{"kfparticleConfigurations.maxDcaProDeu", 1000., "Maximum geometrical distance between proton and deuteron at the SV in 3D with KFParticle"};
    Configurable<float> maxDcaProPi{"kfparticleConfigurations.maxDcaProPi", 1000., "Maximum geometrical distance between proton and pion at the SV in 3D with KFParticle"};
    Configurable<float> maxDcaPiDe{"kfparticleConfigurations.maxDcaPiDe", 1000., "Maximum geometrical distance between pion and deuteron at the SV in 3D with KFParticle"};
    Configurable<float> maxDcaXYSVDau{"kfparticleConfigurations.maxDcaXYSVDau", 1.0, "Maximum geometrical distance of daughter tracks from the SV in XY with KFParticle"};
    Configurable<float> maxRapidityHt{"kfparticleConfigurations.maxRapidityHt", 1., "Maximum rapidity for Hypertriton candidates with KFParticle"};
    Configurable<float> minPtHt{"kfparticleConfigurations.minPtHt", 0., "Minimum momentum for Hypertriton candidates with KFParticle"};
    Configurable<float> maxPtHt{"kfparticleConfigurations.maxPtHt", 36., "Maximum momentum for Hypertriton candidates with KFParticle"};
    Configurable<float> minMassHt{"kfparticleConfigurations.minMassHt", 2.96, "Minimum candidate mass with KFParticle"};
    Configurable<float> maxMassHt{"kfparticleConfigurations.maxMassHt", 3.05, "Maximum candidate mass with KFParticle"};
    Configurable<float> maxctauHt{"kfparticleConfigurations.maxctauHt", 40., "Maximum candidate ctau with KFParticle before topological constraint"};
    Configurable<float> maxChi2geo{"kfparticleConfigurations.maxChi2geo", 1000., "Maximum chi2 geometrical with KFParticle"};
    Configurable<float> minCosPA{"kfparticleConfigurations.minCosPA", 0.5, "Minimum cosine pointing angle with KFParticle"};
    Configurable<float> minCosPAxy{"kfparticleConfigurations.minCosPAxy", 0.5, "Minimum cosine pointing angle in xy with KFParticle"};
    Configurable<bool> applyTopoSel{"kfparticleConfigurations.applyTopoSel", false, "Apply selection constraining the mother to the PV with KFParticle"};
    Configurable<float> maxChi2topo{"kfparticleConfigurations.maxChi2topo", 1000., "Maximum chi2 topological with KFParticle"};
  } kfparticleConfigurations;

  //------------------------------------------------------------------
  // Sets for event mixing
  struct : ConfigurableGroup {
    Configurable<int> nUseMixedEvent{"nUseMixedEvent", 5, "nUseMixedEvent"};
    Configurable<bool> em_event_sel8_selection{"em_event_sel8_selection", true, "event selection count post sel8 cut"};
    Configurable<float> etacut{"etacut", 0.9, "etacut"};
    Configurable<float> minProtonPt{"minProtonPt", 0.3, "minProtonPt"};
    Configurable<float> maxProtonPt{"maxProtonPt", 5, "maxProtonPt"};
    Configurable<float> minPionPt{"minPionPt", 0.1, "minPionPt"};
    Configurable<float> maxPionPt{"maxPionPt", 1.2, "maxPionPt"};
    Configurable<float> minDeuteronPt{"minDeuteronPt", 0.6, "minDeuteronPt"};
    Configurable<float> maxDeuteronPt{"maxDeuteronPt", 10, "maxDeuteronPt"};
    Configurable<int> mintpcNClsproton{"mintpcNClsproton", 90, "min tpc Nclusters for proton"};
    Configurable<int> mintpcNClspion{"mintpcNClspion", 70, "min tpc Nclusters for pion"};
    Configurable<int> mintpcNClsbachelor{"mintpcNClsbachelor", 100, "min tpc Nclusters for bachelor"};
    Configurable<float> emTpcPidNsigmaCut{"emTpcPidNsigmaCut", 5, "emTpcPidNsigmaCut"};
  } EMTrackSel;

  Preslice<TrackExtPIDIUwithEvTimes> tracksperCol = aod::track::collisionId;
  SliceCache cache;
  ConfigurableAxis axisPosZ{"axisPosZ", {40, -10, 10}, "Mixing bins - posZ"};
  ConfigurableAxis axisCentrality{"axisCentrality", {10, 0, 100}, "Mixing bins - centrality"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  // Filters and slices
  // Filter collisionFilter = (aod::evsel::sel8 == true && nabs(aod::collision::posZ) < 10.f);
  Preslice<aod::Decay3Bodys> perCollision = o2::aod::decay3body::collisionId;

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::vertexing::DCAFitterN<3> fitter3body;
  o2::pid::tof::TOFResoParamsV2 mRespParamsV2;

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  // could be changed later
    maxStep = 2.00f; // could be changed later

    // set hypothesis corresponds to hyp3body, tpcpid to be implemented
    switch (motherhyp) {
      case hyp3body::kH3L:
        bachelorcharge = 1;
        bachelorTOFPID.SetPidType(o2::track::PID::Deuteron);
        break;
      case hyp3body::kH4L:
        bachelorcharge = 1;
        bachelorTOFPID.SetPidType(o2::track::PID::Triton);
        break;
      case hyp3body::kHe4L:
        bachelorcharge = 2;
        bachelorTOFPID.SetPidType(o2::track::PID::Helium3);
        break;
      case hyp3body::kHe5L:
        bachelorcharge = 2;
        bachelorTOFPID.SetPidType(o2::track::PID::Alpha);
        break;
      default:
        LOG(fatal) << "Wrong hypothesis for decay3body";
        return;
    }

    fitter3body.setPropagateToPCA(true);
    fitter3body.setMaxR(200.); //->maxRIni3body
    fitter3body.setMinParamChange(1e-3);
    fitter3body.setMinRelChi2Change(0.9);
    fitter3body.setMaxDZIni(1e9);
    fitter3body.setMaxChi2(1e9);
    fitter3body.setUseAbsDCA(d_UseAbsDCA);

    // Material correction in the DCA fitter
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (useMatCorrType == 1) {
      LOGF(info, "TGeo correction requested, loading geometry");
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
      }
    }
    if (useMatCorrType == 2) {
      LOGF(info, "LUT correction requested, loading LUT");
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    }

    // Material correction in the DCA fitter
    if (useMatCorrType == 1)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    if (useMatCorrType == 2)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

    fitter3body.setMatCorrType(matCorr);

    // Add histograms separately for different process functions
    if (doprocessRun3 == true || doprocessRun3EM == true || doprocessRun3EMLikeSign == true) {
      registry.add("hEventCounter", "hEventCounter", HistType::kTH1F, {{1, 0.0f, 1.0f}});
      auto hVtx3BodyCounter = registry.add<TH1>("hVtx3BodyCounter", "hVtx3BodyCounter", HistType::kTH1F, {{6, 0.0f, 6.0f}});
      hVtx3BodyCounter->GetXaxis()->SetBinLabel(1, "Total");
      hVtx3BodyCounter->GetXaxis()->SetBinLabel(2, "TPCNcls");
      hVtx3BodyCounter->GetXaxis()->SetBinLabel(3, "PIDCut");
      hVtx3BodyCounter->GetXaxis()->SetBinLabel(4, "HasSV");
      hVtx3BodyCounter->GetXaxis()->SetBinLabel(5, "DcaDau");
      hVtx3BodyCounter->GetXaxis()->SetBinLabel(6, "CosPA");
      registry.add("hBachelorTOFNSigmaDe", "", HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}});
    }

    if (doprocessRun3withKFParticle == true) {
      auto hEventCounterKFParticle = registry.add<TH1>("hEventCounterKFParticle", "hEventCounterKFParticle", HistType::kTH1F, {{4, 0.0f, 4.0f}});
      hEventCounterKFParticle->GetXaxis()->SetBinLabel(1, "total");
      hEventCounterKFParticle->GetXaxis()->SetBinLabel(2, "sel8");
      hEventCounterKFParticle->GetXaxis()->SetBinLabel(3, "vertexZ");
      hEventCounterKFParticle->GetXaxis()->SetBinLabel(4, "has candidate");
      hEventCounterKFParticle->LabelsOption("v");
      auto hVtx3BodyCounterKFParticle = registry.add<TH1>("hVtx3BodyCounterKFParticle", "hVtx3BodyCounterKFParticle", HistType::kTH1F, {{22, 0.0f, 22.0f}});
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(1, "Total");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(2, "CollIds");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(3, "Charge");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(4, "Eta");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(5, "TPCNcls");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(6, "TPCRows");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(7, "TPCpid");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(8, "DCAxyPV");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(9, "DCAzPV");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(10, "V0MassConst");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(11, "HasSV");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(12, "DcaDau");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(13, "DCADauVtx");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(14, "DauPt");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(15, "Rapidity");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(16, "Pt");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(17, "Mass");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(18, "CosPA");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(19, "CosPAXY");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(20, "Chi2geo");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(21, "TopoConstr");
      hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(22, "Chi2topo");
      hVtx3BodyCounterKFParticle->LabelsOption("v");

      registry.add("QA/Tracks/hTrackPosTPCNcls", "hTrackPosTPCNcls", HistType::kTH1F, {{152, 0, 152, "# TPC clusters"}});
      registry.add("QA/Tracks/hTrackNegTPCNcls", "hTrackNegTPCNcls", HistType::kTH1F, {{152, 0, 152, "# TPC clusters"}});
      registry.add("QA/Tracks/hTrackBachTPCNcls", "hTrackBachTPCNcls", HistType::kTH1F, {{152, 0, 152, "# TPC clusters"}});
      registry.add("QA/Tracks/hTrackPosHasTPC", "hTrackPosHasTPC", HistType::kTH1F, {{2, -0.5, 1.5, "has TPC"}});
      registry.add("QA/Tracks/hTrackNegHasTPC", "hTrackNegHasTPC", HistType::kTH1F, {{2, -0.5, 1.5, "has TPC"}});
      registry.add("QA/Tracks/hTrackBachHasTPC", "hTrackBachHasTPC", HistType::kTH1F, {{2, -0.5, 1.5, "has TPC"}});
      registry.add("QA/Tracks/hTrackBachITSClusSizes", "hTrackBachITSClusSizes", HistType::kTH1F, {{10, 0., 10., "ITS cluster sizes"}});
      registry.add("QA/Tracks/hTrackProtonTPCPID", "hTrackProtonTPCPID", HistType::kTH2F, {{100, -10.0f, 10.0f, "p/z (GeV/c)"}, {100, -10.0f, 10.0f, "TPC n#sigma"}});
      registry.add("QA/Tracks/hTrackPionTPCPID", "hTrackPionTPCPID", HistType::kTH2F, {{100, -10.0f, 10.0f, "p/z (GeV/c)"}, {100, -10.0f, 10.0f, "TPC n#sigma"}});
      registry.add("QA/Tracks/hTrackBachTPCPID", "hTrackBachTPCPID", HistType::kTH2F, {{100, -10.0f, 10.0f, "p/z (GeV/c)"}, {100, -10.0f, 10.0f, "TPC n#sigma"}});
      registry.add("QA/Tracks/hTrackProtonPt", "hTrackProtonPt", HistType::kTH1F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}});
      registry.add("QA/Tracks/hTrackPionPt", "hTrackPionPt", HistType::kTH1F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}});
      registry.add("QA/Tracks/hTrackBachPt", "hTrackBachPt", HistType::kTH1F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}});
      registry.add("QA/Event/hVtxXKF", "hVtxXKF", HistType::kTH1F, {{500, -0.1f, 0.1f, "PV X (cm)"}});
      registry.add("QA/Event/hVtxYKF", "hVtxYKF", HistType::kTH1F, {{500, -0.1f, 0.1f, "PV Y (cm)"}});
      registry.add("QA/Event/hVtxZKF", "hVtxZKF", HistType::kTH1F, {{500, -15.0f, 15.0f, "PV Z (cm)"}});
      registry.add("QA/Event/hVtxCovXXKF", "hVtxCovXXKF", HistType::kTH1F, {{200, -0.0001f, 0.0001f, "PV cov(XX) (cm^{2})"}});
      registry.add("QA/Event/hVtxCovYYKF", "hVtxCovYYKF", HistType::kTH1F, {{200, -0.0001f, 0.0001f, "PV cov(YY) (cm^{2})"}});
      registry.add("QA/Event/hVtxCovZZKF", "hVtxCovZZKF", HistType::kTH1F, {{200, -0.0001f, 0.0001f, "PV cov(ZZ) (cm^{2})"}});
      registry.add("QA/Event/hVtxCovXYKF", "hVtxCovXYKF", HistType::kTH1F, {{200, -0.0001f, 0.0001f, "PV cov(XY) (cm^{2})"}});
      registry.add("QA/Event/hVtxCovXZKF", "hVtxCovXZKF", HistType::kTH1F, {{200, -0.0001f, 0.0001f, "PV cov(XZ) (cm^{2})"}});
      registry.add("QA/Event/hVtxCovYZKF", "hVtxCovYZKF", HistType::kTH1F, {{200, -0.0001f, 0.0001f, "PV cov(YZ) (cm^{2})"}});
      registry.add("QA/Event/hVtxX", "hVtxX", HistType::kTH1F, {{500, -0.1f, 0.1f, "PV X (cm)"}});
      registry.add("QA/Event/hVtxY", "hVtxY", HistType::kTH1F, {{500, -0.1f, 0.1f, "PV Y (cm)"}});
      registry.add("QA/Event/hVtxZ", "hVtxZ", HistType::kTH1F, {{500, -15.0f, 15.0f, "PV Z (cm)"}});
      registry.add("QA/Event/hVtxCovXX", "hVtxCovXX", HistType::kTH1F, {{200, -0.0001f, 0.0001f, "PV cov(XX) (cm^{2})"}});
      registry.add("QA/Event/hVtxCovYY", "hVtxCovYY", HistType::kTH1F, {{200, -0.0001f, 0.0001f, "PV cov(YY) (cm^{2})"}});
      registry.add("QA/Event/hVtxCovZZ", "hVtxCovZZ", HistType::kTH1F, {{200, -0.0001f, 0.0001f, "PV cov(ZZ) (cm^{2})"}});
      registry.add("QA/Event/hVtxCovXY", "hVtxCovXY", HistType::kTH1F, {{200, -0.0001f, 0.0001f, "PV cov(XY) (cm^{2})"}});
      registry.add("QA/Event/hVtxCovXZ", "hVtxCovXZ", HistType::kTH1F, {{200, -0.0001f, 0.0001f, "PV cov(XZ) (cm^{2})"}});
      registry.add("QA/Event/hVtxCovYZ", "hVtxCovYZ", HistType::kTH1F, {{200, -0.0001f, 0.0001f, "PV cov(YZ) (cm^{2})"}});
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      fitter3body.setBz(d_bz);
#ifdef HomogeneousField
      KFParticle::SetField(d_bz);
#endif
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      // d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      d_bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
    // Set magnetic field value once known
    fitter3body.setBz(d_bz);
// Set magnetic field for KF vertexing
#ifdef HomogeneousField
    KFParticle::SetField(d_bz);
#endif

    if (useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized
      // (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }

    // Initial TOF PID Paras, copied from PIDTOF.h
    timestamp.value = bc.timestamp();
    ccdb->setTimestamp(timestamp.value);
    // Not later than now objects
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    // TODO: implement the automatic pass name detection from metadata
    if (passName.value == "") {
      passName.value = "unanchored"; // temporary default
      LOG(warning) << "Passed autodetect mode for pass, not implemented yet, waiting for metadata. Taking '" << passName.value << "'";
    }
    LOG(info) << "Using parameter collection, starting from pass '" << passName.value << "'";

    const std::string fname = paramFileName.value;
    if (!fname.empty()) { // Loading the parametrization from file
      LOG(info) << "Loading exp. sigma parametrization from file " << fname << ", using param: " << parametrizationPath.value;
      if (1) {
        o2::tof::ParameterCollection paramCollection;
        paramCollection.loadParamFromFile(fname, parametrizationPath.value);
        LOG(info) << "+++ Loaded parameter collection from file +++";
        if (!paramCollection.retrieveParameters(mRespParamsV2, passName.value)) {
          if (fatalOnPassNotAvailable) {
            LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          } else {
            LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          }
        } else {
          mRespParamsV2.setShiftParameters(paramCollection.getPars(passName.value));
          mRespParamsV2.printShiftParameters();
        }
      } else {
        mRespParamsV2.loadParamFromFile(fname.data(), parametrizationPath.value);
      }
    } else if (loadResponseFromCCDB) { // Loading it from CCDB
      LOG(info) << "Loading exp. sigma parametrization from CCDB, using path: " << parametrizationPath.value << " for timestamp " << timestamp.value;
      o2::tof::ParameterCollection* paramCollection = ccdb->getForTimeStamp<o2::tof::ParameterCollection>(parametrizationPath.value, timestamp.value);
      paramCollection->print();
      if (!paramCollection->retrieveParameters(mRespParamsV2, passName.value)) { // Attempt at loading the parameters with the pass defined
        if (fatalOnPassNotAvailable) {
          LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
        } else {
          LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
        }
      } else { // Pass is available, load non standard parameters
        mRespParamsV2.setShiftParameters(paramCollection->getPars(passName.value));
        mRespParamsV2.printShiftParameters();
      }
    }
    mRespParamsV2.print();
    if (timeShiftCCDBPath.value != "") {
      if (timeShiftCCDBPath.value.find(".root") != std::string::npos) {
        mRespParamsV2.setTimeShiftParameters(timeShiftCCDBPath.value, "gmean_Pos", true);
        mRespParamsV2.setTimeShiftParameters(timeShiftCCDBPath.value, "gmean_Neg", false);
      } else {
        mRespParamsV2.setTimeShiftParameters(ccdb->getForTimeStamp<TGraph>(Form("%s/pos", timeShiftCCDBPath.value.c_str()), timestamp.value), true);
        mRespParamsV2.setTimeShiftParameters(ccdb->getForTimeStamp<TGraph>(Form("%s/neg", timeShiftCCDBPath.value.c_str()), timestamp.value), false);
      }
    }

    bachelorTOFPID.SetParams(mRespParamsV2);
  }

  //------------------------------------------------------------------
  // Select decay3body candidate based on daughter track PID
  template <typename TTrack>
  bool checkPID(TTrack const& trackProton, TTrack const& trackPion, TTrack const& trackBachelor, const double& tofNSigmaBach)
  {
    if ((tofNSigmaBach < TofPidNsigmaMin || tofNSigmaBach > TofPidNsigmaMax) && trackBachelor.p() > minBachPUseTOF) {
      return false;
    }
    if (std::abs(trackProton.tpcNSigmaPr()) > TpcPidNsigmaCut) {
      return false;
    }
    if (std::abs(trackPion.tpcNSigmaPi()) > TpcPidNsigmaCut) {
      return false;
    }
    return true;
  }
  // PID check for H3L
  template <typename TTrack>
  bool checkPIDH3L(TTrack const& trackProton, TTrack const& trackPion, TTrack const& trackBachelor, const double& tofNSigmaBach)
  {
    if ((std::abs(trackBachelor.tpcNSigmaDe()) > TpcPidNsigmaCut) || !checkPID(trackProton, trackPion, trackBachelor, tofNSigmaBach)) {
      return false;
    }
    return true;
  }

  //------------------------------------------------------------------
  // function to select daughter track PID
  template <typename TTrack>
  bool selectTPCPID(TTrack const& trackProton, TTrack const& trackPion, TTrack const& trackDeuteron)
  {
    if (std::abs(trackProton.tpcNSigmaPr()) > kfparticleConfigurations.maxtpcnSigma) {
      return false;
    }
    if (std::abs(trackDeuteron.tpcNSigmaDe()) > kfparticleConfigurations.maxtpcnSigma) {
      return false;
    }
    if (kfparticleConfigurations.useTPCforPion && std::abs(trackPion.tpcNSigmaPi()) > kfparticleConfigurations.maxtpcnSigma) {
      return false;
    }
    return true;
  }

  //------------------------------------------------------------------
  // 3body candidate builder
  template <class TCollisionClass, typename TCollisionTable, typename TTrackTable>
  void fillVtxCand(TCollisionTable const& collision, TTrackTable const& t0, TTrackTable const& t1, TTrackTable const& t2, int64_t decay3bodyId, int bachelorcharge = 1)
  {

    registry.fill(HIST("hVtx3BodyCounter"), kVtxAll);

    if (t0.tpcNClsFound() < mintpcNCls || t1.tpcNClsFound() < mintpcNCls || t2.tpcNClsFound() < mintpcNCls) {
      return;
    }
    registry.fill(HIST("hVtx3BodyCounter"), kVtxTPCNcls);

    // Recalculate the TOF PID
    double tofNSigmaBach = -999;
    if (t2.has_collision() && t2.hasTOF()) {
      if (decay3bodyId == -1) {
        // for event-mixing, the collisionId of tracks not equal to global index
        tofNSigmaBach = bachelorTOFPID.GetTOFNSigma(t2, collision, collision);
      } else {
        auto originalcol = t2.template collision_as<TCollisionClass>();
        tofNSigmaBach = bachelorTOFPID.GetTOFNSigma(t2, originalcol, collision);
      }
    }

    if (enablePidCut) {
      if (t2.sign() > 0) {
        if (!checkPIDH3L(t0, t1, t2, tofNSigmaBach))
          return;
      } else {
        if (!checkPIDH3L(t1, t0, t2, tofNSigmaBach))
          return;
      }
    }

    registry.fill(HIST("hVtx3BodyCounter"), kVtxPIDCut);

    // Calculate DCA with respect to the collision associated to the V0, not individual tracks
    gpu::gpustd::array<float, 2> dcaInfo;

    auto Track0Par = getTrackPar(t0);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, Track0Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
    auto Track0dcaXY = dcaInfo[0];
    auto Track0dca = std::sqrt(Track0dcaXY * Track0dcaXY + dcaInfo[1] * dcaInfo[1]);

    auto Track1Par = getTrackPar(t1);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, Track1Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
    auto Track1dcaXY = dcaInfo[0];
    auto Track1dca = std::sqrt(Track1dcaXY * Track1dcaXY + dcaInfo[1] * dcaInfo[1]);

    auto Track2Par = getTrackPar(t2);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, Track2Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
    auto Track2dcaXY = dcaInfo[0];
    auto Track2dca = std::sqrt(Track2dcaXY * Track2dcaXY + dcaInfo[1] * dcaInfo[1]);

    auto Track0 = getTrackParCov(t0);
    auto Track1 = getTrackParCov(t1);
    auto Track2 = getTrackParCov(t2);
    int n3bodyVtx = fitter3body.process(Track0, Track1, Track2);
    if (n3bodyVtx == 0) { // discard this pair
      return;
    }
    registry.fill(HIST("hVtx3BodyCounter"), kVtxhasSV);

    std::array<float, 3> pos = {0.};
    const auto& vtxXYZ = fitter3body.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      pos[i] = vtxXYZ[i];
    }

    std::array<float, 3> p0 = {0.}, p1 = {0.}, p2{0.};
    const auto& propagatedTrack0 = fitter3body.getTrack(0);
    const auto& propagatedTrack1 = fitter3body.getTrack(1);
    const auto& propagatedTrack2 = fitter3body.getTrack(2);
    propagatedTrack0.getPxPyPzGlo(p0);
    propagatedTrack1.getPxPyPzGlo(p1);
    propagatedTrack2.getPxPyPzGlo(p2);
    for (int i = 0; i < 3; i++) {
      p2[i] *= bachelorcharge;
    }
    std::array<float, 3> p3B = {p0[0] + p1[0] + p2[0], p0[1] + p1[1] + p2[1], p0[2] + p1[2] + p2[2]};

    if (fitter3body.getChi2AtPCACandidate() > dcavtxdau) {
      return;
    }
    registry.fill(HIST("hVtx3BodyCounter"), kVtxDcaDau);

    float VtxcosPA = RecoDecay::cpa(array{collision.posX(), collision.posY(), collision.posZ()}, array{pos[0], pos[1], pos[2]}, array{p3B[0], p3B[1], p3B[2]});
    if (VtxcosPA < minCosPA3body) {
      return;
    }
    registry.fill(HIST("hVtx3BodyCounter"), kVtxCosPA);

    registry.fill(HIST("hBachelorTOFNSigmaDe"), t2.sign() * t2.p(), tofNSigmaBach);

    vtxCandidate candVtx;
    candVtx.track0Id = t0.globalIndex();
    candVtx.track1Id = t1.globalIndex();
    candVtx.track2Id = t2.globalIndex();
    candVtx.collisionId = collision.globalIndex();
    candVtx.decay3bodyId = decay3bodyId;
    candVtx.vtxPos[0] = pos[0];
    candVtx.vtxPos[1] = pos[1];
    candVtx.vtxPos[2] = pos[2];
    candVtx.track0P[0] = p0[0];
    candVtx.track0P[1] = p0[1];
    candVtx.track0P[2] = p0[2];
    candVtx.track1P[0] = p1[0];
    candVtx.track1P[1] = p1[1];
    candVtx.track1P[2] = p1[2];
    candVtx.track2P[0] = p2[0];
    candVtx.track2P[1] = p2[1];
    candVtx.track2P[2] = p2[2];
    candVtx.dcadaughters = fitter3body.getChi2AtPCACandidate();
    candVtx.daudcaxytopv[0] = Track0dcaXY;
    candVtx.daudcaxytopv[1] = Track1dcaXY;
    candVtx.daudcaxytopv[2] = Track2dcaXY;
    candVtx.daudcatopv[0] = Track0dca;
    candVtx.daudcatopv[1] = Track1dca;
    candVtx.daudcatopv[2] = Track2dca;
    candVtx.bachelortofNsigma = tofNSigmaBach;
    vtxCandidates.push_back(candVtx);
  }
  //------------------------------------------------------------------
  // fill the StoredVtx3BodyDatas table
  void fillVtx3BodyTable(vtxCandidate const& candVtx)
  {
    vtx3bodydata(
      candVtx.track0Id, candVtx.track1Id, candVtx.track2Id, candVtx.collisionId, candVtx.decay3bodyId,
      candVtx.vtxPos[0], candVtx.vtxPos[1], candVtx.vtxPos[2],
      candVtx.track0P[0], candVtx.track0P[1], candVtx.track0P[2], candVtx.track1P[0], candVtx.track1P[1], candVtx.track1P[2], candVtx.track2P[0], candVtx.track2P[1], candVtx.track2P[2],
      candVtx.dcadaughters,
      candVtx.daudcaxytopv[0], candVtx.daudcaxytopv[1], candVtx.daudcaxytopv[2],
      candVtx.daudcatopv[0], candVtx.daudcatopv[1], candVtx.daudcatopv[2],
      candVtx.bachelortofNsigma);
  }

  //------------------------------------------------------------------
  // 3body candidate builder with KFParticle
  template <class TTrackTo, class TCollisionTo, typename TCollision>
  void buildVtx3BodyDataTableKFParticle(TCollision const& collision, aod::Decay3Bodys const& decay3bodys, int bachelorcharge = 1)
  {
    LOG(debug) << "buildVtx3BodyDataTableKFParticle called.";

    // initialise KF primary vertex
    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle kfpv(kfpVertex);
    LOG(debug) << "Created KF PV.";

    // fill event QA histograms
    if (kfparticleConfigurations.doVertexQA) {
      registry.fill(HIST("QA/Event/hVtxXKF"), kfpv.GetX());
      registry.fill(HIST("QA/Event/hVtxYKF"), kfpv.GetY());
      registry.fill(HIST("QA/Event/hVtxZKF"), kfpv.GetZ());
      registry.fill(HIST("QA/Event/hVtxCovXXKF"), kfpv.GetCovariance(0));
      registry.fill(HIST("QA/Event/hVtxCovYYKF"), kfpv.GetCovariance(2));
      registry.fill(HIST("QA/Event/hVtxCovZZKF"), kfpv.GetCovariance(5));
      registry.fill(HIST("QA/Event/hVtxCovXYKF"), kfpv.GetCovariance(1));
      registry.fill(HIST("QA/Event/hVtxCovXZKF"), kfpv.GetCovariance(3));
      registry.fill(HIST("QA/Event/hVtxCovYZKF"), kfpv.GetCovariance(4));
      registry.fill(HIST("QA/Event/hVtxX"), collision.posX());
      registry.fill(HIST("QA/Event/hVtxY"), collision.posY());
      registry.fill(HIST("QA/Event/hVtxZ"), collision.posZ());
      registry.fill(HIST("QA/Event/hVtxCovXX"), collision.covXX());
      registry.fill(HIST("QA/Event/hVtxCovYY"), collision.covYY());
      registry.fill(HIST("QA/Event/hVtxCovZZ"), collision.covZZ());
      registry.fill(HIST("QA/Event/hVtxCovXY"), collision.covXY());
      registry.fill(HIST("QA/Event/hVtxCovXZ"), collision.covXZ());
      registry.fill(HIST("QA/Event/hVtxCovYZ"), collision.covYZ());
    }

    for (auto& vtx3body : decay3bodys) {
      LOG(debug) << "Entered decay3bodys loop.";

      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxAll);

      auto trackPos = vtx3body.template track0_as<TTrackTo>();
      auto trackNeg = vtx3body.template track1_as<TTrackTo>();
      auto trackBach = vtx3body.template track2_as<TTrackTo>();
      auto trackParCovPos = getTrackParCov(trackPos);
      auto trackParCovNeg = getTrackParCov(trackNeg);
      auto trackParCovBach = getTrackParCov(trackBach);
      LOG(debug) << "Got all daughter tracks.";

      bool isMatter = trackBach.sign() > 0 ? true : false;

      // ---------- fill track QA histograms ----------
      if (kfparticleConfigurations.doTrackQA) {
        registry.fill(HIST("QA/Tracks/hTrackPosTPCNcls"), trackPos.tpcNClsFound());
        registry.fill(HIST("QA/Tracks/hTrackNegTPCNcls"), trackNeg.tpcNClsFound());
        registry.fill(HIST("QA/Tracks/hTrackBachTPCNcls"), trackBach.tpcNClsFound());
        registry.fill(HIST("QA/Tracks/hTrackPosHasTPC"), trackPos.hasTPC());
        registry.fill(HIST("QA/Tracks/hTrackNegHasTPC"), trackNeg.hasTPC());
        registry.fill(HIST("QA/Tracks/hTrackBachHasTPC"), trackBach.hasTPC());
        registry.fill(HIST("QA/Tracks/hTrackBachITSClusSizes"), trackBach.itsClusterSizes());
        if (isMatter) {
          registry.fill(HIST("QA/Tracks/hTrackProtonTPCPID"), trackPos.sign() * trackPos.tpcInnerParam(), trackPos.tpcNSigmaPr());
          registry.fill(HIST("QA/Tracks/hTrackPionTPCPID"), trackNeg.sign() * trackNeg.tpcInnerParam(), trackNeg.tpcNSigmaPi());
          registry.fill(HIST("QA/Tracks/hTrackProtonPt"), trackPos.pt());
          registry.fill(HIST("QA/Tracks/hTrackPionPt"), trackNeg.pt());
        } else {
          registry.fill(HIST("QA/Tracks/hTrackProtonTPCPID"), trackNeg.sign() * trackNeg.tpcInnerParam(), trackNeg.tpcNSigmaPr());
          registry.fill(HIST("QA/Tracks/hTrackPionTPCPID"), trackPos.sign() * trackPos.tpcInnerParam(), trackPos.tpcNSigmaPi());
          registry.fill(HIST("QA/Tracks/hTrackProtonPt"), trackNeg.pt());
          registry.fill(HIST("QA/Tracks/hTrackPionPt"), trackPos.pt());
        }
        registry.fill(HIST("QA/Tracks/hTrackBachTPCPID"), trackBach.sign() * trackBach.tpcInnerParam(), trackBach.tpcNSigmaDe());
        registry.fill(HIST("QA/Tracks/hTrackBachPt"), trackBach.pt());
      }

      // -------- STEP 1: track selection --------
      // collision ID --> not correct? tracks can have different collisions, but belong to one 3prong vertex!
      // if (trackPos.collisionId() != trackNeg.collisionId() || trackPos.collisionId() != trackBach.collisionId() || trackNeg.collisionId() != trackBach.collisionId()) {
      //   continue;
      // }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxCollIds);
      // track IDs --> already checked in SVertexer!

      // track signs (pos, neg, bach) --> sanity check, should already be in SVertexer
      if (trackPos.sign() != +1 || trackNeg.sign() != -1) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxCharge);

      // track eta
      if (abs(trackPos.eta()) > kfparticleConfigurations.maxEta || abs(trackNeg.eta()) > kfparticleConfigurations.maxEta || abs(trackBach.eta()) > kfparticleConfigurations.maxEta) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxEta);

      // number of TPC clusters
      if (trackBach.tpcNClsFound() <= kfparticleConfigurations.mintpcNClsBach) {
        continue;
      }
      if (isMatter && ((kfparticleConfigurations.useTPCforPion && trackNeg.tpcNClsFound() <= kfparticleConfigurations.mintpcNClsPion) || trackPos.tpcNClsFound() <= kfparticleConfigurations.mintpcNClsProton)) {
        continue;
      } else if (!isMatter && ((kfparticleConfigurations.useTPCforPion && trackPos.tpcNClsFound() <= kfparticleConfigurations.mintpcNClsPion) || trackNeg.tpcNClsFound() <= kfparticleConfigurations.mintpcNClsProton)) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxTPCNcls);

      // number of TPC crossed rows
      if (trackBach.tpcNClsCrossedRows() <= kfparticleConfigurations.mintpcCrossedRows) {
        continue;
      }
      if (isMatter && ((kfparticleConfigurations.useTPCforPion && trackNeg.tpcNClsCrossedRows() <= kfparticleConfigurations.mintpcCrossedRowsPion) || trackPos.tpcNClsCrossedRows() <= kfparticleConfigurations.mintpcCrossedRows)) {
        continue;
      } else if (!isMatter && ((kfparticleConfigurations.useTPCforPion && trackPos.tpcNClsCrossedRows() <= kfparticleConfigurations.mintpcCrossedRowsPion) || trackNeg.tpcNClsCrossedRows() <= kfparticleConfigurations.mintpcCrossedRows)) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxTPCRows);

      // TPC PID
      float tpcNsigmaProton;
      float tpcNsigmaPion;
      float dEdxProton;
      float dEdxPion;
      float tpcNsigmaDeuteron = trackBach.tpcNSigmaDe();
      float tpcNsigmaPionBach = trackBach.tpcNSigmaPi();
      float dEdxDeuteron = trackBach.tpcSignal();
      if (isMatter) { // hypertriton (proton, pi-, deuteron)
        tpcNsigmaProton = trackPos.tpcNSigmaPr();
        tpcNsigmaPion = trackNeg.tpcNSigmaPi();
        dEdxProton = trackPos.tpcSignal();
        dEdxPion = trackNeg.tpcSignal();
        if (!selectTPCPID(trackPos, trackNeg, trackBach)) {
          continue;
        }
      } else if (!isMatter) { // anti-hypertriton (anti-proton, pi+, deuteron)
        tpcNsigmaProton = trackNeg.tpcNSigmaPr();
        tpcNsigmaPion = trackPos.tpcNSigmaPi();
        dEdxProton = trackNeg.tpcSignal();
        dEdxPion = trackPos.tpcSignal();
        if (!selectTPCPID(trackNeg, trackPos, trackBach)) {
          continue;
        }
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxTPCPID);
      LOG(debug) << "Basic track selections done.";

      // TOF PID of deuteron (set motherhyp correctly)
      double tofNSigmaDeuteron = -999;
      if (trackBach.has_collision() && trackBach.hasTOF()) {
        auto originalcol = trackBach.template collision_as<TCollisionTo>();
        tofNSigmaDeuteron = bachelorTOFPID.GetTOFNSigma(trackBach, originalcol, collision);
      }

      // Average ITS cluster size of deuteron track
      double averageClusterSizeDeuteron(0);
      int nCls(0);
      for (int i = 0; i < 7; i++) {
        int clusterSize = trackBach.itsClsSizeInLayer(i);
        averageClusterSizeDeuteron += static_cast<double>(clusterSize);
        if (clusterSize > 0)
          nCls++;
      }
      averageClusterSizeDeuteron = averageClusterSizeDeuteron / static_cast<double>(nCls);

      // track DCAxy and DCAz to PV associated with decay3body
      o2::dataformats::VertexBase mPV;
      o2::dataformats::DCA mDcaInfoCovPos;
      o2::dataformats::DCA mDcaInfoCovNeg;
      o2::dataformats::DCA mDcaInfoCovBach;
      auto trackParCovPVPos = trackParCovPos;
      auto trackParCovPVNeg = trackParCovNeg;
      auto trackParCovPVBach = trackParCovBach;
      mPV.setPos({collision.posX(), collision.posY(), collision.posZ()});
      mPV.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
      o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, trackParCovPVPos, 2.f, matCorr, &mDcaInfoCovPos);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, trackParCovPVNeg, 2.f, matCorr, &mDcaInfoCovNeg);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, trackParCovPVBach, 2.f, matCorr, &mDcaInfoCovBach);
      auto TrackPosDcaXY = mDcaInfoCovPos.getY();
      auto TrackNegDcaXY = mDcaInfoCovNeg.getY();
      auto TrackBachDcaXY = mDcaInfoCovBach.getY();
      auto TrackPosDcaZ = mDcaInfoCovPos.getZ();
      auto TrackNegDcaZ = mDcaInfoCovNeg.getZ();
      auto TrackBachDcaZ = mDcaInfoCovBach.getZ();
      if (isMatter && (fabs(TrackNegDcaXY) <= kfparticleConfigurations.mindcaXYPionPV || fabs(TrackPosDcaXY) <= kfparticleConfigurations.mindcaXYProtonPV)) {
        continue;
      } else if (!isMatter && (fabs(TrackPosDcaXY) <= kfparticleConfigurations.mindcaXYPionPV || fabs(TrackNegDcaXY) <= kfparticleConfigurations.mindcaXYProtonPV)) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxDCAxyPV);
      if (isMatter && (fabs(TrackNegDcaZ) <= kfparticleConfigurations.mindcaZPionPV || fabs(TrackPosDcaZ) <= kfparticleConfigurations.mindcaZProtonPV)) {
        continue;
      } else if (!isMatter && (fabs(TrackPosDcaZ) <= kfparticleConfigurations.mindcaZPionPV || fabs(TrackNegDcaZ) <= kfparticleConfigurations.mindcaZProtonPV)) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxDCAzPV);
      // calculate 3D track DCA
      auto TrackPosDca = std::sqrt(TrackPosDcaXY * TrackPosDcaXY + TrackPosDcaZ * TrackPosDcaZ);
      auto TrackNegDca = std::sqrt(TrackNegDcaXY * TrackNegDcaXY + TrackNegDcaZ * TrackNegDcaZ);
      auto TrackBachDca = std::sqrt(TrackBachDcaXY * TrackBachDcaXY + TrackBachDcaZ * TrackBachDcaZ);

      // daughter track momentum at inner wall of TPC
      float tpcInnerParamProton;
      float tpcInnerParamPion;
      float tpcInnerParamDeuteron = trackBach.tpcInnerParam();
      if (isMatter) { // hypertriton (proton, pi-, deuteron)
        tpcInnerParamProton = trackPos.tpcInnerParam();
        tpcInnerParamPion = trackNeg.tpcInnerParam();
      } else if (!isMatter) { // anti-hypertriton (anti-proton, pi+, deuteron)
        tpcInnerParamProton = trackNeg.tpcInnerParam();
        tpcInnerParamPion = trackPos.tpcInnerParam();
      }

      // -------- STEP 2: fit vertex with proton and pion --------
      // Fit vertex with DCA fitter to find minimization point --> uses material corrections implicitly
      if (kfparticleConfigurations.doDCAFitterPreMinimum) {
        try {
          fitter3body.process(trackParCovPos, trackParCovNeg, trackParCovBach);
        } catch (std::runtime_error& e) {
          LOG(error) << "Exception caught in DCA fitter process call: Not able to fit decay3body vertex!";
          continue;
        }
        // re-acquire tracks at vertex position from DCA fitter
        trackParCovPos = fitter3body.getTrack(0);
        trackParCovNeg = fitter3body.getTrack(1);
        trackParCovBach = fitter3body.getTrack(2);

        LOG(debug) << "Minimum found with DCA fitter for decay3body.";
      }

      // create KFParticle objects from tracks
      KFParticle kfpProton, kfpPion;
      if (isMatter) {
        kfpProton = createKFParticleFromTrackParCov(trackParCovPos, trackPos.sign(), constants::physics::MassProton);
        kfpPion = createKFParticleFromTrackParCov(trackParCovNeg, trackNeg.sign(), constants::physics::MassPionCharged);
      } else if (!isMatter) {
        kfpProton = createKFParticleFromTrackParCov(trackParCovNeg, trackNeg.sign(), constants::physics::MassProton);
        kfpPion = createKFParticleFromTrackParCov(trackParCovPos, trackPos.sign(), constants::physics::MassPionCharged);
      }
      LOG(debug) << "KFParticle objects created from daughter tracks.";

      // Construct V0 as intermediate step
      KFParticle KFV0;
      int nDaughtersV0 = 2;
      const KFParticle* DaughtersV0[2] = {&kfpProton, &kfpPion};
      KFV0.SetConstructMethod(2);
      try {
        KFV0.Construct(DaughtersV0, nDaughtersV0);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to create V0 vertex from daughter tracks." << e.what();
        continue;
      }
      KFV0.TransportToDecayVertex();
      LOG(debug) << "V0 constructed.";

      // check V0 mass and set mass constraint
      float massV0, sigmaMassV0;
      KFV0.GetMass(massV0, sigmaMassV0);
      KFParticle KFV0Mass = KFV0;
      KFV0Mass.SetNonlinearMassConstraint(o2::constants::physics::MassLambda);
      float chi2massV0 = KFV0Mass.GetChi2() / KFV0Mass.GetNDF();
      if (kfparticleConfigurations.useLambdaMassConstraint) {
        LOG(debug) << "V0 mass constraint applied.";
        KFV0 = KFV0Mass;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxV0MassConst);

      // -------- STEP 3: fit three body vertex --------
      // Create KFParticle object from deuteron track
      KFParticle kfpDeuteron;
      kfpDeuteron = createKFParticleFromTrackParCov(trackParCovBach, trackBach.sign() * bachelorcharge, constants::physics::MassDeuteron);
      LOG(debug) << "KFParticle created from deuteron track.";
      // Construct 3body vertex
      int nDaughters3body = 3;
      const KFParticle* Daughters3body[3] = {&kfpProton, &kfpPion, &kfpDeuteron};
      KFParticle KFHt;
      KFHt.SetConstructMethod(2);
      try {
        KFHt.Construct(Daughters3body, nDaughters3body);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to create Hyper triton 3-body vertex." << e.what();
        continue;
      }
      // transport all daughter tracks to hypertriton vertex
      float HtVtx[3] = {0.};
      HtVtx[0] = KFHt.GetX();
      HtVtx[1] = KFHt.GetY();
      HtVtx[2] = KFHt.GetZ();
      kfpProton.TransportToPoint(HtVtx);
      kfpPion.TransportToPoint(HtVtx);
      kfpDeuteron.TransportToPoint(HtVtx);
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxhasSV);
      LOG(debug) << "Hypertriton vertex constructed.";

      // -------- STEP 4: selections after geometrical vertex fit --------
      // daughter DCAs with KF
      if ((kfpProton.GetDistanceFromParticle(kfpPion) >= kfparticleConfigurations.maxDcaProPi) || (kfpProton.GetDistanceFromParticle(kfpDeuteron) >= kfparticleConfigurations.maxDcaProDeu) || (kfpPion.GetDistanceFromParticle(kfpDeuteron) >= kfparticleConfigurations.maxDcaPiDe)) {
        continue;
      }
      float DCAvtxDaughters3D = kfpProton.GetDistanceFromParticle(kfpPion) + kfpProton.GetDistanceFromParticle(kfpDeuteron) + kfpPion.GetDistanceFromParticle(kfpDeuteron);
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxDcaDau);
      LOG(debug) << "DCA selection after vertex fit applied.";

      // daughter DCAs to vertex
      if (kfpProton.GetDistanceFromVertexXY(KFHt) >= kfparticleConfigurations.maxDcaXYSVDau || kfpPion.GetDistanceFromVertexXY(KFHt) >= kfparticleConfigurations.maxDcaXYSVDau || kfpDeuteron.GetDistanceFromVertexXY(KFHt) >= kfparticleConfigurations.maxDcaXYSVDau) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxDcaDauVtx);
      LOG(debug) << "DCA to vertex selection after vertex fit applied.";

      // daughter pT
      if (kfpProton.GetPt() < kfparticleConfigurations.minPtProton || kfpProton.GetPt() > kfparticleConfigurations.maxPtProton || kfpPion.GetPt() < kfparticleConfigurations.minPtPion || kfpPion.GetPt() > kfparticleConfigurations.maxPtPion || kfpDeuteron.GetPt() < kfparticleConfigurations.minPtDeuteron || kfpDeuteron.GetPt() > kfparticleConfigurations.maxPtDeuteron) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxDauPt);
      LOG(debug) << "Daughter pT selection applied.";

      // -------- STEP 5: candidate selection after geometrical vertex fit --------
      // Rapidity
      float rapHt = RecoDecay::y(std::array{KFHt.GetPx(), KFHt.GetPy(), KFHt.GetPz()}, o2::constants::physics::MassHyperTriton);
      if (std::abs(rapHt) > kfparticleConfigurations.maxRapidityHt) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxRap);

      // Pt selection
      if (KFHt.GetPt() <= kfparticleConfigurations.minPtHt || KFHt.GetPt() >= kfparticleConfigurations.maxPtHt) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxPt);

      // Mass window
      float massHt, sigmaMassHt;
      KFHt.GetMass(massHt, sigmaMassHt);
      if (massHt <= kfparticleConfigurations.minMassHt || massHt >= kfparticleConfigurations.maxMassHt) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxMass);

      // cos(PA) to PV
      if (std::abs(cpaFromKF(KFHt, kfpv)) <= kfparticleConfigurations.minCosPA) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxCosPA);

      // cos(PA) xy to PV
      if (std::abs(cpaXYFromKF(KFHt, kfpv)) <= kfparticleConfigurations.minCosPAxy) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxCosPAXY);

      // chi2 geometrical
      float chi2geoNDF = KFHt.GetChi2() / KFHt.GetNDF();
      if (chi2geoNDF >= kfparticleConfigurations.maxChi2geo) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxChi2geo);
      LOG(debug) << "Basic selections after vertex fit done.";

      // ctau before topo constraint
      if (KFHt.GetLifeTime() > kfparticleConfigurations.maxctauHt) {
        return;
      }

      // -------- STEP 6: topological constraint --------
      /// Set vertex constraint and topological selection
      KFParticle KFHtPV = KFHt;
      try {
        KFHtPV.SetProductionVertex(kfpv);
      } catch (std::runtime_error& e) {
        LOG(error) << "Exception caught KFParticle process call: Topological constraint failed";
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxTopoConstr); // to check if topo constraint fails
      // get topological chi2
      float chi2topoNDF = KFHtPV.GetChi2() / KFHtPV.GetNDF();
      KFHtPV.TransportToDecayVertex();
      if (kfparticleConfigurations.applyTopoSel && chi2topoNDF >= kfparticleConfigurations.maxChi2topo) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxChi2topo);

      //------------------------------------------------------------------
      // table filling
      kfvtx3bodydata(
        collision.globalIndex(), trackPos.globalIndex(), trackNeg.globalIndex(), trackBach.globalIndex(), vtx3body.globalIndex(),
        // hypertriton
        massHt,
        KFHt.GetX(), KFHt.GetY(), KFHt.GetZ(),
        KFHt.GetErrX(), KFHt.GetErrY(), KFHt.GetErrZ(),
        KFHt.GetPx(), KFHt.GetPy(), KFHt.GetPz(), KFHt.GetPt(),
        KFHt.GetErrPx(), KFHt.GetErrPy(), KFHt.GetErrPz(), KFHt.GetErrPt(),
        KFHt.GetQ(),
        KFHt.GetDistanceFromVertex(kfpv), KFHt.GetDistanceFromVertexXY(kfpv),
        cpaFromKF(KFHt, kfpv), // before topo constraint
        cpaXYFromKF(KFHt, kfpv),
        cpaFromKF(KFHtPV, kfpv), // after topo constraint
        cpaXYFromKF(KFHtPV, kfpv),
        KFHtPV.GetDecayLength(), KFHtPV.GetDecayLengthXY(),   // decay length defined after topological constraint
        KFHtPV.GetDecayLength() / KFHtPV.GetErrDecayLength(), // ldl
        chi2geoNDF, chi2topoNDF,
        KFHtPV.GetLifeTime(),
        // V0
        massV0, chi2massV0,
        cpaFromKF(KFV0, kfpv),
        // daughter momenta at vertex
        kfpProton.GetPx(), kfpProton.GetPy(), kfpProton.GetPz(),
        kfpPion.GetPx(), kfpPion.GetPy(), kfpPion.GetPz(),
        kfpDeuteron.GetPx(), kfpDeuteron.GetPy(), kfpDeuteron.GetPz(),
        // daughter momenta at inner wall of TPC
        tpcInnerParamProton, tpcInnerParamPion, tpcInnerParamDeuteron,
        // daughter DCAs KF
        kfpProton.GetDistanceFromVertex(kfpv),
        kfpPion.GetDistanceFromVertex(kfpv),
        kfpDeuteron.GetDistanceFromVertex(kfpv),
        kfpProton.GetDistanceFromVertexXY(kfpv),
        kfpPion.GetDistanceFromVertexXY(kfpv),
        kfpDeuteron.GetDistanceFromVertexXY(kfpv),
        kfpProton.GetDistanceFromVertexXY(KFHt),
        kfpPion.GetDistanceFromVertexXY(KFHt),
        kfpDeuteron.GetDistanceFromVertexXY(KFHt),
        kfpProton.GetDistanceFromParticle(kfpPion),
        kfpProton.GetDistanceFromParticle(kfpDeuteron),
        kfpPion.GetDistanceFromParticle(kfpDeuteron),
        DCAvtxDaughters3D,
        // daughter DCAs to PV in XY propagated with material
        TrackPosDcaXY, TrackNegDcaXY, TrackBachDcaXY,
        // daughter DCAs to PV in 3D propagated with material
        TrackPosDca, TrackNegDca, TrackBachDca,
        // daughter signs
        kfpProton.GetQ(),
        kfpPion.GetQ(),
        trackBach.sign(),
        // daughter PID
        tpcNsigmaProton, tpcNsigmaPion, tpcNsigmaDeuteron, tpcNsigmaPionBach,
        dEdxProton, dEdxPion, dEdxDeuteron,
        tofNSigmaDeuteron,
        averageClusterSizeDeuteron,
        trackBach.pidForTracking());

      if (kfparticleConfigurations.fillCandidateLiteTable) {
        kfvtx3bodydatalite(
          // hypertriton
          massHt,
          KFHt.GetX(), KFHt.GetY(), KFHt.GetZ(),
          KFHt.GetPx(), KFHt.GetPy(), KFHt.GetPz(), KFHt.GetPt(),
          KFHt.GetQ(),
          KFHt.GetDistanceFromVertex(kfpv), KFHt.GetDistanceFromVertexXY(kfpv),
          cpaFromKF(KFHt, kfpv), // before topo constraint
          cpaXYFromKF(KFHt, kfpv),
          KFHtPV.GetDecayLength(), KFHtPV.GetDecayLengthXY(),   // decay length defined after topological constraint
          KFHtPV.GetDecayLength() / KFHtPV.GetErrDecayLength(), // ldl
          chi2geoNDF, chi2topoNDF,
          KFHtPV.GetLifeTime(),
          // V0
          massV0, chi2massV0,
          cpaFromKF(KFV0, kfpv),
          // daughter momenta at vertex
          kfpProton.GetPx(), kfpProton.GetPy(), kfpProton.GetPz(),
          kfpPion.GetPx(), kfpPion.GetPy(), kfpPion.GetPz(),
          kfpDeuteron.GetPx(), kfpDeuteron.GetPy(), kfpDeuteron.GetPz(),
          // daughter momenta at inner wall of TPC
          tpcInnerParamProton, tpcInnerParamPion, tpcInnerParamDeuteron,
          // daughter DCAs KF
          kfpProton.GetDistanceFromVertex(kfpv),
          kfpPion.GetDistanceFromVertex(kfpv),
          kfpDeuteron.GetDistanceFromVertex(kfpv),
          kfpProton.GetDistanceFromVertexXY(kfpv),
          kfpPion.GetDistanceFromVertexXY(kfpv),
          kfpDeuteron.GetDistanceFromVertexXY(kfpv),
          kfpProton.GetDistanceFromVertexXY(KFHt),
          kfpPion.GetDistanceFromVertexXY(KFHt),
          kfpDeuteron.GetDistanceFromVertexXY(KFHt),
          kfpProton.GetDistanceFromParticle(kfpPion),
          kfpProton.GetDistanceFromParticle(kfpDeuteron),
          kfpPion.GetDistanceFromParticle(kfpDeuteron),
          DCAvtxDaughters3D,
          // daughter signs
          kfpProton.GetQ(),
          kfpPion.GetQ(),
          trackBach.sign(),
          // daughter PID
          tpcNsigmaProton, tpcNsigmaPion, tpcNsigmaDeuteron, tpcNsigmaPionBach,
          dEdxProton, dEdxPion, dEdxDeuteron,
          tofNSigmaDeuteron,
          averageClusterSizeDeuteron,
          trackBach.pidForTracking());
      }
      LOG(debug) << "Table filled.";

      // fill event counter hist (has selected candidate)
      registry.fill(HIST("hEventCounterKFParticle"), 3.5);
    }
  }

  //------------------------------------------------------------------
  void processRun3(ColwithEvTimes const& collisions, TrackExtPIDIUwithEvTimes const& /*tracksIU*/, aod::Decay3Bodys const& decay3bodys, aod::BCsWithTimestamps const&)
  {
    vtxCandidates.clear();

    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      registry.fill(HIST("hEventCounter"), 0.5);

      const auto& d3bodysInCollision = decay3bodys.sliceBy(perCollision, collision.globalIndex());
      for (auto& d3body : d3bodysInCollision) {
        auto t0 = d3body.template track0_as<TrackExtPIDIUwithEvTimes>();
        auto t1 = d3body.template track1_as<TrackExtPIDIUwithEvTimes>();
        auto t2 = d3body.template track2_as<TrackExtPIDIUwithEvTimes>();
        fillVtxCand<ColwithEvTimes>(collision, t0, t1, t2, d3body.globalIndex(), bachelorcharge);
      }
    }

    for (auto& candVtx : vtxCandidates) {
      fillVtx3BodyTable(candVtx);
    }
  }
  PROCESS_SWITCH(decay3bodyBuilder, processRun3, "Produce DCA fitter decay3body tables", true);

  //------------------------------------------------------------------
  // Event-mixing background
  void processRun3EM(FullCols const& collisions, TrackExtPIDIUwithEvTimes const& tracksIU, aod::BCsWithTimestamps const&)
  {

    vtxCandidates.clear();

    auto tracksTuple = std::make_tuple(tracksIU);
    BinningType binningEvent{{axisPosZ, axisCentrality}, true};
    SameKindPair<FullCols, TrackExtPIDIUwithEvTimes, BinningType> pair{binningEvent, EMTrackSel.nUseMixedEvent, -1, collisions, tracksTuple, &cache};

    Partition<TrackExtPIDIUwithEvTimes> candProtons = aod::track::signed1Pt > 0.f && nabs(aod::track::eta) <= EMTrackSel.etacut && aod::track::pt >= EMTrackSel.minProtonPt && aod::track::pt <= EMTrackSel.maxProtonPt && aod::track::tpcNClsFindable >= (uint8_t)EMTrackSel.mintpcNClsproton && nabs(aod::pidtpc::tpcNSigmaPr) <= EMTrackSel.emTpcPidNsigmaCut;
    Partition<TrackExtPIDIUwithEvTimes> candAntiProtons = aod::track::signed1Pt < 0.f && nabs(aod::track::eta) <= EMTrackSel.etacut && aod::track::pt >= EMTrackSel.minProtonPt && aod::track::pt <= EMTrackSel.maxProtonPt && aod::track::tpcNClsFindable >= (uint8_t)EMTrackSel.mintpcNClsproton && nabs(aod::pidtpc::tpcNSigmaPr) <= EMTrackSel.emTpcPidNsigmaCut;
    Partition<TrackExtPIDIUwithEvTimes> candPionPlus = aod::track::signed1Pt > 0.f && nabs(aod::track::eta) <= EMTrackSel.etacut && aod::track::pt >= EMTrackSel.minPionPt && aod::track::pt <= EMTrackSel.maxPionPt && aod::track::tpcNClsFindable >= (uint8_t)EMTrackSel.mintpcNClspion && nabs(aod::pidtpc::tpcNSigmaPi) <= EMTrackSel.emTpcPidNsigmaCut;
    Partition<TrackExtPIDIUwithEvTimes> candPionMinus = aod::track::signed1Pt < 0.f && nabs(aod::track::eta) <= EMTrackSel.etacut && aod::track::pt >= EMTrackSel.minPionPt && aod::track::pt <= EMTrackSel.maxPionPt && aod::track::tpcNClsFindable >= (uint8_t)EMTrackSel.mintpcNClspion && nabs(aod::pidtpc::tpcNSigmaPi) <= EMTrackSel.emTpcPidNsigmaCut;
    Partition<TrackExtPIDIUwithEvTimes> candBachelors = aod::track::signed1Pt > 0.f && nabs(aod::track::eta) <= EMTrackSel.etacut && aod::track::pt >= EMTrackSel.minDeuteronPt && aod::track::pt <= EMTrackSel.maxDeuteronPt && aod::track::tpcNClsFindable >= (uint8_t)EMTrackSel.mintpcNClsbachelor && nabs(aod::pidtpc::tpcNSigmaDe) <= EMTrackSel.emTpcPidNsigmaCut;
    Partition<TrackExtPIDIUwithEvTimes> candAntiBachelors = aod::track::signed1Pt < 0.f && nabs(aod::track::eta) <= EMTrackSel.etacut && aod::track::pt >= EMTrackSel.minDeuteronPt && aod::track::pt <= EMTrackSel.maxDeuteronPt && aod::track::tpcNClsFindable >= (uint8_t)EMTrackSel.mintpcNClsbachelor && nabs(aod::pidtpc::tpcNSigmaDe) <= EMTrackSel.emTpcPidNsigmaCut;
    candProtons.bindTable(tracksIU);
    candPionPlus.bindTable(tracksIU);
    candAntiProtons.bindTable(tracksIU);
    candPionMinus.bindTable(tracksIU);
    candBachelors.bindTable(tracksIU);
    candAntiBachelors.bindTable(tracksIU);

    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      if (EMTrackSel.em_event_sel8_selection && (!c1.sel8() || !c2.sel8())) {
        continue;
      }
      auto bc = c1.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto protons = candProtons->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
      auto pionsplus = candPionPlus->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
      auto antiprotons = candAntiProtons->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
      auto pionsminus = candPionMinus->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
      auto bachelors = candBachelors->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);
      auto antibachelors = candAntiBachelors->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);

      for (auto const& [tpos, tneg, tbach] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(protons, pionsminus, bachelors))) {
        fillVtxCand<FullCols>(c1, tpos, tneg, tbach, -1, bachelorcharge);
      }
      for (auto const& [tpos, tneg, tbach] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(pionsplus, antiprotons, antibachelors))) {
        fillVtxCand<FullCols>(c1, tpos, tneg, tbach, -1, bachelorcharge);
      }
    }

    // Aviod break of preslice in following workflow
    std::sort(vtxCandidates.begin(), vtxCandidates.end(), [](const vtxCandidate a, const vtxCandidate b) {
      return a.collisionId < b.collisionId;
    });

    for (auto& candVtx : vtxCandidates) {
      fillVtx3BodyTable(candVtx);
    }
  }
  PROCESS_SWITCH(decay3bodyBuilder, processRun3EM, "Produce event-mix background", false);

  //------------------------------------------------------------------
  // Event-mixing background + like-sign (to aviod deuteron with wrong collisionId)
  void processRun3EMLikeSign(FullCols const& collisions, TrackExtPIDIUwithEvTimes const& tracksIU, aod::BCsWithTimestamps const&)
  {

    vtxCandidates.clear();

    auto tracksTuple = std::make_tuple(tracksIU);
    BinningType binningEvent{{axisPosZ, axisCentrality}, true};
    SameKindPair<FullCols, TrackExtPIDIUwithEvTimes, BinningType> pair{binningEvent, 5, -1, collisions, tracksTuple, &cache};

    Partition<TrackExtPIDIUwithEvTimes> candProtons = aod::track::signed1Pt > 0.f && nabs(aod::track::eta) <= EMTrackSel.etacut && aod::track::pt >= EMTrackSel.minProtonPt && aod::track::pt <= EMTrackSel.maxProtonPt && aod::track::tpcNClsFindable >= (uint8_t)EMTrackSel.mintpcNClsproton && nabs(aod::pidtpc::tpcNSigmaPr) <= EMTrackSel.emTpcPidNsigmaCut;
    Partition<TrackExtPIDIUwithEvTimes> candPionPlus = aod::track::signed1Pt > 0.f && nabs(aod::track::eta) <= EMTrackSel.etacut && aod::track::pt >= EMTrackSel.minPionPt && aod::track::pt <= EMTrackSel.maxPionPt && aod::track::tpcNClsFindable >= (uint8_t)EMTrackSel.mintpcNClspion && nabs(aod::pidtpc::tpcNSigmaPi) <= EMTrackSel.emTpcPidNsigmaCut;
    Partition<TrackExtPIDIUwithEvTimes> candAntiProtons = aod::track::signed1Pt < 0.f && nabs(aod::track::eta) <= EMTrackSel.etacut && aod::track::pt >= EMTrackSel.minProtonPt && aod::track::pt <= EMTrackSel.maxProtonPt && aod::track::tpcNClsFindable >= (uint8_t)EMTrackSel.mintpcNClsproton && nabs(aod::pidtpc::tpcNSigmaPr) <= EMTrackSel.emTpcPidNsigmaCut;
    Partition<TrackExtPIDIUwithEvTimes> candPionMinus = aod::track::signed1Pt < 0.f && nabs(aod::track::eta) <= EMTrackSel.etacut && aod::track::pt >= EMTrackSel.minPionPt && aod::track::pt <= EMTrackSel.maxPionPt && aod::track::tpcNClsFindable >= (uint8_t)EMTrackSel.mintpcNClspion && nabs(aod::pidtpc::tpcNSigmaPi) <= EMTrackSel.emTpcPidNsigmaCut;
    Partition<TrackExtPIDIUwithEvTimes> candBachelors = aod::track::signed1Pt > 0.f && nabs(aod::track::eta) <= EMTrackSel.etacut && aod::track::pt >= EMTrackSel.minDeuteronPt && aod::track::pt <= EMTrackSel.maxDeuteronPt && aod::track::tpcNClsFindable >= (uint8_t)EMTrackSel.mintpcNClsbachelor && nabs(aod::pidtpc::tpcNSigmaDe) <= EMTrackSel.emTpcPidNsigmaCut;
    Partition<TrackExtPIDIUwithEvTimes> candAntiBachelors = aod::track::signed1Pt < 0.f && nabs(aod::track::eta) <= EMTrackSel.etacut && aod::track::pt >= EMTrackSel.minDeuteronPt && aod::track::pt <= EMTrackSel.maxDeuteronPt && aod::track::tpcNClsFindable >= (uint8_t)EMTrackSel.mintpcNClsbachelor && nabs(aod::pidtpc::tpcNSigmaDe) <= EMTrackSel.emTpcPidNsigmaCut;
    candProtons.bindTable(tracksIU);
    candPionPlus.bindTable(tracksIU);
    candAntiProtons.bindTable(tracksIU);
    candPionMinus.bindTable(tracksIU);
    candBachelors.bindTable(tracksIU);
    candAntiBachelors.bindTable(tracksIU);

    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      if (EMTrackSel.em_event_sel8_selection && (!c1.sel8() || !c2.sel8())) {
        continue;
      }
      auto bc = c1.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto protons = candProtons->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
      auto pionsplus = candPionPlus->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
      auto antiprotons = candAntiProtons->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
      auto pionsminus = candPionMinus->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
      auto bachelors = candBachelors->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);
      auto antibachelors = candAntiBachelors->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);

      for (auto const& [tpos, tneg, tbach] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(protons, pionsminus, antibachelors))) {
        fillVtxCand<FullCols>(c1, tpos, tneg, tbach, -1, bachelorcharge);
      }
      for (auto const& [tpos, tneg, tbach] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(pionsplus, antiprotons, bachelors))) {
        fillVtxCand<FullCols>(c1, tpos, tneg, tbach, -1, bachelorcharge);
      }
    }

    // Aviod break of preslice in following workflow
    std::sort(vtxCandidates.begin(), vtxCandidates.end(), [](const vtxCandidate a, const vtxCandidate b) {
      return a.collisionId < b.collisionId;
    });

    for (auto& candVtx : vtxCandidates) {
      fillVtx3BodyTable(candVtx);
    }
  }
  PROCESS_SWITCH(decay3bodyBuilder, processRun3EMLikeSign, "Produce event-mix background with like-sign method", false);
  //------------------------------------------------------------------

  void processRun3withKFParticle(ColwithEvTimes const& collisions, TrackExtPIDIUwithEvTimes const&, aod::Decay3Bodys const& decay3bodys, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      // event selection
      registry.fill(HIST("hEventCounterKFParticle"), 0.5);
      if (kfparticleConfigurations.doSel8selection && !collision.sel8()) {
        continue;
      }
      registry.fill(HIST("hEventCounterKFParticle"), 1.5);
      if (kfparticleConfigurations.doPosZselection && abs(collision.posZ()) > 10.f) {
        continue;
      }
      registry.fill(HIST("hEventCounterKFParticle"), 2.5);

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      LOG(debug) << "CCDB initialised.";

      // slice Decay3Body table by collision
      const uint64_t collIdx = collision.globalIndex();
      // LOG(debug) << "Collision index: " << collIdx;
      auto Decay3BodyTable_thisCollision = decay3bodys.sliceBy(perCollision, collIdx);
      // LOG(debug) << "Decay3Body tables sliced per collision. Calling buildVtx3BodyDataTableKFParticle function...";
      buildVtx3BodyDataTableKFParticle<TrackExtPIDIUwithEvTimes, ColwithEvTimes>(collision, Decay3BodyTable_thisCollision, bachelorcharge);
      LOG(debug) << "End of processKFParticle.";
    }
  }
  PROCESS_SWITCH(decay3bodyBuilder, processRun3withKFParticle, "Produce KFParticle decay3body tables", false);
};

// build link from decay3body -> vtx3body
struct decay3bodyDataLinkBuilder {
  Produces<aod::Decay3BodyDataLink> vtxdataLink;

  void init(InitContext const&) {}

  void process(aod::Decay3Bodys const& decay3bodytable, aod::Vtx3BodyDatas const& vtxdatatable)
  {
    std::vector<int> lIndices;
    lIndices.reserve(decay3bodytable.size());
    for (int ii = 0; ii < decay3bodytable.size(); ii++)
      lIndices[ii] = -1;
    for (auto& vtxdata : vtxdatatable) {
      if (vtxdata.decay3bodyId() != -1) {
        lIndices[vtxdata.decay3bodyId()] = vtxdata.globalIndex();
      }
    }
    for (int ii = 0; ii < decay3bodytable.size(); ii++) {
      vtxdataLink(lIndices[ii]);
    }
  }
};

struct kfdecay3bodyDataLinkBuilder {
  Produces<aod::KFDecay3BodyDataLink> kfvtxdataLink;

  void init(InitContext const&) {}

  // build Decay3Body -> KFDecay3BodyData link table
  void process(aod::Decay3Bodys const& decay3bodytable, aod::KFVtx3BodyDatas const& vtxdatatable)
  {
    std::vector<int> lIndices;
    lIndices.reserve(decay3bodytable.size());
    for (int ii = 0; ii < decay3bodytable.size(); ii++)
      lIndices[ii] = -1;
    for (auto& vtxdata : vtxdatatable) {
      lIndices[vtxdata.decay3bodyId()] = vtxdata.globalIndex();
    }
    for (int ii = 0; ii < decay3bodytable.size(); ii++) {
      kfvtxdataLink(lIndices[ii]);
    }
  }
};

struct decay3bodyLabelBuilder {

  Produces<aod::McVtx3BodyLabels> vtxlabels;
  Produces<aod::McFullVtx3BodyLabels> vtxfulllabels;
  Produces<aod::McKFVtx3BodyLabels> kfvtxlabels;
  Produces<aod::McFullKFVtx3BodyLabels> kfvtxfulllabels;

  HistogramRegistry registry{"registry", {}};

  void init(InitContext const&)
  {
    if (doprocessDoNotBuildLabels == false) {
      auto hLabelCounter = registry.add<TH1>("hLabelCounter", "hLabelCounter", HistType::kTH1D, {{3, 0.0f, 3.0f}});
      hLabelCounter->GetXaxis()->SetBinLabel(1, "Total");
      hLabelCounter->GetXaxis()->SetBinLabel(2, "Have Same MotherTrack");
      hLabelCounter->GetXaxis()->SetBinLabel(3, "True H3L");

      registry.add("hHypertritonMCPt", "hHypertritonMCPt", HistType::kTH1F, {{100, 0.0f, 10.0f}});
      registry.add("hAntiHypertritonMCPt", "hAntiHypertritonMCPt", HistType::kTH1F, {{100, 0.0f, 10.0f}});
      registry.add("hHypertritonMCMass", "hHypertritonMCMass", HistType::kTH1F, {{40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}});
      registry.add("hAntiHypertritonMCMass", "hAntiHypertritonMCMass", HistType::kTH1F, {{40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}});
      registry.add("hHypertritonMCLifetime", "hHypertritonMCLifetime", HistType::kTH1F, {{50, 0.0f, 50.0f, "ct(cm)"}});
      registry.add("hAntiHypertritonMCLifetime", "hAntiHypertritonMCLifetime", HistType::kTH1F, {{50, 0.0f, 50.0f, "ct(cm)"}});
    }
  }

  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};

  void processDoNotBuildLabels(aod::Collisions::iterator const&)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(decay3bodyLabelBuilder, processDoNotBuildLabels, "Do not produce MC label tables", true);

  void processBuildLabels(aod::Decay3BodysLinked const& decay3bodys, aod::Vtx3BodyDatas const& vtx3bodydatas, MCLabeledTracksIU const&, aod::McParticles const&)
  {
    std::vector<int> lIndices;
    lIndices.reserve(vtx3bodydatas.size());
    for (int ii = 0; ii < vtx3bodydatas.size(); ii++) {
      lIndices[ii] = -1;
    }

    for (auto& decay3body : decay3bodys) {

      int lLabel = -1;
      int lPDG = -1;
      float lPt = -1;
      double MClifetime = -1;
      bool is3bodyDecay = false;
      int lGlobalIndex = -1;

      auto lTrack0 = decay3body.track0_as<MCLabeledTracksIU>();
      auto lTrack1 = decay3body.track1_as<MCLabeledTracksIU>();
      auto lTrack2 = decay3body.track2_as<MCLabeledTracksIU>();
      registry.fill(HIST("hLabelCounter"), 0.5);

      // Association check
      // There might be smarter ways of doing this in the future
      if (!lTrack0.has_mcParticle() || !lTrack1.has_mcParticle() || !lTrack2.has_mcParticle()) {
        vtxfulllabels(-1);
        continue;
      }
      auto lMCTrack0 = lTrack0.mcParticle_as<aod::McParticles>();
      auto lMCTrack1 = lTrack1.mcParticle_as<aod::McParticles>();
      auto lMCTrack2 = lTrack2.mcParticle_as<aod::McParticles>();
      if (!lMCTrack0.has_mothers() || !lMCTrack1.has_mothers() || !lMCTrack2.has_mothers()) {
        vtxfulllabels(-1);
        continue;
      }

      for (auto& lMother0 : lMCTrack0.mothers_as<aod::McParticles>()) {
        for (auto& lMother1 : lMCTrack1.mothers_as<aod::McParticles>()) {
          for (auto& lMother2 : lMCTrack2.mothers_as<aod::McParticles>()) {
            if (lMother0.globalIndex() == lMother1.globalIndex() && lMother0.globalIndex() == lMother2.globalIndex()) {
              lGlobalIndex = lMother1.globalIndex();
              lPt = lMother1.pt();
              lPDG = lMother1.pdgCode();
              MClifetime = RecoDecay::sqrtSumOfSquares(lMCTrack2.vx() - lMother2.vx(), lMCTrack2.vy() - lMother2.vy(), lMCTrack2.vz() - lMother2.vz()) * o2::constants::physics::MassHyperTriton / lMother2.p(); // only for hypertriton
              is3bodyDecay = true;                                                                                                                                                                               // vtxs with the same mother
            }
          }
        }
      } // end association check
      if (!is3bodyDecay) {
        vtxfulllabels(-1);
        continue;
      }
      registry.fill(HIST("hLabelCounter"), 1.5);

      // Intended for hypertriton cross-checks only
      if (lPDG == 1010010030 && lMCTrack0.pdgCode() == 2212 && lMCTrack1.pdgCode() == -211 && lMCTrack2.pdgCode() == 1000010020) {
        lLabel = lGlobalIndex;
        double hypertritonMCMass = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
        registry.fill(HIST("hLabelCounter"), 2.5);
        registry.fill(HIST("hHypertritonMCPt"), lPt);
        registry.fill(HIST("hHypertritonMCLifetime"), MClifetime);
        registry.fill(HIST("hHypertritonMCMass"), hypertritonMCMass);
      }
      if (lPDG == -1010010030 && lMCTrack0.pdgCode() == 211 && lMCTrack1.pdgCode() == -2212 && lMCTrack2.pdgCode() == -1000010020) {
        lLabel = lGlobalIndex;
        double antiHypertritonMCMass = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron});
        registry.fill(HIST("hLabelCounter"), 2.5);
        registry.fill(HIST("hAntiHypertritonMCPt"), lPt);
        registry.fill(HIST("hAntiHypertritonMCLifetime"), MClifetime);
        registry.fill(HIST("hAntiHypertritonMCMass"), antiHypertritonMCMass);
      }

      // Construct label table, only vtx which corresponds to true mother and true daughters with a specified order is labeled
      // for matter: track0->p, track1->pi, track2->bachelor
      // for antimatter: track0->pi, track1->p, track2->bachelor
      vtxfulllabels(lLabel);
      if (decay3body.vtx3BodyDataId() != -1) {
        lIndices[decay3body.vtx3BodyDataId()] = lLabel;
      }
    }
    for (int ii = 0; ii < vtx3bodydatas.size(); ii++) {
      vtxlabels(lIndices[ii]);
    }
  }
  PROCESS_SWITCH(decay3bodyLabelBuilder, processBuildLabels, "Produce MC label tables", false);

  void processBuildKFLabels(aod::KFDecay3BodysLinked const& decay3bodys, aod::KFVtx3BodyDatas const& vtx3bodydatas, MCLabeledTracksIU const&, aod::McParticles const&)
  {
    std::vector<int> lIndices;
    lIndices.reserve(vtx3bodydatas.size());
    for (int ii = 0; ii < vtx3bodydatas.size(); ii++) {
      lIndices[ii] = -1;
    }

    for (auto& decay3body : decay3bodys) {

      int lLabel = -1;

      auto lTrack0 = decay3body.track0_as<MCLabeledTracksIU>();
      auto lTrack1 = decay3body.track1_as<MCLabeledTracksIU>();
      auto lTrack2 = decay3body.track2_as<MCLabeledTracksIU>();

      // counter total
      registry.fill(HIST("hLabelCounter"), 0.5);

      // Association check
      if (lTrack0.has_mcParticle() && lTrack1.has_mcParticle() && lTrack2.has_mcParticle()) {
        auto lMCTrack0 = lTrack0.mcParticle_as<aod::McParticles>();
        auto lMCTrack1 = lTrack1.mcParticle_as<aod::McParticles>();
        auto lMCTrack2 = lTrack2.mcParticle_as<aod::McParticles>();
        // check if mother is the same
        if (lMCTrack0.has_mothers() && lMCTrack1.has_mothers() && lMCTrack2.has_mothers()) {
          for (auto& lMother0 : lMCTrack0.mothers_as<aod::McParticles>()) {
            for (auto& lMother1 : lMCTrack1.mothers_as<aod::McParticles>()) {
              for (auto& lMother2 : lMCTrack2.mothers_as<aod::McParticles>()) {
                if (lMother0.globalIndex() == lMother1.globalIndex() && lMother0.globalIndex() == lMother2.globalIndex()) {
                  lLabel = lMother1.globalIndex();
                  // fill counter same mother
                  registry.fill(HIST("hLabelCounter"), 1.5);
                } // end same mother conditional
              }
            }
          } // end loop over daughters
        } // end conditional of mothers existing
      } // end association check

      // Construct label table, only vtx which corresponds to true mother and true daughters with a specified order is labeled
      // for matter: track0->p, track1->pi, track2->bachelor
      // for antimatter: track0->pi, track1->p, track2->bachelor
      kfvtxfulllabels(lLabel);
      if (decay3body.kfvtx3BodyDataId() != -1) {
        lIndices[decay3body.kfvtx3BodyDataId()] = lLabel;
      }
    }
    for (int ii = 0; ii < vtx3bodydatas.size(); ii++) {
      kfvtxlabels(lIndices[ii]);
    }
  }
  PROCESS_SWITCH(decay3bodyLabelBuilder, processBuildKFLabels, "Produce MC KF label tables", false);
};

struct decay3bodyInitializer {
  Spawns<aod::Vtx3BodyDatas> vtx3bodydatas;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<decay3bodyBuilder>(cfgc),
    adaptAnalysisTask<decay3bodyDataLinkBuilder>(cfgc),
    adaptAnalysisTask<kfdecay3bodyDataLinkBuilder>(cfgc),
    adaptAnalysisTask<decay3bodyLabelBuilder>(cfgc),
    adaptAnalysisTask<decay3bodyInitializer>(cfgc),
  };
}
