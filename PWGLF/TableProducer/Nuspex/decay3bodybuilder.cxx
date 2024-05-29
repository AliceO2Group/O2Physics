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

using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using FullTracksExtPIDIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe>;
using MCLabeledTracksIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

struct decay3bodyBuilder {

  Produces<aod::StoredVtx3BodyDatas> vtx3bodydata;
  Produces<aod::StoredKFVtx3BodyDatas> kfvtx3bodydata;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // Configurables
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};

  enum vtxstep { kVtxAll = 0,
                 kVtxTPCNcls,
                 kVtxhasSV,
                 kVtxDcaDau,
                 kVtxCosPA,
                 kNVtxSteps };

  enum kfvtxstep {kKfVtxAll = 0,
                  kKfVtxCollIds,
                  kKfVtxCharge,
                  kKfVtxTPCNcls,
                  kKfVtxTPCRows,
                  kKfVtxDCAxyPV,
                  kKfVtxDCAzPV,
                  kKfVtxTPCPID,
                  kKfVtxhasSV,
                  kKfVtxDcaDau,
                  kKfVtxDcaDauVtx,
                  kKfVtxPt,
                  kKfVtxMass,
                  kKfVtxCosPA,
                  kKfVtxCosPAXY,
                  kKfVtxChi2geo,
                  kKfVtxChi2topo,
                  kKfNVtxSteps};

  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
      {"hVtx3BodyCounter", "hVtx3BodyCounter", {HistType::kTH1F, {{5, 0.0f, 5.0f}}}},
      {"hVtx3BodyCounterKFParticle", "hVtx3BodyCounterKFParticle", {HistType::kTH1F, {{17, 0.0f, 17.0f}}}},
      {"hBachelorTOFNSigmaDe", "hBachelorTOFNSigmaDe", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},
    },
  };

  // hypothesis
  Configurable<int> bachelorcharge{"bachelorcharge", 1, "charge of the bachelor track"};
  // Selection criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<int> mintpcNCls{"mintpcNCls", 70, "min tpc Nclusters"};
  Configurable<float> minCosPA3body{"minCosPA3body", 0.9, "minCosPA3body"};
  Configurable<float> dcavtxdau{"dcavtxdau", 1.0, "DCA Vtx Daughters"};

  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  // CCDB TOF PID paras
  Configurable<std::string> paramFileName{"paramFileName", "", "Path to the parametrization object. If empty the parametrization is not taken from file"};
  Configurable<std::string> parametrizationPath{"parametrizationPath", "TOF/Calib/Params", "Path of the TOF parametrization on the CCDB or in the file, if the paramFileName is not empty"};
  Configurable<std::string> passName{"passName", "", "Name of the pass inside of the CCDB parameter collection. If empty, the automatically deceted from metadata (to be implemented!!!)"};
  Configurable<std::string> timeShiftCCDBPath{"timeShiftCCDBPath", "", "Path of the TOF time shift vs eta. If empty none is taken"};
  Configurable<bool> loadResponseFromCCDB{"loadResponseFromCCDB", false, "Flag to load the response from the CCDB"};
  Configurable<bool> enableTimeDependentResponse{"enableTimeDependentResponse", false, "Flag to use the collision timestamp to fetch the PID Response"};
  Configurable<bool> fatalOnPassNotAvailable{"fatalOnPassNotAvailable", true, "Flag to throw a fatal if the pass is not available in the retrieved CCDB object"};
  // for KFParticle reconstruction
  struct : ConfigurableGroup {
    Configurable<bool> kfDoDCAFitterPreMinimum{"kfDoDCAFitterPreMinimum", true, "KF: do DCAFitter pre-optimization before KF fit to include material corrections for decay3body vertex"};
    Configurable<float> mintpcCrossedRows{"mintpcCrossedRows", 70, "Minimum number of TPC crossed rows for daughter tracks"};
    Configurable<float> mindcaXYTrackPV{"mindcaXYTrackPV", 2., "Minimum DCA XY of the daughter tracks to the PV"};
    Configurable<float> mindcaZTrackPV{"mindcaZTrackPV", 10., "Minimum DCA Z of the daughter tracks to the PV"};
    Configurable<float> maxtpcnSigma{"maxtpcnSigma", 5., "Maximum nSigma TPC for daughter tracks"};
    Configurable<float> lambdaMassWindow{"lambdaMassWindow", 0.01, "Window cut around lambda mass for proton-pion vertex with KFParticle"};
    Configurable<float> maxDca3dDau{"maxDca3dDau", 1000., "Maximum geometrical distance between daughter tracks at the SV in 3D with KFParticle"};
    Configurable<float> maxDcaXYSVDau{"maxDcaXYSVDau", 1.0, "Maximum geometrical distance of daughter tracks from the SV in XY with KFParticle"};
    Configurable<float> minPtHt{"minPtHt", 0., "Minimum momentum for Hypertriton candidates with KFParticle"};
    Configurable<float> maxPtHt{"maxPtHt", 36., "Maximum momentum for Hypertriton candidates with KFParticle"};
    Configurable<float> minMassHt{"minMassHt", 2.96, "Minimum candidate mass with KFParticle"};
    Configurable<float> maxMassHt{"maxMassHt", 3.04, "Maximum candidate mass with KFParticle"};
    Configurable<float> maxChi2geo{"maxChi2geo", 1000., "Maximum chi2 geometrical with KFParticle"};
    Configurable<float> minCosPA{"minCosPA", 0.5, "Minimum cosine pointing angle with KFParticle"};
    Configurable<float> minCosPAxy{"minCosPAxy", 0.5, "Minimum cosine pointing angle in xy with KFParticle"};
    Configurable<bool> applyTopoSels{"applyTopoSels", false, "Apply selections constraining the mother to the PV with KFParticle"};
    Configurable<float> maxChi2topo{"maxChi2topo", 1000., "Maximum chi2 topological with KFParticle"};
  } kfparticleConfigurations;
  

  // Filters and slices
  Filter collisionFilter = (aod::evsel::sel8 == true && nabs(aod::collision::posZ) < 10.f);
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

    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(2, "TPCNcls");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(3, "HasSV");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(4, "DcaDau");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(5, "CosPA");

    // Material correction in the DCA fitter
    if (useMatCorrType == 1)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    if (useMatCorrType == 2)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

    fitter3body.setMatCorrType(matCorr);
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
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
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
    ccdb->setTimestamp(bc.timestamp());
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
      LOG(info) << "Loading exp. sigma parametrization from CCDB, using path: " << parametrizationPath.value << " for timestamp " << bc.timestamp();
      o2::tof::ParameterCollection* paramCollection = ccdb->getForTimeStamp<o2::tof::ParameterCollection>(parametrizationPath.value, bc.timestamp());
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
        mRespParamsV2.setTimeShiftParameters(ccdb->getForTimeStamp<TGraph>(Form("%s/pos", timeShiftCCDBPath.value.c_str()), bc.timestamp()), true);
        mRespParamsV2.setTimeShiftParameters(ccdb->getForTimeStamp<TGraph>(Form("%s/neg", timeShiftCCDBPath.value.c_str()), bc.timestamp()), false);
      }
    }
  }

  //------------------------------------------------------------------
  // Recalculate TOF PID for bachelors (deuteron), copied from PIDTOF.h
  template <typename TrackType>
  static float GetExpectedSigma(const o2::pid::tof::TOFResoParamsV2& parameters, const TrackType& track, const float tofSignal, const float collisionTimeRes, double mMassZ)
  {
    const float& mom = track.p();
    if (mom <= 0) {
      return -999.f;
    }
    const float dpp = parameters[9] + parameters[10] * mom + parameters[11] * mMassZ / mom; // mean relative pt resolution;
    const float sigma = dpp * tofSignal / (1. + mom * mom / (mMassZ * mMassZ));
    return std::sqrt(sigma * sigma + parameters[12] * parameters[12] / mom / mom + parameters[4] * parameters[4] + collisionTimeRes * collisionTimeRes);
  }

  //------------------------------------------------------------------
  // 3body candidate builder
  template <class TTrackClass, typename TCollisionTable, typename TTrackTable>
  void buildVtx3BodyDataTable(TCollisionTable const& collision, TTrackTable const& /*tracks*/, aod::Decay3Bodys const& decay3bodys, int bachelorcharge = 1)
  {

    for (auto& vtx3body : decay3bodys) {

      registry.fill(HIST("hVtx3BodyCounter"), kVtxAll);

      auto t0 = vtx3body.template track0_as<TTrackClass>();
      auto t1 = vtx3body.template track1_as<TTrackClass>();
      auto t2 = vtx3body.template track2_as<TTrackClass>();

      if (t0.tpcNClsFound() < mintpcNCls && t1.tpcNClsFound() < mintpcNCls && t2.tpcNClsFound() < mintpcNCls) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), kVtxTPCNcls);

      // Calculate DCA with respect to the collision associated to the V0, not individual tracks
      gpu::gpustd::array<float, 2> dcaInfo;

      auto Track0Par = getTrackPar(t0);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, Track0Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
      auto Track0dcaXY = dcaInfo[0];

      auto Track1Par = getTrackPar(t1);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, Track1Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
      auto Track1dcaXY = dcaInfo[0];

      auto Track2Par = getTrackPar(t2);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, Track2Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
      auto Track2dcaXY = dcaInfo[0];

      auto Track0 = getTrackParCov(t0);
      auto Track1 = getTrackParCov(t1);
      auto Track2 = getTrackParCov(t2);
      int n3bodyVtx = fitter3body.process(Track0, Track1, Track2);
      if (n3bodyVtx == 0) { // discard this pair
        continue;
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
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), kVtxDcaDau);

      float VtxcosPA = RecoDecay::cpa(array{collision.posX(), collision.posY(), collision.posZ()}, array{pos[0], pos[1], pos[2]}, array{p3B[0], p3B[1], p3B[2]});
      if (VtxcosPA < minCosPA3body) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), kVtxCosPA);

      // Recalculate the TOF PID
      double tofNsigmaDe = -999;
      static constexpr float kCSPEED = TMath::C() * 1.0e2f * 1.0e-12f; // c in cm/ps

      if (t2.hasTOF()) {
        double bachExpTime = t2.length() * sqrt((o2::constants::physics::MassDeuteron * o2::constants::physics::MassDeuteron) + (t2.tofExpMom() * t2.tofExpMom())) / (kCSPEED * t2.tofExpMom()); // L*E/(p*c) = L/v
        double tofsignal = t2.trackTime() * 1000 + bachExpTime;
        // double bachtime = t2.trackTime() * 1000 + bachExpTime - collision.collisionTime(); // in ps

        double expSigma = GetExpectedSigma(mRespParamsV2, t2, tofsignal, collision.collisionTimeRes(), o2::constants::physics::MassDeuteron);
        double corrTofMom = t2.tofExpMom() / (1.f + t2.sign() * mRespParamsV2.getShift(t2.eta()));
        double corrSignal = t2.length() * sqrt((o2::constants::physics::MassDeuteron * o2::constants::physics::MassDeuteron) + (corrTofMom * corrTofMom)) / (kCSPEED * corrTofMom) + mRespParamsV2.getTimeShift(t2.eta(), t2.sign());
        tofNsigmaDe = (tofsignal - collision.collisionTime() - corrSignal) / expSigma;
      }

      registry.fill(HIST("hBachelorTOFNSigmaDe"), t2.sign() * t2.p(), tofNsigmaDe);

      vtx3bodydata(
        t0.globalIndex(), t1.globalIndex(), t2.globalIndex(), collision.globalIndex(), vtx3body.globalIndex(),
        pos[0], pos[1], pos[2],
        p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2],
        fitter3body.getChi2AtPCACandidate(),
        Track0dcaXY, Track1dcaXY, Track2dcaXY,
        tofNsigmaDe);
    }
  }

  //------------------------------------------------------------------
  // function to select daughter track PID
  template <typename TTrack>
  bool selectTPCPID(TTrack const& trackProton, TTrack const& trackPion, TTrack const& trackDeuteron)
  {
    if (!(trackProton.hasTPC() && abs(trackProton.tpcNSigmaPr()) < kfparticleConfigurations.maxtpcnSigma)) {
      return false;
    }
    if (!(trackPion.hasTPC() && abs(trackPion.tpcNSigmaPi()) < kfparticleConfigurations.maxtpcnSigma)) {
      return false;
    }
    if (!(trackDeuteron.hasTPC() && abs(trackDeuteron.tpcNSigmaDe()) < kfparticleConfigurations.maxtpcnSigma)) {
      return false;
    }
    return true;
  }

  //------------------------------------------------------------------
  // 3body candidate builder with KFParticle
  template <class TTrackTo, typename TCollision>
  void buildVtx3BodyDataTableKFParticle(TCollision const& collision, aod::Decay3Bodys const& decay3bodys, int bachelorcharge)
  {

    for (auto& vtx3body : decay3bodys) {

      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxAll);

      auto trackPos = vtx3body.template track0_as<TTrackTo>();
      auto trackNeg = vtx3body.template track1_as<TTrackTo>();
      auto trackBach = vtx3body.template track2_as<TTrackTo>();
      auto trackParCovPos = getTrackParCov(trackPos);
      auto trackParCovNeg = getTrackParCov(trackNeg);
      auto trackParCovBach = getTrackParCov(trackBach);

      KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
      KFParticle kfpv(kfpVertex);

      bool isMatter = trackBach.sign() > 0 ? true : false;

      // -------- STEP 1: track selection --------
      // collision ID
      if (trackPos.collisionId() != trackNeg.collisionId() || trackPos.collisionId() != trackBach.collisionId() || trackNeg.collisionId() != trackBach.collisionId()) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxCollIds);
      // track IDs --> already checked in SVertexer!
      // track signs (pos, neg, bach) --> sanity check, should already be in SVertex
      if (!(trackPos.sign() == 1 && trackNeg.sign() == -1 && abs(trackBach.sign() == 1))) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxCharge);
      // number of TPC clusters
      if (trackPos.tpcNClsFound() < mintpcNCls || trackNeg.tpcNClsFound() < mintpcNCls || trackBach.tpcNClsFound() < mintpcNCls) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxTPCNcls);
      // number of TPC crossed rows
      if (trackPos.tpcNClsCrossedRows() < kfparticleConfigurations.mintpcCrossedRows || trackNeg.tpcNClsCrossedRows() < kfparticleConfigurations.mintpcCrossedRows || trackBach.tpcNClsCrossedRows() < kfparticleConfigurations.mintpcCrossedRows) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxTPCRows);
      // track DCAxy and DCAz to PV associated with decay3body
      o2::dataformats::VertexBase mPV;
      o2::dataformats::DCA mDcaInfoCovPos;
      o2::dataformats::DCA mDcaInfoCovNeg;
      o2::dataformats::DCA mDcaInfoCovBach;
      auto trackParCovPVPos = trackParCovPos;
      auto trackParCovPVNeg = trackParCovNeg;
      auto trackParCovPVBach = trackParCovBach;
      mPV.setPos({collision.posX(), collision.posY(), collision.posZ()});
      mPV.setCov(collision.covXX(), collision.covXX(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
      o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, trackParCovPVPos, 2.f, matCorr, &mDcaInfoCovPos);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, trackParCovPVNeg, 2.f, matCorr, &mDcaInfoCovNeg);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, trackParCovPVBach, 2.f, matCorr, &mDcaInfoCovBach);
      auto TrackPosDcaXY = mDcaInfoCovPos.getY();
      auto TrackNegDcaXY = mDcaInfoCovNeg.getY();
      auto TrackBachDcaXY = mDcaInfoCovBach.getY();
      if (TrackPosDcaXY < kfparticleConfigurations.mindcaXYTrackPV || TrackNegDcaXY < kfparticleConfigurations.mindcaXYTrackPV || TrackBachDcaXY < kfparticleConfigurations.mindcaXYTrackPV) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxDCAxyPV);
      if (fabs(mDcaInfoCovPos.getZ()) < kfparticleConfigurations.mindcaZTrackPV || fabs(mDcaInfoCovNeg.getZ()) < kfparticleConfigurations.mindcaZTrackPV || fabs(mDcaInfoCovBach.getZ()) < kfparticleConfigurations.mindcaZTrackPV) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxDCAzPV);
      // TPC PID
      if (isMatter && !selectTPCPID(trackPos, trackNeg, trackBach)) { // hypertriton (proton, pi-, deuteron)
        continue;
      } else if (!isMatter && !selectTPCPID(trackNeg, trackPos, trackBach)) { // anti-hypertriton (pi+, anti-proton, deuteron)
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxTPCPID);

      // -------- STEP 2: fit vertex with proton and pion --------
      // Fit vertex with DCA fitter to find minimization point --> uses material corrections implicitly
      if (kfparticleConfigurations.kfDoDCAFitterPreMinimum) {
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
      // Construct V0
      KFParticle KFV0;
      int nDaughters = 2;
      const KFParticle* Daughters[2] = {&kfpProton, &kfpPion};
      KFV0.SetConstructMethod(2);
      try {
        KFV0.Construct(Daughters, nDaughters);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to create V0 vertex from daughter tracks." << e.what();
        continue;
      }
      KFV0.TransportToDecayVertex();
      // check V0 mass and set mass constraint
      float massV0, sigmaMassV0;
      KFV0.GetMass(massV0, sigmaMassV0);
      if (abs(massV0 - constants::physics::MassLambda) < kfparticleConfigurations.lambdaMassWindow) {
        continue;
      }
      KFParticle KFV0Mass = KFV0;
      KFV0Mass.SetNonlinearMassConstraint(o2::constants::physics::MassLambda);
      float chi2massV0 = KFV0Mass.GetChi2() / KFV0Mass.GetNDF();

      // -------- STEP 3: fit vertex with V0 and deuteron --------
      // Create KFParticle object from deuteron track
      KFParticle kfpDeuteron;
      kfpDeuteron = createKFParticleFromTrackParCov(trackParCovBach, trackBach.sign()*bachelorcharge, constants::physics::MassDeuteron);
      // Add deuteron to V0 vertex
      KFParticle KFHt;
      KFHt = KFV0;
      KFHt.SetConstructMethod(2);
      try {
        KFHt.AddDaughter(kfpDeuteron);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to create Hyper triton from V0 and deuteron." << e.what();
        continue;
      }

      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxhasSV);

      // -------- STEP 4: selections after geometrical vertex fit --------
      // Get updated daughter tracks
      kfpProton.SetProductionVertex(KFHt);
      kfpPion.SetProductionVertex(KFHt);
      kfpDeuteron.SetProductionVertex(KFHt);
      // daughter DCAs
      if ((kfpProton.GetDistanceFromParticle(kfpPion) > kfparticleConfigurations.maxDca3dDau) || (kfpProton.GetDistanceFromParticle(kfpDeuteron) > kfparticleConfigurations.maxDca3dDau) || (kfpPion.GetDistanceFromParticle(kfpDeuteron) > kfparticleConfigurations.maxDca3dDau)) {
        continue;
      }
      float DCAvtxDaughters3D = kfpProton.GetDistanceFromParticle(kfpPion) + kfpProton.GetDistanceFromParticle(kfpDeuteron) + kfpPion.GetDistanceFromParticle(kfpDeuteron);
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kVtxDcaDau);
      // daughter DCAs to vertex
      if (kfpProton.GetDistanceFromVertexXY(KFHt) > kfparticleConfigurations.maxDcaXYSVDau || kfpPion.GetDistanceFromVertexXY(KFHt) > kfparticleConfigurations.maxDcaXYSVDau || kfpDeuteron.GetDistanceFromVertexXY(KFHt) > kfparticleConfigurations.maxDcaXYSVDau) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxDcaDauVtx);

      // -------- STEP 5: candidate selection after geometrical vertex fit --------
      // function to select hyper triton daughters after the geometrical fitting
      // Pt selection
      if (KFHt.GetPt() < kfparticleConfigurations.minPtHt || KFHt.GetPt() > kfparticleConfigurations.maxPtHt) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxPt);
      // Mass window
      float massHt, sigmaMassHt;
      KFHt.GetMass(massHt, sigmaMassHt);
      if (massHt < kfparticleConfigurations.minMassHt || massHt > kfparticleConfigurations.maxMassHt) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxMass);
      // cos(PA) to PV
      if (abs(cpaFromKF(KFHt, kfpv)) < kfparticleConfigurations.minCosPA) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kVtxCosPA);
      // cos(PA) xy to PV
      if (abs(cpaXYFromKF(KFHt, kfpv)) < kfparticleConfigurations.minCosPAxy) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxCosPAXY);
      // chi2 geometrical
      float chi2geoNDF = KFHt.GetChi2() / KFHt.GetNDF();
      if (chi2geoNDF > kfparticleConfigurations.maxChi2geo) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxChi2geo);

      // -------- STEP 6: topological constraint --------
      /// Set vertex constraint and topological selection
      KFParticle KFHtPV = KFHt;
      KFHtPV.SetProductionVertex(kfpv);
      KFHtPV.TransportToDecayVertex();
      float chi2topoNDF = KFHtPV.GetChi2() / KFHtPV.GetNDF();
      if (kfparticleConfigurations.applyTopoSels && chi2topoNDF > kfparticleConfigurations.maxChi2topo) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounterKFParticle"), kKfVtxChi2topo);

      //------------------------------------------------------------------
      // Recalculate the bachelor TOF PID
      double tofNsigmaDe = -999;
      static constexpr float kCSPEED = TMath::C() * 1.0e2f * 1.0e-12f; // c in cm/ps
      if (trackBach.hasTOF()) {
        double bachExpTime = trackBach.length() * sqrt((o2::constants::physics::MassDeuteron * o2::constants::physics::MassDeuteron) + (trackBach.tofExpMom() * trackBach.tofExpMom())) / (kCSPEED * trackBach.tofExpMom()); // L*E/(p*c) = L/v
        double tofsignal = trackBach.trackTime() * 1000 + bachExpTime;
        // double bachtime = trackBach.trackTime() * 1000 + bachExpTime - collision.collisionTime(); // in ps

        double expSigma = GetExpectedSigma(mRespParamsV2, trackBach, tofsignal, collision.collisionTimeRes(), o2::constants::physics::MassDeuteron);
        double corrTofMom = trackBach.tofExpMom() / (1.f + trackBach.sign() * mRespParamsV2.getShift(trackBach.eta()));
        double corrSignal = trackBach.length() * sqrt((o2::constants::physics::MassDeuteron * o2::constants::physics::MassDeuteron) + (corrTofMom * corrTofMom)) / (kCSPEED * corrTofMom) + mRespParamsV2.getTimeShift(trackBach.eta(), trackBach.sign());
        tofNsigmaDe = (tofsignal - collision.collisionTime() - corrSignal) / expSigma;
      }
      registry.fill(HIST("hBachelorTOFNSigmaDe"), trackBach.sign() * trackBach.p(), tofNsigmaDe);

      //------------------------------------------------------------------
      // table filling
      kfvtx3bodydata(
        collision.globalIndex(), trackPos.globalIndex(), trackNeg.globalIndex(), trackBach.globalIndex(), vtx3body.globalIndex(),
        // hypertriton
        massHt,
        KFHt.GetX(), KFHt.GetY(), KFHt.GetZ(),
        KFHt.GetPx(), KFHt.GetPy(), KFHt.GetPz(), KFHt.GetPt(),
        KFHt.GetQ(),
        KFHt.GetDistanceFromVertex(kfpv), KFHt.GetDistanceFromVertexXY(kfpv),
        cpaFromKF(KFHt, kfpv), // before topo constraint
        cpaXYFromKF(KFHt, kfpv),
        cpaFromKF(KFHtPV, kfpv), // after topo constraint
        cpaXYFromKF(KFHtPV, kfpv),
        KFHtPV.GetDecayLength(), KFHtPV.GetDecayLengthXY(), // decay length defined after topological constraint
        KFHtPV.GetDecayLength()/KFHtPV.GetErrDecayLength(), // ldl
        chi2geoNDF, chi2topoNDF,
        // V0
        massV0, chi2massV0,
        // daughter momenta
        kfpProton.GetPx(), kfpProton.GetPy(), kfpProton.GetPz(),
        kfpPion.GetPx(), kfpPion.GetPy(), kfpPion.GetPz(),
        kfpDeuteron.GetPx(), kfpDeuteron.GetPy(), kfpDeuteron.GetPz(),
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
        DCAvtxDaughters3D,
        // daughter DCAs to PV propagated with material
        TrackPosDcaXY, TrackNegDcaXY, TrackBachDcaXY,
        // daughter signs
        trackPos.sign(),
        trackNeg.sign(),
        trackBach.sign(),
        // bachelor TOF PID
        tofNsigmaDe);
    }
  }


  //------------------------------------------------------------------
  void process(aod::Collision const& collision, FullTracksExtIU const& tracksIU, aod::Decay3Bodys const& decay3bodys, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    registry.fill(HIST("hEventCounter"), 0.5);

    buildVtx3BodyDataTable<FullTracksExtIU>(collision, tracksIU, decay3bodys, bachelorcharge);
  }
  PROCESS_SWITCH(decay3bodyBuilder, process, "Produce DCA fitter decay3body tables", true);

  void processKFParticle(soa::Filtered<soa::Join<aod::Collisions, 
                         aod::EvSels>>::iterator const& collision, 
                         FullTracksExtPIDIU const& tracksIU, 
                         aod::Decay3Bodys const& decay3bodys, 
                         aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    registry.fill(HIST("hEventCounter"), 0.5);

    // slice Decay3Body table by collision
    const uint64_t collIdx = collision.globalIndex();
    auto Decay3BodyTable_thisCollision = decay3bodys.sliceBy(perCollision, collIdx);
    buildVtx3BodyDataTableKFParticle<FullTracksExtPIDIU>(collision, decay3bodys, bachelorcharge);
  }
  PROCESS_SWITCH(decay3bodyBuilder, processKFParticle, "Produce KFParticle decay3body tables", false);

};

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
      lIndices[vtxdata.decay3bodyId()] = vtxdata.globalIndex();
    }
    for (int ii = 0; ii < decay3bodytable.size(); ii++) {
      vtxdataLink(lIndices[ii]);
    }
  }
};

struct decay3bodyLabelBuilder {

  Produces<aod::McVtx3BodyLabels> vtxlabels;
  Produces<aod::McFullVtx3BodyLabels> vtxfulllabels;

  // for bookkeeping purposes: how many V0s come from same mother etc
  HistogramRegistry registry{
    "registry",
    {
      {"hLabelCounter", "hLabelCounter", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
      {"hHypertritonMCPt", "hHypertritonMCPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hAntiHypertritonMCPt", "hAntiHypertritonMCPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hHypertritonMCMass", "hHypertritonMCMass", {HistType::kTH1F, {{40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
      {"hAntiHypertritonMCMass", "hAntiHypertritonMCMass", {HistType::kTH1F, {{40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
      {"hHypertritonMCLifetime", "hHypertritonMCLifetime", {HistType::kTH1F, {{50, 0.0f, 50.0f, "ct(cm)"}}}},
      {"hAntiHypertritonMCLifetime", "hAntiHypertritonMCLifetime", {HistType::kTH1F, {{50, 0.0f, 50.0f, "ct(cm)"}}}},
    },
  };

  void init(InitContext const&)
  {
    registry.get<TH1>(HIST("hLabelCounter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hLabelCounter"))->GetXaxis()->SetBinLabel(2, "Have Same MotherTrack");
    registry.get<TH1>(HIST("hLabelCounter"))->GetXaxis()->SetBinLabel(3, "True H3L");
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
    adaptAnalysisTask<decay3bodyLabelBuilder>(cfgc),
    adaptAnalysisTask<decay3bodyInitializer>(cfgc),
  };
}
