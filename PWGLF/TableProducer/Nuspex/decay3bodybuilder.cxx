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
// ========================
/// \file decay3bodybuilder.cxx
/// \brief Builder task for 3-body hypertriton decay reconstruction (proton + pion + deuteron)
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>
/// \author Carolina Reetz <c.reetz@cern.ch>
// ========================

#include "TableHelper.h"

#include "PWGLF/DataModel/LFPIDTOFGenericTables.h"
#include "PWGLF/DataModel/Reduced3BodyTables.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "PWGLF/Utils/decay3bodyBuilderHelper.h"
#include "PWGLF/Utils/pidTOFGeneric.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef HomogeneousField
#define HomogeneousField
#endif

// includes KFParticle
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticle.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

o2::common::core::MetadataHelper metadataInfo;

static constexpr int nParameters = 1;
static const std::vector<std::string> tableNames{
  "Decay3BodyIndices",
  "Vtx3BodyDatas",
  "Vtx3BodyCovs",
  "McVtx3BodyDatas"};

static constexpr int nTablesConst = 4;

static const std::vector<std::string> parameterNames{"enable"};
static const int defaultParameters[nTablesConst][nParameters]{
  {0}, // Decay3BodyIndices
  {0}, // Vtx3BodyDatas
  {0}, // Vtx3BodyCovs
  {0}  // McVtx3BodyDatas
};

using TracksExtPIDIUwithEvTimes = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe, aod::EvTimeTOFFT0ForTrack>;
using TracksExtPIDIUwithEvTimesLabeled = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe, aod::EvTimeTOFFT0ForTrack, aod::McTrackLabels>;

using ColswithEvTimes = o2::soa::Join<aod::Collisions, aod::EvSels, aod::EvTimeTOFFT0>;
using ColswithEvTimesLabeled = o2::soa::Join<aod::Collisions, aod::EvSels, aod::EvTimeTOFFT0, aod::McCollisionLabels>;

struct decay3bodyBuilder {

  // helper object
  o2::pwglf::decay3bodyBuilderHelper helper;

  // table index : match order above
  enum tableIndex { kDecay3BodyIndices = 0,
                    kVtx3BodyDatas,
                    kVtx3BodyCovs,
                    kMcVtx3BodyDatas,
                    nTables };

  struct : ProducesGroup {
    Produces<aod::Decay3BodyIndices> decay3bodyindices;
    Produces<aod::Vtx3BodyDatas> vtx3bodydatas;
    Produces<aod::Vtx3BodyCovs> vtx3bodycovs;
    Produces<aod::McVtx3BodyDatas> mcvtx3bodydatas;
  } products;

  // enablde tables
  Configurable<LabeledArray<int>> enabledTables{"enabledTables",
                                                {defaultParameters[0], nTables, nParameters, tableNames, parameterNames},
                                                "Produce this table: 0 - false, 1 - true"};
  std::vector<int> mEnabledTables; // Vector of enabled tables

  // general options
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
  Configurable<bool> doTrackQA{"doTrackQA", false, "Flag to fill QA histograms for daughter tracks of (selected) decay3body candidates."};
  Configurable<bool> doVertexQA{"doVertexQA", false, "Flag to fill QA histograms for PV of (selected) events."};
  Configurable<bool> doSel8selection{"doSel8selection", true, "flag for sel8 event selection"};
  Configurable<bool> doPosZselection{"doPosZselection", true, "flag for posZ event selection"};

  // data processing options
  Configurable<bool> doSkimmedProcessing{"doSkimmedProcessing", false, "Apply Zoroo counting in case of skimmed data input"};
  Configurable<std::string> triggerList{"triggerList", "fTriggerEventF1Proton, fTrackedOmega, fTrackedXi, fOmegaLargeRadius, fDoubleOmega, fOmegaHighMult, fSingleXiYN, fQuadrupleXi, fDoubleXi, fhadronOmega, fOmegaXi, fTripleXi, fOmega, fGammaVeryLowPtEMCAL, fGammaVeryLowPtDCAL, fGammaHighPtEMCAL, fGammaLowPtEMCAL, fGammaVeryHighPtDCAL, fGammaVeryHighPtEMCAL, fGammaLowPtDCAL, fJetNeutralLowPt, fJetNeutralHighPt, fGammaHighPtDCAL, fJetFullLowPt, fJetFullHighPt, fEMCALReadout, fPCMandEE, fPHOSnbar, fPCMHighPtPhoton, fPHOSPhoton, fLD, fPPPHI, fPD, fLLL, fPLL, fPPL, fPPP, fLeadingPtTrack, fHighFt0cFv0Flat, fHighFt0cFv0Mult, fHighFt0Flat, fHighFt0Mult, fHighMultFv0, fHighTrackMult, fHfSingleNonPromptCharm3P, fHfSingleNonPromptCharm2P, fHfSingleCharm3P, fHfPhotonCharm3P, fHfHighPt2P, fHfSigmaC0K0, fHfDoubleCharm2P, fHfBeauty3P, fHfFemto3P, fHfFemto2P, fHfHighPt3P, fHfSigmaCPPK, fHfDoubleCharm3P, fHfDoubleCharmMix, fHfPhotonCharm2P, fHfV0Charm2P, fHfBeauty4P, fHfV0Charm3P, fHfSingleCharm2P, fHfCharmBarToXiBach, fSingleMuHigh, fSingleMuLow, fLMeeHMR, fDiMuon, fDiElectron, fLMeeIMR, fSingleE, fTrackHighPt, fTrackLowPt, fJetChHighPt, fJetChLowPt, fUDdiffLarge, fUDdiffSmall, fITSextremeIonisation, fITSmildIonisation, fH3L3Body, fHe, fH2", "List of triggers used to select events"};
  Configurable<bool> onlyKeepInterestedTrigger{"onlyKeepInterestedTrigger", false, "Flag to keep only interested trigger"};
  Configurable<bool> doLikeSign{"doLikeSign", false, "Flag to produce like-sign background. If true, require the sign of pion is as same as deuteron but not proton."};

  // CCDB options
  struct : ConfigurableGroup {
    std::string prefix = "ccdb";
    Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  } ccdbConfigurations;

  // Decay3body building options
  struct : ConfigurableGroup {
    std::string prefix = "decay3bodyBuilderOpts";
    // building options
    Configurable<bool> useKFParticle{"useKFParticle", false, "Use KFParticle for decay3body building"};
    Configurable<bool> kfSetTopologicalConstraint{"kfSetTopologicalConstraint", false, "Set topological vertex constraint in case of KFParticle reconstruction"};
    Configurable<bool> buildOnlyTracked{"buildOnlyTracked", false, "Build only tracked decay3bodys"};
    Configurable<bool> useSelections{"useSelections", true, "Apply selections during decay3body building"};
    Configurable<bool> useChi2Selection{"useChi2Selection", true, "Apply chi2 selection during decay3body building"};
    Configurable<bool> useTPCforPion{"useTPCforPion", false, "Flag to ask for TPC info for pion track (PID, nClusters), false: pion track can be ITS only"};
    Configurable<bool> acceptTPCOnly{"acceptTPCOnly", false, "Accept TPC only tracks as daughters"};
    Configurable<bool> askOnlyITSMatch{"askOnlyITSMatch", true, "ask only ITS match to distinguish TPC only tracks"};
    Configurable<bool> calculateCovariance{"calculateCovariance", true, "Calculate candidate and daughter covariance matrices"};
    // daughter track selections
    Configurable<float> maxEtaDaughters{"maxEtaDaughters", 0.9, "Max eta of daughters"};
    Configurable<int> minTPCNClProton{"minTPCNClProton", 90, "Min TPC NClusters of proton daughter"};
    Configurable<int> minTPCNClPion{"minTPCNClPion", 70, "Min TPC NClusters of pion daughter"};
    Configurable<int> minTPCNClDeuteron{"minTPCNClDeuteron", 100, "Min TPC NClusters of deuteron daughter"};
    Configurable<float> minDCAProtonToPV{"minDCAProtonToPV", 0.1, "Min DCA of proton to PV"};
    Configurable<float> minDCAPionToPV{"minDCAPionToPV", 0.1, "Min DCA of pion to PV"};
    Configurable<float> minDCADeuteronToPV{"minDCADeuteronToPV", 0.1, "Min DCA of deuteron to PV"};
    Configurable<float> minPtProton{"minPtProton", 0.3, "Min Pt of proton daughter"};
    Configurable<float> minPtPion{"minPtPion", 0.1, "Min Pt of pion daughter"};
    Configurable<float> minPtDeuteron{"minPtDeuteron", 0.6, "Min Pt of deuteron daughter"};
    Configurable<float> maxPtProton{"maxPtProton", 5.0, "Max Pt of proton daughter"};
    Configurable<float> maxPtPion{"maxPtPion", 1.2, "Max Pt of pion daughter"};
    Configurable<float> maxPtDeuteron{"maxPtDeuteron", 10.0, "Max Pt of deuteron daughter"};
    Configurable<float> maxTPCnSigma{"maxTPCnSigma", 5.0, "Min/max TPC nSigma of daughter tracks"};
    Configurable<float> minTOFnSigmaDeuteron{"minTOFnSigmaDeuteron", -5.0, "Min TOF nSigma of deuteron daughter"};
    Configurable<float> maxTOFnSigmaDeuteron{"maxTOFnSigmaDeuteron", 5.0, "Max TOF nSigma of deuteron daughter"};
    Configurable<float> minPDeuteronUseTOF{"minPDeuteronUseTOF", 1.0, "Min P of deuteron to use TOF PID"};
    Configurable<float> maxDCADauToSVaverage{"maxDCADauToSVaverage", 0.5, "Max DCA of daughters to SV (quadratic sum of daughter DCAs to SV / 3)"};
    // candidate selections
    Configurable<float> maxRapidity{"maxRapidity", 1.0, "Max rapidity of decay3body vertex"};
    Configurable<float> minPt{"minPt", 2.0, "Min Pt of decay3body candidate"};
    Configurable<float> maxPt{"maxPt", 5.0, "Max Pt of decay3body candidate"};
    Configurable<float> minMass{"minMass", 2.96, "Min mass of decay3body candidate"};
    Configurable<float> maxMass{"maxMass", 3.04, "Max mass of decay3body candidate"};
    Configurable<float> minCtau{"minCtau", 0.0, "Min ctau of decay3body candidate"};
    Configurable<float> maxCtau{"maxCtau", 100.0, "Max ctau of decay3body candidate"};
    Configurable<float> minCosPA{"minCosPA", 0.9, "Min cosPA of decay3body candidate"};
    Configurable<float> maxChi2{"maxChi2", 100.0, "Max chi2 of decay3body candidate"};
  } decay3bodyBuilderOpts;

  struct : ConfigurableGroup {
    std::string prefix = "mixingOpts";
    Configurable<int> n3bodyMixing{"n3bodyMixing", 0, "Number of decay3bodys to mix: 0 - value set to maximum bin entry in hDecay3BodyRadiusPhi, > 0 - manual setting"};
    Configurable<int> mixingType{"mixingType", 0, "0: mix V0 from one event with bachelor from another, 1: mix pion and bachelor from one event with proton from another, 1: mix proton and bachelor from one event with pion from another "};
    ConfigurableAxis bins3BodyRadius{"mixingOpts.bins3BodyRadius", {VARIABLE_WIDTH, 0.0f, 2.0f, 4.0f, 7.0f, 10.0f, 14.0f, 18.0f, 22.0f, 30.0f, 40.0f}, "Mixing bins - 3body radius"};
    ConfigurableAxis bins3BodyPhi{"mixingOpts.bins3BodyPhi", {VARIABLE_WIDTH, -180 * o2::constants::math::Deg2Rad, -120 * o2::constants::math::Deg2Rad, -60 * o2::constants::math::Deg2Rad, 0, 60 * o2::constants::math::Deg2Rad, 120 * o2::constants::math::Deg2Rad, 180 * o2::constants::math::Deg2Rad}, "Mixing bins - 3body phi (rad)"};
    ConfigurableAxis bins3BodyPhiDegree{"mixingOpts.bins3BodyPhiDegree", {VARIABLE_WIDTH, -180, -120, -60, 0, 60, 120, 180}, "Mixing bins - 3body phi (degree)"};
    ConfigurableAxis bins3BodyPosZ{"mixingOpts.bins3BodyPosZ", {VARIABLE_WIDTH, -500.0f, -200.0f, -100.0f, -70.0f, -60.0f, -50.0f, -40.0f, -35.0f, -30.0f, -25.0f, -20.0f, -15.0f, -13.0f, -10.0f, -8.0f, -6.0f, -4.0f, -2.0f, 0.0f, 2.0f, 4.0f, 6.0f, 8.0f, 10.0f, 13.0f, 15.0f, 20.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f, 60.0f, 70.0f, 100.0f, 200.0f, 500.0f}, "3body SV z position"};
    Configurable<bool> selectPVPosZ3bodyMixing{"selectPVPosZ3bodyMixing", true, "Select same pvPosZ events in case of 3body mixing"};
    Configurable<float> maxDeltaPVPosZ3bodyMixing{"maxDeltaPVPosZ3bodyMixing", 1., "max difference between PV z position in case of 3body mixing"};
    // SVertexer selections
    Configurable<bool> doApplySVertexerCuts{"doApplySVertexerCuts", false, "Apply SVertexer selections during event mixing"};
    Configurable<float> minPt2V0{"minPt2V0", 0.5, "Min Pt squared of V0"};
    Configurable<float> maxTgl2V0{"maxTgl2V0", 4, "Max tgl squared of V0"};
    Configurable<float> maxDCAXY2ToMeanVertex3bodyV0{"maxDCAXY2ToMeanVertex3bodyV0", 4, "Max DCA XY squared of V0 to mean vertex"};
    Configurable<float> minCosPAXYMeanVertex3bodyV0{"minCosPAXYMeanVertex3bodyV0", 0.9, "Min cosPA XY of V0 to mean vertex"};
    Configurable<float> minCosPA3bodyV0{"minCosPA3bodyV0", 0.8, "Min cosPA of V0"};
    Configurable<float> maxRDiffV03body{"maxRDiffV03body", 3, "Max RDiff of V0 to 3body"};
    Configurable<float> minPt3Body{"minPt3Body", 0.5, "Min Pt of 3body"};
    Configurable<float> maxTgl3Body{"maxTgl3Body", 0.01, "Max tgl of 3body"};
    Configurable<float> maxDCAXY3Body{"maxDCAXY3Body", 0.5, "Max DCA XY of 3body"};
    Configurable<float> maxDCAZ3Body{"maxDCAZ3Body", 1.0, "Max DCA Z of 3body"};
  } mixingOpts;

  // Helper struct to contain MC information prior to filling
  struct mc3Bodyinfo {
    int label;
    std::array<float, 3> genDecVtx{0.f};
    std::array<float, 3> genMomentum{0.f};
    float genCt;
    float genPhi;
    float genEta;
    float genRapidity;
    float genMomProton;
    float genMomPion;
    float genMomDeuteron;
    float genPtProton;
    float genPtPion;
    float genPtDeuteron;
    bool isTrueH3L;
    bool isTrueAntiH3L;
    bool isReco;
    int motherPdgCode;
    int daughterPrPdgCode;
    int daughterPiPdgCode;
    int daughterDePdgCode;
    bool isDeuteronPrimary;
    bool survivedEventSel;
  };
  mc3Bodyinfo this3BodyMCInfo;

  // CCDB and magnetic field
  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  std::unordered_map<int, float> ccdbCache; // Maps runNumber -> d_bz
  o2::base::MatLayerCylSet* lut = nullptr;

  // histogram registry
  HistogramRegistry registry{"Registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  // bachelor TOF PID
  o2::aod::pidtofgeneric::TofPidNewCollision<TracksExtPIDIUwithEvTimes::iterator> bachelorTOFPID;               // to be updated in Init based on the hypothesis
  o2::aod::pidtofgeneric::TofPidNewCollision<TracksExtPIDIUwithEvTimesLabeled::iterator> bachelorTOFPIDLabeled; // to be updated in Init based on the hypothesis
  // TOF response and input parameters
  o2::pid::tof::TOFResoParamsV3 mRespParamsV3;
  o2::aod::pidtofgeneric::TOFCalibConfig mTOFCalibConfig; // TOF Calib configuration

  // 3body mixing
  using Binning3BodyKF = ColumnBinningPolicy<aod::reduceddecay3body::RadiusKF, aod::reduceddecay3body::PhiKF>;
  using Binning3BodyDCAfitter = ColumnBinningPolicy<aod::reduceddecay3body::RadiusDCA, aod::reduceddecay3body::PhiDCA>;

  // skimmed processing
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // tracked cluster size
  std::vector<int> fTrackedClSizeVector;

  // trigger info
  std::vector<bool> isTriggeredCollision;
  // MC info
  std::vector<bool> isGoodCollision;

  void init(InitContext& initContext)
  {
    zorroSummary.setObject(zorro.getZorroSummary());

    mRunNumber = 0;
    d_bz = 0;

    mEnabledTables.resize(nTables, 0);

    // CCDB options
    ccdb->setURL(ccdbConfigurations.ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // TOF PID parameters initialization
    if (doprocessRealData == true || doprocessMonteCarlo == true) {
      mTOFCalibConfig.metadataInfo = metadataInfo;
      mTOFCalibConfig.inheritFromBaseTask(initContext);
      mTOFCalibConfig.initSetup(mRespParamsV3, ccdb);
    }

    // Set material correction
    if (useMatCorrType == 1) {
      LOGF(info, "TGeo correction requested, loading geometry");
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(ccdbConfigurations.geoPath);
      }
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    }
    if (useMatCorrType == 2) {
      LOGF(info, "LUT correction requested, loading LUT");
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbConfigurations.lutPath));
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    }

    helper.fitterV0.setMatCorrType(matCorr);
    helper.fitter3body.setMatCorrType(matCorr);

    // set bachelor PID
    bachelorTOFPID.SetPidType(o2::track::PID::Deuteron);
    bachelorTOFPIDLabeled.SetPidType(o2::track::PID::Deuteron);

    // set decay3body parameters in the helper
    helper.decay3bodyselections.maxEtaDaughters = decay3bodyBuilderOpts.maxEtaDaughters;
    helper.decay3bodyselections.minTPCNClProton = decay3bodyBuilderOpts.minTPCNClProton;
    helper.decay3bodyselections.minTPCNClPion = decay3bodyBuilderOpts.minTPCNClPion;
    helper.decay3bodyselections.minTPCNClDeuteron = decay3bodyBuilderOpts.minTPCNClDeuteron;
    helper.decay3bodyselections.minDCAProtonToPV = decay3bodyBuilderOpts.minDCAProtonToPV;
    helper.decay3bodyselections.minDCAPionToPV = decay3bodyBuilderOpts.minDCAPionToPV;
    helper.decay3bodyselections.minDCADeuteronToPV = decay3bodyBuilderOpts.minDCADeuteronToPV;
    helper.decay3bodyselections.minPtProton = decay3bodyBuilderOpts.minPtProton;
    helper.decay3bodyselections.minPtPion = decay3bodyBuilderOpts.minPtPion;
    helper.decay3bodyselections.minPtDeuteron = decay3bodyBuilderOpts.minPtDeuteron;
    helper.decay3bodyselections.maxPtProton = decay3bodyBuilderOpts.maxPtProton;
    helper.decay3bodyselections.maxPtPion = decay3bodyBuilderOpts.maxPtPion;
    helper.decay3bodyselections.maxPtDeuteron = decay3bodyBuilderOpts.maxPtDeuteron;
    helper.decay3bodyselections.maxTPCnSigma = decay3bodyBuilderOpts.maxTPCnSigma;
    helper.decay3bodyselections.minTOFnSigmaDeuteron = decay3bodyBuilderOpts.minTOFnSigmaDeuteron;
    helper.decay3bodyselections.maxTOFnSigmaDeuteron = decay3bodyBuilderOpts.maxTOFnSigmaDeuteron;
    helper.decay3bodyselections.minPDeuteronUseTOF = decay3bodyBuilderOpts.minPDeuteronUseTOF;
    helper.decay3bodyselections.maxDCADauToSVaverage = decay3bodyBuilderOpts.maxDCADauToSVaverage;
    helper.decay3bodyselections.maxRapidity = decay3bodyBuilderOpts.maxRapidity;
    helper.decay3bodyselections.minPt = decay3bodyBuilderOpts.minPt;
    helper.decay3bodyselections.maxPt = decay3bodyBuilderOpts.maxPt;
    helper.decay3bodyselections.minMass = decay3bodyBuilderOpts.minMass;
    helper.decay3bodyselections.maxMass = decay3bodyBuilderOpts.maxMass;
    helper.decay3bodyselections.minCtau = decay3bodyBuilderOpts.minCtau;
    helper.decay3bodyselections.maxCtau = decay3bodyBuilderOpts.maxCtau;
    helper.decay3bodyselections.minCosPA = decay3bodyBuilderOpts.minCosPA;
    helper.decay3bodyselections.maxChi2 = decay3bodyBuilderOpts.maxChi2;

    // set SVertexer selection parameters in the helper
    helper.svertexerselections.minPt2V0 = mixingOpts.minPt2V0;
    helper.svertexerselections.maxTgl2V0 = mixingOpts.maxTgl2V0;
    helper.svertexerselections.maxDCAXY2ToMeanVertex3bodyV0 = mixingOpts.maxDCAXY2ToMeanVertex3bodyV0;
    helper.svertexerselections.minCosPAXYMeanVertex3bodyV0 = mixingOpts.minCosPAXYMeanVertex3bodyV0;
    helper.svertexerselections.minCosPA3bodyV0 = mixingOpts.minCosPA3bodyV0;
    helper.svertexerselections.maxRDiffV03body = mixingOpts.maxRDiffV03body;
    helper.svertexerselections.minPt3Body = mixingOpts.minPt3Body;
    helper.svertexerselections.maxTgl3Body = mixingOpts.maxTgl3Body;
    helper.svertexerselections.maxDCAXY3Body = mixingOpts.maxDCAXY3Body;
    helper.svertexerselections.maxDCAZ3Body = mixingOpts.maxDCAZ3Body;

    // list enabled process functions
    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
    LOGF(info, " Decay3body builder: basic configuration listing");
    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");

    if (doprocessRealData) {
      LOGF(info, " ===> process function enabled: processRealData");
    }
    if (doprocessRealDataReduced) {
      LOGF(info, " ===> process function enabled: processRealDataReduced");
    }
    if (doprocessRealDataReduced3bodyMixing) {
      LOGF(info, " ===> process function enabled: processRealDataReduced3bodyMixing");
    }
    if (doprocessMonteCarlo) {
      LOGF(info, " ===> process function enabled: processMonteCarlo");
    }

    // list enabled tables
    for (int i = 0; i < nTables; i++) {
      if (mEnabledTables[i]) {
        LOGF(info, " -~> Table enabled: %s", tableNames[i]);
      }
    }

    // print base cuts
    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
    LOGF(info, "-~> max daughter eta ..............: %f", decay3bodyBuilderOpts.maxEtaDaughters.value);
    LOGF(info, "-~> min TPC ncls proton ...........: %i", decay3bodyBuilderOpts.minTPCNClProton.value);
    LOGF(info, "-~> min TPC ncls pion .............: %i", decay3bodyBuilderOpts.minTPCNClPion.value);
    LOGF(info, "-~> min TPC ncls bach .............: %i", decay3bodyBuilderOpts.minTPCNClDeuteron.value);
    LOGF(info, "-~> min DCA proton to PV ..........: %f", decay3bodyBuilderOpts.minDCAProtonToPV.value);
    LOGF(info, "-~> min DCA pion to PV ............: %f", decay3bodyBuilderOpts.minDCAPionToPV.value);
    LOGF(info, "-~> min DCA bach to PV ............: %f", decay3bodyBuilderOpts.minDCADeuteronToPV.value);
    LOGF(info, "-~> min pT proton .................: %f", decay3bodyBuilderOpts.minPtProton.value);
    LOGF(info, "-~> min pT pion ...................: %f", decay3bodyBuilderOpts.minPtPion.value);
    LOGF(info, "-~> min pT bach ...................: %f", decay3bodyBuilderOpts.minPtDeuteron.value);
    LOGF(info, "-~> max pT proton .................: %f", decay3bodyBuilderOpts.maxPtProton.value);
    LOGF(info, "-~> max pT pion ...................: %f", decay3bodyBuilderOpts.maxPtPion.value);
    LOGF(info, "-~> max pT bach ...................: %f", decay3bodyBuilderOpts.maxPtDeuteron.value);
    LOGF(info, "-~> max TPC nSigma ...............: %f", decay3bodyBuilderOpts.maxTPCnSigma.value);
    LOGF(info, "-~> min TOF nSigma deuteron ......: %f", decay3bodyBuilderOpts.minTOFnSigmaDeuteron.value);
    LOGF(info, "-~> max TOF nSigma deuteron ......: %f", decay3bodyBuilderOpts.maxTOFnSigmaDeuteron.value);
    LOGF(info, "-~> min p bach use TOF ...........: %f", decay3bodyBuilderOpts.minPDeuteronUseTOF.value);
    LOGF(info, "-~> max DCA dau at SV ............: %f", decay3bodyBuilderOpts.maxDCADauToSVaverage.value);
    LOGF(info, "-~> max rapidity .................: %f", decay3bodyBuilderOpts.maxRapidity.value);
    LOGF(info, "-~> min pT .......................: %f", decay3bodyBuilderOpts.minPt.value);
    LOGF(info, "-~> max pT .......................: %f", decay3bodyBuilderOpts.maxPt.value);
    LOGF(info, "-~> min mass .....................: %f", decay3bodyBuilderOpts.minMass.value);
    LOGF(info, "-~> max mass .....................: %f", decay3bodyBuilderOpts.maxMass.value);
    LOGF(info, "-~> min ctau .....................: %f", decay3bodyBuilderOpts.minCtau.value);
    LOGF(info, "-~> max ctau .....................: %f", decay3bodyBuilderOpts.maxCtau.value);
    LOGF(info, "-~> min cosPA ....................: %f", decay3bodyBuilderOpts.minCosPA.value);
    LOGF(info, "-~> max chi2 .....................: %f", decay3bodyBuilderOpts.maxChi2.value);
    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");

    // bookkeeping histograms
    auto h = registry.add<TH1>("Counters/hTableBuildingStatistics", "hTableBuildingStatistics", kTH1D, {{nTablesConst, -0.5f, static_cast<float>(nTablesConst)}});
    auto h2 = registry.add<TH1>("Counters/hInputStatistics", "hInputStatistics", kTH1D, {{nTablesConst, -0.5f, static_cast<float>(nTablesConst)}});
    h2->SetTitle("Input table sizes");

    // configure tables to generate
    for (int i = 0; i < nTables; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, tableNames[i].c_str());
      h2->GetXaxis()->SetBinLabel(i + 1, tableNames[i].c_str());
      h->SetBinContent(i + 1, 0); // mark all as disabled to start

      int f = enabledTables->get(tableNames[i].c_str(), "enable");
      if (f == 1) {
        mEnabledTables[i] = 1;
        h->SetBinContent(i + 1, 1); // mark enabled
      }
    }

    if (mEnabledTables[kVtx3BodyDatas] && mEnabledTables[kMcVtx3BodyDatas]) {
      LOG(fatal) << "Tables Vtx3BodyDatas and McVtx3BodyDatas cannot both be enabled at the same time. Choose one!";
    }

    // Add histograms separately for different process functions
    if (doprocessRealData == true || doprocessMonteCarlo == true) {
      auto hEventCounter = registry.add<TH1>("Counters/hEventCounter", "hEventCounter", HistType::kTH1D, {{3, 0.0f, 3.0f}});
      hEventCounter->GetXaxis()->SetBinLabel(1, "total");
      hEventCounter->GetXaxis()->SetBinLabel(2, "sel8");
      hEventCounter->GetXaxis()->SetBinLabel(3, "vertexZ");
      hEventCounter->LabelsOption("v");
    }

    if (doprocessRealData == true || doprocessRealDataReduced == true || doprocessMonteCarlo == true) {
      if (doTrackQA) { // histograms for all daughter tracks of (selected) 3body candidates
        registry.add("QA/Tracks/hTrackProtonTPCNcls", "hTrackProtonTPCNcls", HistType::kTH1F, {{152, 0, 152, "# TPC clusters"}});
        registry.add("QA/Tracks/hTrackPionTPCNcls", "hTrackPionTPCNcls", HistType::kTH1F, {{152, 0, 152, "# TPC clusters"}});
        registry.add("QA/Tracks/hTrackDeuteronTPCNcls", "hTrackDeuteronTPCNcls", HistType::kTH1F, {{152, 0, 152, "# TPC clusters"}});
        registry.add("QA/Tracks/hTrackProtonHasTPC", "hTrackProtonHasTPC", HistType::kTH1F, {{2, -0.5, 1.5, "has TPC"}});
        registry.add("QA/Tracks/hTrackPionHasTPC", "hTrackPionHasTPC", HistType::kTH1F, {{2, -0.5, 1.5, "has TPC"}});
        registry.add("QA/Tracks/hTrackDeuteronHasTPC", "hTrackDeuteronHasTPC", HistType::kTH1F, {{2, -0.5, 1.5, "has TPC"}});
        registry.add("QA/Tracks/hTrackDeuteronITSClusSizes", "hTrackDeuteronITSClusSizes", HistType::kTH1F, {{10, 0., 10., "ITS cluster sizes"}});
        registry.add("QA/Tracks/hTrackProtonTPCPID", "hTrackProtonTPCPID", HistType::kTH2F, {{100, -10.0f, 10.0f, "p/z (GeV/c)"}, {100, -10.0f, 10.0f, "TPC n#sigma"}});
        registry.add("QA/Tracks/hTrackPionTPCPID", "hTrackPionTPCPID", HistType::kTH2F, {{100, -10.0f, 10.0f, "p/z (GeV/c)"}, {100, -10.0f, 10.0f, "TPC n#sigma"}});
        registry.add("QA/Tracks/hTrackDeuteronTPCPID", "hTrackDeuteronTPCPID", HistType::kTH2F, {{100, -10.0f, 10.0f, "p/z (GeV/c)"}, {100, -10.0f, 10.0f, "TPC n#sigma"}});
        registry.add("QA/Tracks/hTrackProtonPt", "hTrackProtonPt", HistType::kTH1F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}});
        registry.add("QA/Tracks/hTrackPionPt", "hTrackPionPt", HistType::kTH1F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}});
        registry.add("QA/Tracks/hTrackDeuteronPt", "hTrackDeuteronPt", HistType::kTH1F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}});
      }
      if (doVertexQA) {
        registry.add("QA/Event/hAllSelEventsVtxZ", "hAllSelEventsVtxZ", HistType::kTH1F, {{500, -15.0f, 15.0f, "PV Z (cm)"}});
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

    if (doprocessRealDataReduced3bodyMixing == true) {
      auto h3bodyCombinationCounter = registry.add<TH1>("Mixing/h3bodyCombinationCounter", "h3bodyCombinationCounter", HistType::kTH1D, {{4, 0.0f, 4.0f}});
      h3bodyCombinationCounter->GetXaxis()->SetBinLabel(1, "total");
      h3bodyCombinationCounter->GetXaxis()->SetBinLabel(2, "not same collision");
      h3bodyCombinationCounter->GetXaxis()->SetBinLabel(3, "collision VtxZ");
      h3bodyCombinationCounter->GetXaxis()->SetBinLabel(4, "bach sign/ID");
      h3bodyCombinationCounter->LabelsOption("v");
      registry.add("Mixing/hDecay3BodyRadiusPhi", "hDecay3BodyRadiusPhi", HistType::kTH2F, {mixingOpts.bins3BodyRadius, mixingOpts.bins3BodyPhi});
      registry.add("Mixing/hDecay3BodyPosZ", "hDecay3BodyPosZ", HistType::kTH1F, {mixingOpts.bins3BodyPosZ});
    }
  }

  template <typename TCollisions>
  bool initCCDB(aod::BCsWithTimestamps const& bcs, TCollisions const& collisions)
  {
    auto bc = collisions.size() ? collisions.begin().template bc_as<aod::BCsWithTimestamps>() : bcs.begin();
    if (!bcs.size()) {
      LOGF(warn, "No BC found, skipping this DF.");
      return false; // signal to skip this DF
    }

    if (mRunNumber == bc.runNumber()) {
      return true;
    }

    if (doSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerList);
      zorro.populateHistRegistry(registry, bc.runNumber());
    }

    auto timestamp = bc.timestamp();
    o2::parameters::GRPMagField* grpmag = 0x0;
    ccdb->clearCache(ccdbConfigurations.grpmagPath);
    grpmag = ccdb->getSpecific<o2::parameters::GRPMagField>(ccdbConfigurations.grpmagPath, timestamp);
    if (!grpmag) {
      LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfigurations.grpmagPath << " of object GRPMagField for timestamp " << timestamp;
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);
    // Fetch magnetic field from ccdb for current collision
    auto d_bz = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp << " with magnetic field of " << d_bz << " kG";

    // set magnetic field value for DCA fitter
    helper.fitterV0.setBz(d_bz);
    helper.fitter3body.setBz(d_bz);
// Set magnetic field for KF vertexing
#ifdef HomogeneousField
    KFParticle::SetField(d_bz);
#endif

    if (useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized
      // (setMatLUT has implicit and problematic init field call if not)
      LOG(info) << "Loading material look-up table for timestamp: " << timestamp;
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(ccdbConfigurations.lutPath, timestamp));
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }

    // mark run as configured
    mRunNumber = bc.runNumber();

    mTOFCalibConfig.processSetup(mRespParamsV3, ccdb, bc);

    return true;
  }

  float getMagFieldFromRunNumber(int runNumber)
  {
    float magField;
    // Check if the CCDB data for this run is already cached
    if (ccdbCache.find(runNumber) != ccdbCache.end()) {
      LOG(debug) << "CCDB data already cached for run " << runNumber;
      magField = ccdbCache[runNumber];
      // if not, retrieve it from CCDB
    } else {
      std::shared_ptr<o2::parameters::GRPMagField> grpmag = std::make_shared<o2::parameters::GRPMagField>(*ccdb->getForRun<o2::parameters::GRPMagField>(ccdbConfigurations.grpmagPath, runNumber));
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfigurations.grpmagPath << " of object GRPMagField and " << ccdbConfigurations.grpPath << " of object GRPObject for run number " << runNumber;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag.get());
      // Fetch magnetic field from ccdb for current collision
      magField = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << "Retrieved GRP for run number " << runNumber << " with magnetic field of " << d_bz << " kZG";

      // cache magnetic field info
      ccdbCache[runNumber] = magField;
    }
    return magField;
  }

  void initFittersWithMagField(int runNumber, float magField)
  {
    // set magnetic field only when run number changes
    if (mRunNumber == runNumber) {
      LOG(debug) << "CCDB initialized for run " << mRunNumber;
      return;
    }
    mRunNumber = runNumber; // Update the last run number

    // update propagator
    o2::base::Propagator::Instance()->setNominalBz(magField);

    // Set magnetic field for KF vertexing
#ifdef HomogeneousField
    KFParticle::SetField(magField);
#endif
    // Set field for DCAfitter
    helper.fitterV0.setBz(magField);
    helper.fitter3body.setBz(magField);

    if (useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
  }

  // ______________________________________________________________
  // function to build decay3body candidates
  template <class TTracksTo, typename TBCs, typename TCollisions, typename T3Bodys, typename TMCParticles, typename TMCCollisions>
  void buildCandidates(TBCs const&,
                       TCollisions const& collisions,
                       T3Bodys const& decay3bodys,
                       TMCParticles const& mcParticles,
                       TMCCollisions const& mcCollisions)
  {
    if (!(mEnabledTables[kVtx3BodyDatas] || mEnabledTables[kMcVtx3BodyDatas])) {
      LOG(info) << "No request for candidate analysis table in place, skipping candidate building." << std::endl;
      return; // don't do if no request for decay3bodys in place
    }

    // prepare MC container (not necessarily used)
    std::vector<bool> mcParticleIsReco;

    if constexpr (soa::is_table<TBCs>) {
      isTriggeredCollision.clear();
      isTriggeredCollision.resize(collisions.size(), false);
    }
    // clear and reserve size for MC info vectors
    if constexpr (soa::is_table<TMCParticles>) {
      isGoodCollision.clear();
      mcParticleIsReco.clear();
      isGoodCollision.resize(mcCollisions.size(), false);
      mcParticleIsReco.resize(mcParticles.size(), false);
    }

    // Loop over collisions for vertex QA
    for (const auto& collision : collisions) {
      if constexpr (soa::is_table<TBCs>) { // only do if NOT running over reduced data (already done in reducedCreator)
        // Zorro event counting
        bool isZorroSelected = false;
        if (doSkimmedProcessing) {
          isZorroSelected = zorro.isSelected(collision.template bc_as<TBCs>().globalBC());
          if (!isZorroSelected && onlyKeepInterestedTrigger) {
            continue;
          }
        }

        isTriggeredCollision[collision.globalIndex()] = true;
        // event counting
        registry.fill(HIST("Counters/hEventCounter"), 0.5);
        if (doSel8selection && !collision.sel8()) {
          continue;
        }
        registry.fill(HIST("Counters/hEventCounter"), 1.5);
        if (doPosZselection && (collision.posZ() >= 10.0f || collision.posZ() <= -10.0f)) {
          continue;
        }
        registry.fill(HIST("Counters/hEventCounter"), 2.5);
      }

      // vertex QA and counting
      if (doVertexQA) {
        registry.fill(HIST("QA/Event/hAllSelEventsVtxZ"), collision.posZ());
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

      // In case of MC: reco collision survived event selection filter --> fill value for MC collision if collision is "true" MC collision
      if constexpr (soa::is_table<TMCParticles>) {
        if (collision.mcCollisionId() >= 0) {
          isGoodCollision[collision.mcCollisionId()] = true;
        }
      }
    } // loop over collisions

    // Loop over all decay3bodys in same time frame
    registry.fill(HIST("Counters/hInputStatistics"), kVtx3BodyDatas, decay3bodys.size());
    int lastRunNumber = -1;
    for (const auto& decay3body : decay3bodys) {
      // only build tracked decay3body if aksed
      if (decay3bodyBuilderOpts.buildOnlyTracked && fTrackedClSizeVector[decay3body.globalIndex()] == 0) {
        continue;
      }

      // skip decay3body without assigned collision
      /// TODO: do we want this??
      if (decay3body.collisionId() < 0) {
        continue;
      }

      // aquire collision
      auto const& collision = collisions.rawIteratorAt(decay3body.collisionId());

      // initialise CCDB from run number saved in reduced collisions table when running over reduced data
      if constexpr (!soa::is_table<TBCs>) { // only do if running over reduced data (otherwise CCDB is initialised in process function)
        if (collision.runNumber() != lastRunNumber) {
          initFittersWithMagField(collision.runNumber(), getMagFieldFromRunNumber(collision.runNumber()));
          lastRunNumber = collision.runNumber(); // Update the last run number
          LOG(debug) << "CCDB initialized for run " << lastRunNumber;
        }
      }

      // event selection
      if constexpr (soa::is_table<TBCs>) { // only when NOT running over reduced data
        if (doSel8selection && !collision.sel8()) {
          continue;
        }
        if (onlyKeepInterestedTrigger && !isTriggeredCollision[collision.globalIndex()]) {
          continue;
        }
      }
      if (doPosZselection && (collision.posZ() >= 10.0f || collision.posZ() <= -10.0f)) {
        continue;
      }

      // aquire tracks
      auto trackPos = decay3body.template track0_as<TTracksTo>();
      auto trackNeg = decay3body.template track1_as<TTracksTo>();
      auto trackDeuteron = decay3body.template track2_as<TTracksTo>();
      int protonSign = doLikeSign ? -trackDeuteron.sign() : trackDeuteron.sign();
      auto trackProton = protonSign > 0 ? trackPos : trackNeg;
      auto trackPion = protonSign > 0 ? trackNeg : trackPos;

      // get deuteron TOF PID
      float tofNSigmaDeuteron;
      if constexpr (!soa::is_table<TBCs>) { // running over derived data
        tofNSigmaDeuteron = trackDeuteron.tofNSigmaDe();
      } else if constexpr (soa::is_table<TBCs>) {    // running over AO2Ds
        if constexpr (soa::is_table<TMCParticles>) { // running over MC (track table with labels)
          tofNSigmaDeuteron = getTOFnSigma<true /*isMC*/, TCollisions>(mRespParamsV3, collision, trackDeuteron);
        } else { // running over real data
          tofNSigmaDeuteron = getTOFnSigma<false /*isMC*/, TCollisions>(mRespParamsV3, collision, trackDeuteron);
        }
      }

      /// build Decay3body candidate
      if (!helper.buildDecay3BodyCandidate(collision,
                                           trackProton,
                                           trackPion,
                                           trackDeuteron,
                                           decay3body.globalIndex(),
                                           tofNSigmaDeuteron,
                                           fTrackedClSizeVector[decay3body.globalIndex()],
                                           decay3bodyBuilderOpts.useKFParticle,
                                           decay3bodyBuilderOpts.kfSetTopologicalConstraint,
                                           decay3bodyBuilderOpts.useSelections,
                                           decay3bodyBuilderOpts.useChi2Selection,
                                           decay3bodyBuilderOpts.useTPCforPion,
                                           decay3bodyBuilderOpts.acceptTPCOnly,
                                           decay3bodyBuilderOpts.askOnlyITSMatch,
                                           decay3bodyBuilderOpts.calculateCovariance,
                                           false /*isEventMixing*/,
                                           false /*applySVertexerCuts*/)) {
        continue;
      }

      // fill QA histograms
      if (doTrackQA) { // histograms filled for daughter tracks of (selected) 3body candidates
        registry.fill(HIST("QA/Tracks/hTrackProtonTPCNcls"), trackProton.tpcNClsFound());
        registry.fill(HIST("QA/Tracks/hTrackPionTPCNcls"), trackPion.tpcNClsFound());
        registry.fill(HIST("QA/Tracks/hTrackDeuteronTPCNcls"), trackDeuteron.tpcNClsFound());
        registry.fill(HIST("QA/Tracks/hTrackProtonHasTPC"), trackProton.hasTPC());
        registry.fill(HIST("QA/Tracks/hTrackPionHasTPC"), trackPion.hasTPC());
        registry.fill(HIST("QA/Tracks/hTrackDeuteronHasTPC"), trackDeuteron.hasTPC());
        registry.fill(HIST("QA/Tracks/hTrackDeuteronITSClusSizes"), trackDeuteron.itsClusterSizes());
        registry.fill(HIST("QA/Tracks/hTrackProtonTPCPID"), trackProton.sign() * trackProton.tpcInnerParam(), trackProton.tpcNSigmaPr());
        registry.fill(HIST("QA/Tracks/hTrackPionTPCPID"), trackPion.sign() * trackPion.tpcInnerParam(), trackPion.tpcNSigmaPi());
        registry.fill(HIST("QA/Tracks/hTrackDeuteronTPCPID"), trackDeuteron.sign() * trackDeuteron.tpcInnerParam(), trackDeuteron.tpcNSigmaDe());
        registry.fill(HIST("QA/Tracks/hTrackProtonPt"), trackProton.pt());
        registry.fill(HIST("QA/Tracks/hTrackPionPt"), trackPion.pt());
        registry.fill(HIST("QA/Tracks/hTrackDeuteronPt"), trackDeuteron.pt());
      }

      // generate analysis tables with current candidate (only Vtx3BodyDatas is filled here, McVtx3BodyDatas table is filled later)
      if (!mEnabledTables[kMcVtx3BodyDatas]) {
        fillAnalysisTables();
      }

      // ___________________________________________________________
      // MC handling part: matching of reconstructed candidates
      // ___________________________________________________________
      // fill MC table with reco MC candidate information and gen information if matched to MC particle
      if constexpr (soa::is_table<TMCParticles>) {
        // MC info
        resetMCInfo(this3BodyMCInfo);
        this3BodyMCInfo.isReco = true;

        // set flag if selected reco collision has matched gen collision
        if (collision.mcCollisionId() >= 0) { // reco collision is matched to gen collision
          this3BodyMCInfo.survivedEventSel = isGoodCollision[collision.mcCollisionId()];
        } else {
          this3BodyMCInfo.survivedEventSel = false; // false if reco collision not matched to gen collision
        }

        // check if daughters have MC particle
        if (!trackProton.has_mcParticle() || !trackPion.has_mcParticle() || !trackDeuteron.has_mcParticle()) {
          continue;
        }

        // get MC daughter particles
        auto mcTrackProton = trackProton.template mcParticle_as<aod::McParticles>();
        auto mcTrackPion = trackPion.template mcParticle_as<aod::McParticles>();
        auto mcTrackDeuteron = trackDeuteron.template mcParticle_as<aod::McParticles>();

        // set daughter MC info (also for non-matched candidates)
        this3BodyMCInfo.daughterPrPdgCode = mcTrackProton.pdgCode();
        this3BodyMCInfo.daughterPiPdgCode = mcTrackPion.pdgCode();
        this3BodyMCInfo.daughterDePdgCode = mcTrackDeuteron.pdgCode();
        this3BodyMCInfo.isDeuteronPrimary = mcTrackDeuteron.isPhysicalPrimary();
        this3BodyMCInfo.genMomProton = mcTrackProton.p();
        this3BodyMCInfo.genMomPion = mcTrackPion.p();
        this3BodyMCInfo.genMomDeuteron = mcTrackDeuteron.p();
        this3BodyMCInfo.genPtProton = mcTrackProton.pt();
        this3BodyMCInfo.genPtPion = mcTrackPion.pt();
        this3BodyMCInfo.genPtDeuteron = mcTrackDeuteron.pt();

        // check if reco mother is true H3L/Anti-H3l
        bool isMuonReco;
        int motherID = checkH3LTruth(mcTrackProton, mcTrackPion, mcTrackDeuteron, isMuonReco);

        // get generated mother MC info
        if (motherID > 0) {
          auto mcTrackH3L = mcParticles.rawIteratorAt(motherID);
          this3BodyMCInfo.motherPdgCode = mcTrackH3L.pdgCode();
          this3BodyMCInfo.label = motherID;
          this3BodyMCInfo.genMomentum = {mcTrackH3L.px(), mcTrackH3L.py(), mcTrackH3L.pz()};
          this3BodyMCInfo.genDecVtx = {mcTrackProton.vx(), mcTrackProton.vy(), mcTrackProton.vz()};
          this3BodyMCInfo.genCt = RecoDecay::sqrtSumOfSquares(mcTrackProton.vx() - mcTrackH3L.vx(), mcTrackProton.vy() - mcTrackH3L.vy(), mcTrackProton.vz() - mcTrackH3L.vz()) * o2::constants::physics::MassHyperTriton / mcTrackH3L.p();
          this3BodyMCInfo.genPhi = mcTrackH3L.phi();
          this3BodyMCInfo.genEta = mcTrackH3L.eta();
          this3BodyMCInfo.genRapidity = mcTrackH3L.y();
          this3BodyMCInfo.isTrueH3L = this3BodyMCInfo.motherPdgCode > 0 ? true : false;
          this3BodyMCInfo.isTrueAntiH3L = this3BodyMCInfo.motherPdgCode < 0 ? true : false;
        }

        // fill analysis tables (only McVtx3BodyDatas is filled here)
        fillAnalysisTables();

        // mark mcParticle as reconstructed
        if (this3BodyMCInfo.label > -1) {
          mcParticleIsReco[this3BodyMCInfo.label] = true;
        }
      } // constexpr requires mcParticles check
    } // decay3body loop

    // ____________________________________________________________________
    // MC handling part: generated information of non-reco candidates
    // ____________________________________________________________________
    if constexpr (soa::is_table<TMCParticles>) {
      for (const auto& mcparticle : mcParticles) {
        // MC info
        resetMCInfo(this3BodyMCInfo);

        // skip MC particle if reconstructed and already filled previously
        if (mcParticleIsReco[mcparticle.globalIndex()] == true) {
          continue;
        }
        this3BodyMCInfo.isReco = false;

        // set flag if corresponding MC collision has matched reconstructed collision which passed event selection
        this3BodyMCInfo.survivedEventSel = isGoodCollision[mcparticle.mcCollisionId()];

        // check if MC particle is hypertriton
        if (std::abs(mcparticle.pdgCode()) != o2::constants::physics::Pdg::kHyperTriton) {
          continue;
        }

        // check daughter identities
        bool haveProton = false, havePion = false, haveDeuteron = false;
        bool haveAntiProton = false, haveAntiPion = false, haveAntiDeuteron = false;
        for (const auto& mcparticleDaughter : mcparticle.template daughters_as<TMCParticles>()) {
          if (mcparticleDaughter.pdgCode() == PDG_t::kProton)
            haveProton = true;
          if (mcparticleDaughter.pdgCode() == PDG_t::kProtonBar)
            haveAntiProton = true;
          if (mcparticleDaughter.pdgCode() == PDG_t::kPiPlus)
            havePion = true;
          if (mcparticleDaughter.pdgCode() == PDG_t::kPiMinus)
            haveAntiPion = true;
          if (mcparticleDaughter.pdgCode() == o2::constants::physics::Pdg::kDeuteron)
            haveDeuteron = true;
          if (mcparticleDaughter.pdgCode() == -o2::constants::physics::Pdg::kDeuteron)
            haveAntiDeuteron = true;
        }

        // check if hypertriton decayed via 3-body decay and is particle or anti-particle
        if ((haveProton && haveAntiPion && haveDeuteron && !(haveAntiProton || havePion || haveAntiDeuteron)) || (haveAntiProton && havePion && haveAntiDeuteron && !(haveProton || haveAntiPion || haveDeuteron))) {
          if (mcparticle.pdgCode() > 0) {
            this3BodyMCInfo.isTrueH3L = true;
          } else if (mcparticle.pdgCode() < 0) {
            this3BodyMCInfo.isTrueAntiH3L = true;
          }
          // get daughters
          for (const auto& mcparticleDaughter : mcparticle.template daughters_as<aod::McParticles>()) {
            if (std::abs(mcparticleDaughter.pdgCode()) == PDG_t::kProton) { // proton
              this3BodyMCInfo.genMomProton = mcparticleDaughter.p();
              this3BodyMCInfo.genPtProton = mcparticleDaughter.pt();
              this3BodyMCInfo.daughterPrPdgCode = mcparticleDaughter.pdgCode();
              this3BodyMCInfo.genDecVtx = {mcparticleDaughter.vx(), mcparticleDaughter.vy(), mcparticleDaughter.vz()};
            } else if (std::abs(mcparticleDaughter.pdgCode()) == PDG_t::kPiPlus) { // pion
              this3BodyMCInfo.genMomPion = mcparticleDaughter.p();
              this3BodyMCInfo.genPtPion = mcparticleDaughter.pt();
              this3BodyMCInfo.daughterPiPdgCode = mcparticleDaughter.pdgCode();
            } else if (std::abs(mcparticleDaughter.pdgCode()) == o2::constants::physics::Pdg::kDeuteron) { // deuteron
              this3BodyMCInfo.genMomDeuteron = mcparticleDaughter.p();
              this3BodyMCInfo.genPtDeuteron = mcparticleDaughter.pt();
              this3BodyMCInfo.daughterDePdgCode = mcparticleDaughter.pdgCode();
              this3BodyMCInfo.isDeuteronPrimary = mcparticleDaughter.isPhysicalPrimary();
            }
          }
        } else {
          continue; // stop if particle is not decayed via 3-body decay
        }

        // calculate ctau
        this3BodyMCInfo.genCt = RecoDecay::sqrtSumOfSquares(this3BodyMCInfo.genDecVtx[0] - mcparticle.vx(), this3BodyMCInfo.genDecVtx[1] - mcparticle.vy(), this3BodyMCInfo.genDecVtx[2] - mcparticle.vz()) * o2::constants::physics::MassHyperTriton / mcparticle.p();

        // fill MCDecay3BodyCores table if requested
        if (mEnabledTables[kMcVtx3BodyDatas]) {
          products.mcvtx3bodydatas(-1,                 // sign
                                   -1., -1.,           // mass, massV0
                                   -1., -1., -1.,      // position
                                   -1., -1., -1.,      // momentum
                                   -1.,                // chi2
                                   -1.,                // trackedClSize
                                   -1., -1., -1.,      // momProton
                                   -1., -1., -1.,      // momPion
                                   -1., -1., -1.,      // momDeuteron
                                   -1., -1., -1.,      // trackDCAxyToPV: 0 - proton, 1 - pion, 2 - deuteron
                                   -1., -1., -1.,      // trackDCAToPV: 0 - proton, 1 - pion, 2 - deuteron
                                   -1., -1., -1.,      // trackDCAxyToPVprop: 0 - proton, 1 - pion, 2 - deuteron
                                   -1., -1., -1.,      // trackDCAToPVprop: 0 - proton, 1 - pion, 2 - deuteron
                                   -1., -1., -1.,      // daughterDCAtoSV: 0 - proton, 1 - pion, 2 - deuteron
                                   -1.,                // daughterDCAtoSVaverage
                                   -1., -1.,           // cosPA, ctau
                                   -1., -1., -1., -1., // tpcNsigma: 0 - proton, 1 - pion, 2 - deuteron, 3 - bach with pion hyp
                                   -1.,                // tofNsigmaDeuteron
                                   -1., -1., -1.,      // average ITS cluster sizes: proton, pion, deuteron
                                   -1., -1., -1.,      // TPCNCl: proton, pion, deuteron
                                   -1.,                // pidForTrackingDeuteron
                                   // MC information
                                   mcparticle.px(), mcparticle.py(), mcparticle.pz(),
                                   this3BodyMCInfo.genDecVtx[0], this3BodyMCInfo.genDecVtx[1], this3BodyMCInfo.genDecVtx[2],
                                   this3BodyMCInfo.genCt,
                                   mcparticle.phi(), mcparticle.eta(), mcparticle.y(),
                                   this3BodyMCInfo.genMomProton, this3BodyMCInfo.genMomPion, this3BodyMCInfo.genMomDeuteron,
                                   this3BodyMCInfo.genPtProton, this3BodyMCInfo.genPtPion, this3BodyMCInfo.genPtDeuteron,
                                   this3BodyMCInfo.isTrueH3L, this3BodyMCInfo.isTrueAntiH3L,
                                   this3BodyMCInfo.isReco,
                                   mcparticle.pdgCode(),
                                   this3BodyMCInfo.daughterPrPdgCode, this3BodyMCInfo.daughterPiPdgCode, this3BodyMCInfo.daughterDePdgCode,
                                   this3BodyMCInfo.isDeuteronPrimary,
                                   this3BodyMCInfo.survivedEventSel);
        } // enabled table check
      } // mcParticles loop
    } // constexpr requires mcParticles check
  }

  // ______________________________________________________________
  // function to build mixed decay3body candidates
  template <class TRedCollisions, class TRedTracks, typename TRedDecay3Bodys, typename TBinningType>
  void buildMixedCandidates(TRedDecay3Bodys const& decay3bodys, TBinningType const& binningType)
  {
    if (!mEnabledTables[kVtx3BodyDatas]) {
      return; // don't do if no request for decay3bodys in place
    }

    // Strictly upper index policy for decay3body objects binned by radius, phi
    for (const auto& [decay3body0, decay3body1] : selfPairCombinations(binningType, mixingOpts.n3bodyMixing, -1, decay3bodys)) {
      auto trackPos0 = decay3body0.template track0_as<TRedTracks>();
      auto trackNeg0 = decay3body0.template track1_as<TRedTracks>();
      auto trackDeuteron0 = decay3body0.template track2_as<TRedTracks>();
      auto trackPos1 = decay3body1.template track0_as<TRedTracks>();
      auto trackNeg1 = decay3body1.template track1_as<TRedTracks>();
      auto trackDeuteron1 = decay3body1.template track2_as<TRedTracks>();

      // assign tracks
      auto trackProton0 = trackPos0;
      auto trackPion0 = trackNeg0;
      auto trackProton1 = trackPos1;
      auto trackPion1 = trackNeg1;
      if (trackDeuteron0.sign() < 0) {
        trackProton0 = trackNeg0;
        trackPion0 = trackPos0;
      }
      if (trackDeuteron1.sign() < 0) {
        trackProton1 = trackNeg1;
        trackPion1 = trackPos1;
      }

      registry.fill(HIST("Mixing/h3bodyCombinationCounter"), 0.5);

      // only combine if from different event
      if (decay3body0.collisionId() == decay3body1.collisionId()) {
        continue;
      }
      registry.fill(HIST("Mixing/h3bodyCombinationCounter"), 1.5);

      // collision vertex selection
      auto collision0 = decay3body0.template collision_as<TRedCollisions>();
      auto collision1 = decay3body1.template collision_as<TRedCollisions>();

      // get b_z value for each collision (from CCDB or cache) and cache it for that run number
      float magFieldCol0 = getMagFieldFromRunNumber(collision0.runNumber());
      float magFieldCol1 = getMagFieldFromRunNumber(collision1.runNumber());

      // only combine if collision similar in VtxZ
      if (mixingOpts.selectPVPosZ3bodyMixing && std::abs(collision0.posZ() - collision1.posZ()) > mixingOpts.maxDeltaPVPosZ3bodyMixing) {
        continue;
      }
      registry.fill(HIST("Mixing/h3bodyCombinationCounter"), 2.5);

      // Charge selections
      // same magnetic fields --> mix matter with matter
      if ((magFieldCol0 / std::abs(magFieldCol0)) == (magFieldCol1 / std::abs(magFieldCol1))) {
        if (trackDeuteron0.sign() != trackDeuteron1.sign()) {
          continue;
        }
      }
      // opposite magnetic fields --> mix matter with anti-matter
      if ((magFieldCol0 / std::abs(magFieldCol0)) != (magFieldCol1 / std::abs(magFieldCol1))) {
        if (trackDeuteron0.sign() == trackDeuteron1.sign()) {
          continue;
        }
      }

      // don't mix 3body with itself
      if ((trackDeuteron0.globalIndex() == trackDeuteron1.globalIndex()) || (trackProton0.globalIndex() == trackProton1.globalIndex()) || (trackPion0.globalIndex() == trackPion1.globalIndex())) {
        continue;
      }
      registry.fill(HIST("Mixing/h3bodyCombinationCounter"), 3.5);

      // candidate analysis
      // mix deuteron
      if (mixingOpts.mixingType == 0) {
        doMixing(collision0, trackProton0, trackPion0, trackDeuteron1, magFieldCol0);
        doMixing(collision1, trackProton1, trackPion1, trackDeuteron0, magFieldCol1);
      }
      // mix proton
      if (mixingOpts.mixingType == 1) {
        doMixing(collision0, trackProton1, trackPion0, trackDeuteron0, magFieldCol0);
        doMixing(collision1, trackProton0, trackPion1, trackDeuteron1, magFieldCol1);
      }
      // mix pion
      if (mixingOpts.mixingType == 2) {
        doMixing(collision0, trackProton0, trackPion1, trackDeuteron0, magFieldCol0);
        doMixing(collision1, trackProton1, trackPion0, trackDeuteron1, magFieldCol1);
      }
    } // end decay3body combinations loop
  }

  // ______________________________________________________________
  // function to calculate correct TOF nSigma for deuteron track
  template <bool isMC, class TCollisionTo, typename TCollision, typename TTrack>
  double getTOFnSigma(o2::pid::tof::TOFResoParamsV3 const& parameters, TCollision const& collision, TTrack const& track)
  {
    // TOF PID of deuteron
    if (track.has_collision() && track.hasTOF()) {
      auto originalcol = track.template collision_as<TCollisionTo>();
      if constexpr (isMC) {
        return bachelorTOFPIDLabeled.GetTOFNSigma(parameters, track, originalcol, collision);
      } else {
        return bachelorTOFPID.GetTOFNSigma(parameters, track, originalcol, collision);
      }
    }
    return -999;
  }

  // ______________________________________________________________
  // function to fill analysis tables
  void fillAnalysisTables()
  {
    // generate analysis tables
    if (mEnabledTables[kDecay3BodyIndices]) {
      products.decay3bodyindices(helper.decay3body.decay3bodyID,
                                 helper.decay3body.protonID, helper.decay3body.pionID, helper.decay3body.deuteronID,
                                 helper.decay3body.collisionID);
      registry.fill(HIST("Counters/hTableBuildingStatistics"), kDecay3BodyIndices);
    }
    if (mEnabledTables[kVtx3BodyDatas]) {
      products.vtx3bodydatas(helper.decay3body.sign,
                             helper.decay3body.mass, helper.decay3body.massV0,
                             helper.decay3body.position[0], helper.decay3body.position[1], helper.decay3body.position[2],
                             helper.decay3body.momentum[0], helper.decay3body.momentum[1], helper.decay3body.momentum[2],
                             helper.decay3body.chi2,
                             helper.decay3body.trackedClSize,
                             helper.decay3body.momProton[0], helper.decay3body.momProton[1], helper.decay3body.momProton[2],
                             helper.decay3body.momPion[0], helper.decay3body.momPion[1], helper.decay3body.momPion[2],
                             helper.decay3body.momDeuteron[0], helper.decay3body.momDeuteron[1], helper.decay3body.momDeuteron[2],
                             helper.decay3body.trackDCAxyToPV[0], helper.decay3body.trackDCAxyToPV[1], helper.decay3body.trackDCAxyToPV[2],             // 0 - proton, 1 - pion, 2 - deuteron
                             helper.decay3body.trackDCAToPV[0], helper.decay3body.trackDCAToPV[1], helper.decay3body.trackDCAToPV[2],                   // 0 - proton, 1 - pion, 2 - deuteron
                             helper.decay3body.trackDCAxyToPVprop[0], helper.decay3body.trackDCAxyToPVprop[1], helper.decay3body.trackDCAxyToPVprop[2], // 0 - proton, 1 - pion, 2 - deuteron
                             helper.decay3body.trackDCAToPVprop[0], helper.decay3body.trackDCAToPVprop[1], helper.decay3body.trackDCAToPVprop[2],       // 0 - proton, 1 - pion, 2 - deuteron
                             helper.decay3body.daughterDCAtoSV[0], helper.decay3body.daughterDCAtoSV[1], helper.decay3body.daughterDCAtoSV[2],          // 0 - proton, 1 - pion, 2 - deuteron
                             helper.decay3body.daughterDCAtoSVaverage,
                             helper.decay3body.cosPA, helper.decay3body.ctau,
                             helper.decay3body.tpcNsigma[0], helper.decay3body.tpcNsigma[1], helper.decay3body.tpcNsigma[2], helper.decay3body.tpcNsigma[2], // 0 - proton, 1 - pion, 2 - deuteron, 3 - bach with pion hyp
                             helper.decay3body.tofNsigmaDeuteron,
                             helper.decay3body.averageITSClSize[0], helper.decay3body.averageITSClSize[1], helper.decay3body.averageITSClSize[2], // 0 - proton, 1 - pion, 2 - deuteron
                             helper.decay3body.tpcNCl[0], helper.decay3body.tpcNCl[1], helper.decay3body.tpcNCl[2],                               // 0 - proton, 1 - pion, 2 - deuteron
                             helper.decay3body.pidForTrackingDeuteron);
      registry.fill(HIST("Counters/hTableBuildingStatistics"), kVtx3BodyDatas);
    }
    if (mEnabledTables[kVtx3BodyCovs]) {
      products.vtx3bodycovs(helper.decay3body.covProton,
                            helper.decay3body.covPion,
                            helper.decay3body.covDeuteron,
                            helper.decay3body.covariance);
      registry.fill(HIST("Counters/hTableBuildingStatistics"), kVtx3BodyCovs);
    }
    if (mEnabledTables[kMcVtx3BodyDatas]) {
      products.mcvtx3bodydatas(helper.decay3body.sign,
                               helper.decay3body.mass, helper.decay3body.massV0,
                               helper.decay3body.position[0], helper.decay3body.position[1], helper.decay3body.position[2],
                               helper.decay3body.momentum[0], helper.decay3body.momentum[1], helper.decay3body.momentum[2],
                               helper.decay3body.chi2,
                               helper.decay3body.trackedClSize,
                               helper.decay3body.momProton[0], helper.decay3body.momProton[1], helper.decay3body.momProton[2],
                               helper.decay3body.momPion[0], helper.decay3body.momPion[1], helper.decay3body.momPion[2],
                               helper.decay3body.momDeuteron[0], helper.decay3body.momDeuteron[1], helper.decay3body.momDeuteron[2],
                               helper.decay3body.trackDCAxyToPV[0], helper.decay3body.trackDCAxyToPV[1], helper.decay3body.trackDCAxyToPV[2],             // 0 - proton, 1 - pion, 2 - deuteron
                               helper.decay3body.trackDCAToPV[0], helper.decay3body.trackDCAToPV[1], helper.decay3body.trackDCAToPV[2],                   // 0 - proton, 1 - pion, 2 - deuteron
                               helper.decay3body.trackDCAxyToPVprop[0], helper.decay3body.trackDCAxyToPVprop[1], helper.decay3body.trackDCAxyToPVprop[2], // 0 - proton, 1 - pion, 2 - deuteron
                               helper.decay3body.trackDCAToPVprop[0], helper.decay3body.trackDCAToPVprop[1], helper.decay3body.trackDCAToPVprop[2],       // 0 - proton, 1 - pion, 2 - deuteron
                               helper.decay3body.daughterDCAtoSV[0], helper.decay3body.daughterDCAtoSV[1], helper.decay3body.daughterDCAtoSV[2],          // 0 - proton, 1 - pion, 2 - deuteron
                               helper.decay3body.daughterDCAtoSVaverage,
                               helper.decay3body.cosPA, helper.decay3body.ctau,
                               helper.decay3body.tpcNsigma[0], helper.decay3body.tpcNsigma[1], helper.decay3body.tpcNsigma[2], helper.decay3body.tpcNsigma[2], // 0 - proton, 1 - pion, 2 - deuteron, 3 - bach with pion hyp
                               helper.decay3body.tofNsigmaDeuteron,
                               helper.decay3body.averageITSClSize[0], helper.decay3body.averageITSClSize[1], helper.decay3body.averageITSClSize[2], // 0 - proton, 1 - pion, 2 - deuteron
                               helper.decay3body.tpcNCl[0], helper.decay3body.tpcNCl[1], helper.decay3body.tpcNCl[2],                               // 0 - proton, 1 - pion, 2 - deuteron
                               helper.decay3body.pidForTrackingDeuteron,
                               // MC information
                               this3BodyMCInfo.genMomentum[0], this3BodyMCInfo.genMomentum[1], this3BodyMCInfo.genMomentum[2],
                               this3BodyMCInfo.genDecVtx[0], this3BodyMCInfo.genDecVtx[1], this3BodyMCInfo.genDecVtx[2],
                               this3BodyMCInfo.genCt,
                               this3BodyMCInfo.genPhi, this3BodyMCInfo.genEta, this3BodyMCInfo.genRapidity,
                               this3BodyMCInfo.genMomProton, this3BodyMCInfo.genMomPion, this3BodyMCInfo.genMomDeuteron,
                               this3BodyMCInfo.genPtProton, this3BodyMCInfo.genPtPion, this3BodyMCInfo.genPtDeuteron,
                               this3BodyMCInfo.isTrueH3L, this3BodyMCInfo.isTrueAntiH3L,
                               this3BodyMCInfo.isReco,
                               this3BodyMCInfo.motherPdgCode,
                               this3BodyMCInfo.daughterPrPdgCode, this3BodyMCInfo.daughterPiPdgCode, this3BodyMCInfo.daughterDePdgCode,
                               this3BodyMCInfo.isDeuteronPrimary,
                               this3BodyMCInfo.survivedEventSel);
      registry.fill(HIST("Counters/hTableBuildingStatistics"), kMcVtx3BodyDatas);
    }
  }

  // ______________________________________________________________
  // function to build mixed 3body candidate from selected tracks
  template <typename TCollision, typename TTrack>
  void doMixing(TCollision const& collision, TTrack const& trackProton, TTrack const& trackPion, TTrack const& trackDeuteron, float magField)
  {
    // set vertexers and propagator with correct mag field of this collision (only if run number changed compared to previous candidate build)
    initFittersWithMagField(collision.runNumber(), magField);
    if (helper.buildDecay3BodyCandidate(collision, trackProton, trackPion, trackDeuteron,
                                        -1 /*decay3bodyIndex*/,
                                        trackDeuteron.tofNSigmaDe(),
                                        0 /*trackedClSize*/,
                                        decay3bodyBuilderOpts.useKFParticle,
                                        decay3bodyBuilderOpts.kfSetTopologicalConstraint,
                                        decay3bodyBuilderOpts.useSelections,
                                        decay3bodyBuilderOpts.useChi2Selection,
                                        decay3bodyBuilderOpts.useTPCforPion,
                                        decay3bodyBuilderOpts.acceptTPCOnly,
                                        decay3bodyBuilderOpts.askOnlyITSMatch,
                                        decay3bodyBuilderOpts.calculateCovariance,
                                        true, /*isEventMixing*/
                                        mixingOpts.doApplySVertexerCuts /*applySVertexerCuts*/)) {
      // fill analysis tables with built candidate
      fillAnalysisTables();
      return;
    } else {
      return;
    }
  }

  // ______________________________________________________________
  // function to check if a reconstructed mother is a true H3L/Anti-H3L (returns -1 if not)
  template <typename MCTrack3B>
  int checkH3LTruth(MCTrack3B const& mcParticlePr, MCTrack3B const& mcParticlePi, MCTrack3B const& mcParticleDe, bool& isMuonReco)
  {
    if (std::abs(mcParticlePr.pdgCode()) != PDG_t::kProton || std::abs(mcParticleDe.pdgCode()) != o2::constants::physics::Pdg::kDeuteron) {
      return -1;
    }
    // check proton and deuteron mother
    int prDeMomID = -1;
    for (const auto& motherPr : mcParticlePr.template mothers_as<aod::McParticles>()) {
      for (const auto& motherDe : mcParticleDe.template mothers_as<aod::McParticles>()) {
        if (motherPr.globalIndex() == motherDe.globalIndex() && std::abs(motherPr.pdgCode()) == o2::constants::physics::Pdg::kHyperTriton) {
          prDeMomID = motherPr.globalIndex();
          break;
        }
      }
    }
    if (prDeMomID == -1) {
      return -1;
    }
    if (std::abs(mcParticlePi.pdgCode()) != PDG_t::kPiPlus && std::abs(mcParticlePi.pdgCode()) != PDG_t::kMuonMinus) {
      return -1;
    }
    // check if the pion track is a muon coming from a pi -> mu + vu decay, if yes, take the mother pi
    auto mcParticlePiTmp = mcParticlePi;
    if (std::abs(mcParticlePiTmp.pdgCode()) == PDG_t::kMuonMinus) {
      for (const auto& motherPi : mcParticlePiTmp.template mothers_as<aod::McParticles>()) {
        if (std::abs(motherPi.pdgCode()) == PDG_t::kPiPlus) {
          mcParticlePiTmp = motherPi;
          isMuonReco = true;
          break;
        }
      }
    }
    // now loop over the pion mother
    for (const auto& motherPi : mcParticlePiTmp.template mothers_as<aod::McParticles>()) {
      if (motherPi.globalIndex() == prDeMomID) {
        return motherPi.globalIndex();
      }
    }
    return -1;
  }

  // ______________________________________________________________
  // function to reset MCInfo
  void resetMCInfo(mc3Bodyinfo& mcInfo)
  {
    mcInfo.label = -1;
    mcInfo.genMomentum[0] = -1., mcInfo.genMomentum[1] = -1., mcInfo.genMomentum[2] = -1.;
    mcInfo.genDecVtx[0] = -1., mcInfo.genDecVtx[1] = -1., mcInfo.genDecVtx[2] = -1.;
    mcInfo.genCt = -1.;
    mcInfo.genPhi = -1., mcInfo.genEta = -1., mcInfo.genRapidity = -1.;
    mcInfo.genMomProton = -1., mcInfo.genMomPion = -1., mcInfo.genMomDeuteron = -1.;
    mcInfo.genPtProton = -1., mcInfo.genPtPion = -1., mcInfo.genPtDeuteron = -1.;
    mcInfo.isTrueH3L = false, mcInfo.isTrueAntiH3L = false;
    mcInfo.isReco = false;
    mcInfo.motherPdgCode = -1;
    mcInfo.daughterPrPdgCode = -1, mcInfo.daughterPiPdgCode = -1, mcInfo.daughterDePdgCode = -1;
    mcInfo.isDeuteronPrimary = false;
    return;
  }

  // ______________________________________________________________
  // process functions
  void processRealData(ColswithEvTimes const& collisions,
                       aod::Decay3Bodys const& decay3bodys,
                       aod::Tracked3Bodys const& tracked3bodys,
                       TracksExtPIDIUwithEvTimes const&,
                       aod::BCsWithTimestamps const& bcs)
  {
    // initialise CCDB from BCs
    if (!initCCDB(bcs, collisions)) {
      LOG(info) << "CCDB initialisation failed, skipping candidate building." << std::endl;
      return;
    }

    // get tracked cluster size info
    fTrackedClSizeVector.clear();
    fTrackedClSizeVector.resize(decay3bodys.size(), 0);
    for (const auto& tvtx3body : tracked3bodys) {
      fTrackedClSizeVector[tvtx3body.decay3BodyId()] = tvtx3body.itsClsSize();
    }

    // do candidate analysis without MC processing
    buildCandidates<TracksExtPIDIUwithEvTimes>(bcs,                             // bc table
                                               collisions,                      // collision table
                                               decay3bodys,                     // decay3body table
                                               static_cast<TObject*>(nullptr),  // MC particle table
                                               static_cast<TObject*>(nullptr)); // MC collision table
  }

  void processRealDataReduced(aod::RedCollisions const& collisions,
                              soa::Join<aod::RedDecay3Bodys, aod::Red3BodyInfo> const& decay3bodys,
                              aod::RedIUTracks const&)
  {
    // get tracked cluster size info (saved in aod::Red3BodyInfo)
    fTrackedClSizeVector.clear();
    fTrackedClSizeVector.resize(decay3bodys.size(), 0);
    for (const auto& vtx3body : decay3bodys) {
      fTrackedClSizeVector[vtx3body.globalIndex()] = vtx3body.trackedClSize();
    }

    // do candidate analysis without MC processing
    buildCandidates<aod::RedIUTracks>(static_cast<TObject*>(nullptr),  // bc table
                                      collisions,                      // collision table
                                      decay3bodys,                     // decay3body table
                                      static_cast<TObject*>(nullptr),  // MC particle table
                                      static_cast<TObject*>(nullptr)); // MC collision table
  }

  void processRealDataReduced3bodyMixing(aod::RedCollisions const&,
                                         soa::Join<aod::RedDecay3Bodys, aod::Red3BodyInfo> const& decay3bodys,
                                         aod::RedIUTracks const&)
  {
    auto xAxis = registry.get<TH2>(HIST("Mixing/hDecay3BodyRadiusPhi"))->GetXaxis();
    auto yAxis = registry.get<TH2>(HIST("Mixing/hDecay3BodyRadiusPhi"))->GetYaxis();

    for (const auto& decay3body : decay3bodys) {
      int bin_Radius, bin_Phi;
      if (decay3bodyBuilderOpts.useKFParticle) {
        bin_Radius = xAxis->FindBin(decay3body.radiusKF());
        bin_Phi = yAxis->FindBin(decay3body.phiKF());
        registry.fill(HIST("Mixing/hDecay3BodyPosZ"), decay3body.poszKF());
      } else {
        bin_Radius = xAxis->FindBin(decay3body.radiusDCA());
        bin_Phi = yAxis->FindBin(decay3body.phiDCA());
        registry.fill(HIST("Mixing/hDecay3BodyPosZ"), decay3body.poszDCA());
      }
      registry.fill(HIST("Mixing/hDecay3BodyRadiusPhi"), xAxis->GetBinCenter(bin_Radius), yAxis->GetBinCenter(bin_Phi));
    }

    if (decay3bodyBuilderOpts.useKFParticle) {
      Binning3BodyKF binningOnRadPhiKF{{mixingOpts.bins3BodyRadius, mixingOpts.bins3BodyPhi}, true};
      buildMixedCandidates<aod::RedCollisions, aod::RedIUTracks>(decay3bodys, binningOnRadPhiKF);
    } else {
      Binning3BodyDCAfitter binningOnRadPhiDCA{{mixingOpts.bins3BodyRadius, mixingOpts.bins3BodyPhi}, true};
      buildMixedCandidates<aod::RedCollisions, aod::RedIUTracks>(decay3bodys, binningOnRadPhiDCA);
    }
  }

  void processMonteCarlo(ColswithEvTimesLabeled const& collisions,
                         aod::Decay3Bodys const& decay3bodys,
                         aod::Tracked3Bodys const& tracked3bodys,
                         TracksExtPIDIUwithEvTimesLabeled const&,
                         aod::BCsWithTimestamps const& bcs,
                         aod::McParticles const& mcParticles,
                         aod::McCollisions const& mcCollisions)
  {
    // initialise CCDB from BCs
    if (!initCCDB(bcs, collisions)) {
      LOG(info) << "CCDB initialisation failed, skipping candidate building." << std::endl;
      return;
    }

    // get tracked cluster size info
    fTrackedClSizeVector.clear();
    fTrackedClSizeVector.resize(decay3bodys.size(), 0);
    for (const auto& tvtx3body : tracked3bodys) {
      fTrackedClSizeVector[tvtx3body.decay3BodyId()] = tvtx3body.itsClsSize();
    }

    // do candidate analysis with MC processing
    buildCandidates<TracksExtPIDIUwithEvTimesLabeled>(bcs,           // bc table
                                                      collisions,    // collision table
                                                      decay3bodys,   // decay3body table
                                                      mcParticles,   // MC particle table
                                                      mcCollisions); // MC collision table
  }

  PROCESS_SWITCH(decay3bodyBuilder, processRealData, "process real data", true);
  PROCESS_SWITCH(decay3bodyBuilder, processRealDataReduced, "process real reduced data", false);
  PROCESS_SWITCH(decay3bodyBuilder, processRealDataReduced3bodyMixing, "process real reduced data", false);
  PROCESS_SWITCH(decay3bodyBuilder, processMonteCarlo, "process monte carlo", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{
    adaptAnalysisTask<decay3bodyBuilder>(cfgc)};
}
