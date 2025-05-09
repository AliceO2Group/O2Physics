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
/// \brief Builder task for 3-body decay reconstruction (p + pion + bachelor)
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>
/// \author Carolina Reetz <c.reetz@cern.ch>
// ========================

#include <cmath>
#include <array>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <algorithm>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "DetectorsVertexing/SVertexHypothesis.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/Reduced3BodyTables.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "PWGLF/DataModel/pidTOFGeneric.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/PID/PIDTOF.h"
#include "TableHelper.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DataFormatsCalibration/MeanVertexObject.h"

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

static constexpr int nParameters = 1;
static const std::vector<std::string> tableNames {
  "Decay3BodyIndices",
  "StoredDecay3BodyCores",
  "Decay3BodyCovs",
  "Decay3BodyDaughterCovs",
  "McDecay3BodyLabels"
};

static constexpr int nTablesConst = 5;

static const std::vector<std::string> parameterNames{"enable"};
static const int defaultParameters[nTablesConst][nParameters]{
  {0}, // Decay3BodyIndices
  {0}, // StoredDecay3BodyCores
  {0}, // Decay3BodyCovs
  {0}, // Decay3BodyDaughterCovs
  {0}  // McDecay3BodyLabels
};

using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using FullTracksExtPIDIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe>;
using FullTracksExtIULabeled = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::McTrackLabels>;
using FullTracksExtIUwithEvTimes = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::EvTimeTOFFT0ForTrack>;
using FullTracksExtPIDIUwithEvTimes = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe, aod::EvTimeTOFFT0ForTrack>;

using ColwithEvTimes = o2::soa::Join<aod::Collisions, aod::EvSels, aod::EvTimeTOFFT0>;
using ColwithEvTimesMults = o2::soa::Join<aod::Collisions, aod::EvSels, aod::EvTimeTOFFT0, aod::Mults>;

using ReducedCollisionsMults = soa::Join<aod::RedCollisions, aod::RedPVMults>;
using ReducedCollisionsMultsCents = soa::Join<aod::RedCollisions, aod::RedPVMults, aod::RedCentFT0Cs>;

// namespace
// {
// const float pidCutsLambda[o2::vertexing::SVertexHypothesis::NPIDParams] = {0., 20, 0., 5.0, 0.0, 1.09004e-03, 2.62291e-04, 8.93179e-03, 2.83121}; // Lambda
// } // namespace

// struct VtxCandidate {
//   int track0Id;
//   int track1Id;
//   int track2Id;
//   int collisionId;
//   int decay3bodyId;
//   float vtxPos[3];
//   float track0P[3];
//   float track1P[3];
//   float track2P[3];
//   float dcadaughters;
//   float daudcaxytopv[3]; // 0 - proton, 1 - pion, 2 - bachelor
//   float daudcatopv[3];   // 0 - proton, 1 - pion, 2 - bachelor
//   float bachelortofNsigma;
// };

// struct kfCandidate {
//   // hypertriton
//   int collisionID;
//   int trackPosID;
//   int trackNegID;
//   int trackBachID;
//   int decay3bodyID;
//   float mass;
//   float pos[3];
//   float posErr[3];
//   float mom[4];
//   float momErr[4];
//   float charge;
//   float dcaToPV[2];     // 3D, xy
//   float cpaToPV[2];     // 3D, xy
//   float cpaToPVtopo[2]; // 3D, xy
//   float decLen[2];      // 3D, xy
//   float ldl;
//   float chi2geoNDF;
//   float chi2topoNDF;
//   float ctau;
//   float trackedClSize;
//   float DeltaPhiRotDeuteron;
//   float DeltaPhiRotProton;
//   // V0
//   float massV0;
//   float chi2massV0;
//   float cpaV0ToPV;
//   // daughter momenta
//   float protonMom[3];
//   float pionMom[3];
//   float deuteronMom[3];
//   float tpcInnerParam[3]; // proton, pion, deuteron
//   // daughter track quality
//   int tpcNClDaughters[3]; // proton, pion, deuteron
//   float tpcChi2NClDeuteron;
//   // daughter DCAs KF
//   float DCAdaughterToPV[3];   // proton, pion, deuteron
//   float DCAdaughterToPVxy[3]; // proton, pion, deuteron
//   float DCAdaughterToSVxy[3]; // proton, pion, deuteron
//   float DCAprotonToPion;
//   float DCAprotonToDeuteron;
//   float DCApionToDeuteron;
//   float DCAvtxDaughters3D;
//   // daughter DCAs to PV propagated with material
//   float trackDCAxy[3]; // pos, neg, bach
//   float trackDCA[3];   // pos, neg, bach
//   // daughter signs
//   float daughterCharge[3]; // proton, pion, deuteron
//   // daughter PID
//   float tpcNsigma[4]; // proton, pion, deuteron, bach with pion hyp
//   float tpcdEdx[3];   // proton, pion, deuteron
//   float tofNsigmaDeuteron;
//   float averageClusterSizeDeuteron;
//   float pidForTrackingDeuteron;
// };

struct decay3bodyBuilder {

  // helper object
  o2::pwglf::decay3bodyBuilderHelper helper;

  // table index : match order above
  enum tableIndex { kDecay3BodyIndices = 0,
                    kStoredDecay3BodyCores,
                    kDecay3BodyCovs,
                    kDecay3BodyDaughterCovs,
                    kMcDecay3BodyLabels,
                    kMcDecay3BodyCores,
                    kMcDecay3BodyCollRefs,
                    nTables };

  struct : ProducesGroup {
    Produces<aod::Decay3BodyIndices> decay3bodyindices;
    Produces<aod::StoredDecay3BodyCores> decay3bodycores;
    Produces<aod::Decay3BodyCovs> decay3bodycovs;
    Produces<aod::Decay3BodyDaughterCovs> decay3bodydaugcovs;
    Produces<aod::McDecay3BodyLabels> mcdecay3bodylabels;
  } products;

  Configurable<LabelArray<int>> enabledTables{"enabledTables", {defaultParameters[0], nTables, nParameters, tableNames, parameterNames}, "Produce this table: 0 - false, 1 - true"}
  std::vector<int> mEnabledTables; // Vector of enabled tables

  // general options
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
  Configurable<bool> doTrackQA{"doTrackQA", false, "Flag to fill QA histograms for daughter tracks."};
  Configurable<bool> doVertexQA{"doVertexQA", false, "Flag to fill QA histograms for PV."};

  // data processing options
  Configurable<bool> doSkimmedProcessing{"doSkimmedProcessing", false, "Apply Zoroo counting in case of skimmed data input"};
  Configurable<std::string> triggerList{"triggerList", "fTriggerEventF1Proton, fTrackedOmega, fTrackedXi, fOmegaLargeRadius, fDoubleOmega, fOmegaHighMult, fSingleXiYN, fQuadrupleXi, fDoubleXi, fhadronOmega, fOmegaXi, fTripleXi, fOmega, fGammaVeryLowPtEMCAL, fGammaVeryLowPtDCAL, fGammaHighPtEMCAL, fGammaLowPtEMCAL, fGammaVeryHighPtDCAL, fGammaVeryHighPtEMCAL, fGammaLowPtDCAL, fJetNeutralLowPt, fJetNeutralHighPt, fGammaHighPtDCAL, fJetFullLowPt, fJetFullHighPt, fEMCALReadout, fPCMandEE, fPHOSnbar, fPCMHighPtPhoton, fPHOSPhoton, fLD, fPPPHI, fPD, fLLL, fPLL, fPPL, fPPP, fLeadingPtTrack, fHighFt0cFv0Flat, fHighFt0cFv0Mult, fHighFt0Flat, fHighFt0Mult, fHighMultFv0, fHighTrackMult, fHfSingleNonPromptCharm3P, fHfSingleNonPromptCharm2P, fHfSingleCharm3P, fHfPhotonCharm3P, fHfHighPt2P, fHfSigmaC0K0, fHfDoubleCharm2P, fHfBeauty3P, fHfFemto3P, fHfFemto2P, fHfHighPt3P, fHfSigmaCPPK, fHfDoubleCharm3P, fHfDoubleCharmMix, fHfPhotonCharm2P, fHfV0Charm2P, fHfBeauty4P, fHfV0Charm3P, fHfSingleCharm2P, fHfCharmBarToXiBach, fSingleMuHigh, fSingleMuLow, fLMeeHMR, fDiMuon, fDiElectron, fLMeeIMR, fSingleE, fTrackHighPt, fTrackLowPt, fJetChHighPt, fJetChLowPt, fUDdiffLarge, fUDdiffSmall, fITSextremeIonisation, fITSmildIonisation, fH3L3Body, fHe, fH2", "List of triggers used to select events"};

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
    Configurable<bool> useSelections{"useSelections", true, "Apply selections during decay3body building"};
    Configurable<bool> acceptTPCOnly{"acceptTPCOnly", false, "Accept TPC only tracks as daughters"};
    Configurable<bool> calculateCovariance{"calculateCovariance", true, "Calculate candidate and daughter covariance matrices"};
    // daughter track selections
    Configurable<float> maxEtaDaughters{"maxEtaDaughters", 0.9, "Max eta of daughters"};
    Configurable<int> minTPCNClProton{"minTPCNClProton", 90, "Min TPC NClusters of proton daughter"};
    Configurable<int> minTPCNClPion{"minTPCNClPion", 70, "Min TPC NClusters of pion daughter"};
    Configurable<int> minTPCNClBach{"minTPCNClBach", 100, "Min TPC NClusters of bachelor daughter"};
    Configurable<float> minDCAProtonToPV{"minDCAProtonToPV", 0.1, "Min DCA of proton to PV"};
    Configurable<float> minDCAPionToPV{"minDCAPionToPV", 0.1, "Min DCA of pion to PV"};
    Configurable<float> minDCABachToPV{"minDCABachToPV", 0.1, "Min DCA of bachelor to PV"};
    Configurable<float> minPtProton{"minPtProton", 0.3, "Min Pt of proton daughter"};
    Configurable<float> minPtPion{"minPtPion", 0.1, "Min Pt of pion daughter"};
    Configurable<float> minPtBach{"minPtBach", 0.6, "Min Pt of bachelor daughter"};
    Configurable<float> maxPtProton{"maxPtProton", 5.0, "Max Pt of proton daughter"};
    Configurable<float> maxPtPion{"maxPtPion", 1.2, "Max Pt of pion daughter"};
    Configurable<float> maxPtBach{"maxPtBach", 10.0, "Max Pt of bachelor daughter"};
    Configurable<float> maxTPCnSigma{"maxTPCnSigma", 5.0, "Min/max TPC nSigma of daughter tracks"};
    Configurable<float> minTOFnSigmaDeuteron{"minTOFnSigmaDeuteron", -5.0, "Min TOF nSigma of deuteron daughter"};
    Configurable<float> maxTOFnSigmaDeuteron{"maxTOFnSigmaDeuteron", 5.0, "Max TOF nSigma of deuteron daughter"};
    Configurable<float> minPBachUseTOF{"minPBaChUseTOF", 1.0, "Min P of bachelor to use TOF PID"};
    Configurable<float> maxDCADauAtSV{"maxDCADauAtSV", 0.5, "Max DCA of daughters at SV (quadratic sum of daughter DCAs between each other)"};
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
    std::string prefix = "eventMixingOpts";
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
  } eventMixingOpts;

  struct : ConfigurableGroup {
    std::string prefix = "tofPIDOpts";
    Configurable<int64_t> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};
    Configurable<std::string> paramFileName{"paramFileName", "", "Path to the parametrization object. If empty the parametrization is not taken from file"};
    Configurable<std::string> parametrizationPath{"parametrizationPath", "TOF/Calib/Params", "Path of the TOF parametrization on the CCDB or in the file, if the paramFileName is not empty"};
    Configurable<std::string> passName{"passName", "", "Name of the pass inside of the CCDB parameter collection. If empty, the automatically deceted from metadata (to be implemented!!!)"};
    Configurable<std::string> timeShiftCCDBPath{"timeShiftCCDBPath", "", "Path of the TOF time shift vs eta. If empty none is taken"};
    Configurable<bool> loadResponseFromCCDB{"loadResponseFromCCDB", false, "Flag to load the response from the CCDB"};
    Configurable<bool> fatalOnPassNotAvailable{"fatalOnPassNotAvailable", true, "Flag to throw a fatal if the pass is not available in the retrieved CCDB object"};
  } tofPIDOpts;
  

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  int mRunNumber;
  o2::base::MatLayerCylSel* lut = nullptr;

  HistogramRegistry registry{"Registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::vector<std::size_t> sorted_decay3body;

  // bachelor TOF PID
  o2::aod::pidtofgeneric::TofPidNewCollision<TrackExtPIDIUwithEvTimes::iterator> bachelorTOFPID; // to be updated in Init base on the hypothesis
  o2::pid::tof::TOFResoParamsV2 mRespParamsV2;

  // event mixing selection
  std::array<o2::vertexing::SVertexHypothesis, 2> mV0Hyps; // 0 - Lambda, 1 - AntiLambda
  bool doUpdateGRPMagField = false;                        // if initialize magnetic field for each bc

  // skimmed processing
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};


  // ==========================
  std::vector<VtxCandidate> VtxCandidates;

  std::unordered_map<int, float> ccdbCache;                                          // Maps runNumber -> d_bz
  std::unordered_map<int, std::shared_ptr<o2::parameters::GRPMagField>> grpMagCache; // Maps runNumber -> grpmap


  std::vector<int> fTrackedClSizeVector;

  // Configurables
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};

  enum Hyp3Body { kH3L = 0,
                  kH4L,
                  kHe4L,
                  kHe5L,
                  kNHyp3body };

  enum VtxStep { kVtxAll = 0,
                 kVtxTPCNcls,
                 kVtxPIDCut,
                 kVtxhasSV,
                 kVtxDcaDau,
                 kVtxCosPA,
                 kNVtxSteps };

  enum kfvtxstep { kKfVtxAll = 0,
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



  // for KFParticle reconstruction
  struct : ConfigurableGroup {
    Configurable<bool> cfgOnlyKeepInterestedTrigger{"kfparticleConfigurations.cfgOnlyKeepInterestedTrigger", false, "Flag to keep only interested trigger"};
    Configurable<bool> fillCandidateFullTable{"kfparticleConfigurations.fillCandidateFullTable", false, "Switch to fill full table with candidate properties"};
    Configurable<bool> doSel8selection{"kfparticleConfigurations.doSel8selection", true, "flag for sel8 event selection"};
    Configurable<bool> doPosZselection{"kfparticleConfigurations.doPosZselection", true, "flag for posZ event selection"};
  
    Configurable<bool> applyTopoSel{"kfparticleConfigurations.applyTopoSel", false, "Apply selection constraining the mother to the PV with KFParticle"};
    Configurable<float> maxChi2topo{"kfparticleConfigurations.maxChi2topo", 1000., "Maximum chi2 topological with KFParticle"};
    // 3body mixing
    Configurable<int> mixingType{"kfparticleConfigurations.mixingType", 0, "0: mix V0 from one event with bachelor from another, 1: mix pion and bachelor from one event with proton from another "};
    ConfigurableAxis bins3BodyRadius{"kfparticleConfigurations.bins3BodyRadius", {VARIABLE_WIDTH, 0.0f, 0.5f, 1.0f, 1.5f, 2.0f, 3.0f, 4.0f, 6.0f, 8.0f, 10.0f, 12.0f, 14.0f, 16.0f, 18.0f, 20.0f, 30.0f, 1000.0}, "Mixing bins - 3body radius"};
    // ConfigurableAxis bins3BodyPhi{"kfparticleConfigurations.bins3BodyPhi", {VARIABLE_WIDTH, -180.0f*TMath::Pi()/180, -170.0f*TMath::Pi()/180, -160.0f*TMath::Pi()/180, -150.0f*TMath::Pi()/180, -140.0f*TMath::Pi()/180, -130.0f*TMath::Pi()/180, -120.0f*TMath::Pi()/180, -110.0f*TMath::Pi()/180, -100.0f*TMath::Pi()/180, -90.0f*TMath::Pi()/180, -80.0f*TMath::Pi()/180, -70.0f*TMath::Pi()/180, -60.0f*TMath::Pi()/180, -50.0f*TMath::Pi()/180, -40.0f*TMath::Pi()/180, -30.0f*TMath::Pi()/180, -20.0f*TMath::Pi()/180, -10.0f*TMath::Pi()/180, 0.0f, 10.0f*TMath::Pi()/180, 20.0f*TMath::Pi()/180, 30.0f*TMath::Pi()/180, 40.0f*TMath::Pi()/180, 50.0f*TMath::Pi()/180, 60.0f*TMath::Pi()/180, 70.0f*TMath::Pi()/180, 80.0f*TMath::Pi()/180, 90.0f*TMath::Pi()/180, 100.0f*TMath::Pi()/180, 110.0f*TMath::Pi()/180, 120.0f*TMath::Pi()/180, 130.0f*TMath::Pi()/180, 140.0f*TMath::Pi()/180, 150.0f*TMath::Pi()/180, 160.0f*TMath::Pi()/180, 170.0f*TMath::Pi()/180, 180.0f*TMath::Pi()/180}, "Mixing bins - 3body phi"};
    ConfigurableAxis bins3BodyPhi{"kfparticleConfigurations.bins3BodyPhi", {VARIABLE_WIDTH, -180.0f * TMath::Pi() / 180, -160.0f * TMath::Pi() / 180, -140.0f * TMath::Pi() / 180, -120.0f * TMath::Pi() / 180, -100.0f * TMath::Pi() / 180, -80.0f * TMath::Pi() / 180, -60.0f * TMath::Pi() / 180, -40.0f * TMath::Pi() / 180, -20.0f * TMath::Pi() / 180, 0.0f, 20.0f * TMath::Pi() / 180, 40.0f * TMath::Pi() / 180, 60.0f * TMath::Pi() / 180, 80.0f * TMath::Pi() / 180, 100.0f * TMath::Pi() / 180, 120.0f * TMath::Pi() / 180, 140.0f * TMath::Pi() / 180, 160.0f * TMath::Pi() / 180, 180.0f * TMath::Pi() / 180}, "Mixing bins - 3body phi"};
    ConfigurableAxis bins3BodyPosZ{"kfparticleConfigurations.bins3BodyPosZ", {VARIABLE_WIDTH, -300.0f, -42.0f, -13.0f, -6.0f, -4.0f, -2.0f, 0.0f, 2.0f, 4.0f, 6.0f, 13.0f, 42.0f, 300.0f}, "Mixing bins - 3body z position"};
    Configurable<bool> selectVtxZ3bodyMixing{"kfparticleConfigurations.selectVtxZ3bodyMixing", true, "Select same VtxZ events in case of 3body mixing"};
    Configurable<float> VtxZBin3bodyMixing{"kfparticleConfigurations.VtxZBin3bodyMixing", 1., "Bin width for event vtx z position in case of 3body mixing"};
  } kfparticleConfigurations;

  //------------------------------------------------------------------
  // Sets for DCAFitter event mixing
  struct : ConfigurableGroup {
    // Binning for mixing events
    ConfigurableAxis binsVtxZ{"dcaFitterEMSel.binsVtxZ", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis binsMultiplicity{"dcaFitterEMSel.binsMultiplicity", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - multiplicity"};
    Configurable<float> maxDeltaRadiusColMixing{"dcaFitterEMSel.maxDeltaRadiusColMixing", 2., "max difference between pv z position in case of collision mixing"};
    Configurable<float> maxDeltaPhiColMixing{"dcaFitterEMSel.maxDeltaPhiColMixing", 30., "max difference between Phi of monther particle in case of collision mixing (degree)"};
    // Configurations for mixing decay3bodys
    // Configurable<bool> cfgUseDCAFitterInfo{"dcaFitterEMSel.cfgUseDCAFitterInfo", true, ""}; // if use information from dcatFitter while mixing reduced 3bodys
    Configurable<int> cfgMix3BodyMethod{"dcaFitterEMSel.cfgMix3BodyMethod", 0, ""}; // 0: bachelor, 1: pion, 2: proton
    Configurable<bool> cfgApplyV0Cut{"dcaFitterEMSel.cfgApplyV0Cut", true, "if apply V0 cut while performing event-mixing"};
    ConfigurableAxis bins3BodyRadius{"dcaFitterEMSel.bins3BodyRadius", {VARIABLE_WIDTH, 0.0f, 2.0f, 4.0f, 7.0f, 10.0f, 14.0f, 18.0f, 22.0f, 30.0f, 40.0f}, "Mixing bins - 3body radius"};
    ConfigurableAxis bins3BodyPhi{"dcaFitterEMSel.bins3BodyPhi", {VARIABLE_WIDTH, -3.15, -2.15, -1, 0, 1, 2.15, 3.15}, "Mixing bins - 3body phi"};
    ConfigurableAxis bins3BodyPhiDegree{"dcaFitterEMSel.bins3BodyPhiDegree", {VARIABLE_WIDTH, -180, -120, -60, 0, 60, 120, 180}, "Mixing bins - 3body phi"};
    ConfigurableAxis bins3BodyPosZ{"dcaFitterEMSel.bins3BodyPosZ", {VARIABLE_WIDTH, -500.0f, -200.0f, -100.0f, -70.0f, -60.0f, -50.0f, -40.0f, -35.0f, -30.0f, -25.0f, -20.0f, -15.0f, -13.0f, -10.0f, -8.0f, -6.0f, -4.0f, -2.0f, 0.0f, 2.0f, 4.0f, 6.0f, 8.0f, 10.0f, 13.0f, 15.0f, 20.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f, 60.0f, 70.0f, 100.0f, 200.0f, 500.0f}, "3body SV z position"};
    Configurable<bool> selectPVPosZ3bodyMixing{"dcaFitterEMSel.selectPVPosZ3bodyMixing", true, "Select same pvPosZ events in case of 3body mixing"};
    Configurable<float> maxDeltaPVPosZ3bodyMixing{"dcaFitterEMSel.maxDeltaPVPosZ3bodyMixing", 1., "max difference between pv z position in case of 3body mixing"};
  } dcaFitterEMSel;

  SliceCache cache;
  using BinningTypeColEM = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultNTracksPV>;
  using Binning3BodyDCAFitter = ColumnBinningPolicy<aod::dcafittersvinfo::SVRadius, aod::dcafittersvinfo::MomPhi>;
  using Binning3BodyKFInfo = ColumnBinningPolicy<aod::reduceddecay3body::Radius, aod::reduceddecay3body::Phi>;

  // KF event mixing
  using BinningTypeKF = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultNTracksPV>;

  // 3body mixing
  using Binning3Body = ColumnBinningPolicy<aod::reduceddecay3body::Radius, aod::reduceddecay3body::Phi>;

  // Filters and slices
  Preslice<aod::Decay3Bodys> perCollision = o2::aod::decay3body::collisionId;
  Preslice<aod::RedDecay3Bodys> perReducedCollision = o2::aod::reduceddecay3body::collisionId;

  int mRunNumber;
  float d_bz;

  o2::pid::tof::TOFResoParamsV2 mRespParamsV2;
  std::array<o2::vertexing::SVertexHypothesis, 2> mV0Hyps; // 0 - Lambda, 1 - AntiLambda
  bool doUpdateGRPMagField = false;                        // if initialize magnetic field for each bc
  o2::dataformats::VertexBase mMeanVertex{{0., 0., 0.}, {0.1 * 0.1, 0., 0.1 * 0.1, 0., 0., 6. * 6.}};


  // ==========================





  void init(InitContext&)
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

    fitterV0.setMatCorrType(matCorr);
    fitter3body.setMatCorrType(matCorr);
    
    // set bachelor PID
    bachelorTOFPID.SetPidType(o2::track::PID::Deuteron);

    // set decay3body parameters in the helper
    helper.decay3bodyselections.maxEtaDaughters = decay3bodyBuilderOpts.maxEtaDaughters;
    helper.decay3bodyselections.minTPCNClProton = decay3bodyBuilderOpts.minTPCNClProton;
    helper.decay3bodyselections.minTPCNClPion = decay3bodyBuilderOpts.minTPCNClPion;
    helper.decay3bodyselections.minTPCNClBach = decay3bodyBuilderOpts.minTPCNClBach;
    helper.decay3bodyselections.minDCAProtonToPV = decay3bodyBuilderOpts.minDCAProtonToPV;
    helper.decay3bodyselections.minDCAPionToPV = decay3bodyBuilderOpts.minDCAPionToPV;
    helper.decay3bodyselections.minDCABachToPV = decay3bodyBuilderOpts.minDCABachToPV;
    helper.decay3bodyselections.minPtProton = decay3bodyBuilderOpts.minPtProton;
    helper.decay3bodyselections.minPtPion = decay3bodyBuilderOpts.minPtPion;
    helper.decay3bodyselections.minPtBach = decay3bodyBuilderOpts.minPtBach;
    helper.decay3bodyselections.maxPtProton = decay3bodyBuilderOpts.maxPtProton;
    helper.decay3bodyselections.maxPtPion = decay3bodyBuilderOpts.maxPtPion;
    helper.decay3bodyselections.maxPtBach = decay3bodyBuilderOpts.maxPtBach;
    helper.decay3bodyselections.maxTPCnSigma = decay3bodyBuilderOpts.maxTPCnSigma;
    helper.decay3bodyselections.minTOFnSigmaDeuteron = decay3bodyBuilderOpts.minTOFnSigmaDeuteron;
    helper.decay3bodyselections.maxTOFnSigmaDeuteron = decay3bodyBuilderOpts.maxTOFnSigmaDeuteron;
    helper.decay3bodyselections.minPBachUseTOF = decay3bodyBuilderOpts.minPBachUseTOF;
    helper.decay3bodyselections.maxDCADauAtSV = decay3bodyBuilderOpts.maxDCADauAtSV;
    helper.decay3bodylelections.maxRapidity = decay3bodyBuilderOpts.maxRapidity;
    helper.decay3bodyselections.minPt = decay3bodyBuilderOpts.minPt;
    helper.decay3bodyselections.maxPt = decay3bodyBuilderOpts.maxPt;
    helper.decay3bodyselections.minMass = decay3bodyBuilderOpts.minMass;
    helper.decay3bodyselections.maxMass = decay3bodyBuilderOpts.maxMass;
    helper.decay3bodyselections.minCtau = decay3bodyBuilderOpts.minCtau;
    helper.decay3bodyselections.maxCtau = decay3bodyBuilderOpts.maxCtau;
    helper.decay3bodyselections.minCosPA = decay3bodyBuilderOpts.minCosPA;
    helper.decay3bodyselections.maxChi2 = decay3bodyBuilderOpts.maxChi2;

    // set SVertexer selection parameters in the helper
    helper.svselections.minPt2V0 = eventMixingOpts.minPt2V0;
    helper.svselections.maxTgl2V0 = eventMixingOpts.maxTgl2V0;
    helper.svselections.maxDCAXY2ToMeanVertex3bodyV0 = eventMixingOpts.maxDCAXY2ToMeanVertex3bodyV0;
    helper.svselections.minCosPAXYMeanVertex3bodyV0 = eventMixingOpts.minCosPAXYMeanVertex3bodyV0;
    helper.svselections.minCosPA3bodyV0 = eventMixingOpts.minCosPA3bodyV0;
    helper.svselections.maxRDiffV03body = eventMixingOpts.maxRDiffV03body;
    helper.svselections.minPt3Body = eventMixingOpts.minPt3Body;
    helper.svselections.maxTgl3Body = eventMixingOpts.maxTgl3Body;
    helper.svselections.maxDCAXY3Body = eventMixingOpts.maxDCAXY3Body;
    helper.svselections.maxDCAZ3Body = eventMixingOpts.maxDCAZ3Body;

    // configure tables to generate
    // setup bookkeeping histograms
    auto h = registry.add<TH1>("hTableBuildingStatistics", "hTableBuildingStatistics", kTH1D, {{nTablesConst, -0.5f, static_cast<float>(nTablesConst)}});
    auto h2 = registry.add<TH1>("hInputStatistics", "hInputStatistics", kTH1D, {{nTablesConst, -0.5f, static_cast<float>(nTablesConst)}});
    h2->SetTitle("Input table sizes");

    for (int i = 0; i < nTables; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, tableNames[i].c_str());
      h->SetBinContent(i + 1, -1); // mark all as disabled to start

      int f = enabledTables->get(tableNames[i].c_str(), "enable");
      if (f == 1) {
        mEnabledTables[i] = 1;
        h->SetBinContent(i + 1, 0); // mark enabled
      }
    }

    // list enabled process functions
    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
    LOGF(info, " Decay3body builder: basic configuration listing");
    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");

    if (doprocessRealData) {
      LOGF(info, " ===> process function enabled: processRealData");
    }
    if (doprocessRealDataReduced) {
      LOGF(info, " ===> process function enabled: processRealData");
    }
    if (doprocessMonteCarlo) {
      LOGF(info, " ===> process function enabled: processMonteCarlo");
    }
    if (doprocessEventMixing) {
      LOGF(info, " ===> process function enabled: processEventMixing");
    }

    // list enabled tables
    for (int i = 0; i < nTables; i++) {
      // printout to be improved in the future
      if (mEnabledTables[i]) {
        LOGF(info, " -~> Table enabled: %s", tableNames[i]);
      }
    }
    // print base cuts
    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
    LOGF(info, "-~> max daughter eta ..............: %i", decay3bodyBuilderOpts.maxEtaDaughters.value);
    LOGF(info, "-~> min TPC ncls proton ...........: %i", decay3bodyBuilderOpts.minTPCNClProton.value);
    LOGF(info, "-~> min TPC ncls pion .............: %i", decay3bodyBuilderOpts.minTPCNClPion.value);
    LOGF(info, "-~> min TPC ncls bach .............: %i", decay3bodyBuilderOpts.minTPCNClBach.value);
    LOGF(info, "-~> min DCA proton to PV ..........: %f", decay3bodyBuilderOpts.minDCAProtonToPV.value);
    LOGF(info, "-~> min DCA pion to PV ............: %f", decay3bodyBuilderOpts.minDCAPionToPV.value);
    LOGF(info, "-~> min DCA bach to PV ............: %f", decay3bodyBuilderOpts.minDCABachToPV.value);
    LOGF(info, "-~> min pT proton .................: %f", decay3bodyBuilderOpts.minPtProton.value);
    LOGF(info, "-~> min pT pion ...................: %f", decay3bodyBuilderOpts.minPtPion.value);
    LOGF(info, "-~> min pT bach ...................: %f", decay3bodyBuilderOpts.minPtBach.value);
    LOGF(info, "-~> max pT proton .................: %f", decay3bodyBuilderOpts.maxPtProton.value);
    LOGF(info, "-~> max pT pion ...................: %f", decay3bodyBuilderOpts.maxPtPion.value);
    LOGF(info, "-~> max pT bach ...................: %f", decay3bodyBuilderOpts.maxPtBach.value);
    LOGF(info, "-~> max TPC nSigma ...............: %f", decay3bodyBuilderOpts.maxTPCnSigma.value);
    LOGF(info, "-~> min TOF nSigma deuteron ......: %f", decay3bodyBuilderOpts.minTOFnSigmaDeuteron.value);
    LOGF(info, "-~> max TOF nSigma deuteron ......: %f", decay3bodyBuilderOpts.maxTOFnSigmaDeuteron.value);
    LOGF(info, "-~> min p bach use TOF ...........: %f", decay3bodyBuilderOpts.minPBachUseTOF.value);
    LOGF(info, "-~> max DCA dau at SV ............: %f", decay3bodyBuilderOpts.maxDCADauAtSV.value);
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

    // Add histograms separately for different process functions
    if (doprocessRealData == true || doprocessRealDataReduced == true) {
      registry.add("hEventCounter", "Counters/hEventCounter", HistType::kTH1F, {{1, 0.0f, 1.0f}});
    }

    if (doTrackQA) {
      registry.add("QA/Tracks/hTrackPosTPCNcls", "hTrackPosTPCNcls", HistType::kTH1F, {{152, 0, 152, "# TPC clusters"}});
      registry.add("QA/Tracks/hTrackNegTPCNcls", "hTrackNegTPCNcls", HistType::kTH1F, {{152, 0, 152, "# TPC clusters"}});
      registry.add("QA/Tracks/hTrackBachTPCNcls", "hTrackBachTPCNcls", HistType::kTH1F, {{152, 0, 152, "# TPC clusters"}});
      registry.add("QA/Tracks/hTrackPosHasTPC", "hTrackPosHasTPC", HistType::kTH1F, {{2, -0.5, 1.5, "has TPC"}});
      registry.add("QA/Tracks/hTrackNegHasTPC", "hTrackNegHasTPC", HistType::kTH1F, {{2, -0.5, 1.5, "has TPC"}});
      registry.add("QA/Tracks/hTrackBachHasTPC", "hTrackBachHasTPC", HistType::kTH1F, {{2, -0.5, 1.5, "has TPC"}});
      registry.add("QA/Tracks/hTrackProtonTPCPID", "hTrackProtonTPCPID", HistType::kTH2F, {{100, -10.0f, 10.0f, "p/z (GeV/c)"}, {100, -10.0f, 10.0f, "TPC n#sigma"}});
      registry.add("QA/Tracks/hTrackPionTPCPID", "hTrackPionTPCPID", HistType::kTH2F, {{100, -10.0f, 10.0f, "p/z (GeV/c)"}, {100, -10.0f, 10.0f, "TPC n#sigma"}});
      registry.add("QA/Tracks/hTrackBachTPCPID", "hTrackBachTPCPID", HistType::kTH2F, {{100, -10.0f, 10.0f, "p/z (GeV/c)"}, {100, -10.0f, 10.0f, "TPC n#sigma"}});
      registry.add("QA/Tracks/hTrackProtonPt", "hTrackProtonPt", HistType::kTH1F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}});
      registry.add("QA/Tracks/hTrackPionPt", "hTrackPionPt", HistType::kTH1F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}});
      registry.add("QA/Tracks/hTrackBachPt", "hTrackBachPt", HistType::kTH1F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}});
    }
    if (doEventQA) {
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

    // if (doprocessRun3 == true || doprocessRun3Reduced || doprocessRun3ReducedEM == true || doprocessRun3Reduced3bodyMixing == true || doprocessRun3Reduced3bodyMixingKFInfo == true) {
    //   auto hVtx3BodyCounter = registry.add<TH1>("hVtx3BodyCounter", "hVtx3BodyCounter", HistType::kTH1D, {{6, 0.0f, 6.0f}});
    //   hVtx3BodyCounter->GetXaxis()->SetBinLabel(1, "Total");
    //   hVtx3BodyCounter->GetXaxis()->SetBinLabel(2, "TPCNcls");
    //   hVtx3BodyCounter->GetXaxis()->SetBinLabel(3, "PIDCut");
    //   hVtx3BodyCounter->GetXaxis()->SetBinLabel(4, "HasSV");
    //   hVtx3BodyCounter->GetXaxis()->SetBinLabel(5, "DcaDau");
    //   hVtx3BodyCounter->GetXaxis()->SetBinLabel(6, "CosPA");
    //   registry.add("hBachelorTOFNSigmaDe", "", HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}});
    // }

    // if (doprocessRun3ReducedEM == true) {
    //   registry.add("hEventCount", "hEventCount", HistType::kTH2F, {dcaFitterEMSel.binsVtxZ, dcaFitterEMSel.binsMultiplicity});
    //   registry.add("hEventPairs", "hEventPairs", HistType::kTH2F, {dcaFitterEMSel.binsVtxZ, dcaFitterEMSel.binsMultiplicity});
    //   registry.add("hDecay3BodyPairsBeforeCut", "hDecay3BodyPairsBeforeCut", HistType::kTH2F, {dcaFitterEMSel.binsVtxZ, dcaFitterEMSel.binsMultiplicity});
    //   registry.add("hDecay3BodyPairsAfterCut", "hDecay3BodyPairsAfterCut", HistType::kTH2F, {dcaFitterEMSel.binsVtxZ, dcaFitterEMSel.binsMultiplicity});
    //   registry.add("hRadius0", "hRadius0", HistType::kTH1F, {{200, 0.0f, 20.0f, "Radius (cm)"}});
    //   registry.add("hRadius1", "hRadius1", HistType::kTH1F, {{200, 0.0f, 20.0f, "Radius (cm)"}});
    //   registry.add("hDeltaRadius", "hDeltaRadius", HistType::kTH1F, {{400, -20.0f, 20.0f, "#Delta Radius (cm)"}});
    //   registry.add("hPhi0", "hPhi0", HistType::kTH1F, {{360, -180.0f, 180.0f, "#phi (degree)"}});
    //   registry.add("hPhi1", "hPhi1", HistType::kTH1F, {{360, -180.0f, 180.0f, "#phi (degree)"}});
    //   registry.add("hDeltaPhi", "hDeltaPhi", HistType::kTH1F, {{360, -180.0f, 180.0f, "#Delta #phi (degree)"}});
    // }

    // if (doprocessRun3Reduced3bodyMixing == true || doprocessRun3Reduced3bodyMixingKFInfo == true) {
    //   registry.add("hDecay3BodyRadiusPhi", "hDecay3BodyRadiusPhi", HistType::kTH2F, {dcaFitterEMSel.bins3BodyRadius, dcaFitterEMSel.bins3BodyPhi});
    //   registry.add("hDecay3BodyPosZ", "hDecay3BodyPosZ", HistType::kTH1F, {dcaFitterEMSel.bins3BodyPosZ});
    //   auto h3bodyCombinationCounter = registry.add<TH1>("h3bodyCombinationCounter", "h3bodyCombinationCounter", HistType::kTH1D, {{4, 0.0f, 4.0f}});
    //   h3bodyCombinationCounter->GetXaxis()->SetBinLabel(1, "total");
    //   h3bodyCombinationCounter->GetXaxis()->SetBinLabel(2, "bach sign/ID");
    //   h3bodyCombinationCounter->GetXaxis()->SetBinLabel(3, "not same collision");
    //   h3bodyCombinationCounter->GetXaxis()->SetBinLabel(4, "collision VtxZ");
    // }

    // if (doprocessRun3ReducedEM == true || doprocessRun3Reduced3bodyMixing == true || doprocessRun3Reduced3bodyMixingKFInfo == true) {
    //   doUpdateGRPMagField = true;
    //   registry.add("h3bodyEMCutCounter", "h3bodyEMCutCounter", HistType::kTH1D, {{14, 0.0f, 14.0f}});
    // }

    // if (doprocessRun3withKFParticle == true || doprocessRun3withKFParticleStrangenessTracking == true || doprocessRun3withKFParticleReduced == true || doprocessRun3withKFParticleReducedEM == true || doprocessRun3withKFParticleReduced3bodyMixing == true) {
    //   auto hEventCounterZorro = registry.add<TH1>("Counters/hEventCounterZorro", "hEventCounterZorro", HistType::kTH1D, {{2, -0.5, 1.5}});
    //   hEventCounterZorro->GetXaxis()->SetBinLabel(1, "Zorro before evsel");
    //   hEventCounterZorro->GetXaxis()->SetBinLabel(2, "Zorro after evsel");
    //   auto hEventCounterKFParticle = registry.add<TH1>("Counters/hEventCounterKFParticle", "hEventCounterKFParticle", HistType::kTH1D, {{4, 0.0f, 4.0f}});
    //   hEventCounterKFParticle->GetXaxis()->SetBinLabel(1, "total");
    //   hEventCounterKFParticle->GetXaxis()->SetBinLabel(2, "sel8");
    //   hEventCounterKFParticle->GetXaxis()->SetBinLabel(3, "vertexZ");
    //   hEventCounterKFParticle->GetXaxis()->SetBinLabel(4, "has candidate");
    //   hEventCounterKFParticle->LabelsOption("v");
    //   auto hVtx3BodyCounterKFParticle = registry.add<TH1>("Counters/hVtx3BodyCounterKFParticle", "hVtx3BodyCounterKFParticle", HistType::kTH1D, {{21, 0.0f, 21.0f}});
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(1, "Total");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(2, "Charge");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(3, "Eta");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(4, "TPCNcls");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(5, "TPCRows");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(6, "TPCpid");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(7, "DCAxyPV");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(8, "DCAzPV");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(9, "V0MassConst");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(10, "HasSV");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(11, "DcaDau");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(12, "DCADauVtx");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(13, "DauPt");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(14, "Rapidity");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(15, "Pt");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(16, "Mass");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(17, "CosPA");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(18, "CosPAXY");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(19, "Chi2geo");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(20, "TopoConstr");
    //   hVtx3BodyCounterKFParticle->GetXaxis()->SetBinLabel(21, "Chi2topo");
    //   hVtx3BodyCounterKFParticle->LabelsOption("v");
    // }

    // if (doprocessRun3withKFParticleReducedEM == true) {
    //   registry.add("QA/EM/hPairCounterMixing", "hPairCounterMixing", HistType::kTH1D, {{1, 0.0f, 1.0f}});
    //   auto hCombinationCounterMixing = registry.add<TH1>("QA/EM/hCombinationCounterMixing", "hCombinationCounterMixing", HistType::kTH1D, {{3, 0.0f, 3.0f}});
    //   hCombinationCounterMixing->GetXaxis()->SetBinLabel(1, "total");
    //   hCombinationCounterMixing->GetXaxis()->SetBinLabel(2, "bach sign/ID");
    //   hCombinationCounterMixing->GetXaxis()->SetBinLabel(3, "radius, phi");
    //   hCombinationCounterMixing->LabelsOption("v");

    //   registry.add("QA/EM/hEventBinCounts", "hEventBinCounts", HistType::kTH2D, {{10, 0, 10, "bins VtxZ"}, {13, 0, 13, "bins mult"}});
    //   registry.add("QA/EM/h3bodyBinCounts", "h3bodyBinCounts", HistType::kTH2D, {{10, 0, 10, "bins VtxZ"}, {13, 0, 13, "bins mult"}});
    //   registry.add("QA/EM/hPairBinCounts", "hPairBinCounts", HistType::kTH2D, {{10, 0, 10, "bins VtxZ"}, {13, 0, 13, "bins mult"}});

    //   registry.add("QA/EM/hRadius1", "hRadius1", HistType::kTH1F, {{200, 0.0f, 20.0f, "Radius (cm)"}});
    //   registry.add("QA/EM/hRadius2", "hRadius2", HistType::kTH1F, {{200, 0.0f, 20.0f, "Radius (cm)"}});
    //   registry.add("QA/EM/hPhi1", "hPhi1", HistType::kTH1F, {{360, 0.0f, 360.0f, "#phi (degree)"}});
    //   registry.add("QA/EM/hPhi2", "hPhi2", HistType::kTH1F, {{360, 0.0f, 360.0f, "#phi (degree)"}});
    //   registry.add("QA/EM/hDeltaRadius", "hDeltaRadius", HistType::kTH1F, {{200, 0.0f, 10.0f, "#Delta Radius (cm)"}});
    //   registry.add("QA/EM/hDeltaPhi", "hDeltaPhi", HistType::kTH1F, {{360, 0.0f, 360.0f, "#Delta #phi (degree)"}});
    // }

    // if (doprocessRun3withKFParticleReduced3bodyMixing == true) {
    //   auto h3bodyCombinationCounter = registry.add<TH1>("QA/EM/h3bodyCombinationCounter", "h3bodyCombinationCounter", HistType::kTH1D, {{4, 0.0f, 4.0f}});
    //   h3bodyCombinationCounter->GetXaxis()->SetBinLabel(1, "total");
    //   h3bodyCombinationCounter->GetXaxis()->SetBinLabel(2, "not same collision");
    //   h3bodyCombinationCounter->GetXaxis()->SetBinLabel(3, "collision VtxZ");
    //   h3bodyCombinationCounter->GetXaxis()->SetBinLabel(4, "bach sign/ID");
    //   h3bodyCombinationCounter->LabelsOption("v");
    //   // registry.add("QA/EM/h3bodyBinCounts", "h3bodyBinCounts", HistType::kTH3D, {{16, 0, 16, "bins radius"}, {36, 0, 36, "bins phi"}, {12, 0, 12, "bins pos Z"}});
    //   registry.add("QA/EM/h3bodyBinCounts", "h3bodyBinCounts", HistType::kTH2D, {{16, 0, 16, "bins radius"}, {18, 0, 18, "bins phi"}});

    //   AxisSpec radiusAxis = {kfparticleConfigurations.bins3BodyRadius, "Radius (cm)"};
    //   AxisSpec phiAxis = {kfparticleConfigurations.bins3BodyPhi, "#phi (degree)"};
    //   AxisSpec posZAxis = {kfparticleConfigurations.bins3BodyPosZ, "position in z (cm)"};

    //   registry.add("QA/EM/hRadius", "hRadius", HistType::kTH1F, {radiusAxis});
    //   registry.add("QA/EM/hPhi", "hPhi", HistType::kTH1F, {phiAxis});
    //   registry.add("QA/EM/hPosZ", "hPosZ", HistType::kTH1F, {posZAxis});
    // }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    if (doSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerList);
      zorro.populateHistRegistry(registry, bc.runNumber());
    }

    auto timestamp = bc.timestamp();
    o2::parameters::GRPMagField* grpmag = 0x0;
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbConfigurations.grpmagPath, timestamp);
    if (!grpmag) {
      LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfigurations.grpmagPath << " of object GRPMagField for timestamp " << timestamp;
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);
    // Fetch magnetic field from ccdb for current collision
    auto d_bz = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp << " with magnetic field of " << magneticField << " kG";

    // set magnetic field value for DCA fitter
    fitterV0.setBz(d_bz);
    fitter3body.setBz(d_bz);
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
      helper.lut = lut;
    }

    // mark run as configured
    mRunNumber = bc.runNumber();

    // Initial TOF PID Paras, copied from PIDTOF.h
    tofPIDOpts.timestamp.value = bc.timestamp();
    ccdb->setTimestamp(tofPIDOpts.timestamp.value);
    // Not later than now objects
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    // TODO: implement the automatic pass name detection from metadata
    if (tofPIDOpts.passName.value == "") {
      tofPIDOpts.passName.value = "unanchored"; // temporary default
      LOG(warning) << "Passed autodetect mode for pass, not implemented yet, waiting for metadata. Taking '" << tofPIDOpts.passName.value << "'";
    }
    LOG(info) << "Using parameter collection, starting from pass '" << tofPIDOpts.passName.value << "'";

    const std::string fname = tofPIDOpts.paramFileName.value;
    if (!fname.empty()) { // Loading the parametrization from file
      LOG(info) << "Loading exp. sigma parametrization from file " << fname << ", using param: " << tofPIDOpts.parametrizationPath.value;
      if (1) {
        o2::tof::ParameterCollection paramCollection;
        paramCollection.loadParamFromFile(fname, tofPIDOpts.parametrizationPath.value);
        LOG(info) << "+++ Loaded parameter collection from file +++";
        if (!paramCollection.retrieveParameters(mRespParamsV2, tofPIDOpts.passName.value)) {
          if (tofPIDOpts.fatalOnPassNotAvailable) {
            LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", tofPIDOpts.passName.value.data());
          } else {
            LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", tofPIDOpts.passName.value.data());
          }
        } else {
          mRespParamsV2.setShiftParameters(paramCollection.getPars(tofPIDOpts.passName.value));
          mRespParamsV2.printShiftParameters();
        }
      } else {
        mRespParamsV2.loadParamFromFile(fname.data(), tofPIDOpts.parametrizationPath.value);
      }
    } else if (tofPIDOpts.loadResponseFromCCDB) { // Loading it from CCDB
      LOG(info) << "Loading exp. sigma parametrization from CCDB, using path: " << tofPIDOpts.parametrizationPath.value << " for timestamp " << timestamp.value;
      o2::tof::ParameterCollection* paramCollection = ccdb->getForTimeStamp<o2::tof::ParameterCollection>(tofPIDOpts.parametrizationPath.value, timestamp.value);
      paramCollection->print();
      if (!paramCollection->retrieveParameters(mRespParamsV2, tofPIDOpts.passName.value)) { // Attempt at loading the parameters with the pass defined
        if (tofPIDOpts.fatalOnPassNotAvailable) {
          LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", tofPIDOpts.passName.value.data());
        } else {
          LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", tofPIDOpts.passName.value.data());
        }
      } else { // Pass is available, load non standard parameters
        mRespParamsV2.setShiftParameters(paramCollection->getPars(tofPIDOpts.passName.value));
        mRespParamsV2.printShiftParameters();
      }
    }
    mRespParamsV2.print();
    if (tofPIDOpts.timeShiftCCDBPath.value != "") {
      if (tofPIDOpts.timeShiftCCDBPath.value.find(".root") != std::string::npos) {
        mRespParamsV2.setTimeShiftParameters(tofPIDOpts.timeShiftCCDBPath.value, "gmean_Pos", true);
        mRespParamsV2.setTimeShiftParameters(tofPIDOpts.timeShiftCCDBPath.value, "gmean_Neg", false);
      } else {
        mRespParamsV2.setTimeShiftParameters(ccdb->getForTimeStamp<TGraph>(Form("%s/pos", tofPIDOpts.timeShiftCCDBPath.value.c_str()), tofPIDOpts.timestamp.value), true);
        mRespParamsV2.setTimeShiftParameters(ccdb->getForTimeStamp<TGraph>(Form("%s/neg", tofPIDOpts.timeShiftCCDBPath.value.c_str()), tofPIDOpts.timestamp.value), false);
      }
    }

    bachelorTOFPID.SetParams(mRespParamsV2);
  }

  void initCCDBfromRunNumber(int runNumber)
  {
    // set magnetic field only when run number changes
    if (mRunNumber == runNumber) {
      LOG(debug) << "CCDB initialized for run " << mRunNumber;
      return;
    }
    mRunNumber = runNumber; // Update the last run number

    // Check if the CCDB data for this run is already cached
    if (ccdbCache.find(runNumber) != ccdbCache.end()) {
      LOG(debug) << "CCDB data already cached for run " << runNumber;
      d_bz = ccdbCache[runNumber];
      if (doUpdateGRPMagField == true) {
        o2::base::Propagator::initFieldFromGRP(grpMagCache[runNumber].get());
      }
    } else {
      std::shared_ptr<o2::parameters::GRPMagField> grpmag = std::make_shared<o2::parameters::GRPMagField>(*ccdb->getForRun<o2::parameters::GRPMagField>(ccdbConfigurations.grpmagPath, runNumber));
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfigurations.grpmagPath << " of object GRPMagField and " << ccdbConfigurations.grpPath << " of object GRPObject for run number " << runNumber;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag.get());
      // Fetch magnetic field from ccdb for current collision
      d_bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << "Retrieved GRP for run number " << runNumber << " with magnetic field of " << d_bz << " kZG";

      ccdbCache[runNumber] = d_bz;
      grpMagCache[runNumber] = grpmag;
    }

    // Set magnetic field for KF vertexing
#ifdef HomogeneousField
    KFParticle::SetField(d_bz);
#endif
    // Set field for DCAfitter
    fitterV0.setBz(d_bz);
    fitter3body.setBz(d_bz);

    mV0Hyps[0].set(o2::track::PID::Lambda, o2::track::PID::Proton, o2::track::PID::Pion, pidCutsLambda, d_bz);
    mV0Hyps[1].set(o2::track::PID::Lambda, o2::track::PID::Pion, o2::track::PID::Proton, pidCutsLambda, d_bz);

    if (useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }

    // cache magnetic field info
    ccdbCache[runNumber] = d_bz;
  }




  //------------------------------------------------------------------
  // event mixing
  template <class TCollision, class TTrack, typename TMixed3bodys, typename TBinningType>
  void doMixed3Body(TMixed3bodys decay3bodys, TBinningType binningType)
  {
    // Strictly upper index policy for decay3body objects binned by radius, phi
    for (const auto& [decay3body0, decay3body1] : selfCombinations(binningType, dcaFitterEMSel.nUseMixed, -1, decay3bodys, decay3bodys)) {
      auto tpos0 = decay3body0.template track0_as<TTrack>();
      auto tneg0 = decay3body0.template track1_as<TTrack>();
      auto tbach0 = decay3body0.template track2_as<TTrack>();
      auto tpos1 = decay3body1.template track0_as<TTrack>();
      auto tneg1 = decay3body1.template track1_as<TTrack>();
      auto tbach1 = decay3body1.template track2_as<TTrack>();

      registry.fill(HIST("h3bodyCombinationCounter"), 0.5);

      // ---------- selections ----------
      if ((tbach0.sign() > 0 && !(tbach1.sign() > 0)) || (tbach0.sign() < 0 && !(tbach1.sign() < 0)) || tbach0.globalIndex() == tbach1.globalIndex()) { // only combine if tbach1 has correct sign and is not same as tbach0
        continue;
      }
      registry.fill(HIST("h3bodyCombinationCounter"), 1.5);

      if (decay3body0.collisionId() == decay3body1.collisionId()) { // only combine if from different event
        continue;
      }
      registry.fill(HIST("h3bodyCombinationCounter"), 2.5);

      auto c0 = decay3body0.template collision_as<TCollision>();
      auto c1 = decay3body1.template collision_as<TCollision>();

      if (dcaFitterEMSel.selectPVPosZ3bodyMixing && std::abs(c0.posZ() - c1.posZ()) > dcaFitterEMSel.maxDeltaPVPosZ3bodyMixing) { // only combine if collision similar in PV posZ
        continue;
      }
      registry.fill(HIST("h3bodyCombinationCounter"), 3.5);

      initCCDBfromRunNumber(c0.runNumber());

      if (dcaFitterEMSel.cfgMix3BodyMethod == 0) { // mix bachelor (deuteron)
        fillVtxCand(c0, tpos0, tneg0, tbach1, -1, bachelorcharge, tbach1.tofNSigmaDe());
        fillVtxCand(c1, tpos1, tneg1, tbach0, -1, bachelorcharge, tbach0.tofNSigmaDe());
      } else if ((dcaFitterEMSel.cfgMix3BodyMethod == 1 && tbach0.sign() > 0) || (dcaFitterEMSel.cfgMix3BodyMethod == 2 && tbach0.sign() < 0)) { // mix piMinus or proton
        fillVtxCand(c0, tpos0, tneg1, tbach0, -1, bachelorcharge, tbach0.tofNSigmaDe());
        fillVtxCand(c1, tpos1, tneg0, tbach1, -1, bachelorcharge, tbach1.tofNSigmaDe());
      } else if ((dcaFitterEMSel.cfgMix3BodyMethod == 1 && tbach0.sign() < 0) || (dcaFitterEMSel.cfgMix3BodyMethod == 2 && tbach0.sign() > 0)) { // mix piPlus or anti-proton
        fillVtxCand(c0, tpos1, tneg0, tbach0, -1, bachelorcharge, tbach0.tofNSigmaDe());
        fillVtxCand(c1, tpos0, tneg1, tbach1, -1, bachelorcharge, tbach1.tofNSigmaDe());
      }

      VtxCandidates.clear();
    } // end decay3body combinations loop
  }
  //------------------------------------------------------------------
  // fill the StoredVtx3BodyDatas table
  void fillVtx3BodyTable(VtxCandidate const& candVtx)
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
  //-------------------- KFParticle reconstruction -------------------
  //------------------------------------------------------------------
  
  template <class TCollisionTo, typename TCollision, typename TTrack>
  double getTOFnSigma(TCollision const& collision, TTrack const& track, bool isEventMixing)
  {
    // TOF PID of deuteron (set motherhyp correctly)
    double tofNSigmaDeuteron = -999;
    if (track.has_collision() && track.hasTOF()) {
      if (isEventMixing) {
        tofNSigmaDeuteron = bachelorTOFPID.GetTOFNSigma(track, collision, collision);
      } else {
        auto originalcol = track.template collision_as<TCollisionTo>();
        tofNSigmaDeuteron = bachelorTOFPID.GetTOFNSigma(track, originalcol, collision);
      }
    }
    return tofNSigmaDeuteron;
  }


  // ______________________________________________________________
  // function to build decay3body candidates
  template <typename TCollision, typename T3Bodies, typename TTracks>
  void buildDecay3Bodies(TCollision const& collisions, T3Bodies const& decay3bodies, TTrack const& tracks)
  {
    if (!mEnabledTables[kStoredDecay3BodyCores]) {
      return; // don't do if no request for decay3bodies in place
    }

    int nDecay3Bodies = 0;
    // Loop over all decay3bodies in same time frame
    resgistry.fill(HIST("hInputStatistics"), kStoredDecay3BodyCores, decay3bodies.size());
    for (size_t i3body = 0; i3body < decay3bodies.size(); i3body++) {
      // Get tracks and generate candidate
      auto const& decay3body = decay3bodies[sorted_decay3body[i3body]]; /// TODO: FIX!!
      // if collisionId positive: get vertex, negative: origin
      float pvX = 0.0f, pvY = 0.0f, pvZ = 0.0f;
      if (decay3body.collisionId >= 0) {
        auto const& collision = collisions.rawIteratorAt(decay3body.collisionID);
        pvX = collision.posX();
        pvY = collision.posY();
        pvZ = collision.posZ();
      }

      // Aquire tracks
      auto const& trackPos = tracks.rawIteratorAt(decay3body.trackPosID);
      auto const& trackNeg = tracks.rawIteratorAt(decay3body.trackNegID);
      auto const& trackBach = tracks.rawIteratorAt(decay3body.traclBachID);

      // Calculate TOD nSigma
      /// TODO: Calculate TOD nSigma

      // generate analysis tables
      if (mEnabledTables[kDecay3BodyIndices]) {
        products.decay3bodyindices(decay3body.)
        registry.fill(HIST("hTableBuildingStatistics"), kDecay3BodyIndices);
      }
      if (mEnabledTables[kStoredDecay3BodyCores]) {
        products.decay3bodycores();
        registry.fill(HIST("hTableBuildingStatistics"), kStoredDecay3BodyCores);
      }
      if (mEnabledTables[kDecay3BodyCovs]) {
        products.decay3bodycovs();
        registry.fill(HIST("hTableBuildingStatistics"), kDecay3BodyCovs);
      }

      /// TODO:
      // ___________________________________________________________
      // MC handling part
      // prepare MC containers (not necessarily used)
      std::vector<mcCascinfo> mcCascinfos; // Decay3bodyMCCore information
      std::vector<bool> mcParticleIsReco;

      if constexpr (soa::is_table<TMCParticles>) {
        // do this if provided with a mcParticle table as well
        mcParticleIsReco.resize(mcParticles.size(), false);

        if ((mEnabledTables[kMcDEcay3BodyCores] || mEnabledTables[kMcDEcay3BodyLabels] || mEnabledTables[kMcDecay3BodyCollRefs])) {
          extractMonteCarloProperties();

          // generate label table (joinable with Decay3BodyCores)
          if (mEnabledTables[kMcDecay3BodyLabels]) {
            products.mcDecay3bodylabels();
            registry.fill(HIST("hTableBuildingStatistics"), kMcDecay3BodyLabels);
          }

          // mark mcParticle as reconstructed
          if (thisDecay3bodyInfo.label > -1) {
            mcParticleIsReco[thisDecay3bodyInfo.label] = true;
          }

          // generate MC cores table (joinable with Decay3BodyCores)
          if (mEnabledTables[kMcDecay3BodyCores]) {
            products.mcDecay3bodycores();
            registry.fill(HIST("hTableBuildingStatistics"), kMcDecay3BodyCores);
          }
          if (mEnabledTables[kMcDecay3BodyCollRefs]) {
            products.mcDecay3bodycollrefs(thisDecay3bodyInfo.mcCollision);
            registry.fill(HIST("hTableBuildingStatistics"), kMcDecay3BodyCollRefs);
          }
        } // enabled tables check
      } // constexpr requires mcParticles check
    } // decay3body loop

    /// TODO:
    // ____________________________________________________________
    // MC handling part
    if constexpr (soa::is_table<TMCParticles>) {
      if ((mEnabledTables[kMcDEcay3BodyCores] || mEnabledTables[kMcDEcay3BodyLabels] || mEnabledTables[kMcDecay3BodyCollRefs])) {

      } // enabled tables check
    } // constexpr requires mcParticles check
  }

  // ______________________________________________________________
  // function to get MC properties
  template <typename TTrack, typename TMMParticles>
  void extractMonteCarloProperties(TTrack const& trackPos, TTrack const& trackNeg, TTrack const& trackBach, TMMParticles const& mcParticles)
  {
    // encapsulates acquisition of MC properties from MC
    thisDecay3bodyInfo.pdgCode = -1, thisDecay3bodyInfo.pdgCodeMother = -1;
    thisDecay3bodyInfo.pdgCodePositive = -1, thisDecay3bodyInfo.pdgCodeNegative = -1;
    thisDecay3bodyInfo.pdgCodeBachelor = -1, thisDecay3bodyInfo.pdgCodeV0 = -1;
    thisDecay3bodyInfo.isPhysicalPrimary = false;
    thisDecay3bodyInfo.xyz[0] = -999.0f, thisDecay3bodyInfo.xyz[1] = -999.0f, thisDecay3bodyInfo.xyz[2] = -999.0f;
    thisDecay3bodyInfo.lxyz[0] = -999.0f, thisDecay3bodyInfo.lxyz[1] = -999.0f, thisDecay3bodyInfo.lxyz[2] = -999.0f;
    thisDecay3bodyInfo.posP[0] = -999.0f, thisDecay3bodyInfo.posP[1] = -999.0f, thisDecay3bodyInfo.posP[2] = -999.0f;
    thisDecay3bodyInfo.negP[0] = -999.0f, thisDecay3bodyInfo.negP[1] = -999.0f, thisDecay3bodyInfo.negP[2] = -999.0f;
    thisDecay3bodyInfo.bachP[0] = -999.0f, thisDecay3bodyInfo.bachP[1] = -999.0f, thisDecay3bodyInfo.bachP[2] = -999.0f;
    thisDecay3bodyInfo.momentum[0] = -999.0f, thisDecay3bodyInfo.momentum[1] = -999.0f, thisDecay3bodyInfo.momentum[2] = -999.0f;
    thisDecay3bodyInfo.label = -1, thisDecay3bodyInfo.motherLabel = -1;
    thisDecay3bodyInfo.mcParticlePositive = -1;
    thisDecay3bodyInfo.mcParticleNegative = -1;
    thisDecay3bodyInfo.mcParticleBachelor = -1;

    // association check
    if (trackPos.has_mcParticle() && trackNeg.has_mcParticle() && trackBach.has_mcParticle()) {
      auto lMCBachTrack = trackPos.template mcParticle_as<aod::McParticles>();
      auto lMCNegTrack = trackNeg.template mcParticle_as<aod::McParticles>();
      auto lMCPosTrack = trackBach.template mcParticle_as<aod::McParticles>();

      thisDecay3bodyInfo.mcParticlePositive = lMCPosTrack.globalIndex();
      thisDecay3bodyInfo.mcParticleNegative = lMCNegTrack.globalIndex();
      thisDecay3bodyInfo.mcParticleBachelor = lMCBachTrack.globalIndex();
      thisDecay3bodyInfo.pdgCodePositive = lMCPosTrack.pdgCode();
      thisDecay3bodyInfo.pdgCodeNegative = lMCNegTrack.pdgCode();
      thisDecay3bodyInfo.pdgCodeBachelor = lMCBachTrack.pdgCode();
      thisDecay3bodyInfo.posP[0] = lMCPosTrack.px();
      thisDecay3bodyInfo.posP[1] = lMCPosTrack.py();
      thisDecay3bodyInfo.posP[2] = lMCPosTrack.pz();
      thisDecay3bodyInfo.negP[0] = lMCNegTrack.px();
      thisDecay3bodyInfo.negP[1] = lMCNegTrack.py();
      thisDecay3bodyInfo.negP[2] = lMCNegTrack.pz();
      thisDecay3bodyInfo.bachP[0] = lMCBachTrack.px();
      thisDecay3bodyInfo.bachP[1] = lMCBachTrack.py();
      thisDecay3bodyInfo.bachP[2] = lMCBachTrack.pz();
    }

    // treat pi -> mu + antineutrino decays
    

  }

  //------------------------------------------------------------------
  // 3body candidate builder with KFParticle
  template <typename TCollision, typename TTrack>
  void buildVtx3BodyDataTableKFParticle(TCollision const& collision, TTrack const& trackPos, TTrack const& trackNeg, TTrack const& trackBach, int64_t decay3bodyID, int bachelorcharge, double tofNSigmaDeuteron)
  {
    gROOT->SetBatch(true);
    gRandom->SetSeed(42);

    // initialise KF primary vertex
    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle kfpv(kfpVertex);
    LOG(debug) << "Created KF PV.";

    // fill event QA histograms --> only for events with a decay3body!
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

    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxAll);

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
    // track IDs --> already checked in SVertexer!

    // track signs (pos, neg, bach) --> sanity check, should already be in SVertexer
    if (trackPos.sign() != +1 || trackNeg.sign() != -1) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxCharge);

    // track eta
    if (std::abs(trackPos.eta()) > kfparticleConfigurations.maxEta || std::abs(trackNeg.eta()) > kfparticleConfigurations.maxEta || std::abs(trackBach.eta()) > kfparticleConfigurations.maxEtaDeuteron) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxEta);

    // number of TPC clusters
    if (trackBach.tpcNClsFound() <= kfparticleConfigurations.mintpcNClsBach) {
      return;
    }
    if (isMatter && ((kfparticleConfigurations.useTPCforPion && trackNeg.tpcNClsFound() <= kfparticleConfigurations.mintpcNClsPion) || trackPos.tpcNClsFound() <= kfparticleConfigurations.mintpcNClsProton)) {
      return;
    } else if (!isMatter && ((kfparticleConfigurations.useTPCforPion && trackPos.tpcNClsFound() <= kfparticleConfigurations.mintpcNClsPion) || trackNeg.tpcNClsFound() <= kfparticleConfigurations.mintpcNClsProton)) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxTPCNcls);

    // number of TPC crossed rows
    if (trackBach.tpcNClsCrossedRows() <= kfparticleConfigurations.mintpcCrossedRows) {
      return;
    }
    if (isMatter && ((kfparticleConfigurations.useTPCforPion && trackNeg.tpcNClsCrossedRows() <= kfparticleConfigurations.mintpcCrossedRowsPion) || trackPos.tpcNClsCrossedRows() <= kfparticleConfigurations.mintpcCrossedRows)) {
      return;
    } else if (!isMatter && ((kfparticleConfigurations.useTPCforPion && trackPos.tpcNClsCrossedRows() <= kfparticleConfigurations.mintpcCrossedRowsPion) || trackNeg.tpcNClsCrossedRows() <= kfparticleConfigurations.mintpcCrossedRows)) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxTPCRows);

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
        return;
      }
    } else if (!isMatter) { // anti-hypertriton (anti-proton, pi+, deuteron)
      tpcNsigmaProton = trackNeg.tpcNSigmaPr();
      tpcNsigmaPion = trackPos.tpcNSigmaPi();
      dEdxProton = trackNeg.tpcSignal();
      dEdxPion = trackPos.tpcSignal();
      if (!selectTPCPID(trackNeg, trackPos, trackBach)) {
        return;
      }
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxTPCPID);
    LOG(debug) << "Basic track selections done.";

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
    // calculate 3D track DCA
    auto TrackPosDca = std::sqrt(TrackPosDcaXY * TrackPosDcaXY + TrackPosDcaZ * TrackPosDcaZ);
    auto TrackNegDca = std::sqrt(TrackNegDcaXY * TrackNegDcaXY + TrackNegDcaZ * TrackNegDcaZ);
    auto TrackBachDca = std::sqrt(TrackBachDcaXY * TrackBachDcaXY + TrackBachDcaZ * TrackBachDcaZ);
    // selection
    if (kfparticleConfigurations.doDCAPreSel && isMatter && (std::fabs(TrackNegDcaXY) <= kfparticleConfigurations.mindcaXYPionPV || std::fabs(TrackPosDcaXY) <= kfparticleConfigurations.mindcaXYProtonPV)) {
      return;
    } else if (kfparticleConfigurations.doDCAPreSel && !isMatter && (std::fabs(TrackPosDcaXY) <= kfparticleConfigurations.mindcaXYPionPV || std::fabs(TrackNegDcaXY) <= kfparticleConfigurations.mindcaXYProtonPV)) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxDCAxyPV);
    if (kfparticleConfigurations.doDCAPreSel && isMatter && (std::fabs(TrackNegDcaZ) <= kfparticleConfigurations.mindcaZPionPV || std::fabs(TrackPosDcaZ) <= kfparticleConfigurations.mindcaZProtonPV)) {
      return;
    } else if (kfparticleConfigurations.doDCAPreSel && !isMatter && (std::fabs(TrackPosDcaZ) <= kfparticleConfigurations.mindcaZPionPV || std::fabs(TrackNegDcaZ) <= kfparticleConfigurations.mindcaZProtonPV)) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxDCAzPV);

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
        return;
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
      return;
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
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxV0MassConst);

    // apply virtual V0 cuts used in SVertexer in case of 3body mixing with proton track
    if (kfparticleConfigurations.mixingType == 1 && kfparticleConfigurations.applySVertexerV0Cuts) {
      // V0 radius
      if (std::sqrt(KFV0.GetX() * KFV0.GetX() + KFV0.GetY() * KFV0.GetY()) <= 0.5) {
        return;
      }
      // pT
      if (KFV0.GetPt() <= 0.01) {
        return;
      }
      // pz/pT
      if (KFV0.GetPz() / KFV0.GetPt() >= 2) {
        return;
      }
      // cos(PA)
      if (cpaXYFromKF(KFV0, kfpv) <= 0.9 || cpaFromKF(KFV0, kfpv) <= 0.8) {
        return;
      }
    }

    // -------- STEP 3: fit three body vertex --------
    // Create KFParticle object from deuteron track
    KFParticle kfpDeuteron;
    kfpDeuteron = createKFParticleFromTrackParCov(trackParCovBach, trackBach.sign() * bachelorcharge, constants::physics::MassDeuteron);
    LOG(debug) << "KFParticle created from deuteron track.";

    // Construct vertex
    KFParticle KFHt;
    fit3bodyVertex(kfpProton, kfpPion, kfpDeuteron, KFHt);
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxhasSV);

    // -------- STEP 4: daughter selections after geometrical vertex fit --------
    // daughter DCAs with KF
    if ((kfpProton.GetDistanceFromParticle(kfpPion) >= kfparticleConfigurations.maxDcaProPi) || (kfpProton.GetDistanceFromParticle(kfpDeuteron) >= kfparticleConfigurations.maxDcaProDeu) || (kfpPion.GetDistanceFromParticle(kfpDeuteron) >= kfparticleConfigurations.maxDcaPiDe)) {
      return;
    }
    float DCAvtxDaughters3D = kfpProton.GetDistanceFromParticle(kfpPion) + kfpProton.GetDistanceFromParticle(kfpDeuteron) + kfpPion.GetDistanceFromParticle(kfpDeuteron);
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxDcaDau);
    LOG(debug) << "DCA selection after vertex fit applied.";

    // daughter DCAs to vertex
    if (kfpProton.GetDistanceFromVertexXY(KFHt) >= kfparticleConfigurations.maxDcaXYSVDau || kfpPion.GetDistanceFromVertexXY(KFHt) >= kfparticleConfigurations.maxDcaXYSVDau || kfpDeuteron.GetDistanceFromVertexXY(KFHt) >= kfparticleConfigurations.maxDcaXYSVDau) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxDcaDauVtx);
    LOG(debug) << "DCA to vertex selection after vertex fit applied.";

    // daughter pT
    if (kfpProton.GetPt() < kfparticleConfigurations.minPtProton || kfpProton.GetPt() > kfparticleConfigurations.maxPtProton || kfpPion.GetPt() < kfparticleConfigurations.minPtPion || kfpPion.GetPt() > kfparticleConfigurations.maxPtPion || kfpDeuteron.GetPt() < kfparticleConfigurations.minPtDeuteron || kfpDeuteron.GetPt() > kfparticleConfigurations.maxPtDeuteron) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxDauPt);
    LOG(debug) << "Daughter pT selection applied.";

    // -------- STEP 5: candidate selection and constraint after geometrical vertex fit --------
    // Rapidity
    float rapHt = RecoDecay::y(std::array{KFHt.GetPx(), KFHt.GetPy(), KFHt.GetPz()}, o2::constants::physics::MassHyperTriton);
    if (std::abs(rapHt) > kfparticleConfigurations.maxRapidityHt) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxRap);

    // Pt selection
    if (KFHt.GetPt() <= kfparticleConfigurations.minPtHt || KFHt.GetPt() >= kfparticleConfigurations.maxPtHt) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxPt);

    // Mass window
    float massHt, sigmaMassHt;
    KFHt.GetMass(massHt, sigmaMassHt);
    if (massHt <= kfparticleConfigurations.minMassHt || massHt >= kfparticleConfigurations.maxMassHt) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxMass);

    // cos(PA) to PV
    if (std::abs(cpaFromKF(KFHt, kfpv)) <= kfparticleConfigurations.minCosPA) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxCosPA);

    // cos(PA) xy to PV
    if (std::abs(cpaXYFromKF(KFHt, kfpv)) <= kfparticleConfigurations.minCosPAxy) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxCosPAXY);

    // chi2 geometrical
    float chi2geoNDF = KFHt.GetChi2() / KFHt.GetNDF();
    if (kfparticleConfigurations.applyTopoSel && chi2geoNDF >= kfparticleConfigurations.maxChi2geo) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxChi2geo);
    LOG(debug) << "Basic selections after vertex fit done.";

    // ctau before topo constraint
    if (KFHt.GetLifeTime() > kfparticleConfigurations.maxctauHt) {
      return;
    }

    // Set vertex constraint and topological selection
    KFParticle KFHtPV = KFHt;
    try {
      KFHtPV.SetProductionVertex(kfpv);
    } catch (std::runtime_error& e) {
      LOG(error) << "Exception caught KFParticle process call: Topological constraint failed";
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxTopoConstr); // to check if topo constraint fails
    // get topological chi2
    float chi2topoNDF = KFHtPV.GetChi2() / KFHtPV.GetNDF();
    KFHtPV.TransportToDecayVertex();
    if (kfparticleConfigurations.applyTopoSel && chi2topoNDF >= kfparticleConfigurations.maxChi2topo) {
      return;
    }
    registry.fill(HIST("Counters/hVtx3BodyCounterKFParticle"), kKfVtxChi2topo);

    // -------- STEP 6: collect and fill candidate info --------
    // get cluster size of strangeness tracked 3bodies
    float trackedClSize;
    if (decay3bodyID == -1) {
      trackedClSize = 0;
    } else {
      trackedClSize = !fTrackedClSizeVector.empty() ? fTrackedClSizeVector[decay3bodyID] : 0;
    }

    

    //------------------------------------------------------------------
    // table filling
    fillCandidateTable(candidate);
    LOG(debug) << "Table filled.";

    // fill event counter hist (has selected candidate) --> only filled once per vertex
    registry.fill(HIST("Counters/hEventCounterKFParticle"), 3.5);
  } // end buildVtx3BodyDataTableKFParticle

  //------------------------------------------------------------------
  void processRun3(ColwithEvTimes const& collisions, aod::Decay3Bodys const& decay3bodys, TrackExtPIDIUwithEvTimes const&, aod::BCsWithTimestamps const&)
  {
    VtxCandidates.clear();

    registry.fill(HIST("hEventCounter"), 0.5, collisions.size());

    for (const auto& d3body : decay3bodys) {
      auto t0 = d3body.track0_as<TrackExtPIDIUwithEvTimes>();
      auto t1 = d3body.track1_as<TrackExtPIDIUwithEvTimes>();
      auto t2 = d3body.track2_as<TrackExtPIDIUwithEvTimes>();
      auto collision = d3body.collision_as<ColwithEvTimes>();
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      // Recalculate the TOF PID
      double tofNSigmaBach = -999;
      if (t2.has_collision() && t2.hasTOF()) {
        auto originalcol = t2.template collision_as<ColwithEvTimes>();
        tofNSigmaBach = bachelorTOFPID.GetTOFNSigma(t2, originalcol, collision);
      }

      fillVtxCand(collision, t0, t1, t2, d3body.globalIndex(), bachelorcharge, tofNSigmaBach);
    }
  }
  PROCESS_SWITCH(decay3bodyBuilder, processRun3, "Produce DCA fitter decay3body tables", true);

  //------------------------------------------------------------------
  void processRun3Reduced(aod::RedCollisions const& collisions, aod::RedDecay3Bodys const& decay3bodys, aod::RedIUTracks const&)
  {
    VtxCandidates.clear();

    registry.fill(HIST("hEventCounter"), 0.5, collisions.size());

    for (const auto& d3body : decay3bodys) {
      auto t0 = d3body.track0_as<aod::RedIUTracks>();
      auto t1 = d3body.track1_as<aod::RedIUTracks>();
      auto t2 = d3body.track2_as<aod::RedIUTracks>();
      auto collision = d3body.collision_as<aod::RedCollisions>();

      initCCDBfromRunNumber(collision.runNumber());
      fillVtxCand(collision, t0, t1, t2, d3body.globalIndex(), bachelorcharge, t2.tofNSigmaDe());
    }
  }
  PROCESS_SWITCH(decay3bodyBuilder, processRun3Reduced, "Produce DCA fitter decay3body tables with reduced data", false);


  void processRun3Reduced3bodyMixing(ReducedCollisionsMults const&, aod::RedIUTracks const&, soa::Join<aod::RedDecay3Bodys, aod::DCAFitterSVInfo> const& decay3bodys)
  {
    VtxCandidates.clear();

    auto xAxis = registry.get<TH2>(HIST("hDecay3BodyRadiusPhi"))->GetXaxis();
    auto yAxis = registry.get<TH2>(HIST("hDecay3BodyRadiusPhi"))->GetYaxis();

    for (const auto& decay3body : decay3bodys) {
      int bin_Radius = xAxis->FindBin(decay3body.svRadius());
      int bin_Phi = yAxis->FindBin(decay3body.momPhi());
      registry.fill(HIST("hDecay3BodyRadiusPhi"), xAxis->GetBinCenter(bin_Radius), yAxis->GetBinCenter(bin_Phi));
      registry.fill(HIST("hDecay3BodyPosZ"), decay3body.svPosZ());
    }

    Binning3BodyDCAFitter binningOnRadiusPhi{{dcaFitterEMSel.bins3BodyRadius, dcaFitterEMSel.bins3BodyPhiDegree}, true};
    doMixed3Body<ReducedCollisionsMults, aod::RedIUTracks>(decay3bodys, binningOnRadiusPhi);
  }
  PROCESS_SWITCH(decay3bodyBuilder, processRun3Reduced3bodyMixing, "Produce mixing background directly from mixed decay3bodys based on DCAFitter Info", false);

  void processRun3Reduced3bodyMixingKFInfo(ReducedCollisionsMults const&, aod::RedIUTracks const&, soa::Join<aod::RedDecay3Bodys, aod::Red3BodyInfo> const& decay3bodys)
  {
    VtxCandidates.clear();

    auto xAxis = registry.get<TH2>(HIST("hDecay3BodyRadiusPhi"))->GetXaxis();
    auto yAxis = registry.get<TH2>(HIST("hDecay3BodyRadiusPhi"))->GetYaxis();

    for (const auto& decay3body : decay3bodys) {
      int bin_Radius = xAxis->FindBin(decay3body.radius());
      int bin_Phi = yAxis->FindBin(decay3body.phi());
      registry.fill(HIST("hDecay3BodyRadiusPhi"), xAxis->GetBinCenter(bin_Radius), yAxis->GetBinCenter(bin_Phi));
      registry.fill(HIST("hDecay3BodyPosZ"), decay3body.posz());
    }

    Binning3BodyKFInfo binningOnRadiusPhi{{dcaFitterEMSel.bins3BodyRadius, dcaFitterEMSel.bins3BodyPhi}, true};
    doMixed3Body<ReducedCollisionsMults, aod::RedIUTracks>(decay3bodys, binningOnRadiusPhi);
  }
  PROCESS_SWITCH(decay3bodyBuilder, processRun3Reduced3bodyMixingKFInfo, "Produce mixing background directly from mixed decay3bodys based on KF Info", false);

  //------------------------------------------------------------------
  void processRun3withKFParticle(ColwithEvTimes const& collisions, TrackExtPIDIUwithEvTimes const&, aod::Decay3Bodys const& decay3bodys, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      LOG(debug) << "CCDB initialised.";

      // Zorro event counting
      bool isZorroSelected = false;
      if (kfparticleConfigurations.cfgSkimmedProcessing) {
        isZorroSelected = zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC());
        if (isZorroSelected) {
          registry.fill(HIST("Counters/hEventCounterZorro"), 0.);
        } else {
          if (kfparticleConfigurations.cfgOnlyKeepInterestedTrigger) {
            continue;
          }
        }
      }

      // event selection
      registry.fill(HIST("Counters/hEventCounterKFParticle"), 0.5);
      if (kfparticleConfigurations.doSel8selection && !collision.sel8()) {
        continue;
      }
      registry.fill(HIST("Counters/hEventCounterKFParticle"), 1.5);
      if (kfparticleConfigurations.doPosZselection && (collision.posZ() >= 10.0f || collision.posZ() <= -10.0f)) {
        continue;
      }
      registry.fill(HIST("Counters/hEventCounterKFParticle"), 2.5);
      registry.fill(HIST("QA/Event/hAllSelEventsVtxZ"), collision.posZ());

      if (isZorroSelected) {
        registry.fill(HIST("Counters/hEventCounterZorro"), 1.);
      }

      // slice Decay3Body table by collision
      const uint64_t collIdx = collision.globalIndex();
      auto Decay3BodyTable_thisCollision = decay3bodys.sliceBy(perCollision, collIdx);
      for (auto& vtx3body : Decay3BodyTable_thisCollision) {
        auto trackPos = vtx3body.template track0_as<TrackExtPIDIUwithEvTimes>();
        auto trackNeg = vtx3body.template track1_as<TrackExtPIDIUwithEvTimes>();
        auto trackBach = vtx3body.template track2_as<TrackExtPIDIUwithEvTimes>();
        buildVtx3BodyDataTableKFParticle(collision, trackPos, trackNeg, trackBach, vtx3body.globalIndex(), bachelorcharge, getTOFnSigma<ColwithEvTimes>(collision, trackBach, false /*isEventMixing*/));
        LOG(debug) << "End of processKFParticle.";
      }
    }
  }
  PROCESS_SWITCH(decay3bodyBuilder, processRun3withKFParticle, "Produce KFParticle decay3body tables", false);

  void processRun3withKFParticleStrangenessTracking(ColwithEvTimes const& collisions, TrackExtPIDIUwithEvTimes const& tracks, aod::Decay3Bodys const& decay3bodys, aod::Tracked3Bodys const& tracked3bodys, aod::BCsWithTimestamps const& bcs)
  {
    fTrackedClSizeVector.clear();
    fTrackedClSizeVector.resize(decay3bodys.size(), 0);
    for (const auto& tvtx3body : tracked3bodys) {
      fTrackedClSizeVector[tvtx3body.decay3BodyId()] = tvtx3body.itsClsSize();
    }
    processRun3withKFParticle(collisions, tracks, decay3bodys, bcs);
  }
  PROCESS_SWITCH(decay3bodyBuilder, processRun3withKFParticleStrangenessTracking, "Produce KFParticle strangeness tracked decay3body tables", false);

  void processRun3withKFParticleReduced(aod::RedCollisions const& collisions, aod::RedIUTracks const&, aod::RedDecay3Bodys const& decay3bodys)
  {
    int lastRunNumber = -1;

    for (const auto& collision : collisions) {
      // set magnetic field only when run number changes
      if (collision.runNumber() != lastRunNumber) {
        initCCDBfromRunNumber(collision.runNumber());
        lastRunNumber = collision.runNumber(); // Update the last run number
        LOG(debug) << "CCDB initialized for run " << lastRunNumber;
      }

      // event selection
      registry.fill(HIST("Counters/hEventCounterKFParticle"), 2.5);
      registry.fill(HIST("QA/Event/hAllSelEventsVtxZ"), collision.posZ());

      // slice Decay3Body table by collision
      const uint64_t collIdx = collision.globalIndex();
      auto Decay3BodyTable_thisCollision = decay3bodys.sliceBy(perReducedCollision, collIdx);
      for (auto& vtx3body : Decay3BodyTable_thisCollision) {
        auto trackPos = vtx3body.template track0_as<aod::RedIUTracks>();
        auto trackNeg = vtx3body.template track1_as<aod::RedIUTracks>();
        auto trackBach = vtx3body.template track2_as<aod::RedIUTracks>();
        buildVtx3BodyDataTableKFParticle(collision, trackPos, trackNeg, trackBach, vtx3body.globalIndex(), bachelorcharge, trackBach.tofNSigmaDe());
      }
      LOG(debug) << "End of processKFParticleDerived.";
    }
  }
  PROCESS_SWITCH(decay3bodyBuilder, processRun3withKFParticleReduced, "Produce KFParticle decay3body tables from derived decay3body data", false);

  void processRun3withKFParticleReduced3bodyMixing(ReducedCollisionsMults const&, aod::RedIUTracks const&, soa::Join<aod::RedDecay3Bodys, aod::Red3BodyInfo> const& decay3bodys)
  {
    // Define a 2D array to count 3bodies per bin (radius, phi, posZ)
    std::vector<std::vector<int>> bin3bodyCounts(16, std::vector<int>(36, 0));

    // Function to find bin index (returns -1 if out of range)
    auto findBin = [](float value, const std::vector<float>& binEdges) -> int {
      for (size_t i = 0; i < binEdges.size() - 1; ++i) {
        if (value > binEdges[i] && value <= binEdges[i + 1]) {
          return i;
        }
      }
      return -1; // Out of range
    };

    int counter = 0;
    // Loop over all collisions to count them in bins
    for (auto& decay3body : decay3bodys) {
      counter++;
      float radius = decay3body.radius();
      float phi = decay3body.phi();
      float posZ = decay3body.posz();

      registry.fill(HIST("QA/EM/hRadius"), radius);
      registry.fill(HIST("QA/EM/hPhi"), phi);
      registry.fill(HIST("QA/EM/hPosZ"), posZ);

      // float degToRad = TMath::Pi()/180;

      // Determine bin indices
      int radiusBin = findBin(radius, {0.0f, 0.5f, 1.0f, 1.5f, 2.0f, 3.0f, 4.0f, 6.0f, 8.0f, 10.0f, 12.0f, 14.0f, 16.0f, 18.0f, 20.0f, 30.0f, 1000.0});
      int phiBin = findBin(phi, {-180.0f * TMath::Pi() / 180, -160.0f * TMath::Pi() / 180, -140.0f * TMath::Pi() / 180, -120.0f * TMath::Pi() / 180, -100.0f * TMath::Pi() / 180, -80.0f * TMath::Pi() / 180, -60.0f * TMath::Pi() / 180, -40.0f * TMath::Pi() / 180, -20.0f * TMath::Pi() / 180, 0.0f, 20.0f * TMath::Pi() / 180, 40.0f * TMath::Pi() / 180, 60.0f * TMath::Pi() / 180, 80.0f * TMath::Pi() / 180, 100.0f * TMath::Pi() / 180, 120.0f * TMath::Pi() / 180, 140.0f * TMath::Pi() / 180, 160.0f * TMath::Pi() / 180, 180.0f * TMath::Pi() / 180});
      if (radiusBin >= 0 && phiBin >= 0) {   // && posZBin >= 0) {
        bin3bodyCounts[radiusBin][phiBin]++; //[posZBin]++;
      }
    }
    LOG(info) << "3body counter: " << counter;

    // Print out the number of 3-body decays per bin
    LOG(info) << "3body count per bin (radius, phi, posZ):";
    for (size_t i = 0; i < bin3bodyCounts.size(); ++i) {
      for (size_t j = 0; j < bin3bodyCounts[i].size(); ++j) {
        LOG(info) << "Bin (" << i << ", " << j << "): " << bin3bodyCounts[i][j] << " 3bodies";
      }
    }
    // Fill 3D histogram with numbers per bin
    for (size_t i = 0; i < bin3bodyCounts.size(); ++i) {
      for (size_t j = 0; j < bin3bodyCounts[i].size(); ++j) {
        registry.fill(HIST("QA/EM/h3bodyBinCounts"), i, j, bin3bodyCounts[i][j]);
      }
    }
    LOG(info) << "Integral of h3bodyBinCounts: " << registry.get<TH2>(HIST("QA/EM/h3bodyBinCounts"))->Integral();

    Binning3Body binningOnRadPhi{{kfparticleConfigurations.bins3BodyRadius, kfparticleConfigurations.bins3BodyPhi}, true};

    // Strictly upper index policy for decay3body objects binned by radius, phi and z position
    for (auto& [decay3body1, decay3body2] : selfPairCombinations(binningOnRadPhi, kfparticleConfigurations.nEvtMixing, -1, decay3bodys)) {
      auto trackPos1 = decay3body1.template track0_as<aod::RedIUTracks>();
      auto trackNeg1 = decay3body1.template track1_as<aod::RedIUTracks>();
      auto trackBach1 = decay3body1.template track2_as<aod::RedIUTracks>();
      auto trackPos2 = decay3body2.template track0_as<aod::RedIUTracks>();
      auto trackNeg2 = decay3body2.template track1_as<aod::RedIUTracks>();
      auto trackBach2 = decay3body2.template track2_as<aod::RedIUTracks>();

      registry.fill(HIST("QA/EM/h3bodyCombinationCounter"), 0.5);

      // collision vertex selections
      auto collision1 = decay3body1.template collision_as<ReducedCollisionsMults>();
      auto collision2 = decay3body2.template collision_as<ReducedCollisionsMults>();
      initCCDBfromRunNumber(collision2.runNumber());
      initCCDBfromRunNumber(collision1.runNumber());

      if (decay3body1.collisionId() == decay3body2.collisionId()) { // only combine if from different event
        continue;
      }
      registry.fill(HIST("QA/EM/h3bodyCombinationCounter"), 1.5);
      if (kfparticleConfigurations.selectVtxZ3bodyMixing && std::abs(collision1.posZ() - collision2.posZ()) > kfparticleConfigurations.VtxZBin3bodyMixing) { // only combine if collision similar in VtxZ
        continue;
      }
      registry.fill(HIST("QA/EM/h3bodyCombinationCounter"), 2.5);

      // ---------- selections ----------
      if ((trackBach1.sign() > 0 && !(trackBach2.sign() > 0)) || (trackBach1.sign() < 0 && !(trackBach2.sign() < 0)) || trackBach1.globalIndex() == trackBach2.globalIndex()) { // only combine if trackBach2 has correct sign and is not same as trackBach1
        continue;
      }
      registry.fill(HIST("QA/EM/h3bodyCombinationCounter"), 3.5);

      // ---------- do candidate analysis ----------
      bool isMatter1 = false;
      if (trackBach1.sign() > 0) {
        isMatter1 = true;
      }
      if (kfparticleConfigurations.mixingType == 0) { // mix deuteron
        buildVtx3BodyDataTableKFParticle(collision1, trackPos1, trackNeg1, trackBach2, -1 /*vtx3bodyID*/, bachelorcharge, trackBach2.tofNSigmaDe());
        buildVtx3BodyDataTableKFParticle(collision2, trackPos2, trackNeg2, trackBach1, -1 /*vtx3bodyID*/, bachelorcharge, trackBach1.tofNSigmaDe());
      } else if (kfparticleConfigurations.mixingType == 1) { // mix proton
        if (isMatter1 == true) {
          buildVtx3BodyDataTableKFParticle(collision1, trackPos2, trackNeg1, trackBach1, -1 /*vtx3bodyID*/, bachelorcharge, trackBach1.tofNSigmaDe());
          buildVtx3BodyDataTableKFParticle(collision2, trackPos1, trackNeg2, trackBach2, -1 /*vtx3bodyID*/, bachelorcharge, trackBach2.tofNSigmaDe());
        } else if (isMatter1 == false) {
          buildVtx3BodyDataTableKFParticle(collision1, trackPos1, trackNeg2, trackBach1, -1 /*vtx3bodyID*/, bachelorcharge, trackBach1.tofNSigmaDe());
          buildVtx3BodyDataTableKFParticle(collision2, trackPos2, trackNeg1, trackBach2, -1 /*vtx3bodyID*/, bachelorcharge, trackBach2.tofNSigmaDe());
        }
      }
    } // end decay3body combinations loop
  }
  PROCESS_SWITCH(decay3bodyBuilder, processRun3withKFParticleReduced3bodyMixing, "Produce KFParticle mixed decay3body tables from derived decay3body data", false);
};





// build link from decay3body -> vtx3body
struct decay3bodyDataLinkBuilder {
  Produces<aod::Decay3BodyDataLink> VtxDataLink;

  void init(InitContext const&) {}

  template <typename TDecay3Bodys, typename TVtx3BodyDatas>
  void buildDecay3BodyLabel(TDecay3Bodys const& decay3bodytable, TVtx3BodyDatas const& vtxdatatable)
  {
    std::vector<int> lIndices;
    lIndices.reserve(decay3bodytable.size());
    for (int ii = 0; ii < decay3bodytable.size(); ii++)
      lIndices[ii] = -1;
    for (const auto& vtxdata : vtxdatatable) {
      if (vtxdata.decay3bodyId() != -1) {
        lIndices[vtxdata.decay3bodyId()] = vtxdata.globalIndex();
      }
    }
    for (int ii = 0; ii < decay3bodytable.size(); ii++) {
      VtxDataLink(lIndices[ii]);
    }
  }

  void processStandard(aod::Decay3Bodys const& decay3bodytable, aod::Vtx3BodyDatas const& vtxdatatable)
  {
    buildDecay3BodyLabel(decay3bodytable, vtxdatatable);
  }
  PROCESS_SWITCH(decay3bodyDataLinkBuilder, processStandard, "Produce label from decay3body to vtx3body", true);

  void processReduced(aod::RedDecay3Bodys const& decay3bodytable, aod::Vtx3BodyDatas const& vtxdatatable)
  {
    buildDecay3BodyLabel(decay3bodytable, vtxdatatable);
  }
  PROCESS_SWITCH(decay3bodyDataLinkBuilder, processReduced, "Produce label from reducedDecay3body to vtx3body", false);
};

struct kfdecay3bodyDataLinkBuilder {
  Produces<aod::KFDecay3BodyDataLink> kfvtxdataLink;

  void init(InitContext const&) {}

  template <typename TDecay3Bodys, typename TVtx3BodyDatas>
  void buildDataLink(TDecay3Bodys const& decay3bodytable, TVtx3BodyDatas const& vtxdatatable)
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
      kfvtxdataLink(lIndices[ii]);
    }
  }

  void processStandard(aod::Decay3Bodys const& decay3bodytable, aod::KFVtx3BodyDatas const& vtxdatatable)
  {
    buildDataLink(decay3bodytable, vtxdatatable); // build Decay3Body -> KFDecay3BodyData link table
  }
  PROCESS_SWITCH(kfdecay3bodyDataLinkBuilder, processStandard, "Build data link table.", true);

  void processReduced(aod::RedDecay3Bodys const& decay3bodytable, aod::KFVtx3BodyDatas const& vtxdatatable)
  {
    buildDataLink(decay3bodytable, vtxdatatable); // build ReducedDecay3Body -> KFDecay3BodyData link table
  }
  PROCESS_SWITCH(kfdecay3bodyDataLinkBuilder, processReduced, "Build data link table for reduced data.", false);
};

struct decay3bodyLabelBuilder {

  Produces<aod::McVtx3BodyLabels> vtxlabels;
  Produces<aod::McFullVtx3BodyLabels> vtxfulllabels;

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

  void processDoNotBuildLabels(aod::Decay3BodyDataLink const&) // is it possible to have none parameter?
  {
    // dummy process function - should not be required in the future
  };
  PROCESS_SWITCH(decay3bodyLabelBuilder, processDoNotBuildLabels, "Do not produce MC label tables", true);

  void processBuildLabels(aod::Decay3BodysLinked const& decay3bodys, aod::Vtx3BodyDatas const& vtx3bodydatas, MCLabeledTracksIU const&, aod::McParticles const&)
  {
    std::vector<int> lIndices;
    lIndices.reserve(vtx3bodydatas.size());
    for (int ii = 0; ii < vtx3bodydatas.size(); ii++) {
      lIndices[ii] = -1;
    }

    for (const auto& decay3body : decay3bodys) {

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

      for (const auto& lMother0 : lMCTrack0.mothers_as<aod::McParticles>()) {
        for (const auto& lMother1 : lMCTrack1.mothers_as<aod::McParticles>()) {
          for (const auto& lMother2 : lMCTrack2.mothers_as<aod::McParticles>()) {
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
        double hypertritonMCMass = RecoDecay::m(std::array{std::array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, std::array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, std::array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
        registry.fill(HIST("hLabelCounter"), 2.5);
        registry.fill(HIST("hHypertritonMCPt"), lPt);
        registry.fill(HIST("hHypertritonMCLifetime"), MClifetime);
        registry.fill(HIST("hHypertritonMCMass"), hypertritonMCMass);
      }
      if (lPDG == -1010010030 && lMCTrack0.pdgCode() == 211 && lMCTrack1.pdgCode() == -2212 && lMCTrack2.pdgCode() == -1000010020) {
        lLabel = lGlobalIndex;
        double antiHypertritonMCMass = RecoDecay::m(std::array{std::array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, std::array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, std::array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron});
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
    adaptAnalysisTask<kfdecay3bodyDataLinkBuilder>(cfgc),
    adaptAnalysisTask<decay3bodyInitializer>(cfgc),
  };
}
