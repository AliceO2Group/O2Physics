// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file trHeAnalysis.cxx
///
/// \brief triton and helion analysis on Run 3 pp data
///
/// \author Esther Bartsch <esther.bartsch@cern.ch>, Goethe University Frankfurt

#include "MetadataHelper.h"

#include "PWGLF/DataModel/LFNucleiTables.h"
#include "PWGLF/DataModel/LFPIDTOFGenericTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/Utils/pidTOFGeneric.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/BetheBlochAleph.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include "TRandom3.h"

#include <limits>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

using CollisionsFull =
  soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs,
            aod::CentFT0Ms, aod::CentFV0As>;
using CollisionsFullMC = soa::Join<CollisionsFull, aod::McCollisionLabels>;
using TracksFull =
  soa::Join<aod::TracksIU, aod::TracksCovIU, o2::aod::TracksDCA,
            aod::TracksExtra, aod::TrackSelectionExtension, aod::TracksDCA,
            aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta,
            o2::aod::EvTimeTOFFT0ForTrack>;
using TracksFullMC = soa::Join<TracksFull, aod::McTrackLabels>;
using TracksFullPid =
  soa::Join<TracksFull, aod::pidTPCFullPr, aod::pidTPCFullDe,
            aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl>;
using TracksFullPidMC = soa::Join<TracksFullPid, aod::McTrackLabels>;
using TracksFullLfPid =
  soa::Join<TracksFull, aod::pidTPCLfFullPr, aod::pidTPCLfFullDe,
            aod::pidTPCLfFullTr, aod::pidTPCLfFullHe, aod::pidTPCLfFullAl>;
using TracksFullLfPidMC = soa::Join<TracksFullLfPid, aod::McTrackLabels>;

namespace o2::aod
{
namespace trhe_data
{
DECLARE_SOA_COLUMN(Species, species, uint8_t);
DECLARE_SOA_COLUMN(Charge, charge, int8_t);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Rigidity, rigidity, float);
DECLARE_SOA_COLUMN(DeDx, deDx, float);
DECLARE_SOA_COLUMN(NSigmaTpc, nSigmaTpc, float);
DECLARE_SOA_COLUMN(NSigmaIts, nSigmaIts, float);
DECLARE_SOA_COLUMN(TofMass2, tofMass2, float);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);
DECLARE_SOA_COLUMN(SigmaYX, sigmaYX, float);
DECLARE_SOA_COLUMN(SigmaXYZ, sigmaXYZ, float);
DECLARE_SOA_COLUMN(SigmaZ, sigmaZ, float);
DECLARE_SOA_COLUMN(NTpcCluster, nTpcCluster, uint8_t);
DECLARE_SOA_COLUMN(NItsCluster, nItsCluster, uint8_t);
DECLARE_SOA_COLUMN(TpcChi2NCl, tpcChi2NCl, float);
DECLARE_SOA_COLUMN(ItsChi2NCl, itsChi2NCl, float);
DECLARE_SOA_COLUMN(ItsClusterSize, itsClusterSize, float);
DECLARE_SOA_COLUMN(DetectorMap, detectorMap, uint8_t);
DECLARE_SOA_COLUMN(Centrality, centrality, int);
DECLARE_SOA_COLUMN(Occupancy, occupancy, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(McTrue, mcTrue, bool);
DECLARE_SOA_COLUMN(McCollTrue, mcCollTrue, bool);
DECLARE_SOA_COLUMN(IsPhysPrimary, isPhysPrimary, bool);
DECLARE_SOA_COLUMN(IsReconstructed, isReconstructed, bool);
DECLARE_SOA_COLUMN(YGen, yGen, float);
DECLARE_SOA_COLUMN(ChargeGen, chargeGen, float);
DECLARE_SOA_COLUMN(PtGen, ptGen, float);
DECLARE_SOA_COLUMN(EtaGen, etaGen, float);
DECLARE_SOA_COLUMN(PhiGen, phiGen, float);
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);
} // namespace trhe_data
DECLARE_SOA_TABLE(TrHeData, "AOD", "trhe_data", trhe_data::Species,
                  trhe_data::Charge, trhe_data::Y, trhe_data::Pt,
                  trhe_data::Eta, trhe_data::Phi, trhe_data::Rigidity,
                  trhe_data::DeDx, trhe_data::NSigmaTpc, trhe_data::NSigmaIts,
                  trhe_data::TofMass2, trhe_data::DcaXY, trhe_data::DcaZ,
                  trhe_data::SigmaYX, trhe_data::SigmaXYZ, trhe_data::SigmaZ,
                  trhe_data::NTpcCluster, trhe_data::NItsCluster,
                  trhe_data::TpcChi2NCl, trhe_data::ItsChi2NCl,
                  trhe_data::ItsClusterSize, trhe_data::DetectorMap,
                  trhe_data::Centrality, trhe_data::Occupancy,
                  trhe_data::RunNumber);
DECLARE_SOA_TABLE(TrHeMcRec, "AOD", "trhe_mc_rec", trhe_data::Species,
                  trhe_data::Charge, trhe_data::Y, trhe_data::Pt,
                  trhe_data::Eta, trhe_data::Phi, trhe_data::Rigidity,
                  trhe_data::DeDx, trhe_data::NSigmaTpc, trhe_data::NSigmaIts,
                  trhe_data::TofMass2, trhe_data::DcaXY, trhe_data::DcaZ,
                  trhe_data::SigmaYX, trhe_data::SigmaXYZ, trhe_data::SigmaZ,
                  trhe_data::NTpcCluster, trhe_data::NItsCluster,
                  trhe_data::TpcChi2NCl, trhe_data::ItsChi2NCl,
                  trhe_data::ItsClusterSize, trhe_data::DetectorMap,
                  trhe_data::Centrality, trhe_data::Occupancy,
                  trhe_data::RunNumber, trhe_data::McTrue,
                  trhe_data::IsPhysPrimary, trhe_data::PdgCode, trhe_data::McCollTrue);
DECLARE_SOA_TABLE(TrHeMcGen, "AOD", "trhe_mc_gen", trhe_data::Species,
                  trhe_data::ChargeGen, trhe_data::YGen, trhe_data::PtGen,
                  trhe_data::EtaGen, trhe_data::PhiGen,
                  trhe_data::IsPhysPrimary, trhe_data::PdgCode);
DECLARE_SOA_TABLE(TrHeMc, "AOD", "trhe_mc", trhe_data::Species, trhe_data::YGen,
                  trhe_data::PtGen, trhe_data::EtaGen, trhe_data::PhiGen,
                  trhe_data::Charge, trhe_data::Y, trhe_data::Pt,
                  trhe_data::Eta, trhe_data::Phi, trhe_data::Rigidity,
                  trhe_data::DeDx, trhe_data::NSigmaTpc, trhe_data::NSigmaIts,
                  trhe_data::TofMass2, trhe_data::DcaXY, trhe_data::DcaZ,
                  trhe_data::SigmaYX, trhe_data::SigmaXYZ, trhe_data::SigmaZ,
                  trhe_data::NTpcCluster, trhe_data::NItsCluster,
                  trhe_data::TpcChi2NCl, trhe_data::ItsChi2NCl,
                  trhe_data::ItsClusterSize, trhe_data::DetectorMap,
                  trhe_data::Centrality, trhe_data::Occupancy,
                  trhe_data::RunNumber, trhe_data::IsPhysPrimary,
                  trhe_data::IsReconstructed, trhe_data::PdgCode, trhe_data::McCollTrue);
} // namespace o2::aod

namespace
{
static const int nParticles = 5;
enum Species { kProton,
               kDeuteron,
               kTriton,
               kHe3,
               kAlpha };
static const std::vector<std::string> particleNames{
  "proton", "deuteron", "triton", "helion", "alpha"};
static const std::vector<int> particlePdgCodes{
  2212, o2::constants::physics::kDeuteron, o2::constants::physics::kTriton,
  o2::constants::physics::kHelium3, o2::constants::physics::kAlpha};
static const std::vector<double> particleMasses{
  o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron,
  o2::constants::physics::MassTriton, o2::constants::physics::MassHelium3,
  o2::constants::physics::MassAlpha};
static const std::vector<int> particleCharge{1, 1, 1, 2, 2};
enum MassMethod { kMassFromTrack,
                  kMassFromBeta,
                  kMassFromTime };
const int nBetheParams = 8;
enum BBPAR { kP0,
             kP1,
             kP2,
             kP3,
             kP4,
             kResolution,
             kMip,
             kExp };
static const std::vector<std::string> betheBlochParNames{
  "p0", "p1", "p2", "p3", "p4", "resolution", "mip", "exp"};
constexpr float BetheBlochDefault[nParticles][nBetheParams]{
  {0.248753, 3.58634, 0.0167065, 2.29194, 0.774344, 0.07, 50., 1.}, // proton
  {0.248753, 3.58634, 0.0167065, 2.29194, 0.774344, 0.07, 50.,
   1.},                                                             // deuteron
  {0.248753, 3.58634, 0.0167065, 2.29194, 0.774344, 0.07, 50., 1.}, // triton
  {0.0274556, 18.3054, 3.99987e-05, 3.17219, 11.1775, 0.07, 50.,
   2.55}, // helion
  {0.0274556, 18.3054, 3.99987e-05, 3.17219, 11.1775, 0.07, 50.,
   2.55}}; // alpha

const int nTrkSettings = 12;
enum TPCPIDMETHOD {
  kSkipParticle = -1,
  kNone = 0,
  kParamBB = 1,
  kCentral = 2,
  kMC = 3
};
enum TRACKPIDSETTINGS {
  kPIDmethodTPC,
  kMinRigidity,
  kMaxRigidity,
  kMaxTPCnSigma,
  kMaxITSnSigma,
  kTOFrequiredabove,
  kMinTOFmass,
  kMaxTOFmass,
  kMinITSmeanClsSize,
  kMaxITSmeanClsSize,
  kMaxTPCchi2,
  kMaxITSchi2
};
static const std::vector<std::string> trackPIDsettingsNames{
  "PIDmethodTPC", "minRigidity", "maxRigidity", "maxTPCnSigma",
  "maxITSnSigma", "TOFrequiredabove", "minTOFmass2", "maxTOFmass2",
  "minITSclsSize", "maxITSclsSize", "maxTPCchi2", "maxITSchi2"};
constexpr float TrackPIDsettings[nParticles][nTrkSettings]{
  {-1, 0, 10, 3.0, -1., 100, 0, 10, 0., 100, 100, 10000},  // proton
  {-1, 0, 10, 3.0, -1., 100, 0, 10, 0., 100, 100, 10000},  // deuteron
  {2, 0, 10, 3.0, -1., 100, 0, 10, 0., 100, 100, 10000},   // triton
  {2, 0, 10, 3.0, -1., 100, 0, 10, 0., 100, 100, 10000},   // helion
  {-1, 0, 10, 3.0, -1., 100, 0, 10, 0., 100, 100, 10000}}; // alpha
} // end namespace

struct TrHeAnalysis {
  Produces<o2::aod::TrHeData> outData;
  Produces<o2::aod::TrHeMcRec> outMcRec;
  Produces<o2::aod::TrHeMcGen> outMcGen;
  Produces<o2::aod::TrHeMc> outMc;

  Preslice<TracksFull> perCollision = aod::track::collisionId;

  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};
  o2::aod::ITSResponse itsResponse;

  // configurables
  Configurable<float> cfgLowMultCut{"cfgLowMultCut", 0.0f,
                                    "Accepted multiplicity percentage lower limit"};
  Configurable<float> cfgHighMultCut{"cfgHighMultCut", 100.0f,
                                     "Accepted multiplicity percentage higher limit"};
  Configurable<bool> cfgRigidityCorrection{"cfgRigidityCorrection", false,
                                           "Enable Rigidity correction"};
  Configurable<float> cfgVtxCutZ{"cfgVtxCutZ", 10.0f,
                                 "Accepted z-vertex range"};
  Configurable<int> cfgMassMethod{"cfgMassMethod", 0,
                                  "0: Using built in 1: mass calculated with beta 2: mass calculated with "
                                  "the event time"};
  ConfigurableAxis binsVtxZ{"binsVtxZ", {100, -20.f, 20.f}, ""};
  ConfigurableAxis binsRigidity{"binsRigidity", {3000, -10.f, 10.f}, ""};
  ConfigurableAxis binsTpcSignal{"binsTpcSignal", {3000, 0.f, 3000.f}, ""};
  ConfigurableAxis binsPt{"binsPt", {20, 0.f, 10.f}, ""};
  struct : ConfigurableGroup {
    Configurable<float> cfgMinPt{"cfgMinPt", 0.f,
                                 "Min value of the pt selection"};
    Configurable<float> cfgMaxPt{"cfgMaxPt", 15.f,
                                 "Max value of the pt selection"};
    Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8f,
                                  "Value of the eta selection"};
    Configurable<float> cfgMaxRapidity{"cfgMaxRapidity", 1.0f,
                                       "Value of rapidity selection"};
    Configurable<int> cfgMinTpcCls{"cfgMinTpcCls", 0,
                                   "Minimum numbers of TPC clusters"};
    Configurable<int> cfgMinItsCls{"cfgMinItsCls", 0,
                                   "Minimum numbers of ITS clusters"};
    Configurable<float> cfgCutTpcXRows{"cfgCutTpcXRows", -1.f,
                                       "Minimum number of crossed TPC rows"};
    Configurable<float> cfgCutTpcCrRowToFindableCl{"cfgCutTpcCrRowToFindableCl",
                                                   0.8f, "Minimum ratio of crossed rows to findable cluster in TPC"};
    Configurable<bool> cfgCutTpcRefit{"cfgCutTpcRefit", 1, "TPC refit required"};
    Configurable<bool> cfgCutItsRefit{"cfgCutItsRefit", 1, "ITS refit required"};
    Configurable<float> cfgMinDCAXY{"cfgMinDCAXY", 0.f, "Minimum DCA to PV in XY"};
    Configurable<float> cfgMaxDCAXY{"cfgMaxDCAXY", 10000.f, "Maximum DCA to PV in Z"};
    Configurable<float> cfgMinDCAZ{"cfgMinDCAZ", 0.f, "Minimum DCA to PV in XY"};
    Configurable<float> cfgMaxDCAZ{"cfgMaxDCAZ", 10000.f, "Maximum DCA to PV in Z"};
    Configurable<int> cfgTrackSign{"cfgTrackSign", 0, "1: positive only, -1: negative only, 0: all tracks"};  
  } trackCuts;

  Configurable<LabeledArray<float>> cfgBetheBlochParams{"cfgBetheBlochParams",
                                                        {BetheBlochDefault[0], nParticles, nBetheParams, particleNames,
                                                         betheBlochParNames},
                                                        "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<float>> cfgTrackPIDsettings{"cfgTrackPIDsettings",
                                                        {TrackPIDsettings[0], nParticles, nTrkSettings, particleNames,
                                                         trackPIDsettingsNames},
                                                        "track PID criteria"};

  std::vector<std::shared_ptr<TH1>> histCuts;
  std::vector<int64_t> recoMcs;
  std::vector<int> goodEvents;
  //----------------------------------------------------------------------------
  void init(o2::framework::InitContext& context)
  {
    bool isMC = doprocessMC || doprocessMCCentralPid || doprocessMCLfPid;
    o2::aod::ITSResponse::setParameters(context, isMC);

    // define histogram axes
    const AxisSpec axisNev{3, 0., 3., "Number of events"};
    const AxisSpec axisCent{100, 0, 100, "centrality"};
    const AxisSpec axisOccupancy{5000, 0, 50000, "occupancy"};
    const AxisSpec axisVtxZ{binsVtxZ, "#it{z} (cm)"};
    const AxisSpec axisRigidity{binsRigidity, "#it{p/z} (GeV/#it{c})"};
    const AxisSpec axisdEdx{binsTpcSignal, "d#it{E}/d#it{x} (arb. u.)"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} GeV/#it{c}"};
    // create histograms
    histos.add("event/histVtxZ", "histVtxZ", kTH1F, {axisVtxZ});
    histos.add("event/histCentFT0A", "histCentFT0A", kTH1F, {axisCent});
    histos.add("event/histCentFT0C", "histCentFT0C", kTH1F, {axisCent});
    histos.add("event/histCentFT0M", "histCentFT0M", kTH1F, {axisCent});
    histos.add("event/histEvents", "histEvents", kTH2F,
               {axisCent, axisOccupancy});
    histos.add<TH1>("event/eventSelection", "eventSelection", HistType::kTH1D,
                    {{8, -0.5, 7.5}});
    auto h = histos.get<TH1>(HIST("event/eventSelection"));
    h->GetXaxis()->SetBinLabel(1, "Total");
    h->GetXaxis()->SetBinLabel(2, "TVX trigger cut");
    h->GetXaxis()->SetBinLabel(3, "TF border cut");
    h->GetXaxis()->SetBinLabel(4, "ITS ROF cut");
    h->GetXaxis()->SetBinLabel(5, "TVX + TF + ITS ROF");
    h->GetXaxis()->SetBinLabel(6, "Sel8 cut");
    h->GetXaxis()->SetBinLabel(7, "Z-vert Cut");
    h->GetXaxis()->SetBinLabel(7, "Centrality Cut");
    histos.add<TH2>("PID/histdEdx",
                    ";p_{TPC}/z (GeV/#it{c}); d#it{E}/d#it{x} (arb. u.)",
                    HistType::kTH2F, {axisRigidity, axisdEdx});
    histCuts.resize(nParticles);
    for (unsigned int species = 0; species < nParticles; species++) {
      int tpcMethod =
        static_cast<int>(cfgTrackPIDsettings->get(species, "PIDmethodTPC"));
      if (tpcMethod == kSkipParticle)
        continue;
      auto histName = Form("histCuts_%s", particleNames.at(species).data());
      histCuts.at(species) =
        histos.add<TH2>(Form("cuts/%s", histName), histName, HistType::kTH2F,
                        {{15, -0.5, 14.5}, axisPt});
      histCuts.at(species)->GetXaxis()->SetBinLabel(1, "TPC PID");
      histCuts.at(species)->GetXaxis()->SetBinLabel(2, "ITS PID");
      histCuts.at(species)->GetXaxis()->SetBinLabel(3, "TOF PID");
      histCuts.at(species)->GetXaxis()->SetBinLabel(4, "MeanItsClsSize");
      histCuts.at(species)->GetXaxis()->SetBinLabel(5, "eta");
      histCuts.at(species)->GetXaxis()->SetBinLabel(6, "nTPCcls");
      histCuts.at(species)->GetXaxis()->SetBinLabel(7, "nITScls");
      histCuts.at(species)->GetXaxis()->SetBinLabel(8, "ClsCrossedRows");
      histCuts.at(species)->GetXaxis()->SetBinLabel(9, "CrRowToFind");
      histCuts.at(species)->GetXaxis()->SetBinLabel(10, "TPC refit");
      histCuts.at(species)->GetXaxis()->SetBinLabel(11, "ITS refit");
      histCuts.at(species)->GetXaxis()->SetBinLabel(12, "TPC chi2");
      histCuts.at(species)->GetXaxis()->SetBinLabel(13, "ITS chi2");
      histCuts.at(species)->GetXaxis()->SetBinLabel(14, "DCA XY");
      histCuts.at(species)->GetXaxis()->SetBinLabel(15, "DCA Z");
    }
  }
  //----------------------------------------------------------------------------
  void processData(CollisionsFull const& collisions, TracksFull const& tracks,
                   aod::BCsWithTimestamps const& bcs)
  {
    fillTree<false, false>(collisions, tracks, true, bcs);
  }
  PROCESS_SWITCH(TrHeAnalysis, processData, "data analysis", false);
  //----------------------------------------------------------------------------
  void processDataCentralPid(CollisionsFull const& collisions,
                             TracksFullPid const& tracks,
                             aod::BCsWithTimestamps const& bcs)
  {
    fillTree<false, true>(collisions, tracks, true, bcs);
  }
  PROCESS_SWITCH(TrHeAnalysis, processDataCentralPid,
                 "data analysis using central PID", false);
  //----------------------------------------------------------------------------
  void processDataLfPid(CollisionsFull const& collisions,
                        TracksFullLfPid const& tracks,
                        aod::BCsWithTimestamps const& bcs)
  {
    fillTree<false, true>(collisions, tracks, true, bcs);
  }
  PROCESS_SWITCH(TrHeAnalysis, processDataLfPid, "data analysis using LF PID",
                 false);
  //----------------------------------------------------------------------------
  void processMC(CollisionsFullMC const& collisions, TracksFullMC const& tracks,
                 aod::BCsWithTimestamps const& bcs,
                 aod::McParticles const& particlesMC)
  {
    fillTree<true, false>(collisions, tracks, particlesMC, bcs);
  }
  PROCESS_SWITCH(TrHeAnalysis, processMC, "Monte Carlo analysis", true);
  //----------------------------------------------------------------------------
  void processMCCentralPid(CollisionsFullMC const& collisions,
                           TracksFullPidMC const& tracks,
                           aod::BCsWithTimestamps const& bcs,
                           aod::McParticles const& particlesMC)
  {
    fillTree<true, true>(collisions, tracks, particlesMC, bcs);
  }
  PROCESS_SWITCH(TrHeAnalysis, processMCCentralPid,
                 "Monte Carlo analysis using central PID", false);
  //----------------------------------------------------------------------------
  void processMCLfPid(CollisionsFullMC const& collisions,
                      TracksFullLfPidMC const& tracks,
                      aod::BCsWithTimestamps const& bcs,
                      aod::McParticles const& particlesMC)
  {
    fillTree<true, true>(collisions, tracks, particlesMC, bcs);
  }
  PROCESS_SWITCH(TrHeAnalysis, processMCLfPid,
                 "Monte Carlo analysis using LF PID", false);
  //----------------------------------------------------------------------------
  template <bool IsMC, bool UseExtPid, typename C, typename T, typename P>
  void fillTree(const C& collision, const T& tracks, const P& particles,
                aod::BCsWithTimestamps const& bcs)
  {
    recoMcs.clear();
    goodEvents.clear();
    // event loop
    for (const auto& collision : collision) {
      const auto& bc = bcs.rawIteratorAt(collision.bcId());
      // event selection
      histos.fill(HIST("event/eventSelection"), 0);
      if ((collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) &&
          (collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) &&
          (collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
        histos.fill(HIST("event/eventSelection"), 4);
      }
      if (!collision.sel8())
        continue;
      histos.fill(HIST("event/eventSelection"), 5);
      if (std::abs(collision.posZ()) > cfgVtxCutZ)
        continue;
      histos.fill(HIST("event/eventSelection"), 6);
      float centrality = collision.centFT0C();
      if (centrality < cfgLowMultCut || centrality > cfgHighMultCut)
        continue;
      histos.fill(HIST("event/eventSelection"), 7);
      float occupancy = collision.trackOccupancyInTimeRange();
      histos.fill(HIST("event/histVtxZ"), collision.posZ());
      histos.fill(HIST("event/histCentFT0A"), collision.centFT0A());
      histos.fill(HIST("event/histCentFT0C"), collision.centFT0C());
      histos.fill(HIST("event/histCentFT0M"), collision.centFT0M());
      histos.fill(HIST("event/histEvents"), collision.centFT0C(), occupancy);
      if constexpr (IsMC) {
        goodEvents.push_back(collision.mcCollisionId());
      }
      auto tracksByEvent =
        tracks.sliceBy(perCollision, collision.globalIndex());

      // track loop
      for (const auto& track : tracksByEvent) {
        if constexpr (IsMC) {
          if (!track.has_mcParticle())
            continue;
        }
        if (track.sign() * trackCuts.cfgTrackSign < 0)
          continue;
        float rigidity = getRigidity(track);
        histos.fill(HIST("PID/histdEdx"), track.sign() * rigidity,
                    track.tpcSignal());
        for (int species = 0; species < nParticles; species++) {
          int tpcMethod = static_cast<int>(
            cfgTrackPIDsettings->get(species, "PIDmethodTPC"));
          if (tpcMethod == kSkipParticle)
            continue;
          if (rigidity < cfgTrackPIDsettings->get(species, "minRigidity") ||
              rigidity > cfgTrackPIDsettings->get(species, "maxRigidity"))
            continue;
          float pt = particleCharge.at(species) * track.pt();
          if (pt < trackCuts.cfgMinPt || pt > trackCuts.cfgMaxPt)
            continue;

          float rapidity =
            RecoDecayPtEtaPhi::y(pt, track.eta(), particleMasses.at(species));
          if (std::abs(rapidity) > trackCuts.cfgMaxRapidity)
            continue;

          // TPC PID
          float tpcNsigma;
          switch (tpcMethod) {
            case kNone:
              tpcNsigma = 0;
              break;

            case kParamBB:
              tpcNsigma = getTpcNsigmaBB(track, species);
              break;

            case kCentral:
              if constexpr (UseExtPid) {
                tpcNsigma = getTpcNsigmaExt(track, species);
              } else {
                tpcNsigma = -999;
              }
              break;

            case kMC:
              if constexpr (IsMC) {
                tpcNsigma = (std::abs(track.mcParticle().pdgCode()) ==
                             particlePdgCodes[species])
                              ? 0
                              : -999;
              } else {
                tpcNsigma = -999;
              }
              break;

            default:
              tpcNsigma = -999;
          }

          if (std::abs(tpcNsigma) >
              cfgTrackPIDsettings->get(species, "maxTPCnSigma"))
            continue;
          histCuts.at(species)->Fill(0., pt);

          // ITS PID
          float itsNsigma = getItsNsigma(track, species);
          if (cfgTrackPIDsettings->get(species, "maxITSnSigma") > 0 &&
              std::abs(itsNsigma) >
                cfgTrackPIDsettings->get(species, "maxITSnSigma"))
            continue;
          histCuts.at(species)->Fill(1., pt);

          // TOF PID
          float tofMass2 = -1.f;
          if (track.hasTOF())
            tofMass2 = getMass2(track);
          if (pt > cfgTrackPIDsettings->get(species, "TOFrequiredabove") &&
              (tofMass2 < cfgTrackPIDsettings->get(species, "minTOFmass2") ||
               tofMass2 > cfgTrackPIDsettings->get(species, "maxTOFmass2")))
            continue;
          histCuts.at(species)->Fill(2., pt);

          // apply selection criteria
          const float itsMeanClsSize = getMeanItsClsSize(track);
          if (itsMeanClsSize <
                cfgTrackPIDsettings->get(species, "minITSclsSize") ||
              itsMeanClsSize >
                cfgTrackPIDsettings->get(species, "maxITSclsSize"))
            continue;
          histCuts.at(species)->Fill(3., pt);

          if (std::abs(track.eta()) > trackCuts.cfgMaxEta)
            continue;
          histCuts.at(species)->Fill(4., pt);

          if (track.tpcNClsFound() < trackCuts.cfgMinTpcCls)
            continue;
          histCuts.at(species)->Fill(5., pt);

          if (track.itsNCls() < trackCuts.cfgMinItsCls)
            continue;
          histCuts.at(species)->Fill(6., pt);

          if (track.tpcNClsCrossedRows() < trackCuts.cfgCutTpcXRows)
            continue;
          histCuts.at(species)->Fill(7., pt);

          if (track.tpcCrossedRowsOverFindableCls() <=
              trackCuts.cfgCutTpcCrRowToFindableCl)
            continue;
          histCuts.at(species)->Fill(8., pt);

          if (trackCuts.cfgCutTpcRefit && !track.passedTPCRefit())
            continue;
          histCuts.at(species)->Fill(9., pt);

          if (trackCuts.cfgCutItsRefit && !track.passedITSRefit())
            continue;
          histCuts.at(species)->Fill(10., pt);

          if (track.tpcChi2NCl() >
              cfgTrackPIDsettings->get(species, "maxTPCchi2"))
            continue;
          histCuts.at(species)->Fill(11., pt);

          if (track.itsChi2NCl() >
              cfgTrackPIDsettings->get(species, "maxITSchi2"))
            continue;
          histCuts.at(species)->Fill(12., pt);

          if (std::abs(track.dcaXY()) < trackCuts.cfgMinDCAXY ||
              std::abs(track.dcaXY()) > trackCuts.cfgMaxDCAXY)
            continue;
          histCuts.at(species)->Fill(13., pt);

          if (std::abs(track.dcaZ()) < trackCuts.cfgMinDCAZ ||
              std::abs(track.dcaZ()) > trackCuts.cfgMaxDCAZ)
            continue;
          histCuts.at(species)->Fill(14., pt);

          // write output tables
          if constexpr (IsMC) {
            const auto& mcPart = particles.rawIteratorAt(track.mcParticleId());
            const bool isMcTrue =
              mcPart.pdgCode() == particlePdgCodes.at(species) * track.sign();
            const bool mcTrueColl = mcPart.mcCollisionId() == collision.mcCollisionId();
            // mc reconstructed tree
            outMcRec(species, particleCharge.at(species) * track.sign(),
                     rapidity, pt, track.eta(), track.phi(), rigidity,
                     track.tpcSignal(), tpcNsigma, itsNsigma, tofMass2,
                     track.dcaXY(), track.dcaZ(), track.sigmaY(),
                     track.sigmaSnp(), track.sigmaZ(), track.tpcNClsFound(),
                     track.itsNCls(), track.tpcChi2NCl(), track.itsChi2NCl(),
                     itsMeanClsSize, track.detectorMap(), centrality,
                     collision.trackOccupancyInTimeRange(), bc.runNumber(),
                     isMcTrue, mcPart.isPhysicalPrimary(), mcPart.pdgCode(), mcTrueColl);
            // mc generated & reconstructed tree (mcTrue particles only)
            if (!isMcTrue)
              continue;
            outMc(species, mcPart.y(), mcPart.pt(), mcPart.eta(), mcPart.phi(),
                  particleCharge.at(species) * track.sign(), rapidity, pt,
                  track.eta(), track.phi(), rigidity, track.tpcSignal(),
                  tpcNsigma, itsNsigma, tofMass2, track.dcaXY(), track.dcaZ(),
                  track.sigmaY(), track.sigmaSnp(), track.sigmaZ(),
                  track.tpcNClsFound(), track.itsNCls(), track.tpcChi2NCl(),
                  track.itsChi2NCl(), itsMeanClsSize, track.detectorMap(),
                  centrality, collision.trackOccupancyInTimeRange(),
                  bc.runNumber(), mcPart.isPhysicalPrimary(), true, mcPart.pdgCode(), mcTrueColl);
            recoMcs.push_back(mcPart.globalIndex());
          } else { // data tree
            outData(species, particleCharge.at(species) * track.sign(),
                    rapidity, pt, track.eta(), track.phi(), rigidity,
                    track.tpcSignal(), tpcNsigma, itsNsigma, tofMass2,
                    track.dcaXY(), track.dcaZ(), track.sigmaY(),
                    track.sigmaSnp(), track.sigmaZ(), track.tpcNClsFound(),
                    track.itsNCls(), track.tpcChi2NCl(), track.itsChi2NCl(),
                    itsMeanClsSize, track.detectorMap(), centrality,
                    collision.trackOccupancyInTimeRange(), bc.runNumber());
          }
        } // end species loop
      } // end track loop
    } // end event loop

    // MC generated
    if constexpr (IsMC) {
      for (const auto& mcPart : particles) {
        if (!isInVector<int>(goodEvents, mcPart.mcCollisionId()))
          continue;
        for (int species = 0; species < nParticles; species++) {
          int tpcMethod = static_cast<int>(
            cfgTrackPIDsettings->get(species, "PIDmethodTPC"));
          if (tpcMethod == kSkipParticle)
            continue;
          if (std::abs(mcPart.pdgCode()) != particlePdgCodes.at(species))
            continue;
          if (std::abs(mcPart.y()) > trackCuts.cfgMaxRapidity)
            continue;
          float charge =
            particleCharge[species] * (mcPart.pdgCode() < 0 ? -1 : +1);
          // MC generated tree
          outMcGen(species, charge, mcPart.y(), mcPart.pt(), mcPart.eta(),
                   mcPart.phi(), mcPart.isPhysicalPrimary(), mcPart.pdgCode());
          // mc generated & reconstructed tree (non-reconstructed particles
          // only)
          if (isInVector<int64_t>(recoMcs, mcPart.globalIndex()))
            continue;
          outMc(species, mcPart.y(), mcPart.pt(), mcPart.eta(), mcPart.phi(),
                charge, -1.f, -1.f, -1.f, -1.f, -1.f, -1.f, -999.f, -999.f,
                -1.f, -1.f, -1.f, -1.f, -1.f, -1.f, -1, -1, -1.f, -1.f, -1.f, 0,
                0, 0, 0, mcPart.isPhysicalPrimary(), false, mcPart.pdgCode(), false);
        } // end species loop
      } // end mc particle loop
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  float getRigidity(T const& track)
  {
    if (!cfgRigidityCorrection)
      return track.tpcInnerParam();
    bool hePID = track.pidForTracking() == o2::track::PID::Helium3 ||
                 track.pidForTracking() == o2::track::PID::Alpha;
    return hePID ? track.tpcInnerParam() / 2 : track.tpcInnerParam();
  }
  //----------------------------------------------------------------------------
  template <class T>
  float getTpcNsigmaExt(T const& track, int species)
  {
    switch (species) {
      case Species::kProton:
        return track.tpcNSigmaPr();
      case Species::kDeuteron:
        return track.tpcNSigmaDe();
      case Species::kTriton:
        return track.tpcNSigmaTr();
      case Species::kHe3:
        return track.tpcNSigmaHe();
      case Species::kAlpha:
        return track.tpcNSigmaAl();
      default:
        return -999;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  float getTpcNsigmaBB(T const& track, int species)
  {
    const float mass = particleMasses.at(species);
    const float charge = particleCharge.at(species);
    const float mip = cfgBetheBlochParams->get(species, "mip");
    const float exp =
      std::pow(charge, cfgBetheBlochParams->get(species, "exp"));
    const float bg = charge * getRigidity(track) / mass;
    const float expBethe =
      mip * exp *
      o2::common::BetheBlochAleph(bg, cfgBetheBlochParams->get(species, "p0"),
                                  cfgBetheBlochParams->get(species, "p1"),
                                  cfgBetheBlochParams->get(species, "p2"),
                                  cfgBetheBlochParams->get(species, "p3"),
                                  cfgBetheBlochParams->get(species, "p4"));
    const float expSigma =
      expBethe * cfgBetheBlochParams->get(species, "resolution");
    return (track.tpcSignal() - expBethe) / expSigma;
  }
  //----------------------------------------------------------------------------
  template <class T>
  float getItsNsigma(T const& track, int species)
  {
    switch (species) {
      case Species::kProton:
        return itsResponse.nSigmaITS<o2::track::PID::Proton>(track);
      case Species::kDeuteron:
        return itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track);
      case Species::kTriton:
        return itsResponse.nSigmaITS<o2::track::PID::Triton>(track);
      case Species::kHe3:
        return itsResponse.nSigmaITS<o2::track::PID::Helium3>(track);
      case Species::kAlpha:
        return itsResponse.nSigmaITS<o2::track::PID::Alpha>(track);
      default:
        return -999;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  float getMass2(const T& track)
  {
    if (cfgMassMethod == kMassFromTrack) {
      const float& m = track.mass();
      return m * m;
    }
    if (cfgMassMethod == kMassFromBeta) {
      const float& beta = track.beta();
      const float& p = track.p();
      float gamma = 1.f / std::sqrt(1.f - beta * beta);
      float mass = p / std::sqrt(gamma * gamma - 1.f);
      return mass * mass;
    }
    if (cfgMassMethod == kMassFromTime) {
      const float p = track.p();
      const float& tofStartTime = track.evTimeForTrack();
      const float& tofTime = track.tofSignal();
      constexpr float CInCmPs = 2.99792458e-2f;
      const float& length = track.length();
      const float time = tofTime - tofStartTime;
      if (time > 0.f && length > 0.f) {
        const float beta = length / (CInCmPs * time);
        const float gamma = 1.f / std::sqrt(1.f - beta * beta);
        const float mass = p / std::sqrt(gamma * gamma - 1.f);
        return mass * mass;
      }
      return -1.f;
    }
    return -1.f;
  }
  //----------------------------------------------------------------------------
  template <class T>
  float getMeanItsClsSize(T const& track)
  {
    constexpr int NLayers = 8;
    constexpr int NBitsPerLayer = 4;
    constexpr int BitMask = (1 << NBitsPerLayer) - 1;

    int sum = 0, n = 0;
    for (int i = 0; i < NLayers; i++) {
      int clsSize = (track.itsClusterSizes() >> (NBitsPerLayer * i)) & BitMask;
      sum += clsSize;
      if (clsSize) {
        n++;
      }
    }
    return n > 0 ? static_cast<float>(sum) / n : 0.f;
  }
  //----------------------------------------------------------------------------
  template <typename T>
  bool isInVector(std::vector<T> const& vec, T const& val)
  {
    return std::find(vec.begin(), vec.end(), val) != vec.end();
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TrHeAnalysis>(cfgc)};
}
