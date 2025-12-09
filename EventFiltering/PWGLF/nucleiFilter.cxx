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
// O2 includes

#include "../filterTables.h"

#include "PWGLF/DataModel/LFPIDTOFGenericTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "PWGLF/Utils/pidTOFGeneric.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTOF/ParameterContainers.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector4D.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

o2::common::core::MetadataHelper metadataInfo;

namespace
{

static constexpr int nNuclei{3};
static constexpr int nCutsPID{5};
static constexpr std::array<float, nNuclei> masses{
  constants::physics::MassDeuteron, constants::physics::MassTriton,
  constants::physics::MassHelium3};
static constexpr std::array<int, nNuclei> charges{1, 1, 2};
static const std::vector<std::string> matterOrNot{"Matter", "Antimatter"};
static const std::vector<std::string> nucleiNames{"H2", "H3", "Helium"};
static const std::vector<std::string> hypernucleiNames{"H3L"}; // 3-body decay case
static const std::vector<std::string> columnsNames{o2::aod::filtering::H2::columnLabel(), o2::aod::filtering::He::columnLabel(), o2::aod::filtering::HeV0::columnLabel(), o2::aod::filtering::TritonFemto::columnLabel(), o2::aod::filtering::H3L3Body::columnLabel(), o2::aod::filtering::Tracked3Body::columnLabel(), o2::aod::filtering::ITSmildIonisation::columnLabel(), o2::aod::filtering::ITSextremeIonisation::columnLabel()};
static const std::vector<std::string> cutsNames{
  "TPCnSigmaMin", "TPCnSigmaMax", "TOFnSigmaMin", "TOFnSigmaMax", "TOFpidStartPt"};
constexpr double betheBlochDefault[nNuclei][6]{
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static constexpr float cutsPID[nNuclei][nCutsPID]{
  {-3.f, +400.f, -4.f, +400.f, 1.0f},  /*H2*/
  {-3.f, +3.f, -4.f, +4.f, 1.6f},      /*H3*/
  {-5.f, +400.f, -4.f, +4.f, 14000.f}, /*He3*/
};
constexpr double bbMomScalingDefault[nNuclei][2]{
  {1., 1.},
  {1., 1.},
  {1., 1.}};
constexpr double minTPCmom[nNuclei][2]{
  {0.8, 0.},
  {0., 0.},
  {0.8, 0.}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};

std::shared_ptr<TH2> h2TPCsignal[nNuclei];
std::shared_ptr<TH2> h2TPCnSigma[nNuclei];

} // namespace

struct nucleiFilter {

  Produces<aod::NucleiFilters> tags;

  // configurable for nuclei
  Configurable<float> cfgCutVertex{"cfgCutVertex", 12.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 1.f, "Eta range for tracks"};

  Configurable<float> cfgCutNclusITS{"cfgCutNclusITS", 2, "Minimum number of ITS clusters"};
  Configurable<int> cfgCutNclusExtremeIonisationITS{"cfgCutNclusExtremeIonisationITS", 5, "Minimum number of ITS clusters for the extreme ionisation trigger"};
  Configurable<float> cfgMomentumCutExtremeIonisation{"cfgMomentumCutExtremeIonisation", 1.2, "Minimum momentum for the extreme ionisation trigger"};
  Configurable<float> cfgCutClsSizeExtremeIonisation{"cfgCutClsSizeExtremeIonisation", 8, "Minimum average size of ITS clusters for the extreme ionisation trigger"};
  Configurable<float> cfgCutClsSizeMildIonisation{"cfgCutClsSizeMildIonisation", 5, "Minimum average size of ITS clusters for the mild ionisation trigger"};
  Configurable<float> cfgCutNclusTPC{"cfgCutNclusTPC", 80, "Minimum number of TPC clusters"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 3, "Max DCAxy"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 10, "Max DCAz"};
  Configurable<float> cfgCutKstar{"cfgCutKstar", 1.f, "Kstar cut for triton femto trigger"};
  Configurable<double> cfgCutCosPAheV0{"cfgCutCosPAheV0", 0.99, "CosPA cut for HeV0"};

  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], nNuclei, 6, nucleiNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgMomentumScalingBetheBloch{"cfgMomentumScalingBetheBloch", {bbMomScalingDefault[0], nNuclei, 2, nucleiNames, matterOrNot}, "TPC Bethe-Bloch momentum scaling for light nuclei"};
  Configurable<LabeledArray<double>> cfgMinTPCmom{"cfgMinTPCmom", {minTPCmom[0], nNuclei, 2, nucleiNames, matterOrNot}, "Minimum TPC p/Z for nuclei PID"};

  Configurable<LabeledArray<float>> cfgCutsPID{"nucleiCutsPID", {cutsPID[0], nNuclei, nCutsPID, nucleiNames, cutsNames}, "Nuclei PID selections"};
  Configurable<bool> cfgFixTPCinnerParam{"cfgFixTPCinnerParam", false, "Fix TPC inner param"};

  // variable/tool for hypertriton 3body decay
  int mRunNumber;
  float mBz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TrackSelection, aod::TracksDCA, aod::EvTimeTOFFT0ForTrack, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl>; // FIXME: positio has been changed
  o2::aod::pidtofgeneric::TofPidNewCollision<TrackCandidates::iterator> bachelorTOFPID;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::vertexing::DCAFitterN<2> fitter2body;
  o2::vertexing::DCAFitterN<3> fitter3body;
  // TOF response and input parameters
  o2::pid::tof::TOFResoParamsV3 mRespParamsV3;
  o2::aod::pidtofgeneric::TOFCalibConfig mTOFCalibConfig; // TOF Calib configuration
  // configurable for hypertriton 3body decay
  struct : ConfigurableGroup {
    Configurable<double> bFieldInput{"trgH3L3Body.mBz", -999, "bz field, -999 is automatic"};
    Configurable<float> minCosPA3body{"trgH3L3Body.minCosPA3body", 0.9995, "minCosPA3body"};
    Configurable<float> dcavtxdau{"trgH3L3Body.dcavtxdau", 0.15, "meen DCA among Daughters"};
    Configurable<float> dcapiontopv{"trgH3L3Body.dcapiontopv", 0.05, "DCA Pion To PV"};
    Configurable<float> tofPIDNSigmaMin{"trgH3L3Body.tofPIDNSigmaMin", -5, "tofPIDNSigmaMin"};
    Configurable<float> tofPIDNSigmaMax{"trgH3L3Body.tofPIDNSigmaMax", 5, "tofPIDNSigmaMax"};
    Configurable<float> tpcPIDNSigmaCut{"trgH3L3Body.tpcPIDNSigmaCut", 5, "tpcPIDNSigmaCut"};
    Configurable<float> lifetimecut{"trgH3L3Body.lifetimecut", 40., "lifetimecut"};
    Configurable<float> minDaughtersEta{"trgH3L3Body.minDaughtersEta", 1.f, "minDaughtersEta"};
    Configurable<float> minProtonPt{"trgH3L3Body.minProtonPt", 0.3, "minProtonPt"};
    Configurable<float> maxProtonPt{"trgH3L3Body.maxProtonPt", 5, "maxProtonPt"};
    Configurable<float> minPionPt{"trgH3L3Body.minPionPt", 0.1, "minPionPt"};
    Configurable<float> maxPionPt{"trgH3L3Body.maxPionPt", 1.2, "maxPionPt"};
    Configurable<float> minDeuteronPt{"trgH3L3Body.minDeuteronPt", 0.6, "minDeuteronPt"};
    Configurable<float> maxDeuteronPt{"trgH3L3Body.maxDeuteronPt", 10, "maxDeuteronPt"};
    Configurable<float> minDeuteronPUseTOF{"trgH3L3Body.minDeuteronPUseTOF", 999, "minDeuteronPt Enable TOF PID"};
    Configurable<float> h3LMassLowerlimit{"trgH3L3Body.h3LMassLowerlimit", 2.96, "Hypertriton mass lower limit"};
    Configurable<float> h3LMassUpperlimit{"trgH3L3Body.h3LMassUpperlimit", 3.04, "Hypertriton mass upper limit"};
    Configurable<float> minP3Body{"trgH3L3Body.minP3Body", 1.5, "min P3Body"};
    Configurable<int> mintpcNClsproton{"trgH3L3Body.mintpcNClsproton", 90, "min tpc Nclusters for proton"};
    Configurable<int> mintpcNClspion{"trgH3L3Body.mintpcNClspion", 70, "min tpc Nclusters for pion"};
    Configurable<int> mintpcNClsdeuteron{"trgH3L3Body.mintpcNClsdeuteron", 100, "min tpc Nclusters for deuteron"};

    Configurable<int> useMatCorrType{"trgH3L3Body.useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
    // CCDB options
    Configurable<std::string> ccdburl{"trgH3L3Body.ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"trgH3L3Body.grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"trgH3L3Body.grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"trgH3L3Body.lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"trgH3L3Body.geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  } trgH3L3Body;

  HistogramRegistry qaHists{"qaHists", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", ";;Number of filtered events", kNtriggers + 1, -0.5, static_cast<double>(kNtriggers) + 0.5)};

  void init(InitContext& initContext)
  {
    std::vector<double> ptBinning = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5.};

    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    qaHists.add("fCollZpos", "collision z position", HistType::kTH1F, {{600, -20., +20., "z position (cm)"}});
    qaHists.add("fTPCsignalAll", "Specific energy loss (before filter)", HistType::kTH2F, {{1200, -6, 6, "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    qaHists.add("fTPCsignal", "Specific energy loss", HistType::kTH2F, {{1200, -6, 6, "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    qaHists.add("fDeuTOFNsigma", "Deuteron TOF Nsigma distribution", HistType::kTH2F, {{1200, -6, 6, "#it{p} (GeV/#it{c})"}, {2000, -100, 100, "TOF n#sigma"}});
    qaHists.add("fBachDeuTOFNsigma", "Bachelor Deuteron TOF Nsigma distribution", HistType::kTH2F, {{1200, -6, 6, "#it{p} (GeV/#it{c})"}, {2000, -100, 100, "TOF n#sigma"}});
    qaHists.add("fH3LMassVsPt", "Hypertrion mass Vs pT", HistType::kTH2F, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {80, 2.96, 3.04, "Inv. Mass (GeV/c^{2})"}});
    qaHists.add("fH3LDcaVsPt", "DCA vs pT", HistType::kTH2F, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {100, 0, 0.05, "DCA (cm)"}});
    qaHists.add("fH3LCosPAVsPt", "CosPA vs pT", HistType::kTH2F, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {100, 0.999, 1.0, "CosPA"}});
    qaHists.add("fExtremeIonisationITS", "ITS clusters for extreme ionisation trigger", HistType::kTH3F, {{4, 3.5, 7.5, "Number of ITS clusters"}, {150, 0, 15, "Average cluster size in ITS x cos#lambda"}, {100, 0.1, 10, "#it{p} (GeV/#it{c})"}});

    for (int iN{0}; iN < nNuclei; ++iN) {
      h2TPCsignal[iN] = qaHists.add<TH2>(Form("fTPCsignal_%s", nucleiNames[iN].data()), "Specific energy loss", HistType::kTH2F, {{1200, -6, 6., "#it{p}/Z (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
      h2TPCnSigma[iN] = qaHists.add<TH2>(Form("fTPCcounts_%s", nucleiNames[iN].data()), "n-sigma TPC", HistType::kTH2F, {{100, -5, 5, "#it{p} /Z (GeV/#it{c})"}, {200, -10., +10., "n#sigma_{He} (a. u.)"}});
    }

    hProcessedEvents->GetXaxis()->SetBinLabel(1, "Processed events");
    for (uint32_t iS{0}; iS < columnsNames.size(); ++iS) {
      hProcessedEvents->GetXaxis()->SetBinLabel(iS + 2, columnsNames[iS].data());
    }

    // for fH3L3Body
    bachelorTOFPID.SetPidType(o2::track::PID::Deuteron);
    fitter3body.setPropagateToPCA(true);
    fitter3body.setMaxR(200.);
    fitter3body.setMinParamChange(1e-3);
    fitter3body.setMinRelChi2Change(0.9);
    fitter3body.setMaxDZIni(1e9);
    fitter3body.setMaxChi2(1e9);
    fitter3body.setUseAbsDCA(true);

    fitter2body.setPropagateToPCA(true);
    fitter2body.setMaxR(200.);
    fitter2body.setMinParamChange(1e-3);
    fitter2body.setMinRelChi2Change(0.9);
    fitter2body.setMaxDZIni(1e9);
    fitter2body.setMaxChi2(1e9);
    fitter2body.setUseAbsDCA(true);

    ccdb->setURL(trgH3L3Body.ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // Initialization of TOF PID parameters for fH3L3Body
    mTOFCalibConfig.metadataInfo = metadataInfo;
    mTOFCalibConfig.inheritFromBaseTask(initContext);
    mTOFCalibConfig.initSetup(mRespParamsV3, ccdb); // Getting the parametrization parameters
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (trgH3L3Body.bFieldInput > -990) {
      mBz = trgH3L3Body.bFieldInput;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(mBz) > 1e-5) {
        grpmag.setL3Current(30000.f / (mBz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(trgH3L3Body.grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      mBz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << mBz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(trgH3L3Body.grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << trgH3L3Body.grpmagPath << " of object GRPMagField and " << trgH3L3Body.grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      // mBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      mBz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << mBz << " kZG";
    }
    mRunNumber = bc.runNumber();
    // Set magnetic field value once known
    fitter2body.setBz(mBz);
    fitter3body.setBz(mBz);

    if (trgH3L3Body.useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized
      // (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }

    mTOFCalibConfig.processSetup(mRespParamsV3, ccdb, bc);
  }

  enum {
    kH2 = 0,
    kHe,
    kHeV0,
    kTritonFemto,
    kH3L3Body,
    kTracked3Body,
    kITSmildIonisation,
    kITSextremeIonisation,
    kNtriggers
  } TriggerType;
  // void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::Vtx3BodyDatas const& vtx3bodydatas, TrackCandidates const& tracks)
  using ColWithEvTime = soa::Join<aod::Collisions, aod::EvSels, aod::EvTimeTOFFT0>;
  void process(ColWithEvTime::iterator const& collision, aod::Decay3Bodys const& decay3bodys, TrackCandidates const& tracks, aod::AssignedTracked3Bodys const& tracked3Bodys, aod::V0s const& v0s, aod::BCsWithTimestamps const&)
  {
    // collision process loop
    std::array<bool, kNtriggers> keepEvent{false};
    //
    qaHists.fill(HIST("fCollZpos"), collision.posZ());
    hProcessedEvents->Fill(0);
    //
    if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      tags(keepEvent[kH2], keepEvent[kHe], keepEvent[kHeV0], keepEvent[kTritonFemto], keepEvent[kH3L3Body], keepEvent[kTracked3Body], keepEvent[kITSmildIonisation], keepEvent[kITSextremeIonisation]);
      return;
    }

    //
    const double bgScalings[nNuclei][2]{
      {charges[0] * cfgMomentumScalingBetheBloch->get(0u, 0u) / masses[0], charges[0] * cfgMomentumScalingBetheBloch->get(0u, 1u) / masses[0]},
      {charges[1] * cfgMomentumScalingBetheBloch->get(1u, 0u) / masses[1], charges[1] * cfgMomentumScalingBetheBloch->get(1u, 1u) / masses[1]},
      {charges[2] * cfgMomentumScalingBetheBloch->get(2u, 0u) / masses[2], charges[2] * cfgMomentumScalingBetheBloch->get(2u, 1u) / masses[2]}};

    constexpr int nucleusIndex[nNuclei]{kH2, -1, kHe}; /// remap for nuclei triggers
    std::vector<int> h3indices;
    std::vector<ROOT::Math::PtEtaPhiMVector> h3vectors;

    auto getNsigma = [&](const auto& track, int iN, int iC) {
      float fixTPCrigidity{(cfgFixTPCinnerParam && (track.pidForTracking() == track::PID::Helium3 || track.pidForTracking() == track::PID::Alpha)) ? 0.5f : 1.f};
      double expBethe{tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() * fixTPCrigidity * bgScalings[iN][iC]), cfgBetheBlochParams->get(iN, 0u), cfgBetheBlochParams->get(iN, 1u), cfgBetheBlochParams->get(iN, 2u), cfgBetheBlochParams->get(iN, 3u), cfgBetheBlochParams->get(iN, 4u))};
      double expSigma{expBethe * cfgBetheBlochParams->get(iN, 5u)};
      return static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
    };

    for (const auto& track : tracks) { // start loop over tracks
      if (track.itsNCls() >= cfgCutNclusExtremeIonisationITS) {
        double avgClsSize{0.};
        double cosL{std::sqrt(1. / (1. + track.tgl() * track.tgl()))};
        for (int iC{0}; iC < 7; ++iC) {
          avgClsSize += track.itsClsSizeInLayer(iC);
        }
        avgClsSize = avgClsSize * cosL / track.itsNCls();
        qaHists.fill(HIST("fExtremeIonisationITS"), track.itsNCls(), avgClsSize, track.p());
        keepEvent[kITSmildIonisation] = track.p() > cfgMomentumCutExtremeIonisation && avgClsSize > cfgCutClsSizeMildIonisation;
        keepEvent[kITSextremeIonisation] = track.p() > cfgMomentumCutExtremeIonisation && avgClsSize > cfgCutClsSizeExtremeIonisation;
      }
      if (track.itsNCls() < cfgCutNclusITS ||
          track.tpcNClsFound() < cfgCutNclusTPC) {
        continue;
      }

      if (std::abs(track.tpcNSigmaDe()) < 5) {
        qaHists.fill(HIST("fDeuTOFNsigma"), track.p() * track.sign(), track.tofNSigmaDe());
      }

      bool passesDCAselection{(track.sign() < 0 || (std::abs(track.dcaXY()) < cfgCutDCAxy &&
                                                    std::abs(track.dcaZ()) < cfgCutDCAz))};

      float nSigmaTPC[nNuclei]{
        track.tpcNSigmaDe(), track.tpcNSigmaTr(), track.tpcNSigmaHe()};
      const float nSigmaTOF[nNuclei]{
        track.tofNSigmaDe(), track.tofNSigmaTr(), track.tofNSigmaHe()};
      const int iC{track.sign() < 0};

      float fixTPCrigidity{(cfgFixTPCinnerParam && (track.pidForTracking() == track::PID::Helium3 || track.pidForTracking() == track::PID::Alpha)) ? 0.5f : 1.f};

      // fill QA hist: dEdx for all charged tracks
      qaHists.fill(HIST("fTPCsignalAll"), track.sign() * track.tpcInnerParam() * fixTPCrigidity, track.tpcSignal());

      for (int iN{0}; iN < nNuclei; ++iN) {
        /// Cheap checks first
        if (track.tpcInnerParam() * fixTPCrigidity < cfgMinTPCmom->get(iN, iC)) {
          continue;
        }

        if (cfgBetheBlochParams->get(iN, 5u) > 0.f) {
          nSigmaTPC[iN] = getNsigma(track, iN, iC);
        }
        h2TPCnSigma[iN]->Fill(track.sign() * track.tpcInnerParam() * fixTPCrigidity, nSigmaTPC[iN]);
        if (nSigmaTPC[iN] < cfgCutsPID->get(iN, 0u) || nSigmaTPC[iN] > cfgCutsPID->get(iN, 1u)) {
          continue;
        }
        if (track.p() > cfgCutsPID->get(iN, 4u) && (nSigmaTOF[iN] < cfgCutsPID->get(iN, 2u) || nSigmaTOF[iN] > cfgCutsPID->get(iN, 3u))) {
          continue;
        }
        if (iN == 1 && passesDCAselection) {
          h3indices.push_back(track.globalIndex());
          h3vectors.emplace_back(track.pt(), track.eta(), track.phi(), masses[iN]);
        }
        if (nucleusIndex[iN] < 0) {
          continue;
        }
        keepEvent[nucleusIndex[iN]] = passesDCAselection;
        if (keepEvent[nucleusIndex[iN]]) {
          h2TPCsignal[iN]->Fill(track.sign() * track.tpcInnerParam() * fixTPCrigidity, track.tpcSignal());
        }
      }
      //
      // fill QA histograms
      //
      qaHists.fill(HIST("fTPCsignal"), track.sign() * track.tpcInnerParam() * fixTPCrigidity, track.tpcSignal());

    } // end loop over tracks

    for (const auto& track : tracks) {
      if (track.itsNCls() < cfgCutNclusITS ||
          track.tpcNClsFound() < cfgCutNclusTPC ||
          std::abs(track.dcaXY()) > cfgCutDCAxy ||
          std::abs(track.dcaZ()) > cfgCutDCAz ||
          std::abs(track.eta()) > cfgCutEta) {
        continue;
      }
      const ROOT::Math::PtEtaPhiMVector trackVector(track.pt(), track.eta(), track.phi(), constants::physics::MassPiMinus);
      for (size_t iH3{0}; iH3 < h3vectors.size(); ++iH3) {
        if (h3indices[iH3] == track.globalIndex()) {
          continue;
        }
        const auto& h3vector = h3vectors[iH3];
        auto pivector = trackVector;
        auto cm = h3vector + trackVector;
        const ROOT::Math::Boost boost(cm.BoostToCM());
        boost(pivector);
        if (pivector.P() < cfgCutKstar) {
          keepEvent[kTritonFemto] = true;
          break;
        }
      }
    }

    for (const auto& v0 : v0s) {
      const auto& posTrack = v0.posTrack_as<TrackCandidates>();
      const auto& negTrack = v0.negTrack_as<TrackCandidates>();
      if ((posTrack.itsNCls() < cfgCutNclusITS || posTrack.tpcNClsFound() < cfgCutNclusTPC) &&
          (negTrack.itsNCls() < cfgCutNclusITS || negTrack.tpcNClsFound() < cfgCutNclusTPC)) {
        continue;
      }
      float nSigmas[2]{
        cfgBetheBlochParams->get(2, 5u) > 0.f ? getNsigma(posTrack, 2, 0) : posTrack.tpcNSigmaHe(),
        cfgBetheBlochParams->get(2, 5u) > 0.f ? getNsigma(negTrack, 2, 1) : negTrack.tpcNSigmaHe()};

      bool isHe3 = nSigmas[0] > cfgCutsPID->get(2, 0u) && nSigmas[0] < cfgCutsPID->get(2, 1u);
      bool isAntiHe3 = nSigmas[1] > cfgCutsPID->get(2, 0u) && nSigmas[1] < cfgCutsPID->get(2, 1u);
      if (!isHe3 && !isAntiHe3) {
        continue;
      }
      auto& he3Track = isHe3 ? posTrack : negTrack;
      auto& piTrack = isHe3 ? negTrack : posTrack;

      int n2bodyVtx = fitter2body.process(getTrackParCov(he3Track), getTrackParCov(piTrack));
      if (n2bodyVtx == 0) {
        continue;
      }
      auto vtxXYZ = fitter2body.getPCACandidate();
      vtxXYZ[0] -= collision.posX();
      vtxXYZ[1] -= collision.posY();
      vtxXYZ[2] -= collision.posZ();

      std::array<float, 3> momHe3 = {0.};
      std::array<float, 3> momPi = {0.};
      std::array<float, 3> momTot = {0.};
      auto& hePropTrack = fitter2body.getTrack(0);
      auto& piPropTrack = fitter2body.getTrack(1);
      hePropTrack.getPxPyPzGlo(momHe3);
      piPropTrack.getPxPyPzGlo(momPi);
      for (int i = 0; i < 3; ++i) {
        momHe3[i] *= 2;
        momTot[i] = momHe3[i] + momPi[i];
      }
      double cosPA = (vtxXYZ[0] * momTot[0] + vtxXYZ[1] * momTot[1] + vtxXYZ[2] * momTot[2]) /
                     std::sqrt((vtxXYZ[0] * vtxXYZ[0] + vtxXYZ[1] * vtxXYZ[1] + vtxXYZ[2] * vtxXYZ[2]) *
                               (momTot[0] * momTot[0] + momTot[1] * momTot[1] + momTot[2] * momTot[2]));
      if (cosPA < cfgCutCosPAheV0) {
        continue;
      }
      keepEvent[kHeV0] = true;
      break;
    }

    // fH3L3Body trigger
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    for (const auto& decay3body : decay3bodys) {
      auto track0 = decay3body.track0_as<TrackCandidates>();
      auto track1 = decay3body.track1_as<TrackCandidates>();
      auto track2 = decay3body.track2_as<TrackCandidates>();

      // track selection
      // keep like-sign triplets, do not check sign of deuteron, use same cut for p and pi
      if (track0.tpcNClsFound() < trgH3L3Body.mintpcNClspion || track1.tpcNClsFound() < trgH3L3Body.mintpcNClspion || track2.tpcNClsFound() < trgH3L3Body.mintpcNClsdeuteron) {
        continue;
      }

      if (std::abs(track0.eta()) > trgH3L3Body.minDaughtersEta || std::abs(track1.eta()) > trgH3L3Body.minDaughtersEta || std::abs(track2.eta()) > trgH3L3Body.minDaughtersEta) {
        continue;
      }

      bool isProton = false, isPion = false, isAntiProton = false, isAntiPion = false;
      if (std::abs(track0.tpcNSigmaPr()) < std::abs(track0.tpcNSigmaPi())) {
        if (track0.p() >= trgH3L3Body.minProtonPt && track0.p() <= trgH3L3Body.maxProtonPt) {
          if (track0.tpcNClsFound() >= trgH3L3Body.mintpcNClsproton) {
            isProton = true;
          }
        }
      }
      if (std::abs(track0.tpcNSigmaPi()) < std::abs(track0.tpcNSigmaPr())) {
        if (track0.p() >= trgH3L3Body.minPionPt && track0.p() <= trgH3L3Body.maxPionPt) {
          isPion = true;
        }
      }
      if (std::abs(track1.tpcNSigmaPr()) < std::abs(track1.tpcNSigmaPi())) {
        if (track1.p() >= trgH3L3Body.minProtonPt && track1.p() <= trgH3L3Body.maxProtonPt) {
          if (track1.tpcNClsFound() >= trgH3L3Body.mintpcNClsproton) {
            isAntiProton = true;
          }
        }
      }
      if (std::abs(track1.tpcNSigmaPi()) < std::abs(track1.tpcNSigmaPr())) {
        if (track1.p() >= trgH3L3Body.minPionPt && track1.p() <= trgH3L3Body.maxPionPt) {
          isAntiPion = true;
        }
      }

      if (!(isProton && isAntiPion) && !(isAntiProton && isPion)) {
        continue;
      }

      if (std::abs(track2.tpcNSigmaDe()) > trgH3L3Body.tpcPIDNSigmaCut || track2.p() < trgH3L3Body.minDeuteronPt || track2.p() > trgH3L3Body.maxDeuteronPt) {
        continue;
      }

      float tofNSigmaDeuteron = -999;
      if (track2.has_collision() && track2.hasTOF()) {
        auto originalcol = track2.collision_as<ColWithEvTime>();
        tofNSigmaDeuteron = bachelorTOFPID.GetTOFNSigma(mRespParamsV3, track2, originalcol, collision);
      }
      if (track2.p() > trgH3L3Body.minDeuteronPUseTOF && (tofNSigmaDeuteron < trgH3L3Body.tofPIDNSigmaMin || tofNSigmaDeuteron > trgH3L3Body.tofPIDNSigmaMax)) {
        continue;
      }

      // reconstruct the secondary vertex
      auto Track0 = getTrackParCov(track0);
      auto Track1 = getTrackParCov(track1);
      auto Track2 = getTrackParCov(track2);
      int n3bodyVtx = fitter3body.process(Track0, Track1, Track2);
      if (n3bodyVtx == 0) { // discard this pair
        continue;
      }

      std::array<float, 3> pos = {0.};
      const auto& vtxXYZ = fitter3body.getPCACandidate();
      for (int i = 0; i < 3; i++) {
        pos[i] = vtxXYZ[i];
      }

      // Calculate DCA with respect to the collision associated to the SV, not individual tracks
      std::array<float, 2> dcaInfo;

      auto track0Par = getTrackPar(track0);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track0Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
      auto track0dcaXY = dcaInfo[0];
      auto track0dca = std::sqrt(track0dcaXY * track0dcaXY + dcaInfo[1] * dcaInfo[1]);

      auto track1Par = getTrackPar(track1);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track1Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
      auto track1dcaXY = dcaInfo[0];
      auto track1dca = std::sqrt(track1dcaXY * track1dcaXY + dcaInfo[1] * dcaInfo[1]);

      std::array<float, 3> p0 = {0.}, p1 = {0.}, p2{0.};
      const auto& propagatedTrack0 = fitter3body.getTrack(0);
      const auto& propagatedTrack1 = fitter3body.getTrack(1);
      const auto& propagatedTrack2 = fitter3body.getTrack(2);
      propagatedTrack0.getPxPyPzGlo(p0);
      propagatedTrack1.getPxPyPzGlo(p1);
      propagatedTrack2.getPxPyPzGlo(p2);
      std::array<float, 3> p3BXYZ = {p0[0] + p1[0] + p2[0], p0[1] + p1[1] + p2[1], p0[2] + p1[2] + p2[2]};
      float sqpt3B = p3BXYZ[0] * p3BXYZ[0] + p3BXYZ[1] * p3BXYZ[1];
      float sqp3B = sqpt3B + p3BXYZ[2] * p3BXYZ[2];
      float pt3B = std::sqrt(sqpt3B), p3B = std::sqrt(sqp3B);

      if (p3B < trgH3L3Body.minP3Body) {
        continue;
      }

      float dcaDaughters = std::sqrt(fitter3body.getChi2AtPCACandidate());
      if (dcaDaughters > trgH3L3Body.dcavtxdau) {
        continue;
      }

      float vtxCosPA = RecoDecay::cpa(std::array{collision.posX(), collision.posY(), collision.posZ()}, std::array{pos[0], pos[1], pos[2]}, std::array{p3BXYZ[0], p3BXYZ[1], p3BXYZ[2]});
      if (vtxCosPA < trgH3L3Body.minCosPA3body) {
        continue;
      }
      float ct = std::sqrt(std::pow(pos[0] - collision.posX(), 2) + std::pow(pos[1] - collision.posY(), 2) + std::pow(pos[2] - collision.posZ(), 2)) / (p3B + 1E-10) * constants::physics::MassHyperTriton;
      if (ct > trgH3L3Body.lifetimecut) {
        continue;
      }

      float invmassH3L = RecoDecay::m(std::array{std::array{p0[0], p0[1], p0[2]}, std::array{p1[0], p1[1], p1[2]}, std::array{p2[0], p2[1], p2[2]}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
      float invmassAntiH3L = RecoDecay::m(std::array{std::array{p0[0], p0[1], p0[2]}, std::array{p1[0], p1[1], p1[2]}, std::array{p2[0], p2[1], p2[2]}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron});

      if (invmassH3L >= trgH3L3Body.h3LMassLowerlimit && invmassH3L <= trgH3L3Body.h3LMassUpperlimit) {
        // Hypertriton hypothesis
        if (isProton && isAntiPion && std::abs(track1dca) >= trgH3L3Body.dcapiontopv) {
          qaHists.fill(HIST("fH3LMassVsPt"), pt3B, invmassH3L);
          qaHists.fill(HIST("fBachDeuTOFNsigma"), track2.p() * track2.sign(), tofNSigmaDeuteron);
          qaHists.fill(HIST("fH3LDcaVsPt"), pt3B, dcaDaughters);
          qaHists.fill(HIST("fH3LCosPAVsPt"), pt3B, vtxCosPA);
          keepEvent[kH3L3Body] = true;
        }
      }
      if (invmassAntiH3L >= trgH3L3Body.h3LMassLowerlimit && invmassAntiH3L <= trgH3L3Body.h3LMassUpperlimit) {
        // Anti-Hypertriton hypothesis
        if (isAntiProton && isPion && std::abs(track0dca) >= trgH3L3Body.dcapiontopv) {
          qaHists.fill(HIST("fH3LMassVsPt"), pt3B, invmassAntiH3L);
          qaHists.fill(HIST("fBachDeuTOFNsigma"), track2.p() * track2.sign(), tofNSigmaDeuteron);
          qaHists.fill(HIST("fH3LDcaVsPt"), pt3B, dcaDaughters);
          qaHists.fill(HIST("fH3LCosPAVsPt"), pt3B, vtxCosPA);
          keepEvent[kH3L3Body] = true;
        }
      }
    }

    keepEvent[kTracked3Body] = tracked3Bodys.size() > 0;

    for (int iDecision{0}; iDecision < kNtriggers; ++iDecision) {
      if (keepEvent[iDecision]) {
        hProcessedEvents->Fill(iDecision + 1);
      }
    }

    tags(keepEvent[kH2], keepEvent[kHe], keepEvent[kHeV0], keepEvent[kTritonFemto], keepEvent[kH3L3Body], keepEvent[kTracked3Body], keepEvent[kITSmildIonisation], keepEvent[kITSextremeIonisation]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  metadataInfo.initMetadata(cfg);
  return WorkflowSpec{
    adaptAnalysisTask<nucleiFilter>(cfg)};
}
