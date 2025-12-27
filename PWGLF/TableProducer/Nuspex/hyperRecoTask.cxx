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
// Build hypertriton candidates from V0s and tracks

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFHypernucleiTables.h"
#include "PWGLF/Utils/svPoolCreator.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "CCDB/BasicCCDBManager.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "MathUtils/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using CollBracket = o2::math_utils::Bracket<int>;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TOFSignal, aod::TOFEvTime>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms>;

using CollisionsFullWithFlow = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::EPCalibrationTables>;

namespace
{
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleName{"He3"};
std::shared_ptr<TH1> hEvents;
std::shared_ptr<TH1> hEventsZorro;
std::shared_ptr<TH1> hZvtx;
std::shared_ptr<TH1> hCentFT0A;
std::shared_ptr<TH1> hCentFT0C;
std::shared_ptr<TH1> hCentFT0M;
std::shared_ptr<TH2> hNsigma3HeSel;
std::shared_ptr<TH2> hDeDx3HeSel;
std::shared_ptr<TH2> hDeDxTot;
std::shared_ptr<TH1> hH3LMassBefSel;
std::shared_ptr<TH1> hH3LMassTracked;
std::shared_ptr<TH1> hH4LMassBefSel;
std::shared_ptr<TH1> hH4LMassTracked;
std::shared_ptr<TH1> hDecayChannel;
std::shared_ptr<TH1> hIsMatterGen;
std::shared_ptr<TH1> hIsMatterGenTwoBody;
} // namespace

struct hyperCandidate {
  float recoPtHe3() const { return std::hypot(momHe3[0], momHe3[1]); }
  float recoPhiHe3() const { return std::atan2(momHe3[1], momHe3[0]); }
  float recoEtaHe3() const { return std::asinh(momHe3[2] / recoPtHe3()); }
  float recoPtPi() const { return std::hypot(momPi[0], momPi[1]); }
  float recoPhiPi() const { return std::atan2(momPi[1], momPi[0]); }
  float recoEtaPi() const { return std::asinh(momPi[2] / recoPtPi()); }
  float genPt() const { return std::hypot(gMom[0], gMom[1]); }
  float genPtHe3() const { return std::hypot(gMomHe3[0], gMomHe3[1]); }
  float genPhi() const { return std::atan2(gMom[1], gMom[0]); }
  float genEta() const { return std::asinh(gMom[2] / genPt()); }

  int v0ID = -1;
  int heTrackID;
  int piTrackID;
  float dcaV0dau = -10;
  float cosPA = -10;
  float nSigmaHe3 = -10;
  float he3DCAXY = -10;
  float piDCAXY = -10;
  float momHe3TPC = -10.f;
  float momPiTPC = -10.f;
  std::array<float, 3> momHe3;
  std::array<float, 3> momPi;
  std::array<float, 3> decVtx;
  std::array<float, 3> gMom;
  std::array<float, 3> gMomHe3;
  std::array<float, 3> gDecVtx;
  uint16_t tpcSignalHe3 = 0u;
  uint16_t tpcSignalPi = 0u;
  float tpcChi2He3 = 0.f;
  float itsChi2He3 = 0.f;
  float itsChi2Pi = 0.f;
  float massTOFHe3 = 0.f;
  uint8_t nTPCClustersHe3 = 0u;
  uint8_t nTPCClustersPi = 0u;
  uint8_t nTPCpidClusHe3 = 0u;
  uint8_t nTPCpidClusPi = 0u;
  uint32_t clusterSizeITSHe3 = 0u;
  uint32_t clusterSizeITSPi = 0u;

  // collision information
  int64_t collisionID = 0;

  bool isMatter = false;
  bool isSignal = false;           // true MC signal
  bool isReco = false;             // true if the candidate is actually reconstructed
  uint8_t isFakeHeOnITSLayer = 0u; // bit map for fake He on ITS layers

  bool isRecoMCCollision = false; // true if the corresponding MC collision has been reconstructed
  bool isSurvEvSelection = false; // true if the corresponding event passed the event selection
  int pdgCode = 0;                // PDG code of the hypernucleus
  uint8_t flags = 0u;             // flags for dughter particles
};

struct hyperRecoTask {

  Produces<aod::DataHypCands> outputDataTable;
  Produces<aod::DataHypCandsFlow> outputDataTableWithFlow;
  Produces<aod::MCHypCands> outputMCTable;
  Produces<aod::DataHypCandsWColl> outputDataTableWithCollID;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // PDG codes
  Configurable<int> hyperPdg{"hyperPDG", 1010010030, "PDG code of the hyper-mother (could be 3LamH or 4LamH)"};
  Configurable<int> heDauPdg{"heDauPDG", 1000020030, "PDG code of the helium (could be 3He or 4He)"};

  // Selection criteria
  Configurable<double> v0cospacut{"hypcospa", 0.95, "V0 CosPA"};
  Configurable<float> masswidth{"hypmasswidth", 0.06, "Mass width (GeV/c^2)"};
  Configurable<float> dcaToPvPion{"dcapvPi", 0., "DCA to PV pion"};
  Configurable<float> dcaToPvHe{"dcapvHe", 0., "DCA to PV helium"};
  Configurable<float> dcav0dau{"hypdcaDau", 1.0, "DCA V0 Daughters"};
  Configurable<float> ptMin{"ptMin", 0.5, "Minimum pT of the hypercandidate"};
  Configurable<float> TPCRigidityMinHe{"TPCRigidityMinHe", 0.2, "Minimum rigidity of the helium candidate"};
  Configurable<float> etaMax{"eta", 1., "eta daughter"};
  Configurable<float> nSigmaMaxHe{"nSigmaMaxHe", 5, "helium dEdx cut (n sigma)"};
  Configurable<float> nTPCClusMinHe{"nTPCClusMinHe", 70, "helium NTPC clusters cut"};
  Configurable<float> nTPCClusMinPi{"nTPCClusMinPi", -1., "pion NTPC clusters cut"};
  Configurable<bool> mcSignalOnly{"mcSignalOnly", true, "If true, save only signal in MC"};
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};
  Configurable<bool> isEventUsedForEPCalibration{"isEventUsedForEPCalibration", 1, "Event is used for EP calibration"};

  // Define o2 fitter, 2-prong, active memory (no need to redefine per event)
  o2::vertexing::DCAFitterN<2> fitter;
  svPoolCreator svCreator{heDauPdg, PDG_t::kPiPlus};

  // daughter masses
  float he3Mass = o2::constants::physics::MassHelium3;
  float he4Mass = o2::constants::physics::MassAlpha;
  float piMass = o2::constants::physics::MassPionCharged;

  Configurable<bool> useCustomVertexer{"useCustomVertexer", false, "Use custom vertexer"};
  Configurable<bool> skipAmbiTracks{"skipAmbiTracks", false, "Skip ambiguous tracks"};
  Configurable<bool> disableITSROFCut{"disableITSROFCut", false, "Disable ITS ROC cut for event selection"};
  Configurable<float> customVertexerTimeMargin{"customVertexerTimeMargin", 800, "Time margin for custom vertexer (ns)"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleName, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};
  Configurable<bool> cfgCompensatePIDinTracking{"cfgCompensatePIDinTracking", true, "If true, divide tpcInnerParam by the electric charge"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};

  // CCDB options
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> pidPath{"pidPath", "", "Path to the PID response object"};

  // histogram axes
  ConfigurableAxis rigidityBins{"rigidityBins", {200, -10.f, 10.f}, "Binning for rigidity #it{p}^{TPC}/#it{z}"};
  ConfigurableAxis dedxBins{"dedxBins", {1000, 0.f, 1000.f}, "Binning for dE/dx"};
  ConfigurableAxis nSigmaBins{"nSigmaBins", {200, -5.f, 5.f}, "Binning for n sigma"};
  ConfigurableAxis zVtxBins{"zVtxBins", {100, -20.f, 20.f}, "Binning for n sigma"};
  ConfigurableAxis centBins{"centBins", {100, 0.f, 100.f}, "Binning for centrality"};

  // std vector of candidates
  std::vector<hyperCandidate> hyperCandidates;
  // vector to keep track of MC mothers already filled
  std::vector<unsigned int> filledMothers;
  // vector to keep track of the collisions passing the event selection in the MC
  std::vector<int> recoCollisionIds;
  std::vector<bool> isSurvEvSelCollision;
  std::vector<bool> goodCollision;
  std::vector<int> trackedClSize;

  Preslice<aod::V0s> perCollision = o2::aod::v0::collisionId;

  HistogramRegistry qaRegistry{"QA", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  float d_bz;
  std::array<float, 6> mBBparamsHe;

  void init(InitContext const&)
  {

    zorroSummary.setObject(zorro.getZorroSummary());

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    int mat{static_cast<int>(cfgMaterialCorrection)};
    fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));

    svCreator.setTimeMargin(customVertexerTimeMargin);
    if (skipAmbiTracks) {
      svCreator.setSkipAmbiTracks();
    }

    const AxisSpec rigidityAxis{rigidityBins, "#it{p}^{TPC}/#it{z}"};
    const AxisSpec dedxAxis{dedxBins, "d#it{E}/d#it{x}"};
    const AxisSpec nSigma3HeAxis{nSigmaBins, "n_{#sigma}({}^{3}He)"};
    const AxisSpec zVtxAxis{zVtxBins, "z_{vtx} (cm)"};
    const AxisSpec centAxis{centBins, "Centrality"};

    hNsigma3HeSel = qaRegistry.add<TH2>("hNsigma3HeSel", "; p_{TPC}/z (GeV/#it{c}); n_{#sigma} ({}^{3}He)", HistType::kTH2F, {rigidityAxis, nSigma3HeAxis});
    hDeDx3HeSel = qaRegistry.add<TH2>("hDeDx3HeSel", ";p_{TPC}/z (GeV/#it{c}); dE/dx", HistType::kTH2F, {rigidityAxis, dedxAxis});
    hDeDxTot = qaRegistry.add<TH2>("hDeDxTot", ";p_{TPC}/z (GeV/#it{c}); dE/dx", HistType::kTH2F, {rigidityAxis, dedxAxis});
    hH3LMassBefSel = qaRegistry.add<TH1>("hH3LMassBefSel", ";M (GeV/#it{c}^{2}); ", HistType::kTH1D, {{60, 2.96, 3.04}});
    hH3LMassTracked = qaRegistry.add<TH1>("hH3LMassTracked", ";M (GeV/#it{c}^{2}); ", HistType::kTH1D, {{60, 2.96, 3.04}});
    hH4LMassBefSel = qaRegistry.add<TH1>("hH4LMassBefSel", ";M (GeV/#it{c}^{2}); ", HistType::kTH1D, {{60, 3.76, 3.84}});
    hH4LMassTracked = qaRegistry.add<TH1>("hH4LMassTracked", ";M (GeV/#it{c}^{2}); ", HistType::kTH1D, {{60, 3.76, 3.84}});

    hEvents = qaRegistry.add<TH1>("hEvents", ";Events; ", HistType::kTH1D, {{2, -0.5, 1.5}});
    hEvents->GetXaxis()->SetBinLabel(1, "All");
    hEvents->GetXaxis()->SetBinLabel(2, "Selected");

    hEventsZorro = qaRegistry.add<TH1>("hEventsZorro", ";Events; ", HistType::kTH1D, {{2, -0.5, 1.5}});
    hEventsZorro->GetXaxis()->SetBinLabel(1, "Zorro before evsel");
    hEventsZorro->GetXaxis()->SetBinLabel(2, "Zorro after evsel");

    if (doprocessMC || doprocessMCTracked) {
      hDecayChannel = qaRegistry.add<TH1>("hDecayChannel", ";Decay channel; ", HistType::kTH1D, {{2, -0.5, 1.5}});
      hDecayChannel->GetXaxis()->SetBinLabel(1, "2-body");
      hDecayChannel->GetXaxis()->SetBinLabel(2, "3-body");
      hIsMatterGen = qaRegistry.add<TH1>("hIsMatterGen", ";; ", HistType::kTH1D, {{2, -0.5, 1.5}});
      hIsMatterGen->GetXaxis()->SetBinLabel(1, "Matter");
      hIsMatterGen->GetXaxis()->SetBinLabel(2, "Antimatter");
      hIsMatterGenTwoBody = qaRegistry.add<TH1>("hIsMatterGenTwoBody", ";; ", HistType::kTH1D, {{2, -0.5, 1.5}});
      hIsMatterGenTwoBody->GetXaxis()->SetBinLabel(1, "Matter");
      hIsMatterGenTwoBody->GetXaxis()->SetBinLabel(2, "Antimatter");
    }
    hZvtx = qaRegistry.add<TH1>("hZvtx", ";z_{vtx} (cm); ", HistType::kTH1D, {{100, -20, 20}});
    hCentFT0A = qaRegistry.add<TH1>("hCentFT0A", ";Centrality; ", HistType::kTH1D, {{100, 0, 100}});
    hCentFT0C = qaRegistry.add<TH1>("hCentFT0C", ";Centrality; ", HistType::kTH1D, {{100, 0, 100}});
    hCentFT0M = qaRegistry.add<TH1>("hCentFT0M", ";Centrality; ", HistType::kTH1D, {{100, 0, 100}});
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fHe");
      zorro.populateHistRegistry(qaRegistry, bc.runNumber());
    }
    auto run3grp_timestamp = bc.timestamp();

    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
      } else {
        d_bz = d_bz_input;
      }
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
      } else {
        d_bz = d_bz_input;
      }
    }
    if (!pidPath.value.empty()) {
      auto he3pid = ccdb->getForTimeStamp<std::array<float, 6>>(pidPath.value + "_He3", run3grp_timestamp);
      std::copy(he3pid->begin(), he3pid->end(), mBBparamsHe.begin());
    } else {
      for (int i = 0; i < 5; i++) {
        mBBparamsHe[i] = cfgBetheBlochParams->get("He3", Form("p%i", i));
      }
      mBBparamsHe[5] = cfgBetheBlochParams->get("He3", "resolution");
    }
    fitter.setBz(d_bz);
    mRunNumber = bc.runNumber();
  }

  template <typename T>
  float computeNSigmaHe3(const T& candidate)
  {
    bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParam = (heliumPID && cfgCompensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
    float expTPCSignal = o2::common::BetheBlochAleph(static_cast<float>(correctedTPCinnerParam * 2.f / constants::physics::MassHelium3), mBBparamsHe[0], mBBparamsHe[1], mBBparamsHe[2], mBBparamsHe[3], mBBparamsHe[4]);
    double resoTPC{expTPCSignal * mBBparamsHe[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <class Tcoll>
  void selectGoodCollisions(const Tcoll& collisions)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      hEvents->Fill(0.);

      if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder) && !disableITSROFCut) {
        continue;
      }

      bool zorroSelected = false;
      if (cfgSkimmedProcessing) {
        // accounting done after ITS border cut, to properly correct with the MC
        zorroSelected = zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC());
        if (zorroSelected) {
          hEventsZorro->Fill(0.);
        }
      }

      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || std::abs(collision.posZ()) > 10) {
        continue;
      }

      if (zorroSelected) {
        hEventsZorro->Fill(1.);
      }

      goodCollision[collision.globalIndex()] = true;
      hEvents->Fill(1.);
      hZvtx->Fill(collision.posZ());
      hCentFT0A->Fill(collision.centFT0A());
      hCentFT0C->Fill(collision.centFT0C());
      hCentFT0M->Fill(collision.centFT0M());
    }
  }

  template <class Tcoll>
  void selectGoodCollisionsMC(const Tcoll& collisions)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      hEvents->Fill(0.);
      if (collision.has_mcCollision()) {
        recoCollisionIds[collision.mcCollisionId()] = collision.globalIndex();
      }
      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || std::abs(collision.posZ()) > 10)
        continue;

      if (collision.has_mcCollision()) {
        isSurvEvSelCollision[collision.mcCollisionId()] = true;
      }
      goodCollision[collision.globalIndex()] = true;
      hEvents->Fill(1.);
      hZvtx->Fill(collision.posZ());
      hCentFT0A->Fill(collision.centFT0A());
      hCentFT0C->Fill(collision.centFT0C());
      hCentFT0M->Fill(collision.centFT0M());
    }
  }

  template <class Ttrack, class Tcolls>
  void fillHyperCand(Ttrack& heTrack, Ttrack& piTrack, CollBracket collBracket, const Tcolls& collisions, hyperCandidate& hypCand)
  {

    hypCand.isMatter = heTrack.sign() > 0;
    hypCand.nSigmaHe3 = computeNSigmaHe3(heTrack);
    hypCand.nTPCClustersHe3 = heTrack.tpcNClsFound();
    hypCand.tpcSignalHe3 = heTrack.tpcSignal();
    hypCand.nTPCpidClusHe3 = static_cast<int16_t>(heTrack.tpcNClsFindable()) - heTrack.tpcNClsFindableMinusPID();
    hypCand.clusterSizeITSHe3 = heTrack.itsClusterSizes();
    hypCand.nTPCClustersPi = piTrack.tpcNClsFound();
    hypCand.nTPCpidClusPi = static_cast<int16_t>(piTrack.tpcNClsFindable()) - piTrack.tpcNClsFindableMinusPID();
    hypCand.tpcSignalPi = piTrack.tpcSignal();
    hypCand.tpcChi2He3 = heTrack.tpcChi2NCl();
    hypCand.itsChi2He3 = heTrack.itsChi2NCl();
    hypCand.itsChi2Pi = piTrack.itsChi2NCl();
    hypCand.clusterSizeITSPi = piTrack.itsClusterSizes();
    bool heliumPID = heTrack.pidForTracking() == o2::track::PID::Helium3 || heTrack.pidForTracking() == o2::track::PID::Alpha;
    hypCand.momHe3TPC = (heliumPID && cfgCompensatePIDinTracking) ? heTrack.tpcInnerParam() / 2 : heTrack.tpcInnerParam();
    if (hypCand.momHe3TPC < TPCRigidityMinHe)
      return;
    hypCand.momPiTPC = piTrack.tpcInnerParam();
    hDeDxTot->Fill(hypCand.momHe3TPC * heTrack.sign(), heTrack.tpcSignal());
    hDeDxTot->Fill(hypCand.momPiTPC * piTrack.sign(), piTrack.tpcSignal());
    hypCand.flags |= static_cast<uint8_t>((heTrack.pidForTracking() & 0xF) << 4);
    hypCand.flags |= static_cast<uint8_t>(piTrack.pidForTracking() & 0xF);
    auto heTrackCov = getTrackParCov(heTrack);
    auto piTrackCov = getTrackParCov(piTrack);

    int nCand = 0;
    try {
      nCand = fitter.process(heTrackCov, piTrackCov);
    } catch (...) {
      LOG(error) << "Exception caught in DCA fitter process call!";
      return;
    }
    if (nCand == 0) {
      return;
    }

    auto& hePropTrack = fitter.getTrack(0);
    auto& piPropTrack = fitter.getTrack(1);
    hePropTrack.getPxPyPzGlo(hypCand.momHe3);
    piPropTrack.getPxPyPzGlo(hypCand.momPi);
    // the momentum has to be multiplied by 2 (charge)
    for (int i = 0; i < 3; i++) {
      hypCand.momHe3[i] *= 2;
    }
    float heP2 = hypCand.momHe3[0] * hypCand.momHe3[0] + hypCand.momHe3[1] * hypCand.momHe3[1] + hypCand.momHe3[2] * hypCand.momHe3[2];
    float piP2 = hypCand.momPi[0] * hypCand.momPi[0] + hypCand.momPi[1] * hypCand.momPi[1] + hypCand.momPi[2] * hypCand.momPi[2];
    float he3E = std::sqrt(heP2 + he3Mass * he3Mass);
    float he4E = std::sqrt(heP2 + he4Mass * he4Mass);
    float piE = std::sqrt(piP2 + piMass * piMass);
    float h3lE = he3E + piE;
    float h4lE = he4E + piE;
    std::array<float, 3> hypMom;
    const auto& vtx = fitter.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      hypCand.decVtx[i] = vtx[i];
      hypMom[i] = hypCand.momHe3[i] + hypCand.momPi[i];
    }
    float hypPt = std::hypot(hypMom[0], hypMom[1]);
    if (hypPt < ptMin)
      return;
    float massH3L = std::sqrt(h3lE * h3lE - hypMom[0] * hypMom[0] - hypMom[1] * hypMom[1] - hypMom[2] * hypMom[2]);
    float massH4L = std::sqrt(h4lE * h4lE - hypMom[0] * hypMom[0] - hypMom[1] * hypMom[1] - hypMom[2] * hypMom[2]);
    bool isHypMass = false;
    if (massH3L > o2::constants::physics::MassHyperTriton - masswidth && massH3L < o2::constants::physics::MassHyperTriton + masswidth)
      isHypMass = true;
    if (massH4L > o2::constants::physics::MassHyperhydrog4 - masswidth && massH4L < o2::constants::physics::MassHyperhydrog4 + masswidth)
      isHypMass = true;
    if (!isHypMass)
      return;

    hH3LMassBefSel->Fill(massH3L);
    hH4LMassBefSel->Fill(massH4L);
    if (!trackedClSize.empty() && trackedClSize[hypCand.v0ID] > 0) {
      hH3LMassTracked->Fill(massH3L);
      hH4LMassTracked->Fill(massH4L);
    }

    hypCand.dcaV0dau = std::sqrt(fitter.getChi2AtPCACandidate());
    if (hypCand.dcaV0dau > dcav0dau) {
      return;
    }

    // select the collision that maximizes the cosine of the pointing angle
    double cosPAmax = -2;
    unsigned int collIDmax = 0;

    for (int collID = collBracket.getMin(); collID <= collBracket.getMax(); collID++) {
      auto collision = collisions.rawIteratorAt(collID);
      std::array<float, 3> collVtx = {collision.posX(), collision.posY(), collision.posZ()};
      double cosPA = RecoDecay::cpa(collVtx, hypCand.decVtx, hypMom);
      if (cosPA > cosPAmax) {
        cosPAmax = cosPA;
        collIDmax = collID;
      }
    }

    if (!goodCollision[collIDmax]) {
      return;
    }

    if (cosPAmax < v0cospacut) {
      return;
    }

    auto collision = collisions.rawIteratorAt(collIDmax);
    std::array<float, 3> primVtx = {collision.posX(), collision.posY(), collision.posZ()};
    for (int i = 0; i < 3; i++) {
      hypCand.decVtx[i] = hypCand.decVtx[i] - primVtx[i];
    }

    // if survived all selections, propagate decay daughters to PV
    std::array<float, 2> dcaInfo;
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, heTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
    hypCand.he3DCAXY = dcaInfo[0];
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, piTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
    hypCand.piDCAXY = dcaInfo[0];

    if (std::abs(hypCand.piDCAXY) < dcaToPvPion || std::abs(hypCand.he3DCAXY) < dcaToPvHe) {
      return;
    }

    // finally, fill collision info and push back the candidate
    hypCand.isReco = true;
    hypCand.heTrackID = heTrack.globalIndex();
    hypCand.piTrackID = piTrack.globalIndex();
    hypCand.collisionID = collision.globalIndex();

    if (heTrack.hasTOF()) {
      float beta = o2::pid::tof::Beta::GetBeta(heTrack);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      hypCand.massTOFHe3 = hypCand.momHe3TPC * 2.f * std::sqrt(1.f / (beta * beta) - 1.f);
    }

    hDeDx3HeSel->Fill(heTrack.sign() * hypCand.momHe3TPC, heTrack.tpcSignal());
    hNsigma3HeSel->Fill(heTrack.sign() * hypCand.momHe3TPC, hypCand.nSigmaHe3);
    hyperCandidates.push_back(hypCand);
  }

  template <class Tcolls, class Ttracks>
  void fillV0s(const Tcolls& collisions, const Ttracks& tracks, aod::V0s const& V0s)
  {
    if (mBBparamsHe[5] < 0) {
      LOG(fatal) << "Bethe-Bloch parameters for He3 not set, please check your CCDB and configuration";
    }
    for (const auto& v0 : V0s) {
      // if(v0.isStandardV0())
      //   continue;
      auto posTrack = tracks.rawIteratorAt(v0.posTrackId());
      auto negTrack = tracks.rawIteratorAt(v0.negTrackId());

      if (std::abs(posTrack.eta()) > etaMax || std::abs(negTrack.eta()) > etaMax)
        continue;

      // temporary fix: tpcInnerParam() returns the momentum in all the software tags before: https://github.com/AliceO2Group/AliceO2/pull/12521
      auto nSigmaTPCpos = computeNSigmaHe3(posTrack);
      auto nSigmaTPCneg = computeNSigmaHe3(negTrack);
      // ITS only tracks do not have TPC information. TPCnSigma: only lower cut to allow for both hypertriton and hyperhydrogen4 reconstruction
      bool isHe = posTrack.hasTPC() && nSigmaTPCpos > -1 * nSigmaMaxHe;
      bool isAntiHe = negTrack.hasTPC() && nSigmaTPCneg > -1 * nSigmaMaxHe;
      if (!isHe && !isAntiHe) {
        continue;
      }
      auto& heTrack = isHe ? posTrack : negTrack;
      auto& piTrack = isHe ? negTrack : posTrack;
      if (heTrack.tpcNClsFound() < nTPCClusMinHe || piTrack.tpcNClsFound() < nTPCClusMinPi) {
        continue;
      }

      hyperCandidate hypCand;
      hypCand.v0ID = v0.globalIndex();

      int collIdx = v0.collision_as<Tcolls>().globalIndex();
      CollBracket collBracket{collIdx, collIdx};

      fillHyperCand(heTrack, piTrack, collBracket, collisions, hypCand);
    }
  }

  template <class Tcolls, class Ttracks>
  void fillCustomV0s(const Tcolls& collisions, const Ttracks& tracks, aod::AmbiguousTracks const& ambiguousTracks, aod::BCsWithTimestamps const& bcs)
  {

    svCreator.clearPools();
    svCreator.fillBC2Coll(collisions, bcs);

    for (const auto& track : tracks) {

      if (std::abs(track.eta()) > etaMax)
        continue;

      if (!track.hasITS())
        continue;

      auto nSigmaHe = computeNSigmaHe3(track);
      bool isHe = nSigmaHe > -1 * nSigmaMaxHe;
      int pdgHypo = isHe ? heDauPdg : PDG_t::kPiPlus;
      // LOG(info) << "ncls found: " << track.tpcNClsFound();
      if (isHe && track.tpcNClsFound() < nTPCClusMinHe)
        continue;
      if (!isHe && track.tpcNClsFound() < nTPCClusMinPi)
        continue;

      svCreator.appendTrackCand(track, collisions, pdgHypo, ambiguousTracks, bcs);
    }
    auto& svPool = svCreator.getSVCandPool(collisions);
    LOG(debug) << "SV pool size: " << svPool.size();

    for (const auto& svCand : svPool) {
      auto heTrack = tracks.rawIteratorAt(svCand.tr0Idx);
      auto piTrack = tracks.rawIteratorAt(svCand.tr1Idx);
      auto collIdxs = svCand.collBracket;
      hyperCandidate hypCand;
      hypCand.v0ID = -1;
      fillHyperCand(heTrack, piTrack, collIdxs, collisions, hypCand);
    }
  }
  void fillMCinfo(aod::McTrackLabels const& trackLabels, aod::McParticles const&)
  {
    for (auto& hypCand : hyperCandidates) {
      auto mcLabHe = trackLabels.rawIteratorAt(hypCand.heTrackID);
      auto mcLabPi = trackLabels.rawIteratorAt(hypCand.piTrackID);

      if (mcLabHe.has_mcParticle() && mcLabPi.has_mcParticle()) {
        auto mcTrackHe = mcLabHe.mcParticle_as<aod::McParticles>();
        auto mcTrackPi = mcLabPi.mcParticle_as<aod::McParticles>();
        if (mcTrackHe.has_mothers() && mcTrackPi.has_mothers()) {
          for (const auto& heMother : mcTrackHe.mothers_as<aod::McParticles>()) {
            for (const auto& piMother : mcTrackPi.mothers_as<aod::McParticles>()) {
              if (heMother.globalIndex() != piMother.globalIndex())
                continue;
              if (std::abs(mcTrackHe.pdgCode()) != heDauPdg || std::abs(mcTrackPi.pdgCode()) != PDG_t::kPiPlus)
                continue;
              if (std::abs(heMother.pdgCode()) != hyperPdg)
                continue;

              auto primVtx = std::array<float, 3>{heMother.vx(), heMother.vy(), heMother.vz()};
              auto secVtx = std::array<float, 3>{mcTrackHe.vx(), mcTrackHe.vy(), mcTrackHe.vz()};
              hypCand.gMom = std::array<float, 3>{heMother.px(), heMother.py(), heMother.pz()};
              hypCand.gMomHe3 = std::array<float, 3>{mcTrackHe.px(), mcTrackHe.py(), mcTrackHe.pz()};
              for (int i = 0; i < 3; i++) {
                hypCand.gDecVtx[i] = secVtx[i] - primVtx[i];
              }
              hypCand.isSignal = true;
              hypCand.isFakeHeOnITSLayer = mcLabHe.mcMask() & 0x7F; // check if any of the first 7 bits is set
              hypCand.pdgCode = heMother.pdgCode();
              hypCand.isRecoMCCollision = recoCollisionIds[heMother.mcCollisionId()] > 0;
              hypCand.isSurvEvSelection = isSurvEvSelCollision[heMother.mcCollisionId()];
              filledMothers.push_back(heMother.globalIndex());
            }
          }
        }
      }
    }
  }

  void processDataTracked(CollisionsFull const& collisions, aod::V0s const& V0s, aod::TrackedV0s const& tV0s, TracksFull const& tracks, aod::AmbiguousTracks const& ambiTracks, aod::BCsWithTimestamps const& bcs)
  {
    trackedClSize.clear();
    trackedClSize.resize(V0s.size(), 0);
    for (const auto& tV0 : tV0s) {
      trackedClSize[tV0.v0Id()] = tV0.itsClsSize();
    }
    processData(collisions, V0s, tracks, ambiTracks, bcs);
  }
  PROCESS_SWITCH(hyperRecoTask, processDataTracked, "Data analysis wit tracked V0s information", false);

  void processData(CollisionsFull const& collisions, aod::V0s const& V0s, TracksFull const& tracks, aod::AmbiguousTracks const& ambiTracks, aod::BCsWithTimestamps const& bcs)
  {
    goodCollision.clear();
    goodCollision.resize(collisions.size(), false);
    hyperCandidates.clear();

    selectGoodCollisions(collisions);
    useCustomVertexer ? fillCustomV0s(collisions, tracks, ambiTracks, bcs) : fillV0s(collisions, tracks, V0s);

    for (const auto& hypCand : hyperCandidates) {
      auto collision = collisions.rawIteratorAt(hypCand.collisionID);
      float trackedHypClSize = !trackedClSize.empty() ? trackedClSize[hypCand.v0ID] : 0;
      outputDataTable(collision.centFT0A(), collision.centFT0C(), collision.centFT0M(),
                      collision.posX(), collision.posY(), collision.posZ(),
                      mRunNumber, hypCand.isMatter,
                      hypCand.recoPtHe3(), hypCand.recoPhiHe3(), hypCand.recoEtaHe3(),
                      hypCand.recoPtPi(), hypCand.recoPhiPi(), hypCand.recoEtaPi(),
                      hypCand.decVtx[0], hypCand.decVtx[1], hypCand.decVtx[2],
                      hypCand.dcaV0dau, hypCand.he3DCAXY, hypCand.piDCAXY,
                      hypCand.nSigmaHe3, hypCand.nTPCClustersHe3, hypCand.nTPCClustersPi,
                      hypCand.nTPCpidClusHe3, hypCand.nTPCpidClusPi,
                      hypCand.momHe3TPC, hypCand.momPiTPC, hypCand.tpcSignalHe3, hypCand.tpcSignalPi, hypCand.tpcChi2He3, hypCand.itsChi2He3, hypCand.itsChi2Pi,
                      hypCand.massTOFHe3,
                      hypCand.clusterSizeITSHe3, hypCand.clusterSizeITSPi, hypCand.flags, trackedHypClSize);
    }
  }
  PROCESS_SWITCH(hyperRecoTask, processData, "Data analysis", true);

  void processDataWithFlow(CollisionsFullWithFlow const& collisions, aod::V0s const& V0s, TracksFull const& tracks, aod::AmbiguousTracks const& ambiTracks, aod::BCsWithTimestamps const& bcs)
  {

    goodCollision.clear();
    goodCollision.resize(collisions.size(), false);
    hyperCandidates.clear();

    selectGoodCollisions(collisions);
    useCustomVertexer ? fillCustomV0s(collisions, tracks, ambiTracks, bcs) : fillV0s(collisions, tracks, V0s);

    for (const auto& hypCand : hyperCandidates) {
      auto collision = collisions.rawIteratorAt(hypCand.collisionID);
      if (isEventUsedForEPCalibration && !collision.triggereventep()) {
        return;
      }
      float trackedHypClSize = !trackedClSize.empty() ? trackedClSize[hypCand.v0ID] : 0;
      outputDataTableWithFlow(collision.centFT0A(), collision.centFT0C(), collision.centFT0M(),
                              collision.psiFT0A(), collision.multFT0A(),
                              collision.psiFT0C(), collision.multFT0C(), collision.qFT0C(),
                              collision.psiTPC(), collision.multTPC(),
                              collision.posX(), collision.posY(), collision.posZ(),
                              mRunNumber, hypCand.isMatter,
                              hypCand.recoPtHe3(), hypCand.recoPhiHe3(), hypCand.recoEtaHe3(),
                              hypCand.recoPtPi(), hypCand.recoPhiPi(), hypCand.recoEtaPi(),
                              hypCand.decVtx[0], hypCand.decVtx[1], hypCand.decVtx[2],
                              hypCand.dcaV0dau, hypCand.he3DCAXY, hypCand.piDCAXY,
                              hypCand.nSigmaHe3, hypCand.nTPCClustersHe3, hypCand.nTPCClustersPi,
                              hypCand.nTPCpidClusHe3, hypCand.nTPCpidClusPi,
                              hypCand.momHe3TPC, hypCand.momPiTPC, hypCand.tpcSignalHe3, hypCand.tpcSignalPi, hypCand.tpcChi2He3, hypCand.itsChi2He3, hypCand.itsChi2Pi,
                              hypCand.massTOFHe3,
                              hypCand.clusterSizeITSHe3, hypCand.clusterSizeITSPi, hypCand.flags, trackedHypClSize);
    }
  }
  PROCESS_SWITCH(hyperRecoTask, processDataWithFlow, "Data analysis with flow", false);

  void processDataWithCollID(CollisionsFull const& collisions, aod::V0s const& V0s, TracksFull const& tracks, aod::AmbiguousTracks const& ambiTracks, aod::BCsWithTimestamps const& bcs)
  {
    goodCollision.clear();
    goodCollision.resize(collisions.size(), false);
    hyperCandidates.clear();

    selectGoodCollisions(collisions);
    useCustomVertexer ? fillCustomV0s(collisions, tracks, ambiTracks, bcs) : fillV0s(collisions, tracks, V0s);

    for (const auto& hypCand : hyperCandidates) {
      auto collision = collisions.rawIteratorAt(hypCand.collisionID);
      float trackedHypClSize = !trackedClSize.empty() ? trackedClSize[hypCand.v0ID] : 0;
      outputDataTableWithCollID(hypCand.collisionID, collision.centFT0A(), collision.centFT0C(), collision.centFT0M(),
                                collision.posX(), collision.posY(), collision.posZ(),
                                mRunNumber, hypCand.isMatter,
                                hypCand.recoPtHe3(), hypCand.recoPhiHe3(), hypCand.recoEtaHe3(),
                                hypCand.recoPtPi(), hypCand.recoPhiPi(), hypCand.recoEtaPi(),
                                hypCand.decVtx[0], hypCand.decVtx[1], hypCand.decVtx[2],
                                hypCand.dcaV0dau, hypCand.he3DCAXY, hypCand.piDCAXY,
                                hypCand.nSigmaHe3, hypCand.nTPCClustersHe3, hypCand.nTPCClustersPi,
                                hypCand.nTPCpidClusHe3, hypCand.nTPCpidClusPi,
                                hypCand.momHe3TPC, hypCand.momPiTPC, hypCand.tpcSignalHe3, hypCand.tpcSignalPi, hypCand.tpcChi2He3, hypCand.itsChi2He3, hypCand.itsChi2Pi,
                                hypCand.massTOFHe3,
                                hypCand.clusterSizeITSHe3, hypCand.clusterSizeITSPi, hypCand.flags, trackedHypClSize);
    }
  }
  PROCESS_SWITCH(hyperRecoTask, processDataWithCollID, "Data analysis with collision ID", false);

  void processMC(CollisionsFullMC const& collisions, aod::McCollisions const& mcCollisions, aod::V0s const& V0s, TracksFull const& tracks, aod::AmbiguousTracks const& ambiTracks, aod::BCsWithTimestamps const& bcs, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC)
  {
    filledMothers.clear();
    recoCollisionIds.clear();
    recoCollisionIds.resize(mcCollisions.size(), -1);
    isSurvEvSelCollision.clear();
    isSurvEvSelCollision.resize(mcCollisions.size(), false);
    goodCollision.clear();
    goodCollision.resize(collisions.size(), false);
    hyperCandidates.clear();

    selectGoodCollisionsMC(collisions);
    useCustomVertexer ? fillCustomV0s(collisions, tracks, ambiTracks, bcs) : fillV0s(collisions, tracks, V0s);
    fillMCinfo(trackLabelsMC, particlesMC);
    for (const auto& hypCand : hyperCandidates) {
      auto collision = collisions.rawIteratorAt(hypCand.collisionID);
      if (!hypCand.isSignal && mcSignalOnly)
        continue;
      int chargeFactor = -1 + 2 * (hypCand.pdgCode > 0);
      float trackedHypClSize = !trackedClSize.empty() ? trackedClSize[hypCand.v0ID] : 0;
      outputMCTable(collision.centFT0A(), collision.centFT0C(), collision.centFT0M(),
                    collision.posX(), collision.posY(), collision.posZ(),
                    mRunNumber, hypCand.isMatter,
                    hypCand.recoPtHe3(), hypCand.recoPhiHe3(), hypCand.recoEtaHe3(),
                    hypCand.recoPtPi(), hypCand.recoPhiPi(), hypCand.recoEtaPi(),
                    hypCand.decVtx[0], hypCand.decVtx[1], hypCand.decVtx[2],
                    hypCand.dcaV0dau, hypCand.he3DCAXY, hypCand.piDCAXY,
                    hypCand.nSigmaHe3, hypCand.nTPCClustersHe3, hypCand.nTPCClustersPi, hypCand.nTPCpidClusHe3, hypCand.nTPCpidClusPi,
                    hypCand.momHe3TPC, hypCand.momPiTPC, hypCand.tpcSignalHe3, hypCand.tpcSignalPi, hypCand.tpcChi2He3, hypCand.itsChi2He3, hypCand.itsChi2Pi,
                    hypCand.massTOFHe3,
                    hypCand.clusterSizeITSHe3, hypCand.clusterSizeITSPi, hypCand.flags, trackedHypClSize,
                    chargeFactor * hypCand.genPt(), hypCand.genPhi(), hypCand.genEta(), hypCand.genPtHe3(),
                    hypCand.gDecVtx[0], hypCand.gDecVtx[1], hypCand.gDecVtx[2],
                    hypCand.isReco, hypCand.isFakeHeOnITSLayer, hypCand.isSignal, hypCand.isRecoMCCollision, hypCand.isSurvEvSelection, 1, 0);
    }

    // now we fill only the signal candidates that were not reconstructed
    for (const auto& mcPart : particlesMC) {

      if (std::abs(mcPart.pdgCode()) != hyperPdg)
        continue;
      std::array<float, 3> secVtx;
      std::array<float, 3> primVtx = {mcPart.vx(), mcPart.vy(), mcPart.vz()};
      std::array<float, 3> momMother = {mcPart.px(), mcPart.py(), mcPart.pz()};
      std::array<float, 3> momHe3;
      bool isHeFound = false;
      int mcProcess = {0};
      for (const auto& mcDaught : mcPart.daughters_as<aod::McParticles>()) {
        if (std::abs(mcDaught.pdgCode()) == heDauPdg) {
          secVtx = {mcDaught.vx(), mcDaught.vy(), mcDaught.vz()};
          momHe3 = {mcDaught.px(), mcDaught.py(), mcDaught.pz()};
          isHeFound = true;
        }
        if (mcDaught.pdgCode() != PDG_t::kElectron) { // we do not care about delta electrons
          mcProcess = mcDaught.getProcess();
        }
      }
      if (mcPart.pdgCode() > 0) {
        hIsMatterGen->Fill(0.);
      } else {
        hIsMatterGen->Fill(1.);
      }
      if (!isHeFound) {
        hDecayChannel->Fill(1.);
      }
      hDecayChannel->Fill(0.);
      if (mcPart.pdgCode() > 0) {
        hIsMatterGenTwoBody->Fill(0.);
      } else {
        hIsMatterGenTwoBody->Fill(1.);
      }
      if (std::find(filledMothers.begin(), filledMothers.end(), mcPart.globalIndex()) != std::end(filledMothers)) {
        continue;
      }
      hyperCandidate hypCand;
      hypCand.pdgCode = mcPart.pdgCode();
      hypCand.isRecoMCCollision = recoCollisionIds[mcPart.mcCollisionId()] > 0;
      hypCand.isSurvEvSelection = isSurvEvSelCollision[mcPart.mcCollisionId()];
      int chargeFactor = -1 + 2 * (hypCand.pdgCode > 0);
      for (int i = 0; i < 3; i++) {
        hypCand.gDecVtx[i] = secVtx[i] - primVtx[i];
        hypCand.gMom[i] = momMother[i];
        hypCand.gMomHe3[i] = momHe3[i];
      }
      hypCand.heTrackID = -1;
      hypCand.piTrackID = -1;
      hypCand.isSignal = true;

      float centFT0A = -1, centFT0C = -1, centFT0M = -1;
      if (hypCand.isRecoMCCollision) {
        auto recoCollision = collisions.rawIteratorAt(recoCollisionIds[mcPart.mcCollisionId()]);
        centFT0A = recoCollision.centFT0A();
        centFT0C = recoCollision.centFT0C();
        centFT0M = recoCollision.centFT0M();
      }

      outputMCTable(centFT0A, centFT0C, centFT0M,
                    mRunNumber, -1, -1, -1,
                    0,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, 0, 0, 0, 0,
                    -1, -1, -1, false,
                    chargeFactor * hypCand.genPt(), hypCand.genPhi(), hypCand.genEta(), hypCand.genPtHe3(),
                    hypCand.gDecVtx[0], hypCand.gDecVtx[1], hypCand.gDecVtx[2],
                    hypCand.isReco, -1, hypCand.isSignal, hypCand.isRecoMCCollision, hypCand.isSurvEvSelection, isHeFound, mcProcess);
    }
  }
  PROCESS_SWITCH(hyperRecoTask, processMC, "MC analysis", false);

  void processMCTracked(CollisionsFullMC const& collisions, aod::McCollisions const& mcCollisions, aod::V0s const& V0s, aod::TrackedV0s const& tV0s, TracksFull const& tracks, aod::AmbiguousTracks const& ambiTracks, aod::BCsWithTimestamps const& bcs, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC)
  {
    trackedClSize.clear();
    trackedClSize.resize(V0s.size(), 0);
    for (const auto& tV0 : tV0s) {
      trackedClSize[tV0.v0Id()] = tV0.itsClsSize();
    }
    processMC(collisions, mcCollisions, V0s, tracks, ambiTracks, bcs, trackLabelsMC, particlesMC);
  }
  PROCESS_SWITCH(hyperRecoTask, processMCTracked, "MC analysis with tracked V0s", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hyperRecoTask>(cfgc)};
}
