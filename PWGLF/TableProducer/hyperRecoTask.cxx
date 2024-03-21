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
// ==============================================================================

#include <array>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DCAFitter/DCAFitterN.h"

#include "PWGLF/DataModel/LFHypernucleiTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;

using CollisionsFullWithFlow = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::EPCalibrationTables>;

namespace
{
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNames{"He3"};
std::shared_ptr<TH1> hEvents;
std::shared_ptr<TH1> hZvtx;
std::shared_ptr<TH1> hCentFT0A;
std::shared_ptr<TH1> hCentFT0C;
std::shared_ptr<TH1> hCentFT0M;
std::shared_ptr<TH1> hCentFV0A;
std::shared_ptr<TH2> hNsigma3HeSel;
std::shared_ptr<TH2> hDeDx3HeSel;
std::shared_ptr<TH2> hDeDxTot;
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

  int posTrackID;
  int negTrackID;
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
  uint8_t nTPCClustersHe3 = 0u;
  uint8_t nTPCClustersPi = 0u;
  uint32_t clusterSizeITSHe3 = 0u;
  uint32_t clusterSizeITSPi = 0u;
  bool isMatter = false;
  bool isSignal = false;          // true MC signal
  bool isReco = false;            // true if the candidate is actually reconstructed
  bool isRecoMCCollision = false; // true if the corresponding MC collision has been reconstructed
  bool isSurvEvSelection = false; // true if the corresponding event passed the event selection
  int pdgCode = 0;                // PDG code of the hypernucleus
  uint8_t flags = 0u;             // flags for dughter particles
};

struct hyperRecoTask {

  Produces<aod::DataHypCands> outputDataTable;
  Produces<aod::DataHypCandsFlow> outputDataTableWithFlow;
  Produces<aod::MCHypCands> outputMCTable;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Selection criteria
  Configurable<double> v0cospa{"hypcospa", 0.95, "V0 CosPA"};
  Configurable<float> masswidth{"hypmasswidth", 0.06, "Mass width (GeV/c^2)"};
  Configurable<float> dcav0dau{"hypdcaDau", 1.0, "DCA V0 Daughters"};
  Configurable<float> ptMin{"ptMin", 0.5, "Minimum pT of the hypercandidate"};
  Configurable<float> TPCRigidityMinHe{"TPCRigidityMinHe", 0.2, "Minimum rigidity of the helium candidate"};
  Configurable<float> etaMax{"eta", 1., "eta daughter"};
  Configurable<float> nSigmaMaxHe{"nSigmaMaxHe", 5, "helium dEdx cut (n sigma)"};
  Configurable<float> nTPCClusMinHe{"nTPCClusMinHe", 80, "helium NTPC clusters cut"};
  Configurable<bool> mcSignalOnly{"mcSignalOnly", true, "If true, save only signal in MC"};

  // Define o2 fitter, 2-prong, active memory (no need to redefine per event)
  o2::vertexing::DCAFitterN<2> fitter;

  // daughter masses
  float he3Mass = o2::constants::physics::MassHelium3;
  float he4Mass = o2::constants::physics::MassAlpha;
  float piMass = o2::constants::physics::MassPionCharged;

  // bethe bloch parameters
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};

  // CCDB options
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> pidPath{"pidPath", "", "Path to the PID response object"};

  // PDG codes
  Configurable<int> hyperPdg{"hyperPDG", 1010010030, "PDG code of the hyper-mother (could be 3LamH or 4LamH)"};
  Configurable<int> heDauPdg{"heDauPDG", 1000020030, "PDG code of the helium (could be 3He or 4He)"};

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
  std::vector<bool> isRecoCollision;
  std::vector<bool> isSurvEvSelCollision;

  Preslice<aod::V0s> perCollision = o2::aod::v0::collisionId;

  HistogramRegistry qaRegistry{"QA", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  float d_bz;
  std::array<float, 6> mBBparamsHe;

  void init(InitContext const&)
  {

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

    const AxisSpec rigidityAxis{rigidityBins, "#it{p}^{TPC}/#it{z}"};
    const AxisSpec dedxAxis{dedxBins, "d#it{E}/d#it{x}"};
    const AxisSpec nSigma3HeAxis{nSigmaBins, "n_{#sigma}({}^{3}He)"};
    const AxisSpec zVtxAxis{zVtxBins, "z_{vtx} (cm)"};
    const AxisSpec centAxis{centBins, "Centrality"};

    hNsigma3HeSel = qaRegistry.add<TH2>("hNsigma3HeSel", "; p_{TPC}/z (GeV/#it{c}); n_{#sigma} ({}^{3}He)", HistType::kTH2F, {rigidityAxis, nSigma3HeAxis});
    hDeDx3HeSel = qaRegistry.add<TH2>("hDeDx3HeSel", ";p_{TPC}/z (GeV/#it{c}); dE/dx", HistType::kTH2F, {rigidityAxis, dedxAxis});
    hDeDxTot = qaRegistry.add<TH2>("hDeDxTot", ";p_{TPC}/z (GeV/#it{c}); dE/dx", HistType::kTH2F, {rigidityAxis, dedxAxis});
    hEvents = qaRegistry.add<TH1>("hEvents", ";Events; ", HistType::kTH1D, {{2, -0.5, 1.5}});
    hEvents->GetXaxis()->SetBinLabel(1, "All");
    hEvents->GetXaxis()->SetBinLabel(2, "Selected");
    if (doprocessMC) {
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
    hCentFV0A = qaRegistry.add<TH1>("hCentFV0A", ";Centrality; ", HistType::kTH1D, {{100, 0, 100}});
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
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

  template <class Tcoll>
  void fillCandidateData(Tcoll const& collision, aod::V0s const& V0s)
  {
    if (mBBparamsHe[5] < 0) {
      LOG(fatal) << "Bethe-Bloch parameters for He3 not set, please check your CCDB and configuration";
    }
    for (auto& v0 : V0s) {

      auto posTrack = v0.posTrack_as<TracksFull>();
      auto negTrack = v0.negTrack_as<TracksFull>();

      if (std::abs(posTrack.eta()) > etaMax || std::abs(negTrack.eta()) > etaMax)
        continue;

      // temporary fix: tpcInnerParam() returns the momentum in all the software tags before: https://github.com/AliceO2Group/AliceO2/pull/12521
      bool posHeliumPID = posTrack.pidForTracking() == o2::track::PID::Helium3 || posTrack.pidForTracking() == o2::track::PID::Alpha;
      bool negHeliumPID = negTrack.pidForTracking() == o2::track::PID::Helium3 || negTrack.pidForTracking() == o2::track::PID::Alpha;
      float posRigidity = posHeliumPID ? posTrack.tpcInnerParam() / 2 : posTrack.tpcInnerParam();
      float negRigidity = negHeliumPID ? negTrack.tpcInnerParam() / 2 : negTrack.tpcInnerParam();

      hDeDxTot->Fill(posRigidity, posTrack.tpcSignal());
      hDeDxTot->Fill(-negRigidity, negTrack.tpcSignal());

      double expBethePos{tpc::BetheBlochAleph(static_cast<float>(posRigidity * 2 / constants::physics::MassHelium3), mBBparamsHe[0], mBBparamsHe[1], mBBparamsHe[2], mBBparamsHe[3], mBBparamsHe[4])};
      double expBetheNeg{tpc::BetheBlochAleph(static_cast<float>(negRigidity * 2 / constants::physics::MassHelium3), mBBparamsHe[0], mBBparamsHe[1], mBBparamsHe[2], mBBparamsHe[3], mBBparamsHe[4])};
      double expSigmaPos{expBethePos * mBBparamsHe[5]};
      double expSigmaNeg{expBetheNeg * mBBparamsHe[5]};
      auto nSigmaTPCpos = static_cast<float>((posTrack.tpcSignal() - expBethePos) / expSigmaPos);
      auto nSigmaTPCneg = static_cast<float>((negTrack.tpcSignal() - expBetheNeg) / expSigmaNeg);

      // ITS only tracks do not have TPC information. TPCnSigma: only lower cut to allow for both hypertriton and hyperhydrogen4 reconstruction
      bool isHe = posTrack.hasTPC() && nSigmaTPCpos > -1 * nSigmaMaxHe;
      bool isAntiHe = negTrack.hasTPC() && nSigmaTPCneg > -1 * nSigmaMaxHe;

      if (!isHe && !isAntiHe)
        continue;

      hyperCandidate hypCand;
      hypCand.isMatter = isHe && isAntiHe ? std::abs(nSigmaTPCpos) < std::abs(nSigmaTPCneg) : isHe;
      auto& he3track = hypCand.isMatter ? posTrack : negTrack;
      auto& he3Rigidity = hypCand.isMatter ? posRigidity : negRigidity;
      if (he3track.tpcNClsFound() < nTPCClusMinHe || he3Rigidity < TPCRigidityMinHe)
        continue;

      hypCand.nSigmaHe3 = hypCand.isMatter ? nSigmaTPCpos : nSigmaTPCneg;
      hypCand.nTPCClustersHe3 = hypCand.isMatter ? posTrack.tpcNClsFound() : negTrack.tpcNClsFound();
      hypCand.tpcSignalHe3 = hypCand.isMatter ? posTrack.tpcSignal() : negTrack.tpcSignal();
      hypCand.clusterSizeITSHe3 = hypCand.isMatter ? posTrack.itsClusterSizes() : negTrack.itsClusterSizes();
      hypCand.nTPCClustersPi = !hypCand.isMatter ? posTrack.tpcNClsFound() : negTrack.tpcNClsFound();
      hypCand.tpcSignalPi = !hypCand.isMatter ? posTrack.tpcSignal() : negTrack.tpcSignal();
      hypCand.clusterSizeITSPi = !hypCand.isMatter ? posTrack.itsClusterSizes() : negTrack.itsClusterSizes();
      hypCand.momHe3TPC = hypCand.isMatter ? posRigidity : negRigidity;
      hypCand.momPiTPC = !hypCand.isMatter ? posRigidity : negRigidity;

      hypCand.flags |= hypCand.isMatter ? static_cast<uint8_t>((posTrack.pidForTracking() & 0xF) << 4) : static_cast<uint8_t>((negTrack.pidForTracking() & 0xF) << 4);
      hypCand.flags |= hypCand.isMatter ? static_cast<uint8_t>(negTrack.pidForTracking() & 0xF) : static_cast<uint8_t>(posTrack.pidForTracking() & 0xF);

      auto posTrackCov = getTrackParCov(posTrack);
      auto negTrackCov = getTrackParCov(negTrack);

      int nCand = 0;
      try {
        nCand = fitter.process(posTrackCov, negTrackCov);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        continue;
      }
      if (nCand == 0) {
        continue;
      }

      auto& hePropTrack = hypCand.isMatter ? fitter.getTrack(0) : fitter.getTrack(1);
      auto& piPropTrack = hypCand.isMatter ? fitter.getTrack(1) : fitter.getTrack(0);
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
        continue;

      float massH3L = std::sqrt(h3lE * h3lE - hypMom[0] * hypMom[0] - hypMom[1] * hypMom[1] - hypMom[2] * hypMom[2]);
      float massH4L = std::sqrt(h4lE * h4lE - hypMom[0] * hypMom[0] - hypMom[1] * hypMom[1] - hypMom[2] * hypMom[2]);
      bool isHypMass = false;
      if (massH3L > o2::constants::physics::MassHyperTriton - masswidth && massH3L < o2::constants::physics::MassHyperTriton + masswidth)
        isHypMass = true;
      if (massH4L > o2::constants::physics::MassHyperhydrog4 - masswidth && massH4L < o2::constants::physics::MassHyperhydrog4 + masswidth)
        isHypMass = true;
      if (!isHypMass)
        continue;

      hypCand.dcaV0dau = std::sqrt(fitter.getChi2AtPCACandidate());
      if (hypCand.dcaV0dau > dcav0dau) {
        continue;
      }

      std::array<float, 3> primVtx = {collision.posX(), collision.posY(), collision.posZ()};

      double cosPA = RecoDecay::cpa(primVtx, hypCand.decVtx, hypMom);
      if (cosPA < v0cospa) {
        continue;
      }

      for (int i = 0; i < 3; i++) {
        hypCand.decVtx[i] = hypCand.decVtx[i] - primVtx[i];
      }

      // if survived all selections, propagate decay daughters to PV
      gpu::gpustd::array<float, 2> dcaInfo;

      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, posTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
      hypCand.isMatter ? hypCand.he3DCAXY = dcaInfo[0] : hypCand.piDCAXY = dcaInfo[0];

      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, negTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
      hypCand.isMatter ? hypCand.piDCAXY = dcaInfo[0] : hypCand.he3DCAXY = dcaInfo[0];

      // finally, push back the candidate
      hypCand.isReco = true;
      hypCand.posTrackID = posTrack.globalIndex();
      hypCand.negTrackID = negTrack.globalIndex();

      int chargeFactor = -1 + 2 * hypCand.isMatter;
      hDeDx3HeSel->Fill(chargeFactor * hypCand.momHe3TPC, he3track.tpcSignal());
      hNsigma3HeSel->Fill(chargeFactor * hypCand.momHe3TPC, hypCand.nSigmaHe3);

      hyperCandidates.push_back(hypCand);
    }
  }

  void fillMCinfo(aod::McTrackLabels const& trackLabels, aod::McParticles const& particlesMC)
  {
    for (auto& hypCand : hyperCandidates) {
      auto mcLabPos = trackLabels.rawIteratorAt(hypCand.posTrackID);
      auto mcLabNeg = trackLabels.rawIteratorAt(hypCand.negTrackID);

      if (mcLabPos.has_mcParticle() && mcLabNeg.has_mcParticle()) {
        auto mcTrackPos = mcLabPos.mcParticle_as<aod::McParticles>();
        auto mcTrackNeg = mcLabNeg.mcParticle_as<aod::McParticles>();
        if (mcTrackPos.has_mothers() && mcTrackNeg.has_mothers()) {
          for (auto& negMother : mcTrackNeg.mothers_as<aod::McParticles>()) {
            for (auto& posMother : mcTrackPos.mothers_as<aod::McParticles>()) {
              if (posMother.globalIndex() != negMother.globalIndex())
                continue;
              if (!((mcTrackPos.pdgCode() == heDauPdg && mcTrackNeg.pdgCode() == -211) || (mcTrackPos.pdgCode() == 211 && mcTrackNeg.pdgCode() == -1 * heDauPdg))) // TODO: check warning for constant comparison
                continue;
              if (std::abs(posMother.pdgCode()) != hyperPdg)
                continue;

              auto posPrimVtx = array{posMother.vx(), posMother.vy(), posMother.vz()};
              auto secVtx = array{mcTrackPos.vx(), mcTrackPos.vy(), mcTrackPos.vz()};
              hypCand.gMom = array{posMother.px(), posMother.py(), posMother.pz()};
              hypCand.gMomHe3 = mcTrackPos.pdgCode() == heDauPdg ? array{mcTrackPos.px(), mcTrackPos.py(), mcTrackPos.pz()} : array{mcTrackNeg.px(), mcTrackNeg.py(), mcTrackNeg.pz()};
              for (int i = 0; i < 3; i++) {
                hypCand.gDecVtx[i] = secVtx[i] - posPrimVtx[i];
              }
              hypCand.isSignal = true;
              hypCand.pdgCode = posMother.pdgCode();
              hypCand.isRecoMCCollision = isRecoCollision[posMother.mcCollisionId()];
              hypCand.isSurvEvSelection = isSurvEvSelCollision[posMother.mcCollisionId()];
              filledMothers.push_back(posMother.globalIndex());
            }
          }
        }
      }
    }
  }

  void processData(CollisionsFull const& collisions, aod::V0s const& V0s, TracksFull const& tracks, aod::BCsWithTimestamps const&)
  {

    for (const auto& collision : collisions) {
      hyperCandidates.clear();

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      hEvents->Fill(0.);
      if (!collision.sel8() || std::abs(collision.posZ()) > 10 || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        continue;
      hEvents->Fill(1.);
      hZvtx->Fill(collision.posZ());
      hCentFT0A->Fill(collision.centFT0A());
      hCentFT0C->Fill(collision.centFT0C());
      hCentFT0M->Fill(collision.centFT0M());
      hCentFV0A->Fill(collision.centFV0A());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollision, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      fillCandidateData(collision, V0Table_thisCollision);

      for (auto& hypCand : hyperCandidates) {
        outputDataTable(collision.centFT0A(), collision.centFT0C(), collision.centFT0M(),
                        collision.posX(), collision.posY(), collision.posZ(),
                        hypCand.isMatter,
                        hypCand.recoPtHe3(), hypCand.recoPhiHe3(), hypCand.recoEtaHe3(),
                        hypCand.recoPtPi(), hypCand.recoPhiPi(), hypCand.recoEtaPi(),
                        hypCand.decVtx[0], hypCand.decVtx[1], hypCand.decVtx[2],
                        hypCand.dcaV0dau, hypCand.he3DCAXY, hypCand.piDCAXY,
                        hypCand.nSigmaHe3, hypCand.nTPCClustersHe3, hypCand.nTPCClustersPi,
                        hypCand.momHe3TPC, hypCand.momPiTPC, hypCand.tpcSignalHe3, hypCand.tpcSignalPi,
                        hypCand.clusterSizeITSHe3, hypCand.clusterSizeITSPi, hypCand.flags);
      }
    }
  }
  PROCESS_SWITCH(hyperRecoTask, processData, "Data analysis", true);

  void processDataWithFlow(CollisionsFullWithFlow const& collisions, aod::V0s const& V0s, TracksFull const& tracks, aod::BCsWithTimestamps const&)
  {

    for (const auto& collision : collisions) {
      hyperCandidates.clear();

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      hEvents->Fill(0.);
      if (!collision.sel8() || std::abs(collision.posZ()) > 10 || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        continue;
      hEvents->Fill(1.);
      hZvtx->Fill(collision.posZ());
      hCentFT0A->Fill(collision.centFT0A());
      hCentFT0C->Fill(collision.centFT0C());
      hCentFT0M->Fill(collision.centFT0M());
      hCentFV0A->Fill(collision.centFV0A());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollision, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      fillCandidateData(collision, V0Table_thisCollision);

      for (auto& hypCand : hyperCandidates) {
        outputDataTableWithFlow(collision.centFT0A(), collision.centFT0C(), collision.centFT0M(),
                                collision.psiFT0A(), collision.multFT0A(),
                                collision.psiFT0C(), collision.multFT0C(),
                                collision.psiTPC(), collision.multTPC(),
                                collision.posX(), collision.posY(), collision.posZ(),
                                hypCand.isMatter,
                                hypCand.recoPtHe3(), hypCand.recoPhiHe3(), hypCand.recoEtaHe3(),
                                hypCand.recoPtPi(), hypCand.recoPhiPi(), hypCand.recoEtaPi(),
                                hypCand.decVtx[0], hypCand.decVtx[1], hypCand.decVtx[2],
                                hypCand.dcaV0dau, hypCand.he3DCAXY, hypCand.piDCAXY,
                                hypCand.nSigmaHe3, hypCand.nTPCClustersHe3, hypCand.nTPCClustersPi,
                                hypCand.momHe3TPC, hypCand.momPiTPC, hypCand.tpcSignalHe3, hypCand.tpcSignalPi,
                                hypCand.clusterSizeITSHe3, hypCand.clusterSizeITSPi, hypCand.flags);
      }
    }
  }
  PROCESS_SWITCH(hyperRecoTask, processDataWithFlow, "Data analysis with flow", false);

  void processMC(CollisionsFullMC const& collisions, aod::McCollisions const& mcCollisions, aod::V0s const& V0s, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC)
  {
    filledMothers.clear();
    isRecoCollision.clear();
    isSurvEvSelCollision.clear();
    isRecoCollision.resize(mcCollisions.size(), false);
    isSurvEvSelCollision.resize(mcCollisions.size(), false);

    for (const auto& collision : collisions) {
      hyperCandidates.clear();
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      hEvents->Fill(0.);

      if (collision.has_mcCollision()) {
        isRecoCollision[collision.mcCollisionId()] = true;
      }

      if (!collision.sel8() || std::abs(collision.posZ()) > 10)
        continue;
      hEvents->Fill(1.);
      hZvtx->Fill(collision.posZ());
      hCentFT0A->Fill(collision.centFT0A());
      hCentFT0C->Fill(collision.centFT0C());
      hCentFT0M->Fill(collision.centFT0M());
      hCentFV0A->Fill(collision.centFV0A());

      if (collision.has_mcCollision()) {
        isSurvEvSelCollision[collision.mcCollisionId()] = true;
      }

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollision, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      fillCandidateData(collision, V0Table_thisCollision);
      fillMCinfo(trackLabelsMC, particlesMC);

      for (auto& hypCand : hyperCandidates) {
        if (!hypCand.isSignal && mcSignalOnly)
          continue;
        int chargeFactor = -1 + 2 * (hypCand.pdgCode > 0);
        outputMCTable(collision.centFT0A(), collision.centFT0C(), collision.centFT0M(),
                      collision.posX(), collision.posY(), collision.posZ(),
                      hypCand.isMatter,
                      hypCand.recoPtHe3(), hypCand.recoPhiHe3(), hypCand.recoEtaHe3(),
                      hypCand.recoPtPi(), hypCand.recoPhiPi(), hypCand.recoEtaPi(),
                      hypCand.decVtx[0], hypCand.decVtx[1], hypCand.decVtx[2],
                      hypCand.dcaV0dau, hypCand.he3DCAXY, hypCand.piDCAXY,
                      hypCand.nSigmaHe3, hypCand.nTPCClustersHe3, hypCand.nTPCClustersPi,
                      hypCand.momHe3TPC, hypCand.momPiTPC, hypCand.tpcSignalHe3, hypCand.tpcSignalPi,
                      hypCand.clusterSizeITSHe3, hypCand.clusterSizeITSPi, hypCand.flags,
                      chargeFactor * hypCand.genPt(), hypCand.genPhi(), hypCand.genEta(), hypCand.genPtHe3(),
                      hypCand.gDecVtx[0], hypCand.gDecVtx[1], hypCand.gDecVtx[2],
                      hypCand.isReco, hypCand.isSignal, hypCand.isRecoMCCollision, hypCand.isSurvEvSelection);
      }
    }

    // now we fill only the signal candidates that were not reconstructed
    for (auto& mcPart : particlesMC) {

      if (std::abs(mcPart.pdgCode()) != hyperPdg)
        continue;
      std::array<float, 3> secVtx;
      std::array<float, 3> primVtx = {mcPart.vx(), mcPart.vy(), mcPart.vz()};
      std::array<float, 3> momMother = {mcPart.px(), mcPart.py(), mcPart.pz()};
      std::array<float, 3> momHe3;
      bool isHeFound = false;
      for (auto& mcDaught : mcPart.daughters_as<aod::McParticles>()) {
        if (std::abs(mcDaught.pdgCode()) == heDauPdg) {
          secVtx = {mcDaught.vx(), mcDaught.vy(), mcDaught.vz()};
          momHe3 = {mcDaught.px(), mcDaught.py(), mcDaught.pz()};
          isHeFound = true;
          break;
        }
      }
      if (mcPart.pdgCode() > 0) {
        hIsMatterGen->Fill(0.);
      } else {
        hIsMatterGen->Fill(1.);
      }
      if (!isHeFound) {
        hDecayChannel->Fill(1.);
        continue;
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
      hypCand.isRecoMCCollision = isRecoCollision[mcPart.mcCollisionId()];
      hypCand.isSurvEvSelection = isSurvEvSelCollision[mcPart.mcCollisionId()];
      int chargeFactor = -1 + 2 * (hypCand.pdgCode > 0);
      for (int i = 0; i < 3; i++) {
        hypCand.gDecVtx[i] = secVtx[i] - primVtx[i];
        hypCand.gMom[i] = momMother[i];
        hypCand.gMomHe3[i] = momHe3[i];
      }
      hypCand.posTrackID = -1;
      hypCand.negTrackID = -1;
      hypCand.isSignal = true;
      outputMCTable(-1, -1, -1,
                    -1, -1, -1,
                    0,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1, -1,
                    -1, -1, -1,
                    chargeFactor * hypCand.genPt(), hypCand.genPhi(), hypCand.genEta(), hypCand.genPtHe3(),
                    hypCand.gDecVtx[0], hypCand.gDecVtx[1], hypCand.gDecVtx[2],
                    hypCand.isReco, hypCand.isSignal, hypCand.isRecoMCCollision, hypCand.isSurvEvSelection);
    }
  }
  PROCESS_SWITCH(hyperRecoTask, processMC, "MC analysis", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hyperRecoTask>(cfgc)};
}
