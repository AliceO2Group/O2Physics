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
// Build \Lambda-n-n candidates from V0s and tracks
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

#include "PWGLF/DataModel/LFLnnTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;

namespace
{
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNames{"3H"};

std::shared_ptr<TH1> hEvents;
std::shared_ptr<TH1> hZvtx;
std::shared_ptr<TH1> hCentFT0A;
std::shared_ptr<TH1> hCentFT0C;
std::shared_ptr<TH1> hCentFT0M;
std::shared_ptr<TH1> hCentFV0A;
std::shared_ptr<TH2> hNsigma3HSel;
std::shared_ptr<TH2> hdEdx3HSel;
std::shared_ptr<TH2> hdEdx3HTPCMom;
std::shared_ptr<TH2> hdEdxTot;
std::shared_ptr<TH1> hDecayChannel;
std::shared_ptr<TH1> hIsMatterGen;
std::shared_ptr<TH1> hIsMatterGenTwoBody;

} // namespace

struct lnnCandidate {
  float recoPt3H() const { return std::hypot(mom3H[0], mom3H[1]); }
  float recoPhi3H() const { return std::atan2(mom3H[1], mom3H[0]); }
  float recoEta3H() const { return std::asinh(mom3H[2] / recoPt3H()); }
  float recoPtPi() const { return std::hypot(momPi[0], momPi[1]); }
  float recoPhiPi() const { return std::atan2(momPi[1], momPi[0]); }
  float recoEtaPi() const { return std::asinh(momPi[2] / recoPtPi()); }
  float genPt() const { return std::hypot(gMom[0], gMom[1]); }
  float genPt3H() const { return std::hypot(gMom3H[0], gMom3H[1]); }
  float genPhi() const { return std::atan2(gMom[1], gMom[0]); }
  float genEta() const { return std::asinh(gMom[2] / genPt()); }

  int posTrackID;
  int negTrackID;
  float dcaV0dau = -10;
  float cosPA = -10;
  float nSigma3H = -10;
  float h3DCAXY = -10;
  float piDCAXY = -10;
  float mom3HTPC = -10.f;
  float momPiTPC = -10.f;
  std::array<float, 3> mom3H;
  std::array<float, 3> momPi;
  std::array<float, 3> decVtx;
  std::array<float, 3> gMom;
  std::array<float, 3> gMom3H;
  std::array<float, 3> gDecVtx;
  uint16_t tpcSignal3H = 0u;
  uint16_t tpcSignalPi = 0u;
  uint8_t nTPCClusters3H = 0u;
  uint8_t nTPCClustersPi = 0u;
  uint32_t clusterSizeITS3H = 0u;
  uint32_t clusterSizeITSPi = 0u;
  bool isMatter = false;
  bool isSignal = false;        // true MC signal
  bool isReco = false;          // true if the candidate is actually reconstructed
  bool survEvSelection = false; // true if the corresponding event passed the event selection
  int pdgCode = 0;              // PDG code of the hypernucleus
  uint8_t flags = 0u;           // flags for dughter particles
};

struct lnnRecoTask {

  Produces<aod::DataLnnCands> outputDataTable;
  Produces<aod::MCLnnCands> outputMCTable;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Selection criteria
  Configurable<double> v0cospa{"lnncospa", 0.95, "V0 CosPA"};
  Configurable<float> masswidth{"lnnmasswidth", 0.006, "Mass width (GeV/c^2)"};
  Configurable<float> dcav0dau{"lnndcaDau", 1.0, "DCA V0 Daughters"};
  Configurable<float> ptMin{"ptMin", 0.5, "Minimum pT of the lnncandidate"};
  Configurable<float> TPCRigidityMin3H{"TPCRigidityMin3H", 1, "Minimum rigidity of the triton candidate"};
  Configurable<float> etaMax{"eta", 1., "eta daughter"};
  Configurable<float> nSigmaMax3H{"nSigmaMax3H", 5, "triton dEdx cut (n sigma)"};
  Configurable<float> nTPCClusMin3H{"nTPCClusMin3H", 80, "triton NTPC clusters cut"};
  Configurable<bool> mcSignalOnly{"mcSignalOnly", true, "If true, save only signal in MC"};

  // Define o2 fitter, 2-prong, active memory (no need to redefine per event)
  o2::vertexing::DCAFitterN<2> fitter;

  // daughter masses
  float h3Mass = o2::constants::physics::MassTriton;
  float piMass = o2::constants::physics::MassPionCharged;

  // bethe bloch parameters
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for 3H"};
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
  Configurable<int> h3DauPdg{"h3DauPdg", 1000010030, "PDG Triton"}; // PDG Triton
  Configurable<int> lnnPdg{"lnnPdg", 1010000030, "PDG Lnn"};        // PDG Lnn

  // histogram axes
  ConfigurableAxis rigidityBins{"rigidityBins", {200, -6.f, 6.f}, "Binning for rigidity #it{p}^{TPC}/#it{z}"};
  ConfigurableAxis dEdxBins{"dEdxBins", {1000, 0.f, 1000.f}, "Binning for dE/dx"};
  ConfigurableAxis nSigmaBins{"nSigmaBins", {200, -5.f, 5.f}, "Binning for n sigma"};
  ConfigurableAxis zVtxBins{"zVtxBins", {100, -20.f, 20.f}, "Binning for n sigma"};
  ConfigurableAxis centBins{"centBins", {100, 0.f, 100.f}, "Binning for centrality"};
  ConfigurableAxis TritMomBins{"TritMom", {100, 0.f, 6.f}, "Binning for Triton TPC momentum"};

  // std vector of candidates
  std::vector<lnnCandidate> lnnCandidates;
  // vector to keep track of MC mothers already filled
  std::vector<unsigned int> filledMothers;
  // vector to keep track of the collisions passing the event selection in the MC
  std::vector<bool> isGoodCollision;
  // vector to armazenade h3Track

  Preslice<aod::V0s> perCollision = o2::aod::v0::collisionId;

  HistogramRegistry qaRegistry{"QA", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  float d_bz;
  std::array<float, 6> mBBparams3H;

  // Definiton of histograms to real data [hNsigma3HSelected, hdEdx3HSelected, dEdxtotal, hEVents, hCentFT0(A/C/M) and hCentFV0A] and MC [hDecayChannel, hIsMatterGen, hIsMatterGenTwoBody]
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
    const AxisSpec dEdxAxis{dEdxBins, "d#it{E}/d#it{x}"};
    const AxisSpec nSigma3HAxis{nSigmaBins, "n_{#sigma}({}^{3}H)"};
    const AxisSpec zVtxAxis{zVtxBins, "z_{vtx} (cm)"};
    const AxisSpec centAxis{centBins, "Centrality"};
    const AxisSpec TritMomAxis{TritMomBins, "#it{p}^{TPC}({}^{3}H)"};

    hNsigma3HSel = qaRegistry.add<TH2>("hNsigma3HSel", "; p_{TPC}/z (GeV/#it{c}); n_{#sigma} ({}^{3}H)", HistType::kTH2F, {rigidityAxis, nSigma3HAxis});
    hdEdx3HSel = qaRegistry.add<TH2>("hdEdx3HSel", ";p_{TPC}/z (GeV/#it{c}); dE/dx", HistType::kTH2F, {rigidityAxis, dEdxAxis});
    hdEdx3HTPCMom = qaRegistry.add<TH2>("hdEdx3HTPCMom", "; #it{p}^{TPC}({}^{3}H); dE/dx", HistType::kTH2F, {TritMomAxis, dEdxAxis});
    hdEdxTot = qaRegistry.add<TH2>("hdEdxTot", ";p_{TPC}/z (GeV/#it{c}); dE/dx", HistType::kTH2F, {rigidityAxis, dEdxAxis});
    hEvents = qaRegistry.add<TH1>("hEvents", ";Events; ", HistType::kTH1D, {{2, -0.5, 1.5}});

    hEvents->GetXaxis()->SetBinLabel(1, "All");
    hEvents->GetXaxis()->SetBinLabel(2, "sel8");
    if (doprocessMC) {
      hDecayChannel = qaRegistry.add<TH1>("hDecayChannel", ";Decay channel; ", HistType::kTH1D, {{2, -0.5, 1.5}});
      hDecayChannel->GetXaxis()->SetBinLabel(1, "2-body");
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

  // group BCs
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
      auto h3pid = ccdb->getForTimeStamp<std::array<float, 6>>(pidPath.value + "_3H", run3grp_timestamp);
      std::copy(h3pid->begin(), h3pid->end(), mBBparams3H.begin());
    } else {
      for (int i = 0; i < 5; i++) {
        mBBparams3H[i] = cfgBetheBlochParams->get("3H", Form("p%i", i));
      }
      mBBparams3H[5] = cfgBetheBlochParams->get("3H", "resolution");
    }
    fitter.setBz(d_bz);
    mRunNumber = bc.runNumber();
  }

  // Template
  template <class Tcoll>
  void fillCandidateData(Tcoll const& collision, aod::V0s const& V0s)
  {
    if (mBBparams3H[5] < 0) {
      LOG(fatal) << "Bethe-Bloch parameters for 3H not set, please check your CCDB and configuration";
    }
    for (auto& v0 : V0s) {

      auto posTrack = v0.posTrack_as<TracksFull>();
      auto negTrack = v0.negTrack_as<TracksFull>();

      if (std::abs(posTrack.eta()) > etaMax || std::abs(negTrack.eta()) > etaMax) {
        continue;
      }

      float posRigidity = posTrack.tpcInnerParam();
      float negRigidity = negTrack.tpcInnerParam();

      hdEdxTot->Fill(posRigidity, posTrack.tpcSignal());
      hdEdxTot->Fill(-negRigidity, negTrack.tpcSignal());

      // Bethe-Bloch calcution for 3H & nSigma calculation
      double expBethePos{tpc::BetheBlochAleph(static_cast<float>(posRigidity / constants::physics::MassTriton), mBBparams3H[0], mBBparams3H[1], mBBparams3H[2], mBBparams3H[3], mBBparams3H[4])};
      double expBetheNeg{tpc::BetheBlochAleph(static_cast<float>(negRigidity / constants::physics::MassTriton), mBBparams3H[0], mBBparams3H[1], mBBparams3H[2], mBBparams3H[3], mBBparams3H[4])};
      double expSigmaPos{expBethePos * mBBparams3H[5]};
      double expSigmaNeg{expBetheNeg * mBBparams3H[5]};
      auto nSigmaTPCpos = static_cast<float>((posTrack.tpcSignal() - expBethePos) / expSigmaPos);
      auto nSigmaTPCneg = static_cast<float>((negTrack.tpcSignal() - expBetheNeg) / expSigmaNeg);

      // ITS only tracks do not have TPC information. TPCnSigma: only lower cut to allow for triton reconstruction
      bool is3H = posTrack.hasTPC() && nSigmaTPCpos > -1 * nSigmaMax3H;
      bool isAnti3H = negTrack.hasTPC() && nSigmaTPCneg > -1 * nSigmaMax3H;

      if (!is3H && !isAnti3H)
        continue;

      // Describing lnn as matter candidate
      lnnCandidate lnnCand;
      lnnCand.isMatter = is3H && isAnti3H ? std::abs(nSigmaTPCpos) < std::abs(nSigmaTPCneg) : is3H;
      auto& h3track = lnnCand.isMatter ? posTrack : negTrack;
      auto& h3Rigidity = lnnCand.isMatter ? posRigidity : negRigidity;
      if (h3track.tpcNClsFound() < nTPCClusMin3H || h3Rigidity < TPCRigidityMin3H) {
        continue;
      }

      lnnCand.nSigma3H = lnnCand.isMatter ? nSigmaTPCpos : nSigmaTPCneg;
      lnnCand.nTPCClusters3H = lnnCand.isMatter ? posTrack.tpcNClsFound() : negTrack.tpcNClsFound();
      lnnCand.tpcSignal3H = lnnCand.isMatter ? posTrack.tpcSignal() : negTrack.tpcSignal();
      lnnCand.clusterSizeITS3H = lnnCand.isMatter ? posTrack.itsClusterSizes() : negTrack.itsClusterSizes();
      lnnCand.nTPCClustersPi = !lnnCand.isMatter ? posTrack.tpcNClsFound() : negTrack.tpcNClsFound();
      lnnCand.tpcSignalPi = !lnnCand.isMatter ? posTrack.tpcSignal() : negTrack.tpcSignal();
      lnnCand.clusterSizeITSPi = !lnnCand.isMatter ? posTrack.itsClusterSizes() : negTrack.itsClusterSizes();
      lnnCand.mom3HTPC = lnnCand.isMatter ? posRigidity : negRigidity;
      lnnCand.momPiTPC = !lnnCand.isMatter ? posRigidity : negRigidity;

      lnnCand.flags |= lnnCand.isMatter ? static_cast<uint8_t>((posTrack.pidForTracking() & 0xF) << 4) : static_cast<uint8_t>((negTrack.pidForTracking() & 0xF) << 4);
      lnnCand.flags |= lnnCand.isMatter ? static_cast<uint8_t>(negTrack.pidForTracking() & 0xF) : static_cast<uint8_t>(posTrack.pidForTracking() & 0xF);

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

      auto& h3PropTrack = lnnCand.isMatter ? fitter.getTrack(0) : fitter.getTrack(1);
      auto& piPropTrack = lnnCand.isMatter ? fitter.getTrack(1) : fitter.getTrack(0);
      h3PropTrack.getPxPyPzGlo(lnnCand.mom3H);
      piPropTrack.getPxPyPzGlo(lnnCand.momPi);

      // Definition of relativistic momentum and energy to triton and pion and total energy
      float h3P2 = lnnCand.mom3H[0] * lnnCand.mom3H[0] + lnnCand.mom3H[1] * lnnCand.mom3H[1] + lnnCand.mom3H[2] * lnnCand.mom3H[2];
      float piP2 = lnnCand.momPi[0] * lnnCand.momPi[0] + lnnCand.momPi[1] * lnnCand.momPi[1] + lnnCand.momPi[2] * lnnCand.momPi[2];
      float h3E = std::sqrt(h3P2 + h3Mass * h3Mass);
      float piE = std::sqrt(piP2 + piMass * piMass);
      float h3lE = h3E + piE;

      // Building the mother particle: lnn
      std::array<float, 3> lnnMom;

      const auto& vtx = fitter.getPCACandidate();
      for (int i = 0; i < 3; i++) {
        lnnCand.decVtx[i] = vtx[i];
        lnnMom[i] = lnnCand.mom3H[i] + lnnCand.momPi[i];
      }

      float lnnPt = std::hypot(lnnMom[0], lnnMom[1]);
      if (lnnPt < ptMin) {
        continue;
      }

      // Definition of lnn mass
      float mLNN_HypHI = 2.994; // value in GeV, but 2993.7 MeV/c**2

      float massLNNL = std::sqrt(h3lE * h3lE - lnnMom[0] * lnnMom[0] - lnnMom[1] * lnnMom[1] - lnnMom[2] * lnnMom[2]);
      bool isLNNMass = false;
      if (massLNNL > mLNN_HypHI - masswidth && massLNNL < mLNN_HypHI + masswidth) {
        isLNNMass = true;
      }
      if (!isLNNMass) {
        continue;
      }

      // V0, primary vertex and poiting angle
      lnnCand.dcaV0dau = std::sqrt(fitter.getChi2AtPCACandidate());
      if (lnnCand.dcaV0dau > dcav0dau) {
        continue;
      }

      std::array<float, 3> primVtx = {collision.posX(), collision.posY(), collision.posZ()};

      double cosPA = RecoDecay::cpa(primVtx, lnnCand.decVtx, lnnMom);
      if (cosPA < v0cospa) {
        continue;
      }

      for (int i = 0; i < 3; i++) {
        lnnCand.decVtx[i] = lnnCand.decVtx[i] - primVtx[i];
      }

      // if survived all selections, propagate decay daughters to PV
      gpu::gpustd::array<float, 2> dcaInfo;

      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, posTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
      lnnCand.isMatter ? lnnCand.h3DCAXY = dcaInfo[0] : lnnCand.piDCAXY = dcaInfo[0];

      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, negTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
      lnnCand.isMatter ? lnnCand.piDCAXY = dcaInfo[0] : lnnCand.h3DCAXY = dcaInfo[0];

      // finally, push back the candidate
      lnnCand.isReco = true;
      lnnCand.posTrackID = posTrack.globalIndex();
      lnnCand.negTrackID = negTrack.globalIndex();

      int chargeFactor = -1 + 2 * lnnCand.isMatter;
      hdEdx3HSel->Fill(chargeFactor * lnnCand.mom3HTPC, h3track.tpcSignal());
      hNsigma3HSel->Fill(chargeFactor * lnnCand.mom3HTPC, lnnCand.nSigma3H);
      lnnCandidates.push_back(lnnCand);

      if (is3H) {
        hdEdx3HTPCMom->Fill(lnnCand.mom3HTPC, h3track.tpcSignal());
      }
    }
  }

  // Monte Carlo information
  void fillMCinfo(aod::McTrackLabels const& trackLabels, aod::McParticles const&)
  {
    for (auto& lnnCand : lnnCandidates) {
      auto mcLabPos = trackLabels.rawIteratorAt(lnnCand.posTrackID);
      auto mcLabNeg = trackLabels.rawIteratorAt(lnnCand.negTrackID);

      // Checking lnn, tritons and pions with MC simulations
      if (mcLabPos.has_mcParticle() && mcLabNeg.has_mcParticle()) {
        auto mcTrackPos = mcLabPos.mcParticle_as<aod::McParticles>();
        auto mcTrackNeg = mcLabNeg.mcParticle_as<aod::McParticles>();
        if (mcTrackPos.has_mothers() && mcTrackNeg.has_mothers()) {
          for (auto& negMother : mcTrackNeg.mothers_as<aod::McParticles>()) {
            for (auto& posMother : mcTrackPos.mothers_as<aod::McParticles>()) {
              if (posMother.globalIndex() != negMother.globalIndex())
                continue;
              if (!((mcTrackPos.pdgCode() == h3DauPdg && mcTrackNeg.pdgCode() == -211) || (mcTrackPos.pdgCode() == 211 && mcTrackNeg.pdgCode() == -1 * h3DauPdg)))
                continue;
              if (std::abs(posMother.pdgCode()) != lnnPdg)
                continue;

              // Checking primary and second vertex with MC simulations
              std::array<float, 3> posPrimVtx = {posMother.vx(), posMother.vy(), posMother.vz()};

              std::array<float, 3> secVtx = {mcTrackPos.vx(), mcTrackPos.vy(), mcTrackPos.vz()};

              lnnCand.gMom = posMother.pVector();

              lnnCand.gMom3H = mcTrackPos.pdgCode() == h3DauPdg ? mcTrackPos.pVector() : mcTrackNeg.pVector();

              for (int i = 0; i < 3; i++) {
                lnnCand.gDecVtx[i] = secVtx[i] - posPrimVtx[i];
              }
              lnnCand.isSignal = true;
              lnnCand.pdgCode = posMother.pdgCode();
              lnnCand.survEvSelection = isGoodCollision[posMother.mcCollisionId()];

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
      lnnCandidates.clear();

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      hEvents->Fill(0.);
      if ((!collision.sel8()) || std::abs(collision.posZ()) > 10) {
        continue;
      }
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

      for (auto& lnnCand : lnnCandidates) {
        outputDataTable(collision.centFT0A(), collision.centFT0C(), collision.centFT0M(),
                        collision.posX(), collision.posY(), collision.posZ(),
                        lnnCand.isMatter,
                        lnnCand.recoPt3H(), lnnCand.recoPhi3H(), lnnCand.recoEta3H(),
                        lnnCand.recoPtPi(), lnnCand.recoPhiPi(), lnnCand.recoEtaPi(),
                        lnnCand.decVtx[0], lnnCand.decVtx[1], lnnCand.decVtx[2],
                        lnnCand.dcaV0dau, lnnCand.h3DCAXY, lnnCand.piDCAXY,
                        lnnCand.nSigma3H, lnnCand.nTPCClusters3H, lnnCand.nTPCClustersPi,
                        lnnCand.mom3HTPC, lnnCand.momPiTPC, lnnCand.tpcSignal3H, lnnCand.tpcSignalPi,
                        lnnCand.clusterSizeITS3H, lnnCand.clusterSizeITSPi, lnnCand.flags);
      }
    }
  }
  PROCESS_SWITCH(lnnRecoTask, processData, "Data analysis", true);

  // MC process
  void processMC(CollisionsFullMC const& collisions, aod::McCollisions const& mcCollisions, aod::V0s const& V0s, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC)
  {
    filledMothers.clear();

    isGoodCollision.clear();
    isGoodCollision.resize(mcCollisions.size(), false);

    for (const auto& collision : collisions) {
      lnnCandidates.clear();
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      hEvents->Fill(0.);

      if ((collision.posZ()) > 10) {
        continue;
      }
      hEvents->Fill(1.);
      hZvtx->Fill(collision.posZ());
      hCentFT0A->Fill(collision.centFT0A());
      hCentFT0C->Fill(collision.centFT0C());
      hCentFT0M->Fill(collision.centFT0M());
      hCentFV0A->Fill(collision.centFV0A());

      if (collision.has_mcCollision()) {
        isGoodCollision[collision.mcCollisionId()] = true;
      }

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollision, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      fillCandidateData(collision, V0Table_thisCollision);
      fillMCinfo(trackLabelsMC, particlesMC);

      for (auto& lnnCand : lnnCandidates) {
        if (!lnnCand.isSignal && mcSignalOnly) {
          continue;
        }
        int chargeFactor = -1 + 2 * (lnnCand.pdgCode > 0);
        outputMCTable(collision.centFT0A(), collision.centFT0C(), collision.centFT0M(),
                      collision.posX(), collision.posY(), collision.posZ(),
                      lnnCand.isMatter,
                      lnnCand.recoPt3H(), lnnCand.recoPhi3H(), lnnCand.recoEta3H(),
                      lnnCand.recoPtPi(), lnnCand.recoPhiPi(), lnnCand.recoEtaPi(),
                      lnnCand.decVtx[0], lnnCand.decVtx[1], lnnCand.decVtx[2],
                      lnnCand.dcaV0dau, lnnCand.h3DCAXY, lnnCand.piDCAXY,
                      lnnCand.nSigma3H, lnnCand.nTPCClusters3H, lnnCand.nTPCClustersPi,
                      lnnCand.mom3HTPC, lnnCand.momPiTPC, lnnCand.tpcSignal3H, lnnCand.tpcSignalPi,
                      lnnCand.clusterSizeITS3H, lnnCand.clusterSizeITSPi, lnnCand.flags,
                      chargeFactor * lnnCand.genPt(), lnnCand.genPhi(), lnnCand.genEta(), lnnCand.genPt3H(),
                      lnnCand.gDecVtx[0], lnnCand.gDecVtx[1], lnnCand.gDecVtx[2], lnnCand.isReco, lnnCand.isSignal, lnnCand.survEvSelection);
      }
    }

    // now we fill only the signal candidates that were not reconstructed
    for (auto& mcPart : particlesMC) {

      if (std::abs(mcPart.pdgCode()) != lnnPdg) {
        continue;
      }
      std::array<float, 3> secVtx;
      std::array<float, 3> primVtx = {mcPart.vx(), mcPart.vy(), mcPart.vz()};

      std::array<float, 3> momMother = mcPart.pVector();

      std::array<float, 3> mom3H;
      bool is3HFound = false;
      for (auto& mcDaught : mcPart.daughters_as<aod::McParticles>()) {
        if (std::abs(mcDaught.pdgCode()) == h3DauPdg) {
          secVtx = {mcDaught.vx(), mcDaught.vy(), mcDaught.vz()};

          mom3H = mcDaught.pVector();

          is3HFound = true;
          break;
        }
      }
      if (mcPart.pdgCode() > 0) {
        hIsMatterGen->Fill(0.);
      } else {
        hIsMatterGen->Fill(1.);
      }
      if (!is3HFound) {
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

      lnnCandidate lnnCand;
      lnnCand.pdgCode = mcPart.pdgCode();
      lnnCand.survEvSelection = isGoodCollision[mcPart.mcCollisionId()];
      int chargeFactor = -1 + 2 * (lnnCand.pdgCode > 0);
      for (int i = 0; i < 3; i++) {
        lnnCand.gDecVtx[i] = secVtx[i] - primVtx[i];
        lnnCand.gMom[i] = momMother[i];
        lnnCand.gMom3H[i] = mom3H[i];
      }
      lnnCand.posTrackID = -1;
      lnnCand.negTrackID = -1;
      lnnCand.isSignal = true;
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
                    chargeFactor * lnnCand.genPt(), lnnCand.genPhi(), lnnCand.genEta(), lnnCand.genPt3H(),
                    lnnCand.gDecVtx[0], lnnCand.gDecVtx[1], lnnCand.gDecVtx[2], lnnCand.isReco, lnnCand.isSignal, lnnCand.survEvSelection);
    }
  }
  PROCESS_SWITCH(lnnRecoTask, processMC, "MC analysis", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lnnRecoTask>(cfgc)};
}
