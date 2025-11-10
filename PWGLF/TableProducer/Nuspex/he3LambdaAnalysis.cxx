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

#include "PWGLF/DataModel/LFSlimHeLambda.h"

#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <CCDB/BasicCCDBManager.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DataFormatsTPC/BetheBlochAleph.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <Math/Vector4D.h>
#include <TRandom3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
namespace
{
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleName{"He3"};
o2::base::MatLayerCylSet* matLUT = nullptr;

float alphaAP(std::array<float, 3> const& momA, std::array<float, 3> const& momB, std::array<float, 3> const& momC)
{
  const float lQlPos = (momB[0] * momA[0] + momB[1] * momA[1] + momB[2] * momA[2]);
  const float lQlNeg = (momC[0] * momA[0] + momC[1] * momA[1] + momC[2] * momA[2]);
  return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

float qtAP(std::array<float, 3> const& momA, std::array<float, 3> const& momB)
{
  const float dp = momA[0] * momB[0] + momA[1] * momB[1] + momA[2] * momB[2];
  const float p2A = momA[0] * momA[0] + momA[1] * momA[1] + momA[2] * momA[2];
  const float p2B = momB[0] * momB[0] + momB[1] * momB[1] + momB[2] * momB[2];
  return std::sqrt(p2B - dp * dp / p2A);
}

std::shared_ptr<TH2> hTPCsignalAll;
std::shared_ptr<TH2> hTPCsignalHe3;
std::shared_ptr<TH2> hTPCnSigmaAll;
std::shared_ptr<TH2> hTPCnSigmaHe3;
std::shared_ptr<TH2> hArmenterosPodolanskiAll;
std::shared_ptr<TH2> hArmenterosPodolanskiSelected;
std::shared_ptr<TH2> hInvariantMassUS;
std::shared_ptr<TH2> hInvariantMassLS;

}; // namespace

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPr, aod::pidTPCFullPi>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms>;

struct he3LambdaAnalysis {

  // Services
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  o2::vertexing::DCAFitterN<2> fitter;

  Produces<o2::aod::LFEvents> lfHe3V0Collision;
  Produces<o2::aod::LFHe3_001> lfHe3;
  Produces<o2::aod::LFLambda_001> lfLambda;

  // Configurables for event selection
  struct : ConfigurableGroup {
    std::string prefix = "cfgEventSelection";
    Configurable<float> zVertexMax{"zVertexMax", 10.0f, "Accepted z-vertex range"};
    Configurable<bool> useSel8{"useSel8", true, "Use Sel8 event selection"};
    Configurable<bool> skimmedProcessing{"skimmedProcessing", false, "Skimmed dataset processing"};
  } cfgEventSelection;

  // He3 selection criteria
  struct : ConfigurableGroup {
    std::string prefix = "cfgHe3";
    Configurable<float> ptMin{"ptMin", 1.0f, "Minimum He3 pT"};
    Configurable<float> ptMax{"ptMax", 10.0f, "Maximum He3 pT"};
    Configurable<float> etaMax{"etaMax", 0.9f, "Maximum He3 pseudorapidity"};
    Configurable<float> minTPCrigidity{"minTPCrigidity", 0.5f, "Minimum He3 rigidity"};
    Configurable<float> nSigmaTPCMax{"nSigmaTPCMax", 4.0f, "Maximum He3 TPC nSigma"};
    Configurable<float> dcaxyMax{"dcaxyMax", 0.5f, "Maximum He3 DCA xy"};
    Configurable<float> dcazMax{"dcazMax", 0.5f, "Maximum He3 DCA z"};
    Configurable<int> tpcClusMin{"tpcClusMin", 100, "Minimum He3 TPC clusters"};
    Configurable<int> itsClusMin{"itsClusMin", 5, "Minimum He3 ITS clusters"};
    Configurable<LabeledArray<double>> betheBlochParams{"betheBlochParams", {betheBlochDefault[0], 1, 6, particleName, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};
  } cfgHe3;

  // Lambda selection criteria
  struct : ConfigurableGroup {
    std::string prefix = "cfgLambda";
    Configurable<float> ptMin{"ptMin", 0.5f, "Minimum Lambda pT"};
    Configurable<float> ptMax{"ptMax", 10.0f, "Maximum Lambda pT"};
    Configurable<float> massWindow{"massWindow", 0.015f, "Lambda mass window"};
    Configurable<float> cosPAMin{"cosPAMin", 0.99f, "Minimum Lambda cosPA"};
    Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.5f, "Maximum Lambda DCA V0 daughters"};
    Configurable<float> v0RadiusMin{"v0RadiusMin", 0.5f, "Minimum Lambda V0 radius"};
    Configurable<float> v0RadiusMax{"v0RadiusMax", 35.0f, "Maximum Lambda V0 radius"};
    Configurable<int> tpcNClsMin{"tpcNClsMin", 70, "Minimum TPC clusters for Lambda daughters"};
    Configurable<float> protonNSigmaTPCMax{"protonNSigmaTPCMax", 4.0f, "Maximum proton TPC nSigma"};
    Configurable<float> pionNSigmaTPCMax{"pionNSigmaTPCMax", 4.0f, "Maximum pion TPC nSigma"};
  } cfgLambda;

  // Pair selection criteria
  struct : ConfigurableGroup {
    std::string prefix = "cfgPair";
    Configurable<float> ptMin{"PtMin", 1.0f, "Minimum pair pT"};
    Configurable<float> ptMax{"PtMax", 20.0f, "Maximum pair pT"};
    Configurable<float> rapidityMax{"RapidityMax", 0.5f, "Maximum pair rapidity"};
  } cfgPair;

  // CCDB options
  struct : ConfigurableGroup {
    std::string prefix = "ccdb";
    Configurable<std::string> url{"url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  } ccdbOptions;

  std::array<double, 6> mBBparamsHe;
  float mBz = 0.0f; // Magnetic field in T
  HistogramRegistry mRegistry{"He3LambdaAnalysis"};
  int mRunNumber = 0; // Current run number

  void init(InitContext const&)
  {
    // Initialize CCDB
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(true);

    for (int i = 0; i < 5; i++) {
      mBBparamsHe[i] = cfgHe3.betheBlochParams->get("He3", Form("p%i", i));
    }
    mBBparamsHe[5] = cfgHe3.betheBlochParams->get("He3", "resolution");
    matLUT = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    fitter.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrLUT);

    zorroSummary.setObject(zorro.getZorroSummary());

    mRegistry.add("hEventSelection", "Event Selection", {HistType::kTH1L, {{6, -.5, 5.5}}});
    std::vector<std::string> labels{"Total Events", "Sel8 Events", "Z-Vertex OK", "Additional Event Selections", "He3 Candidates Found", "He3 and Lambda Candidates Found"};
    for (size_t i = 1; i <= labels.size(); ++i) {
      mRegistry.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(i, labels[i - 1].c_str());
    }

    mRegistry.add("hCentralityAll", "Centrality All", {HistType::kTH1L, {{100, 0., 100.}}});
    mRegistry.add("hCentralitySelected", "Centrality Selected", {HistType::kTH1L, {{100, 0., 100.}}});

    hTPCsignalAll = mRegistry.add<TH2>("hTPCsignalAll", "TPC Signal All", {HistType::kTH2D, {{400, -10, 10}, {1000, 0, 2000}}});
    hTPCsignalHe3 = mRegistry.add<TH2>("hTPCsignalHe3", "TPC Signal He3", {HistType::kTH2D, {{400, -10, 10}, {1000, 0, 2000}}});

    hTPCnSigmaAll = mRegistry.add<TH2>("hTPCnSigmaAll", "TPC nSigma All", {HistType::kTH2D, {{400, -10, 10}, {100, -5., 5.}}});
    hTPCnSigmaHe3 = mRegistry.add<TH2>("hTPCnSigmaHe3", "TPC nSigma He3", {HistType::kTH2D, {{400, -10, 10}, {100, -5., 5.}}});

    hArmenterosPodolanskiAll = mRegistry.add<TH2>("hArmenterosPodolanskiAll", "Armenteros-Podolanski All", {HistType::kTH2D, {{100, -1., 1.}, {100, 0., 0.5}}});
    hArmenterosPodolanskiSelected = mRegistry.add<TH2>("hArmenterosPodolanskiSelected", "Armenteros-Podolanski Selected", {HistType::kTH2D, {{100, -1., 1.}, {100, 0., 0.5}}});

    constexpr double ConstituentsMass = o2::constants::physics::MassProton + o2::constants::physics::MassNeutron * 2 + o2::constants::physics::MassSigmaPlus;
    hInvariantMassUS = mRegistry.add<TH2>("hInvariantMassUS", "Invariant Mass", {HistType::kTH2D, {{45, 1., 10}, {100, ConstituentsMass - 0.05, ConstituentsMass + 0.05}}});
    hInvariantMassLS = mRegistry.add<TH2>("hInvariantMassLS", "Invariant Mass", {HistType::kTH2D, {{45, 1., 10}, {100, ConstituentsMass - 0.05, ConstituentsMass + 0.05}}});

    LOGF(info, "He3-Lambda analysis initialized");
  }

  void initCCDB(const auto& bc)
  {
    int runNumber = bc.runNumber();
    if (runNumber == mRunNumber) {
      return; // Already initialized for this run
    }
    mRunNumber = runNumber;
    if (cfgEventSelection.skimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fHe");
      zorro.populateHistRegistry(mRegistry, bc.runNumber());
    }
    o2::parameters::GRPMagField* grpmag = ccdb->getForRun<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", runNumber);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    mBz = static_cast<float>(grpmag->getNominalL3Field());
    fitter.setBz(mBz);
    o2::base::Propagator::Instance()->setMatLUT(matLUT);
  }

  void processData(CollisionsFull::iterator const& collision,
                   TracksFull const& tracks,
                   aod::V0s const& v0s,
                   aod::BCsWithTimestamps const&)
  {
    const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    mRegistry.get<TH1>(HIST("hEventSelection"))->Fill(0); // Total events
    mRegistry.get<TH1>(HIST("hCentralityAll"))->Fill(collision.centFT0C());
    if (cfgEventSelection.useSel8 && !collision.sel8()) {
      return; // Skip events not passing Sel8 selection
    }
    mRegistry.get<TH1>(HIST("hEventSelection"))->Fill(1); // Sel8 events
    if (std::abs(collision.posZ()) > cfgEventSelection.zVertexMax) {
      return; // Skip events with z-vertex outside range
    }
    mRegistry.get<TH1>(HIST("hEventSelection"))->Fill(2); // Z-vertex OK

    // Additional event selections not implemented, but can be added here
    if (cfgEventSelection.skimmedProcessing) {
      if (!zorro.isSelected(bc.globalBC())) {
        return; // Skip events not passing Zorro selection
      }
    }
    mRegistry.get<TH1>(HIST("hEventSelection"))->Fill(3); // Additional event selections

    // Process He3 candidates
    std::vector<he3Candidate> he3Candidates;
    o2::track::TrackParCov trackParCov;
    trackParCov.setPID(o2::track::PID::Helium3);
    const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};

    for (auto const& track : tracks) {
      if (track.tpcNClsFound() < cfgHe3.tpcClusMin || track.itsNCls() < cfgHe3.itsClusMin) {
        continue; // Skip tracks with insufficient clusters
      }
      hTPCsignalAll->Fill(track.tpcInnerParam() * track.sign(), track.tpcSignal());
      const float pt = track.pt() * 2.0f;
      float expTPCSignal = o2::tpc::BetheBlochAleph(track.tpcInnerParam() * 2.0f / constants::physics::MassHelium3, mBBparamsHe[0], mBBparamsHe[1], mBBparamsHe[2], mBBparamsHe[3], mBBparamsHe[4]);
      double nSigmaTPC = (track.tpcSignal() - expTPCSignal) / (expTPCSignal * mBBparamsHe[5]);
      hTPCnSigmaAll->Fill(track.tpcInnerParam() * track.sign(), nSigmaTPC);
      if (pt < cfgHe3.ptMin || pt > cfgHe3.ptMax || std::abs(track.eta()) > cfgHe3.etaMax || track.tpcInnerParam() < cfgHe3.minTPCrigidity || std::abs(nSigmaTPC) > cfgHe3.nSigmaTPCMax) {
        continue; // Skip tracks outside He3 PID+kinematics selection criteria
      }
      setTrackParCov(track, trackParCov);
      std::array<float, 2> dcaInfo;
      o2::base::Propagator::Instance()->propagateToDCA(collVtx, trackParCov, mBz, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrLUT, &dcaInfo);
      if (std::abs(dcaInfo[0]) > cfgHe3.dcaxyMax || std::abs(dcaInfo[1]) > cfgHe3.dcazMax) {
        continue; // Skip tracks with DCA outside range
      }
      hTPCsignalHe3->Fill(track.tpcInnerParam() * track.sign(), track.tpcSignal());
      hTPCnSigmaHe3->Fill(track.tpcInnerParam() * track.sign(), nSigmaTPC);
      he3Candidate candidate;
      candidate.momentum.SetCoordinates(track.pt() * 2.0f, track.eta(), track.phi(), o2::constants::physics::MassHelium3);
      candidate.nSigmaTPC = nSigmaTPC;
      candidate.dcaXY = dcaInfo[0];
      candidate.dcaZ = dcaInfo[1];
      candidate.tpcNClsFound = track.tpcNClsFound();
      candidate.tpcNClsPID = track.tpcNClsPID();
      candidate.itsNCls = track.itsNCls();
      candidate.itsClusterSizes = track.itsClusterSizes();
      candidate.sign = track.sign() > 0 ? 1 : -1;
      he3Candidates.push_back(candidate);
    }
    if (he3Candidates.empty()) {
      return; // No valid He3 candidates found
    }
    mRegistry.get<TH1>(HIST("hEventSelection"))->Fill(4); // He3 candidates found

    // Process Lambda candidates
    std::vector<lambdaCandidate> lambdaCandidates;
    for (auto const& v0 : v0s) {
      if (v0.v0Type() != 1) {
        continue;
      }
      const auto posTrack = v0.posTrack_as<TracksFull>();
      const auto negTrack = v0.negTrack_as<TracksFull>();

      if (posTrack.tpcNClsFound() < cfgLambda.tpcNClsMin || negTrack.tpcNClsFound() < cfgLambda.tpcNClsMin) {
        continue; // Skip V0s with insufficient TPC clusters
      }
      auto trackParPos = getTrackParCov(posTrack);
      auto trackParNeg = getTrackParCov(negTrack);
      int nCand = 0;
      try {
        nCand = fitter.process(trackParPos, trackParNeg);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        return;
      }
      if (nCand == 0) {
        continue;
      }
      auto& propParPos = fitter.getTrack(0);
      auto& propParNeg = fitter.getTrack(1);
      std::array<float, 3> momPos, momNeg;
      propParPos.getPxPyPzGlo(momPos);
      propParNeg.getPxPyPzGlo(momNeg);
      const std::array<float, 3> momV0{momPos[0] + momNeg[0], momPos[1] + momNeg[1], momPos[2] + momNeg[2]};
      float alpha = alphaAP(momV0, momPos, momNeg);
      float qt = qtAP(momV0, momPos);
      hArmenterosPodolanskiAll->Fill(alpha, qt);

      bool matter = alpha > 0;
      const auto& protonTrack = matter ? posTrack : negTrack;
      const auto& pionTrack = matter ? negTrack : posTrack;
      const auto& protonMom = matter ? momPos : momNeg;
      const auto& pionMom = matter ? momNeg : momPos;

      if (std::abs(protonTrack.tpcNSigmaPr()) > cfgLambda.protonNSigmaTPCMax ||
          std::abs(pionTrack.tpcNSigmaPi()) > cfgLambda.pionNSigmaTPCMax) {
        continue; // Skip V0s with TPC nSigma outside range
      }
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>> protonMom4D(protonMom[0], protonMom[1], protonMom[2], o2::constants::physics::MassProton);
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>> pionMom4D(pionMom[0], pionMom[1], pionMom[2], o2::constants::physics::MassPionCharged);
      auto lambdaMom4D = protonMom4D + pionMom4D;
      float massLambda = lambdaMom4D.M();

      if (std::abs(massLambda - o2::constants::physics::MassLambda0) > cfgLambda.massWindow) {
        continue; // Skip V0s outside mass window
      }
      hArmenterosPodolanskiSelected->Fill(alpha, qt);

      std::array<float, 2> dcaInfoProton, dcaInfoPion;
      o2::base::Propagator::Instance()->propagateToDCA(collVtx, matter ? trackParPos : trackParNeg, mBz, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrLUT, &dcaInfoProton);
      o2::base::Propagator::Instance()->propagateToDCA(collVtx, matter ? trackParNeg : trackParPos, mBz, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrLUT, &dcaInfoPion);

      const auto sv = fitter.getPCACandidate(0);

      lambdaCandidate candidate;
      candidate.momentum.SetCoordinates(lambdaMom4D.Pt(), lambdaMom4D.Eta(), lambdaMom4D.Phi(), o2::constants::physics::MassLambda0);
      candidate.mass = massLambda;
      candidate.cosPA = (sv[0] - collVtx.x()) * lambdaMom4D.Px() +
                        (sv[1] - collVtx.y()) * lambdaMom4D.Py() +
                        (sv[2] - collVtx.z()) * lambdaMom4D.Pz();
      candidate.cosPA /= std::hypot(sv[0] - collVtx.x(), sv[1] - collVtx.y(), sv[2] - collVtx.z()) * lambdaMom4D.P();
      candidate.dcaV0Daughters = std::sqrt(fitter.getChi2AtPCACandidate(0));
      candidate.dcaProtonToPV = std::hypot(dcaInfoProton[0], dcaInfoProton[1]);
      candidate.dcaPionToPV = std::hypot(dcaInfoPion[0], dcaInfoPion[1]);
      candidate.v0Radius = std::hypot(sv[0], sv[1]);
      candidate.protonNSigmaTPC = protonTrack.tpcNSigmaPr();
      candidate.pionNSigmaTPC = pionTrack.tpcNSigmaPi();
      candidate.sign = matter ? 1 : -1; // Positive sign for Lambda, negative for anti-Lambda
      lambdaCandidates.push_back(candidate);
    }
    if (lambdaCandidates.empty()) {
      return; // No valid Lambda candidates found
    }
    mRegistry.get<TH1>(HIST("hEventSelection"))->Fill(5); // He3 and Lambda candidates found
    mRegistry.get<TH1>(HIST("hCentralitySelected"))->Fill(collision.centFT0C());

    // Fill output tables
    lfHe3V0Collision(collision.posZ(), collision.centFT0C());
    for (const auto& he3 : he3Candidates) {
      lfHe3(lfHe3V0Collision.lastIndex(), he3.momentum.Pt(), he3.momentum.Eta(), he3.momentum.Phi(),
            he3.dcaXY, he3.dcaZ, he3.tpcNClsFound, he3.tpcNClsPID, he3.itsClusterSizes, he3.nSigmaTPC, he3.sign);
    }
    for (const auto& lambda : lambdaCandidates) {
      lfLambda(lfHe3V0Collision.lastIndex(), lambda.momentum.Pt(), lambda.momentum.Eta(), lambda.momentum.Phi(),
               lambda.mass, lambda.cosPA, lambda.dcaV0Daughters, lambda.dcaProtonToPV, lambda.dcaPionToPV, lambda.v0Radius, lambda.protonNSigmaTPC, lambda.pionNSigmaTPC, lambda.sign);
    }

    for (const auto& he3 : he3Candidates) {
      for (const auto& lambda : lambdaCandidates) {
        auto pairMomentum = lambda.momentum + he3.momentum; // Calculate invariant mass
        (he3.sign * lambda.sign > 0 ? hInvariantMassLS : hInvariantMassUS)->Fill(pairMomentum.Pt(), pairMomentum.M());
      }
    }
  }
  PROCESS_SWITCH(he3LambdaAnalysis, processData, "Process data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<he3LambdaAnalysis>(cfgc)};
}
