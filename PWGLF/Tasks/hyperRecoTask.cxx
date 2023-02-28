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
// =====================
//
//
//
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DCAFitter/DCAFitterN.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;

namespace
{
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNames{"He3"};
// float kHyperMass = 2.99131;
// float kHyperPDG = 1010010030;

} // namespace

struct hyperCandidate {
  int posTrackID;
  int negTrackID;
  bool isMatter;
  std::array<float, 3> mom;
  std::array<float, 3> decVtx;
  float dcaV0dau = -1;
  float cosPA;
  float rapidity;
  float nSigmaHe3;
  float nTPCClustersHe3;
  float he3DCAXY;
  float piDCAXY;

  float gCt;
  std::array<float, 3> gMom;

  bool isSignal=false; // true MC signal
  bool isReco=false;  // true if the candidate is actually reconstructed

};

struct hyperRecoTask {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Selection criteria
  Configurable<double> v0cospa{"hypcospa", 0.95, "V0 CosPA"}; // very open, please
  Configurable<float> dcav0dau{"hypdcaDau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dca3hetopv{"hedcatopv", -1, "DCA helium To PV"};
  Configurable<float> dcapiontopv{"pidcatopv", .1, "DCA pion To PV"};
  Configurable<float> v0radius{"hypradius", 1.0, "hyp radius"};
  Configurable<float> rapidity{"rapidity", 0.8, "rapidity"};
  Configurable<float> heliumNsigmaMax{"heliumNsigmaMax", 8, "helium dEdx cut (n sigma)"};

  // Define o2 fitter, 2-prong, active memory (no need to redefine per event)
  o2::vertexing::DCAFitterN<2> fitter;

  // bethe bloch parameters
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};

  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};

  o2::track::TrackPar lHyTrack;

  // CCDB options
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // std vector of candidates
  std::vector<hyperCandidate> hyperCandidates;

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
    // fitter.setWeightedFinalPCA(d_UseWeightedPCA);
    int mat{static_cast<int>(cfgMaterialCorrection)};
    fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));
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
    fitter.setBz(d_bz);
    mRunNumber = bc.runNumber();
  }

  int mRunNumber;
  float d_bz;

  void fillCandidateData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0s const& V0s, TracksFull const& tracks)
  {
    for (auto& v0 : V0s) {

      LOG(info) << "collision: " << collision.globalIndex();
      LOG(info) << "V0 collision ID: " << v0.collisionId();
      auto posTrack = v0.posTrack_as<TracksFull>();
      auto negTrack = v0.negTrack_as<TracksFull>();
      LOG(info) << "pos track: "
                << " pt: " << posTrack.pt();

      LOG(info) << "neg track: "
                << " pt: " << negTrack.pt();

      double expBethePos{tpc::BetheBlochAleph(static_cast<double>(posTrack.tpcInnerParam()), cfgBetheBlochParams->get("He3", 0u), cfgBetheBlochParams->get("He3", 1u), cfgBetheBlochParams->get("He3", 2u), cfgBetheBlochParams->get("He3", 3u), cfgBetheBlochParams->get("He3", 4u))};

      double expBetheNeg{tpc::BetheBlochAleph(static_cast<double>(negTrack.tpcInnerParam()), cfgBetheBlochParams->get("He3", 0u), cfgBetheBlochParams->get("He3", 1u), cfgBetheBlochParams->get("He3", 2u), cfgBetheBlochParams->get("He3", 3u), cfgBetheBlochParams->get("He3", 4u))};

      double expSigmaPos{expBethePos * cfgBetheBlochParams->get("He3", 5u)};
      double expSigmaNeg{expBetheNeg * cfgBetheBlochParams->get("He3", 5u)};
      auto nSigmaTPCpos = static_cast<float>((posTrack.tpcSignal() - expBethePos) / expSigmaPos);
      auto nSigmaTPCneg = static_cast<float>((negTrack.tpcSignal() - expBetheNeg) / expSigmaNeg);
      if (abs(nSigmaTPCpos) > heliumNsigmaMax || abs(nSigmaTPCneg) > heliumNsigmaMax)
        continue;

      hyperCandidate hyperCand;
      hyperCand.isMatter = abs(nSigmaTPCpos) < abs(nSigmaTPCneg);

      hyperCand.nSigmaHe3 = hyperCand.isMatter ? nSigmaTPCpos : nSigmaTPCneg;
      hyperCand.nTPCClustersHe3 = hyperCand.isMatter ? posTrack.tpcNClsFindable() : negTrack.tpcNClsFindable();

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

      auto& propPosTrack = fitter.getTrack(0);
      auto& propNegTrack = fitter.getTrack(1);

      std::array<float, 3> posTrackP;
      std::array<float, 3> negTrackP;

      propPosTrack.getPxPyPzGlo(posTrackP);
      propNegTrack.getPxPyPzGlo(negTrackP);

      unsigned int posAbsCharge = hyperCand.isMatter ? 2 : 1;
      unsigned int negAbsCharge = !hyperCand.isMatter ? 2 : 1;

      posTrackP[0] *= posAbsCharge, posTrackP[1] *= posAbsCharge, posTrackP[2] *= posAbsCharge;
      negTrackP[0] *= negAbsCharge, negTrackP[1] *= negAbsCharge, negTrackP[2] *= negAbsCharge;

      // get decay vertex coordinates
      const auto& vtx = fitter.getPCACandidate();
      for (int i = 0; i < 3; i++) {
        hyperCand.decVtx[i] = vtx[i];
        hyperCand.mom[i] = posTrackP[i] + negTrackP[i];
      }

      hyperCand.dcaV0dau = TMath::Sqrt(fitter.getChi2AtPCACandidate());

      // Apply selections so a skimmed table is created only
      if (hyperCand.dcaV0dau > dcav0dau) {
        continue;
      }

      hyperCand.cosPA = RecoDecay::cpa(array{collision.posX(), collision.posY(), collision.posZ()}, array{hyperCand.decVtx[0], hyperCand.decVtx[1], hyperCand.decVtx[2]}, array{hyperCand.mom[0], hyperCand.mom[1], hyperCand.mom[2]});
      if (hyperCand.cosPA < v0cospa) {
        continue;
      }

      // if survived all selections, propagate decay daughters to PV

      gpu::gpustd::array<float, 2> dcaInfo;

      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, posTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
      hyperCand.isMatter ? hyperCand.he3DCAXY = dcaInfo[0] : hyperCand.piDCAXY = dcaInfo[0];

      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, negTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
      hyperCand.isMatter ? hyperCand.piDCAXY = dcaInfo[0] : hyperCand.he3DCAXY = dcaInfo[0];

      // finally, push back the candidate
      hyperCand.isReco = true;
      hyperCandidates.push_back(hyperCand);
    }
  }

  void processRealData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0s const& V0s, TracksFull const& tracks, aod::BCsWithTimestamps const&)
  {
    /* check the previous run number */
    if (!collision.sel8())
      return;

    if (abs(collision.posZ()) > 10.f) {
      return;
    
    fillCandidateData(collision, V0s, tracks);

    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
  }
  PROCESS_SWITCH(hyperRecoTask, processRealData, "Regular analysis", true);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hyperRecoTask>(cfgc)};
}
