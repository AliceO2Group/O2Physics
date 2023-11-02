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

#include <array>
#include <random>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/PIDResponse.h"
#include "DCAFitter/DCAFitterN.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi>;

namespace
{
std::shared_ptr<TH2> hPiRec;
std::shared_ptr<TH3> hPiRecMass;
std::shared_ptr<TH2> hTagCuts;
} // namespace

struct efficiencyQA {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::vertexing::DCAFitterN<2> fitter;
  std::mt19937 gen32;

  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};

  Configurable<bool> debugFlag{"debugFlag", false, "debug flag"};
  Configurable<float> zVtxMax{"zVtxMax", 10.f, "maximum z position of the primary vertex"};

  Configurable<float> yMax{"yMax", .5f, "maximum track rapidity"};
  Configurable<float> etaMax{"etaMax", .8f, "maximum track eta"};
  Configurable<int> tagNtpcClsMin{"tagNtpcClsMin", 80, "minimum number of TPC clusters of the tag track"};
  Configurable<int> tagNitsClsMin{"tagNitsClsMin", 7, "minimum number of ITS clusters of the tag track"};
  Configurable<float> tagNsigmaTpcMax{"tagNsigmaTpcMax", 4.f, "maximum TPC nsigma of the tag track"};
  Configurable<float> tagDcaMin{"tagDcaXYZMin", 0.1f, "tagDCAXYZMin"};

  Configurable<float> massWidth{"massWidth", 0.1f, "invariant mass width"};
  Configurable<float> massMax{"massMax", 0.02f, "maximum invariant mass deviation"};
  Configurable<float> v0radiusMax{"v0radiusMax", 2.2f, "maximum V0 decay radius"};
  Configurable<float> dcaV0dauMax{"dcaV0dauMax", 1.f, "maximum DCA between V0 daughters"};
  Configurable<float> v0cosPaMin{"v0cosPaMin", 0.99f, "minimum cosine of pointing angle of V0"};

  // CCDB options
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> pidPath{"pidPath", "", "Path to the PID response object"};

  ConfigurableAxis massLambdaAxis{"massLambdaAxis", {400, 0.497648f - 0.1f, 0.497648f + 0.1f}, "binning for the Lambda invariant-mass"};
  ConfigurableAxis zVtxAxis{"zVtxAxis", {200, -10.f, 10.f}, "binning for the z coordinate of the primary vertex"};
  ConfigurableAxis recAxis{"recAxis", {5, 0.f, 5.f}, "binning for the probe reconstruction flag"};
  ConfigurableAxis recAxisTag{"recAxisTag", {8, 0.f, 8.f}, "binning for the tag reconstruction flag"};
  ConfigurableAxis ptAxis{"ptAxis", {200, -10.f, 10.f}, "binning for the pt of V0 daughter tracks"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  float d_bz;

  Preslice<aod::V0s> perCollision = o2::aod::v0::collisionId;

  void init(InitContext const&)
  {
    uint32_t randomSeed = static_cast<uint32_t>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    gen32.seed(randomSeed);

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

    histos.add<TH1>("zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});
    histos.add<TH1>("massV0", ";#it{M}(#pi^{+} + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH1F, {massLambdaAxis});
    hPiRec = histos.add<TH2>("piRec", ";;#it{p}_{T} (GeV/#it{c});Entries", HistType::kTH2F, {recAxis, ptAxis});
    hPiRecMass = histos.add<TH3>("piRecMass", ";;#it{p}_{T} (GeV/#it{c});Entries", HistType::kTH3F, {recAxis, ptAxis, massLambdaAxis});
    hTagCuts = histos.add<TH2>("tagCuts", ";;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {recAxisTag, ptAxis});
    std::string binLabels[]{"All", "ITS only", "TPC only", "ITS+TPC", "TPC+TOF"};
    for (int iB{0}; iB < 5; ++iB) {
      hPiRec->GetXaxis()->SetBinLabel(iB + 1, binLabels[iB].data());
      hPiRecMass->GetXaxis()->SetBinLabel(iB + 1, binLabels[iB].data());
    }
    std::string binLabelsTag[]{"hasITS && hasTPC", "tracking", "PID", "v0 mass", "dcaV0dau", "cosPA", "dcaXYZ", "V0radius"};
    for (int iB{0}; iB < 8; ++iB) {
      hTagCuts->GetXaxis()->SetBinLabel(iB + 1, binLabelsTag[iB].data());
    }
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

  template <class T, class Hist>
  void fillHistV0Daugh(T const& track, std::shared_ptr<Hist> hist, float const& pt, float const& mass = 1)
  {
    auto trackPt = track.sign() * pt;
    hist->Fill(0., trackPt, mass);
    bool itsAccept = !(track.itsChi2NCl() > 36.);
    bool tpcAccept = !(track.tpcCrossedRowsOverFindableCls() < 0.8 || track.tpcNClsCrossedRows() < 70 || track.tpcChi2NCl() > 4.);
    if (track.hasITS() && itsAccept && !track.hasTPC()) {
      hist->Fill(1., trackPt, mass);
    }
    if (track.hasTPC() && tpcAccept) {
      if (!track.hasITS()) {
        hist->Fill(2., trackPt, mass);
      }
      if (track.hasITS() && itsAccept) {
        hist->Fill(3., trackPt, mass);
      }
      if (track.hasTOF() && !track.hasITS()) {
        hist->Fill(4., trackPt, mass);
      }
    }
  }

  template <class T>
  void fillCandidateData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0s const& V0s)
  {
    for (auto& v0 : V0s) {
      auto posTrack = v0.posTrack_as<T>();
      auto negTrack = v0.negTrack_as<T>();

      if (std::abs(posTrack.eta()) > etaMax || std::abs(negTrack.eta()) > etaMax) {
        continue;
      }

      bool isPosTag = posTrack.hasTPC() && posTrack.hasITS();
      bool isNegTag = negTrack.hasTPC() && negTrack.hasITS();

      if (!isPosTag && !isNegTag) {
        continue;
      }

      float rnd = static_cast<float>(gen32()) / static_cast<float>(gen32.max());
      bool flagTagProbe = (isPosTag && !isNegTag) || (isPosTag && isNegTag && rnd > 0.5);

      auto& tagTrack = flagTagProbe ? posTrack : negTrack;
      auto& probeTrack = flagTagProbe ? negTrack : posTrack;

      histos.fill(HIST("tagCuts"), 0., tagTrack.sign() * tagTrack.pt());

      // track selections
      bool itsRejectTag = tagTrack.itsNCls() < tagNitsClsMin || tagTrack.itsChi2NCl() > 36.;
      bool tpcRejectTag = tagTrack.tpcNClsFound() < tagNtpcClsMin || tagTrack.tpcCrossedRowsOverFindableCls() < 0.8 || tagTrack.tpcNClsCrossedRows() < 70 || tagTrack.tpcChi2NCl() > 4.;
      if (itsRejectTag || tpcRejectTag) {
        continue;
      }

      histos.fill(HIST("tagCuts"), 1., tagTrack.sign() * tagTrack.pt());

      // pid selections
      if (std::abs(tagTrack.tpcNSigmaPi()) > tagNsigmaTpcMax) {
        continue;
      }

      histos.fill(HIST("tagCuts"), 2., tagTrack.sign() * tagTrack.pt());

      auto tagTrackCov = getTrackParCov(tagTrack);
      auto probeTrackCov = getTrackParCov(probeTrack);

      int nCand = 0;
      try {
        nCand = fitter.process(tagTrackCov, probeTrackCov);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        continue;
      }
      if (nCand == 0) {
        continue;
      }

      std::array<float, 3> momTag;
      std::array<float, 3> momProbe;

      auto& propTrackTag = fitter.getTrack(0);
      auto& propTrackProbe = fitter.getTrack(1);
      propTrackTag.getPxPyPzGlo(momTag);
      propTrackProbe.getPxPyPzGlo(momProbe);

      float piMass = o2::constants::physics::MassPionCharged;
      float tagP2 = momTag[0] * momTag[0] + momTag[1] * momTag[1] + momTag[2] * momTag[2];
      float probeP2 = momProbe[0] * momProbe[0] + momProbe[1] * momProbe[1] + momProbe[2] * momProbe[2];
      float tagE = std::sqrt(tagP2 + piMass * piMass);
      float probeE = std::sqrt(probeP2 + piMass * piMass);
      float v0E = tagE + probeE;

      std::array<float, 3> primVtx = array{collision.posX(), collision.posY(), collision.posZ()};
      std::array<float, 3> v0Mom;
      const auto& vtx = fitter.getPCACandidate();
      for (int i = 0; i < 3; i++) {
        v0Mom[i] = momTag[i] + momProbe[i];
      }

      float massV0 = std::sqrt(v0E * v0E - v0Mom[0] * v0Mom[0] - v0Mom[1] * v0Mom[1] - v0Mom[2] * v0Mom[2]);
      bool isMassV0 = false;
      if (std::abs(massV0 - o2::constants::physics::MassKaonNeutral) < massWidth)
        isMassV0 = true;
      if (!isMassV0) {
        continue;
      }

      histos.fill(HIST("tagCuts"), 3., tagTrack.sign() * tagTrack.pt());

      auto dcaV0dau = std::sqrt(fitter.getChi2AtPCACandidate());
      if (dcaV0dau > dcaV0dauMax) {
        continue;
      }

      histos.fill(HIST("tagCuts"), 4., tagTrack.sign() * tagTrack.pt());

      double cosPA = RecoDecay::cpa(primVtx, vtx, v0Mom);
      if (cosPA < v0cosPaMin) {
        continue;
      }

      histos.fill(HIST("tagCuts"), 5., tagTrack.sign() * tagTrack.pt());

      // if survived all selections, propagate decay daughters to PV
      gpu::gpustd::array<float, 2> dcaInfo;

      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, tagTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
      float tagDcaXYZ = dcaInfo[0];

      if (std::abs(tagDcaXYZ) < tagDcaMin) {
        continue;
      }

      histos.fill(HIST("tagCuts"), 6., tagTrack.sign() * tagTrack.pt());

      float v0radius = std::hypot(vtx[0], vtx[1]);
      if (v0radius > v0radiusMax) {
        continue;
      }

      histos.fill(HIST("tagCuts"), 7., tagTrack.sign() * tagTrack.pt());

      histos.fill(HIST("massV0"), massV0);

      float propPt = std::hypot(momProbe[0], momProbe[1]);
      fillHistV0Daugh(probeTrack, hPiRecMass, propPt, massV0);
      if (std::abs(massV0 - o2::constants::physics::MassKaonNeutral) < massMax) {
        fillHistV0Daugh(probeTrack, hPiRec, propPt);
      }
    }
  }

  void processData(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::V0s const& V0s, TracksFull const& tracks, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.sel8())
        continue;

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      histos.fill(HIST("zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollision, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      fillCandidateData<TracksFull>(collision, V0Table_thisCollision);
    }
  }
  PROCESS_SWITCH(efficiencyQA, processData, "Data analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<efficiencyQA>(cfgc)};
}
