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
#include <vector>
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
const float piMass = o2::constants::physics::MassPionCharged;
std::shared_ptr<TH2> hPiRec;
std::shared_ptr<TH3> hPiRecMass;
std::shared_ptr<TH2> hTagCuts;
std::shared_ptr<TH2> hTpcSegment;
std::shared_ptr<TH3> hPtRes;
std::shared_ptr<TH3> hEtaRes;
std::shared_ptr<TH3> hPhiRes;
float invMass2Body(std::array<float, 3>& momA, std::array<float, 3> const& momB, std::array<float, 3> const& momC, float const& massB, float const& massC)
{
  float p2B = momB[0] * momB[0] + momB[1] * momB[1] + momB[2] * momB[2];
  float p2C = momC[0] * momC[0] + momC[1] * momC[1] + momC[2] * momC[2];
  for (int i = 0; i < 3; ++i) {
    momA[i] = momB[i] + momC[i];
  }
  float eB = std::sqrt(p2B + massB * massB);
  float eC = std::sqrt(p2C + massC * massC);
  float eA = eB + eC;
  float massA = std::sqrt(eA * eA - momA[0] * momA[0] - momA[1] * momA[1] - momA[2] * momA[2]);
  return massA;
}
} // namespace

struct ProbeTrack {
  ProbeTrack(uint64_t const& idx, uint64_t const& idxTpc, uint64_t const& idxTag, float const& p, float const& pt, float const& eta, float const& phi, float const& mass, uint8_t const& map, float const& xv, float const& yv, float const& zv) : globalIndex{idx},
                                                                                                                                                                                                                                                   globalIndexTpc{idxTpc},
                                                                                                                                                                                                                                                   p{p},
                                                                                                                                                                                                                                                   pt{pt},
                                                                                                                                                                                                                                                   eta{eta},
                                                                                                                                                                                                                                                   phi{phi},
                                                                                                                                                                                                                                                   massTagProbe{mass},
                                                                                                                                                                                                                                                   vtx0{xv},
                                                                                                                                                                                                                                                   vtx1{yv},
                                                                                                                                                                                                                                                   vtx2{zv},
                                                                                                                                                                                                                                                   detectorMap{map},
                                                                                                                                                                                                                                                   globalIndexTag{idxTag}
  {
  }
  uint64_t globalIndex;
  uint64_t globalIndexTpc;
  float p;
  float pt;
  float eta;
  float phi;
  float massTagProbe;
  float vtx0;
  float vtx1;
  float vtx2;
  uint8_t detectorMap;
  uint64_t globalIndexTag;
};

struct efficiencyQA {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::vertexing::DCAFitterN<2> fitter;
  std::mt19937 gen32;

  std::vector<ProbeTrack> probeTracks;

  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};

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

  Configurable<bool> findTpcLeg{"findTpcLeg", false, "toggle search of missing tpc segment"};
  Configurable<bool> refitVertex{"refitVertex", false, "toggle refit of decay vertex using tag and tpc prolongation"};
  Configurable<float> ptWindow{"ptWindow", 0.05f, "pt window to search tpc segment"};
  Configurable<float> etaWindow{"etaWindow", 0.3f, "eta window to search tpc segment"};
  Configurable<float> phiWindow{"phiWindow", 0.2f, "phi window to search tpc segment"};
  Configurable<float> massWindow{"massWindow", 0.03f, "mass window to search tpc segment"};
  Configurable<float> cosPaWindow{"cosPaWindow", 0.8f, "cosPa window to search tpc segment"};

  // CCDB options
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> pidPath{"pidPath", "", "Path to the PID response object"};

  ConfigurableAxis massK0sAxis{"massK0sAxis", {400, 0.497648f - 0.1f, 0.497648f + 0.1f}, "binning for the K0s invariant-mass"};
  ConfigurableAxis massResAxis{"massResAxis", {400, -0.1f, 0.1f}, "binning for the invariant-mass resolution"};
  ConfigurableAxis zVtxAxis{"zVtxAxis", {200, -10.f, 10.f}, "binning for the z coordinate of the primary vertex"};
  ConfigurableAxis recAxis{"recAxis", {8, 0.f, 8.f}, "binning for the reconstruction flag"};
  ConfigurableAxis tpcAxis{"tpcAxis", {3, 0.f, 3.f}, "binning for the matching of tpc segment"};
  ConfigurableAxis ptAxis{"ptAxis", {100, -5.f, 5.f}, "binning for the pt of V0 daughter tracks"};
  ConfigurableAxis ptResAxis{"ptResAxis", {200, -1.f, 1.f}, "binning for the pt resolution of V0 daughter tracks"};
  ConfigurableAxis etaAxis{"etaAxis", {900, -.9f, .9f}, "binning for the eta of V0 daughter tracks"};
  ConfigurableAxis etaResAxis{"etaResAxis", {500, -.5f, .5f}, "binning for the eta resolution of V0 daughter tracks"};
  ConfigurableAxis phiAxis{"phiAxis", {630, 0.f, 6.3f}, "binning for the phi of V0 daughter tracks"};
  ConfigurableAxis phiResAxis{"phiResAxis", {800, -4.f, 4.f}, "binning for the phi resolution of V0 daughter tracks"};
  ConfigurableAxis cosPaAxis{"cosPaAxis", {1000, -1.f, 1.f}, "binning for the cosine of pointing angle"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  float d_bz;

  Preslice<aod::V0s> perCollisionV0s = o2::aod::v0::collisionId;
  Preslice<TracksFull> perCollisionTracks = o2::aod::track::collisionId;

  void init(InitContext const&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    histos.add<TH1>("zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});
    hPiRec = histos.add<TH2>("piRec", ";;#it{p}_{T} (GeV/#it{c});Entries", HistType::kTH2F, {recAxis, ptAxis});
    std::string binLabels[]{"Decays", "ITS", "ITS only", "TPC", "TPC only", "ITS+TPC", "TPC+TOF", " "};
    for (int iB{0}; iB < 8; ++iB) {
      hPiRec->GetXaxis()->SetBinLabel(iB + 1, binLabels[iB].data());
    }
    if (doprocessMcTracks) {
      hPiRec->GetXaxis()->SetBinLabel(1, "Generated");
      hPtRes = histos.add<TH3>("ptRes", ";;#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec} - #it{p}_{T}^{MC} (GeV/#it{c})", HistType::kTH3F, {recAxis, ptAxis, ptResAxis});
      hEtaRes = histos.add<TH3>("etaRes", ";;#eta^{rec};#eta^{rec} - #eta^{MC} (rad)", HistType::kTH3F, {recAxis, etaAxis, etaResAxis});
      hPhiRes = histos.add<TH3>("phiRes", ";;#phi^{rec} (rad);#phi^{rec} - #phi^{MC} (rad)", HistType::kTH3F, {recAxis, phiAxis, phiResAxis});
      for (int iB{1}; iB < 8; ++iB) {
        hPtRes->GetXaxis()->SetBinLabel(iB + 1, binLabels[iB].data());
        hEtaRes->GetXaxis()->SetBinLabel(iB + 1, binLabels[iB].data());
        hPhiRes->GetXaxis()->SetBinLabel(iB + 1, binLabels[iB].data());
      }
      hPtRes->GetXaxis()->SetBinLabel(1, " ");
      hEtaRes->GetXaxis()->SetBinLabel(1, " ");
      hPhiRes->GetXaxis()->SetBinLabel(1, " ");
    } else if (doprocessTagAndProbe || doprocessTagAndProbeMC) {
      hPiRec->GetXaxis()->SetBinLabel(8, "ITS w/ TPC leg");
      uint32_t randomSeed = static_cast<uint32_t>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      gen32.seed(randomSeed);

      mRunNumber = 0;
      d_bz = 0;

      fitter.setPropagateToPCA(true);
      fitter.setMaxR(200.);
      fitter.setMinParamChange(1e-3);
      fitter.setMinRelChi2Change(0.9);
      fitter.setMaxDZIni(1e9);
      fitter.setMaxChi2(1e9);
      fitter.setUseAbsDCA(true);
      int mat{static_cast<int>(cfgMaterialCorrection)};
      fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));

      histos.add<TH1>("massV0", ";#it{M}(#pi^{+} + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH1F, {massK0sAxis});
      hPiRecMass = histos.add<TH3>("piRecMass", ";;#it{p}_{T} (GeV/#it{c});#it{M}(#pi^{+} + #pi^{-}) (GeV/#it{c}^{2})", HistType::kTH3F, {recAxis, ptAxis, massK0sAxis});
      hTagCuts = histos.add<TH2>("tagCuts", ";;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {recAxis, ptAxis});

      std::string binLabelsTag[]{"hasITS && hasTPC", "tracking", "PID", "v0 mass", "dcaV0dau", "cosPA", "dcaXYZ", "V0radius"};
      for (int iB{0}; iB < 8; ++iB) {
        hTagCuts->GetXaxis()->SetBinLabel(iB + 1, binLabelsTag[iB].data());
      }
      for (int iB{0}; iB < 7; ++iB) {
        hPiRecMass->GetXaxis()->SetBinLabel(iB + 1, binLabels[iB].data());
      }
      hPiRecMass->GetXaxis()->SetBinLabel(8, "ITS w/ TPC leg");

      if (doprocessTagAndProbeMC) {
        std::string binLabelsTpc[]{"hasTPCsegment", "foundTPCsegment", "allFoundTPCsegment"};
        hTpcSegment = histos.add<TH2>("tpcSegment", ";;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {tpcAxis, ptAxis});
        for (int iB{0}; iB < 3; ++iB) {
          hTpcSegment->GetXaxis()->SetBinLabel(iB + 1, binLabelsTpc[iB].data());
        }
        histos.add<TH2>("pTpcIts", ";#it{p}^{ITS} (GeV/#it{c});#it{p}^{TPC} - #it{p}^{ITS} (GeV/#it{c});Entries", HistType::kTH2F, {ptAxis, ptResAxis});
        histos.add<TH2>("ptTpcIts", ";#it{p}_{T}^{ITS} (GeV/#it{c});#it{p}^{TPC}_{T} - #it{p}^{ITS}_{T} (GeV/#it{c});Entries", HistType::kTH2F, {ptAxis, ptResAxis});
        histos.add<TH2>("etaTpcIts", ";#it{p}^{ITS} (GeV/#it{c});#eta^{TPC} - #eta^{ITS};Entries", HistType::kTH2F, {ptAxis, etaResAxis});
        histos.add<TH2>("phiTpcIts", ";#it{p}^{ITS} (GeV/#it{c});#phi^{TPC} - #phi^{ITS} (rad);Entries", HistType::kTH2F, {ptAxis, phiResAxis});
        histos.add<TH2>("massTpc", ";#it{p}^{ITS} (GeV/#it{c});#it{M}^{TPC} (GeV/#it{c}^{2});Entries", HistType::kTH2F, {ptAxis, massK0sAxis});
        histos.add<TH2>("massTpcIts", ";#it{p}^{ITS} (GeV/#it{c});#it{M}^{TPC} - #it{M}^{ITS} (GeV/#it{c}^{2});Entries", HistType::kTH2F, {ptAxis, massResAxis});
        histos.add<TH2>("cosPaTpc", ";#it{p}^{ITS} (GeV/#it{c});cos#theta_{p}", HistType::kTH2F, {ptAxis, cosPaAxis});
        histos.add<TH3>("ptEtaPhiTpcIts", ";#it{p}^{TPC}_{T} - #it{p}^{ITS}_{T} (GeV/#it{c});#eta^{TPC} - #eta^{ITS};#phi^{TPC} - #phi^{ITS} (rad)", HistType::kTH3F, {ptResAxis, etaResAxis, phiResAxis});
      }
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
  void fillHistTrack(T const& track, std::shared_ptr<Hist> hist, float const& y, float const& z = 1)
  {
    bool itsAccept = !(track.itsChi2NCl() > 36. || track.itsNCls() < 4);
    bool tpcAccept = !(track.tpcCrossedRowsOverFindableCls() < 0.8 || track.tpcNClsCrossedRows() < 70 || track.tpcChi2NCl() > 4. || track.tpcNClsFound() < 90);
    if (track.hasITS()) {
      hist->Fill(1., y, z);
    }
    if (track.hasITS() && itsAccept && !track.hasTPC()) {
      hist->Fill(2., y, z);
    }
    if (track.hasTPC() && tpcAccept) {
      hist->Fill(3., y, z);
      if (!track.hasITS()) {
        hist->Fill(4., y, z);
      }
      if (track.hasITS() && itsAccept) {
        hist->Fill(5., y, z);
      }
      if (track.hasTOF() && !track.hasITS()) {
        hist->Fill(6., y, z);
      }
    }
  }

  template <class T>
  void fillTagAndProbe(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0s const& V0s, TracksFull const& tracks)
  {
    auto tpcTracks = tracks.sliceBy(perCollisionTracks, collision.globalIndex());
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

      auto propTrackTag = fitter.getTrack(0);
      auto propTrackProbe = fitter.getTrack(1);
      propTrackTag.getPxPyPzGlo(momTag);
      propTrackProbe.getPxPyPzGlo(momProbe);

      std::array<float, 3> v0Mom;
      float massV0 = invMass2Body(v0Mom, momTag, momProbe, piMass, piMass);
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

      std::array<float, 3> primVtx = array{collision.posX(), collision.posY(), collision.posZ()};
      auto vtx = fitter.getPCACandidate();
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

      uint64_t idxTpc{0};
      uint8_t map{0};
      if (probeTrack.hasITS() && !probeTrack.hasTPC() && findTpcLeg) {
        auto acceptIts = !(probeTrack.itsChi2NCl() > 36. || probeTrack.itsNCls() < 4);
        if (!acceptIts) {
          continue;
        }
        map = probeTrack.detectorMap();

        std::array<float, 3> momTpc;
        for (auto& tpcTrack : tpcTracks) {
          if (std::abs(tpcTrack.eta()) > etaMax) {
            continue;
          }

          if (!tpcTrack.hasTPC() || tpcTrack.hasITS()) {
            continue;
          }

          bool acceptTpc = !(tpcTrack.tpcCrossedRowsOverFindableCls() < 0.8 || tpcTrack.tpcNClsCrossedRows() < 70 || tpcTrack.tpcChi2NCl() > 4. || tpcTrack.tpcNClsFound() < 90);
          bool acceptCharge = tpcTrack.sign() == probeTrack.sign();
          if (!acceptTpc || !acceptCharge) {
            continue;
          }

          gpu::gpustd::array<float, 2> dcaInfo;
          auto tpcTrackCov = getTrackParCov(tpcTrack);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({static_cast<float>(vtx[0]), static_cast<float>(vtx[1]), static_cast<float>(vtx[2])}, tpcTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);

          bool acceptTrackPt = std::abs(tpcTrackCov.getPt() - propTrackProbe.getPt()) < ptWindow;
          bool acceptTrackEta = std::abs(tpcTrackCov.getEta() - propTrackProbe.getEta()) < etaWindow;
          bool acceptTrackPhi = std::abs(tpcTrackCov.getPhi() - propTrackProbe.getPhi()) < phiWindow;
          bool acceptTpcTrack = acceptTrackPt && acceptTrackEta && acceptTrackPhi;

          // LOGF(debug, "idx = %lld, Dpt = %f, Deta = %f, Dphi = %f", tpcTrack.globalIndex(), std::abs(tpcTrackCov.getPt() - propTrackProbe.getPt()), std::abs(tpcTrackCov.getEta() - propTrackProbe.getEta()), std::abs(tpcTrackCov.getPhi() - propTrackProbe.getPhi()));

          if (acceptTpcTrack) {
            tpcTrackCov.getPxPyPzGlo(momTpc);
            auto massTpcLeg = invMass2Body(v0Mom, momTag, momTpc, piMass, piMass);
            auto massDiff = std::abs(massTpcLeg - o2::constants::physics::MassKaonNeutral);

            LOGF(debug, "idx = %lld, p = %f, Dpt = %f, Deta = %f, Dphi = %f, mass = %f", tpcTrack.globalIndex(), propTrackProbe.getP(), std::abs(tpcTrackCov.getPt() - propTrackProbe.getPt()), std::abs(tpcTrackCov.getEta() - propTrackProbe.getEta()), std::abs(tpcTrackCov.getPhi() - propTrackProbe.getPhi()), massDiff);
            if (massDiff < massWindow) {
              if (refitVertex) {
                tpcTrackCov = getTrackParCov(tpcTrack);
                tagTrackCov = getTrackParCov(tagTrack);
                int nCand = 0;
                try {
                  nCand = fitter.process(tagTrackCov, tpcTrackCov);
                } catch (...) {
                  LOG(error) << "Exception caught in DCA fitter process call!";
                  continue;
                }
                if (nCand == 0) {
                  continue;
                }

                std::array<float, 3> momTpcRefit, momTagRefit, v0MomRefit;
                auto& propTrackTagRefit = fitter.getTrack(0);
                auto& propTrackTpcRefit = fitter.getTrack(1);
                propTrackTagRefit.getPxPyPzGlo(momTagRefit);
                propTrackTpcRefit.getPxPyPzGlo(momTpcRefit);
                auto massTpcLeg2 = invMass2Body(v0MomRefit, momTagRefit, momTpcRefit, piMass, piMass);
                auto massDiff2 = std::abs(massTpcLeg2 - o2::constants::physics::MassKaonNeutral);

                const auto& vtxRefit = fitter.getPCACandidate();
                double cosPa = RecoDecay::cpa(primVtx, vtxRefit, v0MomRefit);

                LOGF(debug, "idx = %lld, p = %f, Dpt = %f, Deta = %f, Dphi = %f, mass = %f, mass2 = %f, cosPA = %f", tpcTrack.globalIndex(), propTrackProbe.getP(), std::abs(propTrackTpcRefit.getPt() - propTrackProbe.getPt()), std::abs(propTrackTpcRefit.getEta() - propTrackProbe.getEta()), std::abs(propTrackTpcRefit.getPhi() - propTrackProbe.getPhi()), massDiff, massDiff2, cosPa);

                if (cosPa < cosPaWindow) {
                  continue;
                }
              }
              idxTpc = tpcTrack.globalIndex();
            }
          }
        }
      }

      auto trackPt = probeTrack.sign() * std::hypot(momProbe[0], momProbe[1]);
      auto trackP = probeTrack.sign() * std::hypot(trackPt, momProbe[2]);

      hPiRecMass->Fill(0., trackPt, massV0);
      if (idxTpc > 0) {
        hPiRecMass->Fill(7., trackPt, massV0);
      }
      fillHistTrack(probeTrack, hPiRecMass, trackPt, massV0);
      if (std::abs(massV0 - o2::constants::physics::MassKaonNeutral) < massMax) {
        probeTracks.emplace_back(probeTrack.globalIndex(), idxTpc, tagTrack.globalIndex(), trackP, trackPt, propTrackProbe.getEta(), propTrackProbe.getPhi(), massV0, map, vtx[0], vtx[1], vtx[2]);
        hPiRec->Fill(0., trackPt);
        fillHistTrack(probeTrack, hPiRec, trackPt);
        if (idxTpc > 0) {
          hPiRec->Fill(7., trackPt);
        }
      }
    }
  }

  void fillProbeMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, TracksFull const& tracks, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC)
  {
    auto tracks_thisEvent = tracks.sliceBy(perCollisionTracks, collision.globalIndex());

    for (const auto& probeTrack : probeTracks) {
      auto mcLab = trackLabelsMC.rawIteratorAt(probeTrack.globalIndex);
      if (mcLab.mcParticleId() < -1 || mcLab.mcParticleId() >= particlesMC.size()) {
        continue;
      }
      bool hasITS = (probeTrack.detectorMap & o2::aod::track::ITS) == o2::aod::track::ITS;
      bool hasTPC = (probeTrack.detectorMap & o2::aod::track::TPC) == o2::aod::track::TPC;
      if (hasITS && !hasTPC) {
        auto tagTrack = tracks.rawIteratorAt(probeTrack.globalIndexTag);
        auto tagTrackCov = getTrackParCov(tagTrack);
        gpu::gpustd::array<float, 2> dcaInfo;
        o2::base::Propagator::Instance()->propagateToDCABxByBz({probeTrack.vtx0, probeTrack.vtx1, probeTrack.vtx2}, tagTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
        std::array<float, 3> momTag;
        tagTrackCov.getPxPyPzGlo(momTag);

        for (const auto& tpcTrack : tracks_thisEvent) {
          if (std::abs(tpcTrack.eta()) > etaMax) {
            continue;
          }

          if (!tpcTrack.hasTPC() || tpcTrack.hasITS() || tpcTrack.sign() != probeTrack.pt / std::abs(probeTrack.pt)) {
            continue;
          }
          if (tpcTrack.tpcCrossedRowsOverFindableCls() < 0.8 || tpcTrack.tpcNClsCrossedRows() < 70 || tpcTrack.tpcChi2NCl() > 4. || tpcTrack.tpcNClsFound() < 90) {
            continue;
          }
          auto mcLabTpc = trackLabelsMC.rawIteratorAt(tpcTrack.globalIndex());
          if (mcLabTpc.mcParticleId() < -1 || mcLabTpc.mcParticleId() >= particlesMC.size()) {
            continue;
          }
          if (mcLabTpc.mcParticleId() == mcLab.mcParticleId()) {
            hTpcSegment->Fill(0., probeTrack.pt);

            auto trackCov = getTrackParCov(tpcTrack);
            gpu::gpustd::array<float, 2> dcaInfo;
            o2::base::Propagator::Instance()->propagateToDCABxByBz({probeTrack.vtx0, probeTrack.vtx1, probeTrack.vtx2}, trackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
            std::array<float, 3> v0Mom, momTpc;
            trackCov.getPxPyPzGlo(momTpc);
            auto massTpcLeg = invMass2Body(v0Mom, momTag, momTpc, piMass, piMass);

            histos.fill(HIST("pTpcIts"), probeTrack.p, trackCov.getP() - probeTrack.p * probeTrack.p / std::abs(probeTrack.p));
            histos.fill(HIST("ptTpcIts"), probeTrack.p, trackCov.getPt() - probeTrack.pt * probeTrack.pt / std::abs(probeTrack.pt));
            histos.fill(HIST("etaTpcIts"), probeTrack.p, trackCov.getEta() - probeTrack.eta);
            histos.fill(HIST("phiTpcIts"), probeTrack.p, trackCov.getPhi() - probeTrack.phi);
            histos.fill(HIST("massTpc"), probeTrack.p, massTpcLeg);
            histos.fill(HIST("massTpcIts"), probeTrack.p, massTpcLeg - probeTrack.massTagProbe);
            histos.fill(HIST("ptEtaPhiTpcIts"), trackCov.getPt() - probeTrack.pt * probeTrack.pt / std::abs(probeTrack.pt), trackCov.getEta() - probeTrack.eta, trackCov.getPhi() - probeTrack.phi);

            LOGF(debug, "MC::DEBUG: idx = %lld, p = %f, Dpt = %f, Deta = %f, Dphi = %f", tpcTrack.globalIndex(), probeTrack.pt, std::abs(trackCov.getPt() - probeTrack.pt * probeTrack.pt / std::abs(probeTrack.pt)), std::abs(trackCov.getEta() - probeTrack.eta), std::abs(trackCov.getPhi() - probeTrack.phi));

            int nCand = 0;
            try {
              nCand = fitter.process(tagTrackCov, trackCov);
            } catch (...) {
              LOG(error) << "Exception caught in DCA fitter process call!";
              continue;
            }
            if (nCand == 0) {
              continue;
            }

            auto& propTrackTagRefit = fitter.getTrack(0);
            auto& propTrackTpcRefit = fitter.getTrack(1);
            propTrackTagRefit.getPxPyPzGlo(momTag);
            propTrackTpcRefit.getPxPyPzGlo(momTpc);
            invMass2Body(v0Mom, momTag, momTpc, piMass, piMass);

            const auto& vtxRefit = fitter.getPCACandidate();
            double cosPa = RecoDecay::cpa(std::array<float, 3>{collision.posX(), collision.posY(), collision.posZ()}, vtxRefit, v0Mom);

            histos.fill(HIST("cosPaTpc"), probeTrack.p, cosPa);

            break;
          }
        }
        if (probeTrack.globalIndexTpc > 0) {
          hTpcSegment->Fill(2., probeTrack.pt);

          auto mcLabTpc = trackLabelsMC.rawIteratorAt(probeTrack.globalIndexTpc);
          if (mcLabTpc.mcParticleId() < -1 || mcLabTpc.mcParticleId() >= particlesMC.size()) {
            continue;
          }
          if (mcLabTpc.mcParticleId() == mcLab.mcParticleId()) {
            hTpcSegment->Fill(1., probeTrack.pt);
          }
        }
      }
    }
  }

  void processTagAndProbe(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::V0s const& V0s, TracksFull const& tracks, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      probeTracks.clear();

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.sel8())
        continue;

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      histos.fill(HIST("zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0s, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);
      fillTagAndProbe<TracksFull>(collision, V0Table_thisCollision, tracks);
    }
  }
  PROCESS_SWITCH(efficiencyQA, processTagAndProbe, "Tag and probe analysis", true);

  void processTagAndProbeMC(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::V0s const& V0s, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC)
  {
    for (const auto& collision : collisions) {
      probeTracks.clear();

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.sel8())
        continue;

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      histos.fill(HIST("zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0s, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);
      fillTagAndProbe<TracksFull>(collision, V0Table_thisCollision, tracks);

      fillProbeMC(collision, tracks, trackLabelsMC, particlesMC);
    }
  }
  PROCESS_SWITCH(efficiencyQA, processTagAndProbeMC, "Tag and probe analysis on MC", false);

  void processMcTracks(soa::Join<aod::Collisions, aod::EvSels> const& collisions, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC)
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
      auto TrackTable_thisCollision = tracks.sliceBy(perCollisionTracks, collIdx);

      for (auto& track : TrackTable_thisCollision) {
        auto mcLab = trackLabelsMC.rawIteratorAt(track.globalIndex());
        if (mcLab.mcParticleId() < -1 || mcLab.mcParticleId() >= particlesMC.size()) {
          continue;
        }
        if (mcLab.has_mcParticle()) {
          auto mcTrack = mcLab.template mcParticle_as<aod::McParticles>();
          if (std::abs(track.eta()) > etaMax || std::abs(track.rapidity(o2::constants::physics::MassPionCharged)) > yMax) {
            continue;
          }
          if (std::abs(mcTrack.y()) > yMax) {
            continue;
          }
          if (std::abs(mcTrack.pdgCode()) != 211) {
            continue;
          }
          if (!mcTrack.isPhysicalPrimary()) {
            continue;
          }
          const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};

          auto trackParCov = getTrackParCov(track);
          gpu::gpustd::array<float, 2> dcaInfo;
          o2::base::Propagator::Instance()->propagateToDCA(collVtx, trackParCov, d_bz, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfo);

          auto trackPt = track.sign() * trackParCov.getPt();
          fillHistTrack(track, hPiRec, trackPt);

          fillHistTrack(track, hPtRes, track.sign() * trackParCov.getPt(), trackParCov.getPt() - mcTrack.pt());
          fillHistTrack(track, hEtaRes, track.sign() * trackParCov.getPt(), trackParCov.getEta() - mcTrack.eta());
          fillHistTrack(track, hPhiRes, track.sign() * trackParCov.getPt(), trackParCov.getPhi() - mcTrack.phi());
        }
      }
    }

    for (auto& partMC : particlesMC) {
      auto pdgCode = partMC.pdgCode();
      if (std::abs(partMC.y()) > yMax) {
        continue;
      }
      if (std::abs(pdgCode) != 211) {
        continue;
      }
      if (!partMC.isPhysicalPrimary()) {
        continue;
      }
      hPiRec->Fill(0., pdgCode / std::abs(pdgCode) * partMC.pt());
    }
  }
  PROCESS_SWITCH(efficiencyQA, processMcTracks, "MC tracks analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<efficiencyQA>(cfgc)};
}
