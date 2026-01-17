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

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include "CCDB/BasicCCDBManager.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <array>
#include <random>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi>;
using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;

namespace
{
const float piMass = o2::constants::physics::MassPionCharged;
std::shared_ptr<TH3> hPiRec;
std::shared_ptr<THnSparse> hPiRecMass;
std::shared_ptr<TH2> hTagCuts;
std::shared_ptr<TH2> hTpcSegment;
std::shared_ptr<TH3> hPtRes;
std::shared_ptr<TH3> hEtaRes;
std::shared_ptr<TH3> hPhiRes;
std::shared_ptr<TH3> hPiRecSec;
constexpr double deltaPtPar[1][3]{{1.30187e-02, 2.06381e-02, 4.46466e-04}};
constexpr double deltaEtaPar[1][2]{{5.06792e-03, -8.93511e-01}};
constexpr double deltaPhiPar[1][3]{{9.39809e-03, 1.08688e-01, 4.75977e+00}};
static const std::vector<std::string> deltaPt{"deltaPt"};
static const std::vector<std::string> deltaEta{"deltaEta"};
static const std::vector<std::string> deltaPhi{"deltaPhi"};
static const std::vector<std::string> parNames2{"p0", "p1"};
static const std::vector<std::string> parNames3{"p0", "p1", "p2"};
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
  uint64_t globalIndex;
  uint64_t globalIndexTpc;
  int32_t collIndex;
  int64_t bcIndex;
  float p;
  float pt;
  float pProp;
  float ptProp;
  float etaProp;
  float phiProp;
  float massTagProbe;
  float vtx0;
  float vtx1;
  float vtx2;
  float time;
  float timeRes;
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
  Configurable<bool> useTpcTracksFromSameColl{"useTpcTracksFromSameColl", true, "toggle post-matching to tpc segment associated to collision"};
  Configurable<bool> useCollisionWindow{"useCollisionWindow", false, "toogle collision window in re-matching"};
  Configurable<bool> propToTPCinnerWall{"propToTPCinnerWall", false, "toggle propagation of tracks to the TPC inner wall"};
  Configurable<bool> refitVertex{"refitVertex", false, "toggle refit of decay vertex using tag and tpc prolongation"};
  Configurable<float> ptWindow{"ptWindow", 0.05f, "pt window to search tpc segment"};
  Configurable<float> etaWindow{"etaWindow", 0.3f, "eta window to search tpc segment"};
  Configurable<float> phiWindow{"phiWindow", 0.2f, "phi window to search tpc segment"};
  Configurable<float> massWindow{"massWindow", 0.03f, "mass window to search tpc segment"};
  Configurable<float> cosPaWindow{"cosPaWindow", 0.8f, "cosPa window to search tpc segment"};
  Configurable<int> collIdWindow{"collIdWindow", 6, "collision index window to search tpc segment"};

  Configurable<float> trackTimingCut{"trackTimingCut", 3.f, "track timing cut, number of sigmas"};

  Configurable<float> nSigmaDeltaPt{"nSigmaDeltaPt", 3.f, "pt window number of sigmas"};
  Configurable<float> nSigmaDeltaEta{"nSigmaDeltaEta", 3.f, "eta window number of sigmas"};
  Configurable<float> nSigmaDeltaPhi{"nSigmaDeltaPhi", 3.f, "phi window number of sigmas"};

  Configurable<LabeledArray<double>> cfgDeltaPt{"cfgDeltaPt", {deltaPtPar[0], 1, 3, deltaPt, parNames3}, "parameterisation for delta pt standard deviation with respect to momentum"};
  Configurable<LabeledArray<double>> cfgDeltaEta{"cfgDeltaEta", {deltaEtaPar[0], 1, 2, deltaEta, parNames2}, "parameterisation for delta eta standard deviation with respect to momentum"};
  Configurable<LabeledArray<double>> cfgDeltaPhi{"cfgDeltaPhi", {deltaPhiPar[0], 1, 3, deltaPhi, parNames3}, "parameterisation for delta phi standard deviation with respect to momentum"};

  Configurable<bool> setMomFlatWindows{"setMomFlatWindows", true, "toggle momentum-flat windows"};

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
  ConfigurableAxis tpcAxis{"tpcAxis", {4, 0.f, 4.f}, "binning for the matching of tpc segment"};
  ConfigurableAxis ptAxis{"ptAxis", {100, -5.f, 5.f}, "binning for the pt of V0 daughter tracks"};
  ConfigurableAxis ptResAxis{"ptResAxis", {200, -1.f, 1.f}, "binning for the pt resolution of V0 daughter tracks"};
  ConfigurableAxis etaAxis{"etaAxis", {900, -.9f, .9f}, "binning for the eta of V0 daughter tracks"};
  ConfigurableAxis etaResAxis{"etaResAxis", {500, -.5f, .5f}, "binning for the eta resolution of V0 daughter tracks"};
  ConfigurableAxis phiAxis{"phiAxis", {630, 0.f, 6.3f}, "binning for the phi of V0 daughter tracks"};
  ConfigurableAxis phiResAxis{"phiResAxis", {800, -4.f, 4.f}, "binning for the phi resolution of V0 daughter tracks"};
  ConfigurableAxis cosPaAxis{"cosPaAxis", {1000, -1.f, 1.f}, "binning for the cosine of pointing angle"};
  ConfigurableAxis collIdResAxis{"collIdResAxis", {1.e2, -50., 50.}, "binning for the collision ID resolution"};
  ConfigurableAxis timeResAxis{"timeResAxis", {1000, -50., 50.}, "binning for the time difference (normalized by time resolution)"};
  ConfigurableAxis timeResAxisNoNorm{"timeResAxisNoNorm", {500, -20000., 20000.}, "binning for the time difference (not normalized) (ns)"};

  ConfigurableAxis nGenRecAxis{"nGenRecAxis", {20, 0, 20}, "binning for the detector response matrix axis"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  float d_bz;

  Preslice<aod::V0s> perCollisionV0s = o2::aod::v0::collisionId;
  Preslice<TracksFull> perCollisionTracks = o2::aod::track::collisionId;

  double LHCRFFreq = 400.789e6;
  double LHCBunchSpacingNS = 10 * 1.e9 / LHCRFFreq;

  void init(InitContext const&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    histos.add<TH1>("zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});
    hPiRec = histos.add<TH3>("piRec", ";;#it{p}_{T} (GeV/#it{c});#it{#eta}", HistType::kTH3F, {recAxis, ptAxis, etaAxis});
    std::string binLabels[]{"Decays", "ITS", "ITS only", "TPC", "TPC only", "ITS+TPC", "TPC+TOF", " "};
    for (int iB{0}; iB < hPiRec->GetNbinsX(); ++iB) {
      hPiRec->GetXaxis()->SetBinLabel(iB + 1, binLabels[iB].data());
    }
    if (doprocessMcTracks) {
      hPiRecSec = histos.add<TH3>("piRecSec", ";;#it{p}_{T} (GeV/#it{c});#it{#eta}", HistType::kTH3F, {recAxis, ptAxis, etaAxis});
      hPiRec->GetXaxis()->SetBinLabel(1, "Generated");
      hPiRecSec->GetXaxis()->SetBinLabel(1, "Generated");
      hPtRes = histos.add<TH3>("ptRes", ";;#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec} - #it{p}_{T}^{MC} (GeV/#it{c})", HistType::kTH3F, {recAxis, ptAxis, ptResAxis});
      hEtaRes = histos.add<TH3>("etaRes", ";;#it{p}_{T}^{rec} (GeV/#it{c});#eta^{rec} - #eta^{MC} (rad)", HistType::kTH3F, {recAxis, ptAxis, etaResAxis});
      hPhiRes = histos.add<TH3>("phiRes", ";;#it{p}_{T}^{rec} (GeV/#it{c});#phi^{rec} - #phi^{MC} (rad)", HistType::kTH3F, {recAxis, ptAxis, phiResAxis});
      for (int iB{1}; iB < hPtRes->GetNbinsX(); ++iB) {
        hPiRecSec->GetXaxis()->SetBinLabel(iB + 1, binLabels[iB].data());
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
      hPiRecMass = histos.add<THnSparse>("piRecMass", ";;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{M}(#pi^{+} + #pi^{-}) (GeV/#it{c}^{2})", HistType::kTHnSparseF, {recAxis, ptAxis, etaAxis, massK0sAxis});
      hTagCuts = histos.add<TH2>("tagCuts", ";;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {recAxis, ptAxis});

      histos.add<TH2>("massTagTpc", ";#it{p} (GeV/#it{c});#it{M}(#pi^{+} + #pi^{-}) (GeV/#it{c}^{2})", HistType::kTH2F, {ptAxis, massK0sAxis});
      histos.add<TH2>("ptItsTpcSel", ";#it{p} (GeV/#it{c});#it{p}_{T}^{TPC} - #it{p}_{T}^{ITS} (GeV/#it{c})", HistType::kTH2F, {ptAxis, ptResAxis});
      histos.add<TH2>("etaItsTpcSel", ";#it{p} (GeV/#it{c});#eta^{TPC} - #eta^{ITS}", HistType::kTH2F, {ptAxis, etaResAxis});
      histos.add<TH2>("phiItsTpcSel", ";#it{p} (GeV/#it{c});#phi^{TPC} - #phi^{ITS} (rad)", HistType::kTH2F, {ptAxis, phiResAxis});

      std::string binLabelsTag[]{"hasITS && hasTPC", "tracking", "PID", "v0 mass", "dcaV0dau", "cosPA", "dcaXYZ", "V0radius"};
      for (int iB{0}; iB < hTagCuts->GetNbinsX(); ++iB) {
        hTagCuts->GetXaxis()->SetBinLabel(iB + 1, binLabelsTag[iB].data());
      }
      for (int iB{0}; iB < hPiRecMass->GetAxis(0)->GetNbins(); ++iB) {
        hPiRecMass->GetAxis(0)->SetBinLabel(iB + 1, binLabels[iB].data());
      }
      hPiRecMass->GetAxis(0)->SetBinLabel(8, "ITS w/ TPC leg");

      if (doprocessTagAndProbeMC) {
        std::string binLabelsTpc[]{"hasTPCsegment", "foundTPCsegment", "allFoundTPCsegment", "foundTPCsegment (w/ fake)"};
        hTpcSegment = histos.add<TH2>("tpcSegment", ";;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {tpcAxis, ptAxis});
        for (int iB{0}; iB < hTpcSegment->GetNbinsX(); ++iB) {
          hTpcSegment->GetXaxis()->SetBinLabel(iB + 1, binLabelsTpc[iB].data());
        }
        histos.add<TH2>("pTpcIts", ";#it{p}^{ITS} (GeV/#it{c});#it{p}^{TPC} - #it{p}^{ITS} (GeV/#it{c});Entries", HistType::kTH2F, {ptAxis, ptResAxis});
        histos.add<TH2>("ptTpcIts", ";#it{p}^{ITS} (GeV/#it{c});#it{p}^{TPC}_{T} - #it{p}^{ITS}_{T} (GeV/#it{c});Entries", HistType::kTH2F, {ptAxis, ptResAxis});
        histos.add<TH2>("etaTpcIts", ";#it{p}^{ITS} (GeV/#it{c});#eta^{TPC} - #eta^{ITS};Entries", HistType::kTH2F, {ptAxis, etaResAxis});
        histos.add<TH2>("phiTpcIts", ";#it{p}^{ITS} (GeV/#it{c});#phi^{TPC} - #phi^{ITS} (rad);Entries", HistType::kTH2F, {ptAxis, phiResAxis});
        histos.add<TH2>("massTpc", ";#it{p}^{ITS} (GeV/#it{c});#it{M}^{TPC} (GeV/#it{c}^{2});Entries", HistType::kTH2F, {ptAxis, massK0sAxis});
        histos.add<TH2>("massTpcIts", ";#it{p}^{ITS} (GeV/#it{c});#it{M}^{TPC} - #it{M}^{ITS} (GeV/#it{c}^{2});Entries", HistType::kTH2F, {ptAxis, massResAxis});
        histos.add<TH2>("cosPaTpc", ";#it{p}^{ITS} (GeV/#it{c});cos#theta_{p}", HistType::kTH2F, {ptAxis, cosPaAxis});
        histos.add<TH3>("ptEtaPhiTpcIts", ";#it{p}^{TPC}_{T} - #it{p}^{ITS}_{T} (GeV/#it{c});#eta^{TPC} - #eta^{ITS};#phi^{TPC} - #phi^{ITS} (rad)", HistType::kTH3F, {ptResAxis, etaResAxis, phiResAxis});
        histos.add<TH1>("collTpcIts", ";ID_{coll}^{TPC} - ID_{coll}^{ITS};Entries", HistType::kTH1F, {collIdResAxis});
        histos.add<TH1>("collTpcV0", ";ID_{coll}^{TPC} - ID_{coll}^{V0};Entries", HistType::kTH1F, {collIdResAxis});
        histos.add<TH1>("timeTpcIts", ";(#it{t}^{TPC} - #it{t}^{ITS}) / #sigma (a.u.);Entries", HistType::kTH1F, {timeResAxis});
        histos.add<TH1>("timeTpcItsNoNorm", ";(#it{t}^{TPC} - #it{t}^{ITS}) (ns);Entries", HistType::kTH1F, {timeResAxisNoNorm});
      }

      histos.add<TH2>("detRespMatrix", ";#it{N}_{gen};#it{N}_{rec}", HistType::kTH2F, {nGenRecAxis, nGenRecAxis});
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
  void fillHistTrack(T const& track, std::shared_ptr<THnBase> hist, float const& y, float const& z = 1, float const& t = 1)
  {
    bool itsAccept = !(track.itsChi2NCl() > 36. || track.itsNCls() < 4);
    bool tpcAccept = !(track.tpcCrossedRowsOverFindableCls() < 0.8 || track.tpcNClsCrossedRows() < 70 || track.tpcChi2NCl() > 4. || track.tpcNClsFound() < 90);
    if (track.hasITS()) {
      hist->Fill(std::array<double, 4>{1., y, z, t}.data());
    }
    if (track.hasITS() && itsAccept && !track.hasTPC()) {
      hist->Fill(std::array<double, 4>{2., y, z, t}.data());
    }
    if (track.hasTPC() && tpcAccept) {
      hist->Fill(std::array<double, 4>{3., y, z, t}.data());
      if (!track.hasITS()) {
        hist->Fill(std::array<double, 4>{4., y, z, t}.data());
      }
      if (track.hasITS() && itsAccept) {
        hist->Fill(std::array<double, 4>{5., y, z, t}.data());
      }
      if (track.hasTOF() && !track.hasITS()) {
        hist->Fill(std::array<double, 4>{6., y, z, t}.data());
      }
    }
  }

  template <class T>
  void fillTagAndProbe(SelCollisions::iterator const& collision, aod::V0s const& V0s, TracksFull const& tracks, SelCollisions const&)
  {
    float nGenRec[]{0.f, 0.f};
    auto tpcTracks = useTpcTracksFromSameColl ? tracks.sliceBy(perCollisionTracks, collision.globalIndex()) : tracks;
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
      std::array<float, 2> dcaInfo;

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

      nGenRec[0] += 1.f;

      auto trackPt = probeTrack.sign() * std::hypot(momProbe[0], momProbe[1]);
      auto trackP = probeTrack.sign() * std::hypot(trackPt, momProbe[2]);
      auto trackEta = probeTrackCov.getEta();

      ProbeTrack probe;
      probe.globalIndex = probeTrack.globalIndex();
      probe.globalIndexTpc = 0;
      probe.globalIndexTag = tagTrack.globalIndex();
      probe.detectorMap = 0;
      probe.p = trackP;
      probe.pt = trackPt;
      probe.massTagProbe = massV0;
      probe.vtx0 = vtx[0];
      probe.vtx1 = vtx[1];
      probe.vtx2 = vtx[2];
      probe.time = probeTrack.trackTime();
      probe.timeRes = probeTrack.trackTimeRes();
      probe.collIndex = probeTrack.collisionId();
      if (probeTrack.has_collision()) {
        auto collisionIts = probeTrack.template collision_as<SelCollisions>();
        probe.bcIndex = int64_t(collisionIts.bcId());
      } else {
        continue; // TODO: check ambiguous tracks (?)
      }

      if (probeTrack.hasITS() && !probeTrack.hasTPC() && findTpcLeg) {
        auto acceptIts = !(probeTrack.itsChi2NCl() > 36. || probeTrack.itsNCls() < 4);
        if (!acceptIts) {
          continue;
        }
        probe.detectorMap = probeTrack.detectorMap();

        std::array<float, 3> momTpc;
        for (auto& tpcTrack : tpcTracks) {
          if (std::abs(tpcTrack.collisionId() - probeTrack.collisionId()) > collIdWindow) {
            continue;
          }
          if (!useCollisionWindow) {
            if (!tpcTrack.has_collision()) {
              continue;
            }
            auto collisionTpc = tpcTrack.template collision_as<SelCollisions>();
            auto bcTpc = int64_t(collisionTpc.bcId());
            float tdiff = (bcTpc - probe.bcIndex) * LHCBunchSpacingNS + tpcTrack.trackTime() - probe.time;
            float nsigmaT = tdiff / std::sqrt(std::pow(tpcTrack.trackTimeRes(), 2) + std::pow(probe.timeRes, 2));
            if (nsigmaT > trackTimingCut) {
              continue;
            }
          }
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

          std::array<float, 2> dcaInfo;
          auto tpcTrackCov = getTrackParCov(tpcTrack);
          if (propToTPCinnerWall) {
            o2::base::Propagator::Instance()->PropagateToXBxByBz(probeTrackCov, 70.f, 1.f, 2.f, fitter.getMatCorrType());
            o2::base::Propagator::Instance()->PropagateToXBxByBz(tpcTrackCov, 70.f, 1.f, 2.f, fitter.getMatCorrType());
          } else {
            o2::base::Propagator::Instance()->propagateToDCABxByBz({static_cast<float>(vtx[0]), static_cast<float>(vtx[1]), static_cast<float>(vtx[2])}, tpcTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
          }
          // tpcTrackCov.update({probeTrack.y(), probeTrack.z()}, {std::pow(probeTrack.sigmaY(), 2), probeTrack.rhoZY() * probeTrack.sigmaY() * probeTrack.sigmaZ(), std::pow(probeTrack.sigmaZ(), 2)});
          // o2::base::Propagator::Instance()->propagateToDCABxByBz({static_cast<float>(vtx[0]), static_cast<float>(vtx[1]), static_cast<float>(vtx[2])}, tpcTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);

          probe.pProp = probeTrackCov.getP();
          probe.ptProp = probeTrackCov.getPt();
          probe.etaProp = probeTrackCov.getEta();
          probe.phiProp = probeTrackCov.getPhi();

          bool acceptTrackPt;
          bool acceptTrackEta;
          bool acceptTrackPhi;
          if (setMomFlatWindows) {
            acceptTrackPt = std::abs(tpcTrackCov.getPt() - probeTrackCov.getPt()) < ptWindow;
            acceptTrackEta = std::abs(tpcTrackCov.getEta() - probeTrackCov.getEta()) < etaWindow;
            acceptTrackPhi = std::abs(tpcTrackCov.getPhi() - probeTrackCov.getPhi()) < phiWindow;
          } else {
            double ptWindowCustom = cfgDeltaPt->get("p0") + cfgDeltaPt->get("p1") * probe.pProp + cfgDeltaPt->get("p2") * std::pow(probe.pProp, 2);
            double etaWindowCustom = cfgDeltaEta->get("p0") * std::pow(probe.pProp, cfgDeltaEta->get("p1"));
            double phiWindowCustom = cfgDeltaPhi->get("p0") + cfgDeltaPhi->get("p1") * std::exp(cfgDeltaPhi->get("p2") * probe.pProp);

            acceptTrackPt = std::abs(tpcTrackCov.getPt() - probeTrackCov.getPt()) < nSigmaDeltaPt * ptWindowCustom;
            acceptTrackEta = std::abs(tpcTrackCov.getEta() - probeTrackCov.getEta()) < nSigmaDeltaEta * etaWindowCustom;
            acceptTrackPhi = std::abs(tpcTrackCov.getPhi() - probeTrackCov.getPhi()) < nSigmaDeltaPhi * phiWindowCustom;
          }
          bool acceptTpcTrack = acceptTrackPt && acceptTrackEta && acceptTrackPhi;

          // LOGF(debug, "idx = %lld, Dpt = %f, Deta = %f, Dphi = %f", tpcTrack.globalIndex(), std::abs(tpcTrackCov.getPt() - propTrackProbe.getPt()), std::abs(tpcTrackCov.getEta() - propTrackProbe.getEta()), std::abs(tpcTrackCov.getPhi() - propTrackProbe.getPhi()));

          if (acceptTpcTrack) {
            histos.fill(HIST("ptItsTpcSel"), trackP, tpcTrackCov.getPt() - probeTrackCov.getPt());
            histos.fill(HIST("etaItsTpcSel"), trackP, tpcTrackCov.getEta() - probeTrackCov.getEta());
            histos.fill(HIST("phiItsTpcSel"), trackP, tpcTrackCov.getPhi() - probeTrackCov.getPhi());

            if (propToTPCinnerWall) {
              o2::base::Propagator::Instance()->propagateToDCABxByBz({static_cast<float>(vtx[0]), static_cast<float>(vtx[1]), static_cast<float>(vtx[2])}, tpcTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
            }
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

              histos.fill(HIST("massTagTpc"), trackP, massTpcLeg);
              probe.globalIndexTpc = tpcTrack.globalIndex();
              break;
            }
          }
        }
      }

      hPiRecMass->Fill(0., trackPt, trackEta, massV0);
      if (probe.globalIndexTpc > 0) {
        hPiRecMass->Fill(7., trackPt, trackEta, massV0);
      }
      fillHistTrack(probeTrack, hPiRecMass, trackPt, trackEta, massV0);
      if (std::abs(massV0 - o2::constants::physics::MassKaonNeutral) < massMax) {
        probeTracks.push_back(probe);
        hPiRec->Fill(0., trackPt, trackEta);
        fillHistTrack(probeTrack, hPiRec, trackPt, trackEta);
        if (probe.globalIndexTpc > 0) {
          hPiRec->Fill(7., trackPt, trackEta);
        }
        bool itsAccept = !(probeTrack.itsChi2NCl() > 36. || probeTrack.itsNCls() < 4);
        bool tpcAccept = !(probeTrack.tpcCrossedRowsOverFindableCls() < 0.8 || probeTrack.tpcNClsCrossedRows() < 70 || probeTrack.tpcChi2NCl() > 4. || probeTrack.tpcNClsFound() < 90);
        if (probeTrack.hasITS() && probeTrack.hasTPC() && itsAccept && tpcAccept) {
          nGenRec[1] += 1.f;
        }
      }
    }
    histos.fill(HIST("detRespMatrix"), nGenRec[0], nGenRec[1]);
  }

  void fillProbeMC(SelCollisions::iterator const& collision, TracksFull const& tracks, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC, SelCollisions const&)
  {
    auto tracks_thisEvent = useTpcTracksFromSameColl ? tracks.sliceBy(perCollisionTracks, collision.globalIndex()) : tracks;

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
        std::array<float, 2> dcaInfo;
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
            if (probeTrack.globalIndexTpc > 0) {
              hTpcSegment->Fill(3., probeTrack.pt);
            }

            LOGF(debug, "globalID = %lld, probeCollId = %d, tpcCollId = %d", collision.globalIndex(), probeTrack.collIndex, tpcTrack.collisionId());
            histos.fill(HIST("collTpcIts"), tpcTrack.collisionId() - probeTrack.collIndex);
            histos.fill(HIST("collTpcV0"), tpcTrack.collisionId() - collision.globalIndex());

            auto collisionTpc = tpcTrack.template collision_as<SelCollisions>();
            float tdiff = (int64_t(collisionTpc.bcId()) - probeTrack.bcIndex) * LHCBunchSpacingNS + tpcTrack.trackTime() - probeTrack.time;
            float nsigmaT = tdiff / std::sqrt(std::pow(tpcTrack.trackTimeRes(), 2) + std::pow(probeTrack.timeRes, 2));
            histos.fill(HIST("timeTpcIts"), nsigmaT);
            histos.fill(HIST("timeTpcItsNoNorm"), tdiff);

            auto trackCov = getTrackParCov(tpcTrack);
            std::array<float, 2> dcaInfo;
            if (propToTPCinnerWall) {
              o2::base::Propagator::Instance()->PropagateToXBxByBz(trackCov, 70.f, 1.f, 2.f, fitter.getMatCorrType());
            } else {
              o2::base::Propagator::Instance()->propagateToDCABxByBz({probeTrack.vtx0, probeTrack.vtx1, probeTrack.vtx2}, trackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
            }

            histos.fill(HIST("pTpcIts"), probeTrack.p, trackCov.getP() - probeTrack.pProp);
            histos.fill(HIST("ptTpcIts"), probeTrack.p, trackCov.getPt() - probeTrack.ptProp);
            histos.fill(HIST("etaTpcIts"), probeTrack.p, trackCov.getEta() - probeTrack.etaProp);
            histos.fill(HIST("phiTpcIts"), probeTrack.p, trackCov.getPhi() - probeTrack.phiProp);
            histos.fill(HIST("ptEtaPhiTpcIts"), trackCov.getPt() - probeTrack.ptProp, trackCov.getEta() - probeTrack.etaProp, trackCov.getPhi() - probeTrack.phiProp);

            if (propToTPCinnerWall) {
              o2::base::Propagator::Instance()->propagateToDCABxByBz({probeTrack.vtx0, probeTrack.vtx1, probeTrack.vtx2}, trackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
            }
            std::array<float, 3> v0Mom, momTpc;
            trackCov.getPxPyPzGlo(momTpc);
            auto massTpcLeg = invMass2Body(v0Mom, momTag, momTpc, piMass, piMass);

            histos.fill(HIST("massTpc"), probeTrack.p, massTpcLeg);
            histos.fill(HIST("massTpcIts"), probeTrack.p, massTpcLeg - probeTrack.massTagProbe);

            LOGF(debug, "MC::DEBUG: idx = %lld, p = %f, Dpt = %f, Deta = %f, Dphi = %f", tpcTrack.globalIndex(), probeTrack.pt, std::abs(trackCov.getPt() - probeTrack.pt * probeTrack.pt / std::abs(probeTrack.pt)), std::abs(trackCov.getEta() - probeTrack.etaProp), std::abs(trackCov.getPhi() - probeTrack.phiProp));

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

  void processTagAndProbe(SelCollisions const& collisions, aod::V0s const& V0s, TracksFull const& tracks, aod::BCsWithTimestamps const&)
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
      fillTagAndProbe<TracksFull>(collision, V0Table_thisCollision, tracks, collisions);
    }
  }
  PROCESS_SWITCH(efficiencyQA, processTagAndProbe, "Tag and probe analysis", true);

  void processTagAndProbeMC(SelCollisions const& collisions, aod::V0s const& V0s, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC)
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
      fillTagAndProbe<TracksFull>(collision, V0Table_thisCollision, tracks, collisions);

      fillProbeMC(collision, tracks, trackLabelsMC, particlesMC, collisions);
    }
  }
  PROCESS_SWITCH(efficiencyQA, processTagAndProbeMC, "Tag and probe analysis on MC", false);

  void processMcTracks(SelCollisions const& collisions, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC)
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
          const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};

          auto trackParCov = getTrackParCov(track);
          std::array<float, 2> dcaInfo;
          o2::base::Propagator::Instance()->propagateToDCA(collVtx, trackParCov, d_bz, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfo);

          auto trackPt = track.sign() * trackParCov.getPt();
          auto trackEta = trackParCov.getEta();
          if (mcTrack.isPhysicalPrimary()) {
            fillHistTrack(track, hPiRec, trackPt, trackEta);
            fillHistTrack(track, hPtRes, track.sign() * trackParCov.getPt(), trackParCov.getPt() - mcTrack.pt());
            fillHistTrack(track, hEtaRes, track.sign() * trackParCov.getPt(), trackParCov.getEta() - mcTrack.eta());
            fillHistTrack(track, hPhiRes, track.sign() * trackParCov.getPt(), trackParCov.getPhi() - mcTrack.phi());
          } else {
            for (auto& mother : mcTrack.template mothers_as<aod::McParticles>()) {
              if (mother.pdgCode() != 310) {
                continue;
              }
              auto radius = std::hypot(mcTrack.vx(), mcTrack.vy());
              if (radius > v0radiusMax) {
                continue;
              }
              fillHistTrack(track, hPiRecSec, trackPt, trackEta);
              break;
            }
          }
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
      if (partMC.isPhysicalPrimary()) {
        hPiRec->Fill(0., pdgCode / std::abs(pdgCode) * partMC.pt(), partMC.eta());
      } else {
        for (auto& mother : partMC.template mothers_as<aod::McParticles>()) {
          if (mother.pdgCode() != 310) {
            continue;
          }
          auto radius = std::hypot(partMC.vx(), partMC.vy());
          if (radius > v0radiusMax) {
            continue;
          }
          hPiRecSec->Fill(0., pdgCode / std::abs(pdgCode) * partMC.pt(), partMC.eta());
          break;
        }
      }
    }
  }
  PROCESS_SWITCH(efficiencyQA, processMcTracks, "MC tracks analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<efficiencyQA>(cfgc)};
}
