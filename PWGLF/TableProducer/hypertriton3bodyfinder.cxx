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

// This 3-body method is not recommended due to high cost of computing resources
// author: yuanzhe.wang@cern.ch

#include <cmath>
#include <array>
#include <cstdlib>
#include <iterator>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "EventFiltering/filterTables.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe>;
using FullTracksExtMCIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

using LabeledTracks = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

namespace o2::aod
{
namespace v0goodpostrack
{
DECLARE_SOA_INDEX_COLUMN_FULL(GoodTrack, goodTrack, int, Tracks, "_GoodTrack");
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);
} // namespace v0goodpostrack
DECLARE_SOA_TABLE(V0GoodPosTracks, "AOD", "V0GOODPOSTRACKS", o2::soa::Index<>, v0goodpostrack::GoodTrackId, v0goodpostrack::CollisionId, v0goodpostrack::DCAXY);
namespace v0goodnegtrack
{
DECLARE_SOA_INDEX_COLUMN_FULL(GoodTrack, goodTrack, int, Tracks, "_GoodTrack");
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);
} // namespace v0goodnegtrack
DECLARE_SOA_TABLE(V0GoodNegTracks, "AOD", "V0GOODNEGTRACKS", o2::soa::Index<>, v0goodnegtrack::GoodTrackId, v0goodnegtrack::CollisionId, v0goodnegtrack::DCAXY);
namespace v0goodtrack
{
DECLARE_SOA_INDEX_COLUMN_FULL(GoodTrack, goodTrack, int, Tracks, "_GoodTrack");
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);
} // namespace v0goodtrack
DECLARE_SOA_TABLE(V0GoodTracks, "AOD", "V0GOODTRACKS", o2::soa::Index<>, v0goodtrack::GoodTrackId, v0goodtrack::CollisionId, v0goodtrack::DCAXY);
} // namespace o2::aod

struct trackprefilter {
  HistogramRegistry registry{
    "registry",
    {
      {"hCrossedRows", "hCrossedRows", {HistType::kTH1F, {{50, 0.0f, 200.0f}}}},
      {"hGoodTrackCount", "hGoodTrackCount", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
      {"hGoodPosTrackCount", "hGoodPosTrackCount", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
      {"hGoodNegTrackCount", "hGoodNegTrackCount", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
    },
  };

  // change the dca cut for helium3
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<int> tpcrefit{"tpcrefit", 0, "demand TPC refit"};

  Produces<aod::V0GoodPosTracks> v0GoodPosTracks;
  Produces<aod::V0GoodNegTracks> v0GoodNegTracks;
  Produces<aod::V0GoodTracks> v0GoodTracks;

  void process(aod::Collision const& collision,
               FullTracksExtIU const& tracks)
  {
    for (auto& t0 : tracks) {
      registry.fill(HIST("hGoodTrackCount"), 0.5);
      registry.fill(HIST("hCrossedRows"), t0.tpcNClsCrossedRows());
      if (tpcrefit) {
        if (!(t0.trackType() & o2::aod::track::TPCrefit)) {
          continue; // TPC refit
        }
      }
      registry.fill(HIST("hGoodTrackCount"), 1.5);
      if (t0.tpcNClsCrossedRows() < mincrossedrows) {
        continue;
      }
      registry.fill(HIST("hGoodTrackCount"), 2.5);
      if (t0.signed1Pt() > 0.0f) {
        v0GoodPosTracks(t0.globalIndex(), collision.globalIndex(), t0.dcaXY());
        registry.fill(HIST("hGoodPosTrackCount"), 0.5);
        registry.fill(HIST("hGoodTrackCount"), 3.5);
      }
      if (t0.signed1Pt() < 0.0f) {
        v0GoodNegTracks(t0.globalIndex(), collision.globalIndex(), t0.dcaXY());
        registry.fill(HIST("hGoodNegTrackCount"), 0.5);
        registry.fill(HIST("hGoodTrackCount"), 3.5);
      }
      v0GoodTracks(t0.globalIndex(), collision.globalIndex(), t0.dcaXY());
    }
  }
};

template <class TMCTrackTo, typename TMCParticle>
bool is3bodyDecayedH3L(TMCParticle const& particle)
{
  if (particle.pdgCode() != 1010010030 && particle.pdgCode() != -1010010030) {
    return false;
  }
  bool haveProton = false, havePion = false, haveDeuteron = false;
  bool haveAntiProton = false, haveAntiPion = false, haveAntiDeuteron = false;
  for (auto& mcparticleDaughter : particle.template daughters_as<aod::McParticles>()) {
    if (mcparticleDaughter.pdgCode() == 2212)
      haveProton = true;
    if (mcparticleDaughter.pdgCode() == -2212)
      haveAntiProton = true;
    if (mcparticleDaughter.pdgCode() == 211)
      havePion = true;
    if (mcparticleDaughter.pdgCode() == -211)
      haveAntiPion = true;
    if (mcparticleDaughter.pdgCode() == 1000010020)
      haveDeuteron = true;
    if (mcparticleDaughter.pdgCode() == -1000010020)
      haveAntiDeuteron = true;
  }
  if (haveProton && haveAntiPion && haveDeuteron && particle.pdgCode() == 1010010030) {
    return true;
  } else if (haveAntiProton && havePion && haveAntiDeuteron && particle.pdgCode() == -1010010030) {
    return true;
  }
  return false;
}

struct hypertriton3bodyFinder {

  Produces<aod::StoredVtx3BodyDatas> vtx3bodydata;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Configurables
  Configurable<bool> UseCFFilter{"UseCFFilter", true, "Reject event without CF LD trigger"};
  Configurable<bool> RejectBkgInMC{"RejectBkgInMC", false, "Reject fake 3-body pairs in MC check"};
  Configurable<bool> KeepSignalInMC{"KeepSignalInMC", false, "Reject fake 3-body pairs in MC check"};

  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};

  // Selection criteria
  Configurable<float> minRToMeanVertex = {"minRToMeanVertex", 0.5, ""};     ///< min radial distance of V0 from beam line (mean vertex)
                                                                            // Configurable<float> causalityRTolerance = {"causalityRTolerance", 1., ""}; ///< V0 radius cannot exceed its contributors minR by more than this value
  Configurable<float> maxV0ToProngsRDiff = {"maxV0ToProngsRDiff", 50., ""}; ///< V0 radius cannot be lower than this ammount wrt minR of contributors
  Configurable<float> minPtV0 = {"minPtV0", 0.5, ""};                       ///< v0 minimum pT
  Configurable<float> maxTglV0 = {"maxTglV0", 2., ""};                      ///< maximum tgLambda of V0
  Configurable<float> maxDCAXY2ToMeanVertex3bodyV0 = {"maxDCAXY2ToMeanVertex3bodyV0", 2 * 2, ""};
  Configurable<float> minCosPAXYMeanVertex3bodyV0 = {"minCosPAXYMeanVertex3bodyV0", 0.9, ""}; ///< min cos of PA to beam line (mean vertex) in tr. plane for 3body V0 cand.
  Configurable<float> minCosPA3bodyV0 = {"minCosPA3bodyV0", 0.8, ""};                         // min cos of PA to PV for 3body V0

  // for 3 body reconstructed Vertex
  Configurable<float> minbachPt = {"minbachPt", 0.6, ""};           ///< Minimum bachelor Pt
  Configurable<float> maxRDiff3bodyV0 = {"maxRDiff3bodyV0", 3, ""}; ///< Maximum difference between V0 and 3body radii
  Configurable<float> minPt3Body = {"minPt3Body", 0.01, ""};        // minimum pT of 3body Vertex
  Configurable<float> maxTgl3Body = {"maxTgl3Body", 2, ""};         // maximum tgLambda of 3body Vertex
  Configurable<float> minCosPA3body = {"minCosPA3body", 0.8, ""};   // min cos of PA to PV for 3body Vertex

  // for DCA
  Configurable<float> dcavtxdau{"dcavtxdau", 2.0, "DCA Vtx Daughters"};
  Configurable<bool> d_UseH3LDCACut{"d_UseH3LDCACut", true, "Use Cuts for H3L DCA to PV"};
  Configurable<float> maxDCAXY3Body{"maxDCAXY3Body", 0.5, "DCAXY H3L to PV"}; // max DCA of 3 body decay to PV in XY
  Configurable<float> maxDCAZ3Body{"maxDCAZ3Body", 1.0, "DCAZ H3L to PV"};    // max DCA of 3 body decay to PV in Z

  Configurable<int> useMatCorrType{"useMatCorrType", 2, "0: none, 1: TGeo, 2: LUT"};
  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Preslice<aod::V0GoodPosTracks> perCollisionGoodPosTracks = o2::aod::v0goodpostrack::collisionId;
  Preslice<aod::V0GoodNegTracks> perCollisionGoodNegTracks = o2::aod::v0goodnegtrack::collisionId;
  Preslice<aod::V0GoodTracks> perCollisionGoodTracks = o2::aod::v0goodtrack::collisionId;
  Preslice<aod::V0s> perCollisionV0s = o2::aod::v0::collisionId;
  Preslice<aod::V0Datas> perCollisionV0Datas = o2::aod::v0data::collisionId;

  // Helper struct to pass V0 information
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hDauTrackCounter", "hDauTrackCounter", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
      {"hV0Counter", "hV0Counter", {HistType::kTH1F, {{8, -0.5f, 7.5f}}}},
      {"hTrueV0Counter", "hTrueV0Counter", {HistType::kTH1F, {{8, -0.5f, 7.5f}}}},
      {"hVtx3BodyCounter", "hVtx3BodyCounter", {HistType::kTH1F, {{9, -0.5f, 8.5f}}}},
      {"hTrueVtx3BodyCounter", "hTrueVtx3BodyCounter", {HistType::kTH1F, {{9, -0.5f, 8.5f}}}},
      {"hVirtLambaCounter", "hVirtualLambaCounter", {HistType::kTH1F, {{6, -0.5f, 5.5f}}}},
      {"hCFFilteredVirtLambaCounter", "hCFFilteredVirtLambaCounter", {HistType::kTH1F, {{6, -0.5f, 5.5f}}}},
    },
  };

  //------------------------------------------------------------------
  // Fill stats histograms
  enum v0step { kV0All = 0,
                kV0hasSV,
                kV0Radius,
                kV0Pt,
                kV0TgLamda,
                kV0InvMass,
                kV0DcaXY,
                kV0CosPA,
                kNV0Steps };
  enum vtxstep { kVtxAll = 0,
                 kVtxbachPt,
                 kVtxhasSV,
                 kVtxRadius,
                 kVtxPt,
                 kVtxTgLamda,
                 kVtxCosPA,
                 kVtxDcaDau,
                 kVtxDcaH3L,
                 kNVtxSteps };

  // Helper struct to do bookkeeping of building parameters
  struct {
    std::array<int32_t, kNV0Steps> v0stats;
    std::array<int32_t, kNV0Steps> truev0stats;
    std::array<int32_t, kNVtxSteps> vtxstats;
    std::array<int32_t, kNVtxSteps> truevtxstats;
    std::array<int32_t, 12> virtLambdastats;
  } statisticsRegistry;

  void resetHistos()
  {
    for (Int_t ii = 0; ii < kNV0Steps; ii++) {
      statisticsRegistry.v0stats[ii] = 0;
      statisticsRegistry.truev0stats[ii] = 0;
    }
    for (Int_t ii = 0; ii < kNVtxSteps; ii++) {
      statisticsRegistry.vtxstats[ii] = 0;
      statisticsRegistry.truevtxstats[ii] = 0;
    }
    for (Int_t ii = 0; ii < 12; ii++) {
      statisticsRegistry.virtLambdastats[ii] = 0;
    }
  }

  void fillHistos()
  {
    for (Int_t ii = 0; ii < kNV0Steps; ii++) {
      registry.fill(HIST("hV0Counter"), ii, statisticsRegistry.v0stats[ii]);
      registry.fill(HIST("hTrueV0Counter"), ii, statisticsRegistry.truev0stats[ii]);
    }
    for (Int_t ii = 0; ii < kNVtxSteps; ii++) {
      registry.fill(HIST("hVtx3BodyCounter"), ii, statisticsRegistry.vtxstats[ii]);
      registry.fill(HIST("hTrueVtx3BodyCounter"), ii, statisticsRegistry.truevtxstats[ii]);
    }
    for (Int_t ii = 0; ii < 3; ii++) {
      registry.fill(HIST("hVirtLambaCounter"), ii, statisticsRegistry.virtLambdastats[ii]);
      registry.fill(HIST("hVirtLambaCounter"), ii + 3, statisticsRegistry.virtLambdastats[ii + 3]);
      registry.fill(HIST("hCFFilteredVirtLambaCounter"), ii, statisticsRegistry.virtLambdastats[ii + 6]);
      registry.fill(HIST("hCFFilteredVirtLambaCounter"), ii + 3, statisticsRegistry.virtLambdastats[ii + 9]);
    }
  }

  // v0, vtx, and virtual Lambda statiscs
  void FillV0Counter(int kn, bool istrue = false)
  {
    statisticsRegistry.v0stats[kn]++;
    if (istrue) {
      statisticsRegistry.truev0stats[kn]++;
    }
  }
  void FillVtxCounter(int kn, bool istrue = false)
  {
    statisticsRegistry.vtxstats[kn]++;
    if (istrue) {
      statisticsRegistry.truevtxstats[kn]++;
    }
  }
  //------------------------------------------------------------------

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::vertexing::DCAFitterN<2> fitter;
  o2::vertexing::DCAFitterN<3> fitter3body;

  void init(InitContext& context)
  {
    resetHistos();
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  // could be changed later
    maxStep = 2.00f; // could be changed later

    TString DauCounterbinLabel[3] = {"Proton", "Pion", "Deuteron"};
    TString V0CounterbinLabel[8] = {"Total", "hasSV", "V0R", "V0Pt", "TgLambda", "V0Mass", "DcaXY", "CosPA"};
    TString VtxCounterbinLabel[9] = {"Total", "bachPt", "hasSV", "VtxR", "VtxPt", "TgLambda", "CosPA", "DcaDau", "DcaH3L"};
    for (int i{0}; i < 3; i++) {
      registry.get<TH1>(HIST("hDauTrackCounter"))->GetXaxis()->SetBinLabel(i + 1, DauCounterbinLabel[i]);
    }
    for (int i{0}; i < kNV0Steps; i++) {
      registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(i + 1, V0CounterbinLabel[i]);
      registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(i + 1, V0CounterbinLabel[i]);
    }
    for (int i{0}; i < kNVtxSteps; i++) {
      registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(i + 1, VtxCounterbinLabel[i]);
      registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(i + 1, VtxCounterbinLabel[i]);
    }

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // Set 2-body fitter and 3-body fitter3body
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.); //->maxRIni3body
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(d_UseAbsDCA);
    fitter3body.setPropagateToPCA(true);
    fitter3body.setMaxR(200.); //->maxRIni3body
    fitter3body.setMinParamChange(1e-3);
    fitter3body.setMinRelChi2Change(0.9);
    fitter3body.setMaxDZIni(1e9);
    fitter3body.setMaxChi2(1e9);
    fitter3body.setUseAbsDCA(d_UseAbsDCA);

    // Material correction in the DCA fitter
    o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
    if (useMatCorrType == 1) {
      LOGF(info, "TGeo correction requested, loading geometry");
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
      }
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    }
    if (useMatCorrType == 2) {
      LOGF(info, "LUT correction requested, loading LUT");
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    }
    fitter.setMatCorrType(matCorr);
    fitter3body.setMatCorrType(matCorr);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      fitter.setBz(d_bz);
      fitter3body.setBz(d_bz);
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
    // Set magnetic field value once known
    fitter.setBz(d_bz);
    fitter3body.setBz(d_bz);

    if (useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized
      // (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
  }

  //------------------------------------------------------------------
  // Check the info of good tracks
  template <class TTrackTo, typename TGoodTrackTable>
  void CheckGoodTracks(TGoodTrackTable const& dGoodtracks, aod::McParticles const& mcparticles)
  {
    for (auto& goodtrackid : dGoodtracks) {
      auto goodtrack = goodtrackid.template goodTrack_as<TTrackTo>();
      if (!goodtrack.has_mcParticle()) {
        continue;
      }
      auto mcgoodtrack = goodtrack.template mcParticle_as<aod::McParticles>();
      if (!mcgoodtrack.has_mothers()) {
        continue;
      }
      bool flag_H3L = false;
      for (auto& mothertrack : mcgoodtrack.template mothers_as<aod::McParticles>()) {
        if (is3bodyDecayedH3L<aod::McParticles>(mothertrack)) {
          flag_H3L = true;
        }
      }
      if (flag_H3L && std::abs(mcgoodtrack.pdgCode()) == 2212) {
        registry.fill(HIST("hDauTrackCounter"), 0.5);
      }
      if (flag_H3L && std::abs(mcgoodtrack.pdgCode()) == 211) {
        registry.fill(HIST("hDauTrackCounter"), 1.5);
      }
      if (flag_H3L && std::abs(mcgoodtrack.pdgCode()) == 1000010020) {
        registry.fill(HIST("hDauTrackCounter"), 2.5);
      }
    }
  }

  o2::dataformats::VertexBase mMeanVertex{{0., 0., 0.}, {0.1 * 0.1, 0., 0.1 * 0.1, 0., 0., 6. * 6.}};
  //------------------------------------------------------------------
  // 3body decay finder
  template <class TTrackTo, typename TCollisionTable, typename TPosTrackTable, typename TNegTrackTable, typename TGoodTrackTable>
  void DecayFinder(TCollisionTable const& dCollision, TPosTrackTable const& dPtracks, TNegTrackTable const& dNtracks, TGoodTrackTable const& dGoodtracks)
  {
    for (auto& t0id : dPtracks) { // FIXME: turn into combination(...)
      auto t0 = t0id.template goodTrack_as<TTrackTo>();
      auto Track0 = getTrackParCov(t0);

      for (auto& t1id : dNtracks) {

        FillV0Counter(kV0All);
        auto t1 = t1id.template goodTrack_as<TTrackTo>();
        auto Track1 = getTrackParCov(t1);
        int nCand = fitter.process(Track0, Track1);
        if (nCand == 0) {
          continue;
        }
        FillV0Counter(kV0hasSV);

        // validate V0 radial position
        // First check closeness to the beam-line as same as SVertexer
        const auto& v0XYZ = fitter.getPCACandidate();
        float dxv0 = v0XYZ[0] - mMeanVertex.getX(), dyv0 = v0XYZ[1] - mMeanVertex.getY(), r2v0 = dxv0 * dxv0 + dyv0 * dyv0;
        float rv0 = std::sqrt(r2v0);
        if (rv0 < minRToMeanVertex) {
          continue;
        }
        FillV0Counter(kV0Radius);

        // Not involved: Get minR with same way in SVertexer
        // float drv0P = rv0 - Track0minR, drv0N = rv0 - Track1minR;

        // check: if the process function finish the propagation
        if (!fitter.isPropagateTracksToVertexDone() && !fitter.propagateTracksToVertex()) {
          continue;
        }

        int cand = 0;
        auto& trPProp = fitter.getTrack(0, cand);
        auto& trNProp = fitter.getTrack(1, cand);
        std::array<float, 3> pP, pN;
        trPProp.getPxPyPzGlo(pP);
        trNProp.getPxPyPzGlo(pN);
        // estimate DCA of neutral V0 track to beamline: straight line with parametric equation
        // x = X0 + pV0[0]*t, y = Y0 + pV0[1]*t reaches DCA to beamline (Xv, Yv) at
        // t = -[ (x0-Xv)*pV0[0] + (y0-Yv)*pV0[1]) ] / ( pT(pV0)^2 )
        // Similar equation for 3D distance involving pV0[2]
        std::array<float, 3> pV0 = {pP[0] + pN[0], pP[1] + pN[1], pP[2] + pN[2]};
        float pt2V0 = pV0[0] * pV0[0] + pV0[1] * pV0[1], prodXYv0 = dxv0 * pV0[0] + dyv0 * pV0[1], tDCAXY = prodXYv0 / pt2V0;
        float p2V0 = pt2V0 + pV0[2] * pV0[2], ptV0 = std::sqrt(pt2V0);
        if (ptV0 < minPtV0) { // pt cut
          continue;
        }
        FillV0Counter(kV0Pt);

        if (pV0[2] / ptV0 > maxTglV0) { // tgLambda cut
          continue;
        }
        FillV0Counter(kV0TgLamda);

        // apply mass selections
        float massV0LambdaHyp = RecoDecay::m(array{array{pP[0], pP[1], pP[2]}, array{pN[0], pN[1], pN[2]}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
        float massV0AntiLambdaHyp = RecoDecay::m(array{array{pP[0], pP[1], pP[2]}, array{pN[0], pN[1], pN[2]}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
        float massMargin = 20 * (0.001 * (1. + 0.5 * ptV0)) + 0.07;
        if (massV0LambdaHyp - o2::constants::physics::MassLambda > massMargin && massV0AntiLambdaHyp - o2::constants::physics::MassLambda > massMargin) {
          continue;
        }
        FillV0Counter(kV0InvMass);

        float dcaX = dxv0 - pV0[0] * tDCAXY, dcaY = dyv0 - pV0[1] * tDCAXY, dca2 = dcaX * dcaX + dcaY * dcaY;
        if (dca2 > maxDCAXY2ToMeanVertex3bodyV0) {
          continue;
        }
        FillV0Counter(kV0DcaXY);

        float cosPAXY = prodXYv0 / std::sqrt(r2v0 * pt2V0);
        if (cosPAXY < minCosPAXYMeanVertex3bodyV0) {
          continue;
        }
        float dx = v0XYZ[0] - dCollision.posX(), dy = v0XYZ[1] - dCollision.posY(), dz = v0XYZ[2] - dCollision.posZ(), prodXYZv0 = dx * pV0[0] + dy * pV0[1] + dz * pV0[2];
        float cosPA = prodXYZv0 / std::sqrt((dx * dx + dy * dy + dz * dz) * p2V0);
        if (cosPA < minCosPA3bodyV0) {
          continue;
        }
        FillV0Counter(kV0CosPA);

        for (auto& t2id : dGoodtracks) {
          if (t2id.globalIndex() == t0id.globalIndex()) {
            continue; // skip the track used by V0
          }
          FillVtxCounter(kVtxAll);

          auto t2 = t2id.template goodTrack_as<TTrackTo>();
          auto track0 = getTrackParCov(t0);
          auto track1 = getTrackParCov(t1);
          auto bach = getTrackParCov(t2);

          if (bach.getPt() < 0.6) {
            continue;
          }
          FillVtxCounter(kVtxbachPt);

          int n3bodyVtx = fitter3body.process(track0, track1, bach);
          if (n3bodyVtx == 0) { // discard this pair
            continue;
          }
          FillVtxCounter(kVtxhasSV);

          int cand3B = 0;
          const auto& vertexXYZ = fitter3body.getPCACandidatePos(cand3B);
          // make sure the cascade radius is smaller than that of the vertex
          float dxc = vertexXYZ[0] - dCollision.posX(), dyc = vertexXYZ[1] - dCollision.posY(), dzc = vertexXYZ[2] - dCollision.posZ(), r2vertex = dxc * dxc + dyc * dyc;
          float rvertex = std::sqrt(r2vertex);
          if (std::abs(rv0 - rvertex) > maxRDiff3bodyV0 || rvertex < minRToMeanVertex) {
            continue;
          }
          FillVtxCounter(kVtxRadius);

          // Not involved: bach.minR - rveretx check

          if (!fitter3body.isPropagateTracksToVertexDone() && !fitter3body.propagateTracksToVertex()) {
            continue;
          }

          auto& tr0 = fitter3body.getTrack(0, cand3B);
          auto& tr1 = fitter3body.getTrack(1, cand3B);
          auto& tr2 = fitter3body.getTrack(2, cand3B);
          std::array<float, 3> p0, p1, p2;
          tr0.getPxPyPzGlo(p0);
          tr1.getPxPyPzGlo(p1);
          tr2.getPxPyPzGlo(p2);
          std::array<float, 3> p3B = {p0[0] + p1[0] + p2[0], p0[1] + p1[1] + p2[1], p0[2] + p1[2] + p2[2]};

          float pt2 = p3B[0] * p3B[0] + p3B[1] * p3B[1], p2candidate = pt2 + p3B[2] * p3B[2];
          float pt = std::sqrt(pt2);
          if (pt < minPt3Body) { // pt cut
            continue;
          }
          FillVtxCounter(kVtxPt);

          if (p3B[2] / pt > maxTgl3Body) { // tgLambda cut
            continue;
          }
          FillVtxCounter(kVtxTgLamda);

          float cosPA = (p3B[0] * dxc + p3B[1] * dyc + p3B[2] * dzc) / std::sqrt(p2candidate * (r2vertex + dzc * dzc));
          if (cosPA < minCosPA3body) {
            continue;
          }
          FillVtxCounter(kVtxCosPA);

          if (fitter3body.getChi2AtPCACandidate() > dcavtxdau) {
            continue;
          }
          FillVtxCounter(kVtxDcaDau);

          // Calculate DCA with respect to the collision associated to the V0, not individual tracks
          gpu::gpustd::array<float, 2> dcaInfo;

          auto Track0Par = getTrackPar(t0);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({dCollision.posX(), dCollision.posY(), dCollision.posZ()}, Track0Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
          auto Track0dcaXY = dcaInfo[0];

          auto Track1Par = getTrackPar(t1);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({dCollision.posX(), dCollision.posY(), dCollision.posZ()}, Track1Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
          auto Track1dcaXY = dcaInfo[0];

          auto Track2Par = getTrackPar(t2);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({dCollision.posX(), dCollision.posY(), dCollision.posZ()}, Track2Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
          auto Track2dcaXY = dcaInfo[0];

          //  Not involved: H3L DCA Check
          vtx3bodydata(
            t0.globalIndex(), t1.globalIndex(), t2.globalIndex(), dCollision.globalIndex(), 0,
            vertexXYZ[0], vertexXYZ[1], vertexXYZ[2],
            p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2],
            fitter3body.getChi2AtPCACandidate(),
            Track0dcaXY, Track1dcaXY, Track2dcaXY);
        }
      }
    }
    fillHistos();
    resetHistos();
  }
  //------------------------------------------------------------------
  // MC 3body decay finder
  template <class TTrackTo, typename TCollisionTable, typename TPosTrackTable, typename TNegTrackTable, typename TGoodTrackTable>
  void DecayFinderMC(TCollisionTable const& dCollision, TPosTrackTable const& dPtracks, TNegTrackTable const& dNtracks, TGoodTrackTable const& dGoodtracks)
  {
    for (auto& t0id : dPtracks) { // FIXME: turn into combination(...)
      auto t0 = t0id.template goodTrack_as<TTrackTo>();
      auto Track0 = getTrackParCov(t0);
      for (auto& t1id : dNtracks) {

        if (t0id.collisionId() != t1id.collisionId()) {
          continue;
        }
        auto t1 = t1id.template goodTrack_as<TTrackTo>();
        auto Track1 = getTrackParCov(t1);

        bool isTrue3bodyV0 = false;
        if (t0.has_mcParticle() && t1.has_mcParticle()) {
          auto t0mc = t0.template mcParticle_as<aod::McParticles>();
          auto t1mc = t1.template mcParticle_as<aod::McParticles>();
          if ((t0mc.pdgCode() == 2212 && t1mc.pdgCode() == -211) || (t0mc.pdgCode() == 211 && t1mc.pdgCode() == -2212)) {
            if (t0mc.has_mothers() && t1mc.has_mothers()) {
              for (auto& t0mother : t0mc.template mothers_as<aod::McParticles>()) {
                for (auto& t1mother : t1mc.template mothers_as<aod::McParticles>()) {
                  if (t0mother.globalIndex() == t1mother.globalIndex() && std::abs(t0mother.pdgCode()) == 1010010030) {
                    isTrue3bodyV0 = true;
                  }
                }
              }
            }
          }
        }

        FillV0Counter(kV0All, isTrue3bodyV0);
        if (!isTrue3bodyV0 && RejectBkgInMC) {
          continue;
        }

        int nCand = fitter.process(Track0, Track1);
        if (nCand == 0) {
          continue;
        }
        FillV0Counter(kV0hasSV, isTrue3bodyV0);

        // validate V0 radial position
        // First check closeness to the beam-line as same as SVertexer
        const auto& v0XYZ = fitter.getPCACandidate();
        float dxv0 = v0XYZ[0] - mMeanVertex.getX(), dyv0 = v0XYZ[1] - mMeanVertex.getY(), r2v0 = dxv0 * dxv0 + dyv0 * dyv0;
        float rv0 = std::sqrt(r2v0);
        if (rv0 < minRToMeanVertex) {
          continue;
        }
        FillV0Counter(kV0Radius, isTrue3bodyV0);

        // Not involved: Get minR with same way in SVertexer
        // float drv0P = rv0 - Track0minR, drv0N = rv0 - Track1minR;

        // check: if the process function finish the propagation
        if (!fitter.isPropagateTracksToVertexDone() && !fitter.propagateTracksToVertex()) {
          continue;
        }

        int cand = 0;
        auto& trPProp = fitter.getTrack(0, cand);
        auto& trNProp = fitter.getTrack(1, cand);
        std::array<float, 3> pP, pN;
        trPProp.getPxPyPzGlo(pP);
        trNProp.getPxPyPzGlo(pN);
        // estimate DCA of neutral V0 track to beamline: straight line with parametric equation
        // x = X0 + pV0[0]*t, y = Y0 + pV0[1]*t reaches DCA to beamline (Xv, Yv) at
        // t = -[ (x0-Xv)*pV0[0] + (y0-Yv)*pV0[1]) ] / ( pT(pV0)^2 )
        // Similar equation for 3D distance involving pV0[2]
        std::array<float, 3> pV0 = {pP[0] + pN[0], pP[1] + pN[1], pP[2] + pN[2]};
        float pt2V0 = pV0[0] * pV0[0] + pV0[1] * pV0[1], prodXYv0 = dxv0 * pV0[0] + dyv0 * pV0[1], tDCAXY = prodXYv0 / pt2V0;
        float p2V0 = pt2V0 + pV0[2] * pV0[2], ptV0 = std::sqrt(pt2V0);
        if (ptV0 < minPtV0) { // pt cut
          continue;
        }
        FillV0Counter(kV0Pt, isTrue3bodyV0);

        if (pV0[2] / ptV0 > maxTglV0) { // tgLambda cut
          continue;
        }
        FillV0Counter(kV0TgLamda, isTrue3bodyV0);

        // apply mass selections
        float massV0LambdaHyp = RecoDecay::m(array{array{pP[0], pP[1], pP[2]}, array{pN[0], pN[1], pN[2]}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
        float massV0AntiLambdaHyp = RecoDecay::m(array{array{pP[0], pP[1], pP[2]}, array{pN[0], pN[1], pN[2]}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
        float massMargin = 20 * (0.001 * (1. + 0.5 * ptV0)) + 0.07;
        if (massV0LambdaHyp - o2::constants::physics::MassLambda > massMargin && massV0AntiLambdaHyp - o2::constants::physics::MassLambda > massMargin) {
          continue;
        }
        FillV0Counter(kV0InvMass, isTrue3bodyV0);

        float dcaX = dxv0 - pV0[0] * tDCAXY, dcaY = dyv0 - pV0[1] * tDCAXY, dca2 = dcaX * dcaX + dcaY * dcaY;
        float cosPAXY = prodXYv0 / std::sqrt(r2v0 * pt2V0);
        if (dca2 > maxDCAXY2ToMeanVertex3bodyV0) {
          continue;
        }
        FillV0Counter(kV0DcaXY, isTrue3bodyV0);

        if (cosPAXY < minCosPAXYMeanVertex3bodyV0) {
          continue;
        }
        float dx = v0XYZ[0] - dCollision.posX(), dy = v0XYZ[1] - dCollision.posY(), dz = v0XYZ[2] - dCollision.posZ(), prodXYZv0 = dx * pV0[0] + dy * pV0[1] + dz * pV0[2];
        float cosPA = prodXYZv0 / std::sqrt((dx * dx + dy * dy + dz * dz) * p2V0);
        if (cosPA < minCosPA3bodyV0) {
          continue;
        }
        FillV0Counter(kV0CosPA, isTrue3bodyV0);

        for (auto& t2id : dGoodtracks) {
          if (t2id.globalIndex() == t0id.globalIndex()) {
            continue; // skip the track used by V0
          }
          auto t2 = t2id.template goodTrack_as<TTrackTo>();
          auto track0 = getTrackParCov(t0);
          auto track1 = getTrackParCov(t1);
          auto bach = getTrackParCov(t2);

          bool isTrue3bodyVtx = false;
          if (t0.has_mcParticle() && t1.has_mcParticle() && t2.has_mcParticle()) {
            auto t0mc = t0.template mcParticle_as<aod::McParticles>();
            auto t1mc = t1.template mcParticle_as<aod::McParticles>();
            auto t2mc = t2.template mcParticle_as<aod::McParticles>();
            if ((t0mc.pdgCode() == 2212 && t1mc.pdgCode() == -211 && t2mc.pdgCode() == 1000010020) || (t0mc.pdgCode() == 211 && t1mc.pdgCode() == -2212 && t2mc.pdgCode() == -1000010020)) {
              if (t0mc.has_mothers() && t1mc.has_mothers() && t2mc.has_mothers()) {
                for (auto& t0mother : t0mc.template mothers_as<aod::McParticles>()) {
                  for (auto& t1mother : t1mc.template mothers_as<aod::McParticles>()) {
                    for (auto& t2mother : t2mc.template mothers_as<aod::McParticles>()) {
                      if (t0mother.globalIndex() == t1mother.globalIndex() && t0mother.globalIndex() == t2mother.globalIndex() && std::abs(t0mother.pdgCode()) == 1010010030) {
                        isTrue3bodyVtx = true;
                      }
                    }
                  }
                }
              }
            }
          }

          FillVtxCounter(kVtxAll, isTrue3bodyVtx);
          if (!isTrue3bodyVtx && RejectBkgInMC) {
            continue;
          }

          if (bach.getPt() < minbachPt) {
            continue;
          }
          FillVtxCounter(kVtxbachPt, isTrue3bodyVtx);

          int n3bodyVtx = fitter3body.process(track0, track1, bach);
          if (n3bodyVtx == 0) { // discard this pair
            continue;
          }
          FillVtxCounter(kVtxhasSV, isTrue3bodyVtx);

          int cand3B = 0;
          const auto& vertexXYZ = fitter3body.getPCACandidatePos(cand3B);
          // make sure the cascade radius is smaller than that of the vertex
          float dxc = vertexXYZ[0] - dCollision.posX(), dyc = vertexXYZ[1] - dCollision.posY(), dzc = vertexXYZ[2] - dCollision.posZ(), r2vertex = dxc * dxc + dyc * dyc;
          float rvertex = std::sqrt(r2vertex);
          if (std::abs(rv0 - rvertex) > maxRDiff3bodyV0 || rvertex < minRToMeanVertex) {
            continue;
          }
          FillVtxCounter(kVtxRadius, isTrue3bodyVtx);

          // Not involved: bach.minR - rveretx check

          // check: if the process function finish the propagation
          if (!fitter3body.isPropagateTracksToVertexDone() && !fitter3body.propagateTracksToVertex()) {
            continue;
          }

          auto& tr0 = fitter3body.getTrack(0, cand3B);
          auto& tr1 = fitter3body.getTrack(1, cand3B);
          auto& tr2 = fitter3body.getTrack(2, cand3B);
          std::array<float, 3> p0, p1, p2;
          tr0.getPxPyPzGlo(p0);
          tr1.getPxPyPzGlo(p1);
          tr2.getPxPyPzGlo(p2);
          std::array<float, 3> p3B = {p0[0] + p1[0] + p2[0], p0[1] + p1[1] + p2[1], p0[2] + p1[2] + p2[2]};

          float pt2 = p3B[0] * p3B[0] + p3B[1] * p3B[1], p2candidate = pt2 + p3B[2] * p3B[2];
          float pt = std::sqrt(pt2);
          if (pt < minPt3Body) { // pt cut
            continue;
          }
          FillVtxCounter(kVtxPt, isTrue3bodyVtx);

          if (p3B[2] / pt > maxTgl3Body) { // tgLambda cut
            continue;
          }
          FillVtxCounter(kVtxTgLamda, isTrue3bodyVtx);

          float cosPA = (p3B[0] * dxc + p3B[1] * dyc + p3B[2] * dzc) / std::sqrt(p2candidate * (r2vertex + dzc * dzc));
          if (cosPA < minCosPA3body) {
            continue;
          }
          FillVtxCounter(kVtxCosPA, isTrue3bodyVtx);

          if (fitter3body.getChi2AtPCACandidate() > dcavtxdau) {
            continue;
          }
          FillVtxCounter(kVtxDcaDau, isTrue3bodyVtx);

          // Calculate DCA with respect to the collision associated to the V0, not individual tracks
          gpu::gpustd::array<float, 2> dcaInfo;

          auto Track0Par = getTrackPar(t0);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({dCollision.posX(), dCollision.posY(), dCollision.posZ()}, Track0Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
          auto Track0dcaXY = dcaInfo[0];

          auto Track1Par = getTrackPar(t1);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({dCollision.posX(), dCollision.posY(), dCollision.posZ()}, Track1Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
          auto Track1dcaXY = dcaInfo[0];

          auto Track2Par = getTrackPar(t2);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({dCollision.posX(), dCollision.posY(), dCollision.posZ()}, Track2Par, 2.f, fitter3body.getMatCorrType(), &dcaInfo);
          auto Track2dcaXY = dcaInfo[0];

          //  Not involved: H3L DCA Check
          // auto track3B = o2::track::TrackParCov(vertexXYZ, p3B, fitter3body.calcPCACovMatrixFlat(cand3B), t2.sign());
          auto track3B = o2::track::TrackParCov(vertexXYZ, p3B, t2.sign());
          o2::dataformats::DCA dca;
          if (d_UseH3LDCACut && (!track3B.propagateToDCA({{dCollision.posX(), dCollision.posY(), dCollision.posZ()}, {dCollision.covXX(), dCollision.covXY(), dCollision.covYY(), dCollision.covXZ(), dCollision.covYZ(), dCollision.covZZ()}}, fitter3body.getBz(), &dca, 5.) ||
                                 std::abs(dca.getY()) > maxDCAXY3Body || std::abs(dca.getZ()) > maxDCAZ3Body)) {
            continue;
          }
          FillVtxCounter(kVtxDcaH3L, isTrue3bodyVtx);

          vtx3bodydata(
            t0.globalIndex(), t1.globalIndex(), t2.globalIndex(), dCollision.globalIndex(), 0,
            vertexXYZ[0], vertexXYZ[1], vertexXYZ[2],
            p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2],
            fitter3body.getChi2AtPCACandidate(),
            Track0dcaXY, Track1dcaXY, Track2dcaXY);
        }
      }
    }
    fillHistos();
    resetHistos();
  }

  //------------------------------------------------------------------
  // MC virtual lambda check
  template <class TTrackTo, typename TCollisionTable, typename TV0DataTable>
  void VirtualLambdaCheck(TCollisionTable const& dCollision, TV0DataTable const& fullV0s, int bin)
  {
    for (auto& v0 : fullV0s) {
      statisticsRegistry.virtLambdastats[bin]++;
      auto postrack = v0.template posTrack_as<TTrackTo>();
      auto negtrack = v0.template negTrack_as<TTrackTo>();
      if (postrack.has_mcParticle() && negtrack.has_mcParticle()) {
        auto postrackmc = postrack.template mcParticle_as<aod::McParticles>();
        auto negtrackmc = negtrack.template mcParticle_as<aod::McParticles>();

        if ((postrackmc.pdgCode() == 2212 && negtrackmc.pdgCode() == -211) || (postrackmc.pdgCode() == 211 && negtrackmc.pdgCode() == -2212)) {
          if (postrackmc.has_mothers() && negtrackmc.has_mothers()) {
            for (auto& posmother : postrackmc.template mothers_as<aod::McParticles>()) {
              for (auto& negmother : negtrackmc.template mothers_as<aod::McParticles>()) {
                if (posmother.globalIndex() == negmother.globalIndex()) {
                  if (posmother.pdgCode() == 1010010030)
                    statisticsRegistry.virtLambdastats[bin + 1]++;
                  else if (posmother.pdgCode() == -1010010030)
                    statisticsRegistry.virtLambdastats[bin + 2]++;
                }
              }
            }
          }
        }
      }
    }
    fillHistos();
    resetHistos();
  }

  //------------------------------------------------------------------
  // Process Function
  void processData(aod::Collision const& collision, aod::V0GoodPosTracks const& ptracks, aod::V0GoodNegTracks const& ntracks, aod::V0GoodTracks const& goodtracks, FullTracksExtIU const&, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    registry.fill(HIST("hEventCounter"), 0.5);

    DecayFinder<FullTracksExtIU>(collision, ptracks, ntracks, goodtracks);
  }
  PROCESS_SWITCH(hypertriton3bodyFinder, processData, "Produce StoredVtx3BodyDatas with data", false);

  void processCFFilteredData(aod::Collisions const& collisions, aod::CFFilters const& cffilters, aod::V0GoodPosTracks const& Ptracks, aod::V0GoodNegTracks const& Ntracks, aod::V0GoodTracks const& Goodtracks, FullTracksExtIU const&, aod::BCsWithTimestamps const&)
  {
    for (int i{0}; i < collisions.size(); i++) {
      auto collision = collisions.iteratorAt(i);
      auto cffilter = cffilters.iteratorAt(i);
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      registry.fill(HIST("hEventCounter"), 0.5);

      auto ptracks = Ptracks.sliceBy(perCollisionGoodPosTracks, collision.globalIndex());
      auto ntracks = Ntracks.sliceBy(perCollisionGoodNegTracks, collision.globalIndex());
      auto goodtracks = Goodtracks.sliceBy(perCollisionGoodTracks, collision.globalIndex());

      if (!cffilter.hasLD() && UseCFFilter) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 1.5);

      DecayFinder<FullTracksExtIU>(collision, ptracks, ntracks, goodtracks);
    }
  }
  PROCESS_SWITCH(hypertriton3bodyFinder, processCFFilteredData, "Produce StoredVtx3BodyDatas with data using CFtriggers", true);

  void processMC(aod::Collision const& collision, aod::V0GoodPosTracks const& ptracks, aod::V0GoodNegTracks const& ntracks, aod::V0GoodTracks const& goodtracks, aod::McParticles const& mcparticles, FullTracksExtMCIU const&, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    registry.fill(HIST("hEventCounter"), 0.5);

    CheckGoodTracks<FullTracksExtMCIU>(goodtracks, mcparticles);
    DecayFinderMC<FullTracksExtMCIU>(collision, ptracks, ntracks, goodtracks);
  }
  PROCESS_SWITCH(hypertriton3bodyFinder, processMC, "Produce StoredVtx3BodyDatas with MC", false);

  void processCFFilteredMC(aod::Collisions const& collisions, aod::CFFilters const& cffilters, aod::V0GoodPosTracks const& Ptracks, aod::V0GoodNegTracks const& Ntracks, aod::V0GoodTracks const& Goodtracks, aod::V0s const& V0s, aod::V0Datas const& fullV0s, aod::McParticles const& mcparticles, FullTracksExtMCIU const&, aod::BCsWithTimestamps const&)
  {
    for (int i{0}; i < collisions.size(); i++) {
      auto collision = collisions.iteratorAt(i);
      auto cffilter = cffilters.iteratorAt(i);
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      registry.fill(HIST("hEventCounter"), 0.5);

      auto goodtracks = Goodtracks.sliceBy(perCollisionGoodTracks, collision.globalIndex());
      CheckGoodTracks<FullTracksExtMCIU>(goodtracks, mcparticles);
      auto v0s = V0s.sliceBy(perCollisionV0s, collision.globalIndex());
      auto fullv0s = fullV0s.sliceBy(perCollisionV0Datas, collision.globalIndex());
      VirtualLambdaCheck<FullTracksExtMCIU>(collision, v0s, 0);
      VirtualLambdaCheck<FullTracksExtMCIU>(collision, fullv0s, 3);

      if (!cffilter.hasLD() && UseCFFilter) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 1.5);

      auto ptracks = Ptracks.sliceBy(perCollisionGoodPosTracks, collision.globalIndex());
      auto ntracks = Ntracks.sliceBy(perCollisionGoodNegTracks, collision.globalIndex());

      VirtualLambdaCheck<FullTracksExtMCIU>(collision, v0s, 6);
      VirtualLambdaCheck<FullTracksExtMCIU>(collision, fullv0s, 9);
      DecayFinderMC<FullTracksExtMCIU>(collision, ptracks, ntracks, goodtracks);
    }
  }
  PROCESS_SWITCH(hypertriton3bodyFinder, processCFFilteredMC, "Produce StoredVtx3BodyDatas with MC using CFtriggers", false);
};

struct hypertriton3bodyLabelBuilder {

  Produces<aod::McVtx3BodyLabels> vtxlabels;

  // for bookkeeping purposes: how many V0s come from same mother etc
  HistogramRegistry registry{
    "registry",
    {
      {"hLabelCounter", "hLabelCounter", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
      {"hHypertritonCounter", "hHypertritonCounter", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
      {"hPIDCounter", "hPIDCounter", {HistType::kTH1F, {{6, 0.0f, 6.0f}}}},
      {"hHypertritonMCPt", "hHypertritonMCPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hAntiHypertritonMCPt", "hAntiHypertritonMCPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hHypertritonMCMass", "hHypertritonMCMass", {HistType::kTH1F, {{40, 2.95f, 3.05f}}}},
      {"hAntiHypertritonMCMass", "hAntiHypertritonMCMass", {HistType::kTH1F, {{40, 2.95f, 3.05f}}}},
      {"h3dTotalTrueHypertriton", "h3dTotalTrueHypertriton", {HistType::kTH3F, {{50, 0, 50, "ct(cm)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
    },
  };

  void init(InitContext const&)
  {
    registry.get<TH1>(HIST("hLabelCounter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hLabelCounter"))->GetXaxis()->SetBinLabel(2, "Same MotherParticle");
    registry.get<TH1>(HIST("hLabelCounter"))->GetXaxis()->SetBinLabel(3, "True H3L");
    registry.get<TH1>(HIST("hHypertritonCounter"))->GetXaxis()->SetBinLabel(1, "H3L");
    registry.get<TH1>(HIST("hHypertritonCounter"))->GetXaxis()->SetBinLabel(2, "H3L daughters pass PID");
    registry.get<TH1>(HIST("hHypertritonCounter"))->GetXaxis()->SetBinLabel(3, "#bar{H3L}");
    registry.get<TH1>(HIST("hHypertritonCounter"))->GetXaxis()->SetBinLabel(4, "#bar{H3L} daughters pass PID");
    registry.get<TH1>(HIST("hPIDCounter"))->GetXaxis()->SetBinLabel(1, "H3L Proton PID > 5");
    registry.get<TH1>(HIST("hPIDCounter"))->GetXaxis()->SetBinLabel(2, "H3L Pion PID > 5");
    registry.get<TH1>(HIST("hPIDCounter"))->GetXaxis()->SetBinLabel(3, "H3L Deuteron PID > 5");
    registry.get<TH1>(HIST("hPIDCounter"))->GetXaxis()->SetBinLabel(4, "#bar{H3L} Proton PID > 5");
    registry.get<TH1>(HIST("hPIDCounter"))->GetXaxis()->SetBinLabel(5, "#bar{H3L} Pion PID > 5");
    registry.get<TH1>(HIST("hPIDCounter"))->GetXaxis()->SetBinLabel(6, "#bar{H3L} Deuteron PID > 5");
  }

  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};

  void processDoNotBuildLabels(aod::Collisions::iterator const& collision)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(hypertriton3bodyLabelBuilder, processDoNotBuildLabels, "Do not produce MC label tables", true);

  void processBuildLabels(aod::Vtx3BodyDatas const& vtx3bodydatas, LabeledTracks const&, aod::McParticles const& particlesMC)
  {
    std::vector<int> lIndices;
    lIndices.reserve(vtx3bodydatas.size());
    for (int ii = 0; ii < vtx3bodydatas.size(); ii++) {
      lIndices[ii] = -1;
    }

    for (auto& vtx3body : vtx3bodydatas) {

      int lLabel = -1;
      int lPDG = -1;
      float lPt = -1;
      double MClifetime = -1;
      bool is3bodyDecay = false;
      int lGlobalIndex = -1;

      auto lTrack0 = vtx3body.track0_as<LabeledTracks>();
      auto lTrack1 = vtx3body.track1_as<LabeledTracks>();
      auto lTrack2 = vtx3body.track2_as<LabeledTracks>();
      registry.fill(HIST("hLabelCounter"), 0.5);

      // Association check
      // There might be smarter ways of doing this in the future
      if (!lTrack0.has_mcParticle() || !lTrack1.has_mcParticle() || !lTrack2.has_mcParticle()) {
        vtxlabels(-1);
        continue;
      }
      auto lMCTrack0 = lTrack0.mcParticle_as<aod::McParticles>();
      auto lMCTrack1 = lTrack1.mcParticle_as<aod::McParticles>();
      auto lMCTrack2 = lTrack2.mcParticle_as<aod::McParticles>();
      if (!lMCTrack0.has_mothers() || !lMCTrack1.has_mothers() || !lMCTrack2.has_mothers()) {
        vtxlabels(-1);
        continue;
      }

      for (auto& lMother0 : lMCTrack0.mothers_as<aod::McParticles>()) {
        for (auto& lMother1 : lMCTrack1.mothers_as<aod::McParticles>()) {
          for (auto& lMother2 : lMCTrack2.mothers_as<aod::McParticles>()) {
            if (lMother0.globalIndex() == lMother1.globalIndex() && lMother0.globalIndex() == lMother2.globalIndex()) {
              lGlobalIndex = lMother1.globalIndex();
              lPt = lMother1.pt();
              lPDG = lMother1.pdgCode();
              MClifetime = RecoDecay::sqrtSumOfSquares(lMCTrack2.vx() - lMother2.vx(), lMCTrack2.vy() - lMother2.vy(), lMCTrack2.vz() - lMother2.vz()) * o2::constants::physics::MassHyperTriton / lMother2.p();
              is3bodyDecay = true; // vtxs with the same mother
            }
          }
        }
      } // end association check
      if (!is3bodyDecay) {
        vtxlabels(-1);
        continue;
      }
      registry.fill(HIST("hLabelCounter"), 1.5);

      // Intended for cross-checks only
      // N.B. no rapidity cut!
      if (lPDG == 1010010030 && lMCTrack0.pdgCode() == 2212 && lMCTrack1.pdgCode() == -211 && lMCTrack2.pdgCode() == 1000010020) {
        lLabel = lGlobalIndex;
        double hypertritonMCMass = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
        registry.fill(HIST("hLabelCounter"), 2.5);
        registry.fill(HIST("hHypertritonCounter"), 0.5);
        registry.fill(HIST("hHypertritonMCPt"), lPt);
        registry.fill(HIST("hHypertritonMCMass"), hypertritonMCMass);
        registry.fill(HIST("h3dTotalTrueHypertriton"), MClifetime, lPt, hypertritonMCMass);
        if (TMath::Abs(lTrack0.tpcNSigmaPr()) > TpcPidNsigmaCut) {
          registry.fill(HIST("hPIDCounter"), 0.5);
        }
        if (TMath::Abs(lTrack1.tpcNSigmaPi()) > TpcPidNsigmaCut) {
          registry.fill(HIST("hPIDCounter"), 1.5);
        }
        if (TMath::Abs(lTrack2.tpcNSigmaDe()) > TpcPidNsigmaCut) {
          registry.fill(HIST("hPIDCounter"), 2.5);
        }
        if (TMath::Abs(lTrack0.tpcNSigmaPr()) < TpcPidNsigmaCut && TMath::Abs(lTrack1.tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(lTrack2.tpcNSigmaDe()) < TpcPidNsigmaCut) {
          registry.fill(HIST("hHypertritonCounter"), 1.5);
        }
      }
      if (lPDG == -1010010030 && lMCTrack0.pdgCode() == 211 && lMCTrack1.pdgCode() == -2212 && lMCTrack2.pdgCode() == -1000010020) {
        lLabel = lGlobalIndex;
        double antiHypertritonMCMass = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron});
        registry.fill(HIST("hLabelCounter"), 2.5);
        registry.fill(HIST("hHypertritonCounter"), 2.5);
        registry.fill(HIST("hAntiHypertritonMCPt"), lPt);
        registry.fill(HIST("hAntiHypertritonMCMass"), antiHypertritonMCMass);
        registry.fill(HIST("h3dTotalTrueHypertriton"), MClifetime, lPt, antiHypertritonMCMass);
        if (TMath::Abs(lTrack0.tpcNSigmaPi()) > TpcPidNsigmaCut) {
          registry.fill(HIST("hPIDCounter"), 4.5);
        }
        if (TMath::Abs(lTrack1.tpcNSigmaPr()) > TpcPidNsigmaCut) {
          registry.fill(HIST("hPIDCounter"), 3.5);
        }
        if (TMath::Abs(lTrack2.tpcNSigmaDe()) > TpcPidNsigmaCut) {
          registry.fill(HIST("hPIDCounter"), 5.5);
        }
        if (TMath::Abs(lTrack0.tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(lTrack1.tpcNSigmaPr()) < TpcPidNsigmaCut && TMath::Abs(lTrack2.tpcNSigmaDe()) < TpcPidNsigmaCut) {
          registry.fill(HIST("hHypertritonCounter"), 3.5);
        }
      }

      // Construct label table, only true hypertriton and true daughters with a specified order is labeled
      // for matter: track0->p, track1->pi, track2->d
      // for antimatter: track0->pi, track1->p, track2->d
      vtxlabels(lLabel);
    }
  }
  PROCESS_SWITCH(hypertriton3bodyLabelBuilder, processBuildLabels, "Produce MC label tables", false);
};

struct hypertriton3bodyInitializer {
  Spawns<aod::Vtx3BodyDatas> vtx3bodydatas;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<trackprefilter>(cfgc),
    adaptAnalysisTask<hypertriton3bodyFinder>(cfgc),
    adaptAnalysisTask<hypertriton3bodyLabelBuilder>(cfgc),
    adaptAnalysisTask<hypertriton3bodyInitializer>(cfgc),
  };
}
