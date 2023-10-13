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

struct hypertriton3bodyFinder {

  Produces<aod::StoredVtx3BodyDatas> vtx3bodydata;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Configurables
  Configurable<bool> UseCFFilterInMC{"UseCFFilterInMC", false, "Use CFFilter with MC check"};

  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};

  // Selection criteria
  Configurable<float> minR2ToMeanVertex = {"minR2ToMeanVertex", 0.5 * 0.5, ""}; ///< min radial distance of V0 from beam line (mean vertex) //Q: Does this cut need to be removed?

  Configurable<float> causalityRTolerance = {"causalityRTolerance", 1., ""}; ///< V0 radius cannot exceed its contributors minR by more than this value
  Configurable<float> maxV0ToProngsRDiff = {"maxV0ToProngsRDiff", 50., ""};  ///< V0 radius cannot be lower than this ammount wrt minR of contributors
  Configurable<float> minPt2V0 = {"minPt2V0", 0.5 * 0.5, ""};                ///< v0 minimum pT
  Configurable<float> maxTgl2V0 = {"maxTgl2V0", 2. * 2., ""};                ///< maximum tgLambda of V0
  Configurable<float> maxDCAXY2ToMeanVertex3bodyV0 = {"maxDCAXY2ToMeanVertex3bodyV0", 2 * 2, ""};
  Configurable<float> minCosPAXYMeanVertex3bodyV0 = {"minCosPAXYMeanVertex3bodyV0", 0.9, ""}; ///< min cos of PA to beam line (mean vertex) in tr. plane for 3body V0 cand.
  Configurable<float> minCosPA3bodyV0 = {"minCosPA3bodyV0", 0.8, ""};                         // min cos of PA to PV for 3body V0

  // for 3 body reconstructed Vertex
  Configurable<float> maxRDiff3bodyV0 = {"maxRDiffV03bodyV0", 3, ""}; ///< Maximum difference between V0 and 3body radii
  Configurable<float> minPt23Body = {"minPt23Body", 0.01 * 0.01, ""}; // minimum pT of 3body Vertex
  Configurable<float> maxTgl23Body = {"maxTgl23Body", 2. * 2., ""};   // maximum tgLambda of 3body Vertex
  Configurable<float> minCosPA3body = {"minCosPA3body", 0.8, ""};     // min cos of PA to PV for 3body Vertex

  // for DCA
  // Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};

  // for track cut in SVertexer, Can we use it in the production of goodtrack table?
  // float maxDCAXY3Body = 0.3; // max DCA of 3 body decay to PV in XY?
  // float maxDCAZ3Body = 0.3;  // max DCA of 3 body decay to PV in Z

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

  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hDauTrackCounter", "hDauTrackCounter", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
      {"hV0Counter", "hV0Counter", {HistType::kTH1F, {{10, 0.0f, 10.0f}}}},
      {"hVtx3BodyCounter", "hVtx3BodyCounter", {HistType::kTH1F, {{7, 0.0f, 7.0f}}}},
      {"hTrueV0Counter", "hTrueV0Counter", {HistType::kTH1F, {{10, 0.0f, 10.0f}}}},
      {"hTrueVtx3BodyCounter", "hTrueVtx3BodyCounter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},
    },
  };

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::vertexing::DCAFitterN<2> fitter;
  o2::vertexing::DCAFitterN<3> fitter3body;

  void init(InitContext& context)
  {
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  // could be changed later
    maxStep = 2.00f; // could be changed later

    registry.get<TH1>(HIST("hDauTrackCounter"))->GetXaxis()->SetBinLabel(1, "Proton");
    registry.get<TH1>(HIST("hDauTrackCounter"))->GetXaxis()->SetBinLabel(2, "Pion");
    registry.get<TH1>(HIST("hDauTrackCounter"))->GetXaxis()->SetBinLabel(3, "Deuteron");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(2, "hasSV");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(3, "V0R");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(4, "TrackR-V0R(off)");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(5, "IfPropragated");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(6, "V0Pt");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(7, "tgLambda");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(8, "V0Mass");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(9, "dcaXY&CosPAXY");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(10, "CosPA");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(2, "hasSV");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(3, "V0R");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(4, "TrackRi-V0R(off)");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(5, "IfPropragated");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(6, "V0Pt");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(7, "tgLambda");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(8, "V0Mass");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(9, "dcaXY&CosPAXY");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(10, "CosPA");

    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(2, "bachPt");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(3, "hasSV");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(4, "VtxR");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(5, "BachR-VtxR(off)");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(6, "IfPropragated");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(7, "VtxPt");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(8, "tgLambda");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(9, "CosPA");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(2, "bachPt");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(3, "hasSV");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(4, "VtxR");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(5, "BachR-vtxR(off)");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(6, "IfPropragated");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(7, "VtxPt");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(8, "tgLambda");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(9, "CosPA");

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
      bool isHypertriton = false;
      for (auto& mothertrack : mcgoodtrack.template mothers_as<aod::McParticles>()) {
        if (std::abs(mothertrack.pdgCode()) != 1010010030) {
          continue;
        }
        bool haveProton = false, havePion = false, haveDeuteron = false;
        bool haveAntiProton = false, haveAntiPion = false, haveAntiDeuteron = false;
        for (auto& mcparticleDaughter : mothertrack.template daughters_as<aod::McParticles>()) {
          if (mcparticleDaughter.pdgCode() == 2212) {
            haveProton = true;
          }
          if (mcparticleDaughter.pdgCode() == -2212) {
            haveAntiProton = true;
          }
          if (mcparticleDaughter.pdgCode() == 211) {
            havePion = true;
          }
          if (mcparticleDaughter.pdgCode() == -211) {
            haveAntiPion = true;
          }
          if (mcparticleDaughter.pdgCode() == 1000010020) {
            haveDeuteron = true;
          }
          if (mcparticleDaughter.pdgCode() == -1000010020) {
            haveAntiDeuteron = true;
          }
        }
        if ((haveProton && haveAntiPion && haveDeuteron && mothertrack.pdgCode() == 1010010030) || (haveAntiProton && havePion && haveAntiDeuteron && mothertrack.pdgCode() == -1010010030)) {
          isHypertriton = true;
          break;
        }
      }
      if (isHypertriton && std::abs(mcgoodtrack.pdgCode()) == 2212) {
        registry.fill(HIST("hDauTrackCounter"), 0.5);
      }
      if (isHypertriton && std::abs(mcgoodtrack.pdgCode()) == 211) {
        registry.fill(HIST("hDauTrackCounter"), 1.5);
      }
      if (isHypertriton && std::abs(mcgoodtrack.pdgCode()) == 1000010020) {
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
      registry.fill(HIST("hV0Counter"), 0.5);
      auto t0 = t0id.template goodTrack_as<TTrackTo>();
      auto Track0 = getTrackParCov(t0);

      for (auto& t1id : dNtracks) {
        auto t1 = t1id.template goodTrack_as<TTrackTo>();
        auto Track1 = getTrackParCov(t1);
        int nCand = fitter.process(Track0, Track1);
        if (nCand == 0) {
          continue;
        }
        registry.fill(HIST("hV0Counter"), 1.5);
        // validate V0 radial position
        // First check closeness to the beam-line as same as SVertexer
        const auto& v0XYZ = fitter.getPCACandidate();
        float dxv0 = v0XYZ[0] - mMeanVertex.getX(), dyv0 = v0XYZ[1] - mMeanVertex.getY(), r2v0 = dxv0 * dxv0 + dyv0 * dyv0;
        if (r2v0 < minR2ToMeanVertex) {
          continue;
        }
        registry.fill(HIST("hV0Counter"), 2.5);

        float rv0 = std::sqrt(r2v0);
        // Not involved: Get minR with same way in SVertexer
        // float Track0minR =  RecoDecay::sqrtSumOfSquares(t0.x(),  t0.y()), Track1minR =  RecoDecay::sqrtSumOfSquares(t1.x(),  t1.y());
        // float drv0P = rv0 - Track0minR, drv0N = rv0 - Track1minR;
        /*if (drv0P > causalityRTolerance || drv0P < -maxV0ToProngsRDiff ||
          drv0N > causalityRTolerance || drv0N < -maxV0ToProngsRDiff) {
          continue;
          }*/
        registry.fill(HIST("hV0Counter"), 3.5);

        // the process function should finish the propagation
        if (!fitter.isPropagateTracksToVertexDone() && !fitter.propagateTracksToVertex()) {
          continue;
        }
        registry.fill(HIST("hV0Counter"), 4.5);

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
        if (pt2V0 < minPt2V0) { // pt cut
          continue;
        }
        registry.fill(HIST("hV0Counter"), 5.5);

        if (pV0[2] * pV0[2] / pt2V0 > maxTgl2V0) { // tgLambda cut
          continue;
        }
        registry.fill(HIST("hV0Counter"), 6.5);
        float p2V0 = pt2V0 + pV0[2] * pV0[2], ptV0 = std::sqrt(pt2V0);
        // apply mass selections
        float massV0LambdaHyp = RecoDecay::m(array{array{pP[0], pP[1], pP[2]}, array{pN[0], pN[1], pN[2]}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
        float massV0AntiLambdaHyp = RecoDecay::m(array{array{pP[0], pP[1], pP[2]}, array{pN[0], pN[1], pN[2]}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
        float massMargin = 20 * (0.001 * (1. + 0.5 * ptV0)) + 0.07;
        if (massV0LambdaHyp - o2::constants::physics::MassLambda > massMargin && massV0AntiLambdaHyp - o2::constants::physics::MassLambda > massMargin) {
          continue;
        }
        registry.fill(HIST("hV0Counter"), 7.5);

        float dcaX = dxv0 - pV0[0] * tDCAXY, dcaY = dyv0 - pV0[1] * tDCAXY, dca2 = dcaX * dcaX + dcaY * dcaY;
        float cosPAXY = prodXYv0 / std::sqrt(r2v0 * pt2V0);
        if (dca2 > maxDCAXY2ToMeanVertex3bodyV0 || cosPAXY < minCosPAXYMeanVertex3bodyV0) {
          continue;
        }
        registry.fill(HIST("hV0Counter"), 8.5);

        float dx = v0XYZ[0] - dCollision.posX(), dy = v0XYZ[1] - dCollision.posY(), dz = v0XYZ[2] - dCollision.posZ(), prodXYZv0 = dx * pV0[0] + dy * pV0[1] + dz * pV0[2];
        float cosPA = prodXYZv0 / std::sqrt((dx * dx + dy * dy + dz * dz) * p2V0);
        if (cosPA < minCosPA3bodyV0) {
          continue;
        }
        registry.fill(HIST("hV0Counter"), 9.5);

        for (auto& t2id : dGoodtracks) {
          if (t2id.globalIndex() == t0id.globalIndex()) {
            continue; // skip the track used by V0
          }
          registry.fill(HIST("hVtx3BodyCounter"), 0.5);

          auto t2 = t2id.template goodTrack_as<TTrackTo>();
          auto track0 = getTrackParCov(t0);
          auto track1 = getTrackParCov(t1);
          auto bach = getTrackParCov(t2);

          if (bach.getPt() < 0.6) {
            continue;
          }
          registry.fill(HIST("hVtx3BodyCounter"), 1.5);

          int n3bodyVtx = fitter3body.process(track0, track1, bach);
          if (n3bodyVtx == 0) { // discard this pair
            continue;
          }
          registry.fill(HIST("hVtx3BodyCounter"), 2.5);

          int cand3B = 0;
          const auto& vertexXYZ = fitter3body.getPCACandidatePos(cand3B);
          // make sure the cascade radius is smaller than that of the vertex
          float dxc = vertexXYZ[0] - dCollision.posX(), dyc = vertexXYZ[1] - dCollision.posY(), dzc = vertexXYZ[2] - dCollision.posZ(), r2vertex = dxc * dxc + dyc * dyc;
          if (std::abs(rv0 - std::sqrt(r2vertex)) > maxRDiff3bodyV0 || r2vertex < minR2ToMeanVertex) {
            continue;
          }
          registry.fill(HIST("hVtx3BodyCounter"), 3.5);

          // Not involved: bach.minR - rveretx check
          registry.fill(HIST("hVtx3BodyCounter"), 4.5);

          if (!fitter3body.isPropagateTracksToVertexDone() && !fitter3body.propagateTracksToVertex()) {
            continue;
          }
          registry.fill(HIST("hVtx3BodyCounter"), 5.5);

          auto& tr0 = fitter3body.getTrack(0, cand3B);
          auto& tr1 = fitter3body.getTrack(1, cand3B);
          auto& tr2 = fitter3body.getTrack(2, cand3B);
          std::array<float, 3> p0, p1, p2;
          tr0.getPxPyPzGlo(p0);
          tr1.getPxPyPzGlo(p1);
          tr2.getPxPyPzGlo(p2);
          std::array<float, 3> p3B = {p0[0] + p1[0] + p2[0], p0[1] + p1[1] + p2[1], p0[2] + p1[2] + p2[2]};

          float pt2 = p3B[0] * p3B[0] + p3B[1] * p3B[1], p2candidate = pt2 + p3B[2] * p3B[2];
          if (pt2 < minPt23Body) { // pt cut
            continue;
          }
          registry.fill(HIST("hVtx3BodyCounter"), 6.5);

          if (p3B[2] * p3B[2] / pt2 > maxTgl23Body) { // tgLambda cut
            continue;
          }
          registry.fill(HIST("hVtx3BodyCounter"), 7.5);

          float cosPA = (p3B[0] * dxc + p3B[1] * dyc + p3B[2] * dzc) / std::sqrt(p2candidate * (r2vertex + dzc * dzc));
          if (cosPA < minCosPA3body) {
            continue;
          }
          registry.fill(HIST("hVtx3BodyCounter"), 8.5);

          //  Not involved: H3L DCA Check
          vtx3bodydata(
            t0.globalIndex(), t1.globalIndex(), t2.globalIndex(), dCollision.globalIndex(), 0,
            vertexXYZ[0], vertexXYZ[1], vertexXYZ[2],
            p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2],
            fitter3body.getChi2AtPCACandidate(),
            t0.dcaXY(), t1.dcaXY(), t2.dcaXY());
        }
      }
    }
  }

  void processData(aod::Collision const& collision, aod::V0GoodPosTracks const& ptracks, aod::V0GoodNegTracks const& ntracks, aod::V0GoodTracks const& goodtracks, FullTracksExtIU const&, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    registry.fill(HIST("hEventCounter"), 0.5);
    registry.fill(HIST("hEventCounter"), 1.5);

    DecayFinder<FullTracksExtIU>(collision, ptracks, ntracks, goodtracks);
  }
  PROCESS_SWITCH(hypertriton3bodyFinder, processData, "Produce StoredVtx3BodyDatas without MC check", false);

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

      if (!cffilter.hasLD()) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 1.5);

      DecayFinder<FullTracksExtIU>(collision, ptracks, ntracks, goodtracks);
    }
  }
  PROCESS_SWITCH(hypertriton3bodyFinder, processCFFilteredData, "Produce StoredVtx3BodyDatas without MC check under cffilters", true);

  void processMC(aod::Collisions const& collisions, aod::CFFilters const& cffilters, aod::V0GoodPosTracks const& Ptracks, aod::V0GoodNegTracks const& Ntracks, aod::V0GoodTracks const& Goodtracks, aod::McParticles const& mcparticles, FullTracksExtMCIU const&, aod::BCsWithTimestamps const&)
  {
    for (int i{0}; i < collisions.size(); i++) {
      auto collision = collisions.iteratorAt(i);
      auto cffilter = cffilters.iteratorAt(i);
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      registry.fill(HIST("hEventCounter"), 0.5);

      auto goodtracks = Goodtracks.sliceBy(perCollisionGoodTracks, collision.globalIndex());
      CheckGoodTracks<FullTracksExtMCIU>(goodtracks, mcparticles);

      if (!cffilter.hasLD() && UseCFFilterInMC) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 1.5);

      auto ptracks = Ptracks.sliceBy(perCollisionGoodPosTracks, collision.globalIndex());
      auto ntracks = Ntracks.sliceBy(perCollisionGoodNegTracks, collision.globalIndex());

      for (auto& t0id : ptracks) { // FIXME: turn into combination(...)
        auto t0 = t0id.goodTrack_as<FullTracksExtMCIU>();
        auto Track0 = getTrackParCov(t0);
        for (auto& t1id : ntracks) {

          if (t0id.collisionId() != t1id.collisionId()) {
            continue;
          }
          auto t1 = t1id.goodTrack_as<FullTracksExtMCIU>();
          auto Track1 = getTrackParCov(t1);

          bool isTrue3bodyV0 = false;
          if (t0.has_mcParticle() && t1.has_mcParticle()) {
            auto t0mc = t0.mcParticle_as<aod::McParticles>();
            auto t1mc = t1.mcParticle_as<aod::McParticles>();
            if ((t0mc.pdgCode() == 2212 && t1mc.pdgCode() == -211) || (t0mc.pdgCode() == 211 && t1mc.pdgCode() == -2212)) {
              if (t0mc.has_mothers() && t1mc.has_mothers()) {
                for (auto& t0mother : t0mc.mothers_as<aod::McParticles>()) {
                  for (auto& t1mother : t1mc.mothers_as<aod::McParticles>()) {
                    if (t0mother.globalIndex() == t1mother.globalIndex() && std::abs(t0mother.pdgCode()) == 1010010030) {
                      isTrue3bodyV0 = true;
                    }
                  }
                }
              }
            }
          }

          if (!isTrue3bodyV0) {
            continue;
          } // save time

          registry.fill(HIST("hV0Counter"), 0.5);
          if (isTrue3bodyV0) {
            registry.fill(HIST("hTrueV0Counter"), 0.5);
          }

          int nCand = fitter.process(Track0, Track1);
          if (nCand == 0) {
            continue;
          }
          registry.fill(HIST("hV0Counter"), 1.5);
          if (isTrue3bodyV0) {
            registry.fill(HIST("hTrueV0Counter"), 1.5);
          }
          // validate V0 radial position
          // First check closeness to the beam-line as same as SVertexer
          const auto& v0XYZ = fitter.getPCACandidate();
          float dxv0 = v0XYZ[0] - mMeanVertex.getX(), dyv0 = v0XYZ[1] - mMeanVertex.getY(), r2v0 = dxv0 * dxv0 + dyv0 * dyv0;
          if (r2v0 < minR2ToMeanVertex) {
            continue;
          }
          registry.fill(HIST("hV0Counter"), 2.5);
          if (isTrue3bodyV0) {
            registry.fill(HIST("hTrueV0Counter"), 2.5);
          }

          float rv0 = std::sqrt(r2v0);
          // Not involved: Get minR with same way in SVertexer
          // float Track0minR =  RecoDecay::sqrtSumOfSquares(t0.x(),  t0.y()), Track1minR =  RecoDecay::sqrtSumOfSquares(t1.x(),  t1.y());
          // float drv0P = rv0 - Track0minR, drv0N = rv0 - Track1minR;
          // if (drv0P > causalityRTolerance || drv0P < -maxV0ToProngsRDiff ||
          //     drv0N > causalityRTolerance || drv0N < -maxV0ToProngsRDiff) {
          //     continue;
          // }
          registry.fill(HIST("hV0Counter"), 3.5);
          if (isTrue3bodyV0) {
            registry.fill(HIST("hTrueV0Counter"), 3.5);
          }

          // the process function should finish the propagation
          if (!fitter.isPropagateTracksToVertexDone() && !fitter.propagateTracksToVertex()) {
            continue;
          }
          registry.fill(HIST("hV0Counter"), 4.5);
          if (isTrue3bodyV0) {
            registry.fill(HIST("hTrueV0Counter"), 4.5);
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
          if (pt2V0 < minPt2V0) { // pt cut
            continue;
          }
          registry.fill(HIST("hV0Counter"), 5.5);
          if (isTrue3bodyV0) {
            registry.fill(HIST("hTrueV0Counter"), 5.5);
          }

          if (pV0[2] * pV0[2] / pt2V0 > maxTgl2V0) { // tgLambda cut
            continue;
          }
          registry.fill(HIST("hV0Counter"), 6.5);
          if (isTrue3bodyV0) {
            registry.fill(HIST("hTrueV0Counter"), 6.5);
          }
          float p2V0 = pt2V0 + pV0[2] * pV0[2], ptV0 = std::sqrt(pt2V0);
          // apply mass selections
          float massV0LambdaHyp = RecoDecay::m(array{array{pP[0], pP[1], pP[2]}, array{pN[0], pN[1], pN[2]}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
          float massV0AntiLambdaHyp = RecoDecay::m(array{array{pP[0], pP[1], pP[2]}, array{pN[0], pN[1], pN[2]}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
          float massMargin = 20 * (0.001 * (1. + 0.5 * ptV0)) + 0.07;
          if (massV0LambdaHyp - o2::constants::physics::MassLambda > massMargin && massV0AntiLambdaHyp - o2::constants::physics::MassLambda > massMargin) {
            continue;
          }
          registry.fill(HIST("hV0Counter"), 7.5);
          if (isTrue3bodyV0) {
            registry.fill(HIST("hTrueV0Counter"), 7.5);
          }

          float dcaX = dxv0 - pV0[0] * tDCAXY, dcaY = dyv0 - pV0[1] * tDCAXY, dca2 = dcaX * dcaX + dcaY * dcaY;
          float cosPAXY = prodXYv0 / std::sqrt(r2v0 * pt2V0);
          if (dca2 > maxDCAXY2ToMeanVertex3bodyV0 || cosPAXY < minCosPAXYMeanVertex3bodyV0) {
            continue;
          }
          registry.fill(HIST("hV0Counter"), 8.5);
          if (isTrue3bodyV0) {
            registry.fill(HIST("hTrueV0Counter"), 8.5);
          }

          float dx = v0XYZ[0] - collision.posX(), dy = v0XYZ[1] - collision.posY(), dz = v0XYZ[2] - collision.posZ(), prodXYZv0 = dx * pV0[0] + dy * pV0[1] + dz * pV0[2];
          float cosPA = prodXYZv0 / std::sqrt((dx * dx + dy * dy + dz * dz) * p2V0);
          if (cosPA < minCosPA3bodyV0) {
            continue;
          }
          registry.fill(HIST("hV0Counter"), 9.5);
          if (isTrue3bodyV0) {
            registry.fill(HIST("hTrueV0Counter"), 9.5);
          }

          for (auto& t2id : goodtracks) {
            if (t2id.globalIndex() == t0id.globalIndex()) {
              continue; // skip the track used by V0
            }
            auto t2 = t2id.goodTrack_as<FullTracksExtMCIU>();
            auto track0 = getTrackParCov(t0);
            auto track1 = getTrackParCov(t1);
            auto bach = getTrackParCov(t2);

            bool isTrue3bodyVtx = false;
            if (t0.has_mcParticle() && t1.has_mcParticle() && t2.has_mcParticle()) {
              auto t0mc = t0.mcParticle_as<aod::McParticles>();
              auto t1mc = t1.mcParticle_as<aod::McParticles>();
              auto t2mc = t2.mcParticle_as<aod::McParticles>();
              if ((t0mc.pdgCode() == 2212 && t1mc.pdgCode() == -211 && t2mc.pdgCode() == 1000010020) || (t0mc.pdgCode() == 211 && t1mc.pdgCode() == -2212 && t2mc.pdgCode() == -1000010020)) {
                if (t0mc.has_mothers() && t1mc.has_mothers() && t2mc.has_mothers()) {
                  for (auto& t0mother : t0mc.mothers_as<aod::McParticles>()) {
                    for (auto& t1mother : t1mc.mothers_as<aod::McParticles>()) {
                      for (auto& t2mother : t2mc.mothers_as<aod::McParticles>()) {
                        if (t0mother.globalIndex() == t1mother.globalIndex() && t0mother.globalIndex() == t2mother.globalIndex() && std::abs(t0mother.pdgCode()) == 1010010030) {
                          isTrue3bodyVtx = true;
                        }
                      }
                    }
                  }
                }
              }
            }

            // if (!isTrue3bodyVtx){ continue;} // save time
            registry.fill(HIST("hVtx3BodyCounter"), 0.5);
            if (isTrue3bodyVtx) {
              registry.fill(HIST("hTrueVtx3BodyCounter"), 0.5);
            }

            if (bach.getPt() < 0.6) {
              continue;
            }
            registry.fill(HIST("hVtx3BodyCounter"), 1.5);
            if (isTrue3bodyVtx) {
              registry.fill(HIST("hTrueVtx3BodyCounter"), 1.5);
            }

            int n3bodyVtx = fitter3body.process(track0, track1, bach);
            if (n3bodyVtx == 0) { // discard this pair
              continue;
            }
            registry.fill(HIST("hVtx3BodyCounter"), 2.5);
            if (isTrue3bodyVtx) {
              registry.fill(HIST("hTrueVtx3BodyCounter"), 2.5);
            }

            int cand3B = 0;
            const auto& vertexXYZ = fitter3body.getPCACandidatePos(cand3B);
            // make sure the cascade radius is smaller than that of the vertex
            float dxc = vertexXYZ[0] - collision.posX(), dyc = vertexXYZ[1] - collision.posY(), dzc = vertexXYZ[2] - collision.posZ(), r2vertex = dxc * dxc + dyc * dyc;
            if (std::abs(rv0 - std::sqrt(r2vertex)) > maxRDiff3bodyV0 || r2vertex < minR2ToMeanVertex) {
              continue;
            }
            registry.fill(HIST("hVtx3BodyCounter"), 3.5);
            if (isTrue3bodyVtx) {
              registry.fill(HIST("hTrueVtx3BodyCounter"), 3.5);
            }

            // Not involved: bach.minR - rveretx check
            registry.fill(HIST("hVtx3BodyCounter"), 4.5);
            if (isTrue3bodyVtx) {
              registry.fill(HIST("hTrueVtx3BodyCounter"), 4.5);
            }

            if (!fitter3body.isPropagateTracksToVertexDone() && !fitter3body.propagateTracksToVertex()) {
              continue;
            }
            registry.fill(HIST("hVtx3BodyCounter"), 5.5);
            if (isTrue3bodyVtx) {
              registry.fill(HIST("hTrueVtx3BodyCounter"), 5.5);
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
            if (pt2 < minPt23Body) { // pt cut
              continue;
            }
            registry.fill(HIST("hVtx3BodyCounter"), 6.5);
            if (isTrue3bodyVtx) {
              registry.fill(HIST("hTrueVtx3BodyCounter"), 6.5);
            }

            if (p3B[2] * p3B[2] / pt2 > maxTgl23Body) { // tgLambda cut
              continue;
            }
            registry.fill(HIST("hVtx3BodyCounter"), 7.5);
            if (isTrue3bodyVtx) {
              registry.fill(HIST("hTrueVtx3BodyCounter"), 7.5);
            }

            float cosPA = (p3B[0] * dxc + p3B[1] * dyc + p3B[2] * dzc) / std::sqrt(p2candidate * (r2vertex + dzc * dzc));
            if (cosPA < minCosPA3body) {
              continue;
            }
            registry.fill(HIST("hVtx3BodyCounter"), 8.5);
            if (isTrue3bodyVtx) {
              registry.fill(HIST("hTrueVtx3BodyCounter"), 8.5);
            }

            //  Not involved: H3L DCA Check
            vtx3bodydata(
              t0.globalIndex(), t1.globalIndex(), t2.globalIndex(), collision.globalIndex(), 0,
              vertexXYZ[0], vertexXYZ[1], vertexXYZ[2],
              p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2],
              fitter3body.getChi2AtPCACandidate(),
              t0.dcaXY(), t1.dcaXY(), t2.dcaXY());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(hypertriton3bodyFinder, processMC, "Produce StoredVtx3BodyDatas with MC check", false);
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

  void processBuildLabels(aod::Collisions::iterator const& collision, aod::Vtx3BodyDatas const& vtx3bodydatas, LabeledTracks const&, aod::McParticles const& particlesMC)
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
