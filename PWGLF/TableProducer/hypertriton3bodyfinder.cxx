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


#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include <CCDB/BasicCCDBManager.h>

#include "TPCBase/ParameterGas.h"
#include "../DataModel/Vtx3BodyTables.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace ROOT::Math;
using namespace o2::vertexing;

using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullHe>;
using FullTracksExtMCIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;
//using MyTracksIU = FullTracksExtIU;
using MyTracksIU = FullTracksExtMCIU;

namespace o2::aod
{
  namespace v0goodpostracks
  {
    DECLARE_SOA_INDEX_COLUMN_FULL(GoodTrack, goodTrack, int, Tracks, "_GoodTrack");
    DECLARE_SOA_INDEX_COLUMN(Collision, collision);
    DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);
  } // namespace v0goodpostracks
  DECLARE_SOA_TABLE(V0GoodPosTracks, "AOD", "V0GOODPOSTRACKS", o2::soa::Index<>, v0goodpostracks::GoodTrackId, v0goodpostracks::CollisionId, v0goodpostracks::DCAXY);
  namespace v0goodnegtracks
  {
    DECLARE_SOA_INDEX_COLUMN_FULL(GoodTrack, goodTrack, int, Tracks, "_GoodTrack");
    DECLARE_SOA_INDEX_COLUMN(Collision, collision);
    DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);
  } // namespace v0goodnegtracks
  DECLARE_SOA_TABLE(V0GoodNegTracks, "AOD", "V0GOODNEGTRACKS", o2::soa::Index<>, v0goodnegtracks::GoodTrackId, v0goodnegtracks::CollisionId, v0goodnegtracks::DCAXY);
  namespace v0goodtracks
  {
    DECLARE_SOA_INDEX_COLUMN_FULL(GoodTrack, goodTrack, int, Tracks, "_GoodTrack");
    DECLARE_SOA_INDEX_COLUMN(Collision, collision);
    DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);
  } // namespace v0goodnegtracks
  DECLARE_SOA_TABLE(V0GoodTracks, "AOD", "V0GOODTRACKS", o2::soa::Index<>, v0goodtracks::GoodTrackId, v0goodtracks::CollisionId, v0goodtracks::DCAXY);
} // namespace o2::aod

struct trackprefilter {
  HistogramRegistry registry{
    "registry",
      {
        {"hCrossedRows", "hCrossedRows", {HistType::kTH1F, {{50, 0.0f, 200.0f}}}},
        {"hGoodTrackCount", "hGoodTrackCount",{HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
        {"hGoodPosTrackCount", "hGoodPosTrackCount",{HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
        {"hGoodNegTrackCount", "hGoodNegTrackCount",{HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
      },
  };

  //change the dca cut for helium3
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<int> tpcrefit{"tpcrefit", 0, "demand TPC refit"};

  Produces<aod::V0GoodPosTracks> v0GoodPosTracks;
  Produces<aod::V0GoodNegTracks> v0GoodNegTracks;
  Produces<aod::V0GoodTracks> v0GoodTracks;

  void process(aod::Collision const& collision,
      MyTracksIU const& tracks)
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
        v0GoodPosTracks(t0.globalIndex(), t0.collisionId(), t0.dcaXY());
        registry.fill(HIST("hGoodPosTrackCount"), 0.5);
        registry.fill(HIST("hGoodTrackCount"), 3.5);
      }
      if (t0.signed1Pt() < 0.0f) {
        v0GoodNegTracks(t0.globalIndex(), t0.collisionId(), t0.dcaXY());
        registry.fill(HIST("hGoodNegTrackCount"), 0.5);
        registry.fill(HIST("hGoodTrackCount"), 3.5);
      }
      v0GoodTracks(t0.globalIndex(), t0.collisionId(), t0.dcaXY());
    }
  }
};

struct hypertriton3bodyFinder{

  Produces<aod::StoredVtx3BodyDatas> vtx3bodydata;
  //Produces<aod::Decays3Body> decay3bodyv0;

  // Configurables
  Configurable<int> d_UseAbsDCA{"d_UseAbsDCA", kTRUE, "Use Abs DCAs"};
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};

  // Selection criteria
  Configurable<float> MinR2ToMeanVertex = {"MinR2ToMeanVertex", 0.5*0.5, ""};///< min radial distance of V0 from beam line (mean vertex) //Q: Does this cut need to be removed?

  Configurable<float> causalityRTolerance = {"causalityRTolerance", 1., ""}; ///< V0 radius cannot exceed its contributors minR by more than this value
  Configurable<float> maxV0ToProngsRDiff = {"maxV0ToProngsRDiff", 50., ""}; ///< V0 radius cannot be lower than this ammount wrt minR of contributors
  Configurable<float> MinPt2V0 = {"MinPt2V0", 1e-4, ""};   ///< v0 minimum pT
  Configurable<float> MaxTgl2V0 = {"MaxTgl2V0", 2. * 2., ""};///< maximum tgLambda of V0
  float maxDCAXYToMeanVertex3bodyV0 = 0.5;///< max DCA of V0 from beam line (mean vertex) for 3body V0 candidates
  Configurable<float> MaxDCAXY2ToMeanVertex3bodyV0 = {"MaxDCAXY2ToMeanVertex3bodyV0", 0.5*0.5, ""};
  float minCosPAXYMeanVertex3bodyV0 = 0.7;///< min cos of PA to beam line (mean vertex) in tr. plane for 3body V0 cand.
  float minCosPA3body = 0.7; // min cos of PA to PV for 3body V0

  float maxRDiffV03body = 3; ///< Maximum difference between V0 and 3body radii
  float MaxR2Diff3bodyV0 = maxRDiffV03body*maxRDiffV03body;

  float minPt3Body = 0.01;  // minimum pT of 3body V0
  float minPt23Body = 0.01*0.01;
  float maxTgl3Body = 2.;    // maximum tgLambda of 3body V0
  float maxTgl23Body = 2.*2.;
  //for 3 body reconstructed V0

  //for DCA
  //Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};

  // for TPC Track
  // Configurable<float> mMUS2TPCBin = 1.f / (8 * o2::constants::lhc::LHCBunchSpacingMUS);
  // Configurable<float> mTPCBin2Z = 0;

  //for track cut in SVertexer, Can we use it in the production of goodtrack table?
  float maxDCAXY3Body = 0.3; // max DCA of 3 body decay to PV in XY // TODO RS: shall we use real chi2 to vertex?
  float maxDCAZ3Body = 0.3;  // max DCA of 3 body decay to PV in Z


  HistogramRegistry registry{
    "registry",
      {
        {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
        {"hDauTrackCounter", "hDauTrackCounter", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
        {"hV0Counter", "hV0Counter", {HistType::kTH1F, {{8, 0.0f, 8.0f}}}},
        {"hVtx3BodyCounter", "hVtx3BodyCounter", {HistType::kTH1F, {{7, 0.0f, 7.0f}}}},
        {"hTrueV0Counter", "hTrueV0Counter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},
        {"hTrueVtx3BodyCounter", "hTrueVtx3BodyCounter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},
      },
  };


  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
  int mRunNumber;
  float d_bz;
  float maxSnp;  //max sine phi for propagation
  float maxStep; //max step size (cm) for propagation
  void init(InitContext& context)
  {
    // using namespace analysis::lambdakzerobuilder;
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  //could be changed later
    maxStep = 2.00f; //could be changed later

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    auto lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>("GLO/Config/GeometryAligned");
      // it seems this is needed at this level for the material LUT to work properly
      // but what happens if the run changes while doing the processing?  
      constexpr long run3grp_timestamp = (1619781650000 + 1619781529000) / 2;

      o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", run3grp_timestamp);
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }

    registry.get<TH1>(HIST("hDauTrackCounter"))->GetXaxis()->SetBinLabel(1, "Proton");
    registry.get<TH1>(HIST("hDauTrackCounter"))->GetXaxis()->SetBinLabel(2, "Pion");
    registry.get<TH1>(HIST("hDauTrackCounter"))->GetXaxis()->SetBinLabel(3, "Deuteron");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(1, "hasSV");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(2, "V0R");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(3, "TrackR");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(4, "IfPropragated");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(5, "V0Pt");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(6, "tgLambda");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(7, "dcaXY&CosPAXY");
    registry.get<TH1>(HIST("hV0Counter"))->GetXaxis()->SetBinLabel(8, "CosPA");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(2, "hasSV");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(3, "V0R");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(4, "TrackR");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(5, "IfPropragated");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(6, "V0Pt");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(7, "tgLambda");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(8, "dcaXY&CosPAXY");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(9, "CosPA");

    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(1, "hasSV");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(2, "VtxR");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(3, "IfPropragated");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(4, "ProdPos");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(5, "VtxPt");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(6, "tgLambda");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(7, "CosPA");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(2, "bachR");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(3, "hasSV");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(4, "VtxR");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(5, "IfPropragated");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(6, "ProdPos");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(7, "VtxPt");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(8, "tgLambda");
    registry.get<TH1>(HIST("hTrueVtx3BodyCounter"))->GetXaxis()->SetBinLabel(9, "CosPA");
  }

  float getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPObject* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    float output = grpo->getNominalL3Field();
    return output;
  }
  void CheckAndUpdate(Int_t lRunNumber, uint64_t lTimeStamp)
  {
    if (lRunNumber != mRunNumber) {
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = getMagneticField(lTimeStamp);
      } else {
        d_bz = d_bz_input;
      }
      mRunNumber = lRunNumber;
    }
  }

  //__________________________________________________________________

  o2::dataformats::VertexBase mMeanVertex{{0., 0., 0.}, {0.1 * 0.1, 0., 0.1 * 0.1, 0., 0., 6. * 6.}};
  void process( aod::Collision const& collision, MyTracksIU const& tracks, aod::V0GoodPosTracks const& ptracks, aod::V0GoodNegTracks const& ntracks, aod::V0GoodTracks const& goodtracks, aod::McParticles const& mcparticles, aod::BCsWithTimestamps const&){ 

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    CheckAndUpdate(bc.runNumber(), bc.timestamp());

    // Define o2 fitter, 2-prong
    o2::vertexing::DCAFitterN<2> fitter;
    fitter.setBz(d_bz);
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(d_UseAbsDCA);

    // check if V0 is 3-body decay
    o2::vertexing::DCAFitterN<3> fitter3body;
    fitter3body.setBz(d_bz);
    fitter3body.setPropagateToPCA(true);
    fitter3body.setMaxR(200.);//->maxRIni3body
    fitter3body.setMinParamChange(1e-3);
    fitter3body.setMinRelChi2Change(0.9);
    fitter3body.setMaxDZIni(1e9);
    fitter3body.setMaxChi2(1e9);
    fitter3body.setUseAbsDCA(d_UseAbsDCA);

    for (auto& goodtrackid : goodtracks) {
      auto goodtrack = goodtrackid.goodTrack_as<MyTracksIU>();
      if (!goodtrack.has_mcParticle()){
        continue;
      }
      auto mcgoodtrack = goodtrack.mcParticle_as<aod::McParticles>();
      if (!mcgoodtrack.has_mothers()){
        continue;
      }
      bool isHypertriton = false;
      for (auto& mothertrack : mcgoodtrack.mothers_as<aod::McParticles>()){
        if (std::abs(mothertrack.pdgCode()) != 1010010030){
          continue;
        }
        bool haveProton = false, havePion = false, haveDeuteron = false;
        bool haveAntiProton = false, haveAntiPion = false, haveAntiDeuteron = false;
        for (auto& mcparticleDaughter : mothertrack.daughters_as<aod::McParticles>()) {
          if (mcparticleDaughter.pdgCode() == 2212 ){ haveProton = true;}
          if (mcparticleDaughter.pdgCode() == -2212){ haveAntiProton = true;}
          if (mcparticleDaughter.pdgCode() == 211 ){ havePion = true;}
          if (mcparticleDaughter.pdgCode() == -211){ haveAntiPion = true;}
          if (mcparticleDaughter.pdgCode() == 1000010020 ){ haveDeuteron = true;}
          if (mcparticleDaughter.pdgCode() == -1000010020){ haveAntiDeuteron = true;}
        }
        if ( (haveProton && haveAntiPion && haveDeuteron && mothertrack.pdgCode() == 1010010030) || (haveAntiProton && havePion && haveAntiDeuteron && mothertrack.pdgCode() == -1010010030) ){
          isHypertriton = true;
          break;
        }
      }
      if (isHypertriton && std::abs(mcgoodtrack.pdgCode()) == 2212){
        registry.fill(HIST("hDauTrackCounter"), 0.5);
      }
      if (isHypertriton && std::abs(mcgoodtrack.pdgCode()) == 211){
        registry.fill(HIST("hDauTrackCounter"), 1.5);
      }
      if (isHypertriton && std::abs(mcgoodtrack.pdgCode()) == 1000010020){
        registry.fill(HIST("hDauTrackCounter"), 2.5);
      }
    }


    for (auto& t0id : ptracks) { // FIXME: turn into combination(...)
        auto t0 = t0id.goodTrack_as<MyTracksIU>();
        auto Track0 = getTrackParCov(t0);
      for (auto& t1id : ntracks) {

        if (t0id.collisionId() != t1id.collisionId()) {
          continue;
        }
        auto t1 = t1id.goodTrack_as<MyTracksIU>();
        auto Track1 = getTrackParCov(t1);

        bool isTrue3bodyV0 = false;
        if (t0.has_mcParticle() && t1.has_mcParticle()){
          auto t0mc = t0.mcParticle_as<aod::McParticles>();
          auto t1mc = t1.mcParticle_as<aod::McParticles>();
          if ( (t0mc.pdgCode() == 2212 && t1mc.pdgCode() == -211) || (t0mc.pdgCode() == 211 && t1mc.pdgCode() == -2212) ){
            if (t0mc.has_mothers() && t1mc.has_mothers()){
              for (auto& t0mother : t0mc.mothers_as<aod::McParticles>()){
                for (auto& t1mother : t1mc.mothers_as<aod::McParticles>()){
                  if (t0mother.globalIndex() == t1mother.globalIndex() && std::abs(t0mother.pdgCode()) == 1010010030){
                    isTrue3bodyV0 = true;
                  }
                }
              }
            }
          }
        }

        if(!isTrue3bodyV0){ continue;} //save time

        if(isTrue3bodyV0){
          registry.fill(HIST("hTrueV0Counter"), 0.5);
        }

        int nCand = fitter.process(Track0, Track1);
        if (nCand == 0) {
          continue;
        }
        registry.fill(HIST("hV0Counter"), 0.5);
        if(isTrue3bodyV0){
          registry.fill(HIST("hTrueV0Counter"), 1.5);
        }
        // validate V0 radial position
        // First check closeness to the beam-line as same as SVertexer
        const auto& v0XYZ = fitter.getPCACandidate();
        float dxv0 = v0XYZ[0] - mMeanVertex.getX(), dyv0 = v0XYZ[1] - mMeanVertex.getY(), r2v0 = dxv0 * dxv0 + dyv0 * dyv0;
        if (r2v0 < MinR2ToMeanVertex) {
          //continue;
        }
        registry.fill(HIST("hV0Counter"), 1.5);
        if(isTrue3bodyV0){
          registry.fill(HIST("hTrueV0Counter"), 2.5);
        }

        float Track0minR =  RecoDecay::sqrtSumOfSquares(t0.x(),  t0.y()), Track1minR =  RecoDecay::sqrtSumOfSquares(t1.x(),  t1.y());
        float rv0 = std::sqrt(r2v0), drv0P = rv0 - Track0minR, drv0N = rv0 - Track1minR;
        if (drv0P > causalityRTolerance || drv0P < -maxV0ToProngsRDiff ||
            drv0N > causalityRTolerance || drv0N < -maxV0ToProngsRDiff) {
          LOG(debug) << "RejCausality " << drv0P << " " << drv0N;
          //continue;
        }
        registry.fill(HIST("hV0Counter"), 2.5);
        if(isTrue3bodyV0){
          registry.fill(HIST("hTrueV0Counter"), 3.5);
        }

        //the process function should finish the propagation
        if (!fitter.isPropagateTracksToVertexDone() && !fitter.propagateTracksToVertex()) {
          continue;
        }
        registry.fill(HIST("hV0Counter"), 3.5);
        if(isTrue3bodyV0){
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
        if (pt2V0 < MinPt2V0) { // pt cut
          LOG(debug) << "RejPt2 " << pt2V0;
          //continue;
        }
        registry.fill(HIST("hV0Counter"), 4.5);
        if(isTrue3bodyV0){
          registry.fill(HIST("hTrueV0Counter"), 5.5);
        }

        if (pV0[2] * pV0[2] / pt2V0 > MaxTgl2V0) { // tgLambda cut
          LOG(debug) << "RejTgL " << pV0[2] * pV0[2] / pt2V0;
          //continue;
        }
        registry.fill(HIST("hV0Counter"), 5.5);
        if(isTrue3bodyV0){
          registry.fill(HIST("hTrueV0Counter"), 6.5);
        }
        float p2V0 = pt2V0 + pV0[2] * pV0[2], ptV0 = std::sqrt(pt2V0);
        // apply mass selections
        //float p2Pos = pP[0] * pP[0] + pP[1] * pP[1] + pP[2] * pP[2], p2Neg = pN[0] * pN[0] + pN[1] * pN[1] + pN[2] * pN[2];

        float dcaX = dxv0 - pV0[0] * tDCAXY, dcaY = dyv0 - pV0[1] * tDCAXY, dca2 = dcaX * dcaX + dcaY * dcaY;
        float cosPAXY = prodXYv0 / std::sqrt(r2v0 * pt2V0);
        if (dca2 > MaxDCAXY2ToMeanVertex3bodyV0 || cosPAXY < minCosPAXYMeanVertex3bodyV0) {
          //continue;
        }
        registry.fill(HIST("hV0Counter"), 6.5);
        if(isTrue3bodyV0){
          registry.fill(HIST("hTrueV0Counter"), 7.5);
        }

        float dx = v0XYZ[0] - collision.posX(), dy = v0XYZ[1] - collision.posY(), dz = v0XYZ[2] - collision.posZ(), prodXYZv0 = dx * pV0[0] + dy * pV0[1] + dz * pV0[2];
        float cosPA = prodXYZv0 / std::sqrt((dx * dx + dy * dy + dz * dz) * p2V0);
        if (cosPA < minCosPA3body) {
          //continue;
        }
        registry.fill(HIST("hV0Counter"), 7.5);
        if(isTrue3bodyV0){
          registry.fill(HIST("hTrueV0Counter"), 8.5);
        }

        for (auto& t2id : goodtracks) {
          if (t2id.globalIndex() == t0id.globalIndex()) {
            continue; // skip the track used by V0
          }
          auto t2 = t2id.goodTrack_as<MyTracksIU>();
          auto track0 = getTrackParCov(t0);
          auto track1 = getTrackParCov(t1);
          auto bach = getTrackParCov(t2);

          bool isTrue3bodyVtx = false;
          if (t0.has_mcParticle() && t1.has_mcParticle() && t2.has_mcParticle()){
            auto t0mc = t0.mcParticle_as<aod::McParticles>();
            auto t1mc = t1.mcParticle_as<aod::McParticles>();
            auto t2mc = t2.mcParticle_as<aod::McParticles>();
            if ( (t0mc.pdgCode() == 2212 && t1mc.pdgCode() == -211 && t2mc.pdgCode() == 1000010020) || (t0mc.pdgCode() == 211 && t1mc.pdgCode() == -2212 && t2mc.pdgCode() == -1000010020) ){
              if (t0mc.has_mothers() && t1mc.has_mothers() && t2mc.has_mothers()){
                for (auto& t0mother : t0mc.mothers_as<aod::McParticles>()){
                  for (auto& t1mother : t1mc.mothers_as<aod::McParticles>()){
                    for (auto& t2mother : t2mc.mothers_as<aod::McParticles>()){
                      if (t0mother.globalIndex() == t1mother.globalIndex() && t0mother.globalIndex() == t2mother.globalIndex() && std::abs(t0mother.pdgCode()) == 1010010030){
                        isTrue3bodyVtx = true;
                      }
                    }
                  }

                }
              }
            }
          }

          if (!isTrue3bodyVtx){ continue;} // save time
          if (isTrue3bodyVtx){
            registry.fill(HIST("hTrueVtx3BodyCounter"), 0.5);
          }

          float bachminR = RecoDecay::sqrtSumOfSquares(t2.x(),  t2.y());
          if (bachminR > rv0 + causalityRTolerance) {
            //continue;
          }
          if (isTrue3bodyVtx){
            registry.fill(HIST("hTrueVtx3BodyCounter"), 1.5);
          }

          int n3bodyVtx = fitter3body.process(track0, track1, bach);
          if (n3bodyVtx == 0) { // discard this pair
            continue;
          }
          registry.fill(HIST("hVtx3BodyCounter"), 0.5);
          if (isTrue3bodyVtx){
            registry.fill(HIST("hTrueVtx3BodyCounter"), 2.5);
          }

          int cand3B = 0;
          const auto& vertexXYZ = fitter3body.getPCACandidatePos(cand3B);
          // make sure the cascade radius is smaller than that of the vertex
          float dxc = vertexXYZ[0] - collision.posX(), dyc = vertexXYZ[1] - collision.posY(), dzc = vertexXYZ[2] - collision.posZ(), r2vertex = dxc * dxc + dyc * dyc;
          /*if (std::abs(rv0 * rv0 - r2vertex) > MaxR2Diff3bodyV0 || r2vertex < MinR2ToMeanVertex) {
            continue;
          }*/
          registry.fill(HIST("hVtx3BodyCounter"), 1.5);
          if (isTrue3bodyVtx){
            registry.fill(HIST("hTrueVtx3BodyCounter"), 3.5);
          }

          if (!fitter3body.isPropagateTracksToVertexDone() && !fitter3body.propagateTracksToVertex()) {
            continue;
          }
          registry.fill(HIST("hVtx3BodyCounter"), 2.5);
          if (isTrue3bodyVtx){
            registry.fill(HIST("hTrueVtx3BodyCounter"), 4.5);
          }

          auto& tr0 = fitter3body.getTrack(0, cand3B);
          auto& tr1 = fitter3body.getTrack(1, cand3B);
          auto& tr2 = fitter3body.getTrack(2, cand3B);
          std::array<float, 3> p0, p1, p2;
          tr0.getPxPyPzGlo(p0);
          tr1.getPxPyPzGlo(p1);
          tr2.getPxPyPzGlo(p2);
          std::array<float, 3> p3B = {p0[0] + p1[0] + p2[0], p0[1] + p1[1] + p2[1], p0[2] + p1[2] + p2[2]};
          /*auto prodPPos = pV0[0] * dxc + pV0[1] * dyc + pV0[2] * dzc;
          if (prodPPos < 0.) { // causality cut
            continue;
          }*/
          registry.fill(HIST("hVtx3BodyCounter"), 3.5);
          if (isTrue3bodyVtx){
            registry.fill(HIST("hTrueVtx3BodyCounter"), 5.5);
          }

          float pt2 = p3B[0] * p3B[0] + p3B[1] * p3B[1], p2candidate = pt2 + p3B[2] * p3B[2];
          if (pt2 < minPt23Body) { // pt cut
            continue;
          }
          registry.fill(HIST("hVtx3BodyCounter"), 4.5);
          if (isTrue3bodyVtx){
            registry.fill(HIST("hTrueVtx3BodyCounter"), 6.5);
          }

          if (p3B[2] * p3B[2] / pt2 > maxTgl23Body) { // tgLambda cut
            continue;
          }
          registry.fill(HIST("hVtx3BodyCounter"), 5.5);
          if (isTrue3bodyVtx){
            registry.fill(HIST("hTrueVtx3BodyCounter"), 7.5);
          }

          float cosPA = (p3B[0] * dxc + p3B[1] * dyc + p3B[2] * dzc) / std::sqrt(p2candidate * (r2vertex + dzc * dzc));
          if (cosPA < minCosPA3body) {
            continue;
          }
          registry.fill(HIST("hVtx3BodyCounter"), 6.5);
          if (isTrue3bodyVtx){
            registry.fill(HIST("hTrueVtx3BodyCounter"), 8.5);
          }

          //Fix: Daughters DCA Check
          vtx3bodydata( 
              t0.globalIndex(), t1.globalIndex(), t2.globalIndex(), collision.globalIndex(), 
              vertexXYZ[0], vertexXYZ[1], vertexXYZ[2], 
              p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2],
              fitter3body.getChi2AtPCACandidate(),
              t0.dcaXY(), t1.dcaXY(), t2.dcaXY());
        }

      }
    }
  }



  //__________________________________________________________________
  //Track is accepted when X is limited and propagate is fine
  /*bool acceptTrack(GIndex gid, const o2::track::TrackParCov& trc) const
    {
  // DCA to mean vertex
  o2::track::TrackPar trp(trc);
  std::array<float, 2> dca;
  auto* prop = o2::base::Propagator::Instance();
  if (mSVParams->usePropagator) {
  if (trp.getX() > mSVParams->minRFor3DField && !prop->PropagateToXBxByBz(trp, mSVParams->minRFor3DField, mSVParams->maxSnp, mSVParams->maxStep, o2::base::Propagator::MatCorrType(mSVParams->matCorr))) {
  return true; // we don't need actually to propagate to the beam-line
  }
  if (!prop->propagateToDCA(mMeanVertex.getXYZ(), trp, prop->getNominalBz(), mSVParams->maxStep, o2::base::Propagator::MatCorrType(mSVParams->matCorr), &dca)) {
  return true;
  }
  } else {
  if (!trp.propagateParamToDCA(mMeanVertex.getXYZ(), prop->getNominalBz(), &dca)) {
  return true;
  }
  }
  }



  //______________________________________________
  bool processTPCTrack(const o2::tpc::TrackTPC& trTPC, GIndex gid, int vtxid)
  {
  // if TPC trackis unconstrained, try to create in the tracks pool a clone constrained to vtxid vertex time.
  if (trTPC.hasBothSidesClusters()) { // this is effectively constrained track
  return false;                     // let it be processed as such
  }
  const auto& vtx = mPVertices[vtxid];
  auto twe = vtx.getTimeStamp();
  int posneg = trTPC.getSign() < 0 ? 1 : 0;
  auto trLoc = mTracksPool[posneg].emplace_back(TrackCand{trTPC, gid, {vtxid, vtxid}, 0.});
  auto err = correctTPCTrack(trLoc, trTPC, twe.getTimeStamp(), twe.getTimeStampError());
  if (err < 0) {
  mTracksPool[posneg].pop_back(); // discard
  }
  trLoc.minR = std::sqrt(trLoc.getX() * trLoc.getX() + trLoc.getY() * trLoc.getY());
  return true;
  }

  //______________________________________________
  float correctTPCTrack(o2::track::TrackParCov& trc, const o2::tpc::TrackTPC tTPC, float tmus, float tmusErr) const
  {
  // Correct the track copy trc of the TPC track for the assumed interaction time
  // return extra uncertainty in Z due to the interaction time uncertainty
  // TODO: at the moment, apply simple shift, but with Z-dependent calibration we may
  // need to do corrections on TPC cluster level and refit
  // This is a clone of MatchTPCITS::correctTPCTrack
  float dDrift = (tmus * mMUS2TPCBin - tTPC.getTime0()) * mTPCBin2Z;
  float driftErr = tmusErr * mMUS2TPCBin * mTPCBin2Z;
  // eventually should be refitted, at the moment we simply shift...
  trc.setZ(tTPC.getZ() + (tTPC.hasASideClustersOnly() ? dDrift : -dDrift));
  trc.setCov(trc.getSigmaZ2() + driftErr * driftErr, o2::track::kSigZ2);

  return driftErr;
  }*/
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
      adaptAnalysisTask<hypertriton3bodyInitializer>(cfgc),
  };
}
