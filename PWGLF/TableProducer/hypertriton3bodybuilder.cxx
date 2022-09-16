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
#include "Common/DataModel/PIDResponse.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include <CCDB/BasicCCDBManager.h>
#include "DataFormatsTPC/BetheBlochAleph.h"

#include "TPCBase/ParameterGas.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"

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
#include "PWGHF/Utils/UtilsDebugLcK0Sp.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace ROOT::Math;
using namespace o2::vertexing;

using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe>;
using FullTracksExtMCIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;
using MyTracksIU = FullTracksExtIU;

// in case requested
using LabeledTracks = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

inline float GetTPCNSigmaProton(float p, float TPCSignal)
{
  float bg = p/RecoDecay::getMassPDG(kProton);
  return  (TPCSignal - o2::tpc::BetheBlochAleph(bg, -1.80365f, -20.6495f, 16.9085f, 2.15735f, -3.93285f)) / (TPCSignal*0.0922);
}
inline float GetTPCNSigmaPion(float p, float TPCSignal)
{
  float bg = p/RecoDecay::getMassPDG(kPiPlus);
  return  (TPCSignal - o2::tpc::BetheBlochAleph(bg, -2.18229f, -15.4068f, 27.3612f, 2.07923f, -3.94621f)) / (TPCSignal*0.1023);
}

struct hypertriton3bodybuilder{

  Produces<aod::StoredVtx3BodyDatas> vtx3bodydata;

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

  //for track cut in SVertexer, Can we use it in the production of goodtrack table?
  float maxDCAXY3Body = 0.3; // max DCA of 3 body decay to PV in XY // TODO RS: shall we use real chi2 to vertex?
  float maxDCAZ3Body = 0.3;  // max DCA of 3 body decay to PV in Z


  HistogramRegistry registry{
    "registry",
      {
        {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
        {"hVtx3BodyCounter", "hVtx3BodyCounter", {HistType::kTH1F, {{10, 0.0f, 10.0f}}}},
      },
  };

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  int mRunNumber;
  float d_bz;
  float maxSnp;  //max sine phi for propagation
  float maxStep; //max step size (cm) for propagation
  o2::base::MatLayerCylSet* lut = nullptr;

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

    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }

    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(2, "IfSameCollision");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(3, "hasSV");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(4, "VtxR");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(5, "TrackR");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(6, "DiffRR");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(7, "IfPropragated");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(8, "VtxPt");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(9, "tgLambda");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(10, "CosPA");
  }
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto run3grp_timestamp = bc.timestamp();

    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!grpo) {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
    }
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
      o2::base::Propagator::initFieldFromGRP(grpmag);
    }
    o2::base::Propagator::Instance()->setMatLUT(lut);
    mRunNumber = bc.runNumber();
  }
  //__________________________________________________________________

  o2::dataformats::VertexBase mMeanVertex{{0., 0., 0.}, {0.1 * 0.1, 0., 0.1 * 0.1, 0., 0., 6. * 6.}};
  void process( aod::Collision const& collision, MyTracksIU const& tracks, aod::Decays3Body const& decays3body, aod::BCsWithTimestamps const&) {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    registry.fill(HIST("hEventCounter"), 0.5);

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

    for (auto& vtx3body : decays3body) { // FIXME: turn into combination(...)

      registry.fill(HIST("hVtx3BodyCounter"), 0.5);
      auto t0 = vtx3body.track0_as<MyTracksIU>();
      auto t1 = vtx3body.track1_as<MyTracksIU>();
      auto t2 = vtx3body.track2_as<MyTracksIU>();
      auto Track0 = getTrackParCov(t0);
      auto Track1 = getTrackParCov(t1);
      auto Track2 = getTrackParCov(t2);
      if (t0.collisionId() != t1.collisionId() || t0.collisionId() != t2.collisionId() ) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 1.5);

      int n3bodyVtx = fitter3body.process(Track0, Track1, Track2);
      if (n3bodyVtx == 0) { // discard this pair
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 2.5);
      // validate V0 radial position
      // First check closeness to the beam-line as same as SVertexer
      const auto& v0XYZ = fitter3body.getPCACandidate();
      float dxv0 = v0XYZ[0] - mMeanVertex.getX(), dyv0 = v0XYZ[1] - mMeanVertex.getY(), r2v0 = dxv0 * dxv0 + dyv0 * dyv0;
      if (r2v0 < MinR2ToMeanVertex) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 3.5);

      float Track0minR =  RecoDecay::sqrtSumOfSquares(t0.x(), t0.y()), Track1minR =  RecoDecay::sqrtSumOfSquares(t1.x(), t1.y()), Track2minR = RecoDecay::sqrtSumOfSquares(t2.x(), t2.y());
      float rv0 = std::sqrt(r2v0), drv0P = rv0 - Track0minR, drv0N = rv0 - Track1minR, drv0Bach = rv0 - Track2minR;
      if (drv0P > causalityRTolerance || drv0P < -maxV0ToProngsRDiff ||
          drv0N > causalityRTolerance || drv0N < -maxV0ToProngsRDiff ||
          drv0Bach > causalityRTolerance || drv0Bach < -maxV0ToProngsRDiff) {
        LOG(debug) << "RejCausality " << drv0P << " " << drv0N;
        //continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 4.5);


      int cand3B = 0;
      const auto& vertexXYZ = fitter3body.getPCACandidatePos(cand3B);
      // make sure the cascade radius is smaller than that of the vertex
      float dxc = vertexXYZ[0] - collision.posX(), dyc = vertexXYZ[1] - collision.posY(), dzc = vertexXYZ[2] - collision.posZ(), r2vertex = dxc * dxc + dyc * dyc;
      if (std::abs(rv0 * rv0 - r2vertex) > MaxR2Diff3bodyV0 || r2vertex < MinR2ToMeanVertex) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 5.5);

      if (!fitter3body.isPropagateTracksToVertexDone() && !fitter3body.propagateTracksToVertex()) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 6.5);

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
      registry.fill(HIST("hVtx3BodyCounter"), 7.5);

      if (p3B[2] * p3B[2] / pt2 > maxTgl23Body) { // tgLambda cut
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 8.5);

      float cosPA = (p3B[0] * dxc + p3B[1] * dyc + p3B[2] * dzc) / std::sqrt(p2candidate * (r2vertex + dzc * dzc));
      if (cosPA < minCosPA3body) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 9.5);

      //Fix: Daughters DCA Check
      vtx3bodydata( 
          t0.globalIndex(), t1.globalIndex(), t2.globalIndex(), collision.globalIndex(), 
          vertexXYZ[0], vertexXYZ[1], vertexXYZ[2], 
          p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2],
          fitter3body.getChi2AtPCACandidate(),
          t0.dcaXY(), t1.dcaXY(), t2.dcaXY()
          );

    }
  }
};


struct hypertriton3bodyLabelBuilder {

  Produces<aod::McV0Labels> v0labels;

  // for bookkeeping purposes: how many V0s come from same mother etc
  HistogramRegistry registry{
    "registry",
      {
        {"hLabelCounter", "hLabelCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
        {"hHypertritonCounter", "hHypertritonCounter", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
        {"hPIDCounter", "hPIDCounter", {HistType::kTH1F, {{6, 0.0f, 6.0f}}}},
        {"hHypertriton", "hHypertriton", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hAntiHypertriton", "hAntiHypertriton", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hHypertritonMass", "hHypertritonMass", {HistType::kTH1F, {{40, 2.95f, 3.05f}}}},
        {"hAntiHypertritonMass", "hAntiHypertritonMass", {HistType::kTH1F, {{40, 2.95f, 3.05f}}}},
      },
  };

  void init(InitContext const&) {
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

  void processBuildLabels(aod::Collisions::iterator const& collision, aod::Decays3Body const& decays3body, LabeledTracks const&, aod::McParticles const& particlesMC)
  {
    for (auto& vtx : decays3body) {

      int lLabel = -1;
      int lPDG = -1;
      float lPt = -1;
      bool is3bodyDecay = false; // all considered V0s

      auto lTrack0 = vtx.track0_as<LabeledTracks>();
      auto lTrack1 = vtx.track1_as<LabeledTracks>();
      auto lTrack2 = vtx.track2_as<LabeledTracks>();
      registry.fill(HIST("hLabelCounter"), 0.5);

      // Association check
      // There might be smarter ways of doing this in the future
      if (!lTrack0.has_mcParticle() || !lTrack1.has_mcParticle() || !lTrack2.has_mcParticle()) {
        continue;
      }
      auto lMCTrack0 = lTrack0.mcParticle_as<aod::McParticles>();
      auto lMCTrack1 = lTrack1.mcParticle_as<aod::McParticles>();
      auto lMCTrack2 = lTrack2.mcParticle_as<aod::McParticles>();
      if (!lMCTrack0.has_mothers() || !lMCTrack1.has_mothers() || !lMCTrack2.has_mothers()) {
        continue;
      }

      for (auto& lMother0 : lMCTrack0.mothers_as<aod::McParticles>()) {
        for (auto& lMother1 : lMCTrack1.mothers_as<aod::McParticles>()) {
          for (auto& lMother2 : lMCTrack2.mothers_as<aod::McParticles>()) {
            if (lMother0.globalIndex() == lMother1.globalIndex() && lMother0.globalIndex() == lMother2.globalIndex()) {
              lLabel = lMother1.globalIndex();
              lPt = lMother1.pt();
              lPDG = lMother1.pdgCode();
              is3bodyDecay = true; // vtxs with the same mother
            }
          }
        }
      } // end association check
      if (!is3bodyDecay){
        continue;
      } 
      registry.fill(HIST("hLabelCounter"), 1.5);

      // Intended for cross-checks only
      // N.B. no rapidity cut!
      if (lPDG == 1010010030 && lMCTrack0.pdgCode() == 2212 && lMCTrack1.pdgCode() == -211 && lMCTrack2.pdgCode() == 1000010020)
      {
        double hypertritonMCMass = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, array{RecoDecay::getMassPDG(kProton), RecoDecay::getMassPDG(kPiPlus), 1.87561}); 
        registry.fill(HIST("hHypertritonCounter"), 0.5);
        registry.fill(HIST("hHypertriton"), lPt);
        registry.fill(HIST("hHypertritonMass"), hypertritonMCMass);
        if (TMath::Abs( lTrack0.tpcNSigmaPr()) > TpcPidNsigmaCut) {
        //if ( TMath::Abs(GetTPCNSigmaProton(lTrack0.p(), lTrack0.tpcSignal()) ) > TpcPidNsigmaCut ){
          registry.fill(HIST("hPIDCounter"), 0.5);
        }
        if( TMath::Abs(lTrack1.tpcNSigmaPi()) > TpcPidNsigmaCut ){
        //if ( TMath::Abs(GetTPCNSigmaPion(lTrack1.p(), lTrack1.tpcSignal()) ) > TpcPidNsigmaCut ){
          registry.fill(HIST("hPIDCounter"), 1.5);
        }
        if( TMath::Abs( lTrack2.tpcNSigmaDe()) > TpcPidNsigmaCut ) {
          registry.fill(HIST("hPIDCounter"), 2.5);
        }
        if (TMath::Abs( lTrack0.tpcNSigmaPr())  < TpcPidNsigmaCut && TMath::Abs(lTrack1.tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs( lTrack2.tpcNSigmaDe()) < TpcPidNsigmaCut ) {
        //if (TMath::Abs( GetTPCNSigmaProton(lTrack0.p(), lTrack0.tpcSignal()) )  < TpcPidNsigmaCut && TMath::Abs(GetTPCNSigmaPion(lTrack1.p(), lTrack1.tpcSignal() ) ) < TpcPidNsigmaCut && TMath::Abs( lTrack2.tpcNSigmaDe()) < TpcPidNsigmaCut ) {
          registry.fill(HIST("hHypertritonCounter"), 1.5);
        }
      }
      if (lPDG == -1010010030 && lMCTrack0.pdgCode() == 211 && lMCTrack1.pdgCode() == -2212 && lMCTrack2.pdgCode() == -1000010020)
      {
        double antiHypertritonMCMass = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, array{RecoDecay::getMassPDG(kPiPlus), RecoDecay::getMassPDG(kProton), 1.87561}); 
        registry.fill(HIST("hHypertritonCounter"), 2.5);
        registry.fill(HIST("hAntiHypertriton"), lPt);
        registry.fill(HIST("hAntiHypertritonMass"), antiHypertritonMCMass);
        if (TMath::Abs( lTrack0.tpcNSigmaPi()) > TpcPidNsigmaCut) {
        //if ( TMath::Abs(GetTPCNSigmaPion(lTrack0.p(), lTrack0.tpcSignal()) ) > TpcPidNsigmaCut ){
          registry.fill(HIST("hPIDCounter"), 4.5);
        }
        if( TMath::Abs(lTrack1.tpcNSigmaPr()) > TpcPidNsigmaCut ){
        //if ( TMath::Abs(GetTPCNSigmaProton(lTrack1.p(), lTrack1.tpcSignal()) ) > TpcPidNsigmaCut ){
          registry.fill(HIST("hPIDCounter"), 3.5);
        }
        if( TMath::Abs( lTrack2.tpcNSigmaDe()) > TpcPidNsigmaCut ) {
          registry.fill(HIST("hPIDCounter"), 5.5);
        }
        if (TMath::Abs( lTrack0.tpcNSigmaPi())  < TpcPidNsigmaCut && TMath::Abs(lTrack1.tpcNSigmaPr()) < TpcPidNsigmaCut && TMath::Abs( lTrack2.tpcNSigmaDe()) < TpcPidNsigmaCut ) {
        //if (TMath::Abs( GetTPCNSigmaPion(lTrack0.p(), lTrack0.tpcSignal()) )  < TpcPidNsigmaCut && TMath::Abs(GetTPCNSigmaProton(lTrack1.p(), lTrack1.tpcSignal() ) ) < TpcPidNsigmaCut && TMath::Abs( lTrack2.tpcNSigmaDe()) < TpcPidNsigmaCut ) {
          registry.fill(HIST("hHypertritonCounter"), 3.5);
        }
      }

      // Construct label table (note: this will be joinable with V0Datas)
      v0labels(
          lLabel);
    }
  }
  PROCESS_SWITCH(hypertriton3bodyLabelBuilder, processBuildLabels, "Produce MC label tables", false);
};

struct hypertriton3bodyinitializer {
  Spawns<aod::Vtx3BodyDatas> vtx3bodydatas;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertriton3bodyinitializer>(cfgc),
      adaptAnalysisTask<hypertriton3bodyLabelBuilder>(cfgc),
      adaptAnalysisTask<hypertriton3bodybuilder>(cfgc),
  };
}
