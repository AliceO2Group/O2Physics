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
#include "Common/DataModel/StrangenessTables.h"
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
using MyTracksIU = FullTracksExtIU;


struct hypertriton3bodybuilder{

  Produces<aod::StoredVtx3BodyDatas> vtx3bodydata;

  // Configurables
  Configurable<double> d_UseAbsDCA{"d_UseAbsDCA", kTRUE, "Use Abs DCAs"};
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};

  // Selection criteria
  Configurable<float> MinR2ToMeanVertex = {"MinR2ToMeanVertex", 0.5*0.5, ""};///< min radial distance of V0 from beam line (mean vertex) //Q: Does this cut need to be removed?

  Configurable<float> causalityRTolerance = {"causalityRTolerance", 1., ""}; ///< V0 radius cannot exceed its contributors minR by more than this value
  Configurable<float> maxV0ToProngsRDiff = {"maxV0ToProngsRDiff", 50., ""}; ///< V0 radius cannot be lower than this ammount wrt minR of contributors
  Configurable<float> MinPt2V0 = {"MinPt2V0", 1e-4, ""};   ///< v0 minimum pT
  Configurable<float> MaxTgl2V0 = {"MaxTgl2V0", 2. * 2., ""};///< maximum tgLambda of V0
  float maxDCAXYToMeanVertex3bodyV0 = 0.5;///< max DCA of V0 from beam line (mean vertex) for 3body V0 candidates
  Configurable<float> MaxDCAXY2ToMeanVertex3bodyV0 = {"MaxDCAXY2ToMeanVertex3bodyV0", 0.5*0.5, ""};
  float minCosPAXYMeanVertex3bodyV0 = 0.8;///< min cos of PA to beam line (mean vertex) in tr. plane for 3body V0 cand.
  float minCosPA3body = 0.7; // min cos of PA to PV for 3body V0

  float maxRDiffV03body = 0.2; ///< Maximum difference between V0 and 3body radii
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
      /* it seems this is needed at this level for the material LUT to work properly */
      /* but what happens if the run changes while doing the processing?             */
      constexpr long run3grp_timestamp = (1619781650000 + 1619781529000) / 2;

      o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", run3grp_timestamp);
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }

    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(1, "hasSV");
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
  void process( aod::Collision const& collision, MyTracksIU const& tracks, aod::Decays3Body const& decays3body, aod::BCsWithTimestamps const&) {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    CheckAndUpdate(bc.runNumber(), bc.timestamp());

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
        registry.fill(HIST("hV0Counter"), 3.5);

      float Track0minR =  RecoDecay::sqrtSumOfSquares(t0.x(), t0.y()), Track1minR =  RecoDecay::sqrtSumOfSquares(t1.x(), t1.y()), Track2minR = RecoDecay::sqrtSumOfSquares(t2.x(), t2.y());
      float rv0 = std::sqrt(r2v0), drv0P = rv0 - Track0minR, drv0N = rv0 - Track1minR, drv0Bach = rv0 - Track2minR;
      if (drv0P > causalityRTolerance || drv0P < -maxV0ToProngsRDiff ||
          drv0N > causalityRTolerance || drv0N < -maxV0ToProngsRDiff ||
          drv0Bach > causalityRTolerance || drv0Bach < -maxV0ToProngsRDiff) {
        LOG(debug) << "RejCausality " << drv0P << " " << drv0N;
        continue;
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
      //Here is a mass check with Hyps and a cut for dca of candidatesin check3bodyDecays
      vtx3bodydata( 
          t0.globalIndex(), t1.globalIndex(), t2.globalIndex(), collision.globalIndex(), 
          vertexXYZ[0], vertexXYZ[1], vertexXYZ[2], 
          t0.px(), t0.py(), t0.pz(), t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(),
          t0.dcaXY(), t1.dcaXY(), t2.dcaXY()
          );

    }
  }
};

struct hypertriton3bodyinitializer {
  Spawns<aod::Vtx3BodyDatas> vtx3bodydatas;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertriton3bodyinitializer>(cfgc),
      adaptAnalysisTask<hypertriton3bodybuilder>(cfgc),
  };
}
