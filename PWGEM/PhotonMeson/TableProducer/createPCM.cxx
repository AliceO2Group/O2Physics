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
// ========================
//
// This code produces photon data tables.
//    Please write to: daiki.sekihata@cern.ch

#include <array>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

using FullTracksExt = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTPCFullEl, aod::pidTPCFullPi,
                                aod::pidTPCFullKa, aod::pidTPCFullPr>;

using FullTrackExt = FullTracksExt::iterator;

namespace o2::aod
{

namespace ambtrackinfo
{
// Columns to store the information about ambiguous tracks joinable with the track table
DECLARE_SOA_COLUMN(IsAmbiguousTrack, isAmbiguousTrack, bool);                              //!
DECLARE_SOA_COLUMN(NumBC, numBC, unsigned int);                                            //!
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(AmbiguousCollisionIndices, ambiguousCollisionIndices); //!
} // namespace ambtrackinfo

DECLARE_SOA_TABLE(AmbiguousTrackInfos, "AOD", "AMBTRACKINFOS", //!
                  ambtrackinfo::IsAmbiguousTrack,
                  ambtrackinfo::NumBC,
                  ambtrackinfo::AmbiguousCollisionIndicesIds);

using AmbiguousTrackInfo = AmbiguousTrackInfos::iterator;

} // namespace o2::aod

using MyTracks = soa::Join<FullTracksExt, o2::aod::AmbiguousTrackInfos>;
using MyTrack = MyTracks::iterator;

struct AddAmbiguousTrackInfo {
  Produces<o2::aod::AmbiguousTrackInfos> trackWithAmbiguousInfo;
  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
      {"hNonAmbTrackNbc", "Number of BCs for non-ambtrack", {HistType::kTH1F, {{101, -0.5f, 100.5f}}}},
      {"hAmbTrackNbc", "Number of BCs for ambtrack", {HistType::kTH1F, {{101, -0.5f, 100.5f}}}},
      {"hNonAmbTrackNcoll", "Number of collisions for non-ambtrack", {HistType::kTH1F, {{101, -0.5f, 100.5f}}}},
      {"hAmbTrackNcoll", "Number of collisions for ambtrack", {HistType::kTH1F, {{101, -0.5f, 100.5f}}}},
    },
  };

  // Configurables
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max. chi2/NclsTPC"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 5.0, "max. TPC n sigma for electron"};
  Configurable<float> dcamin{"dcamin", 0.0, "dcamin"};
  Configurable<float> dcamax{"dcamax", 1e+10, "dcamax"};

  TrackSelection mysel;
  void init(InitContext& context)
  {
    mysel = CreateTPCOnlyTrackCuts(maxeta, mincrossedrows, maxchi2tpc);
  }

  void process(FullTracksExt const& tracks, aod::AmbiguousTracks const& ambtracks, aod::Collisions const& collisions, aod::BCs const&)
  {
    registry.fill(HIST("hEventCounter"), 1);
    // loop over ambiguous tracks
    std::vector<int> trIndices{};
    std::vector<int> atrIndices{};
    for (auto& ambtrack : ambtracks) {
      auto track = ambtrack.track_as<FullTracksExt>(); // Obtain the corresponding track
      if (mysel.IsSelected(track) && abs(track.tpcNSigmaEl()) < maxTPCNsigmaEl) {
        trIndices.push_back(track.globalIndex());
        atrIndices.push_back(ambtrack.globalIndex());
      }
    }

    for (auto& track : tracks) { // don't use filter here in order to assign isamb to all tracks.
      bool isamb = false;
      unsigned int nBC = 1;
      std::vector<int> collIndices{};
      if (mysel.IsSelected(track) && abs(track.tpcNSigmaEl()) < maxTPCNsigmaEl && (dcamin < abs(track.dcaXY()) && abs(track.dcaXY()) < dcamax)) {
        auto trackIdx = track.globalIndex();
        auto iter = std::find(trIndices.begin(), trIndices.end(), trackIdx);
        if (iter != trIndices.end()) {
          isamb = true;
          auto ambtrack = ambtracks.rawIteratorAt(atrIndices[std::distance(trIndices.begin(), iter)]);
          nBC = ambtrack.bc().size();

          // for (auto& bc : ambtrack.bc()) {
          //   collIndices.push_back(bc.globalIndex());
          // }

          // for (auto& collision : collisions) {
          //     uint64_t mostProbableBC = collision.bc().globalBC();
          //     for (auto& bc : ambtrack.bc()) {
          //         if (bc.globalBC() == mostProbableBC) {
          //             collIndices.push_back(collision.globalIndex());
          //             break;
          //         }
          //     } // end of bc loop in ambtrack
          // } // end of collision loop
        }
      }
      trackWithAmbiguousInfo(isamb, nBC, collIndices);
      if (isamb) {
        registry.fill(HIST("hAmbTrackNbc"), nBC);
        registry.fill(HIST("hAmbTrackNcoll"), collIndices.size());
      } else {
        registry.fill(HIST("hNonAmbTrackNbc"), nBC);
        registry.fill(HIST("hNonAmbTrackNcoll"), collIndices.size());
      }

      collIndices.clear();
      std::vector<int>().swap(collIndices);
    } // end of track loop
  }   // end of process
};

struct createPCM {
  Produces<o2::aod::V0Photons> v0photon;
  Produces<o2::aod::V0Legs> v0leg;

  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
    },
  };

  // Configurables
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Configurable<float> v0max_mee{"v0max_mee", 0.1, "maximum mee for photon conversion in GeV"};
  Configurable<float> maxpsipair{"maxpsipair", 1.6, "maximum psipair for v0"};
  Configurable<float> minv0cospa{"minv0cospa", 0.9, "minimum V0 CosPA"};
  Configurable<float> maxdcav0dau{"maxdcav0dau", 1.5, "max DCA between V0 Daughters"};
  Configurable<float> v0Rmin{"v0Rmin", 0.0, "v0Rmin"};
  Configurable<float> v0Rmax{"v0Rmax", 60.0, "v0Rmax"};
  Configurable<float> dcamin{"dcamin", 0.0, "dcamin"};
  Configurable<float> dcamax{"dcamax", 1e+10, "dcamax"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance for single track"};
  Configurable<float> v0maxeta{"v0maxeta", 0.9, "eta acceptance for v0"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max chi2/NclsTPC"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 5.0, "max. TPC n sigma for electron"};

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  TrackSelection mysel;
  void init(InitContext& context)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }
    mysel = CreateTPCOnlyTrackCuts(maxeta, mincrossedrows, maxchi2tpc);
  }

  float getMagneticField(uint64_t timestamp)
  {
    o2::base::MatLayerCylSet* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    static o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, timestamp);
    if (grpo != nullptr) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
    } else {
      LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
      return 0;
    }
    LOGF(info, "Retrieved GRP for timestamp %llu with L3 ", timestamp, grpo->getL3Current());
    float bz = std::lround(5.f * grpo->getL3Current() / 30000.f); // in kG
    LOGF(info, "magnetig field = %f kG", bz);
    return bz;
  }

  void CheckAndUpdate(Int_t lRunNumber, uint64_t lTimeStamp)
  {
    if (lRunNumber != mRunNumber) {
      if (abs(d_bz_input) > 99) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = getMagneticField(lTimeStamp);
      } else {
        d_bz = d_bz_input;
      }
      mRunNumber = lRunNumber;
    }
  }

  float alphav0(const array<float, 3>& ppos, const array<float, 3>& pneg)
  {
    std::array<float, 3> pv0 = {ppos[0] + pneg[0], ppos[1] + pneg[1], ppos[2] + pneg[2]};
    float momTot = RecoDecay::p(pv0);
    float lQlNeg = RecoDecay::dotProd(pneg, pv0) / momTot;
    float lQlPos = RecoDecay::dotProd(ppos, pv0) / momTot;
    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg); // longitudinal momentum asymmetry
  }

  float qtarmv0(const array<float, 3>& ppos, const array<float, 3>& pneg)
  {
    std::array<float, 3> pv0 = {ppos[0] + pneg[0], ppos[1] + pneg[1], ppos[2] + pneg[2]};
    float momTot2 = RecoDecay::p2(pv0);
    float dp = RecoDecay::dotProd(pneg, pv0);
    return std::sqrt(RecoDecay::p2(pneg) - dp * dp / momTot2); // qtarm
  }

  float psipairv0(const array<float, 3>& ppos, const array<float, 3>& pneg)
  {
    // Following idea to use opening of colinear pairs in magnetic field from e.g. PHENIX to ID conversions.
    float deltat = TMath::ATan(pneg[2] / (TMath::Sqrt(pneg[0] * pneg[0] + pneg[1] * pneg[1]))) - TMath::ATan(ppos[2] / (TMath::Sqrt(ppos[0] * ppos[0] + ppos[1] * ppos[1]))); // difference of angles of the two daughter tracks with z-axis
    float pEle = RecoDecay::p(pneg);                                                                                                                                          // absolute momentum val
    float pPos = RecoDecay::p(ppos);                                                                                                                                          // absolute momentum val
    float chipair = TMath::ACos(RecoDecay::dotProd(ppos, pneg) / (pEle * pPos));                                                                                              // Angle between daughter tracks
    return TMath::Abs(TMath::ASin(deltat / chipair));                                                                                                                         // psipair in [0,pi/2]
  }

  Filter trackFilter = nabs(o2::aod::track::eta) < maxeta && dcamin < nabs(o2::aod::track::dcaXY) && nabs(o2::aod::track::dcaXY) < dcamax;
  using MyFilteredTracks = soa::Filtered<MyTracks>;

  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;
  void process(MyFilteredTracks const& tracks, soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::BCsWithTimestamps const&)
  {
    // printf("number of collisions = %ld , negTracks = %ld, posTracks = %ld\n", collisions.size(), negTracks.size(), posTracks.size());
    for (auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      CheckAndUpdate(bc.runNumber(), bc.timestamp());

      auto groupEle = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());
      auto groupPos = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());
      // printf("collision.globalIndex() = %d , negTracks = %ld , posTracks = %ld\n", collision.globalIndex(), groupEle.size(), groupPos.size());

      // Define o2 fitter, 2-prong, this has to be defined after setting magnetic field!
      o2::vertexing::DCAFitterN<2> fitter;
      fitter.setBz(d_bz); // in kG
      fitter.setPropagateToPCA(true);
      fitter.setMaxR(200.);
      fitter.setMinParamChange(1e-3);
      fitter.setMinRelChi2Change(0.9);
      fitter.setMaxDZIni(1e9);
      fitter.setMaxChi2(1e9);
      fitter.setUseAbsDCA(true);

      array<float, 3> pVtx = {collision.posX(), collision.posY(), collision.posZ()};
      array<float, 3> svpos = {0.}; // secondary vertex position
      array<float, 3> pvec0 = {0.};
      array<float, 3> pvec1 = {0.};

      // create V0 pairs
      for (auto& [ele, pos] : combinations(CombinationsFullIndexPolicy(groupEle, groupPos))) {
        if (!mysel.IsSelected(ele) || abs(ele.tpcNSigmaEl()) > maxTPCNsigmaEl) {
          continue;
        }
        if (!mysel.IsSelected(pos) || abs(pos.tpcNSigmaEl()) > maxTPCNsigmaEl) {
          continue;
        }

        auto pTrack = getTrackParCov(pos); // positive
        auto nTrack = getTrackParCov(ele); // negative

        int nCand = fitter.process(pTrack, nTrack);
        if (nCand != 0) {
          fitter.propagateTracksToVertex();
          const auto& vtx = fitter.getPCACandidate();
          for (int i = 0; i < 3; i++) {
            svpos[i] = vtx[i];
          }
          fitter.getTrack(0).getPxPyPzGlo(pvec0); // positive
          fitter.getTrack(1).getPxPyPzGlo(pvec1); // negative
        } else {
          continue;
        }

        float px = pvec0[0] + pvec1[0];
        float py = pvec0[1] + pvec1[1];
        float pz = pvec0[2] + pvec1[2];
        // float v0pt = RecoDecay::sqrtSumOfSquares(px, py);
        float v0eta = RecoDecay::eta(array{px, py, pz});
        // float v0phi = RecoDecay::phi(px, py);

        if (abs(v0eta) > v0maxeta) {
          continue;
        }

        float v0dca = fitter.getChi2AtPCACandidate(); // distance between 2 legs.
        float v0CosinePA = RecoDecay::cpa(pVtx, array{svpos[0], svpos[1], svpos[2]}, array{px, py, pz});
        float v0radius = RecoDecay::sqrtSumOfSquares(svpos[0], svpos[1]);

        if (v0dca > maxdcav0dau) {
          continue;
        }
        if (v0radius < v0Rmin || v0Rmax < v0radius) {
          continue;
        }
        if (v0CosinePA < minv0cospa) {
          continue;
        }

        float alpha = alphav0(pvec0, pvec1);
        float qtarm = qtarmv0(pvec0, pvec1);
        if (!checkAP(alpha, qtarm))
          continue;

        float psipair = psipairv0(pvec0, pvec1);
        if (psipair > maxpsipair)
          continue;

        float mGamma = RecoDecay::m(array{pvec0, pvec1}, array{RecoDecay::getMassPDG(kElectron), RecoDecay::getMassPDG(kElectron)});
        if (mGamma > v0max_mee)
          continue;

        v0photon(collision.globalIndex(), v0leg.lastIndex() + 1, v0leg.lastIndex() + 2,
                 svpos[0], svpos[1], svpos[2],
                 pvec0[0], pvec0[1], pvec0[2],
                 pvec1[0], pvec1[1], pvec1[2], v0CosinePA, v0dca); // if v0leg is empty, lastIndex = -1.
        v0leg(collision.globalIndex(),
              pos.globalIndex(), pos.sign(), pos.isAmbiguousTrack(),
              pos.pt(), pos.eta(), pos.phi(), pos.dcaXY(), pos.dcaZ(),
              pos.tpcNClsFindable(), pos.tpcNClsFindableMinusFound(), pos.tpcNClsFindableMinusCrossedRows(),
              pos.tpcChi2NCl(), pos.tpcInnerParam(), pos.tpcSignal(),
              pos.tpcNSigmaEl(), pos.tpcNSigmaPi(), pos.tpcNSigmaKa(), pos.tpcNSigmaPr(),
              pos.itsClusterMap(), pos.itsChi2NCl(), pos.detectorMap());
        v0leg(collision.globalIndex(),
              ele.globalIndex(), ele.sign(), ele.isAmbiguousTrack(),
              ele.pt(), ele.eta(), ele.phi(), ele.dcaXY(), ele.dcaZ(),
              ele.tpcNClsFindable(), ele.tpcNClsFindableMinusFound(), ele.tpcNClsFindableMinusCrossedRows(),
              ele.tpcChi2NCl(), ele.tpcInnerParam(), ele.tpcSignal(),
              ele.tpcNSigmaEl(), ele.tpcNSigmaPi(), ele.tpcNSigmaKa(), ele.tpcNSigmaPr(),
              ele.itsClusterMap(), ele.itsChi2NCl(), ele.detectorMap());

      } // end of pairing loop

    } // end of collision loop

  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AddAmbiguousTrackInfo>(cfgc, TaskName{"add-ambtrack-info"}),
    adaptAnalysisTask<createPCM>(cfgc, TaskName{"create-pcm"}),
  };
}
