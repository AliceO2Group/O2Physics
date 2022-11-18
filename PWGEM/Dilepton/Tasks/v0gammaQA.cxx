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
// This code loops over a ambiguous track table and produces some standard analysis output.
//    Please write to: daiki.sekihata@cern.ch

#include <array>
#include <map>
#include "Math/Vector4D.h"
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
#include "DetectorsVertexing/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

using FullTracksExt = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTPCFullEl, aod::pidTPCFullPi,
                                aod::pidTPCFullKa, aod::pidTPCFullPr,
                                aod::pidTOFFullEl, aod::pidTOFFullPi,
                                aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;

using FullTrackExt = FullTracksExt::iterator;

namespace o2::aod
{

namespace v0gamma
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //!
DECLARE_SOA_COLUMN(v0Pt, pt, float);                                    //!
DECLARE_SOA_COLUMN(v0Eta, eta, float);                                  //!
DECLARE_SOA_COLUMN(v0Phi, phi, float);                                  //!
DECLARE_SOA_COLUMN(v0Mass, mass, float);                                //!
} // namespace v0gamma

// basic track information
DECLARE_SOA_TABLE(V0Gammas, "AOD", "AMBIGUOUSV0", //!
                  o2::soa::Index<>, v0gamma::CollisionId, v0gamma::PosTrackId, v0gamma::NegTrackId, v0gamma::v0Pt, v0gamma::v0Eta, v0gamma::v0Phi, v0gamma::v0Mass);

// iterators
using V0Gamma = V0Gammas::iterator;

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

TrackSelection CreateTPCOnlyTrackCuts(float max_eta, int min_tpc_ncr, float max_chi2_tpc)
{
  TrackSelection selectedTracks;
  selectedTracks.SetPtRange(0.01f, 1e10f);
  selectedTracks.SetEtaRange(-max_eta, max_eta);
  selectedTracks.SetMinNCrossedRowsTPC(min_tpc_ncr);
  selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(0.6f);
  selectedTracks.SetMaxChi2PerClusterTPC(max_chi2_tpc);
  selectedTracks.SetMaxDcaXY(2.4f);
  selectedTracks.SetMaxDcaZ(3.2f);
  return selectedTracks;
}

struct AddAmbiguousTrackInfo {
  Produces<o2::aod::AmbiguousTrackInfos> trackWithAmbiguousInfo;
  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
      {"hNonAmbTrackNbc", "Number of BCs for non-ambtrack", {HistType::kTH1F, {{101, -0.5f, 100.5f}}}},
      {"hAmbTrackNbc", "Number of BCs for ambtrack", {HistType::kTH1F, {{101, -0.5f, 100.5f}}}},
    },
  };

  // Configurables
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max. chi2/NclsTPC"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 5.0, "max. TPC n sigma for electron"};

  TrackSelection mysel;
  void init(InitContext& context)
  {
    mysel = CreateTPCOnlyTrackCuts(maxeta, mincrossedrows, maxchi2tpc);
  }

  void process(FullTracksExt const& tracks, aod::AmbiguousTracks const& ambtracks, aod::Collisions const& collisions, aod::BCs const&)
  // void process(FullTracksExt const& tracks, aod::AmbiguousTracks const& ambtracks, aod::BCs const&)
  {
    registry.fill(HIST("hEventCounter"), 1);
    // loop over ambiguous tracks
    std::vector<int> trIndices{};
    std::vector<int> atrIndices{};
    for (auto& ambtrack : ambtracks) {
      if (ambtrack.bc().size() == 1) { // only for protection. Even if nBC = 1, a track is sometimes in ambiguous track table.
        continue;
      }

      auto track = ambtrack.track_as<FullTracksExt>(); // Obtain the corresponding track
      if (mysel.IsSelected(track) && abs(track.tpcNSigmaEl()) < maxTPCNsigmaEl) {
        trIndices.push_back(track.globalIndex());
        atrIndices.push_back(ambtrack.globalIndex());
      }
    }

    for (auto& track : tracks) { // don't use filter here in order to assign isamb to all tracks.
      bool isamb = false;
      // std::size_t nBC = 1;
      unsigned int nBC = 1;
      std::vector<int> collIndices{};
      if (mysel.IsSelected(track) && abs(track.tpcNSigmaEl()) < maxTPCNsigmaEl) {
        auto trackIdx = track.globalIndex();
        auto iter = std::find(trIndices.begin(), trIndices.end(), trackIdx);
        if (iter != trIndices.end()) {
          isamb = true;
          auto ambtrack = ambtracks.rawIteratorAt(atrIndices[std::distance(trIndices.begin(), iter)]);
          nBC = ambtrack.bc().size();

          for (auto& bc : ambtrack.bc()) {
            collIndices.push_back(bc.globalIndex());
          }
        }
        if (nBC == 1) { // only for protection. Even if nBC = 1, a track is sometimes in ambiguous track table.
          isamb = false;
          nBC = 1;
          collIndices.clear();
        }
      }
      trackWithAmbiguousInfo(isamb, nBC, collIndices);

      if (isamb) {
        registry.fill(HIST("hAmbTrackNbc"), nBC);
      } else {
        registry.fill(HIST("hNonAmbTrackNbc"), nBC);
      }
    } // end of track loop
  }   // end of process
};

struct ProducePCMPhoton {
  Produces<o2::aod::V0Gammas> v0gamma;

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

  float psipairv0(const array<float, 3>& ppos, const array<float, 3>& pneg, const float bz)
  {
    // Following idea to use opening of colinear pairs in magnetic field from e.g. PHENIX to ID conversions.
    float deltat = TMath::ATan(pneg[2] / (TMath::Sqrt(pneg[0] * pneg[0] + pneg[1] * pneg[1]))) - TMath::ATan(ppos[2] / (TMath::Sqrt(ppos[0] * ppos[0] + ppos[1] * ppos[1]))); // difference of angles of the two daughter tracks with z-axis
    float pEle = RecoDecay::p(pneg);                                                                                                                                          // absolute momentum val
    float pPos = RecoDecay::p(ppos);                                                                                                                                          // absolute momentum val
    float chipair = TMath::ACos(RecoDecay::dotProd(ppos, pneg) / (pEle * pPos));                                                                                              // Angle between daughter tracks
    return TMath::Abs(TMath::ASin(deltat / chipair));                                                                                                                         // psipair in [0,pi/2]
  }

  bool checkAP(float alpha, float qt)
  {
    const float alpha_max = 0.95;
    const float qt_max = 0.05;
    float ellipse = pow(alpha / alpha_max, 2) + pow(qt / qt_max, 2);
    if (ellipse < 1.0) {
      return true;
    } else {
      return false;
    }
  }

  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
      {"hUnAssignedAmbTrackPt", "unassigned ambiguous track pT;is amb;p_{T} (GeV/c)", {HistType::kTH2F, {{2, -0.5, 1.5}, {1000, 0., 10}}}},
    },
  };

  // Configurables
  Configurable<double> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<double> minv0cospa{"minv0cospa", 0.9, "minimum V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
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
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  TrackSelection mysel;
  void init(InitContext& context)
  {
    addhistograms();
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

  void addhistograms()
  {
    const TString tracktype[2] = {"NonAmb", "Amb"};
    const TString pairtype[3] = {"NonAmbNonAmb", "NonAmbAmb", "AmbAmb"};

    // for single tracks
    for (int i = 0; i < 2; i++) {
      registry.add(Form("h%sTrackPt", tracktype[i].Data()), "pT", HistType::kTH1F, {{1000, 0.0f, 10}});
      registry.add(Form("h%sTrackEtaPhi", tracktype[i].Data()), "#eta vs. #varphi", HistType::kTH2F, {{180, 0, TMath::TwoPi()}, {40, -2.0f, 2.0f}});
      registry.add(Form("h%sTrackDCAxyz", tracktype[i].Data()), "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -5.0f, 5.0f}});
      registry.add(Form("h%sTrackNclsTPC", tracktype[i].Data()), "number of TPC clusters", HistType::kTH1F, {{161, -0.5, 160.5}});
      registry.add(Form("h%sTrackNcrTPC", tracktype[i].Data()), "number of TPC crossed rows", HistType::kTH1F, {{161, -0.5, 160.5}});
      registry.add(Form("h%sTrackChi2TPC", tracktype[i].Data()), "chi2/number of TPC clusters", HistType::kTH1F, {{100, 0, 10}});
      registry.add(Form("h%sTrackTPCdEdx", tracktype[i].Data()), "TPC dE/dx", HistType::kTH2F, {{1000, 0, 10}, {200, 0, 200}});
      registry.add(Form("h%sTrackTPCNsigmaEl", tracktype[i].Data()), "TPC n sigma el", HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}});
      registry.add(Form("h%sTrackTPCNsigmaPi", tracktype[i].Data()), "TPC n sigma pi", HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}});
      registry.add(Form("h%sTrackTPCNcr2Nf", tracktype[i].Data()), "TPC Ncr/Nfindable", HistType::kTH1F, {{200, 0, 2}});
      registry.add(Form("h%sTrackNclsITS", tracktype[i].Data()), "number of ITS clusters", HistType::kTH1F, {{8, -0.5, 7.5}});
      registry.add(Form("h%sTrackChi2ITS", tracktype[i].Data()), "chi2/number of ITS clusters", HistType::kTH1F, {{36, 0, 36}});
      registry.add(Form("h%sTrackTOFbeta", tracktype[i].Data()), "TOF beta", HistType::kTH2F, {{1000, 0, 10}, {120, 0, 1.2}});
      registry.add(Form("h%sTrackTOFNsigmaEl", tracktype[i].Data()), "TOF n sigma el", HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}});
      registry.add(Form("h%sTrackTOFNsigmaPi", tracktype[i].Data()), "TOF n sigma pi", HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}});
      registry.add(Form("h%sTrackNbc", tracktype[i].Data()), "number of bcs", HistType::kTH1F, {{101, -0.5, 100.5}});
      registry.add(Form("h%sTrackCharge", tracktype[i].Data()), "charge", HistType::kTH1F, {{3, -1.5, +1.5}});
    }

    // for V0s
    for (int i = 0; i < 3; i++) {
      registry.add(Form("h%sV0Pt", pairtype[i].Data()), "pT", HistType::kTH1F, {{1000, 0.0f, 10}});
      registry.add(Form("h%sV0EtaPhi", pairtype[i].Data()), "#eta vs. #varphi", HistType::kTH2F, {{180, 0, TMath::TwoPi()}, {40, -2.0f, 2.0f}});
      registry.add(Form("h%sV0Radius", pairtype[i].Data()), "hV0Radius; radius in Z (cm);radius in XY (cm)", HistType::kTH2F, {{500, -250, 250}, {2500, 0.0f, 250.0f}});
      registry.add(Form("h%sV0CosPA", pairtype[i].Data()), "hV0CosPA", HistType::kTH1F, {{100, 0.9f, 1.0f}});
      registry.add(Form("h%sDCAxyPosToPV", pairtype[i].Data()), "hDCAxyPosToPV;DCA_{xy} (cm);", HistType::kTH1F, {{100, -5.0f, 5.0f}});
      registry.add(Form("h%sDCAxyNegToPV", pairtype[i].Data()), "hDCAxyNegToPV;DCA_{xy} (cm);", HistType::kTH1F, {{100, -5.0f, 5.0f}});
      registry.add(Form("h%sDCAzPosToPV", pairtype[i].Data()), "hDCAzPosToPV;DCA_{z} (cm);", HistType::kTH1F, {{100, -5.0f, 5.0f}});
      registry.add(Form("h%sDCAzNegToPV", pairtype[i].Data()), "hDCAzNegToPV;DCA_{z} (cm);", HistType::kTH1F, {{100, -5.0f, 5.0f}});
      registry.add(Form("h%sDCAV0Dau", pairtype[i].Data()), "hDCAV0Dau", HistType::kTH1F, {{100, 0.0f, 10.0f}});
      registry.add(Form("h%sV0APplot", pairtype[i].Data()), "hV0APplot", HistType::kTH2F, {{200, -1.0f, +1.0f}, {250, 0.0f, 0.25f}});
      registry.add(Form("h%sGammaPsiPair", pairtype[i].Data()), "#psi_{pair} for photon conversion", HistType::kTH2F, {{160, 0.0, TMath::PiOver2()}, {100, 0.0f, 0.1f}});
      registry.add(Form("h%sMassGamma", pairtype[i].Data()), "hMassGamma", HistType::kTH2F, {{900, 0.0f, 90.0f}, {100, 0.0f, 0.1f}});
    }
  }

  Filter trackFilter = nabs(o2::aod::track::eta) < 0.9f && o2::aod::track::tpcChi2NCl < 4.f;
  using MyFilteredTracks = soa::Filtered<MyTracks>;

  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;
  void processQA(MyFilteredTracks const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const&)
  {

    for (auto& track : tracks) {
      if (!mysel.IsSelected(track) || abs(track.tpcNSigmaEl()) > maxTPCNsigmaEl) {
        continue;
      }

      if (abs(track.dcaXY()) < dcamin || dcamax < abs(track.dcaXY())) {
        continue;
      }

      if (track.isAmbiguousTrack()) {
        registry.fill(HIST("hAmbTrackNbc"), track.numBC());
        registry.fill(HIST("hAmbTrackPt"), track.pt());
        registry.fill(HIST("hAmbTrackEtaPhi"), track.phi(), track.eta());
        registry.fill(HIST("hAmbTrackDCAxyz"), track.dcaXY(), track.dcaZ());
        registry.fill(HIST("hAmbTrackNclsTPC"), track.tpcNClsFound());
        registry.fill(HIST("hAmbTrackNclsITS"), track.itsNCls());
        registry.fill(HIST("hAmbTrackNcrTPC"), track.tpcNClsCrossedRows());
        registry.fill(HIST("hAmbTrackTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
        registry.fill(HIST("hAmbTrackChi2TPC"), track.tpcChi2NCl());
        registry.fill(HIST("hAmbTrackChi2ITS"), track.itsChi2NCl());
        registry.fill(HIST("hAmbTrackTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("hAmbTrackTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
        registry.fill(HIST("hAmbTrackTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
        registry.fill(HIST("hAmbTrackTOFbeta"), track.tpcInnerParam(), track.beta());
        registry.fill(HIST("hAmbTrackTOFNsigmaEl"), track.tpcInnerParam(), track.tofNSigmaEl());
        registry.fill(HIST("hAmbTrackTOFNsigmaPi"), track.tpcInnerParam(), track.tofNSigmaPi());
        if (track.sign() < 0) {
          registry.fill(HIST("hAmbTrackCharge"), -1);
        } else {
          registry.fill(HIST("hAmbTrackCharge"), +1);
        }
      } else {
        registry.fill(HIST("hNonAmbTrackNbc"), track.numBC());
        registry.fill(HIST("hNonAmbTrackPt"), track.pt());
        registry.fill(HIST("hNonAmbTrackEtaPhi"), track.phi(), track.eta());
        registry.fill(HIST("hNonAmbTrackDCAxyz"), track.dcaXY(), track.dcaZ());
        registry.fill(HIST("hNonAmbTrackNclsTPC"), track.tpcNClsFound());
        registry.fill(HIST("hNonAmbTrackNclsITS"), track.itsNCls());
        registry.fill(HIST("hNonAmbTrackNcrTPC"), track.tpcNClsCrossedRows());
        registry.fill(HIST("hNonAmbTrackTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
        registry.fill(HIST("hNonAmbTrackChi2TPC"), track.tpcChi2NCl());
        registry.fill(HIST("hNonAmbTrackChi2ITS"), track.itsChi2NCl());
        registry.fill(HIST("hNonAmbTrackTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("hNonAmbTrackTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
        registry.fill(HIST("hNonAmbTrackTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
        registry.fill(HIST("hNonAmbTrackTOFbeta"), track.tpcInnerParam(), track.beta());
        registry.fill(HIST("hNonAmbTrackTOFNsigmaEl"), track.tpcInnerParam(), track.tofNSigmaEl());
        registry.fill(HIST("hNonAmbTrackTOFNsigmaPi"), track.tpcInnerParam(), track.tofNSigmaPi());
        if (track.sign() < 0) {
          registry.fill(HIST("hNonAmbTrackCharge"), -1);
        } else {
          registry.fill(HIST("hNonAmbTrackCharge"), +1);
        }
      }

      if (track.collisionId() < 0) {
        registry.fill(HIST("hUnAssignedAmbTrackPt"), track.isAmbiguousTrack(), track.pt());
      }

    } // end of track loop

    // Define o2 fitter, 2-prong
    o2::vertexing::DCAFitterN<2> fitter;
    fitter.setBz(d_bz); // in kG
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true); // use d_UseAbsDCA once we want to use the weighted DCA

    //printf("number of collisions = %ld , negTracks = %ld, posTracks = %ld\n", collisions.size(), negTracks.size(), posTracks.size());
    for (auto& collision : collisions) {
      auto groupEle = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());
      auto groupPos = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());
      // printf("collision.globalIndex() = %d , negTracks = %ld , posTracks = %ld\n", collision.globalIndex(), groupEle.size(), groupPos.size());
      registry.fill(HIST("hEventCounter"), 1);

      // auto const& collision = ele.collision();
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      CheckAndUpdate(bc.runNumber(), bc.timestamp());
      std::array<float, 3> pVtx = {collision.posX(), collision.posY(), collision.posZ()};
      std::array<float, 3> svpos = {0.}; // secondary vertex position
      std::array<float, 3> pvec0 = {0.};
      std::array<float, 3> pvec1 = {0.};

      // create V0 pairs
      // for (auto& [ele, pos] : combinations(CombinationsFullIndexPolicy(negTracks, posTracks))) {
      for (auto& [ele, pos] : combinations(CombinationsFullIndexPolicy(groupEle, groupPos))) {
        // if (ele.collisionId() != pos.collisionId()) {
        //     continue;
        // }
        // if (ele.collisionId() < 0 || pos.collisionId() < 0) {
        //     continue;
        // }

        if (!mysel.IsSelected(ele) || abs(ele.tpcNSigmaEl()) > maxTPCNsigmaEl) {
          continue;
        }
        if (!mysel.IsSelected(pos) || abs(pos.tpcNSigmaEl()) > maxTPCNsigmaEl) {
          continue;
        }

        if (abs(ele.dcaXY()) < dcamin || dcamax < abs(ele.dcaXY())) {
          continue;
        }
        if (abs(pos.dcaXY()) < dcamin || dcamax < abs(pos.dcaXY())) {
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
        float v0pt = RecoDecay::sqrtSumOfSquares(px, py);
        float v0eta = RecoDecay::eta(array{px, py, pz});
        float v0phi = RecoDecay::phi(px, py);

        if (v0pt < 0.05) {
          continue;
        }
        if (abs(v0eta) > v0maxeta) {
          continue;
        }

        float v0dca = fitter.getChi2AtPCACandidate(); // distance between 2 legs.
        float v0CosinePA = RecoDecay::cpa(pVtx, array{svpos[0], svpos[1], svpos[2]}, array{px, py, pz});
        float v0radius = RecoDecay::sqrtSumOfSquares(svpos[0], svpos[1]);
        float v0z = svpos[2];

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
        float psipair = psipairv0(pvec0, pvec1, d_bz);
        if (!checkAP(alpha, qtarm)) {
          continue;
        }
        if (psipair > 0.3) {
          continue;
        }

        float mGamma = RecoDecay::m(array{pvec0, pvec1}, array{RecoDecay::getMassPDG(kElectron), RecoDecay::getMassPDG(kElectron)});

        if (!ele.isAmbiguousTrack() && !pos.isAmbiguousTrack()) {
          registry.fill(HIST("hNonAmbNonAmbV0Pt"), v0pt);
          registry.fill(HIST("hNonAmbNonAmbV0EtaPhi"), v0phi, v0eta);
          registry.fill(HIST("hNonAmbNonAmbV0CosPA"), v0CosinePA);
          registry.fill(HIST("hNonAmbNonAmbDCAV0Dau"), v0dca);

          if (ele.sign() > 0) {
            registry.fill(HIST("hNonAmbNonAmbDCAxyPosToPV"), ele.dcaXY());
            registry.fill(HIST("hNonAmbNonAmbDCAzPosToPV"), ele.dcaZ());
            registry.fill(HIST("hNonAmbNonAmbDCAxyNegToPV"), pos.dcaXY());
            registry.fill(HIST("hNonAmbNonAmbDCAzNegToPV"), pos.dcaZ());
          } else {
            registry.fill(HIST("hNonAmbNonAmbDCAxyPosToPV"), pos.dcaXY());
            registry.fill(HIST("hNonAmbNonAmbDCAzPosToPV"), pos.dcaZ());
            registry.fill(HIST("hNonAmbNonAmbDCAxyNegToPV"), ele.dcaXY());
            registry.fill(HIST("hNonAmbNonAmbDCAzNegToPV"), ele.dcaZ());
          }
          registry.fill(HIST("hNonAmbNonAmbV0Radius"), v0z, v0radius);
          registry.fill(HIST("hNonAmbNonAmbV0APplot"), alpha, qtarm);
          registry.fill(HIST("hNonAmbNonAmbMassGamma"), v0radius, mGamma);
          registry.fill(HIST("hNonAmbNonAmbGammaPsiPair"), psipair, mGamma);
        } else if (ele.isAmbiguousTrack() && pos.isAmbiguousTrack()) {
          registry.fill(HIST("hAmbAmbV0Pt"), v0pt);
          registry.fill(HIST("hAmbAmbV0EtaPhi"), v0phi, v0eta);
          registry.fill(HIST("hAmbAmbV0CosPA"), v0CosinePA);
          registry.fill(HIST("hAmbAmbDCAV0Dau"), v0dca);

          if (ele.sign() > 0) {
            registry.fill(HIST("hAmbAmbDCAxyPosToPV"), ele.dcaXY());
            registry.fill(HIST("hAmbAmbDCAzPosToPV"), ele.dcaZ());
            registry.fill(HIST("hAmbAmbDCAxyNegToPV"), pos.dcaXY());
            registry.fill(HIST("hAmbAmbDCAzNegToPV"), pos.dcaZ());
          } else {
            registry.fill(HIST("hAmbAmbDCAxyPosToPV"), pos.dcaXY());
            registry.fill(HIST("hAmbAmbDCAzPosToPV"), pos.dcaZ());
            registry.fill(HIST("hAmbAmbDCAxyNegToPV"), ele.dcaXY());
            registry.fill(HIST("hAmbAmbDCAzNegToPV"), ele.dcaZ());
          }
          registry.fill(HIST("hAmbAmbV0Radius"), v0z, v0radius);
          registry.fill(HIST("hAmbAmbV0APplot"), alpha, qtarm);
          registry.fill(HIST("hAmbAmbMassGamma"), v0radius, mGamma);
          registry.fill(HIST("hAmbAmbGammaPsiPair"), psipair, mGamma);
        } else {
          registry.fill(HIST("hNonAmbAmbV0Pt"), v0pt);
          registry.fill(HIST("hNonAmbAmbV0EtaPhi"), v0phi, v0eta);
          registry.fill(HIST("hNonAmbAmbV0CosPA"), v0CosinePA);
          registry.fill(HIST("hNonAmbAmbDCAV0Dau"), v0dca);

          if (ele.sign() > 0) {
            registry.fill(HIST("hNonAmbAmbDCAxyPosToPV"), ele.dcaXY());
            registry.fill(HIST("hNonAmbAmbDCAzPosToPV"), ele.dcaZ());
            registry.fill(HIST("hNonAmbAmbDCAxyNegToPV"), pos.dcaXY());
            registry.fill(HIST("hNonAmbAmbDCAzNegToPV"), pos.dcaZ());
          } else {
            registry.fill(HIST("hNonAmbAmbDCAxyPosToPV"), pos.dcaXY());
            registry.fill(HIST("hNonAmbAmbDCAzPosToPV"), pos.dcaZ());
            registry.fill(HIST("hNonAmbAmbDCAxyNegToPV"), ele.dcaXY());
            registry.fill(HIST("hNonAmbAmbDCAzNegToPV"), ele.dcaZ());
          }
          registry.fill(HIST("hNonAmbAmbV0Radius"), v0z, v0radius);
          registry.fill(HIST("hNonAmbAmbV0APplot"), alpha, qtarm);
          registry.fill(HIST("hNonAmbAmbMassGamma"), v0radius, mGamma);
          registry.fill(HIST("hNonAmbAmbGammaPsiPair"), psipair, mGamma);
        }

        // v0gamma(pos.collisionId(), v0pt, v0eta, v0phi, mGamma);
        v0gamma(collision.globalIndex(), pos.globalIndex(), ele.globalIndex(), v0pt, v0eta, v0phi, mGamma);
      }
    }

  } // end of process
  void processDummy(soa::Join<FullTracksExt, o2::aod::AmbiguousTrackInfos> const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const&)
  {
    // do nothing
  }

  PROCESS_SWITCH(ProducePCMPhoton, processQA, "Run calo QA", false);
  PROCESS_SWITCH(ProducePCMPhoton, processDummy, "Dummy function", true);
};

struct PhotonPairing {

  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
      {"hV0Pt", "pT", {HistType::kTH1F, {{1000, 0.0, 10}}}},
      {"hNgamma", "number of photon conversion per event", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
      {"hV0EtaPhi", "#eta vs. #varphi", {HistType::kTH2F, {{180, 0, TMath::TwoPi()}, {40, -2.0f, 2.0f}}}},
      {"h2MggPt", "M_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", {HistType::kTH2F, {{400, 0.0, 0.8}, {100, 0.0, 10.}}}},
    },
  };

  void init(InitContext& context)
  {
    addhistograms();
  }

  void addhistograms()
  {
    const TString pairtype[2] = {"NonAmb", "Amb"};
    for (int i = 0; i < 2; i++) {
      registry.add(Form("h2MggPt_%s", pairtype[i].Data()), Form("M_{#gamma#gamma} vs. p_{T} %s;m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", pairtype[i].Data()), HistType::kTH2F, {{400, 0.0, 0.8}, {100, 0.0, 10.}}, true);
    }
  }

  void processNM(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Gammas const& PCMPhotons, MyTracks const&)
  {
    registry.fill(HIST("hEventCounter"), 1.0); // all
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 2.0); // FT0VX i.e. FT0and

    if (collision.numContrib() < 0.5) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 3.0); // Ncontrib > 0

    if (abs(collision.posZ()) > 10.0) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 4.0); //|Zvtx| < 10 cm

    registry.fill(HIST("hNgamma"), PCMPhotons.size());
    for (auto& g : PCMPhotons) {
      registry.fill(HIST("hV0Pt"), g.pt());
      registry.fill(HIST("hV0EtaPhi"), g.phi(), g.eta());
    }

    for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(PCMPhotons, PCMPhotons))) {
      ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), g1.mass());
      ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mass());
      ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      registry.fill(HIST("h2MggPt"), v12.M(), v12.Pt());

      auto ele1 = g1.negTrack_as<MyTracks>(); // Obtain the corresponding track
      auto pos1 = g1.posTrack_as<MyTracks>(); // Obtain the corresponding track
      auto ele2 = g2.negTrack_as<MyTracks>(); // Obtain the corresponding track
      auto pos2 = g2.posTrack_as<MyTracks>(); // Obtain the corresponding track

      if (
        (!ele1.isAmbiguousTrack() && !pos1.isAmbiguousTrack()) && (!ele2.isAmbiguousTrack() && !pos2.isAmbiguousTrack())) {
        registry.fill(HIST("h2MggPt_NonAmb"), v12.M(), v12.Pt());
      } else {
        registry.fill(HIST("h2MggPt_Amb"), v12.M(), v12.Pt());
      }

    } // end of combination

  } // end of process

  void processDummy(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Gammas const& PCMPhotons)
  {
    // do nothing
  }

  PROCESS_SWITCH(PhotonPairing, processNM, "Run gamma->ee QA", false);
  PROCESS_SWITCH(PhotonPairing, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AddAmbiguousTrackInfo>(cfgc, TaskName{"add-ambtrack-info"}),
    adaptAnalysisTask<ProducePCMPhoton>(cfgc, TaskName{"produce-pcm-photon"}),
    adaptAnalysisTask<PhotonPairing>(cfgc, TaskName{"photon-pairing"})};
}
