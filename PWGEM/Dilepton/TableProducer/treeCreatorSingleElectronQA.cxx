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

/// \brief write relevant information for dalitz ee analysis to an AO2D.root file. This file is then the only necessary input to perform pcm analysis.
/// \author daiki.sekihata@cern.ch

#include <random>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using MyCollisions_Cent = soa::Join<MyCollisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>; // centrality table has dependency on multiplicity table.

using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;
using MyCollisionsMC_Cent = soa::Join<MyCollisionsMC, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>; // centrality table has dependency on multiplicity table.

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov,
                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyTrack = MyTracks::iterator;
using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels>;
using MyTrackMC = MyTracksMC::iterator;

struct TreeCreatorSingleElectronQA {

  Produces<o2::aod::EMEvents> event;
  Produces<o2::aod::EMEventsMult> event_mult;
  Produces<o2::aod::EMEventsCent> event_cent;
  Produces<o2::aod::EMPrimaryElectrons> emprimaryelectrons;
  Produces<o2::aod::EMPrimaryElectronEMEventIds> prmeleventid;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};

  // Operation and minimisation criteria
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 0, "min ncluster tpc"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  Configurable<int> minitsncls{"minitsncls", 4, "min. number of ITS clusters"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max. chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 6.0, "max. chi2/NclsITS"};
  Configurable<float> minpt{"minpt", 0.1, "min pt for track"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<float> dca_xy_max{"dca_xy_max", 1.0f, "max DCAxy in cm"};
  Configurable<float> dca_z_max{"dca_z_max", 1.0f, "max DCAz in cm"};
  Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -1e+10, "min. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", +1e+10, "max. TPC n sigma for electron inclusion"};
  Configurable<float> minTPCNsigmaPi{"minTPCNsigmaPi", 0.0, "min. TPC n sigma for pion exclusion"};
  Configurable<float> maxTPCNsigmaPi{"maxTPCNsigmaPi", 0.0, "max. TPC n sigma for pion exclusion"};
  Configurable<float> down_scaling{"down_scaling", 1e-3, "down scaling factor to store charged particles"};

  std::pair<int8_t, std::set<uint8_t>> itsRequirement = {1, {0, 1, 2}}; // any hits on 3 ITS ib layers.

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  std::mt19937 engine;
  std::uniform_real_distribution<float> dist01;

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    addHistograms();
    std::random_device seed_gen;
    engine = std::mt19937(seed_gen());
    dist01 = std::uniform_real_distribution<float>(0.0f, 1.0f);
  }

  void addHistograms()
  {
    fRegistry.add("hCollisionCounter", "collision counter", kTH1F, {{2, -0.5f, 1.5}}, false);
    // for track info
    fRegistry.add("Track/positive/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/positive/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/positive/hPtPhi", "pT vs. #varphi;#varphi (rad.);p_{T} (GeV/c)", kTH2F, {{360, 0, 2 * M_PI}, {1000, 0.0f, 10.f}}, false);
    fRegistry.add("Track/positive/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{360, 0, 2 * M_PI}, {40, -2.0f, 2.0f}}, false);
    fRegistry.add("Track/positive/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/positive/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/positive/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{1000, 0, 10}, {100, 0., 1000}}, false);
    fRegistry.add("Track/positive/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{1000, 0, 10}, {100, 0., 1000}}, false);
    fRegistry.add("Track/positive/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/positive/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/positive/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/positive/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/positive/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTOFbeta", "TOF beta;p_{in} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {600, 0, 1.2}}, false);
    fRegistry.add("Track/positive/h1OverTOFbeta", "1/TOF beta;p_{in} (GeV/c);1/#beta", kTH2F, {{1000, 0, 10}, {1000, 0.8, 1.8}}, false);
    fRegistry.add("Track/positive/hTOFNsigmaEl", "TOF n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTOFNsigmaMu", "TOF n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTOFNsigmaPi", "TOF n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTOFNsigmaKa", "TOF n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTOFNsigmaPr", "TOF n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/positive/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/positive/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("Track/positive/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/positive/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.addClone("Track/positive/", "Track/negative/");
  }

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
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
  }

  template <bool isMC, typename TTrack>
  bool checkTrack(TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    if (track.pt() < minpt || abs(track.eta()) > maxeta) {
      return false;
    }

    if (track.tpcChi2NCl() > maxchi2tpc) {
      return false;
    }

    if (track.itsChi2NCl() > maxchi2its) {
      return false;
    }

    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.itsNCls() < minitsncls) {
      return false;
    }

    auto hits = std::count_if(itsRequirement.second.begin(), itsRequirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
    if (hits < itsRequirement.first) {
      return false;
    }

    if (track.tpcNClsFound() < min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < min_tpc_cr_findable_ratio) {
      return false;
    }

    if (track.tpcNSigmaEl() < minTPCNsigmaEl || maxTPCNsigmaEl < track.tpcNSigmaEl()) {
      return false;
    }
    if (minTPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < maxTPCNsigmaPi) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  void fillTrackHistograms(TTrack const& track)
  {
    if (track.sign() > 0) {
      fRegistry.fill(HIST("Track/positive/hPt"), track.pt());
      fRegistry.fill(HIST("Track/positive/hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/positive/hEtaPhi"), track.phi(), track.eta());
      fRegistry.fill(HIST("Track/positive/hPtPhi"), track.phi(), track.pt());
      fRegistry.fill(HIST("Track/positive/hDCAxyz"), track.dcaXY(), track.dcaZ());
      fRegistry.fill(HIST("Track/positive/hDCAxyzSigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
      fRegistry.fill(HIST("Track/positive/hDCAxyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("Track/positive/hDCAzRes_Pt"), track.pt(), sqrt(track.cZZ()) * 1e+4);  // convert cm to um
      fRegistry.fill(HIST("Track/positive/hNclsITS"), track.itsNCls());
      fRegistry.fill(HIST("Track/positive/hNclsTPC"), track.tpcNClsFound());
      fRegistry.fill(HIST("Track/positive/hNcrTPC"), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("Track/positive/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
      fRegistry.fill(HIST("Track/positive/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
      fRegistry.fill(HIST("Track/positive/hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("Track/positive/hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/positive/hITSClusterMap"), track.itsClusterMap());
      fRegistry.fill(HIST("Track/positive/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
      fRegistry.fill(HIST("Track/positive/hTOFbeta"), track.tpcInnerParam(), track.beta());
      fRegistry.fill(HIST("Track/positive/h1OverTOFbeta"), track.tpcInnerParam(), 1.f / track.beta());
      fRegistry.fill(HIST("Track/positive/hTOFNsigmaEl"), track.tpcInnerParam(), track.tofNSigmaEl());
      fRegistry.fill(HIST("Track/positive/hTOFNsigmaMu"), track.tpcInnerParam(), track.tofNSigmaMu());
      fRegistry.fill(HIST("Track/positive/hTOFNsigmaPi"), track.tpcInnerParam(), track.tofNSigmaPi());
      fRegistry.fill(HIST("Track/positive/hTOFNsigmaKa"), track.tpcInnerParam(), track.tofNSigmaKa());
      fRegistry.fill(HIST("Track/positive/hTOFNsigmaPr"), track.tpcInnerParam(), track.tofNSigmaPr());
    } else {
      fRegistry.fill(HIST("Track/negative/hPt"), track.pt());
      fRegistry.fill(HIST("Track/negative/hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/negative/hEtaPhi"), track.phi(), track.eta());
      fRegistry.fill(HIST("Track/negative/hPtPhi"), track.phi(), track.pt());
      fRegistry.fill(HIST("Track/negative/hDCAxyz"), track.dcaXY(), track.dcaZ());
      fRegistry.fill(HIST("Track/negative/hDCAxyzSigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
      fRegistry.fill(HIST("Track/negative/hDCAxyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("Track/negative/hDCAzRes_Pt"), track.pt(), sqrt(track.cZZ()) * 1e+4);  // convert cm to um
      fRegistry.fill(HIST("Track/negative/hNclsITS"), track.itsNCls());
      fRegistry.fill(HIST("Track/negative/hNclsTPC"), track.tpcNClsFound());
      fRegistry.fill(HIST("Track/negative/hNcrTPC"), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("Track/negative/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
      fRegistry.fill(HIST("Track/negative/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
      fRegistry.fill(HIST("Track/negative/hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("Track/negative/hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/negative/hITSClusterMap"), track.itsClusterMap());
      fRegistry.fill(HIST("Track/negative/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
      fRegistry.fill(HIST("Track/negative/hTOFbeta"), track.tpcInnerParam(), track.beta());
      fRegistry.fill(HIST("Track/negative/h1OverTOFbeta"), track.tpcInnerParam(), 1.f / track.beta());
      fRegistry.fill(HIST("Track/negative/hTOFNsigmaEl"), track.tpcInnerParam(), track.tofNSigmaEl());
      fRegistry.fill(HIST("Track/negative/hTOFNsigmaMu"), track.tpcInnerParam(), track.tofNSigmaMu());
      fRegistry.fill(HIST("Track/negative/hTOFNsigmaPi"), track.tpcInnerParam(), track.tofNSigmaPi());
      fRegistry.fill(HIST("Track/negative/hTOFNsigmaKa"), track.tpcInnerParam(), track.tofNSigmaKa());
      fRegistry.fill(HIST("Track/negative/hTOFNsigmaPr"), track.tpcInnerParam(), track.tofNSigmaPr());
    }
  }

  template <typename TTrack>
  void fillTrackTable(TTrack const& track)
  {
    emprimaryelectrons(track.collisionId(), track.globalIndex(), track.sign(),
                       track.pt(), track.eta(), track.phi(), track.dcaXY(), track.dcaZ(),
                       track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(), track.tpcNClsShared(),
                       track.tpcChi2NCl(), track.tpcInnerParam(),
                       track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                       track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                       track.itsClusterSizes(), track.itsChi2NCl(), track.detectorMap(), track.x(), track.alpha(), track.y(), track.z(), track.snp(), track.tgl(), true);
  }

  SliceCache cache;
  Preslice<aod::Tracks> perCollision_track = o2::aod::track::collisionId;

  Filter trackFilter = o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& nabs(o2::aod::track::dcaXY) < dca_xy_max&& nabs(o2::aod::track::dcaZ) < dca_z_max&& o2::aod::track::tpcChi2NCl < maxchi2tpc&& o2::aod::track::itsChi2NCl < maxchi2its&& ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  Filter pidFilter = (minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl) && (o2::aod::pidtpc::tpcNSigmaPi < minTPCNsigmaPi || maxTPCNsigmaPi < o2::aod::pidtpc::tpcNSigmaPi);
  using MyFilteredTracks = soa::Filtered<MyTracks>;

  // ---------- for data ----------
  void processRec(MyCollisions_Cent const& collisions, MyBCs const&, MyFilteredTracks const& tracks)
  {
    for (auto& collision : collisions) {
      auto bc = collision.bc_as<MyBCs>();
      initCCDB(bc);

      fRegistry.fill(HIST("hCollisionCounter"), 0.f);
      if (dist01(engine) > down_scaling) {
        continue;
      }
      fRegistry.fill(HIST("hCollisionCounter"), 1.f);

      event(collision.globalIndex(), bc.runNumber(), bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(),
            collision.posX(), collision.posY(), collision.posZ(),
            collision.numContrib(), collision.trackOccupancyInTimeRange());
      event_mult(collision.multFT0A(), collision.multFT0C(), collision.multNTracksPV(), collision.multNTracksPVeta1(), collision.multNTracksPVetaHalf());
      event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C());

      auto tracks_per_collision = tracks.sliceBy(perCollision_track, collision.globalIndex());

      for (auto& track : tracks_per_collision) {
        if (!checkTrack<false>(track)) {
          continue;
        }
        fillTrackTable(track);
        prmeleventid(event.lastIndex());
        fillTrackHistograms(track);
      } // end of track loop

    } // end of collision loop
  }
  PROCESS_SWITCH(TreeCreatorSingleElectronQA, processRec, "process reconstructed info only", true); // standalone

  // ---------- for MC ----------
  using MyFilteredTracksMC = soa::Filtered<MyTracksMC>;
  void processMC(MyCollisionsMC_Cent const& collisions, aod::McCollisions const&, MyBCs const&, MyFilteredTracksMC const& tracks)
  {
    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.bc_as<MyBCs>();
      initCCDB(bc);

      fRegistry.fill(HIST("hCollisionCounter"), 0.f);
      if (dist01(engine) > down_scaling) {
        continue;
      }
      fRegistry.fill(HIST("hCollisionCounter"), 1.f);

      event(collision.globalIndex(), bc.runNumber(), bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(),
            collision.posX(), collision.posY(), collision.posZ(),
            collision.numContrib(), collision.trackOccupancyInTimeRange());
      event_mult(collision.multFT0A(), collision.multFT0C(), collision.multNTracksPV(), collision.multNTracksPVeta1(), collision.multNTracksPVetaHalf());
      event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C());

      auto tracks_per_collision = tracks.sliceBy(perCollision_track, collision.globalIndex());

      for (auto& track : tracks_per_collision) {
        if (!checkTrack<true>(track)) {
          continue;
        }
        fillTrackTable(track);
        prmeleventid(event.lastIndex());
        fillTrackHistograms(track);
      } // end of track loop

    } // end of collision loop
  }
  PROCESS_SWITCH(TreeCreatorSingleElectronQA, processMC, "process reconstructed and MC info ", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TreeCreatorSingleElectronQA>(cfgc, TaskName{"tree-creator-single-electron-qa"}),
  };
}
