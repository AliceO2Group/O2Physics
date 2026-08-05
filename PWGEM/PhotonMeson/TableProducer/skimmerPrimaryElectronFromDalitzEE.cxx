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

/// \file skimmerPrimaryElectronFromDalitzEE.cxx
/// \brief write relevant information about primary electrons.
/// \author daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Concepts.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <PID/PIDTOFParamService.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/PID.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels>;
using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerBitsTMP>;

using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov,
                           aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                           aod::pidTOFbeta, aod::TOFSignal, aod::TOFEvTime>;
using MyTrack = MyTracks::iterator;
using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels, aod::mcTPCTuneOnData>;
using MyTrackMC = MyTracksMC::iterator;

namespace o2::aod
{
namespace pwgem::pm::recalculatedtofpid
{
DECLARE_SOA_COLUMN(BetaRecalculated, betaRecalculated, float);
DECLARE_SOA_COLUMN(TOFNSigmaElRecalculated, tofNSigmaElRecalculated, float);
} // namespace pwgem::pm::recalculatedtofpid

DECLARE_SOA_TABLE(EMTOFNSigmas, "AOD", "EMTOFNSIGMA", // make std::map in your tasks later. // Don't store this table in the derived data.
                  o2::aod::emprimaryelectron::CollisionId, o2::aod::emprimaryelectron::TrackId,
                  o2::aod::pwgem::pm::recalculatedtofpid::BetaRecalculated, o2::aod::pwgem::pm::recalculatedtofpid::TOFNSigmaElRecalculated);

using EMTOFNSigma = EMTOFNSigmas::iterator;
} // namespace o2::aod

struct skimmerPrimaryElectronFromDalitzEE {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = o2::aod::track::collisionId;
  PresliceOptional<MyTracksMC> perTracksCollision = aod::track::collisionId;
  Preslice<aod::V0PhotonsKF> perCol_pcm = o2::aod::v0photonkf::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Produces<aod::EMPrimaryElectronsFromDalitz> emprimaryelectrons;
  Produces<aod::EMPrimaryElectronsDeDxMC> emprimaryelectronsDeDxMC;
  Service<o2::pid::tof::TOFResponse> mTOFResponse{};

  Produces<aod::EMTOFNSigmas> emtofs;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};

  // Operation and minimisation criteria
  Configurable<float> dBzInput{"dBzInput", -999, "bz field in kG, -999 is automatic"};
  Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 0, "min ncluster tpc"}; // o2-linter: disable=name/function-variable (renaming configs would mess up hyperloop)
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  Configurable<float> max_frac_shared_clusters_tpc{"max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
  Configurable<int> min_ncluster_its{"min_ncluster_its", 4, "min ncluster its"};
  Configurable<int> min_ncluster_itsib{"min_ncluster_itsib", 1, "min ncluster itsib"};
  Configurable<float> minchi2tpc{"minchi2tpc", 0.0, "min. chi2/NclsTPC"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max. chi2/NclsTPC"};
  Configurable<float> minchi2its{"minchi2its", 0.0, "min. chi2/NclsITS"};
  Configurable<float> maxchi2its{"maxchi2its", 36.0, "max. chi2/NclsITS"};
  Configurable<float> minpt{"minpt", 0.05, "min pt for ITS-TPC track"};
  Configurable<float> maxeta{"maxeta", 2.0, "max eta acceptance"};
  Configurable<float> dca_xy_max{"dca_xy_max", 1, "max DCAxy in cm"};                 // o2-linter: disable=name/function-variable (renaming configs would mess up hyperloop)
  Configurable<float> dca_z_max{"dca_z_max", 1, "max DCAz in cm"};                    // o2-linter: disable=name/function-variable (renaming configs would mess up hyperloop)
  Configurable<float> dca_3d_sigma_max{"dca_3d_sigma_max", 2, "max DCA 3D in sigma"}; // o2-linter: disable=name/function-variable (renaming configs would mess up hyperloop)
  Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -2.5, "min. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", +3.5, "max. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaPi{"maxTPCNsigmaPi", 0.0, "max. TPC n sigma for pion exclusion"};
  Configurable<float> minTPCNsigmaPi{"minTPCNsigmaPi", 0.0, "min. TPC n sigma for pion exclusion"};
  Configurable<float> minTOFNsigmaEl{"minTOFNsigmaEl", -3.5, "min. TOF n sigma for electron inclusion"};
  Configurable<float> maxTOFNsigmaEl{"maxTOFNsigmaEl", +3.5, "max. TOF n sigma for electron inclusion"};
  Configurable<float> minTPCNsigmaKa{"minTPCNsigmaKa", -2.5, "min. TPC n sigma for kaon exclusion"};
  Configurable<float> maxTPCNsigmaKa{"maxTPCNsigmaKa", 2.5, "max. TPC n sigma for kaon exclusion"};
  Configurable<float> minTPCNsigmaPr{"minTPCNsigmaPr", -2.5, "min. TPC n sigma for proton exclusion"};
  Configurable<float> maxTPCNsigmaPr{"maxTPCNsigmaPr", 2.5, "max. TPC n sigma for proton exclusion"};
  Configurable<bool> requireTOF{"requireTOF", false, "require TOF hit"};
  Configurable<float> min_pin_for_pion_rejection{"min_pin_for_pion_rejection", 0.0, "pion rejection is applied above this pin"}; // this is used only in TOFreq
  Configurable<float> max_pin_for_pion_rejection{"max_pin_for_pion_rejection", 0.5, "pion rejection is applied below this pin"};
  Configurable<float> maxMee{"maxMee", 0.04, "max. mee to store dalitz ee pairs"};
  Configurable<bool> fillLS{"fillLS", true, "flag to fill LS histograms for QA"};
  Configurable<bool> fillWithPairs{"fillWithPairs", false, "flag to fill table based on pair information"};
  Configurable<bool> includeITSsa{"includeITSsa", false, "Flag to include ITSsa tracks"};
  Configurable<float> maxpt_itssa{"maxpt_itssa", 0.15, "max pt for ITSsa track"}; // o2-linter: disable=name/function-variable (renaming configs would mess up hyperloop)
  Configurable<float> maxMeanITSClusterSize{"maxMeanITSClusterSize", 16, "max <ITS cluster size> x cos(lambda)"};
  Configurable<float> slope{"slope", 0.0185, "slope for m vs. phiv"};
  Configurable<float> intercept{"intercept", -0.0380, "intercept for m vs. phiv"};
  Configurable<bool> useTOFNSigmaDeltaBC{"useTOFNSigmaDeltaBC", false, "Flag to shift delta BC for TOF n sigma (only with TTCA)"};
  Configurable<bool> storeOnlyTrueElectronMC{"storeOnlyTrueElectronMC", false, "Flag to store only true electron in MC"};

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::array<std::string_view, 3> DileptonSigns = {"uls/", "lspp/", "lsmm/"};

  int mRunNumber = 0;
  float dBz = 0.;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::dataformats::VertexBase mVtx;

  void init(InitContext& initContext)
  {
    // if (doprocessRec && doprocessRec_SWT) {
    //   LOGF(fatal, "Cannot enable doprocessRec and doprocessRec_SWT at the same time. Please choose one.");
    // }

    mRunNumber = 0;
    dBz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    LOGF(info, "before TOF initSetup");
    mTOFResponse->initSetup(ccdb, initContext);
    LOGF(info, "after TOF initSetup");

    fRegistry.add("Track/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, o2::constants::math::TwoPI}, {400, -2.0f, 2.0f}}, false);
    fRegistry.add("Track/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/hRelDeltaPt", "pT resolution;p_{T} (GeV/c);#Deltap_{T}/p_{T}", kTH2F, {{1000, 0, 10}, {100, 0, 0.1}}, false);
    fRegistry.add("Track/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/hDCAxy_Pt", "DCA_{xy} vs. pT;p_{T} (GeV/c);DCA_{xy} (cm)", kTH2F, {{200, 0, 10}, {200, -1, 1}}, false);
    fRegistry.add("Track/hDCAz_Pt", "DCA_{z} vs. pT;p_{T} (GeV/c);DCA_{z} (cm)", kTH2F, {{200, 0, 10}, {200, -1, 1}}, false);
    fRegistry.add("Track/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{200, 0, 10}, {500, 0., 500}}, false);
    fRegistry.add("Track/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{200, 0, 10}, {500, 0., 500}}, false);

    // TPC
    fRegistry.add("Track/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
    fRegistry.add("Track/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/hTPCdEdxMC", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);

    // ITS
    fRegistry.add("Track/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("Track/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{400, 0, 40}}, false);
    fRegistry.add("Track/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("Track/hMeanClusterSizeITS", "mean cluster size ITS;p_{pv} (GeV/c);<ITS cluster size> #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {150, 0, 15}}, false);
    fRegistry.add("Track/hMeanClusterSizeITSib", "mean cluster size ITSib;p_{pv} (GeV/c);<ITSib cluster size> #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {150, 0, 15}}, false);
    fRegistry.add("Track/hMeanClusterSizeITSob", "mean cluster size ITSob;p_{pv} (GeV/c);<ITSob cluster size> #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {150, 0, 15}}, false);

    // TOF
    fRegistry.add("Track/hChi2TOF", "chi2 of TOF", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/hTOFbeta", "TOF beta;p_{pv} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);
    fRegistry.add("Track/hTOFNsigmaEl", "TOF n sigma el;p_{pv} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTOFNsigmaPi", "TOF n sigma pi;p_{pv} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);

    // pair
    fRegistry.add("Pair/uls/hTrackMvsPt", "m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {{100, 0, 0.1}, {200, 0, 2}}, false);
    fRegistry.add("Pair/uls/hCheckEMvsPt", "m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {{100, 0, 0.1}, {200, 0, 2}}, false);
    fRegistry.add("Pair/uls/hMvsPt", "m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {{100, 0, 0.1}, {200, 0, 2}}, false);
    fRegistry.add("Pair/uls/hMCutMvsPt", "m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {{100, 0, 0.1}, {200, 0, 2}}, false);
    fRegistry.add("Pair/uls/hMPhiCutMvsPt", "m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {{100, 0, 0.1}, {200, 0, 2}}, false);

    fRegistry.add("Pair/uls/hTrackMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{180, 0, o2::constants::math::PI}, {100, 0, 0.1}}, false);
    fRegistry.add("Pair/uls/hCheckEMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{180, 0, o2::constants::math::PI}, {100, 0, 0.1}}, false);
    fRegistry.add("Pair/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{180, 0, o2::constants::math::PI}, {100, 0, 0.1}}, false);
    fRegistry.add("Pair/uls/hMCutMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{180, 0, o2::constants::math::PI}, {100, 0, 0.1}}, false);
    fRegistry.add("Pair/uls/hMPhiCutMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{180, 0, o2::constants::math::PI}, {100, 0, 0.1}}, false);

    fRegistry.addClone("Pair/uls/", "Pair/lspp/");
    fRegistry.addClone("Pair/uls/", "Pair/lsmm/");
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (dBzInput > -990) { // o2-linter: disable=magic-number (check against some default number)
      dBz = dBzInput;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(dBz) > 1e-5) {                   // o2-linter: disable=magic-number (check against some default number)
        grpmag.setL3Current(30000.f / (dBz / 5.0f)); // o2-linter: disable=magic-number (values to calculate the magnetic field)
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grpTimestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = nullptr;
    o2::parameters::GRPMagField* grpmag = nullptr;
    if (!skipGRPOquery) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grpTimestamp);
    }
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      dBz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grpTimestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grpTimestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      dBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
    }
    mRunNumber = bc.runNumber();
  }

  template <bool withTTCA, o2::soa::is_table TCollisions, o2::soa::is_table TBCs, o2::soa::is_table TTracks, o2::soa::is_table TTrackAssoc>
  void calculateTOFNSigmaWithReassociation(TCollisions const& collisions, TBCs const&, TTracks const& tracks, TTrackAssoc const& trackIndices)
  {
    if (useTOFNSigmaDeltaBC) {
      if constexpr (withTTCA) {
        for (const auto& collision : collisions) {
          if (mapCollisionTime.find(collision.globalIndex()) == mapCollisionTime.end()) {
            continue;
          }
          auto bcCollision = collision.template bc_as<TBCs>();
          auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
          for (const auto& trackId : trackIdsThisCollision) {
            auto track = trackId.template track_as<TTracks>();
            if (!track.hasITS() || !track.hasTPC()) { // apply only minimal cut
              continue;
            }

            if (track.hasTOF() && track.has_collision()) { // TTCA may use orphan tracks.
              auto bcTrack = track.template collision_as<TCollisions>().template bc_as<TBCs>();
              float tofNSigmaEl = mTOFResponse->nSigma<o2::track::PID::Electron>(track.tofSignalInAnotherBC(bcTrack.globalBC(), bcCollision.globalBC()), track.tofExpMom(), track.length(), track.p(), track.eta(), mapCollisionTime[collision.globalIndex()], mapCollisionTimeError[collision.globalIndex()]);
              float beta = track.length() / (track.tofSignalInAnotherBC(bcTrack.globalBC(), bcCollision.globalBC()) - mapCollisionTime[collision.globalIndex()]) / (o2::constants::physics::LightSpeedCm2S * 1e-12);
              mapTOFNsigmaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = tofNSigmaEl;
              mapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = beta;
              emtofs(collision.globalIndex(), track.globalIndex(), mapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())], mapTOFNsigmaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())]);
            } else {
              mapTOFNsigmaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaEl();
              mapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.beta();
              emtofs(collision.globalIndex(), track.globalIndex(), mapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())], mapTOFNsigmaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())]);
            }
          } // end of track loop
        } // end of collision loop
      } else {
        for (const auto& collision : collisions) {
          auto tracks_per_coll = tracks.sliceBy(perCol, collision.globalIndex());
          for (const auto& track : tracks_per_coll) {
            if (!track.hasITS() || !track.hasTPC()) { // apply only minimal cut
              continue;
            }
            mapTOFNsigmaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaEl();
            mapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.beta();
            emtofs(collision.globalIndex(), track.globalIndex(), mapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())], mapTOFNsigmaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())]);
          }
        } // end of track loop
      } // end of collision loop
    } else {
      if constexpr (withTTCA) {
        for (const auto& collision : collisions) {
          auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
          for (const auto& trackId : trackIdsThisCollision) {
            auto track = trackId.template track_as<TTracks>();
            if (!track.hasITS() || !track.hasTPC()) { // apply only minimal cut
              continue;
            }
            mapTOFNsigmaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaEl();
            mapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.beta();
            emtofs(collision.globalIndex(), track.globalIndex(), mapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())], mapTOFNsigmaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())]);
          } // end of track loop
        } // end of collision loop
      } else {
        for (const auto& collision : collisions) {
          auto tracks_per_coll = tracks.sliceBy(perCol, collision.globalIndex());
          for (const auto& track : tracks_per_coll) {
            if (!track.hasITS() || !track.hasTPC()) { // apply only minimal cut
              continue;
            }
            mapTOFNsigmaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaEl();
            mapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.beta();
            emtofs(collision.globalIndex(), track.globalIndex(), mapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())], mapTOFNsigmaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())]);
          }
        } // end of track loop
      } // end of collision loop
    }
  }

  template <bool isMC, o2::soa::is_iterator TCollision, o2::soa::is_iterator TTrack>
  bool checkTrack(TCollision const& collision, TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
      if (storeOnlyTrueElectronMC) {
        const auto& mcParticle = track.template mcParticle_as<aod::McParticles>();
        if (std::abs(mcParticle.pdgCode()) != 11) {
          return false;
        }
      }
    }

    float tofNSigmaEl = mapTOFNsigmaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())];
    if (requireTOF && !(track.hasTOF() && std::fabs(tofNSigmaEl) < maxTOFNsigmaEl)) {
      return false;
    }

    if (!track.hasITS()) {
      return false;
    }

    if (track.itsChi2NCl() < minchi2its || maxchi2its < track.itsChi2NCl()) {
      return false;
    }

    if (track.itsNCls() < min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < min_ncluster_itsib) {
      return false;
    }

    if (!includeITSsa && !track.hasTPC()) {
      return false;
    }

    if (track.hasTPC()) {
      if (track.tpcChi2NCl() < minchi2tpc || maxchi2tpc < track.tpcChi2NCl()) {
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

      if (track.tpcFractionSharedCls() > max_frac_shared_clusters_tpc) {
        return false;
      }
    }

    o2::dataformats::DCA mDcaInfoCov; // for Track association
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    trackParCov.setPID(o2::track::PID::Electron);
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    if (!isPropOK) {
      return false;
    }
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();

    if (std::fabs(dcaXY) > dca_xy_max || std::fabs(dcaZ) > dca_z_max) {
      return false;
    }

    float dca3D = 999.f;
    float det = trackParCov.getSigmaY2() * trackParCov.getSigmaZ2() - trackParCov.getSigmaZY() * trackParCov.getSigmaZY();
    if (det < 0) {
      dca3D = 999.f;
    } else {
      float chi2 = (dcaXY * dcaXY * trackParCov.getSigmaZ2() + dcaZ * dcaZ * trackParCov.getSigmaY2() - 2. * dcaXY * dcaZ * trackParCov.getSigmaZY()) / det;
      dca3D = std::sqrt(std::fabs(chi2) / 2.);
    }
    if (dca3D > dca_3d_sigma_max) {
      return false;
    }

    if (trackParCov.getPt() < minpt || std::fabs(trackParCov.getEta()) > maxeta) {
      return false;
    }

    int total_cluster_size = 0, nl = 0;
    for (unsigned int layer = 0; layer < 7; layer++) {
      int cluster_size_per_layer = track.itsClsSizeInLayer(layer);
      if (cluster_size_per_layer > 0) {
        nl++;
      }
      total_cluster_size += cluster_size_per_layer;
    }

    if (maxMeanITSClusterSize < static_cast<float>(total_cluster_size) / static_cast<float>(nl) * std::cos(std::atan(trackParCov.getTgl()))) {
      return false;
    }

    if ((track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF()) && maxpt_itssa < track.pt()) {
      return false;
    }

    // if ((track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF()) && maxpt_itssa < track.pt()) {
    //   return false;
    // }

    return true;
  }

  template <o2::soa::is_iterator TCollision, o2::soa::is_iterator TTrack>
  bool isElectron(TCollision const& collision, TTrack const& track)
  {
    if (includeITSsa && (track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF())) [[unlikely]] {
      int total_cluster_size = 0, nl = 0;
      for (unsigned int layer = 0; layer < 7; layer++) {
        int cluster_size_per_layer = track.itsClsSizeInLayer(layer);
        if (cluster_size_per_layer > 0) {
          nl++;
        }
        total_cluster_size += cluster_size_per_layer;
      }

      return (maxMeanITSClusterSize > static_cast<float>(total_cluster_size) / static_cast<float>(nl) * std::cos(std::atan(track.tgl())));
    }
    return isElectron_TPChadrej(collision, track) || isElectron_TOFreq(collision, track);
  }

  template <o2::soa::is_iterator TCollision, o2::soa::is_iterator TTrack>
  bool isElectron_TPChadrej(TCollision const& collision, TTrack const& track)
  {
    float tofNSigmaEl = mapTOFNsigmaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())];

    if (track.tpcNSigmaEl() < minTPCNsigmaEl || maxTPCNsigmaEl < track.tpcNSigmaEl()) {
      return false;
    }
    if (minTPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < maxTPCNsigmaPi && track.tpcInnerParam() < max_pin_for_pion_rejection) {
      return false;
    }
    if (minTPCNsigmaKa < track.tpcNSigmaKa() && track.tpcNSigmaKa() < maxTPCNsigmaKa) {
      return false;
    }
    if (minTPCNsigmaPr < track.tpcNSigmaPr() && track.tpcNSigmaPr() < maxTPCNsigmaPr) {
      return false;
    }
    if (track.hasTOF() && (maxTOFNsigmaEl < std::fabs(tofNSigmaEl))) {
      return false;
    }
    return true;
  }

  template <o2::soa::is_iterator TCollision, o2::soa::is_iterator TTrack>
  bool isElectron_TOFreq(TCollision const& collision, TTrack const& track)
  {
    float tofNSigmaEl = mapTOFNsigmaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())];

    if (minTPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < maxTPCNsigmaPi && (min_pin_for_pion_rejection < track.tpcInnerParam() && track.tpcInnerParam() < max_pin_for_pion_rejection)) {
      return false;
    }
    return minTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < maxTPCNsigmaEl && std::fabs(tofNSigmaEl) < maxTOFNsigmaEl;
  }

  template <bool isMC, o2::soa::is_iterator TCollision, o2::soa::is_table TTracks>
  void fillTrackInfo(TCollision const& collision, TTracks const& tracks)
  {

    for (const auto& track : tracks) {
      if (!checkTrack<isMC>(collision, track)) {
        continue;
      }

      if (!isElectron(collision, track)) {
        continue;
      }

      fillTrackHistograms<isMC>(track);
      fillTrackTable<isMC>(collision, track);
    }
  }

  template <bool isMC, int pairtype, o2::soa::is_iterator TCollision, o2::soa::is_table TTracks1, o2::soa::is_table TTracks2>
  void fillPairInfo(TCollision const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    if constexpr (pairtype == 0) { // ULS
      for (const auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!checkTrack<isMC>(collision, t1) || !checkTrack<isMC>(collision, t2)) {
          continue;
        }

        if (!isElectron(collision, t1) || !isElectron(collision, t2)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), dBz);

        fRegistry.fill(HIST("Pair/") + HIST(DileptonSigns[pairtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/") + HIST(DileptonSigns[pairtype]) + HIST("hMvsPhiV"), phiv, v12.M());

        if (v12.M() > maxMee) { // don't store
          continue;
        }

        fRegistry.fill(HIST("Pair/") + HIST(DileptonSigns[pairtype]) + HIST("hMCutMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/") + HIST(DileptonSigns[pairtype]) + HIST("hMCutMvsPhiV"), phiv, v12.M());

        if (v12.M() < slope * phiv + intercept) {
          continue;
        }

        fRegistry.fill(HIST("Pair/") + HIST(DileptonSigns[pairtype]) + HIST("hMPhiCutMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/") + HIST(DileptonSigns[pairtype]) + HIST("hMPhiCutMvsPhiV"), phiv, v12.M());

        if (t1.sign() > 0) { // for positron
          if (std::find(acceptedPosTrackIds_per_collision.begin(), acceptedPosTrackIds_per_collision.end(), t1.globalIndex()) == acceptedPosTrackIds_per_collision.end()) {
            fillTrackHistograms<isMC>(t1);
            acceptedPosTrackIds_per_collision.emplace_back(t1.globalIndex());
          }
        } else { // for electron
          if (std::find(acceptedNegTrackIds_per_collision.begin(), acceptedNegTrackIds_per_collision.end(), t1.globalIndex()) == acceptedNegTrackIds_per_collision.end()) {
            fillTrackHistograms<isMC>(t1);
            acceptedNegTrackIds_per_collision.emplace_back(t1.globalIndex());
          }
        }

        if (t2.sign() > 0) { // for positron
          if (std::find(acceptedPosTrackIds_per_collision.begin(), acceptedPosTrackIds_per_collision.end(), t2.globalIndex()) == acceptedPosTrackIds_per_collision.end()) {
            fillTrackHistograms<isMC>(t2);
            acceptedPosTrackIds_per_collision.emplace_back(t2.globalIndex());
          }
        } else { // for electron
          if (std::find(acceptedNegTrackIds_per_collision.begin(), acceptedNegTrackIds_per_collision.end(), t2.globalIndex()) == acceptedNegTrackIds_per_collision.end()) {
            fillTrackHistograms<isMC>(t2);
            acceptedNegTrackIds_per_collision.emplace_back(t2.globalIndex());
          }
        }
      } // end of ULS pairing
    } else { // LS
      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks1, tracks2))) {
        if (!checkTrack<isMC>(collision, t1) || !checkTrack<isMC>(collision, t2)) {
          continue;
        }
        if (!isElectron(collision, t1) || !isElectron(collision, t2)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), dBz);
        fRegistry.fill(HIST("Pair/") + HIST(DileptonSigns[pairtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/") + HIST(DileptonSigns[pairtype]) + HIST("hMvsPhiV"), phiv, v12.M());
      } // end of LS pairing
    }
  }

  template <bool isMC, o2::soa::is_iterator TTrack>
  void fillTrackHistograms(TTrack const& track)
  {
    float mcTunedTPCSignal = 0.f;
    if constexpr (isMC) {
      mcTunedTPCSignal = track.mcTunedTPCSignal();
    }

    fRegistry.fill(HIST("Track/hPt"), track.pt());
    fRegistry.fill(HIST("Track/hEtaPhi"), track.phi(), track.eta());
    fRegistry.fill(HIST("Track/hQoverPt"), track.sign() / track.pt());
    fRegistry.fill(HIST("Track/hRelDeltaPt"), track.pt(), track.sigma1Pt() * track.pt());
    fRegistry.fill(HIST("Track/hDCAxyz"), track.dcaXY(), track.dcaZ());
    fRegistry.fill(HIST("Track/hDCAxy_Pt"), track.pt(), track.dcaXY());
    fRegistry.fill(HIST("Track/hDCAz_Pt"), track.pt(), track.dcaZ());
    fRegistry.fill(HIST("Track/hDCAxyzSigma"), track.dcaXY() / std::sqrt(track.cYY()), track.dcaZ() / std::sqrt(track.cZZ()));
    fRegistry.fill(HIST("Track/hDCAxyRes_Pt"), track.pt(), std::sqrt(track.cYY()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Track/hDCAzRes_Pt"), track.pt(), std::sqrt(track.cZZ()) * 1e+4);  // convert cm to um

    fRegistry.fill(HIST("Track/hNclsTPC"), track.tpcNClsFound());
    fRegistry.fill(HIST("Track/hNcrTPC"), track.tpcNClsCrossedRows());
    fRegistry.fill(HIST("Track/hChi2TPC"), track.tpcChi2NCl());
    fRegistry.fill(HIST("Track/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
    fRegistry.fill(HIST("Track/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
    fRegistry.fill(HIST("Track/hTPCNclsShared"), track.pt(), track.tpcFractionSharedCls());
    fRegistry.fill(HIST("Track/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
    fRegistry.fill(HIST("Track/hTPCdEdxMC"), track.tpcInnerParam(), mcTunedTPCSignal);
    fRegistry.fill(HIST("Track/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
    fRegistry.fill(HIST("Track/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());

    fRegistry.fill(HIST("Track/hChi2TOF"), track.tofChi2());
    fRegistry.fill(HIST("Track/hTOFbeta"), track.p(), track.beta());
    fRegistry.fill(HIST("Track/hTOFNsigmaEl"), track.p(), track.tofNSigmaEl());
    fRegistry.fill(HIST("Track/hTOFNsigmaPi"), track.p(), track.tofNSigmaPi());

    fRegistry.fill(HIST("Track/hNclsITS"), track.itsNCls());
    fRegistry.fill(HIST("Track/hChi2ITS"), track.itsChi2NCl());
    fRegistry.fill(HIST("Track/hITSClusterMap"), track.itsClusterMap());

    int totalClusterSize = 0, nl = 0;
    for (unsigned int layer = 0; layer < 7; layer++) { // o2-linter: disable=magic-number (number of ITS layers)
      int clusterSizePerLayer = track.itsClsSizeInLayer(layer);
      if (clusterSizePerLayer > 0) {
        nl++;
      }
      totalClusterSize += clusterSizePerLayer;
    }

    int totalClusterSizeIb = 0, nl_ib = 0;
    for (unsigned int layer = 0; layer < 3; layer++) { // o2-linter: disable=magic-number (number of inner barrel layers)
      int clusterSizePerLayer = track.itsClsSizeInLayer(layer);
      if (clusterSizePerLayer > 0) {
        nl_ib++;
      }
      totalClusterSizeIb += clusterSizePerLayer;
    }

    int totalClusterSizeOb = 0, nlOb = 0;
    for (unsigned int layer = 3; layer < 7; layer++) { // o2-linter: disable=magic-number (number of outer barrel layers)
      int clusterSizePerLayer = track.itsClsSizeInLayer(layer);
      if (clusterSizePerLayer > 0) {
        nlOb++;
      }
      totalClusterSizeOb += clusterSizePerLayer;
    }
    fRegistry.fill(HIST("Track/hMeanClusterSizeITS"), track.p(), static_cast<float>(totalClusterSize) / static_cast<float>(nl) * std::cos(std::atan(track.tgl())));
    fRegistry.fill(HIST("Track/hMeanClusterSizeITSib"), track.p(), static_cast<float>(totalClusterSizeIb) / static_cast<float>(nl_ib) * std::cos(std::atan(track.tgl())));
    fRegistry.fill(HIST("Track/hMeanClusterSizeITSob"), track.p(), static_cast<float>(totalClusterSizeOb) / static_cast<float>(nlOb) * std::cos(std::atan(track.tgl())));
  }

  template <bool isMC, o2::soa::is_iterator TCollision, o2::soa::is_iterator TTrack>
  void fillTrackTable(TCollision const& collision, TTrack const& track)
  {
    emprimaryelectrons(collision.globalIndex(), track.globalIndex(), track.sign(),
                       track.pt(), track.eta(), track.phi(), track.dcaXY(), track.dcaZ(), track.cYY(), track.cZY(), track.cZZ(),
                       track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(), track.tpcNClsShared(),
                       track.tpcChi2NCl(), track.tpcInnerParam(),
                       track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaPi(),
                       track.beta(), track.tofNSigmaEl(),
                       track.itsClusterSizes(), track.itsChi2NCl(), track.tofChi2(), track.detectorMap());

    if constexpr (isMC) {
      emprimaryelectronsDeDxMC(track.mcTunedTPCSignal());
    }
  }

  std::vector<int> acceptedPosTrackIds_per_collision;
  std::vector<int> acceptedNegTrackIds_per_collision;
  std::vector<int> acceptedTrackIds_per_collision;
  // Filter trackFilter = minpt < o2::aod::track::pt && nabs(o2::aod::track::eta) < maxeta && o2::aod::track::itsChi2NCl < maxchi2its && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && nabs(o2::aod::track::dcaXY) < dca_xy_max && nabs(o2::aod::track::dcaZ) < dca_z_max;
  // Filter trackFilter
  // using MyFilteredTracks = soa::Filtered<MyTracks>;
  Partition<MyTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  std::unordered_map<int, double> mapCollisionTime;
  std::unordered_map<int, double> mapCollisionTimeError;

  std::map<std::pair<int, int>, float> mapTOFNsigmaReassociated; // map pair(collisionId, trackId) -> tof n sigma
  std::map<std::pair<int, int>, float> mapTOFBetaReassociated;   // map pair(collisionId, trackId) -> tof beta

  // ---------- for data ----------
  void processRec(MyCollisions const& collisions, aod::BCsWithTimestamps const& bcs, MyTracks const& tracks, aod::V0PhotonsKF const& v0photons, aod::TrackAssoc const& trackIndices)
  {
    initCCDB(bcs.iteratorAt(0));
    mTOFResponse->processSetup(bcs.iteratorAt(0));

    for (const auto& track : tracks) {
      mapCollisionTime.try_emplace(track.collisionId(), track.tofEvTime());
      mapCollisionTimeError.try_emplace(track.collisionId(), track.tofEvTimeErr());
    }
    calculateTOFNSigmaWithReassociation<true>(collisions, bcs, tracks, trackIndices);

    for (const auto& collision : collisions) {
      if (!collision.isSelected()) {
        continue;
      }

      const auto& v0photons_per_coll = v0photons.sliceBy(perCol_pcm, collision.globalIndex());
      const auto& posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      const auto& negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      const auto& slicedTracks = tracks.sliceBy(perTracksCollision, collision.globalIndex());
      acceptedPosTrackIds_per_collision.reserve(posTracks_per_coll.size());
      acceptedNegTrackIds_per_collision.reserve(negTracks_per_coll.size());

      if (!fillWithPairs) {
        fillTrackInfo<false>(collision, slicedTracks);
      } else {
        fillPairInfo<false, 0>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
        if (fillLS) {
          fillPairInfo<false, 1>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
          fillPairInfo<false, 2>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
        }

        if ((v0photons_per_coll.size() >= 1 && !acceptedPosTrackIds_per_collision.empty() && !acceptedNegTrackIds_per_collision.empty()) || (acceptedPosTrackIds_per_collision.size() >= 2 && acceptedNegTrackIds_per_collision.size() >= 2)) {
          for (const auto& posId : acceptedPosTrackIds_per_collision) {
            const auto& pos = tracks.rawIteratorAt(posId);
            fillTrackTable<false>(collision, pos);
          }
          for (const auto& eleId : acceptedNegTrackIds_per_collision) {
            const auto& ele = tracks.rawIteratorAt(eleId);
            fillTrackTable<false>(collision, ele);
          }
        }
      }

      acceptedPosTrackIds_per_collision.clear();
      acceptedNegTrackIds_per_collision.clear();

    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerPrimaryElectronFromDalitzEE, processRec, "process reconstructed info only", false); // standalone

  // void processRec_SWT(MyCollisionsWithSWT const& collisions, aod::BCsWithTimestamps const& bcs, MyTracks const& tracks, aod::V0PhotonsKF const& v0photons, aod::TrackAssoc const& trackIndices)
  // {
  //   initCCDB(bcs.iteratorAt(0));

  //   mTOFResponse->processSetup(bcs.iteratorAt(0));
  //   for (const auto& track : tracks) {
  //     mapCollisionTime.try_emplace(track.collisionId(), track.tofEvTime());
  //     mapCollisionTimeError.try_emplace(track.collisionId(), track.tofEvTimeErr());
  //   }
  //   calculateTOFNSigmaWithReassociation<true>(collisions, bcs, tracks, trackIndices);

  //   for (const auto& collision : collisions) {
  //     if (!collision.isSelected()) {
  //       continue;
  //     }

  //     if (collision.swtaliastmp_raw() == 0) {
  //       continue;
  //     }

  //     const auto& v0photons_per_coll = v0photons.sliceBy(perCol_pcm, collision.globalIndex());
  //     const auto& posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
  //     const auto& negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
  //     const auto& slicedTracks = tracks.sliceBy(perTracksCollision, collision.globalIndex());
  //     acceptedPosTrackIds_per_collision.reserve(posTracks_per_coll.size());
  //     acceptedNegTrackIds_per_collision.reserve(negTracks_per_coll.size());

  //     if(!fillWithPairs){
  //       fillTrackInfo<false>(collision, slicedTracks);
  //     }else{
  //       fillPairInfo<false, 0>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
  //       if (fillLS) {
  //         fillPairInfo<false, 1>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
  //         fillPairInfo<false, 2>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
  //       }

  //       if ((v0photons_per_coll.size() >= 1 &&  !acceptedPosTrackIds_per_collision.empty() && !acceptedNegTrackIds_per_collision.empty()) || (acceptedPosTrackIds_per_collision.size() >= 2 && acceptedNegTrackIds_per_collision.size() >= 2)) {
  //         // LOGF(info, "v0photons_per_coll.size() = %d, acceptedPosTrackIds_per_collision.size() = %d, acceptedNegTrackIds_per_collision.size() = %d", v0photons_per_coll.size(), acceptedPosTrackIds_per_collision.size(), acceptedNegTrackIds_per_collision.size());
  //         for (const auto& posId : acceptedPosTrackIds_per_collision) {
  //           const auto& pos = tracks.rawIteratorAt(posId);
  //           fillTrackTable<false>(collision, pos);
  //         }
  //         for (const auto& eleId : acceptedNegTrackIds_per_collision) {
  //           const auto& ele = tracks.rawIteratorAt(eleId);
  //           fillTrackTable<false>(collision, ele);
  //         }
  //       }
  //     }

  //     acceptedPosTrackIds_per_collision.clear();
  //     acceptedNegTrackIds_per_collision.clear();

  //   } // end of collision loop
  // }
  // PROCESS_SWITCH(skimmerPrimaryElectronFromDalitzEE, processRec_SWT, "process reconstructed info with CEFP", false); // with cefp

  // using MyFilteredTracksMC = soa::Filtered<MyTracksMC>;
  Partition<MyTracksMC> posTracksMC = o2::aod::track::signed1Pt > 0.f;
  Partition<MyTracksMC> negTracksMC = o2::aod::track::signed1Pt < 0.f;
  // ---------- for MC ----------
  void processMC(MyCollisionsMC const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const& bcs, MyTracksMC const& tracks, aod::V0PhotonsKF const& v0photons, aod::TrackAssoc const& trackIndices)
  {
    uint64_t nCollisSel = 0;
    uint64_t nNoMcColl = 0;
    uint64_t nProcessedCollisions = 0;
    uint64_t nCheckTrack = 0;

    initCCDB(bcs.iteratorAt(0));

    mTOFResponse->processSetup(bcs.iteratorAt(0));
    for (const auto& track : tracks) {
      mapCollisionTime.try_emplace(track.collisionId(), track.tofEvTime());
      mapCollisionTimeError.try_emplace(track.collisionId(), track.tofEvTimeErr());
    }
    calculateTOFNSigmaWithReassociation<true>(collisions, bcs, tracks, trackIndices);

    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        nNoMcColl++;
        continue;
      }

      if (!collision.isSelected()) {
        nCollisSel++;
        continue;
      }
      nProcessedCollisions++;

      const auto& v0photons_per_coll = v0photons.sliceBy(perCol_pcm, collision.globalIndex());
      const auto& posTracks_per_coll = posTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      const auto& negTracks_per_coll = negTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      const auto& slicedTracks = tracks.sliceBy(perTracksCollision, collision.globalIndex());
      acceptedPosTrackIds_per_collision.reserve(posTracks_per_coll.size());
      acceptedNegTrackIds_per_collision.reserve(negTracks_per_coll.size());
      acceptedTrackIds_per_collision.reserve(2 * (negTracks_per_coll.size()));

      if (!fillWithPairs) {
        fillTrackInfo<true>(collision, slicedTracks);
      } else {
        fillPairInfo<true, 0>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
        if (fillLS) {
          fillPairInfo<true, 1>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
          fillPairInfo<true, 2>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
        }
        if ((acceptedPosTrackIds_per_collision.empty() && acceptedNegTrackIds_per_collision.empty()) || (acceptedPosTrackIds_per_collision.size() >= 2 && acceptedNegTrackIds_per_collision.size() >= 2)) { // v0photons_per_coll.size() >= 1 &&
          for (const auto& posId : acceptedPosTrackIds_per_collision) {
            const auto& pos = tracks.rawIteratorAt(posId);
            fillTrackTable<true>(collision, pos);
          }
          for (const auto& eleId : acceptedNegTrackIds_per_collision) {
            const auto& ele = tracks.rawIteratorAt(eleId);
            fillTrackTable<true>(collision, ele);
          }
        }
      } // end of fill loop

      acceptedPosTrackIds_per_collision.clear();
      acceptedNegTrackIds_per_collision.clear();

    } // end of collision loop

    LOG(info) << "Total collisions: " << collisions.size();
    LOG(info) << "No MC collision: " << nNoMcColl;
    LOG(info) << "Rejected by selection: " << nCollisSel;
    LOG(info) << "Processed collisions: " << nProcessedCollisions;
    LOG(info) << "Invalid track: " << nCheckTrack;
  }
  PROCESS_SWITCH(skimmerPrimaryElectronFromDalitzEE, processMC, "process reconstructed and MC info ", true);
};

// WorkflowSpec defineDataProcessing(ConfigContext const& context)
// {
//   return WorkflowSpec{adaptAnalysisTask<skimmerPrimaryElectronFromDalitzEE>(context, TaskName{"skimmer-primary-electron-from-dalitzee"})};
// }
WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  o2::pid::tof::TOFResponseImpl::metadataInfo.initMetadata(context);

  return WorkflowSpec{
    adaptAnalysisTask<skimmerPrimaryElectronFromDalitzEE>(context, TaskName{"skimmer-primary-electron-from-dalitzee"})};
}
