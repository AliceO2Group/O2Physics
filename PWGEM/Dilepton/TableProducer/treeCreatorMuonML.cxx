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
// This code will create data table for inputs to machine learning for electrons.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/lmeeMLTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/CollisionTypeHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MCHTracking/TrackExtrap.h>
#include <MathUtils/Utils.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <TGeoGlobalMagField.h>
#include <TH1.h>

#include <cmath>
#include <map>
#include <random>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <math.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using namespace o2::aod::fwdtrackutils;

struct TreeCreatorMuonML {
  using MyCollisionsMC = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::McCollisionLabels>;
  using MyCollisionMC = MyCollisionsMC::iterator;

  using MyFwdTracksMC = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;
  using MyFwdTrackMC = MyFwdTracksMC::iterator;

  using MyMFTTracksMC = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
  using MyMFTTrackMC = MyMFTTracksMC::iterator;

  Produces<aod::EMFwdTracksForML> trackTable;
  Produces<aod::EMFwdTrackErrsForML> trackErrTable;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // for z shift for propagation
  Configurable<bool> cfgApplyZShiftFromCCDB{"cfgApplyZShiftFromCCDB", false, "flag to apply z shift"};
  Configurable<std::string> cfgZShiftPath{"cfgZShiftPath", "Users/m/mcoquet/ZShift", "CCDB path for z shift to apply to forward tracks"};
  Configurable<float> cfgManualZShift{"cfgManualZShift", 0, "manual z-shift for propagation of global muon to PV"};
  Configurable<float> cfgDownSampling{"cfgDownSampling", 1.1, "down sampling for fake matches"};
  Configurable<float> matchingZ{"matchingZ", -77.5, "z position where matching is performed"};
  Configurable<float> maxMatchingChi2MCHMFT{"maxMatchingChi2MCHMFT", 100.f, "max. chi2 for MCH-MFT matching"};

  struct : ConfigurableGroup {
    std::string prefix = "eventCutGroup";
    Configurable<float> cfgMinZvtx{"cfgMinZvtx", -10.f, "min. Zvtx of collision"};
    Configurable<float> cfgMaxZvtx{"cfgMaxZvtx", +10.f, "max. Zvtx of collision"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    // for RCT
    o2::framework::Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", false, "require good detector flag in run condtion table"};
    o2::framework::Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_muon_glo", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
    o2::framework::Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for AA"};
    o2::framework::Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  } eventCutGroup;

  struct : ConfigurableGroup {
    std::string prefix = "glMuonCutGroup";
    // Configurable<float> minEta{"minEta", -3.6, "min. eta acceptance for MFT-MCH-MID"};
    // Configurable<float> maxEta{"maxEta", -2.5, "max. eta acceptance for MFT-MCH-MID"};
    // Configurable<float> maxChi2{"maxChi2", 1e+10, "max. chi2 for global muon"};
    Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "flag to refit global muon"};
  } glMuonCutGroup;

  struct : ConfigurableGroup {
    std::string prefix = "MFTCutGroup";
    Configurable<float> minPt{"minPt", 0.04f, "min. pT for MFTsa to reject crazy tracks"};
    Configurable<float> minEta{"minEta", -4.6f, "min. eta acceptance for MFTsa to reject crazy tracks"};
    Configurable<float> maxEta{"maxEta", -1.5f, "max. eta acceptance for MFTsa to reject crazy tracks"};
  } MFTCutGroup;

  struct : ConfigurableGroup {
    std::string prefix = "MCHMIDCutGroup";
    Configurable<float> minEta{"minEta", -4.6f, "min. eta acceptance for MCHMID to reject crazy tracks"};
    Configurable<float> maxEta{"maxEta", -1.5f, "max. eta acceptance for MCHMID to reject crazy tracks"};
  } MCHMIDCutGroup;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  std::mt19937 engine;
  std::uniform_real_distribution<float> dist01;

  o2::aod::rctsel::RCTFlagsChecker rctChecker;
  ctpRateFetcher mRateFetcher;
  std::string irSourceForCptFetcher{""};

  int mRunNumber = 0;
  float mBz = 0;
  float mZShift = 0;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view muon_types[5] = {"MFTMCHMID/", "MFTMCHMIDOtherMatch/", "MFTMCH/", "MCHMID/", "MCH/"};

  void init(InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);
    rctChecker.init(eventCutGroup.cfgRCTLabel.value, eventCutGroup.cfgCheckZDC.value, eventCutGroup.cfgTreatLimitedAcceptanceAsBad.value);

    addHistograms();

    mRunNumber = 0;
    mBz = 0;
    mZShift = 0;

    std::random_device seed_gen;
    engine = std::mt19937(seed_gen());
    dist01 = std::uniform_real_distribution<float>(0.0f, 1.0f);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();

    std::map<std::string, std::string> metadata;
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
    auto ts = soreor.first;
    auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }
    o2::mch::TrackExtrap::setField();
    const double centerMFT[3] = {0, 0, -61.4};
    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    mBz = field->getBz(centerMFT); // Get field at centre of MFT
    LOGF(info, "Bz at center of MFT = %f kZG", mBz);

    if (cfgApplyZShiftFromCCDB) {
      auto* zShift = ccdb->getForTimeStamp<std::vector<float>>(cfgZShiftPath, bc.timestamp());
      if (zShift != nullptr && !zShift->empty()) {
        LOGF(info, "reading z shift %f from %s", (*zShift)[0], cfgZShiftPath.value);
        mZShift = (*zShift)[0];
      } else {
        LOGF(info, "z shift is not found in ccdb path %s. set to 0 cm", cfgZShiftPath.value);
        mZShift = 0;
      }
    } else {
      LOGF(info, "z shift is manually set to %f cm", cfgManualZShift.value);
      mZShift = cfgManualZShift;
    }

    o2::parameters::GRPLHCIFData* grplhcif = ccdb.service->getSpecificForRun<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", bc.runNumber());
    auto collsys = o2::common::core::CollisionSystemType::getCollisionTypeFromGrp(grplhcif);
    if (collsys == o2::common::core::CollisionSystemType::kCollSyspp || collsys == o2::common::core::CollisionSystemType::kCollSysOO) {
      irSourceForCptFetcher = std::string("T0VTX");
    } else {
      irSourceForCptFetcher = std::string("ZNC hadronic");
    }
    LOGF(info, "irSourceForCptFetcher = %s", irSourceForCptFetcher.data());
  }

  void addHistograms()
  {
    fRegistry.add("Event/hCorrFT0CvsMFT", "mult. corr. between FT0C and MFT;multFT0C;N_{MFTsa}", kTH2F, {{600, 0, 60000}, {100, 0, 10000}}, false);

    auto hMuonType = fRegistry.add<TH1>("hMuonType", "muon type", kTH1F, {{5, -0.5f, 4.5f}}, false);
    hMuonType->GetXaxis()->SetBinLabel(1, "MFT-MCH-MID (global muon)");
    hMuonType->GetXaxis()->SetBinLabel(2, "MFT-MCH-MID (global muon other match)");
    hMuonType->GetXaxis()->SetBinLabel(3, "MFT-MCH");
    hMuonType->GetXaxis()->SetBinLabel(4, "MCH-MID");
    hMuonType->GetXaxis()->SetBinLabel(5, "MCH standalone");

    // information at matching plane
    fRegistry.add("MCHMID/correct/hP", "p;p (GeV/c)", kTH1F, {{1000, 0.0f, 100}}, false);
    fRegistry.add("MCHMID/correct/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("MCHMID/correct/hPz", "pz;p_{z} (GeV/c)", kTH1F, {{1000, -100, 0}}, false);
    fRegistry.add("MCHMID/correct/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {80, -5.f, -1.f}}, false);
    fRegistry.add("MCHMID/correct/hPEta", "p vs. #eta;#eta;p (GeV/c)", kTH2F, {{80, -5.f, -1.f}, {1000, 0, 100}}, false);
    fRegistry.add("MCHMID/correct/hPtEta", "pT vs. #eta;#eta;p_{T} (GeV/c)", kTH2F, {{80, -5.f, -1.f}, {1000, 0, 10}}, false);
    fRegistry.add("MCHMID/correct/hPzEta", "pz vs. #eta;#eta;p_{z} (GeV/c)", kTH2F, {{80, -5.f, -1.f}, {1000, -100, 0}}, false);
    fRegistry.addClone("MCHMID/correct/", "MCHMID/wrong/");
    fRegistry.addClone("MCHMID/", "MFT/");
  }

  template <typename TCollision>
  bool isSelectedCollision(TCollision const& collision)
  {
    if (!(eventCutGroup.cfgMinZvtx < collision.posZ() && collision.posZ() < eventCutGroup.cfgMaxZvtx)) {
      return false;
    }
    if (eventCutGroup.cfgRequireFT0AND && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (eventCutGroup.cfgRequireNoTFB && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (eventCutGroup.cfgRequireNoITSROFB && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (eventCutGroup.cfgRequireNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (eventCutGroup.cfgRequireGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    return true;
  }

  SliceCache cache;
  Preslice<aod::FwdTracks> perCollision = o2::aod::fwdtrack::collisionId;
  Preslice<aod::MFTTracks> perCollision_MFT = o2::aod::fwdtrack::collisionId;
  PresliceUnsorted<aod::FwdTracks> fwdtracksPerMCHTrack = aod::fwdtrack::matchMCHTrackId;
  Partition<MyFwdTracksMC> glMuons = o2::aod::fwdtrack::trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack);
  Partition<MyFwdTracksMC> saMuons = o2::aod::fwdtrack::trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack);

  std::unordered_map<int, int> map_mfttrackcovs;

  void processWithMFTCov(MyCollisionsMC const& collisions, aod::BCsWithTimestamps const&, MyFwdTracksMC const&, MyMFTTracksMC const& mfttracks, aod::MFTTracksCov const& mftCovs, aod::McParticles const&, aod::McCollisions const&)
  {
    for (const auto& mfttrackCov : mftCovs) {
      map_mfttrackcovs[mfttrackCov.matchMFTTrackId()] = mfttrackCov.globalIndex();
    }

    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.has_mcCollision()) {
        continue;
      }

      if (!isSelectedCollision(collision)) {
        continue;
      }

      if (eventCutGroup.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }

      float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSourceForCptFetcher) * 1.e-3; // kHz

      auto mfttracks_per_collision = mfttracks.sliceBy(perCollision_MFT, collision.globalIndex());
      uint16_t multMFT = mfttracks_per_collision.size();
      fRegistry.fill(HIST("Event/hCorrFT0CvsMFT"), collision.multFT0C(), mfttracks_per_collision.size());

      auto saMuons_per_coll = saMuons->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      for (const auto& mchtrack : saMuons_per_coll) {
        if (!mchtrack.has_mcParticle()) {
          continue;
        }
        if (mchtrack.chi2() < 0.f) {
          continue;
        }
        if (mchtrack.chi2MatchMCHMID() < 0.f) {
          continue;
        }

        if (dist01(engine) > cfgDownSampling) { // random sampling, if necessary
          continue;
        }

        auto mcParticle_MCHMID = mchtrack.template mcParticle_as<aod::McParticles>(); // this is identical to mcParticle_MFTMCHMID
        bool isPrimary_MCHMID = mcParticle_MCHMID.isPhysicalPrimary() || mcParticle_MCHMID.producedByGenerator();

        auto muonAtMP = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToMatchingPlane, matchingZ, mBz, mZShift); // propagated to matching plane
        float signed1PtMCHMIDatMP = muonAtMP.getInvQPt();
        float tglMCHMIDatMP = muonAtMP.getTanl();
        float phiMCHMIDatMP = RecoDecay::constrainAngle(muonAtMP.getPhi(), 0, 1U);
        float xMCHMIDatMP = muonAtMP.getX();
        float yMCHMIDatMP = muonAtMP.getY();

        float tglErrMCHMIDatMP = std::sqrt(muonAtMP.getSigma2Tanl());
        float phiErrMCHMIDatMP = std::sqrt(muonAtMP.getSigma2Phi());
        float xErrMCHMIDatMP = std::sqrt(muonAtMP.getSigma2X());
        float yErrMCHMIDatMP = std::sqrt(muonAtMP.getSigma2Y());
        if (muonAtMP.getEta() < MCHMIDCutGroup.minEta || MCHMIDCutGroup.maxEta < muonAtMP.getEta()) {
          continue;
        }

        // auto glMuons_per_MCHMID = fwdtracks.sliceBy(fwdtracksPerMCHTrack, mchtrack.globalIndex());
        auto glMuons_per_MCHMID = glMuons->sliceByCachedUnsorted(o2::aod::fwdtrack::matchMCHTrackId, mchtrack.globalIndex(), cache);
        // LOGF(info, "glMuons_per_MCHMID.size() = %d", glMuons_per_MCHMID.size());

        for (const auto& fwdtrack : glMuons_per_MCHMID) {
          // LOGF(info, "mchtrack.globalIndex() = %d, fwdtrack.globalIndex() = %d, mchtrack.collisionId() = %d, fwdtrack.collisionId() = %d, fwdtrack.trackType() = %d", mchtrack.globalIndex(), fwdtrack.globalIndex(), mchtrack.collisionId(), fwdtrack.collisionId(), fwdtrack.trackType());
          if (!fwdtrack.has_mcParticle()) {
            continue;
          }
          // auto mcParticle_MFTMCHMID = fwdtrack.template mcParticle_as<aod::McParticles>(); // this is identical to mcParticle_MCHMID
          if (fwdtrack.chi2MatchMCHMFT() < 0.f || maxMatchingChi2MCHMFT < fwdtrack.chi2MatchMCHMFT()) {
            continue;
          }
          if (fwdtrack.chi2MatchMCHMID() < 0.f) {
            continue;
          }

          auto mfttrack = fwdtrack.template matchMFTTrack_as<MyMFTTracksMC>(); // MFTsa
          if (mfttrack.eta() < MFTCutGroup.minEta || MFTCutGroup.maxEta < mfttrack.eta() || mfttrack.pt() < MFTCutGroup.minPt || mfttrack.chi2() < 0.f) {
            continue;
          }
          if (!mfttrack.has_mcParticle()) {
            continue;
          }
          auto mcParticle_MFT = mfttrack.template mcParticle_as<aod::McParticles>();
          bool isPrimary_MFT = mcParticle_MFT.isPhysicalPrimary() || mcParticle_MFT.producedByGenerator();

          bool isMatched = (mcParticle_MFT.globalIndex() == mcParticle_MCHMID.globalIndex()) && (mcParticle_MFT.mcCollisionId() == mcParticle_MCHMID.mcCollisionId());

          auto mfttrackcov = mftCovs.rawIteratorAt(map_mfttrackcovs[mfttrack.globalIndex()]);
          o2::track::TrackParCovFwd mftsaAtMP = getTrackParCovFwdShift(mfttrack, mZShift, mfttrackcov); // values at innermost update
          mftsaAtMP.propagateToZhelix(matchingZ, mBz);                                                  // propagated to matching plane

          float signed1PtMFTatMP = mftsaAtMP.getInvQPt();
          float tglMFTatMP = mftsaAtMP.getTanl();
          float phiMFTatMP = RecoDecay::constrainAngle(mftsaAtMP.getPhi(), 0, 1U);
          float xMFTatMP = mftsaAtMP.getX();
          float yMFTatMP = mftsaAtMP.getY();

          float tglErrMFTatMP = std::sqrt(mftsaAtMP.getSigma2Tanl());
          float phiErrMFTatMP = std::sqrt(mftsaAtMP.getSigma2Phi());
          float xErrMFTatMP = std::sqrt(mftsaAtMP.getSigma2X());
          float yErrMFTatMP = std::sqrt(mftsaAtMP.getSigma2Y());

          trackTable(collision.posZ(), collision.multFT0C(), collision.ft0cOccupancyInTimeRange(), hadronicRate, multMFT,
                     signed1PtMFTatMP, tglMFTatMP, phiMFTatMP, xMFTatMP, yMFTatMP,
                     signed1PtMCHMIDatMP, tglMCHMIDatMP, phiMCHMIDatMP, xMCHMIDatMP, yMCHMIDatMP,
                     fwdtrack.chi2MatchMCHMFT(),
                     mcParticle_MFT.pdgCode(), isPrimary_MFT,
                     mcParticle_MCHMID.pdgCode(), isPrimary_MCHMID,
                     isMatched);

          trackErrTable(tglErrMFTatMP, phiErrMFTatMP, xErrMFTatMP, yErrMFTatMP,
                        tglErrMCHMIDatMP, phiErrMCHMIDatMP, xErrMCHMIDatMP, yErrMCHMIDatMP);

          if (isMatched) {
            fRegistry.fill(HIST("MCHMID/correct/hPt"), muonAtMP.getPt());
            fRegistry.fill(HIST("MCHMID/correct/hPz"), muonAtMP.getPz());
            fRegistry.fill(HIST("MCHMID/correct/hP"), muonAtMP.getP());
            fRegistry.fill(HIST("MCHMID/correct/hEtaPhi"), RecoDecay::constrainAngle(muonAtMP.getPhi(), 0, 1U), muonAtMP.getEta());
            fRegistry.fill(HIST("MCHMID/correct/hPEta"), muonAtMP.getEta(), muonAtMP.getP());
            fRegistry.fill(HIST("MCHMID/correct/hPtEta"), muonAtMP.getEta(), muonAtMP.getPt());
            fRegistry.fill(HIST("MCHMID/correct/hPzEta"), muonAtMP.getEta(), muonAtMP.getPz());

            fRegistry.fill(HIST("MFT/correct/hPt"), mftsaAtMP.getPt());
            fRegistry.fill(HIST("MFT/correct/hPz"), mftsaAtMP.getPz());
            fRegistry.fill(HIST("MFT/correct/hP"), mftsaAtMP.getP());
            fRegistry.fill(HIST("MFT/correct/hEtaPhi"), RecoDecay::constrainAngle(mftsaAtMP.getPhi(), 0, 1U), mftsaAtMP.getEta());
            fRegistry.fill(HIST("MFT/correct/hPEta"), mftsaAtMP.getEta(), mftsaAtMP.getP());
            fRegistry.fill(HIST("MFT/correct/hPtEta"), mftsaAtMP.getEta(), mftsaAtMP.getPt());
            fRegistry.fill(HIST("MFT/correct/hPzEta"), mftsaAtMP.getEta(), mftsaAtMP.getPz());
          } else {
            fRegistry.fill(HIST("MCHMID/wrong/hPt"), muonAtMP.getPt());
            fRegistry.fill(HIST("MCHMID/wrong/hPz"), muonAtMP.getPz());
            fRegistry.fill(HIST("MCHMID/wrong/hP"), muonAtMP.getP());
            fRegistry.fill(HIST("MCHMID/wrong/hEtaPhi"), RecoDecay::constrainAngle(muonAtMP.getPhi(), 0, 1U), muonAtMP.getEta());
            fRegistry.fill(HIST("MCHMID/wrong/hPEta"), muonAtMP.getEta(), muonAtMP.getP());
            fRegistry.fill(HIST("MCHMID/wrong/hPtEta"), muonAtMP.getEta(), muonAtMP.getPt());
            fRegistry.fill(HIST("MCHMID/wrong/hPzEta"), muonAtMP.getEta(), muonAtMP.getPz());

            fRegistry.fill(HIST("MFT/wrong/hPt"), mftsaAtMP.getPt());
            fRegistry.fill(HIST("MFT/wrong/hPz"), mftsaAtMP.getPz());
            fRegistry.fill(HIST("MFT/wrong/hP"), mftsaAtMP.getP());
            fRegistry.fill(HIST("MFT/wrong/hEtaPhi"), RecoDecay::constrainAngle(mftsaAtMP.getPhi(), 0, 1U), mftsaAtMP.getEta());
            fRegistry.fill(HIST("MFT/wrong/hPEta"), mftsaAtMP.getEta(), mftsaAtMP.getP());
            fRegistry.fill(HIST("MFT/wrong/hPtEta"), mftsaAtMP.getEta(), mftsaAtMP.getPt());
            fRegistry.fill(HIST("MFT/wrong/hPzEta"), mftsaAtMP.getEta(), mftsaAtMP.getPz());
          }
        }
      } // end of MCHMID track

    } // end of collision loop

    map_mfttrackcovs.clear();
  }
  PROCESS_SWITCH(TreeCreatorMuonML, processWithMFTCov, "produce ML input for single track level", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TreeCreatorMuonML>(cfgc, TaskName{"tree-creator-muon-ml"})};
}
