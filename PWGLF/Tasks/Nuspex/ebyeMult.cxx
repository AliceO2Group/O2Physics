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

#include <vector>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
// #include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "TDatabasePDG.h"
#include "TFormula.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>;
using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;

namespace
{
constexpr float dcaSels[3]{10., 10., 10.};
static const std::vector<std::string> dcaSelsNames{"dcaxy", "dcaz", "dca"};
static const std::vector<std::string> particleName{"tracks"};
} // namespace

struct CandidateTrack {
  float pt = -999.f;
  float eta = -999.f;
  float dcapv = 0;
  float dcaxypv = 0;
  float dcazpv = 0;
  float genpt = -999.f;
  float geneta = -999.f;
  int pdgcode = -999;
  bool isreco = 0;
  int64_t mcIndex = -999;
  int64_t globalIndex = -999;
};

struct CandidateEvent {
  int nTrkRec = -1;
  int nTklRec = -1;
};

struct tagRun2V0MCalibration {
  bool mCalibrationStored = false;
  TH1* mhVtxAmpCorrV0A = nullptr;
  TH1* mhVtxAmpCorrV0C = nullptr;
  TH1* mhMultSelCalib = nullptr;
  float mMCScalePars[6] = {0.0};
  TFormula* mMCScale = nullptr;
} Run2V0MInfo;

enum PartTypes {
  kPi = 0,
  kKa = 1,
  kPr = 2,
  kEl = 3,
  kMu = 4,
  kSig = 5,
  kXi = 6,
  kOm = 7,
  kOther = 8
};

struct ebyeMult {
  std::vector<CandidateTrack> candidateTracks;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  CandidateEvent candidateEvent;

  int mRunNumber;
  float d_bz;
  uint8_t nTrackletsColl;

  ConfigurableAxis centAxis{"centAxis", {106, 0, 106}, "binning for the centrality"};
  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  ConfigurableAxis multAxis{"multAxis", {100, 0, 10000}, "Binning for the multiplicity axis"};
  ConfigurableAxis multFt0Axis{"multFt0Axis", {100, 0, 100000}, "Binning for the ft0 multiplicity axis"};
  Configurable<std::string> genName{"genname", "", "Genearator name: HIJING, PYTHIA8, ... Default: \"\""};

  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.8f, "maximum eta"};

  Configurable<float> ptMin{"ptMin", 0.4f, "minimum pT (GeV/c)"};
  Configurable<float> ptMax{"ptMax", 4.f, "maximum pT (GeV/c)"};

  Configurable<float> trackNcrossedRows{"trackNcrossedRows", 70, "Minimum number of crossed TPC rows"};
  Configurable<float> trackNclusItsCut{"trackNclusITScut", 2, "Minimum number of ITS clusters"};
  Configurable<float> trackNclusTpcCut{"trackNclusTPCcut", 60, "Minimum number of TPC clusters"};
  Configurable<float> trackChi2Cut{"trackChi2Cut", 4.f, "Maximum chi2/ncls in TPC"};
  Configurable<LabeledArray<float>> cfgDcaSels{"cfgDcaSels", {dcaSels, 1, 3, particleName, dcaSelsNames}, "DCA selections"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Preslice<TracksFull> perCollisionTracksFull = o2::aod::track::collisionId;
  Preslice<aod::McParticles> perCollisionMcParts = o2::aod::mcparticle::mcCollisionId;

  // TODO: add function to extract the particle type based on the pdg code
  int getPartType(int const pdgCode)
  {
    switch (std::abs(pdgCode)) {
      case 211:
        return PartTypes::kPi;
      case 321:
        return PartTypes::kKa;
      case 2212:
        return PartTypes::kPr;
      case 11:
        return PartTypes::kEl;
      case 13:
        return PartTypes::kMu;
      case 3222:
        return PartTypes::kSig;
      case 3112:
        return PartTypes::kSig;
      case 3312:
        return PartTypes::kXi;
      case 3334:
        return PartTypes::kOm;
      default:
        return PartTypes::kOther;
    }
  }

  template <class T>
  bool selectTrack(T const& track)
  {
    if (std::abs(track.eta()) > etaMax) {
      return false;
    }
    if (!(track.itsClusterMap() & 0x01) && !(track.itsClusterMap() & 0x02)) {
      return false;
    }
    if (track.itsNCls() < trackNclusItsCut ||
        track.tpcNClsFound() < trackNclusTpcCut ||
        track.tpcNClsCrossedRows() < trackNcrossedRows ||
        track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
        track.tpcChi2NCl() > trackChi2Cut ||
        track.itsChi2NCl() > 36.f) {
      return false;
    }
    if (doprocessRun2 || doprocessMcRun2) {
      if (!(track.trackType() & o2::aod::track::Run2Track) ||
          !(track.flags() & o2::aod::track::TPCrefit) ||
          !(track.flags() & o2::aod::track::ITSrefit)) {
        return false;
      }
    }
    return true;
  }

  template <class Bc>
  void initCCDB(Bc const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    auto timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (doprocessRun2 || doprocessMcRun2) {
      auto grpPath{"GLO/GRP/GRP"};
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (!grpo) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpPath << " of object GRPObject for timestamp " << timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      TList* callst = ccdb->getForTimeStamp<TList>("Centrality/Estimators", bc.timestamp());
      auto getccdb = [callst](const char* ccdbhname) {
        TH1* h = reinterpret_cast<TH1*>(callst->FindObject(ccdbhname));
        return h;
      };
      auto getformulaccdb = [callst](const char* ccdbhname) {
        TFormula* f = reinterpret_cast<TFormula*>(callst->FindObject(ccdbhname));
        return f;
      };
      Run2V0MInfo.mhVtxAmpCorrV0A = getccdb("hVtx_fAmplitude_V0A_Normalized");
      Run2V0MInfo.mhVtxAmpCorrV0C = getccdb("hVtx_fAmplitude_V0C_Normalized");
      Run2V0MInfo.mhMultSelCalib = getccdb("hMultSelCalib_V0M");
      Run2V0MInfo.mMCScale = getformulaccdb(TString::Format("%s-V0M", genName->c_str()).Data());
      if ((Run2V0MInfo.mhVtxAmpCorrV0A != nullptr) && (Run2V0MInfo.mhVtxAmpCorrV0C != nullptr) && (Run2V0MInfo.mhMultSelCalib != nullptr)) {
        if (genName->length() != 0) {
          if (Run2V0MInfo.mMCScale != nullptr) {
            for (int ixpar = 0; ixpar < 6; ++ixpar) {
              Run2V0MInfo.mMCScalePars[ixpar] = Run2V0MInfo.mMCScale->GetParameter(ixpar);
            }
          } else {
            LOGF(fatal, "MC Scale information from V0M for run %d not available", bc.runNumber());
          }
        }
        Run2V0MInfo.mCalibrationStored = true;
      } else {
        LOGF(fatal, "Calibration information from V0M for run %d corrupted", bc.runNumber());
      }
    } else {
      auto grpmagPath{"GLO/Config/GRPMagField"};
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField for timestamp " << timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
    }
    // Fetch magnetic field from ccdb for current collision
    d_bz = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp << " with magnetic field of " << d_bz << " kG";
    mRunNumber = bc.runNumber();
  }

  // float getV0M(int64_t const id, float const zvtx, aod::FV0As const& fv0as, aod::FV0Cs const& fv0cs)
  // {
  //   auto fv0a = fv0as.rawIteratorAt(id);
  //   auto fv0c = fv0cs.rawIteratorAt(id);
  //   float multFV0A = 0;
  //   float multFV0C = 0;
  //   for (float amplitude : fv0a.amplitude()) {
  //     multFV0A += amplitude;
  //   }

  //   for (float amplitude : fv0c.amplitude()) {
  //     multFV0C += amplitude;
  //   }

  //   float v0m = -1;
  //   auto scaleMC = [](float x, float pars[6]) {
  //     return pow(((pars[0] + pars[1] * pow(x, pars[2])) - pars[3]) / pars[4], 1.0f / pars[5]);
  //   };

  //   if (Run2V0MInfo.mMCScale != nullptr) {
  //     float multFV0M = multFV0A + multFV0C;
  //     v0m = scaleMC(multFV0M, Run2V0MInfo.mMCScalePars);
  //     LOGF(debug, "Unscaled v0m: %f, scaled v0m: %f", multFV0M, v0m);
  //   } else {
  //     v0m = multFV0A * Run2V0MInfo.mhVtxAmpCorrV0A->GetBinContent(Run2V0MInfo.mhVtxAmpCorrV0A->FindFixBin(zvtx)) +
  //           multFV0C * Run2V0MInfo.mhVtxAmpCorrV0C->GetBinContent(Run2V0MInfo.mhVtxAmpCorrV0C->FindFixBin(zvtx));
  //   }
  //   return v0m;
  // }

  void init(o2::framework::InitContext&)
  {

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // event QA
    histos.add<TH1>("QA/zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});
    histos.add<TH2>("QA/V0MvsCL0", ";Centrality CL0 (%);Centrality V0M (%)", HistType::kTH2F, {centAxis, centAxis});
    histos.add<TH2>("QA/trackletsVsV0M", ";Centrality CL0 (%);Centrality V0M (%)", HistType::kTH2F, {centAxis, multAxis});
    histos.add<TH2>("QA/nTrklCorrelation", ";Tracklets |#eta| > 0.7; Tracklets |#eta| < 0.6", HistType::kTH2D, {{201, -0.5, 200.5}, {201, -0.5, 200.5}});
    histos.add<TH1>("QA/TrklEta", ";Tracklets #eta; Entries", HistType::kTH1D, {{100, -3., 3.}});

    // rec tracks
    histos.add<TH3>("RecTracks", ";Tracklets |#eta| > 0.7;#it{p}_{T} (GeV/#it{c});DCA_{#it{xy}} (cm)", HistType::kTH3D, {{201, -0.5, 200.5}, {100, -5., 5.}, {200, -1., 1.}});

    // rec  tracks and tracklets distribution
    histos.add<TH2>("TracksDistr", ";Tracklets |#eta| > 0.7;#it{N}_{trk}", HistType::kTH2D, {{201, -0.5, 200.5}, {201, -0.5, 200.5}});
    histos.add<TH2>("TrackletsDistr", ";Tracklets |#eta| > 0.7;#it{N}_{tkl}", HistType::kTH2D, {{201, -0.5, 200.5}, {201, -0.5, 200.5}});

    if (doprocessMcRun2) {
      // rec & gen particles (per species)
      histos.add<TH3>("RecPart", ";Tracklets |#eta| > 0.7;#it{p}_{T} (GeV/#it{c});Species", HistType::kTH3D, {{201, -0.5, 200.5}, {100, -5., 5.}, {10, 0, 10}});
      histos.add<TH3>("GenPart", ";Tracklets |#eta| > 0.7;#it{p}_{T} (GeV/#it{c});Species", HistType::kTH3D, {{201, -0.5, 200.5}, {100, -5., 5.}, {10, 0, 10}});

      // dca_xy templates
      histos.add<TH3>("PrimTracks", ";Tracklets |#eta| > 0.7;#it{p}_{T} (GeV/#it{c});DCA_{#it{xy}} (cm)", HistType::kTH3D, {{201, -0.5, 200.5}, {100, -5., 5.}, {200, -1., 1.}});
      histos.add<TH3>("SecWDTracks", ";Tracklets |#eta| > 0.7;#it{p}_{T} (GeV/#it{c});DCA_{#it{xy}} (cm)", HistType::kTH3D, {{201, -0.5, 200.5}, {100, -5., 5.}, {200, -1., 1.}});
      histos.add<TH3>("SecTracks", ";Tracklets |#eta| > 0.7;#it{p}_{T} (GeV/#it{c});DCA_{#it{xy}} (cm)", HistType::kTH3D, {{201, -0.5, 200.5}, {100, -5., 5.}, {200, -1., 1.}});

      // response
      histos.add<TH3>("GenRecTracks", ";Tracklets |#eta| > 0.7#it;#it{N}_{trk};#it{N}_{gen}", HistType::kTH3D, {{201, -0.5, 200.5}, {201, -0.5, 200.5}, {201, -0.5, 200.5}});
      histos.add<TH3>("GenRecTracklets", ";Tracklets |#eta| > 0.7;#it{N}_{tkl};#it{N}_{gen}", HistType::kTH3D, {{201, -0.5, 200.5}, {201, -0.5, 200.5}, {201, -0.5, 200.5}});
    }
  }

  template <class C, class T>
  void fillRecoEvent(C const& collision, T const& tracksAll /* , float const& centrality */)
  {
    auto tracks = tracksAll.sliceBy(perCollisionTracksFull, collision.globalIndex());
    candidateTracks.clear();

    gpu::gpustd::array<float, 2> dcaInfo;
    int nTracklets[2]{0, 0};
    int nTracks{0};
    for (const auto& track : tracks) {

      if (track.trackType() == 255 && std::abs(track.eta()) < 1.2) { // tracklet
        if (std::abs(track.eta()) < 0.6)
          nTracklets[0]++;
        else if (std::abs(track.eta()) > 0.7)
          nTracklets[1]++;
      }

      if (!selectTrack(track)) {
        continue;
      }

      auto trackParCov = getTrackParCov(track);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCov, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrNONE, &dcaInfo);
      auto dca = std::hypot(dcaInfo[0], dcaInfo[1]);
      auto trackPt = trackParCov.getPt();
      auto trackEta = trackParCov.getEta();
      if (dca > cfgDcaSels->get("dca")) { // dca
        continue;
      }
      if (std::abs(dcaInfo[1]) > cfgDcaSels->get("dcaz")) { // dcaz
        continue;
      }

      CandidateTrack candTrack;
      candTrack.pt = track.sign() > 0. ? trackPt : -trackPt;
      candTrack.eta = trackEta;
      candTrack.dcapv = dca;
      candTrack.dcaxypv = dcaInfo[0];
      candTrack.dcazpv = dcaInfo[1];
      candTrack.globalIndex = track.globalIndex();
      candidateTracks.push_back(candTrack);

      if (std::abs(dcaInfo[0]) < cfgDcaSels->get("dcaxy")) { // dcaxy
        ++nTracks;
      }
    }

    histos.fill(HIST("QA/nTrklCorrelation"), nTracklets[0], nTracklets[1]);
    nTrackletsColl = nTracklets[1];

    candidateEvent.nTklRec = nTracklets[0];
    histos.fill(HIST("TracksDistr"), nTracklets[1], nTracks);
    histos.fill(HIST("TrackletsDistr"), nTracklets[1], nTracklets[0]);
  }

  template <class C, class T>
  void fillMcEvent(C const& collision, T const& tracks /* , float const& centrality */, aod::McParticles const&, aod::McTrackLabels const& mcLabels)
  {
    fillRecoEvent<C, T>(collision, tracks /* , centrality */);

    int nTracks{0};
    for (auto& candidateTrack : candidateTracks) {
      candidateTrack.isreco = true;

      auto mcLab = mcLabels.rawIteratorAt(candidateTrack.globalIndex);
      if (mcLab.has_mcParticle()) {
        auto mcTrack = mcLab.template mcParticle_as<aod::McParticles>();
        if (((mcTrack.flags() & 0x8) && (doprocessMcRun2)) || (mcTrack.flags() & 0x2))
          continue;
        if (!mcTrack.isPhysicalPrimary()) {
          if (mcTrack.has_mothers()) { // sec WD
            histos.fill(HIST("SecWDTracks"), nTrackletsColl, candidateTrack.pt, candidateTrack.dcaxypv);
          } else { // from material
            histos.fill(HIST("SecTracks"), nTrackletsColl, candidateTrack.pt, candidateTrack.dcaxypv);
          }
        }
        if ((mcTrack.flags() & 0x1))
          continue;
        if (mcTrack.isPhysicalPrimary()) { // primary
          histos.fill(HIST("PrimTracks"), nTrackletsColl, candidateTrack.pt, candidateTrack.dcaxypv);
        }

        if (std::abs(candidateTrack.dcaxypv) > cfgDcaSels->get("dcaxy"))
          continue;
        int partType = getPartType(mcTrack.pdgCode());
        if (mcTrack.isPhysicalPrimary()) { // primary
          histos.fill(HIST("RecPart"), nTrackletsColl, candidateTrack.pt, partType);
          if (partType < PartTypes::kOther) {
            ++nTracks;
          }
        }
        auto genPt = std::hypot(mcTrack.px(), mcTrack.py());
        candidateTrack.pdgcode = mcTrack.pdgCode();
        candidateTrack.genpt = genPt;
        candidateTrack.geneta = mcTrack.eta();
        candidateTrack.mcIndex = mcTrack.globalIndex();
      }
    }
    candidateEvent.nTrkRec = nTracks;
  }

  void fillMcGen(aod::McParticles const& mcParticles, aod::McTrackLabels const& /*mcLab*/, uint64_t const& collisionId)
  {
    int nParticles = 0;
    auto mcParticles_thisCollision = mcParticles.sliceBy(perCollisionMcParts, collisionId);
    for (auto& mcPart : mcParticles_thisCollision) {
      auto genEta = mcPart.eta();
      if (std::abs(genEta) > etaMax) {
        continue;
      }
      if (((mcPart.flags() & 0x8) && (doprocessMcRun2)) || (mcPart.flags() & 0x2) || (mcPart.flags() & 0x1))
        continue;
      if (!mcPart.isPhysicalPrimary() /* && !mcPart.has_mothers() */)
        continue;
      auto genPt = std::hypot(mcPart.px(), mcPart.py());
      CandidateTrack candTrack;
      candTrack.genpt = genPt;
      candTrack.geneta = mcPart.eta();
      candTrack.pdgcode = mcPart.pdgCode();

      int partType = getPartType(mcPart.pdgCode());
      if (partType < PartTypes::kOther) {
        ++nParticles;
      }
      histos.fill(HIST("GenPart"), nTrackletsColl, mcPart.pdgCode() > 0 ? genPt : -genPt, partType);

      auto it = find_if(candidateTracks.begin(), candidateTracks.end(), [&](CandidateTrack trk) { return trk.mcIndex == mcPart.globalIndex(); });
      if (it != candidateTracks.end()) {
        continue;
      } else {
        candidateTracks.emplace_back(candTrack);
      }
    }
    histos.fill(HIST("GenRecTracks"), nTrackletsColl, candidateEvent.nTrkRec, nParticles);
    histos.fill(HIST("GenRecTracklets"), nTrackletsColl, candidateEvent.nTklRec, nParticles);
  }

  void processRun2(soa::Join<aod::Collisions, aod::EvSels> const& collisions, TracksFull const& tracks /* , aod::FV0As const& fv0as, aod::FV0Cs const& fv0cs */, BCsWithRun2Info const&)
  {

    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<BCsWithRun2Info>();
      initCCDB(bc);

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      if (!collision.alias_bit(kINT7))
        continue;

      if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted)))
        continue;

      // float v0m = getV0M(bc.globalIndex(), collision.posZ(), fv0as, fv0cs);
      // float cV0M = Run2V0MInfo.mhMultSelCalib->GetBinContent(Run2V0MInfo.mhMultSelCalib->FindFixBin(v0m));

      histos.fill(HIST("QA/zVtx"), collision.posZ());

      fillRecoEvent(collision, tracks /* , cV0M */);

      for (auto t : candidateTracks) {
        histos.fill(HIST("RecTracks"), nTrackletsColl, t.pt, t.dcaxypv);
      }
    }
  }
  PROCESS_SWITCH(ebyeMult, processRun2, "process (Run 2)", false);

  void processMcRun2(soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions, aod::McCollisions const& /*mcCollisions*/, TracksFull const& tracks /* , aod::FV0As const& fv0as, aod::FV0Cs const& fv0cs */, aod::McParticles const& mcParticles, aod::McTrackLabels const& mcLab, BCsWithRun2Info const&)
  {

    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<BCsWithRun2Info>();
      initCCDB(bc);

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted)))
        continue;

      // float v0m = getV0M(bc.globalIndex(), collision.posZ(), fv0as, fv0cs);
      // float cV0M = Run2V0MInfo.mhMultSelCalib->GetBinContent(Run2V0MInfo.mhMultSelCalib->FindFixBin(v0m));

      histos.fill(HIST("QA/zVtx"), collision.posZ());

      fillMcEvent(collision, tracks /* , cV0M */, mcParticles, mcLab);
      fillMcGen(mcParticles, mcLab, collision.mcCollisionId());
    }
  }
  PROCESS_SWITCH(ebyeMult, processMcRun2, "process mc (Run 2)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ebyeMult>(cfgc)};
}
