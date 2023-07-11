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

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "CCDB/BasicCCDBManager.h"
#include "DetectorsBase/Propagator.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
  static constexpr int nParticles{3};
  static constexpr int nCutsPID{2};
  static constexpr std::array<float, nParticles> masses{
    constants::physics::MassKaonCharged,
    constants::physics::MassProton,
    constants::physics::MassPionCharged};
  static constexpr std::array<int, nParticles> charges{-1, 1, -1};
  static const std::vector<std::string> matterOrNot{"Matter", "Antimatter"};
  static const std::vector<std::string> particlesNames{"K", "Pr", "Pi"};
  static const std::vector<std::string> cutsNames{"TPCnSigmaMin", "TPCnSigmaMax"};
  static constexpr float cutsPID[nParticles][nCutsPID]{
    {-4.f, +4.f}, /*K*/
    {-4.f, +4.f}, /*Pr*/
    {-4.f, +4.f}, /*Pi*/
  };
  std::shared_ptr<TH2> h2TPCsignal[nParticles];
  std::shared_ptr<TH2> h2TPCnSigma[nParticles];
} // namespace

struct NonPromptOmegaTask {
  using TracksExt = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::McTrackLabels, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTPCFullPr>;

  Configurable<double> bz{"bz", -50., "magnetic field"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> cfgGRPmagPath{"cfgGRPmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> cfgGRPpath{"cfgGRPpath", "GLO/GRP/GRP", "Path of the grp file"};

  Configurable<float> cfgCutNclusTPC{"cfgCutNclusTPC", 70, "Minimum number of TPC clusters"};
  Configurable<LabeledArray<float>> cfgCutsPID{"particlesCutsPID", {cutsPID[0], nParticles, nCutsPID, particlesNames, cutsNames}, "Nuclei PID selections"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  float mBz = 0.f;

  HistogramRegistry registry{
    "registry",
    {
      {"h_dca", "DCA;DCA (cm)", {HistType::kTH1D, {{200, 0., .5}}}},
      {"h_dcaxy", "DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"h_dcaz", "DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"h_bachdcaxyM", "Bachelor DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcaxyAM", "Bachelor DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcazM", "Bachelor DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcazAM", "Bachelor DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_dcavspt", "DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_bachdcavspt", "Bachelor DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_bachdcavsr", "Bachelor DCA vs R (cm);DCA (cm);R (cm)", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 30.}}}},
      {"h_ntrackdcavspt", "N track DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_ptrackdcavspt", "P track DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_dcavsr", "DCA vs R;DCA (cm);R (cm)", {HistType::kTH2D, {{200, -.5, .5}, {200, 0., 10.}}}},
      {"h_massvspt", "Mass vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, 0., 10.}, {200, 0., 10.}}}},
      {"h_buildermassvspt", "Mass (from builder) vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, 0., 10.}, {200, 0., 10.}}}},
      {"h_massvsmass", "Mass vs mass;Mass (GeV/#it{c}^{2});Mass (GeV/#it{c}^{2})", {HistType::kTH2D, {{200, 0., 10.}, {200, 0., 10.}}}},
      {"h_bachelorsign", "Bachelor sign;Sign;Counts", {HistType::kTH1D, {{10, -5., 5.}}}},
      {"h_ptmassdcaxyM", ";p_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});DCA_{xy} (cm)", {HistType::kTH3D, {{200, 0., 10.},{200, 0., 10.},{200, -.5, .5}}}},
      {"h_ptmassdcaxyAM", ";p_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});DCA_{xy} (cm)", {HistType::kTH3D, {{200, 0., 10.},{200, 0., 10.},{200, -.5, .5}}}},
      {"h_ptmassdcazM", ";p_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});DCA_{x} (cm)", {HistType::kTH3D, {{200, 0., 10.},{200, 0., 10.},{200, -.5, .5}}}},
      {"h_ptmassdcazAM", ";p_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});DCA_{x} (cm)", {HistType::kTH3D, {{200, 0., 10.},{200, 0., 10.},{200, -.5, .5}}}},
    }};

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto run3grp_timestamp = bc.timestamp();
    mRunNumber = bc.runNumber();

    if (o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(cfgGRPpath, run3grp_timestamp)) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      LOG(debug) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << mBz << " kZG";
    } else if (o2::parameters::GRPMagField* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(cfgGRPmagPath, run3grp_timestamp)) {
      o2::base::Propagator::initFieldFromGRP(grpmag);
      LOG(debug) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << mBz << " kZG";
    } else {
      LOG(fatal) << "Got nullptr from CCDB for path " << cfgGRPpath << " of object GRPMagField and " << cfgGRPmagPath << " of object GRPObject for timestamp " << run3grp_timestamp;
    }
  }

  void init(InitContext const&)
  {
    if (static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value) == o2::base::Propagator::MatCorrType::USEMatCorrLUT) {
      auto* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
      o2::base::Propagator::Instance(true)->setMatLUT(lut);
    }

    std::vector<double> ptBinning = {0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    registry.add("fTPCsignal", "Specific energy loss", HistType::kTH2F, {{1200, -6, 6, "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    registry.add("fCollZpos", "collision z position", HistType::kTH1F, {{600, -20., +20., "z position (cm)"}});

    for (int iN{0}; iN < nParticles; ++iN) {
      h2TPCsignal[iN] = registry.add<TH2>(Form("fTPCsignal_%s", particlesNames[iN].data()), "Specific energy loss", HistType::kTH2F, {{1200, -6, 6., "#it{p}/Z (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
      h2TPCnSigma[iN] = registry.add<TH2>(Form("fTPCcounts_%s", particlesNames[iN].data()), "n-sigma TPC", HistType::kTH2F, {{100, -5, 5, "#it{p} /Z (GeV/#it{c})"}, {200, -10., +10., Form("n#sigma_%s (a. u.)", particlesNames[iN].data())}});
    }
    auto scalers{std::get<std::shared_ptr<TH1>>(registry.add("fProcessedEvents", ";;Number of filtered events", HistType::kTH1F, {{nParticles + 1, -0.5, nParticles + 0.5}}))};
    scalers->GetXaxis()->SetBinLabel(1, "Processed events");
    for (uint32_t iS{0}; iS < particlesNames.size(); ++iS) {
      scalers->GetXaxis()->SetBinLabel(iS + 2, particlesNames[iS].data());
    }
  }

  void processTrackedCascades(aod::Collision const& collision,
                              aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& cascades,
                              aod::V0s const& v0s, TracksExt const& tracks, aod::McParticles const& mcParticles,
                              aod::BCsWithTimestamps const&)
  {
    bool keepEvent[nParticles]{false};

    registry.fill(HIST("fProcessedEvents"), 0);

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    const auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value);

    const auto primaryVertex = getPrimaryVertex(collision);

    for (const auto& trackedCascade : trackedCascades) {
      const auto track = trackedCascade.track_as<TracksExt>();

      const auto& casc = trackedCascade.cascade();
      const auto& bachelor = casc.bachelor_as<TracksExt>();
      const auto& v0 = casc.v0();
      const auto& ptrack = v0.posTrack_as<TracksExt>();
      const auto& ntrack = v0.negTrack_as<TracksExt>();

      float nSigmaTPC[nParticles]{
        bachelor.tpcNSigmaKa(), ptrack.tpcNSigmaPr(), ntrack.tpcNSigmaPi()};

      if (!bachelor.hasTOF() && !ptrack.hasTOF() && !ntrack.hasTOF() ) {
        LOG(debug)<< "no TOF: "<<bachelor.hasTOF()<<"/"<<ptrack.hasTOF()<<"/"<<ntrack.hasTOF();
        continue;
      }

      if (bachelor.tpcNClsFound() < cfgCutNclusTPC || ptrack.tpcNClsFound() < cfgCutNclusTPC || ntrack.tpcNClsFound() < cfgCutNclusTPC) {
        LOG(debug)<< "no tpcNClsFound: "<<bachelor.tpcNClsFound()<<"/"<<ptrack.tpcNClsFound()<<"/"<<ntrack.tpcNClsFound();
        continue;
      }

      h2TPCnSigma[0]->Fill(bachelor.sign() * bachelor.tpcInnerParam(), nSigmaTPC[0]);
      registry.fill(HIST("fTPCsignal"), bachelor.sign() * bachelor.tpcInnerParam(), bachelor.tpcSignal());
      LOG(debug) << "TPCSignal bachelor " << bachelor.sign() << "/"<<bachelor.tpcInnerParam()<<"/"<<bachelor.tpcSignal();
      if (nSigmaTPC[0] < cfgCutsPID->get(0u, 0u) || nSigmaTPC[0] > cfgCutsPID->get(0u, 1u)) {
        continue;
      }
      keepEvent[0] = true;
      if (keepEvent[0]) {
        h2TPCsignal[0]->Fill(bachelor.sign() * bachelor.tpcInnerParam(), bachelor.tpcSignal());
      }

      h2TPCnSigma[1]->Fill(ptrack.sign() * ptrack.tpcInnerParam(), nSigmaTPC[1]);
      registry.fill(HIST("fTPCsignal"), ptrack.sign() * ptrack.tpcInnerParam(), ptrack.tpcSignal());
      LOG(debug) << "TPCSignal ptrack " << ptrack.sign() << "/"<<ptrack.tpcInnerParam()<<"/"<<ptrack.tpcSignal();
      if (nSigmaTPC[1] < cfgCutsPID->get(1u, 0u) || nSigmaTPC[1] > cfgCutsPID->get(1u, 1u)) {
        continue;
      }
      keepEvent[1] = true;
      if (keepEvent[1]) {
        h2TPCsignal[1]->Fill(ptrack.sign() * ptrack.tpcInnerParam(), ptrack.tpcSignal());
      }

      h2TPCnSigma[2]->Fill(ntrack.sign() * ntrack.tpcInnerParam(), nSigmaTPC[2]);
      registry.fill(HIST("fTPCsignal"), ntrack.sign() * ntrack.tpcInnerParam(), ntrack.tpcSignal());
      LOG(debug) << "TPCSignal ntrack " << ntrack.sign() << "/"<<ntrack.tpcInnerParam()<<"/"<<ntrack.tpcSignal();
      if (nSigmaTPC[2] < cfgCutsPID->get(2u, 0u) || nSigmaTPC[2] > cfgCutsPID->get(2u, 1u)) {
        continue;
      }
      keepEvent[2] = true;
      if (keepEvent[2]) {
        h2TPCsignal[2]->Fill(ntrack.sign() * ntrack.tpcInnerParam(), ntrack.tpcSignal());
      }

      registry.fill(HIST("fCollZpos"), collision.posZ());
  
      for (int iDecision{0}; iDecision < 3; ++iDecision) {
        if (keepEvent[iDecision]) {
          registry.fill(HIST("fProcessedEvents"), iDecision + 1);
        }
      }

      auto trackCovTrk = getTrackParCov(track);
      o2::dataformats::DCA impactParameterTrk;
      if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovTrk, mBz, 2.f, matCorr, &impactParameterTrk)) {
        registry.fill(HIST("h_dca"), TMath::Sqrt(impactParameterTrk.getR2()));
        registry.fill(HIST("h_dcaxy"), impactParameterTrk.getY());
        registry.fill(HIST("h_dcaz"), impactParameterTrk.getZ());
        registry.fill(HIST("h_dcavspt"), impactParameterTrk.getY(), track.pt());
        registry.fill(HIST("h_dcavsr"), impactParameterTrk.getY(), std::hypot(track.x(), track.y()));
        registry.fill(HIST("h_massvspt"), trackedCascade.omegaMass(), track.pt());
        if (bachelor.sign()<0){
          registry.fill(HIST("h_ptmassdcazM"),track.pt(), trackedCascade.omegaMass(),  impactParameterTrk.getZ());
          registry.fill(HIST("h_ptmassdcaxyM"),track.pt(), trackedCascade.omegaMass(),  impactParameterTrk.getY());
        }
        else if (bachelor.sign()>0){
          registry.fill(HIST("h_ptmassdcazAM"),track.pt(), trackedCascade.omegaMass(),  impactParameterTrk.getZ());
          registry.fill(HIST("h_ptmassdcaxyAM"),track.pt(), trackedCascade.omegaMass(),  impactParameterTrk.getY());
        }
      }

      LOGF(debug, "ptrack (id: %d, pdg: %d) has mother %d", ptrack.mcParticleId(),
           ptrack.mcParticle().pdgCode(), ptrack.mcParticle().has_mothers() ? ptrack.mcParticle().mothersIds()[0] : -1);
      LOGF(debug, "ntrack (id: %d, pdg: %d) has mother %d", ntrack.mcParticleId(),
           ntrack.mcParticle().pdgCode(), ntrack.mcParticle().has_mothers() ? ntrack.mcParticle().mothersIds()[0] : -1);

      LOG(debug) << "bachelor with PDG code: " << bachelor.mcParticle().pdgCode() << ". Charge: " << bachelor.sign();
      if (ptrack.mcParticle().has_mothers() && ntrack.mcParticle().has_mothers() &&
          ptrack.mcParticle().mothersIds()[0] == ntrack.mcParticle().mothersIds()[0]) {
        const auto v0part = ptrack.mcParticle().mothers_as<aod::McParticles>()[0];
        LOG(debug) << "v0 with PDG code: " << v0part.pdgCode();
        if (v0part.has_mothers() && bachelor.mcParticle().has_mothers() &&
            v0part.mothersIds()[0] == bachelor.mcParticle().mothersIds()[0]) {
          LOG(debug) << "cascade with PDG code: " << v0part.mothers_as<aod::McParticles>()[0].pdgCode();
        } else {
          LOG(debug) << "rejecting particle.";
          continue;
        }
      }
      auto trackCovBach = getTrackParCov(bachelor);
      o2::dataformats::DCA impactParameterBach;
      if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovBach, mBz, 2.f, matCorr, &impactParameterBach)) {
        if (bachelor.sign() < 0) {
          registry.fill(HIST("h_bachdcaxyM"), impactParameterBach.getY());
          registry.fill(HIST("h_bachdcazM"), impactParameterBach.getZ());
        } else if (bachelor.sign() > 0) {
          registry.fill(HIST("h_bachdcaxyAM"), impactParameterBach.getY());
          registry.fill(HIST("h_bachdcazAM"), impactParameterBach.getZ());
        }
        registry.fill(HIST("h_bachdcavspt"), impactParameterBach.getY(), bachelor.pt());
        registry.fill(HIST("h_bachdcavsr"), impactParameterBach.getY(), std::hypot(trackedCascade.decayX(), trackedCascade.decayY()));
        registry.fill(HIST("h_bachelorsign"), bachelor.sign());
      }

      auto trackCovNtrack = getTrackParCov(ntrack);
      o2::dataformats::DCA impactParameterNtrack;
      if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovNtrack, mBz, 2.f, matCorr, &impactParameterNtrack)) {
        registry.fill(HIST("h_ntrackdcavspt"), impactParameterNtrack.getY(), ntrack.pt());
      }

      auto trackCovPtrack = getTrackParCov(ptrack);
      o2::dataformats::DCA impactParameterPtrack;
      if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovPtrack, mBz, 2.f, matCorr, &impactParameterPtrack)) {
        registry.fill(HIST("h_ptrackdcavspt"), impactParameterPtrack.getY(), ptrack.pt());
      }
    }
  }

  PROCESS_SWITCH(NonPromptOmegaTask, processTrackedCascades, "process cascades from strangeness tracking", true);

  void processCascades(aod::Collision const& collision, aod::TrackedCascades const& trackedCascades, aod::Cascades const& cascades, aod::V0s const& v0s,
                       soa::Join<aod::TraCascDatas, aod::McTraCascLabels> const& trackedcascdata, TracksExt const& tracks, aod::McParticles const& mcParticles)
  {
    for (const auto& trackedCascadeData : trackedcascdata) {
      registry.fill(HIST("h_buildermassvspt"), trackedCascadeData.mOmega(), trackedCascadeData.pt());
    }

    for (const auto& trackedCascade : trackedCascades) {
      const auto& casc = trackedCascade.cascade();
      const auto& bachelor = casc.bachelor_as<TracksExt>();
      const auto& v0 = casc.v0();
      const auto& ptrack = v0.posTrack_as<TracksExt>();
      const auto& ntrack = v0.negTrack_as<TracksExt>();
      if (ptrack.mcParticle().has_mothers() && ntrack.mcParticle().has_mothers() &&
          ptrack.mcParticle().mothersIds()[0] == ntrack.mcParticle().mothersIds()[0]) {
        const auto v0part = ptrack.mcParticle().mothers_as<aod::McParticles>()[0];
        if (v0part.has_mothers() && bachelor.mcParticle().has_mothers() &&
            v0part.mothersIds()[0] == bachelor.mcParticle().mothersIds()[0]) {
          int mcid = v0part.mothersIds()[0];
          for (const auto& trackedCascadeData : trackedcascdata) {
            if (trackedCascadeData.mcParticleId() == mcid) {
              registry.fill(HIST("h_massvsmass"), trackedCascade.omegaMass(), trackedCascadeData.mOmega());

              break;
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(NonPromptOmegaTask, processCascades, "process cascades from builder", true);



  
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NonPromptOmegaTask>(cfgc)};
}
