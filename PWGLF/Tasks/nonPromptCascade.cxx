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
#include "Common/DataModel/Centrality.h"
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
static constexpr int nParticles{4};
static constexpr int nCutsPID{2};
static const std::vector<std::string> matterOrNot{"Matter", "Antimatter"};
static const std::vector<std::string> particlesNames{"K-bachelor", "Pi-bachelor", "Pr", "Pi"};
static const std::vector<std::string> cutsNames{"TPCnSigmaMin", "TPCnSigmaMax"};
static constexpr float cutsPID[nParticles][nCutsPID]{
  {-4.f, +4.f}, /*K bachelor*/
  {-4.f, +4.f}, /*Pi bachelor*/
  {-4.f, +4.f}, /*Pr*/
  {-4.f, +4.f}, /*Pi*/
};
std::shared_ptr<TH2> h2TPCsignal[nParticles];
std::shared_ptr<TH2> h2TPCnSigma[nParticles];

std::shared_ptr<TH1> invMassBCOmega;
std::shared_ptr<TH1> invMassACOmega;
std::shared_ptr<TH1> invMassBCXi;
std::shared_ptr<TH1> invMassACXi;

} // namespace

struct NonPromptCascadeTask {

  using TracksExtData = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  using TracksExtMC = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::McTrackLabels, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  using CollisionCandidatesRun3 = soa::Join<aod::Collisions, aod::EvSels>::iterator;

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
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
      {"h_dca_Omega", "DCA;DCA (cm)", {HistType::kTH1D, {{200, 0., .5}}}},
      {"h_dcaxy_Omega", "DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"h_dcaz_Omega", "DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"h_bachdcaxyM_Omega", "Bachelor DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcaxyAM_Omega", "Bachelor DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcazM_Omega", "Bachelor DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcazAM_Omega", "Bachelor DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_dcavspt_Omega", "DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_bachdcavspt_Omega", "Bachelor DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_bachdcavsr_Omega", "Bachelor DCA vs R (cm);DCA (cm);R (cm)", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 30.}}}},
      {"h_ntrackdcavspt_Omega", "N track DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_ptrackdcavspt_Omega", "P track DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_dcavsr_Omega", "DCA vs R;DCA (cm);R (cm)", {HistType::kTH2D, {{200, -.5, .5}, {200, 0., 10.}}}},
      {"h_massvspt_Omega", "Mass vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, 0., 10.}, {200, 0., 10.}}}},
      {"h_buildermassvspt_Omega", "Mass (from builder) vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, 0., 10.}, {200, 0., 10.}}}},
      {"h_massvsmass_Omega", "Mass vs mass;Mass (GeV/#it{c}^{2});Mass (GeV/#it{c}^{2})", {HistType::kTH2D, {{200, 0., 10.}, {200, 0., 10.}}}},
      {"h_bachelorsign_Omega", "Bachelor sign;Sign;Counts", {HistType::kTH1D, {{10, -5., 5.}}}},
      {"h_ptmassdcaxyM_Omega", ";p_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});DCA_{xy} (cm)", {HistType::kTH3D, {{200, 0., 10.}, {200, 0., 10.}, {200, -.5, .5}}}},
      {"h_ptmassdcaxyAM_Omega", ";p_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});DCA_{xy} (cm)", {HistType::kTH3D, {{200, 0., 10.}, {200, 0., 10.}, {200, -.5, .5}}}},
      {"h_ptmassdcazM_Omega", ";p_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});DCA_{x} (cm)", {HistType::kTH3D, {{200, 0., 10.}, {200, 0., 10.}, {200, -.5, .5}}}},
      {"h_ptmassdcazAM_Omega", ";p_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});DCA_{x} (cm)", {HistType::kTH3D, {{200, 0., 10.}, {200, 0., 10.}, {200, -.5, .5}}}},

      {"h_dca_Xi", "DCA;DCA (cm)", {HistType::kTH1D, {{200, 0., .5}}}},
      {"h_dcaxy_Xi", "DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"h_dcaz_Xi", "DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"h_bachdcaxyM_Xi", "Bachelor DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcaxyAM_Xi", "Bachelor DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcazM_Xi", "Bachelor DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcazAM_Xi", "Bachelor DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_dcavspt_Xi", "DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_bachdcavspt_Xi", "Bachelor DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_bachdcavsr_Xi", "Bachelor DCA vs R (cm);DCA (cm);R (cm)", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 30.}}}},
      {"h_ntrackdcavspt_Xi", "N track DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_ptrackdcavspt_Xi", "P track DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_dcavsr_Xi", "DCA vs R;DCA (cm);R (cm)", {HistType::kTH2D, {{200, -.5, .5}, {200, 0., 10.}}}},
      {"h_massvspt_Xi", "Mass vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, 0., 10.}, {200, 0., 10.}}}},
      {"h_buildermassvspt_Xi", "Mass (from builder) vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, 0., 10.}, {200, 0., 10.}}}},
      {"h_massvsmass_Xi", "Mass vs mass;Mass (GeV/#it{c}^{2});Mass (GeV/#it{c}^{2})", {HistType::kTH2D, {{200, 0., 10.}, {200, 0., 10.}}}},
      {"h_bachelorsign_Xi", "Bachelor sign;Sign;Counts", {HistType::kTH1D, {{10, -5., 5.}}}},
      {"h_ptmassdcaxyM_Xi", ";p_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});DCA_{xy} (cm)", {HistType::kTH3D, {{200, 0., 10.}, {200, 0., 10.}, {200, -.5, .5}}}},
      {"h_ptmassdcaxyAM_Xi", ";p_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});DCA_{xy} (cm)", {HistType::kTH3D, {{200, 0., 10.}, {200, 0., 10.}, {200, -.5, .5}}}},
      {"h_ptmassdcazM_Xi", ";p_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});DCA_{x} (cm)", {HistType::kTH3D, {{200, 0., 10.}, {200, 0., 10.}, {200, -.5, .5}}}},
      {"h_ptmassdcazAM_Xi", ";p_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});DCA_{x} (cm)", {HistType::kTH3D, {{200, 0., 10.}, {200, 0., 10.}, {200, -.5, .5}}}},

    }};

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {

    if (mRunNumber == bc.runNumber()) {
      return;
    }
    LOG(debug) << "Run number: " << mRunNumber << "   bc runNumber: " << bc.runNumber();

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
    return;
  }

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setFatalWhenNull(false);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

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
      h2TPCnSigma[iN] = registry.add<TH2>(Form("fTPCnSigma_%s", particlesNames[iN].data()), "n-sigma TPC", HistType::kTH2F, {{100, -5, 5, "#it{p} /Z (GeV/#it{c})"}, {200, -10., +10., Form("n#sigma_%s (a. u.)", particlesNames[iN].data())}});
    }
    auto scalers{std::get<std::shared_ptr<TH1>>(registry.add("fProcessedCascades", ";;Number of filtered cascades", HistType::kTH1D, {{nParticles + 1, -0.5, nParticles + 0.5}}))};
    scalers->GetXaxis()->SetBinLabel(1, "Processed cascades");
    for (uint32_t iS{0}; iS < particlesNames.size(); ++iS) {
      scalers->GetXaxis()->SetBinLabel(iS + 2, particlesNames[iS].data());
    }

    auto cutsOmega{std::get<std::shared_ptr<TH2>>(registry.add("h_PIDcutsOmega", ";;Invariant mass (GeV/#it{c}^{2})", HistType::kTH2D, {{6, -0.5, 5.5}, {200, 1.550, 1.750}}))};
    cutsOmega->GetXaxis()->SetBinLabel(1, "Tot #Omega");
    cutsOmega->GetXaxis()->SetBinLabel(2, "hasTof");
    cutsOmega->GetXaxis()->SetBinLabel(3, "nClusTPC");
    cutsOmega->GetXaxis()->SetBinLabel(4, "nSigmaTPCbach");
    cutsOmega->GetXaxis()->SetBinLabel(5, "nSigmaTPCprotontrack");
    cutsOmega->GetXaxis()->SetBinLabel(6, "nSigmaTPCpiontrack");

    auto cutsXi{std::get<std::shared_ptr<TH2>>(registry.add("h_PIDcutsXi", ";;Invariant mass (GeV/#it{c}^{2})", HistType::kTH2D, {{6, -0.5, 5.5}, {200, 1.250, 1.450}}))};
    cutsXi->GetXaxis()->SetBinLabel(1, "Tot #Xi");
    cutsXi->GetXaxis()->SetBinLabel(2, "hasTof");
    cutsXi->GetXaxis()->SetBinLabel(3, "nClusTPC");
    cutsXi->GetXaxis()->SetBinLabel(4, "nSigmaTPCbach");
    cutsXi->GetXaxis()->SetBinLabel(5, "nSigmaTPCprotontrack");
    cutsXi->GetXaxis()->SetBinLabel(6, "nSigmaTPCpiontrack");

    auto nClusPerTrack{std::get<std::shared_ptr<TH2>>(registry.add("h_nClusPerTrack", ";;# Clusters", HistType::kTH2D, {{3, -0.5, 2.5}, {210, 0., 210.}}))};
    nClusPerTrack->GetXaxis()->SetBinLabel(1, "Bachelor");
    nClusPerTrack->GetXaxis()->SetBinLabel(2, "Proton track");
    nClusPerTrack->GetXaxis()->SetBinLabel(3, "Pion track");

    invMassBCOmega = registry.add<TH1>("h_invariantmass_beforeCuts_Omega", "Invariant Mass (GeV/#it{c}^{2})", HistType::kTH1D, {{200, 1.609, 1.735, "Invariant Mass (GeV/#it{c}^{2})"}});
    invMassACOmega = registry.add<TH1>("h_invariantmass_afterCuts_Omega", "Invariant Mass (GeV/#it{c}^{2})", HistType::kTH1D, {{200, 1.609, 1.735, "Invariant Mass (GeV/#it{c}^{2})"}});
    invMassBCXi = registry.add<TH1>("h_invariantmass_beforeCuts_Xi", "Invariant Mass (GeV/#it{c}^{2})", HistType::kTH1D, {{200, 1.25, 1.392, "Invariant Mass (GeV/#it{c}^{2})"}});
    invMassACXi = registry.add<TH1>("h_invariantmass_afterCuts_Xi", "Invariant Mass (GeV/#it{c}^{2})", HistType::kTH1D, {{200, 1.25, 1.392, "Invariant Mass (GeV/#it{c}^{2})"}});
  }

  template <typename B, typename PR, typename PI>
  void fillQAPID(float nSigmaTPC[4], B const bachelor, PR const protonTrack, PI const pionTrack)
  {
    registry.fill(HIST("h_nClusPerTrack"), 0, bachelor.tpcNClsFound());
    registry.fill(HIST("h_nClusPerTrack"), 1, protonTrack.tpcNClsFound());
    registry.fill(HIST("h_nClusPerTrack"), 2, pionTrack.tpcNClsFound());

    registry.fill(HIST("fTPCsignal"), bachelor.sign() * bachelor.tpcInnerParam(), bachelor.tpcSignal());
    registry.fill(HIST("fTPCsignal"), protonTrack.sign() * protonTrack.tpcInnerParam(), protonTrack.tpcSignal());
    registry.fill(HIST("fTPCsignal"), pionTrack.sign() * pionTrack.tpcInnerParam(), pionTrack.tpcSignal());

    h2TPCnSigma[0]->Fill(bachelor.sign() * bachelor.tpcInnerParam(), nSigmaTPC[0]);
    h2TPCnSigma[1]->Fill(bachelor.sign() * bachelor.tpcInnerParam(), nSigmaTPC[1]);
    h2TPCnSigma[2]->Fill(protonTrack.sign() * protonTrack.tpcInnerParam(), nSigmaTPC[2]);
    h2TPCnSigma[3]->Fill(pionTrack.sign() * pionTrack.tpcInnerParam(), nSigmaTPC[3]);

    h2TPCsignal[1]->Fill(bachelor.sign() * bachelor.tpcInnerParam(), bachelor.tpcSignal());
    h2TPCsignal[2]->Fill(protonTrack.sign() * protonTrack.tpcInnerParam(), protonTrack.tpcSignal());
    h2TPCsignal[3]->Fill(pionTrack.sign() * pionTrack.tpcInnerParam(), pionTrack.tpcSignal());
  }

  template <typename TC, typename T, typename B>
  void fillCascadeDCA(TC const& trackedCascade, T const track, B const bachelor, o2::dataformats::VertexBase primaryVertex, bool isOmega)
  {
    const auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value);

    auto trackCovTrk = getTrackParCov(track);
    o2::dataformats::DCA impactParameterTrk;
    if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovTrk, mBz, 2.f, matCorr, &impactParameterTrk)) {

      if (isOmega) {
        registry.fill(HIST("h_dca_Omega"), TMath::Sqrt(impactParameterTrk.getR2()));
        registry.fill(HIST("h_dcaxy_Omega"), impactParameterTrk.getY());
        registry.fill(HIST("h_dcaz_Omega"), impactParameterTrk.getZ());
        registry.fill(HIST("h_dcavspt_Omega"), impactParameterTrk.getY(), track.pt());
        registry.fill(HIST("h_dcavsr_Omega"), impactParameterTrk.getY(), std::hypot(track.x(), track.y()));
        registry.fill(HIST("h_massvspt_Omega"), trackedCascade.omegaMass(), track.pt());
        if (bachelor.sign() < 0) {
          registry.fill(HIST("h_ptmassdcazM_Omega"), track.pt(), trackedCascade.omegaMass(), impactParameterTrk.getZ());
          registry.fill(HIST("h_ptmassdcaxyM_Omega"), track.pt(), trackedCascade.omegaMass(), impactParameterTrk.getY());
        } else if (bachelor.sign() > 0) {
          registry.fill(HIST("h_ptmassdcazAM_Omega"), track.pt(), trackedCascade.omegaMass(), impactParameterTrk.getZ());
          registry.fill(HIST("h_ptmassdcaxyAM_Omega"), track.pt(), trackedCascade.omegaMass(), impactParameterTrk.getY());
        }
      }
      registry.fill(HIST("h_dca_Xi"), TMath::Sqrt(impactParameterTrk.getR2()));
      registry.fill(HIST("h_dcaxy_Xi"), impactParameterTrk.getY());
      registry.fill(HIST("h_dcaz_Xi"), impactParameterTrk.getZ());
      registry.fill(HIST("h_dcavspt_Xi"), impactParameterTrk.getY(), track.pt());
      registry.fill(HIST("h_dcavsr_Xi"), impactParameterTrk.getY(), std::hypot(track.x(), track.y()));
      registry.fill(HIST("h_massvspt_Xi"), trackedCascade.omegaMass(), track.pt());
      if (bachelor.sign() < 0) {
        registry.fill(HIST("h_ptmassdcazM_Xi"), track.pt(), trackedCascade.omegaMass(), impactParameterTrk.getZ());
        registry.fill(HIST("h_ptmassdcaxyM_Xi"), track.pt(), trackedCascade.omegaMass(), impactParameterTrk.getY());
      } else if (bachelor.sign() > 0) {
        registry.fill(HIST("h_ptmassdcazAM_Xi"), track.pt(), trackedCascade.omegaMass(), impactParameterTrk.getZ());
        registry.fill(HIST("h_ptmassdcaxyAM_Xi"), track.pt(), trackedCascade.omegaMass(), impactParameterTrk.getY());
      }
    }
  }

  template <typename TC, typename B, typename PR, typename PI>
  void fillDauDCA(TC const& trackedCascade, B const& bachelor, PR const& protonTrack, PI const& pionTrack, o2::dataformats::VertexBase primaryVertex, bool isOmega)
  {

    const auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value);

    auto trackCovBach = getTrackParCov(bachelor);
    o2::dataformats::DCA impactParameterBach;
    if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovBach, mBz, 2.f, matCorr, &impactParameterBach)) {
      if (isOmega) {
        if (bachelor.sign() < 0) {
          registry.fill(HIST("h_bachdcaxyM_Omega"), impactParameterBach.getY());
          registry.fill(HIST("h_bachdcazM_Omega"), impactParameterBach.getZ());
        } else if (bachelor.sign() > 0) {
          registry.fill(HIST("h_bachdcaxyAM_Omega"), impactParameterBach.getY());
          registry.fill(HIST("h_bachdcazAM_Omega"), impactParameterBach.getZ());
        }
        registry.fill(HIST("h_bachdcavspt_Omega"), impactParameterBach.getY(), bachelor.pt());
        registry.fill(HIST("h_bachdcavsr_Omega"), impactParameterBach.getY(), std::hypot(trackedCascade.decayX(), trackedCascade.decayY()));
        registry.fill(HIST("h_bachelorsign_Omega"), bachelor.sign());
      }
      if (bachelor.sign() < 0) {
        registry.fill(HIST("h_bachdcaxyM_Xi"), impactParameterBach.getY());
        registry.fill(HIST("h_bachdcazM_Xi"), impactParameterBach.getZ());
      } else if (bachelor.sign() > 0) {
        registry.fill(HIST("h_bachdcaxyAM_Xi"), impactParameterBach.getY());
        registry.fill(HIST("h_bachdcazAM_Xi"), impactParameterBach.getZ());
      }
      registry.fill(HIST("h_bachdcavspt_Xi"), impactParameterBach.getY(), bachelor.pt());
      registry.fill(HIST("h_bachdcavsr_Xi"), impactParameterBach.getY(), std::hypot(trackedCascade.decayX(), trackedCascade.decayY()));
      registry.fill(HIST("h_bachelorsign_Xi"), bachelor.sign());
    }

    auto trackCovNtrack = getTrackParCov(pionTrack);
    o2::dataformats::DCA impactParameterNtrack;
    if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovNtrack, mBz, 2.f, matCorr, &impactParameterNtrack)) {
      if (isOmega) {
        registry.fill(HIST("h_ntrackdcavspt_Omega"), impactParameterNtrack.getY(), pionTrack.pt());
      }
      registry.fill(HIST("h_ntrackdcavspt_Xi"), impactParameterNtrack.getY(), pionTrack.pt());
    }

    auto trackCovPtrack = getTrackParCov(protonTrack);
    o2::dataformats::DCA impactParameterPtrack;
    if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovPtrack, mBz, 2.f, matCorr, &impactParameterPtrack)) {
      if (isOmega) {
        registry.fill(HIST("h_ptrackdcavspt_Omega"), impactParameterPtrack.getY(), protonTrack.pt());
      }
      registry.fill(HIST("h_ptrackdcavspt_Xi"), impactParameterPtrack.getY(), protonTrack.pt());
    }
  }

  void processTrackedCascadesMC(CollisionCandidatesRun3 const& collision,
                                aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& cascades,
                                aod::V0s const& v0s, TracksExtMC const& tracks,
                                soa::Join<aod::TraCascDatas, aod::McTraCascLabels> const& trackedcascdata,
                                aod::McParticles const& mcParticles, aod::BCsWithTimestamps const&)
  {
    bool keepEvent[nParticles]{false};
    bool isOmega{false};

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    registry.fill(HIST("fCollZpos"), collision.posZ());
    const auto primaryVertex = getPrimaryVertex(collision);

    for (const auto& trackedCascadeData : trackedcascdata) {
      registry.fill(HIST("h_buildermassvspt_Omega"), trackedCascadeData.mOmega(), trackedCascadeData.pt());
      registry.fill(HIST("h_buildermassvspt_Xi"), trackedCascadeData.mXi(), trackedCascadeData.pt());
    }

    for (const auto& trackedCascade : trackedCascades) {
      registry.fill(HIST("fProcessedCascades"), 0);

      isOmega = false;

      const auto track = trackedCascade.track_as<TracksExtMC>();
      const auto& casc = trackedCascade.cascade();
      const auto& bachelor = casc.bachelor_as<TracksExtMC>();
      const auto& v0 = casc.v0();
      const auto& ptrack = v0.posTrack_as<TracksExtMC>();
      const auto& ntrack = v0.negTrack_as<TracksExtMC>();
      const auto& protonTrack = bachelor.sign() > 0 ? ntrack : ptrack;
      const auto& pionTrack = bachelor.sign() > 0 ? ptrack : ntrack;

      ////Omega hypohesis -> rejecting Xi
      if (TMath::Abs(trackedCascade.xiMass() - constants::physics::MassXiMinus) > 0.005) {
        isOmega = true;
        invMassBCOmega->Fill(trackedCascade.omegaMass());
      }

      invMassBCXi->Fill(trackedCascade.xiMass());

      registry.fill(HIST("h_PIDcutsOmega"), 0, trackedCascade.omegaMass());
      registry.fill(HIST("h_PIDcutsXi"), 0, trackedCascade.xiMass());

      // if (!bachelor.hasTOF() && !ptrack.hasTOF() && !ntrack.hasTOF()) {
      //   LOG(debug) << "no TOF: " << bachelor.hasTOF() << "/" << ptrack.hasTOF() << "/" << ntrack.hasTOF();
      //   continue;
      // }

      registry.fill(HIST("h_PIDcutsOmega"), 1, trackedCascade.omegaMass());
      registry.fill(HIST("h_PIDcutsXi"), 1, trackedCascade.xiMass());

      if (protonTrack.tpcNClsFound() < cfgCutNclusTPC || pionTrack.tpcNClsFound() < cfgCutNclusTPC) {
        LOG(debug) << "no tpcNClsFound: " << bachelor.tpcNClsFound() << "/" << protonTrack.tpcNClsFound() << "/" << pionTrack.tpcNClsFound();
        continue;
      }

      registry.fill(HIST("h_PIDcutsOmega"), 2, trackedCascade.omegaMass());
      registry.fill(HIST("h_PIDcutsXi"), 2, trackedCascade.xiMass());

      // QA PID
      float nSigmaTPC[nParticles]{bachelor.tpcNSigmaKa(), bachelor.tpcNSigmaPi(), protonTrack.tpcNSigmaPr(), pionTrack.tpcNSigmaPi()};
      fillQAPID(nSigmaTPC, bachelor, protonTrack, pionTrack);

      if (isOmega) {
        h2TPCsignal[0]->Fill(bachelor.sign() * bachelor.tpcInnerParam(), bachelor.tpcSignal());
        if (bachelor.hasTPC()) {
          LOG(debug) << "TPCSignal bachelor " << bachelor.sign() << "/" << bachelor.tpcInnerParam() << "/" << bachelor.tpcSignal();
          if (nSigmaTPC[0] < cfgCutsPID->get(0u, 0u) || nSigmaTPC[0] > cfgCutsPID->get(0u, 1u)) {
            continue;
          }
        }
        keepEvent[0] = true;
        registry.fill(HIST("h_PIDcutsOmega"), 3, trackedCascade.omegaMass());
      }

      if (bachelor.hasTPC()) {
        LOG(debug) << "TPCSignal bachelor " << bachelor.sign() << "/" << bachelor.tpcInnerParam() << "/" << bachelor.tpcSignal();
        if (nSigmaTPC[1] < cfgCutsPID->get(1u, 0u) || nSigmaTPC[1] > cfgCutsPID->get(1u, 1u)) {
          continue;
        }
      }
      keepEvent[1] = true;
      registry.fill(HIST("h_PIDcutsXi"), 3, trackedCascade.xiMass());

      LOG(debug) << "TPCSignal protonTrack " << protonTrack.sign() << "/" << protonTrack.tpcInnerParam() << "/" << protonTrack.tpcSignal();
      if (nSigmaTPC[2] < cfgCutsPID->get(2u, 0u) || nSigmaTPC[2] > cfgCutsPID->get(2u, 1u)) {
        continue;
      }
      keepEvent[2] = true;

      registry.fill(HIST("h_PIDcutsOmega"), 4, trackedCascade.omegaMass());
      registry.fill(HIST("h_PIDcutsXi"), 4, trackedCascade.xiMass());

      LOG(debug) << "TPCSignal ntrack " << pionTrack.sign() << "/" << pionTrack.tpcInnerParam() << "/" << pionTrack.tpcSignal();
      if (nSigmaTPC[3] < cfgCutsPID->get(3u, 0u) || nSigmaTPC[3] > cfgCutsPID->get(3u, 1u)) {
        continue;
      }
      keepEvent[3] = true;

      registry.fill(HIST("h_PIDcutsXi"), 5, trackedCascade.xiMass());
      registry.fill(HIST("h_PIDcutsOmega"), 5, trackedCascade.omegaMass());

      invMassACXi->Fill(trackedCascade.xiMass());
      invMassACOmega->Fill(trackedCascade.omegaMass());

      for (int iDecision{0}; iDecision < 4; ++iDecision) {
        if (keepEvent[iDecision]) {
          registry.fill(HIST("fProcessedCascades"), iDecision + 1);
        }
      }

      fillCascadeDCA(trackedCascade, track, bachelor, primaryVertex, isOmega);

      LOGF(debug, "protonTrack (id: %d, pdg: %d) has mother %d", protonTrack.mcParticleId(),
           protonTrack.mcParticle().pdgCode(), protonTrack.mcParticle().has_mothers() ? protonTrack.mcParticle().mothersIds()[0] : -1);
      LOGF(debug, "pionTrack (id: %d, pdg: %d) has mother %d", pionTrack.mcParticleId(),
           pionTrack.mcParticle().pdgCode(), pionTrack.mcParticle().has_mothers() ? pionTrack.mcParticle().mothersIds()[0] : -1);

      LOG(debug) << "bachelor with PDG code: " << bachelor.mcParticle().pdgCode() << ". Charge: " << bachelor.sign();
      if (ptrack.mcParticle().has_mothers() && pionTrack.mcParticle().has_mothers() &&
          protonTrack.mcParticle().mothersIds()[0] == pionTrack.mcParticle().mothersIds()[0]) {
        const auto v0part = protonTrack.mcParticle().mothers_as<aod::McParticles>()[0];
        LOG(debug) << "v0 with PDG code: " << v0part.pdgCode();
        if (v0part.has_mothers() && bachelor.mcParticle().has_mothers() &&
            v0part.mothersIds()[0] == bachelor.mcParticle().mothersIds()[0]) {
          int mcid = v0part.mothersIds()[0];
          for (const auto& trackedCascadeData : trackedcascdata) {
            if (trackedCascadeData.mcParticleId() == mcid) {
              if (isOmega) {
                registry.fill(HIST("h_massvsmass_Omega"), trackedCascade.omegaMass(), trackedCascadeData.mOmega());
              } else {
                registry.fill(HIST("h_massvsmass_Xi"), trackedCascade.omegaMass(), trackedCascadeData.mOmega());
              }
              break;
            }
            LOG(debug) << "cascade with PDG code: " << v0part.mothers_as<aod::McParticles>()[0].pdgCode();
          }
        } else {
          LOG(debug) << "Rejecting particle.";
          continue;
        }
      }
      fillDauDCA(trackedCascade, bachelor, protonTrack, pionTrack, primaryVertex, isOmega);
    }
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processTrackedCascadesMC, "process cascades from strangeness tracking: MC analysis", true);

  void processTrackedCascadesData(CollisionCandidatesRun3 const& collision,
                                  aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& cascades,
                                  aod::V0s const& v0s, TracksExtData const& tracks,
                                  aod::BCsWithTimestamps const&)
  {

    bool keepEvent[nParticles]{false};

    bool isOmega{false};

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    registry.fill(HIST("fCollZpos"), collision.posZ());
    const auto primaryVertex = getPrimaryVertex(collision);

    for (const auto& trackedCascade : trackedCascades) {
      registry.fill(HIST("fProcessedCascades"), 0);

      isOmega = false;

      const auto track = trackedCascade.track_as<TracksExtData>();
      const auto& casc = trackedCascade.cascade();
      const auto& bachelor = casc.bachelor_as<TracksExtData>();
      const auto& v0 = casc.v0();
      const auto& ptrack = v0.posTrack_as<TracksExtData>();
      const auto& ntrack = v0.negTrack_as<TracksExtData>();
      const auto& protonTrack = bachelor.sign() > 0 ? ntrack : ptrack;
      const auto& pionTrack = bachelor.sign() > 0 ? ptrack : ntrack;

      ////Omega hypohesis -> rejecting Xi
      if (TMath::Abs(trackedCascade.xiMass() - constants::physics::MassXiMinus) > 0.005) {
        isOmega = true;
        invMassBCOmega->Fill(trackedCascade.omegaMass());
      }

      invMassBCXi->Fill(trackedCascade.xiMass());

      registry.fill(HIST("h_PIDcutsXi"), 0, trackedCascade.xiMass());
      registry.fill(HIST("h_PIDcutsOmega"), 0, trackedCascade.omegaMass());

      // if (!bachelor.hasTOF() && !ptrack.hasTOF() && !ntrack.hasTOF() ) {
      //   LOG(debug)<< "no TOF: "<<bachelor.hasTOF()<<"/"<<ptrack.hasTOF()<<"/"<<ntrack.hasTOF();
      //   continue;
      // }

      registry.fill(HIST("h_PIDcutsXi"), 1, trackedCascade.xiMass());
      registry.fill(HIST("h_PIDcutsOmega"), 1, trackedCascade.omegaMass());

      if (protonTrack.tpcNClsFound() < cfgCutNclusTPC || pionTrack.tpcNClsFound() < cfgCutNclusTPC) {
        LOG(debug) << "no tpcNClsFound: " << bachelor.tpcNClsFound() << "/" << protonTrack.tpcNClsFound() << "/" << pionTrack.tpcNClsFound();
        continue;
      }
      registry.fill(HIST("h_PIDcutsXi"), 2, trackedCascade.xiMass());
      registry.fill(HIST("h_PIDcutsOmega"), 2, trackedCascade.omegaMass());

      // QA PID
      float nSigmaTPC[nParticles]{bachelor.tpcNSigmaKa(), bachelor.tpcNSigmaPi(), protonTrack.tpcNSigmaPr(), pionTrack.tpcNSigmaPi()};
      fillQAPID(nSigmaTPC, bachelor, protonTrack, pionTrack);

      if (isOmega) {
        h2TPCsignal[0]->Fill(bachelor.sign() * bachelor.tpcInnerParam(), bachelor.tpcSignal());
        if (bachelor.hasTPC()) {
          LOG(debug) << "TPCSignal bachelor " << bachelor.sign() << "/" << bachelor.tpcInnerParam() << "/" << bachelor.tpcSignal();
          if (nSigmaTPC[0] < cfgCutsPID->get(0u, 0u) || nSigmaTPC[0] > cfgCutsPID->get(0u, 1u)) {
            continue;
          }
        }
        keepEvent[0] = true;
        registry.fill(HIST("h_PIDcutsOmega"), 3, trackedCascade.omegaMass());
      }

      if (bachelor.hasTPC()) {
        LOG(debug) << "TPCSignal bachelor " << bachelor.sign() << "/" << bachelor.tpcInnerParam() << "/" << bachelor.tpcSignal();
        if (nSigmaTPC[1] < cfgCutsPID->get(1u, 0u) || nSigmaTPC[1] > cfgCutsPID->get(1u, 1u)) {
          continue;
        }
      }
      keepEvent[1] = true;
      registry.fill(HIST("h_PIDcutsXi"), 3, trackedCascade.xiMass());

      LOG(debug) << "TPCSignal protonTrack " << protonTrack.sign() << "/" << protonTrack.tpcInnerParam() << "/" << protonTrack.tpcSignal();
      if (nSigmaTPC[2] < cfgCutsPID->get(2u, 0u) || nSigmaTPC[2] > cfgCutsPID->get(2u, 1u)) {
        continue;
      }
      keepEvent[2] = true;

      registry.fill(HIST("h_PIDcutsXi"), 4, trackedCascade.xiMass());
      registry.fill(HIST("h_PIDcutsOmega"), 4, trackedCascade.omegaMass());

      LOG(debug) << "TPCSignal ntrack " << pionTrack.sign() << "/" << pionTrack.tpcInnerParam() << "/" << pionTrack.tpcSignal();
      if (nSigmaTPC[3] < cfgCutsPID->get(3u, 0u) || nSigmaTPC[3] > cfgCutsPID->get(3u, 1u)) {
        continue;
      }
      keepEvent[3] = true;

      registry.fill(HIST("h_PIDcutsXi"), 5, trackedCascade.xiMass());
      registry.fill(HIST("h_PIDcutsOmega"), 5, trackedCascade.omegaMass());

      invMassACXi->Fill(trackedCascade.xiMass());
      invMassACOmega->Fill(trackedCascade.omegaMass());

      for (int iDecision{0}; iDecision < 4; ++iDecision) {
        if (keepEvent[iDecision]) {
          registry.fill(HIST("fProcessedCascades"), iDecision + 1);
        }
      }

      fillCascadeDCA(trackedCascade, track, bachelor, primaryVertex, isOmega);
      fillDauDCA(trackedCascade, bachelor, protonTrack, pionTrack, primaryVertex, isOmega);
    }
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processTrackedCascadesData, "process cascades from strangeness tracking: Data analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NonPromptCascadeTask>(cfgc)};
}
