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

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
// #include "PWGHF/Core/PDG.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "PWGLF/DataModel/LFNonPromptCascadeTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
struct NPCascCandidate {
  int64_t mcParticleId;
  int64_t trackGlobID;
  int64_t trackITSID;
  int64_t collisionID;
  float matchingChi2;
  float deltaPt;
  float itsClusSize;
  bool hasReassociatedCluster;
  bool isGoodMatch;
  bool isGoodCascade;
  int pdgCodeMom;
  int pdgCodeITStrack;
  bool isFromBeauty;
  bool isFromCharm;
  uint16_t pvContributors;
  float pvTimeResolution;
  float pvX;
  float pvY;
  float pvZ;
  float cascPt;
  float cascEta;
  float cascPhi;
  float protonPt;
  float protonEta;
  float pionPt;
  float pionEta;
  float bachPt;
  float bachEta;
  float cascDCAxy;
  float cascDCAz;
  float protonDCAxy;
  float protonDCAz;
  float pionDCAxy;
  float pionDCAz;
  float bachDCAxy;
  float bachDCAz;
  float casccosPA;
  float v0cosPA;
  double massXi;
  double massOmega;
  double massV0;
  float cascRadius;
  float v0radius;
  float cascLength;
  float v0length;
  int cascNClusITS;
  int protonNClusITS;
  int pionNClusITS;
  int bachKaonNClusITS;
  int bachPionNClusITS;
  int protonNClusTPC;
  int pionNClusTPC;
  int bachKaonNClusTPC;
  int bachPionNClusTPC;
  float protonTPCNSigma;
  float pionTPCNSigma;
  float bachKaonTPCNSigma;
  float bachPionTPCNSigma;
  bool protonHasTOF;
  bool pionHasTOF;
  bool bachKaonHasTOF;
  bool bachPionHasTOF;
  float protonTOFNSigma;
  float pionTOFNSigma;
  float bachKaonTOFNSigma;
  float bachPionTOFNSigma;
};

struct motherDCA {
  float DCAxy;
  float DCAz;
};

struct daughtersDCA {
  float bachDCAxy;
  float bachDCAz;
  float protonDCAxy;
  float protonDCAz;
  float pionDCAxy;
  float pionDCAz;
};

std::array<bool, 2> isFromHF(auto& particle)
{
  bool fromBeauty = false;
  bool fromCharm = false;
  if (particle.has_mothers()) {
    auto mom = particle.template mothers_as<aod::McParticles>()[0];
    int pdgCodeMom = mom.pdgCode();
    fromBeauty = std::abs(pdgCodeMom) / 5000 == 1 || std::abs(pdgCodeMom) / 500 == 1 || std::abs(pdgCodeMom) == 5;
    fromCharm = std::abs(pdgCodeMom) / 4000 == 1 || std::abs(pdgCodeMom) / 400 == 1 || std::abs(pdgCodeMom) == 4;
    while (mom.has_mothers()) {
      const auto grandma = mom.template mothers_as<aod::McParticles>()[0];
      int pdgCodeGrandma = std::abs(grandma.pdgCode());
      fromBeauty = fromBeauty || (pdgCodeGrandma / 5000 == 1 || pdgCodeGrandma / 500 == 1 || pdgCodeGrandma == 5);
      fromCharm = fromCharm || (pdgCodeGrandma / 4000 == 1 || pdgCodeGrandma / 400 == 1 || pdgCodeGrandma == 4);
      mom = grandma;
    }
  }
  return {fromBeauty, fromCharm};
}

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

std::shared_ptr<TH1> invMassBCV0;
std::shared_ptr<TH1> invMassACV0;

std::vector<NPCascCandidate> candidates;

} // namespace

struct NonPromptCascadeTask {

  Produces<o2::aod::NPCascTable> NPCTable;
  Produces<o2::aod::NPCascTableMC> NPCTableMC;
  Produces<o2::aod::NPCascTableGen> NPCTableGen;

  using TracksExtData = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullKa, aod::pidTOFFullPi, aod::pidTOFFullPr>;
  using TracksExtMC = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::McTrackLabels, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullKa, aod::pidTOFFullPi, aod::pidTOFFullPr>;
  using CollisionCandidatesRun3 = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionCandidatesRun3MC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> propToDCA{"propToDCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> cfgGRPmagPath{"cfgGRPmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Configurable<float> cfgCutNclusTPC{"cfgCutNclusTPC", 70, "Minimum number of TPC clusters"};
  Configurable<LabeledArray<float>> cfgCutsPID{"particlesCutsPID", {cutsPID[0], nParticles, nCutsPID, particlesNames, cutsNames}, "Nuclei PID selections"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  float mBz = 0.f;

  HistogramRegistry registry{
    "registry",
    {
      {"h_PV_x", "Primary vertex x;x (cm)", {HistType::kTH1D, {{100, -1., 1.}}}},
      {"h_PV_y", "Primary vertex y;y (cm)", {HistType::kTH1D, {{100, -1., 1.}}}},
      {"h_PV_z", "Primary vertex z;z (cm)", {HistType::kTH1D, {{100, -1., 1.}}}},

      {"h_dca_Omega", "DCA;DCA (cm)", {HistType::kTH1D, {{200, 0., .5}}}},
      {"h_dcaxy_Omega", "DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"h_dcaz_Omega", "DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"h_bachdcaxyM_Omega", "Bachelor DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcaxyAM_Omega", "Bachelor DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcazM_Omega", "Bachelor DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcazAM_Omega", "Bachelor DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_dcavspt_Omega", "DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{100, -0.1, 0.1}, {50, 0., 10.}}}},
      {"h_bachdcavspt_Omega", "Bachelor DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {50, 0., 10.}}}},
      {"h_bachdcavsr_Omega", "Bachelor DCA vs R (cm);DCA (cm);R (cm)", {HistType::kTH2D, {{200, -1., 1.}, {50, 0., 30.}}}},
      {"h_ntrackdcavspt_Omega", "N track DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {50, 0., 10.}}}},
      {"h_ptrackdcavspt_Omega", "P track DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {50, 0., 10.}}}},
      {"h_dcavsr_Omega", "DCA vs R;DCA (cm);R (cm)", {HistType::kTH2D, {{200, -.5, .5}, {200, 0., 5.}}}},
      {"h_massvspt_Omega", "Mass vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{125, 1.650, 1.700}, {50, 0., 10.}}}},
      {"h_buildermassvspt_Omega", "Mass (from builder) vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{125, 1.650, 1.700}, {50, 0., 10.}}}},
      {"h_massvsmass_Omega", "Mass vs mass;Mass (GeV/#it{c}^{2});Mass (GeV/#it{c}^{2})", {HistType::kTH2D, {{125, 1.650, 1.700}, {125, 1.650, 1.700}}}},
      {"h_bachelorsign_Omega", "Bachelor sign;Sign;Counts", {HistType::kTH1D, {{6, -3., 3.}}}},

      {"h_dca_Xi", "DCA;DCA (cm)", {HistType::kTH1D, {{200, 0., .5}}}},
      {"h_dcaxy_Xi", "DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"h_dcaz_Xi", "DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"h_bachdcaxyM_Xi", "Bachelor DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcaxyAM_Xi", "Bachelor DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcazM_Xi", "Bachelor DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcazAM_Xi", "Bachelor DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_dcavspt_Xi", "DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{100, -0.1, 0.1}, {50, 0., 10.}}}},
      {"h_bachdcavspt_Xi", "Bachelor DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {50, 0., 10.}}}},
      {"h_bachdcavsr_Xi", "Bachelor DCA vs R (cm);DCA (cm);R (cm)", {HistType::kTH2D, {{200, -1., 1.}, {50, 0., 30.}}}},
      {"h_ntrackdcavspt_Xi", "N track DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {50, 0., 10.}}}},
      {"h_ptrackdcavspt_Xi", "P track DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {50, 0., 10.}}}},
      {"h_dcavsr_Xi", "DCA vs R;DCA (cm);R (cm)", {HistType::kTH2D, {{200, -.5, .5}, {200, 0., 5.}}}},
      {"h_massvspt_Xi", "Mass vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{125, 1.296, 1.346}, {50, 0., 10.}}}},
      {"h_buildermassvspt_Xi", "Mass (from builder) vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{125, 1.296, 1.346}, {50, 0., 10.}}}},
      {"h_massvsmass_Xi", "Mass vs mass;Mass (GeV/#it{c}^{2});Mass (GeV/#it{c}^{2})", {HistType::kTH2D, {{125, 1.296, 1.346}, {125, 1.296, 1.346}}}},
      {"h_bachelorsign_Xi", "Bachelor sign;Sign;Counts", {HistType::kTH1D, {{6, -3., 3.}}}},

      {"h_massvspt_V0", "Mass vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{125, 1.090, 1.140}, {50, 0., 10.}}}},

    }};

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();

    if (o2::parameters::GRPMagField* grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(cfgGRPmagPath, mRunNumber)) {
      o2::base::Propagator::initFieldFromGRP(grpmag);
      mBz = static_cast<float>(grpmag->getNominalL3Field());
    }

    if (static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value) == o2::base::Propagator::MatCorrType::USEMatCorrLUT) {
      auto* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForRun<o2::base::MatLayerCylSet>("GLO/Param/MatLUT", mRunNumber));
      o2::base::Propagator::Instance(true)->setMatLUT(lut);
    }
  }

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setFatalWhenNull(true);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    std::vector<double> ptBinning = {0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    auto cutsOmega{std::get<std::shared_ptr<TH2>>(registry.add("h_PIDcutsOmega", ";;Invariant mass (GeV/#it{c}^{2})", HistType::kTH2D, {{6, -0.5, 5.5}, {125, 1.650, 1.700}}))};
    cutsOmega->GetXaxis()->SetBinLabel(1, "Tot #Omega");
    cutsOmega->GetXaxis()->SetBinLabel(2, "hasTof");
    cutsOmega->GetXaxis()->SetBinLabel(3, "nClusTPC");
    cutsOmega->GetXaxis()->SetBinLabel(4, "nSigmaTPCbach");
    cutsOmega->GetXaxis()->SetBinLabel(5, "nSigmaTPCprotontrack");
    cutsOmega->GetXaxis()->SetBinLabel(6, "nSigmaTPCpiontrack");

    auto cutsXi{std::get<std::shared_ptr<TH2>>(registry.add("h_PIDcutsXi", ";;Invariant mass (GeV/#it{c}^{2})", HistType::kTH2D, {{6, -0.5, 5.5}, {125, 1.296, 1.346}}))};
    cutsXi->GetXaxis()->SetBinLabel(1, "Tot #Xi");
    cutsXi->GetXaxis()->SetBinLabel(2, "hasTof");
    cutsXi->GetXaxis()->SetBinLabel(3, "nClusTPC");
    cutsXi->GetXaxis()->SetBinLabel(4, "nSigmaTPCbach");
    cutsXi->GetXaxis()->SetBinLabel(5, "nSigmaTPCprotontrack");
    cutsXi->GetXaxis()->SetBinLabel(6, "nSigmaTPCpiontrack");

    invMassBCV0 = registry.add<TH1>("h_invariantmass_beforeCuts_V0", "Invariant Mass (GeV/#it{c}^{2})", HistType::kTH1D, {{125, 1.090, 1.140, "Invariant Mass (GeV/#it{c}^{2})"}});
    invMassACV0 = registry.add<TH1>("h_invariantmass_afterCuts_V0", "Invariant Mass (GeV/#it{c}^{2})", HistType::kTH1D, {{125, 1.090, 1.140, "Invariant Mass (GeV/#it{c}^{2})"}});
  }

  template <typename T, typename PR, typename PI>
  void fillCascadeDCA(T const track, PR const& protonTrack, PI const& pionTrack, o2::dataformats::VertexBase primaryVertex, bool isOmega, motherDCA& motherDCA)
  {
    const auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value);
    auto trackCovTrk = getTrackParCov(track);
    o2::dataformats::DCA impactParameterTrk{-999.f, -999.f};

    if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovTrk, mBz, 2.f, matCorr, &impactParameterTrk)) {
      if (protonTrack.hasTPC() && pionTrack.hasTPC()) {
        if (isOmega) {
          registry.fill(HIST("h_dca_Omega"), std::sqrt(impactParameterTrk.getR2()));
          registry.fill(HIST("h_dcaxy_Omega"), impactParameterTrk.getY());
          registry.fill(HIST("h_dcaz_Omega"), impactParameterTrk.getZ());
          registry.fill(HIST("h_dcavspt_Omega"), impactParameterTrk.getY(), track.pt());
          registry.fill(HIST("h_dcavsr_Omega"), impactParameterTrk.getY(), std::hypot(track.x(), track.y()));
        }
      }

      if (protonTrack.hasTPC() && pionTrack.hasTPC()) {
        registry.fill(HIST("h_dca_Xi"), std::sqrt(impactParameterTrk.getR2()));
        registry.fill(HIST("h_dcaxy_Xi"), impactParameterTrk.getY());
        registry.fill(HIST("h_dcaz_Xi"), impactParameterTrk.getZ());
        registry.fill(HIST("h_dcavspt_Xi"), impactParameterTrk.getY(), track.pt());
        registry.fill(HIST("h_dcavsr_Xi"), impactParameterTrk.getY(), std::hypot(track.x(), track.y()));
      }
    }
    motherDCA.DCAxy = impactParameterTrk.getY();
    motherDCA.DCAz = impactParameterTrk.getZ();
  }

  template <typename TC, typename B, typename PR, typename PI>
  void fillDauDCA(TC const& trackedCascade, B const& bachelor, PR const& protonTrack, PI const& pionTrack, o2::dataformats::VertexBase primaryVertex, bool isOmega, daughtersDCA& dDCA)
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
    o2::dataformats::DCA impactParameterPiontrack;
    if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovNtrack, mBz, 2.f, matCorr, &impactParameterPiontrack)) {
      if (isOmega) {
        registry.fill(HIST("h_ntrackdcavspt_Omega"), impactParameterPiontrack.getY(), pionTrack.pt());
      }
      registry.fill(HIST("h_ntrackdcavspt_Xi"), impactParameterPiontrack.getY(), pionTrack.pt());
    }

    auto trackCovPtrack = getTrackParCov(protonTrack);
    o2::dataformats::DCA impactParameterProtontrack;
    if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovPtrack, mBz, 2.f, matCorr, &impactParameterProtontrack)) {
      if (isOmega) {
        registry.fill(HIST("h_ptrackdcavspt_Omega"), impactParameterProtontrack.getY(), protonTrack.pt());
      }
      registry.fill(HIST("h_ptrackdcavspt_Xi"), impactParameterProtontrack.getY(), protonTrack.pt());
    }

    dDCA.bachDCAxy = impactParameterBach.getY();
    dDCA.bachDCAz = impactParameterBach.getZ();
    dDCA.protonDCAxy = impactParameterProtontrack.getY();
    dDCA.protonDCAz = impactParameterProtontrack.getZ();
    dDCA.pionDCAxy = impactParameterPiontrack.getY();
    dDCA.pionDCAz = impactParameterPiontrack.getZ();
  }

  template <typename TrackType, typename CollisionType>
  void fillCandidatesVector(CollisionType const&, auto const& trackedCascades)
  {
    candidates.clear();
    for (const auto& trackedCascade : trackedCascades) {

      auto collision = trackedCascade.template collision_as<CollisionType>();
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      const auto primaryVertex = getPrimaryVertex(collision);

      o2::vertexing::DCAFitterN<2> df2;
      df2.setBz(mBz);
      df2.setPropagateToPCA(propToDCA);
      df2.setMaxR(maxR);
      df2.setMaxDZIni(maxDZIni);
      df2.setMinParamChange(minParamChange);
      df2.setMinRelChi2Change(minRelChi2Change);
      df2.setUseAbsDCA(useAbsDCA);

      const auto& track = trackedCascade.template track_as<TrackType>();
      const auto& ITStrack = trackedCascade.template itsTrack_as<TrackType>();
      const auto& casc = trackedCascade.cascade();
      const auto& bachelor = casc.template bachelor_as<TrackType>();
      const auto& v0 = casc.v0();
      const auto& ptrack = v0.template posTrack_as<TrackType>();
      const auto& ntrack = v0.template negTrack_as<TrackType>();
      const auto& protonTrack = bachelor.sign() > 0 ? ntrack : ptrack;
      const auto& pionTrack = bachelor.sign() > 0 ? ptrack : ntrack;
      bool hasReassociatedClusters = (track.itsNCls() != ITStrack.itsNCls());

      std::array<std::array<float, 3>, 2> momenta;
      std::array<float, 3> v0Pos;
      std::array<double, 2> masses;

      // track propagation
      o2::track::TrackParCov trackParCovV0;
      o2::track::TrackPar trackParV0;
      o2::track::TrackPar trackParBachelor;
      std::array<float, 3> cascadeMomentum;

      float cascCpa = -1;
      float v0Cpa = -1;
      if (df2.process(getTrackParCov(pionTrack), getTrackParCov(protonTrack))) {
        trackParCovV0 = df2.createParentTrackParCov(0); // V0 track retrieved from p and pi daughters
        v0Pos = {trackParCovV0.getX(), trackParCovV0.getY(), trackParCovV0.getZ()};
        if (df2.process(trackParCovV0, getTrackParCov(bachelor))) {
          trackParV0 = df2.getTrackParamAtPCA(0);
          trackParBachelor = df2.getTrackParamAtPCA(1);
          trackParV0.getPxPyPzGlo(momenta[0]);       // getting the V0 momentum
          trackParBachelor.getPxPyPzGlo(momenta[1]); // getting the bachelor momentum
          df2.createParentTrackParCov().getPxPyPzGlo(cascadeMomentum);
          std::array<float, 3> pvPos = {primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()};
          cascCpa = RecoDecay::cpa(pvPos, df2.getPCACandidate(), cascadeMomentum);
          v0Cpa = RecoDecay::cpa(pvPos, df2.getPCACandidate(), momenta[0]);
        } else {
          continue;
        }
      } else {
        continue;
      }
      float deltaPtITSCascade = std::hypot(cascadeMomentum[0], cascadeMomentum[1]) - ITStrack.pt();

      // PV
      registry.fill(HIST("h_PV_x"), primaryVertex.getX());
      registry.fill(HIST("h_PV_y"), primaryVertex.getY());
      registry.fill(HIST("h_PV_z"), primaryVertex.getZ());

      // Omega
      masses = {o2::constants::physics::MassLambda0, o2::constants::physics::MassKPlus};
      const auto massOmega = RecoDecay::m(momenta, masses);

      // Xi
      masses = {o2::constants::physics::MassLambda0, o2::constants::physics::MassPiPlus};
      const auto massXi = RecoDecay::m(momenta, masses);

      // Lambda
      masses = {o2::constants::physics::MassProton, o2::constants::physics::MassPiMinus};
      momenta[0] = {protonTrack.px(), protonTrack.py(), protonTrack.pz()};
      momenta[1] = {pionTrack.px(), pionTrack.py(), pionTrack.pz()};
      const auto v0mass = RecoDecay::m(momenta, masses);

      //// Omega hypohesis -> rejecting Xi, we don't do it in the MC as we can identify the particle with the MC truth
      bool isOmega{std::abs(massXi - constants::physics::MassXiMinus) > 0.005};

      std::array<bool, 2> fromHF{false, false};
      bool isGoodMatch{false}, isGoodCascade{false};
      int itsTrackPDG{0}, pdgCodeMom{0};
      int64_t mcParticleID{-1};

      if constexpr (TrackType::template contains<aod::McTrackLabels>()) {
        if (protonTrack.mcParticle().has_mothers() && pionTrack.mcParticle().has_mothers() && bachelor.mcParticle().has_mothers()) {
          if (protonTrack.mcParticle().mothersIds()[0] == pionTrack.mcParticle().mothersIds()[0]) {
            const auto v0part = protonTrack.mcParticle().template mothers_first_as<aod::McParticles>();
            if (std::abs(v0part.pdgCode()) == 3122 && v0part.has_mothers()) {
              const auto motherV0 = v0part.template mothers_as<aod::McParticles>()[0];
              if (v0part.mothersIds()[0] == bachelor.mcParticle().mothersIds()[0]) {
                if (std::abs(motherV0.pdgCode()) == 3312 || std::abs(motherV0.pdgCode()) == 3334) {
                  isGoodCascade = true;
                  isOmega = (std::abs(motherV0.pdgCode()) == 3334);
                  fromHF = isFromHF(motherV0);
                  mcParticleID = v0part.mothersIds()[0];
                }
              }
            }
          }
        }
        isGoodMatch = ((mcParticleID == ITStrack.mcParticleId())) ? true : false;

        if (isGoodMatch) {
          pdgCodeMom = track.mcParticle().has_mothers() ? track.mcParticle().template mothers_as<aod::McParticles>()[0].pdgCode() : 0;
        }
        itsTrackPDG = ITStrack.has_mcParticle() ? ITStrack.mcParticle().pdgCode() : 0;
      }

      invMassBCV0->Fill(v0mass);

      registry.fill(HIST("h_PIDcutsXi"), 0, massXi);
      registry.fill(HIST("h_PIDcutsOmega"), 0, massOmega);

      int bachKaonNClusTPC = -1;
      int bachPionNClusTPC = -1;
      int bachKaonNClusITS = -1;
      int bachPionNClusITS = -1;
      if (isOmega) {
        bachKaonNClusTPC = bachelor.tpcNClsFound();
        bachKaonNClusITS = bachelor.itsNCls();
      }
      bachPionNClusTPC = bachelor.tpcNClsFound(); /// by default cascade = Xi
      bachPionNClusITS = bachelor.itsNCls();      /// by default cascade = Xi

      bool bachKaonHasTOF = 0;
      bool bachPionHasTOF = 0;
      if (isOmega) {
        bachKaonHasTOF = bachelor.hasTOF();
      }
      bachPionHasTOF = bachelor.hasTOF();

      registry.fill(HIST("h_PIDcutsXi"), 1, massXi);
      registry.fill(HIST("h_PIDcutsOmega"), 1, massOmega);

      if (protonTrack.tpcNClsFound() < cfgCutNclusTPC || pionTrack.tpcNClsFound() < cfgCutNclusTPC) {
        LOG(debug) << "no tpcNClsFound: " << bachelor.tpcNClsFound() << "/" << protonTrack.tpcNClsFound() << "/" << pionTrack.tpcNClsFound();
        continue;
      }
      registry.fill(HIST("h_PIDcutsXi"), 2, massXi);
      registry.fill(HIST("h_PIDcutsOmega"), 2, massOmega);

      // QA PID
      float nSigmaTPC[nParticles]{bachelor.tpcNSigmaKa(), bachelor.tpcNSigmaPi(), protonTrack.tpcNSigmaPr(), pionTrack.tpcNSigmaPi()};

      bool isBachelorSurvived = false;
      if (isOmega) {
        if (bachelor.hasTPC()) {
          LOG(debug) << "TPCSignal bachelor " << bachelor.sign() << "/" << bachelor.tpcInnerParam() << "/" << bachelor.tpcSignal();
          if (nSigmaTPC[0] > cfgCutsPID->get(0u, 0u) && nSigmaTPC[0] < cfgCutsPID->get(0u, 1u)) {
            registry.fill(HIST("h_PIDcutsOmega"), 3, massOmega);
            isBachelorSurvived = true;
          }
        }
      }

      if (bachelor.hasTPC()) {
        LOG(debug) << "TPCSignal bachelor " << bachelor.sign() << "/" << bachelor.tpcInnerParam() << "/" << bachelor.tpcSignal();
        if (nSigmaTPC[1] > cfgCutsPID->get(1u, 0u) && nSigmaTPC[1] < cfgCutsPID->get(1u, 1u)) {
          registry.fill(HIST("h_PIDcutsXi"), 3, massXi);
          isBachelorSurvived = true;
        }
      }

      if (!isBachelorSurvived) {
        continue;
      }

      LOG(debug) << "TPCSignal protonTrack " << protonTrack.sign() << "/" << protonTrack.tpcInnerParam() << "/" << protonTrack.tpcSignal();
      if (nSigmaTPC[2] < cfgCutsPID->get(2u, 0u) || nSigmaTPC[2] > cfgCutsPID->get(2u, 1u)) {
        continue;
      }

      if (isOmega) {
        registry.fill(HIST("h_PIDcutsOmega"), 4, massOmega);
      }
      registry.fill(HIST("h_PIDcutsXi"), 4, massXi);

      LOG(debug) << "TPCSignal ntrack " << pionTrack.sign() << "/" << pionTrack.tpcInnerParam() << "/" << pionTrack.tpcSignal();
      if (nSigmaTPC[3] < cfgCutsPID->get(3u, 0u) || nSigmaTPC[3] > cfgCutsPID->get(3u, 1u)) {
        continue;
      }

      if (isOmega) {
        registry.fill(HIST("h_PIDcutsOmega"), 5, massOmega);
        registry.fill(HIST("h_massvspt_Omega"), massOmega, track.pt());
      }

      registry.fill(HIST("h_PIDcutsXi"), 5, massXi);
      registry.fill(HIST("h_massvspt_Xi"), massXi, track.pt());

      invMassACV0->Fill(v0mass);
      registry.fill(HIST("h_massvspt_V0"), v0mass, track.pt());

      motherDCA motherDCA;
      fillCascadeDCA(track, protonTrack, pionTrack, primaryVertex, isOmega, motherDCA);
      daughtersDCA dDCA;
      fillDauDCA(trackedCascade, bachelor, protonTrack, pionTrack, primaryVertex, isOmega, dDCA);

      candidates.emplace_back(NPCascCandidate{mcParticleID, track.globalIndex(), ITStrack.globalIndex(), trackedCascade.collisionId(), trackedCascade.matchingChi2(), deltaPtITSCascade, trackedCascade.itsClsSize(), hasReassociatedClusters, isGoodMatch, isGoodCascade, pdgCodeMom, itsTrackPDG, fromHF[0], fromHF[1],
                                              collision.numContrib(), collision.collisionTimeRes(), primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                                              track.pt(), track.eta(), track.phi(),
                                              protonTrack.pt(), protonTrack.eta(), pionTrack.pt(), pionTrack.eta(), bachelor.pt(), bachelor.eta(),
                                              motherDCA.DCAxy, motherDCA.DCAz, dDCA.protonDCAxy, dDCA.protonDCAz, dDCA.pionDCAxy, dDCA.pionDCAz, dDCA.bachDCAxy, dDCA.bachDCAz,
                                              cascCpa, v0Cpa,
                                              massXi, massOmega, v0mass,
                                              std::hypot(trackedCascade.decayX(), trackedCascade.decayY()), std::hypot(v0Pos[0], v0Pos[1]), std::hypot(trackedCascade.decayX(), trackedCascade.decayY(), trackedCascade.decayZ()), std::hypot(v0Pos[0], v0Pos[1], v0Pos[2]),
                                              track.itsNCls(), protonTrack.itsNCls(), pionTrack.itsNCls(), bachKaonNClusITS, bachPionNClusITS, protonTrack.tpcNClsFound(), pionTrack.tpcNClsFound(), bachKaonNClusTPC, bachPionNClusTPC,
                                              protonTrack.tpcNSigmaPr(), pionTrack.tpcNSigmaPi(), bachelor.tpcNSigmaKa(), bachelor.tpcNSigmaPi(),
                                              protonTrack.hasTOF(), pionTrack.hasTOF(), bachKaonHasTOF, bachPionHasTOF,
                                              protonTrack.tofNSigmaPr(), pionTrack.tofNSigmaPi(), bachelor.tofNSigmaKa(), bachelor.tofNSigmaPi()});
    }
  }

  void processTrackedCascadesMC(CollisionCandidatesRun3MC const& collisions,
                                aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& /*cascades*/,
                                aod::V0s const& /*v0s*/, TracksExtMC const& /*tracks*/,
                                aod::McParticles const& mcParticles, aod::McCollisions const&, aod::BCsWithTimestamps const&)
  {

    fillCandidatesVector<TracksExtMC>(collisions, trackedCascades);

    for (size_t i = 0; i < candidates.size(); ++i) {

      if (candidates[i].mcParticleId < 0) {
        continue;
      }
      auto particle = mcParticles.iteratorAt(candidates[i].mcParticleId);
      auto& c = candidates[i];
      auto mcCollision = particle.mcCollision_as<aod::McCollisions>();
      auto label = collisions.iteratorAt(c.collisionID);

      NPCTableMC(c.matchingChi2, c.deltaPt, c.itsClusSize, c.hasReassociatedCluster, c.isGoodMatch, c.isGoodCascade, c.pdgCodeMom, c.pdgCodeITStrack, c.isFromBeauty, c.isFromCharm,
                 c.pvContributors, c.pvTimeResolution, c.pvX, c.pvY, c.pvZ,
                 c.cascPt, c.cascEta, c.cascPhi,
                 c.protonPt, c.protonEta, c.pionPt, c.pionEta, c.bachPt, c.bachEta,
                 c.cascDCAxy, c.cascDCAz, c.protonDCAxy, c.protonDCAz, c.pionDCAxy, c.pionDCAz, c.bachDCAxy, c.bachDCAz,
                 c.casccosPA, c.v0cosPA,
                 c.massXi, c.massOmega, c.massV0,
                 c.cascRadius, c.v0radius, c.cascLength, c.v0length,
                 c.cascNClusITS, c.protonNClusITS, c.pionNClusITS, c.bachKaonNClusITS, c.bachPionNClusITS, c.protonNClusTPC, c.pionNClusTPC, c.bachKaonNClusTPC, c.bachPionNClusTPC,
                 c.protonTPCNSigma, c.pionTPCNSigma, c.bachKaonTPCNSigma, c.bachPionTPCNSigma,
                 c.protonHasTOF, c.pionHasTOF, c.bachKaonHasTOF, c.bachPionHasTOF,
                 c.protonTOFNSigma, c.pionTOFNSigma, c.bachKaonTOFNSigma, c.bachPionTOFNSigma,
                 particle.pt(), particle.eta(), particle.phi(), particle.pdgCode(), mcCollision.posX() - particle.vx(), mcCollision.posY() - particle.vy(), mcCollision.posZ() - particle.vz(), mcCollision.globalIndex() == label.mcCollisionId());
    }

    for (const auto& p : mcParticles) {
      auto absCode = std::abs(p.pdgCode());
      if (absCode != 3312 && absCode != 3334) {
        continue;
      }
      auto fromHF = isFromHF(p);
      int pdgCodeMom = p.has_mothers() ? p.mothers_as<aod::McParticles>()[0].pdgCode() : 0;
      auto mcCollision = p.mcCollision_as<aod::McCollisions>();

      NPCTableGen(p.pt(), p.eta(), p.phi(), p.pdgCode(), pdgCodeMom, mcCollision.posX() - p.vx(), mcCollision.posY() - p.vy(), mcCollision.posZ() - p.vz(), fromHF[0], fromHF[1]);
    }
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processTrackedCascadesMC, "process cascades from strangeness tracking: MC analysis", true);

  void processTrackedCascadesData(CollisionCandidatesRun3 const& collisions,
                                  aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& /*cascades*/,
                                  aod::V0s const& /*v0s*/, TracksExtData const& /*tracks*/,
                                  aod::BCsWithTimestamps const&)
  {
    fillCandidatesVector<TracksExtData>(collisions, trackedCascades);
    for (const auto& c : candidates) {
      NPCTable(c.matchingChi2, c.deltaPt, c.itsClusSize, c.hasReassociatedCluster,
               c.pvContributors, c.pvTimeResolution, c.pvX, c.pvY, c.pvZ,
               c.cascPt, c.cascEta, c.cascPhi,
               c.protonPt, c.protonEta, c.pionPt, c.pionEta, c.bachPt, c.bachEta,
               c.cascDCAxy, c.cascDCAz, c.protonDCAxy, c.protonDCAz, c.pionDCAxy, c.pionDCAz, c.bachDCAxy, c.bachDCAz,
               c.casccosPA, c.v0cosPA,
               c.massXi, c.massOmega, c.massV0,
               c.cascRadius, c.v0radius, c.cascLength, c.v0length,
               c.cascNClusITS, c.protonNClusITS, c.pionNClusITS, c.bachKaonNClusITS, c.bachPionNClusITS, c.protonNClusTPC, c.pionNClusTPC, c.bachKaonNClusTPC, c.bachPionNClusTPC,
               c.protonTPCNSigma, c.pionTPCNSigma, c.bachKaonTPCNSigma, c.bachPionTPCNSigma,
               c.protonHasTOF, c.pionHasTOF, c.bachKaonHasTOF, c.bachPionHasTOF,
               c.protonTOFNSigma, c.pionTOFNSigma, c.bachKaonTOFNSigma, c.bachPionTOFNSigma);
    }
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processTrackedCascadesData, "process cascades from strangeness tracking: Data analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NonPromptCascadeTask>(cfgc)};
}
