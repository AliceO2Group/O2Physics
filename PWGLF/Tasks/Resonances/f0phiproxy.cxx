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

/// \file f0phiproxy.cxx
/// \brief Dipion-based half-momentum kaon proxies with kaon/phi and pion/K*0 bachelors, in data and MC.
/// \author Sushanta Tripathy

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <initializer_list>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct F0phiproxy {
  using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>;
  using Tracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                           aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi,
                           aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using CollisionsMC = soa::Join<Collisions, aod::McCollisionLabels>;
  using TracksMC = soa::Join<Tracks, aod::McTrackLabels>;
  using FourVector = ROOT::Math::PxPyPzMVector;
  Preslice<Tracks> perCollision = aod::track::collisionId;
  Preslice<TracksMC> perCollisionMC = aod::track::collisionId;
  Service<o2::framework::O2DatabasePDG> pdgDB{};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> vertexMax{"vertexMax", 10.f, "Maximum absolute vertex z (cm)"};
  Configurable<bool> requireSel8{"requireSel8", true, "Require sel8"};
  Configurable<bool> requireINELgt0{"requireINELgt0", true, "Require at least one PV track in |eta|<1; MC gen uses charged physical primaries"};
  Configurable<bool> requireNoSameBunchPileup{"requireNoSameBunchPileup", false, "Require kNoSameBunchPileup, as in Lambda-star trigger mode 4"};
  Configurable<bool> requireGoodZvtxFT0vsPV{"requireGoodZvtxFT0vsPV", false, "Require kIsGoodZvtxFT0vsPV, as in Lambda-star trigger mode 4"};
  Configurable<bool> useFT0M{"useFT0M", true, "Use calibrated FT0M percentiles for reconstructed data/MC; false retains legacy NPV axes and mixing"};
  ConfigurableAxis ft0mAxis{"ft0mAxis", {0., 0., 1., 5., 10., 15., 20., 30., 40., 50., 70., 100.}, "FT0M percentile classes; first zero selects variable-width bins; must cover [0,100]"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15f, "Minimum track pT (GeV/c)"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.8f, "Maximum track |eta|"};
  Configurable<bool> requireGlobalTrack{"requireGlobalTrack", true, "Require global track selection"};
  Configurable<bool> requirePVContributor{"requirePVContributor", true, "Require PV contributor"};
  Configurable<float> pionTPC{"pionTPC", 3.f, "Pion TPC nSigma cut"};
  Configurable<float> pionTOF{"pionTOF", 3.f, "Pion TOF veto when matched"};
  Configurable<float> kaonTPC{"kaonTPC", 3.f, "Kaon TPC nSigma cut"};
  Configurable<float> kaonTOF{"kaonTOF", 3.f, "Kaon TOF veto when matched"};
  Configurable<bool> requireTOF{"requireTOF", false, "Require TOF for both species instead of veto-only"};
  Configurable<float> dipionYMax{"dipionYMax", 0.5f, "Maximum dipion |y|"};
  Configurable<float> proxyYMax{"proxyYMax", 0.5f, "Maximum proxy+bachelor and standard-control |y|"};
  Configurable<float> dipionMassMin{"dipionMassMin", 0.6f, "Lower dipion mass for proxy construction"};
  Configurable<float> dipionMassMax{"dipionMassMax", 1.4f, "Upper dipion mass for proxy construction"};
  Configurable<int> mixDepth{"mixDepth", 5, "Previous events per mixing bin; zero disables ME booking and mixing"};
  Configurable<float> mixZWidth{"mixZWidth", 2.f, "Vertex bin width (cm)"};
  Configurable<float> mixMultWidth{"mixMultWidth", 20.f, "PV multiplicity bin width"};
  Configurable<float> mixFT0MWidth{"mixFT0MWidth", 1.f, "FT0M mixing width in percentage points; pools also respect analysis-class boundaries"};

  Configurable<bool> enablePhiProxy{"enablePhiProxy", true, "Dipion/2 as on-shell K plus bachelor K"};
  Configurable<bool> enableKstarProxy{"enableKstarProxy", false, "Dipion/2 as on-shell K plus bachelor pion; use the dedicated Kstar configuration"};
  Configurable<bool> enablePhiControls{"enablePhiControls", true, "Standard KK SE/ME OS/LS spectra"};
  Configurable<bool> enableKstarControls{"enableKstarControls", false, "Standard K pi SE/ME OS/LS spectra"};
  Configurable<bool> enableCompactProxy{"enableCompactProxy", true, "Proxy maps retaining mass, pT, dipion mass, grouped multiplicity, charges and TOF"};
  Configurable<bool> enableFullProxySparse{"enableFullProxySparse", false, "Diagnostic 7D maps including multiplicity, bachelor charge and TOF; expensive to merge"};
  Configurable<bool> enablePidQA{"enablePidQA", false, "Book/fill the four PID QA histograms"};
  Configurable<int> multiplicityBins{"multiplicityBins", 100, "Multiplicity bins in dipion/control and full diagnostic maps; identical settings required for merging"};
  Configurable<int> proxyDipionMassBins{"proxyDipionMassBins", 300, "Dipion mass bins in compact proxy maps; default preserves original 5 MeV resolution"};
  ConfigurableAxis proxyMultiplicityAxis{"proxyMultiplicityAxis", {0., 0., 20., 50., 200., 1000000.}, "Compact proxy NPV intervals; first zero selects variable-width bins"};
  ConfigurableAxis proxyPtAxis{"proxyPtAxis", {0., 0., 0.3, 1., 2., 3., 4., 6., 10.}, "Compact proxy pT intervals covering both the data and Kstar-study selections"};
  Configurable<int> kstarMassBins{"kstarMassBins", 200, "K pi and Kstar-proxy mass bins over 0.6-1.4 GeV; 200 gives the study's 4 MeV bins"};
  Configurable<bool> mcPrimaryOnly{"mcPrimaryOnly", true, "Generated daughters must be physical primaries; strong-decay daughters are included"};
  Configurable<int> mcF0PdgCode{"mcF0PdgCode", 0, "0 resolves f0_980 from the O2 PDG service; override only for a verified generator convention"};

  static constexpr size_t Phi = 0;
  static constexpr size_t Kstar = 1;
  static constexpr size_t NChannels = 2;
  static constexpr size_t NCharges = 3;
  static constexpr size_t NMixClasses = 2;
  static constexpr size_t NTruthClasses = 4;
  static constexpr double ProxyMomentumScale = 0.5;
  static constexpr int PdgPion = static_cast<int>(PDG_t::kPiPlus);
  static constexpr int PdgKaon = static_cast<int>(PDG_t::kKPlus);
  static constexpr int PdgPhi = static_cast<int>(constants::physics::Pdg::kPhi);
  static constexpr int PdgKstar = static_cast<int>(constants::physics::Pdg::kK0Star892);
  static constexpr size_t TwoBodyDaughters = 2;
  static constexpr int IonCodeThreshold = 1000000000;
  static constexpr int IonZDivisor = 10000;
  static constexpr int IonZModulo = 1000;
  static constexpr float FT0MPercentileMax = 100.f;
  static constexpr size_t MinAxisSpecificationSize = 3;
  static constexpr double MaxUniformFT0MBins = 10000.;
  int resolvedF0Code = 0;
  std::vector<double> ft0mEdges{};

  struct TruthInfo {
    int64_t id = -1;
    int collisionId = -1;
    int pdg = 0;
    double px = 0., py = 0., pz = 0.;
    bool primary = false;
    std::vector<std::pair<int64_t, int>> parents{};
  };
  struct TrackCandidate {
    int64_t id = -1;
    int sign = 0;
    double px = 0., py = 0., pz = 0.;
    bool hasTOF = false;
    TruthInfo truth{};
    [[nodiscard]] FourVector vector(double mass) const { return {px, py, pz, mass}; }
  };
  struct Event {
    std::vector<TrackCandidate> pions{}, kaons{};
  };

  using MixingKey = std::tuple<int, int, int, int>;
  enum EventStage : int { Received = 1,
                          VertexAccepted,
                          Sel8Accepted,
                          InelAccepted,
                          QualityAccepted,
                          ActivityAccepted,
                          NEventStages = ActivityAccepted };
  std::shared_ptr<TH1> eventSelection{}, eventStages{}, selectedEvents{}, ft0mStatusQA{}, mixPartners{};
  std::shared_ptr<TH2> eventStagesByActivity{}, activityVsNPV{}, activityVsVertex{}, mixPartnersByActivity{};
  std::array<std::shared_ptr<TH2>, 4> pidQA{};
  std::array<std::shared_ptr<TH3>, 4> pidQAByActivity{};

  std::array<std::shared_ptr<THnSparse>, NMixClasses> dipionSpectra{};
  std::array<std::array<std::shared_ptr<THnSparse>, NMixClasses>, NChannels> compactProxy{}, fullProxy{};
  std::array<std::array<std::array<std::shared_ptr<THnSparse>, NCharges>, NMixClasses>, NChannels> controls{};
  std::array<std::array<std::shared_ptr<TH2>, NTruthClasses>, NChannels> recoComponents{};
  std::array<std::shared_ptr<TH2>, NChannels> recoTrueParent{}, recoPtResponse{}, proxyPtResponse{};
  std::shared_ptr<TH2> recoTrueF0{};
  std::array<std::array<std::shared_ptr<THnSparse>, NTruthClasses>, NChannels> recoComponentsByActivity{};
  std::array<std::shared_ptr<THnSparse>, NChannels> recoTrueParentByActivity{}, recoPtResponseByActivity{}, proxyPtResponseByActivity{};
  std::shared_ptr<THnSparse> recoTrueF0ByActivity{};
  enum GenSpectrum : size_t { GenDipion,
                              GenF0,
                              GenPhi,
                              GenPhiTrue,
                              GenKstar,
                              GenKstarTrue,
                              GenProxyPhi,
                              GenProxyPhiF0,
                              GenProxyKstar,
                              GenProxyKstarF0,
                              NGenSpectra };
  std::array<std::shared_ptr<TH2>, NGenSpectra> generated{};

  bool mixingEnabled() const { return mixDepth.value > 0; }
  bool channelEnabled(size_t channel) const { return channel == Phi ? enablePhiProxy.value : enableKstarProxy.value; }
  bool controlEnabled(size_t channel) const { return channel == Phi ? enablePhiControls.value : enableKstarControls.value; }
  static double bachelorMass(size_t channel) { return channel == Phi ? constants::physics::MassKaonCharged : constants::physics::MassPionCharged; }
  static int bachelorPdg(size_t channel) { return channel == Phi ? PdgKaon : PdgPion; }

  std::string recoName(std::string const& name) const
  {
    return useFT0M.value ? std::string("FT0M/").append(name) : name;
  }

  static int ft0mStatus(float percentile)
  {
    if (!std::isfinite(percentile)) {
      return 1;
    }
    if (percentile < 0.f) {
      return 2;
    }
    return percentile < FT0MPercentileMax ? 0 : 3;
  }

  static std::vector<double> resolveFT0MEdges(std::vector<double> const& specification)
  {
    if (specification.size() < MinAxisSpecificationSize ||
        !std::all_of(specification.begin(), specification.end(), [](double x) { return std::isfinite(x); })) {
      throw std::runtime_error("Invalid FT0M axis specification");
    }
    std::vector<double> edges;
    if (specification.front() == 0.) {
      edges.assign(specification.begin() + 1, specification.end());
    } else {
      const double bins = specification.front();
      if (specification.size() != MinAxisSpecificationSize || bins < 1. || bins > MaxUniformFT0MBins || std::floor(bins) != bins) {
        throw std::runtime_error("Invalid uniform FT0M axis specification");
      }
      for (int i = 0; i <= static_cast<int>(bins); ++i) {
        edges.push_back(specification[1] + (specification[2] - specification[1]) * i / bins);
      }
    }
    if (edges.front() != 0. || edges.back() != FT0MPercentileMax ||
        std::adjacent_find(edges.begin(), edges.end(), [](double a, double b) { return a >= b; }) != edges.end()) {
      throw std::runtime_error("FT0M classes must have increasing edges covering [0,100]");
    }
    return edges;
  }

  AxisSpec activityAxis(bool compact = false) const
  {
    if (useFT0M.value) {
      return AxisSpec{ft0mEdges, "FT0M percentile (MB>0), calibrated"};
    }
    return compact ? AxisSpec{proxyMultiplicityAxis, "N_{PV}(|#eta|<1), grouped"} : AxisSpec{multiplicityBins.value, 0., 200., "N_{PV}(|#eta|<1)"};
  }

  int ft0mClass(float percentile) const
  {
    if (ft0mStatus(percentile) != 0) {
      return -1;
    }
    return static_cast<int>(std::upper_bound(ft0mEdges.begin(), ft0mEdges.end(), percentile) - ft0mEdges.begin()) - 1;
  }

  MixingKey mixingKey(int run, float vertex, float mult) const
  {
    const int activityClass = useFT0M.value ? ft0mClass(mult) : 0;
    if (useFT0M.value && activityClass < 0) {
      throw std::runtime_error("Invalid FT0M event cannot enter a mixing pool");
    }
    const float width = useFT0M.value ? mixFT0MWidth.value : mixMultWidth.value;
    auto bin = [](float value, float binWidth) {
      const double index = std::floor(value / binWidth);
      if (!std::isfinite(index) || index < std::numeric_limits<int>::min() || index > std::numeric_limits<int>::max()) {
        throw std::runtime_error("Mixing coordinate cannot be represented by the configured bin width");
      }
      return static_cast<int>(index);
    };
    return {run, bin(vertex, mixZWidth.value), activityClass, bin(mult, width)};
  }

  void init(InitContext const&)
  {
    if (doprocess && doprocessMCReco) {
      throw std::runtime_error("Enable either process (data) or processMCReco, never both");
    }
    if (mixDepth.value < 0 || !std::isfinite(mixZWidth.value) || mixZWidth.value <= 0.f ||
        !std::isfinite(mixMultWidth.value) || mixMultWidth.value <= 0.f ||
        !std::isfinite(mixFT0MWidth.value) || mixFT0MWidth.value <= 0.f ||
        vertexMax.value <= 0.f || trackPtMin.value < 0.f || trackEtaMax.value <= 0.f ||
        pionTPC.value <= 0.f || pionTOF.value <= 0.f || kaonTPC.value <= 0.f || kaonTOF.value <= 0.f ||
        dipionMassMin.value >= dipionMassMax.value || dipionYMax.value <= 0.f || proxyYMax.value <= 0.f ||
        multiplicityBins.value <= 0 || proxyDipionMassBins.value <= 0 || kstarMassBins.value <= 0 || mcF0PdgCode.value < 0) {
      throw std::runtime_error("Invalid F0PhiProxy configuration");
    }
    if (useFT0M.value) {
      ft0mEdges = resolveFT0MEdges(ft0mAxis.value);
    }
    if (doprocess || doprocessMCReco) {
      if ((enablePhiProxy.value || enableKstarProxy.value) && !enableCompactProxy.value && !enableFullProxySparse.value) {
        throw std::runtime_error("Enable compact or full proxy output when a proxy channel is enabled");
      }
      bookReco();
    }
    if (doprocessMCReco || doprocessMCGen) {
      resolvedF0Code = mcF0PdgCode.value;
      if (resolvedF0Code == 0) {
        auto const* particle = pdgDB->GetParticle("f0_980");
        if (!particle) {
          throw std::runtime_error("f0_980 not found in O2 PDG service; configure mcF0PdgCode explicitly");
        }
        resolvedF0Code = particle->PdgCode();
      }
      bookMC();
    }
  }

  void bookReco()
  {
    const AxisSpec massPi{300, 0.2, 1.7, "m_{#pi#pi} (GeV/c^{2})"};
    const AxisSpec compactMassPi{proxyDipionMassBins.value, 0.2, 1.7, "m_{#pi#pi} (GeV/c^{2})"};
    const AxisSpec pt{100, 0., 10., "p_{T} (GeV/c)"};
    const AxisSpec mult = activityAxis();
    const AxisSpec compactMult = activityAxis(true);
    const AxisSpec compactPt{proxyPtAxis, "proxy pair p_{T} (GeV/c), grouped"};
    const AxisSpec charge{3, -0.5, 2.5, "dipion: 0=unlike, 1=++, 2=--"};
    const AxisSpec bachelorCharge{2, -1.5, 1.5, "bachelor charge"};
    const AxisSpec tof{2, -0.5, 1.5, "bachelor has TOF"};
    eventSelection = histos.add<TH1>(recoName("hEventSelection").c_str(), "Received; selected including activity validity", HistType::kTH1D, {{2, 0.5, 2.5}});
    eventSelection->GetXaxis()->SetBinLabel(1, "Received");
    eventSelection->GetXaxis()->SetBinLabel(2, "Selected");
    const AxisSpec stages{NEventStages, 0.5, static_cast<double>(NEventStages) + 0.5, "Event selection stage"};
    eventStages = histos.add<TH1>(recoName("hEventSelectionSteps").c_str(), "Reconstructed event cut flow", HistType::kTH1D, {stages});
    eventStagesByActivity = histos.add<TH2>(recoName("hEventSelectionVsMultiplicity").c_str(), "Valid-activity collisions;activity;selection stage", HistType::kTH2D, {compactMult, stages});
    const std::array<const char*, NEventStages> stageNames{"Received", "Vertex", "sel8 or disabled", "INEL>0 or disabled", "Quality flags or disabled", "Valid activity selected"};
    for (int i = 0; i < NEventStages; ++i) {
      eventStages->GetXaxis()->SetBinLabel(i + 1, stageNames[i]);
      eventStagesByActivity->GetYaxis()->SetBinLabel(i + 1, stageNames[i]);
    }
    selectedEvents = histos.add<TH1>(recoName("hSelectedEvents").c_str(), "Selected collisions counted once, including zero-candidate events;activity;Events", HistType::kTH1D, {compactMult});
    activityVsNPV = histos.add<TH2>(recoName("hMultiplicityVsNPV").c_str(), "Selected collisions;activity;N_{PV}(|#eta|<1), overflow retained", HistType::kTH2D, {compactMult, {400, 0., 400.}});
    activityVsVertex = histos.add<TH2>(recoName("hMultiplicityVsVertex").c_str(), "Selected collisions;activity;vertex z (cm)", HistType::kTH2D, {compactMult, {100, -vertexMax.value, vertexMax.value}});
    auto estimator = histos.add<TH1>("hMultiplicityEstimator", "Configured reconstructed estimator;mode;Workers", HistType::kTH1D, {{2, -0.5, 1.5}});
    estimator->GetXaxis()->SetBinLabel(1, "NPV");
    estimator->GetXaxis()->SetBinLabel(2, "FT0M");
    estimator->Fill(useFT0M.value ? 1. : 0.);
    if (useFT0M.value) {
      ft0mStatusQA = histos.add<TH1>("FT0M/hPercentileStatus", "After base event selection;FT0M status;Events", HistType::kTH1D, {{4, -0.5, 3.5}});
      const std::array<const char*, 4> statusNames{"Valid [0,100)", "Nonfinite", "Below zero", "At or above 100"};
      for (size_t i = 0; i < statusNames.size(); ++i) {
        ft0mStatusQA->GetXaxis()->SetBinLabel(static_cast<int>(i + 1), statusNames[i]);
      }
    }
    if (mixingEnabled()) {
      mixPartnersByActivity = histos.add<TH2>(recoName("hMixPartnersVsMultiplicity").c_str(), "Available previous events;activity;mixing partners", HistType::kTH2D, {compactMult, {101, -0.5, 100.5}});
      mixPartners = histos.add<TH1>(recoName("hMixPartners").c_str(), "Mixing partners;N", HistType::kTH1D, {{101, -0.5, 100.5}});
    }
    if (enablePidQA.value) {
      const std::array<const char*, 4> names{"hPionTPC", "hPionTOF", "hKaonTPC", "hKaonTOF"};
      for (size_t i = 0; i < names.size(); ++i) {
        pidQA[i] = histos.add<TH2>(recoName(names[i]).c_str(), "Quality-selected tracks, before PID;pT;nSigma", HistType::kTH2F, {pt, {100, -10., 10.}});
        if (useFT0M.value) {
          const std::string name = std::string(names[i]).append("ByMultiplicity");
          pidQAByActivity[i] = histos.add<TH3>(recoName(name).c_str(), "Quality-selected tracks, before PID;pT;nSigma;FT0M percentile", HistType::kTH3F, {pt, {100, -10., 10.}, compactMult});
        }
      }
    }
    const std::array<const char*, NMixClasses> dipionNames{"hDipionSE", "hDipionME"};
    for (size_t mix = 0; mix < dipionNames.size(); ++mix) {
      if (mix != 0 && !mixingEnabled()) {
        continue;
      }
      dipionSpectra[mix] = histos.add<THnSparse>(recoName(dipionNames[mix]).c_str(), "Dipion SE or raw ME", HistType::kTHnSparseF, {massPi, pt, mult, charge});
    }
    const std::array<std::array<const char*, NMixClasses>, NChannels> compactNames{{{"hProxyCompactSE", "hProxyCompactME"}, {"hKstarProxyCompactSE", "hKstarProxyCompactME"}}};
    const std::array<std::array<const char*, NMixClasses>, NChannels> fullNames{{{"hProxySE", "hProxyME"}, {"hKstarProxySE", "hKstarProxyME"}}};
    const std::array<std::array<std::array<const char*, NCharges>, NMixClasses>, NChannels> controlNames{{{{{"hPhiKK", "hPhiKKLikePP", "hPhiKKLikeMM"}, {"hPhiKKMixed", "hPhiKKMixedLikePP", "hPhiKKMixedLikeMM"}}},
                                                                                                          {{{"hKstarKPi", "hKstarKPiLikePP", "hKstarKPiLikeMM"}, {"hKstarKPiMixed", "hKstarKPiMixedLikePP", "hKstarKPiMixedLikeMM"}}}}};
    for (size_t channel = 0; channel < NChannels; ++channel) {
      const AxisSpec mass = channel == Phi ? AxisSpec{600, 0.9, 1.5, "m_{K,proxy K} (GeV/c^{2})"} : AxisSpec{kstarMassBins.value, 0.6, 1.4, "m_{K,proxy #pi} (GeV/c^{2})"};
      for (size_t mix = 0; mix < NMixClasses; ++mix) {
        if (mix != 0 && !mixingEnabled()) {
          continue;
        }
        if (channelEnabled(channel)) {
          // The measured dipion mass is retained; a mass window is not an f0 tag.
          if (enableCompactProxy.value) {
            compactProxy[channel][mix] = histos.add<THnSparse>(recoName(compactNames[channel][mix]).c_str(), "Grouped-pT/multiplicity proxy; no constituent charge is measured", HistType::kTHnSparseF, {mass, compactPt, compactMassPi, compactMult, charge, bachelorCharge, tof});
          }
          if (enableFullProxySparse.value) {
            fullProxy[channel][mix] = histos.add<THnSparse>(recoName(fullNames[channel][mix]).c_str(), "Diagnostic full proxy", HistType::kTHnSparseF, {mass, pt, massPi, mult, charge, bachelorCharge, tof});
          }
        }
        if (controlEnabled(channel)) {
          for (size_t sign = 0; sign < NCharges; ++sign) {
            controls[channel][mix][sign] = histos.add<THnSparse>(recoName(controlNames[channel][mix][sign]).c_str(), "Standard control SE or raw ME", HistType::kTHnSparseF, {mass, pt, mult});
          }
        }
      }
    }
  }

  void bookMC()
  {
    const AxisSpec pt{100, 0., 10., "pair pT (GeV/c)"};
    const AxisSpec pionMass{300, 0.2, 1.7, "dipion mass (GeV/c2)"};
    const std::array<AxisSpec, NChannels> mass{AxisSpec{600, 0.9, 1.5, "KK mass (GeV/c2)"}, AxisSpec{kstarMassBins.value, 0.6, 1.4, "K pi mass (GeV/c2)"}};
    if (doprocessMCGen) {
      histos.add("hMCGenEventSelection", "Generated events;1=all,2=vertex,3=truth INEL>0", HistType::kTH1D, {{3, 0.5, 3.5}});
      const std::array<const char*, NGenSpectra> names{"hMCGenDipion", "hMCGenTrueF0", "hMCGenPhiKK", "hMCGenTruePhi", "hMCGenKstarKPi", "hMCGenTrueKstar", "hMCGenProxyPhi", "hMCGenProxyPhiTrueF0", "hMCGenProxyKstar", "hMCGenProxyKstarTrueF0"};
      for (size_t i = 0; i < names.size(); ++i) {
        const AxisSpec axis = i <= GenF0 ? pionMass : (i == GenKstar || i == GenKstarTrue || i >= GenProxyKstar ? mass[Kstar] : mass[Phi]);
        generated[i] = histos.add<TH2>(names[i], "Generated unweighted fiducial spectrum; no detector FT0M selection", HistType::kTH2D, {axis, pt}, true);
      }
    }
    if (doprocessMCReco) {
      recoTrueF0 = histos.add<TH2>(recoName("hMCRecoTrueF0").c_str(), "Selected common-parent f0 -> pi pi", HistType::kTH2D, {pionMass, pt}, true);
      if (useFT0M.value) {
        recoTrueF0ByActivity = histos.add<THnSparse>("FT0M/hMCRecoTrueF0ByMultiplicity", "Selected common-parent f0;dipion mass;pair pT;FT0M percentile", HistType::kTHnSparseD, {pionMass, pt, activityAxis(true)}, true);
      }
      const std::array<const char*, NChannels> channelNames{"Phi", "Kstar"};
      const std::array<const char*, NTruthClasses> classNames{"TrueF0", "Combinatorial", "MisID", "Unmatched"};
      for (size_t channel = 0; channel < NChannels; ++channel) {
        if (channelEnabled(channel)) {
          for (size_t category = 0; category < classNames.size(); ++category) {
            const std::string name = std::string("hMCRecoProxy").append(channelNames[channel]).append(classNames[category]);
            recoComponents[channel][category] = histos.add<TH2>(recoName(name).c_str(), "Exclusive proxy truth component; not a coalescence ancestry tag", HistType::kTH2D, {mass[channel], pt}, true);
            if (useFT0M.value) {
              recoComponentsByActivity[channel][category] = histos.add<THnSparse>(recoName(std::string(name).append("ByMultiplicity")).c_str(), "Exclusive proxy truth component; not coalescence ancestry", HistType::kTHnSparseD, {mass[channel], pt, activityAxis(true)}, true);
            }
          }
          const std::string name = std::string("hMCRecoProxy").append(channelNames[channel]).append("PtResponse");
          proxyPtResponse[channel] = histos.add<TH2>(recoName(name).c_str(), "Primary true f0 daughters and true bachelor;generated proxy pair pT;reconstructed proxy pair pT", HistType::kTH2D, {pt, pt}, true);
          if (useFT0M.value) {
            proxyPtResponseByActivity[channel] = histos.add<THnSparse>(recoName(std::string(name).append("ByMultiplicity")).c_str(), "Selected true-f0 proxy response;generated pair pT;reconstructed pair pT;FT0M percentile", HistType::kTHnSparseD, {pt, pt, activityAxis(true)}, true);
          }
        }
        if (controlEnabled(channel)) {
          const std::string truthName = std::string("hMCRecoTrue").append(channelNames[channel]);
          const std::string responseName = std::string("hMCReco").append(channelNames[channel]).append("PtResponse");
          recoTrueParent[channel] = histos.add<TH2>(recoName(truthName).c_str(), "Selected exact common-parent standard control", HistType::kTH2D, {mass[channel], pt}, true);
          recoPtResponse[channel] = histos.add<TH2>(recoName(responseName).c_str(), "Selected standard control;generated pair pT;reconstructed pair pT", HistType::kTH2D, {pt, pt}, true);
          if (useFT0M.value) {
            recoTrueParentByActivity[channel] = histos.add<THnSparse>(recoName(std::string(truthName).append("ByMultiplicity")).c_str(), "Selected exact common-parent standard control", HistType::kTHnSparseD, {mass[channel], pt, activityAxis(true)}, true);
            recoPtResponseByActivity[channel] = histos.add<THnSparse>(recoName(std::string(responseName).append("ByMultiplicity")).c_str(), "Selected standard response;generated pair pT;reconstructed pair pT;FT0M percentile", HistType::kTHnSparseD, {pt, pt, activityAxis(true)}, true);
          }
        }
      }
    }
  }

  static void fillSparse(std::shared_ptr<THnSparse> const& histogram, std::initializer_list<double> coordinates)
  {
    if (histogram) {
      histogram->Fill(coordinates.begin());
    }
  }

  bool passPID(float tpc, float tof, bool matched, float tpcCut, float tofCut) const
  {
    return std::isfinite(tpc) && std::abs(tpc) < tpcCut && (!requireTOF.value || matched) &&
           (!matched || (std::isfinite(tof) && std::abs(tof) < tofCut));
  }
  static int chargeCategory(int a, int b) { return a != b ? 0 : (a > 0 ? 1 : 2); }
  static FourVector constituent(FourVector const& dipion)
  {
    return {ProxyMomentumScale * dipion.Px(), ProxyMomentumScale * dipion.Py(), ProxyMomentumScale * dipion.Pz(), constants::physics::MassKaonCharged};
  }
  int chargeFromPdg(int code) const
  {
    const int absolute = std::abs(code);
    if (absolute >= IonCodeThreshold) {
      const int z = (absolute / IonZDivisor) % IonZModulo;
      return code > 0 ? z : -z;
    }
    if (auto const* particle = pdgDB->GetParticle(code)) {
      return static_cast<int>(std::lround(particle->Charge() / 3.));
    }
    return 0;
  }
  template <typename TParticle>
  TruthInfo particleTruth(TParticle const& particle) const
  {
    TruthInfo truth{.id = particle.globalIndex(), .collisionId = particle.mcCollisionId(), .pdg = particle.pdgCode(), .px = particle.px(), .py = particle.py(), .pz = particle.pz(), .primary = particle.isPhysicalPrimary(), .parents = {}};
    for (auto const& mother : particle.template mothers_as<aod::McParticles>()) {
      const int code = std::abs(mother.pdgCode());
      if (mother.mcCollisionId() == truth.collisionId && (code == resolvedF0Code || code == PdgPhi || code == PdgKstar)) {
        // The complementary-channel study uses exact two-body parent tags.
        if (mother.template daughters_as<aod::McParticles>().size() == TwoBodyDaughters) {
          truth.parents.emplace_back(mother.globalIndex(), code);
        }
      }
    }
    return truth;
  }
  template <typename TTrack>
  TruthInfo trackTruth(TTrack const& track, int collisionId) const
  {
    if (collisionId < 0 || !track.has_mcParticle()) {
      return {};
    }
    auto particle = track.template mcParticle_as<aod::McParticles>();
    return particle.mcCollisionId() == collisionId ? particleTruth(particle) : TruthInfo{};
  }
  static bool commonParent(TruthInfo const& first, TruthInfo const& second, int parentCode)
  {
    if (first.id < 0 || second.id < 0 || first.id == second.id || first.collisionId != second.collisionId) {
      return false;
    }
    return std::any_of(first.parents.begin(), first.parents.end(), [&second, parentCode](auto const& parent) {
      return parent.second == parentCode && std::find(second.parents.begin(), second.parents.end(), parent) != second.parents.end();
    });
  }
  static bool identified(TrackCandidate const& candidate, int positiveCode)
  {
    return candidate.truth.id >= 0 && candidate.truth.pdg == candidate.sign * positiveCode;
  }
  bool trueF0(TrackCandidate const& first, TrackCandidate const& second) const
  {
    return first.sign != second.sign && identified(first, PdgPion) && identified(second, PdgPion) && commonParent(first.truth, second.truth, resolvedF0Code);
  }
  static double truthPairPt(TruthInfo const& first, TruthInfo const& second)
  {
    return std::hypot(first.px + second.px, first.py + second.py);
  }
  size_t proxyTruthClass(TrackCandidate const& first, TrackCandidate const& second, TrackCandidate const& bachelor, size_t channel) const
  {
    if (first.truth.id < 0 || second.truth.id < 0 || bachelor.truth.id < 0 ||
        first.truth.collisionId != second.truth.collisionId || first.truth.collisionId != bachelor.truth.collisionId) {
      return NTruthClasses - 1;
    }
    if (!identified(first, PdgPion) || !identified(second, PdgPion) || !identified(bachelor, bachelorPdg(channel))) {
      return 2;
    }
    const bool distinctBachelor = bachelor.truth.id != first.truth.id && bachelor.truth.id != second.truth.id;
    return distinctBachelor && trueF0(first, second) ? 0 : 1;
  }

  void fillProxy(FourVector const& dipion, TrackCandidate const& first, TrackCandidate const& second,
                 std::vector<TrackCandidate> const& bachelors, float mult, bool mixed, size_t channel, bool mc = false)
  {
    if (!channelEnabled(channel)) {
      return;
    }
    const FourVector proxy = constituent(dipion);
    const int category = chargeCategory(first.sign, second.sign);
    for (auto const& bachelor : bachelors) {
      if (!mixed && (bachelor.id == first.id || bachelor.id == second.id)) {
        continue;
      }
      const auto pair = proxy + bachelor.vector(bachelorMass(channel));
      if (std::abs(pair.Rapidity()) >= proxyYMax.value) {
        continue;
      }
      // Neutral dipion: one entry per bachelor. There is no measured proxy-kaon charge.
      fillSparse(compactProxy[channel][mixed], {pair.M(), pair.Pt(), dipion.M(), mult, static_cast<double>(category), static_cast<double>(bachelor.sign), static_cast<double>(bachelor.hasTOF)});
      fillSparse(fullProxy[channel][mixed], {pair.M(), pair.Pt(), dipion.M(), mult, static_cast<double>(category), static_cast<double>(bachelor.sign), static_cast<double>(bachelor.hasTOF)});
      if (mc && !mixed) {
        const size_t truthClass = proxyTruthClass(first, second, bachelor, channel);
        recoComponents[channel][truthClass]->Fill(pair.M(), pair.Pt());
        fillSparse(recoComponentsByActivity[channel][truthClass], {pair.M(), pair.Pt(), mult});
        if (truthClass == 0 && first.truth.primary && second.truth.primary && bachelor.truth.primary &&
            bachelor.truth.id != first.truth.id && bachelor.truth.id != second.truth.id) {
          const double px = ProxyMomentumScale * (first.truth.px + second.truth.px) + bachelor.truth.px;
          const double py = ProxyMomentumScale * (first.truth.py + second.truth.py) + bachelor.truth.py;
          proxyPtResponse[channel]->Fill(std::hypot(px, py), pair.Pt());
          fillSparse(proxyPtResponseByActivity[channel], {std::hypot(px, py), pair.Pt(), mult});
        }
      }
    }
  }

  void fillControlPair(TrackCandidate const& first, TrackCandidate const& second, float mult, bool mixed, size_t channel, bool mc = false)
  {
    if ((!mixed && first.id == second.id) || !controlEnabled(channel)) {
      return;
    }
    const auto pair = first.vector(constants::physics::MassKaonCharged) + second.vector(bachelorMass(channel));
    if (std::abs(pair.Rapidity()) >= proxyYMax.value) {
      return;
    }
    fillSparse(controls[channel][mixed][chargeCategory(first.sign, second.sign)], {pair.M(), pair.Pt(), mult});
    const int parent = channel == Phi ? PdgPhi : PdgKstar;
    if (mc && !mixed && first.sign != second.sign && identified(first, PdgKaon) && identified(second, bachelorPdg(channel)) && commonParent(first.truth, second.truth, parent)) {
      recoTrueParent[channel]->Fill(pair.M(), pair.Pt());
      recoPtResponse[channel]->Fill(truthPairPt(first.truth, second.truth), pair.Pt());
      fillSparse(recoTrueParentByActivity[channel], {pair.M(), pair.Pt(), mult});
      fillSparse(recoPtResponseByActivity[channel], {truthPairPt(first.truth, second.truth), pair.Pt(), mult});
    }
  }
  void fillControls(Event const& event, std::deque<Event> const& pool, float mult, bool mc)
  {
    for (size_t i = 0; i < event.kaons.size(); ++i) {
      auto const& kaon = event.kaons[i];
      if (enablePhiControls.value) {
        for (size_t j = i + 1; j < event.kaons.size(); ++j) {
          fillControlPair(kaon, event.kaons[j], mult, false, Phi, mc);
        }
      }
      if (enableKstarControls.value) {
        for (auto const& pion : event.pions) {
          fillControlPair(kaon, pion, mult, false, Kstar, mc);
        }
      }
      for (auto const& previous : pool) {
        if (enablePhiControls.value) {
          for (auto const& other : previous.kaons) {
            fillControlPair(kaon, other, mult, true, Phi);
          }
        }
        if (enableKstarControls.value) {
          for (auto const& pion : previous.pions) {
            fillControlPair(kaon, pion, mult, true, Kstar);
          }
        }
      }
    }
    if (enableKstarControls.value) {
      // Both K(current)-pi(previous) and K(previous)-pi(current) are included.
      for (auto const& previous : pool) {
        for (auto const& kaon : previous.kaons) {
          for (auto const& pion : event.pions) {
            fillControlPair(kaon, pion, mult, true, Kstar);
          }
        }
      }
    }
  }

  template <bool IsMC, typename TCollisions, typename TTracks, typename TSlice>
  void processReco(TCollisions const& collisions, TTracks const& tracks, TSlice const& slice)
  {
    // Pools reset at each dataframe; the key also separates runs, z and multiplicity.
    std::map<MixingKey, std::deque<Event>> pools;
    const std::deque<Event> emptyPool;
    for (auto const& collision : collisions) {
      const float npv = collision.multNTracksPVeta1();
      const float mult = useFT0M.value ? collision.centFT0M() : npv;
      const bool validActivity = useFT0M.value ? ft0mStatus(mult) == 0 : std::isfinite(mult) && mult >= 0.f;
      auto stage = [&](int step) {
        eventStages->Fill(step);
        if (validActivity) {
          eventStagesByActivity->Fill(mult, step);
        }
      };
      eventSelection->Fill(1.);
      stage(Received);
      if (!std::isfinite(collision.posZ()) || std::abs(collision.posZ()) >= vertexMax.value) {
        continue;
      }
      stage(VertexAccepted);
      if (requireSel8.value && !collision.sel8()) {
        continue;
      }
      stage(Sel8Accepted);
      if (requireINELgt0.value && npv < 1.f) {
        continue;
      }
      stage(InelAccepted);
      if ((requireNoSameBunchPileup.value && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) ||
          (requireGoodZvtxFT0vsPV.value && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
        continue;
      }
      stage(QualityAccepted);
      if (useFT0M.value) {
        // Missing calibration (105) and non-INEL markers (100.5) are not percentiles.
        // Never substitute NPV in an FT0M histogram or pool.
        ft0mStatusQA->Fill(ft0mStatus(mult));
      }
      if (!validActivity) {
        continue;
      }
      stage(ActivityAccepted);
      eventSelection->Fill(2.);
      selectedEvents->Fill(mult);
      activityVsNPV->Fill(mult, npv);
      activityVsVertex->Fill(mult, collision.posZ());
      std::deque<Event>* activePool = nullptr;
      if (mixingEnabled()) {
        const int run = collision.template bc_as<aod::BCsWithTimestamps>().runNumber();
        activePool = &pools[mixingKey(run, collision.posZ(), mult)];
        mixPartners->Fill(activePool->size());
        mixPartnersByActivity->Fill(mult, activePool->size());
      }
      auto const& pool = activePool ? *activePool : emptyPool;
      Event event;
      auto selected = tracks.sliceBy(slice, collision.globalIndex());
      for (auto const& track : selected) {
        if (track.sign() == 0 || track.pt() < trackPtMin.value || std::abs(track.eta()) >= trackEtaMax.value ||
            (requireGlobalTrack.value && !track.isGlobalTrack()) || (requirePVContributor.value && !track.isPVContributor())) {
          continue;
        }
        const bool tof = track.hasTOF();
        if (enablePidQA.value) {
          auto fillPID = [&](size_t index, float value) {
            pidQA[index]->Fill(track.pt(), value);
            if (pidQAByActivity[index]) {
              pidQAByActivity[index]->Fill(track.pt(), value, mult);
            }
          };
          fillPID(0, track.tpcNSigmaPi());
          fillPID(2, track.tpcNSigmaKa());
          if (tof) {
            fillPID(1, track.tofNSigmaPi());
            fillPID(3, track.tofNSigmaKa());
          }
        }
        TrackCandidate candidate{.id = track.globalIndex(), .sign = track.sign(), .px = track.px(), .py = track.py(), .pz = track.pz(), .hasTOF = tof, .truth = {}};
        if constexpr (IsMC) {
          candidate.truth = trackTruth(track, collision.mcCollisionId());
        }
        if (passPID(track.tpcNSigmaPi(), track.tofNSigmaPi(), tof, pionTPC.value, pionTOF.value)) {
          event.pions.push_back(candidate);
        }
        if (passPID(track.tpcNSigmaKa(), track.tofNSigmaKa(), tof, kaonTPC.value, kaonTOF.value)) {
          event.kaons.push_back(candidate);
        }
      }
      fillDipions(event, pool, mult, IsMC);
      fillControls(event, pool, mult, IsMC);
      if (activePool) {
        // Mixed events never use ancestry; release that storage before caching.
        if constexpr (IsMC) {
          for (auto const& collection : {&event.pions, &event.kaons}) {
            for (size_t i = 0; i < collection->size(); ++i) {
              (*collection)[i].truth = {};
            }
          }
        }
        activePool->push_back(std::move(event));
        if (activePool->size() > static_cast<size_t>(mixDepth.value)) {
          activePool->pop_front();
        }
      }
    }
  }

  void fillDipions(Event const& event, std::deque<Event> const& pool, float mult, bool mc)
  {
    for (size_t i = 0; i < event.pions.size(); ++i) {
      auto const& first = event.pions[i];
      for (size_t j = i + 1; j < event.pions.size(); ++j) {
        auto const& second = event.pions[j];
        const auto pair = first.vector(constants::physics::MassPionCharged) + second.vector(constants::physics::MassPionCharged);
        if (std::abs(pair.Rapidity()) >= dipionYMax.value) {
          continue;
        }
        fillSparse(dipionSpectra[0], {pair.M(), pair.Pt(), mult, static_cast<double>(chargeCategory(first.sign, second.sign))});
        if (mc && trueF0(first, second)) {
          recoTrueF0->Fill(pair.M(), pair.Pt());
          fillSparse(recoTrueF0ByActivity, {pair.M(), pair.Pt(), mult});
        }
        if (pair.M() < dipionMassMin.value || pair.M() >= dipionMassMax.value) {
          continue;
        }
        fillProxy(pair, first, second, event.kaons, mult, false, Phi, mc);
        fillProxy(pair, first, second, event.pions, mult, false, Kstar, mc);
        for (auto const& previous : pool) {
          fillProxy(pair, first, second, previous.kaons, mult, true, Phi);
          fillProxy(pair, first, second, previous.pions, mult, true, Kstar);
        }
      }
      for (auto const& previous : pool) {
        for (auto const& second : previous.pions) {
          const auto pair = first.vector(constants::physics::MassPionCharged) + second.vector(constants::physics::MassPionCharged);
          if (std::abs(pair.Rapidity()) < dipionYMax.value) {
            fillSparse(dipionSpectra[1], {pair.M(), pair.Pt(), mult, static_cast<double>(chargeCategory(first.sign, second.sign))});
          }
        }
      }
    }
  }

  void fillGenerated(Event const& event)
  {
    for (size_t i = 0; i < event.pions.size(); ++i) {
      auto const& first = event.pions[i];
      for (size_t j = i + 1; j < event.pions.size(); ++j) {
        auto const& second = event.pions[j];
        if (first.sign == second.sign) {
          continue;
        }
        const auto pair = first.vector(constants::physics::MassPionCharged) + second.vector(constants::physics::MassPionCharged);
        if (std::abs(pair.Rapidity()) >= dipionYMax.value) {
          continue;
        }
        const bool tagged = trueF0(first, second);
        generated[GenDipion]->Fill(pair.M(), pair.Pt());
        if (tagged) {
          generated[GenF0]->Fill(pair.M(), pair.Pt());
        }
        if (pair.M() < dipionMassMin.value || pair.M() >= dipionMassMax.value) {
          continue;
        }
        const FourVector proxy = constituent(pair);
        for (size_t channel = 0; channel < NChannels; ++channel) {
          if (!channelEnabled(channel)) {
            continue;
          }
          auto const& bachelors = channel == Phi ? event.kaons : event.pions;
          for (auto const& bachelor : bachelors) {
            if (bachelor.id == first.id || bachelor.id == second.id) {
              continue;
            }
            const auto candidate = proxy + bachelor.vector(bachelorMass(channel));
            if (std::abs(candidate.Rapidity()) >= proxyYMax.value) {
              continue;
            }
            generated[channel == Phi ? GenProxyPhi : GenProxyKstar]->Fill(candidate.M(), candidate.Pt());
            if (tagged) {
              generated[channel == Phi ? GenProxyPhiF0 : GenProxyKstarF0]->Fill(candidate.M(), candidate.Pt());
            }
          }
        }
      }
    }
    for (size_t channel = 0; channel < NChannels; ++channel) {
      if (!controlEnabled(channel)) {
        continue;
      }
      for (size_t i = 0; i < event.kaons.size(); ++i) {
        auto const& first = event.kaons[i];
        auto const& others = channel == Phi ? event.kaons : event.pions;
        const size_t start = channel == Phi ? i + 1 : 0;
        for (size_t j = start; j < others.size(); ++j) {
          auto const& second = others[j];
          if (first.id == second.id || first.sign == second.sign) {
            continue;
          }
          const auto pair = first.vector(constants::physics::MassKaonCharged) + second.vector(bachelorMass(channel));
          if (std::abs(pair.Rapidity()) >= proxyYMax.value) {
            continue;
          }
          generated[channel == Phi ? GenPhi : GenKstar]->Fill(pair.M(), pair.Pt());
          if (commonParent(first.truth, second.truth, channel == Phi ? PdgPhi : PdgKstar)) {
            generated[channel == Phi ? GenPhiTrue : GenKstarTrue]->Fill(pair.M(), pair.Pt());
          }
        }
      }
    }
  }

  void processData(Collisions const& collisions, Tracks const& tracks, aod::BCsWithTimestamps const&)
  {
    processReco<false>(collisions, tracks, perCollision);
  }
  // Keep the original JSON option "process"; no unconditional process() method.
  PROCESS_SWITCH_FULL(F0phiproxy, processData, process, "Dipion phi/Kstar proxies in data", true);

  void processMCReco(CollisionsMC const& collisions, TracksMC const& tracks, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    processReco<true>(collisions, tracks, perCollisionMC);
  }
  PROCESS_SWITCH(F0phiproxy, processMCReco, "Data-identical reconstruction with truth components and responses", false);

  void processMCGen(aod::McCollision const& collision, aod::McParticles const& particles)
  {
    histos.fill(HIST("hMCGenEventSelection"), 1.);
    if (std::abs(collision.posZ()) >= vertexMax.value) {
      return;
    }
    histos.fill(HIST("hMCGenEventSelection"), 2.);
    bool inel = false;
    for (auto const& particle : particles) {
      // O2 ASoA row iterators and their end sentinel are not STL input iterators.
      // cppcheck-suppress useStlAlgorithm
      if (particle.isPhysicalPrimary() && std::abs(particle.eta()) < 1.f && chargeFromPdg(particle.pdgCode()) != 0) {
        inel = true;
        break;
      }
    }
    if (requireINELgt0.value && !inel) {
      return;
    }
    histos.fill(HIST("hMCGenEventSelection"), 3.);
    Event event;
    for (auto const& particle : particles) {
      const int absolute = std::abs(particle.pdgCode());
      if ((absolute != PdgPion && absolute != PdgKaon) || particle.pt() < trackPtMin.value ||
          std::abs(particle.eta()) >= trackEtaMax.value || (mcPrimaryOnly.value && !particle.isPhysicalPrimary())) {
        continue;
      }
      TrackCandidate candidate{.id = particle.globalIndex(), .sign = particle.pdgCode() > 0 ? 1 : -1, .px = particle.px(), .py = particle.py(), .pz = particle.pz(), .hasTOF = false, .truth = particleTruth(particle)};
      (absolute == PdgPion ? event.pions : event.kaons).push_back(std::move(candidate));
    }
    // Unweighted primary-daughter fiducial spectra, with no detector selection.
    // Native f0 truth is not a resonance-fed coalescence ancestry tag.
    fillGenerated(event);
  }
  PROCESS_SWITCH(F0phiproxy, processMCGen, "Independent generated fiducial spectra and exact common-parent tags", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<F0phiproxy>(context)};
}
