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

/// \file uccZdc.cxx
///
/// \brief task for analysis of UCC with the ZDC
/// \author Omar Vazquez (omar.vazquez.rueda@cern.ch)
/// \since January 29, 2025

#include <cmath>
#include <cstdint>
#include <cstdlib>

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/ZDCConstants.h"
#include "Framework/ASoAHelpers.h" // required for Filter op.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"
#include "TPDGCode.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::evsel;
using namespace o2::constants::physics;
using namespace o2::constants::math;

namespace o2::aod
{
using ColEvSels =
  soa::Join<aod::Collisions, aod::EvSels, aod::FT0MultZeqs,
            o2::aod::CentFT0Cs, aod::TPCMults, o2::aod::BarrelMults>;
using BCsRun3 =
  soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using TracksSel =
  soa::Join<aod::FullTracks, aod::TrackSelection, aod::TracksDCA>;
using SimCollisions = soa::Join<aod::Collisions, aod::EvSels,
                                aod::McCollisionLabels, o2::aod::CentFT0Cs>;
using SimTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra,
                            aod::TracksDCA, aod::McTrackLabels>;
} // namespace o2::aod

struct UccZdc {
  // Event selection
  Configurable<float> posZcut{"posZcut", +10.0, "z-vertex position cut"};
  Configurable<float> minT0CcentCut{"minT0CcentCut", 0.0, "Min T0C Cent. cut"};
  Configurable<float> maxT0CcentCut{"maxT0CcentCut", 90.0, "Max T0C Cent. cut"};

  // Track selection settings
  // Configurable<int> minItsNclusters{"minItsNclusters", 5,
  //                                   "minimum number of ITS clusters"};
  // Configurable<int> minTpcNclusters{"minTpcNclusters", 70,
  //                                   "minimum number of TPC clusters"};
  // Configurable<int> minTpcNcrossedRows{
  //     "minTpcNcrossedRows", 70, "minimum number of TPC crossed pad rows"};
  // Configurable<double> maxChiSquareTpc{"maxChiSquareTpc", 4.0,
  //                                      "maximum TPC chi^2/Ncls"};
  // Configurable<double> maxChiSquareIts{"maxChiSquareIts", 36.0,
  //                                      "maximum ITS chi^2/Ncls"};
  Configurable<double> minPt{"minPt", 0.1, "minimum pt of the tracks"};
  Configurable<double> maxPt{"maxPt", 50., "maximum pt of the tracks"};
  Configurable<double> minEta{"minEta", -0.8, "minimum eta"};
  Configurable<double> maxEta{"maxEta", +0.8, "maximum eta"};
  // Configurable<double> maxDcaxy{"maxDcaxy", 0.05, "Maximum DCAxy"};
  // Configurable<double> maxDcaz{"maxDcaz", 0.05, "Maximum DCAz"};
  // Configurable<bool> setDCAselectionPtDep{"setDCAselectionPtDep", true,
  //                                         "require pt dependent selection"};
  // Configurable<double> par0{"par0", 0.0105, "par 0"};
  // Configurable<double> par1{"par1", 0.035, "par 1"};
  // Configurables, binning
  Configurable<int> nBinsAmpFV0{"nBinsAmpFV0", 1000, "N bins FV0 amp"};
  Configurable<float> maxAmpFV0{"maxAmpFV0", 3000, "Max FV0 amp"};
  Configurable<int> nBinsAmpFT0{"nBinsAmpFT0", 1000, "N bins FT0 amp"};
  Configurable<float> maxAmpFT0{"maxAmpFT0", 3000, "Max FT0 amp"};
  Configurable<int> nBinsNch{"nBinsNch", 2500, "N bins Nch (|eta|<0.8)"};
  Configurable<float> maxNch{"maxNch", 2500, "Max Nch (|eta|<0.8)"};
  Configurable<int> nBinsZDC{"nBinsZDC", 1025, "nBinsZDC"};
  Configurable<int> nBinsZEM{"nBinsZEM", 100, "nBinsZEM"};
  Configurable<float> maxZN{"maxZN", 4099.5, "Max ZN signal"};
  Configurable<float> maxZP{"maxZP", 3099.5, "Max ZP signal"};
  Configurable<float> maxZEM{"maxZEM", 3099.5, "Max ZEM signal"};
  Configurable<int> nBinsTDC{"nBinsTDC", 480, "nbinsTDC"};
  Configurable<float> minTdc{"minTdc", -15.0, "minimum TDC"};
  Configurable<float> maxTdc{"maxTdc", 15.0, "maximum TDC"};
  Configurable<float> minMeanpT{"minMeanpT", 0.5, "minimum [pT]"};
  Configurable<float> maxMeanpT{"maxMeanpT", 1.1, "maximum [pT]"};
  Configurable<int> nBinsMeanpT{"nBinsMeanpT", 160, "# bins [pT]"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0., 0.1, 0.25, 0.5, 1., 2., 4., 6., 8., 10., 20.}, "pT binning"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.}, "T0C binning"};

  // Configurable event selectiond and flags ZDC
  Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", true,
                                            "Enable SameBunchPileup cut"};
  Configurable<bool> isApplyGoodZvtxFT0vsPV{"isApplyGoodZvtxFT0vsPV", true,
                                            "Enable GoodZvtxFT0vsPV cut"};
  Configurable<bool> isApplyVertexITSTPC{"isApplyVertexITSTPC", true,
                                         "Enable VertexITSTPC cut"};
  Configurable<bool> isApplyVertexTOFmatched{"isApplyVertexTOFmatched", true,
                                             "Enable VertexTOFmatched cut"};
  Configurable<bool> isAmpZDC{"isAmpZDC", false, "Use amplitude ZDC?"};
  Configurable<bool> isCommPMT{"isCommPMT", false, "Use common PMT ZDC?"};
  Configurable<bool> isSumTowers{"isSumTowers", false, "Use sum of Tow ZDC?"};
  Configurable<bool> isTDCcut{"isTDCcut", false, "Use TDC cut?"};
  Configurable<bool> isZEMcut{"isZEMcut", true, "Use ZEM cut?"};
  Configurable<float> zemCut{"zemCut", 1000.0, "ZEM cut"};
  Configurable<float> tdcCut{"tdcCut", 1.0, "TDC cut"};

  enum EvCutLabel { All = 1,
                    SelEigth,
                    NoSameBunchPileup,
                    IsGoodZvtxFT0vsPV,
                    IsVertexITSTPC,
                    IsVertexTOFmatched,
                    Centrality,
                    VtxZ,
                    CentralityCut,
                    Zdc,
                    TZero,
                    Tdc,
                    Zem };

  // Filters
  Filter collFilter = (nabs(aod::collision::posZ) < posZcut);
  Filter trackFilter = (requireGlobalTrackInFilter());

  // Apply Filters
  // using TheFilteredCollisions = soa::Filtered<o2::aod::ColEvSels>;
  // using TheFilteredCollision = TheFilteredCollisions::iterator;
  using TheFilteredTracks = soa::Filtered<o2::aod::TracksSel>;
  // using TheFilteredTrack = TheFilteredTracks::iterator;

  using TheFilteredSimCollisions = soa::Filtered<o2::aod::SimCollisions>;
  using TheFilteredSimTracks = soa::Filtered<o2::aod::SimTracks>;

  // Histograms: Data
  HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisZpos{48, -12., 12., "Vtx_{z} (cm)"};
    const AxisSpec axisEvent{14, 0.5, 14.5, ""};
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisDeltaPt{100, -1.0, +1.0, "#Delta(p_{T})"};
    const AxisSpec axisCent{binsCent, "T0C centrality"};
    const AxisSpec axisAmpCh{250, 0., 2500., "Amplitude of non-zero channels"};
    const AxisSpec axisEneCh{300, 0., 300., "Energy of non-zero channels"};

    registry.add("zPos", ";;Entries;", kTH1F, {axisZpos});
    registry.add("hEventCounter", ";;Events", kTH1F, {axisEvent});
    auto hstat = registry.get<TH1>(HIST("hEventCounter"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All events");
    x->SetBinLabel(2, "sel8");
    x->SetBinLabel(3, "kNoSameBunchPileup");
    x->SetBinLabel(4, "kIsGoodZvtxFT0vsPV");
    x->SetBinLabel(5, "kIsVertexITSTPC");
    x->SetBinLabel(6, "kIsVertexTOFmatched");
    x->SetBinLabel(7, "Centrality");
    x->SetBinLabel(8, "VtxZ");
    x->SetBinLabel(9, "Centrality cut");
    x->SetBinLabel(10, "has ZDC?");
    x->SetBinLabel(11, "has T0?");
    x->SetBinLabel(12, "inside TDC cut?");
    x->SetBinLabel(13, "within ZEM cut?");

    //  Histograms: paritcle-level info
    if (doprocessZdcCollAss) {
      registry.add("EtaVsPhi", ";#eta;#varphi", kTH2F,
                   {{{axisEta}, {100, -0.1 * PI, +2.1 * PI}}});
      registry.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
      registry.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
      registry.add("dcaXYvspT", "", kTH2F, {{{50, -1., 1.}, {axisPt}}});
      registry.add("T0Ccent", ";T0C centrality;Entries", kTH1F, {axisCent});
      registry.add("Nch", ";Nch (|#eta|<0.8);", kTH1F,
                   {{nBinsNch, -0.5, maxNch}});
      registry.add("ZNA", "", kTH1F, {{nBinsZDC, -0.5, maxZN}});
      registry.add("ZPA", "", kTH1F, {{nBinsZDC, -0.5, maxZP}});
      registry.add("ZNC", "", kTH1F, {{nBinsZDC, -0.5, maxZN}});
      registry.add("ZPC", "", kTH1F, {{nBinsZDC, -0.5, maxZP}});
      registry.add("ZNvsZEM", "ZNvsZEM; ZEM; ZNA+ZNC", kTH2F,
                   {{{nBinsZDC, -0.5, maxZEM}, {nBinsZEM, -0.5, maxZN}}});
      registry.add("ZNAvsZNC", "ZNAvsZNC; ZNC; ZNA", kTH2F,
                   {{{nBinsZDC, -0.5, maxZN}, {nBinsZDC, -0.5, maxZN}}});
      registry.add("ZPAvsZPC", "ZPAvsZPC; ZPA; ZPC", kTH2F,
                   {{{nBinsZDC, -0.5, maxZP}, {nBinsZDC, -0.5, maxZP}}});
      registry.add("ZNAvsZPA", "ZNAvsZPA; ZPA; ZNA", kTH2F,
                   {{{nBinsZDC, -0.5, maxZP}, {nBinsZDC, -0.5, maxZN}}});
      registry.add("ZNCvsZPC", "ZNCvsZPC; ZPC; ZNC", kTH2F,
                   {{{nBinsZDC, -0.5, maxZP}, {nBinsZDC, -0.5, maxZN}}});
      registry.add("ZNCvstdc", "ZNCvstdc; time ZNC; ZNC", kTH2F,
                   {{{nBinsTDC, minTdc, maxTdc}, {nBinsZDC, -0.5, maxZN}}});
      registry.add("ZNAvstdc", "ZNAvstdc; time ZNA; ZNA", kTH2F,
                   {{{nBinsTDC, minTdc, maxTdc}, {nBinsZDC, -0.5, maxZN}}});
      registry.add("ZPCvstdc", "ZPCvstdc; time ZPC; ZPC", kTH2F,
                   {{{nBinsTDC, minTdc, maxTdc}, {nBinsZDC, -0.5, maxZP}}});
      registry.add("ZPAvstdc", "ZPAvstdc; time ZPA; ZPA", kTH2F,
                   {{{nBinsTDC, minTdc, maxTdc}, {nBinsZDC, -0.5, maxZP}}});
      registry.add("ZEM1vstdc", "ZEM1vstdc; time ZEM1; ZEM1", kTH2F,
                   {{{nBinsTDC, minTdc, maxTdc}, {nBinsZEM, -0.5, maxZEM}}});
      registry.add("ZEM2vstdc", "ZEM2vstdc; time ZEM2; ZEM2", kTH2F,
                   {{{nBinsTDC, minTdc, maxTdc}, {nBinsZEM, -0.5, maxZEM}}});
      registry.add("debunch", ";t_{ZDC}-t_{ZDA};t_{ZDC}+t_{ZDA}", kTH2F,
                   {{{nBinsTDC, minTdc, maxTdc}, {nBinsTDC, minTdc, maxTdc}}});
      registry.add("NchvsFT0C", ";T0C;N_{ch} (|#eta|<0.8);", kTH2F,
                   {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsNch, -0.5, maxNch}}});
      registry.add("NchvsFT0A", ";T0A;N_{ch} (|#eta|<0.8);", kTH2F,
                   {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsNch, -0.5, maxNch}}});
      registry.add("NchvsFV0A", ";V0A;N_{ch} (|#eta|<0.8);", kTH2F,
                   {{{nBinsAmpFV0, 0., maxAmpFV0}, {nBinsNch, -0.5, maxNch}}});
      registry.add("NchvsNPV", ";NPVTracks (|#eta|<1);N_{ch} (|#eta|<0.8);",
                   kTH2F,
                   {{{nBinsNch, -0.5, maxNch}, {nBinsNch, -0.5, maxNch}}});
      registry.add("ZNCvsNch", ";Nch (|#eta|<0.8);ZNC", kTH2F,
                   {{{nBinsNch, -0.5, maxNch}, {nBinsZDC, -0.5, maxZN}}});
      registry.add("ZNAvsNch", ";Nch (|#eta|<0.8);ZNA", kTH2F,
                   {{{nBinsNch, -0.5, maxNch}, {nBinsZDC, -0.5, maxZN}}});
      registry.add("ZNCvsNchvspT", ";Nch (|#eta|<0.8);ZNC;[p_{T}]", kTH3F,
                   {{{nBinsNch, -0.5, maxNch},
                     {nBinsZDC, -0.5, maxZN},
                     {nBinsMeanpT, minMeanpT, maxMeanpT}}});
      registry.add("ZNAvsNchvspT", ";Nch (|#eta|<0.8);ZNA;[p_{T}]", kTH3F,
                   {{{nBinsNch, -0.5, maxNch},
                     {nBinsZDC, -0.5, maxZN},
                     {nBinsMeanpT, minMeanpT, maxMeanpT}}});
    }

    // MC Histograms
    if (doprocessMC) {
      registry.add("hT0C_cent_rec", ";T0C centrality;Entries", kTH1F,
                   {axisCent});
      registry.add("Pt_MC_rec_all_ch", ";;Entries;", kTH1F, {axisPt});
      registry.add("Pt_MC_rec_ch", ";;Entries;", kTH1F, {axisPt});
      registry.add("Pt_MC_rec_pi", ";;Entries;", kTH1F, {axisPt});
      registry.add("Pt_MC_rec_ka", ";;Entries;", kTH1F, {axisPt});
      registry.add("Pt_MC_rec_pr", ";;Entries;", kTH1F, {axisPt});
      registry.add("Pt_MC_rec_sigpos", ";;Entries;", kTH1F, {axisPt});
      registry.add("Pt_MC_rec_signeg", ";;Entries;", kTH1F, {axisPt});
      registry.add("Pt_MC_rec_re", ";;Entries;", kTH1F, {axisPt});
      registry.add("EtaVsPhi_MC_rec", ";;#varphi;", kTH2F,
                   {{{axisEta}, {100, -0.1 * PI, +2.1 * PI}}});

      // registry.add("numberOfRecoCollisions", "",
      //              {HistType::kTH1F, {{6, -0.5, 5.5}}});
      registry.add("hEventCounter_MC", "Event counter", kTH1F, {axisEvent});
      registry.add("zPos_MC", ";;Entries;", kTH1F, {axisZpos});
      registry.add("aV0Avsb", ";V0A amplitude;Impact parameter", kTH2F,
                   {{{nBinsAmpFV0, 0., maxAmpFV0}, {19, 0., 18.}}});
      registry.add("aT0Avsb", ";T0A amplitude; Impact parameter", kTH2F,
                   {{{nBinsAmpFT0, 0., maxAmpFT0}, {19, 0., 18.}}});
      registry.add("aT0Cvsb", ";T0C amplitude; Impact parameter", kTH2F,
                   {{{nBinsAmpFT0, 0., maxAmpFT0}, {19, -0.5, 18.5}}});
      registry.add("Pt_MC_tru_ch", ";p_{T};Entries;", kTH1F, {axisPt});
      registry.add("Pt_MC_tru_pi", ";p_{T};Entries;", kTH1F, {axisPt});
      registry.add("Pt_MC_tru_ka", ";p_{T};Entries;", kTH1F, {axisPt});
      registry.add("Pt_MC_tru_pr", ";p_{T};Entries;", kTH1F, {axisPt});
      registry.add("Pt_MC_tru_sigpos", "#Sigma^{+};p_{T};Entries;", kTH1F,
                   {axisPt});
      registry.add("Pt_MC_tru_signeg", "#Sigma^{-};p_{T};Entries;", kTH1F,
                   {axisPt});
      registry.add("Pt_MC_tru_re", "Remaining ch particles;p_{T};Entries;",
                   kTH1F, {axisPt});
    }
  }

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    registry.fill(HIST("hEventCounter"), EvCutLabel::All);
    if (!col.sel8()) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::SelEigth);

    if (isApplySameBunchPileup &&
        !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::NoSameBunchPileup);

    if (isApplyGoodZvtxFT0vsPV &&
        !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::IsGoodZvtxFT0vsPV);

    if (isApplyVertexITSTPC &&
        !col.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::IsVertexITSTPC);

    if (isApplyVertexTOFmatched &&
        !col.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::IsVertexTOFmatched);

    // if (isApplyVertexTRDmatched &&
    //     !col.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
    //   return false;
    // }
    // histos.fill(HIST("EventHist"), 7);

    if (col.centFT0C() < 0. || col.centFT0C() > 100.) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::Centrality);

    // Z-vertex position cut
    if (std::fabs(col.posZ()) > posZcut) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::VtxZ);

    // T0C centrality cut
    if (col.centFT0C() < minT0CcentCut || col.centFT0C() > maxT0CcentCut) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::CentralityCut);

    // if (isApplyExtraCorrCut && col.multNTracksPV() > npvTracksCut &&
    //     col.multFT0C() < (10 * col.multNTracksPV() - ft0cCut)) {
    //   return false;
    // }
    // histos.fill(HIST("EventHist"), 9);
    //
    // if (isApplyNoCollInTimeRangeStandard &&
    //     !col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
    //   return false;
    // }
    // histos.fill(HIST("EventHist"), 10);
    //
    // if (isApplyNoCollInRofStandard &&
    //     !col.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
    //   return false;
    // }
    // histos.fill(HIST("EventHist"), 11);
    //
    // if (isApplyNoHighMultCollInPrevRof &&
    //     !col.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
    //   return false;
    // }
    // histos.fill(HIST("EventHist"), 12);
    return true;
  }

  void processZdcCollAss(o2::aod::ColEvSels::iterator const& collision,
                         o2::aod::BCsRun3 const& /*bcs*/,
                         aod::Zdcs const& /*zdcs*/, aod::FV0As const& /*fv0as*/,
                         aod::FT0s const& /*ft0s*/,
                         TheFilteredTracks const& tracks)
  {
    if (!isEventSelected(collision)) {
      return;
    }

    const auto& foundBC = collision.foundBC_as<o2::aod::BCsRun3>();
    if (!foundBC.has_zdc()) { // has ZDC?
      return;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::Zdc);

    float aT0A{0.0};
    float aT0C{0.0};
    float aV0A{0.0};
    float znA{0.0};
    float znC{0.0};
    float zpA{0.0};
    float zpC{0.0};
    float aZEM1{0.0};
    float aZEM2{0.0};
    float tZEM1{0.0};
    float tZEM2{0.0};
    float tZNA{0.0};
    float tZNC{0.0};
    float tZPA{0.0};
    float tZPC{0.0};
    float sumZNs{0.0};
    float sumZEMs{0.0};
    float tZDCdif{0.0};
    float tZDCsum{0.0};

    aZEM1 = foundBC.zdc().amplitudeZEM1();
    aZEM2 = foundBC.zdc().amplitudeZEM2();
    tZEM1 = foundBC.zdc().timeZEM1();
    tZEM2 = foundBC.zdc().timeZEM2();
    tZNA = foundBC.zdc().timeZNA();
    tZNC = foundBC.zdc().timeZNC();
    tZPA = foundBC.zdc().timeZPA();
    tZPC = foundBC.zdc().timeZPC();
    tZDCdif = tZNC + tZPC - tZNA - tZPA;
    tZDCsum = tZNC + tZPC + tZNA + tZPA;

    if (isAmpZDC) {
      znA = foundBC.zdc().amplitudeZNA();
      znC = foundBC.zdc().amplitudeZNC();
      zpA = foundBC.zdc().amplitudeZPA();
      zpC = foundBC.zdc().amplitudeZPC();
    } else if (isCommPMT) {
      znA = foundBC.zdc().energyCommonZNA();
      znC = foundBC.zdc().energyCommonZNC();
      zpA = foundBC.zdc().energyCommonZPA();
      zpC = foundBC.zdc().energyCommonZPC();
    } else if (isSumTowers) {
      for (const auto& eZNA : foundBC.zdc().energySectorZNA())
        znA += eZNA;
      for (const auto& eZNC : foundBC.zdc().energySectorZNC())
        znC += eZNC;
      for (const auto& eZPA : foundBC.zdc().energySectorZPA())
        zpA += eZPA;
      for (const auto& eZPC : foundBC.zdc().energySectorZPC())
        zpC += eZPC;
    } else {
      znA = -999.;
      znC = -999.;
      zpA = -999.;
      zpC = -999.;
    }
    sumZNs = znA + znC;
    sumZEMs = aZEM1 + aZEM2;

    if (foundBC.has_ft0()) {
      for (const auto& amplitude : foundBC.ft0().amplitudeA()) {
        aT0A += amplitude;
      }
      for (const auto& amplitude : foundBC.ft0().amplitudeC()) {
        aT0C += amplitude;
      }
    } else {
      return;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::TZero);

    if (foundBC.has_fv0a()) {
      for (const auto& amplitude : foundBC.fv0a().amplitude()) {
        aV0A += amplitude;
      }
    } else {
      aV0A = -999.;
    }

    // TDC cut
    if (isTDCcut) {
      if (std::sqrt(std::pow(tZDCdif, 2.) + std::pow(tZDCsum, 2.)) > tdcCut) {
        return;
      }
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::Tdc);

    // ZEM cut
    if (isZEMcut) {
      if (sumZEMs < zemCut) {
        return;
      }
    }

    registry.fill(HIST("hEventCounter"), EvCutLabel::Zem);
    registry.fill(HIST("zPos"), collision.posZ());
    registry.fill(HIST("T0Ccent"), collision.centFT0C());
    registry.fill(HIST("ZNA"), znA);
    registry.fill(HIST("ZNC"), znC);
    registry.fill(HIST("ZPA"), zpA);
    registry.fill(HIST("ZPC"), zpC);
    registry.fill(HIST("ZNAvsZNC"), znC, znA);
    registry.fill(HIST("ZNAvsZPA"), zpA, znA);
    registry.fill(HIST("ZNCvsZPC"), zpC, znC);
    registry.fill(HIST("ZPAvsZPC"), zpC, zpA);
    registry.fill(HIST("ZNvsZEM"), sumZEMs, sumZNs);
    registry.fill(HIST("ZNCvstdc"), tZNC, znC);
    registry.fill(HIST("ZNAvstdc"), tZNA, znA);
    registry.fill(HIST("ZPCvstdc"), tZPC, zpC);
    registry.fill(HIST("ZPAvstdc"), tZPA, zpA);
    registry.fill(HIST("ZEM1vstdc"), tZEM1, aZEM1);
    registry.fill(HIST("ZEM2vstdc"), tZEM2, aZEM2);
    registry.fill(HIST("debunch"), tZDCdif, tZDCsum);

    float p1{0.0};
    // float p2{0.0};
    float oneParCorr{0.0};
    // float twoParCorr{0.0};
    const int64_t nch{tracks.size()};
    for (const auto& track : tracks) {
      // Track Selection
      if (!track.isGlobalTrack()) {
        continue;
      }
      // if (track.eta() < minEta || track.eta() > maxEta) {
      //   continue;
      // }
      // if (track.pt() < minPt || track.pt() > maxPt) {
      //   continue;
      // }

      float pt{track.pt()};
      p1 += pt;
      // p2 += std::pow(pt, 2.);

      registry.fill(HIST("EtaVsPhi"), track.eta(), track.phi());
      registry.fill(HIST("etaHistogram"), track.eta());
      registry.fill(HIST("ptHistogram"), track.pt());
      registry.fill(HIST("dcaXYvspT"), track.dcaXY(), track.pt());
    }

    oneParCorr = p1;
    // twoParCorr = std::pow(p1, 2.) - p2;
    // std::cout << "twoParCorr= " << twoParCorr << '\n';

    if (nch > 0) {
      oneParCorr /= nch;
    }

    registry.fill(HIST("NchvsFV0A"), aV0A / 100., nch);
    registry.fill(HIST("NchvsFT0A"), aT0A / 100., nch);
    registry.fill(HIST("NchvsFT0C"), aT0C / 100., nch);
    registry.fill(HIST("NchvsNPV"), collision.multNTracksPVeta1(), nch);
    registry.fill(HIST("Nch"), nch);
    registry.fill(HIST("ZNAvsNch"), nch, znA);
    registry.fill(HIST("ZNCvsNch"), nch, znC);
    registry.fill(HIST("ZNCvsNchvspT"), nch, znC, oneParCorr);
    registry.fill(HIST("ZNAvsNchvspT"), nch, znA, oneParCorr);
  }
  PROCESS_SWITCH(UccZdc, processZdcCollAss, "Process ZDC W/Coll Ass.", true);

  Preslice<aod::McParticles> perMCCollision = aod::mcparticle::mcCollisionId;
  Preslice<TheFilteredSimTracks> perCollision = aod::track::collisionId;
  void processMC(aod::McCollisions const& mcCollisions,
                 o2::aod::BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcs*/,
                 aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/,
                 o2::aod::SimCollisions const& collisions,
                 aod::McParticles const& mcParticles,
                 TheFilteredSimTracks const& simTracks)
  {
    // Generated MC
    for (const auto& mccollision : mcCollisions) {
      registry.fill(HIST("hEventCounter_MC"), EvCutLabel::All);
      // Z-vtx position cut
      if (std::fabs(mccollision.posZ()) > posZcut) {
        continue;
      }
      registry.fill(HIST("zPos_MC"), mccollision.posZ());
      registry.fill(HIST("hEventCounter_MC"), EvCutLabel::VtxZ);

      auto mcParticlesPerColl =
        mcParticles.sliceBy(perMCCollision, mccollision.globalIndex());

      for (const auto& particle : mcParticlesPerColl) {
        if (particle.eta() < minEta || particle.eta() > maxEta) {
          continue;
        }
        if (particle.pt() < minPt || particle.pt() > maxPt) {
          continue;
        }
        if (!particle.isPhysicalPrimary()) {
          continue;
        }
        registry.fill(HIST("Pt_MC_tru_ch"), particle.pt());
        if (particle.pdgCode() == PDG_t::kPiPlus ||
            particle.pdgCode() == PDG_t::kPiMinus) { // pion
          registry.fill(HIST("Pt_MC_tru_pi"), particle.pt());
        } else if (particle.pdgCode() == PDG_t::kKPlus ||
                   particle.pdgCode() == PDG_t::kKMinus) { // kaon
          registry.fill(HIST("Pt_MC_tru_ka"), particle.pt());
        } else if (particle.pdgCode() == PDG_t::kProton ||
                   particle.pdgCode() == PDG_t::kProtonBar) { // proton
          registry.fill(HIST("Pt_MC_tru_pr"), particle.pt());
        } else if (particle.pdgCode() == PDG_t::kSigmaPlus ||
                   particle.pdgCode() ==
                     PDG_t::kSigmaBarMinus) { // positive sigma
          registry.fill(HIST("Pt_MC_tru_sigpos"), particle.pt());
        } else if (particle.pdgCode() == PDG_t::kSigmaMinus ||
                   particle.pdgCode() ==
                     PDG_t::kSigmaBarPlus) { // negative sigma
          registry.fill(HIST("Pt_MC_tru_signeg"), particle.pt());
        } else { // rest
          registry.fill(HIST("Pt_MC_tru_re"), particle.pt());
        }
      }
    }
    // registry.fill(HIST("numberOfRecoCollisions"), collisions.size());
    // if (collisions.size() == 0 || collisions.size() > 1) {
    //   return;
    // }
    //----- MC reconstructed -----//
    for (const auto& collision : collisions) {
      // Event selection
      if (!isEventSelected(collision)) {
        continue;
      }
      // MC collision?
      if (!collision.has_mcCollision()) {
        continue;
      }
      const auto& mccollision = collision.mcCollision_as<aod::McCollisions>();
      const auto& foundBC = collision.foundBC_as<o2::aod::BCsRun3>();
      registry.fill(HIST("zPos"), collision.posZ());

      float aT0A{0.0};
      float aT0C{0.0};
      float aV0A{0.0};
      float b{mccollision.impactParameter()};
      if (foundBC.has_ft0()) {
        for (const auto& amplitude : foundBC.ft0().amplitudeA()) {
          aT0A += amplitude;
        }
        for (const auto& amplitude : foundBC.ft0().amplitudeC()) {
          aT0C += amplitude;
        }
      } else {
        aT0A = aT0C = -999;
      }
      if (foundBC.has_fv0a()) {
        for (const auto& amplitude : foundBC.fv0a().amplitude()) {
          aV0A += amplitude;
        }
      } else {
        aV0A = -999;
      }

      registry.fill(HIST("hT0C_cent_rec"), collision.centFT0C());
      registry.fill(HIST("aT0Avsb"), aT0A / 100., b);
      registry.fill(HIST("aT0Cvsb"), aT0C / 100., b);
      registry.fill(HIST("aV0Avsb"), aV0A / 100., b);

      const auto groupedTracks =
        simTracks.sliceBy(perCollision, collision.globalIndex());

      for (const auto& track : groupedTracks) {
        if (!track.has_mcParticle()) {
          continue;
        }
        const auto particle = track.mcParticle();

        registry.fill(HIST("Pt_MC_rec_all_ch"), track.pt());
        if (!particle.isPhysicalPrimary()) {
          continue;
        }

        registry.fill(HIST("EtaVsPhi_MC_rec"), track.eta(), track.phi());
        registry.fill(HIST("Pt_MC_rec_ch"), track.pt());
        if (particle.pdgCode() == PDG_t::kPiPlus ||
            particle.pdgCode() == PDG_t::kPiMinus) {
          registry.fill(HIST("Pt_MC_rec_pi"), track.pt());
        } else if (particle.pdgCode() == PDG_t::kKPlus ||
                   particle.pdgCode() == PDG_t::kKMinus) {
          registry.fill(HIST("Pt_MC_rec_ka"), track.pt());
        } else if (particle.pdgCode() == PDG_t::kProton ||
                   particle.pdgCode() == PDG_t::kProtonBar) {
          registry.fill(HIST("Pt_MC_rec_pr"), track.pt());
        } else if (particle.pdgCode() == PDG_t::kSigmaPlus ||
                   particle.pdgCode() == PDG_t::kSigmaBarMinus) {
          registry.fill(HIST("Pt_MC_rec_sigpos"), track.pt());
        } else if (particle.pdgCode() == PDG_t::kSigmaMinus ||
                   particle.pdgCode() == PDG_t::kSigmaBarPlus) {
          registry.fill(HIST("Pt_MC_rec_signeg"), track.pt());
        } else {
          registry.fill(HIST("Pt_MC_rec_re"), track.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(UccZdc, processMC, "Process MC", false);

  // Single-Track Selection
  // template <typename T2>
  // bool passedTrackSelection(const T2& track) {
  //   if (track.eta() < minEta || track.eta() > maxEta) return false;
  //   if (track.pt() < minPt) return false;
  //
  //   if (!track.hasITS()) return false;
  //   if (track.itsNCls() < minItsNclusters) return false;
  //   if (!track.hasTPC()) return false;
  //   if (track.tpcNClsFound() < minTpcNclusters) return false;
  //   if (track.tpcNClsCrossedRows() < minTpcNcrossedRows) return false;
  //   if (track.tpcChi2NCl() > maxChiSquareTpc) return false;
  //   if (track.itsChi2NCl() > maxChiSquareIts) return false;
  //   // pt-dependent selection
  //   if (setDCAselectionPtDep) {
  //     if (std::fabs(track.dcaXY()) > (par0 + par1 / track.pt())) return
  //     false; if (std::fabs(track.dcaZ()) > (par0 + par1 / track.pt()))
  //     return false;
  //   }
  //   // standard selection
  //   if (!setDCAselectionPtDep) {
  //     if (std::fabs(track.dcaXY()) > maxDcaxy) return false;
  //     if (std::fabs(track.dcaZ()) > maxDcaz) return false;
  //   }
  //   return true;
  // }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UccZdc>(cfgc)};
}
