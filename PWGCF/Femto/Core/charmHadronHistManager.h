// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file charmHadronHistManager.h
/// \brief histogram manager for charm hadron histograms
/// \author Igor Ptak, WUT, igor.ptak.stud@pw.edu.pl

#ifndef PWGCF_FEMTO_CORE_CHARMHADRONHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_CHARMHADRONHISTMANAGER_H_

#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/trackHistManager.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <TPDGCode.h>

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto::charmhadronhistmanager 
{
enum CharmHadronHist {
    kPt,
    kEta,
    kPhi,
    kMass,
    kSign,
    kPtVsMass,

    kCharmHadronHistLast
};

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define CHARMHADRON_DEFAULT_BINNING(defaultMassMin, defaultMassMax)                                           \
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};                                     \
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};                             \
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};   \
  o2::framework::ConfigurableAxis mass{"mass", {{200, (defaultMassMin), (defaultMassMax)}}, "Mass"}; \
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};                            \
  o2::framework::ConfigurableAxis charmHadrons{"charmHadrons", {{8001, -4000.5, 4000.5}}, "MC ONLY: CharmHadrons codes of reconstructed D0s"};

template <auto& Prefix>
struct ConfD0Binning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  CHARMHADRON_DEFAULT_BINNING(1.7, 2.0)
};

#undef CHARMHADRON_DEFAULT_BINNING

constexpr const char PrefixD0Binning1[] = "D0Binning1";
using ConfD0Binning1 = ConfD0Binning<PrefixD0Binning1>;


// must be in sync with enum CharmHadronHist
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<CharmHadronHist>, kCharmHadronHistLast> HistTable = {
  {{kPt, o2::framework::HistType::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::HistType::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::HistType::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::HistType::kTH1F, "hMass", "Invariant Mass; m_{Inv} (GeV/#it{c}^{2}); Entries"},
   {kSign, o2::framework::HistType::kTH1F, "hSign", "Sign (-1 -> D0bar, +1 -> D0); sign; Entries"},
   {kPtVsMass, o2::framework::HistType::kTH2F, "hPtVsMass", "Transverse momentum vs invariant mass; p_{T} (GeV/#it{c}); m_{Inv} (GeV/#it{c}^{2})"}},
};

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define CHARMHADRON_HIST_ANALYSIS_MAP(conf) \
  {kPt, {(conf).pt}},                       \
    {kEta, {(conf).eta}},                   \
    {kPhi, {(conf).phi}},                   \
    {kMass, {(conf).mass}},                 \
    {kSign, {(conf).sign}},                 \
    {kPtVsMass, {(conf).pt, (conf).mass}},

template <typename T>
auto makeD0HistSpecMap(const T& confBinningAnalysis)
{
  return std::map<CharmHadronHist, std::vector<o2::framework::AxisSpec>>{
    CHARMHADRON_HIST_ANALYSIS_MAP(confBinningAnalysis)};
}

#undef CHARMHADRON_HIST_ANALYSIS_MAP

// prefixes for the output directories in the histogram registry
constexpr char PrefixD01[] = "D01/";
constexpr char PrefixD02[] = "D02/";

constexpr std::string_view AnalysisDir = "Analysis/";

/// \class CharmHadronHistManager
/// \brief Class for histogramming charm hadron properties
template <auto& charmHadronPrefix,
          auto& prong0Prefix,
          auto& prong1Prefix,
          modes::CharmHadron charmHadron>
class CharmHadronHistManager
{
 public:
  CharmHadronHistManager() = default;
  ~CharmHadronHistManager() = default;

  // init for analysis
  template <modes::Mode mode, typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CharmHadronHist, std::vector<o2::framework::AxisSpec>> const& CharmHadronSpecs,
            T const& ConfCharmHadronSelection,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& Prong0Specs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& Prong1Specs)
  {
    mHistogramRegistry = registry;
    mPdgCode = std::abs(ConfCharmHadronSelection.pdgCodeAbs.value);

    // in PWGHF the prong charge is fixed by the reconstruction: prong0 is always the positive
    // daughter, prong1 the negative one. The D0/D0bar hypothesis only swaps which prong is the
    // pion and which is the kaon, not their charge.
    int prong0PdgCodeAbs = 0;
    int prong1PdgCodeAbs = 0;
    const int absCharge = 1;
    const int signPlus = 1;
    const int signMinus = -1;

    constexpr int PdgD0 = 421; // not defined in ROOT's TPDGCode.h
    if (mPdgCode == PdgD0) {
      if (ConfCharmHadronSelection.sign.value > 0) {
        // D0 -> pi+ K-
        prong0PdgCodeAbs = std::abs(PDG_t::kPiPlus);
        prong1PdgCodeAbs = std::abs(PDG_t::kKMinus);
      } else {
        // D0bar -> K+ pi-
        mPdgCode = -1 * mPdgCode; // switch sign for D0bar
        prong0PdgCodeAbs = std::abs(PDG_t::kKPlus);
        prong1PdgCodeAbs = std::abs(PDG_t::kPiMinus);
      }
    } else {
      LOG(fatal) << "PDG code for charm hadron has to be D0 (421)";
    }

    mProng0Manager.template init<mode>(registry, Prong0Specs, absCharge, signPlus, prong0PdgCodeAbs);
    mProng1Manager.template init<mode>(registry, Prong1Specs, absCharge, signMinus, prong1PdgCodeAbs);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->initAnalysis(CharmHadronSpecs);
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fill(T1 const& charmHadronCandidate, T2 const& tracks)
  {
    auto prong0 = tracks.rawIteratorAt(charmHadronCandidate.posDauId() - tracks.offset());
    mProng0Manager.template fill<mode>(prong0, tracks);
    auto prong1 = tracks.rawIteratorAt(charmHadronCandidate.negDauId() - tracks.offset());
    mProng1Manager.template fill<mode>(prong1, tracks);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->fillAnalysis(charmHadronCandidate);
    }
  }

 private:
  void initAnalysis(std::map<CharmHadronHist, std::vector<o2::framework::AxisSpec>> const& CharmHadronSpecs)
  {
    std::string analysisDir = std::string(charmHadronPrefix) + std::string(AnalysisDir);
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {CharmHadronSpecs.at(kPt)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {CharmHadronSpecs.at(kEta)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {CharmHadronSpecs.at(kPhi)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kMass, HistTable), getHistDesc(kMass, HistTable), getHistType(kMass, HistTable), {CharmHadronSpecs.at(kMass)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kSign, HistTable), getHistDesc(kSign, HistTable), getHistType(kSign, HistTable), {CharmHadronSpecs.at(kSign)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPtVsMass, HistTable), getHistDesc(kPtVsMass, HistTable), getHistType(kPtVsMass, HistTable), {CharmHadronSpecs.at(kPtVsMass)});
  }

  template <typename T>
  void fillAnalysis(T const& charmHadronCandidate)
  {
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPt, HistTable)), charmHadronCandidate.pt());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kEta, HistTable)), charmHadronCandidate.eta());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPhi, HistTable)), charmHadronCandidate.phi());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kMass, HistTable)), charmHadronCandidate.mass());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), charmHadronCandidate.sign());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPtVsMass, HistTable)), charmHadronCandidate.pt(), charmHadronCandidate.mass());
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  int mPdgCode = 0;

  trackhistmanager::TrackHistManager<prong0Prefix> mProng0Manager;
  trackhistmanager::TrackHistManager<prong1Prefix> mProng1Manager;
};
}; // namespace o2::analysis::femto::charmhadronhistmanager
#endif // PWGCF_FEMTO_CORE_CHARMHADRONHISTMANAGER_H_
