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

/// \file mcParticleHistManager.h
/// \brief histogram manager for mc particle histograms
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTO_CORE_MCPARTICLEHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_MCPARTICLEHISTMANAGER_H_

#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <TH1.h>

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto
{
namespace mcparticlehistmanager
{

// enum for mc particle histograms
enum McParticleHist {
  // kinemtics
  kPt,
  kEta,
  kPhi,
  kSign,
  // pdg codes
  kPdg,
  kPdgMother,
  kPdgPartonicMother,
  // mother kinematics
  kMotherPt,
  kMotherEta,
  kMotherPhi,
  kMotherSign,
  // origin (kept in sync with TrackHist, even though some categories
  // are not expected to be populated for pure MC truth particles)
  kOrigin,
  kNoMcParticle,
  kPrimary,
  kFromWrongCollision,
  kFromMaterial,
  kMissidentified,
  kSecondary1,
  kSecondary2,
  kSecondary3,
  kSecondaryOther,

  kMcParticleHistLast
};

constexpr std::size_t MaxSecondary = 3;

template <const char* Prefix>
struct ConfMcParticleBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};
  o2::framework::ConfigurableAxis pdgCodes{"pdgCodes", {{8001, -4000.5, 4000.5}}, "PDG codes of selected particles"};
  o2::framework::Configurable<bool> plotOrigins{"plotOrigins", true, "Plot pt distributions for different particle origins"};
  o2::framework::Configurable<std::vector<int>> pdgCodesForMothersOfSecondary{"pdgCodesForMothersOfSecondary", {3122}, "PDG codes of mothers of secondaries (Max 3 will be considered)"};
};

constexpr const char PrefixMcParticleBinning1[] = "McParticleBinning1";
constexpr const char PrefixMcParticleBinning2[] = "McParticleBinning2";

using ConfMcParticleBinning1 = ConfMcParticleBinning<PrefixMcParticleBinning1>;
using ConfMcParticleBinning2 = ConfMcParticleBinning<PrefixMcParticleBinning2>;

// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<McParticleHist>, kMcParticleHistLast>
  HistTable = {
    {
      {kPt, o2::framework::HistType::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
      {kEta, o2::framework::HistType::kTH1F, "hEta", "Pseudorapidity; #eta; Entries"},
      {kPhi, o2::framework::HistType::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
      {kSign, o2::framework::HistType::kTH1F, "hSign", "Sign of charge; Sign; Entries"},
      {kPdg, o2::framework::HistType::kTH1F, "hPdg", "PDG code of the particle; PDG Code; Entries"},
      {kPdgMother, o2::framework::HistType::kTH1F, "hPdgMother", "PDG code of the mother; PDG Code; Entries"},
      {kPdgPartonicMother, o2::framework::HistType::kTH1F, "hPdgPartonicMother", "PDG code of the partonic mother; PDG Code; Entries"},
      {kMotherPt, o2::framework::HistType::kTH1F, "hMotherPt", "Mother transverse momentum; p_{T} (GeV/#it{c}); Entries"},
      {kMotherEta, o2::framework::HistType::kTH1F, "hMotherEta", "Mother pseudorapidity; #eta; Entries"},
      {kMotherPhi, o2::framework::HistType::kTH1F, "hMotherPhi", "Mother azimuthal angle; #varphi; Entries"},
      {kMotherSign, o2::framework::HistType::kTH1F, "hMotherSign", "Sign of mother charge; Sign; Entries"},
      {kOrigin, o2::framework::HistType::kTH1F, "hOrigin", "Status Codes (=Origin); Status Code; Entries"},
      {kNoMcParticle, o2::framework::HistType::kTH1F, "hNoMcParticle", "Particles with no associated MC information; p_{T} (GeV/#it{c}); Entries"},
      {kPrimary, o2::framework::HistType::kTH1F, "hPrimary", "Primary particles; p_{T} (GeV/#it{c}); Entries"},
      {kFromWrongCollision, o2::framework::HistType::kTH1F, "hFromWrongCollision", "Particles associated to wrong collision; p_{T} (GeV/#it{c}); Entries"},
      {kFromMaterial, o2::framework::HistType::kTH1F, "hFromMaterial", "Particles from material; p_{T} (GeV/#it{c}); Entries"},
      {kMissidentified, o2::framework::HistType::kTH1F, "hMissidentified", "Missidentified particles (fake/wrong PDG code); p_{T} (GeV/#it{c}); Entries"},
      {kSecondary1, o2::framework::HistType::kTH1F, "hFromSecondary1", "Particles from secondary decay; p_{T} (GeV/#it{c}); Entries"},
      {kSecondary2, o2::framework::HistType::kTH1F, "hFromSecondary2", "Particles from secondary decay; p_{T} (GeV/#it{c}); Entries"},
      {kSecondary3, o2::framework::HistType::kTH1F, "hFromSecondary3", "Particles from secondary decay; p_{T} (GeV/#it{c}); Entries"},
      {kSecondaryOther, o2::framework::HistType::kTH1F, "hFromSecondaryOther", "Particles from every other secondary decay; p_{T} (GeV/#it{c}); Entries"},
    }};

template <typename T>
auto makeMcParticleHistSpecMap(const T& confBinning)
{
  return std::map<McParticleHist, std::vector<o2::framework::AxisSpec>>{
    {kPt, {confBinning.pt}},
    {kEta, {confBinning.eta}},
    {kPhi, {confBinning.phi}},
    {kSign, {confBinning.sign}},
    {kPdg, {confBinning.pdgCodes}},
    {kPdgMother, {confBinning.pdgCodes}},
    {kPdgPartonicMother, {confBinning.pdgCodes}},
    {kMotherPt, {confBinning.pt}},
    {kMotherEta, {confBinning.eta}},
    {kMotherPhi, {confBinning.phi}},
    {kMotherSign, {confBinning.sign}},
    {kNoMcParticle, {confBinning.pt}},
    {kPrimary, {confBinning.pt}},
    {kFromWrongCollision, {confBinning.pt}},
    {kFromMaterial, {confBinning.pt}},
    {kMissidentified, {confBinning.pt}},
    {kSecondary1, {confBinning.pt}},
    {kSecondary2, {confBinning.pt}},
    {kSecondary3, {confBinning.pt}},
    {kSecondaryOther, {confBinning.pt}},
  };
}

constexpr char PrefixMcParticle1[] = "McParticle1/";
constexpr char PrefixMcParticle2[] = "McParticle2/";

constexpr std::string_view McDir = "MC/";

template <const char* prefix>
class McParticleHistManager
{
 public:
  McParticleHistManager() = default;
  ~McParticleHistManager() = default;

  // init with origin histograms enabled, controlled via confBinning.plotOrigins
  template <typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<McParticleHist, std::vector<o2::framework::AxisSpec>> const& Specs,
            T const& confBinning)
  {
    mHistogramRegistry = registry;

    std::string mcDir = std::string(prefix) + std::string(McDir);
    mHistogramRegistry->add(mcDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {Specs.at(kPt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {Specs.at(kEta)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {Specs.at(kPhi)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kSign, HistTable), getHistDesc(kSign, HistTable), getHistType(kSign, HistTable), {Specs.at(kSign)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdg, HistTable), getHistDesc(kPdg, HistTable), getHistType(kPdg, HistTable), {Specs.at(kPdg)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgMother, HistTable), getHistDesc(kPdgMother, HistTable), getHistType(kPdgMother, HistTable), {Specs.at(kPdgMother)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgPartonicMother, HistTable), getHistDesc(kPdgPartonicMother, HistTable), getHistType(kPdgPartonicMother, HistTable), {Specs.at(kPdgPartonicMother)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kMotherPt, HistTable), getHistDesc(kMotherPt, HistTable), getHistType(kMotherPt, HistTable), {Specs.at(kMotherPt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kMotherEta, HistTable), getHistDesc(kMotherEta, HistTable), getHistType(kMotherEta, HistTable), {Specs.at(kMotherEta)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kMotherPhi, HistTable), getHistDesc(kMotherPhi, HistTable), getHistType(kMotherPhi, HistTable), {Specs.at(kMotherPhi)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kMotherSign, HistTable), getHistDesc(kMotherSign, HistTable), getHistType(kMotherSign, HistTable), {Specs.at(kMotherSign)});

    this->enableOptionalHistograms(confBinning);
    this->initOrigin(Specs);
  }

  template <typename T1, typename T2, typename T3>
  void fill(T1 const& mcParticle, T2 const& /*mcMothers*/, T3 const& /*mcPartonicMothers*/)
  {
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPt, HistTable)), mcParticle.pt());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kEta, HistTable)), mcParticle.eta());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPhi, HistTable)), mcParticle.phi());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kSign, HistTable)), mcParticle.sign());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), mcParticle.pdgCode());

    bool hasMother = mcParticle.has_fMcMother();
    if (hasMother) {
      auto mother = mcParticle.template fMcMother_as<T2>();
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), mother.pdgCode());
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kMotherPt, HistTable)), mother.pt());
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kMotherEta, HistTable)), mother.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kMotherPhi, HistTable)), mother.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kMotherSign, HistTable)), mother.sign());
    } else {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), 0);
    }

    if (mcParticle.has_fMcPartMoth()) {
      auto partonicMother = mcParticle.template fMcPartMoth_as<T3>();
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), partonicMother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), 0);
    }

    if (mPlotOrigins) {
      // NOTE: kNoMcParticle and kMissidentified are kept here only for enum/histogram
      // symmetry with TrackHist — they describe reco-vs-truth mismatches that cannot
      // occur for a pure MC-truth particle, so these bins are expected to stay empty
      // unless the producer assigns those origin codes for some other reason
      // (e.g. pileup-associated particles tagged as kFromWrongCollision).
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), mcParticle.origin());
      switch (static_cast<modes::McOrigin>(mcParticle.origin())) {
        case modes::McOrigin::kNoMcParticle:
          mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kNoMcParticle, HistTable)), mcParticle.pt());
          break;
        case modes::McOrigin::kPhysicalPrimary:
          mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPrimary, HistTable)), mcParticle.pt());
          break;
        case modes::McOrigin::kFromWrongCollision:
          mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kFromWrongCollision, HistTable)), mcParticle.pt());
          break;
        case modes::McOrigin::kFromMaterial:
          mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kFromMaterial, HistTable)), mcParticle.pt());
          break;
        case modes::McOrigin::kMissidentified:
          mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kMissidentified, HistTable)), mcParticle.pt());
          break;
        case modes::McOrigin::kFromSecondaryDecay:
          if (hasMother) {
            auto mother = mcParticle.template fMcMother_as<T2>();
            int motherPdgCode = std::abs(mother.pdgCode());
            if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel1 && motherPdgCode == mPdgCodesSecondaryMother[0]) {
              mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kSecondary1, HistTable)), mcParticle.pt());
            } else if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel2 && motherPdgCode == mPdgCodesSecondaryMother[1]) {
              mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kSecondary2, HistTable)), mcParticle.pt());
            } else if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel3 && motherPdgCode == mPdgCodesSecondaryMother[2]) {
              mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kSecondary3, HistTable)), mcParticle.pt());
            } else {
              mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kSecondaryOther, HistTable)), mcParticle.pt());
            }
          }
          break;
        default:
          LOG(warn) << "Encountered mc particle with unknown origin!";
          break;
      }
    }
  }

 private:
  template <typename T>
  void enableOptionalHistograms(T const& confBinning)
  {
    mPlotOrigins = confBinning.plotOrigins.value;
    mPlotNSecondaries = confBinning.pdgCodesForMothersOfSecondary.value.size();

    for (std::size_t i = 0; i < MaxSecondary; i++) {
      if (i < confBinning.pdgCodesForMothersOfSecondary.value.size()) {
        mPdgCodesSecondaryMother.at(i) = std::abs(confBinning.pdgCodesForMothersOfSecondary.value.at(i));
      } else {
        mPdgCodesSecondaryMother.at(i) = 0;
      }
    }
  }

  void initOrigin(std::map<McParticleHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    if (!mPlotOrigins) {
      return;
    }
    std::string mcDir = std::string(prefix) + std::string(McDir);

    // mc origin axis can be configured here
    const o2::framework::AxisSpec axisOrigin = {static_cast<int>(modes::McOrigin::kMcOriginLast), -0.5, static_cast<double>(modes::McOrigin::kMcOriginLast) - 0.5};
    mHistogramRegistry->add(mcDir + getHistNameV2(kOrigin, HistTable), getHistDesc(kOrigin, HistTable), getHistType(kOrigin, HistTable), {axisOrigin});
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kNoMcParticle), modes::mcOriginToString(modes::McOrigin::kNoMcParticle));
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromWrongCollision), modes::mcOriginToString(modes::McOrigin::kFromWrongCollision));
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kPhysicalPrimary), modes::mcOriginToString(modes::McOrigin::kPhysicalPrimary));
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromSecondaryDecay), modes::mcOriginToString(modes::McOrigin::kFromSecondaryDecay));
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromMaterial), modes::mcOriginToString(modes::McOrigin::kFromMaterial));
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kMissidentified), modes::mcOriginToString(modes::McOrigin::kMissidentified));

    mHistogramRegistry->add(mcDir + getHistNameV2(kNoMcParticle, HistTable), getHistDesc(kNoMcParticle, HistTable), getHistType(kNoMcParticle, HistTable), {Specs.at(kNoMcParticle)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPrimary, HistTable), getHistDesc(kPrimary, HistTable), getHistType(kPrimary, HistTable), {Specs.at(kPrimary)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kFromWrongCollision, HistTable), getHistDesc(kFromWrongCollision, HistTable), getHistType(kFromWrongCollision, HistTable), {Specs.at(kFromWrongCollision)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kFromMaterial, HistTable), getHistDesc(kFromMaterial, HistTable), getHistType(kFromMaterial, HistTable), {Specs.at(kFromMaterial)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kMissidentified, HistTable), getHistDesc(kMissidentified, HistTable), getHistType(kMissidentified, HistTable), {Specs.at(kMissidentified)});

    if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel1) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary1, HistTable), getHistDesc(kSecondary1, HistTable), getHistType(kSecondary1, HistTable), {Specs.at(kSecondary1)});
    }
    if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel2) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary2, HistTable), getHistDesc(kSecondary2, HistTable), getHistType(kSecondary2, HistTable), {Specs.at(kSecondary2)});
    }
    if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel3) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary3, HistTable), getHistDesc(kSecondary3, HistTable), getHistType(kSecondary3, HistTable), {Specs.at(kSecondary3)});
    }
    mHistogramRegistry->add(mcDir + getHistNameV2(kSecondaryOther, HistTable), getHistDesc(kSecondaryOther, HistTable), getHistType(kSecondaryOther, HistTable), {Specs.at(kSecondaryOther)});
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  bool mPlotOrigins = true;
  int mPlotNSecondaries = 0;
  std::array<int, MaxSecondary> mPdgCodesSecondaryMother = {0};
};
} // namespace mcparticlehistmanager
} // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_MCPARTICLEHISTMANAGER_H_
