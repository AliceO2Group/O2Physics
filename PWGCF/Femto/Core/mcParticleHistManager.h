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

  void init(o2::framework::HistogramRegistry* registry,
            std::map<McParticleHist, std::vector<o2::framework::AxisSpec>> const& Specs)
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
  }

  /// T1 = mc particle row — must carry FMcMotherId/FMcPartMothId, i.e. this should be
  ///      a row from o2::soa::Join<FMcParticles, FMcMotherLabels>, not bare FMcParticles
  /// T2 = FMcMothers (now with kinematics: Origin, PdgCode, SignedPt, Eta, Phi)
  /// T3 = FMcPartMoths
  template <typename T1, typename T2, typename T3>
  void fill(T1 const& mcParticle, T2 const& /*mcMothers*/, T3 const& /*mcPartonicMothers*/)
  {
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPt, HistTable)), mcParticle.pt());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kEta, HistTable)), mcParticle.eta());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPhi, HistTable)), mcParticle.phi());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kSign, HistTable)), mcParticle.sign());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), mcParticle.pdgCode());

    if (mcParticle.has_fMcMother()) {
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
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
};
}; // namespace mcparticlehistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_MCPARTICLEHISTMANAGER_H_
