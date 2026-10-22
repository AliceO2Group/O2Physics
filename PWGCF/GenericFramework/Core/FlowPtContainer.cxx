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

/// \file FlowPtContainer.cxx
/// \brief Class to handle angular and transverse momentum correlations
/// \author Emil Gorm Nielsen, NBI, emil.gorm.nielsen@cern.ch

#include "FlowPtContainer.h"

#include "PWGCF/GenericFramework/Core/BootstrapProfile.h"
#include "PWGCF/GenericFramework/Core/GFWConfig.h"

#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <TCollection.h>
#include <TH1.h>
#include <TList.h>
#include <TNamed.h>
#include <TString.h>

#include <RtypesCore.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <charconv>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

using namespace o2::analysis::genericframework;
using namespace o2::analysis::genericframework::eventweight;
using o2::analysis::genericframework::flowptcontainer::PtCorrConfig;
using o2::analysis::genericframework::flowptcontainer::PtMomentConfig;
using o2::analysis::genericframework::flowptcontainer::PtMomentGroup;

namespace
{
constexpr std::size_t TwoParticleOrder = 2;
constexpr int ParityDivisor = 2;

bool parseInteger(std::string_view token, int& value)
{
  if (!token.empty() && token.front() == '+') {
    token.remove_prefix(1);
  }
  const auto result = std::from_chars(token.data(), token.data() + token.size(), value);
  return !token.empty() && result.ec == std::errc{} && result.ptr == token.data() + token.size();
}

std::optional<std::vector<std::vector<bool>>> replacePtTerms(std::string& expression, const std::string& head)
{
  std::vector<std::vector<bool>> isPt;
  size_t cursor = 0;
  while (true) {
    auto openBrace = expression.find('{', cursor);
    if (openBrace == std::string::npos) {
      break;
    }
    if (openBrace > 0 && !std::isspace(static_cast<unsigned char>(expression[openBrace - 1]))) {
      expression.insert(openBrace++, 1, ' ');
    }
    const auto closeBrace = expression.find('}', openBrace + 1);
    if (closeBrace == std::string::npos || expression.find('{', openBrace + 1) < closeBrace) {
      LOGF(error, "Invalid pT correlation %s: unmatched brace", head.c_str());
      return std::nullopt;
    }

    std::string terms = expression.substr(openBrace + 1, closeBrace - openBrace - 1);
    std::replace(terms.begin(), terms.end(), ',', ' ');
    std::replace(terms.begin(), terms.end(), ';', ' ');
    std::istringstream stream(terms);
    std::string replacement;
    std::string token;
    std::vector<bool> groupIsPt;
    while (stream >> token) {
      int harmonic = 0;
      const bool tokenIsPt = token == "pt";
      if (!tokenIsPt && !parseInteger(token, harmonic)) {
        LOGF(error, "Invalid pT correlation %s: unknown observable %s", head.c_str(), token.c_str());
        return std::nullopt;
      }
      if (!replacement.empty()) {
        replacement += ' ';
      }
      replacement += tokenIsPt ? "0" : token;
      groupIsPt.push_back(tokenIsPt);
    }
    if (groupIsPt.empty()) {
      LOGF(error, "Invalid pT correlation %s: empty observable group", head.c_str());
      return std::nullopt;
    }

    expression.replace(openBrace + 1, closeBrace - openBrace - 1, replacement);
    cursor = openBrace + replacement.size() + 2;
    isPt.push_back(std::move(groupIsPt));
  }

  if (isPt.empty()) {
    LOGF(error, "Invalid pT correlation %s: no observable groups", head.c_str());
    return std::nullopt;
  }
  return isPt;
}

std::string makeFullTopologyName(const std::string& head, std::uint64_t activeFactors, std::uint64_t allFactors, bool hasFlow, int ptOrder)
{
  if (activeFactors == allFactors) {
    return head;
  }

  std::string name = head + "_full";
  if (hasFlow && (activeFactors & std::uint64_t{1})) {
    name += "_flow";
  }
  const int ptOffset = hasFlow ? 1 : 0;
  for (int i = 0; i < ptOrder; ++i) {
    if (activeFactors & (std::uint64_t{1} << (ptOffset + i))) {
      name += "_pt" + std::to_string(i + 1);
    }
  }
  return name;
}

double binomialCoefficient(int n, int k)
{
  k = std::min(k, n - k);
  double result = 1.;
  for (int i = 1; i <= k; ++i) {
    result *= static_cast<double>(n - k + i) / i;
  }
  return result;
}

std::string makePostprocessName(std::string_view observable, const char* head, int ind)
{
  return std::string{head} + "_" + std::string{observable} + "_" + std::to_string(ind);
}

bool sameMomentConfigs(const std::vector<PtMomentConfig>& first, const std::vector<PtMomentConfig>& second)
{
  if (first.size() != second.size()) {
    return false;
  }
  for (size_t config = 0; config < first.size(); ++config) {
    if (first[config].head != second[config].head || first[config].groups.size() != second[config].groups.size()) {
      return false;
    }
    for (size_t group = 0; group < first[config].groups.size(); ++group) {
      const auto& lhs = first[config].groups[group];
      const auto& rhs = second[config].groups[group];
      if (lhs.region != rhs.region || lhs.order != rhs.order || lhs.meanProfile != rhs.meanProfile) {
        return false;
      }
    }
  }
  return true;
}
} // namespace

FlowPtContainer::FlowPtContainer() : fCMTermList(nullptr),
                                     fCorrList(nullptr),
                                     fCovList(nullptr),
                                     fSubList(nullptr),
                                     fSubCMList(nullptr),
                                     fCumulantList(nullptr),
                                     fCentralMomentList(nullptr),
                                     mpar(0),
                                     nSubevents(0),
                                     fillCounter(0),
                                     fEventWeight(EventWeight::UnityWeight),
                                     fUseCentralMoments(true),
                                     fUseGap(false),
                                     fPtCentralMomentList(nullptr),
                                     fPtCumulantList(nullptr)
{
}
FlowPtContainer::~FlowPtContainer()
{
  delete fCMTermList;
  delete fCorrList;
  delete fCovList;
  delete fSubList;
  delete fSubCMList;
  delete fCumulantList;
  delete fCentralMomentList;
  delete fPtCentralMomentList;
  delete fPtCumulantList;
};
FlowPtContainer::FlowPtContainer(const char* name) : TNamed(name, name),
                                                     fCMTermList(nullptr),
                                                     fCorrList(nullptr),
                                                     fCovList(nullptr),
                                                     fSubList(nullptr),
                                                     fSubCMList(nullptr),
                                                     fCumulantList(nullptr),
                                                     fCentralMomentList(nullptr),
                                                     mpar(0),
                                                     nSubevents(0),
                                                     fillCounter(0),
                                                     fEventWeight(EventWeight::UnityWeight),
                                                     fUseCentralMoments(true),
                                                     fUseGap(false),
                                                     fPtCentralMomentList(nullptr),
                                                     fPtCumulantList(nullptr)
{
}
FlowPtContainer::FlowPtContainer(const char* name, const char* title) : TNamed(name, title),
                                                                        fCMTermList(nullptr),
                                                                        fCorrList(nullptr),
                                                                        fCovList(nullptr),
                                                                        fSubList(nullptr),
                                                                        fSubCMList(nullptr),
                                                                        fCumulantList(nullptr),
                                                                        fCentralMomentList(nullptr),
                                                                        mpar(0),
                                                                        nSubevents(0),
                                                                        fillCounter(0),
                                                                        fEventWeight(EventWeight::UnityWeight),
                                                                        fUseCentralMoments(true),
                                                                        fUseGap(false),
                                                                        fPtCentralMomentList(nullptr),
                                                                        fPtCumulantList(nullptr)
{
}
void FlowPtContainer::initialise(const o2::framework::AxisSpec& axis, const GFWRegions& regions, const GFWCorrConfigs& configs, const int& nsub)
{
  const auto numberOfRegions = regions.GetNames().size();
  const auto numberOfConfigs = configs.GetCorrs().size();
  if (regions.GetEtaMin().size() != numberOfRegions || regions.GetEtaMax().size() != numberOfRegions ||
      regions.GetpTDifs().size() != numberOfRegions || regions.GetBitmasks().size() != numberOfRegions ||
      configs.GetHeads().size() != numberOfConfigs || configs.GetpTDifs().size() != numberOfConfigs) {
    LOGF(error, "Cannot initialize %s: GFW configuration vectors have different sizes", GetName());
    return;
  }
  if (nsub < 0) {
    LOGF(error, "Cannot initialize %s: number of subsamples cannot be negative", GetName());
    return;
  }

  std::vector<double> multiBins = axis.binEdges;
  int numberOfBins = axis.nBins.value_or(0);
  const bool isVariableBinning = numberOfBins <= 0;
  if (isVariableBinning) {
    numberOfBins = multiBins.size() - 1;
  }
  const size_t expectedEdges = isVariableBinning ? static_cast<size_t>(numberOfBins + 1) : 2u;
  if (numberOfBins <= 0 || multiBins.size() != expectedEdges) {
    LOGF(error, "Cannot initialize %s: invalid multiplicity axis", GetName());
    return;
  }

  GFW parser;
  for (size_t i = 0; i < numberOfRegions; ++i) {
    parser.AddRegion(regions.GetNames()[i], regions.GetEtaMin()[i], regions.GetEtaMax()[i], regions.GetpTDifs()[i] ? 2 : 1, regions.GetBitmasks()[i]);
  }
  if (parser.fRegions.size() != numberOfRegions) {
    LOGF(error, "Cannot initialize %s: invalid GFW region", GetName());
    return;
  }

  struct ParsedPtRequest {
    GFW::CorrConfig config;
    std::vector<std::vector<bool>> isPt;
    std::vector<std::vector<int>> factors;
    std::vector<std::pair<size_t, size_t>> ptTerms;
    int ptOrder{};
    bool hasFlow{};
    std::uint64_t allFactors{};
  };

  std::vector<ParsedPtRequest> requests;
  requests.reserve(numberOfConfigs);
  std::unordered_set<std::string> correlationHeads;
  std::vector<int> regionMaxOrders(numberOfRegions, 0);
  int maxPtOrder = 0;
  for (size_t i = 0; i < numberOfConfigs; ++i) {
    if (configs.GetpTDifs()[i] != 0) {
      LOGF(error, "Invalid pT correlation %s: differential configurations are not supported", configs.GetHeads()[i].c_str());
      return;
    }
    std::string expression = configs.GetCorrs()[i];
    auto isPt = replacePtTerms(expression, configs.GetHeads()[i]);
    if (!isPt) {
      return;
    }
    auto config = parser.GetCorrelatorConfig(expression, configs.GetHeads()[i], configs.GetpTDifs()[i] != 0);
    if (config.Regs.size() != isPt->size() || config.Hars.size() != isPt->size() || config.Overlap.size() != isPt->size() || config.ptInd.size() != isPt->size()) {
      LOGF(error, "Invalid pT correlation %s", config.Head.c_str());
      return;
    }
    if (std::any_of(config.ptInd.begin(), config.ptInd.end(), [](int ptIndex) { return ptIndex >= 0; })) {
      LOGF(error, "Invalid pT correlation %s: explicit pT-bin indices are not supported", config.Head.c_str());
      return;
    }
    int ptOrder = 0;
    bool hasFlow = false;
    std::vector<std::pair<size_t, size_t>> ptTerms;
    std::vector<std::vector<int>> factors(isPt->size());
    for (size_t group = 0; group < isPt->size(); ++group) {
      if (config.Regs[group].empty() || config.Hars[group].size() != isPt->at(group).size()) {
        LOGF(error, "Invalid pT correlation %s", config.Head.c_str());
        return;
      }
      const bool groupHasPt = std::any_of(isPt->at(group).begin(), isPt->at(group).end(), [](bool value) { return value; });
      const bool groupHasFlow = std::any_of(isPt->at(group).begin(), isPt->at(group).end(), [](bool value) { return !value; });
      if (groupHasPt && groupHasFlow) {
        LOGF(error, "Invalid pT correlation %s: mixed flow and pt terms in one group are not supported", config.Head.c_str());
        return;
      }
      if (groupHasPt && (config.Regs[group].size() != 1 || config.Overlap[group] >= 0)) {
        LOGF(error, "Invalid pT correlation %s: pt groups must contain exactly one region without an explicit overlap", config.Head.c_str());
        return;
      }
      factors[group].resize(isPt->at(group).size());
      for (size_t term = 0; term < isPt->at(group).size(); ++term) {
        if (isPt->at(group)[term]) {
          ptTerms.emplace_back(group, term);
          ++ptOrder;
        } else {
          hasFlow = true;
        }
      }
    }
    if (ptOrder == 0) {
      LOGF(error, "Invalid pT correlation %s: configuration does not contain pt", config.Head.c_str());
      return;
    }
    if (config.Head.empty() || !correlationHeads.insert(config.Head).second) {
      LOGF(error, "Cannot initialize %s: correlation heads must be non-empty and unique", GetName());
      return;
    }

    const int factorCount = ptOrder + (hasFlow ? 1 : 0);
    if (factorCount >= std::numeric_limits<std::uint64_t>::digits) {
      LOGF(error, "Invalid pT correlation %s: too many observable factors", config.Head.c_str());
      return;
    }
    const int ptOffset = hasFlow ? 1 : 0;
    for (int ptIndex = 0; ptIndex < ptOrder; ++ptIndex) {
      const auto [group, term] = ptTerms[ptIndex];
      factors[group][term] = ptOffset + ptIndex;
    }
    if (hasFlow) {
      for (size_t group = 0; group < isPt->size(); ++group) {
        if (!isPt->at(group).front()) {
          std::fill(factors[group].begin(), factors[group].end(), 0);
        }
      }
    }
    const std::uint64_t allFactors = (std::uint64_t{1} << factorCount) - 1;
    for (size_t group = 0; group < isPt->size(); ++group) {
      if (!isPt->at(group).front()) {
        continue;
      }
      const int order = isPt->at(group).size();
      for (const int& region : config.Regs[group]) {
        regionMaxOrders[region] = std::max(regionMaxOrders[region], order);
      }
      if (config.Overlap[group] >= 0) {
        regionMaxOrders[config.Overlap[group]] = std::max(regionMaxOrders[config.Overlap[group]], order);
      }
    }

    maxPtOrder = std::max(maxPtOrder, ptOrder);
    requests.push_back({std::move(config), std::move(*isPt), std::move(factors), std::move(ptTerms), ptOrder, hasFlow, allFactors});
  }

  struct NaturalMean {
    int region{};
    int ptIndex{};
    bool pTDif{};
    std::string profile;
  };
  std::vector<NaturalMean> naturalMeans;
  auto getPtRegion = [](const ParsedPtRequest& request, int ptIndex) {
    const auto group = request.ptTerms[ptIndex].first;
    return request.config.Regs[group].front();
  };
  auto findNaturalMean = [&](const ParsedPtRequest& request, int ptIndex) {
    const auto group = request.ptTerms[ptIndex].first;
    const int region = getPtRegion(request, ptIndex);
    return std::find_if(naturalMeans.begin(), naturalMeans.end(), [&](const auto& mean) {
      return mean.region == region && mean.ptIndex == request.config.ptInd[group] && mean.pTDif == request.config.pTDif;
    });
  };

  // Prefer explicitly requested one-particle observables as the shared native means,
  // irrespective of their order in the configuration list.
  for (const auto& request : requests) {
    if (request.hasFlow || request.ptOrder != 1 || findNaturalMean(request, 0) != naturalMeans.end()) {
      continue;
    }
    const auto group = request.ptTerms.front().first;
    naturalMeans.push_back({getPtRegion(request, 0), request.config.ptInd[group], request.config.pTDif, request.config.Head});
  }

  std::vector<PtCorrConfig> parsedConfigs;
  parsedConfigs.reserve(requests.size());
  std::vector<PtMomentConfig> momentConfigs;
  momentConfigs.reserve(requests.size());
  for (const auto& request : requests) {
    const int ptOffset = request.hasFlow ? 1 : 0;
    PtCorrConfig parsed{request.config.Head, request.allFactors, request.hasFlow, {}, std::vector<std::string>(request.ptOrder)};

    // Factors in one single-region pT group are exchangeable. Their full-topology
    // terms therefore depend only on how many factors are active.
    std::vector<std::vector<int>> factorClasses;
    if (request.hasFlow) {
      factorClasses.push_back({0});
    }
    for (size_t group = 0; group < request.isPt.size(); ++group) {
      if (!request.isPt[group].front()) {
        continue;
      }
      factorClasses.push_back(request.factors[group]);
    }

    std::vector<std::uint64_t> activeMasks{0};
    for (const auto& factors : factorClasses) {
      const size_t previousSize = activeMasks.size();
      for (size_t i = 0; i < previousSize; ++i) {
        std::uint64_t addedFactors = 0;
        for (const int& factor : factors) {
          addedFactors |= std::uint64_t{1} << factor;
          activeMasks.push_back(activeMasks[i] | addedFactors);
        }
      }
    }
    std::sort(activeMasks.begin(), activeMasks.end());
    parsed.profiles.reserve(activeMasks.size() + request.ptOrder);

    // Inactive factors stay in the tuple as harmonic-zero multiplicity insertions,
    // so every subset retains the event weight of the full correlation.
    auto addFullTopologyProfile = [&](std::uint64_t activeFactors) {
      auto profileConfig = request.config;
      if (request.hasFlow && !(activeFactors & std::uint64_t{1})) {
        for (size_t group = 0; group < request.isPt.size(); ++group) {
          if (!request.isPt[group].front()) {
            std::fill(profileConfig.Hars[group].begin(), profileConfig.Hars[group].end(), 0);
          }
        }
      }
      profileConfig.Head = makeFullTopologyName(request.config.Head, activeFactors, request.allFactors, request.hasFlow, request.ptOrder);
      parsed.profiles.push_back({std::move(profileConfig), request.isPt, request.factors, activeFactors, true});
    };

    addFullTopologyProfile(request.allFactors);
    for (const auto& activeFactors : activeMasks) {
      if (activeFactors != 0 && activeFactors != request.allFactors) {
        addFullTopologyProfile(activeFactors);
      }
    }

    // Each distinct pT source has one naturally weighted mean. Explicit one-particle
    // requests provide that mean directly; otherwise the first consumer books it.
    for (int ptIndex = 0; ptIndex < request.ptOrder; ++ptIndex) {
      const std::uint64_t factor = std::uint64_t{1} << (ptOffset + ptIndex);
      if (factor == request.allFactors) {
        parsed.meanProfiles[ptIndex] = request.config.Head;
        continue;
      }
      auto existingMean = findNaturalMean(request, ptIndex);
      if (existingMean != naturalMeans.end()) {
        parsed.meanProfiles[ptIndex] = existingMean->profile;
        continue;
      }

      const auto group = request.ptTerms[ptIndex].first;
      GFW::CorrConfig meanConfig;
      meanConfig.Regs = {{getPtRegion(request, ptIndex)}};
      meanConfig.Hars = {{0}};
      meanConfig.Overlap = {-1};
      meanConfig.ptInd = {request.config.ptInd[group]};
      meanConfig.pTDif = request.config.pTDif;
      meanConfig.Head = request.config.Head + "_mean_pt" + std::to_string(ptIndex + 1);
      parsed.meanProfiles[ptIndex] = meanConfig.Head;
      naturalMeans.push_back({meanConfig.Regs.front().front(), meanConfig.ptInd.front(), meanConfig.pTDif, meanConfig.Head});
      parsed.profiles.push_back({std::move(meanConfig), {{true}}, {{ptOffset + ptIndex}}, factor, false});
    }

    if (!request.hasFlow) {
      PtMomentConfig momentConfig{request.config.Head, {}};
      momentConfig.groups.reserve(request.isPt.size());
      for (size_t group = 0; group < request.isPt.size(); ++group) {
        const int firstFactor = request.factors[group].front();
        const int region = request.config.Regs[group].front();
        momentConfig.groups.push_back({regions.GetNames()[region], static_cast<int>(request.factors[group].size()), parsed.meanProfiles[firstFactor]});
      }
      momentConfigs.push_back(std::move(momentConfig));
    }

    parsedConfigs.push_back(std::move(parsed));
  }

  auto* profiles = new TList();
  profiles->SetOwner(kTRUE);
  std::unordered_set<std::string> profileNames;
  for (const auto& parsed : parsedConfigs) {
    for (const auto& profile : parsed.profiles) {
      const auto& head = profile.config.Head;
      if (!profileNames.insert(head).second) {
        LOGF(error, "Cannot initialize %s: generated profile name %s is not unique", GetName(), head.c_str());
        delete profiles;
        return;
      }
      if (isVariableBinning) {
        profiles->Add(new BootstrapProfile(head.c_str(), head.c_str(), numberOfBins, multiBins.data()));
      } else {
        profiles->Add(new BootstrapProfile(head.c_str(), head.c_str(), numberOfBins, multiBins.front(), multiBins.back()));
      }
    }
  }
  if (nsub > 0) {
    for (int i = 0; i < profiles->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(profiles->At(i))->InitializeSubsamples(nsub);
    }
  }

  delete fCorrList;
  fCorrList = profiles;
  fPtRegions = regions;
  fPtCorrConfigs = std::move(parsedConfigs);
  fPtMomentConfigs = std::move(momentConfigs);
  clearPtPostprocessLists();
  fPtRegionMoments.clear();
  fPtRegionMoments.reserve(regionMaxOrders.size());
  for (const int& order : regionMaxOrders) {
    fPtRegionMoments.push_back({order, std::vector<double>((order + 1) * (order + 1), 0.)});
  }
  mpar = maxPtOrder;
  nSubevents = static_cast<int>(numberOfRegions);
  LOGF(info, "Container %s initialized with %zu pT correlations, %i profiles and %i subsamples", GetName(), numberOfConfigs, profiles->GetEntries(), nsub);
}
void FlowPtContainer::initialise(const o2::framework::AxisSpec& axis, const int& maxOrder, const GFWCorrConfigs& configs, const int& nsub)
{
  arr.resize(3 * 3 * 3 * 3);
  warr.resize(3 * 3 * 3 * 3);
  if (mpar == 0) {
    mpar = maxOrder;
  }
  std::vector<double> multiBins = axis.binEdges;
  int nMultiBins = axis.nBins.value_or(0);
  if (nMultiBins <= 0) {
    nMultiBins = multiBins.size() - 1;
  }
  if (nMultiBins <= 0) {
    LOGF(warning, "Multiplicity axis does not exist");
    return;
  }
  delete fCMTermList;
  fCMTermList = new TList();
  fCMTermList->SetOwner(kTRUE);
  delete fCorrList;
  fCorrList = new TList();
  fCorrList->SetOwner(kTRUE);
  delete fCovList;
  fCovList = new TList();
  fCovList->SetOwner(kTRUE);
  for (int m = 0; m < mpar; ++m) {
    fCorrList->Add(new BootstrapProfile(Form("mpt%i", m + 1), Form("mpt%i", m + 1), nMultiBins, multiBins.data()));
  }
  for (int m = 0; m < centralMomentMaxOrder; ++m) {
    for (int i = 0; i <= m; ++i) {
      fCMTermList->Add(new BootstrapProfile(Form("cm%i_Mpt%i", m + 1, i), Form("cm%i_Mpt%i", m + 1, i), nMultiBins, multiBins.data()));
    }
  }
  if (fUseGap) {
    int obsIndex = 0;
    fCovFirstIndex.resize(configs.GetSize(), 0);
    for (int i = 0; i < configs.GetSize(); ++i) {
      fCovFirstIndex[i] = obsIndex;
      for (auto m(1); m <= mpar; ++m) {
        if ((configs.GetpTCorrMasks()[i] & (1 << (m - 1))) == 0) {
          continue;
        }
        if (fUseCentralMoments) {
          for (auto j = 0; j <= m; ++j) {
            fCovList->Add(new BootstrapProfile(Form("%spt%i_Mpt%i", configs.GetHeads()[i].c_str(), m, j), Form("%spt%i_Mpt%i", configs.GetHeads()[i].c_str(), m, j), nMultiBins, multiBins.data()));
            obsIndex++;
          }
        } else {
          fCovList->Add(new BootstrapProfile(Form("%spt%i", configs.GetHeads()[i].c_str(), m), Form("%spt%i", configs.GetHeads()[i].c_str(), m), nMultiBins, multiBins.data()));
          obsIndex++;
        }
      }
    }
  } else {
    if (fUseCentralMoments) {
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt0", "ChFull24pt2_Mpt0", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt1", "ChFull24pt2_Mpt1", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt2", "ChFull24pt2_Mpt2", nMultiBins, multiBins.data()));

      fCovList->Add(new BootstrapProfile("ChFull24pt1_Mpt0", "ChFull24pt1_Mpt0", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull24pt1_Mpt1", "ChFull24pt1_Mpt1", nMultiBins, multiBins.data()));

      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt0", "ChFull22pt2_Mpt0", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt1", "ChFull22pt2_Mpt1", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt2", "ChFull22pt2_Mpt2", nMultiBins, multiBins.data()));

      fCovList->Add(new BootstrapProfile("ChFull22pt1_Mpt0", "ChFull22pt1_Mpt0", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull22pt1_Mpt1", "ChFull22pt1_Mpt1", nMultiBins, multiBins.data()));

      fCovList->Add(new BootstrapProfile("ChFull22pt3_Mpt0", "ChFull22pt3_Mpt0", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull22pt3_Mpt1", "ChFull22pt3_Mpt1", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull22pt3_Mpt2", "ChFull22pt3_Mpt2", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull22pt3_Mpt3", "ChFull22pt3_Mpt3", nMultiBins, multiBins.data()));

      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt0", "ChFull22pt4_Mpt0", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt1", "ChFull22pt4_Mpt1", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt2", "ChFull22pt4_Mpt2", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt3", "ChFull22pt4_Mpt3", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt4", "ChFull22pt4_Mpt4", nMultiBins, multiBins.data()));
    } else {
      fCovList->Add(new BootstrapProfile("ChFull24pt2", "ChFull24pt2", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull24pt1", "ChFull24pt1", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull22pt2", "ChFull22pt2", nMultiBins, multiBins.data()));
      fCovList->Add(new BootstrapProfile("ChFull22pt1", "ChFull22pt1", nMultiBins, multiBins.data()));
    }
  }

  if (nsub != 0) {
    for (int i = 0; i < fCorrList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fCorrList->At(i))->InitializeSubsamples(nsub);
    }
    for (int i = 0; i < fCMTermList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fCMTermList->At(i))->InitializeSubsamples(nsub);
    }
    for (int i = 0; i < fCovList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fCovList->At(i))->InitializeSubsamples(nsub);
    }
  }
  LOGF(info, "Container %s initialized with m = %i\n and %i subsamples", this->GetName(), mpar, nsub);
};
void FlowPtContainer::initialise(int nbinsx, double* xbins, const int& maxOrder, const GFWCorrConfigs& configs, const int& nsub)
{
  arr.resize(3 * 3 * 5 * 5);
  warr.resize(3 * 3 * 5 * 5);
  if (mpar == 0) {
    mpar = maxOrder;
  }
  delete fCMTermList;
  fCMTermList = new TList();
  fCMTermList->SetOwner(kTRUE);
  delete fCorrList;
  fCorrList = new TList();
  fCorrList->SetOwner(kTRUE);
  delete fCovList;
  fCovList = new TList();
  fCovList->SetOwner(kTRUE);
  for (int m = 0; m < mpar; ++m) {
    fCorrList->Add(new BootstrapProfile(Form("mpt%i", m + 1), Form("mpt%i", m + 1), nbinsx, xbins));
  }
  for (int m = 0; m < centralMomentMaxOrder; ++m) {
    for (int i = 0; i <= m; ++i) {
      fCMTermList->Add(new BootstrapProfile(Form("cm%i_Mpt%i", m + 1, i), Form("cm%i_Mpt%i", m + 1, i), nbinsx, xbins));
    }
  }
  if (fUseGap) {
    int obsIndex = 0;
    fCovFirstIndex.resize(configs.GetSize(), 0);
    for (int i = 0; i < configs.GetSize(); ++i) {
      fCovFirstIndex[i] = obsIndex;
      for (auto m(1); m <= mpar; ++m) {
        if ((configs.GetpTCorrMasks()[i] & (1 << (m - 1))) == 0) {
          continue;
        }
        if (fUseCentralMoments) {
          for (auto j = 0; j <= m; ++j) {
            fCovList->Add(new BootstrapProfile(Form("%spt%i_Mpt%i", configs.GetHeads()[i].c_str(), m, j), Form("%spt%i_Mpt%i", configs.GetHeads()[i].c_str(), m, j), nbinsx, xbins));
            obsIndex++;
          }
        } else {
          fCovList->Add(new BootstrapProfile(Form("%spt%i", configs.GetHeads()[i].c_str(), m), Form("%spt%i", configs.GetHeads()[i].c_str(), m), nbinsx, xbins));
          obsIndex++;
        }
      }
    }
  } else {
    if (fUseCentralMoments) {
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt0", "ChFull24pt2_Mpt0", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt1", "ChFull24pt2_Mpt1", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt2", "ChFull24pt2_Mpt2", nbinsx, xbins));

      fCovList->Add(new BootstrapProfile("ChFull24pt1_Mpt0", "ChFull24pt1_Mpt0", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull24pt1_Mpt1", "ChFull24pt1_Mpt1", nbinsx, xbins));

      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt0", "ChFull22pt2_Mpt0", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt1", "ChFull22pt2_Mpt1", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt2", "ChFull22pt2_Mpt2", nbinsx, xbins));

      fCovList->Add(new BootstrapProfile("ChFull22pt1_Mpt0", "ChFull22pt1_Mpt0", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt1_Mpt1", "ChFull22pt1_Mpt1", nbinsx, xbins));

      fCovList->Add(new BootstrapProfile("ChFull22pt3_Mpt0", "ChFull22pt3_Mpt0", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt3_Mpt1", "ChFull22pt3_Mpt1", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt3_Mpt2", "ChFull22pt3_Mpt2", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt3_Mpt3", "ChFull22pt3_Mpt3", nbinsx, xbins));

      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt0", "ChFull22pt4_Mpt0", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt1", "ChFull22pt4_Mpt1", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt2", "ChFull22pt4_Mpt2", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt3", "ChFull22pt4_Mpt3", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt4", "ChFull22pt4_Mpt4", nbinsx, xbins));
    } else {
      fCovList->Add(new BootstrapProfile("ChFull24pt2", "ChFull24pt2", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull24pt1", "ChFull24pt1", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt2", "ChFull22pt2", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt1", "ChFull22pt1", nbinsx, xbins));
    }
  }
  if (nsub != 0) {
    for (int i = 0; i < fCorrList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fCorrList->At(i))->InitializeSubsamples(nsub);
    }
    for (int i = 0; i < fCMTermList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fCMTermList->At(i))->InitializeSubsamples(nsub);
    }
    for (int i = 0; i < fCovList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fCovList->At(i))->InitializeSubsamples(nsub);
    }
  }
  LOGF(info, "Container %s initialized with m = %i\n", this->GetName(), mpar);
};
void FlowPtContainer::initialise(int nbinsx, double xlow, double xhigh, const int& maxOrder, const GFWCorrConfigs& configs, const int& nsub)
{
  arr.resize(3 * 3 * 5 * 5);
  warr.resize(3 * 3 * 5 * 5);
  if (mpar == 0) {
    mpar = maxOrder;
  }
  delete fCMTermList;
  fCMTermList = new TList();
  fCMTermList->SetOwner(kTRUE);
  delete fCorrList;
  fCorrList = new TList();
  fCorrList->SetOwner(kTRUE);
  delete fCovList;
  fCovList = new TList();
  fCovList->SetOwner(kTRUE);
  for (int m = 0; m < mpar; ++m) {
    fCorrList->Add(new BootstrapProfile(Form("mpt%i", m + 1), Form("mpt%i", m + 1), nbinsx, xlow, xhigh));
  }
  for (int m = 0; m < centralMomentMaxOrder; ++m) {
    for (int i = 0; i <= m; ++i) {
      fCMTermList->Add(new BootstrapProfile(Form("cm%i_Mpt%i", m + 1, i), Form("cm%i_Mpt%i", m + 1, i), nbinsx, xlow, xhigh));
    }
  }
  if (fUseGap) {
    int obsIndex = 0;
    fCovFirstIndex.resize(configs.GetSize(), 0);
    for (int i = 0; i < configs.GetSize(); ++i) {
      fCovFirstIndex[i] = obsIndex;
      for (auto m(1); m <= mpar; ++m) {
        if ((configs.GetpTCorrMasks()[i] & (1 << (m - 1))) == 0) {
          continue;
        }
        if (fUseCentralMoments) {
          for (auto j = 0; j <= m; ++j) {
            fCovList->Add(new BootstrapProfile(Form("%spt%i_Mpt%i", configs.GetHeads()[i].c_str(), m, j), Form("%spt%i_Mpt%i", configs.GetHeads()[i].c_str(), m, j), nbinsx, xlow, xhigh));
            obsIndex++;
          }
        } else {
          fCovList->Add(new BootstrapProfile(Form("%spt%i", configs.GetHeads()[i].c_str(), m), Form("%spt%i", configs.GetHeads()[i].c_str(), m), nbinsx, xlow, xhigh));
          obsIndex++;
        }
      }
    }
  } else {
    if (fUseCentralMoments) {
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt0", "ChFull24pt2_Mpt0", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt1", "ChFull24pt2_Mpt1", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt2", "ChFull24pt2_Mpt2", nbinsx, xlow, xhigh));

      fCovList->Add(new BootstrapProfile("ChFull24pt1_Mpt0", "ChFull24pt1_Mpt0", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull24pt1_Mpt1", "ChFull24pt1_Mpt1", nbinsx, xlow, xhigh));

      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt0", "ChFull22pt2_Mpt0", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt1", "ChFull22pt2_Mpt1", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt2", "ChFull22pt2_Mpt2", nbinsx, xlow, xhigh));

      fCovList->Add(new BootstrapProfile("ChFull22pt1_Mpt0", "ChFull22pt1_Mpt0", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt1_Mpt1", "ChFull22pt1_Mpt1", nbinsx, xlow, xhigh));

      fCovList->Add(new BootstrapProfile("ChFull22pt3_Mpt0", "ChFull22pt3_Mpt0", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt3_Mpt1", "ChFull22pt3_Mpt1", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt3_Mpt2", "ChFull22pt3_Mpt2", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt3_Mpt3", "ChFull22pt3_Mpt3", nbinsx, xlow, xhigh));

      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt0", "ChFull22pt4_Mpt0", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt1", "ChFull22pt4_Mpt1", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt2", "ChFull22pt4_Mpt2", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt3", "ChFull22pt4_Mpt3", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt4_Mpt4", "ChFull22pt4_Mpt4", nbinsx, xlow, xhigh));
    } else {
      fCovList->Add(new BootstrapProfile("ChFull24pt2", "ChFull24pt2", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull24pt1", "ChFull24pt1", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt2", "ChFull22pt2", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt1", "ChFull22pt1", nbinsx, xlow, xhigh));
    }
  }
  if (nsub != 0) {
    for (int i = 0; i < fCorrList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fCorrList->At(i))->InitializeSubsamples(nsub);
    }
    for (int i = 0; i < fCMTermList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fCMTermList->At(i))->InitializeSubsamples(nsub);
    }
    for (int i = 0; i < fCovList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fCovList->At(i))->InitializeSubsamples(nsub);
    }
  }
  LOGF(info, "Container %s initialized with m = %i\n", this->GetName(), mpar);
};
void FlowPtContainer::initialiseSubevent(const o2::framework::AxisSpec& axis, const int& maxOrder, const int& nsubev, const int& nsub)
{
  if (nsubev < 1) {
    LOGF(fatal, "Need at least one subevent");
    return;
  }
  nSubevents = nsubev;
  if (mpar == 0) {
    mpar = maxOrder;
  }
  std::vector<double> multiBins = axis.binEdges;
  int nMultiBins = axis.nBins.value_or(0);
  if (nMultiBins <= 0) {
    nMultiBins = multiBins.size() - 1;
  }
  if (nMultiBins <= 0) {
    LOGF(warning, "Multiplicity axis does not exist");
    return;
  }
  delete fSubList;
  fSubList = new TList();
  fSubList->SetOwner(kTRUE);

  // Get all possible subevent combinations given m particles and nsubev subevents - also considering not using all m particles, e.g. all lower orders
  std::vector<int> current;
  getSubevents(mpar, nsubev + 1, current, subevents);
  // remove unused "extra" subevent
  for (auto& subevent : subevents) { // o2-linter: disable=const-ref-in-for-loop (modified through pop_back())
    subevent.pop_back();
  }
  subevents.erase(subevents.begin(), subevents.begin() + 1);

  std::vector<std::string> histnames;
  for (const auto& subevent : subevents) {
    std::string outstr = "ptptsub";
    int index = 0;
    for (const auto& particles : subevent) {
      outstr += "_" + std::to_string(particles) +
                "sub" + std::to_string(index + 1);
      ++index;
    }
    histnames.push_back(outstr);
  }
  for (const auto& name : histnames) {
    fSubList->Add(new BootstrapProfile(name.c_str(), this->GetTitle(), nMultiBins, multiBins.data()));
  }
  delete fSubCMList;
  fSubCMList = new TList();
  fSubCMList->SetOwner(kTRUE);
  const int maxSubEv = 2;
  for (int subEv = 0; subEv < maxSubEv; ++subEv) {
    for (int m = 0; m < centralMomentMaxOrder; ++m) {
      for (int i = 0; i <= m; ++i) {
        fSubCMList->Add(new BootstrapProfile(Form("cm%i_sub%i_Mpt%i", m + 1, subEv + 1, i), this->GetTitle(), nMultiBins, multiBins.data()));
      }
    }
  }
  for (int m = 2; m <= centralMomentMaxOrder; ++m) {
    for (int first = 1; first < m; ++first) {
      for (int second = first; second < m; ++second) {
        if (first > second) {
          continue;
        }
        int fourth = m - second;
        for (int third = 1; third < m; ++third) {
          if (third > fourth) {
            continue;
          }
          fSubCMList->Add(new BootstrapProfile(Form("cm%i_%i%isub1_%i%isub2", m, first, second, third, fourth), this->GetTitle(), nMultiBins, multiBins.data()));
        }
      }
    }
  }

  if (nsub != 0) {
    for (int i = 0; i < fSubList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fSubList->At(i))->InitializeSubsamples(nsub);
    }
    for (int i = 0; i < fSubCMList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(i))->InitializeSubsamples(nsub);
    }
  }
  LOGF(info, "Container %s initialized Subevents and %i subsamples", this->GetName(), nsub);
}
void FlowPtContainer::initialiseSubevent(int nbinsx, double* xbins, const int& maxOrder, const int& nsubev, const int& nsub)
{
  if (nsubev < 1) {
    LOGF(fatal, "Need at least one subevent");
    return;
  }
  nSubevents = nsubev;
  if (mpar == 0) {
    mpar = maxOrder;
  }
  delete fSubList;
  fSubList = new TList();
  fSubList->SetOwner(kTRUE);

  // Get all possible subevent combinations given m particles and nsubev subevents - also considering not using all m particles, e.g. all lower orders
  std::vector<int> current;
  getSubevents(mpar, nsubev + 1, current, subevents);
  // remove unused "extra" subevent
  for (auto& subevent : subevents) { // o2-linter: disable=const-ref-in-for-loop (modified through pop_back())
    subevent.pop_back();
  }
  subevents.erase(subevents.begin(), subevents.begin() + 1);

  std::vector<std::string> histnames;
  for (const auto& subevent : subevents) {
    std::string outstr = "ptptsub";
    int index = 0;
    for (const auto& particles : subevent) {
      outstr += "_" + std::to_string(particles) +
                "sub" + std::to_string(index + 1);
      ++index;
    }
    histnames.push_back(outstr);
  }
  for (const auto& name : histnames) {
    fSubList->Add(new BootstrapProfile(name.c_str(), this->GetTitle(), nbinsx, xbins));
  }
  delete fSubCMList;
  fSubCMList = new TList();
  fSubCMList->SetOwner(kTRUE);
  const int maxSubEv = 2;
  for (int subEv = 0; subEv < maxSubEv; ++subEv) {
    for (int m = 0; m < centralMomentMaxOrder; ++m) {
      for (int i = 0; i <= m; ++i) {
        fSubCMList->Add(new BootstrapProfile(Form("cm%i_sub%i_Mpt%i", m + 1, subEv + 1, i), this->GetTitle(), nbinsx, xbins));
      }
    }
  }
  for (int m = 2; m <= centralMomentMaxOrder; ++m) {
    for (int first = 1; first < m; ++first) {
      for (int second = first; second < m; ++second) {
        if (first > second) {
          continue;
        }
        int fourth = m - second;
        for (int third = 1; third < m; ++third) {
          if (third > fourth) {
            continue;
          }
          fSubCMList->Add(new BootstrapProfile(Form("cm%i_%i%isub1_%i%isub2", m, first, second, third, fourth), this->GetTitle(), nbinsx, xbins));
        }
      }
    }
  }

  if (nsub != 0) {
    for (int i = 0; i < fSubList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fSubList->At(i))->InitializeSubsamples(nsub);
    }
    for (int i = 0; i < fSubCMList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(i))->InitializeSubsamples(nsub);
    }
  }
  LOGF(info, "Container %s initialized Subevents and %i subsamples", this->GetName(), nsub);
}
void FlowPtContainer::initialiseSubevent(int nbinsx, double xlow, double xhigh, const int& maxOrder, const int& nsubev, const int& nsub)
{
  if (nsubev < 1) {
    LOGF(fatal, "Need at least one subevent");
    return;
  }
  nSubevents = nsubev;
  if (mpar == 0) {
    mpar = maxOrder;
  }
  delete fSubList;
  fSubList = new TList();
  fSubList->SetOwner(kTRUE);

  // Get all possible subevent combinations given m particles and nsubev subevents - also considering not using all m particles, e.g. all lower orders
  std::vector<int> current;
  getSubevents(mpar, nsubev + 1, current, subevents);
  // remove unused "extra" subevent
  for (auto& subevent : subevents) { // o2-linter: disable=const-ref-in-for-loop (modified through pop_back())
    subevent.pop_back();
  }
  subevents.erase(subevents.begin(), subevents.begin() + 1);

  std::vector<std::string> histnames;
  for (const auto& subevent : subevents) {
    std::string outstr = "ptptsub";
    int index = 0;
    for (const auto& particles : subevent) {
      outstr += "_" + std::to_string(particles) +
                "sub" + std::to_string(index + 1);
      ++index;
    }
    histnames.push_back(outstr);
  }
  for (const auto& name : histnames) {
    fSubList->Add(new BootstrapProfile(name.c_str(), this->GetTitle(), nbinsx, xlow, xhigh));
  }
  delete fSubCMList;
  fSubCMList = new TList();
  fSubCMList->SetOwner(kTRUE);
  const int maxSubEv = 2;
  for (int subEv = 0; subEv < maxSubEv; ++subEv) {
    for (int m = 0; m < centralMomentMaxOrder; ++m) {
      for (int i = 0; i <= m; ++i) {
        fSubCMList->Add(new BootstrapProfile(Form("cm%i_sub%i_Mpt%i", m + 1, subEv + 1, i), this->GetTitle(), nbinsx, xlow, xhigh));
      }
    }
  }
  for (int m = 2; m <= centralMomentMaxOrder; ++m) {
    for (int first = 1; first < m; ++first) {
      for (int second = first; second < m; ++second) {
        if (first > second) {
          continue;
        }
        int fourth = m - second;
        for (int third = 1; third < m; ++third) {
          if (third > fourth) {
            continue;
          }
          fSubCMList->Add(new BootstrapProfile(Form("cm%i_%i%isub1_%i%isub2", m, first, second, third, fourth), this->GetTitle(), nbinsx, xlow, xhigh));
        }
      }
    }
  }
  if (nsub != 0) {
    for (int i = 0; i < fSubList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fSubList->At(i))->InitializeSubsamples(nsub);
    }
    for (int i = 0; i < fSubCMList->GetEntries(); ++i) {
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(i))->InitializeSubsamples(nsub);
    }
  }
  LOGF(info, "Container %s initialized Subevents and %i subsamples", this->GetName(), nsub);
}
void FlowPtContainer::fill(const double& w, const double& pt)
{
  for (size_t i = 0; i < sumP.size(); ++i) {
    sumP[i] += std::pow(w, i % (mpar + 1)) * std::pow(pt, i / (mpar + 1));
  }
}
void FlowPtContainer::fill(const double& eta, const double& w, const double& pt, int mask)
{
  int maxOrder = 0;
  for (const auto& moments : fPtRegionMoments) {
    maxOrder = std::max(maxOrder, moments.maxOrder);
  }
  if (maxOrder == 0) {
    return;
  }

  std::vector<double> weightPowers(maxOrder + 1, 1.);
  std::vector<double> ptPowers(maxOrder + 1, 1.);
  for (int power = 1; power <= maxOrder; ++power) {
    weightPowers[power] = weightPowers[power - 1] * w;
    ptPowers[power] = ptPowers[power - 1] * pt;
  }

  for (size_t region = 0; region < fPtRegionMoments.size(); ++region) {
    auto& moments = fPtRegionMoments[region];
    if (moments.maxOrder == 0 || !(fPtRegions.GetEtaMin()[region] < eta && eta < fPtRegions.GetEtaMax()[region]) ||
        !(fPtRegions.GetBitmasks()[region] & mask)) {
      continue;
    }
    for (int weightPower = 1; weightPower <= moments.maxOrder; ++weightPower) {
      for (int ptPower = 0; ptPower <= weightPower; ++ptPower) {
        moments.values[weightPower * (moments.maxOrder + 1) + ptPower] += weightPowers[weightPower] * ptPowers[ptPower];
      }
    }
  }
}
void FlowPtContainer::clearPtRegionMoments()
{
  for (auto& moments : fPtRegionMoments) { // o2-linter: disable=const-ref-in-for-loop (contents are reset)
    std::fill(moments.values.begin(), moments.values.end(), 0.);
  }
}
void FlowPtContainer::clearPtPostprocessLists()
{
  delete fPtCentralMomentList;
  fPtCentralMomentList = nullptr;
  delete fPtCumulantList;
  fPtCumulantList = nullptr;
}
double FlowPtContainer::getPtMoment(int region, int weightPower, int ptPower) const
{
  if (region < 0 || region >= static_cast<int>(fPtRegionMoments.size())) {
    return 0.;
  }
  const auto& moments = fPtRegionMoments[region];
  if (weightPower < 1 || weightPower > moments.maxOrder || ptPower < 0 || ptPower > weightPower) {
    return 0.;
  }
  return moments.values[weightPower * (moments.maxOrder + 1) + ptPower];
}
double FlowPtContainer::recursivePtCorrelation(int poi, int ref, int overlap, std::vector<int>& weightPowers, std::vector<int>& ptPowers) const
{
  if (weightPowers.empty() || weightPowers.size() != ptPowers.size()) {
    return 0.;
  }
  if (weightPowers.front() != 1 && overlap >= 0) {
    poi = overlap;
  }
  if (weightPowers.size() == 1) {
    return getPtMoment(poi, weightPowers.front(), ptPowers.front());
  }
  if (weightPowers.size() == TwoParticleOrder) {
    const double product = getPtMoment(poi, weightPowers[0], ptPowers[0]) * getPtMoment(ref, weightPowers[1], ptPowers[1]);
    const double selfCorrelation = overlap >= 0 ? getPtMoment(overlap, weightPowers[0] + weightPowers[1], ptPowers[0] + ptPowers[1]) : 0.;
    return product - selfCorrelation;
  }

  const int lastWeightPower = weightPowers.back();
  const int lastPtPower = ptPowers.back();
  weightPowers.pop_back();
  ptPowers.pop_back();
  double result = recursivePtCorrelation(poi, ref, overlap, weightPowers, ptPowers) * getPtMoment(ref, lastWeightPower, lastPtPower);
  for (int term = static_cast<int>(weightPowers.size()) - 1; term >= 0; --term) {
    weightPowers[term] += lastWeightPower;
    ptPowers[term] += lastPtPower;
    result -= recursivePtCorrelation(poi, ref, overlap, weightPowers, ptPowers);
    weightPowers[term] -= lastWeightPower;
    ptPowers[term] -= lastPtPower;
  }
  weightPowers.push_back(lastWeightPower);
  ptPowers.push_back(lastPtPower);
  return result;
}
void FlowPtContainer::fillSub(const double& w, const double& pt, int subIndex)
{
  for (size_t i = 0; i < insub[subIndex].size(); ++i) {
    insub[subIndex][i] += std::pow(w, i % (mpar + 1)) * std::pow(pt, i / (mpar + 1));
  }
}
void FlowPtContainer::calculateCorrelations()
{
  corrNum.clear();
  corrNum.resize(mpar + 1, 0);
  corrNum[0] = 1.0;
  corrDen.clear();
  corrDen.resize(mpar + 1, 0);
  corrDen[0] = 1.0;
  double sumNum = 0;
  double sumDenum = 0;
  std::vector<double> valNum;
  std::vector<double> valDenum;
  for (int m(1); m <= mpar; ++m) {
    for (int k(1); k <= m; ++k) {
      valNum.push_back(SignArray[k - 1] * corrNum[m - k] * (FactorialArray[m - 1] / FactorialArray[m - k]) * sumP[getVectorIndex(k, k)]);
      valDenum.push_back(SignArray[k - 1] * corrDen[m - k] * (FactorialArray[m - 1] / FactorialArray[m - k]) * sumP[getVectorIndex(k, 0)]);
    }
    sumNum = orderedAddition(valNum);
    sumDenum = orderedAddition(valDenum);
    valNum.clear();
    valDenum.clear();

    corrNum[m] = sumNum;
    corrDen[m] = sumDenum;
  }
}
void FlowPtContainer::calculateSubeventCorrelations()
{
  corrNumSub.clear();
  corrNumSub.resize(nSubevents, std::vector<double>(mpar + 1, 0));
  for (auto& corrnum : corrNumSub) { // o2-linter: disable=const-ref-in-for-loop (assigned a value)
    corrnum[0] = 1.0;
  }
  corrDenSub.resize(nSubevents, std::vector<double>(mpar + 1, 0));
  for (auto& corrden : corrDenSub) { // o2-linter: disable=const-ref-in-for-loop (assigned a value)
    corrden[0] = 1.0;
  }

  for (int subIndex = 0; subIndex < nSubevents; ++subIndex) {
    double sumNum = 0.0;
    double sumDenum = 0.0;
    std::vector<double> valNum;
    std::vector<double> valDenum;
    for (int m(1); m <= mpar; ++m) {
      for (int k(1); k <= m; ++k) {
        valNum.push_back(SignArray[k - 1] * corrNumSub[subIndex][m - k] * (FactorialArray[m - 1] / FactorialArray[m - k]) * insub[subIndex][getVectorIndex(k, k)]);
        valDenum.push_back(SignArray[k - 1] * corrDenSub[subIndex][m - k] * (FactorialArray[m - 1] / FactorialArray[m - k]) * insub[subIndex][getVectorIndex(k, 0)]);
      }
      sumNum = orderedAddition(valNum);
      sumDenum = orderedAddition(valDenum);
      valNum.clear();
      valDenum.clear();
      corrNumSub[subIndex][m] = sumNum;
      corrDenSub[subIndex][m] = sumDenum;
    }
  }
}
void FlowPtContainer::fillPtProfiles(const double& centmult, const double& rn)
{
  for (int m = 1; m <= mpar; ++m) {
    if (corrDen[m] != 0) {
      dynamic_cast<BootstrapProfile*>(fCorrList->At(m - 1))->FillProfile(centmult, corrNum[m] / corrDen[m], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : corrDen[m], rn);
    }
  }
}
void FlowPtContainer::fillPtProfiles(GFW& gfw, const double& centmult, const double& rn)
{
  if (!fCorrList) {
    return;
  }
  clearPtPostprocessLists();
  for (const auto& correlation : fPtCorrConfigs) {
    for (const auto& profileConfig : correlation.profiles) {
      std::complex<double> numerator{1., 0.};
      double denominator = 1.;
      for (size_t group = 0; group < profileConfig.config.Regs.size(); ++group) {
        if (profileConfig.isPt[group].front()) {
          const int poi = profileConfig.config.Regs[group].front();
          const int ref = profileConfig.config.Regs[group].size() > 1 ? profileConfig.config.Regs[group][1] : poi;
          int overlap = profileConfig.config.Overlap[group];
          if (overlap < 0 && poi == ref) {
            overlap = ref;
          }
          std::vector<int> weightPowers(profileConfig.isPt[group].size(), 1);
          std::vector<int> ptPowers;
          ptPowers.reserve(profileConfig.isPt[group].size());
          for (const int& factor : profileConfig.factors[group]) {
            ptPowers.push_back((profileConfig.activeFactors & (std::uint64_t{1} << factor)) ? 1 : 0);
          }
          auto denominatorPtPowers = std::vector<int>(ptPowers.size(), 0);
          denominator *= recursivePtCorrelation(poi, ref, overlap, weightPowers, denominatorPtPowers);
          numerator *= recursivePtCorrelation(poi, ref, overlap, weightPowers, ptPowers);
        } else {
          GFW::CorrConfig groupConfig;
          groupConfig.Regs = {profileConfig.config.Regs[group]};
          groupConfig.Hars = {profileConfig.config.Hars[group]};
          groupConfig.Overlap = {profileConfig.config.Overlap[group]};
          groupConfig.ptInd = {profileConfig.config.ptInd[group]};
          groupConfig.pTDif = profileConfig.config.pTDif;
          numerator *= gfw.Calculate(groupConfig, 0, false);
          denominator *= gfw.Calculate(groupConfig, 0, true).real();
        }
      }
      if (!std::isfinite(numerator.real()) || !std::isfinite(numerator.imag()) || !std::isfinite(denominator) || denominator <= 0.) {
        continue;
      }
      auto* profile = dynamic_cast<BootstrapProfile*>(fCorrList->FindObject(profileConfig.config.Head.c_str()));
      if (!profile) {
        LOGF(error, "pT profile %s has not been booked", profileConfig.config.Head.c_str());
        continue;
      }
      const double eventWeight = fEventWeight == EventWeight::UnityWeight ? 1. : denominator;
      profile->FillProfile(centmult, numerator.real() / denominator, eventWeight, rn);
    }
  }
}
bool FlowPtContainer::addPtProfile(const char* name, int observableOrder)
{
  if (!fCorrList || !name || (name[0] == 0) || observableOrder < 1 || observableOrder > mpar) {
    LOGF(error, "Cannot add pT profile %s for order %d", name ? name : "(null)", observableOrder);
    return false;
  }
  const std::string profileName{name};
  if (fCorrList->FindObject(profileName.c_str())) {
    LOGF(error, "pT profile %s already exists", profileName.c_str());
    return false;
  }
  auto* original = dynamic_cast<BootstrapProfile*>(fCorrList->At(observableOrder - 1));
  const auto* axis = original->GetXaxis();
  BootstrapProfile* profile = nullptr;
  if (axis->GetXbins()->GetSize() != 0) {
    profile = new BootstrapProfile(profileName.c_str(), profileName.c_str(), axis->GetNbins(), axis->GetXbins()->GetArray());
  } else {
    profile = new BootstrapProfile(profileName.c_str(), profileName.c_str(), axis->GetNbins(), axis->GetXmin(), axis->GetXmax());
  }
  if (original->fListOfEntries) {
    profile->InitializeSubsamples(original->fListOfEntries->GetEntries());
  }
  fCorrList->Add(profile);
  return true;
}
bool FlowPtContainer::fillPtProfile(const char* name, int observableOrder, double mult, double eventWeight, double rn)
{
  if (!fCorrList || !name || observableOrder < 1 || observableOrder > mpar ||
      static_cast<size_t>(observableOrder) >= corrDen.size() || corrDen[observableOrder] == 0. || eventWeight == 0.) {
    return false;
  }
  auto* profile = dynamic_cast<BootstrapProfile*>(fCorrList->FindObject(name));
  if (!profile) {
    LOGF(error, "pT profile %s has not been booked", name);
    return false;
  }
  profile->FillProfile(mult, corrNum[observableOrder] / corrDen[observableOrder], eventWeight, rn);
  return true;
}
void FlowPtContainer::fillSubeventPtProfiles(const double& centmult, const double& rn)
{
  int histCounter = 0;
  for (const auto& subevent : subevents) {
    double val = 1.0;
    double dn = 1.0;
    int subIndex = 0;
    bool valid = true;
    for (const auto& m : subevent) {
      if (corrDenSub[subIndex][m] == 0) {
        valid = false;
        break;
      }
      val *= corrNumSub[subIndex][m] / corrDenSub[subIndex][m];
      dn *= corrDenSub[subIndex][m];

      ++subIndex;
    }
    if (valid) {
      dynamic_cast<BootstrapProfile*>(fSubList->At(histCounter))->FillProfile(centmult, val, (fEventWeight == EventWeight::UnityWeight) ? 1.0 : dn, rn);
    }
    ++histCounter;
  }
}
void FlowPtContainer::fillVnPtCorrProfiles(const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask)
{
  if (mask == 0u) {
    return;
  }
  for (auto m(1); m <= mpar; ++m) {
    if ((mask & (1 << (m - 1))) == 0) {
      continue;
    }
    if (corrDen[m] != 0) {
      dynamic_cast<BootstrapProfile*>(fCovList->At(fillCounter))->FillProfile(centmult, flowval * corrNum[m] / corrDen[m], (fEventWeight == UnityWeight) ? 1.0 : flowtuples * corrDen[m], rn);
    }
    ++fillCounter;
  }
}
void FlowPtContainer::fillVnDeltaPtProfiles(const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask)
{
  if (mask == 0u) {
    return;
  }
  for (auto m(1); m <= mpar; ++m) {
    if ((mask & (1 << (m - 1))) == 0) {
      continue;
    }
    for (auto i = 0; i <= m; ++i) {
      if (cmDen[m] != 0) {
        dynamic_cast<BootstrapProfile*>(fCovList->At(fillCounter))->FillProfile(centmult, flowval * ((i == m) ? cmVal[0] : cmVal[m * (m - 1) / 2 + i + 1]), (fEventWeight == UnityWeight) ? 1.0 : flowtuples * cmDen[m], rn);
      }
      ++fillCounter;
    }
  }
}
void FlowPtContainer::fillVnPtCorrProfiles(const int configIndex, const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask)
{
  if (mask == 0u) {
    return;
  }
  int startIndex = fCovFirstIndex[configIndex];
  for (auto m(1); m <= mpar; ++m) {
    if ((mask & (1 << (m - 1))) == 0) {
      continue;
    }
    if (corrDen[m] != 0) {
      dynamic_cast<BootstrapProfile*>(fCovList->At(startIndex))->FillProfile(centmult, flowval * corrNum[m] / corrDen[m], (fEventWeight == UnityWeight) ? 1.0 : flowtuples * corrDen[m], rn);
    }
    ++startIndex;
  }
}
void FlowPtContainer::fillVnDeltaPtProfiles(const int configIndex, const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask)
{
  if (mask == 0u) {
    return;
  }

  if (configIndex < 0 || static_cast<size_t>(configIndex) >= fCovFirstIndex.size()) {
    LOGF(error, "Invalid configuration index %d", configIndex);
    return;
  }
  int profileIndex = fCovFirstIndex[configIndex];
  for (auto m(1); m <= mpar; ++m) {
    if ((mask & (1 << (m - 1))) == 0) {
      continue;
    }
    for (auto i = 0; i <= m; ++i, ++profileIndex) {
      const size_t cmIndex = (i == m) ? 0u : static_cast<size_t>(m * (m - 1) / 2 + i + 1);

      if (static_cast<size_t>(m) >= cmDen.size() || cmIndex >= cmVal.size() || cmDen[m] == 0.) {
        continue;
      }

      if (!fCovList || profileIndex >= fCovList->GetEntries()) {
        LOGF(error, "Profile index %d out of range (entries=%d, config=%d, m=%d)", profileIndex, fCovList ? fCovList->GetEntries() : 0, configIndex, m);
        continue;
      }

      auto* profile = dynamic_cast<BootstrapProfile*>(fCovList->At(profileIndex));

      if (!profile) {
        LOGF(error, "Missing BootstrapProfile at index %d for config %d", profileIndex, configIndex);
        continue;
      }
      profile->FillProfile(centmult, flowval * cmVal[cmIndex], (fEventWeight == UnityWeight) ? 1.0 : flowtuples * cmDen[m], rn);
    }
  }
}
void FlowPtContainer::fillVnPtCorrStdProfiles(const double& centmult, const double& rn)
{
  double wAABBCC = getStdAABBCC(warr);
  if (wAABBCC != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(0))->FillProfile(centmult, getStdAABBCC(arr) / wAABBCC, (fEventWeight == UnityWeight) ? 1.0 : wAABBCC, rn);
  }
  double wAABBC = getStdAABBC(warr);
  if (wAABBC != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(1))->FillProfile(centmult, getStdAABBCC(arr) / wAABBC, (fEventWeight == UnityWeight) ? 1.0 : wAABBC, rn);
  }
  double wABCC = getStdAABBC(warr);
  if (wABCC != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(2))->FillProfile(centmult, getStdABCC(arr) / wABCC, (fEventWeight == UnityWeight) ? 1.0 : wABCC, rn);
  }
  double wABC = getStdABC(warr);
  if (wABC != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(3))->FillProfile(centmult, getStdABC(arr) / wABC, (fEventWeight == UnityWeight) ? 1.0 : wABC, rn);
  }
}
void FlowPtContainer::fillVnDeltaPtStdProfiles(const double& centmult, const double& rn)
{
  double wAABBCC = getStdAABBCC(warr);
  if (wAABBCC != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(0))->FillProfile(centmult, getStdAABBCC(arr) / wAABBCC, (fEventWeight == UnityWeight) ? 1.0 : wAABBCC, rn);
  }
  double wAABBCD = getStdAABBCD(warr);
  if (wAABBCD != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(1))->FillProfile(centmult, getStdAABBCD(arr) / wAABBCD, (fEventWeight == UnityWeight) ? 1.0 : wAABBCD, rn);
  }
  double wAABBDD = getStdAABBDD(warr);
  if (wAABBDD != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(2))->FillProfile(centmult, getStdAABBDD(arr) / wAABBDD, (fEventWeight == UnityWeight) ? 1.0 : wAABBDD, rn);
  }

  double wAABBC = getStdAABBC(warr);
  if (wAABBC != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(3))->FillProfile(centmult, getStdAABBC(arr) / wAABBC, (fEventWeight == UnityWeight) ? 1.0 : wAABBC, rn);
  }
  double wAABBD = getStdAABBD(warr);
  if (wAABBD != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(4))->FillProfile(centmult, getStdAABBD(arr) / wAABBD, (fEventWeight == UnityWeight) ? 1.0 : wAABBD, rn);
  }

  double wABCC = getStdABCC(warr);
  if (wABCC != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(5))->FillProfile(centmult, getStdABCC(arr) / wABCC, (fEventWeight == UnityWeight) ? 1.0 : wABCC, rn);
  }
  double wABCD = getStdABCD(warr);
  if (wABCD != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(6))->FillProfile(centmult, getStdABCD(arr) / wABCD, (fEventWeight == UnityWeight) ? 1.0 : wABCD, rn);
  }
  double wABDD = getStdABDD(warr);
  if (wABDD != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(7))->FillProfile(centmult, getStdABDD(arr) / wABDD, (fEventWeight == UnityWeight) ? 1.0 : wABDD, rn);
  }

  double wABC = getStdABC(warr);
  if (wABC != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(8))->FillProfile(centmult, getStdABC(arr) / wABC, (fEventWeight == UnityWeight) ? 1.0 : wABC, rn);
  }
  double wABD = getStdABD(warr);
  if (wABD != 0) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(9))->FillProfile(centmult, getStdABD(arr) / wABD, (fEventWeight == UnityWeight) ? 1.0 : wABD, rn);
  }
  double wABCCCC = getStdABCCCC(warr);
  if (wABCCCC != 0.) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(14))->FillProfile(centmult, getStdABCCCC(arr) / wABCCCC, (fEventWeight == UnityWeight) ? 1. : wABCCCC, rn);
  }
  double wABCCCD = getStdABCCCD(warr);
  if (wABCCCD != 0.) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(15))->FillProfile(centmult, getStdABCCCD(arr) / wABCCCD, (fEventWeight == UnityWeight) ? 1. : wABCCCD, rn);
  }
  double wABCCDD = getStdABCCDD(warr);
  if (wABCCDD != 0.) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(16))->FillProfile(centmult, getStdABCCDD(arr) / wABCCDD, (fEventWeight == UnityWeight) ? 1. : wABCCDD, rn);
  }
  double wABCDDD = getStdABCDDD(warr);
  if (wABCDDD != 0.) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(17))->FillProfile(centmult, getStdABCDDD(arr) / wABCDDD, (fEventWeight == UnityWeight) ? 1. : wABCDDD, rn);
  }
  double wABDDDD = getStdABDDDD(warr);
  if (wABDDDD != 0.) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(18))->FillProfile(centmult, getStdABDDDD(arr) / wABDDDD, (fEventWeight == UnityWeight) ? 1. : wABDDDD, rn);
  }
  double wABCCC = getStdABCCC(warr);
  if (wABCCC != 0.) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(10))->FillProfile(centmult, getStdABCCC(arr) / wABCCC, (fEventWeight == UnityWeight) ? 1. : wABCCC, rn);
  }
  double wABCCD = getStdABCCD(warr);
  if (wABCCD != 0.) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(11))->FillProfile(centmult, getStdABCCD(arr) / wABCCD, (fEventWeight == UnityWeight) ? 1. : wABCCD, rn);
  }
  double wABCDD = getStdABCDD(warr);
  if (wABCDD != 0.) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(12))->FillProfile(centmult, getStdABCDD(arr) / wABCDD, (fEventWeight == UnityWeight) ? 1. : wABCDD, rn);
  }
  double wABDDD = getStdABDDD(warr);
  if (wABDDD != 0.) {
    dynamic_cast<BootstrapProfile*>(fCovList->At(13))->FillProfile(centmult, getStdABDDD(arr) / wABDDD, (fEventWeight == UnityWeight) ? 1. : wABDDD, rn);
  }
}
void FlowPtContainer::fillCMProfiles(const double& centmult, const double& rn)
{
  if (sumP[getVectorIndex(0, 0)] == 0) {
    return;
  }
  // 0th order correlation
  cmDen.push_back(1.);
  cmVal.push_back(1.);

  cmDen.push_back(sumP[getVectorIndex(1, 0)]);
  cmDen.push_back(sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 0)] - sumP[getVectorIndex(2, 0)]);
  cmDen.push_back(sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 0)] - 3 * sumP[getVectorIndex(2, 0)] * sumP[getVectorIndex(1, 0)] + 2 * sumP[getVectorIndex(3, 0)]);
  cmDen.push_back(sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 0)] - 6 * sumP[getVectorIndex(2, 0)] * sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 0)] + 8 * sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(3, 0)] + 3 * sumP[getVectorIndex(2, 0)] * sumP[getVectorIndex(2, 0)] - 6 * sumP[getVectorIndex(4, 0)]);
  if (mpar < 1 || cmDen[1] == 0) {
    return;
  }
  cmVal.push_back(sumP[getVectorIndex(1, 1)] / cmDen[1]);
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(0))->FillProfile(centmult, cmVal[1], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDen[1], rn);
  if (mpar < 2 || sumP[getVectorIndex(2, 0)] == 0 || cmDen[2] == 0) { // o2-linter: disable=magic-number (less than order 2)
    return;
  }
  cmVal.push_back(1 / cmDen[2] * (sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 1)] - sumP[getVectorIndex(2, 2)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(1))->FillProfile(centmult, cmVal[2], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDen[2], rn);
  cmVal.push_back(-2 * 1 / cmDen[2] * (sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 1)] - sumP[getVectorIndex(2, 1)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(2))->FillProfile(centmult, cmVal[3], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDen[2], rn);
  if (mpar < 3 || sumP[getVectorIndex(3, 0)] == 0 || cmDen[3] == 0) { // o2-linter: disable=magic-number (less than order 3)
    return;
  }
  cmVal.push_back(1 / cmDen[3] * (sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 1)] - 3 * sumP[getVectorIndex(2, 2)] * sumP[getVectorIndex(1, 1)] + 2 * sumP[getVectorIndex(3, 3)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(3))->FillProfile(centmult, cmVal[4], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDen[3], rn);
  cmVal.push_back(-3 * 1 / cmDen[3] * (sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 0)] - 2 * sumP[getVectorIndex(2, 1)] * sumP[getVectorIndex(1, 1)] + 2 * sumP[getVectorIndex(3, 2)] - sumP[getVectorIndex(2, 2)] * sumP[getVectorIndex(1, 0)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(4))->FillProfile(centmult, cmVal[5], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDen[3], rn);
  cmVal.push_back(3 * 1 / cmDen[3] * (sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 0)] - 2 * sumP[getVectorIndex(2, 1)] * sumP[getVectorIndex(1, 0)] + 2 * sumP[getVectorIndex(3, 1)] - sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(2, 0)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(5))->FillProfile(centmult, cmVal[6], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDen[3], rn);
  if (mpar < 4 || sumP[getVectorIndex(4, 0)] == 0 || cmDen[4] == 0) { // o2-linter: disable=magic-number (less than order 4)
    return;
  }
  cmVal.push_back(1 / cmDen[4] * (sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 1)] - 6 * sumP[getVectorIndex(2, 2)] * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 1)] + 3 * sumP[getVectorIndex(2, 2)] * sumP[getVectorIndex(2, 2)] + 8 * sumP[getVectorIndex(3, 3)] * sumP[getVectorIndex(1, 1)] - 6 * sumP[getVectorIndex(4, 4)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(6))->FillProfile(centmult, cmVal[7], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDen[4], rn);
  cmVal.push_back(-4 * 1 / cmDen[4] * (sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 0)] - 3 * sumP[getVectorIndex(2, 2)] * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 0)] - 3 * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(2, 1)] + 3 * sumP[getVectorIndex(2, 2)] * sumP[getVectorIndex(2, 1)] + 2 * sumP[getVectorIndex(3, 3)] * sumP[getVectorIndex(1, 0)] + 6 * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(3, 2)] - 6 * sumP[getVectorIndex(4, 3)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(7))->FillProfile(centmult, cmVal[8], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDen[4], rn);
  cmVal.push_back(6 * 1 / cmDen[4] * (sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 0)] - sumP[getVectorIndex(2, 2)] * sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 0)] - sumP[getVectorIndex(2, 0)] * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 1)] + sumP[getVectorIndex(2, 0)] * sumP[getVectorIndex(2, 2)] - 4 * sumP[getVectorIndex(2, 1)] * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 0)] + 4 * sumP[getVectorIndex(3, 2)] * sumP[getVectorIndex(1, 0)] + 4 * sumP[getVectorIndex(3, 1)] * sumP[getVectorIndex(1, 1)] + 2 * sumP[getVectorIndex(2, 1)] * sumP[getVectorIndex(2, 1)] - 6 * sumP[getVectorIndex(4, 2)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(8))->FillProfile(centmult, cmVal[9], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDen[4], rn);
  cmVal.push_back(-4 * 1 / cmDen[4] * (sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 0)] - 3 * sumP[getVectorIndex(2, 1)] * sumP[getVectorIndex(1, 0)] * sumP[getVectorIndex(1, 0)] - 3 * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(2, 0)] * sumP[getVectorIndex(1, 0)] + 3 * sumP[getVectorIndex(2, 1)] * sumP[getVectorIndex(2, 0)] + 2 * sumP[getVectorIndex(1, 1)] * sumP[getVectorIndex(3, 0)] + 6 * sumP[getVectorIndex(3, 1)] * sumP[getVectorIndex(1, 0)] - 6 * sumP[getVectorIndex(4, 1)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(9))->FillProfile(centmult, cmVal[10], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDen[4], rn);
}
void FlowPtContainer::fillCMSubeventProfiles(const double& centmult, const double& rn)
{
  // do I need to add an extra return statement here to match fillCMProfiles?
  if (mpar < 1) {
    return;
  }
  const int minSubevents = 2;
  if (nSubevents < minSubevents) {
    return;
  }

  int indOffset = 0;
  for (int im = 1; im <= mpar; im++) {
    indOffset += im;
  }

  std::vector<int> validMpar = {0, 0, 0};
  // 0th order correlation
  cmDenSub[0].push_back(1.);
  cmValSub[0].push_back(1.);
  cmDenSub[nSubevents - 1].push_back(1.);
  cmValSub[nSubevents - 1].push_back(1.);

  cmDenSub[0].push_back(insub[0][getVectorIndex(1, 0)]);                                                                                                                                                                                                                                                                                                                                                                                            // cmDensub[0][1] = sum w
  cmDenSub[0].push_back(insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 0)] - insub[0][getVectorIndex(2, 0)]);                                                                                                                                                                                                                                                                                                                          // cmDensub[0][2] = sum wi wj
  cmDenSub[0].push_back(insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 0)] - 3 * insub[0][getVectorIndex(2, 0)] * insub[0][getVectorIndex(1, 0)] + 2 * insub[0][getVectorIndex(3, 0)]);                                                                                                                                                                                                               // cmDensub[0][3] = sum wi wj wk
  cmDenSub[0].push_back(insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 0)] - 6 * insub[0][getVectorIndex(2, 0)] * insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 0)] + 8 * insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(3, 0)] + 3 * insub[0][getVectorIndex(2, 0)] * insub[0][getVectorIndex(2, 0)] - 6 * insub[0][getVectorIndex(4, 0)]); // cmDensub[0][4] = sum wi wj wk wl

  cmDenSub[nSubevents - 1].push_back(insub[nSubevents - 1][getVectorIndex(1, 0)]);
  cmDenSub[nSubevents - 1].push_back(insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] - insub[nSubevents - 1][getVectorIndex(2, 0)]);
  cmDenSub[nSubevents - 1].push_back(insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] - 3 * insub[nSubevents - 1][getVectorIndex(2, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] + 2 * insub[nSubevents - 1][getVectorIndex(3, 0)]);
  cmDenSub[nSubevents - 1].push_back(insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] - 6 * insub[nSubevents - 1][getVectorIndex(2, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] + 8 * insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(3, 0)] + 3 * insub[nSubevents - 1][getVectorIndex(2, 0)] * insub[nSubevents - 1][getVectorIndex(2, 0)] - 6 * insub[nSubevents - 1][getVectorIndex(4, 0)]);

  if (cmDenSub[0][1] != 0) {
    cmValSub[0].push_back(insub[0][getVectorIndex(1, 1)] / cmDenSub[0][1]);
    dynamic_cast<BootstrapProfile*>(fSubCMList->At(0))->FillProfile(centmult, cmValSub[0][1], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[0][1], rn);
  }
  if (cmDenSub[nSubevents - 1][1] != 0) {
    cmValSub[nSubevents - 1].push_back(insub[nSubevents - 1][getVectorIndex(1, 1)] / cmDenSub[nSubevents - 1][1]);
    dynamic_cast<BootstrapProfile*>(fSubCMList->At(indOffset + 0))->FillProfile(centmult, cmValSub[nSubevents - 1][1], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[nSubevents - 1][1], rn);
  }

  if (mpar >= 2) { // o2-linter: disable=magic-number (greater than order 2)
    if (insub[0][getVectorIndex(2, 0)] != 0 && cmDenSub[0][2] != 0) {
      cmValSub[0].push_back(1 / cmDenSub[0][2] * (insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 1)] - insub[0][getVectorIndex(2, 2)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(1))->FillProfile(centmult, cmValSub[0][2], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[0][2], rn);
      cmValSub[0].push_back(1 / cmDenSub[0][2] * (insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 1)] - insub[0][getVectorIndex(2, 1)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(2))->FillProfile(centmult, cmValSub[0][3], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[0][2], rn);
    }
    if (insub[nSubevents - 1][getVectorIndex(2, 0)] != 0 && cmDenSub[nSubevents - 1][2] != 0) {
      cmValSub[nSubevents - 1].push_back(1 / cmDenSub[nSubevents - 1][2] * (insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] - insub[nSubevents - 1][getVectorIndex(2, 2)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(indOffset + 1))->FillProfile(centmult, cmValSub[nSubevents - 1][2], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[nSubevents - 1][2], rn);
      cmValSub[nSubevents - 1].push_back(1 / cmDenSub[nSubevents - 1][2] * (insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 1)] - insub[nSubevents - 1][getVectorIndex(2, 1)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(indOffset + 2))->FillProfile(centmult, cmValSub[nSubevents - 1][3], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[nSubevents - 1][2], rn);
    }
    validMpar[0] = 1;
  }

  if (mpar >= 3) { // o2-linter: disable=magic-number (greater than order 3)
    if (insub[0][getVectorIndex(3, 0)] != 0 && cmDenSub[0][3] != 0) {
      cmValSub[0].push_back(1 / cmDenSub[0][3] * (insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 1)] - 3 * insub[0][getVectorIndex(2, 2)] * insub[0][getVectorIndex(1, 1)] + 2 * insub[0][getVectorIndex(3, 3)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(3))->FillProfile(centmult, cmValSub[0][4], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[0][3], rn);
      cmValSub[0].push_back(1 / cmDenSub[0][3] * (insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 0)] - 2 * insub[0][getVectorIndex(2, 1)] * insub[0][getVectorIndex(1, 1)] + 2 * insub[0][getVectorIndex(3, 2)] - insub[0][getVectorIndex(2, 2)] * insub[0][getVectorIndex(1, 0)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(4))->FillProfile(centmult, cmValSub[0][5], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[0][3], rn);
      cmValSub[0].push_back(1 / cmDenSub[0][3] * (insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 0)] - 2 * insub[0][getVectorIndex(2, 1)] * insub[0][getVectorIndex(1, 0)] + 2 * insub[0][getVectorIndex(3, 1)] - insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(2, 0)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(5))->FillProfile(centmult, cmValSub[0][6], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[0][3], rn);
    }
    if (insub[nSubevents - 1][getVectorIndex(3, 0)] != 0 && cmDenSub[nSubevents - 1][3] != 0) {
      cmValSub[nSubevents - 1].push_back(1 / cmDenSub[nSubevents - 1][3] * (insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] - 3 * insub[nSubevents - 1][getVectorIndex(2, 2)] * insub[nSubevents - 1][getVectorIndex(1, 1)] + 2 * insub[nSubevents - 1][getVectorIndex(3, 3)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(indOffset + 3))->FillProfile(centmult, cmValSub[nSubevents - 1][4], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[nSubevents - 1][3], rn);
      cmValSub[nSubevents - 1].push_back(1 / cmDenSub[nSubevents - 1][3] * (insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 0)] - 2 * insub[nSubevents - 1][getVectorIndex(2, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] + 2 * insub[nSubevents - 1][getVectorIndex(3, 2)] - insub[nSubevents - 1][getVectorIndex(2, 2)] * insub[nSubevents - 1][getVectorIndex(1, 0)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(indOffset + 4))->FillProfile(centmult, cmValSub[nSubevents - 1][5], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[nSubevents - 1][3], rn);
      cmValSub[nSubevents - 1].push_back(1 / cmDenSub[nSubevents - 1][3] * (insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] - 2 * insub[nSubevents - 1][getVectorIndex(2, 1)] * insub[nSubevents - 1][getVectorIndex(1, 0)] + 2 * insub[nSubevents - 1][getVectorIndex(3, 1)] - insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(2, 0)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(indOffset + 5))->FillProfile(centmult, cmValSub[nSubevents - 1][6], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[nSubevents - 1][3], rn);
    }
    validMpar[1] = 1;
  }
  if (mpar >= 4) { // o2-linter: disable=magic-number (greater than order 4)
    if (insub[0][getVectorIndex(4, 0)] != 0 && cmDenSub[0][4] != 0) {
      cmValSub[0].push_back(1 / cmDenSub[0][4] * (insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 1)] - 6 * insub[0][getVectorIndex(2, 2)] * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 1)] + 3 * insub[0][getVectorIndex(2, 2)] * insub[0][getVectorIndex(2, 2)] + 8 * insub[0][getVectorIndex(3, 3)] * insub[0][getVectorIndex(1, 1)] - 6 * insub[0][getVectorIndex(4, 4)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(6))->FillProfile(centmult, cmValSub[0][7], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[0][4], rn);
      cmValSub[0].push_back(1 / cmDenSub[0][4] * (insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 0)] - 3 * insub[0][getVectorIndex(2, 2)] * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 0)] - 3 * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(2, 1)] + 3 * insub[0][getVectorIndex(2, 2)] * insub[0][getVectorIndex(2, 1)] + 6 * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(3, 2)] - 6 * insub[0][getVectorIndex(4, 3)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(7))->FillProfile(centmult, cmValSub[0][8], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[0][4], rn);
      cmValSub[0].push_back(1 / cmDenSub[0][4] * (insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 0)] - insub[0][getVectorIndex(2, 2)] * insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 0)] - insub[0][getVectorIndex(2, 0)] * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 1)] + insub[0][getVectorIndex(2, 0)] * insub[0][getVectorIndex(2, 2)] - 4 * insub[0][getVectorIndex(2, 1)] * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 0)] + 4 * insub[0][getVectorIndex(3, 2)] * insub[0][getVectorIndex(1, 0)] + 4 * insub[0][getVectorIndex(3, 1)] * insub[0][getVectorIndex(1, 1)] + 2 * insub[0][getVectorIndex(2, 1)] * insub[0][getVectorIndex(2, 1)] - 6 * insub[0][getVectorIndex(4, 2)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(8))->FillProfile(centmult, cmValSub[0][9], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[0][4], rn);
      cmValSub[0].push_back(1 / cmDenSub[0][4] * (insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 0)] - 3 * insub[0][getVectorIndex(2, 1)] * insub[0][getVectorIndex(1, 0)] * insub[0][getVectorIndex(1, 0)] - 3 * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(2, 0)] * insub[0][getVectorIndex(1, 0)] + 3 * insub[0][getVectorIndex(2, 1)] * insub[0][getVectorIndex(2, 0)] + 2 * insub[0][getVectorIndex(1, 1)] * insub[0][getVectorIndex(3, 0)] + 6 * insub[0][getVectorIndex(3, 1)] * insub[0][getVectorIndex(1, 0)] - 6 * insub[0][getVectorIndex(4, 1)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(9))->FillProfile(centmult, cmValSub[0][10], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[0][4], rn);
    }
    if (insub[nSubevents - 1][getVectorIndex(4, 0)] != 0 && cmDenSub[nSubevents - 1][4] != 0) {
      cmValSub[nSubevents - 1].push_back(1 / cmDenSub[nSubevents - 1][4] * (insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] - 6 * insub[nSubevents - 1][getVectorIndex(2, 2)] * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] + 3 * insub[nSubevents - 1][getVectorIndex(2, 2)] * insub[nSubevents - 1][getVectorIndex(2, 2)] + 8 * insub[nSubevents - 1][getVectorIndex(3, 3)] * insub[nSubevents - 1][getVectorIndex(1, 1)] - 6 * insub[nSubevents - 1][getVectorIndex(4, 4)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(indOffset + 6))->FillProfile(centmult, cmValSub[nSubevents - 1][7], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[nSubevents - 1][4], rn);
      cmValSub[nSubevents - 1].push_back(1 / cmDenSub[nSubevents - 1][4] * (insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 0)] - 3 * insub[nSubevents - 1][getVectorIndex(2, 2)] * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 0)] - 3 * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(2, 1)] + 3 * insub[nSubevents - 1][getVectorIndex(2, 2)] * insub[nSubevents - 1][getVectorIndex(2, 1)] + 6 * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(3, 2)] - 6 * insub[nSubevents - 1][getVectorIndex(4, 3)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(indOffset + 7))->FillProfile(centmult, cmValSub[nSubevents - 1][8], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[nSubevents - 1][4], rn);
      cmValSub[nSubevents - 1].push_back(1 / cmDenSub[nSubevents - 1][4] * (insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] - insub[nSubevents - 1][getVectorIndex(2, 2)] * insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] - insub[nSubevents - 1][getVectorIndex(2, 0)] * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] + insub[nSubevents - 1][getVectorIndex(2, 0)] * insub[nSubevents - 1][getVectorIndex(2, 2)] - 4 * insub[nSubevents - 1][getVectorIndex(2, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 0)] + 4 * insub[nSubevents - 1][getVectorIndex(3, 2)] * insub[nSubevents - 1][getVectorIndex(1, 0)] + 4 * insub[nSubevents - 1][getVectorIndex(3, 1)] * insub[nSubevents - 1][getVectorIndex(1, 1)] + 2 * insub[nSubevents - 1][getVectorIndex(2, 1)] * insub[nSubevents - 1][getVectorIndex(2, 1)] - 6 * insub[nSubevents - 1][getVectorIndex(4, 2)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(indOffset + 8))->FillProfile(centmult, cmValSub[nSubevents - 1][9], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[nSubevents - 1][4], rn);
      cmValSub[nSubevents - 1].push_back(1 / cmDenSub[nSubevents - 1][4] * (insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] - 3 * insub[nSubevents - 1][getVectorIndex(2, 1)] * insub[nSubevents - 1][getVectorIndex(1, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] - 3 * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(2, 0)] * insub[nSubevents - 1][getVectorIndex(1, 0)] + 3 * insub[nSubevents - 1][getVectorIndex(2, 1)] * insub[nSubevents - 1][getVectorIndex(2, 0)] + 2 * insub[nSubevents - 1][getVectorIndex(1, 1)] * insub[nSubevents - 1][getVectorIndex(3, 0)] + 6 * insub[nSubevents - 1][getVectorIndex(3, 1)] * insub[nSubevents - 1][getVectorIndex(1, 0)] - 6 * insub[nSubevents - 1][getVectorIndex(4, 1)]));
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(indOffset + 9))->FillProfile(centmult, cmValSub[nSubevents - 1][10], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[nSubevents - 1][4], rn);
    }
    validMpar[2] = 1;
  }
  // Fill cross terms
  for (int m = 2; m <= centralMomentMaxOrder; ++m) {
    if (validMpar[m - 2] == 0) {
      continue;
    }
    for (int first = 1; first < m; ++first) {
      for (int second = first; second < m; ++second) {
        if (first > second) {
          continue;
        }
        int fourth = m - second;
        for (int third = 1; third < m; ++third) {
          if (third > fourth) {
            continue;
          }
          if (insub[0][getVectorIndex(m, 0)] != 0 && insub[nSubevents - 1][getVectorIndex(m, 0)] != 0 && cmDenSub[0][m] * cmDenSub[nSubevents - 1][m] != 0) {
            dynamic_cast<BootstrapProfile*>(fSubCMList->FindObject(Form("cm%i_%i%isub1_%i%isub2", m, first, second, third, fourth)))->FillProfile(centmult, cmValSub[0][second * (second - 1) / 2 + second - first + 1] * cmValSub[nSubevents - 1][fourth * (fourth - 1) / 2 + fourth - third + 1], (fEventWeight == EventWeight::UnityWeight) ? 1.0 : cmDenSub[0][m] * cmDenSub[nSubevents - 1][m], rn);
          }
        }
      }
    }
  }
}
void FlowPtContainer::fillArray(FillType a, FillType b, double c, double d)
{
  const int nTerms = 3 * 3 * 5 * 5;
  for (int idx = 0; idx < nTerms; ++idx) {
    int i = idx % 3;
    int j = ((idx - i) / 3) % 3;
    int k = ((idx - j * 3 - i) / 9) % 5;
    int l = ((idx - k * 9 - j * 3 - i) / 45) % 5;
    if (std::holds_alternative<std::complex<double>>(a) && std::holds_alternative<std::complex<double>>(b)) {
      arr[idx] += std::pow(std::get<0>(a), i) * std::pow(std::get<0>(b), j) * std::pow(c, k) * std::pow(d, l);
    } else if (std::holds_alternative<double>(a) && std::holds_alternative<double>(b)) {
      warr[idx] += std::pow(std::get<1>(a), i) * std::pow(std::get<1>(b), j) * std::pow(c, k) * std::pow(d, l);
    } else {
      LOGF(error, "FillType variant should hold same type for a and b during single function c");
    }
  }
}
template <typename T>
double FlowPtContainer::getStdAABBCC(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> c = inarr[getVectorIndex(0, 0, 1, 0)];
  std::complex<double> aa = inarr[getVectorIndex(2, 0, 0, 0)];
  std::complex<double> bb = inarr[getVectorIndex(0, 2, 0, 0)];
  std::complex<double> cc = inarr[getVectorIndex(0, 0, 2, 0)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ac = inarr[getVectorIndex(1, 0, 1, 0)];
  std::complex<double> bc = inarr[getVectorIndex(0, 1, 1, 0)];
  std::complex<double> aab = inarr[getVectorIndex(2, 1, 0, 0)];
  std::complex<double> aac = inarr[getVectorIndex(2, 0, 1, 0)];
  std::complex<double> abb = inarr[getVectorIndex(1, 2, 0, 0)];
  std::complex<double> acc = inarr[getVectorIndex(1, 0, 2, 0)];
  std::complex<double> abc = inarr[getVectorIndex(1, 1, 1, 0)];
  std::complex<double> bbc = inarr[getVectorIndex(0, 2, 1, 0)];
  std::complex<double> bcc = inarr[getVectorIndex(0, 1, 2, 0)];
  std::complex<double> aabb = inarr[getVectorIndex(2, 2, 0, 0)];
  std::complex<double> aacc = inarr[getVectorIndex(2, 0, 2, 0)];
  std::complex<double> aabc = inarr[getVectorIndex(2, 1, 1, 0)];
  std::complex<double> abbc = inarr[getVectorIndex(1, 2, 1, 0)];
  std::complex<double> abcc = inarr[getVectorIndex(1, 1, 2, 0)];
  std::complex<double> bbcc = inarr[getVectorIndex(0, 2, 2, 0)];
  std::complex<double> aabbc = inarr[getVectorIndex(2, 2, 1, 0)];
  std::complex<double> aabcc = inarr[getVectorIndex(2, 1, 2, 0)];
  std::complex<double> abbcc = inarr[getVectorIndex(1, 2, 2, 0)];
  std::complex<double> aabbcc = inarr[getVectorIndex(2, 2, 2, 0)];
  return (a * a * b * b * c * c - aa * b * b * c * c - a * a * bb * c * c - a * a * b * b * cc - 4. * a * ab * b * c * c -
          4. * a * ac * b * b * c - 4. * a * a * b * bc * c + 4. * aab * b * c * c + 4. * aac * b * b * c +
          4. * a * abb * c * c + 4. * a * acc * b * b + 4. * a * a * bbc * c + 4. * a * a * b * bcc +
          16. * a * abc * b * c + aa * bb * c * c + aa * b * b * cc + a * a * bb * cc + 2. * ab * ab * c * c +
          2. * ac * ac * b * b + 2. * a * a * bc * bc + 4. * aa * b * bc * c + 4. * a * ac * bb * c +
          4. * a * ab * b * cc + 8. * ab * ac * b * c + 8. * a * ab * bc * c + 8. * a * ac * b * bc - 6. * aabb * c * c -
          24. * aabc * b * c - 6. * aacc * b * b - 24. * abbc * a * c - 24. * abcc * a * b - 6. * bbcc * a * a -
          8. * aab * bc * c - 8. * aac * b * bc - 4. * aac * bb * c - 4. * aab * b * cc - 8. * abb * ac * c -
          4. * abb * a * cc - 8. * acc * ab * b - 4. * acc * a * bb - 8. * bbc * a * ac - 4. * bbc * aa * c -
          8. * bcc * a * ab - 4. * bcc * aa * b - 16. * abc * ab * c - 16. * abc * ac * b - 16. * abc * a * bc -
          aa * bb * cc - 2. * ab * ab * cc - 2. * ac * ac * bb - 2. * bc * bc * aa - 8. * ab * ac * bc +
          48. * aabbc * c + 48. * aabcc * b + 48. * abbcc * a + 6. * aabb * cc + 6. * aacc * bb +
          6. * bbcc * aa + 24. * aabc * bc + 24. * abbc * ac + 24. * abcc * ab + 8. * aab * bcc +
          8. * aac * bbc + 8. * abb * acc + 16. * abc * abc - 120. * aabbcc)
    .real();
}
template <typename T>
double FlowPtContainer::getStdAABBCD(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> c = inarr[getVectorIndex(0, 0, 1, 0)];
  std::complex<double> d = inarr[getVectorIndex(0, 0, 0, 1)];
  std::complex<double> aa = inarr[getVectorIndex(2, 0, 0, 0)];
  std::complex<double> bb = inarr[getVectorIndex(0, 2, 0, 0)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ac = inarr[getVectorIndex(1, 0, 1, 0)];
  std::complex<double> ad = inarr[getVectorIndex(1, 0, 0, 1)];
  std::complex<double> bc = inarr[getVectorIndex(0, 1, 1, 0)];
  std::complex<double> bd = inarr[getVectorIndex(0, 1, 0, 1)];
  std::complex<double> cd = inarr[getVectorIndex(0, 0, 1, 1)];
  std::complex<double> aab = inarr[getVectorIndex(2, 1, 0, 0)];
  std::complex<double> aac = inarr[getVectorIndex(2, 0, 1, 0)];
  std::complex<double> aad = inarr[getVectorIndex(2, 0, 0, 1)];
  std::complex<double> abb = inarr[getVectorIndex(1, 2, 0, 0)];
  std::complex<double> abc = inarr[getVectorIndex(1, 1, 1, 0)];
  std::complex<double> abd = inarr[getVectorIndex(1, 1, 0, 1)];
  std::complex<double> acd = inarr[getVectorIndex(1, 0, 1, 1)];
  std::complex<double> bbc = inarr[getVectorIndex(0, 2, 1, 0)];
  std::complex<double> bbd = inarr[getVectorIndex(0, 2, 0, 1)];
  std::complex<double> bcd = inarr[getVectorIndex(0, 1, 1, 1)];
  std::complex<double> aabb = inarr[getVectorIndex(2, 2, 0, 0)];
  std::complex<double> aabc = inarr[getVectorIndex(2, 1, 1, 0)];
  std::complex<double> aabd = inarr[getVectorIndex(2, 1, 0, 1)];
  std::complex<double> aacd = inarr[getVectorIndex(2, 0, 1, 1)];
  std::complex<double> abbc = inarr[getVectorIndex(1, 2, 1, 0)];
  std::complex<double> abbd = inarr[getVectorIndex(1, 2, 0, 1)];
  std::complex<double> abcd = inarr[getVectorIndex(1, 1, 1, 1)];
  std::complex<double> bbcd = inarr[getVectorIndex(0, 2, 1, 1)];
  std::complex<double> aabbc = inarr[getVectorIndex(2, 2, 1, 0)];
  std::complex<double> aabbd = inarr[getVectorIndex(2, 2, 0, 1)];
  std::complex<double> aabcd = inarr[getVectorIndex(2, 1, 1, 1)];
  std::complex<double> abbcd = inarr[getVectorIndex(1, 2, 1, 1)];
  std::complex<double> aabbcd = inarr[getVectorIndex(2, 2, 1, 1)];
  return (-120. * aabbcd + 48. * a * abbcd + 24. * ab * abcd + 16. * abc * abd + 12. * abbd * ac +
          8. * abb * acd + 12. * abbc * ad + 48. * aabcd * b - 24. * a * abcd * b - 8. * abd * ac * b -
          8. * ab * acd * b - 8. * abc * ad * b - 6. * aacd * b * b + 4. * a * acd * b * b + 2. * ac * ad * b * b +
          6. * aacd * bb - 4. * a * acd * bb - 2. * ac * ad * bb + 4. * aad * bbc - 4. * a * ad * bbc -
          6. * a * a * bbcd + 6. * aa * bbcd + 4. * aac * bbd - 4. * a * ac * bbd + 12. * aabd * bc -
          8. * a * abd * bc - 4. * ab * ad * bc - 4. * aad * b * bc + 4. * a * ad * b * bc + 8. * aab * bcd -
          8. * a * ab * bcd + 4. * a * a * b * bcd - 4. * aa * b * bcd + 12. * aabc * bd - 8. * a * abc * bd -
          4. * ab * ac * bd - 4. * aac * b * bd + 4. * a * ac * b * bd + 2. * a * a * bc * bd - 2. * aa * bc * bd +
          24. * aabbd * c - 12. * a * abbd * c - 8. * ab * abd * c - 4. * abb * ad * c - 12. * aabd * b * c +
          8. * a * abd * b * c + 4. * ab * ad * b * c + 2. * aad * b * b * c - 2. * a * ad * b * b * c -
          2. * aad * bb * c + 2. * a * ad * bb * c + 2. * a * a * bbd * c - 2. * aa * bbd * c - 4. * aab * bd * c +
          4. * a * ab * bd * c - 2. * a * a * b * bd * c + 2. * aa * b * bd * c + 6. * aabb * cd - 2. * ab * ab * cd -
          4. * a * abb * cd - 4. * aab * b * cd + 4. * a * ab * b * cd - a * a * b * b * cd + aa * b * b * cd +
          a * a * bb * cd - aa * bb * cd + 24. * aabbc * d - 12. * a * abbc * d - 8. * ab * abc * d -
          4. * abb * ac * d - 12. * aabc * b * d + 8. * a * abc * b * d + 4. * ab * ac * b * d + 2. * aac * b * b * d -
          2. * a * ac * b * b * d - 2. * aac * bb * d + 2. * a * ac * bb * d + 2. * a * a * bbc * d - 2. * aa * bbc * d -
          4. * aab * bc * d + 4. * a * ab * bc * d - 2. * a * a * b * bc * d + 2. * aa * b * bc * d - 6. * aabb * c * d +
          2. * ab * ab * c * d + 4. * a * abb * c * d + 4. * aab * b * c * d - 4. * a * ab * b * c * d +
          a * a * b * b * c * d - aa * b * b * c * d - a * a * bb * c * d + aa * bb * c * d)
    .real();
}
template <typename T>
double FlowPtContainer::getStdAABBDD(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> d = inarr[getVectorIndex(0, 0, 0, 1)];
  std::complex<double> aa = inarr[getVectorIndex(2, 0, 0, 0)];
  std::complex<double> bb = inarr[getVectorIndex(0, 2, 0, 0)];
  std::complex<double> dd = inarr[getVectorIndex(0, 0, 0, 2)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ad = inarr[getVectorIndex(1, 0, 0, 1)];
  std::complex<double> bd = inarr[getVectorIndex(0, 1, 0, 1)];
  std::complex<double> aab = inarr[getVectorIndex(2, 1, 0, 0)];
  std::complex<double> aad = inarr[getVectorIndex(2, 0, 0, 1)];
  std::complex<double> abb = inarr[getVectorIndex(1, 2, 0, 0)];
  std::complex<double> add = inarr[getVectorIndex(1, 0, 0, 2)];
  std::complex<double> abd = inarr[getVectorIndex(1, 1, 0, 1)];
  std::complex<double> bbd = inarr[getVectorIndex(0, 2, 0, 1)];
  std::complex<double> bdd = inarr[getVectorIndex(0, 1, 0, 2)];
  std::complex<double> aabb = inarr[getVectorIndex(2, 2, 0, 0)];
  std::complex<double> aadd = inarr[getVectorIndex(2, 0, 0, 2)];
  std::complex<double> aabd = inarr[getVectorIndex(2, 1, 0, 1)];
  std::complex<double> abbd = inarr[getVectorIndex(1, 2, 0, 1)];
  std::complex<double> abdd = inarr[getVectorIndex(1, 1, 0, 2)];
  std::complex<double> bbdd = inarr[getVectorIndex(0, 2, 0, 2)];
  std::complex<double> aabbd = inarr[getVectorIndex(2, 2, 0, 1)];
  std::complex<double> aabdd = inarr[getVectorIndex(2, 1, 0, 2)];
  std::complex<double> abbdd = inarr[getVectorIndex(1, 2, 0, 2)];
  std::complex<double> aabbdd = inarr[getVectorIndex(2, 2, 0, 2)];
  return (-120. * aabbdd + 48. * a * abbdd + 16. * abd * abd + 24. * ab * abdd + 24. * abbd * ad +
          8. * abb * add + 48. * aabdd * b - 24. * a * abdd * b - 16. * abd * ad * b - 8. * ab * add * b -
          6. * aadd * b * b + 2. * ad * ad * b * b + 4. * a * add * b * b + 6. * aadd * bb - 2. * ad * ad * bb -
          4. * a * add * bb + 8. * aad * bbd - 8. * a * ad * bbd - 6. * a * a * bbdd + 6. * aa * bbdd +
          24. * aabd * bd - 16. * a * abd * bd - 8. * ab * ad * bd - 8. * aad * b * bd + 8. * a * ad * b * bd +
          2. * a * a * bd * bd - 2. * aa * bd * bd + 8. * aab * bdd - 8. * a * ab * bdd + 4. * a * a * b * bdd -
          4. * aa * b * bdd + 48. * aabbd * d - 24. * a * abbd * d - 16. * ab * abd * d - 8. * abb * ad * d -
          24. * aabd * b * d + 16. * a * abd * b * d + 8. * ab * ad * b * d + 4. * aad * b * b * d -
          4. * a * ad * b * b * d - 4. * aad * bb * d + 4. * a * ad * bb * d + 4. * a * a * bbd * d - 4. * aa * bbd * d -
          8. * aab * bd * d + 8. * a * ab * bd * d - 4. * a * a * b * bd * d + 4. * aa * b * bd * d - 6. * aabb * d * d +
          2. * ab * ab * d * d + 4. * a * abb * d * d + 4. * aab * b * d * d - 4. * a * ab * b * d * d +
          a * a * b * b * d * d - aa * b * b * d * d - a * a * bb * d * d + aa * bb * d * d + 6. * aabb * dd -
          2. * ab * ab * dd - 4. * a * abb * dd - 4. * aab * b * dd + 4. * a * ab * b * dd - a * a * b * b * dd +
          aa * b * b * dd + a * a * bb * dd - aa * bb * dd)
    .real();
}
template <typename T>
double FlowPtContainer::getStdAABBC(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> c = inarr[getVectorIndex(0, 0, 1, 0)];
  std::complex<double> aa = inarr[getVectorIndex(2, 0, 0, 0)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ac = inarr[getVectorIndex(1, 0, 1, 0)];
  std::complex<double> bb = inarr[getVectorIndex(0, 2, 0, 0)];
  std::complex<double> bc = inarr[getVectorIndex(0, 1, 1, 0)];
  std::complex<double> aab = inarr[getVectorIndex(2, 1, 0, 0)];
  std::complex<double> aac = inarr[getVectorIndex(2, 0, 1, 0)];
  std::complex<double> abb = inarr[getVectorIndex(1, 2, 0, 0)];
  std::complex<double> abc = inarr[getVectorIndex(1, 1, 1, 0)];
  std::complex<double> bbc = inarr[getVectorIndex(0, 2, 1, 0)];
  std::complex<double> aabb = inarr[getVectorIndex(2, 2, 0, 0)];
  std::complex<double> aabc = inarr[getVectorIndex(2, 1, 1, 0)];
  std::complex<double> abbc = inarr[getVectorIndex(1, 2, 1, 0)];
  std::complex<double> aabbc = inarr[getVectorIndex(2, 2, 1, 0)];
  return (a * a * b * b * c - aa * b * b * c - a * a * bb * c - 4. * ab * a * b * c - 2. * a * ac * b * b - 2. * a * a * bc * b + 2. * ab * ab * c + 4. * ab * ac * b + 4. * ab * bc * a + 8. * abc * a * b + 4. * aab * b * c + 2. * aac * b * b + 4. * abb * a * c + 2. * bbc * a * a + aa * bb * c + 2. * aa * b * bc + 2. * bb * a * ac - 12. * aabc * b - 12. * abbc * a - 6. * aabb * c - 8. * abc * ab - 2. * bbc * aa - 2. * aac * bb - 4. * aab * bc - 4. * abb * ac + 24. * aabbc).real();
}
template <typename T>
double FlowPtContainer::getStdAABBD(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> d = inarr[getVectorIndex(0, 0, 0, 1)];
  std::complex<double> aa = inarr[getVectorIndex(2, 0, 0, 0)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ad = inarr[getVectorIndex(1, 0, 0, 1)];
  std::complex<double> bb = inarr[getVectorIndex(0, 2, 0, 0)];
  std::complex<double> bd = inarr[getVectorIndex(0, 1, 0, 1)];
  std::complex<double> aab = inarr[getVectorIndex(2, 1, 0, 0)];
  std::complex<double> aad = inarr[getVectorIndex(2, 0, 0, 1)];
  std::complex<double> abb = inarr[getVectorIndex(1, 2, 0, 0)];
  std::complex<double> abd = inarr[getVectorIndex(1, 1, 0, 1)];
  std::complex<double> bbd = inarr[getVectorIndex(0, 2, 0, 1)];
  std::complex<double> aabb = inarr[getVectorIndex(2, 2, 0, 0)];
  std::complex<double> aabd = inarr[getVectorIndex(2, 1, 0, 1)];
  std::complex<double> abbd = inarr[getVectorIndex(1, 2, 0, 1)];
  std::complex<double> aabbd = inarr[getVectorIndex(2, 2, 0, 1)];
  return (a * a * b * b * d - aa * b * b * d - a * a * bb * d - 4. * ab * a * b * d - 2. * a * ad * b * b - 2. * a * a * bd * b + 2. * ab * ab * d + 4. * ab * ad * b + 4. * ab * bd * a + 8. * abd * a * b + 4. * aab * b * d + 2. * aad * b * b + 4. * abb * a * d + 2. * bbd * a * a + aa * bb * d + 2. * aa * b * bd + 2. * bb * a * ad - 12. * aabd * b - 12. * abbd * a - 6. * aabb * d - 8. * abd * ab - 2. * bbd * aa - 2. * aad * bb - 4. * aab * bd - 4. * abb * ad + 24. * aabbd).real();
}
template <typename T>
double FlowPtContainer::getStdABCC(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> c = inarr[getVectorIndex(0, 0, 1, 0)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ac = inarr[getVectorIndex(1, 0, 1, 0)];
  std::complex<double> bc = inarr[getVectorIndex(0, 1, 1, 0)];
  std::complex<double> cc = inarr[getVectorIndex(0, 0, 2, 0)];
  std::complex<double> abc = inarr[getVectorIndex(1, 1, 1, 0)];
  std::complex<double> acc = inarr[getVectorIndex(1, 0, 2, 0)];
  std::complex<double> bcc = inarr[getVectorIndex(0, 1, 2, 0)];
  std::complex<double> abcc = inarr[getVectorIndex(1, 1, 2, 0)];
  return (a * b * c * c - a * b * cc - 2. * a * bc * c - 2. * ac * b * c - ab * c * c + 2. * acc * b + 2. * a * bcc + 4. * abc * c + ab * cc + 2. * ac * bc - 6. * abcc).real();
}
template <typename T>
double FlowPtContainer::getStdABCD(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> c = inarr[getVectorIndex(0, 0, 1, 0)];
  std::complex<double> d = inarr[getVectorIndex(0, 0, 0, 1)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ac = inarr[getVectorIndex(1, 0, 1, 0)];
  std::complex<double> ad = inarr[getVectorIndex(1, 0, 0, 1)];
  std::complex<double> bc = inarr[getVectorIndex(0, 1, 1, 0)];
  std::complex<double> bd = inarr[getVectorIndex(0, 1, 0, 1)];
  std::complex<double> cd = inarr[getVectorIndex(0, 0, 1, 1)];
  std::complex<double> abc = inarr[getVectorIndex(1, 1, 1, 0)];
  std::complex<double> abd = inarr[getVectorIndex(1, 1, 0, 1)];
  std::complex<double> acd = inarr[getVectorIndex(1, 0, 1, 1)];
  std::complex<double> bcd = inarr[getVectorIndex(0, 1, 1, 1)];
  std::complex<double> abcd = inarr[getVectorIndex(1, 1, 0, 1)];
  return (-6. * abcd + 2. * acd * b + ad * bc + 2. * a * bcd + ac * bd + 2. * abd * c - ad * b * c -
          a * bd * c + ab * cd - a * b * cd + 2. * abc * d - ac * b * d - a * bc * d - ab * c * d +
          a * b * c * d)
    .real();
}
template <typename T>
double FlowPtContainer::getStdABDD(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> d = inarr[getVectorIndex(0, 0, 0, 1)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ad = inarr[getVectorIndex(1, 0, 0, 1)];
  std::complex<double> bd = inarr[getVectorIndex(0, 1, 0, 1)];
  std::complex<double> dd = inarr[getVectorIndex(0, 0, 0, 2)];
  std::complex<double> abd = inarr[getVectorIndex(1, 1, 0, 1)];
  std::complex<double> add = inarr[getVectorIndex(1, 0, 0, 2)];
  std::complex<double> bdd = inarr[getVectorIndex(0, 1, 0, 2)];
  std::complex<double> abdd = inarr[getVectorIndex(1, 1, 0, 2)];
  return (a * b * d * d - a * b * dd - 2. * a * bd * d - 2. * ad * b * d - ab * d * d + 2. * add * b + 2. * a * bdd + 4. * abd * d + ab * dd + 2. * ad * bd - 6. * abdd).real();
}
template <typename T>
double FlowPtContainer::getStdABC(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> c = inarr[getVectorIndex(0, 0, 1, 0)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ac = inarr[getVectorIndex(1, 0, 1, 0)];
  std::complex<double> bc = inarr[getVectorIndex(0, 1, 1, 0)];
  std::complex<double> abc = inarr[getVectorIndex(1, 1, 1, 0)];
  return (a * b * c - ab * c - ac * b - a * bc + 2. * abc).real();
}
template <typename T>
double FlowPtContainer::getStdABD(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> d = inarr[getVectorIndex(0, 0, 0, 1)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ad = inarr[getVectorIndex(1, 0, 0, 1)];
  std::complex<double> bd = inarr[getVectorIndex(0, 1, 0, 1)];
  std::complex<double> abd = inarr[getVectorIndex(1, 1, 0, 1)];
  return (a * b * d - ab * d - ad * b - a * bd + 2. * abd).real();
}
template <typename T>
double FlowPtContainer::getStdABCCCC(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> c = inarr[getVectorIndex(0, 0, 1, 0)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ac = inarr[getVectorIndex(1, 0, 1, 0)];
  std::complex<double> bc = inarr[getVectorIndex(0, 1, 1, 0)];
  std::complex<double> cc = inarr[getVectorIndex(0, 0, 2, 0)];
  std::complex<double> abc = inarr[getVectorIndex(1, 1, 1, 0)];
  std::complex<double> acc = inarr[getVectorIndex(1, 0, 2, 0)];
  std::complex<double> bcc = inarr[getVectorIndex(0, 1, 2, 0)];
  std::complex<double> ccc = inarr[getVectorIndex(0, 0, 3, 0)];
  std::complex<double> abcc = inarr[getVectorIndex(1, 1, 2, 0)];
  std::complex<double> accc = inarr[getVectorIndex(1, 0, 3, 0)];
  std::complex<double> bccc = inarr[getVectorIndex(0, 1, 3, 0)];
  std::complex<double> cccc = inarr[getVectorIndex(0, 0, 4, 0)];
  std::complex<double> abccc = inarr[getVectorIndex(1, 1, 3, 0)];
  std::complex<double> acccc = inarr[getVectorIndex(1, 0, 4, 0)];
  std::complex<double> bcccc = inarr[getVectorIndex(0, 1, 4, 0)];
  std::complex<double> abcccc = inarr[getVectorIndex(1, 1, 4, 0)];
  return (-120. * abcccc + 24. * acccc * b + 24. * accc * bc + 24. * acc * bcc +
          24. * ac * bccc + 24. * a * bcccc + 96. * abccc * c - 24. * accc * b * c - 24. * acc * bc * c -
          24. * ac * bcc * c - 24. * a * bccc * c - 36. * abcc * c * c + 12. * acc * b * c * c +
          12. * ac * bc * c * c + 12. * a * bcc * c * c + 8. * abc * c * c * c - 4. * ac * b * c * c * c - 4. * a * bc * c * c * c -
          ab * c * c * c * c + a * b * c * c * c * c + 36. * abcc * cc - 12. * acc * b * cc - 12. * ac * bc * cc - 12. * a * bcc * cc - 24. * abc * c * cc + 12. * ac * b * c * cc + 12. * a * bc * c * cc + 6. * ab * c * c * cc - 6. * a * b * c * c * cc - 3. * ab * cc * cc + 3. * a * b * cc * cc + 16. * abc * ccc - 8. * ac * b * ccc - 8. * a * bc * ccc - 8. * ab * c * ccc + 8. * a * b * c * ccc + 6. * ab * cccc - 6. * a * b * cccc)
    .real();
}
template <typename T>
double FlowPtContainer::getStdABCCCD(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> c = inarr[getVectorIndex(0, 0, 1, 0)];
  std::complex<double> d = inarr[getVectorIndex(0, 0, 0, 1)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ac = inarr[getVectorIndex(1, 0, 1, 0)];
  std::complex<double> ad = inarr[getVectorIndex(1, 0, 0, 1)];
  std::complex<double> bc = inarr[getVectorIndex(0, 1, 1, 0)];
  std::complex<double> bd = inarr[getVectorIndex(0, 1, 0, 1)];
  std::complex<double> cc = inarr[getVectorIndex(0, 0, 2, 0)];
  std::complex<double> cd = inarr[getVectorIndex(0, 0, 1, 1)];
  std::complex<double> abc = inarr[getVectorIndex(1, 1, 1, 0)];
  std::complex<double> abd = inarr[getVectorIndex(1, 1, 0, 1)];
  std::complex<double> acc = inarr[getVectorIndex(1, 0, 2, 0)];
  std::complex<double> acd = inarr[getVectorIndex(1, 0, 1, 1)];
  std::complex<double> bcc = inarr[getVectorIndex(0, 1, 2, 0)];
  std::complex<double> bcd = inarr[getVectorIndex(0, 1, 1, 1)];
  std::complex<double> ccc = inarr[getVectorIndex(0, 0, 3, 0)];
  std::complex<double> ccd = inarr[getVectorIndex(0, 0, 2, 1)];
  std::complex<double> abcc = inarr[getVectorIndex(1, 1, 2, 0)];
  std::complex<double> abcd = inarr[getVectorIndex(1, 1, 1, 1)];
  std::complex<double> accc = inarr[getVectorIndex(1, 0, 3, 0)];
  std::complex<double> accd = inarr[getVectorIndex(1, 0, 2, 1)];
  std::complex<double> bccc = inarr[getVectorIndex(0, 1, 3, 0)];
  std::complex<double> bccd = inarr[getVectorIndex(0, 1, 2, 1)];
  std::complex<double> cccd = inarr[getVectorIndex(0, 0, 3, 1)];
  std::complex<double> abccc = inarr[getVectorIndex(1, 1, 3, 0)];
  std::complex<double> abccd = inarr[getVectorIndex(1, 1, 2, 1)];
  std::complex<double> acccd = inarr[getVectorIndex(1, 0, 3, 1)];
  std::complex<double> bcccd = inarr[getVectorIndex(0, 1, 3, 1)];
  std::complex<double> abcccd = inarr[getVectorIndex(1, 1, 3, 1)];
  return (-120. * abcccd + 24. * acccd * b + 18. * accd * bc + 12. * acd * bcc + 6. * ad * bccc +
          24. * a * bcccd + 18. * ac * bccd + 12. * acc * bcd + 6. * accc * bd + 72. * abccd * c -
          18. * accd * b * c - 12. * acd * bc * c - 6. * ad * bcc * c - 18. * a * bccd * c - 12. * ac * bcd * c -
          6. * acc * bd * c - 18. * abcd * c * c + 6. * acd * b * c * c + 3. * ad * bc * c * c + 6. * a * bcd * c * c +
          3. * ac * bd * c * c + 2. * abd * c * c * c - ad * b * c * c * c - a * bd * c * c * c + 18. * abcd * cc -
          6. * acd * b * cc - 3. * ad * bc * cc - 6. * a * bcd * cc - 3. * ac * bd * cc - 6. * abd * c * cc +
          3. * ad * b * c * cc + 3. * a * bd * c * cc + 4. * abd * ccc - 2. * ad * b * ccc - 2. * a * bd * ccc +
          6. * ab * cccd - 6. * a * b * cccd + 12. * abc * ccd - 6. * ac * b * ccd - 6. * a * bc * ccd -
          6. * ab * c * ccd + 6. * a * b * c * ccd + 18. * abcc * cd - 6. * acc * b * cd - 6. * ac * bc * cd -
          6. * a * bcc * cd - 12. * abc * c * cd + 6. * ac * b * c * cd + 6. * a * bc * c * cd + 3. * ab * c * c * cd -
          3. * a * b * c * c * cd - 3. * ab * cc * cd + 3. * a * b * cc * cd + 24. * abccc * d - 6. * accc * b * d -
          6. * acc * bc * d - 6. * ac * bcc * d - 6. * a * bccc * d - 18. * abcc * c * d + 6. * acc * b * c * d +
          6. * ac * bc * c * d + 6. * a * bcc * c * d + 6. * abc * c * c * d - 3. * ac * b * c * c * d -
          3. * a * bc * c * c * d - ab * c * c * c * d + a * b * c * c * c * d - 6. * abc * cc * d + 3. * ac * b * cc * d +
          3. * a * bc * cc * d + 3. * ab * c * cc * d - 3. * a * b * c * cc * d - 2. * ab * ccc * d + 2. * a * b * ccc * d)
    .real();
}
template <typename T>
double FlowPtContainer::getStdABCCDD(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> c = inarr[getVectorIndex(0, 0, 1, 0)];
  std::complex<double> d = inarr[getVectorIndex(0, 0, 0, 1)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ac = inarr[getVectorIndex(1, 0, 1, 0)];
  std::complex<double> ad = inarr[getVectorIndex(1, 0, 0, 1)];
  std::complex<double> bc = inarr[getVectorIndex(0, 1, 1, 0)];
  std::complex<double> bd = inarr[getVectorIndex(0, 1, 0, 1)];
  std::complex<double> cc = inarr[getVectorIndex(0, 0, 2, 0)];
  std::complex<double> cd = inarr[getVectorIndex(0, 0, 1, 1)];
  std::complex<double> dd = inarr[getVectorIndex(0, 0, 0, 2)];
  std::complex<double> abc = inarr[getVectorIndex(1, 1, 1, 0)];
  std::complex<double> abd = inarr[getVectorIndex(1, 1, 0, 1)];
  std::complex<double> acc = inarr[getVectorIndex(1, 0, 2, 0)];
  std::complex<double> acd = inarr[getVectorIndex(1, 0, 1, 1)];
  std::complex<double> add = inarr[getVectorIndex(1, 0, 0, 2)];
  std::complex<double> bcc = inarr[getVectorIndex(0, 1, 2, 0)];
  std::complex<double> bcd = inarr[getVectorIndex(0, 1, 1, 1)];
  std::complex<double> bdd = inarr[getVectorIndex(0, 1, 0, 2)];
  std::complex<double> ccd = inarr[getVectorIndex(0, 0, 2, 1)];
  std::complex<double> cdd = inarr[getVectorIndex(0, 0, 1, 2)];
  std::complex<double> abcc = inarr[getVectorIndex(1, 1, 2, 0)];
  std::complex<double> abcd = inarr[getVectorIndex(1, 1, 1, 1)];
  std::complex<double> abdd = inarr[getVectorIndex(1, 1, 0, 2)];
  std::complex<double> accd = inarr[getVectorIndex(1, 0, 2, 1)];
  std::complex<double> acdd = inarr[getVectorIndex(1, 0, 1, 2)];
  std::complex<double> bccd = inarr[getVectorIndex(0, 1, 2, 1)];
  std::complex<double> bcdd = inarr[getVectorIndex(0, 1, 1, 2)];
  std::complex<double> ccdd = inarr[getVectorIndex(0, 0, 2, 2)];
  std::complex<double> abccd = inarr[getVectorIndex(1, 1, 2, 1)];
  std::complex<double> abcdd = inarr[getVectorIndex(1, 1, 1, 2)];
  std::complex<double> accdd = inarr[getVectorIndex(1, 0, 2, 2)];
  std::complex<double> bccdd = inarr[getVectorIndex(0, 1, 2, 2)];
  std::complex<double> abccdd = inarr[getVectorIndex(1, 1, 2, 2)];
  return (-120. * abccdd + 24. * accdd * b + 12. * acdd * bc + 4. * add * bcc + 12. * ad * bccd +
          24. * a * bccdd + 16. * acd * bcd + 12. * ac * bcdd + 12. * accd * bd + 4. * acc * bdd +
          48. * abcdd * c - 12. * acdd * b * c - 4. * add * bc * c - 8. * ad * bcd * c - 12. * a * bcdd * c -
          8. * acd * bd * c - 4. * ac * bdd * c - 6. * abdd * c * c + 2. * add * b * c * c + 2. * ad * bd * c * c +
          2. * a * bdd * c * c + 6. * abdd * cc - 2. * add * b * cc - 2. * ad * bd * cc - 2. * a * bdd * cc +
          8. * abd * ccd - 4. * ad * b * ccd - 4. * a * bd * ccd + 6. * ab * ccdd - 6. * a * b * ccdd +
          24. * abcd * cd - 8. * acd * b * cd - 4. * ad * bc * cd - 8. * a * bcd * cd - 4. * ac * bd * cd -
          8. * abd * c * cd + 4. * ad * b * c * cd + 4. * a * bd * c * cd - 2. * ab * cd * d + 2. * a * b * cd * d +
          8. * abc * cdd - 4. * ac * b * cdd - 4. * a * bc * cdd - 4. * ab * c * cdd + 4. * a * b * c * cdd +
          48. * abccd * d - 12. * accd * b * d - 8. * acd * bc * d - 4. * ad * bcc * d - 12. * a * bccd * d -
          8. * ac * bcd * d - 4. * acc * bd * d - 24. * abcd * c * d + 8. * acd * b * c * d + 4. * ad * bc * c * d +
          8. * a * bcd * c * d + 4. * ac * bd * c * d + 4. * abd * c * c * d - 2. * ad * b * c * c * d -
          2. * a * bd * c * c * d - 4. * abd * cc * d + 2. * ad * b * cc * d + 2. * a * bd * cc * d - 4. * ab * ccd * d +
          4. * a * b * ccd * d - 8. * abc * cd * d + 4. * ac * b * cd * d + 4. * a * bc * cd * d + 4. * ab * c * cd * d -
          4. * a * b * c * cd * d - 6. * abcc * d * d + 2. * acc * b * d * d + 2. * ac * bc * d * d +
          2. * a * bcc * d * d + 4. * abc * c * d * d - 2. * ac * b * c * d * d - 2. * a * bc * c * d * d -
          ab * c * c * d * d + a * b * c * c * d * d + ab * cc * d * d - a * b * cc * d * d + 6. * abcc * dd -
          2. * acc * b * dd - 2. * ac * bc * dd - 2. * a * bcc * dd - 4. * abc * c * dd + 2. * ac * b * c * dd +
          2. * a * bc * c * dd + ab * c * c * dd - a * b * c * c * dd - ab * cc * dd + a * b * cc * dd)
    .real();
}
template <typename T>
double FlowPtContainer::getStdABCDDD(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> c = inarr[getVectorIndex(0, 0, 1, 0)];
  std::complex<double> d = inarr[getVectorIndex(0, 0, 0, 1)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ac = inarr[getVectorIndex(1, 0, 1, 0)];
  std::complex<double> ad = inarr[getVectorIndex(1, 0, 0, 1)];
  std::complex<double> bc = inarr[getVectorIndex(0, 1, 1, 0)];
  std::complex<double> bd = inarr[getVectorIndex(0, 1, 0, 1)];
  std::complex<double> cd = inarr[getVectorIndex(0, 0, 1, 1)];
  std::complex<double> dd = inarr[getVectorIndex(0, 0, 0, 2)];
  std::complex<double> abc = inarr[getVectorIndex(1, 1, 1, 0)];
  std::complex<double> abd = inarr[getVectorIndex(1, 1, 0, 1)];
  std::complex<double> acd = inarr[getVectorIndex(1, 0, 1, 1)];
  std::complex<double> add = inarr[getVectorIndex(1, 0, 0, 2)];
  std::complex<double> bcd = inarr[getVectorIndex(0, 1, 1, 1)];
  std::complex<double> bdd = inarr[getVectorIndex(0, 1, 0, 2)];
  std::complex<double> cdd = inarr[getVectorIndex(0, 0, 1, 2)];
  std::complex<double> ddd = inarr[getVectorIndex(0, 0, 0, 3)];
  std::complex<double> abcd = inarr[getVectorIndex(1, 1, 1, 1)];
  std::complex<double> abdd = inarr[getVectorIndex(1, 1, 0, 2)];
  std::complex<double> acdd = inarr[getVectorIndex(1, 0, 1, 2)];
  std::complex<double> addd = inarr[getVectorIndex(1, 0, 0, 3)];
  std::complex<double> bcdd = inarr[getVectorIndex(0, 1, 1, 2)];
  std::complex<double> bddd = inarr[getVectorIndex(0, 1, 0, 3)];
  std::complex<double> cddd = inarr[getVectorIndex(0, 0, 1, 3)];
  std::complex<double> abcdd = inarr[getVectorIndex(1, 1, 1, 2)];
  std::complex<double> abddd = inarr[getVectorIndex(1, 1, 0, 3)];
  std::complex<double> acddd = inarr[getVectorIndex(1, 0, 1, 3)];
  std::complex<double> bcddd = inarr[getVectorIndex(0, 1, 1, 3)];
  std::complex<double> abcddd = inarr[getVectorIndex(1, 1, 1, 3)];
  return (-120. * abcddd + 24. * acddd * b + 6. * addd * bc + 12. * add * bcd + 18. * ad * bcdd +
          24. * a * bcddd + 18. * acdd * bd + 12. * acd * bdd + 6. * ac * bddd + 24. * abddd * c -
          6. * addd * b * c - 6. * add * bd * c - 6. * ad * bdd * c - 6. * a * bddd * c + 18. * abdd * cd -
          6. * add * b * cd - 6. * ad * bd * cd - 6. * a * bdd * cd + 12. * abd * cdd - 6. * ad * b * cdd -
          6. * a * bd * cdd + 6. * ab * cddd - 6. * a * b * cddd + 72. * abcdd * d - 18. * acdd * b * d -
          6. * add * bc * d - 12. * ad * bcd * d - 18. * a * bcdd * d - 12. * acd * bd * d - 6. * ac * bdd * d -
          18. * abdd * c * d + 6. * add * b * c * d + 6. * ad * bd * c * d + 6. * a * bdd * c * d -
          12. * abd * cd * d + 6. * ad * b * cd * d + 6. * a * bd * cd * d - 6. * ab * cdd * d + 6. * a * b * cdd * d -
          18. * abcd * d * d + 6. * acd * b * d * d + 3. * ad * bc * d * d + 6. * a * bcd * d * d +
          3. * ac * bd * d * d + 6. * abd * c * d * d - 3. * ad * b * c * d * d - 3. * a * bd * c * d * d +
          3. * ab * cd * d * d - 3. * a * b * cd * d * d + 2. * abc * d * d * d - ac * b * d * d * d - a * bc * d * d * d -
          ab * c * d * d * d + a * b * c * d * d * d + 18. * abcd * dd - 6. * acd * b * dd - 3. * ad * bc * dd -
          6. * a * bcd * dd - 3. * ac * bd * dd - 6. * abd * c * dd + 3. * ad * b * c * dd + 3. * a * bd * c * dd -
          3. * ab * cd * dd + 3. * a * b * cd * dd - 6. * abc * d * dd + 3. * ac * b * d * dd + 3. * a * bc * d * dd +
          3. * ab * c * d * dd - 3. * a * b * c * d * dd + 4. * abc * ddd - 2. * ac * b * ddd - 2. * a * bc * ddd -
          2. * ab * c * ddd + 2. * a * b * c * ddd)
    .real();
}
template <typename T>
double FlowPtContainer::getStdABDDDD(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> d = inarr[getVectorIndex(0, 0, 0, 1)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ad = inarr[getVectorIndex(1, 0, 0, 1)];
  std::complex<double> bd = inarr[getVectorIndex(0, 1, 0, 1)];
  std::complex<double> dd = inarr[getVectorIndex(0, 0, 0, 2)];
  std::complex<double> abd = inarr[getVectorIndex(1, 1, 0, 1)];
  std::complex<double> add = inarr[getVectorIndex(1, 0, 0, 2)];
  std::complex<double> bdd = inarr[getVectorIndex(0, 1, 0, 2)];
  std::complex<double> ddd = inarr[getVectorIndex(0, 0, 0, 3)];
  std::complex<double> abdd = inarr[getVectorIndex(1, 1, 0, 2)];
  std::complex<double> addd = inarr[getVectorIndex(1, 0, 0, 3)];
  std::complex<double> bddd = inarr[getVectorIndex(0, 1, 0, 3)];
  std::complex<double> dddd = inarr[getVectorIndex(0, 0, 0, 4)];
  std::complex<double> abddd = inarr[getVectorIndex(1, 1, 0, 3)];
  std::complex<double> adddd = inarr[getVectorIndex(1, 0, 0, 4)];
  std::complex<double> bdddd = inarr[getVectorIndex(0, 1, 0, 4)];
  std::complex<double> abdddd = inarr[getVectorIndex(1, 1, 0, 4)];
  return (-120. * abdddd + 24. * adddd * b + 24. * addd * bd + 24. * add * bdd + 24. * ad * bddd +
          24. * a * bdddd + 96. * abddd * d - 24. * addd * b * d - 24. * add * bd * d - 24. * ad * bdd * d -
          24. * a * bddd * d - 36. * abdd * d * d + 12. * add * b * d * d + 12. * ad * bd * d * d +
          12. * a * bdd * d * d + 8. * abd * d * d * d - 4. * ad * b * d * d * d - 4. * a * bd * d * d * d - ab * d * d * d * d +
          a * b * d * d * d * d + 36. * abdd * dd - 12. * add * b * dd - 12. * ad * bd * dd - 12. * a * bdd * dd -
          24. * abd * d * dd + 12. * ad * b * d * dd + 12. * a * bd * d * dd + 6. * ab * d * d * dd -
          6. * a * b * d * d * dd - 3. * ab * dd * d + 3. * a * b * dd * d + 16. * abd * ddd - 8. * ad * b * ddd -
          8. * a * bd * ddd - 8. * ab * d * ddd + 8. * a * b * d * ddd + 6. * ab * dddd - 6. * a * b * dddd)
    .real();
}
template <typename T>
double FlowPtContainer::getStdABCCC(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> c = inarr[getVectorIndex(0, 0, 1, 0)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ac = inarr[getVectorIndex(1, 0, 1, 0)];
  std::complex<double> bc = inarr[getVectorIndex(0, 1, 1, 0)];
  std::complex<double> cc = inarr[getVectorIndex(0, 0, 2, 0)];
  std::complex<double> abc = inarr[getVectorIndex(1, 1, 1, 0)];
  std::complex<double> acc = inarr[getVectorIndex(1, 0, 2, 0)];
  std::complex<double> bcc = inarr[getVectorIndex(0, 1, 2, 0)];
  std::complex<double> ccc = inarr[getVectorIndex(0, 0, 3, 0)];
  std::complex<double> abcc = inarr[getVectorIndex(1, 1, 2, 0)];
  std::complex<double> accc = inarr[getVectorIndex(1, 0, 3, 0)];
  std::complex<double> bccc = inarr[getVectorIndex(0, 1, 3, 0)];
  std::complex<double> abccc = inarr[getVectorIndex(1, 1, 3, 0)];
  return (24. * abccc - 6. * accc * b - 6. * acc * bc - 6. * ac * bcc - 6. * a * bccc - 18. * abcc * c +
          6. * acc * b * c + 6. * ac * bc * c + 6. * a * bcc * c + 6. * abc * c * c - 3. * ac * b * c * c -
          3. * a * bc * c * c - ab * c * c * c + a * b * c * c * c - 6. * abc * cc + 3. * ac * b * cc + 3. * a * bc * cc +
          3. * ab * c * cc - 3. * a * b * c * cc - 2. * ab * ccc + 2. * a * b * ccc)
    .real();
}
template <typename T>
double FlowPtContainer::getStdABCCD(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> c = inarr[getVectorIndex(0, 0, 1, 0)];
  std::complex<double> d = inarr[getVectorIndex(0, 0, 0, 1)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ac = inarr[getVectorIndex(1, 0, 1, 0)];
  std::complex<double> ad = inarr[getVectorIndex(1, 0, 0, 1)];
  std::complex<double> bc = inarr[getVectorIndex(0, 1, 1, 0)];
  std::complex<double> bd = inarr[getVectorIndex(0, 1, 0, 1)];
  std::complex<double> cc = inarr[getVectorIndex(0, 0, 2, 0)];
  std::complex<double> cd = inarr[getVectorIndex(0, 0, 1, 1)];
  std::complex<double> abc = inarr[getVectorIndex(1, 1, 1, 0)];
  std::complex<double> abd = inarr[getVectorIndex(1, 1, 0, 1)];
  std::complex<double> acc = inarr[getVectorIndex(1, 0, 2, 0)];
  std::complex<double> bcc = inarr[getVectorIndex(0, 1, 2, 0)];
  std::complex<double> ccd = inarr[getVectorIndex(0, 0, 2, 1)];
  std::complex<double> acd = inarr[getVectorIndex(1, 0, 1, 1)];
  std::complex<double> bcd = inarr[getVectorIndex(0, 1, 1, 1)];
  std::complex<double> abcc = inarr[getVectorIndex(1, 1, 2, 0)];
  std::complex<double> abcd = inarr[getVectorIndex(1, 1, 1, 1)];
  std::complex<double> accd = inarr[getVectorIndex(1, 0, 2, 1)];
  std::complex<double> bccd = inarr[getVectorIndex(0, 1, 2, 1)];
  std::complex<double> abccd = inarr[getVectorIndex(1, 1, 2, 1)];
  return (24. * abccd - 6. * accd * b - 4. * acd * bc - 2. * ad * bcc - 6. * a * bccd - 4. * ac * bcd -
          2. * acc * bd - 12. * abcd * c + 4. * acd * b * c + 2. * ad * bc * c + 4. * a * bcd * c +
          2. * ac * bd * c + 2. * abd * c * c - ad * b * c * c - a * bd * c * c - 2. * abd * cc + ad * b * cc +
          a * bd * cc - 2. * ab * ccd + 2. * a * b * ccd - 4. * abc * cd + 2. * ac * b * cd + 2. * a * bc * cd +
          2. * ab * c * cd - 2. * a * b * c * cd - 6. * abcc * d + 2. * acc * b * d + 2. * ac * bc * d +
          2. * a * bcc * d + 4. * abc * c * d - 2. * ac * b * c * d - 2. * a * bc * c * d - ab * c * c * d +
          a * b * c * c * d + ab * cc * d - a * b * cc * d)
    .real();
}
template <typename T>
double FlowPtContainer::getStdABCDD(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> c = inarr[getVectorIndex(0, 0, 1, 0)];
  std::complex<double> d = inarr[getVectorIndex(0, 0, 0, 1)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ac = inarr[getVectorIndex(1, 0, 1, 0)];
  std::complex<double> ad = inarr[getVectorIndex(1, 0, 0, 1)];
  std::complex<double> bc = inarr[getVectorIndex(0, 1, 1, 0)];
  std::complex<double> bd = inarr[getVectorIndex(0, 1, 0, 1)];
  std::complex<double> cd = inarr[getVectorIndex(0, 0, 1, 1)];
  std::complex<double> dd = inarr[getVectorIndex(0, 0, 0, 2)];
  std::complex<double> abc = inarr[getVectorIndex(1, 1, 1, 0)];
  std::complex<double> abd = inarr[getVectorIndex(1, 1, 0, 1)];
  std::complex<double> add = inarr[getVectorIndex(1, 0, 0, 2)];
  std::complex<double> bdd = inarr[getVectorIndex(0, 1, 0, 2)];
  std::complex<double> cdd = inarr[getVectorIndex(0, 0, 1, 2)];
  std::complex<double> acd = inarr[getVectorIndex(1, 0, 1, 1)];
  std::complex<double> bcd = inarr[getVectorIndex(0, 1, 1, 1)];
  std::complex<double> abdd = inarr[getVectorIndex(1, 1, 0, 2)];
  std::complex<double> abcd = inarr[getVectorIndex(1, 1, 1, 1)];
  std::complex<double> acdd = inarr[getVectorIndex(1, 0, 1, 2)];
  std::complex<double> bcdd = inarr[getVectorIndex(0, 1, 1, 2)];
  std::complex<double> abcdd = inarr[getVectorIndex(1, 1, 1, 2)];
  return (24. * abcdd - 6. * acdd * b - 2. * add * bc - 4. * ad * bcd - 6. * a * bcdd - 4. * acd * bd -
          2. * ac * bdd - 6. * abdd * c + 2. * add * b * c + 2. * ad * bd * c + 2. * a * bdd * c - 4. * abd * cd +
          2. * ad * b * cd + 2. * a * bd * cd - 2. * ab * cdd + 2. * a * b * cdd - 12. * abcd * d +
          4. * acd * b * d + 2. * ad * bc * d + 4. * a * bcd * d + 2. * ac * bd * d + 4. * abd * c * d -
          2. * ad * b * c * d - 2. * a * bd * c * d + 2. * ab * cd * d - 2. * a * b * cd * d + 2. * abc * d * d -
          ac * b * d * d - a * bc * d * d - ab * c * d * d + a * b * c * d * d - 2. * abc * dd + ac * b * dd +
          a * bc * dd + ab * c * dd - a * b * c * dd)
    .real();
}
template <typename T>
double FlowPtContainer::getStdABDDD(T& inarr)
{
  std::complex<double> a = inarr[getVectorIndex(1, 0, 0, 0)];
  std::complex<double> b = inarr[getVectorIndex(0, 1, 0, 0)];
  std::complex<double> d = inarr[getVectorIndex(0, 0, 0, 1)];
  std::complex<double> ab = inarr[getVectorIndex(1, 1, 0, 0)];
  std::complex<double> ad = inarr[getVectorIndex(1, 0, 0, 1)];
  std::complex<double> bd = inarr[getVectorIndex(0, 1, 0, 1)];
  std::complex<double> dd = inarr[getVectorIndex(0, 0, 0, 2)];
  std::complex<double> abd = inarr[getVectorIndex(1, 1, 0, 1)];
  std::complex<double> add = inarr[getVectorIndex(1, 0, 0, 2)];
  std::complex<double> bdd = inarr[getVectorIndex(0, 1, 0, 2)];
  std::complex<double> ddd = inarr[getVectorIndex(0, 0, 0, 3)];
  std::complex<double> abdd = inarr[getVectorIndex(1, 1, 0, 2)];
  std::complex<double> addd = inarr[getVectorIndex(1, 0, 0, 3)];
  std::complex<double> bddd = inarr[getVectorIndex(0, 1, 0, 3)];
  std::complex<double> abddd = inarr[getVectorIndex(1, 1, 0, 3)];
  return (24. * abddd - 6. * addd * b - 6. * add * bd - 6. * ad * bdd - 6. * a * bddd - 18. * abdd * d +
          6. * add * b * d + 6. * ad * bd * d + 6. * a * bdd * d + 6. * abd * d * d - 3. * ad * b * d * d -
          3. * a * bd * d * d - ab * d * d * d + a * b * d * d * d - 6. * abd * dd + 3. * ad * b * dd + 3. * a * bd * dd +
          3. * ab * d * dd - 3. * a * b * d * dd - 2. * ab * ddd + 2. * a * b * ddd)
    .real();
}
double FlowPtContainer::orderedAddition(std::vector<double> vec)
{
  double sum = 0;
  std::sort(vec.begin(), vec.end());
  for (size_t i = 0; i < vec.size(); i++) {
    sum += vec[i];
  }
  return sum;
}
void FlowPtContainer::rebinMulti(int nbins)
{
  clearPtPostprocessLists();
  if (fCMTermList) {
    for (int i = 0; i < fCMTermList->GetEntries(); i++) {
      dynamic_cast<BootstrapProfile*>(fCMTermList->At(i))->RebinMulti(nbins);
    }
  }
  if (fCorrList) {
    for (int i = 0; i < fCorrList->GetEntries(); i++) {
      dynamic_cast<BootstrapProfile*>(fCorrList->At(i))->RebinMulti(nbins);
    }
  }
  if (fCovList) {
    for (int i = 0; i < fCovList->GetEntries(); i++) {
      dynamic_cast<BootstrapProfile*>(fCovList->At(i))->RebinMulti(nbins);
    }
  }
  if (fSubList) {
    for (int i = 0; i < fSubList->GetEntries(); i++) {
      dynamic_cast<BootstrapProfile*>(fSubList->At(i))->RebinMulti(nbins);
    }
  }
  if (fSubCMList) {
    for (int i = 0; i < fSubCMList->GetEntries(); i++) {
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(i))->RebinMulti(nbins);
    }
  }
}
void FlowPtContainer::rebinMulti(int nbins, double* binedges)
{
  clearPtPostprocessLists();
  if (fCMTermList) {
    for (int i = 0; i < fCMTermList->GetEntries(); i++) {
      dynamic_cast<BootstrapProfile*>(fCMTermList->At(i))->RebinMulti(nbins, binedges);
    }
  }
  if (fCorrList) {
    for (int i = 0; i < fCorrList->GetEntries(); i++) {
      dynamic_cast<BootstrapProfile*>(fCorrList->At(i))->RebinMulti(nbins, binedges);
    }
  }
  if (fCovList) {
    for (int i = 0; i < fCovList->GetEntries(); i++) {
      dynamic_cast<BootstrapProfile*>(fCovList->At(i))->RebinMulti(nbins, binedges);
    }
  }
  if (fSubList) {
    for (int i = 0; i < fSubList->GetEntries(); i++) {
      dynamic_cast<BootstrapProfile*>(fSubList->At(i))->RebinMulti(nbins, binedges);
    }
  }
  if (fSubCMList) {
    for (int i = 0; i < fSubCMList->GetEntries(); i++) {
      dynamic_cast<BootstrapProfile*>(fSubCMList->At(i))->RebinMulti(nbins, binedges);
    }
  }
}
TH1* FlowPtContainer::getCorrHist(int ind, int m)
{
  return dynamic_cast<BootstrapProfile*>(fCorrList->FindObject(Form("mpt%i", m)))->getHist(ind);
}
TH1* FlowPtContainer::getPtCentralMomentHist(const char* head, int ind)
{
  if (!head || head[0] == 0 || !fCorrList) {
    return nullptr;
  }

  const std::string resultName = makePostprocessName("ptcentral", head, ind);
  if (fPtCentralMomentList) {
    if (auto* result = dynamic_cast<TH1*>(fPtCentralMomentList->FindObject(resultName.c_str()))) {
      return result;
    }
  }

  const auto config = std::find_if(fPtMomentConfigs.begin(), fPtMomentConfigs.end(), [&](const auto& candidate) { return candidate.head == head; });
  if (config == fPtMomentConfigs.end()) {
    LOGF(error, "Cannot calculate pT central moment %s: no pure-pT postprocessing configuration was stored", head);
    return nullptr;
  }

  int order = 0;
  for (const auto& group : config->groups) {
    order += group.order;
  }
  if (order < 1 || order >= std::numeric_limits<std::uint64_t>::digits) {
    LOGF(error, "Cannot calculate pT central moment %s: invalid order %i", head, order);
    return nullptr;
  }

  auto getProfileHist = [&](const std::string& name) -> std::unique_ptr<TH1> {
    auto* profile = dynamic_cast<BootstrapProfile*>(fCorrList->FindObject(name.c_str()));
    if (!profile) {
      LOGF(error, "Cannot calculate pT central moment %s: profile %s is missing", head, name.c_str());
      return nullptr;
    }
    std::unique_ptr<TH1> histogram{profile->getHist(ind)};
    if (histogram) {
      histogram->SetDirectory(nullptr);
    }
    return histogram;
  };

  auto shape = getProfileHist(config->head);
  if (!shape) {
    return nullptr;
  }
  auto* result = dynamic_cast<TH1*>(shape->Clone(resultName.c_str()));
  result->SetDirectory(nullptr);
  result->Reset();

  const std::uint64_t allFactors = (std::uint64_t{1} << order) - 1;
  std::vector<int> activeCounts(config->groups.size(), 0);
  bool finished = false;
  while (!finished) {
    int activeOrder = 0;
    int factorOffset = 0;
    double coefficient = 1.;
    std::uint64_t activeFactors = 0;
    std::map<std::string, int> meanPowers;
    for (size_t groupIndex = 0; groupIndex < config->groups.size(); ++groupIndex) {
      const auto& group = config->groups[groupIndex];
      const int active = activeCounts[groupIndex];
      activeOrder += active;
      coefficient *= binomialCoefficient(group.order, active);
      for (int factor = 0; factor < active; ++factor) {
        activeFactors |= std::uint64_t{1} << (factorOffset + factor);
      }
      if (active < group.order) {
        meanPowers[group.meanProfile] += group.order - active;
      }
      factorOffset += group.order;
    }

    std::unique_ptr<TH1> term;
    if (activeOrder == 0) {
      term.reset(dynamic_cast<TH1*>(shape->Clone()));
      term->SetDirectory(nullptr);
      term->Reset();
      for (int bin = 1; bin <= term->GetNbinsX(); ++bin) {
        if (shape->GetBinContent(bin) != 0. || shape->GetBinError(bin) != 0.) {
          term->SetBinContent(bin, 1.);
        }
      }
    } else {
      term = getProfileHist(makeFullTopologyName(config->head, activeFactors, allFactors, false, order));
      if (!term) {
        delete result;
        return nullptr;
      }
      term->SetDirectory(nullptr);
    }

    for (const auto& [meanName, power] : meanPowers) {
      auto mean = getProfileHist(meanName);
      if (!mean) {
        delete result;
        return nullptr;
      }
      mean->SetDirectory(nullptr);
      if (power > 1) {
        for (int bin = 0; bin <= mean->GetNbinsX() + 1; ++bin) {
          const double value = mean->GetBinContent(bin);
          const double error = mean->GetBinError(bin);
          mean->SetBinContent(bin, std::pow(value, power));
          mean->SetBinError(bin, std::abs(power * std::pow(value, power - 1) * error));
        }
      }
      term->Multiply(mean.get());
    }

    if ((order - activeOrder) % ParityDivisor != 0) {
      coefficient = -coefficient;
    }
    result->Add(term.get(), coefficient);

    size_t groupIndex = 0;
    for (; groupIndex < activeCounts.size(); ++groupIndex) {
      if (++activeCounts[groupIndex] <= config->groups[groupIndex].order) {
        break;
      }
      activeCounts[groupIndex] = 0;
    }
    finished = groupIndex == activeCounts.size();
  }

  if (!fPtCentralMomentList) {
    fPtCentralMomentList = new TList();
    fPtCentralMomentList->SetOwner(kTRUE);
  }
  fPtCentralMomentList->Add(result);
  return result;
}
TH1* FlowPtContainer::getPtCumulantHist(const char* head, int ind)
{
  if (!head || head[0] == 0 || !fCorrList) {
    return nullptr;
  }

  const std::string resultName = makePostprocessName("ptcumulant", head, ind);
  if (fPtCumulantList) {
    if (auto* result = dynamic_cast<TH1*>(fPtCumulantList->FindObject(resultName.c_str()))) {
      return result;
    }
  }

  const auto config = std::find_if(fPtMomentConfigs.begin(), fPtMomentConfigs.end(), [&](const auto& candidate) { return candidate.head == head; });
  if (config == fPtMomentConfigs.end()) {
    LOGF(error, "Cannot calculate pT cumulant %s: no pure-pT postprocessing configuration was stored", head);
    return nullptr;
  }

  std::vector<size_t> factorGroups;
  for (size_t groupIndex = 0; groupIndex < config->groups.size(); ++groupIndex) {
    factorGroups.insert(factorGroups.end(), config->groups[groupIndex].order, groupIndex);
  }
  const int order = factorGroups.size();
  if (order < 1) {
    return nullptr;
  }

  if (order == 1) {
    auto* profile = dynamic_cast<BootstrapProfile*>(fCorrList->FindObject(config->groups.front().meanProfile.c_str()));
    if (!profile) {
      LOGF(error, "Cannot calculate pT cumulant %s: profile %s is missing", head, config->groups.front().meanProfile.c_str());
      return nullptr;
    }
    auto* result = profile->getHist(ind);
    if (!result) {
      return nullptr;
    }
    result->SetName(resultName.c_str());
    result->SetDirectory(nullptr);
    if (!fPtCumulantList) {
      fPtCumulantList = new TList();
      fPtCumulantList->SetOwner(kTRUE);
    }
    fPtCumulantList->Add(result);
    return result;
  }

  auto* centralMoment = getPtCentralMomentHist(head, ind);
  if (!centralMoment) {
    return nullptr;
  }
  auto* result = dynamic_cast<TH1*>(centralMoment->Clone(resultName.c_str()));
  result->SetDirectory(nullptr);
  result->Reset();

  using Topology = std::vector<std::pair<std::string, int>>;
  auto topology = [](const PtMomentConfig& candidate, const std::vector<int>* counts = nullptr) {
    Topology value;
    for (size_t groupIndex = 0; groupIndex < candidate.groups.size(); ++groupIndex) {
      const int count = counts ? counts->at(groupIndex) : candidate.groups[groupIndex].order;
      if (count > 0) {
        value.emplace_back(candidate.groups[groupIndex].region, count);
      }
    }
    std::sort(value.begin(), value.end());
    return value;
  };
  auto findMomentConfig = [&](const std::vector<int>& counts) {
    const auto required = topology(*config, &counts);
    return std::find_if(fPtMomentConfigs.begin(), fPtMomentConfigs.end(), [&](const auto& candidate) { return topology(candidate) == required; });
  };
  auto describeTopology = [&](const std::vector<int>& counts) {
    std::string description;
    for (const auto& [region, count] : topology(*config, &counts)) {
      if (!description.empty()) {
        description += " ";
      }
      description += region + "{" + std::to_string(count) + "pt}";
    }
    return description;
  };

  bool failed = false;
  std::vector<int> blockForFactor(order, 0);
  std::function<void(int, int)> addPartitions = [&](int factor, int highestBlock) {
    if (failed) {
      return;
    }
    if (factor < order) {
      for (int block = 0; block <= highestBlock + 1; ++block) {
        blockForFactor[factor] = block;
        addPartitions(factor + 1, std::max(highestBlock, block));
      }
      return;
    }

    const int numberOfBlocks = highestBlock + 1;
    std::vector<std::vector<int>> blockCounts(numberOfBlocks, std::vector<int>(config->groups.size(), 0));
    std::vector<int> blockSizes(numberOfBlocks, 0);
    for (int factorIndex = 0; factorIndex < order; ++factorIndex) {
      const int block = blockForFactor[factorIndex];
      ++blockCounts[block][factorGroups[factorIndex]];
      ++blockSizes[block];
    }
    if (std::any_of(blockSizes.begin(), blockSizes.end(), [](int size) { return size == 1; })) {
      return;
    }

    std::unique_ptr<TH1> term;
    for (int block = 0; block < numberOfBlocks; ++block) {
      const PtMomentConfig* blockConfig = &*config;
      if (numberOfBlocks > 1) {
        const auto found = findMomentConfig(blockCounts[block]);
        if (found == fPtMomentConfigs.end()) {
          LOGF(error, "Cannot calculate pT cumulant %s: missing lower-order correlation %s", head, describeTopology(blockCounts[block]).c_str());
          failed = true;
          return;
        }
        blockConfig = &*found;
      }
      auto* blockMoment = getPtCentralMomentHist(blockConfig->head.c_str(), ind);
      if (!blockMoment) {
        failed = true;
        return;
      }
      if (!term) {
        term.reset(dynamic_cast<TH1*>(blockMoment->Clone()));
        term->SetDirectory(nullptr);
      } else {
        term->Multiply(blockMoment);
      }
    }

    double coefficient = 1.;
    for (int i = 2; i < numberOfBlocks; ++i) {
      coefficient *= i;
    }
    if (numberOfBlocks % ParityDivisor == 0) {
      coefficient = -coefficient;
    }
    result->Add(term.get(), coefficient);
  };

  blockForFactor.front() = 0;
  addPartitions(1, 0);
  if (failed) {
    delete result;
    return nullptr;
  }

  if (!fPtCumulantList) {
    fPtCumulantList = new TList();
    fPtCumulantList->SetOwner(kTRUE);
  }
  fPtCumulantList->Add(result);
  return result;
}
TH1* FlowPtContainer::getCentralMomentHist(int ind, int m)
{
  if (!fCentralMomentList) {
    createCentralMomentList();
  }
  if (!fCentralMomentList) {
    return nullptr;
  }
  if (ind + 1 < fCentralMomentList->GetEntries()) {
    return dynamic_cast<TH1*>(fCentralMomentList->FindObject(Form("cm%i_%i", m, ind)));
  }
  return nullptr;
}
void FlowPtContainer::createCentralMomentList()
{
  delete fCentralMomentList;
  fCentralMomentList = new TList();
  fCentralMomentList->SetOwner();
  for (auto m(1); m <= centralMomentMaxOrder; ++m) {
    for (int i = -1; i < dynamic_cast<BootstrapProfile*>(fCMTermList->At(0))->getNSubs(); ++i) {
      TH1* hMpt = dynamic_cast<BootstrapProfile*>(fCMTermList->At(0))->getHist(i);
      std::vector<TH1*> hTs;
      for (int j = 0; j < m; ++j) {
        dynamic_cast<BootstrapProfile*>(fCMTermList->FindObject(Form("cm%i_Mpt%i", m, j)))->SetErrorOption("g");
        hTs.push_back(dynamic_cast<BootstrapProfile*>(fCMTermList->FindObject(Form("cm%i_Mpt%i", m, j)))->getHist(i));
      }
      calculateCentralMomentHists(hTs, i, m, hMpt);
    }
  }
}
void FlowPtContainer::calculateCentralMomentHists(std::vector<TH1*> inh, int ind, int m, TH1* hMpt)
{
  TH1* reth = dynamic_cast<TH1*>(inh[0]->Clone(Form("cm%i_%i", m, ind)));
  for (auto i(1); i < m; ++i) {
    TH1* mptPow = raiseHistToPower(hMpt, i);
    inh[i]->Multiply(mptPow);
    reth->Add(inh[i]);
  }
  TH1* mptLast = raiseHistToPower(hMpt, m);
  reth->Add(mptLast, ((m % 2) != 0) ? (-1) : 1);
  fCentralMomentList->Add(reth);
}
TH1* FlowPtContainer::getCumulantHist(int ind, int m)
{
  if (!fCumulantList) {
    createCumulantList();
  }
  if (!fCumulantList) {
    return nullptr;
  }
  if (ind + 1 < fCumulantList->GetEntries()) {
    return dynamic_cast<TH1*>(fCumulantList->At((ind + 1) * mpar + m - 1));
  }
  return nullptr;
}
void FlowPtContainer::createCumulantList()
{
  delete fCumulantList;
  fCumulantList = new TList();
  fCumulantList->SetOwner();
  //((BootstrapProfile*)fCorrList->At(0))->PresetWeights((BootstrapProfile*)fCorrList->At(mpar-1));
  for (int i = -1; i < dynamic_cast<BootstrapProfile*>(fCorrList->At(0))->getNSubs(); ++i) {
    std::vector<TH1*> hTs;
    for (int j = 0; j < mpar; ++j) {
      dynamic_cast<BootstrapProfile*>(fCorrList->FindObject(Form("mpt%i", j + 1)))->SetErrorOption("g");
      hTs.push_back(dynamic_cast<BootstrapProfile*>(fCorrList->FindObject(Form("mpt%i", j + 1)))->getHist(i));
    }
    calculateCumulantHists(hTs, i);
  }
  //((BootstrapProfile*)fCorrList->At(0))->PresetWeights(0);
}
void FlowPtContainer::calculateCumulantHists(std::vector<TH1*> inh, int ind)
{
  auto binomial = [&](const int n, const int m) { assert(n >= m); return FactorialArray[n]/(FactorialArray[m]*FactorialArray[n-m]); };
  for (int m = 1; m <= mpar; ++m) {
    TH1* reth = dynamic_cast<TH1*>(inh[m - 1]->Clone(Form("reth%i_%i", m, ind)));
    // TH1* hWeights = (TH1*)inh[m-1]->Clone(Form("hWeights%i_%i",m,ind));
    for (int k = 1; k < m; ++k) {
      TH1* corrh = dynamic_cast<TH1*>(inh[m - k - 1]->Clone(Form("hcorr%i%i_%i", m, k, ind)));
      corrh->Multiply(dynamic_cast<TH1*>(fCumulantList->At((ind + 1) * mpar + k - 1)));
      corrh->Scale(binomial(m - 1, k - 1));
      reth->Add(corrh, -1);
      delete corrh;
    }
    // for(int i=1;i<=hWeights->GetNbinsX();++i) reth->SetBinError(i,hWeights->GetBinError(i));
    // delete hWeights;
    fCumulantList->Add(dynamic_cast<TH1*>(reth->Clone(Form("kappa%i_%i", m, ind))));
  }
}
Long64_t FlowPtContainer::Merge(TCollection* collist) // o2-linter: disable=name/function-variable (Keep name consistent with ROOT streamer merging)
{
  delete fCumulantList;
  fCumulantList = nullptr;
  delete fCentralMomentList;
  fCentralMomentList = nullptr;
  clearPtPostprocessLists();

  Long64_t nmerged = 0;
  TIter allPTC(collist);
  FlowPtContainer* lPTC = nullptr;
  while ((lPTC = dynamic_cast<FlowPtContainer*>(allPTC()))) {
    if (!lPTC->fPtMomentConfigs.empty()) {
      if (fPtMomentConfigs.empty()) {
        fPtMomentConfigs = lPTC->fPtMomentConfigs;
      } else if (!sameMomentConfigs(fPtMomentConfigs, lPTC->fPtMomentConfigs)) {
        LOGF(error, "Cannot merge FlowPtContainer %s: pT postprocessing configurations differ", GetName());
        continue;
      }
    } else if (!fPtMomentConfigs.empty()) {
      LOGF(error, "Cannot merge FlowPtContainer %s: incoming pT postprocessing configuration is missing", GetName());
      continue;
    }

    TList* tCMTerm = lPTC->fCMTermList;
    TList* tCorr = lPTC->fCorrList;
    TList* tCov = lPTC->fCovList;
    TList* tSub = lPTC->fSubList;
    TList* tSubCM = lPTC->fSubCMList;
    bool merged = false;
    if (tCMTerm) {
      if (!fCMTermList) {
        fCMTermList = dynamic_cast<TList*>(tCMTerm->Clone());
      } else {
        mergeBSLists(fCMTermList, tCMTerm);
      }
      merged = true;
    }
    if (tCorr) {
      if (!fCorrList) {
        fCorrList = dynamic_cast<TList*>(tCorr->Clone());
      } else {
        mergeBSLists(fCorrList, tCorr);
      }
      merged = true;
    }
    if (tCov) {
      if (!fCovList) {
        fCovList = dynamic_cast<TList*>(tCov->Clone());
      } else {
        mergeBSLists(fCovList, tCov);
      }
      merged = true;
    }
    if (tSub) {
      if (!fSubList) {
        fSubList = dynamic_cast<TList*>(tSub->Clone());
      } else {
        mergeBSLists(fSubList, tSub);
      }
      merged = true;
    }
    if (tSubCM) {
      if (!fSubCMList) {
        fSubCMList = dynamic_cast<TList*>(tSubCM->Clone());
      } else {
        mergeBSLists(fSubCMList, tSubCM);
      }
      merged = true;
    }
    if (merged) {
      ++nmerged;
    }
  }
  return nmerged;
}
void FlowPtContainer::mergeBSLists(TList* source, TList* target)
{
  if (source->GetEntries() != target->GetEntries()) {
    LOGF(warning, "Number in lists to be merged are not the same, skipping...\n");
    return;
  }
  for (int i = 0; i < source->GetEntries(); i++) {
    auto* lObj = dynamic_cast<BootstrapProfile*>(source->At(i));
    auto* tObj = dynamic_cast<BootstrapProfile*>(target->At(i));
    lObj->MergeBS(tObj);
  }
}
TH1* FlowPtContainer::raiseHistToPower(TH1* inh, double p)
{
  TH1D* reth = dynamic_cast<TH1D*>(inh->Clone("reth"));
  reth->SetName(Form("power%.2f_%s", p, inh->GetName()));
  for (int i = 1; i <= inh->GetNbinsX(); i++) {
    if (inh->GetBinContent(i) >= 0 || std::floor(p) == p) {
      reth->SetBinContent(i, std::pow(inh->GetBinContent(i), p));
      reth->SetBinError(i, p * std::pow(reth->GetBinContent(i), p - 1) * inh->GetBinError(i));
    } else {
      reth->SetBinContent(i, -999);
      reth->SetBinError(i, 0.000000001);
    }
  }
  return reth;
}
void FlowPtContainer::getSubevents(int k, int n, std::vector<int>& current, std::vector<std::vector<int>>& outputSubevents)
{
  if (n == 1) {
    // Last box gets all remaining objects
    current.push_back(k);
    outputSubevents.push_back(current);
    current.pop_back();
    return;
  }

  for (int i = 0; i <= k; ++i) {
    current.push_back(i);
    getSubevents(k - i, n - 1, current, outputSubevents);
    current.pop_back();
  }
}
