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
#ifndef O2_ANALYSIS_TRIGGER_H_
#define O2_ANALYSIS_TRIGGER_H_

#include <array>
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace filtering
{
DECLARE_SOA_COLUMN(H2, hasH2, bool);   //!
DECLARE_SOA_COLUMN(H3, hasH3, bool);   //!
DECLARE_SOA_COLUMN(He3, hasHe3, bool); //!
DECLARE_SOA_COLUMN(He4, hasHe4, bool); //!

// diffraction
DECLARE_SOA_COLUMN(TwoPi, has2pi, bool);  //! Double Gap events, DG, 2 pion
DECLARE_SOA_COLUMN(FourPi, has4pi, bool); //! 4 pion
DECLARE_SOA_COLUMN(TwoK, has2K, bool);    //! 2 K
DECLARE_SOA_COLUMN(FourK, has4K, bool);   //! 4 K

// Dileptons & Quarkonia
DECLARE_SOA_COLUMN(SingleE, hasSingleE, bool);           //! single electron trigger
DECLARE_SOA_COLUMN(SingleMuLow, hasSingleMuLow, bool);   //! single muon with low pT trigger
DECLARE_SOA_COLUMN(SingleMuHigh, hasSingleMuHigh, bool); //! single muon with high pT trigger
DECLARE_SOA_COLUMN(DiElectron, hasDiElectron, bool);     //! dielectron trigger
DECLARE_SOA_COLUMN(DiMuon, hasDiMuon, bool);             //! dimuon trigger with low pT on muons

// heavy flavours
DECLARE_SOA_COLUMN(HfHighPt, hasHfHighPt, bool);           //! high-pT charm hadron
DECLARE_SOA_COLUMN(HfBeauty, hasHfBeauty, bool);           //! beauty hadron
DECLARE_SOA_COLUMN(HfFemto, hasHfFemto, bool);             //! charm-hadron - N pair
DECLARE_SOA_COLUMN(HfDoubleCharm, hasHfDoubleCharm, bool); //! at least two charm-hadron candidates

// CF three body triggers
DECLARE_SOA_COLUMN(PPP, hasPPP, bool); //! has p-p-p triplet
DECLARE_SOA_COLUMN(PPL, hasPPL, bool); //! has p-p-L triplet
DECLARE_SOA_COLUMN(PLL, hasPLL, bool); //! has p-L-L triplet
DECLARE_SOA_COLUMN(LLL, hasLLL, bool); //! has L-L-L tripletD

// jets
DECLARE_SOA_COLUMN(HighPtTrack, hasHighPtTrack, bool); //! event contains high-pT track
DECLARE_SOA_COLUMN(JetChHighPt, hasJetChHighPt, bool); //! high-pT charged jet

// strangeness (lf)
DECLARE_SOA_COLUMN(Omega, hasOmega, bool);             //! at leat 1 Omega
DECLARE_SOA_COLUMN(hadronXi, hashadronXi, bool);       //! at least 1 Xi + high-pt hadron
DECLARE_SOA_COLUMN(DoubleXi, hasDoubleXi, bool);       //! at least 2 Xi
DECLARE_SOA_COLUMN(TripleXi, hasTripleXi, bool);       //! at least 3 Xi
DECLARE_SOA_COLUMN(QuadrupleXi, hasQuadrupleXi, bool); //! at least 4 Xi
DECLARE_SOA_COLUMN(SingleXiYN, hasSingleXiYN, bool);   //! at least 1 Xi with R > 24.39 cm (YN interactions)
// multiplicity
DECLARE_SOA_COLUMN(LeadingPtTrack, hasLeadingPtTrack, bool);         //! event contains leading track
DECLARE_SOA_COLUMN(HighMultFv0, hasHighMultFv0, bool);               //! high FV0 muliplicity
DECLARE_SOA_COLUMN(HighTrackMult, hasHighTrackMult, bool);           //! high muliplicity collision

} // namespace filtering

namespace decision
{

DECLARE_SOA_COLUMN(BCId, hasBCId, int);                     //! Bunch crossing Id
DECLARE_SOA_COLUMN(CollisionTime, hasCollisionTime, float); //! Collision time
DECLARE_SOA_COLUMN(CefpSelected, hasCefpSelected, bool);    //! CEFP decision

} // namespace decision

// nuclei
DECLARE_SOA_TABLE(NucleiFilters, "AOD", "NucleiFilters", //!
                  filtering::H2, filtering::H3, filtering::He3, filtering::He4);
using NucleiFilter = NucleiFilters::iterator;

// diffraction
DECLARE_SOA_TABLE(DiffractionFilters, "AOD", "DiffFilters", //! Diffraction filters
                  filtering::TwoPi, filtering::FourPi, filtering::TwoK, filtering::FourK);
using DiffractionFilter = DiffractionFilters::iterator;

// Dileptons & Quarkonia
DECLARE_SOA_TABLE(DqFilters, "AOD", "DqFilters", //!
                  filtering::SingleE, filtering::DiElectron, filtering::SingleMuLow, filtering::SingleMuHigh, filtering::DiMuon);
using DqFilter = DqFilters::iterator;

// heavy flavours
DECLARE_SOA_TABLE(HfFilters, "AOD", "HfFilters", //!
                  filtering::HfHighPt, filtering::HfBeauty, filtering::HfFemto, filtering::HfDoubleCharm);

using HfFilter = HfFilters::iterator;

// correlations
DECLARE_SOA_TABLE(CFFilters, "AOD", "CFFilters", //!
                  filtering::PPP, filtering::PPL, filtering::PLL, filtering::LLL);
using CfFilter = CFFilters::iterator;

// jets
DECLARE_SOA_TABLE(JetFilters, "AOD", "JetFilters", //!
                  filtering::HighPtTrack, filtering::JetChHighPt);

using JetFilter = JetFilters::iterator;

// strangeness (lf)
DECLARE_SOA_TABLE(StrangenessFilters, "AOD", "LFStrgFilters", //!
                  filtering::Omega, filtering::hadronXi, filtering::DoubleXi, filtering::TripleXi, filtering::QuadrupleXi, filtering::SingleXiYN);

using StrangenessFilter = StrangenessFilters::iterator;

// multiplicity
DECLARE_SOA_TABLE(MultFilters, "AOD", "MultFilters", //!
                  filtering::LeadingPtTrack, filtering::HighMultFv0, filtering::HighTrackMult);
using MultFilter = MultFilters::iterator;

// cefp decision
DECLARE_SOA_TABLE(CefpDecisions, "AOD", "CefpDecision", //!
                  decision::BCId, decision::CollisionTime, decision::CefpSelected);
using CefpDecision = CefpDecisions::iterator;

/// List of the available filters, the description of their tables and the name of the tasks
constexpr int NumberOfFilters{8};
constexpr std::array<char[32], NumberOfFilters> AvailableFilters{"NucleiFilters", "DiffractionFilters", "DqFilters", "HfFilters", "CFFilters", "JetFilters", "StrangenessFilters", "MultFilters"};
constexpr std::array<char[16], NumberOfFilters> FilterDescriptions{"NucleiFilters", "DiffFilters", "DQFilters", "HFFilters", "CFFilters", "JetFilters", "LFStrgFilters", "MultFilters"};
constexpr std::array<char[128], NumberOfFilters> FilteringTaskNames{"o2-analysis-nuclei-filter", "o2-analysis-diffraction-filter", "o2-analysis-dq-filter-pp", "o2-analysis-hf-filter", "o2-analysis-cf-threebodyfemto-filter", "o2-analysis-je-filter", "o2-analysis-lf-strangeness-filter", "o2-analysis-mult-filter"};
constexpr o2::framework::pack<NucleiFilters, DiffractionFilters, DqFilters, HfFilters, CFFilters, JetFilters, StrangenessFilters, MultFilters> FiltersPack;
static_assert(o2::framework::pack_size(FiltersPack) == NumberOfFilters);

template <typename T, typename C>
void addColumnToMap(std::unordered_map<std::string, std::unordered_map<std::string, float>>& map)
{
  map[MetadataTrait<T>::metadata::tableLabel()][C::columnLabel()] = 1.f;
}

template <typename T, typename... C>
void addColumnsToMap(o2::framework::pack<C...>, std::unordered_map<std::string, std::unordered_map<std::string, float>>& map)
{
  (addColumnToMap<T, C>(map), ...);
}

template <typename... T>
void FillFiltersMap(o2::framework::pack<T...>, std::unordered_map<std::string, std::unordered_map<std::string, float>>& map)
{
  (addColumnsToMap<T>(typename T::iterator::persistent_columns_t{}, map), ...);
}

template <typename... C>
static std::vector<std::string> ColumnsNames(o2::framework::pack<C...>)
{
  return {C::columnLabel()...};
}

template <typename T>
unsigned int NumberOfColumns()
{
  return o2::framework::pack_size(typename T::iterator::persistent_columns_t{});
}

} // namespace o2::aod

#endif // O2_ANALYSIS_TRIGGER_H_
