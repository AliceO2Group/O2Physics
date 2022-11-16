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
#ifndef EVENTFILTERING_FILTERTABLES_H_
#define EVENTFILTERING_FILTERTABLES_H_

#include <array>
#include <unordered_map>
#include <string>
#include <vector>
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
DECLARE_SOA_COLUMN(HfHighPt2P, hasHfHighPt2P, bool);             //! high-pT 2-prong charm hadron
DECLARE_SOA_COLUMN(HfHighPt3P, hasHfHighPt3P, bool);             //! high-pT 3-prong charm hadron
DECLARE_SOA_COLUMN(HfBeauty3P, hasHfBeauty3P, bool);             //! 3-prong beauty hadron
DECLARE_SOA_COLUMN(HfBeauty4P, hasHfBeauty4P, bool);             //! 4-prong beauty hadron
DECLARE_SOA_COLUMN(HfFemto2P, hasHfFemto2P, bool);               //! 2-prong charm-hadron - N pair
DECLARE_SOA_COLUMN(HfFemto3P, hasHfFemto3P, bool);               //! 3-prong charm-hadron - N pair
DECLARE_SOA_COLUMN(HfDoubleCharm2P, hasHfDoubleCharm2P, bool);   //! at least two 2-prong charm-hadron candidates
DECLARE_SOA_COLUMN(HfDoubleCharm3P, hasHfDoubleCharm3P, bool);   //! at least two 3-prong charm-hadron candidates
DECLARE_SOA_COLUMN(HfDoubleCharmMix, hasHfDoubleCharmMix, bool); //! at least one 2-prong and one 3-prong charm-hadron candidates

// CF two body triggers
DECLARE_SOA_COLUMN(PD, hasPD, bool); //! has d-p pair
DECLARE_SOA_COLUMN(LD, hasLD, bool); //! has l-d pair

// CF three body triggers
DECLARE_SOA_COLUMN(PPP, hasPPP, bool); //! has p-p-p triplet
DECLARE_SOA_COLUMN(PPL, hasPPL, bool); //! has p-p-L triplet
DECLARE_SOA_COLUMN(PLL, hasPLL, bool); //! has p-L-L triplet
DECLARE_SOA_COLUMN(LLL, hasLLL, bool); //! has L-L-L tripletD

// jets
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(JetChHighPt, hasJetChHighPt, bool); //! high-pT charged jet

// strangeness (lf)
DECLARE_SOA_COLUMN(Omega, hasOmega, bool);             //! at leat 1 Omega
DECLARE_SOA_COLUMN(hadronXi, hashadronXi, bool);       //! at least 1 Xi + high-pt hadron
DECLARE_SOA_COLUMN(DoubleXi, hasDoubleXi, bool);       //! at least 2 Xi
DECLARE_SOA_COLUMN(TripleXi, hasTripleXi, bool);       //! at least 3 Xi
DECLARE_SOA_COLUMN(QuadrupleXi, hasQuadrupleXi, bool); //! at least 4 Xi
DECLARE_SOA_COLUMN(SingleXiYN, hasSingleXiYN, bool);   //! at least 1 Xi with R > 24.39 cm (YN interactions)

// multiplicity
DECLARE_SOA_COLUMN(HighTrackMult, hasHighTrackMult, bool);                 //! high trk muliplicity
DECLARE_SOA_COLUMN(HighFddFt0cMftFv0Mult, hasHighFddFt0cMftFv0Mult, bool); //! high FDD MFT FT0 FV0 multiplicity
DECLARE_SOA_COLUMN(HighFddMftFv0Mult, hasHighFddMftFv0Mult, bool);         //! high FDD MFT FV0 multiplicity
DECLARE_SOA_COLUMN(HighFv0MftMult, hasHighFv0MftMult, bool);               //! high MFT FV0 multiplicity
DECLARE_SOA_COLUMN(HighMultFv0, hasHighMultFv0, bool);                     //! high FV0 muliplicity
DECLARE_SOA_COLUMN(HighMftMult, hasHighMftMult, bool);                     //! high MFT multiplicity
DECLARE_SOA_COLUMN(HighFv0Flat, hasHighFv0Flat, bool);                     //! isotropic event FV0
DECLARE_SOA_COLUMN(HighMftTrkFlat, hasHighMftTrkFlat, bool);               //! isotropic event MFT
DECLARE_SOA_COLUMN(HighMftFv0TrkFlat, hasHighMftFv0TrkFlat, bool);         //! isotropic event MFT FV0 trk
DECLARE_SOA_COLUMN(HighMftFv0Flat, hasHighMftFv0Flat, bool);               //! isotropic event MFT FV0
DECLARE_SOA_COLUMN(HighMftFt0aMult, hasHighMftFt0aMult, bool);             //! high MFT T0A multiplicit
DECLARE_SOA_COLUMN(HighFt0Mult, hasHighFt0Mult, bool);                     //! high FT0 multiplicity
DECLARE_SOA_COLUMN(HighFt0Flat, hasHighFt0Flat, bool);                     //! isotropic event FT0
DECLARE_SOA_COLUMN(HighMftFt0aFlat, hasHighMftFt0aFlat, bool);             //! isotropic event MFT FT0A
DECLARE_SOA_COLUMN(HighFt0cFv0Mult, hasHighFt0cFv0Mult, bool);             //! high FT0C FV0 multiplicity
DECLARE_SOA_COLUMN(HighFt0cFv0Flat, hasHighFt0cFv0Flat, bool);             //! isotropic event FT0C FV0
DECLARE_SOA_COLUMN(LeadingPtTrack, hasLeadingPtTrack, bool);               //! event contains leading track

} // namespace filtering

namespace decision
{

DECLARE_SOA_COLUMN(BCId, hasBCId, int);                     //! Bunch crossing Id
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
                  filtering::HfHighPt2P, filtering::HfHighPt3P, filtering::HfBeauty3P, filtering::HfBeauty4P, filtering::HfFemto2P, filtering::HfFemto3P, filtering::HfDoubleCharm2P, filtering::HfDoubleCharm3P, filtering::HfDoubleCharmMix);

using HfFilter = HfFilters::iterator;

// correlations
DECLARE_SOA_TABLE(CFFiltersTwoN, "AOD", "CFFiltersTwoN", //!
                  filtering::PD, filtering::LD);
using CFFilterTwoN = CFFiltersTwoN::iterator;

DECLARE_SOA_TABLE(CFFilters, "AOD", "CFFilters", //!
                  filtering::PPP, filtering::PPL, filtering::PLL, filtering::LLL);
using CfFilter = CFFilters::iterator;

// jets
DECLARE_SOA_TABLE(JetFilters, "AOD", "JetFilters", //!
                  filtering::CollisionId,
                  filtering::JetChHighPt);

using JetFilter = JetFilters::iterator;

// strangeness (lf)
DECLARE_SOA_TABLE(StrangenessFilters, "AOD", "LFStrgFilters", //!
                  filtering::Omega, filtering::hadronXi, filtering::DoubleXi, filtering::TripleXi, filtering::QuadrupleXi, filtering::SingleXiYN);

using StrangenessFilter = StrangenessFilters::iterator;

// multiplicity
DECLARE_SOA_TABLE(MultFilters, "AOD", "MultFilters", //!
                  filtering::HighTrackMult, filtering::HighFddFt0cMftFv0Mult, filtering::HighFddMftFv0Mult, filtering::HighFv0MftMult, filtering::HighMultFv0, filtering::HighMftMult, filtering::HighFv0Flat, filtering::HighMftTrkFlat, filtering::HighMftFv0TrkFlat, filtering::HighMftFv0Flat, filtering::HighMftFt0aMult, filtering::HighFt0Mult, filtering::HighFt0Flat, filtering::HighMftFt0aFlat, filtering::HighFt0cFv0Mult, filtering::HighFt0cFv0Flat, filtering::LeadingPtTrack);

using MultFilter = MultFilters::iterator;

// cefp decision
DECLARE_SOA_TABLE(CefpDecisions, "AOD", "CefpDecision", //!
                  decision::BCId, decision::CefpSelected);
using CefpDecision = CefpDecisions::iterator;

/// List of the available filters, the description of their tables and the name of the tasks
constexpr int NumberOfFilters{9};
constexpr std::array<char[32], NumberOfFilters> AvailableFilters{"NucleiFilters", "DiffractionFilters", "DqFilters", "HfFilters", "CFFiltersTwoN", "CFFilters", "JetFilters", "StrangenessFilters", "MultFilters"};
constexpr std::array<char[16], NumberOfFilters> FilterDescriptions{"NucleiFilters", "DiffFilters", "DqFilters", "HfFilters", "CFFiltersTwoN", "CFFilters", "JetFilters", "LFStrgFilters", "MultFilters"};
constexpr std::array<char[128], NumberOfFilters> FilteringTaskNames{"o2-analysis-nuclei-filter", "o2-analysis-diffraction-filter", "o2-analysis-dq-filter-pp", "o2-analysis-hf-filter", "o2-analysis-cf-twobodyfemto-filter", "o2-analysis-cf-threebodyfemto-filter", "o2-analysis-je-filter", "o2-analysis-lf-strangeness-filter", "o2-analysis-mult-filter"};
constexpr o2::framework::pack<NucleiFilters, DiffractionFilters, DqFilters, HfFilters, CFFilters, CFFiltersTwoN, JetFilters, StrangenessFilters, MultFilters> FiltersPack;
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

#endif // EVENTFILTERING_FILTERTABLES_H_
