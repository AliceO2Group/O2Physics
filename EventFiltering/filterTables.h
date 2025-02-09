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
#include <cstdint>

namespace o2::aod
{
template <uint32_t T>
struct Hash;
}

#include "Framework/ASoA.h"

namespace o2::soa
{
template <typename T>
  requires(!std::same_as<typename o2::aod::MetadataTrait<std::decay_t<T>>::metadata, void>)
const char* getTableLabel()
{
  return o2::aod::MetadataTrait<std::decay_t<T>>::metadata::tableLabel();
}

template <typename T>
  requires requires { T::ref.label_hash; }
const char* getTableLabel()
{
  return o2::aod::Hash<T::ref.label_hash>::str;
}
} // namespace o2::soa

namespace o2::aod
{
namespace filtering
{
DECLARE_SOA_COLUMN(H2, hasH2, bool);             //! deuteron trigger for the helium normalisation (to be downscaled)
DECLARE_SOA_COLUMN(He, hasHe, bool);             //! helium
DECLARE_SOA_COLUMN(H3L3Body, hasH3L3Body, bool); //! hypertriton 3body
DECLARE_SOA_COLUMN(ITSextremeIonisation, hasITSextremeIonisation, bool); //! ITS extreme ionisation
DECLARE_SOA_COLUMN(ITSmildIonisation, hasITSmildIonisation, bool);       //! ITS mild ionisation (normalisation of the extreme ionisation), to be downscaled

// diffraction
DECLARE_SOA_COLUMN(UDdiffSmall, hasDiffSmall, bool); //! Double Gap events, <= 3 prongs
DECLARE_SOA_COLUMN(UDdiffLarge, hasDiffLarge, bool); //! Double Gap events, >= 4 prongs

DECLARE_SOA_COLUMN(UDdiffBC, hasDiffBC, bool); //! diffractive BC

// Dileptons & Quarkonia
DECLARE_SOA_COLUMN(SingleE, hasSingleE, bool);           //! single electron trigger
DECLARE_SOA_COLUMN(SingleMuLow, hasSingleMuLow, bool);   //! single muon with low pT trigger
DECLARE_SOA_COLUMN(SingleMuHigh, hasSingleMuHigh, bool); //! single muon with high pT trigger
DECLARE_SOA_COLUMN(DiElectron, hasDiElectron, bool);     //! dielectron trigger
DECLARE_SOA_COLUMN(DiMuon, hasDiMuon, bool);             //! dimuon trigger with low pT on muons
// EM dielectrons
DECLARE_SOA_COLUMN(LMeeIMR, hasLMeeIMR, bool); //! dielectron trigger for intermediate mass region
DECLARE_SOA_COLUMN(LMeeHMR, hasLMeeHMR, bool); //! dielectron trigger for high mass region
// Electron-muon pair
DECLARE_SOA_COLUMN(ElectronMuon, hasElectronMuon, bool); //! dimuon trigger with low pT on muons

// heavy flavours
DECLARE_SOA_COLUMN(HfHighPt2P, hasHfHighPt2P, bool);                             //! high-pT 2-prong charm hadron
DECLARE_SOA_COLUMN(HfHighPt3P, hasHfHighPt3P, bool);                             //! high-pT 3-prong charm hadron
DECLARE_SOA_COLUMN(HfBeauty3P, hasHfBeauty3P, bool);                             //! 3-prong beauty hadron
DECLARE_SOA_COLUMN(HfBeauty4P, hasHfBeauty4P, bool);                             //! 4-prong beauty hadron
DECLARE_SOA_COLUMN(HfFemto2P, hasHfFemto2P, bool);                               //! 2-prong charm-hadron - N pair
DECLARE_SOA_COLUMN(HfFemto3P, hasHfFemto3P, bool);                               //! 3-prong charm-hadron - N pair
DECLARE_SOA_COLUMN(HfDoubleCharm2P, hasHfDoubleCharm2P, bool);                   //! at least two 2-prong charm-hadron candidates
DECLARE_SOA_COLUMN(HfDoubleCharm3P, hasHfDoubleCharm3P, bool);                   //! at least two 3-prong charm-hadron candidates
DECLARE_SOA_COLUMN(HfDoubleCharmMix, hasHfDoubleCharmMix, bool);                 //! at least one 2-prong and one 3-prong charm-hadron candidates
DECLARE_SOA_COLUMN(HfV0Charm2P, hasHfV0Charm2P, bool);                           //! V0 with 2-prong charm hadron
DECLARE_SOA_COLUMN(HfV0Charm3P, hasHfV0Charm3P, bool);                           //! V0 with 3-prong charm hadron
DECLARE_SOA_COLUMN(HfCharmBarToXiBach, hasHfCharmBarToXiBach, bool);             //! Charm baryon to Xi + bachelor
DECLARE_SOA_COLUMN(HfCharmBarToXiBachBach, hasHfCharmBarToXiBachBach, bool);     //! Charm baryon to Xi + 2 bachelors
DECLARE_SOA_COLUMN(HfSigmaCPPK, hasHfSigmaCPPK, bool);                           //! SigmaC(2455)++K- and SigmaC(2520)++K- + c.c.
DECLARE_SOA_COLUMN(HfSigmaC0K0, hasHfSigmaC0K0, bool);                           //! SigmaC(2455)0KS0 and SigmaC(2520)0KS0
DECLARE_SOA_COLUMN(HfPhotonCharm2P, hasHfPhotonCharm2P, bool);                   //! photon with 2-prong charm hadron
DECLARE_SOA_COLUMN(HfPhotonCharm3P, hasHfPhotonCharm3P, bool);                   //! photon with 3-prong charm hadron
DECLARE_SOA_COLUMN(HfSingleCharm2P, hasHfSingleCharm2P, bool);                   //! 2-prong charm hadron (for efficiency studies)
DECLARE_SOA_COLUMN(HfSingleCharm3P, hasHfSingleCharm3P, bool);                   //! 3-prong charm hadron (for efficiency studies)
DECLARE_SOA_COLUMN(HfSingleNonPromptCharm2P, hasHfSingleNonPromptCharm2P, bool); //! 2-prong charm hadron (for efficiency studies)
DECLARE_SOA_COLUMN(HfSingleNonPromptCharm3P, hasHfSingleNonPromptCharm3P, bool); //! 3-prong charm hadron (for efficiency studies)

// CF two body triggers
DECLARE_SOA_COLUMN(PD, hasPD, bool); //! has d-p pair
DECLARE_SOA_COLUMN(LD, hasLD, bool); //! has l-d pair

// CF three body triggers
DECLARE_SOA_COLUMN(PPP, hasPPP, bool);     //! has p-p-p triplet
DECLARE_SOA_COLUMN(PPL, hasPPL, bool);     //! has p-p-L triplet
DECLARE_SOA_COLUMN(PLL, hasPLL, bool);     //! has p-L-L triplet
DECLARE_SOA_COLUMN(LLL, hasLLL, bool);     //! has L-L-L tripletD
DECLARE_SOA_COLUMN(PPPHI, hasPPPHI, bool); //! has P-P-PHI triplet

// jets
DECLARE_SOA_COLUMN(JetChLowPt, hasJetChLowPt, bool);   //! low-pT charged jet
DECLARE_SOA_COLUMN(JetChHighPt, hasJetChHighPt, bool); //! high-pT charged jet
DECLARE_SOA_COLUMN(TrackLowPt, hasTrackLowPt, bool);   //! low-pT track
DECLARE_SOA_COLUMN(TrackHighPt, hasTrackHighPt, bool); //! high-pT track

// hf-jets
DECLARE_SOA_COLUMN(JetD0ChLowPt, hasJetD0ChLowPt, bool);   //! low-pT charged D0 jet
DECLARE_SOA_COLUMN(JetD0ChHighPt, hasJetD0ChHighPt, bool); //! high-pT charged D0 jet
DECLARE_SOA_COLUMN(JetLcChLowPt, hasJetLcChLowPt, bool);   //! low-pT charged Lc jet
DECLARE_SOA_COLUMN(JetLcChHighPt, hasJetLcChHighPt, bool); //! high-pT charged Lc jet

// full jets
DECLARE_SOA_COLUMN(EMCALReadout, hasEMCALinReadout, bool);               //! EMCAL readout
DECLARE_SOA_COLUMN(JetFullHighPt, hasJetFullHighPt, bool);               //! high-pT full jet
DECLARE_SOA_COLUMN(JetFullLowPt, hasJetFullLowPt, bool);                 //! low-pT full jet
DECLARE_SOA_COLUMN(JetNeutralHighPt, hasJetNeutralHighPt, bool);         //! high-pT neutral jet
DECLARE_SOA_COLUMN(JetNeutralLowPt, hasJetNeutralLowPt, bool);           //! low-pT neutral jet
DECLARE_SOA_COLUMN(GammaVeryHighPtEMCAL, hasGammaVeryHighPtEMCAL, bool); //! Photon trigger in EMCAL, very high threshold
DECLARE_SOA_COLUMN(GammaVeryHighPtDCAL, hasGammaVeryHighPtDCAL, bool);   //! Photon trigger in DCAL, very high threshold
DECLARE_SOA_COLUMN(GammaHighPtEMCAL, hasGammaHighPtEMCAL, bool);         //! Photon trigger in EMCAL, high threshold
DECLARE_SOA_COLUMN(GammaHighPtDCAL, hasGammaHighPtDCAL, bool);           //! Photon trigger in DCAL, high threshold
DECLARE_SOA_COLUMN(GammaLowPtEMCAL, hasGammaLowPtEMCAL, bool);           //! Photon trigger in EMCAL, low threshold
DECLARE_SOA_COLUMN(GammaLowPtDCAL, hasGammaLowPtDCAL, bool);             //! Photon trigger in DCAL, low threshold
DECLARE_SOA_COLUMN(GammaVeryLowPtEMCAL, hasGammaVeryLowPtEMCAL, bool);   //! Photon trigger in EMCAL, very low threshold
DECLARE_SOA_COLUMN(GammaVeryLowPtDCAL, hasGammaVeryLowPtDCAL, bool);     //! Photon trigger in DCAL, very low threshold

// strangeness (lf)
DECLARE_SOA_COLUMN(Omega, hasOmega, bool);                       //! at leat 1 Omega
DECLARE_SOA_COLUMN(hadronOmega, hashadronOmega, bool);           //! at least 1 Omega + high-pt hadron
DECLARE_SOA_COLUMN(DoubleXi, hasDoubleXi, bool);                 //! at least 2 Xi
DECLARE_SOA_COLUMN(TripleXi, hasTripleXi, bool);                 //! at least 3 Xi
DECLARE_SOA_COLUMN(QuadrupleXi, hasQuadrupleXi, bool);           //! at least 4 Xi
DECLARE_SOA_COLUMN(DoubleOmega, hasDoubleOmega, bool);           //! at least 2 Omega
DECLARE_SOA_COLUMN(OmegaXi, hasOmegaXi, bool);                   //! at least 1 Omega + 1 Xi
DECLARE_SOA_COLUMN(SingleXiYN, hasSingleXiYN, bool);             //! at least 1 Xi with high radius (YN interactions)
DECLARE_SOA_COLUMN(OmegaLargeRadius, hasOmegaLargeRadius, bool); //! at least 1 Omega with high radius
DECLARE_SOA_COLUMN(TrackedCascade, hasTrackedCascade, bool);     //! at least 1 tracked cascade
DECLARE_SOA_COLUMN(TrackedXi, hasTrackedXi, bool);               //! at least 1 tracked Xi
DECLARE_SOA_COLUMN(TrackedOmega, hasTrackedOmega, bool);         //! at least 1 tracked Omega
DECLARE_SOA_COLUMN(Tracked3Body, hasTracked3Body, bool);         //! at least 1 tracked 3Body
DECLARE_SOA_COLUMN(OmegaHighMult, hasOmegaHighMult, bool);       //! at least 1 Omega + high-mult event

// F1-proton
DECLARE_SOA_COLUMN(TriggerEventF1Proton, triggereventf1proton, bool); //! F1 - proton femto trigger event

// multiplicity
DECLARE_SOA_COLUMN(HighTrackMult, hasHighTrackMult, bool);     //! high trk muliplicity
DECLARE_SOA_COLUMN(HighMultFv0, hasHighMultFv0, bool);         //! high FV0 muliplicity
DECLARE_SOA_COLUMN(HighFt0Mult, hasHighFt0Mult, bool);         //! high FT0 multiplicity
DECLARE_SOA_COLUMN(HighFt0Flat, hasHighFt0Flat, bool);         //! isotropic event FT0
DECLARE_SOA_COLUMN(HighFt0cFv0Mult, hasHighFt0cFv0Mult, bool); //! high FT0C FV0 multiplicity
DECLARE_SOA_COLUMN(HighFt0cFv0Flat, hasHighFt0cFv0Flat, bool); //! isotropic event FT0C FV0
DECLARE_SOA_COLUMN(LeadingPtTrack, hasLeadingPtTrack, bool);   //! event contains leading track

// photons
DECLARE_SOA_COLUMN(PHOSPhoton, hasPHOSPhoton, bool); //! PHOS single photons
DECLARE_SOA_COLUMN(PHOSnbar, hasPHOSnbar, bool);     //! PHOS antineutrons
// DECLARE_SOA_COLUMN(PHOSElectron, hasPHOSElectron, bool);       //! PHOS single electron
// DECLARE_SOA_COLUMN(PHOSPair, hasPHOSpair, bool);               //! PHOS photon pair
DECLARE_SOA_COLUMN(PCMHighPtPhoton, hasPCMHighPtPhoton, bool); //! PCM high pT photon
// DECLARE_SOA_COLUMN(PCMMatCalib, hasPCMMatCalib, bool);         //! PCM material budget calibration
// DECLARE_SOA_COLUMN(PCMEtaDalitz, hasPCMEtaDalitz, bool);       //! PCM eta -> ee gamma
// DECLARE_SOA_COLUMN(PCMEtaGG, hasPCMEtaGG, bool);               //! PCM eta -> ee gamma
DECLARE_SOA_COLUMN(PCMandEE, hasPCMandEE, bool); //! PCM and ee
} // namespace filtering

namespace decision
{

DECLARE_SOA_COLUMN(BCId, bcIndex, uint64_t);                   //! Bunch crossing Id
DECLARE_SOA_COLUMN(GlobalBCId, globalBC, uint64_t);            //! Global Bunch crossing Id
DECLARE_SOA_COLUMN(EvSelBC, evSelBC, uint64_t);                //! Global Bunch crossing Id
DECLARE_SOA_COLUMN(CollisionTime, collisionTime, float);       //! Collision time
DECLARE_SOA_COLUMN(CollisionTimeRes, collisionTimeRes, float); //! Collision time resolution
DECLARE_SOA_COLUMN(CefpTriggered0, cefpTriggered0, uint64_t);  //! CEFP triggers before downscalings
DECLARE_SOA_COLUMN(CefpTriggered1, cefpTriggered1, uint64_t);  //! CEFP triggers before downscalings
DECLARE_SOA_COLUMN(CefpSelected0, cefpSelected0, uint64_t);    //! CEFP decision
DECLARE_SOA_COLUMN(CefpSelected1, cefpSelected1, uint64_t);    //! CEFP decision

} // namespace decision

namespace bcrange
{
DECLARE_SOA_COLUMN(BCstart, hasBCstart, uint64_t); //! CEFP triggers before downscalings
DECLARE_SOA_COLUMN(BCend, hasBCend, uint64_t);     //! CEFP bcrange

} // namespace bcrange

// nuclei
DECLARE_SOA_TABLE(NucleiFilters, "AOD", "NucleiFilters", //!
                  filtering::H2, filtering::He, filtering::H3L3Body, filtering::ITSmildIonisation,
                  filtering::ITSextremeIonisation);
using NucleiFilter = NucleiFilters::iterator;

// diffraction
DECLARE_SOA_TABLE(DiffractionFilters, "AOD", "DiffFilters", //! Diffraction filters (Collisions)
                  filtering::UDdiffSmall, filtering::UDdiffLarge);
using DiffractionFilter = DiffractionFilters::iterator;

DECLARE_SOA_TABLE(DiffractionBCFilters, "AOD", "DiffBCFilters", //! Diffraction filters (BCs)
                  filtering::UDdiffBC);
using DiffractionBCFilter = DiffractionBCFilters::iterator;

// Dileptons & Quarkonia
DECLARE_SOA_TABLE(DqFilters, "AOD", "DqFilters", //!
                  filtering::SingleE, filtering::LMeeIMR, filtering::LMeeHMR, filtering::DiElectron, filtering::SingleMuLow, filtering::SingleMuHigh, filtering::DiMuon, filtering::ElectronMuon);
using DqFilter = DqFilters::iterator;

// heavy flavours
DECLARE_SOA_TABLE(HfFilters, "AOD", "HfFilters", //!
                  filtering::HfHighPt2P,
                  filtering::HfHighPt3P,
                  filtering::HfBeauty3P,
                  filtering::HfBeauty4P,
                  filtering::HfFemto2P,
                  filtering::HfFemto3P,
                  filtering::HfDoubleCharm2P,
                  filtering::HfDoubleCharm3P,
                  filtering::HfDoubleCharmMix,
                  filtering::HfV0Charm2P,
                  filtering::HfV0Charm3P,
                  filtering::HfCharmBarToXiBach,
                  filtering::HfSigmaCPPK,
                  filtering::HfSigmaC0K0,
                  filtering::HfPhotonCharm2P,
                  filtering::HfPhotonCharm3P,
                  filtering::HfSingleCharm2P,
                  filtering::HfSingleCharm3P,
                  filtering::HfSingleNonPromptCharm2P,
                  filtering::HfSingleNonPromptCharm3P,
                  filtering::HfCharmBarToXiBachBach);

using HfFilter = HfFilters::iterator;

DECLARE_SOA_TABLE(CFFilters, "AOD", "CFFilters", //!
                  filtering::PPP, filtering::PPL, filtering::PLL, filtering::LLL, filtering::PPPHI, filtering::PD, filtering::LD);
using CfFilter = CFFilters::iterator;

// jets
DECLARE_SOA_TABLE(JetFilters, "AOD", "JetFilters", //!
                  filtering::JetChLowPt,
                  filtering::JetChHighPt,
                  filtering::TrackLowPt,
                  filtering::TrackHighPt);

using JetFilter = JetFilters::iterator;

DECLARE_SOA_TABLE(JetHFFilters, "AOD", "JetHFFilters", //!
                  filtering::JetD0ChLowPt,
                  filtering::JetD0ChHighPt,
                  filtering::JetLcChLowPt,
                  filtering::JetLcChHighPt);

using JetHFFilter = JetHFFilters::iterator;

DECLARE_SOA_TABLE(FullJetFilters, "AOD", "FullJetFilters", //!
                  filtering::EMCALReadout, filtering::JetFullHighPt, filtering::JetFullLowPt, filtering::JetNeutralHighPt, filtering::JetNeutralLowPt, filtering::GammaVeryHighPtEMCAL, filtering::GammaVeryHighPtDCAL, filtering::GammaHighPtEMCAL, filtering::GammaHighPtDCAL, filtering::GammaLowPtEMCAL, filtering::GammaLowPtDCAL, filtering::GammaVeryLowPtEMCAL, filtering::GammaVeryLowPtDCAL);

using FullJetFilter = FullJetFilters::iterator;

// strangeness (lf)
DECLARE_SOA_TABLE(StrangenessFilters, "AOD", "LFStrgFilters", //!
                  filtering::Omega, filtering::hadronOmega, filtering::DoubleXi, filtering::TripleXi, filtering::QuadrupleXi, filtering::SingleXiYN, filtering::OmegaLargeRadius, filtering::TrackedXi, filtering::TrackedOmega, filtering::OmegaHighMult, filtering::DoubleOmega, filtering::OmegaXi);

using StrangenessFilter = StrangenessFilters::iterator;

// F1 proton
DECLARE_SOA_TABLE(F1ProtonFilters, "AOD", "F1ProtonFilters", //!
                  filtering::TriggerEventF1Proton);
using F1ProtonFilter = F1ProtonFilters::iterator;

// multiplicity
DECLARE_SOA_TABLE(MultFilters, "AOD", "MultFilters", //!
                  filtering::HighTrackMult, filtering::HighMultFv0, filtering::HighFt0Mult, filtering::HighFt0Flat, filtering::HighFt0cFv0Mult, filtering::HighFt0cFv0Flat, filtering::LeadingPtTrack);

using MultFilter = MultFilters::iterator;

// photons
DECLARE_SOA_TABLE(PhotonFilters, "AOD", "PhotonFilters", //!
                  filtering::PHOSPhoton, filtering::PHOSnbar, filtering::PCMHighPtPhoton, filtering::PCMandEE);

using PhotonFilter = PhotonFilters::iterator;

// cefp decision
DECLARE_SOA_TABLE(CefpDecisions, "AOD", "CefpDecision", //!
                  decision::BCId, decision::GlobalBCId, decision::EvSelBC, decision::CollisionTime, decision::CollisionTimeRes, decision::CefpTriggered0, decision::CefpTriggered1, decision::CefpSelected0, decision::CefpSelected1);
using CefpDecision = CefpDecisions::iterator;

// cefp decision
DECLARE_SOA_TABLE(BCRanges, "AOD", "BCRanges", //!
                  bcrange::BCstart, bcrange::BCend);
using BCRange = BCRanges::iterator;

/// List of the available filters, the description of their tables and the name of the tasks
constexpr int NumberOfFilters{12};
constexpr std::array<char[32], NumberOfFilters> AvailableFilters{"NucleiFilters", "DiffractionFilters", "DqFilters", "HfFilters", "CFFilters", "JetFilters", "JetHFFilters", "FullJetFilters", "StrangenessFilters", "MultFilters", "PhotonFilters", "F1ProtonFilters"};
constexpr std::array<char[16], NumberOfFilters> FilterDescriptions{"NucleiFilters", "DiffFilters", "DqFilters", "HfFilters", "CFFilters", "JetFilters", "JetHFFilters", "FullJetFilters", "LFStrgFilters", "MultFilters", "PhotonFilters", "F1ProtonFilters"};
constexpr std::array<char[128], NumberOfFilters> FilteringTaskNames{"o2-analysis-nuclei-filter", "o2-analysis-diffraction-filter", "o2-analysis-dq-filter-pp-with-association", "o2-analysis-hf-filter", "o2-analysis-cf-filter", "o2-analysis-je-filter", "o2-analysis-je-hf-filter", "o2-analysis-fje-filter", "o2-analysis-lf-strangeness-filter", "o2-analysis-mult-filter", "o2-analysis-em-photon-filter", "o2-analysis-lf-f1proton-filter"};
constexpr o2::framework::pack<NucleiFilters, DiffractionFilters, DqFilters, HfFilters, CFFilters, JetFilters, JetHFFilters, FullJetFilters, StrangenessFilters, MultFilters, PhotonFilters, F1ProtonFilters> FiltersPack;
static_assert(o2::framework::pack_size(FiltersPack) == NumberOfFilters);

template <typename T, typename C>
void addColumnToMap(std::unordered_map<std::string, std::unordered_map<std::string, float>>& map)
{
  map[o2::soa::getTableLabel<T>()][C::columnLabel()] = 1.f;
}

template <typename T, typename... C>
void addColumnsToMap(o2::framework::pack<C...>, std::unordered_map<std::string, std::unordered_map<std::string, float>>& map)
{
  ([&]() {
    if constexpr (soa::is_persistent_v<C>) {
      addColumnToMap<T, C>(map);
    }
  }(),
   ...);
}

template <typename... T>
void FillFiltersMap(o2::framework::pack<T...>, std::unordered_map<std::string, std::unordered_map<std::string, float>>& map)
{
  (addColumnsToMap<T>(typename T::table_t::persistent_columns_t{}, map), ...);
}

template <typename... C>
static std::vector<std::string> ColumnsNames(o2::framework::pack<C...>)
{
  std::vector<std::string> result;
  ([&]() {
    if constexpr (soa::is_persistent_v<C>) {
      result.push_back(C::columnLabel());
    }
  }(),
   ...);
  return result;
}

template <typename... C>
unsigned int NumberOfColumns(o2::framework::pack<C...>)
{
  unsigned int result = 0;
  ([&]() {
    if constexpr (soa::is_persistent_v<C>) {
      ++result;
    }
  }(),
   ...);
  return result;
}

} // namespace o2::aod

#endif // EVENTFILTERING_FILTERTABLES_H_
