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

/// \file CorrelationTables.h
/// \brief Correlation table definitions.
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#ifndef PWGHF_HFC_DATAMODEL_CORRELATIONTABLES_H_
#define PWGHF_HFC_DATAMODEL_CORRELATIONTABLES_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
// definition of columns and tables for D-Dbar correlation pairs
namespace hf_correlation_d_dbar
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);
DECLARE_SOA_COLUMN(PtD, ptD, float);
DECLARE_SOA_COLUMN(PtDbar, ptDbar, float);
DECLARE_SOA_COLUMN(MD, mD, float);
DECLARE_SOA_COLUMN(MDbar, mDbar, float);
DECLARE_SOA_COLUMN(SignalStatus, signalStatus, int);
} // namespace hf_correlation_d_dbar

DECLARE_SOA_TABLE(DDbarPair, "AOD", "DDBARPAIR",
                  aod::hf_correlation_d_dbar::DeltaPhi,
                  aod::hf_correlation_d_dbar::DeltaEta,
                  aod::hf_correlation_d_dbar::PtD,
                  aod::hf_correlation_d_dbar::PtDbar);

DECLARE_SOA_TABLE(DDbarRecoInfo, "AOD", "DDBARRECOINFO",
                  aod::hf_correlation_d_dbar::MD,
                  aod::hf_correlation_d_dbar::MDbar,
                  aod::hf_correlation_d_dbar::SignalStatus);

// definition of columns and tables for D0-Hadron correlation pairs
namespace hf_correlation_d0_hadron
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);       //! DeltaPhi between D0 and Hadrons
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);       //! DeltaEta between D0 and Hadrons
DECLARE_SOA_COLUMN(PtD, ptD, float);                 //! Transverse momentum of D0
DECLARE_SOA_COLUMN(PtHadron, ptHadron, float);       //! Transverse momentum of Hadron
DECLARE_SOA_COLUMN(MD, mD, float);                   //! Invariant mass of D0
DECLARE_SOA_COLUMN(MDbar, mDbar, float);             //! Invariant mass of D0bar
DECLARE_SOA_COLUMN(SignalStatus, signalStatus, int); //! Tag for D0,D0bar
DECLARE_SOA_COLUMN(PoolBin, poolBin, int);           //! Pool Bin for the MixedEvent

enum ParticleTypeData {
  D0Only = 1,        // Identified as D0
  D0barOnly,         // Identified as D0bar
  D0D0barBoth,       // Identified as both D0 and D0bar
  D0OnlySoftPi = 11, // Identified as D0 with soft pion
  D0barOnlySoftPi,   // Identified as D0bar with soft pion
  D0D0barBothSoftPi  // Identified as both D0 and D0bar with soft pion
};

enum ParticleTypeMcRec {
  D0Sig = 0, // D0 signal
  D0Ref,     // D0 reflection
  D0Bg,      // D0 background
  D0barSig,  // D0bar signal
  D0barRef,  // D0bar reflection
  D0barBg,   // D0bar background
  SoftPi     // pairs including soft pion
};
} // namespace hf_correlation_d0_hadron

DECLARE_SOA_TABLE(DHadronPair, "AOD", "DHADRONPAIR", //! D0-Hadrons pairs Informations
                  aod::hf_correlation_d0_hadron::DeltaPhi,
                  aod::hf_correlation_d0_hadron::DeltaEta,
                  aod::hf_correlation_d0_hadron::PtD,
                  aod::hf_correlation_d0_hadron::PtHadron,
                  aod::hf_correlation_d0_hadron::PoolBin);

DECLARE_SOA_TABLE(DHadronRecoInfo, "AOD", "DHADRONRECOINFO", //! D0-Hadrons pairs Reconstructed Informations
                  aod::hf_correlation_d0_hadron::MD,
                  aod::hf_correlation_d0_hadron::MDbar,
                  aod::hf_correlation_d0_hadron::SignalStatus);

// definition of columns and tables for Ds-Hadron correlation pairs
namespace hf_correlation_ds_hadron
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float); //! DeltaPhi between Ds and Hadrons
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float); //! DeltaEta between Ds and Hadrons
DECLARE_SOA_COLUMN(PtD, ptD, float);           //! Transverse momentum of Ds
DECLARE_SOA_COLUMN(PtHadron, ptHadron, float); //! Transverse momentum of Hadron
DECLARE_SOA_COLUMN(MD, mD, float);             //! Invariant mass of Ds
DECLARE_SOA_COLUMN(PoolBin, poolBin, int);     //! Pool Bin for the MixedEvent
DECLARE_SOA_COLUMN(IsSignal, isSignal, bool);  //! Used in MC-Rec, Ds Signal
DECLARE_SOA_COLUMN(IsPrompt, isPrompt, bool);  //! Used in MC-Rec, Ds Prompt or Non-Prompt
} // namespace hf_correlation_ds_hadron

DECLARE_SOA_TABLE(DsHadronPair, "AOD", "DSHPAIR", //! Ds-Hadrons pairs Informations
                  aod::hf_correlation_ds_hadron::DeltaPhi,
                  aod::hf_correlation_ds_hadron::DeltaEta,
                  aod::hf_correlation_ds_hadron::PtD,
                  aod::hf_correlation_ds_hadron::PtHadron,
                  aod::hf_correlation_ds_hadron::PoolBin);

DECLARE_SOA_TABLE(DsHadronRecoInfo, "AOD", "DSHRECOINFO", //! Ds-Hadrons pairs Reconstructed Informations
                  aod::hf_correlation_ds_hadron::MD,
                  aod::hf_correlation_ds_hadron::IsSignal);

DECLARE_SOA_TABLE(DsHadronGenInfo, "AOD", "DSHGENINFO", //! Ds-Hadrons pairs Generated Informations
                  aod::hf_correlation_ds_hadron::IsPrompt);

// definition of columns and tables for Dplus-Hadron correlation pairs
namespace hf_correlation_dplus_hadron
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);        //! DeltaPhi between D+ and Hadrons
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);        //! DeltaEta between D+ and Hadrons
DECLARE_SOA_COLUMN(PtD, ptD, float);                  //! Transverse momentum of D+
DECLARE_SOA_COLUMN(PtHadron, ptHadron, float);        //! Transverse momentum of Hadron
DECLARE_SOA_COLUMN(MD, mD, float);                    //! Invariant mass of D+
DECLARE_SOA_COLUMN(SignalStatus, signalStatus, bool); //! Used in MC-Rec, D+ Signal
DECLARE_SOA_COLUMN(PoolBin, poolBin, int);            //! Pool Bin of event defined using zvtx and multiplicity
} // namespace hf_correlation_dplus_hadron

DECLARE_SOA_TABLE(DplusHadronPair, "AOD", "DPLUSHPAIR", //! D+-Hadrons pairs Informations
                  aod::hf_correlation_dplus_hadron::DeltaPhi,
                  aod::hf_correlation_dplus_hadron::DeltaEta,
                  aod::hf_correlation_dplus_hadron::PtD,
                  aod::hf_correlation_dplus_hadron::PtHadron,
                  aod::hf_correlation_dplus_hadron::PoolBin);

DECLARE_SOA_TABLE(DplusHadronRecoInfo, "AOD", "DPLUSHRECOINFO", //! D+-Hadrons pairs Reconstructed Informations
                  aod::hf_correlation_dplus_hadron::MD,
                  aod::hf_correlation_dplus_hadron::SignalStatus);

// Table for selection of Dmeson in a collision
namespace hf_selection_dmeson_collision
{
DECLARE_SOA_COLUMN(DmesonSel, dmesonSel, bool); //! Selection flag for D meson in a collision
} // namespace hf_selection_dmeson_collision

DECLARE_SOA_TABLE(DmesonSelection, "AOD", "DINCOLL", // Selection of D meson in collisions
                  aod::hf_selection_dmeson_collision::DmesonSel);
} // namespace o2::aod

#endif // PWGHF_HFC_DATAMODEL_CORRELATIONTABLES_H_
