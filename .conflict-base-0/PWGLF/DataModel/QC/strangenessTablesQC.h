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
///
/// \brief QA task for Cascade analysis using derived data
///
/// \author Roman Nepeivoda (roman.nepeivoda@cern.ch)

#ifndef PWGLF_DATAMODEL_QC_STRANGENESSTABLESQC_H_
#define PWGLF_DATAMODEL_QC_STRANGENESSTABLESQC_H_

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace cascadesQC
{
DECLARE_SOA_COLUMN(CascCosPA, casccosPA, float);                 //!
DECLARE_SOA_COLUMN(V0CosPA, v0cosPA, float);                     //! needs to be changed to double
DECLARE_SOA_COLUMN(CascRadius, cascradius, float);               //!
DECLARE_SOA_COLUMN(V0Radius, v0radius, float);                   //! V0 decay radius (2D, centered at zero)
DECLARE_SOA_COLUMN(Sign, sign, int);                             //!
DECLARE_SOA_COLUMN(YXi, yXi, float);                             //!
DECLARE_SOA_COLUMN(YOmega, yOmega, float);                       //!
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);             //! cascade decay length
DECLARE_SOA_COLUMN(LifetimeXi, lifetimeXi, float);               //!
DECLARE_SOA_COLUMN(LifetimeOmega, lifetimeOmega, float);         //!
DECLARE_SOA_COLUMN(LifetimeV0, lifetimeV0, float);               //!
DECLARE_SOA_COLUMN(DCAV0Daughters, dcaV0daughters, float);       //! DCA between V0 daughters
DECLARE_SOA_COLUMN(DCACascDaughters, dcacascdaughters, float);   //!
DECLARE_SOA_COLUMN(DCAV0toPV, dcav0topv, float);                 //!
DECLARE_SOA_COLUMN(DCAbachtoPV, dcabachtopv, float);             //!
DECLARE_SOA_COLUMN(DCAPosToPV, dcapostopv, float);               //! DCA positive prong to PV
DECLARE_SOA_COLUMN(DCANegToPV, dcanegtopv, float);               //! DCA negative prong to PV
DECLARE_SOA_COLUMN(PosNSigmaV0Pion, posNSigmaV0Pion, float);     //!
DECLARE_SOA_COLUMN(PosNSigmaV0Proton, posNSigmaV0Proton, float); //!
DECLARE_SOA_COLUMN(NegNSigmaV0Pion, negNSigmaV0Pion, float);     //!
DECLARE_SOA_COLUMN(NegNSigmaV0Proton, negNSigmaV0Proton, float); //!
DECLARE_SOA_COLUMN(BachNSigmaV0Pion, bachNSigmaV0Pion, float);   //!
DECLARE_SOA_COLUMN(BachNSigmaV0Kaon, bachNSigmaV0Kaon, float);   //!
DECLARE_SOA_COLUMN(Pt, pt, float);                               //!
DECLARE_SOA_COLUMN(Eta, eta, float);                             //!
DECLARE_SOA_COLUMN(MLambda, mLambda, float);                     //!
DECLARE_SOA_COLUMN(MOmega, mOmega, float);                       //!
DECLARE_SOA_COLUMN(MXi, mXi, float);                             //!

} // namespace cascadesQC

DECLARE_SOA_TABLE(CascadesQC, "AOD", "CASCADESQC", o2::soa::Index<>,
                  cascadesQC::Sign, cascadesQC::YXi, cascadesQC::YOmega,
                  cascadesQC::CascCosPA, cascadesQC::V0CosPA,
                  cascadesQC::CascRadius, cascadesQC::V0Radius,
                  cascadesQC::DecayLength, cascadesQC::LifetimeXi, cascadesQC::LifetimeOmega, cascadesQC::LifetimeV0,
                  cascadesQC::DCAV0Daughters, cascadesQC::DCACascDaughters, cascadesQC::DCAV0toPV,
                  cascadesQC::DCAbachtoPV, cascadesQC::DCAPosToPV, cascadesQC::DCANegToPV,
                  cascadesQC::PosNSigmaV0Pion, cascadesQC::PosNSigmaV0Proton,
                  cascadesQC::NegNSigmaV0Pion, cascadesQC::NegNSigmaV0Proton,
                  cascadesQC::BachNSigmaV0Pion, cascadesQC::BachNSigmaV0Kaon,
                  cascadesQC::Pt, cascadesQC::Eta,
                  cascadesQC::MLambda, cascadesQC::MOmega, cascadesQC::MXi);

namespace vZerosQC
{
DECLARE_SOA_COLUMN(V0CosPA, v0cosPA, float);                     //! needs to be changed to double
DECLARE_SOA_COLUMN(YK0Short, yK0Short, float);                   //! V0 y with K0short hypothesis
DECLARE_SOA_COLUMN(YLambda, yLambda, float);                     //! V0 y with lambda or antilambda hypothesis
DECLARE_SOA_COLUMN(DCAV0Daughters, dcaV0daughters, float);       //! DCA between V0 daughters
DECLARE_SOA_COLUMN(DCAV0toPV, dcav0topv, float);                 //! DCA V0 to PV
DECLARE_SOA_COLUMN(DCAPosToPV, dcapostopv, float);               //! DCA positive prong to PV
DECLARE_SOA_COLUMN(DCANegToPV, dcanegtopv, float);               //! DCA negative prong to PV
DECLARE_SOA_COLUMN(V0Radius, v0radius, float);                   //! V0 decay radius (2D, centered at zero)
DECLARE_SOA_COLUMN(PosNSigmaV0Pion, posNSigmaV0Pion, float);     //! number of TPC sigmas for a pos daughter to be a pion
DECLARE_SOA_COLUMN(PosNSigmaV0Proton, posNSigmaV0Proton, float); //! number of TPC sigmas for a pos daughter to be a proton
DECLARE_SOA_COLUMN(NegNSigmaV0Pion, negNSigmaV0Pion, float);     //! number of TPC sigmas for a neg daughter to be a pion
DECLARE_SOA_COLUMN(NegNSigmaV0Proton, negNSigmaV0Proton, float); //! number of TPC sigmas for a neg daughter to be a proton
DECLARE_SOA_COLUMN(LifetimeLambda, lifetimeLambda, float);       //! lifetime with lambda mass assumption
DECLARE_SOA_COLUMN(LifetimeK0s, lifetimeK0s, float);             //! lifetime with K0s mass assumption
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);             //! V0 decay length
DECLARE_SOA_COLUMN(PosITSNhits, posITSNhits, int);               //! number of ITS hits of pos daughter
DECLARE_SOA_COLUMN(NegITSNhits, negITSNhits, int);               //! number of ITS hits of neg daughter
DECLARE_SOA_COLUMN(Pt, pt, float);                               //! transverse momentum in GeV/c
DECLARE_SOA_COLUMN(Eta, eta, float);                             //! pseudorapidity of V0
DECLARE_SOA_COLUMN(PosEta, poseta, float);                       //! pseudorapidity of V0 pos daughter
DECLARE_SOA_COLUMN(NegEta, negeta, float);                       //! pseudorapidity of V0 neg daughter
DECLARE_SOA_COLUMN(PosPhi, posphi, float);                       //!
DECLARE_SOA_COLUMN(NegPhi, negphi, float);                       //!
DECLARE_SOA_COLUMN(MK0Short, mK0Short, float);                   //! inv mass with K0s assumption
DECLARE_SOA_COLUMN(MLambda, mLambda, float);                     //! inv mass with lambda assumption
DECLARE_SOA_COLUMN(MAntiLambda, mAntiLambda, float);             //! inv mass with anti-lambda assumption
} // namespace vZerosQC

DECLARE_SOA_TABLE(VZerosQC, "AOD", "VZEROSQC", o2::soa::Index<>,
                  vZerosQC::YK0Short, vZerosQC::YLambda,
                  vZerosQC::DCAV0Daughters, vZerosQC::DCAV0toPV,
                  vZerosQC::DCAPosToPV, vZerosQC::DCANegToPV,
                  vZerosQC::V0Radius, vZerosQC::V0CosPA,
                  vZerosQC::PosNSigmaV0Pion, vZerosQC::PosNSigmaV0Proton,
                  vZerosQC::NegNSigmaV0Pion, vZerosQC::NegNSigmaV0Proton,
                  vZerosQC::DecayLength, vZerosQC::LifetimeLambda, vZerosQC::LifetimeK0s,
                  vZerosQC::PosITSNhits, vZerosQC::NegITSNhits,
                  vZerosQC::Pt, vZerosQC::Eta,
                  vZerosQC::MK0Short, vZerosQC::MLambda, vZerosQC::MAntiLambda,
                  vZerosQC::PosEta, vZerosQC::NegEta,
                  vZerosQC::PosPhi, vZerosQC::NegPhi);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_QC_STRANGENESSTABLESQC_H_
