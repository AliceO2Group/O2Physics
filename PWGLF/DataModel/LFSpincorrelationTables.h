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
/// \author Sourav Kundu <sourav.kundu@cern.ch>

#ifndef PWGLF_DATAMODEL_LFSPINCORRELATIONTABLES_H_
#define PWGLF_DATAMODEL_LFSPINCORRELATIONTABLES_H_

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

#include <cmath>

namespace o2::aod
{
namespace lambdaspincorrevent
{
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Posz, posz, float);
} // namespace lambdaspincorrevent
DECLARE_SOA_TABLE(LambdaSpinCorrEvents, "AOD", "LAMBDASPINCOREVENT",
                  o2::soa::Index<>,
                  lambdaspincorrevent::Cent,
                  lambdaspincorrevent::Posz)
using LambdaSpinCorrEvent = LambdaSpinCorrEvents::iterator;

namespace lambdaspincorrpair
{
DECLARE_SOA_INDEX_COLUMN(LambdaSpinCorrEvent, lambdaspincorrevent);
DECLARE_SOA_COLUMN(V0Status, v0Status, float);                     //! Lambda or Anti-Lambda status
DECLARE_SOA_COLUMN(DoubleStatus, doubleStatus, bool);              //! Double status
DECLARE_SOA_COLUMN(V0Cospa, v0Cospa, float);                       //! V0 Cospa
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);                     //! V0 Radius
DECLARE_SOA_COLUMN(DcaPositive, dcaPositive, float);               //! DCA Positive
DECLARE_SOA_COLUMN(DcaNegative, dcaNegative, float);               //! DCA Negative
DECLARE_SOA_COLUMN(DcaBetweenDaughter, dcaBetweenDaughter, float); //! DCA between daughters
DECLARE_SOA_COLUMN(LambdaPt, lambdaPt, float);                     //! Lambda Pt
DECLARE_SOA_COLUMN(LambdaEta, lambdaEta, float);                   //! Lambda Eta
DECLARE_SOA_COLUMN(LambdaPhi, lambdaPhi, float);                   //! Lambda Phi
DECLARE_SOA_COLUMN(LambdaMass, lambdaMass, float);                 //! Lambda Mass
DECLARE_SOA_COLUMN(ProtonPt, protonPt, float);                     //! Proton Pt
DECLARE_SOA_COLUMN(ProtonEta, protonEta, float);                   //! Proton Eta
DECLARE_SOA_COLUMN(ProtonPhi, protonPhi, float);                   //! Proton Phi
DECLARE_SOA_COLUMN(ProtonIndex, protonIndex, int);                 //! Proton index
DECLARE_SOA_COLUMN(PionIndex, pionIndex, int);                     //! Pion index
} // namespace lambdaspincorrpair
DECLARE_SOA_TABLE(LambdaSpinCorrPairs, "AOD", "LAMBDAPAIR",
                  o2::soa::Index<>,
                  lambdaspincorrpair::LambdaSpinCorrEventId,
                  lambdaspincorrpair::V0Status,
                  lambdaspincorrpair::DoubleStatus,
                  lambdaspincorrpair::V0Cospa,
                  lambdaspincorrpair::V0Radius,
                  lambdaspincorrpair::DcaPositive,
                  lambdaspincorrpair::DcaNegative,
                  lambdaspincorrpair::DcaBetweenDaughter,
                  lambdaspincorrpair::LambdaPt,
                  lambdaspincorrpair::LambdaEta,
                  lambdaspincorrpair::LambdaPhi,
                  lambdaspincorrpair::LambdaMass,
                  lambdaspincorrpair::ProtonPt,
                  lambdaspincorrpair::ProtonEta,
                  lambdaspincorrpair::ProtonPhi,
                  lambdaspincorrpair::ProtonIndex,
                  lambdaspincorrpair::PionIndex);

using LambdaSpinCorrPair = LambdaSpinCorrPairs::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFSPINCORRELATIONTABLES_H_
