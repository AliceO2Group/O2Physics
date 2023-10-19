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
/// \file   Qvectors.h
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  Declaration of the table for the (un)corrected Q-vectors for the event plane
/// determination.
///

#ifndef COMMON_DATAMODEL_QVECTORS_H_
#define COMMON_DATAMODEL_QVECTORS_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace qvec
{
DECLARE_SOA_COLUMN(Cent, cent, float);             //! Centrality percentile.
DECLARE_SOA_COLUMN(CentBin, centBin, int);
DECLARE_SOA_COLUMN(QvecUncorRe, qvecUncorRe, float);
DECLARE_SOA_COLUMN(QvecUncorIm, qvecUncorIm, float);
DECLARE_SOA_COLUMN(QvecRectrRe, qvecRectrRe, float);
DECLARE_SOA_COLUMN(QvecRectrIm, qvecRectrIm, float);
DECLARE_SOA_COLUMN(QvecTwistRe, qvecTwistRe, float);
DECLARE_SOA_COLUMN(QvecTwistIm, qvecTwistIm, float);
DECLARE_SOA_COLUMN(QvecFinalRe, qvecFinalRe, float);
DECLARE_SOA_COLUMN(QvecFinalIm, qvecFinalIm, float);
DECLARE_SOA_COLUMN(QvecBPosUncorRe, qvecBPosUncorRe, float);
DECLARE_SOA_COLUMN(QvecBPosUncorIm, qvecBPosUncorIm, float);
DECLARE_SOA_COLUMN(QvecBPosRectrRe, qvecBPosRectrRe, float);
DECLARE_SOA_COLUMN(QvecBPosRectrIm, qvecBPosRectrIm, float);
DECLARE_SOA_COLUMN(QvecBPosTwistRe, qvecBPosTwistRe, float);
DECLARE_SOA_COLUMN(QvecBPosTwistIm, qvecBPosTwistIm, float);
DECLARE_SOA_COLUMN(QvecBPosFinalRe, qvecBPosFinalRe, float);
DECLARE_SOA_COLUMN(QvecBPosFinalIm, qvecBPosFinalIm, float);
DECLARE_SOA_COLUMN(QvecBNegUncorRe, qvecBNegUncorRe, float);
DECLARE_SOA_COLUMN(QvecBNegUncorIm, qvecBNegUncorIm, float);
DECLARE_SOA_COLUMN(QvecBNegRectrRe, qvecBNegRectrRe, float);
DECLARE_SOA_COLUMN(QvecBNegRectrIm, qvecBNegRectrIm, float);
DECLARE_SOA_COLUMN(QvecBNegTwistRe, qvecBNegTwistRe, float);
DECLARE_SOA_COLUMN(QvecBNegTwistIm, qvecBNegTwistIm, float);
DECLARE_SOA_COLUMN(QvecBNegFinalRe, qvecBNegFinalRe, float);
DECLARE_SOA_COLUMN(QvecBNegFinalIm, qvecBNegFinalIm, float);
DECLARE_SOA_COLUMN(NTrkBPos, nTrkBPos, int);
DECLARE_SOA_COLUMN(NTrkBNeg, nTrkBNeg, int);
/// NOTE: Add here Qx,Qy for other systems.
} // namespace qvec

DECLARE_SOA_TABLE(Qvectors, "AOD", "QVECTORS", //! Table with all Qvectors.
                  qvec::Cent, qvec::CentBin,
                  qvec::QvecUncorRe, qvec::QvecUncorIm, qvec::QvecRectrRe, qvec::QvecRectrIm,
                  qvec::QvecTwistRe, qvec::QvecTwistIm, qvec::QvecFinalRe, qvec::QvecFinalIm,
                  qvec::QvecBPosUncorRe, qvec::QvecBPosUncorIm, qvec::QvecBPosRectrRe, qvec::QvecBPosRectrIm,
                  qvec::QvecBPosTwistRe, qvec::QvecBPosTwistIm, qvec::QvecBPosFinalRe, qvec::QvecBPosFinalIm,
                  qvec::QvecBNegUncorRe, qvec::QvecBNegUncorIm, qvec::QvecBNegRectrRe, qvec::QvecBNegRectrIm,
                  qvec::QvecBNegTwistRe, qvec::QvecBNegTwistIm, qvec::QvecBNegFinalRe, qvec::QvecBNegFinalIm,
                  qvec::NTrkBPos, qvec::NTrkBNeg);
using Qvector = Qvectors::iterator;
} // namespace o2::aod

#endif // COMMON_DATAMODEL_QVECTORS_H_
