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
DECLARE_SOA_COLUMN(QvecFT0ARe, qvecFT0ARe, float); //! Real part of Qvec in FT0A.
DECLARE_SOA_COLUMN(QvecFT0AIm, qvecFT0AIm, float); //! Imaginary part for FT0A.
DECLARE_SOA_COLUMN(QvecFT0CRe, qvecFT0CRe, float); //! Real part of Qvec in FT0C.
DECLARE_SOA_COLUMN(QvecFT0CIm, qvecFT0CIm, float); //! Imaginary part for FT0C.
DECLARE_SOA_COLUMN(QvecFV0ARe, qvecFV0ARe, float); //! Real part of Qvec in FV0A.
DECLARE_SOA_COLUMN(QvecFV0AIm, qvecFV0AIm, float); //! Imaginary part for FV0A.
DECLARE_SOA_COLUMN(QvecBPosRe, qvecBPosRe, float);
DECLARE_SOA_COLUMN(QvecBPosIm, qvecBPosIm, float);
DECLARE_SOA_COLUMN(QvecBNegRe, qvecBNegRe, float);
DECLARE_SOA_COLUMN(QvecBNegIm, qvecBNegIm, float);
DECLARE_SOA_COLUMN(QvecFT0CUncorRe, qvecFT0CUncorRe, float);
DECLARE_SOA_COLUMN(QvecFT0CUncorIm, qvecFT0CUncorIm, float);
DECLARE_SOA_COLUMN(QvecFT0CRectrRe, qvecFT0CRectrRe, float);
DECLARE_SOA_COLUMN(QvecFT0CRectrIm, qvecFT0CRectrIm, float);
DECLARE_SOA_COLUMN(QvecFT0CTwistRe, qvecFT0CTwistRe, float);
DECLARE_SOA_COLUMN(QvecFT0CTwistIm, qvecFT0CTwistIm, float);
/// NOTE: Add here Qx,Qy for other systems.
} // namespace qvec

DECLARE_SOA_TABLE(Qvectors, "AOD", "QVECTORS", //! Table with all Qvectors.
                  qvec::Cent,
                  qvec::QvecFT0ARe, qvec::QvecFT0AIm,
                  qvec::QvecFT0CRe, qvec::QvecFT0CIm,
                  qvec::QvecFV0ARe, qvec::QvecFV0AIm,
                  qvec::QvecBPosRe, qvec::QvecBPosIm,
                  qvec::QvecBNegRe, qvec::QvecBNegIm,
                  qvec::QvecFT0CUncorRe, qvec::QvecFT0CUncorIm,
                  qvec::QvecFT0CRectrRe, qvec::QvecFT0CRectrIm,
                  qvec::QvecFT0CTwistRe, qvec::QvecFT0CTwistIm);

using Qvector = Qvectors::iterator;
} // namespace o2::aod

#endif // COMMON_DATAMODEL_QVECTORS_H_
