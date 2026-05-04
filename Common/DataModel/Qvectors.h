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

#include <Framework/AnalysisDataModel.h>

#include <vector>

namespace o2::aod
{
namespace qvec
{
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(IsCalibrated, isCalibrated, bool);

DECLARE_SOA_COLUMN(QvecRe, qvecRe, std::vector<float>);
DECLARE_SOA_COLUMN(QvecIm, qvecIm, std::vector<float>);
DECLARE_SOA_COLUMN(QvecAmp, qvecAmp, std::vector<float>);

DECLARE_SOA_COLUMN(QvecFT0CReVec, qvecFT0CReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFT0CImVec, qvecFT0CImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFT0AReVec, qvecFT0AReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFT0AImVec, qvecFT0AImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFT0MReVec, qvecFT0MReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFT0MImVec, qvecFT0MImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFV0AReVec, qvecFV0AReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFV0AImVec, qvecFV0AImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecTPCposReVec, qvecTPCposReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecTPCposImVec, qvecTPCposImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecTPCnegReVec, qvecTPCnegReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecTPCnegImVec, qvecTPCnegImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecTPCallReVec, qvecTPCallReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecTPCallImVec, qvecTPCallImVec, std::vector<float>);

DECLARE_SOA_COLUMN(QvecFT0CRe, qvecFT0CRe, float);
DECLARE_SOA_COLUMN(QvecFT0CIm, qvecFT0CIm, float);
DECLARE_SOA_COLUMN(QvecFT0ARe, qvecFT0ARe, float);
DECLARE_SOA_COLUMN(QvecFT0AIm, qvecFT0AIm, float);
DECLARE_SOA_COLUMN(QvecFT0MRe, qvecFT0MRe, float);
DECLARE_SOA_COLUMN(QvecFT0MIm, qvecFT0MIm, float);
DECLARE_SOA_COLUMN(QvecFV0ARe, qvecFV0ARe, float);
DECLARE_SOA_COLUMN(QvecFV0AIm, qvecFV0AIm, float);
DECLARE_SOA_COLUMN(QvecTPCposRe, qvecTPCposRe, float);
DECLARE_SOA_COLUMN(QvecTPCposIm, qvecTPCposIm, float);
DECLARE_SOA_COLUMN(QvecTPCnegRe, qvecTPCnegRe, float);
DECLARE_SOA_COLUMN(QvecTPCnegIm, qvecTPCnegIm, float);
DECLARE_SOA_COLUMN(QvecTPCallRe, qvecTPCallRe, float);
DECLARE_SOA_COLUMN(QvecTPCallIm, qvecTPCallIm, float);

DECLARE_SOA_COLUMN(SumAmplFT0C, sumAmplFT0C, float);
DECLARE_SOA_COLUMN(SumAmplFT0A, sumAmplFT0A, float);
DECLARE_SOA_COLUMN(SumAmplFT0M, sumAmplFT0M, float);
DECLARE_SOA_COLUMN(SumAmplFV0A, sumAmplFV0A, float);
DECLARE_SOA_COLUMN(NTrkTPCpos, nTrkTPCpos, int);
DECLARE_SOA_COLUMN(NTrkTPCneg, nTrkTPCneg, int);
DECLARE_SOA_COLUMN(NTrkTPCall, nTrkTPCall, int);
DECLARE_SOA_COLUMN(LabelsTPCpos, labelsTPCpos, std::vector<int>);
DECLARE_SOA_COLUMN(LabelsTPCneg, labelsTPCneg, std::vector<int>);
DECLARE_SOA_COLUMN(LabelsTPCall, labelsTPCall, std::vector<int>);

// Deprecated, will be removed in future after transition time //
DECLARE_SOA_COLUMN(QvecBPosReVec, qvecBPosReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecBPosImVec, qvecBPosImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecBNegReVec, qvecBNegReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecBNegImVec, qvecBNegImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecBTotReVec, qvecBTotReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecBTotImVec, qvecBTotImVec, std::vector<float>);

DECLARE_SOA_COLUMN(QvecBPosRe, qvecBPosRe, float);
DECLARE_SOA_COLUMN(QvecBPosIm, qvecBPosIm, float);
DECLARE_SOA_COLUMN(QvecBNegRe, qvecBNegRe, float);
DECLARE_SOA_COLUMN(QvecBNegIm, qvecBNegIm, float);
DECLARE_SOA_COLUMN(QvecBTotRe, qvecBTotRe, float);
DECLARE_SOA_COLUMN(QvecBTotIm, qvecBTotIm, float);

DECLARE_SOA_COLUMN(NTrkBPos, nTrkBPos, int);
DECLARE_SOA_COLUMN(NTrkBNeg, nTrkBNeg, int);
DECLARE_SOA_COLUMN(NTrkBTot, nTrkBTot, int);
DECLARE_SOA_COLUMN(LabelsBPos, labelsBPos, std::vector<int>);
DECLARE_SOA_COLUMN(LabelsBNeg, labelsBNeg, std::vector<int>);
DECLARE_SOA_COLUMN(LabelsBTot, labelsBTot, std::vector<int>);
/////////////////////////////////////////////////////////////////
} // namespace qvec

DECLARE_SOA_TABLE(Qvectors, "AOD", "QVECTORDEVS", //! Table with all Qvectors.
                  qvec::Cent, qvec::IsCalibrated, qvec::QvecRe, qvec::QvecIm, qvec::QvecAmp);
using Qvector = Qvectors::iterator;

DECLARE_SOA_TABLE(QvectorFT0Cs, "AOD", "QVECTORSFT0C", qvec::IsCalibrated, qvec::QvecFT0CRe, qvec::QvecFT0CIm, qvec::SumAmplFT0C);
DECLARE_SOA_TABLE(QvectorFT0As, "AOD", "QVECTORSFT0A", qvec::IsCalibrated, qvec::QvecFT0ARe, qvec::QvecFT0AIm, qvec::SumAmplFT0A);
DECLARE_SOA_TABLE(QvectorFT0Ms, "AOD", "QVECTORSFT0M", qvec::IsCalibrated, qvec::QvecFT0MRe, qvec::QvecFT0MIm, qvec::SumAmplFT0M);
DECLARE_SOA_TABLE(QvectorFV0As, "AOD", "QVECTORSFV0A", qvec::IsCalibrated, qvec::QvecFV0ARe, qvec::QvecFV0AIm, qvec::SumAmplFV0A);
DECLARE_SOA_TABLE(QvectorTPCposs, "AOD", "QVECTORSTPCPOS", qvec::IsCalibrated, qvec::QvecTPCposRe, qvec::QvecTPCposIm, qvec::NTrkTPCpos, qvec::LabelsTPCpos);
DECLARE_SOA_TABLE(QvectorTPCnegs, "AOD", "QVECTORSTPCNEG", qvec::IsCalibrated, qvec::QvecTPCnegRe, qvec::QvecTPCnegIm, qvec::NTrkTPCneg, qvec::LabelsTPCneg);
DECLARE_SOA_TABLE(QvectorTPCalls, "AOD", "QVECTORSTPCALL", qvec::IsCalibrated, qvec::QvecTPCallRe, qvec::QvecTPCallIm, qvec::NTrkTPCall, qvec::LabelsTPCall);

DECLARE_SOA_TABLE(QvectorFT0CVecs, "AOD", "QVECTORSFT0CVEC", qvec::IsCalibrated, qvec::QvecFT0CReVec, qvec::QvecFT0CImVec, qvec::SumAmplFT0C);
DECLARE_SOA_TABLE(QvectorFT0AVecs, "AOD", "QVECTORSFT0AVEC", qvec::IsCalibrated, qvec::QvecFT0AReVec, qvec::QvecFT0AImVec, qvec::SumAmplFT0A);
DECLARE_SOA_TABLE(QvectorFT0MVecs, "AOD", "QVECTORSFT0MVEC", qvec::IsCalibrated, qvec::QvecFT0MReVec, qvec::QvecFT0MImVec, qvec::SumAmplFT0M);
DECLARE_SOA_TABLE(QvectorFV0AVecs, "AOD", "QVECTORSFV0AVEC", qvec::IsCalibrated, qvec::QvecFV0AReVec, qvec::QvecFV0AImVec, qvec::SumAmplFV0A);
DECLARE_SOA_TABLE(QvectorTPCposVecs, "AOD", "QVECTORSTPCPVEC", qvec::IsCalibrated, qvec::QvecTPCposReVec, qvec::QvecTPCposImVec, qvec::NTrkTPCpos, qvec::LabelsTPCpos);
DECLARE_SOA_TABLE(QvectorTPCnegVecs, "AOD", "QVECTORSTPCNVEC", qvec::IsCalibrated, qvec::QvecTPCnegReVec, qvec::QvecTPCnegImVec, qvec::NTrkTPCneg, qvec::LabelsTPCneg);
DECLARE_SOA_TABLE(QvectorTPCallVecs, "AOD", "QVECTORSTPCAVEC", qvec::IsCalibrated, qvec::QvecTPCallReVec, qvec::QvecTPCallImVec, qvec::NTrkTPCall, qvec::LabelsTPCall);

using QvectorFT0C = QvectorFT0Cs::iterator;
using QvectorFT0A = QvectorFT0As::iterator;
using QvectorFT0M = QvectorFT0Ms::iterator;
using QvectorFV0A = QvectorFV0As::iterator;
using QvectorTPCpos = QvectorTPCposs::iterator;
using QvectorTPCneg = QvectorTPCnegs::iterator;
using QvectorTPCall = QvectorTPCalls::iterator;

using QvectorFT0CVec = QvectorFT0CVecs::iterator;
using QvectorFT0AVec = QvectorFT0AVecs::iterator;
using QvectorFT0MVec = QvectorFT0MVecs::iterator;
using QvectorFV0AVec = QvectorFV0AVecs::iterator;
using QvectorTPCposVec = QvectorTPCposVecs::iterator;
using QvectorTPCnegVec = QvectorTPCnegVecs::iterator;
using QvectorTPCallVec = QvectorTPCallVecs::iterator;

// Deprecated, will be removed in future after transition time //
DECLARE_SOA_TABLE(QvectorBPoss, "AOD", "QVECTORSBPOS", qvec::IsCalibrated, qvec::QvecBPosRe, qvec::QvecBPosIm, qvec::NTrkBPos, qvec::LabelsBPos);
DECLARE_SOA_TABLE(QvectorBNegs, "AOD", "QVECTORSBNEG", qvec::IsCalibrated, qvec::QvecBNegRe, qvec::QvecBNegIm, qvec::NTrkBNeg, qvec::LabelsBNeg);
DECLARE_SOA_TABLE(QvectorBTots, "AOD", "QVECTORSBTOT", qvec::IsCalibrated, qvec::QvecBTotRe, qvec::QvecBTotIm, qvec::NTrkBTot, qvec::LabelsBTot);

DECLARE_SOA_TABLE(QvectorBPosVecs, "AOD", "QVECTORSBPOSVEC", qvec::IsCalibrated, qvec::QvecBPosReVec, qvec::QvecBPosImVec, qvec::NTrkBPos, qvec::LabelsBPos);
DECLARE_SOA_TABLE(QvectorBNegVecs, "AOD", "QVECTORSBNEGVEC", qvec::IsCalibrated, qvec::QvecBNegReVec, qvec::QvecBNegImVec, qvec::NTrkBNeg, qvec::LabelsBNeg);
DECLARE_SOA_TABLE(QvectorBTotVecs, "AOD", "QVECTORSBTOTVEC", qvec::IsCalibrated, qvec::QvecBTotReVec, qvec::QvecBTotImVec, qvec::NTrkBTot, qvec::LabelsBTot);

using QvectorBPos = QvectorBPoss::iterator;
using QvectorBNeg = QvectorBNegs::iterator;
using QvectorBTot = QvectorBTots::iterator;

using QvectorBPosVec = QvectorBPosVecs::iterator;
using QvectorBNegVec = QvectorBNegVecs::iterator;
using QvectorBTotVec = QvectorBTotVecs::iterator;
/////////////////////////////////////////////////////////////////

namespace ese_qvec
{
DECLARE_SOA_COLUMN(IsCalibrated, isCalibrated, bool);

DECLARE_SOA_COLUMN(EseQvecRe, eseQvecRe, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecIm, eseQvecIm, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecAmp, eseQvecAmp, std::vector<float>);

DECLARE_SOA_COLUMN(EseQvecFT0CReVec, eseQvecFT0CReVec, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecFT0CImVec, eseQvecFT0CImVec, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecFT0AReVec, eseQvecFT0AReVec, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecFT0AImVec, eseQvecFT0AImVec, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecFT0MReVec, eseQvecFT0MReVec, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecFT0MImVec, eseQvecFT0MImVec, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecFV0AReVec, eseQvecFV0AReVec, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecFV0AImVec, eseQvecFV0AImVec, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecTPCposReVec, eseQvecTPCposReVec, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecTPCposImVec, eseQvecTPCposImVec, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecTPCnegReVec, eseQvecTPCnegReVec, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecTPCnegImVec, eseQvecTPCnegImVec, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecTPCallReVec, eseQvecTPCallReVec, std::vector<float>);
DECLARE_SOA_COLUMN(EseQvecTPCallImVec, eseQvecTPCallImVec, std::vector<float>);

DECLARE_SOA_COLUMN(EseQvecFT0CRe, eseQvecFT0CRe, float);
DECLARE_SOA_COLUMN(EseQvecFT0CIm, eseQvecFT0CIm, float);
DECLARE_SOA_COLUMN(EseQvecFT0ARe, eseQvecFT0ARe, float);
DECLARE_SOA_COLUMN(EseQvecFT0AIm, eseQvecFT0AIm, float);
DECLARE_SOA_COLUMN(EseQvecFT0MRe, eseQvecFT0MRe, float);
DECLARE_SOA_COLUMN(EseQvecFT0MIm, eseQvecFT0MIm, float);
DECLARE_SOA_COLUMN(EseQvecFV0ARe, eseQvecFV0ARe, float);
DECLARE_SOA_COLUMN(EseQvecFV0AIm, eseQvecFV0AIm, float);
DECLARE_SOA_COLUMN(EseQvecTPCposRe, eseQvecTPCposRe, float);
DECLARE_SOA_COLUMN(EseQvecTPCposIm, eseQvecTPCposIm, float);
DECLARE_SOA_COLUMN(EseQvecTPCnegRe, eseQvecTPCnegRe, float);
DECLARE_SOA_COLUMN(EseQvecTPCnegIm, eseQvecTPCnegIm, float);
DECLARE_SOA_COLUMN(EseQvecTPCallRe, eseQvecTPCallRe, float);
DECLARE_SOA_COLUMN(EseQvecTPCallIm, eseQvecTPCallIm, float);

DECLARE_SOA_COLUMN(EseRedQFT0C, eseRedQFT0C, float);
DECLARE_SOA_COLUMN(EseRedQFT0A, eseRedQFT0A, float);
DECLARE_SOA_COLUMN(EseRedQFT0M, eseRedQFT0M, float);
DECLARE_SOA_COLUMN(EseRedQFV0A, eseRedQFV0A, float);
/////////////////////////////////////////////////////////////////
} // namespace ese_qvec

DECLARE_SOA_TABLE(EseQvectors, "AOD", "ESEQVECTORDEVS", //! Table with all Qvectors.
                  qvec::Cent, ese_qvec::IsCalibrated, ese_qvec::EseQvecRe, ese_qvec::EseQvecIm, ese_qvec::EseQvecAmp);
using Qvector = Qvectors::iterator;

DECLARE_SOA_TABLE(EseQvecFT0Cs, "AOD", "ESEQVECFT0C", ese_qvec::IsCalibrated, ese_qvec::EseQvecFT0CRe, ese_qvec::EseQvecFT0CIm, qvec::SumAmplFT0C);
DECLARE_SOA_TABLE(EseQvecFT0As, "AOD", "ESEQVECFT0A", ese_qvec::IsCalibrated, ese_qvec::EseQvecFT0ARe, ese_qvec::EseQvecFT0AIm, qvec::SumAmplFT0A);
DECLARE_SOA_TABLE(EseQvecFT0Ms, "AOD", "ESEQVECFT0M", ese_qvec::IsCalibrated, ese_qvec::EseQvecFT0MRe, ese_qvec::EseQvecFT0MIm, qvec::SumAmplFT0M);
DECLARE_SOA_TABLE(EseQvecFV0As, "AOD", "ESEQVECFV0A", ese_qvec::IsCalibrated, ese_qvec::EseQvecFV0ARe, ese_qvec::EseQvecFV0AIm, qvec::SumAmplFV0A);
DECLARE_SOA_TABLE(EseQvecTPCposs, "AOD", "ESEQVECTPCPOS", ese_qvec::IsCalibrated, ese_qvec::EseQvecTPCposRe, ese_qvec::EseQvecTPCposIm, qvec::NTrkTPCpos, qvec::LabelsTPCpos);
DECLARE_SOA_TABLE(EseQvecTPCnegs, "AOD", "ESEQVECTPCNEG", ese_qvec::IsCalibrated, ese_qvec::EseQvecTPCnegRe, ese_qvec::EseQvecTPCnegIm, qvec::NTrkTPCneg, qvec::LabelsTPCneg);
DECLARE_SOA_TABLE(EseQvecTPCalls, "AOD", "ESEQVECTPCALL", ese_qvec::IsCalibrated, ese_qvec::EseQvecTPCallRe, ese_qvec::EseQvecTPCallIm, qvec::NTrkTPCall, qvec::LabelsTPCall);

DECLARE_SOA_TABLE(EseQvecFT0CVecs, "AOD", "ESEQVECFT0CVEC", ese_qvec::IsCalibrated, ese_qvec::EseQvecFT0CReVec, ese_qvec::EseQvecFT0CImVec, qvec::SumAmplFT0C);
DECLARE_SOA_TABLE(EseQvecFT0AVecs, "AOD", "ESEQVECFT0AVEC", ese_qvec::IsCalibrated, ese_qvec::EseQvecFT0AReVec, ese_qvec::EseQvecFT0AImVec, qvec::SumAmplFT0A);
DECLARE_SOA_TABLE(EseQvecFT0MVecs, "AOD", "ESEQVECFT0MVEC", ese_qvec::IsCalibrated, ese_qvec::EseQvecFT0MReVec, ese_qvec::EseQvecFT0MImVec, qvec::SumAmplFT0M);
DECLARE_SOA_TABLE(EseQvecFV0AVecs, "AOD", "ESEQVECFV0AVEC", ese_qvec::IsCalibrated, ese_qvec::EseQvecFV0AReVec, ese_qvec::EseQvecFV0AImVec, qvec::SumAmplFV0A);
DECLARE_SOA_TABLE(EseQvecTPCposVecs, "AOD", "ESEQVECTPCPVEC", ese_qvec::IsCalibrated, ese_qvec::EseQvecTPCposReVec, ese_qvec::EseQvecTPCposImVec, qvec::NTrkTPCpos, qvec::LabelsTPCpos);
DECLARE_SOA_TABLE(EseQvecTPCnegVecs, "AOD", "ESEQVECTPCNVEC", ese_qvec::IsCalibrated, ese_qvec::EseQvecTPCnegReVec, ese_qvec::EseQvecTPCnegImVec, qvec::NTrkTPCneg, qvec::LabelsTPCneg);
DECLARE_SOA_TABLE(EseQvecTPCallVecs, "AOD", "ESEQVECTPCAVEC", ese_qvec::IsCalibrated, ese_qvec::EseQvecTPCallReVec, ese_qvec::EseQvecTPCallImVec, qvec::NTrkTPCall, qvec::LabelsTPCall);

DECLARE_SOA_TABLE(EseQvecPercs, "AOD", "ESEQVECPERC", ese_qvec::EseQvecFT0CRe, ese_qvec::EseQvecFT0CIm, qvec::SumAmplFT0C,
                  ese_qvec::EseQvecFT0ARe, ese_qvec::EseQvecFT0AIm, qvec::SumAmplFT0A,
                  ese_qvec::EseQvecFT0MRe, ese_qvec::EseQvecFT0MIm, qvec::SumAmplFT0M,
                  ese_qvec::EseQvecFV0ARe, ese_qvec::EseQvecFV0AIm, qvec::SumAmplFV0A,
                  ese_qvec::EseQvecTPCposRe, ese_qvec::EseQvecTPCposIm, qvec::NTrkTPCpos,
                  ese_qvec::EseQvecTPCnegRe, ese_qvec::EseQvecTPCnegIm, qvec::NTrkTPCneg,
                  ese_qvec::EseQvecTPCallRe, ese_qvec::EseQvecTPCallIm, qvec::NTrkTPCall);

using EseQvectorFT0C = EseQvecFT0Cs::iterator;
using EseQvectorFT0A = EseQvecFT0As::iterator;
using EseQvectorFT0M = EseQvecFT0Ms::iterator;
using EseQvectorFV0A = EseQvecFV0As::iterator;
using EseQvectorTPCpos = EseQvecTPCposs::iterator;
using EseQvectorTPCneg = EseQvecTPCnegs::iterator;
using EseQvectorTPCall = EseQvecTPCalls::iterator;

using EseQvectorFT0CVec = EseQvecFT0CVecs::iterator;
using EseQvectorFT0AVec = EseQvecFT0AVecs::iterator;
using EseQvectorFT0MVec = EseQvecFT0MVecs::iterator;
using EseQvectorFV0AVec = EseQvecFV0AVecs::iterator;
using EseQvectorTPCposVec = EseQvecTPCposVecs::iterator;
using EseQvectorTPCnegVec = EseQvecTPCnegVecs::iterator;
using EseQvectorTPCallVec = EseQvecTPCallVecs::iterator;
/////////////////////////////////////////////////////////////////
} // namespace o2::aod

#endif // COMMON_DATAMODEL_QVECTORS_H_
