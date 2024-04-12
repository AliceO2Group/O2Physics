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

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"

#ifndef PWGLF_DATAMODEL_LFSTRANGENESSMLTABLES_H_
#define PWGLF_DATAMODEL_LFSTRANGENESSMLTABLES_H_

using namespace o2;
using namespace o2::framework;


// Creating output TTree for ML analysis
namespace o2::aod 
{
namespace v0mlcandidates
{
DECLARE_SOA_COLUMN(PosITSCls, posITSCls, int);
DECLARE_SOA_COLUMN(NegITSCls, negITSCls, int);
DECLARE_SOA_COLUMN(PosITSClSize, posITSClSize, int);
DECLARE_SOA_COLUMN(NegITSClSize, negITSClSize, int);
DECLARE_SOA_COLUMN(PosTPCRows, posTPCRows, float);
DECLARE_SOA_COLUMN(NegTPCRows, negTPCRows, float);
DECLARE_SOA_COLUMN(PosTPCSigmaPi, posTPCSigmaPi, float);
DECLARE_SOA_COLUMN(NegTPCSigmaPi, negTPCSigmaPi, float);
DECLARE_SOA_COLUMN(PosTPCSigmaPr, posTPCSigmaPr, float);
DECLARE_SOA_COLUMN(NegTPCSigmaPr, negTPCSigmaPr, float);
DECLARE_SOA_COLUMN(PosTPCSigmaEl, posTPCSigmaEl, float);
DECLARE_SOA_COLUMN(NegTPCSigmaEl, negTPCSigmaEl, float);
DECLARE_SOA_COLUMN(TOFSigmaLaPr, tofSigmaLaPr, float);
DECLARE_SOA_COLUMN(TOFSigmaLaPi, tofSigmaLaPi, float);
DECLARE_SOA_COLUMN(TOFSigmaALaPi, tofSigmaALaPi, float);
DECLARE_SOA_COLUMN(TOFSigmaALaPr, tofSigmaALaPr, float);
DECLARE_SOA_COLUMN(TOFSigmaK0PiPlus, tofSigmaK0PiPlus, float);
DECLARE_SOA_COLUMN(TOFSigmaK0PiMinus, tofSigmaK0PiMinus, float);
DECLARE_SOA_COLUMN(LambdaMass, lambdaMass, float);
DECLARE_SOA_COLUMN(AntiLambdaMass, antiLambdaMass, float);
DECLARE_SOA_COLUMN(GammaMass, gammaMass, float);
DECLARE_SOA_COLUMN(KZeroShortMass, kZeroShortMass, float);
DECLARE_SOA_COLUMN(PT, pT, float);
DECLARE_SOA_COLUMN(Qt, qt, float);
DECLARE_SOA_COLUMN(Alpha, alpha, float);
DECLARE_SOA_COLUMN(PosEta, posEta, float);
DECLARE_SOA_COLUMN(NegEta, negEta, float);
DECLARE_SOA_COLUMN(V0Eta, v0Eta, float);
DECLARE_SOA_COLUMN(Z, z, float);
DECLARE_SOA_COLUMN(V0radius, v0radius, float);
DECLARE_SOA_COLUMN(PA, pa, float);
DECLARE_SOA_COLUMN(DCApostopv, dcapostopv, float);
DECLARE_SOA_COLUMN(DCAnegtopv, dcanegtopv, float);
DECLARE_SOA_COLUMN(DCAV0daughters, dcaV0daughters, float);
DECLARE_SOA_COLUMN(DCAv0topv, dcav0topv, float);
DECLARE_SOA_COLUMN(PsiPair, psiPair, float);
DECLARE_SOA_COLUMN(V0type, v0type, uint8_t);
DECLARE_SOA_COLUMN(IsLambda, isLambda, bool);
DECLARE_SOA_COLUMN(IsAntiLambda, isAntiLambda, bool);
DECLARE_SOA_COLUMN(IsGamma, isGamma, bool);
DECLARE_SOA_COLUMN(IsKZeroShort, isKZeroShort, bool);
} // namespace v0mlcandidates

DECLARE_SOA_TABLE(V0MLCandidates, "AOD", "V0MLCANDIDATES",
				v0mlcandidates::PosITSCls,
				v0mlcandidates::NegITSCls,
				v0mlcandidates::PosITSClSize,
				v0mlcandidates::NegITSClSize,
				v0mlcandidates::PosTPCRows,
				v0mlcandidates::NegTPCRows, 
				v0mlcandidates::PosTPCSigmaPi,
				v0mlcandidates::NegTPCSigmaPi,
				v0mlcandidates::PosTPCSigmaPr,
				v0mlcandidates::NegTPCSigmaPr,
				v0mlcandidates::PosTPCSigmaEl,
				v0mlcandidates::NegTPCSigmaEl,
				v0mlcandidates::TOFSigmaLaPr, 
				v0mlcandidates::TOFSigmaLaPi, 
				v0mlcandidates::TOFSigmaALaPi,
				v0mlcandidates::TOFSigmaALaPr, 
				v0mlcandidates::TOFSigmaK0PiPlus,
				v0mlcandidates::TOFSigmaK0PiMinus, 
				v0mlcandidates::LambdaMass, 
				v0mlcandidates::AntiLambdaMass, 
				v0mlcandidates::GammaMass, 
				v0mlcandidates::KZeroShortMass, 
				v0mlcandidates::PT, 
				v0mlcandidates::Qt, 
				v0mlcandidates::Alpha,
				v0mlcandidates::PosEta, 
				v0mlcandidates::NegEta, 
				v0mlcandidates::V0Eta, 
				v0mlcandidates::Z, 
				v0mlcandidates::V0radius,
				v0mlcandidates::PA, 
				v0mlcandidates::DCApostopv, 
				v0mlcandidates::DCAnegtopv, 
				v0mlcandidates::DCAV0daughters, 
				v0mlcandidates::DCAv0topv, 
				v0mlcandidates::PsiPair, 
				v0mlcandidates::V0type, 
				v0mlcandidates::IsLambda, 
				v0mlcandidates::IsAntiLambda, 
				v0mlcandidates::IsGamma, 
				v0mlcandidates::IsKZeroShort);


namespace V0MLSelection
{
DECLARE_SOA_COLUMN(GammaBDTScore, gammaBDTScore, float);
DECLARE_SOA_COLUMN(LambdaBDTScore, lambdaBDTScore, float);
DECLARE_SOA_COLUMN(AntiLambdaBDTScore, antiLambdaBDTScore, float); 
DECLARE_SOA_COLUMN(K0ShortBDTScore, k0ShortBDTScore, float);
} // namespace V0MLSelection

DECLARE_SOA_TABLE(V0GammaMLScores, "AOD", "V0GaMLScores",
				  V0MLSelection::GammaBDTScore);
DECLARE_SOA_TABLE(V0LambdaMLScores, "AOD", "V0LaMLScores",
				  V0MLSelection::LambdaBDTScore);
DECLARE_SOA_TABLE(V0AntiLambdaMLScores, "AOD", "V0ALaMLScores",
				  V0MLSelection::AntiLambdaBDTScore);
DECLARE_SOA_TABLE(V0K0ShortMLScores, "AOD", "V0K0MLScores",
				  V0MLSelection::K0ShortBDTScore);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSTRANGENESSMLTABLES_H_