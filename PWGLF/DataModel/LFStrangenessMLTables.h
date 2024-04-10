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
DECLARE_SOA_COLUMN(LambdaMass, lambdaMass, float);
DECLARE_SOA_COLUMN(AntiLambdaMass, antilambdaMass, float);
DECLARE_SOA_COLUMN(GammaMass, gammaMass, float);
DECLARE_SOA_COLUMN(KZeroShortMass, kZeroShortMass, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Qt, qt, float);
DECLARE_SOA_COLUMN(Alpha, alpha, float);
DECLARE_SOA_COLUMN(Radius, radius, float);
DECLARE_SOA_COLUMN(CosPA, cosPA, float);
DECLARE_SOA_COLUMN(DCADau, dcaDau, float);
DECLARE_SOA_COLUMN(DCANegPV, dcaNegPV, float);
DECLARE_SOA_COLUMN(DCAPosPV, dcaPosPV, float);
DECLARE_SOA_COLUMN(IsLambda, isLambda, bool);
DECLARE_SOA_COLUMN(IsAntiLambda, isAntiLambda, bool);
DECLARE_SOA_COLUMN(IsGamma, isGamma, bool);
DECLARE_SOA_COLUMN(IsKZeroShort, isKZeroShort, bool);
} // namespace v0mlcandidates

DECLARE_SOA_TABLE(V0MLCandidates, "AOD", "V0MLCANDIDATES",
				v0mlcandidates::LambdaMass,
				v0mlcandidates::AntiLambdaMass,
				v0mlcandidates::GammaMass,
				v0mlcandidates::KZeroShortMass,
				v0mlcandidates::Pt,  
				v0mlcandidates::Qt, 
				v0mlcandidates::Alpha,
				v0mlcandidates::Radius,
				v0mlcandidates::CosPA,
				v0mlcandidates::DCADau, 
				v0mlcandidates::DCANegPV,
				v0mlcandidates::DCAPosPV,
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