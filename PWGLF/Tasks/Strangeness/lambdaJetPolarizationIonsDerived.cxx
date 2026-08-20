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
//
/// \file lambdaJetPolarizationIonsDerived.cxx
/// \brief Lambda and antiLambda polarization analysis task using derived data
/// \author Cicero Domenico Muncinelli <cicero.domenico.muncinelli@cern.ch>, Campinas State University
//
// Jet Polarization Ions task -- Derived data
// ================
//
// This code loops over custom derived data tables defined on
// lambdaJetPolarizationIons.h (JetsRing, LambdaLikeV0sRing).
// From this derived data, calculates polarization on an EbE
// basis (see TProfiles).
// Signal extraction is done out of the framework, based on
// the AnalysisResults of this code.
//
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    cicero.domenico.muncinelli@cern.ch
//

#include "PWGLF/DataModel/lambdaJetPolarizationIons.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/VectorUtil.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom3.h> // For perpendicular jet direction QAs

#include <algorithm> // std::fill, for resetting the Delta Method accumulators
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using ROOT::Math::PtEtaPhiMVector;
using ROOT::Math::XYZVector;

// Declaring constants:
constexpr double ProtonMass = o2::constants::physics::MassProton; // Assumes particle identification for daughter is perfect
constexpr double LambdaWeakDecayConstant = 0.749;                 // DPG 2025 update
constexpr double AntiLambdaWeakDecayConstant = -0.758;            // DPG 2025 update
constexpr double PolPrefactorLambda = 3.0 / LambdaWeakDecayConstant;
constexpr double PolPrefactorAntiLambda = 3.0 / AntiLambdaWeakDecayConstant;

enum CentEstimator {
  kCentFT0C = 0,
  kCentFT0M,
  kCentFV0A
};

// Helper macro to avoid writing the histogram fills 4 times for about 20 histograms:
#define RING_OBSERVABLE_FILL_LIST(X, FOLDER)                                                                               \
  /* Counters */                                                                                                           \
  X(FOLDER "/QA/hDeltaPhi", deltaPhiJet)                                                                                   \
  X(FOLDER "/QA/hDeltaPhiVsDeltaEta", deltaPhiJet, deltaEtaJet)                                                            \
  X(FOLDER "/QA/hDeltaTheta", deltaThetaJet)                                                                               \
  X(FOLDER "/QA/hCosDeltaTheta", cosDeltaThetaJet)                                                                         \
  X(FOLDER "/QA/hIntegrated", 0.)                                                                                          \
  X(FOLDER "/QA/hPtJet", leadingJetPt)                                                                                     \
  /* Lambda pT variation -- Youpeng's proposal */                                                                          \
  X(FOLDER "/QA/hLambdaPt", v0pt)                                                                                          \
  /* Counters */                                                                                                           \
  X(FOLDER "/QA/h2dDeltaPhiVsLambdaPt", deltaPhiJet, v0pt)                                                                 \
  X(FOLDER "/QA/h2dDeltaThetaVsLambdaPt", deltaThetaJet, v0pt)                                                             \
  X(FOLDER "/QA/hDeltaPhiVsLeadJetPhi", deltaPhiJet, leadingJetPhi)                                                        \
  /* Additional plots for instant gratification - 1D Profiles */                                                           \
  X(FOLDER "/hRingObservableCounts", ringObservable)                                                                       \
  X(FOLDER "/pRingObservableDeltaPhi", deltaPhiJet, ringObservable)                                                        \
  X(FOLDER "/pRingObservablePhiJet", leadingJetPhi, ringObservable)                                                        \
  X(FOLDER "/pRingObservablePhiLambda", v0phi, ringObservable)                                                             \
  X(FOLDER "/pRingObservableDeltaTheta", deltaThetaJet, ringObservable)                                                    \
  X(FOLDER "/EtaDependence/pRingObservableEtaLambda", v0eta, ringObservable)                                               \
  X(FOLDER "/EtaDependence/pRingObservableEtaJet", leadingJetEta, ringObservable)                                          \
  X(FOLDER "/pRingObservableIntegrated", 0., ringObservable)                                                               \
  X(FOLDER "/pRingObservableLambdaPt", v0pt, ringObservable)                                                               \
  X(FOLDER "/pRingObservableLeadJetPVz", collisionPVz, ringObservable)                                                     \
  X(FOLDER "/ProxyPtDependence/pRingVsPtJet", leadingJetPt, ringObservable)                                                \
  X(FOLDER "/ProxyPtDependence/pRingVsPtJetVsEtaJet", leadingJetPt, leadingJetEta, ringObservable)                         \
  X(FOLDER "/ProxyPtDependence/pRingVsPtJetVsEtaV0", leadingJetPt, v0eta, ringObservable)                                  \
  /* 2D Profiles */                                                                                                        \
  X(FOLDER "/p2dRingObservableDeltaPhiVsLambdaPt", deltaPhiJet, v0pt, ringObservable)                                      \
  X(FOLDER "/p2dRingObservableDeltaThetaVsLambdaPt", deltaThetaJet, v0pt, ringObservable)                                  \
  X(FOLDER "/p2dRingObservableDeltaPhiVsLeadJetPt", deltaPhiJet, leadingJetPt, ringObservable)                             \
  X(FOLDER "/p2dRingObservableDeltaThetaVsLeadJetPt", deltaThetaJet, leadingJetPt, ringObservable)                         \
  /* 1D Mass */                                                                                                            \
  X(FOLDER "/QA/hMass", v0LambdaLikeMass)                                                                                  \
  X(FOLDER "/QA/hRingObservableNumMass", v0LambdaLikeMass, ringObservable)                                                 \
  X(FOLDER "/hMassSigExtract", v0LambdaLikeMass)                                                                           \
  /* Counters */                                                                                                           \
  X(FOLDER "/QA/h2dDeltaPhiVsMass", deltaPhiJet, v0LambdaLikeMass)                                                         \
  X(FOLDER "/QA/h2dDeltaThetaVsMass", deltaThetaJet, v0LambdaLikeMass)                                                     \
  X(FOLDER "/QA/h3dDeltaPhiVsMassVsLambdaPt", deltaPhiJet, v0LambdaLikeMass, v0pt)                                         \
  X(FOLDER "/QA/h3dDeltaThetaVsMassVsLambdaPt", deltaThetaJet, v0LambdaLikeMass, v0pt)                                     \
  X(FOLDER "/QA/h3dDeltaPhiVsMassVsLeadJetPt", deltaPhiJet, v0LambdaLikeMass, leadingJetPt)                                \
  X(FOLDER "/QA/h3dDeltaThetaVsMassVsLeadJetPt", deltaThetaJet, v0LambdaLikeMass, leadingJetPt)                            \
  X(FOLDER "/QA/h3dDeltaPhiVsMassVsCent", deltaPhiJet, v0LambdaLikeMass, centrality)                                       \
  X(FOLDER "/QA/h3dDeltaThetaVsMassVsCent", deltaThetaJet, v0LambdaLikeMass, centrality)                                   \
  /* TProfile of Ring vs Mass */                                                                                           \
  X(FOLDER "/pRingObservableMass", v0LambdaLikeMass, ringObservable)                                                       \
  /* TProfile of Ring vs Mass -- Leading Particle and 2nd-to-leading jet - QA */                                           \
  X(FOLDER "/pRingObservableLeadPMass", v0LambdaLikeMass, ringObservableLeadP)                                             \
  X(FOLDER "/pRingObservable2ndJetMass", v0LambdaLikeMass, ringObservable2ndJet)                                           \
  /* 2D Profiles: Angle vs Mass */                                                                                         \
  X(FOLDER "/p2dRingObservableDeltaPhiVsMass", deltaPhiJet, v0LambdaLikeMass, ringObservable)                              \
  X(FOLDER "/p2dRingObservableDeltaThetaVsMass", deltaThetaJet, v0LambdaLikeMass, ringObservable)                          \
  /* 2D Profile: Ring vs Eta variables */                                                                                  \
  X(FOLDER "/EtaDependence/hCounterEtaLambdaMinusEtaJet", v0eta - leadingJetEta)                                           \
  X(FOLDER "/EtaDependence/pRingObservableEtaLambdaMinusEtaJet", v0eta - leadingJetEta, ringObservable)                    \
  X(FOLDER "/EtaDependence/p2dRingObservableEtaLambdaVsEtaJet", v0eta, leadingJetEta, ringObservable)                      \
  X(FOLDER "/EtaDependence/h2dCounterEtaLambdaVsEtaJet", v0eta, leadingJetEta)                                             \
  X(FOLDER "/EtaDependence/p2dRingObservableEtaLambdaVsEtaJet_FineBins", v0eta, leadingJetEta, ringObservable)             \
  X(FOLDER "/EtaDependence/h2dCounterEtaLambdaVsEtaJet_FineBins", v0eta, leadingJetEta)                                    \
  /* 3D Profiles: Angle vs Mass vs Lambda pT */                                                                            \
  X(FOLDER "/p3dRingObservableDeltaPhiVsMassVsLambdaPt", deltaPhiJet, v0LambdaLikeMass, v0pt, ringObservable)              \
  X(FOLDER "/p3dRingObservableDeltaThetaVsMassVsLambdaPt", deltaThetaJet, v0LambdaLikeMass, v0pt, ringObservable)          \
  /* 3D Profiles: Angle vs Mass vs Lead Jet pT */                                                                          \
  X(FOLDER "/p3dRingObservableDeltaPhiVsMassVsLeadJetPt", deltaPhiJet, v0LambdaLikeMass, leadingJetPt, ringObservable)     \
  X(FOLDER "/p3dRingObservableDeltaThetaVsMassVsLeadJetPt", deltaThetaJet, v0LambdaLikeMass, leadingJetPt, ringObservable) \
  /* 2D Profile: Mass vs Centrality */                                                                                     \
  X(FOLDER "/p2dRingObservableMassVsCent", v0LambdaLikeMass, centrality, ringObservable)                                   \
  /* 3D Profiles: Angle vs Mass vs Centrality */                                                                           \
  X(FOLDER "/p3dRingObservableDeltaPhiVsMassVsCent", deltaPhiJet, v0LambdaLikeMass, centrality, ringObservable)            \
  X(FOLDER "/p3dRingObservableDeltaThetaVsMassVsCent", deltaThetaJet, v0LambdaLikeMass, centrality, ringObservable)        \
  X(FOLDER "/pRingVsCentrality", centrality, ringObservable)                                                               \
  X(FOLDER "/p2dRingObservableVsPxPy", v0px, v0py, ringObservable)                                                         \
  X(FOLDER "/p2dRingObservableVsPzPx", v0pz, v0px, ringObservable)
// (TODO: add counters for regular TH2Ds about centrality)

// For leading particle
#define RING_OBSERVABLE_LEADP_FILL_LIST(X, FOLDER)                                                      \
  X(FOLDER "/QA/hDeltaPhiLeadP", deltaPhiLeadP)                                                         \
  X(FOLDER "/QA/hDeltaThetaLeadP", deltaThetaLeadP)                                                     \
  X(FOLDER "/QA/hPtLeadP", leadPPt)                                                                     \
  X(FOLDER "/QA/hCosDeltaThetaLeadP", cosDeltaThetaLeadP)                                               \
  X(FOLDER "/hRingObservableLeadPCounts", ringObservableLeadP)                                          \
  X(FOLDER "/pRingObservableLeadPDeltaPhi", deltaPhiLeadP, ringObservableLeadP)                         \
  X(FOLDER "/pRingObservableLeadPDeltaTheta", deltaThetaLeadP, ringObservableLeadP)                     \
  X(FOLDER "/EtaDependence/pRingObservableEtaLambdaLeadP", v0eta, ringObservableLeadP)                  \
  X(FOLDER "/EtaDependence/pRingObservableEtaLeadP", leadPEta, ringObservableLeadP)                     \
  X(FOLDER "/pRingObservableLeadPIntegrated", 0., ringObservableLeadP)                                  \
  X(FOLDER "/pRingObservableLeadPLambdaPt", v0pt, ringObservableLeadP)                                  \
  X(FOLDER "/pRingObservableLeadPPVz", collisionPVz, ringObservableLeadP)                               \
  X(FOLDER "/ProxyPtDependence/pRingVsPtLeadP", leadPPt, ringObservableLeadP)                           \
  X(FOLDER "/ProxyPtDependence/pRingVsPtLeadPVsEtaLeadP", leadPPt, leadPEta, ringObservableLeadP)       \
  X(FOLDER "/ProxyPtDependence/pRingVsPtLeadPVsEtaV0", leadPPt, v0eta, ringObservableLeadP)             \
  X(FOLDER "/EtaDependence/p2dRingObservableEtaLambdaVsEtaLeadP", v0eta, leadPEta, ringObservableLeadP) \
  X(FOLDER "/EtaDependence/h2dCounterEtaLambdaVsEtaLeadP", v0eta, leadPEta)                             \
  X(FOLDER "/p2dRingObservableLeadPVsPxPy", v0px, v0py, ringObservableLeadP)                            \
  X(FOLDER "/p2dRingObservableLeadPVsPzPx", v0pz, v0px, ringObservableLeadP)

// A macro that encapsulates all eta checks for leading particle and V0s, along with the fills
// Parameters:
//   FOLDER       -- histogram folder string (compile-time literal)
//   LEADP_IS_POS -- bool: leadPEtaPos
//   V0_IS_POS    -- bool: lambdaEtaPos
#define RING_OBSERVABLE_LEADP_ETA_SPLIT_FILL_LIST(FOLDER, LEADP_IS_POS, V0_IS_POS)                                      \
  do {                                                                                                                  \
    if (LEADP_IS_POS) {                                                                                                 \
      /* leadP marginal: positive side */                                                                               \
      APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_PosEtaLeadP", leadPPt, ringObservableLeadP)            \
      if (V0_IS_POS) {                                                                                                  \
        /* V0 marginal: positive side */                                                                                \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_PosEtaV0", leadPPt, ringObservableLeadP)             \
        /* Joint: (+leadP, +V0) */                                                                                      \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_PosEtaLeadP_PosEtaV0", leadPPt, ringObservableLeadP) \
      } else {                                                                                                          \
        /* V0 marginal: negative side */                                                                                \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_NegEtaV0", leadPPt, ringObservableLeadP)             \
        /* Joint: (+leadP, -V0) */                                                                                      \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_PosEtaLeadP_NegEtaV0", leadPPt, ringObservableLeadP) \
      }                                                                                                                 \
    } else {                                                                                                            \
      /* leadP marginal: negative side */                                                                               \
      APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_NegEtaLeadP", leadPPt, ringObservableLeadP)            \
      if (V0_IS_POS) {                                                                                                  \
        /* V0 marginal: positive side */                                                                                \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_PosEtaV0", leadPPt, ringObservableLeadP)             \
        /* Joint: (-leadP, +V0) */                                                                                      \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_NegEtaLeadP_PosEtaV0", leadPPt, ringObservableLeadP) \
      } else {                                                                                                          \
        /* V0 marginal: negative side */                                                                                \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_NegEtaV0", leadPPt, ringObservableLeadP)             \
        /* Joint: (-leadP, -V0) */                                                                                      \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_NegEtaLeadP_NegEtaV0", leadPPt, ringObservableLeadP) \
      }                                                                                                                 \
    }                                                                                                                   \
  } while (0)

// For subleading jet:
#define RING_OBSERVABLE_2NDJET_FILL_LIST(X, FOLDER)                                                                  \
  X(FOLDER "/QA/hDeltaPhi2ndJet", deltaPhi2ndJet)                                                                    \
  X(FOLDER "/QA/hDeltaTheta2ndJet", deltaTheta2ndJet)                                                                \
  X(FOLDER "/QA/hCosDeltaTheta2ndJet", cosDeltaTheta2ndJet)                                                          \
  X(FOLDER "/QA/hPt2ndJet", subleadingJetPt)                                                                         \
  X(FOLDER "/hRingObservable2ndJetCounter", ringObservable2ndJet)                                                    \
  X(FOLDER "/pRingObservable2ndJetDeltaPhi", deltaPhi2ndJet, ringObservable2ndJet)                                   \
  X(FOLDER "/pRingObservable2ndJetDeltaTheta", deltaTheta2ndJet, ringObservable2ndJet)                               \
  X(FOLDER "/EtaDependence/pRingObservableEtaLambda2ndJet", v0eta, ringObservable2ndJet)                             \
  X(FOLDER "/EtaDependence/pRingObservableEta2ndJet", subleadingJetEta, ringObservable2ndJet)                        \
  X(FOLDER "/pRingObservable2ndJetIntegrated", 0., ringObservable2ndJet)                                             \
  X(FOLDER "/pRingObservable2ndJetLambdaPt", v0pt, ringObservable2ndJet)                                             \
  X(FOLDER "/pRingObservableSubLeadPVz", collisionPVz, ringObservable2ndJet)                                         \
  X(FOLDER "/ProxyPtDependence/pRingVsPt2ndJet", subleadingJetPt, ringObservable2ndJet)                              \
  X(FOLDER "/ProxyPtDependence/pRingVsPt2ndJetVsEta2ndJet", subleadingJetPt, subleadingJetEta, ringObservable2ndJet) \
  X(FOLDER "/ProxyPtDependence/pRingVsPt2ndJetVsEtaV0", subleadingJetPt, v0eta, ringObservable2ndJet)                \
  X(FOLDER "/EtaDependence/p2dRingObservableEtaLambdaVsEta2ndJet", v0eta, subleadingJetEta, ringObservable2ndJet)    \
  X(FOLDER "/EtaDependence/h2dCounterEtaLambdaVsEta2ndJet", v0eta, subleadingJetEta)

#define POLARIZATION_PROFILE_FILL_LIST(X, FOLDER)                          \
  /* =============================== */                                    \
  /* 1D TProfiles vs v0phi */                                              \
  /* =============================== */                                    \
  X(FOLDER "/QA/pPxStarPhi", v0phiToFillHists, polStarX)                   \
  X(FOLDER "/QA/pPyStarPhi", v0phiToFillHists, polStarY)                   \
  X(FOLDER "/QA/pPzStarPhi", v0phiToFillHists, polStarZ)                   \
  /* =============================== */                                    \
  /* 1D TProfiles vs DeltaPhi_jet */                                       \
  /* =============================== */                                    \
  X(FOLDER "/QA/pPxStarDeltaPhi", deltaPhiJet, polStarX)                   \
  X(FOLDER "/QA/pPyStarDeltaPhi", deltaPhiJet, polStarY)                   \
  X(FOLDER "/QA/pPzStarDeltaPhi", deltaPhiJet, polStarZ)                   \
  /* =============================== */                                    \
  /* 2D TProfiles vs DeltaPhi_jet and Lambda pT */                         \
  /* =============================== */                                    \
  X(FOLDER "/QA/p2dPxStarDeltaPhiVsLambdaPt", deltaPhiJet, v0pt, polStarX) \
  X(FOLDER "/QA/p2dPyStarDeltaPhiVsLambdaPt", deltaPhiJet, v0pt, polStarY) \
  X(FOLDER "/QA/p2dPzStarDeltaPhiVsLambdaPt", deltaPhiJet, v0pt, polStarZ) \
  /* 2D vector-field profiles vs Lambda momentum plane */                  \
  X(FOLDER "/QA/p2dPxStar_vsPxPy", v0px, v0py, polStarX)                   \
  X(FOLDER "/QA/p2dPyStar_vsPxPy", v0px, v0py, polStarY)                   \
  X(FOLDER "/QA/p2dPzStar_vsPxPy", v0px, v0py, polStarZ)                   \
  X(FOLDER "/QA/p2dPxStar_vsPzPx", v0pz, v0px, polStarX)                   \
  X(FOLDER "/QA/p2dPyStar_vsPzPx", v0pz, v0px, polStarY)                   \
  X(FOLDER "/QA/p2dPzStar_vsPzPx", v0pz, v0px, polStarZ)

// Apply the macros (notice I had to include the semicolon (";") after the function, so you don't need to
// write that when calling this APPLY_HISTO_FILL. The code will look weird, but without this the compiler
// would not know to end each statement with a semicolon):
#define APPLY_HISTO_FILL(NAME, ...) histos.fill(HIST(NAME), __VA_ARGS__);

// Delta Method Fill Lists
#define DELTA_INTEGRATED_FILL_LIST(X, FOLDER, r, n)              \
  X(FOLDER "/DeltaMethod/hIntegrated", 0.5, r)                   \
  X(FOLDER "/DeltaMethod/hIntegrated", 1.5, (double)(n))         \
  X(FOLDER "/DeltaMethod/hIntegrated", 2.5, (r) * (r))           \
  X(FOLDER "/DeltaMethod/hIntegrated", 3.5, (double)((n) * (n))) \
  X(FOLDER "/DeltaMethod/hIntegrated", 4.5, (r) * (n))

#define DELTA_2D_FILL_LIST(X, FOLDER, HIST_NAME, center, r, n)          \
  X(FOLDER "/DeltaMethod/" HIST_NAME, center, 0.5, r)                   \
  X(FOLDER "/DeltaMethod/" HIST_NAME, center, 1.5, (double)(n))         \
  X(FOLDER "/DeltaMethod/" HIST_NAME, center, 2.5, (r) * (r))           \
  X(FOLDER "/DeltaMethod/" HIST_NAME, center, 3.5, (double)((n) * (n))) \
  X(FOLDER "/DeltaMethod/" HIST_NAME, center, 4.5, (r) * (n))

// Master flush macro to dump an event tracker into the histograms:
#define FLUSH_DELTA_TRACKER(FOLDER, TRACKER, AXIS_PT, AXIS_MASS, AXIS_DTHETA)                    \
  if ((TRACKER).nInt > 0) {                                                                      \
    DELTA_INTEGRATED_FILL_LIST(APPLY_HISTO_FILL, FOLDER, (TRACKER).rInt, (TRACKER).nInt)         \
  }                                                                                              \
  for (size_t bin = 0; bin < (TRACKER).rPt.size(); ++bin) {                                      \
    int nVal = (TRACKER).nPt[bin];                                                               \
    if (nVal == 0)                                                                               \
      continue;                                                                                  \
    double rVal = (TRACKER).rPt[bin];                                                            \
    double center = (AXIS_PT)->GetBinCenter(bin);                                                \
    DELTA_2D_FILL_LIST(APPLY_HISTO_FILL, FOLDER, "h2dLambdaPtVsDeltaComp", center, rVal, nVal)   \
  }                                                                                              \
  for (size_t bin = 0; bin < (TRACKER).rMass.size(); ++bin) {                                    \
    int nVal = (TRACKER).nMass[bin];                                                             \
    if (nVal == 0)                                                                               \
      continue;                                                                                  \
    double rVal = (TRACKER).rMass[bin];                                                          \
    double center = (AXIS_MASS)->GetBinCenter(bin);                                              \
    DELTA_2D_FILL_LIST(APPLY_HISTO_FILL, FOLDER, "h2dMassVsDeltaComp", center, rVal, nVal)       \
  }                                                                                              \
  for (size_t bin = 0; bin < (TRACKER).rDtheta.size(); ++bin) {                                  \
    int nVal = (TRACKER).nDtheta[bin];                                                           \
    if (nVal == 0)                                                                               \
      continue;                                                                                  \
    double rVal = (TRACKER).rDtheta[bin];                                                        \
    double center = (AXIS_DTHETA)->GetBinCenter(bin);                                            \
    DELTA_2D_FILL_LIST(APPLY_HISTO_FILL, FOLDER, "h2dDeltaThetaVsDeltaComp", center, rVal, nVal) \
  }

struct lambdajetpolarizationionsderived {
  // Define histogram registries:
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Master analysis switches:
  Configurable<bool> analyseLambda{"analyseLambda", true, "process Lambda-like candidates"};
  Configurable<bool> analyseAntiLambda{"analyseAntiLambda", false, "process AntiLambda-like candidates"};
  Configurable<bool> analyseMagField{"analyseMagField", true, "analyse efficiency effects wrt magnetic field"}; // Older DerivedData lacks magField. "if constexpr (requires { collision.magField(); })" would only see the current header definition, so need a flag for retro-comp.
  // Configurable<bool> doPPAnalysis{"doPPAnalysis", false, "if in pp, set to true. Default is HI"};
  Configurable<bool> doJetProxy5dQA{"doJetProxy5dQA", false, "generates expensive THnSparse histograms for joint distribution QA of the jet proxies and collisions"};

  // A very inexpensive "signal extraction" imitation based on v0InMassPeak bool:
  // (Uses a mass interval to remove or include V0s from the final AnalysisResults to take advantage of existing post-processing codes)
  Configurable<bool> excludeOutOfPeakQA{"excludeOutOfPeakQA", false, "removes all V0s outside an approximate +/- 5*sigma window from the mass peak"}; // A naive estimator of signal
  Configurable<bool> excludeInPeakQA{"excludeInPeakQA", false, "uses only the V0s outside an approximate +/- 5*sigma window from the mass peak"};     // A naive estimator of background

  // Per-family histogram switches:
  // (Each family books >100 histograms, so it is necessary to keep some of these switches off to avoid the HistogramRegistry limit)
  struct : ConfigurableGroup {
    std::string prefix = "familySwitches"; // JSON group name
    Configurable<bool> doFamilyRing{"doFamilyRing", true, "Book and fill the 'Ring' family (no additional cuts). Keep this on for most passes."};
    Configurable<bool> doFamilyRingKinematicCuts{"doFamilyRingKinematicCuts", false, "Book and fill the 'RingKinematicCuts' family (Lambda kinematic cuts applied)."};
    Configurable<bool> doFamilyJetKinematicCuts{"doFamilyJetKinematicCuts", true, "Book and fill the 'JetKinematicCuts' family (jet kinematic cuts applied)."};
    Configurable<bool> doFamilyJetAndLambdaKinematicCuts{"doFamilyJetAndLambdaKinematicCuts", false, "Book and fill the 'JetAndLambdaKinematicCuts' family (both cuts applied)."};
  } familySwitches;

  // QA switches:
  struct : ConfigurableGroup {
    std::string prefix = "qaSwitches";                                                                                                           // JSON group name
    Configurable<bool> doFakePolDiagnosticsQA{"doFakePolDiagnosticsQA", true, "Book and fill the EtaStudy/ and HelicityEfficiencyQA/ folders."}; // The largest per-V0 fill cost in this task
  } qaSwitches;

  // Centrality:
  Configurable<int> centralityEstimator{"centralityEstimator", kCentFT0M, "Run 3 centrality estimator (0:CentFT0C, 1:CentFT0M, 2:CentFV0A)"}; // Default is FT0M
  Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10., "max Z vtx position [cm]"};                                                     // An additional post-processing cut after derived data was written. Same default as lambdaJetPolarizationIons.cxx producer

  // QAs that purposefully "break" the analysis
  // -- All of these tests should give us zero signal if the source is truly Lambda Polarization from vortices
  struct : ConfigurableGroup {
    std::string prefix = "fakePolSwitches"; // JSON group name
    Configurable<bool> forcePolSignQA{"forcePolSignQA", false, "force antiLambda decay constant to be positive: should kill all the signal, if any. For QA"};
    Configurable<bool> forcePerpToJet{"forcePerpToJet", false, "force jet direction to be perpendicular to jet estimator. For QA"};
    Configurable<bool> forceJetDirectionSmudge{"forceJetDirectionSmudge", false, "fluctuate jet direction by 10% of R around original axis. For QA (tests sensibility)"};
    Configurable<bool> forceRandJet{"forceRandJet", false, "makes jet direction random. A QA for AEE fake signal and its removal"};
    Configurable<bool> forcePreviousJet{"forcePreviousJet", false, "uses previous event's jet direction instead of a random sample. A baseline for fake signal removal"};
    Configurable<bool> forceDatalikeJet{"forceDatalikeJet", false, "a compromise between forceRandJet and forcePreviousJet. Parameterized distribution from data"};
    Configurable<bool> doMixedEventProxies{"doMixedEventProxies", false, "mix leadP/leadJet/subJet directions between events using (proxy pt, Zvtx, centrality) bins -- three independent mixings, one per proxy"};
    Configurable<int> mixedEventWindowSize{"mixedEventWindowSize", 10, "number of neighbours for doMixedEventProxies: how many similar collisions stay eligible as mixing partners at once (shared by all 3 proxies)."};
    Configurable<int> nProxyResamples{"nProxyResamples", 1, "The amount of resamplings of jet direction per event. Use ONLY for forceRandJet and forceDatalikeJet"};
  } fakePolSwitches;

  // Configurable<float> jetRForSmudging{"jetRForSmudging", 0.4, "QA quantity: the chosen R scale for the jet direction smudge"}; // Superseeded by jetR: kept the same scale in analysis and QA
  Configurable<float> jetR{"jetR", 0.4f, "Radius of the jet"}; // Provide manually, please.
  Configurable<float> minLeadParticlePt{"minLeadParticlePt", 2.0f, "Minimum Pt for a lead track to be considered a valid proxy for a jet (may be more restrictive than TableProducer)"};
  Configurable<float> minLeadJetPt{"minLeadJetPt", 10.0f, "Minimum Pt for leading jet to be considered valid (may be more restrictive than TableProducer)"};
  Configurable<float> minSubLeadJetPt{"minSubLeadJetPt", 5.0f, "Minimum Pt for subleading jet to be considered valid (may be more restrictive than TableProducer)"};

  /////////////////////////
  // Configurable blocks:
  // Histogram axes configuration:
  struct : ConfigurableGroup {
    std::string prefix = "axisConfigurations"; // JSON group name
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
    ConfigurableAxis axisPtCoarseQA{"axisPtCoarseQA", {VARIABLE_WIDTH, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 7.0f, 10.0f, 15.0f}, "pt axis for QA"};
    ConfigurableAxis axisLambdaMass{"axisLambdaMass", {450, 1.08f, 1.15f}, "#Lambda mass in GeV/c"}; // Default is {200, 1.101f, 1.131f}

    // Symmetric momentum-plane axes for the vector-field / ring 2D profiles:
    ConfigurableAxis axisLambdaPx{"axisLambdaPx", {40, -3.0f, 3.0f}, "#Lambda p_{x} (GeV/c)"};
    ConfigurableAxis axisLambdaPy{"axisLambdaPy", {40, -3.0f, 3.0f}, "#Lambda p_{y} (GeV/c)"};
    ConfigurableAxis axisLambdaPz{"axisLambdaPz", {40, -4.0f, 4.0f}, "#Lambda p_{z} (GeV/c)"};

    // Event properties:
    ConfigurableAxis axisPVz{"axisPVz", {60, -15.0f, +15.0f}, "Primary Vertex Z [cm]"};

    // Jet axes:
    // ConfigurableAxis axisLeadingParticlePt{"axisLeadingParticlePt", {100, 0.f, 200.f}, "Leading particle p_{T} (GeV/c)"}; // Simpler version!
    // ConfigurableAxis axisJetPt{"axisJetPt", {50, 0.f, 200.f}, "Jet p_{t} (GeV)"};
    ConfigurableAxis axisJetPt{
      "axisJetPt",
      {VARIABLE_WIDTH,
       0, 2, 4, 6, 8, 10, // 2 GeV bins
       15, 20,            // 5 GeV bins
       30, 40,            // 10 GeV bins
       60, 80,            // 20 GeV bins
       120, 160, 200},    // 40 GeV bins
      "Jet p_{T} (GeV)"};
    ConfigurableAxis axisJetPtSigExtract{"axisJetPtSigExtract", {VARIABLE_WIDTH, 0, 5, 10, 12, 16, 20, 25, 30, 35, 40, 60, 100, 200}, "Jet p_{t} (GeV)"};
    ConfigurableAxis axisEta{"axisEta", {50, -1.0f, 1.0f}, "#eta"};
    ConfigurableAxis axisEtaCoarse{"axisEtaCoarse", {20, -0.9f, 0.9f}, "#eta coarse axis"};
    ConfigurableAxis axisV0Eta{"axisV0Eta", {75, -1.5f, 1.5f}, "V0 #eta"}; // An axis for V0 eta, which can go up to 1.5 given standard producer selections
    ConfigurableAxis axisV0EtaCoarse{"axisV0EtaCoarse", {32, -1.5f, 1.5f}, "V0 #eta coarse"};
    ConfigurableAxis axisDeltaEtaCoarse{"axisDeltaEtaCoarse", {40, -1.8f, 1.8f}, "#Delta#eta coarse axis"};
    ConfigurableAxis axisDeltaTheta{"axisDeltaTheta", {40, 0, constants::math::PI}, "#Delta #theta_{jet}"};
    ConfigurableAxis axisCosTheta{"axisCosTheta", {50, -1, 1}, "cos(#theta)"};
    ConfigurableAxis axisPhi{"axisPhi", {40, 0., constants::math::TwoPI}, "#varphi"};
    ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {40, -constants::math::PI, constants::math::PI}, "#Delta #phi_{jet}"};
    ConfigurableAxis axisRingCounts{"axisRingCounts", {90, -4.5, 4.5}, "<#it{R}>"};
    ConfigurableAxis axisDeltaCollisionIndex{"axisDeltaCollisionIndex", {200, -0.5f, 199.5f}, "#Delta collision index"}; // Always positive: SameKindPair pairs strictly upper

    ConfigurableAxis axisDCAdau{"axisDCAdau", {10, 0., 2.0}, "DCA V0 daughters (cm)"};
    ConfigurableAxis axisDCAdauPV{"axisDCAdauPV", {10, 0., 1.2}, "DCA dauPV (cm)"}; // v0Selections.dcav0dau's default maximum is 1.2f in the TableProducer

    // Coarser axes for signal extraction:
    ConfigurableAxis axisPtSigExtract{"axisPtSigExtract", {VARIABLE_WIDTH, 0.0f, 0.25f, 0.5f, 0.75f, 1.0f, 1.25f, 1.5f, 2.0f, 2.5f, 3.0f, 4.0f, 6.0f, 8.0f, 10.0f, 15.0f, 20.0f, 30.0f, 50.0f}, "pt axis for signal extraction"};
    // ConfigurableAxis axisLambdaMassSigExtract{"axisLambdaMassSigExtract", {175, 1.08f, 1.15f}, "Lambda mass in GeV/c"}; // With a sigma of 0.002 GeV/c, this has about 5 bins per sigma, so that the window is properly grasped.
    // A coarser axis (sigma is still well estimated, with about 8 bins in the peak region)
    ConfigurableAxis axisLambdaMassSigExtract{
      "axisLambdaMassSigExtract",
      {VARIABLE_WIDTH,
       // Left sideband (7 bins, 0.004 width)
       1.0800, 1.0840, 1.0880, 1.0920,
       1.0960, 1.1000, 1.1040, 1.1080,
       // Fine peak region (8 bins, 0.0016 width)
       1.1096, 1.1112, 1.1128, 1.1144,
       1.1160, 1.1176, 1.1192, 1.1208,
       // Right sideband (7 bins, 0.004 width)
       1.1248, 1.1288, 1.1328, 1.1368,
       1.1408, 1.1448, 1.1488},
      "Lambda mass in GeV/c"};
    // ConfigurableAxis axisLeadingParticlePtSigExtract{"axisLeadingParticlePtSigExtract", {VARIABLE_WIDTH, 0, 4, 8, 12, 16, 20, 25, 30, 35, 40, 60, 100, 200}, "Leading particle p_{T} (GeV/c)"}; // Simpler version!

    // (TODO: add a lambdaPt axis that is pre-selected only on the 0.5 to 1.5 Pt region for the Ring observable with lambda cuts to not store a huge histogram with empty bins by construction)

    // ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Centrality"};
    ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 20.0f, 50.0f, 100.0f}, "Centrality"};

    // For the delta method error propagation (slightly better than just SEM error propagation with TProfiles):
    ConfigurableAxis axisDeltaComponents{"axisDeltaComponents", {5, 0.0, 5.0}, "0: r_k, 1: n_k, 2: r_k^2, 3: n_k^2, 4: r_k*n_k"};
  } axisConfigurations;

  // Helper functions:
  // Fast wrapping into [-PI, PI) (restricted to this interval for function speed)
  inline double wrapToPiFast(double phi)
  {
    constexpr double TwoPi = o2::constants::math::TwoPI;
    constexpr double Pi = o2::constants::math::PI;
    if (phi >= Pi)
      phi -= TwoPi;
    else if (phi < -Pi)
      phi += TwoPi;
    return phi;
  }

  // A small tracker struct for convenience -- Accumulates values for the Delta Method error estimator:
  struct EventDeltaTracker {
    double rInt = 0.0; // Ring accumulator
    int nInt = 0;      // Counts accumulator
    std::vector<double> rPt, rMass, rDtheta;
    std::vector<int> nPt, nMass, nDtheta;

    /// \brief Resizes every accumulator. Size includes ROOT's under/overflow bins, so the indices
    ///        returned by TAxis::FindBin() (0 .. nBins+1) can be used directly for dereferencing here.
    void resize(int nBinsPt, int nBinsMass, int nBinsDTheta)
    {
      rPt.assign(nBinsPt + 2, 0.0);
      rMass.assign(nBinsMass + 2, 0.0);
      rDtheta.assign(nBinsDTheta + 2, 0.0);
      nPt.assign(nBinsPt + 2, 0);
      nMass.assign(nBinsMass + 2, 0);
      nDtheta.assign(nBinsDTheta + 2, 0);
    }

    void reset()
    {
      rInt = 0.0;
      nInt = 0;
      std::fill(rPt.begin(), rPt.end(), 0.0);
      std::fill(rMass.begin(), rMass.end(), 0.0);
      std::fill(rDtheta.begin(), rDtheta.end(), 0.0);
      std::fill(nPt.begin(), nPt.end(), 0);
      std::fill(nMass.begin(), nMass.end(), 0);
      std::fill(nDtheta.begin(), nDtheta.end(), 0);
    }

    void addV0(double ringObs, int binPt, int binMass, int binDTheta)
    {
      rInt += ringObs;
      rPt[binPt] += ringObs;
      rMass[binMass] += ringObs;
      rDtheta[binDTheta] += ringObs;
      nInt += 1;
      nPt[binPt] += 1;
      nMass[binMass] += 1;
      nDtheta[binDTheta] += 1;
    }
  };

  // Allocating one tracker per family so the accumulators are allocated only once:
  EventDeltaTracker trackRing, trackRingKinCuts, trackJetKinCuts, trackJetLambdaKinCuts;

  // Axis pointers for Delta Method binning (fetched once in init, declared once here)
  TAxis* mAxisPt = nullptr;
  TAxis* mAxisMass = nullptr;
  TAxis* mAxisDTheta = nullptr;

  void init(InitContext const&)
  {
    // Configuration validation:
    // These combinations would not crash otherwise, so they are rejected at init() time.
    const int nDistortionsOn = static_cast<int>(fakePolSwitches.forcePerpToJet) + static_cast<int>(fakePolSwitches.forceJetDirectionSmudge) +
                               static_cast<int>(fakePolSwitches.forceRandJet) + static_cast<int>(fakePolSwitches.forcePreviousJet) +
                               static_cast<int>(fakePolSwitches.forceDatalikeJet) + static_cast<int>(fakePolSwitches.doMixedEventProxies);
    if (nDistortionsOn > 1) // applyProxyDistortion() is an if/else chain, so extra switches are silently ignored
      LOG(fatal) << "fakePolSwitches: " << nDistortionsOn << " proxy distortions enabled at once. They are mutually exclusive -- "
                 << "applyProxyDistortion() would apply only the first and silently drop the rest.";
    if (fakePolSwitches.nProxyResamples > 1 && (fakePolSwitches.forcePreviousJet || fakePolSwitches.doMixedEventProxies))
      LOG(fatal) << "fakePolSwitches: nProxyResamples > 1 is only meaningful for forceRandJet/forceDatalikeJet. "
                 << "Previous-jet/Mixed-Event proxies do not change between resamplings, so every extra pass would double-count the same proxy.";
    if (excludeOutOfPeakQA && excludeInPeakQA) // Complementary selections
      LOG(fatal) << "excludeOutOfPeakQA and excludeInPeakQA are complementary: enabling both rejects every V0.";
    if (!familySwitches.doFamilyRing) // Todo: think of a smarter way of handling the axis getters for the DeltaMethod
      LOG(fatal) << "doFamilyRing must be on: the Delta Method accumulators take their binning from the Ring/ histograms.";

    // Ring observable histograms:
    // Helper to register one full histogram family (kinematic cut variation of ring observable)
    auto addRingObservableFamily = [&](const std::string& folder) {
      // ===============================
      // QA histograms: angle and pT distributions
      // (No mass dependency -- useful to check kinematic sculpting from cuts)
      // ===============================
      histos.add((folder + "/QA/hDeltaPhi").c_str(), "#Delta#varphi_{jet};#Delta#varphi_{jet};Counts", kTH1D, {axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/QA/hDeltaPhiVsDeltaEta").c_str(), "#Delta#varphi_{jet};#Delta#varphi_{jet}; #eta_{#Lambda}-#eta_{Jet};Counts", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDeltaEtaCoarse});
      histos.add((folder + "/QA/hDeltaPhiVsLeadJetPhi").c_str(), "#Delta#varphi_{jet};#varphi_{Jet};#varphi_{Jet};Counts", kTH2D, {axisConfigurations.axisPhi, axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/QA/hDeltaTheta").c_str(), "#Delta#theta_{jet};#Delta#theta_{jet};Counts", kTH1D, {axisConfigurations.axisDeltaTheta});
      histos.add((folder + "/QA/hCosDeltaTheta").c_str(), "cos(#Delta#theta_{jet});cos(#Delta#theta_{jet});Counts", kTH1D, {axisConfigurations.axisCosTheta}); // Should actually be flat due to the geometry
      histos.add((folder + "/QA/hIntegrated").c_str(), "Integrated counts; ;Counts", kTH1D, {{1, -0.5, 0.5}});

      histos.add((folder + "/QA/hDeltaPhiLeadP").c_str(), "#Delta#varphi_{LeadP};#Delta#varphi_{LeadP};Counts", kTH1D, {axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/QA/hDeltaThetaLeadP").c_str(), "#Delta#theta_{LeadP};#Delta#theta_{LeadP};Counts", kTH1D, {axisConfigurations.axisDeltaTheta});
      histos.add((folder + "/QA/hCosDeltaThetaLeadP").c_str(), "cos(#Delta#theta_{LeadP});cos(#Delta#theta_{LeadP});Counts", kTH1D, {axisConfigurations.axisCosTheta}); // Should actually be flat due to the geometry
      histos.add((folder + "/QA/hDeltaPhi2ndJet").c_str(), "#Delta#varphi_{SubJet};#Delta#varphi_{SubJet};Counts", kTH1D, {axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/QA/hDeltaTheta2ndJet").c_str(), "#Delta#theta_{SubJet};#Delta#theta_{SubJet};Counts", kTH1D, {axisConfigurations.axisDeltaTheta});
      histos.add((folder + "/QA/hCosDeltaTheta2ndJet").c_str(), "cos(#Delta#theta_{SubJet});cos(#Delta#theta_{SubJet});Counts", kTH1D, {axisConfigurations.axisCosTheta}); // Should actually be flat due to the geometry

      // ===============================
      // Lambda pT dependence
      // ===============================
      histos.add((folder + "/QA/hLambdaPt").c_str(), "#Lambda #it{p}_{T};#it{p}_{T}^{#Lambda} (GeV/c);Counts", kTH1D, {axisConfigurations.axisPt});
      histos.add((folder + "/QA/h2dDeltaPhiVsLambdaPt").c_str(), "#Delta#varphi_{jet} vs #Lambda #it{p}_{T};#Delta#varphi_{jet};#it{p}_{T}^{#Lambda} (GeV/c)", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPt});
      histos.add((folder + "/QA/h2dDeltaThetaVsLambdaPt").c_str(), "#Delta#theta_{jet} vs #Lambda #it{p}_{T};#Delta#theta_{jet};#it{p}_{T}^{#Lambda} (GeV/c)", kTH2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisPt});
      // ===============================
      //   Polarization observable QAs
      // (not Ring: actual polarization!)
      // ===============================
      // Will implement these as TProfiles, as polarization is also a measure like P_Lambda = (3/\alpha_Lambda) * <p_{proton}>, so the error is similar
      // ===============================
      // 1D TProfiles
      // ===============================
      histos.add((folder + "/QA/pPxStarPhi").c_str(), "<P_{#Lambda}^{*}>_{x} vs #varphi_{#Lambda};#varphi_{#Lambda};<P_{#Lambda}^{*}>_{x}", kTProfile, {axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/QA/pPyStarPhi").c_str(), "<P_{#Lambda}^{*}>_{y} vs #varphi_{#Lambda};#varphi_{#Lambda};<P_{#Lambda}^{*}>_{y}", kTProfile, {axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/QA/pPzStarPhi").c_str(), "<P_{#Lambda}^{*}>_{z} vs #varphi_{#Lambda};#varphi_{#Lambda};<P_{#Lambda}^{*}>_{z}", kTProfile, {axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/QA/pPxStarDeltaPhi").c_str(), "<P_{#Lambda}^{*}>_{x} vs #Delta#varphi_{jet};#Delta#varphi_{jet};<P_{#Lambda}^{*}>_{x}", kTProfile, {axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/QA/pPyStarDeltaPhi").c_str(), "<P_{#Lambda}^{*}>_{y} vs #Delta#varphi_{jet};#Delta#varphi_{jet};<P_{#Lambda}^{*}>_{y}", kTProfile, {axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/QA/pPzStarDeltaPhi").c_str(), "<P_{#Lambda}^{*}>_{z} vs #Delta#varphi_{jet};#Delta#varphi_{jet};<P_{#Lambda}^{*}>_{z}", kTProfile, {axisConfigurations.axisDeltaPhi});
      // ===============================
      // 2D TProfiles (Lambda correlations)
      // ===============================
      histos.add((folder + "/QA/p2dPxStarDeltaPhiVsLambdaPt").c_str(), "<P_{#Lambda}^{*}>_{x} vs #Delta#varphi_{jet} vs #it{p}_{T}^{#Lambda};#Delta#varphi_{jet};#it{p}_{T}^{#Lambda} (GeV/c);<P_{#Lambda}^{*}>_{x}", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPtSigExtract});
      histos.add((folder + "/QA/p2dPyStarDeltaPhiVsLambdaPt").c_str(), "<P_{#Lambda}^{*}>_{y} vs #Delta#varphi_{jet} vs #it{p}_{T}^{#Lambda};#Delta#varphi_{jet};#it{p}_{T}^{#Lambda} (GeV/c);<P_{#Lambda}^{*}>_{y}", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPtSigExtract});
      histos.add((folder + "/QA/p2dPzStarDeltaPhiVsLambdaPt").c_str(), "<P_{#Lambda}^{*}>_{z} vs #Delta#varphi_{jet} vs #it{p}_{T}^{#Lambda};#Delta#varphi_{jet};#it{p}_{T}^{#Lambda} (GeV/c);<P_{#Lambda}^{*}>_{z}", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPtSigExtract});

      // ===============================
      // Vector-field profiles: <P*_x,y,z> vs Lambda momentum plane:
      // (a companion code plots these as a vector field)
      // ===============================
      histos.add((folder + "/QA/p2dPxStar_vsPxPy").c_str(), "<P_{#Lambda}^{*}>_{x} vs (p_{x}^{#Lambda},p_{y}^{#Lambda});p_{x}^{#Lambda} (GeV/c);p_{y}^{#Lambda} (GeV/c);<P_{#Lambda}^{*}>_{x}", kTProfile2D, {axisConfigurations.axisLambdaPx, axisConfigurations.axisLambdaPy});
      histos.add((folder + "/QA/p2dPyStar_vsPxPy").c_str(), "<P_{#Lambda}^{*}>_{y} vs (p_{x}^{#Lambda},p_{y}^{#Lambda});p_{x}^{#Lambda} (GeV/c);p_{y}^{#Lambda} (GeV/c);<P_{#Lambda}^{*}>_{y}", kTProfile2D, {axisConfigurations.axisLambdaPx, axisConfigurations.axisLambdaPy});
      histos.add((folder + "/QA/p2dPzStar_vsPxPy").c_str(), "<P_{#Lambda}^{*}>_{z} vs (p_{x}^{#Lambda},p_{y}^{#Lambda}) [colormap];p_{x}^{#Lambda} (GeV/c);p_{y}^{#Lambda} (GeV/c);<P_{#Lambda}^{*}>_{z}", kTProfile2D, {axisConfigurations.axisLambdaPx, axisConfigurations.axisLambdaPy});

      histos.add((folder + "/QA/p2dPxStar_vsPzPx").c_str(), "<P_{#Lambda}^{*}>_{x} vs (p_{z}^{#Lambda},p_{x}^{#Lambda});p_{z}^{#Lambda} (GeV/c);p_{x}^{#Lambda} (GeV/c);<P_{#Lambda}^{*}>_{x}", kTProfile2D, {axisConfigurations.axisLambdaPz, axisConfigurations.axisLambdaPx});
      histos.add((folder + "/QA/p2dPyStar_vsPzPx").c_str(), "<P_{#Lambda}^{*}>_{y} vs (p_{z}^{#Lambda},p_{x}^{#Lambda}) [colormap];p_{z}^{#Lambda} (GeV/c);p_{x}^{#Lambda} (GeV/c);<P_{#Lambda}^{*}>_{y}", kTProfile2D, {axisConfigurations.axisLambdaPz, axisConfigurations.axisLambdaPx});
      histos.add((folder + "/QA/p2dPzStar_vsPzPx").c_str(), "<P_{#Lambda}^{*}>_{z} vs (p_{z}^{#Lambda},p_{x}^{#Lambda});p_{z}^{#Lambda} (GeV/c);p_{x}^{#Lambda} (GeV/c);<P_{#Lambda}^{*}>_{z}", kTProfile2D, {axisConfigurations.axisLambdaPz, axisConfigurations.axisLambdaPx});

      // ===============================
      // Ring observable, single scalar 2D profile:
      // ===============================
      histos.add((folder + "/p2dRingObservableVsPxPy").c_str(), "<#it{R}> vs (p_{x}^{#Lambda},p_{y}^{#Lambda});p_{x}^{#Lambda} (GeV/c);p_{y}^{#Lambda} (GeV/c);<#it{R}>", kTProfile2D, {axisConfigurations.axisLambdaPx, axisConfigurations.axisLambdaPy});
      histos.add((folder + "/p2dRingObservableVsPzPx").c_str(), "<#it{R}> vs (p_{z}^{#Lambda},p_{x}^{#Lambda});p_{z}^{#Lambda} (GeV/c);p_{x}^{#Lambda} (GeV/c);<#it{R}>", kTProfile2D, {axisConfigurations.axisLambdaPz, axisConfigurations.axisLambdaPx});
      // For LeadP estimators:
      histos.add((folder + "/p2dRingObservableLeadPVsPxPy").c_str(), "<#it{R}>_{LeadP} vs (p_{x}^{#Lambda},p_{y}^{#Lambda});p_{x}^{#Lambda} (GeV/c);p_{y}^{#Lambda} (GeV/c);<#it{R}>", kTProfile2D, {axisConfigurations.axisLambdaPx, axisConfigurations.axisLambdaPy});
      histos.add((folder + "/p2dRingObservableLeadPVsPzPx").c_str(), "<#it{R}>_{LeadP} vs (p_{z}^{#Lambda},p_{x}^{#Lambda});p_{z}^{#Lambda} (GeV/c);p_{x}^{#Lambda} (GeV/c);<#it{R}>", kTProfile2D, {axisConfigurations.axisLambdaPz, axisConfigurations.axisLambdaPx});

      // TProfiles with correct error bars::
      // -- TProfiles will handle the error estimate of the Ring Observable via the variance, even though
      // they still lack the proper signal extraction and possible efficiency corrections in the current state
      // -- If any efficiency corrections arise, you can fill with the kTH1D as (deltaPhiJet, ringObservable, weight)
      // instead of the simple (deltaPhiJet, ringObservable) --> Notice TProfile knows how to accept 3 entries
      // for a TH1D-like object!
      // -- CAUTION! The TProfile does not utilize unbiased variance estimators with N-1 instead of N in the denominator,
      // so you might get biased errors when counts are too low in higher-dimensional profiles (i.e., kTProfile2Ds)
      // ===============================
      // 1D TProfiles
      // ===============================
      histos.add((folder + "/pRingObservableDeltaPhi").c_str(), "<#it{R}> vs #Delta#varphi_{jet};#Delta#varphi_{jet};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaPhi});
      // To see the actual distribution of counts in data (another differential-like shape of the distribution we are taking an average of):
      histos.add((folder + "/hRingObservableCounts").c_str(), "Counts vs <#it{R}>_{jet};<#it{R}>; Counts", kTH1D, {axisConfigurations.axisRingCounts});
      histos.add((folder + "/pRingObservablePhiJet").c_str(), "<#it{R}> vs #varphi_{jet};#varphi_{jet};<#it{R}>", kTProfile, {axisConfigurations.axisPhi});
      histos.add((folder + "/pRingObservablePhiLambda").c_str(), "<#it{R}> vs #varphi_{#Lambda};#varphi_{#Lambda};<#it{R}>", kTProfile, {axisConfigurations.axisPhi});
      histos.add((folder + "/pRingObservableDeltaTheta").c_str(), "<#it{R}> vs #Delta#theta_{jet};#Delta#theta_{jet};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaTheta});
      histos.add((folder + "/pRingObservableIntegrated").c_str(), "Integrated <#it{R}>; ;<#it{R}>", kTProfile, {{1, -0.5, 0.5}});
      histos.add((folder + "/pRingObservableLambdaPt").c_str(), "<#it{R}> vs #it{p}_{T}^{#Lambda};#it{p}_{T}^{#Lambda} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisPt});

      // Ring vs Jet proxy pT:
      histos.add((folder + "/ProxyPtDependence/pRingVsPtJet").c_str(), "<#it{R}> vs Jet #it{p}_{T};#it{p}_{T}^{Jet} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisJetPt});
      histos.add((folder + "/ProxyPtDependence/pRingVsPtLeadP").c_str(), "<#it{R}> vs LeadP #it{p}_{T};#it{p}_{T}^{LeadP} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisJetPt});
      histos.add((folder + "/ProxyPtDependence/pRingVsPt2ndJet").c_str(), "<#it{R}> vs SubJet #it{p}_{T};#it{p}_{T}^{SubJet} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisJetPt});
      // And some counters to be aware of the amount of Lambdas (and jets) in each pT interval:
      histos.add((folder + "/QA/hPtJet").c_str(), "Jet #it{p}_{T};#it{p}_{T}^{Jet} (GeV/c);Counts", kTH1D, {axisConfigurations.axisJetPt});
      histos.add((folder + "/QA/hPtLeadP").c_str(), "LeadP #it{p}_{T};#it{p}_{T}^{LeadP} (GeV/c);Counts", kTH1D, {axisConfigurations.axisJetPt});
      histos.add((folder + "/QA/hPt2ndJet").c_str(), "SubJet #it{p}_{T};#it{p}_{T}^{SubJet} (GeV/c);Counts", kTH1D, {axisConfigurations.axisJetPt});

      // Splitting into positive and negative eta contributions:
      histos.add((folder + "/ProxyPtDependence/pRingVsPtJetVsEtaJet").c_str(), "<#it{R}> vs Jet #it{p}_{T} vs #eta_{Jet};#it{p}_{T}^{Jet} (GeV/c);#eta_{Jet};<#it{R}>", kTProfile2D, {axisConfigurations.axisJetPt, {2, -0.9, 0.9}});
      histos.add((folder + "/ProxyPtDependence/pRingVsPtLeadPVsEtaLeadP").c_str(), "<#it{R}> vs LeadP #it{p}_{T} vs #eta_{LeadP};#it{p}_{T}^{LeadP} (GeV/c);#eta_{LeadP};<#it{R}>", kTProfile2D, {axisConfigurations.axisJetPt, {2, -0.9, 0.9}});
      histos.add((folder + "/ProxyPtDependence/pRingVsPt2ndJetVsEta2ndJet").c_str(), "<#it{R}> vs SubJet #it{p}_{T} vs #eta_{SubJet};#it{p}_{T}^{SubJet} (GeV/c);#eta_{SubJet};<#it{R}>", kTProfile2D, {axisConfigurations.axisJetPt, {2, -0.9, 0.9}});
      // For each Lambda's eta:
      histos.add((folder + "/ProxyPtDependence/pRingVsPtJetVsEtaV0").c_str(), "<#it{R}> vs Jet #it{p}_{T} vs #eta_{V0};#it{p}_{T}^{Jet} (GeV/c);#eta_{V0};<#it{R}>", kTProfile2D, {axisConfigurations.axisJetPt, {2, -0.9, 0.9}});
      histos.add((folder + "/ProxyPtDependence/pRingVsPtLeadPVsEtaV0").c_str(), "<#it{R}> vs LeadP #it{p}_{T} vs #eta_{V0};#it{p}_{T}^{LeadP} (GeV/c);#eta_{V0};<#it{R}>", kTProfile2D, {axisConfigurations.axisJetPt, {2, -0.9, 0.9}});
      histos.add((folder + "/ProxyPtDependence/pRingVsPt2ndJetVsEtaV0").c_str(), "<#it{R}> vs SubJet #it{p}_{T} vs #eta_{V0};#it{p}_{T}^{SubJet} (GeV/c);#eta_{V0};<#it{R}>", kTProfile2D, {axisConfigurations.axisJetPt, {2, -0.9, 0.9}});

      // Rasterizing, only for LeadP the TProfile2D into two TProfile 1Ds (easier to draw with "same")
      histos.add((folder + "/ProxyPtDependence/pRingVsPtLeadP_PosEtaLeadP").c_str(), "<#it{R}> vs LeadP #it{p}_{T} (#eta_{LeadP}>0);#it{p}_{T}^{LeadP} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisJetPt});
      histos.add((folder + "/ProxyPtDependence/pRingVsPtLeadP_NegEtaLeadP").c_str(), "<#it{R}> vs LeadP #it{p}_{T} (#eta_{LeadP}<0);#it{p}_{T}^{LeadP} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisJetPt});
      // V0 eta:
      histos.add((folder + "/ProxyPtDependence/pRingVsPtLeadP_PosEtaV0").c_str(), "<#it{R}> vs LeadP #it{p}_{T} (#eta_{V0}>0);#it{p}_{T}^{LeadP} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisJetPt});
      histos.add((folder + "/ProxyPtDependence/pRingVsPtLeadP_NegEtaV0").c_str(), "<#it{R}> vs LeadP #it{p}_{T} (#eta_{V0}<0);#it{p}_{T}^{LeadP} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisJetPt});

      // Looking at V0Eta and JetEta combinations (only for LeadP):
      histos.add((folder + "/ProxyPtDependence/pRingVsPtLeadP_PosEtaLeadP_PosEtaV0").c_str(), "<#it{R}> vs LeadP #it{p}_{T} (#eta_{LeadP}>0, #eta_{V0}>0);#it{p}_{T}^{LeadP} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisJetPt});
      histos.add((folder + "/ProxyPtDependence/pRingVsPtLeadP_NegEtaLeadP_PosEtaV0").c_str(), "<#it{R}> vs LeadP #it{p}_{T} (#eta_{LeadP}<0, #eta_{V0}>0);#it{p}_{T}^{LeadP} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisJetPt});
      histos.add((folder + "/ProxyPtDependence/pRingVsPtLeadP_PosEtaLeadP_NegEtaV0").c_str(), "<#it{R}> vs LeadP #it{p}_{T} (#eta_{LeadP}>0, #eta_{V0}<0);#it{p}_{T}^{LeadP} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisJetPt});
      histos.add((folder + "/ProxyPtDependence/pRingVsPtLeadP_NegEtaLeadP_NegEtaV0").c_str(), "<#it{R}> vs LeadP #it{p}_{T} (#eta_{LeadP}<0, #eta_{V0}<0);#it{p}_{T}^{LeadP} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisJetPt});

      // Understanding eta dependence seen in pRingEtaCuts:
      histos.add((folder + "/EtaDependence/pRingObservableEtaLambda").c_str(), "<#it{R}> vs #eta_{#Lambda};#eta_{#Lambda};<#it{R}>", kTProfile, {axisConfigurations.axisV0EtaCoarse});
      histos.add((folder + "/EtaDependence/pRingObservableEtaJet").c_str(), "<#it{R}> vs #eta_{Jet};#eta_{Jet};<#it{R}>", kTProfile, {axisConfigurations.axisEtaCoarse});

      histos.add((folder + "/EtaDependence/pRingObservableEtaLambda2ndJet").c_str(), "<#it{R}> vs #eta_{#Lambda} (SubJet);#eta_{#Lambda};<#it{R}>", kTProfile, {axisConfigurations.axisV0EtaCoarse});
      histos.add((folder + "/EtaDependence/pRingObservableEta2ndJet").c_str(), "<#it{R}> vs #eta_{SubJet};#eta_{SubJet};<#it{R}>", kTProfile, {axisConfigurations.axisEtaCoarse});

      histos.add((folder + "/EtaDependence/pRingObservableEtaLambdaLeadP").c_str(), "<#it{R}> vs #eta_{#Lambda} (LeadP);#eta_{#Lambda};<#it{R}>", kTProfile, {axisConfigurations.axisV0EtaCoarse});
      histos.add((folder + "/EtaDependence/pRingObservableEtaLeadP").c_str(), "<#it{R}> vs #eta_{LeadP};#eta_{LeadP};<#it{R}>", kTProfile, {axisConfigurations.axisEtaCoarse});
      // For the leading particle:
      histos.add((folder + "/hRingObservableLeadPCounts").c_str(), "Counts vs <#it{R}>_{LeadP};<#it{R}>; Counts", kTH1D, {axisConfigurations.axisRingCounts});
      histos.add((folder + "/pRingObservableLeadPDeltaPhi").c_str(), "<#it{R}> vs #Delta#varphi_{LeadP};#Delta#varphi_{LeadP};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/pRingObservableLeadPDeltaTheta").c_str(), "<#it{R}> vs #Delta#theta_{LeadP};#Delta#theta_{LeadP};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaTheta});
      histos.add((folder + "/pRingObservableLeadPIntegrated").c_str(), "Integrated <#it{R}> (LeadP); ;<#it{R}>", kTProfile, {{1, -0.5, 0.5}});
      histos.add((folder + "/pRingObservableLeadPLambdaPt").c_str(), "<#it{R}> vs #it{p}_{T}^{#Lambda} (LeadP);#it{p}_{T}^{#Lambda} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisPt});
      // For the second-to-leading jet:
      histos.add((folder + "/hRingObservable2ndJetCounter").c_str(), "Counts vs <#it{R}>_{SubJet};<#it{R}>; Counts", kTH1D, {axisConfigurations.axisRingCounts});
      histos.add((folder + "/pRingObservable2ndJetDeltaPhi").c_str(), "<#it{R}> vs #Delta#varphi_{SubJet};#Delta#varphi_{SubJet};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/pRingObservable2ndJetDeltaTheta").c_str(), "<#it{R}> vs #Delta#theta_{SubJet};#Delta#theta_{SubJet};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaTheta});
      histos.add((folder + "/pRingObservable2ndJetIntegrated").c_str(), "Integrated <#it{R}> (SubJet); ;<#it{R}>", kTProfile, {{1, -0.5, 0.5}});
      histos.add((folder + "/pRingObservable2ndJetLambdaPt").c_str(), "<#it{R}> vs #it{p}_{T}^{#Lambda} (SubJet);#it{p}_{T}^{#Lambda} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisPt});

      // For the Zvtx dependence:
      histos.add((folder + "/pRingObservableLeadJetPVz").c_str(), "<#it{R}>_{LeadJet} vs PVz;PVz (cm);<#it{R}>", kTProfile, {axisConfigurations.axisPVz});
      histos.add((folder + "/pRingObservableSubLeadPVz").c_str(), "<#it{R}>_{SubLead} vs PVz;PVz (cm);<#it{R}>", kTProfile, {axisConfigurations.axisPVz});
      histos.add((folder + "/pRingObservableLeadPPVz").c_str(), "<#it{R}>_{LeadP} vs PVz;PVz (cm);<#it{R}>", kTProfile, {axisConfigurations.axisPVz});
      // ===============================
      // 2D TProfiles (Lambda correlations)
      // ===============================
      histos.add((folder + "/p2dRingObservableDeltaPhiVsLambdaPt").c_str(), "<#it{R}> vs #Delta#varphi_{jet} vs #it{p}_{T}^{#Lambda};#Delta#varphi_{jet};#it{p}_{T}^{#Lambda} (GeV/c);<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPt});
      histos.add((folder + "/p2dRingObservableDeltaThetaVsLambdaPt").c_str(), "<#it{R}> vs #Delta#theta_{jet} vs #it{p}_{T}^{#Lambda};#Delta#theta_{jet};#it{p}_{T}^{#Lambda} (GeV/c);<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisPt});
      // ===============================
      // 2D TProfiles (Jet correlations)
      // ===============================
      histos.add((folder + "/p2dRingObservableDeltaPhiVsLeadJetPt").c_str(), "<#it{R}> vs #Delta#varphi_{jet} vs Lead Jet #it{p}_{T};#Delta#varphi_{jet};#it{p}_{T}^{LeadJet} (GeV/c);<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisJetPt});
      histos.add((folder + "/p2dRingObservableDeltaThetaVsLeadJetPt").c_str(), "<#it{R}> vs #Delta#theta_{jet} vs Lead Jet #it{p}_{T};#Delta#theta_{jet};#it{p}_{T}^{LeadJet} (GeV/c);<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisJetPt});

      // ===============================
      // Multi-dimensional histograms for signal extraction
      // (Mass-dependent polarization extraction)
      // ===============================
      // Simple invariant mass plot for QA:
      histos.add((folder + "/QA/hMass").c_str(), "#Lambda Mass;m_{p#pi} (GeV/c^{2});Counts", kTH1D, {axisConfigurations.axisLambdaMass});
      histos.add((folder + "/hMassSigExtract").c_str(), "#Lambda Mass (Sig Extract);m_{p#pi} (GeV/c^{2});Counts", kTH1D, {axisConfigurations.axisLambdaMassSigExtract});
      // 1D Mass dependence of observable numerator:
      histos.add((folder + "/QA/hRingObservableNumMass").c_str(), "Ring Observable Numerator vs Mass;m_{p#pi} (GeV/c^{2});Counts", kTH1D, {axisConfigurations.axisLambdaMassSigExtract});
      // --- 2D counters: Angle vs Mass vs ---
      histos.add((folder + "/QA/h2dDeltaPhiVsMass").c_str(), "#Delta#varphi_{jet} vs Mass;#Delta#varphi_{jet};m_{p#pi} (GeV/c^{2})", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract});
      histos.add((folder + "/QA/h2dDeltaThetaVsMass").c_str(), "#Delta#theta_{jet} vs Mass;#Delta#theta_{jet};m_{p#pi} (GeV/c^{2})", kTH2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract});
      // --- 3D counters: Angle vs Mass vs Lambda pT ---
      histos.add((folder + "/QA/h3dDeltaPhiVsMassVsLambdaPt").c_str(), "#Delta#varphi_{jet} vs Mass vs #it{p}_{T}^{#Lambda};#Delta#varphi_{jet};m_{p#pi} (GeV/c^{2});#it{p}_{T}^{#Lambda} (GeV/c)", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
      histos.add((folder + "/QA/h3dDeltaThetaVsMassVsLambdaPt").c_str(), "#Delta#theta_{jet} vs Mass vs #it{p}_{T}^{#Lambda};#Delta#theta_{jet};m_{p#pi} (GeV/c^{2});#it{p}_{T}^{#Lambda} (GeV/c)", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
      // --- 3D counters: Angle vs Mass vs Lead Jet pT ---
      histos.add((folder + "/QA/h3dDeltaPhiVsMassVsLeadJetPt").c_str(), "#Delta#varphi_{jet} vs Mass vs Lead Jet #it{p}_{T};#Delta#varphi_{jet};m_{p#pi} (GeV/c^{2});#it{p}_{T}^{LeadJet} (GeV/c)", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract});
      histos.add((folder + "/QA/h3dDeltaThetaVsMassVsLeadJetPt").c_str(), "#Delta#theta_{jet} vs Mass vs Lead Jet #it{p}_{T};#Delta#theta_{jet};m_{p#pi} (GeV/c^{2});#it{p}_{T}^{LeadJet} (GeV/c)", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract});

      // ===============================
      // TProfiles vs Mass: quick glancing before signal extraction
      // ===============================
      // TProfile of ring vs mass (integrated in all phi, and properly normalized by N_Lambda):
      histos.add((folder + "/pRingObservableMass").c_str(), "<#it{R}> vs Mass;m_{p#pi} (GeV/c^{2});<#it{R}>", kTProfile, {axisConfigurations.axisLambdaMassSigExtract});
      histos.add((folder + "/pRingObservableLeadPMass").c_str(), "<#it{R}> vs Mass (LeadP);m_{p#pi} (GeV/c^{2});<#it{R}>", kTProfile, {axisConfigurations.axisLambdaMassSigExtract});
      histos.add((folder + "/pRingObservable2ndJetMass").c_str(), "<#it{R}> vs Mass (SubJet);m_{p#pi} (GeV/c^{2});<#it{R}>", kTProfile, {axisConfigurations.axisLambdaMassSigExtract});
      // TProfile2D: <R> vs Mass (DeltaPhi)
      histos.add((folder + "/p2dRingObservableDeltaPhiVsMass").c_str(), "<#it{R}> vs #Delta#varphi_{jet} vs Mass;#Delta#varphi_{jet};m_{p#pi} (GeV/c^{2});<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract});
      // TProfile2D: <R> vs Mass (DeltaTheta)
      histos.add((folder + "/p2dRingObservableDeltaThetaVsMass").c_str(), "<#it{R}> vs #Delta#theta_{jet} vs Mass;#Delta#theta_{jet};m_{p#pi} (GeV/c^{2});<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract});
      // TProfile2D: <R> vs Eta Lambda vs Eta Jet (Understanding eta dependence seen in pRingEtaCuts)
      histos.add((folder + "/EtaDependence/hCounterEtaLambdaMinusEtaJet").c_str(), "N_{V0s} vs #eta_{#Lambda} - #eta_{Jet};#eta_{#Lambda} - #eta_{Jet}; N_{V0s}", kTH1D, {axisConfigurations.axisDeltaEtaCoarse});
      histos.add((folder + "/EtaDependence/pRingObservableEtaLambdaMinusEtaJet").c_str(), "<#it{R}> vs #eta_{#Lambda} - #eta_{Jet};#eta_{#Lambda} - #eta_{Jet};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaEtaCoarse});
      histos.add((folder + "/EtaDependence/p2dRingObservableEtaLambdaVsEtaJet").c_str(), "<#it{R}> vs #eta_{#Lambda} vs #eta_{Jet};#eta_{#Lambda};#eta_{Jet};<#it{R}>", kTProfile2D, {axisConfigurations.axisV0EtaCoarse, axisConfigurations.axisEtaCoarse});
      histos.add((folder + "/EtaDependence/p2dRingObservableEtaLambdaVsEtaJet_FineBins").c_str(), "<#it{R}> vs #eta_{#Lambda} vs #eta_{Jet} (fine bins);#eta_{#Lambda};#eta_{Jet};<#it{R}>", kTProfile2D, {axisConfigurations.axisV0Eta, axisConfigurations.axisEta});
      histos.add((folder + "/EtaDependence/p2dRingObservableEtaLambdaVsEtaLeadP").c_str(), "<#it{R}> vs #eta_{#Lambda} vs #eta_{LeadP};#eta_{#Lambda};#eta_{LeadP};<#it{R}>", kTProfile2D, {axisConfigurations.axisV0EtaCoarse, axisConfigurations.axisEtaCoarse});
      histos.add((folder + "/EtaDependence/p2dRingObservableEtaLambdaVsEta2ndJet").c_str(), "<#it{R}> vs #eta_{#Lambda} vs #eta_{SubJet};#eta_{#Lambda};#eta_{SubJet};<#it{R}>", kTProfile2D, {axisConfigurations.axisV0EtaCoarse, axisConfigurations.axisEtaCoarse});
      // Counters for these histograms, instead of only TProfile2Ds:
      histos.add((folder + "/EtaDependence/h2dCounterEtaLambdaVsEtaJet").c_str(), "Counts, #eta_{#Lambda} vs #eta_{Jet};#eta_{#Lambda};#eta_{Jet};Counts", kTH2D, {axisConfigurations.axisV0EtaCoarse, axisConfigurations.axisEtaCoarse});
      histos.add((folder + "/EtaDependence/h2dCounterEtaLambdaVsEtaJet_FineBins").c_str(), "Counts (fine bins), #eta_{#Lambda} vs #eta_{Jet};#eta_{#Lambda};#eta_{Jet};Counts", kTH2D, {axisConfigurations.axisV0Eta, axisConfigurations.axisEta});
      histos.add((folder + "/EtaDependence/h2dCounterEtaLambdaVsEtaLeadP").c_str(), "Counts, #eta_{#Lambda} vs #eta_{LeadP};#eta_{#Lambda};#eta_{LeadP};Counts", kTH2D, {axisConfigurations.axisV0EtaCoarse, axisConfigurations.axisEtaCoarse});
      histos.add((folder + "/EtaDependence/h2dCounterEtaLambdaVsEta2ndJet").c_str(), "Counts, #eta_{#Lambda} vs #eta_{SubJet};#eta_{#Lambda};#eta_{SubJet};Counts", kTH2D, {axisConfigurations.axisV0EtaCoarse, axisConfigurations.axisEtaCoarse});
      // --- TProfile3D: <R> vs DeltaPhi vs Mass vs LambdaPt ---
      histos.add((folder + "/p3dRingObservableDeltaPhiVsMassVsLambdaPt").c_str(), "<#it{R}> vs #Delta#varphi_{jet} vs Mass vs #it{p}_{T}^{#Lambda};#Delta#varphi_{jet};m_{p#pi} (GeV/c^{2});#it{p}_{T}^{#Lambda} (GeV/c)", kTProfile3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
      // --- TProfile3D: <R> vs DeltaTheta vs Mass vs LambdaPt ---
      histos.add((folder + "/p3dRingObservableDeltaThetaVsMassVsLambdaPt").c_str(), "<#it{R}> vs #Delta#theta_{jet} vs Mass vs #it{p}_{T}^{#Lambda};#Delta#theta_{jet};m_{p#pi} (GeV/c^{2});#it{p}_{T}^{#Lambda} (GeV/c)", kTProfile3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
      // --- TProfile3D: <R> vs DeltaPhi vs Mass vs LeadJetPt ---
      histos.add((folder + "/p3dRingObservableDeltaPhiVsMassVsLeadJetPt").c_str(), "<#it{R}> vs #Delta#varphi_{jet} vs Mass vs Lead Jet #it{p}_{T};#Delta#varphi_{jet};m_{p#pi} (GeV/c^{2});#it{p}_{T}^{LeadJet} (GeV/c)", kTProfile3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract});
      // --- TProfile3D: <R> vs DeltaTheta vs Mass vs LeadJetPt ---
      histos.add((folder + "/p3dRingObservableDeltaThetaVsMassVsLeadJetPt").c_str(), "<#it{R}> vs #Delta#theta_{jet} vs Mass vs Lead Jet #it{p}_{T};#Delta#theta_{jet};m_{p#pi} (GeV/c^{2});#it{p}_{T}^{LeadJet} (GeV/c)", kTProfile3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract});

      // ===============================
      // Mass histograms with centrality
      // ===============================
      // Counters
      histos.add((folder + "/QA/h3dDeltaPhiVsMassVsCent").c_str(), "#Delta#varphi_{jet} vs Mass vs Centrality;#Delta#varphi_{jet};m_{p#pi} (GeV/c^{2});Centrality (%)", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
      histos.add((folder + "/QA/h3dDeltaThetaVsMassVsCent").c_str(), "#Delta#theta_{jet} vs Mass vs Centrality;#Delta#theta_{jet};m_{p#pi} (GeV/c^{2});Centrality (%)", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
      // Useful TProfiles:
      // --- TProfile1D: Integrated <R> vs Centrality:
      histos.add((folder + "/pRingVsCentrality").c_str(), "<#it{R}> vs Centrality;Centrality (%);<#it{R}>", kTProfile, {axisConfigurations.axisCentrality});
      // --- TProfile2D: <R> vs Mass vs Centrality ---
      histos.add((folder + "/p2dRingObservableMassVsCent").c_str(), "<#it{R}> vs Mass vs Centrality;m_{p#pi} (GeV/c^{2});Centrality (%);<#it{R}>", kTProfile2D, {axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
      // --- TProfile3D: <R> vs DeltaPhi vs Mass vs Centrality ---
      histos.add((folder + "/p3dRingObservableDeltaPhiVsMassVsCent").c_str(), "<#it{R}> vs #Delta#varphi_{jet} vs Mass vs Centrality;#Delta#varphi_{jet};m_{p#pi} (GeV/c^{2});Centrality (%)", kTProfile3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
      // --- TProfile3D: <R> vs DeltaTheta vs Mass vs Centrality ---
      histos.add((folder + "/p3dRingObservableDeltaThetaVsMassVsCent").c_str(), "<#it{R}> vs #Delta#theta_{jet} vs Mass vs Centrality;#Delta#theta_{jet};m_{p#pi} (GeV/c^{2});Centrality (%)", kTProfile3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});

      // ===============================
      // QA histograms - Useful numbers
      // ===============================
      // (TODO: implement these!)
      // (TODO: implement momentum imbalance checks for jets!)
      // Added to a separate folder for further control (changed the usage of the "folder" string):
      // histos.add(("QA_Numbers/" + folder + "/hValidLeadJets").c_str(), "hValidLeadJets", kTH1D, {{1,0,1}});
      // TODO: Add "frequency of jets per pT" histograms either here or in the TableProducer

      // Estimating error bars with the Delta Method for <R> = A/B:
      // 1D Delta Method for Integrated observable:
      histos.add((folder + "/DeltaMethod/hIntegrated").c_str(), "Delta Method Accumulators Integrated;Component;Counts", kTH1D, {axisConfigurations.axisDeltaComponents});

      // 2D Delta Method for Differentials
      histos.add((folder + "/DeltaMethod/h2dDeltaThetaVsDeltaComp").c_str(), "Delta Method vs #Delta#theta_{jet};#Delta#theta_{jet};Component", kTH2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisDeltaComponents});
      histos.add((folder + "/DeltaMethod/h2dLambdaPtVsDeltaComp").c_str(), "Delta Method vs #Lambda #it{p}_{T};#it{p}_{T}^{#Lambda} (GeV/c);Component", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisDeltaComponents});
      histos.add((folder + "/DeltaMethod/h2dMassVsDeltaComp").c_str(), "Delta Method vs Mass;m_{p#pi} (GeV/c^{2});Component", kTH2D, {axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisDeltaComponents});
    };
    // Execute local lambda to register histogram families:
    if (familySwitches.doFamilyRing)
      addRingObservableFamily("Ring");
    if (familySwitches.doFamilyRingKinematicCuts)
      addRingObservableFamily("RingKinematicCuts");
    if (familySwitches.doFamilyJetKinematicCuts)
      addRingObservableFamily("JetKinematicCuts");
    if (familySwitches.doFamilyJetAndLambdaKinematicCuts)
      addRingObservableFamily("JetAndLambdaKinematicCuts");

    histos.add("IntegratedCuts/pRingCuts", "pRingCuts; ;<#it{R}>", kTProfile, {{4, 0, 4}});
    histos.get<TProfile>(HIST("IntegratedCuts/pRingCuts"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TProfile>(HIST("IntegratedCuts/pRingCuts"))->GetXaxis()->SetBinLabel(2, "p_{T}^{#Lambda}@[0.5,1.5],|y_{#Lambda}|<0.5"); // (v0pt > 0.5 && v0pt < 1.5) && std::abs(lambdaRapidity) < 0.5;
    histos.get<TProfile>(HIST("IntegratedCuts/pRingCuts"))->GetXaxis()->SetBinLabel(3, "|Jet_{#eta}|<0.5");
    histos.get<TProfile>(HIST("IntegratedCuts/pRingCuts"))->GetXaxis()->SetBinLabel(4, "#Lambda + Jet cuts");

    // Same for subleading jet and leading particle:
    histos.add("IntegratedCuts/pRingCutsSubLeadingJet", "pRingCutsSubLeadingJet; ;<#it{R}>", kTProfile, {{4, 0, 4}});
    histos.get<TProfile>(HIST("IntegratedCuts/pRingCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TProfile>(HIST("IntegratedCuts/pRingCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(2, "p_{T,#Lambda}@[0.5,1.5],|y_{#Lambda}|<0.5");
    histos.get<TProfile>(HIST("IntegratedCuts/pRingCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(3, "|SubJet_{#eta}|<0.5");
    histos.get<TProfile>(HIST("IntegratedCuts/pRingCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(4, "#Lambda + SubJet cuts");

    histos.add("IntegratedCuts/pRingCutsLeadingP", "pRingCutsLeadingP; ;<#it{R}>", kTProfile, {{4, 0, 4}});
    histos.get<TProfile>(HIST("IntegratedCuts/pRingCutsLeadingP"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TProfile>(HIST("IntegratedCuts/pRingCutsLeadingP"))->GetXaxis()->SetBinLabel(2, "p_{T}^{#Lambda}@[0.5,1.5],|y_{#Lambda}|<0.5");
    histos.get<TProfile>(HIST("IntegratedCuts/pRingCutsLeadingP"))->GetXaxis()->SetBinLabel(3, "|LeadP_{#eta}|<0.5");
    histos.get<TProfile>(HIST("IntegratedCuts/pRingCutsLeadingP"))->GetXaxis()->SetBinLabel(4, "#Lambda + LeadP cuts");

    // Counters for each case to understand statistics loss:
    histos.add("IntegratedCuts/hCountCuts", "hCountCuts; ;N V0s", kTH1D, {{4, 0, 4}});
    histos.get<TH1>(HIST("IntegratedCuts/hCountCuts"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TH1>(HIST("IntegratedCuts/hCountCuts"))->GetXaxis()->SetBinLabel(2, "p_{T}^{#Lambda}@[0.5,1.5],|y_{#Lambda}|<0.5"); // (v0pt > 0.5 && v0pt < 1.5) && std::abs(lambdaRapidity) < 0.5;
    histos.get<TH1>(HIST("IntegratedCuts/hCountCuts"))->GetXaxis()->SetBinLabel(3, "|Jet_{#eta}|<0.5");
    histos.get<TH1>(HIST("IntegratedCuts/hCountCuts"))->GetXaxis()->SetBinLabel(4, "#Lambda + Jet cuts");

    // Same for subleading jet and leading particle:
    histos.add("IntegratedCuts/hCountCutsSubLeadingJet", "hCountCutsSubLeadingJet; ;N V0s", kTH1D, {{4, 0, 4}});
    histos.get<TH1>(HIST("IntegratedCuts/hCountCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TH1>(HIST("IntegratedCuts/hCountCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(2, "p_{T,#Lambda}@[0.5,1.5],|y_{#Lambda}|<0.5");
    histos.get<TH1>(HIST("IntegratedCuts/hCountCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(3, "|SubJet_{#eta}|<0.5");
    histos.get<TH1>(HIST("IntegratedCuts/hCountCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(4, "#Lambda + SubJet cuts");

    histos.add("IntegratedCuts/hCountCutsLeadingP", "hCountCutsLeadingP; ;N V0s", kTH1D, {{4, 0, 4}});
    histos.get<TH1>(HIST("IntegratedCuts/hCountCutsLeadingP"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TH1>(HIST("IntegratedCuts/hCountCutsLeadingP"))->GetXaxis()->SetBinLabel(2, "p_{T}^{#Lambda}@[0.5,1.5],|y_{#Lambda}|<0.5");
    histos.get<TH1>(HIST("IntegratedCuts/hCountCutsLeadingP"))->GetXaxis()->SetBinLabel(3, "|LeadP_{#eta}|<0.5");
    histos.get<TH1>(HIST("IntegratedCuts/hCountCutsLeadingP"))->GetXaxis()->SetBinLabel(4, "#Lambda + LeadP cuts");

    // Fake-polarization diagnostics:
    if (qaSwitches.doFakePolDiagnosticsQA) {
      // Integrated observable dependent on jet proxy #eta to unfold possible asymmetries in detector:
      histos.add("EtaStudy/pRingEtaCuts", "pRingEtaCuts; ;<#it{R}>", kTProfile, {{15, 0, 15}});
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(10, "#eta_{Jet} > R");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(11, "#eta_{Jet} < -R");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(12, "#eta_{Jet} > R, #eta_{#Lambda} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(13, "#eta_{Jet} > R, #eta_{#Lambda} < 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(14, "#eta_{Jet} < -R, #eta_{#Lambda} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCuts"))->GetXaxis()->SetBinLabel(15, "#eta_{Jet} < -R, #eta_{#Lambda} < 0");

      histos.add("EtaStudy/pRingEtaCutsSubLeadingJet", "pRingEtaCutsSubLeadingJet; ;<#it{R}>", kTProfile, {{15, 0, 15}});
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(2, "#eta_{SubJet} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(3, "#eta_{SubJet} < 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(6, "#eta_{SubJet} #geq 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(7, "#eta_{SubJet} #geq 0, #eta_{#Lambda} < 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(8, "#eta_{SubJet} < 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(9, "#eta_{SubJet} < 0, #eta_{#Lambda} < 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(10, "#eta_{SubJet} > R");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(11, "#eta_{SubJet} < -R");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(12, "#eta_{SubJet} > R, #eta_{#Lambda} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(13, "#eta_{SubJet} > R, #eta_{#Lambda} < 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(14, "#eta_{SubJet} < -R, #eta_{#Lambda} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(15, "#eta_{SubJet} < -R, #eta_{#Lambda} < 0");

      histos.add("EtaStudy/pRingEtaCutsLeadingP", "pRingEtaCutsLeadingP; ;<#it{R}>", kTProfile, {{9, 0, 9}});
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(2, "#eta_{LeadP} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(3, "#eta_{LeadP} < 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(6, "#eta_{LeadP} #geq 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(7, "#eta_{LeadP} #geq 0, #eta_{#Lambda} < 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(8, "#eta_{LeadP} < 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile>(HIST("EtaStudy/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(9, "#eta_{LeadP} < 0, #eta_{#Lambda} < 0");

      // Studying the signal Vs background integral (a very naive estimative of the invariant mass peak)
      histos.add("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground", "pRingEtaCutsLeadingP_MassSignalVsBackground; ; ;<#it{R}>", kTProfile2D, {{9, 0, 9}, {2, 0, 2}});
      histos.get<TProfile2D>(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
      histos.get<TProfile2D>(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"))->GetXaxis()->SetBinLabel(2, "#eta_{LeadP} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"))->GetXaxis()->SetBinLabel(3, "#eta_{LeadP} < 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"))->GetXaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"))->GetXaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"))->GetXaxis()->SetBinLabel(6, "#eta_{LeadP} #geq 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"))->GetXaxis()->SetBinLabel(7, "#eta_{LeadP} #geq 0, #eta_{#Lambda} < 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"))->GetXaxis()->SetBinLabel(8, "#eta_{LeadP} < 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"))->GetXaxis()->SetBinLabel(9, "#eta_{LeadP} < 0, #eta_{#Lambda} < 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"))->GetYaxis()->SetBinLabel(1, "#Lambda out of mass peak");
      histos.get<TProfile2D>(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"))->GetYaxis()->SetBinLabel(2, "#Lambda in mass peak");

      // Fake polarization signal QA
      // --> The "negative helicity problem", where topologies with a proton decaying opposite to the Lambda momentum are enhanced by
      // efficiency of reconstruction. The geometries where the proton moves in the same direction as the boost will have a very small
      // momentum pion, which is not as easily detected as the opposite case! This may insert a fake signal of polarization in the measurement!
      histos.add("EtaStudy/hFakePolCounts", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda};", kTH2D, {axisConfigurations.axisCosTheta, {9, 0, 9}});
      histos.get<TH2>(HIST("EtaStudy/hFakePolCounts"))->GetZaxis()->SetTitle("N_{V0s}");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCounts"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCounts"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCounts"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCounts"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCounts"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCounts"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCounts"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCounts"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCounts"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");

      // The same, but for actual signal instead of counts:
      histos.add("EtaStudy/pFakePolSignalVsCosTheta", "FakePolSignal; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda};", kTProfile2D, {axisConfigurations.axisCosTheta, {9, 0, 9}});
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalVsCosTheta"))->GetZaxis()->SetTitle("<#it{R}>");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");

      // Seeing the dependence between phi* = atan2(p_p_star \cdot (p_Lambda_hat \times (z_hat \cross p_Lambda)), p_p_star \cdot (z_hat \cross p_Lambda))
      // e_z = p_Lambda_hat; // e_x = normalize(z_hat cross p_Lambda); // e_y = e_z cross e_x;
      // phi_star = atan2(p_p_star dot e_y, p_p_star dot e_x);
      histos.add("HelicityEfficiencyQA/hFakePolCounts_CosThetaVsPhiStar", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; #phi^{*}", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaPhi});
      histos.add("HelicityEfficiencyQA/pFakePolSignal_CosThetaVsPhiStar", "FakePolSignal; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; #phi^{*}", kTProfile2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaPhi});
      // Specific counter for when we have leading jets (relates directly to pFakePolSignal_CosThetaVsPhiStar):
      histos.add("HelicityEfficiencyQA/hFakePolCountsJet_CosThetaVsPhiStar", "FakePolCounts - HasValidLeadJet OK; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; #phi^{*}", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaPhi});

      // Similar split, but for AEE instead of HEE:
      histos.add("EtaStudy/hCountsVsPhiStar", "FakePolCounts, AEE dependence; #phi^{*};", kTH2D, {axisConfigurations.axisDeltaPhi, {9, 0, 9}});
      histos.get<TH2>(HIST("EtaStudy/hCountsVsPhiStar"))->GetZaxis()->SetTitle("Counts");
      histos.get<TH2>(HIST("EtaStudy/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
      histos.get<TH2>(HIST("EtaStudy/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
      histos.get<TH2>(HIST("EtaStudy/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
      histos.get<TH2>(HIST("EtaStudy/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
      histos.get<TH2>(HIST("EtaStudy/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");
      // For the ring observable as well:
      // Explicitly, checking the Phi* dependence on a series of #eta cuts.
      // The fake, AEE-induced, signal should be zero when integrating on full solid angle, and then have some shape for each eta slice.
      histos.add("EtaStudy/pFakePolSignalvsPhiStar", "FakePolSignal, AEE dependence; #phi^{*};", kTProfile2D, {axisConfigurations.axisDeltaPhi, {9, 0, 9}});
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiStar"))->GetZaxis()->SetTitle("<#it{R}>");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");

      // For the phi_Lambda - phi_D^* dependency as well:
      histos.add("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar", "FakePolSignal, AEE dependence; #phi_{#Lambda} - #phi_{p}^{*};", kTProfile2D, {axisConfigurations.axisDeltaPhi, {9, 0, 9}});
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetZaxis()->SetTitle("<#it{R}>");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
      histos.get<TProfile2D>(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");

      // More about possible AEE dependencies (should see an invariance with JetEta and <R>/Jz):
      // TODO: think about these error bars: do they still make sense via regular TProfile's SEM error? These are just a quick check, so wouldn't bother much about it.
      histos.add("HelicityEfficiencyQA/pRingVsJetZcomponent", "<#it{R}> vs #hat{t}_{z}; #hat{t}_{z}; <#it{R}>", kTProfile, {{40, -1, 1}}); // Numerically stable and can also show the sign flip (essentially the <#it{R}> vs Eta Jet plot in another scale)
      histos.add("HelicityEfficiencyQA/pRingOverJetZcomponent_VsJetEta", "<#it{R}>/#hat{t}_{z}; #eta_{Jet}; <#it{R}>/#hat{t}_{z}", kTProfile, {axisConfigurations.axisEtaCoarse});
      histos.add("HelicityEfficiencyQA/pRingOverJetZcomponent_VsCosThetaHEE", "<#it{R}>/#hat{t}_{z} Vs cos(#theta) HEE; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; <#it{R}>/#hat{t}_{z}", kTProfile, {axisConfigurations.axisCosTheta});
      histos.add("HelicityEfficiencyQA/pRingOverJetZcomponent_VsPhiStar", "<#it{R}>/#hat{t}_{z} Vs #phi^{*}; #phi^{*} = atan2(#vec{p}^{*}_{p} #cdot [#hat{p}_{#Lambda} #times (#hat{z} #times #hat{p}_{#Lambda})] , #vec{p}^{*}_{p} #cdot (#hat{z} #times #hat{p}_{#Lambda})); <#it{R}>/#hat{t}_{z}", kTProfile, {axisConfigurations.axisDeltaPhi});
      histos.add("HelicityEfficiencyQA/pRingOverJetZcomponent_VsJetEtaVsCosThetaHEE", "<#it{R}>/#hat{t}_{z}; #eta_{Jet}; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; <#it{R}>/#hat{t}_{z}", kTProfile2D, {axisConfigurations.axisEtaCoarse, axisConfigurations.axisCosTheta});
      histos.add("HelicityEfficiencyQA/pRingOverJetZcomponent_VsJetEtaVsPhiStar", "<#it{R}>/#hat{t}_{z}; #eta_{Jet}; #phi^{*}; <#it{R}>/#hat{t}_{z}", kTProfile2D, {axisConfigurations.axisEtaCoarse, axisConfigurations.axisDeltaPhi});

      // Seeing if phi* is indeed influenced by the DCA between daughters:
      histos.add("HelicityEfficiencyQA/hFakePolCountsJet_PhiStarVsDCAdau", "FakePolCounts - HasValidLeadJet OK; #phi^{*}; DCA_{V0 Daughters}", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDCAdau});
      histos.add("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAdau", "FakePolSignal; #phi^{*}; DCA_{V0 Daughters}", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDCAdau});

      // Adding a way to check if the jet eta is positive or negative as well
      // (AEE signal could be closer to zero otherwise: phi^* dependency may not make it fall to zero as we are no longer integrating <R> in full solid angle, yet analyzing this other dependency is also important)
      histos.add("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAdauVsEtaJet", "FakePolSignal; #phi^{*}; DCA_{V0 Daughters}; #eta_{Jet} sign", kTProfile3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDCAdau, {2, -0.9, 0.9}});
      histos.add("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAdauVsEtaLambda", "FakePolSignal; #phi^{*}; DCA_{V0 Daughters}; #eta_{#Lambda} sign", kTProfile3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDCAdau, {2, -0.9, 0.9}});

      // Similar checks for dcaPosToPV and dcaNegToPV, which influence AEE the strongest:
      histos.add("HelicityEfficiencyQA/hFakePolCountsJet_PhiStarVsDCAProLike", "FakePolCounts - HasValidLeadJet OK; #phi^{*}; DCA_{PosPV}", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDCAdauPV});
      histos.add("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAProLike", "FakePolSignal; #phi^{*}; DCA_{PosPV}", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDCAdauPV});
      histos.add("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAProLikeVsEtaJet", "FakePolSignal; #phi^{*}; DCA_{PosPV}; #eta_{Jet} sign", kTProfile3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDCAdauPV, {2, -0.9, 0.9}});
      histos.add("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAProLikeVsEtaLambda", "FakePolSignal; #phi^{*}; DCA_{PosPV}; #eta_{#Lambda} sign", kTProfile3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDCAdauPV, {2, -0.9, 0.9}});

      histos.add("HelicityEfficiencyQA/hFakePolCountsJet_PhiStarVsDCAPiLike", "FakePolCounts - HasValidLeadJet OK; #phi^{*}; DCA_{NegPV}", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDCAdauPV});
      histos.add("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAPiLike", "FakePolSignal; #phi^{*}; DCA_{NegPV}", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDCAdauPV});
      histos.add("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAPiLikeVsEtaJet", "FakePolSignal; #phi^{*}; DCA_{NegPV}; #eta_{Jet} sign", kTProfile3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDCAdauPV, {2, -0.9, 0.9}});
      histos.add("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAPiLikeVsEtaLambda", "FakePolSignal; #phi^{*}; DCA_{NegPV}; #eta_{#Lambda} sign", kTProfile3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDCAdauPV, {2, -0.9, 0.9}});

      // Doing the same HEE study for leading particles:
      // (eta_{Jet} may be a bad estimator!)
      histos.add("EtaStudy/hFakePolCountsLeadP", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda};", kTH2D, {axisConfigurations.axisCosTheta, {9, 0, 9}});
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLeadP"))->GetZaxis()->SetTitle("N_{V0s}");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(2, "#eta_{LeadP} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(3, "#eta_{LeadP} < 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(6, "#eta_{LeadP} #geq 0, #eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(7, "#eta_{LeadP} #geq 0, #eta_{#Lambda} < 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(8, "#eta_{LeadP} < 0, #eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(9, "#eta_{LeadP} < 0, #eta_{#Lambda} < 0");

      // Avoid fake signal by jets boosting the Lambda in its own direction, then modifying efficiency of reconstruction in a similar way:
      histos.add("EtaStudy/hFakePolCountsLambdaPtCut", "FakePol,p_{T}^{#Lambda}#in[0.5,1.5]; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda};", kTH2D, {axisConfigurations.axisCosTheta, {9, 0, 9}});
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtCut"))->GetZaxis()->SetTitle("N_{V0s}");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");

      // Even stricter cut (also demands rapidity cut stricter than jets, so may see different boosting):
      histos.add("EtaStudy/hFakePolCountsLambdaPtYCuts", "FakePol,p_{T}^{#Lambda}#in[0.5,1.5],|y_{#Lambda}|<0.5; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda};", kTH2D, {axisConfigurations.axisCosTheta, {9, 0, 9}});
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"))->GetZaxis()->SetTitle("N_{V0s}");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
      histos.get<TH2>(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");

      // Another useful quantity -- How much is the fake signal related to the jet's momentum (how much the fake signal is correlated with the Jet-Lambda angular separation):
      histos.add("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaJet", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; #Delta#theta_{Jet}", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaTheta});
      histos.add("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaJetPosEta", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; #Delta#theta_{Jet}", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaTheta});
      histos.add("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaJetNegEta", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; #Delta#theta_{Jet}", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaTheta});
      histos.add("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaLeadP", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; #Delta#theta_{LeadP}", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaTheta});
      histos.add("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaLeadPPosEta", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; #Delta#theta_{LeadP}", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaTheta});
      histos.add("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaLeadPNegEta", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; #Delta#theta_{LeadP}", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaTheta});

      // Understanding the dip at the cos = -1 end:
      histos.add("HelicityEfficiencyQA/hFakePolCountsCosThetaVsPtForJets", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; p_{T}^{#Lambda}", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisPtCoarseQA});
      histos.add("HelicityEfficiencyQA/hFakePolCountsCosThetaVsPtForLeadP", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; p_{T}^{#Lambda}", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisPtCoarseQA});

      // Studying the magnetic field dependence of particle reconstruction efficiency (not magnitude, just sign of field):
      // (also for the "negative helicity" problem)
      if (analyseMagField) {
        histos.add("HelicityEfficiencyQA/hLambdaMassDecayGeomRight", "hLambdaMassDecayGeomRight; m_{Inv}; Counts", kTH1D, {axisConfigurations.axisLambdaMass});
        histos.add("HelicityEfficiencyQA/hLambdaMassDecayGeomLeft", "hLambdaMassDecayGeomLeft; m_{Inv}; Counts", kTH1D, {axisConfigurations.axisLambdaMass});
        histos.add("HelicityEfficiencyQA/hAntiLambdaMassDecayGeomRight", "hAntiLambdaMassDecayGeomRight; m_{Inv}; Counts", kTH1D, {axisConfigurations.axisLambdaMass});
        histos.add("HelicityEfficiencyQA/hAntiLambdaMassDecayGeomLeft", "hAntiLambdaMassDecayGeomLeft; m_{Inv}; Counts", kTH1D, {axisConfigurations.axisLambdaMass});
      }

      // Also including an observable that probes spectrum broadening due to the "Azimuthal Efficiency Effect":
      histos.add("HelicityEfficiencyQA/hLambdaMassVsPhiLambdaMinusPhiProtonStar", "m_{#Lambda}, AEE probe; m_{Inv}; #phi_{#Lambda} - #phi_{p}^{*} ; Counts", kTH2D, {axisConfigurations.axisLambdaMass, axisConfigurations.axisDeltaPhi});
      histos.add("HelicityEfficiencyQA/hAntiLambdaMassVsPhiLambdaMinusPhiProtonStar", "m_{#bar{#Lambda}}, AEE probe; m_{Inv}; #phi_{#bar{#Lambda}} - #phi_{p}^{*} ; Counts", kTH2D, {axisConfigurations.axisLambdaMass, axisConfigurations.axisDeltaPhi});
      // Watching the effect on the ring observable as well:
      histos.add("HelicityEfficiencyQA/p2dRing_LambdaMassVsPhiLambdaMinusPhiProtonStar", "<#it{R}>, AEE probe; m_{Inv}; #phi_{#Lambda} - #phi_{p}^{*} ; <#it{R}>", kTProfile2D, {axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisDeltaPhi});
      histos.add("HelicityEfficiencyQA/p2dRing_AntiLambdaMassVsPhiLambdaMinusPhiProtonStar", "<#it{R}>, AEE probe; m_{Inv}; #phi_{#bar{#Lambda}} - #phi_{p}^{*} ; <#it{R}>", kTProfile2D, {axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisDeltaPhi});
    } // end doFakePolDiagnosticsQA bookings

    // Integrated observable for events with NLambda+NAntiLambda V0s per event
    // (an interesting measurement of correlation between <R> and Lambda-like V0s multiplicity. A proxy of covariance)
    // (calculated for leading jets only)
    histos.add("IntegratedCuts/pRingVsNV0s", "pRingVsNV0s; N_{#Lambda}+N_{#bar{#Lambda}};<#it{R}>", kTProfile, {{20, 0, 20}}); // See hNV0sVsCentrality below for the correlation between number of V0s and centrality
    histos.add("hNV0sVsCentrality", "hNV0sVsCentrality; N_{#Lambda}+N_{#bar{#Lambda}};Centrality (%)", kTH2D, {{20, 0, 20}, axisConfigurations.axisCentrality});

    // Proxy Eta QA:
    histos.add("JetKinematicsQA/hLeadJetEta", "hLeadJetEta;#eta;Counts", kTH1D, {axisConfigurations.axisEta});
    histos.add("JetKinematicsQA/hSubLeadJetEta", "hSubLeadJetEta;#eta;Counts", kTH1D, {axisConfigurations.axisEta});
    histos.add("JetKinematicsQA/hLeadPEta", "hLeadPEta;#eta;Counts", kTH1D, {axisConfigurations.axisEta});

    // Proxy Phi QA:
    histos.add("JetKinematicsQA/hLeadJetPhi", "hLeadJetPhi;#varphi;Counts", kTH1D, {axisConfigurations.axisPhi});
    histos.add("JetKinematicsQA/hSubLeadJetPhi", "hSubLeadJetPhi;#varphi;Counts", kTH1D, {axisConfigurations.axisPhi});
    histos.add("JetKinematicsQA/hLeadPPhi", "hLeadPPhi;#varphi;Counts", kTH1D, {axisConfigurations.axisPhi});

    // Counting the number of jets/proxies themselves (these count at most once per event) -- Similar to the FOLDER/QA/hPtJet counters:
    histos.add("JetKinematicsQA/hJetCounterPtJet", "hJetCounterPtJet; p_{T}^{Jet} (GeV/c)", kTH1D, {axisConfigurations.axisJetPt});
    histos.add("JetKinematicsQA/hJetCounterPtLeadP", "hJetCounterPtLeadP; p_{T}^{LeadP} (GeV/c)", kTH1D, {axisConfigurations.axisJetPt});
    histos.add("JetKinematicsQA/hJetCounterPt2ndJet", "hJetCounterPt2ndJet; p_{T}^{SubJet} (GeV/c)", kTH1D, {axisConfigurations.axisJetPt});

    // Proxy Eta vs PVz:
    histos.add("JetKinematicsQA/h2dLeadJetEtaVsPVz", "Lead Jet #eta Vs PVz;#eta;Primary Vertex Z [cm];Counts", kTH2D, {axisConfigurations.axisEta, axisConfigurations.axisPVz});
    histos.add("JetKinematicsQA/h2dSubLeadJetEtaVsPVz", "SubLead Jet #eta Vs PVz;#eta;Primary Vertex Z [cm];Counts", kTH2D, {axisConfigurations.axisEta, axisConfigurations.axisPVz});
    histos.add("JetKinematicsQA/h2dLeadPEtaVsPVz", "Lead Ptc #eta Vs PVz;#eta;Primary Vertex Z [cm];Counts", kTH2D, {axisConfigurations.axisEta, axisConfigurations.axisPVz});

    // For building and event-mixing-like procedure similar to forceDatalikeJet:
    if (doJetProxy5dQA) {
      histos.add("JetKinematicsQA/h5dLeadJetEtaPhiPtPVzCent", "h5dLeadJetEtaPhiPtPVzCent;#eta;#phi;p_{t};Primary Vertex Z [cm];Centrality (%);Counts", kTHnSparseF,
                 {axisConfigurations.axisEtaCoarse, axisConfigurations.axisPhi, axisConfigurations.axisJetPt, axisConfigurations.axisPVz, axisConfigurations.axisCentrality});
      histos.add("JetKinematicsQA/h5dSubLeadJetEtaPhiPtPVzCent", "h5dSubLeadJetEtaPhiPtPVzCent;#eta;#phi;p_{t};Primary Vertex Z [cm];Centrality (%);Counts", kTHnSparseF,
                 {axisConfigurations.axisEtaCoarse, axisConfigurations.axisPhi, axisConfigurations.axisJetPt, axisConfigurations.axisPVz, axisConfigurations.axisCentrality});
      histos.add("JetKinematicsQA/h5dLeadPEtaPhiPtPVzCent", "h5dLeadPEtaPhiPtPVzCent;#eta;#phi;p_{t};Primary Vertex Z [cm];Centrality (%);Counts", kTHnSparseF,
                 {axisConfigurations.axisEtaCoarse, axisConfigurations.axisPhi, axisConfigurations.axisJetPt, axisConfigurations.axisPVz, axisConfigurations.axisCentrality});
    }

    // doMixedEventProxies QA: gauge the size of the "too few collisions per bin" problem (see resonanceMergeDF.cxx)
    // Booked per proxy, since the three mixings are independent and can succeed/fail at different rates.
    if (fakePolSwitches.doMixedEventProxies) {
      histos.add("EventMixingQA/hMixedEventLeadPOutcome", "hMixedEventLeadPOutcome;Outcome (0=skipped, 1=found);Counts", kTH1D, {{2, -0.5f, 1.5f}});
      histos.add("EventMixingQA/hMixedEventLeadJetOutcome", "hMixedEventLeadJetOutcome;Outcome (0=skipped, 1=found);Counts", kTH1D, {{2, -0.5f, 1.5f}});
      histos.add("EventMixingQA/hMixedEventSubJetOutcome", "hMixedEventSubJetOutcome;Outcome (0=skipped, 1=found);Counts", kTH1D, {{2, -0.5f, 1.5f}});

      histos.add("EventMixingQA/hMixedEventLeadPWindowNeighbours", "hMixedEventLeadPWindowNeighbours;Neighbours found in bin window;Counts", kTH1D, {{22, -0.5f, 21.5f}});
      histos.add("EventMixingQA/hMixedEventLeadJetWindowNeighbours", "hMixedEventLeadJetWindowNeighbours;Neighbours found in bin window;Counts", kTH1D, {{22, -0.5f, 21.5f}});
      histos.add("EventMixingQA/hMixedEventSubJetWindowNeighbours", "hMixedEventSubJetWindowNeighbours;Neighbours found in bin window;Counts", kTH1D, {{22, -0.5f, 21.5f}});

      histos.add("EventMixingQA/h2dMixedLeadPEtaVsLeadPEta", "MixedLeadP #eta vs LeadP #eta;MixedLeadP #eta; LeadP #eta;Counts", kTH2D, {axisConfigurations.axisEtaCoarse, axisConfigurations.axisEtaCoarse});
      histos.add("EventMixingQA/h2dMixedLeadPPhiVsLeadPPhi", "MixedLeadP #phi vs LeadP #phi;MixedLeadP #phi; LeadP #phi;Counts", kTH2D, {axisConfigurations.axisPhi, axisConfigurations.axisPhi});
      histos.add("EventMixingQA/h2dMixedLeadJetEtaVsLeadJetEta", "MixedLeadJet #eta vs LeadJet #eta;MixedLeadJet #eta; LeadJet #eta;Counts", kTH2D, {axisConfigurations.axisEtaCoarse, axisConfigurations.axisEtaCoarse});
      histos.add("EventMixingQA/h2dMixedLeadJetPhiVsLeadJetPhi", "MixedLeadJet #phi vs LeadJet #phi;MixedLeadJet #phi; LeadJet #phi;Counts", kTH2D, {axisConfigurations.axisPhi, axisConfigurations.axisPhi});
      histos.add("EventMixingQA/h2dMixedSubJetEtaVsSubJetEta", "MixedSubJet #eta vs SubJet #eta;MixedSubJet #eta; SubJet #eta;Counts", kTH2D, {axisConfigurations.axisEtaCoarse, axisConfigurations.axisEtaCoarse});
      histos.add("EventMixingQA/h2dMixedSubJetPhiVsSubJetPhi", "MixedSubJet #phi vs SubJet #phi;MixedSubJet #phi; SubJet #phi;Counts", kTH2D, {axisConfigurations.axisPhi, axisConfigurations.axisPhi});

      // Collision-index proximity QA:
      // (To understand possible continous readout effects on the choice of event mixing sources -- do notice index and time are both monotonic, but there is a conversion factor)
      // The shape of the whole candidate pool:
      histos.add("EventMixingQA/hMixedEventLeadPDeltaIndexEligible", "hMixedEventLeadPDeltaIndexEligible;#Delta collision index (pair);Counts", kTH1D, {axisConfigurations.axisDeltaCollisionIndex});
      histos.add("EventMixingQA/hMixedEventLeadJetDeltaIndexEligible", "hMixedEventLeadJetDeltaIndexEligible;#Delta collision index (pair);Counts", kTH1D, {axisConfigurations.axisDeltaCollisionIndex});
      histos.add("EventMixingQA/hMixedEventSubJetDeltaIndexEligible", "hMixedEventSubJetDeltaIndexEligible;#Delta collision index (pair);Counts", kTH1D, {axisConfigurations.axisDeltaCollisionIndex});

      // What the reservoir actually picked:
      // (should be equally as narrow as the whole candidate pool, if no preferential mixing exists)
      histos.add("EventMixingQA/hMixedEventLeadPDeltaIndexSelected", "hMixedEventLeadPDeltaIndexSelected;#Delta collision index (selected partner);Counts", kTH1D, {axisConfigurations.axisDeltaCollisionIndex});
      histos.add("EventMixingQA/hMixedEventLeadJetDeltaIndexSelected", "hMixedEventLeadJetDeltaIndexSelected;#Delta collision index (selected partner);Counts", kTH1D, {axisConfigurations.axisDeltaCollisionIndex});
      histos.add("EventMixingQA/hMixedEventSubJetDeltaIndexSelected", "hMixedEventSubJetDeltaIndexSelected;#Delta collision index (selected partner);Counts", kTH1D, {axisConfigurations.axisDeltaCollisionIndex});

      // Event mixing source QA -- How repeated is each collision in the mix:
      // For each source collision for mixing, counts how many times it was used. Flat distribution is ideal.
      histos.add("EventMixingQA/hMixedEventLeadPSourceUsageCount", "hMixedEventLeadPSourceUsageCount;Times collision was used;Counts", kTH1D, {{50, -0.5f, 49.5f}});
      histos.add("EventMixingQA/hMixedEventLeadJetSourceUsageCount", "hMixedEventLeadJetSourceUsageCount;Times collision was used;Counts", kTH1D, {{50, -0.5f, 49.5f}});
      histos.add("EventMixingQA/hMixedEventSubJetSourceUsageCount", "hMixedEventSubJetSourceUsageCount;Times collision was used;Counts", kTH1D, {{50, -0.5f, 49.5f}});
      // Useful TProfiles -- the mean number of times a given source was used, as a function of eta or phi
      // (more convenient than a single TH1, as it gives an average, not the raw counter)
      histos.add("EventMixingQA/pMixedEventLeadPSourceUsageVsEta", "pMixedEventLeadPSourceUsageVsEta;Source #eta;<Times used>", kTProfile, {axisConfigurations.axisEtaCoarse});
      histos.add("EventMixingQA/pMixedEventLeadPSourceUsageVsPhi", "pMixedEventLeadPSourceUsageVsPhi;Source #varphi;<Times used>", kTProfile, {axisConfigurations.axisPhi});
      histos.add("EventMixingQA/pMixedEventLeadJetSourceUsageVsEta", "pMixedEventLeadJetSourceUsageVsEta;Source #eta;<Times used>", kTProfile, {axisConfigurations.axisEtaCoarse});
      histos.add("EventMixingQA/pMixedEventLeadJetSourceUsageVsPhi", "pMixedEventLeadJetSourceUsageVsPhi;Source #varphi;<Times used>", kTProfile, {axisConfigurations.axisPhi});
      histos.add("EventMixingQA/pMixedEventSubJetSourceUsageVsEta", "pMixedEventSubJetSourceUsageVsEta;Source #eta;<Times used>", kTProfile, {axisConfigurations.axisEtaCoarse});
      histos.add("EventMixingQA/pMixedEventSubJetSourceUsageVsPhi", "pMixedEventSubJetSourceUsageVsPhi;Source #varphi;<Times used>", kTProfile, {axisConfigurations.axisPhi});
    }

    // Fetch the X-axes from one of the families (since they all share the same ConfigurableAxis binning)
    mAxisPt = histos.get<TH2>(HIST("Ring/DeltaMethod/h2dLambdaPtVsDeltaComp"))->GetXaxis();
    mAxisMass = histos.get<TH2>(HIST("Ring/DeltaMethod/h2dMassVsDeltaComp"))->GetXaxis();
    mAxisDTheta = histos.get<TH2>(HIST("Ring/DeltaMethod/h2dDeltaThetaVsDeltaComp"))->GetXaxis();
    for (auto const& tracker : {&trackRing, &trackRingKinCuts, &trackJetKinCuts, &trackJetLambdaKinCuts})
      tracker->resize(mAxisPt->GetNbins(), mAxisMass->GetNbins(), mAxisDTheta->GetNbins());
  }

  // Helper to get centrality (same from TableProducer, thanks to templating!):
  template <typename TCollision>
  auto getCentrality(TCollision const& collision)
  {
    if (centralityEstimator == kCentFT0M)
      return collision.centFT0M();
    else if (centralityEstimator == kCentFT0C)
      return collision.centFT0C();
    else if (centralityEstimator == kCentFV0A)
      return collision.centFV0A();
    return -1.f;
  }

  // Initializing a random number generator for the worker (for perpendicular-to-jet direction QAs):
  TRandom3 randomGen{0}; // 0 means we auto-seed from machine entropy. This is called once per device in the pipeline, so we should not see repeated seeds across workers
  std::mt19937 rng{std::random_device{}()};

  // Pre-computed values for helper below:
  const double smearSigma = 0.05 * jetR;
  // Normalized eta weights spanning -0.9 to 0.9 (46 bins) and phi weights for the 50 bins spanning [0, 2pi]:
  static constexpr std::array<double, 46> etaLeadPWeights = {{0.01782505198123039, 0.01826119427561306, 0.01890047073124532, 0.01942224199989093, 0.01993380780273602, 0.02047274597515178, 0.02094135547756474, 0.02140259654932778, 0.02178490245182078, 0.02218346916434517, 0.02252861343224298, 0.02278214932340838, 0.02297395476452691, 0.02311861709583109, 0.02322295246943318, 0.02329274166449468, 0.02335344516182264, 0.02335971904087711, 0.02340163522424806, 0.02352868368676468, 0.02345839093195849, 0.02295391718536531, 0.02306543716698383, 0.02293131780181040, 0.02265098126991631, 0.02318893563623931, 0.02322457978177088, 0.02316188601954564, 0.02308812992419636, 0.02305831097751334, 0.02300695504254397, 0.02296398287895598, 0.02286956017362988, 0.02274321888486162, 0.02254049413132752, 0.02233024817234042, 0.02204811997596179, 0.02170252191737713, 0.02130517220864903, 0.02086349641950970, 0.02036590755841183, 0.01987946074337928, 0.01934550418314029, 0.01879487530409137, 0.01812577211265362, 0.01766945805710696}};
  static constexpr std::array<double, 50> phiLeadPWeights = {{0.01907529231698144, 0.02044679008716434, 0.01948618554713157, 0.02046288887443206, 0.02142057576726765, 0.01961361841611185, 0.02174981627752354, 0.02160945937846856, 0.02027231236207667, 0.02153799273983672, 0.02107609996106984, 0.02001899849606885, 0.02196516947939817, 0.02047654705787587, 0.02059561382369167, 0.02148625027035289, 0.02001510925188416, 0.02183059661361331, 0.02111406548114694, 0.01881826666371129, 0.02112797285609031, 0.02034071592473790, 0.01968993216337670, 0.02126946766383166, 0.02025580366897253, 0.02061136834815962, 0.02083881238552183, 0.01994368379135331, 0.02046280212696892, 0.02131631148368759, 0.01967275608960357, 0.02064965975476278, 0.02155758091535052, 0.02012557837329328, 0.02084718400924722, 0.02065094443305587, 0.01969546187428681, 0.02136531489392276, 0.02084491659074202, 0.01970883494847568, 0.02080349992636310, 0.02098440611808174, 0.02159984658204099, 0.02045819959592145, 0.01952755742791563, 0.02166002908586709, 0.02017512590511229, 0.01932658087972190, 0.02108933627170306, 0.02019288778648293}};
  // Build discrete eta distribution for sampling:
  std::discrete_distribution<int> etaLeadPDist{etaLeadPWeights.begin(), etaLeadPWeights.end()}; // Will be passed as the etaDist variable
  std::discrete_distribution<int> phiLeadPDist{phiLeadPWeights.begin(), phiLeadPWeights.end()};

  /// \brief One jet proxy (leadP, leadJet or subJet). Bound by reference so applyProxyDistortion() edits in-place.
  struct ProxyState {
    bool& hasValidProxy; //! whether the proxy is valid. Re-evaluated against minPtThreshold after distortion
    float& pt;
    float& eta;
    float& phi;
    XYZVector& unitVec; //! the proxy direction as a unit vector
  };

  /// \brief The cache slots shared by forcePreviousJet and doMixedEventProxies, also bound by reference.
  /// \note  Fallback skips this event using ProxyState::hasValidProxy.
  /// \note  forcePreviousJet updates these here for the next collision, in a simplistic event event mixing approach.
  //         doMixedEventProxies instead overwrites the caller's cache fields with this collision's mixed proxy right before the call.
  struct ProxyCacheRef {
    bool& hadProxy;
    float& eta;
    float& phi;
  };

  /// \brief Applies whichever fakePolSwitches distortion is active to a jet-proxy direction, in place. No-op if none are on.
  /// \param proxy input/edited in-place: the proxy's kinematics, overwritten by whichever distortion is active.
  /// \param minPtThreshold pT threshold for re-evaluating hasValidProxy after distortion.
  /// \param cache input/edited in-place: the previous-jet/mixed-event caches.
  /// \param etaDist,phiDist,rng sampling distributions/generator for forceRandJet and forceDatalikeJet.
  /// \note  Shared across leadP/leadJet/subJet: the caller resolves which proxy-specific procedure (e.g., LUT for evtMixing) applies.
  // Helper to modify the jet direction for QA and for spurious signal baseline removal tests:
  inline void applyProxyDistortion(ProxyState proxy, float minPtThreshold, ProxyCacheRef cache,
                                   std::discrete_distribution<int>& etaDist, std::discrete_distribution<int>& phiDist, std::mt19937& rng)
  {
    if (!fakePolSwitches.forcePerpToJet && !fakePolSwitches.forceJetDirectionSmudge && !fakePolSwitches.forceRandJet && !fakePolSwitches.forcePreviousJet && !fakePolSwitches.forceDatalikeJet && !fakePolSwitches.doMixedEventProxies) [[likely]] {
      return; // Skip this function if none of the modifications are actually being executed!
    }

    // QA block -- Purposefully changing the jet direction (should kill signal, if any):
    if (fakePolSwitches.forcePerpToJet) {
      // First, we build a vector perpendicular to the jet by picking an arbitrary vector not parallel to the jet
      XYZVector perpVec;
      if (std::abs(proxy.unitVec.X()) > 0.99) {
        perpVec = XYZVector(-proxy.unitVec.Z(), 0., proxy.unitVec.X()).Unit(); // Cross product with Y-axis (0, 1, 0)
      } else {
        perpVec = XYZVector(0., proxy.unitVec.Z(), -proxy.unitVec.Y()).Unit(); // Cross product with X-axis (1, 0, 0)
      }

      // Now we rotate around the jet axis by a random angle, just to make sure we are not introducing a bias in the QA:
      // We will use Rodrigues' rotation formula (v_rot = v*cos(randomAngle) + (Jet \cross v)*sin(randomAngle))
      const double randomAngle = randomGen.Uniform(0., o2::constants::math::TwoPI);
      proxy.unitVec = perpVec * std::cos(randomAngle) + proxy.unitVec.Cross(perpVec) * std::sin(randomAngle);
    } else if (fakePolSwitches.forceJetDirectionSmudge) {
      // Smear the jet direction by a small random angle to estimate sensitivity to
      // jet axis uncertainty. We rotate the jet axis by angle theta around a uniformly
      // random perpendicular axis -- this is isotropic and coordinate-independent,
      // unlike smearing eta and phi separately (which would break azimuthal symmetry
      // around the jet axis and depend on where in eta the jet sits).

      // 1) We pick a uniformly random axis perpendicular to the jet.
      // (re-using the same Rodrigues formula as in the forcePerpToJet block above)
      XYZVector perpVec;
      if (std::abs(proxy.unitVec.X()) > 0.99) {
        perpVec = XYZVector(-proxy.unitVec.Z(), 0., proxy.unitVec.X()).Unit(); // Cross product with Y-axis (0, 1, 0)
      } else {
        perpVec = XYZVector(0., proxy.unitVec.Z(), -proxy.unitVec.Y()).Unit(); // Cross product with X-axis (1, 0, 0)
      }

      // Rotate perpVec around the jet axis by a uniform random azimuth to get
      // a uniformly distributed random perpendicular direction (the smear axis):
      const double smearAzimuth = randomGen.Uniform(0., o2::constants::math::TwoPI);
      XYZVector smearAxis = perpVec * std::cos(smearAzimuth) + proxy.unitVec.Cross(perpVec) * std::sin(smearAzimuth);

      // 2) draw the smearing polar angle from a Gaussian:
      // sigma = 0.05 * R --> ~68% of events smeared within 5% of R,
      //                      ~95% of events smeared within 10% of R,
      //                       ~5% see a displacement > 0.1*R (a very "badly determined jet", for our QA purposes)
      // std::abs() folds the symmetric Gaussian onto a half-normal ([0, inf))
      // -- R is not really an angle: just gives me a scale for the angular shift I am performing.
      // -- This may pose problems for forward jets: a small displacement in \theta becomes a large displacement in \eta space
      const double smearAngle = std::abs(randomGen.Gaus(0., smearSigma));

      // 3) rotate the jet axis by smearAngle around smearAxis.
      // Rodrigues is v_rot = v*cos(theta) + (k \cross v)*sin(theta) + k*(k \cdot v)*(1-cos(theta))
      // But the last term vanishes because smearAxis is perpendicular to unitVec:
      proxy.unitVec = proxy.unitVec * std::cos(smearAngle) + smearAxis.Cross(proxy.unitVec) * std::sin(smearAngle);
      // Also, rotation preserves the norm, so no re-normalisation is needed for this to be a unit vector.
    } else if (fakePolSwitches.forceRandJet) {
      // This randomization was made different for each proxy (LeadP, LeadJet, SubLeadJet): bear that in mind!
      // 1) Uniformly sample cos(theta) and phi to ensure an isotropic distribution (could also use TRandom::Sphere as well, but may be slower)
      // Notice that uniformly sampling theta would make the distribution non-isotropic, thus we use cos(theta)!
      const double cosTheta = randomGen.Uniform(-1., 1.);
      const double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
      const double randPhi = randomGen.Uniform(0., o2::constants::math::TwoPI);

      // 2) Construct the new random unit vector (there is no need to use the magnitude at all! We only need direction here):
      proxy.unitVec = XYZVector(sinTheta * std::cos(randPhi), sinTheta * std::sin(randPhi), cosTheta);
    } else if (fakePolSwitches.forcePreviousJet) { // Use the jet direction from the immediately preceding collision. The simplest event mixing
      const bool usableProxy = cache.hadProxy;
      if (usableProxy) {
        const double inverseCoshEta = 1.0 / std::cosh(cache.eta);
        const double sinPhi = std::sin(cache.phi);
        const double cosPhi = std::cos(cache.phi);
        proxy.unitVec = XYZVector(cosPhi * inverseCoshEta, sinPhi * inverseCoshEta, std::tanh(cache.eta));
      }

      // Update cache with the current event data:
      cache.hadProxy = proxy.hasValidProxy;
      cache.eta = proxy.eta;
      cache.phi = proxy.phi;

      // Current event cannot use previous-jet mixing if previous event lacked a proxy
      if (!usableProxy)
        proxy.hasValidProxy = false;
    } else if (fakePolSwitches.forceDatalikeJet) { // A compromise between forceRandJet and forcePreviousJet, using data-like weights for sampling jets
      const float etaMin = -0.92f;
      const float etaBinWidth = 0.04f;
      constexpr float phiBinWidth = constants::math::TwoPI / 50.f;

      // Pick one of the 46 bins according to etaWeights:
      const int binEtaIdx = etaDist(rng);
      const int binPhiIdx = phiDist(rng);

      // Uniformly smear inside the chosen bin:
      proxy.eta = etaMin + etaBinWidth * (binEtaIdx + std::generate_canonical<float, 24>(rng));
      proxy.phi = phiBinWidth * (binPhiIdx + std::generate_canonical<float, 24>(rng));

      const double inverseCoshEta = 1.0 / std::cosh(proxy.eta);
      const double sinPhi = std::sin(proxy.phi);
      const double cosPhi = std::cos(proxy.phi);
      proxy.unitVec = XYZVector(cosPhi * inverseCoshEta, sinPhi * inverseCoshEta, std::tanh(proxy.eta));
    } else if (fakePolSwitches.doMixedEventProxies) {
      // Mixes this proxy's direction from a real, similar other collision, rather than sampling a fitted distribution as datalikeJet.
      // Only unitVec is rebuilt, eta/phi are left as this collision's own real values -- same convention forcePreviousJet uses above.
      if (cache.hadProxy) {
        const double inverseCoshEta = 1.0 / std::cosh(cache.eta);
        const double sinPhi = std::sin(cache.phi);
        const double cosPhi = std::cos(cache.phi);
        proxy.unitVec = XYZVector(cosPhi * inverseCoshEta, sinPhi * inverseCoshEta, std::tanh(cache.eta));
      } else {
        // No mixing partner found for this collision (sparse bins!).
        // Same skip that forcePreviousJet does above: an event without a valid mixed proxy cannot be used.
        proxy.hasValidProxy = false;
      }
    }

    // Recalculating pT, phi and eta after distortions:
    // (without this, later kinematic selections make no sense at all! In the forceRandJet case, the ring observable would always sum zero)
    if (proxy.hasValidProxy) { // If you don't check this flag here, the "if (!usableProxy)" change would be silently overwritten
      if (!fakePolSwitches.forcePreviousJet && !fakePolSwitches.forceDatalikeJet) {
        // (In these two cases, we already know the eta and phi variables. No need to recompute)
        // Calculating total jet momentum, which should be preserved, just to recalculate the new jet Pt:
        const double mag = proxy.pt * std::cosh(proxy.eta);

        // For stability (Rho is the projection on the transverse plane, badly behaved for high |eta|):
        const double transverseNorm = std::max(proxy.unitVec.Rho(), 1e-12); // Stability guard

        // Recalculate pT:
        proxy.pt = mag * transverseNorm;

        // Recalculate phi:
        proxy.phi = RecoDecay::constrainAngle(std::atan2(proxy.unitVec.Y(), proxy.unitVec.X()), 0.0f); // atan2 outputs [-PI, PI), and DataModel convention was [0,2PI) as per FastJet's phi() getter
        // proxy.phi = std::atan2(proxy.unitVec.Y(), proxy.unitVec.X());
        // if (proxy.phi < 0.0f)
        // proxy.phi += o2::constants::math::TwoPI;

        // Stable eta computation:
        // Stabler than 0.5 * std::log((1. + cosTheta) / (1. - cosTheta))
        proxy.eta = std::asinh(proxy.unitVec.Z() / transverseNorm);
      }

      // Recalculate the bool after distortion to see if we proceed with this jet proxy:
      proxy.hasValidProxy = proxy.pt > minPtThreshold;
    }
  }

  // Caching the previous collision's jet directions -- An optional feature for forcePreviousJet QA:
  struct PrevJetCache {
    // Leading jet
    bool hadLeadJet;
    float leadJetEta;
    float leadJetPhi;
    // Subleading jet
    bool hadSubJet;
    float subJetEta;
    float subJetPhi;
    // Leading particle
    bool hadLeadP;
    float leadPEta;
    float leadPPhi;
  };

  // Per-collision leading/subleading jet, computed once per dataframe and shared by both the main loop below
  // (instead of re-scanning RingJets once per resampling pass) and the leadJet/subJet event mixing functions:
  struct JetProxyCache {
    bool hasValidLeadingJet = false;
    float leadingJetPt = -1.f;
    float leadingJetEta = 0.f;
    float leadingJetPhi = 0.f;
    bool hasValidSubJet = false;
    float subleadingJetPt = -1.f;
    float subleadingJetEta = 0.f;
    float subleadingJetPhi = 0.f;
  };

  // A simple struct for doMixedEventProxies. Each proxy (leadP/leadJet/subJet) gets its own independent cache.
  struct MixedProxyInfo {
    float eta;
    float phi;
    int64_t sourceCollisionId;
  };

  /// \brief Uniform random pick among each target's candidate window, via reservoir sampling
  ///        (cannibalization of proxies by neighbouring collisions in Continuous Readout is not a worry as ITS hits are being demanded)
  /// \note  Reuses the task's own rng member rather than a separate generator per proxy.
  void reservoirInsert(std::unordered_map<int64_t, int>& candidateCount, std::unordered_map<int64_t, MixedProxyInfo>& lut,
                       int64_t targetId, const MixedProxyInfo& candidate)
  {
    int& nSeen = candidateCount[targetId];
    ++nSeen;
    std::uniform_int_distribution<int> pick(1, nSeen);
    if (pick(rng) == 1) {
      lut[targetId] = candidate;
    }
  }

  /// \brief How many different target collisions each source collision ended up supplying, from a finished LUT.
  static std::unordered_map<int64_t, int> tallySourceUsage(std::unordered_map<int64_t, MixedProxyInfo> const& lut)
  {
    std::unordered_map<int64_t, int> usage;
    for (auto const& kv : lut)
      usage[kv.second.sourceCollisionId]++;
    return usage;
  }

  // Defining filters for events:
  Filter zvtxFilter = (nabs(o2::aod::lambdajetpol::zvtx) < maxZVtxPosition);

  // Preslices for correct collisions association:
  // (tested custom grouping and performs worse here)
  Preslice<aod::RingJets> perColJets = o2::aod::lambdajetpol::ringCollisionId;
  Preslice<aod::RingLaV0s> perColV0s = o2::aod::lambdajetpol::ringCollisionId;
  Preslice<aod::RingLeadPs> perColLeadPs = o2::aod::lambdajetpol::ringCollisionId;
  // For doMixedEventProxies:
  SliceCache mixCache;
  /// \brief Main analysis loop: for each collision, rebuilds the leading jet/particle proxies (with optional fakePolSwitches distortions), then loops over V0s computing the ring observable and polarization-vector profiles for every enabled kinematic-cut family.
  void processPolarizationData(soa::Filtered<o2::aod::RingCollisions> const& collisions, o2::aod::RingJets const& jets, o2::aod::RingLaV0s const& v0s,
                               o2::aod::RingLeadPs const& leadPs)
  {
    // Event mixing, built once per dataframe
    // Simplistic event mixing using previous collision:
    PrevJetCache prevJetCache{}; // Initializing struct before the collisions loop. This zero-initializes, so bools start as false
                                 // This should not be used along nProxyResamples > 1, as it does not apply to that case.

    // Leading/subleading jet per collision, resolved once here instead of inside the (possibly nProxyResamples times resampled) main loop.
    // Both the main loop and the leadJet/subJet mixing functions read from this, so RingJets is scanned exactly once
    std::unordered_map<int64_t, JetProxyCache> jetProxyByCollision;
    for (auto const& collision : collisions) {
      const auto collId = collision.globalIndex();
      JetProxyCache cache;
      // std::optional avoids undefined behaviour from a default-constructed iterator:
      std::optional<o2::aod::RingJets::iterator> leadingJet;
      std::optional<o2::aod::RingJets::iterator> subleadingJet;
      for (auto const& jet : jets.sliceBy(perColJets, collId)) {
        const auto jetpt = jet.jetPt();
        if (jetpt > cache.leadingJetPt) {
          // Current leading becomes subleading:
          cache.subleadingJetPt = cache.leadingJetPt;
          subleadingJet = leadingJet; // may still be std::nullopt on first pass -- that is fine!
          // Now update the leading jet:
          cache.leadingJetPt = jetpt;
          leadingJet = jet;
        } else if (jetpt > cache.subleadingJetPt) { // Update subleading only:
          cache.subleadingJetPt = jetpt;
          subleadingJet = jet;
        }
      }
      // Finer control on jet momentum, further than TableProducer pre-selection:
      cache.hasValidLeadingJet = cache.leadingJetPt > minLeadJetPt;
      cache.hasValidSubJet = cache.subleadingJetPt > minSubLeadJetPt;
      if (cache.hasValidLeadingJet) {
        cache.leadingJetEta = leadingJet->jetEta();
        cache.leadingJetPhi = leadingJet->jetPhi();
      }
      if (cache.hasValidSubJet) {
        cache.subleadingJetEta = subleadingJet->jetEta();
        cache.subleadingJetPhi = subleadingJet->jetPhi();
      }
      jetProxyByCollision[collId] = cache;
    }

    // doMixedEventProxies caches: three independent mixings, one per proxy (a collision may have a valid LeadP, but no SubLeadJet).
    // Mixing is binned on (Zvtx, proxy pt, centrality), varying only eta/phi.
    std::unordered_map<int64_t, MixedProxyInfo> mixedLeadPByCollision;
    std::unordered_map<int64_t, MixedProxyInfo> mixedLeadJetByCollision;
    std::unordered_map<int64_t, MixedProxyInfo> mixedSubJetByCollision;
    // First we build lookup tables based on current dataframe's collisions (connects pairs of jet proxies from similar collisions):
    // (these proxies may come from collisions with no valid Lambdas, by construction, enabling more mixes)
    // (This is performed out of the resampling loop, so nProxyResamples will not resample event mixing candidates)
    if (fakePolSwitches.doMixedEventProxies) {
      auto getMixCentrality = [this](o2::aod::RingCollision const& col) { return this->getCentrality(col); }; // Already filtered even though referencing only RingCollision, not an iterator to Filtered table

      // For the leading particle:
      // Pattern follows MixedEventsLambdaBinning tutorial (captures leadPs table and cache define in task's struct):
      auto getMixLeadPPt = [&leadPs, this](o2::aod::RingCollision const& col) {
        auto rows = leadPs.sliceByCached(o2::aod::lambdajetpol::ringCollisionId, col.globalIndex(), this->mixCache); // As it is cached, grouping happens only once
        return rows.size() > 0 ? rows.begin().leadParticlePt() : -999.f;
      };
      using LeadPBinningType = FlexibleBinningPolicy<std::tuple<decltype(getMixLeadPPt), decltype(getMixCentrality)>,
                                                     o2::aod::lambdajetpol::Zvtx, decltype(getMixLeadPPt), decltype(getMixCentrality)>;
      LeadPBinningType leadPBinning{{getMixLeadPPt, getMixCentrality},
                                    {axisConfigurations.axisPVz, axisConfigurations.axisJetPt, axisConfigurations.axisCentrality},
                                    true}; // Ignore overflows true

      // SameKindPair defaults to CombinationsBlockStrictlyUpperSameIndexPolicy, so no same-event mixing should happen:
      //(Already filtered by Zvtx even though we call by aod::RingCollisions, so no need to access the filtered table)
      SameKindPair<o2::aod::RingCollisions, o2::aod::RingLeadPs, LeadPBinningType> leadPPair{
        leadPBinning, fakePolSwitches.mixedEventWindowSize, -1, collisions, std::make_tuple(leadPs), &mixCache};

      std::unordered_map<int64_t, int> leadPCandidateCount;
      for (auto it = leadPPair.begin(); it != leadPPair.end(); ++it) {
        auto& [c1, leadP1, c2, leadP2] = *it;         // Iterates over collision pairs and leading particle pairs (structured binding)
        if (leadP1.size() > 0 && leadP2.size() > 0) { // There should always be at least one leadP, given the overflow exclusion above
          float eta1 = 0.f, phi1 = 0.f, eta2 = 0.f, phi2 = 0.f;
          for (auto const& lp : leadP1) {
            eta1 = lp.leadParticleEta();
            phi1 = lp.leadParticlePhi();
            break;
          } // Retrieves the first entry
          for (auto const& lp : leadP2) {
            eta2 = lp.leadParticleEta();
            phi2 = lp.leadParticlePhi();
            break;
          }
          // Each side of the pair is one more candidate for the other collision's reservoir:
          reservoirInsert(leadPCandidateCount, mixedLeadPByCollision, c1.globalIndex(), {eta2, phi2, c2.globalIndex()});
          reservoirInsert(leadPCandidateCount, mixedLeadPByCollision, c2.globalIndex(), {eta1, phi1, c1.globalIndex()});
          histos.fill(HIST("EventMixingQA/hMixedEventLeadPDeltaIndexEligible"), c2.globalIndex() - c1.globalIndex());
        }
        if (it.isNewWindow()) { // Count each bin-window once, not once per pair inside it
          histos.fill(HIST("EventMixingQA/hMixedEventLeadPWindowNeighbours"), it.currentWindowNeighbours());
        }
      }

      // Leading jet:
      // pt comes from the jetProxyByCollision pre-pass above.
      // (Collisions whose leading jet failed minLeadJetPt get the -999 sentinel, so they never enter this mixing)
      auto getMixLeadJetPt = [&jetProxyByCollision](o2::aod::RingCollision const& col) {
        auto cacheIt = jetProxyByCollision.find(col.globalIndex());
        if (cacheIt == jetProxyByCollision.end() || !cacheIt->second.hasValidLeadingJet)
          return -999.f;
        return cacheIt->second.leadingJetPt;
      };
      using LeadJetBinningType = FlexibleBinningPolicy<std::tuple<decltype(getMixLeadJetPt), decltype(getMixCentrality)>,
                                                       o2::aod::lambdajetpol::Zvtx, decltype(getMixLeadJetPt), decltype(getMixCentrality)>;
      LeadJetBinningType leadJetBinning{{getMixLeadJetPt, getMixCentrality},
                                        {axisConfigurations.axisPVz, axisConfigurations.axisJetPt, axisConfigurations.axisCentrality},
                                        true};

      // RingJets is still the associated table (SameKindPair requires one), but its sliced content goes unused here:
      // the borrowed direction comes from jetProxyByCollision, which already knows which jet is the leading one.
      SameKindPair<o2::aod::RingCollisions, o2::aod::RingJets, LeadJetBinningType> leadJetPair{
        leadJetBinning, fakePolSwitches.mixedEventWindowSize, -1, collisions, std::make_tuple(jets), &mixCache};

      std::unordered_map<int64_t, int> leadJetCandidateCount;
      for (auto it = leadJetPair.begin(); it != leadJetPair.end(); ++it) {
        auto& [c1, jets1, c2, jets2] = *it; // jets1/jets2 intentionally unused: the leading/subleading jet is already resolved in jetProxyByCollision
        auto cachedLeadJet1 = jetProxyByCollision.find(c1.globalIndex());
        auto cachedLeadJet2 = jetProxyByCollision.find(c2.globalIndex());
        if (cachedLeadJet1 != jetProxyByCollision.end() && cachedLeadJet2 != jetProxyByCollision.end() &&
            cachedLeadJet1->second.hasValidLeadingJet && cachedLeadJet2->second.hasValidLeadingJet) {
          reservoirInsert(leadJetCandidateCount, mixedLeadJetByCollision, c1.globalIndex(),
                          {cachedLeadJet2->second.leadingJetEta, cachedLeadJet2->second.leadingJetPhi, c2.globalIndex()});
          reservoirInsert(leadJetCandidateCount, mixedLeadJetByCollision, c2.globalIndex(),
                          {cachedLeadJet1->second.leadingJetEta, cachedLeadJet1->second.leadingJetPhi, c1.globalIndex()});
          histos.fill(HIST("EventMixingQA/hMixedEventLeadJetDeltaIndexEligible"), c2.globalIndex() - c1.globalIndex());
        }
        if (it.isNewWindow()) {
          histos.fill(HIST("EventMixingQA/hMixedEventLeadJetWindowNeighbours"), it.currentWindowNeighbours());
        }
      }

      // Subleading jet:
      auto getMixSubJetPt = [&jetProxyByCollision](o2::aod::RingCollision const& col) {
        auto cacheIt = jetProxyByCollision.find(col.globalIndex());
        if (cacheIt == jetProxyByCollision.end() || !cacheIt->second.hasValidSubJet)
          return -999.f;
        return cacheIt->second.subleadingJetPt;
      };
      using SubJetBinningType = FlexibleBinningPolicy<std::tuple<decltype(getMixSubJetPt), decltype(getMixCentrality)>,
                                                      o2::aod::lambdajetpol::Zvtx, decltype(getMixSubJetPt), decltype(getMixCentrality)>;
      SubJetBinningType subJetBinning{{getMixSubJetPt, getMixCentrality},
                                      {axisConfigurations.axisPVz, axisConfigurations.axisJetPt, axisConfigurations.axisCentrality},
                                      true};

      SameKindPair<o2::aod::RingCollisions, o2::aod::RingJets, SubJetBinningType> subJetPair{
        subJetBinning, fakePolSwitches.mixedEventWindowSize, -1, collisions, std::make_tuple(jets), &mixCache};

      std::unordered_map<int64_t, int> subJetCandidateCount;
      for (auto it = subJetPair.begin(); it != subJetPair.end(); ++it) {
        auto& [c1, jets1, c2, jets2] = *it; // jets1/jets2 intentionally unused: the leading/subleading jet is already resolved in jetProxyByCollision
        auto cachedSubJet1 = jetProxyByCollision.find(c1.globalIndex());
        auto cachedSubJet2 = jetProxyByCollision.find(c2.globalIndex());
        if (cachedSubJet1 != jetProxyByCollision.end() && cachedSubJet2 != jetProxyByCollision.end() &&
            cachedSubJet1->second.hasValidSubJet && cachedSubJet2->second.hasValidSubJet) {
          reservoirInsert(subJetCandidateCount, mixedSubJetByCollision, c1.globalIndex(),
                          {cachedSubJet2->second.subleadingJetEta, cachedSubJet2->second.subleadingJetPhi, c2.globalIndex()});
          reservoirInsert(subJetCandidateCount, mixedSubJetByCollision, c2.globalIndex(),
                          {cachedSubJet1->second.subleadingJetEta, cachedSubJet1->second.subleadingJetPhi, c1.globalIndex()});
          histos.fill(HIST("EventMixingQA/hMixedEventSubJetDeltaIndexEligible"), c2.globalIndex() - c1.globalIndex());
        }
        if (it.isNewWindow()) {
          histos.fill(HIST("EventMixingQA/hMixedEventSubJetWindowNeighbours"), it.currentWindowNeighbours());
        }
      }

      // Selected-partner proximity: |target - source| for the partner the reservoir actually kept.
      for (auto const& kv : mixedLeadPByCollision)
        histos.fill(HIST("EventMixingQA/hMixedEventLeadPDeltaIndexSelected"), std::abs(kv.first - kv.second.sourceCollisionId));
      for (auto const& kv : mixedLeadJetByCollision)
        histos.fill(HIST("EventMixingQA/hMixedEventLeadJetDeltaIndexSelected"), std::abs(kv.first - kv.second.sourceCollisionId));
      for (auto const& kv : mixedSubJetByCollision)
        histos.fill(HIST("EventMixingQA/hMixedEventSubJetDeltaIndexSelected"), std::abs(kv.first - kv.second.sourceCollisionId));

      // Source-usage QA:
      // (Per-proxy (HIST() needs literal names, so the fills stay unrolled here)
      auto leadPUsage = tallySourceUsage(mixedLeadPByCollision);
      for (auto const& kv : leadPUsage)
        histos.fill(HIST("EventMixingQA/hMixedEventLeadPSourceUsageCount"), kv.second);
      // Eta/phi of a mixed proxy is stored under the collision that receives the mixing -- reuse the first hit found:
      for (auto const& kv : mixedLeadPByCollision) {
        auto usageIt = leadPUsage.find(kv.second.sourceCollisionId);
        if (usageIt == leadPUsage.end()) // already consumed below
          continue;
        histos.fill(HIST("EventMixingQA/pMixedEventLeadPSourceUsageVsEta"), kv.second.eta, usageIt->second);
        histos.fill(HIST("EventMixingQA/pMixedEventLeadPSourceUsageVsPhi"), kv.second.phi, usageIt->second);
        leadPUsage.erase(usageIt); // Fill this source exactly once, not once per target it supplied
      }

      auto leadJetUsage = tallySourceUsage(mixedLeadJetByCollision);
      for (auto const& kv : leadJetUsage)
        histos.fill(HIST("EventMixingQA/hMixedEventLeadJetSourceUsageCount"), kv.second);
      for (auto const& kv : mixedLeadJetByCollision) {
        auto usageIt = leadJetUsage.find(kv.second.sourceCollisionId);
        if (usageIt == leadJetUsage.end())
          continue;
        histos.fill(HIST("EventMixingQA/pMixedEventLeadJetSourceUsageVsEta"), kv.second.eta, usageIt->second);
        histos.fill(HIST("EventMixingQA/pMixedEventLeadJetSourceUsageVsPhi"), kv.second.phi, usageIt->second);
        leadJetUsage.erase(usageIt);
      }

      auto subJetUsage = tallySourceUsage(mixedSubJetByCollision);
      for (auto const& kv : subJetUsage)
        histos.fill(HIST("EventMixingQA/hMixedEventSubJetSourceUsageCount"), kv.second);
      for (auto const& kv : mixedSubJetByCollision) {
        auto usageIt = subJetUsage.find(kv.second.sourceCollisionId);
        if (usageIt == subJetUsage.end())
          continue;
        histos.fill(HIST("EventMixingQA/pMixedEventSubJetSourceUsageVsEta"), kv.second.eta, usageIt->second);
        histos.fill(HIST("EventMixingQA/pMixedEventSubJetSourceUsageVsPhi"), kv.second.phi, usageIt->second);
        subJetUsage.erase(usageIt);
      }
    }

    for (int idxResampling = 0; idxResampling < fakePolSwitches.nProxyResamples; idxResampling++) { // resampling loop for forceRandJet and forceDatalikeJet
      for (auto const& collision : collisions) {
        const float collisionPVz = collision.zvtx();

        const auto collId = collision.globalIndex(); // The self-index accessor
        const float centrality = getCentrality(collision);

        // Used this dummy for backwards compatibility, under the reasonable assumption that the field points always in the same *direction* in the used runs
        // (it is not worth it to fetch and store the magnetic field in the datamodel)
        const float magField = 1.f; // Purely geometric.

        // Slice jets, V0s and leading particle belonging to this collision:
        // (global collision indices repeat a lot, but they are unique to a same TimeFrame (TF) subfolder in the derived data)
        auto v0sInColl = v0s.sliceBy(perColV0s, collId);
        auto leadPsInColl = leadPs.sliceBy(perColLeadPs, collId);

        // Check if there is at least one V0 and one jet in the collision:
        // (in the way I fill the table, there is always at least one V0 in
        //  the stored collision, but the jets table can not be filled for
        //  that collision, and a collision may not be filled when the jets
        //  table is. Be mindful of that!)
        // 1) Require at least one V0:
        const int nLambdaLikeV0s = v0sInColl.size(); // Caching this variable, as it will be reused in the loop
                                                     // In the latest datamodel format, only unambiguous V0s (Lambda XOR antiLambda) are saved,
                                                     // so the number of V0s in the collision table is the number of Lambdas/antiLambda identified
                                                     // by the table producer.
        if (!nLambdaLikeV0s)
          continue;

        // 2) We require at least one leading particle:
        // (The goal is to see how diluted the signal gets with events which don't even have a loose FastJet jet)
        // (The leading particle is built from all tracks that passed the pseudojet
        // selection, so it exists whenever FastJet was run on this collision.
        // Events that have a leading jet always have a leading particle too, but
        // the converse is not true: events can have a leading particle with no jet
        // if no jet survives the pT threshold/the background subtraction)
        // (At least that is the case when minLeadParticlePt = 0)
        float leadPPt = -1.; // pT = -1 means "table entry not found for this collision".
        float leadPEta = 0.;
        float leadPPhi = 0.;
        float leadPPx = 0., leadPPy = 0., leadPPz = 0.;
        for (auto const& lp : leadPsInColl) {
          // Table should contain exactly one entry per collision, but we break immediately to be safe:
          leadPPt = lp.leadParticlePt();
          leadPEta = lp.leadParticleEta();
          leadPPhi = lp.leadParticlePhi();
          // Using dynamic columns to make code cleaner:
          leadPPx = lp.leadParticlePx();
          leadPPy = lp.leadParticlePy();
          leadPPz = lp.leadParticlePz();
        }
        // // Discard events with no leading particle (FastJet didn't even run in these cases!):
        // if (leadPPt < 0.)
        //   continue;

        // Apply minimum pT selection for the leading particle (not necessarily the same as in derived data builder. Can be a stricter cut!):
        bool hasValidLeadingP = leadPPt > minLeadParticlePt;

        // Build leading particle unit vector, outside the V0 loop for performance.
        XYZVector leadPUnitVec(1., 0., 0.); // dummy (overwritten below when hasValidLeadingP)
        if (hasValidLeadingP) {
          leadPUnitVec = XYZVector(leadPPx, leadPPy, leadPPz).Unit();

          // Apply distortion logic:
          const bool hadPreviousProxy = prevJetCache.hadLeadP;

          // doMixedEventProxies: get this collision's mixing partner (if any) from the LUT built above:
          // (this check is performed only if hasValidLeadingP is true for performance, as the LeadPPt matching for the event mixing would also demand a minimum jet pT)
          if (fakePolSwitches.doMixedEventProxies) {
            auto itMix = mixedLeadPByCollision.find(collId);
            prevJetCache.hadLeadP = (itMix != mixedLeadPByCollision.end()); // Check if this event had a valid mixing target. Reusing prevJetCache
            if (prevJetCache.hadLeadP) {
              prevJetCache.leadPEta = itMix->second.eta;
              prevJetCache.leadPPhi = itMix->second.phi;
            }
            histos.fill(HIST("EventMixingQA/hMixedEventLeadPOutcome"), prevJetCache.hadLeadP ? 1 : 0);
            if (prevJetCache.hadLeadP) { // Only fill the comparison histograms on an actual hit:
              histos.fill(HIST("EventMixingQA/h2dMixedLeadPEtaVsLeadPEta"), prevJetCache.leadPEta, leadPEta);
              histos.fill(HIST("EventMixingQA/h2dMixedLeadPPhiVsLeadPPhi"), prevJetCache.leadPPhi, leadPPhi);
            }
          }

          applyProxyDistortion({hasValidLeadingP, leadPPt, leadPEta, leadPPhi, leadPUnitVec},
                               minLeadParticlePt, {prevJetCache.hadLeadP, prevJetCache.leadPEta, prevJetCache.leadPPhi},
                               etaLeadPDist, phiLeadPDist, rng);

          // Fill distorted-proxy QA histograms:
          // Do not gate on the post-distortion hasValidLeadingP (a pT cut) value here!
          if (!fakePolSwitches.forcePreviousJet || hadPreviousProxy) {
            histos.fill(HIST("JetKinematicsQA/hLeadPEta"), leadPEta);
            histos.fill(HIST("JetKinematicsQA/hLeadPPhi"), leadPPhi);
            histos.fill(HIST("JetKinematicsQA/hJetCounterPtLeadP"), leadPPt);

            histos.fill(HIST("JetKinematicsQA/h2dLeadPEtaVsPVz"), leadPEta, collisionPVz);
            if (doJetProxy5dQA)
              histos.fill(HIST("JetKinematicsQA/h5dLeadPEtaPhiPtPVzCent"), leadPEta, leadPPhi, leadPPt, collisionPVz, centrality);
          }
        }

        // 3) Fetching leading jet and subleading jet -- Resolved once per collision in the jetProxyByCollision pre-pass:
        const JetProxyCache& jetProxies = jetProxyByCollision.at(collId); // .at() makes it sure we don't create a new key in the map
        float leadingJetPt = jetProxies.leadingJetPt;
        float subleadingJetPt = jetProxies.subleadingJetPt;

        // Defining local bools that may be changed by applyProxyDistortion:
        bool hasValidLeadingJet = jetProxies.hasValidLeadingJet;
        bool hasValidSubJet = jetProxies.hasValidSubJet;

        // Build jet vectors (only when the corresponding jet exists):
        // Dummy initialisations are safe: all jet-dependent fills are gated on hasValidLeadingJet / hasValidSubJet.
        float leadingJetEta = 0.;
        float leadingJetPhi = 0.;
        XYZVector leadingJetUnitVec(1., 0., 0.); // dummy (overwritten below)
        if (hasValidLeadingJet) {
          leadingJetEta = jetProxies.leadingJetEta;
          leadingJetPhi = jetProxies.leadingJetPhi;
          // Rebuild the direction from the cached eta/phi (cheaper than calling jetPx() internal getters and then normalizing with .Unit()):
          const double inverseCoshEta = 1.0 / std::cosh(leadingJetEta);
          leadingJetUnitVec = XYZVector(std::cos(leadingJetPhi) * inverseCoshEta, std::sin(leadingJetPhi) * inverseCoshEta, std::tanh(leadingJetEta));

          // Apply distortion logic:
          const bool hadPreviousProxy = prevJetCache.hadLeadJet;
          if (fakePolSwitches.doMixedEventProxies) {
            // Get this collision's leading-jet mixing partner (if any) from the LUT built above:
            auto itMix = mixedLeadJetByCollision.find(collId);
            prevJetCache.hadLeadJet = (itMix != mixedLeadJetByCollision.end());
            if (prevJetCache.hadLeadJet) {
              prevJetCache.leadJetEta = itMix->second.eta;
              prevJetCache.leadJetPhi = itMix->second.phi;
            }
            histos.fill(HIST("EventMixingQA/hMixedEventLeadJetOutcome"), prevJetCache.hadLeadJet ? 1 : 0);
            if (prevJetCache.hadLeadJet) {
              histos.fill(HIST("EventMixingQA/h2dMixedLeadJetEtaVsLeadJetEta"), prevJetCache.leadJetEta, leadingJetEta);
              histos.fill(HIST("EventMixingQA/h2dMixedLeadJetPhiVsLeadJetPhi"), prevJetCache.leadJetPhi, leadingJetPhi);
            }
          }
          applyProxyDistortion({hasValidLeadingJet, leadingJetPt, leadingJetEta, leadingJetPhi, leadingJetUnitVec},
                               minLeadJetPt, {prevJetCache.hadLeadJet, prevJetCache.leadJetEta, prevJetCache.leadJetPhi},
                               etaLeadPDist, phiLeadPDist, rng);

          // Fill distorted-proxy QA histograms:
          // Do not gate on the post-distortion hasValidLeadingJet (a pT cut) value here!
          if (!fakePolSwitches.forcePreviousJet || hadPreviousProxy) {
            histos.fill(HIST("JetKinematicsQA/hLeadJetEta"), leadingJetEta);
            histos.fill(HIST("JetKinematicsQA/hLeadJetPhi"), leadingJetPhi);
            histos.fill(HIST("JetKinematicsQA/hJetCounterPtJet"), leadingJetPt);

            histos.fill(HIST("JetKinematicsQA/h2dLeadJetEtaVsPVz"), leadingJetEta, collisionPVz);
            if (doJetProxy5dQA)
              histos.fill(HIST("JetKinematicsQA/h5dLeadJetEtaPhiPtPVzCent"), leadingJetEta, leadingJetPhi, leadingJetPt, collisionPVz, centrality);
          }
        }

        float subleadingJetEta = 0.;
        float subleadingJetPhi = 0.;
        XYZVector subJetUnitVec(1., 0., 0.);
        if (hasValidSubJet) {
          subleadingJetEta = jetProxies.subleadingJetEta;
          subleadingJetPhi = jetProxies.subleadingJetPhi;
          const double inverseCoshEtaSub = 1.0 / std::cosh(subleadingJetEta);
          subJetUnitVec = XYZVector(std::cos(subleadingJetPhi) * inverseCoshEtaSub, std::sin(subleadingJetPhi) * inverseCoshEtaSub, std::tanh(subleadingJetEta));

          // Apply distortion logic:
          const bool hadPreviousProxy = prevJetCache.hadSubJet;
          if (fakePolSwitches.doMixedEventProxies) {
            // Get this collision's subleading-jet mixing partner (if any) from the LUT built above:
            auto itMix = mixedSubJetByCollision.find(collId);
            prevJetCache.hadSubJet = (itMix != mixedSubJetByCollision.end());
            if (prevJetCache.hadSubJet) {
              prevJetCache.subJetEta = itMix->second.eta;
              prevJetCache.subJetPhi = itMix->second.phi;
            }
            histos.fill(HIST("EventMixingQA/hMixedEventSubJetOutcome"), prevJetCache.hadSubJet ? 1 : 0);
            if (prevJetCache.hadSubJet) {
              histos.fill(HIST("EventMixingQA/h2dMixedSubJetEtaVsSubJetEta"), prevJetCache.subJetEta, subleadingJetEta);
              histos.fill(HIST("EventMixingQA/h2dMixedSubJetPhiVsSubJetPhi"), prevJetCache.subJetPhi, subleadingJetPhi);
            }
          }
          applyProxyDistortion({hasValidSubJet, subleadingJetPt, subleadingJetEta, subleadingJetPhi, subJetUnitVec},
                               minSubLeadJetPt, {prevJetCache.hadSubJet, prevJetCache.subJetEta, prevJetCache.subJetPhi},
                               etaLeadPDist, phiLeadPDist, rng);

          // Fill distorted-proxy QA histograms:
          // Do not gate on the post-distortion hasValidSubJet (a pT cut) value here!
          if (!fakePolSwitches.forcePreviousJet || hadPreviousProxy) {
            histos.fill(HIST("JetKinematicsQA/hSubLeadJetEta"), subleadingJetEta);
            histos.fill(HIST("JetKinematicsQA/hSubLeadJetPhi"), subleadingJetPhi);
            histos.fill(HIST("JetKinematicsQA/hJetCounterPt2ndJet"), subleadingJetPt);

            histos.fill(HIST("JetKinematicsQA/h2dSubLeadJetEtaVsPVz"), subleadingJetEta, collisionPVz);
            if (doJetProxy5dQA)
              histos.fill(HIST("JetKinematicsQA/h5dSubLeadJetEtaPhiPtPVzCent"), subleadingJetEta, subleadingJetPhi, subleadingJetPt, collisionPVz, centrality);
          }
        }

        // (jet eta cuts only meaningful when the jet actually exists)
        const bool kinematicJetCheck = hasValidLeadingJet && (std::abs(leadingJetEta) < 0.5);
        const bool kinematic2ndJetCheck = hasValidSubJet && (std::abs(subleadingJetEta) < 0.5);
        const bool kinematicLeadPCheck = hasValidLeadingP && (std::abs(leadPEta) < 0.5);

        // Quick bools that are useful for detector asymmetry QA:
        const bool jetEtaPos = hasValidLeadingJet && (leadingJetEta >= 0.); // Only perform >= check if has validJet
        const bool subJetEtaPos = hasValidSubJet && (subleadingJetEta >= 0.);
        const bool leadPEtaPos = hasValidLeadingP && (leadPEta >= 0.);

        // Stricter QA version of the bools -- Jets have a radius that makes it possible eta_{jet} > 0, yet half its tracks are in eta < 0
        // (This does not apply to leading particles, obviously. They have no substructure in eta)
        const bool jetEtaStrict = hasValidLeadingJet && (std::abs(leadingJetEta) >= jetR);
        const bool subJetEtaStrict = hasValidSubJet && (std::abs(subleadingJetEta) >= jetR);
        // If one was to define bools for each side of the detector (not needed in the current if-else structure on TProfile fills)
        // const bool jetEtaStrictPos = jetEtaPos && jetEtaStrict;
        // const bool jetEtaStrictNeg = !jetEtaPos && jetEtaStrict;
        // const bool subJetEtaStrictPos = subJetEtaPos && subJetEtaStrict;
        // const bool subJetEtaStrictNeg = !subJetEtaPos && subJetEtaStrict;

        // Fetching number of Lambda-like V0s in collision (must be known before full loop, to fill "pRingVsNV0s"):
        // int nLambdaLikeV0s = 0;
        // for (auto const& v0 : v0sInColl) {
        //   if (v0.isLambda() ^ v0.isAntiLambda()){ // XOR (only the non-ambiguous candidates)
        //     nLambdaLikeV0s++;
        //   }
        // }
        // Code above was superseeded: new datamodel does not store ambiguous candidates!
        // The new getter comes at the very start of processPolarizationData() now.

        // Initialize delta method accumulators (reset for new collision):
        for (auto const& tracker : {&trackRing, &trackRingKinCuts, &trackJetKinCuts, &trackJetLambdaKinCuts})
          tracker->reset();
        for (auto const& v0 : v0sInColl) {
          const bool isLambda = v0.isLambda(); // true: is a Lambda. false: is an antiLambda.
          // For now, removing the ambiguous candidates from the analysis. New datamodel does NOT save ambiguous candidates.
          // (From Podolanski-Armenteros plots, the population of ambiguous is ~3.8% without TOF, and without
          //  competing mass rejection. From those, ~99% seem to be K0s, so no real gain in considering the
          //  ambiguous candidates in the analysis)
          // const bool isAntiLambda = v0.isAntiLambda(); // No longer used!
          // if (isLambda && isAntiLambda) continue;
          const float v0pt = v0.v0Pt();
          const float v0eta = v0.v0Eta();
          const float v0phi = v0.v0Phi();
          const float dcaDau = v0.dcaV0Daughters();

          float v0LambdaLikeMass = v0.massV0();
          float protonLikePt = 0;
          float protonLikeEta = 0;
          float protonLikePhi = 0;
          float protonLikeDCADauToPV = 0;
          float pionLikeDCADauToPV = 0;
          if (isLambda) {
            if (!analyseLambda)
              continue;
            protonLikePt = v0.posPt();
            protonLikeEta = v0.posEta();
            protonLikePhi = v0.posPhi();
            protonLikeDCADauToPV = v0.dcaPosToPV();
            pionLikeDCADauToPV = v0.dcaNegToPV();
          } else { // Guaranteed to be an antiLambda candidate, not an ambiguous candidate
            if (!analyseAntiLambda)
              continue;
            protonLikePt = v0.negPt();
            protonLikeEta = v0.negEta();
            protonLikePhi = v0.negPhi();
            protonLikeDCADauToPV = v0.dcaNegToPV();
            pionLikeDCADauToPV = v0.dcaPosToPV();
          }

          PtEtaPhiMVector lambdaLike4Vec(v0pt, v0eta, v0phi, v0LambdaLikeMass);
          PtEtaPhiMVector protonLike4Vec(protonLikePt, protonLikeEta, protonLikePhi, ProtonMass);
          const float lambdaRapidity = lambdaLike4Vec.Rapidity();                                    // For further kinematic selections
          const int v0InMassPeak = (v0LambdaLikeMass >= 1.1020593 && v0LambdaLikeMass <= 1.1288811); // Very naive estimator. Based on signal extractions from outside this code

          // Inexpensive estimates of signal extraction effects on the observable:
          if (excludeOutOfPeakQA && !v0InMassPeak)
            continue;
          else if (excludeInPeakQA && v0InMassPeak)
            continue;

          // Boosting proton into lambda frame:
          XYZVector beta = lambdaLike4Vec.BoostToCM(); // Boost trivector that goes from laboratory frame to Lambda's rest frame (convenient new function, different from TLorentzVector's BoostVector())
          auto protonLike4VecStar = ROOT::Math::VectorUtil::boost(protonLike4Vec, beta);

          // Getting unit vectors and 3-components:
          XYZVector lambdaLike3Vec = lambdaLike4Vec.Vect();
          auto lambdaLikeUnit3Vec = lambdaLike3Vec.Unit();
          XYZVector protonLikeStarUnit3Vec = protonLike4VecStar.Vect().Unit();

          // Lab-frame Lambda momentum components -- Not for polarization, but for actual momenta plotting in XY and ZX planes:
          // (Used for the (px,py) and (pz,px) polarization vector-field / ring 2D profiles)
          const float v0px = lambdaLike3Vec.X();
          const float v0py = lambdaLike3Vec.Y();
          const float v0pz = lambdaLike3Vec.Z();

          // Calculating fake polarization ("negative helicity problem") estimator:
          // (this estimator is calculated outside of any gate, as it does not depend on jet proxy used)
          float cosFakePol = protonLikeStarUnit3Vec.Dot(lambdaLikeUnit3Vec);

          // Calculating the azimuthal angle between the Lambda and the proton:
          float deltaPhiLambdaProtonStar = wrapToPiFast(lambdaLikeUnit3Vec.Phi() - protonLikeStarUnit3Vec.Phi()); // Phi is defined from -PI to PI in ROOT::Math::Cartesian3D, thus kept the wrapping

          // Calculating the phi* angle:
          // e_z = p_Lambda_hat; // e_x = normalize(z_hat cross p_Lambda); // e_y = e_z cross e_x;
          // // phi_star = atan2(p_p_star dot e_y, p_p_star dot e_x);
          // XYZVector e_x(-lambdaLikeUnit3Vec.Y(), lambdaLikeUnit3Vec.X(), 0.); // Same as e_x = zHat.Cross(lambdaLikeUnit3Vec);
          // XYZVector e_y = lambdaLikeUnit3Vec.Cross(e_x); // e_y completes the right-handed coordinate system (e_z is lambdaLikeUnit3Vec)
          // float pX = protonLikeStarUnit3Vec.Dot(e_x);
          // float pY = protonLikeStarUnit3Vec.Dot(e_y);
          // float phiStar = std::atan2(pY, pX);
          // Faster implementation:
          // pX = p_y * L_x - p_x * L_y
          float pX = protonLikeStarUnit3Vec.Y() * lambdaLikeUnit3Vec.X() - protonLikeStarUnit3Vec.X() * lambdaLikeUnit3Vec.Y();
          // pY = p_z - L_z * (p_proton_star dot p_lambda_hat)
          float pY = protonLikeStarUnit3Vec.Z() - lambdaLikeUnit3Vec.Z() * cosFakePol; // (Reusing cosFakePol calculated earlier!)
          float phiStar = std::atan2(pY, pX);                                          // This will give an output from -PI to PI

          if (qaSwitches.doFakePolDiagnosticsQA) {
            // Another reconstruction efficiency measure:
            // (Formula is: p_{Lambda} \cross p_{Daughter}^{*} \cdot B, and B points in Z)
            if (analyseMagField) {
              auto crossGeom = lambdaLike3Vec.Cross(protonLikeStarUnit3Vec);
              const bool positiveGeom = crossGeom.Z() * magField > 0;

              if (isLambda && positiveGeom)
                histos.fill(HIST("HelicityEfficiencyQA/hLambdaMassDecayGeomRight"), v0LambdaLikeMass);
              else if (isLambda && !positiveGeom)
                histos.fill(HIST("HelicityEfficiencyQA/hLambdaMassDecayGeomLeft"), v0LambdaLikeMass);
              else if (!isLambda && positiveGeom)
                histos.fill(HIST("HelicityEfficiencyQA/hAntiLambdaMassDecayGeomRight"), v0LambdaLikeMass);
              else
                histos.fill(HIST("HelicityEfficiencyQA/hAntiLambdaMassDecayGeomLeft"), v0LambdaLikeMass);
            }

            // Measuring the AEE effect differentially (azimuthal efficiency effect, which causes different V0 topologies to be enhanced/suppressed):
            if (isLambda)
              histos.fill(HIST("HelicityEfficiencyQA/hLambdaMassVsPhiLambdaMinusPhiProtonStar"), v0LambdaLikeMass, deltaPhiLambdaProtonStar);
            else
              histos.fill(HIST("HelicityEfficiencyQA/hAntiLambdaMassVsPhiLambdaMinusPhiProtonStar"), v0LambdaLikeMass, deltaPhiLambdaProtonStar);
            // AEE and HEE correlation:
            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCounts_CosThetaVsPhiStar"), cosFakePol, phiStar);
          } // end doFakePolDiagnosticsQA (per-V0 AEE/HEE)

          // Useful kinematic bools:
          const bool lambdaEtaPos = v0eta >= 0.;
          const bool pTLambdaCheck = v0pt > 0.5 && v0pt < 1.5;
          const bool rapidityLambdaCheck = std::abs(lambdaRapidity) < 0.5;
          const bool kinematicLambdaCheck = pTLambdaCheck && rapidityLambdaCheck;

          ////////////////////////////////////////////
          // Ring observable: Leading particle proxy
          // Only computed when a valid leading particle exists (pT > minLeadParticlePt)
          ////////////////////////////////////////////
          float ringObservableLeadP = 0.;
          float deltaPhiLeadP = 0.;
          float deltaThetaLeadP = 0.;
          float cosDeltaThetaLeadP = 0.;
          if (hasValidLeadingP) {
            XYZVector crossLeadP = leadPUnitVec.Cross(lambdaLike3Vec);
            ringObservableLeadP = protonLikeStarUnit3Vec.Dot(crossLeadP) / crossLeadP.R();
            // Adding the prefactor related to the CP-violating decay (decay constants have different signs)
            if (!fakePolSwitches.forcePolSignQA)
              ringObservableLeadP *= (isLambda) ? PolPrefactorLambda : PolPrefactorAntiLambda;
            else
              ringObservableLeadP *= (isLambda) ? PolPrefactorLambda : -1.0 * PolPrefactorAntiLambda;
            // Angular variables
            deltaPhiLeadP = wrapToPiFast(v0phi - leadPPhi); // Wrapped to [-PI, PI), for convenience

            cosDeltaThetaLeadP = leadPUnitVec.Dot(lambdaLikeUnit3Vec); // Uses the pre-calculated unit vectors to avoid recomputation
            deltaThetaLeadP = std::acos(cosDeltaThetaLeadP);           // 3D angular separation. Same as ROOT::Math::VectorUtil::Angle(leadPUnitVec, lambdaLike3Vec);
          }

          //////////////////////////////////////////
          // Ring observable: Leading jet proxy
          // Only computed when a leading jet exists in this collision.
          //////////////////////////////////////////
          float ringObservable = 0.;
          float deltaPhiJet = 0.;
          float deltaEtaJet = 0.;
          float deltaThetaJet = 0.;
          float cosDeltaThetaJet = 0.;
          float Jz = 0.;
          float ringObservableOverJetZ = 0.;
          if (hasValidLeadingJet) {
            XYZVector cross = leadingJetUnitVec.Cross(lambdaLike3Vec);
            ringObservable = protonLikeStarUnit3Vec.Dot(cross) / cross.R();
            // Adding prefactor
            if (!fakePolSwitches.forcePolSignQA)
              ringObservable *= (isLambda) ? PolPrefactorLambda : PolPrefactorAntiLambda;
            else
              ringObservable *= (isLambda) ? PolPrefactorLambda : -1.0 * PolPrefactorAntiLambda;
            // Angular variables
            deltaPhiJet = wrapToPiFast(v0phi - leadingJetPhi);
            deltaEtaJet = v0eta - leadingJetEta;

            cosDeltaThetaJet = leadingJetUnitVec.Dot(lambdaLikeUnit3Vec);
            deltaThetaJet = std::acos(cosDeltaThetaJet);

            // Testing an invariance -- <R>/\hat{t}_z -- Possible source of an artificial (trivial) sign flip in the observable:
            Jz = leadingJetUnitVec.Z();
            if (std::abs(Jz) > 1e-4)
              ringObservableOverJetZ = ringObservable / Jz;
            else
              ringObservableOverJetZ = 0.0; // A simple guard. May not be the best, but works

            // // A second projection schema, where e_x = normalize(t_hat cross p_Lambda), using t_hat instead of z_hat (different from phi*):
            // // (this decomposes in an orthogonal basis related to the jet coordinates)
            // XYZVector ez = lambdaLikeUnit3Vec;
            // XYZVector ex = leadingJetUnitVec.Cross(lambdaLike3Vec);
            // XYZVector ey = ez.Cross(ex);

            // ringObservableExProjection = ringObservable;
            // ringObservableEyProjection =
            // // Assuming that energy can get inside the average in:
            // // P_Lambda \cdot p_Lambda = E_Lambda/m_Lambda * P_Lambda^* \cdot p_Lambda = <E_Lambda/m_Lambda * p_D^*> * p_Lambda
            // ringObservableEzProjection = cosFakePol * lambdaLike4Vec.E()/v0LambdaLikeMass;
          }

          //////////////////////////////////////////
          // Ring observable: Subleading jet proxy
          // Only computed when a subleading jet exists in this collision.
          //////////////////////////////////////////
          float ringObservable2ndJet = 0.;
          float deltaPhi2ndJet = 0.;
          float deltaTheta2ndJet = 0.;
          float cosDeltaTheta2ndJet = 0.;
          if (hasValidSubJet) {
            XYZVector cross2ndJet = subJetUnitVec.Cross(lambdaLike3Vec);
            ringObservable2ndJet = protonLikeStarUnit3Vec.Dot(cross2ndJet) / cross2ndJet.R();
            // Adding prefactor
            if (!fakePolSwitches.forcePolSignQA)
              ringObservable2ndJet *= (isLambda) ? PolPrefactorLambda : PolPrefactorAntiLambda;
            else
              ringObservable2ndJet *= (isLambda) ? PolPrefactorLambda : -1.0 * PolPrefactorAntiLambda;
            // Angular variables
            deltaPhi2ndJet = wrapToPiFast(v0phi - subleadingJetPhi);
            cosDeltaTheta2ndJet = subJetUnitVec.Dot(lambdaLikeUnit3Vec);
            deltaTheta2ndJet = std::acos(cosDeltaTheta2ndJet);
          }

          // Calculating polarization observables (in the Lambda frame, because that is easier -- does not require boosts):
          // To be precise, not actually the polarization, but a part of the summand in P^*_Lambda = (3/\alpha_Lambda) * <p^*_{proton}>
          float polStarX = 0, polStarY = 0, polStarZ = 0; // Dummy initialization: avoid warnings in compile time
          if (isLambda) {                                 // Notice there is no need to check analyseLambda again due to previous checks.
            polStarX = PolPrefactorLambda * protonLikeStarUnit3Vec.X();
            polStarY = PolPrefactorLambda * protonLikeStarUnit3Vec.Y();
            polStarZ = PolPrefactorLambda * protonLikeStarUnit3Vec.Z();
          } else {
            polStarX = PolPrefactorAntiLambda * protonLikeStarUnit3Vec.X();
            polStarY = PolPrefactorAntiLambda * protonLikeStarUnit3Vec.Y();
            polStarZ = PolPrefactorAntiLambda * protonLikeStarUnit3Vec.Z();
          }

          float v0phiToFillHists = wrapToPiFast(v0phi); // A short wrap to reuse some predefined axes

          // Fill ring histograms: (1D, lambda 2D correlations and jet 2D correlations):
          if (hasValidLeadingP) {
            if (familySwitches.doFamilyRing) {
              RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "Ring") // Notice the usage of macros! If you change the variable names, this WILL break the code!
                                                                        // No, there should NOT be any ";" here! Read the macro definition for an explanation
            }
            histos.fill(HIST("IntegratedCuts/pRingCutsLeadingP"), 0, ringObservableLeadP); // First bin of comparison
            histos.fill(HIST("IntegratedCuts/hCountCutsLeadingP"), 0);

            // Filling checks that rely on Eta>0 or Eta<0 checks for V0 and LeadingP eta:
            if (familySwitches.doFamilyRing) {
              RING_OBSERVABLE_LEADP_ETA_SPLIT_FILL_LIST("Ring", leadPEtaPos, lambdaEtaPos);
            }
          }
          if (familySwitches.doFamilyRing) {
            POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "Ring")
          }

          // Binary search using the pre-fetched axes for delta method of error bar estimation:
          int binPt = 0; // Dummy declarations
          int binMass = 0;
          int binDTheta = 0;
          if (hasValidLeadingJet) {
            if (familySwitches.doFamilyRing) {
              RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "Ring")
            }
            histos.fill(HIST("IntegratedCuts/pRingCuts"), 0, ringObservable);
            histos.fill(HIST("IntegratedCuts/hCountCuts"), 0);
            histos.fill(HIST("IntegratedCuts/pRingVsNV0s"), nLambdaLikeV0s, ringObservable);
            histos.fill(HIST("hNV0sVsCentrality"), nLambdaLikeV0s, centrality);

            // Properly fetching values as they are needed:
            binPt = mAxisPt->FindBin(v0pt);
            binMass = mAxisMass->FindBin(v0LambdaLikeMass);
            binDTheta = mAxisDTheta->FindBin(deltaThetaJet);
            if (familySwitches.doFamilyRing) {
              trackRing.addV0(ringObservable, binPt, binMass, binDTheta);
            }

            if (qaSwitches.doFakePolDiagnosticsQA) {
              // Measuring the AEE differentially (azimuthal efficiency effect, which causes different V0 topologies to be enhanced/suppressed):
              if (isLambda)
                histos.fill(HIST("HelicityEfficiencyQA/p2dRing_LambdaMassVsPhiLambdaMinusPhiProtonStar"), v0LambdaLikeMass, deltaPhiLambdaProtonStar, ringObservable);
              else
                histos.fill(HIST("HelicityEfficiencyQA/p2dRing_AntiLambdaMassVsPhiLambdaMinusPhiProtonStar"), v0LambdaLikeMass, deltaPhiLambdaProtonStar, ringObservable);
              // AEE and HEE correlation:
              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsJet_CosThetaVsPhiStar"), cosFakePol, phiStar);
              histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignal_CosThetaVsPhiStar"), cosFakePol, phiStar, ringObservable);

              // AEE and DCA between daughters correlation:
              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsJet_PhiStarVsDCAdau"), phiStar, dcaDau);
              histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAdau"), phiStar, dcaDau, ringObservable);
              histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAdauVsEtaJet"), phiStar, dcaDau, leadingJetEta, ringObservable);
              histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAdauVsEtaLambda"), phiStar, dcaDau, v0eta, ringObservable);

              // DCA dau to PV correlation:
              // For proton-like daughter:
              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsJet_PhiStarVsDCAProLike"), phiStar, protonLikeDCADauToPV);
              histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAProLike"), phiStar, protonLikeDCADauToPV, ringObservable);
              histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAProLikeVsEtaJet"), phiStar, protonLikeDCADauToPV, leadingJetEta, ringObservable);
              histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAProLikeVsEtaLambda"), phiStar, protonLikeDCADauToPV, v0eta, ringObservable);
              // For pion-like daughter:
              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsJet_PhiStarVsDCAPiLike"), phiStar, pionLikeDCADauToPV);
              histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAPiLike"), phiStar, pionLikeDCADauToPV, ringObservable);
              histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAPiLikeVsEtaJet"), phiStar, pionLikeDCADauToPV, leadingJetEta, ringObservable);
              histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalJet_PhiStarVsDCAPiLikeVsEtaLambda"), phiStar, pionLikeDCADauToPV, v0eta, ringObservable);

              histos.fill(HIST("HelicityEfficiencyQA/pRingVsJetZcomponent"), Jz, ringObservable);
              histos.fill(HIST("HelicityEfficiencyQA/pRingOverJetZcomponent_VsJetEta"), leadingJetEta, ringObservableOverJetZ);
              histos.fill(HIST("HelicityEfficiencyQA/pRingOverJetZcomponent_VsCosThetaHEE"), cosFakePol, ringObservableOverJetZ);
              histos.fill(HIST("HelicityEfficiencyQA/pRingOverJetZcomponent_VsPhiStar"), phiStar, ringObservableOverJetZ);
              histos.fill(HIST("HelicityEfficiencyQA/pRingOverJetZcomponent_VsJetEtaVsCosThetaHEE"), leadingJetEta, cosFakePol, ringObservableOverJetZ);
              histos.fill(HIST("HelicityEfficiencyQA/pRingOverJetZcomponent_VsJetEtaVsPhiStar"), leadingJetEta, phiStar, ringObservableOverJetZ);
            } // end doFakePolDiagnosticsQA (leading-jet AEE/HEE)
          }
          if (hasValidSubJet) {
            if (familySwitches.doFamilyRing) {
              RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "Ring")
            }
            histos.fill(HIST("IntegratedCuts/pRingCutsSubLeadingJet"), 0, ringObservable2ndJet);
            histos.fill(HIST("IntegratedCuts/hCountCutsSubLeadingJet"), 0);
          }

          if (qaSwitches.doFakePolDiagnosticsQA) {
            // Filling eta dependence QAs of the result (both for V0 and jet proxy):
            // Defining shared binning which depend on the V0 only:
            const int etaLambdaBin = lambdaEtaPos ? 3 : 4;
            if (hasValidLeadingJet) {
              histos.fill(HIST("EtaStudy/pRingEtaCuts"), 0, ringObservable);
              histos.fill(HIST("EtaStudy/pRingEtaCuts"), etaLambdaBin, ringObservable);

              // Bin indices for this proxy:
              // Bin 0: all, 1/2: proxy #eta sign, 3/4: #Lambda #eta sign, 5-8: joint,
              // 9/10: |#eta_{proxy}| >= R, 11-14: strict joint.
              const int etaProxyBin = jetEtaPos ? 1 : 2;
              const int etaProxyLambdaBin = (jetEtaPos ? 5 : 7) + (lambdaEtaPos ? 0 : 1);
              const int etaProxyStrictBin = jetEtaPos ? 9 : 10;
              const int etaProxyStrictLambdaBin = (jetEtaPos ? 11 : 13) + (lambdaEtaPos ? 0 : 1);

              histos.fill(HIST("EtaStudy/pRingEtaCuts"), etaProxyBin, ringObservable);
              histos.fill(HIST("EtaStudy/pRingEtaCuts"), etaProxyLambdaBin, ringObservable);

              // HEE study (helicity efficiency effect):
              histos.fill(HIST("EtaStudy/hFakePolCounts"), cosFakePol, 0);
              histos.fill(HIST("EtaStudy/hFakePolCounts"), cosFakePol, etaLambdaBin);
              histos.fill(HIST("EtaStudy/hFakePolCounts"), cosFakePol, etaProxyBin);
              histos.fill(HIST("EtaStudy/hFakePolCounts"), cosFakePol, etaProxyLambdaBin);
              // Same for signal:
              histos.fill(HIST("EtaStudy/pFakePolSignalVsCosTheta"), cosFakePol, 0, ringObservable);
              histos.fill(HIST("EtaStudy/pFakePolSignalVsCosTheta"), cosFakePol, etaLambdaBin, ringObservable);
              histos.fill(HIST("EtaStudy/pFakePolSignalVsCosTheta"), cosFakePol, etaProxyBin, ringObservable);
              histos.fill(HIST("EtaStudy/pFakePolSignalVsCosTheta"), cosFakePol, etaProxyLambdaBin, ringObservable);
              // Counter and ring accumulators for AEE study:
              histos.fill(HIST("EtaStudy/hCountsVsPhiStar"), phiStar, 0);
              histos.fill(HIST("EtaStudy/hCountsVsPhiStar"), phiStar, etaLambdaBin);
              histos.fill(HIST("EtaStudy/hCountsVsPhiStar"), phiStar, etaProxyBin);
              histos.fill(HIST("EtaStudy/hCountsVsPhiStar"), phiStar, etaProxyLambdaBin);
              histos.fill(HIST("EtaStudy/pFakePolSignalvsPhiStar"), phiStar, 0, ringObservable);
              histos.fill(HIST("EtaStudy/pFakePolSignalvsPhiStar"), phiStar, etaLambdaBin, ringObservable);
              histos.fill(HIST("EtaStudy/pFakePolSignalvsPhiStar"), phiStar, etaProxyBin, ringObservable);
              histos.fill(HIST("EtaStudy/pFakePolSignalvsPhiStar"), phiStar, etaProxyLambdaBin, ringObservable);

              histos.fill(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"), deltaPhiLambdaProtonStar, 0, ringObservable);
              histos.fill(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"), deltaPhiLambdaProtonStar, etaLambdaBin, ringObservable);
              histos.fill(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"), deltaPhiLambdaProtonStar, etaProxyBin, ringObservable);
              histos.fill(HIST("EtaStudy/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"), deltaPhiLambdaProtonStar, etaProxyLambdaBin, ringObservable);

              // Extra correlations test:
              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaJet"), cosFakePol, deltaThetaJet);
              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsCosThetaVsPtForJets"), cosFakePol, v0pt);
              // Split by proxy #eta sign:
              if (jetEtaPos)
                histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaJetPosEta"), cosFakePol, deltaThetaJet);
              else
                histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaJetNegEta"), cosFakePol, deltaThetaJet);

              if (pTLambdaCheck) {
                histos.fill(HIST("EtaStudy/hFakePolCountsLambdaPtCut"), cosFakePol, 0);
                histos.fill(HIST("EtaStudy/hFakePolCountsLambdaPtCut"), cosFakePol, etaLambdaBin);
                histos.fill(HIST("EtaStudy/hFakePolCountsLambdaPtCut"), cosFakePol, etaProxyBin);
                histos.fill(HIST("EtaStudy/hFakePolCountsLambdaPtCut"), cosFakePol, etaProxyLambdaBin);
                if (rapidityLambdaCheck) { // Stricter check
                  histos.fill(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"), cosFakePol, 0);
                  histos.fill(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"), cosFakePol, etaLambdaBin);
                  histos.fill(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"), cosFakePol, etaProxyBin);
                  histos.fill(HIST("EtaStudy/hFakePolCountsLambdaPtYCuts"), cosFakePol, etaProxyLambdaBin);
                }
              }
              if (jetEtaStrict) { // |eta_{Jet}| >= R
                histos.fill(HIST("EtaStudy/pRingEtaCuts"), etaProxyStrictBin, ringObservable);
                histos.fill(HIST("EtaStudy/pRingEtaCuts"), etaProxyStrictLambdaBin, ringObservable);
              }
            }
            if (hasValidSubJet) {
              // Same bin scheme as the leading jet above:
              const int etaProxyBin = subJetEtaPos ? 1 : 2;
              const int etaProxyLambdaBin = (subJetEtaPos ? 5 : 7) + (lambdaEtaPos ? 0 : 1);
              const int etaProxyStrictBin = subJetEtaPos ? 9 : 10;
              const int etaProxyStrictLambdaBin = (subJetEtaPos ? 11 : 13) + (lambdaEtaPos ? 0 : 1);

              histos.fill(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"), 0, ringObservable2ndJet);
              histos.fill(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"), etaLambdaBin, ringObservable2ndJet);
              histos.fill(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"), etaProxyBin, ringObservable2ndJet);
              histos.fill(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"), etaProxyLambdaBin, ringObservable2ndJet);
              if (subJetEtaStrict) { // |eta_{SubJet}| >= R
                histos.fill(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"), etaProxyStrictBin, ringObservable2ndJet);
                histos.fill(HIST("EtaStudy/pRingEtaCutsSubLeadingJet"), etaProxyStrictLambdaBin, ringObservable2ndJet);
              }
            }
            if (hasValidLeadingP) {
              // Same bin scheme, without the strict variant (this axis stops at 9 bins):
              const int etaProxyBin = leadPEtaPos ? 1 : 2;
              const int etaProxyLambdaBin = (leadPEtaPos ? 5 : 7) + (lambdaEtaPos ? 0 : 1);

              histos.fill(HIST("EtaStudy/pRingEtaCutsLeadingP"), 0, ringObservableLeadP);
              histos.fill(HIST("EtaStudy/pRingEtaCutsLeadingP"), etaLambdaBin, ringObservableLeadP);
              histos.fill(HIST("EtaStudy/pRingEtaCutsLeadingP"), etaProxyBin, ringObservableLeadP);
              histos.fill(HIST("EtaStudy/pRingEtaCutsLeadingP"), etaProxyLambdaBin, ringObservableLeadP);
              histos.fill(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"), 0, v0InMassPeak, ringObservableLeadP);
              histos.fill(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"), etaLambdaBin, v0InMassPeak, ringObservableLeadP);
              histos.fill(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"), etaProxyBin, v0InMassPeak, ringObservableLeadP);
              histos.fill(HIST("EtaStudy/pRingEtaCutsLeadingP_MassSignalVsBackground"), etaProxyLambdaBin, v0InMassPeak, ringObservableLeadP);

              histos.fill(HIST("EtaStudy/hFakePolCountsLeadP"), cosFakePol, 0);
              histos.fill(HIST("EtaStudy/hFakePolCountsLeadP"), cosFakePol, etaLambdaBin);
              histos.fill(HIST("EtaStudy/hFakePolCountsLeadP"), cosFakePol, etaProxyBin);
              histos.fill(HIST("EtaStudy/hFakePolCountsLeadP"), cosFakePol, etaProxyLambdaBin);

              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsCosThetaVsPtForLeadP"), cosFakePol, v0pt); // Understanding the population of events that has a leading particle (even though this does not need one to be calculated!)
              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaLeadP"), cosFakePol, deltaThetaLeadP);
              if (leadPEtaPos)
                histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaLeadPPosEta"), cosFakePol, deltaThetaLeadP);
              else
                histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaLeadPNegEta"), cosFakePol, deltaThetaLeadP);
            }
          } // end doFakePolDiagnosticsQA (eta-dependence block)

          // Extra kinematic criteria for Lambda candidates (removes polarization background):
          if (kinematicLambdaCheck) {
            if (hasValidLeadingP) {
              if (familySwitches.doFamilyRingKinematicCuts) {
                RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
              }
              histos.fill(HIST("IntegratedCuts/pRingCutsLeadingP"), 1, ringObservableLeadP);
              histos.fill(HIST("IntegratedCuts/hCountCutsLeadingP"), 1);

              // Filling checks that rely on Eta>0 or Eta<0 checks for V0 and LeadingP eta:
              if (familySwitches.doFamilyRingKinematicCuts) {
                RING_OBSERVABLE_LEADP_ETA_SPLIT_FILL_LIST("RingKinematicCuts", leadPEtaPos, lambdaEtaPos);
              }
            }
            if (familySwitches.doFamilyRingKinematicCuts) {
              POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
            }
            if (hasValidLeadingJet) {
              if (familySwitches.doFamilyRingKinematicCuts) {
                RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
              }
              histos.fill(HIST("IntegratedCuts/pRingCuts"), 1, ringObservable);
              histos.fill(HIST("IntegratedCuts/hCountCuts"), 1);
              if (familySwitches.doFamilyRingKinematicCuts) {
                trackRingKinCuts.addV0(ringObservable, binPt, binMass, binDTheta);
              }
            }
            if (hasValidSubJet) {
              if (familySwitches.doFamilyRingKinematicCuts) {
                RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
              }
              histos.fill(HIST("IntegratedCuts/pRingCutsSubLeadingJet"), 1, ringObservable2ndJet);
              histos.fill(HIST("IntegratedCuts/hCountCutsSubLeadingJet"), 1);
            }
          }

          // Extra selection criteria on jet candidates:
          // (redundant for jets with R=0.4, but for jets with R<0.4 the leading jet may be farther in eta)
          if (kinematicJetCheck) { // Already includes hasValidLeadingJet in the bool! (no need to check again)
            if (familySwitches.doFamilyJetKinematicCuts) {
              RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
            }
            histos.fill(HIST("IntegratedCuts/pRingCuts"), 2, ringObservable);
            histos.fill(HIST("IntegratedCuts/hCountCuts"), 2);
            if (familySwitches.doFamilyJetKinematicCuts) {
              POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
              trackJetKinCuts.addV0(ringObservable, binPt, binMass, binDTheta);
            }
          }

          // Extra selection criteria on both Lambda and jet candidates:
          if (kinematicLambdaCheck && kinematicJetCheck) {
            if (familySwitches.doFamilyJetAndLambdaKinematicCuts) {
              RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
            }
            histos.fill(HIST("IntegratedCuts/pRingCuts"), 3, ringObservable);
            histos.fill(HIST("IntegratedCuts/hCountCuts"), 3);
            if (familySwitches.doFamilyJetAndLambdaKinematicCuts) {
              POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
              trackJetLambdaKinCuts.addV0(ringObservable, binPt, binMass, binDTheta);
            }
          }

          // Same variations for the leading particle and for the subleading jet:
          // (kinematicLeadPCheck already encodes hasValidLeadingP, so no extra gate needed here)
          if (kinematicLeadPCheck) {
            if (familySwitches.doFamilyJetKinematicCuts) {
              RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
            }
            histos.fill(HIST("IntegratedCuts/pRingCutsLeadingP"), 2, ringObservableLeadP);
            histos.fill(HIST("IntegratedCuts/hCountCutsLeadingP"), 2);

            // Filling checks that rely on Eta>0 or Eta<0 checks for V0 and LeadingP eta:
            if (familySwitches.doFamilyJetKinematicCuts) {
              RING_OBSERVABLE_LEADP_ETA_SPLIT_FILL_LIST("JetKinematicCuts", leadPEtaPos, lambdaEtaPos);
            }
          }
          if (kinematic2ndJetCheck) {
            if (familySwitches.doFamilyJetKinematicCuts) {
              RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
            }
            histos.fill(HIST("IntegratedCuts/pRingCutsSubLeadingJet"), 2, ringObservable2ndJet);
            histos.fill(HIST("IntegratedCuts/hCountCutsSubLeadingJet"), 2);
          }
          if (kinematicLambdaCheck && kinematicLeadPCheck) {
            if (familySwitches.doFamilyJetAndLambdaKinematicCuts) {
              RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
            }
            histos.fill(HIST("IntegratedCuts/pRingCutsLeadingP"), 3, ringObservableLeadP);
            histos.fill(HIST("IntegratedCuts/hCountCutsLeadingP"), 3);

            // Filling checks that rely on Eta>0 or Eta<0 checks for V0 and LeadingP eta:
            if (familySwitches.doFamilyJetAndLambdaKinematicCuts) {
              RING_OBSERVABLE_LEADP_ETA_SPLIT_FILL_LIST("JetAndLambdaKinematicCuts", leadPEtaPos, lambdaEtaPos);
            }
          }
          if (kinematicLambdaCheck && kinematic2ndJetCheck) {
            if (familySwitches.doFamilyJetAndLambdaKinematicCuts) {
              RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
            }
            histos.fill(HIST("IntegratedCuts/pRingCutsSubLeadingJet"), 3, ringObservable2ndJet);
            histos.fill(HIST("IntegratedCuts/hCountCutsSubLeadingJet"), 3);
          }
        } // end v0s loop

        // Flush trackers to the actual O2 histograms (via macros, so that O2 compiles properly):
        if (familySwitches.doFamilyRing) {
          FLUSH_DELTA_TRACKER("Ring", trackRing, mAxisPt, mAxisMass, mAxisDTheta)
        }
        if (familySwitches.doFamilyRingKinematicCuts) {
          FLUSH_DELTA_TRACKER("RingKinematicCuts", trackRingKinCuts, mAxisPt, mAxisMass, mAxisDTheta)
        }
        if (familySwitches.doFamilyJetKinematicCuts) {
          FLUSH_DELTA_TRACKER("JetKinematicCuts", trackJetKinCuts, mAxisPt, mAxisMass, mAxisDTheta)
        }
        if (familySwitches.doFamilyJetAndLambdaKinematicCuts) {
          FLUSH_DELTA_TRACKER("JetAndLambdaKinematicCuts", trackJetLambdaKinCuts, mAxisPt, mAxisMass, mAxisDTheta)
        }
      } // end collisions
    } // end of resampling loop for forceRandJet and forceDatalikeJet
  }

  PROCESS_SWITCH(lambdajetpolarizationionsderived, processPolarizationData, "Process derived data in Run 3 Data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdajetpolarizationionsderived>(cfgc)};
}

// Avoid macro leakage!
#undef APPLY_HISTO_FILL
