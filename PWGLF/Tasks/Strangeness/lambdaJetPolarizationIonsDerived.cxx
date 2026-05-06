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
/// \file lambdajetpolarizationionsderived.cxx
/// \brief Lambda and antiLambda polarization analysis task using derived data
///
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

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/VectorUtil.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TProfile.h>
#include <TRandom3.h> // For perpendicular jet direction QAs

#include <cmath>
#include <optional>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using ROOT::Math::PtEtaPhiMVector;
using ROOT::Math::XYZVector;
// using namespace o2::aod::lambdajetpol; // Used it explicitly along the code for clarity

// Declaring constants:
constexpr double protonMass = o2::constants::physics::MassProton; // Assumes particle identification for daughter is perfect
constexpr double lambdaWeakDecayConstant = 0.749;                 // DPG 2025 update
constexpr double antiLambdaWeakDecayConstant = -0.758;            // DPG 2025 update
constexpr double polPrefactorLambda = 3.0 / lambdaWeakDecayConstant;
constexpr double polPrefactorAntiLambda = 3.0 / antiLambdaWeakDecayConstant;

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
  /* Additional plots for instant gratification - 1D Profiles */                                                           \
  X(FOLDER "/pRingObservableDeltaPhi", deltaPhiJet, ringObservable)                                                        \
  X(FOLDER "/pRingObservablePhiJet", leadingJetPhi, ringObservable)                                                        \
  X(FOLDER "/pRingObservablePhiLambda", v0phi, ringObservable)                                                             \
  X(FOLDER "/pRingObservableDeltaTheta", deltaThetaJet, ringObservable)                                                    \
  X(FOLDER "/EtaDependence/pRingObservableEtaLambda", v0eta, ringObservable)                                               \
  X(FOLDER "/EtaDependence/pRingObservableEtaJet", leadingJetEta, ringObservable)                                          \
  X(FOLDER "/pRingObservableIntegrated", 0., ringObservable)                                                               \
  X(FOLDER "/pRingObservableLambdaPt", v0pt, ringObservable)                                                               \
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
  X(FOLDER "/pRingVsCentrality", centrality, ringObservable)
// (TODO: add counters for regular TH2Ds about centrality)

// For leading particle
#define RING_OBSERVABLE_LEADP_FILL_LIST(X, FOLDER)                                  \
  X(FOLDER "/QA/hDeltaPhiLeadP", deltaPhiLeadP)                                     \
  X(FOLDER "/QA/hDeltaThetaLeadP", deltaThetaLeadP)                                 \
  X(FOLDER "/QA/hPtLeadP", leadPPt)                                                 \
  X(FOLDER "/QA/hCosDeltaThetaLeadP", cosDeltaThetaLeadP)                           \
  X(FOLDER "/pRingObservableLeadPDeltaPhi", deltaPhiLeadP, ringObservableLeadP)     \
  X(FOLDER "/pRingObservableLeadPDeltaTheta", deltaThetaLeadP, ringObservableLeadP) \
  X(FOLDER "/EtaDependence/pRingObservableEtaLambdaLeadP", v0eta, ringObservableLeadP)            \
  X(FOLDER "/EtaDependence/pRingObservableEtaLeadP", leadPEta, ringObservableLeadP)               \
  X(FOLDER "/pRingObservableLeadPIntegrated", 0., ringObservableLeadP)              \
  X(FOLDER "/pRingObservableLeadPLambdaPt", v0pt, ringObservableLeadP)              \
  X(FOLDER "/ProxyPtDependence/pRingVsPtLeadP", leadPPt, ringObservableLeadP)                      \
  X(FOLDER "/ProxyPtDependence/pRingVsPtLeadPVsEtaLeadP", leadPPt, leadPEta, ringObservableLeadP)  \
  X(FOLDER "/ProxyPtDependence/pRingVsPtLeadPVsEtaV0", leadPPt, v0eta, ringObservableLeadP)        \
  X(FOLDER "/EtaDependence/p2dRingObservableEtaLambdaVsEtaLeadP", v0eta, leadPEta, ringObservableLeadP) \
  X(FOLDER "/EtaDependence/h2dCounterEtaLambdaVsEtaLeadP", v0eta, leadPEta)

// A macro that encapsulates all eta checks for leading particle and V0s, along with the fills
// Parameters:
//   FOLDER       -- histogram folder string (compile-time literal)
//   LEADP_IS_POS -- bool: leadPEtaPos
//   V0_IS_POS    -- bool: lambdaEtaPos
#define RING_OBSERVABLE_LEADP_ETA_SPLIT_FILL_LIST(FOLDER, LEADP_IS_POS, V0_IS_POS)       \
  do {                                                                                   \
    if (LEADP_IS_POS) {                                                                  \
      /* leadP marginal: positive side */                                                \
      APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_PosEtaLeadP", leadPPt, ringObservableLeadP) \
      if (V0_IS_POS) {                                                                   \
        /* V0 marginal: positive side */                                                 \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_PosEtaV0", leadPPt, ringObservableLeadP) \
        /* Joint: (+leadP, +V0) */                                                       \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_PosEtaLeadP_PosEtaV0", leadPPt, ringObservableLeadP) \
      } else {                                                                           \
        /* V0 marginal: negative side */                                                 \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_NegEtaV0", leadPPt, ringObservableLeadP) \
        /* Joint: (+leadP, -V0) */                                                       \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_PosEtaLeadP_NegEtaV0", leadPPt, ringObservableLeadP) \
      }                                                                                  \
    } else {                                                                             \
      /* leadP marginal: negative side */                                                \
      APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_NegEtaLeadP", leadPPt, ringObservableLeadP) \
      if (V0_IS_POS) {                                                                   \
        /* V0 marginal: positive side */                                                 \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_PosEtaV0", leadPPt, ringObservableLeadP) \
        /* Joint: (-leadP, +V0) */                                                       \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_NegEtaLeadP_PosEtaV0", leadPPt, ringObservableLeadP) \
      } else {                                                                           \
        /* V0 marginal: negative side */                                                 \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_NegEtaV0", leadPPt, ringObservableLeadP) \
        /* Joint: (-leadP, -V0) */                                                       \
        APPLY_HISTO_FILL(FOLDER "/ProxyPtDependence/pRingVsPtLeadP_NegEtaLeadP_NegEtaV0", leadPPt, ringObservableLeadP) \
      }                                                                                  \
    }                                                                                    \
  } while (0)

// For subleading jet:
#define RING_OBSERVABLE_2NDJET_FILL_LIST(X, FOLDER)                                                   \
  X(FOLDER "/QA/hDeltaPhi2ndJet", deltaPhi2ndJet)                                                     \
  X(FOLDER "/QA/hDeltaTheta2ndJet", deltaTheta2ndJet)                                                 \
  X(FOLDER "/QA/hCosDeltaTheta2ndJet", cosDeltaTheta2ndJet)                                           \
  X(FOLDER "/QA/hPt2ndJet", subleadingJetPt)                                                          \
  X(FOLDER "/pRingObservable2ndJetDeltaPhi", deltaPhi2ndJet, ringObservable2ndJet)                    \
  X(FOLDER "/pRingObservable2ndJetDeltaTheta", deltaTheta2ndJet, ringObservable2ndJet)                \
  X(FOLDER "/EtaDependence/pRingObservableEtaLambda2ndJet", v0eta, ringObservable2ndJet)                            \
  X(FOLDER "/EtaDependence/pRingObservableEta2ndJet", subleadingJetEta, ringObservable2ndJet)                       \
  X(FOLDER "/pRingObservable2ndJetIntegrated", 0., ringObservable2ndJet)                              \
  X(FOLDER "/pRingObservable2ndJetLambdaPt", v0pt, ringObservable2ndJet)                              \
  X(FOLDER "/ProxyPtDependence/pRingVsPt2ndJet", subleadingJetPt, ringObservable2ndJet)                              \
  X(FOLDER "/ProxyPtDependence/pRingVsPt2ndJetVsEta2ndJet", subleadingJetPt, subleadingJetEta, ringObservable2ndJet) \
  X(FOLDER "/ProxyPtDependence/pRingVsPt2ndJetVsEtaV0", subleadingJetPt, v0eta, ringObservable2ndJet)                \
  X(FOLDER "/EtaDependence/p2dRingObservableEtaLambdaVsEta2ndJet", v0eta, subleadingJetEta, ringObservable2ndJet) \
  X(FOLDER "/EtaDependence/h2dCounterEtaLambdaVsEta2ndJet", v0eta, subleadingJetEta) \

#define POLARIZATION_PROFILE_FILL_LIST(X, FOLDER)                          \
  /* =============================== */                                    \
  /* 1D TProfiles vs v0phi */                                              \
  /* =============================== */                                    \
  X(FOLDER "/QA/pPxStarPhi", v0phiToFillHists, PolStarX)                   \
  X(FOLDER "/QA/pPyStarPhi", v0phiToFillHists, PolStarY)                   \
  X(FOLDER "/QA/pPzStarPhi", v0phiToFillHists, PolStarZ)                   \
  /* =============================== */                                    \
  /* 1D TProfiles vs DeltaPhi_jet */                                       \
  /* =============================== */                                    \
  X(FOLDER "/QA/pPxStarDeltaPhi", deltaPhiJet, PolStarX)                   \
  X(FOLDER "/QA/pPyStarDeltaPhi", deltaPhiJet, PolStarY)                   \
  X(FOLDER "/QA/pPzStarDeltaPhi", deltaPhiJet, PolStarZ)                   \
  /* =============================== */                                    \
  /* 2D TProfiles vs DeltaPhi_jet and Lambda pT */                         \
  /* =============================== */                                    \
  X(FOLDER "/QA/p2dPxStarDeltaPhiVsLambdaPt", deltaPhiJet, v0pt, PolStarX) \
  X(FOLDER "/QA/p2dPyStarDeltaPhiVsLambdaPt", deltaPhiJet, v0pt, PolStarY) \
  X(FOLDER "/QA/p2dPzStarDeltaPhiVsLambdaPt", deltaPhiJet, v0pt, PolStarZ)

// Apply the macros (notice I had to include the semicolon (";") after the function, so you don't need to
// write that when calling this APPLY_HISTO_FILL. The code will look weird, but without this the compiler
// would not know to end each statement with a semicolon):
#define APPLY_HISTO_FILL(NAME, ...) histos.fill(HIST(NAME), __VA_ARGS__);


// Delta Method Fill Lists
#define DELTA_INTEGRATED_FILL_LIST(X, FOLDER, r, n) \
  X(FOLDER "/DeltaMethod/hIntegrated", 0.5, r) \
  X(FOLDER "/DeltaMethod/hIntegrated", 1.5, (double)(n)) \
  X(FOLDER "/DeltaMethod/hIntegrated", 2.5, (r)*(r)) \
  X(FOLDER "/DeltaMethod/hIntegrated", 3.5, (double)((n)*(n))) \
  X(FOLDER "/DeltaMethod/hIntegrated", 4.5, (r)*(n))

#define DELTA_2D_FILL_LIST(X, FOLDER, HIST_NAME, center, r, n) \
  X(FOLDER "/DeltaMethod/" HIST_NAME, center, 0.5, r) \
  X(FOLDER "/DeltaMethod/" HIST_NAME, center, 1.5, (double)(n)) \
  X(FOLDER "/DeltaMethod/" HIST_NAME, center, 2.5, (r)*(r)) \
  X(FOLDER "/DeltaMethod/" HIST_NAME, center, 3.5, (double)((n)*(n))) \
  X(FOLDER "/DeltaMethod/" HIST_NAME, center, 4.5, (r)*(n))

// Master flush macro to dump an event tracker into the histograms:
#define FLUSH_DELTA_TRACKER(FOLDER, TRACKER, AXIS_PT, AXIS_MASS, AXIS_DTHETA) \
  if ((TRACKER).n_int > 0) { \
    DELTA_INTEGRATED_FILL_LIST(APPLY_HISTO_FILL, FOLDER, (TRACKER).r_int, (TRACKER).n_int) \
  } \
  for (const auto& [bin, r_val] : (TRACKER).r_pt) { \
    int n_val = (TRACKER).n_pt.at(bin); \
    double center = (AXIS_PT)->GetBinCenter(bin); \
    DELTA_2D_FILL_LIST(APPLY_HISTO_FILL, FOLDER, "h2dLambdaPtVsDeltaComp", center, r_val, n_val) \
  } \
  for (const auto& [bin, r_val] : (TRACKER).r_mass) { \
    int n_val = (TRACKER).n_mass.at(bin); \
    double center = (AXIS_MASS)->GetBinCenter(bin); \
    DELTA_2D_FILL_LIST(APPLY_HISTO_FILL, FOLDER, "h2dMassVsDeltaComp", center, r_val, n_val) \
  } \
  for (const auto& [bin, r_val] : (TRACKER).r_dtheta) { \
    int n_val = (TRACKER).n_dtheta.at(bin); \
    double center = (AXIS_DTHETA)->GetBinCenter(bin); \
    DELTA_2D_FILL_LIST(APPLY_HISTO_FILL, FOLDER, "h2dDeltaThetaVsDeltaComp", center, r_val, n_val) \
  }

struct lambdajetpolarizationionsderived {

  // Define histogram registries:
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Master analysis switches
  Configurable<bool> analyseLambda{"analyseLambda", true, "process Lambda-like candidates"};
  Configurable<bool> analyseAntiLambda{"analyseAntiLambda", false, "process AntiLambda-like candidates"};
  Configurable<bool> analyseMagField{"analyseMagField", true, "analyse efficiency effects wrt magnetic field"}; // Older DerivedData lacks magField. "if constexpr (requires { collision.magField(); })" would only see the current header definition, so need a flag for retro-comp.
  Configurable<bool> doPPAnalysis{"doPPAnalysis", false, "if in pp, set to true. Default is HI"};

  // Centrality:
  Configurable<int> centralityEstimator{"centralityEstimator", kCentFT0M, "Run 3 centrality estimator (0:CentFT0C, 1:CentFT0M, 2:CentFV0A)"}; // Default is FT0M

  // QAs that purposefully break the analysis
  // -- All of these tests should give us zero signal if the source is truly Lambda Polarization from vortices
  Configurable<bool> forcePolSignQA{"forcePolSignQA", false, "force antiLambda decay constant to be positive: should kill all the signal, if any. For QA"};
  Configurable<bool> forcePerpToJet{"forcePerpToJet", false, "force jet direction to be perpendicular to jet estimator. For QA"};
  Configurable<bool> forceJetDirectionSmudge{"forceJetDirectionSmudge", false, "fluctuate jet direction by 10% of R around original axis. For QA (tests sensibility)"};
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
    ConfigurableAxis axisLambdaMass{"axisLambdaMass", {450, 1.08f, 1.15f}, "Lambda mass in GeV/c"}; // Default is {200, 1.101f, 1.131f}

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
    ConfigurableAxis axisDeltaEtaCoarse{"axisDeltaEtaCoarse", {40, -1.8f, 1.8f}, "#Delta#eta coarse axis"};
    ConfigurableAxis axisDeltaTheta{"axisDeltaTheta", {40, 0, constants::math::PI}, "#Delta #theta_{jet}"};
    ConfigurableAxis axisCosTheta{"axisCosTheta", {50, -1, 1}, "cos(#theta)"};
    ConfigurableAxis axisPhi{"axisPhi", {40, 0., constants::math::TwoPI}, "#varphi"};
    ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {40, -constants::math::PI, constants::math::PI}, "#Delta #phi_{jet}"};

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

    ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Centrality"};

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
    double r_int = 0.0; int n_int = 0;
    std::map<int, double> r_pt, r_mass, r_dtheta;
    std::map<int, int> n_pt, n_mass, n_dtheta;

    void addV0(double ringObs, int binPt, int binMass, int binDTheta) {
      r_int += ringObs;               n_int += 1;
      r_pt[binPt] += ringObs;         n_pt[binPt] += 1;
      r_mass[binMass] += ringObs;     n_mass[binMass] += 1;
      r_dtheta[binDTheta] += ringObs; n_dtheta[binDTheta] += 1;
    }
  };

  // Axis pointers for Delta Method binning (fetched once in init, declared once here)
  TAxis* mAxisPt = nullptr;
  TAxis* mAxisMass = nullptr;
  TAxis* mAxisDTheta = nullptr;

  void init(InitContext const&)
  {
    // Ring observable histograms:
    // Helper to register one full histogram family (kinematic cut variation of ring observable)
    auto addRingObservableFamily = [&](const std::string& folder) {
      // ===============================
      // QA histograms: angle and pT distributions
      // (No mass dependency -- useful to check kinematic sculpting from cuts)
      // ===============================
      histos.add((folder + "/QA/hDeltaPhi").c_str(), "#Delta#varphi_{jet};#Delta#varphi_{jet};Counts", kTH1D, {axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/QA/hDeltaPhiVsDeltaEta").c_str(), "#Delta#varphi_{jet};#Delta#varphi_{jet}; #eta_{#Lambda}-#eta_{Jet};Counts", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDeltaEtaCoarse});
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
      histos.add((folder + "/EtaDependence/pRingObservableEtaLambda").c_str(), "<#it{R}> vs #eta_{#Lambda};#eta_{#Lambda};<#it{R}>", kTProfile, {axisConfigurations.axisEtaCoarse});
      histos.add((folder + "/EtaDependence/pRingObservableEtaJet").c_str(), "<#it{R}> vs #eta_{Jet};#eta_{Jet};<#it{R}>", kTProfile, {axisConfigurations.axisEtaCoarse});
      
      histos.add((folder + "/EtaDependence/pRingObservableEtaLambda2ndJet").c_str(), "<#it{R}> vs #eta_{#Lambda} (SubJet);#eta_{#Lambda};<#it{R}>", kTProfile, {axisConfigurations.axisEtaCoarse});
      histos.add((folder + "/EtaDependence/pRingObservableEta2ndJet").c_str(), "<#it{R}> vs #eta_{SubJet};#eta_{SubJet};<#it{R}>", kTProfile, {axisConfigurations.axisEtaCoarse});
      
      histos.add((folder + "/EtaDependence/pRingObservableEtaLambdaLeadP").c_str(), "<#it{R}> vs #eta_{#Lambda} (LeadP);#eta_{#Lambda};<#it{R}>", kTProfile, {axisConfigurations.axisEtaCoarse});
      histos.add((folder + "/EtaDependence/pRingObservableEtaLeadP").c_str(), "<#it{R}> vs #eta_{LeadP};#eta_{LeadP};<#it{R}>", kTProfile, {axisConfigurations.axisEtaCoarse});
      // For the leading particle:
      histos.add((folder + "/pRingObservableLeadPDeltaPhi").c_str(), "<#it{R}> vs #Delta#varphi_{LeadP};#Delta#varphi_{LeadP};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/pRingObservableLeadPDeltaTheta").c_str(), "<#it{R}> vs #Delta#theta_{LeadP};#Delta#theta_{LeadP};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaTheta});
      histos.add((folder + "/pRingObservableLeadPIntegrated").c_str(), "Integrated <#it{R}> (LeadP); ;<#it{R}>", kTProfile, {{1, -0.5, 0.5}});
      histos.add((folder + "/pRingObservableLeadPLambdaPt").c_str(), "<#it{R}> vs #it{p}_{T}^{#Lambda} (LeadP);#it{p}_{T}^{#Lambda} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisPt});
      // For the second-to-leading jet:
      histos.add((folder + "/pRingObservable2ndJetDeltaPhi").c_str(), "<#it{R}> vs #Delta#varphi_{SubJet};#Delta#varphi_{SubJet};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaPhi});
      histos.add((folder + "/pRingObservable2ndJetDeltaTheta").c_str(), "<#it{R}> vs #Delta#theta_{SubJet};#Delta#theta_{SubJet};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaTheta});
      histos.add((folder + "/pRingObservable2ndJetIntegrated").c_str(), "Integrated <#it{R}> (SubJet); ;<#it{R}>", kTProfile, {{1, -0.5, 0.5}});
      histos.add((folder + "/pRingObservable2ndJetLambdaPt").c_str(), "<#it{R}> vs #it{p}_{T}^{#Lambda} (SubJet);#it{p}_{T}^{#Lambda} (GeV/c);<#it{R}>", kTProfile, {axisConfigurations.axisPt});
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
      histos.add((folder + "/EtaDependence/p2dRingObservableEtaLambdaVsEtaJet").c_str(), "<#it{R}> vs #eta_{#Lambda} vs #eta_{Jet};#eta_{#Lambda};#eta_{Jet};<#it{R}>", kTProfile2D, {axisConfigurations.axisEtaCoarse, axisConfigurations.axisEtaCoarse});
      histos.add((folder + "/EtaDependence/p2dRingObservableEtaLambdaVsEtaJet_FineBins").c_str(), "<#it{R}> vs #eta_{#Lambda} vs #eta_{Jet} (fine bins);#eta_{#Lambda};#eta_{Jet};<#it{R}>", kTProfile2D, {axisConfigurations.axisEta, axisConfigurations.axisEta});
      histos.add((folder + "/EtaDependence/p2dRingObservableEtaLambdaVsEtaLeadP").c_str(), "<#it{R}> vs #eta_{#Lambda} vs #eta_{LeadP};#eta_{#Lambda};#eta_{LeadP};<#it{R}>", kTProfile2D, {axisConfigurations.axisEtaCoarse, axisConfigurations.axisEtaCoarse});
      histos.add((folder + "/EtaDependence/p2dRingObservableEtaLambdaVsEta2ndJet").c_str(), "<#it{R}> vs #eta_{#Lambda} vs #eta_{SubJet};#eta_{#Lambda};#eta_{SubJet};<#it{R}>", kTProfile2D, {axisConfigurations.axisEtaCoarse, axisConfigurations.axisEtaCoarse});
        // Counters for these histograms, instead of only TProfile2Ds:
      histos.add((folder + "/EtaDependence/h2dCounterEtaLambdaVsEtaJet").c_str(), "Counts, #eta_{#Lambda} vs #eta_{Jet};#eta_{#Lambda};#eta_{Jet};Counts", kTH2D, {axisConfigurations.axisEtaCoarse, axisConfigurations.axisEtaCoarse});
      histos.add((folder + "/EtaDependence/h2dCounterEtaLambdaVsEtaJet_FineBins").c_str(), "Counts (fine bins), #eta_{#Lambda} vs #eta_{Jet};#eta_{#Lambda};#eta_{Jet};Counts", kTH2D, {axisConfigurations.axisEta, axisConfigurations.axisEta});
      histos.add((folder + "/EtaDependence/h2dCounterEtaLambdaVsEtaLeadP").c_str(), "Counts, #eta_{#Lambda} vs #eta_{LeadP};#eta_{#Lambda};#eta_{LeadP};Counts", kTH2D, {axisConfigurations.axisEtaCoarse, axisConfigurations.axisEtaCoarse});
      histos.add((folder + "/EtaDependence/h2dCounterEtaLambdaVsEta2ndJet").c_str(), "Counts, #eta_{#Lambda} vs #eta_{SubJet};#eta_{#Lambda};#eta_{SubJet};Counts", kTH2D, {axisConfigurations.axisEtaCoarse, axisConfigurations.axisEtaCoarse});
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
    addRingObservableFamily("Ring");
    addRingObservableFamily("RingKinematicCuts");
    addRingObservableFamily("JetKinematicCuts");
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


    // Integrated observable dependent on jet proxy #eta to unfold possible asymmetries in detector:
    histos.add("ProxyEta/pRingEtaCuts", "pRingEtaCuts; ;<#it{R}>", kTProfile, {{15, 0, 15}});
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(1,  "All #Lambda");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(2,  "#eta_{Jet} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(3,  "#eta_{Jet} < 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(4,  "#eta_{#Lambda} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(5,  "#eta_{#Lambda} < 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(6,  "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(7,  "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(8,  "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(9,  "#eta_{Jet} < 0, #eta_{#Lambda} < 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(10, "#eta_{Jet} > R");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(11, "#eta_{Jet} < -R");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(12, "#eta_{Jet} > R, #eta_{#Lambda} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(13, "#eta_{Jet} > R, #eta_{#Lambda} < 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(14, "#eta_{Jet} < -R, #eta_{#Lambda} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCuts"))->GetXaxis()->SetBinLabel(15, "#eta_{Jet} < -R, #eta_{#Lambda} < 0");

    histos.add("ProxyEta/pRingEtaCutsSubLeadingJet", "pRingEtaCutsSubLeadingJet; ;<#it{R}>", kTProfile, {{15, 0, 15}});
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(1,  "All #Lambda");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(2,  "#eta_{SubJet} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(3,  "#eta_{SubJet} < 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(4,  "#eta_{#Lambda} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(5,  "#eta_{#Lambda} < 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(6,  "#eta_{SubJet} #geq 0, #eta_{#Lambda} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(7,  "#eta_{SubJet} #geq 0, #eta_{#Lambda} < 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(8,  "#eta_{SubJet} < 0, #eta_{#Lambda} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(9,  "#eta_{SubJet} < 0, #eta_{#Lambda} < 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(10, "#eta_{SubJet} > R");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(11, "#eta_{SubJet} < -R");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(12, "#eta_{SubJet} > R, #eta_{#Lambda} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(13, "#eta_{SubJet} > R, #eta_{#Lambda} < 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(14, "#eta_{SubJet} < -R, #eta_{#Lambda} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(15, "#eta_{SubJet} < -R, #eta_{#Lambda} < 0");

    histos.add("ProxyEta/pRingEtaCutsLeadingP", "pRingEtaCutsLeadingP; ;<#it{R}>", kTProfile, {{9, 0, 9}});
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(2, "#eta_{LeadP} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(3, "#eta_{LeadP} < 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(6, "#eta_{LeadP} #geq 0, #eta_{#Lambda} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(7, "#eta_{LeadP} #geq 0, #eta_{#Lambda} < 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(8, "#eta_{LeadP} < 0, #eta_{#Lambda} #geq 0");
    histos.get<TProfile>(HIST("ProxyEta/pRingEtaCutsLeadingP"))->GetXaxis()->SetBinLabel(9, "#eta_{LeadP} < 0, #eta_{#Lambda} < 0");

    // Fake polarization signal QA
    // --> The "negative helicity problem", where topologies with a proton decaying opposite to the Lambda momentum are enhanced by
    // efficiency of reconstruction. The geometries where the proton moves in the same direction as the boost will have a very small
    // momentum pion, which is not as easily detected as the opposite case! This may insert a fake signal of polarization in the measurement!
    histos.add("HelicityEfficiencyQA/hFakePolCounts", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda};", kTH2D, {axisConfigurations.axisCosTheta, {9, 0, 9}});
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCounts"))->GetZaxis()->SetTitle("N_{V0s}");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCounts"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCounts"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCounts"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCounts"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCounts"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCounts"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCounts"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCounts"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCounts"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");

      // The same, but for actual signal instead of counts:
    histos.add("HelicityEfficiencyQA/pFakePolSignalVsCosTheta", "FakePolSignal; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda};", kTProfile2D, {axisConfigurations.axisCosTheta, {9, 0, 9}});
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"))->GetZaxis()->SetTitle("<#it{R}>");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");

    // Seeing the dependence between phi* = atan2(p_p_star \cdot (p_Lambda_hat \times (z_hat \cross p_Lambda)), p_p_star \cdot (z_hat \cross p_Lambda))
      // e_z = p_Lambda_hat; // e_x = normalize(z_hat cross p_Lambda); // e_y = e_z cross e_x;
      // phi_star = atan2(p_p_star dot e_y, p_p_star dot e_x);
    histos.add("HelicityEfficiencyQA/hFakePolCounts_CosThetaVsPhiStar", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; #phi^{*}", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaPhi});
    histos.add("HelicityEfficiencyQA/pFakePolSignal_CosThetaVsPhiStar", "FakePolSignal; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; #phi^{*}", kTProfile2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaPhi});
        // Specific counter for when we have leading jets (relates directly to pFakePolSignal_CosThetaVsPhiStar):
    histos.add("HelicityEfficiencyQA/hFakePolCountsJet_CosThetaVsPhiStar", "FakePolCounts - HasValidLeadJet OK; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; #phi^{*}", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaPhi});

      // Similar split, but for AEE instead of HEE:
    histos.add("HelicityEfficiencyQA/hCountsVsPhiStar", "FakePolCounts, AEE dependence; #phi^{*};", kTH2D, {axisConfigurations.axisDeltaPhi, {9, 0, 9}});
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"))->GetZaxis()->SetTitle("Counts");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");
        // For the ring observable as well:
    histos.add("HelicityEfficiencyQA/pFakePolSignalvsPhiStar", "FakePolSignal, AEE dependence; #phi^{*};", kTProfile2D, {axisConfigurations.axisDeltaPhi, {9, 0, 9}});
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"))->GetZaxis()->SetTitle("<#it{R}>");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");

      // For the phi_Lambda - phi_D^* dependency as well:
    histos.add("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar", "FakePolSignal, AEE dependence; #phi_{#Lambda} - #phi_{p}^{*};", kTProfile2D, {axisConfigurations.axisDeltaPhi, {9, 0, 9}});
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetZaxis()->SetTitle("<#it{R}>");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
    histos.get<TProfile2D>(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");

    // More about possible AEE dependencies (should see an invariance with JetEta and <R>/Jz):
      // TODO: think about these error bars: do they still make sense via regular TProfile's SEM error?
    histos.add("HelicityEfficiencyQA/pRingVsJetZcomponent", "<#it{R}> vs #hat{t}_{z}; #hat{t}_{z}; <#it{R}>", kTProfile, {{40, -1, 1}}); // Numerically stable and can also show the sign flip (essentially the <#it{R}> vs Eta Jet plot in another scale)
    histos.add("HelicityEfficiencyQA/pRingOverJetZcomponent_VsJetEta", "<#it{R}>/#hat{t}_{z}; #eta_{Jet}; <#it{R}>/#hat{t}_{z}", kTProfile, {axisConfigurations.axisEtaCoarse});
    histos.add("HelicityEfficiencyQA/pRingOverJetZcomponent_VsCosThetaHEE", "<#it{R}>/#hat{t}_{z} Vs cos(#theta) HEE; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; <#it{R}>/#hat{t}_{z}", kTProfile, {axisConfigurations.axisCosTheta});
    histos.add("HelicityEfficiencyQA/pRingOverJetZcomponent_VsPhiStar", "<#it{R}>/#hat{t}_{z} Vs #phi^{*}; #phi^{*} = atan2(#vec{p}^{*}_{p} #cdot [#hat{p}_{#Lambda} #times (#hat{z} #times #hat{p}_{#Lambda})] , #vec{p}^{*}_{p} #cdot (#hat{z} #times #hat{p}_{#Lambda})); <#it{R}>/#hat{t}_{z}", kTProfile, {axisConfigurations.axisDeltaPhi});
    histos.add("HelicityEfficiencyQA/pRingOverJetZcomponent_VsJetEtaVsCosThetaHEE", "<#it{R}>/#hat{t}_{z}; #eta_{Jet}; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda}; <#it{R}>/#hat{t}_{z}", kTProfile2D, {axisConfigurations.axisEtaCoarse, axisConfigurations.axisCosTheta});
    histos.add("HelicityEfficiencyQA/pRingOverJetZcomponent_VsJetEtaVsPhiStar", "<#it{R}>/#hat{t}_{z}; #eta_{Jet}; #phi^{*}; <#it{R}>/#hat{t}_{z}", kTProfile2D, {axisConfigurations.axisEtaCoarse, axisConfigurations.axisDeltaPhi});

    // Doing the same HEE study for leading particles:
    // (eta_{Jet} may be a bad estimator!)
    histos.add("HelicityEfficiencyQA/hFakePolCountsLeadP", "FakePolCounts; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda};", kTH2D, {axisConfigurations.axisCosTheta, {9, 0, 9}});
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"))->GetZaxis()->SetTitle("N_{V0s}");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(2, "#eta_{LeadP} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(3, "#eta_{LeadP} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(6, "#eta_{LeadP} #geq 0, #eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(7, "#eta_{LeadP} #geq 0, #eta_{#Lambda} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(8, "#eta_{LeadP} < 0, #eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"))->GetYaxis()->SetBinLabel(9, "#eta_{LeadP} < 0, #eta_{#Lambda} < 0");

    // Avoid fake signal by jets boosting the Lambda in its own direction, then modifying efficiency of reconstruction in a similar way:
    histos.add("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut", "FakePol,p_{T}^{#Lambda}#in[0.5,1.5]; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda};", kTH2D, {axisConfigurations.axisCosTheta, {9, 0, 9}});
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"))->GetZaxis()->SetTitle("N_{V0s}");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");

    // Even stricter cut (also demands rapidity cut stricter than jets, so may see different boosting):
    histos.add("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts", "FakePol,p_{T}^{#Lambda}#in[0.5,1.5],|y_{#Lambda}|<0.5; cos(#theta)=#hat{p}^{*}_{D} . #vec{p}_{#Lambda};", kTH2D, {axisConfigurations.axisCosTheta, {9, 0, 9}});
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"))->GetZaxis()->SetTitle("N_{V0s}");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(1, "All #Lambda");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(2, "#eta_{Jet} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(3, "#eta_{Jet} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(4, "#eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(5, "#eta_{#Lambda} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(6, "#eta_{Jet} #geq 0, #eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(7, "#eta_{Jet} #geq 0, #eta_{#Lambda} < 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(8, "#eta_{Jet} < 0, #eta_{#Lambda} #geq 0");
    histos.get<TH2>(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"))->GetYaxis()->SetBinLabel(9, "#eta_{Jet} < 0, #eta_{#Lambda} < 0");

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

    // Integrated observable for events with NLambda+NAntiLambda V0s per event
    // (an interesting measurement of correlation between <R> and Lambda-like V0s multiplicity. A proxy of covariance)
    // (calculated for leading jets only)
    histos.add("IntegratedCuts/pRingVsNV0s", "pRingVsNV0s; N_{#Lambda}+N_{#bar{#Lambda}};<#it{R}>", kTProfile, {{20, 0, 20}});

    // Leading Jet QA:
    histos.add("JetKinematicsQA/hLeadJetEta", "hLeadJetEta", kTH1D, {axisConfigurations.axisEta});
    histos.add("JetKinematicsQA/hSubLeadJetEta", "hSubLeadJetEta", kTH1D, {axisConfigurations.axisEta});
    histos.add("JetKinematicsQA/hLeadPEta", "hLeadPEta", kTH1D, {axisConfigurations.axisEta});

    // Counting the number of jets/proxies themselves (these count at most once per event) -- Similar to the FOLDER/QA/hPtJet counters:
    histos.add("JetKinematicsQA/hJetCounterPtJet", "hJetCounterPtJet; p_{T}^{Jet} (GeV/c)", kTH1D, {axisConfigurations.axisJetPt});
    histos.add("JetKinematicsQA/hJetCounterPtLeadP", "hJetCounterPtLeadP; p_{T}^{LeadP} (GeV/c)", kTH1D, {axisConfigurations.axisJetPt});
    histos.add("JetKinematicsQA/hJetCounterPt2ndJet", "hJetCounterPt2ndJet; p_{T}^{SubJet} (GeV/c)", kTH1D, {axisConfigurations.axisJetPt});

    // Fetch the X-axes from one of the families (since they all share the same ConfigurableAxis binning)
    mAxisPt = histos.get<TH2>(HIST("Ring/DeltaMethod/h2dLambdaPtVsDeltaComp"))->GetXaxis();
    mAxisMass = histos.get<TH2>(HIST("Ring/DeltaMethod/h2dMassVsDeltaComp"))->GetXaxis();
    mAxisDTheta = histos.get<TH2>(HIST("Ring/DeltaMethod/h2dDeltaThetaVsDeltaComp"))->GetXaxis();
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

  // Preslices for correct collisions association:
  // (TODO: test using custom grouping)
  Preslice<aod::RingJets> perColJets = o2::aod::lambdajetpol::ringCollisionId;
  Preslice<aod::RingLaV0s> perColV0s = o2::aod::lambdajetpol::ringCollisionId;
  Preslice<aod::RingLeadPs> perColLeadPs = o2::aod::lambdajetpol::ringCollisionId;
  void processPolarizationData(o2::aod::RingCollisions const& collisions, o2::aod::RingJets const& jets, o2::aod::RingLaV0s const& v0s,
                               o2::aod::RingLeadPs const& leadPs)
  {
    for (auto const& collision : collisions) {
      const auto collId = collision.globalIndex(); // The self-index accessor
      const double centrality = getCentrality(collision);
      
      // Fetch magnetic field only if DataModel is in the latest version:
      float magField = 1.f; // Used this dummy for backwards compatibility, under the reasonable assumption that the field points always in the same *direction* in the used runs
      // if (analyseMagField) 
      //   magField = collision.magField();

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
      if (!v0sInColl.size())
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
        // Table should contain exactly one entry per collision,
        // but we break immediately to be safe:
        leadPPt = lp.leadParticlePt();
        leadPEta = lp.leadParticleEta();
        leadPPhi = lp.leadParticlePhi();
        // Using dynamic columns to make code cleaner:
        leadPPx = lp.leadParticlePx();
        leadPPy = lp.leadParticlePy();
        leadPPz = lp.leadParticlePz();
        break;
      }
      // // Discard events with no leading particle (FastJet didn't even run in these cases!):
      // if (leadPPt < 0.)
      //   continue;

      // Apply minimum pT selection for the leading particle (not necessarily the same as in derived data builder. Can be a stricter cut!):
      const bool hasValidLeadingP = leadPPt > minLeadParticlePt;

      // Build leading particle unit vector, outside the V0 loop for performance.
      XYZVector leadPUnitVec(1., 0., 0.); // dummy (overwritten below when hasValidLeadingP)
      if (hasValidLeadingP) {
        histos.fill(HIST("JetKinematicsQA/hLeadPEta"), leadPEta);
        histos.fill(HIST("JetKinematicsQA/hJetCounterPtLeadP"), leadPPt);

        leadPUnitVec = XYZVector(leadPPx, leadPPy, leadPPz).Unit();
        // QA: same direction-smearing/perp logic as for the leading jet estimator.
        // The hLeadPEta histogram above intentionally uses the unmodified direction.
        if (forcePerpToJet) {
          XYZVector refVec(1., 0., 0.);
          if (std::abs(leadPUnitVec.Dot(refVec)) > 0.99)
            refVec = XYZVector(0., 1., 0.);
          XYZVector perpVec = leadPUnitVec.Cross(refVec).Unit();
          double randomAngle = randomGen.Uniform(0., o2::constants::math::TwoPI);
          leadPUnitVec = perpVec * std::cos(randomAngle) + leadPUnitVec.Cross(perpVec) * std::sin(randomAngle);
        } else if (forceJetDirectionSmudge) {
          XYZVector refVec(1., 0., 0.);
          if (std::abs(leadPUnitVec.Dot(refVec)) > 0.99)
            refVec = XYZVector(0., 1., 0.);
          XYZVector perpVec = leadPUnitVec.Cross(refVec).Unit();
          double smearAzimuth = randomGen.Uniform(0., o2::constants::math::TwoPI);
          XYZVector smearAxis = perpVec * std::cos(smearAzimuth) + leadPUnitVec.Cross(perpVec) * std::sin(smearAzimuth);
          double smearAngle = std::abs(randomGen.Gaus(0., 0.05 * jetR));
          leadPUnitVec = leadPUnitVec * std::cos(smearAngle) + smearAxis.Cross(leadPUnitVec) * std::sin(smearAngle);
        }
      }

      // 3) Checking if the event has a leading jet:
      auto jetsInColl = jets.sliceBy(perColJets, collId);
      float leadingJetPt = -1.;
      float subleadingJetPt = -1.;
      // std::optional avoids undefined behaviour from a default-constructed iterator:
      // (essentially, just protection for when we fetch jetEta() and the such)
      std::optional<o2::aod::RingJets::iterator> leadingJet;
      std::optional<o2::aod::RingJets::iterator> subleadingJet;
      for (auto const& jet : jetsInColl) {
        const auto jetpt = jet.jetPt();
        if (jetpt > leadingJetPt) {
          // Current leading becomes subleading:
          subleadingJetPt = leadingJetPt;
          subleadingJet = leadingJet; // may still be std::nullopt on first pass -- that is fine!
          // Now update the leading jet:
          leadingJetPt = jetpt;
          leadingJet = jet;
        } else if (jetpt > subleadingJetPt) { // Update subleading only:
          subleadingJetPt = jetpt;
          subleadingJet = jet;
        }
      }

      // Some useful bools to check if we have a leading jet and a subleading jet:
      // const bool hasValidLeadingJet = leadingJetPt > 0.;
      const bool hasValidLeadingJet = leadingJetPt > minLeadJetPt; // Finer control on jet momentum
      const bool hasValidSubJet = subleadingJetPt > minSubLeadJetPt;

      // Build jet vectors (only when the corresponding jet exists):
      // Dummy initialisations are safe: all jet-dependent fills are gated on hasValidLeadingJet / hasValidSubJet.
      float leadingJetEta = 0.;
      float leadingJetPhi = 0.;
      XYZVector leadingJetUnitVec(1., 0., 0.); // dummy (overwritten below)
      if (hasValidLeadingJet) {
        leadingJetEta = leadingJet->jetEta();
        leadingJetPhi = leadingJet->jetPhi();
        // Using internal getters to make code cleaner:
        leadingJetUnitVec = XYZVector(leadingJet->jetPx(), leadingJet->jetPy(), leadingJet->jetPz()).Unit();
        histos.fill(HIST("JetKinematicsQA/hLeadJetEta"), leadingJetEta); // This will not be subject to the forcePerpToJet nor the forceJetDirectionSmudge QAs, for simplicity
        histos.fill(HIST("JetKinematicsQA/hJetCounterPtJet"), leadingJetPt);

        // QA block -- Purposefully changing the jet direction (should kill signal, if any):
        if (forcePerpToJet) { // Use modified jet direction (done outside loop to guarantee all V0s inside event use same fake jet)
          // First, we build a vector perpendicular to the jet by picking an arbitrary vector not parallel to the jet
          XYZVector refVec(1., 0., 0.);
          if (std::abs(leadingJetUnitVec.Dot(refVec)) > 0.99)
            refVec = XYZVector(0., 1., 0.);
          // Now we get a perpendicular vector to the jet direction:
          XYZVector perpVec = leadingJetUnitVec.Cross(refVec).Unit();
          // Now we rotate around the jet axis by a random angle, just to make sure we are not introducing a bias in the QA:
          // We will use Rodrigues' rotation formula (v_rot = v*cos(randomAngle) + (Jet \cross v)*sin(randomAngle))
          double randomAngle = randomGen.Uniform(0., o2::constants::math::TwoPI);
          leadingJetUnitVec = perpVec * std::cos(randomAngle) + leadingJetUnitVec.Cross(perpVec) * std::sin(randomAngle);
        } else if (forceJetDirectionSmudge) {
          // Smear the jet direction by a small random angle to estimate sensitivity to
          // jet axis uncertainty. We rotate the jet axis by angle theta around a uniformly
          // random perpendicular axis -- this is isotropic and coordinate-independent,
          // unlike smearing eta and phi separately (which would break azimuthal symmetry
          // around the jet axis and depend on where in eta the jet sits).

          // 1) We pick a uniformly random axis perpendicular to the jet.
          // (re-using the same Rodrigues formula as in the forcePerpToJet block above)
          XYZVector refVec(1., 0., 0.);
          if (std::abs(leadingJetUnitVec.Dot(refVec)) > 0.99)
            refVec = XYZVector(0., 1., 0.);
          XYZVector perpVec = leadingJetUnitVec.Cross(refVec).Unit();
          // Rotate perpVec around the jet axis by a uniform random azimuth to get
          // a uniformly distributed random perpendicular direction (the smear axis):
          double smearAzimuth = randomGen.Uniform(0., o2::constants::math::TwoPI);
          XYZVector smearAxis = perpVec * std::cos(smearAzimuth) + leadingJetUnitVec.Cross(perpVec) * std::sin(smearAzimuth);

          // Step 2: draw the smearing polar angle from a Gaussian:
          // sigma = 0.05 * R --> ~68% of events smeared within 5% of R,
          //                      ~95% of events smeared within 10% of R,
          //                       ~5% see a displacement > 0.1*R (a very "badly determined jet", for our QA purposes)
          // std::abs() folds the symmetric Gaussian onto a half-normal ([0, inf))
          // -- R is not really an angle: just gives me a scale for the angular shift I am performing.
          // -- This may pose problems for forward jets: a small displacement in \theta becomes a large displacement in \eta space
          double smearSigma = 0.05 * jetR;
          double smearAngle = std::abs(randomGen.Gaus(0., smearSigma));

          // Step 3: rotate the jet axis by smearAngle around smearAxis.
          // Rodrigues is v_rot = v*cos(theta) + (k \cross v)*sin(theta) + k*(k \cdot v)*(1-cos(theta))
          // But the last term vanishes because smearAxis is perpendicular to leadingJetUnitVec:
          leadingJetUnitVec = leadingJetUnitVec * std::cos(smearAngle) + smearAxis.Cross(leadingJetUnitVec) * std::sin(smearAngle);
          // Also, rotation preserves the norm, so no re-normalisation is needed for this to be a unit vector.
        }
      }

      float subleadingJetEta = 0.;
      float subleadingJetPhi = 0.;
      XYZVector subJetUnitVec(1., 0., 0.);
      if (hasValidSubJet) {
        subleadingJetEta = subleadingJet->jetEta();
        subleadingJetPhi = subleadingJet->jetPhi();
        subJetUnitVec = XYZVector(subleadingJet->jetPx(), subleadingJet->jetPy(), subleadingJet->jetPz()).Unit();
        histos.fill(HIST("JetKinematicsQA/hSubLeadJetEta"), subleadingJetEta); // Unmodified direction
        histos.fill(HIST("JetKinematicsQA/hJetCounterPt2ndJet"), subleadingJetPt);
        // QA: same direction-smearing/perp logic as for the leading jet estimator.
        if (forcePerpToJet) {
          XYZVector refVec(1., 0., 0.);
          if (std::abs(subJetUnitVec.Dot(refVec)) > 0.99)
            refVec = XYZVector(0., 1., 0.);
          XYZVector perpVec = subJetUnitVec.Cross(refVec).Unit();
          double randomAngle = randomGen.Uniform(0., o2::constants::math::TwoPI);
          subJetUnitVec = perpVec * std::cos(randomAngle) + subJetUnitVec.Cross(perpVec) * std::sin(randomAngle);
        } else if (forceJetDirectionSmudge) {
          XYZVector refVec(1., 0., 0.);
          if (std::abs(subJetUnitVec.Dot(refVec)) > 0.99)
            refVec = XYZVector(0., 1., 0.);
          XYZVector perpVec = subJetUnitVec.Cross(refVec).Unit();
          double smearAzimuth = randomGen.Uniform(0., o2::constants::math::TwoPI);
          XYZVector smearAxis = perpVec * std::cos(smearAzimuth) + subJetUnitVec.Cross(perpVec) * std::sin(smearAzimuth);
          double smearAngle = std::abs(randomGen.Gaus(0., 0.05 * jetR));
          subJetUnitVec = subJetUnitVec * std::cos(smearAngle) + smearAxis.Cross(subJetUnitVec) * std::sin(smearAngle);
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
      int NLambdaLikeV0s = 0;
      for (auto const& v0 : v0sInColl) {
        if (v0.isLambda() ^ v0.isAntiLambda()){ // XOR (only the non-ambiguous candidates)
          NLambdaLikeV0s++;
        }
      }

      // Initialize delta method accumulators:
      EventDeltaTracker trackRing, trackRingKinCuts, trackJetKinCuts, trackJetLambdaKinCuts;
      for (auto const& v0 : v0sInColl) {
        const bool isLambda = v0.isLambda();
        const bool isAntiLambda = v0.isAntiLambda();
        // For now, removing the ambiguous candidates from the analysis. Derived data permits handling both.
        // (From Podolanski-Armenteros plots, the population of ambiguous is ~2% without TOF, and without
        //  competing mass rejection. From those, ~99% seem to be K0s, so no real gain in considering the
        //  ambiguous candidates in the analysis)
        if (isLambda && isAntiLambda)
          continue;
        const float v0pt = v0.v0Pt();
        const float v0eta = v0.v0Eta();
        const float v0phi = v0.v0Phi();

        float v0LambdaLikeMass = 0; // Initialized just to catch any stray behavior
        float protonLikePt = 0;
        float protonLikeEta = 0;
        float protonLikePhi = 0;
        if (isLambda) {
          if (!analyseLambda)
            continue;
          v0LambdaLikeMass = v0.massLambda();
          protonLikePt = v0.posPt();
          protonLikeEta = v0.posEta();
          protonLikePhi = v0.posPhi();
        } else if (isAntiLambda) { // (TODO: add a split histogram where you consider Lambda and AntiLambda polarization separately?)
          if (!analyseAntiLambda)
            continue;
          v0LambdaLikeMass = v0.massAntiLambda();
          protonLikePt = v0.negPt();
          protonLikeEta = v0.negEta();
          protonLikePhi = v0.negPhi();
        }

        PtEtaPhiMVector lambdaLike4Vec(v0pt, v0eta, v0phi, v0LambdaLikeMass);
        PtEtaPhiMVector protonLike4Vec(protonLikePt, protonLikeEta, protonLikePhi, protonMass);
        float lambdaRapidity = lambdaLike4Vec.Rapidity(); // For further kinematic selections

        // Boosting proton into lambda frame:
        XYZVector beta = lambdaLike4Vec.BoostToCM(); // Boost trivector that goes from laboratory frame to Lambda's rest frame (convenient new function, different from TLorentzVector's BoostVector())
        auto protonLike4VecStar = ROOT::Math::VectorUtil::boost(protonLike4Vec, beta);

        // Getting unit vectors and 3-components:
        XYZVector lambdaLike3Vec = lambdaLike4Vec.Vect();
        auto lambdaLikeUnit3Vec = lambdaLike3Vec.Unit();
        XYZVector protonLikeStarUnit3Vec = protonLike4VecStar.Vect().Unit();

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
        float phiStar = std::atan2(pY, pX);

        // Another reconstruction efficiency measure:
        // (Formula is: p_{Lambda} \cross p_{Daughter}^{*} \cdot B, and B points in Z)
        if (analyseMagField) {
          auto crossGeom = lambdaLike3Vec.Cross(protonLikeStarUnit3Vec);
          const bool positiveGeom = crossGeom.Z() * magField > 0;

          if (isLambda && positiveGeom)
            histos.fill(HIST("HelicityEfficiencyQA/hLambdaMassDecayGeomRight"), v0LambdaLikeMass);
          else if (isLambda && !positiveGeom)
            histos.fill(HIST("HelicityEfficiencyQA/hLambdaMassDecayGeomLeft"), v0LambdaLikeMass);
          else if (isAntiLambda && positiveGeom)
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
          // Cross product
          XYZVector crossLeadP = leadPUnitVec.Cross(lambdaLike3Vec);
          ringObservableLeadP = protonLikeStarUnit3Vec.Dot(crossLeadP) / crossLeadP.R();
          // Adding the prefactor related to the CP-violating decay (decay constants have different signs)
          if (!forcePolSignQA)
            ringObservableLeadP *= (isLambda) ? polPrefactorLambda : polPrefactorAntiLambda;
          else
            ringObservableLeadP *= (isLambda) ? polPrefactorLambda : -1.0 * polPrefactorAntiLambda;
          // Angular variables
          deltaPhiLeadP = wrapToPiFast(v0phi - leadPPhi); // Wrapped to [-PI, PI), for convenience
          
          cosDeltaThetaLeadP = leadPUnitVec.Dot(lambdaLikeUnit3Vec); // Uses the pre-calculated unit vectors to avoid recomputation
          deltaThetaLeadP = std::acos(cosDeltaThetaLeadP); // 3D angular separation. Same as ROOT::Math::VectorUtil::Angle(leadPUnitVec, lambdaLike3Vec);
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
          // Cross product
          XYZVector cross = leadingJetUnitVec.Cross(lambdaLike3Vec);
          ringObservable = protonLikeStarUnit3Vec.Dot(cross) / cross.R();
          // Adding prefactor
          if (!forcePolSignQA)
            ringObservable *= (isLambda) ? polPrefactorLambda : polPrefactorAntiLambda;
          else
            ringObservable *= (isLambda) ? polPrefactorLambda : -1.0 * polPrefactorAntiLambda;
          // Angular variables
          deltaPhiJet = wrapToPiFast(v0phi - leadingJetPhi);
          deltaEtaJet = v0eta - leadingJetEta;

          cosDeltaThetaJet = leadingJetUnitVec.Dot(lambdaLikeUnit3Vec);
          deltaThetaJet = std::acos(cosDeltaThetaJet);
          
          // Testing an invariance -- <R>/\hat{t}_z -- Possible source of an artificial (trivial) sign flip in the observable:
          Jz = leadingJetUnitVec.Z();
          if (std::abs(Jz) > 1e-4)
            ringObservableOverJetZ = ringObservable/Jz;
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
          if (!forcePolSignQA)
            ringObservable2ndJet *= (isLambda) ? polPrefactorLambda : polPrefactorAntiLambda;
          else
            ringObservable2ndJet *= (isLambda) ? polPrefactorLambda : -1.0 * polPrefactorAntiLambda;
          // Angular variables
          deltaPhi2ndJet = wrapToPiFast(v0phi - subleadingJetPhi);
          cosDeltaTheta2ndJet = subJetUnitVec.Dot(lambdaLikeUnit3Vec);
          deltaTheta2ndJet = std::acos(cosDeltaTheta2ndJet);
        }

        // Calculating polarization observables (in the Lambda frame, because that is easier -- does not require boosts):
        // To be precise, not actually the polarization, but a part of the summand in P^*_Lambda = (3/\alpha_Lambda) * <p^*_{proton}>
        float PolStarX = 0, PolStarY = 0, PolStarZ = 0; // Dummy initialization: avoid warnings in compile time
        if (isLambda) {                                 // Notice there is no need to check analyseLambda again due to previous checks.
          PolStarX = polPrefactorLambda * protonLikeStarUnit3Vec.X();
          PolStarY = polPrefactorLambda * protonLikeStarUnit3Vec.Y();
          PolStarZ = polPrefactorLambda * protonLikeStarUnit3Vec.Z();
        } else if (isAntiLambda) {
          PolStarX = polPrefactorAntiLambda * protonLikeStarUnit3Vec.X();
          PolStarY = polPrefactorAntiLambda * protonLikeStarUnit3Vec.Y();
          PolStarZ = polPrefactorAntiLambda * protonLikeStarUnit3Vec.Z();
        }

        float v0phiToFillHists = wrapToPiFast(v0phi); // A short wrap to reuse some predefined axes

        // Fill ring histograms: (1D, lambda 2D correlations and jet 2D correlations):
        if (hasValidLeadingP) {
          RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "Ring")       // Notice the usage of macros! If you change the variable names, this WILL break the code!
                                                                          // No, there should NOT be any ";" here! Read the macro definition for an explanation
          histos.fill(HIST("IntegratedCuts/pRingCutsLeadingP"), 0, ringObservableLeadP); // First bin of comparison
          histos.fill(HIST("IntegratedCuts/hCountCutsLeadingP"), 0);

          // Filling checks that rely on Eta>0 or Eta<0 checks for V0 and LeadingP eta:
          RING_OBSERVABLE_LEADP_ETA_SPLIT_FILL_LIST("Ring", leadPEtaPos, lambdaEtaPos);
        }
        POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "Ring")

        // Binary search using the pre-fetched axes for delta method of error bar estimation:
        int binPt = 0; // Dummy declarations
        int binMass = 0;
        int binDTheta = 0;
        if (hasValidLeadingJet) {
          RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "Ring")
          histos.fill(HIST("IntegratedCuts/pRingCuts"), 0, ringObservable);
          histos.fill(HIST("IntegratedCuts/hCountCuts"), 0);
          histos.fill(HIST("IntegratedCuts/pRingVsNV0s"), NLambdaLikeV0s, ringObservable);

          // Properly fetching values as they are needed:
          binPt = mAxisPt->FindBin(v0pt);
          binMass = mAxisMass->FindBin(v0LambdaLikeMass);
          binDTheta = mAxisDTheta->FindBin(deltaThetaJet);
          trackRing.addV0(ringObservable, binPt, binMass, binDTheta);

          // Measuring the AEE differentially (azimuthal efficiency effect, which causes different V0 topologies to be enhanced/suppressed):
          if (isLambda)
            histos.fill(HIST("HelicityEfficiencyQA/p2dRing_LambdaMassVsPhiLambdaMinusPhiProtonStar"), v0LambdaLikeMass, deltaPhiLambdaProtonStar, ringObservable);
          else
            histos.fill(HIST("HelicityEfficiencyQA/p2dRing_AntiLambdaMassVsPhiLambdaMinusPhiProtonStar"), v0LambdaLikeMass, deltaPhiLambdaProtonStar, ringObservable);
          // AEE and HEE correlation:
          histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsJet_CosThetaVsPhiStar"), cosFakePol, phiStar);
          histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignal_CosThetaVsPhiStar"), cosFakePol, phiStar, ringObservable);

          histos.fill(HIST("HelicityEfficiencyQA/pRingVsJetZcomponent"), Jz, ringObservable);
          histos.fill(HIST("HelicityEfficiencyQA/pRingOverJetZcomponent_VsJetEta"), leadingJetEta, ringObservableOverJetZ);
          histos.fill(HIST("HelicityEfficiencyQA/pRingOverJetZcomponent_VsCosThetaHEE"), cosFakePol, ringObservableOverJetZ);
          histos.fill(HIST("HelicityEfficiencyQA/pRingOverJetZcomponent_VsPhiStar"), phiStar, ringObservableOverJetZ);
          histos.fill(HIST("HelicityEfficiencyQA/pRingOverJetZcomponent_VsJetEtaVsCosThetaHEE"), leadingJetEta, cosFakePol, ringObservableOverJetZ);
          histos.fill(HIST("HelicityEfficiencyQA/pRingOverJetZcomponent_VsJetEtaVsPhiStar"), leadingJetEta, phiStar, ringObservableOverJetZ);
        }
        if (hasValidSubJet) {
          RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "Ring")
          histos.fill(HIST("IntegratedCuts/pRingCutsSubLeadingJet"), 0, ringObservable2ndJet);
          histos.fill(HIST("IntegratedCuts/hCountCutsSubLeadingJet"), 0);
        }

        // Filling eta dependence QAs of the result (both for V0 and jet proxy):
        if (hasValidLeadingJet) {
          histos.fill(HIST("ProxyEta/pRingEtaCuts"), 0, ringObservable);
          histos.fill(HIST("ProxyEta/pRingEtaCuts"), lambdaEtaPos ? 3 : 4, ringObservable);

          // HEE study (helicity efficiency effect):
          histos.fill(HIST("HelicityEfficiencyQA/hFakePolCounts"), cosFakePol, 0);
          histos.fill(HIST("HelicityEfficiencyQA/hFakePolCounts"), cosFakePol, lambdaEtaPos ? 3 : 4);
            // Same for signal:
          histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"), cosFakePol, 0, ringObservable);
          histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"), cosFakePol, lambdaEtaPos ? 3 : 4, ringObservable);
          // Counter and ring accumulators for AEE study:
          histos.fill(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"), phiStar, 0);
          histos.fill(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"), phiStar, lambdaEtaPos ? 3 : 4);
          histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"), phiStar, 0, ringObservable);
          histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"), phiStar, lambdaEtaPos ? 3 : 4, ringObservable);

          histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"), deltaPhiLambdaProtonStar, 0, ringObservable);
          histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"), deltaPhiLambdaProtonStar, lambdaEtaPos ? 3 : 4, ringObservable);

          // Extra correlations test:
          histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaJet"), cosFakePol, deltaThetaJet);
          histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsCosThetaVsPtForJets"), cosFakePol, v0pt);

          if (pTLambdaCheck) {
            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"), cosFakePol, 0);
            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"), cosFakePol, lambdaEtaPos ? 3 : 4);
            if (rapidityLambdaCheck){ // Stricter check
              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"), cosFakePol, 0);
              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"), cosFakePol, lambdaEtaPos ? 3 : 4);
            }
          }
          if (jetEtaPos) { // Less readable than "if ( jetEtaPos &&  lambdaEtaPos)", yet more efficient
            histos.fill(HIST("ProxyEta/pRingEtaCuts"), 1, ringObservable);
            histos.fill(HIST("ProxyEta/pRingEtaCuts"), lambdaEtaPos ? 5 : 6, ringObservable);

            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCounts"), cosFakePol, 1);
            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCounts"), cosFakePol, lambdaEtaPos ? 5 : 6);
              // Same for signal:
            histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"), cosFakePol, 1, ringObservable);
            histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"), cosFakePol, lambdaEtaPos ? 5 : 6, ringObservable);
            // Counter and ring accumulators for AEE study:
            histos.fill(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"), phiStar, 1);
            histos.fill(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"), phiStar, lambdaEtaPos ? 5 : 6);
            histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"), phiStar, 1, ringObservable);
            histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"), phiStar, lambdaEtaPos ? 5 : 6, ringObservable);

            histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"), deltaPhiLambdaProtonStar, 1, ringObservable);
            histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"), deltaPhiLambdaProtonStar, lambdaEtaPos ? 5 : 6, ringObservable);

            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaJetPosEta"), cosFakePol, deltaThetaJet);
            if (pTLambdaCheck) {
              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"), cosFakePol, 1);
              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"), cosFakePol, lambdaEtaPos ? 5 : 6);
              if (rapidityLambdaCheck){ // Stricter check
                histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"), cosFakePol, 1);
                histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"), cosFakePol, lambdaEtaPos ? 5 : 6);
              }
            }
            if (jetEtaStrict) { // eta_{Jet} >= R
              histos.fill(HIST("ProxyEta/pRingEtaCuts"),  9, ringObservable);
              histos.fill(HIST("ProxyEta/pRingEtaCuts"), lambdaEtaPos ? 11 : 12, ringObservable);
            }
          } else {
            histos.fill(HIST("ProxyEta/pRingEtaCuts"), 2, ringObservable);
            histos.fill(HIST("ProxyEta/pRingEtaCuts"), lambdaEtaPos ? 7 : 8, ringObservable);

            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCounts"), cosFakePol, 2);
            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCounts"), cosFakePol, lambdaEtaPos ? 7 : 8);
              // Same for signal:
            histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"), cosFakePol, 2, ringObservable);
            histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalVsCosTheta"), cosFakePol, lambdaEtaPos ? 7 : 8, ringObservable);
            // Counter and ring accumulators for AEE study:
            histos.fill(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"), phiStar, 2);
            histos.fill(HIST("HelicityEfficiencyQA/hCountsVsPhiStar"), phiStar, lambdaEtaPos ? 7 : 8);
            histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"), phiStar, 2, ringObservable);
            histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiStar"), phiStar, lambdaEtaPos ? 7 : 8, ringObservable);

            histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"), deltaPhiLambdaProtonStar, 2, ringObservable);
            histos.fill(HIST("HelicityEfficiencyQA/pFakePolSignalvsPhiLambdaMinusPhiProtonStar"), deltaPhiLambdaProtonStar, lambdaEtaPos ? 7 : 8, ringObservable);

            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaJetNegEta"), cosFakePol, deltaThetaJet);

            if (pTLambdaCheck) {
              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"), cosFakePol, 2);
              histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtCut"), cosFakePol, lambdaEtaPos ? 7 : 8);
              if (rapidityLambdaCheck){ // Stricter check
                histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"), cosFakePol, 2);
                histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLambdaPtYCuts"), cosFakePol, lambdaEtaPos ? 7 : 8);
              }
            }
            if (jetEtaStrict) { // eta_{Jet} <= -R
              histos.fill(HIST("ProxyEta/pRingEtaCuts"), 10, ringObservable);
              histos.fill(HIST("ProxyEta/pRingEtaCuts"), lambdaEtaPos ? 13 : 14, ringObservable);
            }
          }
        }
        if (hasValidSubJet) {
          histos.fill(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"), 0, ringObservable2ndJet);
          histos.fill(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"), lambdaEtaPos ? 3 : 4, ringObservable2ndJet);
          if (subJetEtaPos) {
            histos.fill(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"), 1, ringObservable2ndJet);
            histos.fill(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"), lambdaEtaPos ? 5 : 6, ringObservable2ndJet);
            if (subJetEtaStrict) { // eta_{SubJet} >= R
              histos.fill(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"),  9, ringObservable2ndJet);
              histos.fill(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"), lambdaEtaPos ? 11 : 12, ringObservable2ndJet);
            }
          } else {
            histos.fill(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"), 2, ringObservable2ndJet);
            histos.fill(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"), lambdaEtaPos ? 7 : 8, ringObservable2ndJet);
            if (subJetEtaStrict) { // eta_{SubJet} <= -R
              histos.fill(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"), 10, ringObservable2ndJet);
              histos.fill(HIST("ProxyEta/pRingEtaCutsSubLeadingJet"), lambdaEtaPos ? 13 : 14, ringObservable2ndJet);
            }
          }
        }
        if (hasValidLeadingP) {
          histos.fill(HIST("ProxyEta/pRingEtaCutsLeadingP"), 0, ringObservableLeadP);
          histos.fill(HIST("ProxyEta/pRingEtaCutsLeadingP"), lambdaEtaPos ? 3 : 4, ringObservableLeadP);

          histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"), cosFakePol, 0);
          histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"), cosFakePol, lambdaEtaPos ? 3 : 4);

          histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsCosThetaVsPtForLeadP"), cosFakePol, v0pt); // Understanding the population of events that has a leading particle (even though this does not need one to be calculated!)
          histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaLeadP"), cosFakePol, deltaThetaLeadP);

          // Extra correlations test:
          histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"), cosFakePol, deltaThetaLeadP);
          if (leadPEtaPos) {
            histos.fill(HIST("ProxyEta/pRingEtaCutsLeadingP"), 1, ringObservableLeadP);
            histos.fill(HIST("ProxyEta/pRingEtaCutsLeadingP"), lambdaEtaPos ? 5 : 6, ringObservableLeadP);

            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"), cosFakePol, 1);
            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"), cosFakePol, lambdaEtaPos ? 5 : 6);

            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaLeadPPosEta"), cosFakePol, deltaThetaLeadP);
          }
          else {
            histos.fill(HIST("ProxyEta/pRingEtaCutsLeadingP"), 2, ringObservableLeadP);
            histos.fill(HIST("ProxyEta/pRingEtaCutsLeadingP"), lambdaEtaPos ? 7 : 8, ringObservableLeadP);

            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"), cosFakePol, 2);
            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsLeadP"), cosFakePol, lambdaEtaPos ? 7 : 8);

            histos.fill(HIST("HelicityEfficiencyQA/hFakePolCountsVsDeltaThetaLeadPNegEta"), cosFakePol, deltaThetaLeadP);
          }
        }

        // Extra kinematic criteria for Lambda candidates (removes polarization background):
        if (kinematicLambdaCheck) {
          if (hasValidLeadingP) {
            RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
            histos.fill(HIST("IntegratedCuts/pRingCutsLeadingP"), 1, ringObservableLeadP);
            histos.fill(HIST("IntegratedCuts/hCountCutsLeadingP"), 1);

            // Filling checks that rely on Eta>0 or Eta<0 checks for V0 and LeadingP eta:
            RING_OBSERVABLE_LEADP_ETA_SPLIT_FILL_LIST("RingKinematicCuts", leadPEtaPos, lambdaEtaPos);
          }
          POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
          if (hasValidLeadingJet) {
            RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
            histos.fill(HIST("IntegratedCuts/pRingCuts"), 1, ringObservable);
            histos.fill(HIST("IntegratedCuts/hCountCuts"), 1);
            trackRingKinCuts.addV0(ringObservable, binPt, binMass, binDTheta);
          }
          if (hasValidSubJet) {
            RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
            histos.fill(HIST("IntegratedCuts/pRingCutsSubLeadingJet"), 1, ringObservable2ndJet);
            histos.fill(HIST("IntegratedCuts/hCountCutsSubLeadingJet"), 1);
          }
        }

        // Extra selection criteria on jet candidates:
        // (redundant for jets with R=0.4, but for jets with R<0.4 the leading jet may be farther in eta)
        if (kinematicJetCheck) { // Already includes hasValidLeadingJet in the bool! (no need to check again)
          RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
          histos.fill(HIST("IntegratedCuts/pRingCuts"), 2, ringObservable);
          histos.fill(HIST("IntegratedCuts/hCountCuts"), 2);
          POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
          trackJetKinCuts.addV0(ringObservable, binPt, binMass, binDTheta);
        }

        // Extra selection criteria on both Lambda and jet candidates:
        if (kinematicLambdaCheck && kinematicJetCheck) {
          RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
          histos.fill(HIST("IntegratedCuts/pRingCuts"), 3, ringObservable);
          histos.fill(HIST("IntegratedCuts/hCountCuts"), 3);
          POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
          trackJetLambdaKinCuts.addV0(ringObservable, binPt, binMass, binDTheta);
        }

        // Same variations for the leading particle and for the subleading jet:
        // (kinematicLeadPCheck already encodes hasValidLeadingP, so no extra gate needed here)
        if (kinematicLeadPCheck) {
          RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
          histos.fill(HIST("IntegratedCuts/pRingCutsLeadingP"), 2, ringObservableLeadP);
          histos.fill(HIST("IntegratedCuts/hCountCutsLeadingP"), 2);

          // Filling checks that rely on Eta>0 or Eta<0 checks for V0 and LeadingP eta:
          RING_OBSERVABLE_LEADP_ETA_SPLIT_FILL_LIST("JetKinematicCuts", leadPEtaPos, lambdaEtaPos);
        }
        if (kinematic2ndJetCheck) {
          RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
          histos.fill(HIST("IntegratedCuts/pRingCutsSubLeadingJet"), 2, ringObservable2ndJet);
          histos.fill(HIST("IntegratedCuts/hCountCutsSubLeadingJet"), 2);
        }
        if (kinematicLambdaCheck && kinematicLeadPCheck) {
          RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
          histos.fill(HIST("IntegratedCuts/pRingCutsLeadingP"), 3, ringObservableLeadP);
          histos.fill(HIST("IntegratedCuts/hCountCutsLeadingP"), 3);

          // Filling checks that rely on Eta>0 or Eta<0 checks for V0 and LeadingP eta:
          RING_OBSERVABLE_LEADP_ETA_SPLIT_FILL_LIST("JetAndLambdaKinematicCuts", leadPEtaPos, lambdaEtaPos);
        }
        if (kinematicLambdaCheck && kinematic2ndJetCheck) {
          RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
          histos.fill(HIST("IntegratedCuts/pRingCutsSubLeadingJet"), 3, ringObservable2ndJet);
          histos.fill(HIST("IntegratedCuts/hCountCutsSubLeadingJet"), 3);
        }
      } // end v0s loop

      // Flush trackers to the actual O2 histograms (via macros, so that O2 compiles properly):
      FLUSH_DELTA_TRACKER("Ring", trackRing, mAxisPt, mAxisMass, mAxisDTheta)
      FLUSH_DELTA_TRACKER("RingKinematicCuts", trackRingKinCuts, mAxisPt, mAxisMass, mAxisDTheta)
      FLUSH_DELTA_TRACKER("JetKinematicCuts", trackJetKinCuts, mAxisPt, mAxisMass, mAxisDTheta)
      FLUSH_DELTA_TRACKER("JetAndLambdaKinematicCuts", trackJetLambdaKinCuts, mAxisPt, mAxisMass, mAxisDTheta)
    } // end collisions
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
