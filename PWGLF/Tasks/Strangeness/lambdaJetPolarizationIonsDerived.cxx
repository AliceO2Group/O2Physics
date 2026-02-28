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

#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

// Custom data model:
#include "PWGLF/DataModel/lambdaJetPolarizationIons.h"

#include <cmath>
#include <map>
#include <string>
#include <vector>

// #include <TLorentzVector.h>
// #include <TVector3.h>
// New recommended format:
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <TRandom3.h> // For perpendicular jet direction QAs

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using ROOT::Math::XYZVector;
using ROOT::Math::PtEtaPhiMVector;
// using namespace o2::aod::lambdajetpol; // Used it explicitly along the code for clarity

    // Declaring constants:
constexpr double protonMass = o2::constants::physics::MassProton; // Assumes particle identification for daughter is perfect
constexpr double lambdaWeakDecayConstant = 0.749; // DPG 2025 update
constexpr double antiLambdaWeakDecayConstant = -0.758; // DPG 2025 update
constexpr double polPrefactorLambda = 3.0/lambdaWeakDecayConstant;
constexpr double polPrefactorAntiLambda = 3.0/antiLambdaWeakDecayConstant;

// Helper macro to avoid writing the histogram fills 4 times for about 20 histograms:
#define RING_OBSERVABLE_FILL_LIST(X, FOLDER)                                                          \
    /* 1D observable histograms */                                                                    \
    X(FOLDER "/hRingObservableDeltaPhi",                     deltaPhiJet,   ringObservable)           \
    X(FOLDER "/hRingObservableDeltaTheta",                   deltaThetaJet, ringObservable)           \
    X(FOLDER "/hRingObservableIntegrated",                   0.,            ringObservable)           \
    /* Counters */                                                                                    \
    X(FOLDER "/hDeltaPhi",                                   deltaPhiJet)                             \
    X(FOLDER "/hDeltaTheta",                                 deltaThetaJet)                           \
    X(FOLDER "/hIntegrated",                                 0.)                                      \
    /* Lambda pT variation -- Youpeng's proposal */                                                   \
    X(FOLDER "/hRingObservableLambdaPt",                     v0pt,           ringObservable)          \
    X(FOLDER "/hLambdaPt",                                   v0pt)                                    \
    /* 2D Lambda correlations */                                                                      \
    X(FOLDER "/h2dRingObservableDeltaPhiVsLambdaPt",         deltaPhiJet,   v0pt, ringObservable)     \
    X(FOLDER "/h2dRingObservableDeltaThetaVsLambdaPt",       deltaThetaJet, v0pt, ringObservable)     \
    /* Counters */                                                                                    \
    X(FOLDER "/h2dDeltaPhiVsLambdaPt",                       deltaPhiJet,   v0pt)                     \
    X(FOLDER "/h2dDeltaThetaVsLambdaPt",                     deltaThetaJet, v0pt)                     \
    /* 2D Jet correlations */                                                                         \
    X(FOLDER "/h2dRingObservableDeltaPhiVsLeadJetPt",        deltaPhiJet,   leadingJetPt, ringObservable) \
    X(FOLDER "/h2dRingObservableDeltaThetaVsLeadJetPt",      deltaThetaJet, leadingJetPt, ringObservable) \
    /* Counters */                                                                                    \
    X(FOLDER "/h2dDeltaPhiVsLeadJetPt",                      deltaPhiJet,   leadingJetPt)             \
    X(FOLDER "/h2dDeltaThetaVsLeadJetPt",                    deltaThetaJet, leadingJetPt)             \
    /* Additional plots for instant gratification - 1D Profiles */                                    \
    X(FOLDER "/pRingObservableDeltaPhi",                     deltaPhiJet,   ringObservable)           \
    X(FOLDER "/pRingObservableDeltaTheta",                   deltaThetaJet, ringObservable)           \
    X(FOLDER "/pRingObservableIntegrated",                   0.,            ringObservable)           \
    X(FOLDER "/pRingObservableLambdaPt",                     v0pt,          ringObservable)           \
    /* 2D Profiles */                                                                                 \
    X(FOLDER "/p2dRingObservableDeltaPhiVsLambdaPt",         deltaPhiJet,   v0pt, ringObservable)     \
    X(FOLDER "/p2dRingObservableDeltaThetaVsLambdaPt",       deltaThetaJet, v0pt, ringObservable)     \
    X(FOLDER "/p2dRingObservableDeltaPhiVsLeadJetPt",        deltaPhiJet,   leadingJetPt, ringObservable) \
    X(FOLDER "/p2dRingObservableDeltaThetaVsLeadJetPt",      deltaThetaJet, leadingJetPt, ringObservable) \
    /* 1D Mass */                                                                                     \
    X(FOLDER "/hMass", v0LambdaLikeMass)                                                              \
    X(FOLDER "/hRingObservableMass", v0LambdaLikeMass, ringObservable)                                \
    X(FOLDER "/hMassSigExtract", v0LambdaLikeMass)                                                    \
    /* 2D: Observable vs Mass */                                                                      \
    X(FOLDER "/h2dRingObservableDeltaPhiVsMass",      deltaPhiJet,   v0LambdaLikeMass, ringObservable) \
    X(FOLDER "/h2dRingObservableDeltaThetaVsMass",    deltaThetaJet, v0LambdaLikeMass, ringObservable) \
    /* Counters */                                                                                    \
    X(FOLDER "/h2dDeltaPhiVsMass",                    deltaPhiJet,   v0LambdaLikeMass)                \
    X(FOLDER "/h2dDeltaThetaVsMass",                  deltaThetaJet, v0LambdaLikeMass)                \
    /* 3D: Observable vs Mass vs Lambda pT */                                                         \
    X(FOLDER "/h3dRingObservableDeltaPhiVsMassVsLambdaPt",   deltaPhiJet,   v0LambdaLikeMass, v0pt, ringObservable) \
    X(FOLDER "/h3dRingObservableDeltaThetaVsMassVsLambdaPt", deltaThetaJet, v0LambdaLikeMass, v0pt, ringObservable) \
    /* Counters */                                                                                    \
    X(FOLDER "/h3dDeltaPhiVsMassVsLambdaPt",          deltaPhiJet,   v0LambdaLikeMass, v0pt)          \
    X(FOLDER "/h3dDeltaThetaVsMassVsLambdaPt",        deltaThetaJet, v0LambdaLikeMass, v0pt)          \
    /* 3D: Observable vs Mass vs Lead Jet pT */                                                       \
    X(FOLDER "/h3dRingObservableDeltaPhiVsMassVsLeadJetPt",   deltaPhiJet,   v0LambdaLikeMass, leadingJetPt, ringObservable) \
    X(FOLDER "/h3dRingObservableDeltaThetaVsMassVsLeadJetPt", deltaThetaJet, v0LambdaLikeMass, leadingJetPt, ringObservable) \
    /* Counters */                                                                                    \
    X(FOLDER "/h3dDeltaPhiVsMassVsLeadJetPt",                 deltaPhiJet,   v0LambdaLikeMass, leadingJetPt) \
    X(FOLDER "/h3dDeltaThetaVsMassVsLeadJetPt",               deltaThetaJet, v0LambdaLikeMass, leadingJetPt) \
    /* 2D: Observable vs Mass vs Centrality (projected as 2D Mass vs Cent for integrated observable) */ \
    X(FOLDER "/h2dRingObservableMassVsCent", v0LambdaLikeMass, centrality, ringObservable)            \
    /* 3D: Observable vs Mass vs Centrality */                                                        \
    X(FOLDER "/h3dRingObservableDeltaPhiVsMassVsCent",   deltaPhiJet,   v0LambdaLikeMass, centrality, ringObservable) \
    X(FOLDER "/h3dRingObservableDeltaThetaVsMassVsCent", deltaThetaJet, v0LambdaLikeMass, centrality, ringObservable) \
    /* Counters */                                                                                    \
    X(FOLDER "/h3dDeltaPhiVsMassVsCent",   deltaPhiJet,   v0LambdaLikeMass, centrality)               \
    X(FOLDER "/h3dDeltaThetaVsMassVsCent", deltaThetaJet, v0LambdaLikeMass, centrality)               \
    /* TProfile of Ring vs Mass */                                                                    \
    X(FOLDER "/pRingObservableMass", v0LambdaLikeMass, ringObservable)                                \
    /* TProfile of Ring vs Mass -- Leading Particle and 2nd-to-leading jet - QA */                    \
    X(FOLDER "/pRingObservableLeadPMass", v0LambdaLikeMass, ringObservableLeadP)                      \
    X(FOLDER "/pRingObservable2ndJetMass", v0LambdaLikeMass, ringObservable2ndJet)                    \
    /* 2D Profiles: Angle vs Mass */                                                                  \
    X(FOLDER "/p2dRingObservableDeltaPhiVsMass",   deltaPhiJet,   v0LambdaLikeMass, ringObservable)   \
    X(FOLDER "/p2dRingObservableDeltaThetaVsMass", deltaThetaJet, v0LambdaLikeMass, ringObservable)   \
    /* 3D Profiles: Angle vs Mass vs Lambda pT */                                                     \
    X(FOLDER "/p3dRingObservableDeltaPhiVsMassVsLambdaPt",   deltaPhiJet,   v0LambdaLikeMass, v0pt, ringObservable) \
    X(FOLDER "/p3dRingObservableDeltaThetaVsMassVsLambdaPt", deltaThetaJet, v0LambdaLikeMass, v0pt, ringObservable) \
    /* 3D Profiles: Angle vs Mass vs Lead Jet pT */                                                   \
    X(FOLDER "/p3dRingObservableDeltaPhiVsMassVsLeadJetPt",   deltaPhiJet,   v0LambdaLikeMass, leadingJetPt, ringObservable) \
    X(FOLDER "/p3dRingObservableDeltaThetaVsMassVsLeadJetPt", deltaThetaJet, v0LambdaLikeMass, leadingJetPt, ringObservable) \
    /* 2D Profile: Mass vs Centrality */                                                              \
    X(FOLDER "/p2dRingObservableMassVsCent", v0LambdaLikeMass, centrality, ringObservable)            \
    /* 3D Profiles: Angle vs Mass vs Centrality */                                                    \
    X(FOLDER "/p3dRingObservableDeltaPhiVsMassVsCent",   deltaPhiJet,   v0LambdaLikeMass, centrality, ringObservable) \
    X(FOLDER "/p3dRingObservableDeltaThetaVsMassVsCent", deltaThetaJet, v0LambdaLikeMass, centrality, ringObservable)
    // (TODO: add counters for regular TH2Ds about centrality)

// For leading particle
#define RING_OBSERVABLE_LEADP_FILL_LIST(X, FOLDER)                                    \
    X(FOLDER "/pRingObservableLeadPDeltaPhi",   deltaPhiLeadP,   ringObservableLeadP) \
    X(FOLDER "/pRingObservableLeadPDeltaTheta", deltaThetaLeadP, ringObservableLeadP) \
    X(FOLDER "/pRingObservableLeadPIntegrated", 0.,            ringObservableLeadP)   \
    X(FOLDER "/pRingObservableLeadPLambdaPt",   v0pt,          ringObservableLeadP)

// For subleading jet:
#define RING_OBSERVABLE_2NDJET_FILL_LIST(X, FOLDER)                                      \
    X(FOLDER "/pRingObservable2ndJetDeltaPhi",   deltaPhi2ndJet,   ringObservable2ndJet) \
    X(FOLDER "/pRingObservable2ndJetDeltaTheta", deltaTheta2ndJet, ringObservable2ndJet) \
    X(FOLDER "/pRingObservable2ndJetIntegrated", 0.,            ringObservable2ndJet)    \
    X(FOLDER "/pRingObservable2ndJetLambdaPt",   v0pt,          ringObservable2ndJet)


#define POLARIZATION_PROFILE_FILL_LIST(X, FOLDER) \
    /* =============================== */ \
    /* 1D TProfiles vs v0phi */ \
    /* =============================== */ \
    X(FOLDER "/pPxStarPhi",             v0phiToFillHists,           PolStarX) \
    X(FOLDER "/pPyStarPhi",             v0phiToFillHists,           PolStarY) \
    X(FOLDER "/pPzStarPhi",             v0phiToFillHists,           PolStarZ) \
    /* =============================== */ \
    /* 1D TProfiles vs DeltaPhi_jet */ \
    /* =============================== */ \
    X(FOLDER "/pPxStarDeltaPhi",        deltaPhiJet,   PolStarX) \
    X(FOLDER "/pPyStarDeltaPhi",        deltaPhiJet,   PolStarY) \
    X(FOLDER "/pPzStarDeltaPhi",        deltaPhiJet,   PolStarZ) \
    /* =============================== */ \
    /* 2D TProfiles vs DeltaPhi_jet and Lambda pT */ \
    /* =============================== */ \
    X(FOLDER "/p2dPxStarDeltaPhiVsLambdaPt",  deltaPhiJet, v0pt, PolStarX) \
    X(FOLDER "/p2dPyStarDeltaPhiVsLambdaPt",  deltaPhiJet, v0pt, PolStarY) \
    X(FOLDER "/p2dPzStarDeltaPhiVsLambdaPt",  deltaPhiJet, v0pt, PolStarZ)

// Apply the macros (notice I had to include the semicolon (";") after the function, so you don't need to 
// write that when calling this APPLY_HISTO_FILL. The code will look weird, but without this the compiler
// would not know to end each statement with a semicolon):
#define APPLY_HISTO_FILL(NAME, ...) histos.fill(HIST(NAME), __VA_ARGS__);


struct lambdajetpolarizationionsderived {

    // Define histogram registries:
    HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

    // Master analysis switches
    Configurable<bool> analyseLambda{"analyseLambda", true, "process Lambda-like candidates"};
    Configurable<bool> analyseAntiLambda{"analyseAntiLambda", false, "process AntiLambda-like candidates"};
    Configurable<bool> doPPAnalysis{"doPPAnalysis", false, "if in pp, set to true. Default is HI"};
    
    // QAs that purposefully break the analysis
        // -- All of these tests should give us zero signal if the source is truly Lambda Polarization from vortices
    Configurable<bool> forcePolSignQA{"forcePolSignQA", false, "force antiLambda decay constant to be positive: should kill all the signal, if any. For QA"};
    Configurable<bool> forcePerpToJet{"forcePerpToJet", false, "force jet direction to be perpendicular () to jet estimator. For QA"};

    /////////////////////////
    // Configurable blocks:
    // Histogram axes configuration:
    struct : ConfigurableGroup {
        std::string prefix = "axisConfigurations"; // JSON group name
        ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
        ConfigurableAxis axisPtCoarseQA{"axisPtCoarse", {VARIABLE_WIDTH, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 7.0f, 10.0f, 15.0f}, "pt axis for QA"};
        ConfigurableAxis axisLambdaMass{"axisLambdaMass", {450, 1.08f, 1.15f}, "Lambda mass in GeV/c"}; // Default is {200, 1.101f, 1.131f}

        // Jet axes:
        ConfigurableAxis axisLeadingParticlePt{"axisLeadingParticlePt",{100, 0.f, 200.f},"Leading particle p_{T} (GeV/c)"}; // Simpler version!
        ConfigurableAxis axisJetPt{"axisJetPt",{50, 0.f, 200.f},"Jet p_{t} (GeV)"};
        ConfigurableAxis axisDeltaTheta{"axisDeltaTheta", {40, 0, constants::math::PI}, "#Delta #theta_{jet}"};
        ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {40, -constants::math::PI, constants::math::PI}, "#Delta #phi_{jet}"};

        // Coarser axes for signal extraction:
        ConfigurableAxis axisPtSigExtract{"axisPtSigExtract", {VARIABLE_WIDTH, 0.0f, 0.25f, 0.5f, 0.75f, 1.0f, 1.25f, 1.5f, 2.0f, 2.5f, 3.0f, 4.0f, 6.0f, 8.0f, 10.0f, 15.0f, 20.0f, 30.0f, 50.0f}, "pt axis for signal extraction"};
        // ConfigurableAxis axisLambdaMassSigExtract{"axisLambdaMassSigExtract", {175, 1.08f, 1.15f}, "Lambda mass in GeV/c"}; // With a sigma of 0.002 GeV/c, this has about 5 bins per sigma, so that the window is properly grasped.
        // A coarser axis (sigma is still well estimated, with about 8 bins in the peak region)
        ConfigurableAxis axisLambdaMassSigExtract{
            "axisLambdaMassSigExtract", {VARIABLE_WIDTH,
            // Left sideband (7 bins, 0.004 width)
            1.0800, 1.0840, 1.0880, 1.0920,
            1.0960, 1.1000, 1.1040, 1.1080,
            // Fine peak region (8 bins, 0.0016 width)
            1.1096, 1.1112, 1.1128, 1.1144,
            1.1160, 1.1176, 1.1192, 1.1208,
            // Right sideband (7 bins, 0.004 width)
            1.1248, 1.1288, 1.1328, 1.1368,
            1.1408, 1.1448, 1.1488},
            "Lambda mass in GeV/c"
        };
        ConfigurableAxis axisLeadingParticlePtSigExtract{"axisLeadingParticlePtSigExtract", {VARIABLE_WIDTH, 0, 4, 8, 12, 16, 20, 25, 30, 35, 40, 60, 100, 200}, "Leading particle p_{T} (GeV/c)"}; // Simpler version!
        ConfigurableAxis axisJetPtSigExtract{"axisJetPtSigExtract", {VARIABLE_WIDTH, 0, 5, 10, 12, 16, 20, 25, 30, 35, 40, 60, 100, 200},"Jet p_{t} (GeV)"};

        // (TODO: add a lambdaPt axis that is pre-selected only on the 0.5 to 1.5 Pt region for the Ring observable with lambda cuts to not store a huge histogram with empty bins by construction)

        ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Centrality"};
    } axisConfigurations;


    // Helper functions:
    // Fast wrapping into [-PI, PI) (restricted to this interval for function speed)
    inline double wrapToPiFast(double phi){
        constexpr double TwoPi = o2::constants::math::TwoPI;
        constexpr double Pi = o2::constants::math::PI;
        if (phi >= Pi) phi -= TwoPi;
        else if (phi < -Pi) phi += TwoPi;
        return phi;
    }


    void init(InitContext const&){
        // Ring observable histograms:
        // Helper to register one full histogram family (kinematic cut variation of ring observable)
        auto addRingObservableFamily = [&](const std::string& folder){
            // ===============================
            // 1D observable histograms
            // ===============================
            histos.add((folder + "/hRingObservableDeltaPhi").c_str(), "hRingObservableDeltaPhi", kTH1D, {axisConfigurations.axisDeltaPhi}); // Not quite the ring observable itself: this is just the numerator of <R> = \sum R_i / N_\Lambda (error bars WILL be wrong without TProfile)
            histos.add((folder + "/hRingObservableDeltaTheta").c_str(), "hRingObservableDeltaTheta", kTH1D, {axisConfigurations.axisDeltaTheta});
            histos.add((folder + "/hRingObservableIntegrated").c_str(), "hRingObservableIntegrated", kTH1D, {{1, -0.5, 0.5}});
            // Counters (denominators)
            histos.add((folder + "/hDeltaPhi").c_str(), "hDeltaPhi", kTH1D, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/hDeltaTheta").c_str(), "hDeltaTheta", kTH1D, {axisConfigurations.axisDeltaTheta});
            histos.add((folder + "/hIntegrated").c_str(), "hIntegrated", kTH1D, {{1, -0.5, 0.5}});
            // ===============================
            // Lambda pT dependence
            // ===============================
            histos.add((folder + "/hRingObservableLambdaPt").c_str(), "hRingObservableLambdaPt", kTH1D, {axisConfigurations.axisPt});
            histos.add((folder + "/hLambdaPt").c_str(), "hLambdaPt", kTH1D, {axisConfigurations.axisPt});
            // ===============================
            // 2D Lambda correlations
            // ===============================
            histos.add((folder + "/h2dRingObservableDeltaPhiVsLambdaPt").c_str(), "h2dRingObservableDeltaPhiVsLambdaPt", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPt});
            histos.add((folder + "/h2dRingObservableDeltaThetaVsLambdaPt").c_str(), "h2dRingObservableDeltaThetaVsLambdaPt", kTH2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisPt});
            // Counters
            histos.add((folder + "/h2dDeltaPhiVsLambdaPt").c_str(), "h2dDeltaPhiVsLambdaPt", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPt});
            histos.add((folder + "/h2dDeltaThetaVsLambdaPt").c_str(), "h2dDeltaThetaVsLambdaPt", kTH2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisPt});
            // ===============================
            // 2D Jet correlations
            // ===============================
            histos.add((folder + "/h2dRingObservableDeltaPhiVsLeadJetPt").c_str(), "h2dRingObservableDeltaPhiVsLeadJetPt", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisJetPt});
            histos.add((folder + "/h2dRingObservableDeltaThetaVsLeadJetPt").c_str(), "h2dRingObservableDeltaThetaVsLeadJetPt", kTH2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisJetPt});
            // Counters
            histos.add((folder + "/h2dDeltaPhiVsLeadJetPt").c_str(), "h2dDeltaPhiVsLeadJetPt", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisJetPt});
            histos.add((folder + "/h2dDeltaThetaVsLeadJetPt").c_str(), "h2dDeltaThetaVsLeadJetPt", kTH2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisJetPt});
            
            // Additional plots for instant gratification:
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
            histos.add((folder + "/pRingObservableDeltaPhi").c_str(), "pRingObservableDeltaPhi;#Delta#varphi_{jet};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/pRingObservableDeltaTheta").c_str(), "pRingObservableDeltaTheta;#Delta#theta_{jet};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaTheta});
            histos.add((folder + "/pRingObservableIntegrated").c_str(), "pRingObservableIntegrated; ;<#it{R}>", kTProfile, {{1, -0.5, 0.5}});
            histos.add((folder + "/pRingObservableLambdaPt").c_str(), "pRingObservableLambdaPt;#it{p}_{T}^{#Lambda};<#it{R}>", kTProfile, {axisConfigurations.axisPt});

            // For the leading particle:
            histos.add((folder + "/pRingObservableLeadPDeltaPhi").c_str(), "pRingObservableLeadPDeltaPhi;#Delta#varphi_{leadP};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/pRingObservableLeadPDeltaTheta").c_str(), "pRingObservableLeadPDeltaTheta;#Delta#theta_{leadP};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaTheta});
            histos.add((folder + "/pRingObservableLeadPIntegrated").c_str(), "pRingObservableLeadPIntegrated; ;<#it{R}>", kTProfile, {{1, -0.5, 0.5}});
            histos.add((folder + "/pRingObservableLeadPLambdaPt").c_str(), "pRingObservableLeadPLambdaPt;#it{p}_{T}^{#Lambda};<#it{R}>", kTProfile, {axisConfigurations.axisPt});
            // For the second-to-leading jet:
            histos.add((folder + "/pRingObservable2ndJetDeltaPhi").c_str(), "pRingObservable2ndJetDeltaPhi;#Delta#varphi_{2ndJet};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/pRingObservable2ndJetDeltaTheta").c_str(), "pRingObservable2ndJetDeltaTheta;#Delta#theta_{2ndJet};<#it{R}>", kTProfile, {axisConfigurations.axisDeltaTheta});
            histos.add((folder + "/pRingObservable2ndJetIntegrated").c_str(), "pRingObservable2ndJetIntegrated; ;<#it{R}>", kTProfile, {{1, -0.5, 0.5}});
            histos.add((folder + "/pRingObservable2ndJetLambdaPt").c_str(), "pRingObservable2ndJetLambdaPt;#it{p}_{T}^{#Lambda};<#it{R}>", kTProfile, {axisConfigurations.axisPt});
            // ===============================
            // 2D TProfiles (Lambda correlations)
            // ===============================
            histos.add((folder + "/p2dRingObservableDeltaPhiVsLambdaPt").c_str(), "p2dRingObservableDeltaPhiVsLambdaPt;#Delta#varphi_{jet};#it{p}_{T}^{#Lambda};<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPt});
            histos.add((folder + "/p2dRingObservableDeltaThetaVsLambdaPt").c_str(), "p2dRingObservableDeltaThetaVsLambdaPt;#Delta#theta_{jet};#it{p}_{T}^{#Lambda};<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisPt});
            // ===============================
            // 2D TProfiles (Jet correlations)
            // ===============================
            histos.add((folder + "/p2dRingObservableDeltaPhiVsLeadJetPt").c_str(), "p2dRingObservableDeltaPhiVsLeadJetPt;#Delta#varphi_{jet};#it{p}_{T}^{lead jet};<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisJetPt});
            histos.add((folder + "/p2dRingObservableDeltaThetaVsLeadJetPt").c_str(), "p2dRingObservableDeltaThetaVsLeadJetPt;#Delta#theta_{jet};#it{p}_{T}^{lead jet};<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisJetPt});
            
            // ===============================
            // Multi-dimensional histograms for signal extraction
            // (Mass-dependent polarization extraction)
            // ===============================
            // (TODO: possibly remove all TH2Ds that are not counter histograms and deal only with TProfiles!)
            // Simple invariant mass plot for QA:
            histos.add((folder + "/hMass").c_str(), "hMass", kTH1D, {axisConfigurations.axisLambdaMass});
            histos.add((folder + "/hMassSigExtract").c_str(), "hMassSigExtract", kTH1D, {axisConfigurations.axisLambdaMassSigExtract});
            // 1D Mass dependence of observable:
                // Important to know if the signal varies with mass, or if the sideband signal subtraction will probably work well enough!
            histos.add((folder + "/hRingObservableMass").c_str(), "hRingObservableMass", kTH1D, {axisConfigurations.axisLambdaMassSigExtract});
            // --- 2D: Angular observable vs Invariant Mass ---
            // Ring observable weighted with R
            histos.add((folder + "/h2dRingObservableDeltaPhiVsMass").c_str(), "h2dRingObservableDeltaPhiVsMass", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract});
            histos.add((folder + "/h2dRingObservableDeltaThetaVsMass").c_str(), "h2dRingObservableDeltaThetaVsMass", kTH2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract});
            // --- Counters (denominators) ---
            histos.add((folder + "/h2dDeltaPhiVsMass").c_str(), "h2dDeltaPhiVsMass", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract});
            histos.add((folder + "/h2dDeltaThetaVsMass").c_str(), "h2dDeltaThetaVsMass", kTH2D,{axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract});
            // --- 3D: Angular observable vs Mass vs Lambda pT ---
            histos.add((folder + "/h3dRingObservableDeltaPhiVsMassVsLambdaPt").c_str(), "h3dRingObservableDeltaPhiVsMassVsLambdaPt", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
            histos.add((folder + "/h3dRingObservableDeltaThetaVsMassVsLambdaPt").c_str(), "h3dRingObservableDeltaThetaVsMassVsLambdaPt", kTH3D, {axisConfigurations.axisDeltaTheta,axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
            // Counters
            histos.add((folder + "/h3dDeltaPhiVsMassVsLambdaPt").c_str(), "h3dDeltaPhiVsMassVsLambdaPt", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
            histos.add((folder + "/h3dDeltaThetaVsMassVsLambdaPt").c_str(), "h3dDeltaThetaVsMassVsLambdaPt", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
            // --- 3D: Angular observable vs Mass vs Lead Jet pT ---
            histos.add((folder + "/h3dRingObservableDeltaPhiVsMassVsLeadJetPt").c_str(), "h3dRingObservableDeltaPhiVsMassVsLeadJetPt", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract}); 
            histos.add((folder + "/h3dRingObservableDeltaThetaVsMassVsLeadJetPt").c_str(), "h3dRingObservableDeltaThetaVsMassVsLeadJetPt", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract}); 
            // --- Counters ---
            histos.add((folder + "/h3dDeltaPhiVsMassVsLeadJetPt").c_str(), "h3dDeltaPhiVsMassVsLeadJetPt", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract}); 
            histos.add((folder + "/h3dDeltaThetaVsMassVsLeadJetPt").c_str(), "h3dDeltaThetaVsMassVsLeadJetPt", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract});

            /////////////////////////////////////
            /// TProfiles with the proper errors for quick glancing
            /////////////////////////////////////
            // TProfile of ring vs mass (integrated in all phi, and properly normalized by N_\Lambda):
            histos.add((folder + "/pRingObservableMass").c_str(), "pRingObservableMass;m_{p#pi};<#it{R}>", kTProfile, {axisConfigurations.axisLambdaMassSigExtract});
            histos.add((folder + "/pRingObservableLeadPMass").c_str(), "pRingObservableLeadPMass;m_{p#pi};<#it{R}>", kTProfile, {axisConfigurations.axisLambdaMassSigExtract});
            histos.add((folder + "/pRingObservable2ndJetMass").c_str(), "pRingObservableLeadPMass;m_{p#pi};<#it{R}>", kTProfile, {axisConfigurations.axisLambdaMassSigExtract});
                // TProfile2D: <R> vs Mass (DeltaPhi)
            histos.add((folder + "/p2dRingObservableDeltaPhiVsMass").c_str(), "p2dRingObservableDeltaPhiVsMass;#Delta#varphi;m_{p#pi};<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract});
                // TProfile2D: <R> vs Mass (DeltaTheta)
            histos.add((folder + "/p2dRingObservableDeltaThetaVsMass").c_str(), "p2dRingObservableDeltaThetaVsMass;#Delta#theta;m_{p#pi};<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract});

                // --- TProfile3D: <R> vs DeltaPhi vs Mass vs LambdaPt ---
            histos.add((folder + "/p3dRingObservableDeltaPhiVsMassVsLambdaPt").c_str(), "p3dRingObservableDeltaPhiVsMassVsLambdaPt;#Delta#varphi;m_{p#pi};p_{T}^{#Lambda};<#it{R}>", kTProfile3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
                // --- TProfile3D: <R> vs DeltaTheta vs Mass vs LambdaPt ---
            histos.add((folder + "/p3dRingObservableDeltaThetaVsMassVsLambdaPt").c_str(), "p3dRingObservableDeltaThetaVsMassVsLambdaPt;#Delta#theta;m_{p#pi};p_{T}^{#Lambda};<#it{R}>", kTProfile3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
                // --- TProfile3D: <R> vs DeltaPhi vs Mass vs LeadJetPt ---
            histos.add((folder + "/p3dRingObservableDeltaPhiVsMassVsLeadJetPt").c_str(), "p3dRingObservableDeltaPhiVsMassVsLeadJetPt;#Delta#varphi;m_{p#pi};p_{T}^{jet};<#it{R}>", kTProfile3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract});
            // --- TProfile3D: <R> vs DeltaTheta vs Mass vs LeadJetPt ---
            histos.add((folder + "/p3dRingObservableDeltaThetaVsMassVsLeadJetPt").c_str(), "p3dRingObservableDeltaThetaVsMassVsLeadJetPt;#Delta#theta;m_{p#pi};p_{T}^{jet};<#it{R}>", kTProfile3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract});

            // ===============================
            // Mass histograms with centrality
            // ===============================
            // 2D Mass dependence of observable vs Centrality:
                // Important to know if the signal varies with mass, or if the sideband signal subtraction will probably work well enough!
            histos.add((folder + "/h2dRingObservableMassVsCent").c_str(), "h2dRingObservableMassVsCent", kTH2D, {axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
            // --- 3D: Angular observable vs Invariant Mass ---
            histos.add((folder + "/h3dRingObservableDeltaPhiVsMassVsCent").c_str(), "h3dRingObservableDeltaPhiVsMassVsCent", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
            histos.add((folder + "/h3dRingObservableDeltaThetaVsMassVsCent").c_str(), "h3dRingObservableDeltaThetaVsMassVsCent", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
            // --- Counters (denominators) ---
            histos.add((folder + "/h3dDeltaPhiVsMassVsCent").c_str(), "h3dDeltaPhiVsMassVsCent", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
            histos.add((folder + "/h3dDeltaThetaVsMassVsCent").c_str(), "h3dDeltaThetaVsMassVsCent", kTH3D,{axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});

            // Useful TProfiles:
                // --- TProfile2D: <R> vs Mass vs Centrality ---
            histos.add((folder + "/p2dRingObservableMassVsCent").c_str(), "p2dRingObservableMassVsCent;m_{p#pi};Centrality;<#it{R}>", kTProfile2D, {axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
                // --- TProfile3D: <R> vs DeltaPhi vs Mass vs Centrality ---
            histos.add((folder + "/p3dRingObservableDeltaPhiVsMassVsCent").c_str(), "p3dRingObservableDeltaPhiVsMassVsCent;#Delta#varphi;m_{p#pi};Centrality;<#it{R}>", kTProfile3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
                // --- TProfile3D: <R> vs DeltaTheta vs Mass vs Centrality ---
            histos.add((folder + "/p3dRingObservableDeltaThetaVsMassVsCent").c_str(), "p3dRingObservableDeltaThetaVsMassVsCent;#Delta#theta;m_{p#pi};Centrality;<#it{R}>", kTProfile3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});

            // ===============================
            //   Polarization observable QAs 
            // (not Ring: actual polarization!)
            // ===============================
                // Will implement these as TProfiles, as polarization is also a measure like P_\Lambda = (3/\alpha_\Lambda) * <p_{proton}>, so the error is similar
            // ===============================
            // 1D TProfiles
            // ===============================
            histos.add((folder + "/pPxStarPhi").c_str(), "pPxStarPhi;#varphi_{#Lambda};<P_{#Lambda}>_{x}", kTProfile, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/pPyStarPhi").c_str(), "pPyStarPhi;#varphi_{#Lambda};<P_{#Lambda}>_{y}", kTProfile, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/pPzStarPhi").c_str(), "pPzStarPhi;#varphi_{#Lambda};<P_{#Lambda}>_{z}", kTProfile, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/pPxStarDeltaPhi").c_str(), "pPxStarDeltaPhi;#Delta#varphi_{jet};<P_{#Lambda}>_{x}", kTProfile, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/pPyStarDeltaPhi").c_str(), "pPyStarDeltaPhi;#Delta#varphi_{jet};<P_{#Lambda}>_{y}", kTProfile, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/pPzStarDeltaPhi").c_str(), "pPzStarDeltaPhi;#Delta#varphi_{jet};<P_{#Lambda}>_{z}", kTProfile, {axisConfigurations.axisDeltaPhi});
            // ===============================
            // 2D TProfiles (Lambda correlations)
            // ===============================
            histos.add((folder + "/p2dPxStarDeltaPhiVsLambdaPt").c_str(), "p2dPxStarDeltaPhiVsLambdaPt;#Delta#varphi_{jet};#it{p}_{T}^{#Lambda};<P_{#Lambda}>_{x}", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPtSigExtract});
            histos.add((folder + "/p2dPyStarDeltaPhiVsLambdaPt").c_str(), "p2dPyStarDeltaPhiVsLambdaPt;#Delta#varphi_{jet};#it{p}_{T}^{#Lambda};<P_{#Lambda}>_{y}", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPtSigExtract});
            histos.add((folder + "/p2dPzStarDeltaPhiVsLambdaPt").c_str(), "p2dPzStarDeltaPhiVsLambdaPt;#Delta#varphi_{jet};#it{p}_{T}^{#Lambda};<P_{#Lambda}>_{z}", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPtSigExtract});

            // ===============================
            // QA histograms - Useful numbers
            // ===============================
            // (TODO: implement these!)
            // (TODO: implement momentum imbalance checks for jets!)
            // Added to a separate folder for further control (changed the usage of the "folder" string):
            // histos.add(("QA_Numbers/" + folder + "/hEventsWithV0").c_str(), "hEventsWithV0", kTH1D, {{1,0,1}}); // In the current derived data, all saved events have a V0
            // histos.add(("QA_Numbers/" + folder + "/hLambdaCounter").c_str(), "hLambdaCounter", kTH1D, {{1,0,1}});
            // histos.add(("QA_Numbers/" + folder + "/hAntiLambdaCounter").c_str(), "hAntiLambdaCounter", kTH1D, {{1,0,1}});
            // histos.add(("QA_Numbers/" + folder + "/hValidLeadJets").c_str(), "hValidLeadJets", kTH1D, {{1,0,1}});
        };
        // Execute local lambda to register histogram families:
        addRingObservableFamily("Ring");
        addRingObservableFamily("RingKinematicCuts");
        addRingObservableFamily("JetKinematicCuts");
        addRingObservableFamily("JetAndLambdaKinematicCuts");
    }

    // Initializing a random number generator for the worker (for perpendicular-to-jet direction QAs):
    TRandom3 randomGen{0}; // 0 means we auto-seed from machine entropy. This is called once per device in the pipeline, so we should not see repeated seeds across workers

        // Preslices for correct collisions association:
        // (TODO: test using custom grouping)
    Preslice<aod::RingJets> perColJets = o2::aod::lambdajetpol::collisionId;
    Preslice<aod::RingLaV0s> perColV0s = o2::aod::lambdajetpol::collisionId;
    Preslice<aod::RingLeadP> perColLeadPs = o2::aod::lambdajetpol::collisionId;
    void processPolarizationData(o2::aod::RingCollisions const& collisions, o2::aod::RingJets const& jets, o2::aod::RingLaV0s const& v0s,
                                 o2::aod::RingLeadP const& leadPs){
        
        for (auto const& collision : collisions) {
            const auto collId = collision.collisionId();
            const double centrality = collision.centrality();

            // Slice jets and V0s belonging to this collision
                // (global collision indices repeat a lot, but they are unique to a same TimeFrame (TF) subfolder in the derived data)
            auto jetsInColl = jets.sliceBy(perColJets, collId);
            auto v0sInColl  = v0s.sliceBy(perColV0s, collId);

            // Check if there is at least one V0 and one jet in the collision:
            // (in the way I fill the table, there is always at least one V0 in
            //  the stored collision, but the jets table can not be filled for
            //  that collision, and a collision may not be filled when the jets
            //  table is. Be mindful of that!)
            if (!jetsInColl.size() || !v0sInColl.size()) continue;
            // (TODO: either add a separate loop for events that have a leading particle, but no leading jet,
            // or be aware that there are events saved with a leading particle that don't have a single valid jet in them)
            // (by the way, no need to check a leadPsInColl, because if there is a jet in the collision, 
            //  there surely is a leading particle in it)

            // Get leading jet and second to leading jet:
            double leadingJetPt = -1.;
            double subleadingJetPt = -1.;
            o2::aod::RingJets::iterator leadingJet;
            std::optional<o2::aod::RingJets::iterator> subleadingJet;
            // std::optional avoids undefined behaviour from a default-constructed iterator:
            // (will work if subleadingJet was not found!)
            for (auto const& jet : jetsInColl) {
                const auto jetpt = jet.jetPt();
                if (jetpt > leadingJetPt){
                    // Current leading becomes subleading:
                    subleadingJetPt = leadingJetPt;
                    subleadingJet = leadingJet; // may still be std::nullopt on first pass -- that is handled by std::optional!
                    // Now update updating the leading jet:
                    leadingJetPt = jetpt;
                    leadingJet = jet;
                }
                else if (jetpt > subleadingJetPt){ // Update subleading only:
                    subleadingJetPt = jetpt;
                    subleadingJet = jet;
                }
            }

            // For leading jet (always exists):
            const double leadingJetEta = leadingJet.jetEta();
            const double leadingJetPhi = leadingJet.jetPhi();
            // Convert to 3-vector components for inner product:
            const double jetPx = leadingJetPt * std::cos(leadingJetPhi);
            const double jetPy = leadingJetPt * std::sin(leadingJetPhi);
            const double jetPz = leadingJetPt * std::sinh(leadingJetEta);
            // XYZVector leadingJetVec(jetPx, jetPy, jetPz);
            XYZVector leadingJetUnitVec = XYZVector(jetPx, jetPy, jetPz).Unit();

            // QA block -- Purposefully changing the jet direction (should kill signal, if any):
            if (forcePerpToJet) { // Use modified jet direction (done outside loop to guarantee all V0s inside event use same fake jet)
                // First, we build a vector perpendicular to the jet by picking an arbitrary vector not parallel to the jet
                XYZVector refVec(1., 0., 0.);
                if (std::abs(leadingJetUnitVec.Dot(refVec)) > 0.99) refVec = XYZVector(0., 1., 0.);
                // Now we get a perpendicular vector to the jet direction:
                XYZVector perpVec = leadingJetUnitVec.Cross(refVec).Unit();
                
                // Now we rotate around the jet axis by a random angle, just to make sure we are not introducing a bias in the QA:
                // We will use Rodrigues' rotation formula (v_rot = v*cos(randomAngle) + (Jet \cross v)*sin(randomAngle))
                double randomAngle = randomGen.Uniform(0., o2::constants::math::TwoPI);
                XYZVector rotatedPerpVec = perpVec * std::cos(randomAngle) + leadingJetUnitVec.Cross(perpVec) * std::sin(randomAngle);
                leadingJetUnitVec = rotatedPerpVec;
            }

            // -----------------------------------------------------------------------
            // Find the leading particle for this collision.
            // pT = -1 means "not found".
            double leadPPt  = -1.;
            double leadPEta =  0.; // safe dummy values -- only used when leadPPt > 0
            double leadPPhi =  0.;
            // std::optional avoids the unassigned-iterator trap again:
            std::optional<o2::aod::RingLeadP::iterator> leadingParticleOpt;
            auto leadPsInColl = leadPs.sliceBy(perColLeadPs, collId);
            for (auto const& lp : leadPsInColl) {
                // Table should contain exactly one entry per collision,
                // but we break immediately to be safe:
                leadingParticleOpt = lp;
                break;
            }
            // Extract kinematics only if we actually found an entry:
            // (Physically, if there is at least one jet there should always be a
            //  leading particle, but we guard against it anyway)
            if (leadingParticleOpt.has_value()) {
                leadPPt  = leadingParticleOpt->leadParticlePt();
                leadPEta = leadingParticleOpt->leadParticleEta();
                leadPPhi = leadingParticleOpt->leadParticlePhi();
            }

            // Defining some bools to help:
            bool hasValidLeadP = leadPPt > 0.;
            bool hasValidSubJet = subleadingJetPt > 0.;

            // --- Subleading jet (only valid when subleadingJetPt > 0) ---
            // We still build the variables here to keep the V0 loop clean.
            // Their values are irrelevant when subleadingJetPt <= 0.
            double subleadingJetEta = 0.;
            double subleadingJetPhi = 0.;
            XYZVector subJetUnitVec(1., 0., 0.); // dummy unit vector: overwritten below
            if (hasValidSubJet) {
                subleadingJetEta = subleadingJet->jetEta();
                subleadingJetPhi = subleadingJet->jetPhi();
                const double subJetPx = subleadingJetPt * std::cos(subleadingJetPhi);
                const double subJetPy = subleadingJetPt * std::sin(subleadingJetPhi);
                const double subJetPz = subleadingJetPt * std::sinh(subleadingJetEta);
                subJetUnitVec = XYZVector(subJetPx, subJetPy, subJetPz).Unit();
            }

            // --- Leading particle (only valid when leadPPt > 0) ---
            XYZVector leadPUnitVec(1., 0., 0.); // dummy: overwritten below
            if (hasValidLeadP) {
                const double leadPPx = leadPPt * std::cos(leadPPhi);
                const double leadPPy = leadPPt * std::sin(leadPPhi);
                const double leadPPz = leadPPt * std::sinh(leadPEta);
                leadPUnitVec = XYZVector(leadPPx, leadPPy, leadPPz).Unit();
            }

            // Calculating per-event bools only once:
            const bool kinematicJetCheck = std::abs(leadingJetEta) < 0.5;
            const bool kinematic2ndJetCheck = (subleadingJetPt > 0.) && (std::abs(subleadingJetEta) < 0.5);
            const bool kinematicLeadPCheck = (leadPPt > 0.) && (std::abs(leadPEta) < 0.5);

            // TODO: add centrality selection procedure and options (one configurable for no centrality separation at all too!)
            for (auto const& v0 : v0sInColl) {
                const bool isLambda = v0.isLambda();
                const bool isAntiLambda = v0.isAntiLambda();
                // For now, removing the ambiguous candidates from the analysis. Derived data permits handling both.
                // (From Podolanski-Armenteros plots, the population of ambiguous is ~2% without TOF, and without 
                //  competing mass rejection. From those, ~99% seem to be K0s, so no real gain in considering the
                //  ambiguous candidates in the analysis)
                if (isLambda && isAntiLambda) continue;
                const double v0pt = v0.v0Pt();
                const double v0eta = v0.v0Eta();
                const double v0phi = v0.v0Phi();

                double v0LambdaLikeMass = 0; // Initialized just to catch any stray behavior
                double protonLikePt = 0;
                double protonLikeEta = 0;
                double protonLikePhi = 0;
                if (isLambda){
                    if (!analyseLambda) continue;
                    v0LambdaLikeMass = v0.massLambda();
                    protonLikePt = v0.posPt();
                    protonLikeEta = v0.posEta();
                    protonLikePhi = v0.posPhi();
                }
                else if (isAntiLambda){ // (TODO: add a split histogram where you consider Lambda and AntiLambda polarization separately?)
                    if (!analyseAntiLambda) continue;
                    v0LambdaLikeMass = v0.massAntiLambda();
                    protonLikePt = v0.negPt();
                    protonLikeEta = v0.negEta();
                    protonLikePhi = v0.negPhi();
                }
                
                PtEtaPhiMVector lambdaLike4Vec(v0pt, v0eta, v0phi, v0LambdaLikeMass);
                PtEtaPhiMVector protonLike4Vec(protonLikePt, protonLikeEta, protonLikePhi, protonMass);
                double lambdaRapidity = lambdaLike4Vec.Rapidity(); // For further kinematic selections

                // Boosting proton into lambda frame:
                XYZVector beta = -lambdaLike4Vec.BoostToCM(); // Boost trivector that goes from laboratory frame to the rest frame
                auto protonLike4VecStar = ROOT::Math::VectorUtil::boost(protonLike4Vec, beta);

                // Getting unit vectors and 3-components:
                XYZVector lambdaLike3Vec = lambdaLike4Vec.Vect();
                XYZVector protonLikeStarUnit3Vec = protonLike4VecStar.Vect().Unit();

                // Calculating cross product:
                XYZVector cross = leadingJetUnitVec.Cross(lambdaLike3Vec);
                double crossProductNorm = cross.R();

                double ringObservable = protonLikeStarUnit3Vec.Dot(cross) / crossProductNorm;
                    // Adding the prefactor related to the CP-violating decay (decay constants have different signs)
                if (!forcePolSignQA) ringObservable *= (isLambda) ? polPrefactorLambda : polPrefactorAntiLambda;
                else ringObservable *= (isLambda) ? polPrefactorLambda : -1.0*polPrefactorAntiLambda;

                // Angular variables:
                double deltaPhiJet = wrapToPiFast(v0phi - leadingJetPhi); // Wrapped to [-PI, pi), for convenience
                double deltaThetaJet = ROOT::Math::VectorUtil::Angle(leadingJetUnitVec, lambdaLike3Vec); // 3D angular separation

                //////////////////////////////////////////
                // Ring observable: Subleading jet proxy
                //////////////////////////////////////////
                double ringObservable2ndJet = 0.;
                double deltaPhi2ndJet = 0.;
                double deltaTheta2ndJet = 0.;
                if (hasValidSubJet) {
                    // Cross product
                    XYZVector cross2ndJet = subJetUnitVec.Cross(lambdaLike3Vec);
                    double crossProductNorm2ndJet = cross2ndJet.R();
                    double ringObservable2ndJet = protonLikeStarUnit3Vec.Dot(cross2ndJet) / crossProductNorm2ndJet;
                    // Adding prefactor
                    if (!forcePolSignQA) ringObservable2ndJet *= (isLambda) ? polPrefactorLambda : polPrefactorAntiLambda;
                    else ringObservable2ndJet *= (isLambda) ? polPrefactorLambda : -1.0 * polPrefactorAntiLambda;
                    // Angular variables
                    double deltaPhi2ndJet = wrapToPiFast(v0phi - subleadingJetPhi);
                    double deltaTheta2ndJet = ROOT::Math::VectorUtil::Angle(subJetUnitVec, lambdaLike3Vec);
                }

                ////////////////////////////////////////////
                // Ring observable: Leading particle proxy
                ////////////////////////////////////////////
                double ringObservableLeadP = 0.;
                double deltaPhiLeadP = 0.;
                double deltaThetaLeadP = 0.;
                if (hasValidLeadP) {
                    // Cross product
                    XYZVector crossLeadP = leadPUnitVec.Cross(lambdaLike3Vec);
                    double crossProductNormLeadP = crossLeadP.R();
                    double ringObservableLeadP = protonLikeStarUnit3Vec.Dot(crossLeadP) / crossProductNormLeadP;
                    // Adding prefactor
                    if (!forcePolSignQA) ringObservableLeadP *= (isLambda) ? polPrefactorLambda : polPrefactorAntiLambda;
                    else ringObservableLeadP *= (isLambda) ? polPrefactorLambda : -1.0 * polPrefactorAntiLambda;
                    // Angular variables
                    double deltaPhiLeadP = wrapToPiFast(v0phi - leadPPhi);
                    double deltaThetaLeadP = ROOT::Math::VectorUtil::Angle(leadPUnitVec, lambdaLike3Vec);
                }

                // Calculating polarization observables (in the Lambda frame, because that is easier -- does not require boosts):
                    // To be precise, not actually the polarization, but a part of the summand in P^*_\Lambda = (3/\alpha_\Lambda) * <p^*_{proton}>
                double PolStarX, PolStarY, PolStarZ;
                if (isLambda){ // Notice there is no need to check analyseLambda again due to previous checks.
                    PolStarX = polPrefactorLambda * protonLikeStarUnit3Vec.X();
                    PolStarY = polPrefactorLambda * protonLikeStarUnit3Vec.Y();
                    PolStarZ = polPrefactorLambda * protonLikeStarUnit3Vec.Z();
                }
                else if (isAntiLambda){
                    PolStarX = polPrefactorAntiLambda * protonLikeStarUnit3Vec.X();
                    PolStarY = polPrefactorAntiLambda * protonLikeStarUnit3Vec.Y();
                    PolStarZ = polPrefactorAntiLambda * protonLikeStarUnit3Vec.Z();
                }

                double v0phiToFillHists = wrapToPiFast(v0phi); // A short wrap to reuse some predefined axes

                // Fill ring histograms: (1D, lambda 2D correlations and jet 2D correlations):
                RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "Ring") // Notice the usage of macros! If you change the variable names, this WILL break the code!
                                                                    // No, there should NOT be any ";" here! Read the macro definition for an explanation
                if (hasValidLeadP) {RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "Ring")}
                if (hasValidSubJet) {RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "Ring")}
                POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "Ring")

                // Extra kinematic criteria for Lambda candidates (removes polarization background):
                const bool kinematicLambdaCheck = (v0pt > 0.5 && v0pt < 1.5) && std::abs(lambdaRapidity) < 0.5;
                if (kinematicLambdaCheck){
                    RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
                    POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
                    if (hasValidLeadP) {RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")}
                    if (hasValidSubJet) {RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")}
                }
                
                // Extra selection criteria on jet candidates:
                if (kinematicJetCheck){ // This is redundant for jets with R=0.4, but for jets with R<0.4 the leading jet may be farther in eta.
                    RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
                    POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
                }
                
                // Extra selection criteria on both Lambda and jet candidates:
                if (kinematicLambdaCheck && kinematicJetCheck){
                    RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
                    POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
                }

                // Same variations for the leading particle and for the subleading jet:
                if (kinematicLeadPCheck){RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")}
                if (kinematic2ndJetCheck){RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")}
                if (kinematicLambdaCheck && kinematicLeadPCheck){RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")}
                if (kinematicLambdaCheck && kinematic2ndJetCheck){RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")}
            } // end v0s loop
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