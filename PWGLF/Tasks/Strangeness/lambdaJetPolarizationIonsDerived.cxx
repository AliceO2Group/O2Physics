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

enum CentEstimator {
    kCentFT0C = 0,
    kCentFT0M,
    kCentFT0CVariant1,
    kCentMFT,
    kCentNGlobal,
    kCentFV0A
};

// Helper macro to avoid writing the histogram fills 4 times for about 20 histograms:
#define RING_OBSERVABLE_FILL_LIST(X, FOLDER)                                                              \
    /* Counters */                                                                                        \
    X(FOLDER "/QA/hDeltaPhi",                            deltaPhiJet)                                     \
    X(FOLDER "/QA/hDeltaTheta",                          deltaThetaJet)                                   \
    X(FOLDER "/QA/hIntegrated",                          0.)                                              \
    /* Lambda pT variation -- Youpeng's proposal */                                                       \
    X(FOLDER "/QA/hLambdaPt",                            v0pt)                                            \
    /* Counters */                                                                                        \
    X(FOLDER "/QA/h2dDeltaPhiVsLambdaPt",                deltaPhiJet,   v0pt)                             \
    X(FOLDER "/QA/h2dDeltaThetaVsLambdaPt",              deltaThetaJet, v0pt)                             \
    /* Additional plots for instant gratification - 1D Profiles */                                        \
    X(FOLDER "/pRingObservableDeltaPhi",                 deltaPhiJet,   ringObservable)                   \
    X(FOLDER "/pRingObservableDeltaTheta",               deltaThetaJet, ringObservable)                   \
    X(FOLDER "/pRingObservableIntegrated",               0.,            ringObservable)                   \
    X(FOLDER "/pRingObservableLambdaPt",                 v0pt,          ringObservable)                   \
    /* 2D Profiles */                                                                                     \
    X(FOLDER "/p2dRingObservableDeltaPhiVsLambdaPt",     deltaPhiJet,   v0pt, ringObservable)             \
    X(FOLDER "/p2dRingObservableDeltaThetaVsLambdaPt",   deltaThetaJet, v0pt, ringObservable)             \
    X(FOLDER "/p2dRingObservableDeltaPhiVsLeadJetPt",    deltaPhiJet,   leadingJetPt, ringObservable)     \
    X(FOLDER "/p2dRingObservableDeltaThetaVsLeadJetPt",  deltaThetaJet, leadingJetPt, ringObservable)     \
    /* 1D Mass */                                                                                         \
    X(FOLDER "/QA/hMass",                                v0LambdaLikeMass)                                \
    X(FOLDER "/QA/hRingObservableNumMass",               v0LambdaLikeMass, ringObservable)                \
    X(FOLDER "/hMassSigExtract",                         v0LambdaLikeMass)                                \
    /* Counters */                                                                                        \
    X(FOLDER "/QA/h2dDeltaPhiVsMass",                    deltaPhiJet,   v0LambdaLikeMass)                 \
    X(FOLDER "/QA/h2dDeltaThetaVsMass",                  deltaThetaJet, v0LambdaLikeMass)                 \
    X(FOLDER "/QA/h3dDeltaPhiVsMassVsLambdaPt",          deltaPhiJet,   v0LambdaLikeMass, v0pt)           \
    X(FOLDER "/QA/h3dDeltaThetaVsMassVsLambdaPt",        deltaThetaJet, v0LambdaLikeMass, v0pt)           \
    X(FOLDER "/QA/h3dDeltaPhiVsMassVsLeadJetPt",         deltaPhiJet,   v0LambdaLikeMass, leadingJetPt)   \
    X(FOLDER "/QA/h3dDeltaThetaVsMassVsLeadJetPt",       deltaThetaJet, v0LambdaLikeMass, leadingJetPt)   \
    X(FOLDER "/QA/h3dDeltaPhiVsMassVsCent",              deltaPhiJet,   v0LambdaLikeMass, centrality)     \
    X(FOLDER "/QA/h3dDeltaThetaVsMassVsCent",            deltaThetaJet, v0LambdaLikeMass, centrality)     \
    /* TProfile of Ring vs Mass */                                                                        \
    X(FOLDER "/pRingObservableMass",                     v0LambdaLikeMass, ringObservable)                \
    /* TProfile of Ring vs Mass -- Leading Particle and 2nd-to-leading jet - QA */                        \
    X(FOLDER "/pRingObservableLeadPMass",                v0LambdaLikeMass, ringObservableLeadP)           \
    X(FOLDER "/pRingObservable2ndJetMass",               v0LambdaLikeMass, ringObservable2ndJet)          \
    /* 2D Profiles: Angle vs Mass */                                                                      \
    X(FOLDER "/p2dRingObservableDeltaPhiVsMass",         deltaPhiJet,   v0LambdaLikeMass, ringObservable) \
    X(FOLDER "/p2dRingObservableDeltaThetaVsMass",       deltaThetaJet, v0LambdaLikeMass, ringObservable) \
    /* 3D Profiles: Angle vs Mass vs Lambda pT */                                                         \
    X(FOLDER "/p3dRingObservableDeltaPhiVsMassVsLambdaPt",   deltaPhiJet,   v0LambdaLikeMass, v0pt, ringObservable) \
    X(FOLDER "/p3dRingObservableDeltaThetaVsMassVsLambdaPt", deltaThetaJet, v0LambdaLikeMass, v0pt, ringObservable) \
    /* 3D Profiles: Angle vs Mass vs Lead Jet pT */                                                       \
    X(FOLDER "/p3dRingObservableDeltaPhiVsMassVsLeadJetPt",  deltaPhiJet,   v0LambdaLikeMass, leadingJetPt, ringObservable) \
    X(FOLDER "/p3dRingObservableDeltaThetaVsMassVsLeadJetPt",deltaThetaJet, v0LambdaLikeMass, leadingJetPt, ringObservable) \
    /* 2D Profile: Mass vs Centrality */                                                                  \
    X(FOLDER "/p2dRingObservableMassVsCent",             v0LambdaLikeMass, centrality, ringObservable)    \
    /* 3D Profiles: Angle vs Mass vs Centrality */                                                        \
    X(FOLDER "/p3dRingObservableDeltaPhiVsMassVsCent",   deltaPhiJet,   v0LambdaLikeMass, centrality, ringObservable) \
    X(FOLDER "/p3dRingObservableDeltaThetaVsMassVsCent", deltaThetaJet, v0LambdaLikeMass, centrality, ringObservable) \
    X(FOLDER "/pRingIntVsCentrality",                    centrality, ringObservable)
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
    X(FOLDER "/QA/pPxStarPhi",                   v0phiToFillHists,           PolStarX) \
    X(FOLDER "/QA/pPyStarPhi",                   v0phiToFillHists,           PolStarY) \
    X(FOLDER "/QA/pPzStarPhi",                   v0phiToFillHists,           PolStarZ) \
    /* =============================== */ \
    /* 1D TProfiles vs DeltaPhi_jet */ \
    /* =============================== */ \
    X(FOLDER "/QA/pPxStarDeltaPhi",              deltaPhiJet,   PolStarX) \
    X(FOLDER "/QA/pPyStarDeltaPhi",              deltaPhiJet,   PolStarY) \
    X(FOLDER "/QA/pPzStarDeltaPhi",              deltaPhiJet,   PolStarZ) \
    /* =============================== */ \
    /* 2D TProfiles vs DeltaPhi_jet and Lambda pT */ \
    /* =============================== */ \
    X(FOLDER "/QA/p2dPxStarDeltaPhiVsLambdaPt",  deltaPhiJet, v0pt, PolStarX) \
    X(FOLDER "/QA/p2dPyStarDeltaPhiVsLambdaPt",  deltaPhiJet, v0pt, PolStarY) \
    X(FOLDER "/QA/p2dPzStarDeltaPhiVsLambdaPt",  deltaPhiJet, v0pt, PolStarZ)

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

    // Centrality:
    Configurable<int> centralityEstimator{"centralityEstimator", kCentFT0M, "Run 3 centrality estimator (0:CentFT0C, 1:CentFT0M, 2:CentFT0CVariant1, 3:CentMFT, 4:CentNGlobal, 5:CentFV0A)"}; // Default is FT0M
    
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
            // QA histograms: angle and pT distributions
            // (No mass dependency -- useful to check kinematic sculpting from cuts)
            // ===============================
            histos.add((folder + "/QA/hDeltaPhi").c_str(), "hDeltaPhi", kTH1D, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/QA/hDeltaTheta").c_str(), "hDeltaTheta", kTH1D, {axisConfigurations.axisDeltaTheta});
            histos.add((folder + "/QA/hIntegrated").c_str(), "hIntegrated", kTH1D, {{1, -0.5, 0.5}});
            // ===============================
            // Lambda pT dependence
            // ===============================
            histos.add((folder + "/QA/hLambdaPt").c_str(), "hLambdaPt", kTH1D, {axisConfigurations.axisPt});
            histos.add((folder + "/QA/h2dDeltaPhiVsLambdaPt").c_str(), "h2dDeltaPhiVsLambdaPt", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPt});
            histos.add((folder + "/QA/h2dDeltaThetaVsLambdaPt").c_str(), "h2dDeltaThetaVsLambdaPt", kTH2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisPt});
            // ===============================
            //   Polarization observable QAs 
            // (not Ring: actual polarization!)
            // ===============================
                // Will implement these as TProfiles, as polarization is also a measure like P_\Lambda = (3/\alpha_\Lambda) * <p_{proton}>, so the error is similar
            // ===============================
            // 1D TProfiles
            // ===============================
            histos.add((folder + "/QA/pPxStarPhi").c_str(), "pPxStarPhi;#varphi_{#Lambda};<P_{#Lambda}>_{x}", kTProfile, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/QA/pPyStarPhi").c_str(), "pPyStarPhi;#varphi_{#Lambda};<P_{#Lambda}>_{y}", kTProfile, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/QA/pPzStarPhi").c_str(), "pPzStarPhi;#varphi_{#Lambda};<P_{#Lambda}>_{z}", kTProfile, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/QA/pPxStarDeltaPhi").c_str(), "pPxStarDeltaPhi;#Delta#varphi_{jet};<P_{#Lambda}>_{x}", kTProfile, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/QA/pPyStarDeltaPhi").c_str(), "pPyStarDeltaPhi;#Delta#varphi_{jet};<P_{#Lambda}>_{y}", kTProfile, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/QA/pPzStarDeltaPhi").c_str(), "pPzStarDeltaPhi;#Delta#varphi_{jet};<P_{#Lambda}>_{z}", kTProfile, {axisConfigurations.axisDeltaPhi});
            // ===============================
            // 2D TProfiles (Lambda correlations)
            // ===============================
            histos.add((folder + "/QA/p2dPxStarDeltaPhiVsLambdaPt").c_str(), "p2dPxStarDeltaPhiVsLambdaPt;#Delta#varphi_{jet};#it{p}_{T}^{#Lambda};<P_{#Lambda}>_{x}", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPtSigExtract});
            histos.add((folder + "/QA/p2dPyStarDeltaPhiVsLambdaPt").c_str(), "p2dPyStarDeltaPhiVsLambdaPt;#Delta#varphi_{jet};#it{p}_{T}^{#Lambda};<P_{#Lambda}>_{y}", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPtSigExtract});
            histos.add((folder + "/QA/p2dPzStarDeltaPhiVsLambdaPt").c_str(), "p2dPzStarDeltaPhiVsLambdaPt;#Delta#varphi_{jet};#it{p}_{T}^{#Lambda};<P_{#Lambda}>_{z}", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPtSigExtract});
            
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
            // Simple invariant mass plot for QA:
            histos.add((folder + "/QA/hMass").c_str(), "hMass", kTH1D, {axisConfigurations.axisLambdaMass});
            histos.add((folder + "/hMassSigExtract").c_str(), "hMassSigExtract", kTH1D, {axisConfigurations.axisLambdaMassSigExtract});
            // 1D Mass dependence of observable numerator:
            histos.add((folder + "/QA/hRingObservableNumMass").c_str(), "hRingObservableNumMass", kTH1D, {axisConfigurations.axisLambdaMassSigExtract});
            // --- 2D counters: Angle vs Mass vs ---
            histos.add((folder + "/QA/h2dDeltaPhiVsMass").c_str(), "h2dDeltaPhiVsMass", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract});
            histos.add((folder + "/QA/h2dDeltaThetaVsMass").c_str(), "h2dDeltaThetaVsMass", kTH2D,{axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract});
            // --- 3D counters: Angle vs Mass vs Lambda pT ---
            histos.add((folder + "/QA/h3dDeltaPhiVsMassVsLambdaPt").c_str(), "h3dDeltaPhiVsMassVsLambdaPt", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
            histos.add((folder + "/QA/h3dDeltaThetaVsMassVsLambdaPt").c_str(), "h3dDeltaThetaVsMassVsLambdaPt", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
            // --- 3D counters: Angle vs Mass vs Lead Jet pT ---
            histos.add((folder + "/QA/h3dDeltaPhiVsMassVsLeadJetPt").c_str(), "h3dDeltaPhiVsMassVsLeadJetPt", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract}); 
            histos.add((folder + "/QA/h3dDeltaThetaVsMassVsLeadJetPt").c_str(), "h3dDeltaThetaVsMassVsLeadJetPt", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract});

            // ===============================
            // TProfiles vs Mass: quick glancing before signal extraction
            // ===============================
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
            // Counters
            histos.add((folder + "/QA/h3dDeltaPhiVsMassVsCent").c_str(), "h3dDeltaPhiVsMassVsCent", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
            histos.add((folder + "/QA/h3dDeltaThetaVsMassVsCent").c_str(), "h3dDeltaThetaVsMassVsCent", kTH3D,{axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
            // Useful TProfiles:
            // --- TProfile1D: Integrated <R> vs Centrality:
            histos.add((folder + "/pRingIntVsCentrality").c_str(), "pRingIntVsCentrality; Centrality (%);<#it{R}>", kTProfile, {axisConfigurations.axisCentrality});
            // --- TProfile2D: <R> vs Mass vs Centrality ---
            histos.add((folder + "/p2dRingObservableMassVsCent").c_str(), "p2dRingObservableMassVsCent;m_{p#pi};Centrality;<#it{R}>", kTProfile2D, {axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
            // --- TProfile3D: <R> vs DeltaPhi vs Mass vs Centrality ---
            histos.add((folder + "/p3dRingObservableDeltaPhiVsMassVsCent").c_str(), "p3dRingObservableDeltaPhiVsMassVsCent;#Delta#varphi;m_{p#pi};Centrality;<#it{R}>", kTProfile3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
            // --- TProfile3D: <R> vs DeltaTheta vs Mass vs Centrality ---
            histos.add((folder + "/p3dRingObservableDeltaThetaVsMassVsCent").c_str(), "p3dRingObservableDeltaThetaVsMassVsCent;#Delta#theta;m_{p#pi};Centrality;<#it{R}>", kTProfile3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});

            // ===============================
            // QA histograms - Useful numbers
            // ===============================
            // (TODO: implement these!)
            // (TODO: implement momentum imbalance checks for jets!)
            // Added to a separate folder for further control (changed the usage of the "folder" string):
            // histos.add(("QA_Numbers/" + folder + "/hValidLeadJets").c_str(), "hValidLeadJets", kTH1D, {{1,0,1}});
            // TODO: Add "frequency of jets per pT" histograms either here or in the TableProducer
        };
        // Execute local lambda to register histogram families:
        addRingObservableFamily("Ring");
        addRingObservableFamily("RingKinematicCuts");
        addRingObservableFamily("JetKinematicCuts");
        addRingObservableFamily("JetAndLambdaKinematicCuts");

        histos.add("pRingCuts", "pRingCuts; ;<#it{R}>", kTProfile, {{4, 0, 4}});
        histos.get<TProfile>(HIST("pRingCuts"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
        histos.get<TProfile>(HIST("pRingCuts"))->GetXaxis()->SetBinLabel(2, "p_{T}^{#Lambda}@[0.5,1.5],|y_{#Lambda}|<0.5"); // (v0pt > 0.5 && v0pt < 1.5) && std::abs(lambdaRapidity) < 0.5;
        histos.get<TProfile>(HIST("pRingCuts"))->GetXaxis()->SetBinLabel(3, "|Jet_{#eta}|<0.5");
        histos.get<TProfile>(HIST("pRingCuts"))->GetXaxis()->SetBinLabel(4, "#Lambda + Jet cuts");

        // Same for subleading jet and leading particle:
        histos.add("pRingCutsSubLeadingJet", "pRingCutsSubLeadingJet; ;<#it{R}>", kTProfile, {{4, 0, 4}});
        histos.get<TProfile>(HIST("pRingCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
        histos.get<TProfile>(HIST("pRingCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(2, "p_{T,#Lambda}@[0.5,1.5],|y_{#Lambda}|<0.5");
        histos.get<TProfile>(HIST("pRingCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(3, "|SubJet_{#eta}|<0.5");
        histos.get<TProfile>(HIST("pRingCutsSubLeadingJet"))->GetXaxis()->SetBinLabel(4, "#Lambda + SubJet cuts");

        histos.add("pRingCutsLeadingP", "pRingCutsLeadingP; ;<#it{R}>", kTProfile, {{4, 0, 4}});
        histos.get<TProfile>(HIST("pRingCutsLeadingP"))->GetXaxis()->SetBinLabel(1, "All #Lambda");
        histos.get<TProfile>(HIST("pRingCutsLeadingP"))->GetXaxis()->SetBinLabel(2, "p_{T}^{#Lambda}@[0.5,1.5],|y_{#Lambda}|<0.5");
        histos.get<TProfile>(HIST("pRingCutsLeadingP"))->GetXaxis()->SetBinLabel(3, "|LeadP_{#eta}|<0.5");
        histos.get<TProfile>(HIST("pRingCutsLeadingP"))->GetXaxis()->SetBinLabel(4, "#Lambda + LeadP cuts");
    }

    // Helper to get centrality (same from TableProducer, thanks to templating!):
    template <typename TCollision>
    auto getCentrality(TCollision const& collision)
    {
        if (centralityEstimator == kCentFT0M) return collision.centFT0M();
        else if (centralityEstimator == kCentFT0C) return collision.centFT0C();
        else if (centralityEstimator == kCentFT0CVariant1) return collision.centFT0CVariant1();
        else if (centralityEstimator == kCentMFT) return collision.centMFT();
        else if (centralityEstimator == kCentNGlobal) return collision.centNGlobal();
        else if (centralityEstimator == kCentFV0A) return collision.centFV0A();
        return -1.f;
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
            const double centrality = getCentrality(collision);

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
            if (!v0sInColl.size()) continue;

            // 2) We require at least a leading particle, then we get the leading jet only if it exists:
            // (The goal is to see how diluted the signal gets with events which don't even have a loose FastJet jet)
            // (The leading particle is built from all tracks that passed the pseudojet
            // selection, so it exists whenever FastJet was run on this collision.
            // Events that have a leading jet always have a leading particle too, but
            // the converse is not true: events can have a leading particle with no jet
            // if no jet survives the pT threshold/the background subtraction)
            float leadPPt  = -1.; // pT = -1 means "table entry not found for this collision".
            float leadPEta =  0.;
            float leadPPhi =  0.;
            float leadPPx = 0., leadPPy = 0., leadPPz = 0.;
            for (auto const& lp : leadPsInColl) {
                // Table should contain exactly one entry per collision,
                // but we break immediately to be safe:
                leadPPt  = lp.leadParticlePt();
                leadPEta = lp.leadParticleEta();
                leadPPhi = lp.leadParticlePhi();
                // Using dynamic columns to make code cleaner:
                leadPPx = lp.leadParticlePx();
                leadPPy = lp.leadParticlePy();
                leadPPz = lp.leadParticlePz();
                break;
            }
            // Discard events with no leading particle (FastJet didn't even run in these cases!):
            if (leadPPt < 0.) continue;

            // Build leading particle unit vector, outside the V0 loop for performance:
            XYZVector leadPUnitVec = XYZVector(leadPPx, leadPPy, leadPPz).Unit();

            // 3) Checking if the event has a leading jet:
            auto jetsInColl = jets.sliceBy(perColJets, collId);
            float leadingJetPt    = -1.;
            float subleadingJetPt = -1.;
            // std::optional avoids undefined behaviour from a default-constructed iterator:
            // (essentially, just protection for when we fetch jetEta() and the such)
            std::optional<o2::aod::RingJets::iterator> leadingJet;
            std::optional<o2::aod::RingJets::iterator> subleadingJet;
            for (auto const& jet : jetsInColl) {
                const auto jetpt = jet.jetPt();
                if (jetpt > leadingJetPt){
                    // Current leading becomes subleading:
                    subleadingJetPt = leadingJetPt;
                    subleadingJet = leadingJet; // may still be std::nullopt on first pass -- that is fine!
                    // Now update the leading jet:
                    leadingJetPt = jetpt;
                    leadingJet = jet;
                }
                else if (jetpt > subleadingJetPt){ // Update subleading only:
                    subleadingJetPt = jetpt;
                    subleadingJet = jet;
                }
            }

            // Some useful bools to check if we have a leading jet and a subleading jet:
            const bool hasValidLeadingJet = leadingJetPt > 0.;
            const bool hasValidSubJet = subleadingJetPt > 0.;

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
                    leadingJetUnitVec = perpVec * std::cos(randomAngle) + leadingJetUnitVec.Cross(perpVec) * std::sin(randomAngle);
                }
            }

            float subleadingJetEta = 0.;
            float subleadingJetPhi = 0.;
            XYZVector subJetUnitVec(1., 0., 0.);
            if (hasValidSubJet) {
                subleadingJetEta = subleadingJet->jetEta();
                subleadingJetPhi = subleadingJet->jetPhi();
                // Using internal getters to make code cleaner:
                subJetUnitVec = XYZVector(subleadingJet->jetPx(), subleadingJet->jetPy(), subleadingJet->jetPz()).Unit();
            }

            // (jet eta cuts only meaningful when the jet actually exists)
            const bool kinematicJetCheck = hasValidLeadingJet && (std::abs(leadingJetEta) < 0.5);
            const bool kinematic2ndJetCheck = hasValidSubJet && (std::abs(subleadingJetEta) < 0.5);
            const bool kinematicLeadPCheck = std::abs(leadPEta) < 0.5;

            for (auto const& v0 : v0sInColl) {
                const bool isLambda = v0.isLambda();
                const bool isAntiLambda = v0.isAntiLambda();
                // For now, removing the ambiguous candidates from the analysis. Derived data permits handling both.
                // (From Podolanski-Armenteros plots, the population of ambiguous is ~2% without TOF, and without 
                //  competing mass rejection. From those, ~99% seem to be K0s, so no real gain in considering the
                //  ambiguous candidates in the analysis)
                if (isLambda && isAntiLambda) continue;
                const float v0pt = v0.v0Pt();
                const float v0eta = v0.v0Eta();
                const float v0phi = v0.v0Phi();

                float v0LambdaLikeMass = 0; // Initialized just to catch any stray behavior
                float protonLikePt = 0;
                float protonLikeEta = 0;
                float protonLikePhi = 0;
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
                float lambdaRapidity = lambdaLike4Vec.Rapidity(); // For further kinematic selections

                // Boosting proton into lambda frame:
                XYZVector beta = -lambdaLike4Vec.BoostToCM(); // Boost trivector that goes from laboratory frame to the rest frame
                auto protonLike4VecStar = ROOT::Math::VectorUtil::boost(protonLike4Vec, beta);

                // Getting unit vectors and 3-components:
                XYZVector lambdaLike3Vec = lambdaLike4Vec.Vect();
                XYZVector protonLikeStarUnit3Vec = protonLike4VecStar.Vect().Unit();

                ////////////////////////////////////////////
                // Ring observable: Leading particle proxy
                // Always computed -- leading particle existence is guaranteed by the second check above
                ////////////////////////////////////////////
                // Cross product
                XYZVector crossLeadP = leadPUnitVec.Cross(lambdaLike3Vec);
                float ringObservableLeadP = protonLikeStarUnit3Vec.Dot(crossLeadP) / crossLeadP.R();
                // Adding the prefactor related to the CP-violating decay (decay constants have different signs)
                if (!forcePolSignQA) ringObservableLeadP *= (isLambda) ? polPrefactorLambda : polPrefactorAntiLambda;
                else ringObservableLeadP *= (isLambda) ? polPrefactorLambda : -1.0 * polPrefactorAntiLambda;
                // Angular variables
                float deltaPhiLeadP = wrapToPiFast(v0phi - leadPPhi); // Wrapped to [-PI, pi), for convenience
                float deltaThetaLeadP = ROOT::Math::VectorUtil::Angle(leadPUnitVec, lambdaLike3Vec); // 3D angular separation

                //////////////////////////////////////////
                // Ring observable: Leading jet proxy
                // Only computed when a leading jet exists in this collision.
                //////////////////////////////////////////
                float ringObservable = 0.;
                float deltaPhiJet = 0.;
                float deltaThetaJet = 0.;
                if (hasValidLeadingJet) {
                    // Cross product
                    XYZVector cross = leadingJetUnitVec.Cross(lambdaLike3Vec);
                    ringObservable = protonLikeStarUnit3Vec.Dot(cross) / cross.R();
                    // Adding prefactor
                    if (!forcePolSignQA) ringObservable *= (isLambda) ? polPrefactorLambda : polPrefactorAntiLambda;
                    else ringObservable *= (isLambda) ? polPrefactorLambda : -1.0 * polPrefactorAntiLambda;
                    // Angular variables
                    deltaPhiJet = wrapToPiFast(v0phi - leadingJetPhi); 
                    deltaThetaJet = ROOT::Math::VectorUtil::Angle(leadingJetUnitVec, lambdaLike3Vec); 
                }

                //////////////////////////////////////////
                // Ring observable: Subleading jet proxy
                // Only computed when a subleading jet exists in this collision.
                //////////////////////////////////////////
                float ringObservable2ndJet = 0.;
                float deltaPhi2ndJet = 0.;
                float deltaTheta2ndJet = 0.;
                if (hasValidSubJet) {
                    XYZVector cross2ndJet = subJetUnitVec.Cross(lambdaLike3Vec);
                    ringObservable2ndJet = protonLikeStarUnit3Vec.Dot(cross2ndJet) / cross2ndJet.R();
                    // Adding prefactor
                    if (!forcePolSignQA) ringObservable2ndJet *= (isLambda) ? polPrefactorLambda : polPrefactorAntiLambda;
                    else ringObservable2ndJet *= (isLambda) ? polPrefactorLambda : -1.0 * polPrefactorAntiLambda;
                    // Angular variables
                    deltaPhi2ndJet = wrapToPiFast(v0phi - subleadingJetPhi);
                    deltaTheta2ndJet = ROOT::Math::VectorUtil::Angle(subJetUnitVec, lambdaLike3Vec);
                }

                // Calculating polarization observables (in the Lambda frame, because that is easier -- does not require boosts):
                    // To be precise, not actually the polarization, but a part of the summand in P^*_\Lambda = (3/\alpha_\Lambda) * <p^*_{proton}>
                float PolStarX, PolStarY, PolStarZ;
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

                float v0phiToFillHists = wrapToPiFast(v0phi); // A short wrap to reuse some predefined axes

                // Fill ring histograms: (1D, lambda 2D correlations and jet 2D correlations):
                RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "Ring") // Notice the usage of macros! If you change the variable names, this WILL break the code!
                                                                          // No, there should NOT be any ";" here! Read the macro definition for an explanation
                histos.fill(HIST("pRingCutsLeadingP"), 0, ringObservableLeadP); // First bin of comparison
                POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "Ring")

                if (hasValidLeadingJet) {
                    RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "Ring")
                    histos.fill(HIST("pRingCuts"), 0, ringObservable);
                }
                if (hasValidSubJet) {
                    RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "Ring")
                    histos.fill(HIST("pRingCutsSubLeadingJet"), 0, ringObservable2ndJet);
                }

                // Extra kinematic criteria for Lambda candidates (removes polarization background):
                const bool kinematicLambdaCheck = (v0pt > 0.5 && v0pt < 1.5) && std::abs(lambdaRapidity) < 0.5;
                if (kinematicLambdaCheck){
                    RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
                    histos.fill(HIST("pRingCutsLeadingP"), 1, ringObservableLeadP);
                    POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
                    if (hasValidLeadingJet) {
                        RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
                        histos.fill(HIST("pRingCuts"), 1, ringObservable);
                    }
                    if (hasValidSubJet) {
                        RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
                        histos.fill(HIST("pRingCutsSubLeadingJet"), 1, ringObservable2ndJet);
                    }
                }
                
                // Extra selection criteria on jet candidates:
                // (redundant for jets with R=0.4, but for jets with R<0.4 the leading jet may be farther in eta)
                if (kinematicJetCheck){ // Already includes hasValidLeadingJet in the bool! (no need to check again)
                    RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
                    histos.fill(HIST("pRingCuts"), 2, ringObservable);
                    POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
                }

                // Extra selection criteria on both Lambda and jet candidates:
                if (kinematicLambdaCheck && kinematicJetCheck){
                    RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
                    histos.fill(HIST("pRingCuts"), 3, ringObservable);
                    POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
                }

                // Same variations for the leading particle and for the subleading jet:
                if (kinematicLeadPCheck){
                    RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
                    histos.fill(HIST("pRingCutsLeadingP"), 2, ringObservableLeadP);
                }
                if (kinematic2ndJetCheck){
                    RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
                    histos.fill(HIST("pRingCutsSubLeadingJet"), 2, ringObservable2ndJet);
                }
                if (kinematicLambdaCheck && kinematicLeadPCheck){
                    RING_OBSERVABLE_LEADP_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
                    histos.fill(HIST("pRingCutsLeadingP"), 3, ringObservableLeadP);
                }
                if (kinematicLambdaCheck && kinematic2ndJetCheck){
                    RING_OBSERVABLE_2NDJET_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
                    histos.fill(HIST("pRingCutsSubLeadingJet"), 3, ringObservable2ndJet);
                }
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