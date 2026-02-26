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
#define RING_OBSERVABLE_FILL_LIST(X, FOLDER)                                                         \
    /* 1D observable histograms */                                                                     \
    X(FOLDER "/hRingObservableDeltaPhi",                     deltaPhiJet,   ringObservable)           \
    X(FOLDER "/hRingObservableDeltaTheta",                   deltaThetaJet, ringObservable)           \
    X(FOLDER "/hRingObservableIntegrated",                   0.,            ringObservable)           \
    /* Counters */                                                                                     \
    X(FOLDER "/hDeltaPhi",                                   deltaPhiJet)                             \
    X(FOLDER "/hDeltaTheta",                                 deltaThetaJet)                           \
    X(FOLDER "/hIntegrated",                                 0.)                                      \
    /* Lambda pT variation -- Youpeng's proposal */                                                   \
    X(FOLDER "/hRingObservableLambdaPt",                     v0pt,           ringObservable)          \
    X(FOLDER "/hLambdaPt",                                   v0pt)                                    \
    /* 2D Lambda correlations */                                                                       \
    X(FOLDER "/h2dRingObservableDeltaPhiVsLambdaPt",         deltaPhiJet,   v0pt, ringObservable)     \
    X(FOLDER "/h2dRingObservableDeltaThetaVsLambdaPt",       deltaThetaJet, v0pt, ringObservable)     \
    /* Counters */                                                                                     \
    X(FOLDER "/h2dDeltaPhiVsLambdaPt",                       deltaPhiJet,   v0pt)                     \
    X(FOLDER "/h2dDeltaThetaVsLambdaPt",                     deltaThetaJet, v0pt)                     \
    /* 2D Jet correlations */                                                                          \
    X(FOLDER "/h2dRingObservableDeltaPhiVsLeadJetPt",        deltaPhiJet,   leadingJetPt, ringObservable) \
    X(FOLDER "/h2dRingObservableDeltaThetaVsLeadJetPt",      deltaThetaJet, leadingJetPt, ringObservable) \
    /* Counters */                                                                                     \
    X(FOLDER "/h2dDeltaPhiVsLeadJetPt",                      deltaPhiJet,   leadingJetPt)             \
    X(FOLDER "/h2dDeltaThetaVsLeadJetPt",                    deltaThetaJet, leadingJetPt)             \
    /* Additional plots for instant gratification - 1D Profiles */                                    \
    X(FOLDER "/pRingObservableDeltaPhi",                     deltaPhiJet,   ringObservable)           \
    X(FOLDER "/pRingObservableDeltaTheta",                   deltaThetaJet, ringObservable)           \
    X(FOLDER "/pRingObservableIntegrated",                   0.,            ringObservable)           \
    X(FOLDER "/pRingObservableLambdaPt",                     v0pt,          ringObservable)           \
    /* 2D Profiles */                                                                                  \
    X(FOLDER "/p2dRingObservableDeltaPhiVsLambdaPt",         deltaPhiJet,   v0pt, ringObservable)     \
    X(FOLDER "/p2dRingObservableDeltaThetaVsLambdaPt",       deltaThetaJet, v0pt, ringObservable)     \
    X(FOLDER "/p2dRingObservableDeltaPhiVsLeadJetPt",        deltaPhiJet,   leadingJetPt, ringObservable) \
    X(FOLDER "/p2dRingObservableDeltaThetaVsLeadJetPt",      deltaThetaJet, leadingJetPt, ringObservable) \
    /* 1D Mass */ \
    X(FOLDER "/hMass", v0LambdaLikeMass) \
    X(FOLDER "/hRingObservableMass", v0LambdaLikeMass, ringObservable) \
    X(FOLDER "/hMassSigExtract", v0LambdaLikeMass) \
    /* 2D: Observable vs Mass */ \
    X(FOLDER "/h2dRingObservableDeltaPhiVsMass",      deltaPhiJet,   v0LambdaLikeMass, ringObservable) \
    X(FOLDER "/h2dRingObservableDeltaThetaVsMass",    deltaThetaJet, v0LambdaLikeMass, ringObservable) \
    /* Counters */ \
    X(FOLDER "/h2dDeltaPhiVsMass",                    deltaPhiJet,   v0LambdaLikeMass) \
    X(FOLDER "/h2dDeltaThetaVsMass",                  deltaThetaJet, v0LambdaLikeMass) \
    /* 3D: Observable vs Mass vs Lambda pT */ \
    X(FOLDER "/h3dRingObservableDeltaPhiVsMassVsLambdaPt",   deltaPhiJet,   v0LambdaLikeMass, v0pt, ringObservable) \
    X(FOLDER "/h3dRingObservableDeltaThetaVsMassVsLambdaPt", deltaThetaJet, v0LambdaLikeMass, v0pt, ringObservable) \
    /* Counters */ \
    X(FOLDER "/h3dDeltaPhiVsMassVsLambdaPt",          deltaPhiJet,   v0LambdaLikeMass, v0pt) \
    X(FOLDER "/h3dDeltaThetaVsMassVsLambdaPt",        deltaThetaJet, v0LambdaLikeMass, v0pt) \
    /* 3D: Observable vs Mass vs Lead Jet pT */ \
    X(FOLDER "/h3dRingObservableDeltaPhiVsMassVsLeadJetPt",   deltaPhiJet,   v0LambdaLikeMass, leadingJetPt, ringObservable) \
    X(FOLDER "/h3dRingObservableDeltaThetaVsMassVsLeadJetPt", deltaThetaJet, v0LambdaLikeMass, leadingJetPt, ringObservable) \
    /* Counters */ \
    X(FOLDER "/h3dDeltaPhiVsMassVsLeadJetPt",                 deltaPhiJet,   v0LambdaLikeMass, leadingJetPt) \
    X(FOLDER "/h3dDeltaThetaVsMassVsLeadJetPt",               deltaThetaJet, v0LambdaLikeMass, leadingJetPt) \
    /* 2D: Observable vs Mass vs Centrality (projected as 2D Mass vs Cent for integrated observable) */ \
    X(FOLDER "/h2dRingObservableMassVsCent", v0LambdaLikeMass, centrality, ringObservable) \
    /* 3D: Observable vs Mass vs Centrality */ \
    X(FOLDER "/h3dRingObservableDeltaPhiVsMassVsCent",   deltaPhiJet,   v0LambdaLikeMass, centrality, ringObservable) \
    X(FOLDER "/h3dRingObservableDeltaThetaVsMassVsCent", deltaThetaJet, v0LambdaLikeMass, centrality, ringObservable) \
    /* Counters */ \
    X(FOLDER "/h3dDeltaPhiVsMassVsCent",   deltaPhiJet,   v0LambdaLikeMass, centrality) \
    X(FOLDER "/h3dDeltaThetaVsMassVsCent", deltaThetaJet, v0LambdaLikeMass, centrality) \
    /* TProfile of Ring vs Mass */ \
    X(FOLDER "/pRingObservableMass", v0LambdaLikeMass, ringObservable) \
    /* 2D Profiles: Angle vs Mass */ \
    X(FOLDER "/p2dRingObservableDeltaPhiVsMass",   deltaPhiJet,   v0LambdaLikeMass, ringObservable) \
    X(FOLDER "/p2dRingObservableDeltaThetaVsMass", deltaThetaJet, v0LambdaLikeMass, ringObservable) \
    /* 3D Profiles: Angle vs Mass vs Lambda pT */ \
    X(FOLDER "/p3dRingObservableDeltaPhiVsMassVsLambdaPt",   deltaPhiJet,   v0LambdaLikeMass, v0pt, ringObservable) \
    X(FOLDER "/p3dRingObservableDeltaThetaVsMassVsLambdaPt", deltaThetaJet, v0LambdaLikeMass, v0pt, ringObservable) \
    /* 3D Profiles: Angle vs Mass vs Lead Jet pT */ \
    X(FOLDER "/p3dRingObservableDeltaPhiVsMassVsLeadJetPt",   deltaPhiJet,   v0LambdaLikeMass, leadingJetPt, ringObservable) \
    X(FOLDER "/p3dRingObservableDeltaThetaVsMassVsLeadJetPt", deltaThetaJet, v0LambdaLikeMass, leadingJetPt, ringObservable) \
    /* 2D Profile: Mass vs Centrality */ \
    X(FOLDER "/p2dRingObservableMassVsCent", v0LambdaLikeMass, centrality, ringObservable) \
    /* 3D Profiles: Angle vs Mass vs Centrality */ \
    X(FOLDER "/p3dRingObservableDeltaPhiVsMassVsCent",   deltaPhiJet,   v0LambdaLikeMass, centrality, ringObservable) \
    X(FOLDER "/p3dRingObservableDeltaThetaVsMassVsCent", deltaThetaJet, v0LambdaLikeMass, centrality, ringObservable)
    // (TODO: add counters for regular TH2Ds about centrality)


// ======================================================
// Ring Observable SQUARED histogram fill list
// ======================================================
#if 0 // Disabling the whole definition in a cleaner way -- Multiline comments do not work appropriately in this type of macro definition!
#define RING_OBSERVABLE_SQUARED_FILL_LIST(X, FOLDER)                                                \
    /* 1D observable histograms */                                                                    \
    X(FOLDER "/hRingObservableSquaredDeltaPhi",               deltaPhiJet,   ringObservableSquared)  \
    X(FOLDER "/hRingObservableSquaredDeltaTheta",             deltaThetaJet, ringObservableSquared)  \
    X(FOLDER "/hRingObservableSquaredIntegrated",             0.,            ringObservableSquared)  \
    /* Lambda pT variation */                                                                         \
    X(FOLDER "/hRingObservableSquaredLambdaPt",               v0pt,          ringObservableSquared)  \
    /* 2D Lambda correlations */                                                                      \
    X(FOLDER "/h2dRingObservableSquaredDeltaPhiVsLambdaPt",   deltaPhiJet,   v0pt,          ringObservableSquared) \
    X(FOLDER "/h2dRingObservableSquaredDeltaThetaVsLambdaPt", deltaThetaJet, v0pt,          ringObservableSquared) \
    /* 2D Jet correlations */                                                                         \
    X(FOLDER "/h2dRingObservableSquaredDeltaPhiVsLeadJetPt",  deltaPhiJet,   leadingJetPt, ringObservableSquared) \
    X(FOLDER "/h2dRingObservableSquaredDeltaThetaVsLeadJetPt",deltaThetaJet, leadingJetPt, ringObservableSquared) \
    /* 2D - Mass correlations */ \
    X(FOLDER "/h2dRingObservableSquaredDeltaPhiVsMass",      deltaPhiJet,   v0LambdaLikeMass, ringObservableSquared) \
    X(FOLDER "/h2dRingObservableSquaredDeltaThetaVsMass",    deltaThetaJet, v0LambdaLikeMass, ringObservableSquared) \
    /* 3D - LambdaPt */ \
    X(FOLDER "/h3dRingObservableSquaredDeltaPhiVsMassVsLambdaPt",   deltaPhiJet,   v0LambdaLikeMass, v0pt, ringObservableSquared) \
    X(FOLDER "/h3dRingObservableSquaredDeltaThetaVsMassVsLambdaPt", deltaThetaJet, v0LambdaLikeMass, v0pt, ringObservableSquared) \
    /* 3D - LeadJetPt*/ \
    X(FOLDER "/h3dRingObservableSquaredDeltaPhiVsMassVsLeadJetPt",   deltaPhiJet,   v0LambdaLikeMass, leadingJetPt, ringObservableSquared) \
    X(FOLDER "/h3dRingObservableSquaredDeltaThetaVsMassVsLeadJetPt", deltaThetaJet, v0LambdaLikeMass, leadingJetPt, ringObservableSquared) \
    /* 3D: Squared observable vs Mass vs Centrality */ \
    X(FOLDER "/h3dRingObservableSquaredDeltaPhiVsMassVsCent",   deltaPhiJet,   v0LambdaLikeMass, centrality, ringObservableSquared) \
    X(FOLDER "/h3dRingObservableSquaredDeltaThetaVsMassVsCent", deltaThetaJet, v0LambdaLikeMass, centrality, ringObservableSquared)
#endif

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


// // Another macro for the significance histograms expansion: // (Moved into signal extraction post-processing to allow for pipelining)
// #define RING_1DSIGNIFICANCE_LIST(X, FOLDER) \
//     X(FOLDER "/pRingObservableDeltaPhi",        FOLDER "/hRingSignificanceDeltaPhi") \
//     X(FOLDER "/pRingObservableDeltaTheta",      FOLDER "/hRingSignificanceDeltaTheta") \
//     X(FOLDER "/pRingObservableIntegrated",      FOLDER "/hRingSignificanceIntegrated") \
//     X(FOLDER "/pRingObservableLambdaPt",        FOLDER "/hRingSignificanceLambdaPt") \
//     X(FOLDER "/pRingObservableMass",            FOLDER "/hRingSignificanceMass")

// #define APPLY_RING_SIGNIFICANCE(PROFILE, HISTO) \
//     fillSignificance(histos.get<TProfile>(HIST(PROFILE)), histos.get<TH1>(HIST(HISTO)));


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
            // Rewrote the axisLambdaMassSigextract to have 5x coarser bins outside the peak region, and 2x coarser bins in the peak region
            // (this allows for better fits and smaller fluctuations)
        // ConfigurableAxis axisLambdaMassSigExtract{"axisLambdaMassSigExtract", {VARIABLE_WIDTH, 1.080000, 1.082000, 1.084000, 1.086000, 1.088000, 1.090000, 1.092000, 1.094000, 1.096000, 1.098000, 1.100000, 1.102000, 1.104000, 1.106000, 1.108000, 1.109683, 1.110483, 1.111283, 1.112083, 1.112883, 1.113683, 1.114483, 1.115283, 1.116083, 1.116883, 1.117683, 1.118483, 1.119283, 1.120083, 1.120883, 1.121683, 1.123683, 1.125683, 1.127683, 1.129683, 1.131683, 1.133683, 1.135683, 1.137683, 1.139683, 1.141683, 1.143683, 1.145683, 1.147683, 1.149683, 1.150000}, "Lambda mass in GeV/c"};
            // Even coarser axis:
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
            // (TODO: check if all squared observables were actually transformed into a (much better) TProfile version)
            // // ===============================
            // // Squared observable (error propagation)
            // // ===============================
            // histos.add((folder + "/hRingObservableSquaredDeltaPhi").c_str(), "hRingObservableSquaredDeltaPhi", kTH1D, {axisConfigurations.axisDeltaPhi});
            // histos.add((folder + "/hRingObservableSquaredDeltaTheta").c_str(), "hRingObservableSquaredDeltaTheta", kTH1D, {axisConfigurations.axisDeltaTheta});
            // histos.add((folder + "/hRingObservableSquaredIntegrated").c_str(), "hRingObservableSquaredIntegrated", kTH1D, {{1, -0.5, 0.5}});
            // histos.add((folder + "/hRingObservableSquaredLambdaPt").c_str(), "hRingObservableSquaredLambdaPt", kTH1D, {axisConfigurations.axisPt});
            // histos.add((folder + "/h2dRingObservableSquaredDeltaPhiVsLambdaPt").c_str(), "h2dRingObservableSquaredDeltaPhiVsLambdaPt", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPt});
            // histos.add((folder + "/h2dRingObservableSquaredDeltaThetaVsLambdaPt").c_str(), "h2dRingObservableSquaredDeltaThetaVsLambdaPt", kTH2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisPt});
            // histos.add((folder + "/h2dRingObservableSquaredDeltaPhiVsLeadJetPt").c_str(), "h2dRingObservableSquaredDeltaPhiVsLeadJetPt", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisJetPt});
            // histos.add((folder + "/h2dRingObservableSquaredDeltaThetaVsLeadJetPt").c_str(), "h2dRingObservableSquaredDeltaThetaVsLeadJetPt", kTH2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisJetPt});

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
            // // Squared observable (for variance propagation)
            // histos.add((folder + "/h2dRingObservableSquaredDeltaPhiVsMass").c_str(), "h2dRingObservableSquaredDeltaPhiVsMass", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract});
            // histos.add((folder + "/h2dRingObservableSquaredDeltaThetaVsMass").c_str(), "h2dRingObservableSquaredDeltaThetaVsMass", kTH2D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract});
            // --- Counters (denominators) ---
            histos.add((folder + "/h2dDeltaPhiVsMass").c_str(), "h2dDeltaPhiVsMass", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract});
            histos.add((folder + "/h2dDeltaThetaVsMass").c_str(), "h2dDeltaThetaVsMass", kTH2D,{axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract});
            // --- 3D: Angular observable vs Mass vs Lambda pT ---
            histos.add((folder + "/h3dRingObservableDeltaPhiVsMassVsLambdaPt").c_str(), "h3dRingObservableDeltaPhiVsMassVsLambdaPt", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
            histos.add((folder + "/h3dRingObservableDeltaThetaVsMassVsLambdaPt").c_str(), "h3dRingObservableDeltaThetaVsMassVsLambdaPt", kTH3D, {axisConfigurations.axisDeltaTheta,axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
            // // Squared version
            // histos.add((folder + "/h3dRingObservableSquaredDeltaPhiVsMassVsLambdaPt").c_str(), "h3dRingObservableSquaredDeltaPhiVsMassVsLambdaPt", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
            // histos.add((folder + "/h3dRingObservableSquaredDeltaThetaVsMassVsLambdaPt").c_str(), "h3dRingObservableSquaredDeltaThetaVsMassVsLambdaPt", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
            // Counters
            histos.add((folder + "/h3dDeltaPhiVsMassVsLambdaPt").c_str(), "h3dDeltaPhiVsMassVsLambdaPt", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
            histos.add((folder + "/h3dDeltaThetaVsMassVsLambdaPt").c_str(), "h3dDeltaThetaVsMassVsLambdaPt", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisPtSigExtract});
            // --- 3D: Angular observable vs Mass vs Lead Jet pT ---
            histos.add((folder + "/h3dRingObservableDeltaPhiVsMassVsLeadJetPt").c_str(), "h3dRingObservableDeltaPhiVsMassVsLeadJetPt", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract}); 
            histos.add((folder + "/h3dRingObservableDeltaThetaVsMassVsLeadJetPt").c_str(), "h3dRingObservableDeltaThetaVsMassVsLeadJetPt", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract}); 
            // // --- Squared version ---
            // histos.add((folder + "/h3dRingObservableSquaredDeltaPhiVsMassVsLeadJetPt").c_str(), "h3dRingObservableSquaredDeltaPhiVsMassVsLeadJetPt", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract}); 
            // histos.add((folder + "/h3dRingObservableSquaredDeltaThetaVsMassVsLeadJetPt").c_str(), "h3dRingObservableSquaredDeltaThetaVsMassVsLeadJetPt", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract}); 
            // --- Counters ---
            histos.add((folder + "/h3dDeltaPhiVsMassVsLeadJetPt").c_str(), "h3dDeltaPhiVsMassVsLeadJetPt", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract}); 
            histos.add((folder + "/h3dDeltaThetaVsMassVsLeadJetPt").c_str(), "h3dDeltaThetaVsMassVsLeadJetPt", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisJetPtSigExtract});

            /////////////////////////////////////
            /// TProfiles with the proper errors for quick glancing
            /////////////////////////////////////
            // TProfile of ring vs mass (integrated in all phi, and properly normalized by N_\Lambda):
            histos.add((folder + "/pRingObservableMass").c_str(), "pRingObservableMass;m_{p#pi};<#it{R}>", kTProfile, {axisConfigurations.axisLambdaMassSigExtract});
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
            // // Squared observable (for variance propagation)
            // histos.add((folder + "/h3dRingObservableSquaredDeltaPhiVsMassVsCent").c_str(), "h3dRingObservableSquaredDeltaPhiVsMassVsCent", kTH3D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
            // histos.add((folder + "/h3dRingObservableSquaredDeltaThetaVsMassVsCent").c_str(), "h3dRingObservableSquaredDeltaThetaVsMassVsCent", kTH3D, {axisConfigurations.axisDeltaTheta, axisConfigurations.axisLambdaMassSigExtract, axisConfigurations.axisCentrality});
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
            // Added to a separate folder for further control (changed the usage of the "folder" string):
            histos.add(("QA_Numbers/" + folder + "/hEventsWithV0").c_str(), "hEventsWithV0", kTH1D, {{1,0,1}}); // In the current derived data, all saved events have a V0
            histos.add(("QA_Numbers/" + folder + "/hLambdaCounter").c_str(), "hLambdaCounter", kTH1D, {{1,0,1}});
            histos.add(("QA_Numbers/" + folder + "/hAntiLambdaCounter").c_str(), "hAntiLambdaCounter", kTH1D, {{1,0,1}});
            histos.add(("QA_Numbers/" + folder + "/hValidLeadJets").c_str(), "hValidLeadJets", kTH1D, {{1,0,1}});
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
    Preslice<aod::RingJets> perColJets = o2::aod::lambdajetpol::collisionId;
    Preslice<aod::RingLaV0s> perColV0s = o2::aod::lambdajetpol::collisionId;
    void processPolarizationData(o2::aod::RingCollisions const& collisions, o2::aod::RingJets const& jets, o2::aod::RingLaV0s const& v0s){
        // Custom grouping
        // (TODO: test using global index's custom grouping utility as in the sigma0builder)
        // std::vector<std::vector<int>> v0Grouped(collisions.size());
        // for (const auto& v0 : v0s) {v0Grouped[v0.collisionId()].push_back(v0.globalIndex());}
        // std::vector<std::vector<int>> jetsGrouped(collisions.size());
        // for (const auto& jet : jets) {jetsGrouped[jet.collisionId()].push_back(jet.globalIndex());}
            // Only really need the leading jet for now:
            // -1 means "no jet for this collision"
        // std::vector<int> leadingJetIndex(collisions.size(), -1);
        // std::vector<float> leadingJetPt(collisions.size(), -1.f);
        // for (const auto& jet : jets) {
        //     int collId = jet.collisionId();
        //     float pt = jet.pt();  // or whatever pT accessor you use
        //     if (pt > leadingJetPt[collId]) {
        //         leadingJetPt[collId] = pt;
        //         leadingJetIndex[collId] = jet.globalIndex();
        //     }
        // }
        
        for (auto const& collision : collisions) {
            const auto collId = collision.collisionId();
            const double centrality = collision.centrality();

            // Slice jets and V0s belonging to this collision
                // (global collision indices repeat a lot, but they are unique to a same TimeFrame (TF) subfolder in the derived data)
            auto jetsInColl = jets.sliceBy(perColJets, collId);
            auto v0sInColl  = v0s.sliceBy(perColV0s, collId);

            // Alternative custom grouping block:
            // int jetIndex = leadingJetIndex[collision.globalIndex()];
            // if (jetIndex < 0) continue;
            // auto leadingJet = jets.rawIteratorAt(jetIndex);

            // Check if there is at least one V0 and one jet in the collision:
            // (in the way I fill the table, there is always at least one V0 in
            //  the stored collision, but the jets table can not be filled for
            //  that collision, and a collision may not be filled when the jets
            //  table is. Be mindful of that!)
            if (!jetsInColl.size() || !v0sInColl.size()) continue;

            // Get leading jet:
            double leadingJetPt = -1;
            o2::aod::RingJets::iterator leadingJet;
            for (auto const& jet : jetsInColl) {
                const auto jetpt = jet.jetPt();
                if (jetpt > leadingJetPt){
                    leadingJetPt = jetpt;
                    leadingJet = jet;
                }
            }

            // Now you can use:
            const double leadingJetEta = leadingJet.jetEta();
            const double leadingJetPhi = leadingJet.jetPhi();

            // Convert to 3-vector components for inner product:
            const double jetPx = leadingJetPt * std::cos(leadingJetPhi);
            const double jetPy = leadingJetPt * std::sin(leadingJetPhi);
            const double jetPz = leadingJetPt * std::sinh(leadingJetEta);
            XYZVector leadingJetVec(jetPx, jetPy, jetPz);
            XYZVector leadingJetUnitVec = leadingJetVec.Unit();

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

            // TODO: add centrality selection procedure and options (one configurable for no centrality separation at all too!)
            // TODO: add Lambda candidate selection. Think of a statistical method like signal extraction (if possible) for ring polarization
            // TODO: add calculations with second to leading jet too.
            // TODO: Add calculations with leading particle

            // Custom grouping alternative:
            // for (size_t i = 0; i < v0Grouped[coll.globalIndex()].size(); i++) {
            //     auto v0 = v0s.rawIteratorAt(v0Grouped[collision.globalIndex()][i]);
            // }

            for (auto const& v0 : v0sInColl) {
                const bool isLambda = v0.isLambda();
                const bool isAntiLambda = v0.isAntiLambda();
                if (isLambda && isAntiLambda) continue; // For now, removing the ambiguous candidates from the analysis. Derived data permits handling both.
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
                else if (isAntiLambda){ // (TODO: add a split histogram where you consider Lambda and AntiLambda polarization separately)
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

                // // Calculating error bars:
                // double ringObservableSquared = ringObservable*ringObservable;

                // Angular variables:
                double deltaPhiJet = wrapToPiFast(v0phi - leadingJetPhi); // Wrapped to [-PI, pi), for convenience
                double deltaThetaJet = ROOT::Math::VectorUtil::Angle(leadingJetUnitVec, lambdaLike3Vec); // 3D angular separation

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
                // RING_OBSERVABLE_SQUARED_FILL_LIST(APPLY_HISTO_FILL, "Ring") // No, there should NOT be any ";" here! Read the macro definition for an explanation
                POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "Ring")

                // Extra kinematic criteria for Lambda candidates (removes polarization background):
                const bool kinematicLambdaCheck = (v0pt > 0.5 && v0pt < 1.5) && std::abs(lambdaRapidity) < 0.5;
                if (kinematicLambdaCheck){
                    RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
                    // RING_OBSERVABLE_SQUARED_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
                    POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
                }
                
                // Extra selection criteria on jet candidates:
                const bool kinematicJetCheck = std::abs(leadingJetEta) < 0.5;
                if (kinematicJetCheck){ // This is redundant for jets with R=0.4, but for jets with R<0.4 the leading jet may be farther in eta.
                    RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
                    // RING_OBSERVABLE_SQUARED_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
                    POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
                }
                
                // Extra selection criteria on both Lambda and jet candidates:
                if (kinematicLambdaCheck && kinematicJetCheck){
                    RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
                    // RING_OBSERVABLE_SQUARED_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
                    POLARIZATION_PROFILE_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
                }
            } // end v0s loop
        } // end collisions
    }

    // // Filling final histograms for QA based on previous 1D TProfiles, after all processing has been done
    // void finalize(){
    //     RING_1DSIGNIFICANCE_LIST(APPLY_RING_SIGNIFICANCE, "Ring")
    //     RING_1DSIGNIFICANCE_LIST(APPLY_RING_SIGNIFICANCE, "RingKinematicCuts")
    //     RING_1DSIGNIFICANCE_LIST(APPLY_RING_SIGNIFICANCE, "JetKinematicCuts")
    //     RING_1DSIGNIFICANCE_LIST(APPLY_RING_SIGNIFICANCE, "JetAndLambdaKinematicCuts")
    // }
    
    PROCESS_SWITCH(lambdajetpolarizationionsderived, processPolarizationData, "Process derived data in Run 3 Data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdajetpolarizationionsderived>(cfgc)};
}

// Avoid macro leakage!
#undef APPLY_HISTO_FILL