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
    X(FOLDER "/hRingObservableLambdaPt",                     v0pt)                                    \
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
    X(FOLDER "/p2dRingObservableDeltaThetaVsLeadJetPt",      deltaThetaJet, leadingJetPt, ringObservable)
  // Lambda mass correlations (1D + 1D and 2D + 1D): (TODO: signal extraction attempts)

// ======================================================
// Ring Observable SQUARED histogram fill list
// ======================================================
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
    X(FOLDER "/h2dRingObservableSquaredDeltaThetaVsLeadJetPt",deltaThetaJet, leadingJetPt, ringObservableSquared)

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

    /////////////////////////
    // Configurable blocks:
    // Histogram axes configuration:
    struct : ConfigurableGroup {
        std::string prefix = "axisConfigurations"; // JSON group name
        ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
        ConfigurableAxis axisPtCoarse{"axisPtCoarse", {VARIABLE_WIDTH, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 7.0f, 10.0f, 15.0f}, "pt axis for QA"};
        ConfigurableAxis axisLambdaMass{"axisLambdaMass", {450, 1.08f, 1.15f}, "Lambda mass in GeV/c"}; // Default is {200, 1.101f, 1.131f}

        // Jet axes:
        ConfigurableAxis axisLeadingParticlePt{"axisLeadingParticlePt",{100, 0.f, 200.f},"Leading particle p_{T} (GeV/c)"}; // Simpler version!
        ConfigurableAxis axisJetPt{"axisJetPt",{100, 0.f, 200.f},"Jet p_{t} (GeV)"};
        ConfigurableAxis axisCosTheta{"axisDeltaTheta", {100, 0, constants::math::PI}, "#Delta #theta_{jet}"};
        ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {100, -constants::math::PI, constants::math::PI}, "#Delta #phi_{jet}"};
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
            histos.add((folder + "/hRingObservableDeltaPhi").c_str(), "hRingObservableDeltaPhi", kTH1D, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/hRingObservableDeltaTheta").c_str(), "hRingObservableDeltaTheta", kTH1D, {axisConfigurations.axisCosTheta});
            histos.add((folder + "/hRingObservableIntegrated").c_str(), "hRingObservableIntegrated", kTH1D, {{1, -0.5, 0.5}});
            // Counters (denominators)
            histos.add((folder + "/hDeltaPhi").c_str(), "hDeltaPhi", kTH1D, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/hDeltaTheta").c_str(), "hDeltaTheta", kTH1D, {axisConfigurations.axisCosTheta});
            histos.add((folder + "/hIntegrated").c_str(), "hIntegrated", kTH1D, {{1, -0.5, 0.5}});
            // ===============================
            // Lambda pT dependence
            // ===============================
            histos.add((folder + "/hRingObservableLambdaPt").c_str(), "hRingObservableLambdaPt", kTH1D, {axisConfigurations.axisPt});
            // ===============================
            // 2D Lambda correlations
            // ===============================
            histos.add((folder + "/h2dRingObservableDeltaPhiVsLambdaPt").c_str(), "h2dRingObservableDeltaPhiVsLambdaPt", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPt});
            histos.add((folder + "/h2dRingObservableDeltaThetaVsLambdaPt").c_str(), "h2dRingObservableDeltaThetaVsLambdaPt", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisPt});
            // Counters
            histos.add((folder + "/h2dDeltaPhiVsLambdaPt").c_str(), "h2dDeltaPhiVsLambdaPt", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPt});
            histos.add((folder + "/h2dDeltaThetaVsLambdaPt").c_str(), "h2dDeltaThetaVsLambdaPt", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisPt});
            // ===============================
            // 2D Jet correlations
            // ===============================
            histos.add((folder + "/h2dRingObservableDeltaPhiVsLeadJetPt").c_str(), "h2dRingObservableDeltaPhiVsLeadJetPt", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisJetPt});
            histos.add((folder + "/h2dRingObservableDeltaThetaVsLeadJetPt").c_str(), "h2dRingObservableDeltaThetaVsLeadJetPt", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisJetPt});
            // Counters
            histos.add((folder + "/h2dDeltaPhiVsLeadJetPt").c_str(), "h2dDeltaPhiVsLeadJetPt", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisJetPt});
            histos.add((folder + "/h2dDeltaThetaVsLeadJetPt").c_str(), "h2dDeltaThetaVsLeadJetPt", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisJetPt});
            // ===============================
            // Squared observable (error propagation)
            // ===============================
            histos.add((folder + "/hRingObservableSquaredDeltaPhi").c_str(), "hRingObservableSquaredDeltaPhi", kTH1D, {axisConfigurations.axisDeltaPhi});
            histos.add((folder + "/hRingObservableSquaredDeltaTheta").c_str(), "hRingObservableSquaredDeltaTheta", kTH1D, {axisConfigurations.axisCosTheta});
            histos.add((folder + "/hRingObservableSquaredIntegrated").c_str(), "hRingObservableSquaredIntegrated", kTH1D, {{1, -0.5, 0.5}});
            histos.add((folder + "/hRingObservableSquaredLambdaPt").c_str(), "hRingObservableSquaredLambdaPt", kTH1D, {axisConfigurations.axisPt});
            histos.add((folder + "/h2dRingObservableSquaredDeltaPhiVsLambdaPt").c_str(), "h2dRingObservableSquaredDeltaPhiVsLambdaPt", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPt});
            histos.add((folder + "/h2dRingObservableSquaredDeltaThetaVsLambdaPt").c_str(), "h2dRingObservableSquaredDeltaThetaVsLambdaPt", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisPt});
            histos.add((folder + "/h2dRingObservableSquaredDeltaPhiVsLeadJetPt").c_str(), "h2dRingObservableSquaredDeltaPhiVsLeadJetPt", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisJetPt});
            histos.add((folder + "/h2dRingObservableSquaredDeltaThetaVsLeadJetPt").c_str(), "h2dRingObservableSquaredDeltaThetaVsLeadJetPt", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisJetPt});

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
            histos.add((folder + "/pRingObservableDeltaTheta").c_str(), "pRingObservableDeltaTheta;cos#theta_{jet};<#it{R}>", kTProfile, {axisConfigurations.axisCosTheta});
            histos.add((folder + "/pRingObservableIntegrated").c_str(), "pRingObservableIntegrated; ;<#it{R}>", kTProfile, {{1, -0.5, 0.5}});
            histos.add((folder + "/pRingObservableLambdaPt").c_str(), "pRingObservableLambdaPt;#it{p}_{T}^{#Lambda};<#it{R}>", kTProfile, {axisConfigurations.axisPt});
            // ===============================
            // 2D TProfiles (Lambda correlations)
            // ===============================
            histos.add((folder + "/p2dRingObservableDeltaPhiVsLambdaPt").c_str(), "p2dRingObservableDeltaPhiVsLambdaPt;#Delta#varphi_{jet};#it{p}_{T}^{#Lambda};<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisPt});
            histos.add((folder + "/p2dRingObservableDeltaThetaVsLambdaPt").c_str(), "p2dRingObservableDeltaThetaVsLambdaPt;cos#theta_{jet};#it{p}_{T}^{#Lambda};<#it{R}>", kTProfile2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisPt});
            // ===============================
            // 2D TProfiles (Jet correlations)
            // ===============================
            histos.add((folder + "/p2dRingObservableDeltaPhiVsLeadJetPt").c_str(), "p2dRingObservableDeltaPhiVsLeadJetPt;#Delta#varphi_{jet};#it{p}_{T}^{lead jet};<#it{R}>", kTProfile2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisJetPt});
            histos.add((folder + "/p2dRingObservableDeltaThetaVsLeadJetPt").c_str(), "p2dRingObservableDeltaThetaVsLeadJetPt;cos#theta_{jet};#it{p}_{T}^{lead jet};<#it{R}>", kTProfile2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisJetPt});
            // (TODO: add mass histograms for signal extraction)
        };
        // Execute local lambda to register histogram families:
        addRingObservableFamily("Ring");
        addRingObservableFamily("RingKinematicCuts");
        addRingObservableFamily("JetKinematicCuts");
        addRingObservableFamily("JetAndLambdaKinematicCuts");
    }

    ////////////// Fill Ring Observable histograms:
    ///(This block was tranformed into a bunch of #define statements at the top of the code)
    //////////////

        // Preslices for correct collisions association:
    Preslice<aod::RingJets> perColJets = o2::aod::lambdajetpol::collIdx; // Slicing by the key that comes with the index column
    Preslice<aod::RingLaV0s> perColV0s = o2::aod::lambdajetpol::collIdx;
    void processPolarizationData(o2::aod::RingCollisions const& collisions, o2::aod::RingJets const& jets, o2::aod::RingLaV0s const& v0s){
        for (auto const& collision : collisions) {
            const auto collId = collision.collIdx();
            // const double centrality = collision.centrality(); // (TODO: implement centrality!)

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

            // TODO: add centrality selection procedure and options (one configurable for no centrality separation at all too!)
            // TODO: add Lambda candidate selection. Think of a statistical method like signal extraction (if possible) for ring polarization
            // TODO: add calculations with second to leading jet too.
            // TODO: Add calculations with leading particle
            for (auto const& v0 : v0sInColl) {
                const bool isLambda = v0.isLambda();
                const bool isAntiLambda = v0.isAntiLambda();
                if (isLambda && isAntiLambda) continue; // For now, removing the ambiguous candidates from the analysis. Derived data permits handling both.
                const double v0pt = v0.v0Pt();
                const double v0eta = v0.v0Eta();
                const double v0phi = v0.v0Phi();

                double v0LambdaLikeMass;
                double protonLikePt;
                double protonLikeEta;
                double protonLikePhi;
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
                ringObservable *= (isLambda) ? 3./lambdaWeakDecayConstant : 3./antiLambdaWeakDecayConstant;

                // Calculating error bars:
                double ringObservableSquared = ringObservable*ringObservable;

                // Angular variables:
                double deltaPhiJet = wrapToPiFast(v0phi - leadingJetPhi); // Wrapped to [-PI, pi), for convenience
                double deltaThetaJet = ROOT::Math::VectorUtil::Angle(leadingJetUnitVec, lambdaLike3Vec); // 3D angular separation

                // Fill ring histograms: (1D, lambda 2D correlations and jet 2D correlations):
                RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "Ring") // Notice the usage of macros! If you change the variable names, this WILL break the code!
                RING_OBSERVABLE_SQUARED_FILL_LIST(APPLY_HISTO_FILL, "Ring") // No, there should NOT be any ";" here! Read the macro definition for an explanation

                // Extra kinematic criteria for Lambda candidates (removes polarization background):
                const bool kinematicLambdaCheck = (v0pt > 0.5 && v0pt < 1.5) && std::abs(lambdaRapidity) < 0.5;
                if (kinematicLambdaCheck){
                    RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
                    RING_OBSERVABLE_SQUARED_FILL_LIST(APPLY_HISTO_FILL, "RingKinematicCuts")
                }
                
                // Extra selection criteria on jet candidates:
                const bool kinematicJetCheck = std::abs(leadingJetEta) < 0.5;
                if (kinematicJetCheck){ // This is redundant for jets with R=0.4, but for jets with R<0.4 the leading jet may be farther in eta.
                    RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
                    RING_OBSERVABLE_SQUARED_FILL_LIST(APPLY_HISTO_FILL, "JetKinematicCuts")
                }
                
                // Extra selection criteria on both Lambda and jet candidates:
                if (kinematicLambdaCheck && kinematicJetCheck){
                    RING_OBSERVABLE_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
                    RING_OBSERVABLE_SQUARED_FILL_LIST(APPLY_HISTO_FILL, "JetAndLambdaKinematicCuts")
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