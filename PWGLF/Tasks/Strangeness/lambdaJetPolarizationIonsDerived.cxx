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
///
/// \author Cicero Domenico Muncinelli <cicero.domenico.muncinelli@cern.ch>, Campinas State University
//
// Jet Polarization Ions task -- Derived data
// ================
//
// This code loops over custom derived data tables defined on 
// lambdaJetPolarizationIons.h (JetsRing, LambdaLikeV0sRing).
// From this derived data, calculates polarization on an EbE
// basis.
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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace ROOT::Math;
// using namespace o2::aod::lambdajetpol; // Used it explicitly along the code for clarity
struct lambdaJetPolarizationIonsDerived {
    const double protonMass = o2::constants::physics::MassProton; // Assumes particle identification for daughter is perfect
    const double lambdaWeakDecayConstant = 0.747; // Fixed as a constant for this analysis.






    void processPolarizationData(aod::RingCollisions const& collisions, aod::RingJets const& jets, aod::RingLambdaLikeV0s const& v0s){

        for (auto const& collision : collisions) {
            const auto collId = collision.collIdx();
            const double centrality = collision.centrality();

            // Slice jets and V0s belonging to this collision
            auto jetsInColl = jets.sliceBy(o2::aod::lambdajetpol::collIdx, collId);
            auto v0sInColl  = v0s.sliceBy(o2::aod::lambdajetpol::collIdx, collId);

            // Check if there is at least one V0 and one jet in the collision:
            // (in the way I fill the table, there is always at least one V0 in
            //  the stored collision, but the jets table can not be filled for
            //  that collision, and a collision may not be filled when the jets
            //  table is. Be mindful of that!)
            if (jetsInColl.empty() || v0sInColl.empty()) continue;

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

            // TODO: add centrality selection procedure and options (one configurable for no centrality separation at all too!)
            // TODO: add Lambda candidate selection. Think of a statistical method like signal extraction (if possible) for ring polarization
            // TODO: add calculations with second to leading jet too.
            // TODO: Add calculations with leading particle
            for (auto const& v0 : v0sInColl) {
                const double v0pt = v0.V0Pt();
                const double v0eta = v0.v0Eta();
                const double v0phi = v0.Phi();

                double v0LambdaMass;
                double protonLikePt;
                double protonLikeEta;
                double protonLikePhi;
                if (isLambda){
                    v0LambdaMass = v0.MassLambda();
                    protonLikePt = v0.PosPt();
                    protonLikeEta = v0.PosEta();
                    protonLikePhi = v0.PosPhi();
                }
                // (TODO: implement AntiLambda polarization)
                // if (isAntiLambda){
                // }
                
                const double lambdaPx = v0pt * std::cos(v0phi);
                const double lambdaPy = v0pt * std::sin(v0phi);
                const double lambdaPz = v0pt * std::sinh(v0eta);
                const double lambdaMomentumSquared = lambdaPx*lambdaPx + lambdaPy*lambdaPy + lambdaPz*lambdaPz;
                const double lambdaE = std::sqrt(v0LambdaMass*v0LambdaMass + lambdaMomentumSquared);
                // const double lambdaRapidity = 0.5 * std::log((lambdaE + lambdaPz) / (lambdaE - lambdaPz)); // For extra selection criteria later
                const double lambdaRapidity = std::atanh(lambdaPz / lambdaE); // More numerically stable

                const double protonPx = protonLikePt * std::cos(protonLikePhi);
                const double protonPy = protonLikePt * std::sin(protonLikePhi);
                const double protonPz = protonLikePt * std::sinh(protonLikeEta);
                const double protonMomentumSquared = protonPx*protonPx + protonPy*protonPy + protonPz*protonPz;
                const double protonE = std::sqrt(protonMass*protonMass + protonMomentumSquared);

                TLorentzVector lambdaLike4Vec(lambdaPx, lambdaPy, lambdaPz, lambdaE);
                TLorentzVector protonLike4Vec(protonPx, protonPy, protonPz, protonE);

                // Boosting proton into lambda frame:
                TVector3 betaInverse = -lambdaLike4Vec.BoostVector(); // Boost trivector that goes from laboratory frame to the rest frame
                TLorentzVector protonLike4VecStar = protonLike4Vec;
                protonLike4VecStar.Boost(betaInverse);

                TVector3 protonLikeStarUnitVec = (protonLike4VecStar.Vect()).Unit();
                TVector3 jetUnitVec = (TVector3(jetPx, jetPy, jetPz)).Unit();
                TVector3 lambda3Vec = TVector3(lambdaPx, lambdaPy, lambdaPz);

                // Calculating inner product:
                double crossX = jetUnitVec.Y()*lambdaPz - jetUnitVec.Z()*lambdaPy;
                double crossY = jetUnitVec.Z()*lambdaPx - jetUnitVec.X()*lambdaPz;
                double crossZ = jetUnitVec.X()*lambdaPy - jetUnitVec.Y()*lambdaPx;
                double crossProductNorm = std::sqrt(crossX*crossX + crossY*crossY + crossZ*crossZ);

                double ringObservable = 3./(lambdaWeakDecayConstant) * ((protonLikeStarUnitVec.X()*crossX
                                        + protonLikeStarUnitVec.Y()*crossY + protonLikeStarUnitVec.Z()*cross_z)/crossProductNorm);

                // Calculating error bars:
                double ringObservableSquared = ringObservable*ringObservable;

                // Angular variables:
                double deltaPhiJet = v0phi - leadingJetPhi;
                double deltaThetaJet = ; // 3D angular separation

                // Fill ring histograms: (1D, lambda 2D correlations and jet 2D correlations): (TODO)
                histos.fill(HIST("Ring/hRingObservableDeltaPhi"), 0);
                histos.fill(HIST("Ring/hRingObservableIntegrated"), 0);

                histos.fill(HIST("Ring/hRingObservableDeltaPhi"), 0);
                histos.fill(HIST("Ring/hRingObservableIntegrated"), 0);
                
                
                // Extra kinematic criteria for Lambda candidates (removes polarization background):
                const bool kinematicLambdaCheck = (v0pt > 0.5 && v0pt < 1.5) && std::abs(lambdaRapidity) < 0.5;
                if (kinematicLambdaCheck){
                    histos.fill(HIST("RingKinematicCuts/hRingObservableDeltaPhi"), 0);
                    histos.fill(HIST("RingKinematicCuts/hRingObservableIntegrated"), 0);
                }
                
                
                // Extra selection criteria on jet candidates:
                const bool kinematicJetCheck = std::abs(leadingJetEta) < 0.5;
                if (kinematicJetCheck){ // This is redundant for jets with R=0.4, but for jets with R<0.4 the leading jet may be farther in eta.

                }
                

                // Extra selection criteria on Lambda and jet candidates:
                if (kinematicLambdaCheck && kinematicJetCheck){

                }
                

            } // end v0s loop
        } // end collisions
    }
};