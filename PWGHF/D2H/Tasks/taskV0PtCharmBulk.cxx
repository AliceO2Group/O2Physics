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

/// \file taskV0PtCharmBulk.cxx
/// \brief v0 pt for the charm-bulk correlation analysis task
/// \author Wu Chuntai, UNIPD, CCNU, and INFN Padova
/// \author Andrea Rossi, INFN Padova

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>
#include <TF1.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <string>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;

enum DecayChannel {
    D0ToPiK = 0,
    D0ToKPi
};

struct HfTaskV0PtCharmBulk
{
    // General configuration
    Configurable<float> etaAMin{"etaAMin", -0.8, "eta min for A subevent"};
    Configurable<float> etaAMax{"etaAMax", -0.2, "eta max for A subevent"};
    Configurable<float> etaBMin{"etaBMin", 0.2, "eta min for B subevent"};
    Configurable<float> etaBMax{"etaBMax", 0.8, "eta max for B subevent"};
    Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4)"};
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indexes of BDT scores to be stored. Two indexes max."};

    // Track configuration
    Configurable<int> tpcNClsCrossedRowsMin{"tpcNClsCrossedRowsMin", 70, "min. TPC crossed rows for associated tracks"};
    Configurable<float> etaTrkMax{"etaTrkMax", 1., "max. track eta"};
    Configurable<float> ptTrkMin{"ptTrkMin", 0.2, "min. track pT"};
    Configurable<float> ptTrkMax{"ptTrkMax", 5., "max. track pT"};
    Configurable<float> dcaXYTrkMax{"dcaXYTrkMax", 1., "max. track DCA XY"};
    Configurable<float> dcaZTrkMax{"dcaZTrkMax", 1., "max. track DCA Z"};
    Configurable<bool> usePtDiffDcaXYCut{"usePtDiffDcaXYCut", true, "Use pt-differential DCAxy cut for associated tracks"};
    Configurable<float> dcaXYTrkNSigmaMax{"dcaXYTrkNSigmaMax", 7, "Cut on number of sigma deviations from expected DCA in the transverse direction"};
    Configurable<std::string> dcaXYPtPrimTrkFunc{"dcaXYPtPrimTrkFunc", "(0.0026+0.005/(x^1.01))", "Functional form of pt-dependent DCAxy cut"};

    // Candidate configuration
    Configurable<int> selectionFlag{"selectionFlag", 1, "Selection Flag for hadron (ML score tables are required to run the task)"};

    // Event configuration
    Configurable<bool> forceCharmInCollision{"forceCharmInCollision", true, "Flag to force charm in collision"};
    HfEventSelection hfEvSel; // event selection and monitoring

    o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb{};
    SliceCache cache;
    TF1* funcDcaXYPtCutPrimTrk = nullptr;

    // Subcribe and join the tables
    using TracksTable = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra>>;
    using D0CandTable = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>>;
    using CollsWithCentMult = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>;

    // Select tracks and candidates
    Filter filterSelectTracks = (nabs(aod::track::eta) < etaTrkMax) && (aod::track::pt > ptTrkMin) && (aod::track::pt < ptTrkMax) && (nabs(aod::track::dcaXY) < dcaXYTrkMax) && (nabs(aod::track::dcaZ) < dcaZTrkMax);
    Filter filterSelectD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;

    // pre-slice by collision for tracks
    Preslice<TracksTable> tracksTablePerColl = aod::track::collisionId;
    Preslice<D0CandTable> candD0TablePerColl = aod::hf_cand::collisionId;

    // Partitions for selected candidates
    Partition<D0CandTable> selectedD0ToPiK = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag;
    Partition<D0CandTable> selectedD0ToKPi = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;

    // Partition<TracksTable> selectedTracks = (nabs(aod::track::eta) < etaTrkMax) && (aod::track::pt > ptTrkMin) && (aod::track::pt < ptTrkMax) && (nabs(aod::track::dcaXY) < dcaXYTrkMax) && (nabs(aod::track::dcaZ) < dcaZTrkMax);

    // THnSparse configuration
    ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent",  {100, 0, 100}, ""};
    ConfigurableAxis thnConfigAxisCandMass{"thnConfigAxisCandMass", {200, 1.68, 2.08}, ""};
    ConfigurableAxis thnConfigAxisCandPt{"thnConfigAxisCandPt", {10, 0., 10.}, ""};
    ConfigurableAxis thnConfigAxisCandEta{"thnConfigAxisCandEta", {200, -1., 1.}, ""};
    ConfigurableAxis thnConfigAxisMPtTrkA{"thnConfigAxisMPtTrkA", {100, 0., 5.}, ""};
    ConfigurableAxis thnConfigAxisMPtTrkB{"thnConfigAxisMPtTrkB", {100, 0., 5.}, ""};
    ConfigurableAxis thnConfigAxisMlOne{"thnConfigAxisMlOne", {100, 0., 1.}, ""};
    ConfigurableAxis thnConfigAxisMlTwo{"thnConfigAxisMlTwo", {100, 0., 1.}, ""};

    HistogramRegistry registry{"registry", {}};

    void init(InitContext&) {
        std::array<bool, 14> doprocess{doprocessD0WCentMult};
        if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) == 0) {
            LOGP(fatal, "At least one process function should be enabled at a time.");
        }

        hfEvSel.addHistograms(registry); // collision monitoring
        ccdb->setURL(ccdbUrl);
        ccdb->setCaching(true);
        ccdb->setLocalObjectValidityChecking();

        // Define the axes for the THnSparse
        const AxisSpec axisCent = {thnConfigAxisCent, "Centrality"};
        const AxisSpec axisCandMass = {thnConfigAxisCandMass, "Inv. mass (GeV/#it{c}^{2})"};
        const AxisSpec axisCandPt = {thnConfigAxisCandPt, "#it{p}_{T} (GeV/#it{c})"};
        const AxisSpec axisCandEta = {thnConfigAxisCandEta, "Eta"};
        const AxisSpec axisMPtTrkA = {thnConfigAxisMPtTrkA, "Mean pT of tracks in subevent A (GeV/#it{c})"};
        const AxisSpec axisMPtTrkB = {thnConfigAxisMPtTrkB, "Mean pT of tracks in subevent B (GeV/#it{c})"};
        const AxisSpec axisMlOne = {thnConfigAxisMlOne, "ML score 1"};
        const AxisSpec axisMlTwo = {thnConfigAxisMlTwo, "ML score 2"};

        std::vector<AxisSpec> axes = {axisCent, axisCandMass, axisCandPt, axisCandEta, axisMPtTrkA, axisMPtTrkB, axisMlOne, axisMlTwo};
        registry.add("hSparseMeanTrkPtCharm", "THn for mean track pT", HistType::kTHnSparseF, axes);

        registry.add("hMeanTrkPtAVsCent", "Mean track pT A vs centrality", HistType::kTH2F, {axisCent, axisMPtTrkA});
        registry.add("hMeanTrkPtBVsCent", "Mean track pT B vs centrality", HistType::kTH2F, {axisCent, axisMPtTrkB});

        if (usePtDiffDcaXYCut) {
            funcDcaXYPtCutPrimTrk = new TF1("funcDcaXYPtCutPrimTrk", Form("[0]*%s", dcaXYPtPrimTrkFunc.value.data()), 0.001, 100);
            funcDcaXYPtCutPrimTrk->SetParameter(0, dcaXYTrkNSigmaMax);
            LOGF(info, "DCAxy pt-dependence function: %s", Form("[0]*%s", dcaXYPtPrimTrkFunc.value.data()));
        }
    } // End init

    /// Function to calculate mean pT of tracks in a given eta range (subevent) for a specific collision
    /// \param tracks are the tracks to be used for mean pT calculation
    /// \param candidates are the D0 candidates to be analyzed
    /// \return a pair of pairs: {{meanPtA, countA}, {meanPtB, countB}}
    template <typename CandT>
    std::pair<float, float> calculateMeanPt(TracksTable const& tracks, CandT const& candidates)
    {
        float sumPtA = 0.f;
        int countA = 0;
        float sumPtB = 0.f;
        int countB = 0;
        std::vector<int> candProngsA;
        std::vector<int> candProngsB;

        // Gather the global indices of the candidate daughters
        if constexpr (std::is_same_v<CandT, D0CandTable>)
        {
            for (const auto& cand : candidates) {
                if (cand.eta() > etaAMin && cand.eta() < etaAMax) {
                    candProngsA.push_back(cand.prong0Id());
                    candProngsA.push_back(cand.prong1Id());
                } else if (cand.eta() > etaBMin && cand.eta() < etaBMax) {
                    candProngsB.push_back(cand.prong0Id());
                    candProngsB.push_back(cand.prong1Id());
                }
            }
        }

        // Loop over tracks and calculate sum of pT and count for subevent A and B, excluding candidate daughters
        for (const auto& track : tracks)
        {
            float eta = track.eta();
            float pt = track.pt();

            // Select only global tracks with DCA information or with sufficient TPC crossed rows
            if (track.isGlobalTrackWoDCA() || track.tpcNClsCrossedRows() < tpcNClsCrossedRowsMin) {
                continue;
            }

            // Apply DCA cuts
            if (usePtDiffDcaXYCut) 
            {
                float const dcaXYTrkCut = funcDcaXYPtCutPrimTrk->Eval(pt);
                if (std::fabs(track.dcaXY()) > dcaXYTrkCut) 
                {
                    continue;
                }
            }

            int const trackGlobalIndex = track.globalIndex();
            // A side
            if (eta < etaAMax && eta > etaAMin) 
            {
                if (std::find(candProngsB.begin(), candProngsB.end(), trackGlobalIndex) != candProngsB.end()) 
                {
                    continue; // skip tracks that are daughters of the candidate in the opposite subevent
                    /*TODO: Considering remove the overlapped daughter tracks by recalculating the mean pT 
                    excluding the candidate daughters on the fly for each candidate*/
                }
                sumPtA += pt;
                countA++;
            }

            // B side
            if (eta < etaBMax && eta > etaBMin) 
            {
                if (std::find(candProngsA.begin(), candProngsA.end(), trackGlobalIndex) != candProngsA.end()) 
                {
                    continue; // skip tracks that are daughters of the candidate in the opposite subevent
                    /*TODO: Considering remove the overlapped daughter tracks by recalculating the mean pT 
                    excluding the candidate daughters on the fly for each candidate*/
                }
                sumPtB += pt;
                countB++;
            }
        }

        if (countA == 0 && countB == 0) 
        {
            return {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
            // return NaN if no tracks in either subevent
        } else if (countA == 0) 
        {
            return {std::numeric_limits<float>::quiet_NaN(), sumPtB / countB}; // NaN for subevent A if no tracks
        } else if (countB == 0) 
        {
            return {sumPtA / countA, std::numeric_limits<float>::quiet_NaN()}; // NaN for subevent B if no tracks
        }
        return {sumPtA / countA, sumPtB / countB};
    }

    /// Calculate mean pT of tracks in subevent A and B, and fill the THnSparse
    /// \param collision is the collision with the centrality and multiplicity information
    /// \param tracks are the tracks to be used for mean pT calculation
    /// \param candidates are the D0 candidates to be analyzed
    template <DecayChannel Channel, typename CandT>
    void runCharmBulkAnalysis(CandT const& candidates, TracksTable const& tracks, float cent)
    {
        auto [meanPtA, meanPtB] = calculateMeanPt(tracks, candidates);
        // Loop over candidates and fill the THnSparse
        for (const auto& cand : candidates) 
        {
            float invMass = 0.f;
            std::vector<float> outputMl = {-999., -999.};
            if constexpr (std::is_same_v<CandT, D0CandTable>) {
                switch (Channel) {
                    case DecayChannel::D0ToPiK:
                        invMass = HfHelper::invMassD0ToPiK(cand);
                        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
                            outputMl[iclass] = cand.mlProbD0()[classMl->at(iclass)];
                        }
                        break;
                    case DecayChannel::D0ToKPi:
                        invMass = HfHelper::invMassD0barToKPi(cand);
                        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
                            outputMl[iclass] = cand.mlProbD0bar()[classMl->at(iclass)];
                        }
                        break;
                }
            }
            registry.fill(HIST("hSparseMeanTrkPtCharm"), cent, invMass, cand.pt(), cand.eta(), meanPtA, meanPtB, outputMl[0], outputMl[1]);
            registry.fill(HIST("hMeanTrkPtAVsCent"), cent, meanPtA);
            registry.fill(HIST("hMeanTrkPtBVsCent"), cent, meanPtB);
        }
    }

    /// Check event selections for collision and fill event selection histograms
    /// \param collision is the collision
    template <typename Coll>
    bool isSelectedHfCollision(Coll const& collision, float& cent)
    {
        o2::hf_evsel::HfCollisionRejectionMask collRejMask{};
        if (centEstimator == CentralityEstimator::FT0A)
        {
            collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0A, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
        } else if (centEstimator == CentralityEstimator::FT0C)
        {
            collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
        } else if (centEstimator == CentralityEstimator::FT0M)
        {
            collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0M, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
        } else if (centEstimator == CentralityEstimator::FV0A)
        {
            collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FV0A, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
        } else 
        {
            LOG(fatal) << "Centrality estimator not recognized for collision selection";
            std::abort();
        }
            hfEvSel.fillHistograms(collision, collRejMask, cent);
        return collRejMask == 0;
    }

    void processD0WCentMult(CollsWithCentMult::iterator const& collision, TracksTable const& tracks, D0CandTable const& /* D0 candidates */) {
        auto tableD0ToPiK = selectedD0ToPiK->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
        auto tableD0ToKPi = selectedD0ToKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
        if (forceCharmInCollision && tableD0ToPiK.size() < 1 && tableD0ToKPi.size() < 1) {
            return;
        }
        float cent = -1.f;
        if (!isSelectedHfCollision(collision, cent)) {
            return;
        }    
        runCharmBulkAnalysis<DecayChannel::D0ToPiK>(tableD0ToPiK, tracks, cent);
        runCharmBulkAnalysis<DecayChannel::D0ToKPi>(tableD0ToKPi, tracks, cent);
    }
    PROCESS_SWITCH(HfTaskV0PtCharmBulk, processD0WCentMult, "Process D0 candidates for tracks' mean-pT analysis", true);
}; // End struct HfTaskV0PtCharmBulk

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
    return WorkflowSpec{adaptAnalysisTask<HfTaskV0PtCharmBulk>(cfgc)};
}
