#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <TAxis.h>
#include <TEfficiency.h>
#include <THashList.h>
#include <TMathBase.h>
#include <TString.h>

#include <array>
#include <cmath>
#include <memory>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct getMothers {

    HistogramRegistry registry{"registry"};
    
    void init(InitContext&) {
        registry.add("hMother", "Mother particles", kTH1F, {{2, -0.5, 1.5}});
        registry.add("hReco", "Reconstructed D+", kTH1F, {{2, -0.5, 1.5}});
    }
    
    Preslice<o2::aod::McParticles> mcParticlesPerColl = o2::aod::mcparticle::mcCollisionId;
    Preslice<o2::aod::Tracks> tracksPerColl = aod::track::collisionId;


    /// Finds the mother of an MC particle by looking for the expected PDG code in the mother chain.
    /// \tparam acceptFlavourOscillation  switch to accept decays where the mother oscillated (e.g. B0 -> B0bar)
    /// \param particlesMC  table with MC particles
    /// \param particle  MC particle
    /// \param pdgMother  expected mother PDG code
    /// \param acceptAntiParticles  switch to accept the antiparticle of the expected mother
    /// \param sign  antiparticle indicator of the found mother w.r.t. pdgMother; 1 if particle, -1 if antiparticle, 0 if mother not found
    /// \param depthMax  maximum decay tree level to check; Mothers up to this level will be considered. If -1, all levels are considered.
    /// \return index of the mother particle if found, -1 otherwise
    template <bool acceptFlavourOscillation = false, typename T>
    static int getMotherFixed(const T& particlesMC,
                        const typename T::iterator& particle,
                        int pdgMother,
                        bool acceptAntiParticles = false,
                        int8_t* sign = nullptr,
                        int8_t depthMax = -1)
    {
        int8_t sgn = 0;           // 1 if the expected mother is particle, -1 if antiparticle (w.r.t. pdgMother)
        int indexMother = -1;     // index of the final matched mother, if found
        int stage = 0;            // mother tree level (just for debugging)
        bool motherFound = false; // true when the desired mother particle is found in the kine tree
        if (sign) {
        *sign = sgn;
        }

        // vector of vectors with mother indices; each line corresponds to a "stage"
        std::vector<std::vector<int64_t>> arrayIds{};
        std::vector<int64_t> initVec{particle.globalIndex()};
        arrayIds.push_back(initVec); // the first vector contains the index of the original particle

        while (!motherFound && arrayIds[-stage].size() > 0 && (depthMax < 0 || -stage < depthMax)) {
        // vector of mother indices for the current stage
        std::vector<int64_t> arrayIdsStage{};
        for (auto iPart : arrayIds[-stage]) { // check all the particles that were the mothers at the previous stage, o2-linter: disable=const-ref-in-for-loop (int elements)
            auto particleMother = particlesMC.rawIteratorAt(iPart - particlesMC.offset());
            if (particleMother.has_mothers()) {
            for (const auto& mother : particleMother.template mothers_as<aod::McParticles>()) { // loop over the mother particles of the analysed particle
                auto iMother = mother.globalIndex();
                if (std::find(arrayIdsStage.begin(), arrayIdsStage.end(), iMother) != arrayIdsStage.end()) {                       // if a mother is still present in the vector, do not check it again
                continue;
                }
                // Check mother's PDG code.
                auto pdgParticleIMother = mother.pdgCode(); // PDG code of the mother
                // printf("getMother: ");
                // for (int i = stage; i < 0; i++) // Indent to make the tree look nice.
                //   printf(" ");
                // printf("Stage %d: Mother PDG: %d, Index: %d\n", stage, pdgParticleIMother, iMother);
                if (pdgParticleIMother == pdgMother) { // exact PDG match
                sgn = 1;
                indexMother = iMother;
                motherFound = true;
                break;
                } else if (acceptAntiParticles && pdgParticleIMother == -pdgMother) { // antiparticle PDG match
                sgn = -1;
                indexMother = iMother;
                motherFound = true;
                break;
                }
                // add mother index in the vector for the current stage
                arrayIdsStage.push_back(iMother);
            }
            }
        }
        // add vector of mother indices for the current stage
        arrayIds.push_back(arrayIdsStage);
        stage--;
        }
        if (sign) {
        if constexpr (acceptFlavourOscillation) {
            if (std::abs(particle.getGenStatusCode()) == 42) { // take possible flavour oscillation of B0(s) mother into account
            sgn *= -1;                                                                      // select the sign of the mother after oscillation (and not before)
            }
        }
        *sign = sgn;
        }

        return indexMother;
    }


    /// Checks whether the reconstructed decay candidate is the expected decay.
    /// \tparam acceptFlavourOscillation  switch to accept decays where the mother oscillated (e.g. B0 -> B0bar)
    /// \tparam checkProcess  switch to accept only decay daughters by checking the production process of MC particles
    /// \tparam acceptIncompleteReco  switch to accept candidates with only part of the daughters reconstructed
    /// \tparam acceptTrackDecay  switch to accept candidates with daughter tracks of pions and kaons which decayed
    /// \tparam acceptTrackIntWithMaterial switch to accept candidates with final (i.e. p, K, pi) daughter tracks interacting with material
    /// \param particlesMC  table with MC particles
    /// \param arrDaughters  array of candidate daughters
    /// \param pdgMother  expected mother PDG code
    /// \param arrPdgDaughters  array of expected daughter PDG codes
    /// \param acceptAntiParticles  switch to accept the antiparticle version of the expected decay
    /// \param sign  antiparticle indicator of the found mother w.r.t. pdgMother; 1 if particle, -1 if antiparticle, 0 if mother not found
    /// \param depthMax  maximum decay tree level to check; Daughters up to this level will be considered. If -1, all levels are considered.
    /// \param nPiToMu  number of pion prongs decayed to a muon
    /// \param nKaToPi  number of kaon prongs decayed to a pion
    /// \param nInteractionsWithMaterial  number of daughter particles that interacted with material
    /// \return index of the mother particle if the mother and daughters are correct, -1 otherwise
    template <bool acceptFlavourOscillation = false, bool checkProcess = false, bool acceptIncompleteReco = false, bool acceptTrackDecay = false, bool acceptTrackIntWithMaterial = false, std::size_t N, typename T, typename U>
    static int getMatchedMCRecFixed(const T& particlesMC,
                                    const std::array<U, N>& arrDaughters,
                                    int pdgMother,
                                    std::array<int, N> arrPdgDaughters,
                                    bool acceptAntiParticles = false,
                                    int8_t* sign = nullptr,
                                    int depthMax = 1,
                                    int8_t* nPiToMu = nullptr,
                                    int8_t* nKaToPi = nullptr,
                                    int8_t* nInteractionsWithMaterial = nullptr)
    {
        // Printf("MC Rec: Expected mother PDG: %d", pdgMother);
        int8_t coefFlavourOscillation = 1;         // 1 if no B0(s) flavour oscillation occured, -1 else
        int8_t sgn = 0;                            // 1 if the expected mother is particle, -1 if antiparticle (w.r.t. pdgMother)
        int8_t nPiToMuLocal = 0;                   // number of pion prongs decayed to a muon
        int8_t nKaToPiLocal = 0;                   // number of kaon prongs decayed to a pion
        int8_t nInteractionsWithMaterialLocal = 0; // number of interactions with material
        int indexMother = -1;                      // index of the mother particle
        std::vector<int> arrAllDaughtersIndex;     // vector of indices of all daughters of the mother of the first provided daughter
        std::array<int, N> arrDaughtersIndex;      // array of indices of provided daughters
        if (sign) {
        *sign = sgn;
        }
        if constexpr (acceptFlavourOscillation) {
        // Loop over decay candidate prongs to spot possible oscillation decay product
        for (std::size_t iProng = 0; iProng < N; ++iProng) {
            if (!arrDaughters[iProng].has_mcParticle()) {
            return -1;
            }
            auto particleI = arrDaughters[iProng].template mcParticle_as<T>();                 // ith daughter particle
            if (std::abs(particleI.getGenStatusCode()) == 42) { // oscillation decay product spotted
            coefFlavourOscillation = -1;                                                     // select the sign of the mother after oscillation (and not before)
            break;
            }
        }
        }
        // Loop over decay candidate prongs
        for (std::size_t iProng = 0; iProng < N; ++iProng) {
        if (!arrDaughters[iProng].has_mcParticle()) {
            return -1;
        }
        auto particleI = arrDaughters[iProng].template mcParticle_as<T>(); // ith daughter particle
        if constexpr (acceptTrackDecay) {
            // Replace the MC particle associated with the prong by its mother for π → μ and K → π.
            auto motherI = particleI.template mothers_first_as<T>();
            auto pdgI = std::abs(particleI.pdgCode());
            auto pdgMotherI = std::abs(motherI.pdgCode());
            if (pdgI == PDG_t::kMuonMinus && pdgMotherI == PDG_t::kPiPlus) {
            // π → μ
            nPiToMuLocal++;
            particleI = motherI;
            } else if (pdgI == PDG_t::kPiPlus && pdgMotherI == PDG_t::kKPlus) {
            // K → π
            nKaToPiLocal++;
            particleI = motherI;
            }
        }
        if constexpr (acceptTrackIntWithMaterial) {
            // Replace the MC particle associated with the prong by its mother for part → part due to material interactions.
            // It keeps looking at the mother iteratively, until it finds a particle from decay or primary
            auto process = particleI.getProcess();
            auto pdgI = std::abs(particleI.pdgCode());
            auto pdgMotherI = std::abs(particleI.pdgCode());
            while (process != TMCProcess::kPDecay && process != TMCProcess::kPPrimary && pdgI == pdgMotherI) {
            if (!particleI.has_mothers()) {
                break;
            }
            auto motherI = particleI.template mothers_first_as<T>();
            pdgI = std::abs(particleI.pdgCode());
            pdgMotherI = std::abs(motherI.pdgCode());
            if (pdgI == pdgMotherI) {
                particleI = motherI;
                process = particleI.getProcess();
                if (process == TMCProcess::kPDecay || process == TMCProcess::kPPrimary) { // we found the original daughter that interacted with material
                nInteractionsWithMaterialLocal++;
                }
            }
            }
        }
        arrDaughtersIndex[iProng] = particleI.globalIndex();
        // Get the list of daughter indices from the mother of the first prong.
        if (iProng == 0) {
            // Get the mother index and its sign.
            // PDG code of the first daughter's mother determines whether the expected mother is a particle or antiparticle.
            indexMother = getMotherFixed(particlesMC, particleI, pdgMother, acceptAntiParticles, &sgn, depthMax);
            // Check whether mother was found.
            if (indexMother <= -1) {
            // Printf("MC Rec: Rejected: bad mother index or PDG");
            return -1;
            }
            // Printf("MC Rec: Good mother: %d", indexMother);
            auto particleMother = particlesMC.rawIteratorAt(indexMother - particlesMC.offset());
            // Check the daughter indices.
            if (!particleMother.has_daughters()) {
            // Printf("MC Rec: Rejected: bad daughter index range: %d-%d", particleMother.daughtersIds().front(), particleMother.daughtersIds().back());
            return -1;
            }
            // Check that the number of direct daughters is not larger than the number of expected final daughters.
            if constexpr (!acceptIncompleteReco && !checkProcess) {
            if (particleMother.daughtersIds().back() - particleMother.daughtersIds().front() + 1 > static_cast<int>(N)) {
                // Printf("MC Rec: Rejected: too many direct daughters: %d (expected %ld final)", particleMother.daughtersIds().back() - particleMother.daughtersIds().front() + 1, N);
                return -1;
            }
            }
            // Get the list of actual final daughters.
            RecoDecay::getDaughters<checkProcess>(particleMother, &arrAllDaughtersIndex, arrPdgDaughters, depthMax);
            // printf("MC Rec: Mother %d has %d final daughters:", indexMother, arrAllDaughtersIndex.size());
            // for (auto i : arrAllDaughtersIndex) {
            //   printf(" %d", i);
            // }
            // printf("\n");
            //  Check whether the number of actual final daughters is equal to the number of expected final daughters (i.e. the number of provided prongs).
            if (!acceptIncompleteReco && arrAllDaughtersIndex.size() != N) {
            // Printf("MC Rec: Rejected: incorrect number of final daughters: %ld (expected %ld)", arrAllDaughtersIndex.size(), N);
            return -1;
            }
        }
        // Check that the daughter is in the list of final daughters.
        // (Check that the daughter is not a stepdaughter, i.e. particle pointing to the mother while not being its daughter.)
        bool isDaughterFound = false; // Is the index of this prong among the remaining expected indices of daughters?
        for (std::size_t iD = 0; iD < arrAllDaughtersIndex.size(); ++iD) {
            if (arrDaughtersIndex[iProng] == arrAllDaughtersIndex[iD]) {
            arrAllDaughtersIndex[iD] = -1; // Remove this index from the array of expected daughters. (Rejects twin daughters, i.e. particle considered twice as a daughter.)
            isDaughterFound = true;
            break;
            }
        }
        if (!isDaughterFound) {
            // Printf("MC Rec: Rejected: bad daughter index: %d not in the list of final daughters", arrDaughtersIndex[iProng]);
            return -1;
        }
        // Check daughter's PDG code.
        auto pdgParticleI = particleI.pdgCode(); // PDG code of the ith daughter
        // Printf("MC Rec: Daughter %d PDG: %d", iProng, pdgParticleI);
        bool isPdgFound = false; // Is the PDG code of this daughter among the remaining expected PDG codes?
        for (std::size_t iProngCp = 0; iProngCp < N; ++iProngCp) {
            if (pdgParticleI == coefFlavourOscillation * sgn * arrPdgDaughters[iProngCp]) {
            arrPdgDaughters[iProngCp] = 0; // Remove this PDG code from the array of expected ones.
            isPdgFound = true;
            break;
            }
        }
        if (!isPdgFound) {
            // Printf("MC Rec: Rejected: bad daughter PDG: %d", pdgParticleI);
            return -1;
        }
        }
        // Printf("MC Rec: Accepted: m: %d", indexMother);
        if (sign) {
        *sign = sgn;
        }
        if constexpr (acceptTrackDecay) {
        if (nPiToMu) {
            *nPiToMu = nPiToMuLocal;
        }
        if (nKaToPi) {
            *nKaToPi = nKaToPiLocal;
        }
        }
        if constexpr (acceptTrackIntWithMaterial) {
        if (nInteractionsWithMaterial) {
            *nInteractionsWithMaterial = nInteractionsWithMaterialLocal;
        }
        }
        return indexMother;
    }

    std::vector<int64_t> idxPions{}, idxKaons{};
    std::vector<int> pdgPions{}, pdgKaons{};

    void process(o2::aod::McCollisions const& mcCollisions,
                 o2::aod::McParticles const& mcParticles
                //  o2::aod::Collisions const& collisions,
                //  o2::soa::Join<o2::aod::Tracks, aod::McTrackLabels, aod::TrackSelection> const& tracksWithMcLabels
    ) {
        
        for (const auto& collision : mcCollisions) {
            std::cout << "------------------------- Collision " << collision.globalIndex() << "/" << mcCollisions.size() << " -------------------------" << std::endl;
            auto particlesCollision = mcParticles.sliceBy(mcParticlesPerColl, collision.globalIndex());
            for (const auto& particle : particlesCollision) {
                if (!particle.has_mothers()) { // both functions would return -1
                    continue;
                }
                if (std::abs(particle.mothersIds().back() - particle.mothersIds().front()) > 1) {
                    std::cout << "Particle with more than 2 mothers";
                    bool isFaulty = false;
                    for (auto iMother=particle.mothersIds().front(); iMother < particle.mothersIds().back(); ++iMother) {
                        auto mother = mcParticles.rawIteratorAt(iMother - mcParticles.offset());
                        if (particle.globalIndex() < mother.daughtersIds().front() || particle.globalIndex() > mother.daughtersIds().back()) {
                            isFaulty = true;
                        }
                    }
                    if (isFaulty) {
                        std::cout << " is a faulty particle";
                    }
                    std::cout << std::endl;
                }
                std::cout << "Index: " << particle.globalIndex() <<
                             ", PDG: " << particle.pdgCode() << ", Mothers: ";
                for (const auto& iMother : particle.mothersIds()) {
                    std::cout << iMother << " ";
                }
                std::cout << ", Daughters: ";
                for (auto iDaughter = particle.daughtersIds().front(); iDaughter <= particle.daughtersIds().back(); ++iDaughter) {
                    std::cout << iDaughter << " ";
                }
                std::cout << std::endl;
                if (RecoDecay::getMother(mcParticles, particle, 411, true, nullptr, 2, false) != -1) {
                    registry.fill(HIST("hMother"), 0);
                }
                if (getMotherFixed(mcParticles, particle, 411, true, nullptr, 2) != -1) {
                    registry.fill(HIST("hMother"), 1);
                }
            }
        }
        // for (const auto& collision : collisions) {
        //     std::cout << "Processing collision " << collision.globalIndex() << "/" << collisions.size() << std::endl;
        //     const auto groupedTracks = tracksWithMcLabels.sliceBy(tracksPerColl, collision.globalIndex());

        //     idxPions.clear(); pdgPions.clear(); idxKaons.clear(); pdgKaons.clear();
        //     int64_t idx = 0;
        //     for (const auto& track : groupedTracks) {
        //         if (track.has_mcParticle() && track.isGlobalTrackWoDCA()) {
        //             const int pdg = track.mcParticle().pdgCode();
        //             if (std::abs(pdg) == 211) { idxPions.push_back(idx); pdgPions.push_back(pdg); }
        //             else if (std::abs(pdg) == 321) { idxKaons.push_back(idx); pdgKaons.push_back(pdg); }
        //         }
        //         ++idx;
        //     }
        //     if (idxKaons.empty() || idxPions.size() < 2) {
        //         continue;
        //     }
        //     std::cout << "Found " << idxPions.size() << " pions and " << idxKaons.size() << " kaons in collision " << collision.globalIndex() << std::endl;

        //     for (std::size_t a = 0; a < idxPions.size(); ++a) {
        //         const auto trackI = groupedTracks.rawIteratorAt(idxPions[a]);
        //         for (std::size_t b = a + 1; b < idxPions.size(); ++b) {
        //             if (pdgPions[b] != pdgPions[a]) { // the two pions carry the same charge
        //                 continue;
        //             }
        //             const auto trackK = groupedTracks.rawIteratorAt(idxPions[b]);
        //             for (std::size_t c = 0; c < idxKaons.size(); ++c) {
        //                 if (pdgKaons[c] * pdgPions[a] > 0) { // kaon must be opposite-sign
        //                     continue;
        //                 }
        //                 const auto trackJ = groupedTracks.rawIteratorAt(idxKaons[c]);
        //                 const std::array arrayDaughters{trackI, trackJ, trackK};
        //                 // depth 2, as in candidateCreator3Prong: D+ -> K pi pi also proceeds via
        //                 // resonances (K*0bar pi+, ...), which depth 1 would reject outright
        //                 if (RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, 411, {211, -321, 211}, true, nullptr, 2) != -1) {
        //                     registry.fill(HIST("hReco"), 0);
        //                 }
        //                 if (getMatchedMCRecFixed(mcParticles, arrayDaughters, 411, {211, -321, 211}, true, nullptr, 2) != -1) {
        //                     registry.fill(HIST("hReco"), 1);
        //                 }
        //             }
        //         }
        //     }
        // }
    }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<getMothers>(cfgc)}; }