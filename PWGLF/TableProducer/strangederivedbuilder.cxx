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
//__________________________________________________
// this task provides general links between collisions
// and strange objects reconstructed in various ways.
// It is meant to help with providing auxiliary information
// when dealing with derived data.

#include <cmath>
#include <array>
#include <cstdlib>
#include <map>
#include <iterator>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Qvectors.h"
#include "Framework/StaticFor.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using TracksWithExtra = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullHe>;
using TracksCompleteIUMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::McTrackLabels>;
using FullTracksExtIUTOF = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TOFEvTime, aod::TOFSignal>;

// simple checkers
#define bitset(var, nbit) ((var) |= (1 << (nbit)))
#define bitcheck(var, nbit) ((var) & (1 << (nbit)))

struct strangederivedbuilder {
  //__________________________________________________
  // fundamental building blocks of derived data
  Produces<aod::StraCollision> strangeColl;      // characterises collisions
  Produces<aod::StraCents> strangeCents;         // characterises collisions / centrality
  Produces<aod::StraRawCents> strangeRawCents;   // characterises collisions / centrality
  Produces<aod::StraEvSels> strangeEvSels;       // characterises collisions / sel8 selection
  Produces<aod::StraStamps> strangeStamps;       // provides timestamps, run numbers
  Produces<aod::V0CollRefs> v0collref;           // references collisions from V0s
  Produces<aod::CascCollRefs> casccollref;       // references collisions from cascades
  Produces<aod::KFCascCollRefs> kfcasccollref;   // references collisions from KF cascades
  Produces<aod::TraCascCollRefs> tracasccollref; // references collisions from tracked cascades

  //__________________________________________________
  // track extra references
  Produces<aod::DauTrackExtras> dauTrackExtras;   // daughter track detector properties
  Produces<aod::DauTrackTPCPIDs> dauTrackTPCPIDs; // daughter track TPC PID
  Produces<aod::V0Extras> v0Extras;               // references DauTracks from V0s
  Produces<aod::CascExtras> cascExtras;           // references DauTracks from cascades
  Produces<aod::StraTrackExtras> straTrackExtras; // references DauTracks from tracked cascades

  //__________________________________________________
  // cascade interlinks
  Produces<aod::CascToTraRefs> cascToTraRefs; // cascades -> tracked
  Produces<aod::CascToKFRefs> cascToKFRefs;   // cascades -> KF
  Produces<aod::TraToCascRefs> traToCascRefs; // tracked -> cascades
  Produces<aod::KFToCascRefs> kfToCascRefs;   // KF -> cascades

  //__________________________________________________
  // mother information
  Produces<aod::V0MCMothers> v0mothers;       // V0 mother references
  Produces<aod::CascMCMothers> cascmothers;   // casc mother references
  Produces<aod::MotherMCParts> motherMCParts; // mc particles for mothers

  //__________________________________________________
  // raw TOF PID for posterior use if requested
  Produces<aod::V0TOFs> v0tofs;     // V0 part
  Produces<aod::CascTOFs> casctofs; // cascade part

  //__________________________________________________
  // Q-vectors
  Produces<aod::StraFT0AQVs> StraFT0AQVs; // FT0A Q-vector
  Produces<aod::StraFT0CQVs> StraFT0CQVs; // FT0C Q-vector
  Produces<aod::StraFT0MQVs> StraFT0MQVs; // FT0M Q-vector
  Produces<aod::StraFV0AQVs> StraFV0AQVs; // FV0A Q-vector

  // histogram registry for bookkeeping
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr int nSpecies = 14;
  static constexpr int nParameters = 1;
  static const std::vector<std::string> particleNames;
  static const std::vector<int> particlePDGCodes;
  static const std::vector<std::string> parameterNames;
  static const int defaultParameters[nSpecies][nParameters];
  static constexpr std::string_view particleNamesConstExpr[] = {"Gamma", "K0Short", "Lambda", "AntiLambda",
                                                                "Sigma0", "AntiSigma0", "SigmaPlus", "SigmaMinus",
                                                                "Hypertriton", "AntiHypertriton",
                                                                "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};

  uint32_t enabledBits = 0;

  Configurable<LabeledArray<int>> enableGeneratedInfo{"enableGeneratedInfo",
                                                      {defaultParameters[0], nSpecies,
                                                       nParameters, particleNames, parameterNames},
                                                      "Fill generated particle histograms for each species. 0: no, 1: yes"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisCentrality{"axisCentrality", {100, 0.0f, 100.0f}, "Centrality"};
  ConfigurableAxis axisNVertices{"axisNVertices", {10, -0.5f, 9.5f}, "N(vertices)"};

  Configurable<bool> fillEmptyCollisions{"fillEmptyCollisions", false, "fill collision entries without candidates"};

  // round Nsigma variables up to a certain level of precision if requested
  // useful to keep derived data sizes under control
  // variables that are rounded include the DCAs but not the CosPA (precision needed)
  Configurable<bool> roundNSigmaVariables{"roundNSigmaVariables", false, "round NSigma variables"};
  Configurable<float> precisionNSigmas{"precisionNSigmas", 0.1f, "precision to keep NSigmas"};

  Configurable<bool> fillRawFT0A{"fillRawFT0A", false, "Fill raw FT0A information for debug"};
  Configurable<bool> fillRawFT0C{"fillRawFT0C", true, "Fill raw FT0C information for debug"};
  Configurable<bool> fillRawFV0A{"fillRawFV0A", false, "Fill raw FV0A information for debug"};
  Configurable<bool> fillRawNTracksEta1{"fillRawNTracksEta1", true, "Fill raw NTracks |eta|<1 information for debug"};

  // For manual sliceBy
  Preslice<aod::V0Datas> V0perCollision = o2::aod::v0data::collisionId;
  Preslice<aod::CascDatas> CascperCollision = o2::aod::cascdata::collisionId;
  Preslice<aod::KFCascDatas> KFCascperCollision = o2::aod::cascdata::collisionId;
  Preslice<aod::TraCascDatas> TraCascperCollision = o2::aod::cascdata::collisionId;

  int64_t currentCollIdx;

  float roundToPrecision(float number, float step = 0.01)
  {
    // this function rounds a certain number in an axis that is quantized by
    // the variable 'step'; the rounded number is placed halfway between
    // n*step and (n+1)*step such that analysis can be done with absolutely
    // no issue with precision 'step'.
    return step * static_cast<float>(static_cast<int>((number) / step)) + TMath::Sign(1.0f, number) * (0.5f) * step;
  }

  void init(InitContext& context)
  {
    currentCollIdx = -1;
    // setup map for fast checking if enabled
    static_for<0, nSpecies - 1>([&](auto i) {
      constexpr int index = i.value;
      int f = enableGeneratedInfo->get(particleNames[index].c_str(), "Enable");
      if (f == 1) {
        bitset(enabledBits, index);
      }
    });

    // Creation of histograms: MC generated
    for (Int_t i = 0; i < nSpecies; i++) {
      histos.add(Form("hGen%s", particleNames[i].data()), Form("hGen%s", particleNames[i].data()), kTH1D, {axisPt});
      histos.add(Form("h2dGen%s", particleNames[i].data()), Form("h2dGen%s", particleNames[i].data()), kTH2D, {axisCentrality, axisPt});
    }

    histos.add("h2dNVerticesVsCentrality", "h2dNVerticesVsCentrality", kTH2D, {axisCentrality, axisNVertices});
  }

  void processCollisionsV0sOnly(soa::Join<aod::Collisions, aod::FT0Mults, aod::FV0Mults, aod::PVMults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::EvSels> const& collisions, aod::V0Datas const& V0s, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisColl = V0s.sliceBy(V0perCollision, collIdx);
      bool strange = V0Table_thisColl.size() > 0;
      // casc table sliced
      if (strange || fillEmptyCollisions) {
        strangeColl(collision.posX(), collision.posY(), collision.posZ());
        strangeCents(collision.centFT0M(), collision.centFT0A(),
                     collision.centFT0C(), collision.centFV0A());
        strangeEvSels(collision.sel8());
        auto bc = collision.bc_as<aod::BCsWithTimestamps>();
        strangeStamps(bc.runNumber(), bc.timestamp());

        if (fillRawFT0C || fillRawFT0C || fillRawFV0A || fillRawNTracksEta1) {
          strangeRawCents(collision.multFT0A() * static_cast<float>(fillRawFT0A),
                          collision.multFT0C() * static_cast<float>(fillRawFT0C),
                          collision.multFT0A() * static_cast<float>(fillRawFV0A),
                          collision.multNTracksPVeta1() * static_cast<int>(fillRawNTracksEta1));
        }
      }
      for (int i = 0; i < V0Table_thisColl.size(); i++)
        v0collref(strangeColl.lastIndex());
    }
  }

  void processCollisions(soa::Join<aod::Collisions, aod::FT0Mults, aod::FV0Mults, aod::PVMults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::EvSels> const& collisions, aod::V0Datas const& V0s, aod::CascDatas const& Cascades, aod::KFCascDatas const& KFCascades, aod::TraCascDatas const& TraCascades, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisColl = V0s.sliceBy(V0perCollision, collIdx);
      auto CascTable_thisColl = Cascades.sliceBy(CascperCollision, collIdx);
      auto KFCascTable_thisColl = KFCascades.sliceBy(KFCascperCollision, collIdx);
      auto TraCascTable_thisColl = TraCascades.sliceBy(TraCascperCollision, collIdx);
      bool strange = V0Table_thisColl.size() > 0 ||
                     CascTable_thisColl.size() > 0 ||
                     KFCascTable_thisColl.size() > 0 ||
                     TraCascTable_thisColl.size() > 0;
      // casc table sliced
      if (strange || fillEmptyCollisions) {
        strangeColl(collision.posX(), collision.posY(), collision.posZ());
        strangeCents(collision.centFT0M(), collision.centFT0A(),
                     collision.centFT0C(), collision.centFV0A());
        strangeEvSels(collision.sel8());
        auto bc = collision.bc_as<aod::BCsWithTimestamps>();
        strangeStamps(bc.runNumber(), bc.timestamp());

        if (fillRawFT0C || fillRawFT0C || fillRawFV0A || fillRawNTracksEta1) {
          strangeRawCents(collision.multFT0A() * static_cast<float>(fillRawFT0A),
                          collision.multFT0C() * static_cast<float>(fillRawFT0C),
                          collision.multFT0A() * static_cast<float>(fillRawFV0A),
                          collision.multNTracksPVeta1() * static_cast<int>(fillRawNTracksEta1));
        }
      }
      for (int i = 0; i < V0Table_thisColl.size(); i++)
        v0collref(strangeColl.lastIndex());
      for (int i = 0; i < CascTable_thisColl.size(); i++)
        casccollref(strangeColl.lastIndex());
      for (int i = 0; i < KFCascTable_thisColl.size(); i++)
        kfcasccollref(strangeColl.lastIndex());
      for (int i = 0; i < TraCascTable_thisColl.size(); i++)
        tracasccollref(strangeColl.lastIndex());
    }
  }

  void processTrackExtrasV0sOnly(aod::V0Datas const& V0s, TracksWithExtra const& tracksExtra)
  {
    std::vector<int> trackMap(tracksExtra.size(), -1); // index -1: not used

    //__________________________________________________
    // mark tracks that belong to V0s
    for (auto const& v0 : V0s) {
      auto const& posTrack = v0.posTrack_as<TracksWithExtra>();
      auto const& negTrack = v0.negTrack_as<TracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
    }
    //__________________________________________________
    // Figure out the numbering of the new tracks table
    // assume filling per order
    int nTracks = 0;
    for (int i = 0; i < trackMap.size(); i++) {
      if (trackMap[i] >= 0) {
        trackMap[i] = nTracks++;
      }
    }
    //__________________________________________________
    // populate track references
    for (auto const& v0 : V0s) {
      auto const& posTrack = v0.posTrack_as<TracksWithExtra>();
      auto const& negTrack = v0.negTrack_as<TracksWithExtra>();
      v0Extras(trackMap[posTrack.globalIndex()],
               trackMap[negTrack.globalIndex()]); // joinable with V0Datas
    }
    //__________________________________________________
    // circle back and populate actual DauTrackExtra table
    for (auto const& tr : tracksExtra) {
      if (trackMap[tr.globalIndex()] >= 0) {
        dauTrackExtras(tr.detectorMap(), tr.itsClusterSizes(),
                       tr.tpcNClsFound(), tr.tpcNClsCrossedRows());
      }
    }
    // done!
  }

  void processTrackExtras(aod::V0Datas const& V0s, aod::CascDatas const& Cascades, aod::KFCascDatas const& KFCascades, aod::TraCascDatas const& TraCascades, TracksWithExtra const& tracksExtra, aod::V0s const&)
  {
    std::vector<int> trackMap(tracksExtra.size(), -1); // index -1: not used

    //__________________________________________________
    // mark tracks that belong to V0s
    for (auto const& v0 : V0s) {
      auto const& posTrack = v0.posTrack_as<TracksWithExtra>();
      auto const& negTrack = v0.negTrack_as<TracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
    }

    //__________________________________________________
    // index tracks that belong to CascDatas
    for (auto const& casc : Cascades) {
      auto bachTrack = casc.bachelor_as<TracksWithExtra>();
      auto posTrack = casc.posTrack_as<TracksWithExtra>();
      auto negTrack = casc.negTrack_as<TracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
      trackMap[bachTrack.globalIndex()] = 0;
    }
    //__________________________________________________
    // index tracks that belong to KFCascDatas
    for (auto const& casc : KFCascades) {
      auto bachTrack = casc.bachelor_as<TracksWithExtra>();
      auto posTrack = casc.posTrack_as<TracksWithExtra>();
      auto negTrack = casc.negTrack_as<TracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
      trackMap[bachTrack.globalIndex()] = 0;
    }
    //__________________________________________________
    // index tracks that belong to TraCascDatas
    for (auto const& casc : TraCascades) {
      auto bachTrack = casc.bachelor_as<TracksWithExtra>();
      auto posTrack = casc.posTrack_as<TracksWithExtra>();
      auto negTrack = casc.negTrack_as<TracksWithExtra>();
      auto strangeTrack = casc.strangeTrack_as<TracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
      trackMap[bachTrack.globalIndex()] = 0;
      trackMap[strangeTrack.globalIndex()] = 0;
    }
    //__________________________________________________
    // Figure out the numbering of the new tracks table
    // assume filling per order
    int nTracks = 0;
    for (int i = 0; i < trackMap.size(); i++) {
      if (trackMap[i] >= 0) {
        trackMap[i] = nTracks++;
      }
    }
    //__________________________________________________
    // populate track references
    for (auto const& v0 : V0s) {
      auto const& posTrack = v0.posTrack_as<TracksWithExtra>();
      auto const& negTrack = v0.negTrack_as<TracksWithExtra>();
      v0Extras(trackMap[posTrack.globalIndex()],
               trackMap[negTrack.globalIndex()]); // joinable with V0Datas
    }
    //__________________________________________________
    // populate track references
    for (auto const& casc : Cascades) {
      auto bachTrack = casc.bachelor_as<TracksWithExtra>();
      auto posTrack = casc.posTrack_as<TracksWithExtra>();
      auto negTrack = casc.negTrack_as<TracksWithExtra>();
      cascExtras(trackMap[posTrack.globalIndex()],
                 trackMap[negTrack.globalIndex()],
                 trackMap[bachTrack.globalIndex()]); // joinable with CascDatas
    }
    //__________________________________________________
    // populate track references
    for (auto const& casc : TraCascades) {
      auto strangeTrack = casc.strangeTrack_as<TracksWithExtra>();
      straTrackExtras(trackMap[strangeTrack.globalIndex()]); // joinable with TraCascDatas
    }
    //__________________________________________________
    // circle back and populate actual DauTrackExtra table
    for (auto const& tr : tracksExtra) {
      if (trackMap[tr.globalIndex()] >= 0) {
        dauTrackExtras(tr.detectorMap(), tr.itsClusterSizes(),
                       tr.tpcNClsFound(), tr.tpcNClsCrossedRows());

        // round if requested
        if (roundNSigmaVariables) {
          dauTrackTPCPIDs(tr.tpcSignal(),
                          roundToPrecision(tr.tpcNSigmaEl(), precisionNSigmas),
                          roundToPrecision(tr.tpcNSigmaPi(), precisionNSigmas),
                          roundToPrecision(tr.tpcNSigmaKa(), precisionNSigmas),
                          roundToPrecision(tr.tpcNSigmaPr(), precisionNSigmas),
                          roundToPrecision(tr.tpcNSigmaHe(), precisionNSigmas));
        } else {
          dauTrackTPCPIDs(tr.tpcSignal(), tr.tpcNSigmaEl(),
                          tr.tpcNSigmaPi(), tr.tpcNSigmaKa(),
                          tr.tpcNSigmaPr(), tr.tpcNSigmaHe());
        }
      }
    }
    // done!
  }

  void processStrangeMothers(soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s, soa::Join<aod::CascDatas, aod::McCascLabels> const& Cascades, aod::McParticles const& mcParticles)
  {
    std::vector<int> motherReference(mcParticles.size(), -1); // index -1: not used / no reference

    //__________________________________________________
    // mark mcParticles for referencing
    for (auto const& v0 : V0s)
      if (v0.has_mcMotherParticle())
        motherReference[v0.mcMotherParticleId()] = 0;
    for (auto const& ca : Cascades)
      if (ca.has_mcMotherParticle())
        motherReference[ca.mcMotherParticleId()] = 0;
    //__________________________________________________
    // Figure out the numbering of the new mcMother table
    // assume filling per order
    int nParticles = 0;
    for (int i = 0; i < motherReference.size(); i++) {
      if (motherReference[i] >= 0) {
        motherReference[i] = nParticles++; // count particles of interest
      }
    }
    //__________________________________________________
    // populate track references
    for (auto const& v0 : V0s)
      v0mothers(motherReference[v0.mcMotherParticleId()]); // joinable with V0Datas
    for (auto const& ca : Cascades)
      cascmothers(motherReference[ca.mcMotherParticleId()]); // joinable with CascDatas
    //__________________________________________________
    // populate motherMCParticles
    for (auto const& tr : mcParticles) {
      if (motherReference[tr.globalIndex()] >= 0) {
        motherMCParts(tr.px(), tr.py(), tr.pz(), tr.pdgCode(), tr.isPhysicalPrimary());
      }
    }
  }

  using interlinkedCascades = soa::Join<aod::Cascades, aod::CascDataLink, aod::KFCascDataLink, aod::TraCascDataLink>;

  void processCascadeInterlinkTracked(interlinkedCascades const& masterCascades, aod::CascIndices const& Cascades, aod::TraCascIndices const& TraCascades)
  {
    // Standard to tracked
    for (auto const& c : Cascades) {
      int indexTracked = -1;
      if (c.has_cascade()) {
        auto cascade = c.cascade_as<interlinkedCascades>();
        indexTracked = cascade.traCascDataId();
      }
      cascToTraRefs(indexTracked);
    }
    // Tracked to standard
    for (auto const& c : TraCascades) {
      int index = -1;
      if (c.has_cascade()) {
        auto cascade = c.cascade_as<interlinkedCascades>();
        index = cascade.cascDataId();
      }
      traToCascRefs(index);
    }
  }

  void processCascadeInterlinkKF(interlinkedCascades const& masterCascades, aod::CascIndices const& Cascades, aod::KFCascIndices const& KFCascades)
  {
    // Standard to KF
    for (auto const& c : Cascades) {
      int indexKF = -1;
      if (c.has_cascade()) {
        auto cascade = c.cascade_as<interlinkedCascades>();
        indexKF = cascade.kfCascDataId();
      }
      cascToKFRefs(indexKF);
    }
    // KF to standard
    for (auto const& c : KFCascades) {
      int index = -1;
      if (c.has_cascade()) {
        auto cascade = c.cascade_as<interlinkedCascades>();
        index = cascade.cascDataId();
      }
      kfToCascRefs(index);
    }
  }

  void processPureSimulation(aod::McParticles const& mcParticles)
  {
    for (auto& mcp : mcParticles) {
      if (TMath::Abs(mcp.y()) < 0.5) {
        static_for<0, nSpecies - 1>([&](auto i) {
          constexpr int index = i.value;
          if (mcp.pdgCode() == particlePDGCodes[index] && bitcheck(enabledBits, index)) {
            histos.fill(HIST("hGen") + HIST(particleNamesConstExpr[index]), mcp.pt());
          }
        });
      }
    }
  }

  void processReconstructedSimulation(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>> const& collisions, aod::McParticles const& mcParticles)
  {
    // this process function also checks if a given collision was reconstructed and checks explicitly for splitting, etc

    // identify best-of collision
    int biggestNContribs = -1;
    float bestCentrality = 100.5;
    for (auto& collision : collisions) {
      if (biggestNContribs < collision.numContrib()) {
        biggestNContribs = collision.numContrib();
        bestCentrality = collision.centFT0C();
      }
    }
    histos.fill(HIST("h2dNVerticesVsCentrality"), bestCentrality, collisions.size());

    for (auto& mcp : mcParticles) {
      if (TMath::Abs(mcp.y()) < 0.5 && mcp.isPhysicalPrimary()) {
        static_for<0, nSpecies - 1>([&](auto i) {
          constexpr int index = i.value;
          if (mcp.pdgCode() == particlePDGCodes[index] && bitcheck(enabledBits, index)) {
            histos.fill(HIST("h2dGen") + HIST(particleNamesConstExpr[index]), bestCentrality, mcp.pt());
          }
        });
      }
    }
  }

  void processProduceV0TOFs(aod::Collision const& collision, aod::V0Datas const& V0s, FullTracksExtIUTOF const&)
  {
    for (auto const& v0 : V0s) {
      auto const& posTrackRow = v0.posTrack_as<FullTracksExtIUTOF>();
      auto const& negTrackRow = v0.negTrack_as<FullTracksExtIUTOF>();
      v0tofs(posTrackRow.length(), negTrackRow.length(),
             posTrackRow.tofSignal(), negTrackRow.tofSignal(),
             posTrackRow.tofEvTime(), negTrackRow.tofEvTime());
    }
  }
  void processProduceCascTOFs(aod::Collision const& collision, aod::CascDatas const& Cascades, FullTracksExtIUTOF const&)
  {
    for (auto const& cascade : Cascades) {
      auto const& posTrackRow = cascade.posTrack_as<FullTracksExtIUTOF>();
      auto const& negTrackRow = cascade.negTrack_as<FullTracksExtIUTOF>();
      auto const& bachTrackRow = cascade.bachelor_as<FullTracksExtIUTOF>();
      casctofs(posTrackRow.length(), negTrackRow.length(), bachTrackRow.length(),
               posTrackRow.tofSignal(), negTrackRow.tofSignal(), bachTrackRow.tofSignal(),
               posTrackRow.tofEvTime(), negTrackRow.tofEvTime(), bachTrackRow.tofEvTime());
    }
  }

  void processFT0AQVectors(soa::Join<aod::Collisions, aod::QvectorFT0As>::iterator const& collision)
  {
    StraFT0AQVs(collision.qvecFT0ARe(), collision.qvecFT0AIm(), collision.sumAmplFT0A());
  }
  void processFT0CQVectors(soa::Join<aod::Collisions, aod::QvectorFT0Cs>::iterator const& collision)
  {
    StraFT0CQVs(collision.qvecFT0CRe(), collision.qvecFT0CIm(), collision.sumAmplFT0C());
  }
  void processFT0MQVectors(soa::Join<aod::Collisions, aod::QvectorFT0Ms>::iterator const& collision)
  {
    StraFT0MQVs(collision.qvecFT0MRe(), collision.qvecFT0MIm(), collision.sumAmplFT0M());
  }
  void processFV0AQVectors(soa::Join<aod::Collisions, aod::QvectorFV0As>::iterator const& collision)
  {
    StraFV0AQVs(collision.qvecFV0ARe(), collision.qvecFV0AIm(), collision.sumAmplFV0A());
  }

  PROCESS_SWITCH(strangederivedbuilder, processCollisionsV0sOnly, "Produce collisions (V0s only)", true);
  PROCESS_SWITCH(strangederivedbuilder, processCollisions, "Produce collisions (V0s + casc)", true);
  PROCESS_SWITCH(strangederivedbuilder, processTrackExtrasV0sOnly, "Produce track extra information (V0s only)", true);
  PROCESS_SWITCH(strangederivedbuilder, processTrackExtras, "Produce track extra information (V0s + casc)", true);
  PROCESS_SWITCH(strangederivedbuilder, processStrangeMothers, "Produce tables with mother info for V0s + casc", true);
  PROCESS_SWITCH(strangederivedbuilder, processCascadeInterlinkTracked, "Produce tables interconnecting cascades", false);
  PROCESS_SWITCH(strangederivedbuilder, processCascadeInterlinkKF, "Produce tables interconnecting cascades", false);
  PROCESS_SWITCH(strangederivedbuilder, processPureSimulation, "Produce pure simulated information", true);
  PROCESS_SWITCH(strangederivedbuilder, processReconstructedSimulation, "Produce reco-ed simulated information", true);
  PROCESS_SWITCH(strangederivedbuilder, processProduceV0TOFs, "Produce V0TOFs table", true);
  PROCESS_SWITCH(strangederivedbuilder, processProduceCascTOFs, "Produce CascTOFs table", true);
  PROCESS_SWITCH(strangederivedbuilder, processFT0AQVectors, "Produce FT0A Q-vectors table", false);
  PROCESS_SWITCH(strangederivedbuilder, processFT0CQVectors, "Produce FT0C Q-vectors table", false);
  PROCESS_SWITCH(strangederivedbuilder, processFT0MQVectors, "Produce FT0M Q-vectors table", false);
  PROCESS_SWITCH(strangederivedbuilder, processFV0AQVectors, "Produce FV0A Q-vectors table", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangederivedbuilder>(cfgc)};
}

//__________________________________________________
// do not over-populate general namespace, keep scope strangederivedbuilder::
const std::vector<std::string> strangederivedbuilder::particleNames{"Gamma", "K0Short", "Lambda", "AntiLambda",
                                                                    "Sigma0", "AntiSigma0", "SigmaPlus", "SigmaMinus",
                                                                    "Hypertriton", "AntiHypertriton",
                                                                    "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};
const std::vector<int> strangederivedbuilder::particlePDGCodes{22, 310, 3122, -3122, 3212, -3212, 3222, 3112,
                                                               1010010030, -1010010030, 3312, -3312, 3334, -3334};
const std::vector<std::string> strangederivedbuilder::parameterNames{"Enable"};

const int strangederivedbuilder::defaultParameters[strangederivedbuilder::nSpecies][strangederivedbuilder::nParameters] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}};
