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

/// \brief write relevant information for photon conversion analysis to an AO2D.root file. This file is then the only necessary input to perform
/// pcm analysis.
/// dependencies: o2-analysis-lf-lambdakzerobuilder
/// \author stephan.friedrich.stiefelmaier@cern.ch

// runme like: o2-analysis-trackselection -b --aod-file ${sourceFile} --aod-writer-json ${writerFile} | o2-analysis-timestamp -b | o2-analysis-trackextension -b | o2-analysis-lf-lambdakzerobuilder -b | o2-analysis-pid-tpc -b | o2-analysis-em-skimmermc -b

// todo: remove reduantant information in GammaConversionsInfoTrue
#include "gammaTables.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using collisionEvSelIt = soa::Join<aod::Collisions, aod::EvSels>::iterator;
using tracksAndTPCInfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::pidTPCEl, aod::pidTPCPi>;
using tracksAndTPCInfoMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::pidTPCEl, aod::pidTPCPi, aod::McTrackLabels>;

struct skimmerGammaConversions {

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      {"hCollisionZ_MCTrue_all", "hCollisionZ_MCTrue_all;z (cm);counts", {HistType::kTH1F, {{800, -50.f, 50.f}}}},
      {"hCollisionZ_MCTrue", "hCollisionZ_MCTrue (at least one rec. collision);z (cm);counts", {HistType::kTH1F, {{800, -50.f, 50.f}}}},
      {"hCollisionZ_MCRec", "hCollisionZ_MCRec;z (cm);counts", {HistType::kTH1F, {{800, -50.f, 50.f}}}},
      {"hCollisionZ_Rec", "hCollisionZ_Rec;z (cm);counts", {HistType::kTH1F, {{800, -50.f, 50.f}}}},
      {"hMcParticlesSize", "hMcParticlesSize", {HistType::kTH1F, {{100, 0.f, 1000000.f}}}},
      {"hV0Confirmation", "hV0Confirmation", {HistType::kTH1F, {{10, 0.f, 10.f}}}},
    },
  };

  // declare this here in order to be able to access it from a lambda
  std::shared_ptr<TH1> fMotherSizesHisto{};

  enum eV0Confirmation {
    kV0In,
    kNoMcParticleForTrack,
    kNoTrackComesFromMother,
    kOneTrackHasOneMother,
    kNotSameMothers,
    kMotherHasNoDaughter,
    kGoodMcMother
  };

  std::map<size_t, std::string> fV0ConfirmationLabel{
    {kV0In, "kV0In"},
    {kNoMcParticleForTrack, "kNoMcParticleForTrack"},
    {kNoTrackComesFromMother, "kNoTrackComesFromMother"},
    {kOneTrackHasOneMother, "kOneTrackHasOneMother"},
    {kNotSameMothers, "kNotSameMothers"},
    {kMotherHasNoDaughter, "kMotherHasNoDaughter"},
    {kGoodMcMother, "kGoodMcMother"}};

  Produces<aod::V0DaughterTracks> fFuncTableV0DaughterTracks;
  Produces<aod::McGammasTrue> fFuncTableMcGammasFromConfirmedV0s;

  void init(InitContext const&)
  {
    if (doprocessRec && doprocessMc) {
      LOGF(fatal, "Cannot enable doprocessRec and doprocessMc at the same time. Please choose one.");
    }

    fMotherSizesHisto = std::get<std::shared_ptr<TH1>>(fRegistry.add("hMotherSizes", "hMotherSizes", {HistType::kTH1F, {{10, 0.f, 10.f}}}));

    // set axis lables
    TAxis* lXaxis = fRegistry.get<TH1>(HIST("hV0Confirmation"))->GetXaxis();
    for (auto& lPairIt : fV0ConfirmationLabel) {
      lXaxis->SetBinLabel(lPairIt.first + 1, lPairIt.second.data());
    }
  }

  template <typename TV0, typename TTRACK>
  void fillTrackTable(TV0 const& theV0, TTRACK const& theTrack, bool theIsPositive)
  {
    fFuncTableV0DaughterTracks(
      theV0.v0Id(),
      theTrack.dcaXY(),
      theTrack.eta(),
      theTrack.p(),
      theTrack.phi(),
      theTrack.pt(),
      theIsPositive,
      theTrack.tpcCrossedRowsOverFindableCls(),
      theTrack.tpcFoundOverFindableCls(),
      theTrack.tpcNClsCrossedRows(),
      theTrack.tpcNSigmaEl(),
      theTrack.tpcNSigmaPi(),
      theTrack.tpcSignal());
  }

  // ============================ FUNCTION DEFINITIONS ====================================================
  void processRec(aod::Collisions::iterator const& theCollision,
                  aod::V0Datas const& theV0s,
                  tracksAndTPCInfo const& theTracks)
  {
    fRegistry.fill(HIST("hCollisionZ_Rec"), theCollision.posZ());

    for (auto& lV0 : theV0s) {

      auto lTrackPos = lV0.template posTrack_as<tracksAndTPCInfo>(); // positive daughter
      auto lTrackNeg = lV0.template negTrack_as<tracksAndTPCInfo>(); // negative daughter

      fillTrackTable(lV0, lTrackPos, true);
      fillTrackTable(lV0, lTrackNeg, false);
    }
  }
  PROCESS_SWITCH(skimmerGammaConversions, processRec, "process reconstructed info only", true);

  void processMc(aod::McCollision const& theMcCollision,
                 soa::SmallGroups<soa::Join<aod::McCollisionLabels,
                                            aod::Collisions>> const& theCollisions,
                 aod::V0Datas const& theV0s,
                 tracksAndTPCInfoMC const& theTracks,
                 aod::McParticles const& theMcParticles)
  {
    fRegistry.fill(HIST("hCollisionZ_MCTrue_all"), theMcCollision.posZ());

    if (theCollisions.size() == 0) {
      return;
    }

    fRegistry.fill(HIST("hCollisionZ_MCTrue"), theMcCollision.posZ());
    fRegistry.fill(HIST("hMcParticlesSize"), theMcParticles.size());

    for (auto& lCollision : theCollisions) {
      fRegistry.fill(HIST("hCollisionZ_MCRec"), lCollision.posZ());

      // todo: replace by sliceByCached
      auto lGroupedV0s = theV0s.sliceBy(aod::v0data::collisionId, lCollision.globalIndex());
      for (auto& lV0 : lGroupedV0s) {

        auto lTrackPos = lV0.template posTrack_as<tracksAndTPCInfoMC>(); // positive daughter
        auto lTrackNeg = lV0.template negTrack_as<tracksAndTPCInfoMC>(); // negative daughter

        eV0Confirmation lV0Status = isTrueV0(lV0,
                                             lTrackPos,
                                             lTrackNeg);

        fRegistry.get<TH1>(HIST("hV0Confirmation"))->Fill(lV0Status);

        fillTrackTable(lV0, lTrackPos, true);
        fillTrackTable(lV0, lTrackNeg, false);
      }
    }
  }
  PROCESS_SWITCH(skimmerGammaConversions, processMc, "process reconstructed and mc info ", false);

  template <typename TV0, typename TTRACK>
  eV0Confirmation isTrueV0(TV0 const& theV0,
                           TTRACK const& theTrackPos,
                           TTRACK const& theTrackNeg)
  {
    auto getMothersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lMothersIndeces{};
      for (auto& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        lMothersIndeces.push_back(lMother.globalIndex());
      }
      fMotherSizesHisto->Fill(0.5 + (float)lMothersIndeces.size());
      return lMothersIndeces;
    };

    fRegistry.get<TH1>(HIST("hV0Confirmation"))->Fill(kV0In);

    if (!(theTrackPos.has_mcParticle() && theTrackNeg.has_mcParticle())) {
      return kNoMcParticleForTrack;
    }

    // get mcParticles
    auto lMcPos = theTrackPos.mcParticle();
    auto lMcNeg = theTrackNeg.mcParticle();

    // get indeces of mcMother of tracks
    std::vector<int> lMothersIndecesPos = getMothersIndeces(lMcPos);
    std::vector<int> lMothersIndecesNeg = getMothersIndeces(lMcNeg);

    // none of tracks has a mother, has been accounted for in fMotherSizesHisto
    if ((lMothersIndecesPos.size() + lMothersIndecesNeg.size()) == 0) {
      return kNoTrackComesFromMother;
    }

    // exactly one track has one mother
    if ((lMothersIndecesPos.size() + lMothersIndecesNeg.size()) == 1) {
      return kOneTrackHasOneMother;
    }

    // we know now both tracks have at least one mother
    // check if it is the same
    if (lMothersIndecesPos[0] != lMothersIndecesNeg[0]) {
      return kNotSameMothers;
    }

    // both tracks have the same first mother
    // SFS todo: actually no loop required here, for this. Access directly mc particle with lMothersIndecesFlags
    for (auto& lMcMother : lMcPos.template mothers_as<aod::McParticles>()) {

      // get mc daughter in order to compute true conversion point
      auto lDaughters = lMcMother.template daughters_as<aod::McParticles>();
      if (lDaughters.begin() == lDaughters.end()) {
        // mc converted mother has no mc daughters, should never happen
        return kMotherHasNoDaughter;
      }
      // todo: lMcPos instead here should give identical results
      auto lDaughter0 = lDaughters.begin();
      float lDaughter0Vx = lDaughter0.vx();
      float lDaughter0Vy = lDaughter0.vy();
      float lDaughter0Vz = lDaughter0.vz();
      float lV0Radius = sqrt(pow(lDaughter0Vx, 2) + pow(lDaughter0Vy, 2));

      fFuncTableMcGammasFromConfirmedV0s(
        lMcMother.mcCollisionId(),
        lMcMother.globalIndex(),
        theV0.v0Id(),
        lMcMother.pdgCode(), lMcMother.statusCode(), lMcMother.flags(),
        lMcMother.px(), lMcMother.py(), lMcMother.pz(),
        lMcMother.vx(), lMcMother.vy(), lMcMother.vz(), lMcMother.vt(),
        lDaughters.size(),
        lMcMother.eta(), lMcMother.phi(), lMcMother.p(), lMcMother.pt(), lMcMother.y(),
        lDaughter0Vx, lDaughter0Vy, lDaughter0Vz,
        lV0Radius);
      break; // because we only want to look at the first mother. If there are more it will show up in fMotherSizesHisto
    }
    return kGoodMcMother;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerGammaConversions>(cfgc)};
}
