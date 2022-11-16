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

// *****revision history*****:
//
// added recalculation of the conversion point on 08.07.22 by Nikita Philip Tatsch (tatsch@physi.uni-heidelberg.de)
//
// **************************

// *****revision history*****:
//
// added recalculation of the conversion point on 08.07.22 by Nikita Philip Tatsch (tatsch@physi.uni-heidelberg.de)
//
// **************************

// runme like: o2-analysis-trackselection -b --aod-file ${sourceFile} --aod-writer-json ${writerFile} | o2-analysis-timestamp -b | o2-analysis-trackextension -b | o2-analysis-lf-lambdakzerobuilder -b | o2-analysis-pid-tpc -b | o2-analysis-em-skimmermc -b

// todo: remove reduantant information in GammaConversionsInfoTrue
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/gammaConvDefinitions.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

//includes for the R recalculation
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include <CCDB/BasicCCDBManager.h>

#include "DetectorsVertexing/HelixHelper.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Common/Core/trackUtilities.h"

#include <TMath.h> // for ATan2, Cos, Sin, Sqrt
#include "TVector2.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using collisionEvSelIt = soa::Join<aod::Collisions, aod::EvSels>::iterator;
using tracksAndTPCInfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCEl, aod::pidTPCPi>;
using tracksAndTPCInfoMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCEl, aod::pidTPCPi, aod::McTrackLabels>;

struct skimmerGammaConversions {

  //configurables for CCDB access
  Configurable<std::string> path{"ccdb-path", "GLO/GRP/GRP", "path to the ccdb object"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<long> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      gHistoSpec_hCollisionZ_all_MCTrue,
      gHistoSpec_hCollisionZ_MCTrue,
      gHistoSpec_hCollisionZ_MCRec,
      {"hCollisionZ_Rec", "hCollisionZ_Rec;z (cm);counts", {HistType::kTH1F, {gAxis_zColl}}},
      {"hMcParticlesSize", "hMcParticlesSize", {HistType::kTH1F, {{100, 0.f, 1000000.f}}}},
      {"hV0Confirmation", "hV0Confirmation", {HistType::kTH1I, {{10, 0.f, 10.f}}}},
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
  Produces<aod::V0Recalculated> fFuncTableV0Recalculated;
  Produces<aod::V0DaughterMcParticles> fFuncTableMCTrackInformation;
  Produces<aod::MCParticleIndex> fIndexTableMCTrackIndex;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int runNumber = -1;

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

    // This is added in order to access the ccdb

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking(); // no idea wether this is useful or not, there is no documentation
    // Not later than now, will be replaced by the value of the train creation
    // This avoids that users can replace objects **while** a train is running
    ccdb->setCreatedNotAfter(nolaterthan.value); // was like that in the tutorial efficiencyPerRun
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    // if run number matches then magnetic field should also match. Avoids unessecary acceses.
    if (runNumber == bc.runNumber()) {
      return;
    }
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(path.value, bc.timestamp());
    if (!grpo) {
      LOGF(fatal, "Efficiency object not found!");
    }
    o2::base::Propagator::initFieldFromGRP(grpo);
    //o2::base::Propagator::Instance()->setMatLUT(lut);
    runNumber = bc.runNumber();
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

  template <typename TV0>
  void fillV0RecalculatedTable(TV0 const& theV0, float* recalculatedVtx)
  {
    fFuncTableV0Recalculated(
      recalculatedVtx[0],
      recalculatedVtx[1],
      recalculatedVtx[2]);
  }

  template <typename TTRACK>
  void fillfFuncTableMCTrackInformation(TTRACK theTrack, bool sameMother)
  {
    fFuncTableMCTrackInformation(
      theTrack.mcParticle().pdgCode(),
      theTrack.mcParticle().px(),
      theTrack.mcParticle().py(),
      theTrack.mcParticle().pz(),
      sameMother);
  }

  // ============================ FUNCTION DEFINITIONS ====================================================

  void processRec(aod::Collisions::iterator const& theCollision,
                  aod::BCsWithTimestamps const& bcs,
                  aod::V0Datas const& theV0s,
                  tracksAndTPCInfo const& theTracks)
  {
    // skip if bc has no Collisions
    if (theCollision.size() == 0) {
      return;
    }

    initCCDB(bcs.begin());

    fRegistry.fill(HIST("hCollisionZ_Rec"), theCollision.posZ());

    for (auto& lV0 : theV0s) {

      auto lTrackPos = lV0.template posTrack_as<tracksAndTPCInfo>(); // positive daughter
      auto lTrackNeg = lV0.template negTrack_as<tracksAndTPCInfo>(); // negative daughter

      float recalculatedVtx[3];
      Vtx_recalculation(lTrackPos, lTrackNeg, recalculatedVtx);

      fillTrackTable(lV0, lTrackPos, true);
      fillTrackTable(lV0, lTrackNeg, false);
      fillV0RecalculatedTable(lV0, recalculatedVtx);
    }
  }
  PROCESS_SWITCH(skimmerGammaConversions, processRec, "process reconstructed info only", true);

  Preslice<aod::V0Datas> perCollision = aod::v0data::collisionId;

  void processMc(aod::McCollision const& theMcCollision,
                 soa::SmallGroups<soa::Join<aod::McCollisionLabels,
                                            aod::Collisions>> const& theCollisions,
                 aod::BCsWithTimestamps const& bcs,
                 aod::V0Datas const& theV0s,
                 tracksAndTPCInfoMC const& theTracks,
                 aod::McParticles const& theMcParticles)
  {

    initCCDB(bcs.begin());

    fRegistry.fill(HIST("hCollisionZ_all_MCTrue"), theMcCollision.posZ());

    if (theCollisions.size() == 0) {
      return;
    }

    fRegistry.fill(HIST("hCollisionZ_MCTrue"), theMcCollision.posZ());
    fRegistry.fill(HIST("hMcParticlesSize"), theMcParticles.size());

    for (auto& lCollision : theCollisions) {
      fRegistry.fill(HIST("hCollisionZ_MCRec"), lCollision.posZ());

      auto lGroupedV0s = theV0s.sliceBy(perCollision, lCollision.globalIndex());
      for (auto& lV0 : lGroupedV0s) {

        auto lTrackPos = lV0.template posTrack_as<tracksAndTPCInfoMC>(); // positive daughter
        auto lTrackNeg = lV0.template negTrack_as<tracksAndTPCInfoMC>(); // negative daughter

        eV0Confirmation lV0Status = isTrueV0(lV0,
                                             lTrackPos,
                                             lTrackNeg);

        fRegistry.get<TH1>(HIST("hV0Confirmation"))->Fill(lV0Status);

        float recalculatedVtx[3];
        Vtx_recalculation(lTrackPos, lTrackNeg, recalculatedVtx);

        fillTrackTable(lV0, lTrackPos, true);
        fillTrackTable(lV0, lTrackNeg, false);
        fillV0RecalculatedTable(lV0, recalculatedVtx);
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

    int theReturnReason = -1;
    bool hasSameMother = false;
    bool MCTrackInformationHasEntry = false;

    int lPosEntryInMCTrack = -1;
    int lNegEntryInMCTrack = -1;

    // none of tracks has a mother, has been accounted for in fMotherSizesHisto
    if ((lMothersIndecesPos.size() + lMothersIndecesNeg.size()) == 0) {
      theReturnReason = kNoTrackComesFromMother;
    }

    // exactly one track has one mother
    if ((lMothersIndecesPos.size() + lMothersIndecesNeg.size()) == 1) {
      theReturnReason = kOneTrackHasOneMother;
    }

    // we know now both tracks have at least one mother
    // check if it is the same
    if (lMothersIndecesPos[0] != lMothersIndecesNeg[0]) {
      // fill Track Mc true table
      hasSameMother = false;
      MCTrackInformationHasEntry = true;
      theReturnReason = kNotSameMothers;
    } else {
      // fill Track Mc true table
      hasSameMother = true;
      MCTrackInformationHasEntry = true;
    }

    if (MCTrackInformationHasEntry) {
      fillfFuncTableMCTrackInformation(theTrackPos, hasSameMother);
      lPosEntryInMCTrack = fFuncTableMCTrackInformation.lastIndex();
      fillfFuncTableMCTrackInformation(theTrackNeg, hasSameMother);
      lNegEntryInMCTrack = fFuncTableMCTrackInformation.lastIndex();
    }

    fIndexTableMCTrackIndex(lPosEntryInMCTrack);
    fIndexTableMCTrackIndex(lNegEntryInMCTrack);

    if ((theReturnReason == kNoTrackComesFromMother) || (theReturnReason == kOneTrackHasOneMother) || (theReturnReason == kNotSameMothers)) {
      return static_cast<eV0Confirmation>(theReturnReason);
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
        lV0Radius,
        -1, -1);
      break; // because we only want to look at the first mother. If there are more it will show up in fMotherSizesHisto
    }
    return kGoodMcMother;
  }

  template <typename TrackPrecision = float, typename T>
  void Vtx_recalculation(T lTrackPos, T lTrackNeg, float* conversionPosition)
  {
    o2::base::Propagator* prop = o2::base::Propagator::Instance(); //This singleton propagator requires some initialisation of the CCDB object.

    float bz = prop->getNominalBz();

    //*******************************************************

    // o2::track::TrackParametrization<TrackPrecision> = TrackPar, I use the full version to have control over the data type
    o2::track::TrackParametrization<TrackPrecision> trackPosInformation = getTrackPar(lTrackPos); //first get an object that stores Track information (positive)
    o2::track::TrackParametrization<TrackPrecision> trackNegInformation = getTrackPar(lTrackNeg); //first get an object that stores Track information (negative)

    o2::track::TrackAuxPar helixPos(trackPosInformation, bz); //This object is a descendant of a CircleXY and stores cirlce information with respect to the magnetic field. This object uses functions and information of the o2::track::TrackParametrizationWithError<TrackPrecision> object (positive)
    o2::track::TrackAuxPar helixNeg(trackNegInformation, bz); //This object is a descendant of a CircleXY and stores cirlce information with respect to the magnetic field. This object uses functions and information of the o2::track::TrackParametrizationWithError<TrackPrecision> object (negative)

    conversionPosition[0] = (helixPos.xC * helixNeg.rC + helixNeg.xC * helixPos.rC) / (helixPos.rC + helixNeg.rC); //This calculates the coordinates of the conversion point as an weighted average of the two helix centers. xC and yC should be the global coordinates for the helix center as far as I understand. But you can double check the code of trackPosInformation.getCircleParamsLoc
    conversionPosition[1] = (helixPos.yC * helixNeg.rC + helixNeg.yC * helixPos.rC) / (helixPos.rC + helixNeg.rC); //If this calculation doesn't work check if the rotateZ function, because the "documentation" says I get global coordinates but maybe i don't.

    //I am unsure about the Z calculation but this is how it is done in AliPhysics as far as I understand
    o2::track::TrackParametrization<TrackPrecision> trackPosInformationCopy = o2::track::TrackParametrization<TrackPrecision>(trackPosInformation);
    o2::track::TrackParametrization<TrackPrecision> trackNegInformationCopy = o2::track::TrackParametrization<TrackPrecision>(trackNegInformation);

    //I think this calculation gets the closest point on the track to the conversion point
    //This alpha is a different alpha than the usual alpha and I think it is the angle between X axis and conversion point
    Double_t alphaPos = TMath::Pi() + TMath::ATan2(-(conversionPosition[1] - helixPos.yC), (conversionPosition[0] - helixPos.xC));
    Double_t alphaNeg = TMath::Pi() + TMath::ATan2(-(conversionPosition[1] - helixNeg.yC), (conversionPosition[0] - helixNeg.xC));

    Double_t vertexXPos = helixPos.xC + helixPos.rC * TMath::Cos(alphaPos);
    Double_t vertexYPos = helixPos.yC + helixPos.rC * TMath::Sin(alphaPos);
    Double_t vertexXNeg = helixNeg.xC + helixNeg.rC * TMath::Cos(alphaNeg);
    Double_t vertexYNeg = helixNeg.yC + helixNeg.rC * TMath::Sin(alphaNeg);

    TVector2 vertexPos(vertexXPos, vertexYPos);
    TVector2 vertexNeg(vertexXNeg, vertexYNeg);

    // Convert to local coordinate system
    TVector2 vertexPosRot = vertexPos.Rotate(-trackPosInformationCopy.getAlpha());
    TVector2 vertexNegRot = vertexNeg.Rotate(-trackNegInformationCopy.getAlpha());

    prop->propagateToX(trackPosInformationCopy,
                       vertexPosRot.X(),
                       bz,
                       o2::base::PropagatorImpl<TrackPrecision>::MAX_SIN_PHI,
                       o2::base::PropagatorImpl<TrackPrecision>::MAX_STEP,
                       o2::base::PropagatorImpl<TrackPrecision>::MatCorrType::USEMatCorrNONE);

    prop->propagateToX(trackNegInformationCopy,
                       vertexNegRot.X(),
                       bz,
                       o2::base::PropagatorImpl<TrackPrecision>::MAX_SIN_PHI,
                       o2::base::PropagatorImpl<TrackPrecision>::MAX_STEP,
                       o2::base::PropagatorImpl<TrackPrecision>::MatCorrType::USEMatCorrNONE);

    // TODO: This is still off and needs to be checked...
    conversionPosition[2] = (trackPosInformationCopy.getZ() * helixNeg.rC + trackNegInformationCopy.getZ() * helixPos.rC) / (helixPos.rC + helixNeg.rC);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerGammaConversions>(cfgc)};
}
