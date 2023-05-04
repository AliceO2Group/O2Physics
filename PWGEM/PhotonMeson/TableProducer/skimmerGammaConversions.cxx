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
// adding accesing to ccdb objects for 2022 data taking on 30.11.22 by A. Marin (a.marin@gsi.de)
//
// **************************

// runme like: o2-analysis-trackselection -b --aod-file ${sourceFile} --aod-writer-json ${writerFile} | o2-analysis-timestamp -b | o2-analysis-trackextension -b | o2-analysis-lf-lambdakzerobuilder -b | o2-analysis-pid-tpc -b | o2-analysis-em-skimmermc -b

// todo: remove reduantant information in GammaConversionsInfoTrue
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/gammaConvDefinitions.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

// includes for the R recalculation
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "DCAFitter/HelixHelper.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Common/Core/trackUtilities.h"

#include <TMath.h>
#include <TVector2.h>

#include "Tools/KFparticle/KFUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using tracksAndTPCInfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::TracksCov>;
using tracksAndTPCInfo = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::TracksCovIU>;
using tracksAndTPCInfoMC = soa::Join<tracksAndTPCInfo, aod::McTrackLabels>;

struct skimmerGammaConversions {

  // configurables for CCDB access
  Configurable<std::string> ccdbPath{"ccdb-path", "GLO/GRP/GRP", "path to the ccdb object"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "path to the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<float> kfMassConstrain{"KFParticleMassConstrain", 0.f, "mass constrain for the KFParticle mother particle"};

  Configurable<int> mincrossedrows{"mincrossedrows", 10, "min. crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max. chi2/NclsTPC"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 5.0, "max. TPC n sigma for electron"};
  Configurable<float> dcamin{"dcamin", 0.1, "dcamin"};
  Configurable<float> dcamax{"dcamax", 1e+10, "dcamax"};

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

  struct recalculatedVertexParameters {
    float recalculatedConversionPoint[3];
    float KFParticleChi2DividedByNDF;
  };

  Produces<aod::V0Photons> v0photons;
  Produces<aod::V0Legs> v0legs;
  Produces<aod::McGammasTrue> fFuncTableMcGammasFromConfirmedV0s;
  Produces<aod::V0RecalculationAndKF> fFuncTableV0Recalculated;
  Produces<aod::V0DaughterMcParticles> fFuncTableMCTrackInformation;
  Produces<aod::MCParticleIndex> fIndexTableMCTrackIndex;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int runNumber = -1;
  o2::base::MatLayerCylSet* lut = nullptr;

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

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    // if run number matches then magnetic field should also match. Avoids unessecary acceses.
    if (runNumber == bc.runNumber()) {
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = nullptr;

    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << ccdbPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      } else {
        LOG(info) << "Magnetic field initialized from GRPMagField";
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
    }

    o2::base::Propagator::Instance()->setMatLUT(lut);
    runNumber = bc.runNumber();
  }

  template <typename TV0, typename TTRACK>
  void fillTrackTable(TV0 const& theV0, TTRACK const& theTrack)
  {
    v0legs(theTrack.collisionId(),
           theTrack.globalIndex(), theTrack.sign(), false,
           theTrack.pt(), theTrack.eta(), theTrack.phi(), theTrack.p(), theTrack.dcaXY(), theTrack.dcaZ(),
           theTrack.tpcNClsFindable(), theTrack.tpcNClsFindableMinusFound(), theTrack.tpcNClsFindableMinusCrossedRows(),
           theTrack.tpcChi2NCl(), theTrack.tpcInnerParam(), theTrack.tpcSignal(),
           theTrack.tpcNSigmaEl(), theTrack.tpcNSigmaPi(),
           theTrack.itsClusterMap(), theTrack.itsChi2NCl(), theTrack.detectorMap());
  }

  template <typename TV0>
  void fillV0RecalculatedTable(TV0 const& theV0, recalculatedVertexParameters recalculatedVertex)
  {
    fFuncTableV0Recalculated(
      recalculatedVertex.recalculatedConversionPoint[0],
      recalculatedVertex.recalculatedConversionPoint[1],
      recalculatedVertex.recalculatedConversionPoint[2],
      recalculatedVertex.KFParticleChi2DividedByNDF);
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

  template <typename TLeg>
  bool checkV0leg(TLeg const& leg)
  {
    if (abs(leg.eta()) > maxeta) {
      return false;
    }
    if (abs(leg.tpcNSigmaEl()) > maxTPCNsigmaEl) {
      return false;
    }
    if (leg.tpcChi2NCl() > maxchi2tpc) {
      return false;
    }
    if (leg.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }
    if (abs(leg.dcaXY()) < dcamin) {
      return false;
    }
    if (dcamax < abs(leg.dcaXY())) {
      return false;
    }
    return true;
  }

  // ============================ FUNCTION DEFINITIONS ====================================================

  Preslice<aod::V0Datas> perCollision = aod::v0data::collisionId;

  void processRec(aod::Collisions const& collisions,
                  aod::BCsWithTimestamps const& bcs,
                  aod::V0Datas const& V0s,
                  tracksAndTPCInfo const& theTracks)
  {
    for (auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("hCollisionZ_Rec"), collision.posZ());

      auto groupedV0s = V0s.sliceBy(perCollision, collision.globalIndex());
      for (auto& v0 : groupedV0s) {
        if (!checkAP(v0.alpha(), v0.qtarm())) { // store only photon conversions
          continue;
        }
        auto pos = v0.template posTrack_as<tracksAndTPCInfo>(); // positive daughter
        auto ele = v0.template negTrack_as<tracksAndTPCInfo>(); // negative daughter
        if (!checkV0leg(pos) || !checkV0leg(ele)) {
          continue;
        }

        bool flag_closer = true;
        for (auto& v0tmp : groupedV0s) {
          if (!checkAP(v0tmp.alpha(), v0tmp.qtarm())) { // store only photon conversions
            continue;
          }
          auto pos_tmp = v0tmp.template posTrack_as<tracksAndTPCInfo>(); // positive daughter
          auto ele_tmp = v0tmp.template negTrack_as<tracksAndTPCInfo>(); // negative daughter
          if (!checkV0leg(pos_tmp) || !checkV0leg(ele_tmp)) {
            continue;
          }

          if (v0.index() == v0tmp.index()) { // don't check onviously, exactly the same v0.
            // LOGF(info, "don't check the exactly the same 2 V0s");
            continue;
          }
          if ((ele.globalIndex() == ele_tmp.globalIndex() || pos.globalIndex() == pos_tmp.globalIndex()) && v0.dcaV0daughters() > v0tmp.dcaV0daughters()) {
            // LOGF(info, "!reject! | collision id = %d | g1 id = %d , g2 id = %d , posid1 = %d , eleid1 = %d , posid2 = %d , eleid2 = %d , pca1 = %f , pca2 = %f",
            //     collision.globalIndex(), v0.index(), v0tmp.index(), pos.globalIndex(), ele.globalIndex(), pos_tmp.globalIndex(), ele_tmp.globalIndex(), v0.dcaV0daughters(), v0tmp.dcaV0daughters());
            flag_closer = false;
            break;
          }
        } // end of v0tmp loop

        if (!flag_closer) {
          continue;
        }

        recalculatedVertexParameters recalculatedVertex;
        Vtx_recalculation(pos, ele, &recalculatedVertex);

        v0photons(collision.globalIndex(), v0legs.lastIndex() + 1, v0legs.lastIndex() + 2,
                  v0.x(), v0.y(), v0.z(),
                  v0.pxpos(), v0.pypos(), v0.pzpos(),
                  v0.pxneg(), v0.pyneg(), v0.pzneg(),
                  v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()), v0.dcaV0daughters()); // if v0legs is empty, lastIndex = -1.

        fillTrackTable(v0, pos);
        fillTrackTable(v0, ele);
        fillV0RecalculatedTable(v0, recalculatedVertex);
      } // end of v0 loop
    }   // end of collision loop
  }
  PROCESS_SWITCH(skimmerGammaConversions, processRec, "process reconstructed info only", true);

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  void processMc(soa::Join<aod::McCollisionLabels, aod::Collisions> const& collisions,
                 aod::McCollisions const&,
                 aod::BCsWithTimestamps const& bcs,
                 aod::V0Datas const& theV0s,
                 tracksAndTPCInfoMC const& theTracks,
                 aod::McParticles const& mcTracks)
  {
    for (auto& collision : collisions) {

      if (!collision.has_mcCollision()) {
        continue;
      }
      auto mcCollision = collision.mcCollision();
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      fRegistry.fill(HIST("hCollisionZ_all_MCTrue"), mcCollision.posZ());
      fRegistry.fill(HIST("hCollisionZ_MCTrue"), mcCollision.posZ());

      auto groupedMcTracks = mcTracks.sliceBy(perMcCollision, mcCollision.globalIndex());
      fRegistry.fill(HIST("hMcParticlesSize"), groupedMcTracks.size());
      fRegistry.fill(HIST("hCollisionZ_MCRec"), collision.posZ());

      auto lGroupedV0s = theV0s.sliceBy(perCollision, collision.globalIndex());
      for (auto& v0 : lGroupedV0s) {
        if (!checkAP(v0.alpha(), v0.qtarm())) { // store only photon conversions
          continue;
        }
        auto pos = v0.template posTrack_as<tracksAndTPCInfoMC>(); // positive daughter
        auto ele = v0.template negTrack_as<tracksAndTPCInfoMC>(); // negative daughter
        if (!checkV0leg(pos) || !checkV0leg(ele)) {
          continue;
        }

        bool flag_closer = true;
        for (auto& v0tmp : lGroupedV0s) {
          if (!checkAP(v0tmp.alpha(), v0tmp.qtarm())) { // store only photon conversions
            continue;
          }
          auto pos_tmp = v0tmp.template posTrack_as<tracksAndTPCInfoMC>(); // positive daughter
          auto ele_tmp = v0tmp.template negTrack_as<tracksAndTPCInfoMC>(); // negative daughter
          if (!checkV0leg(pos_tmp) || !checkV0leg(ele_tmp)) {
            continue;
          }

          if (v0.index() == v0tmp.index()) { // don't check onviously, exactly the same v0.
            // LOGF(info, "don't check the exactly the same 2 V0s");
            continue;
          }
          if ((ele.globalIndex() == ele_tmp.globalIndex() || pos.globalIndex() == pos_tmp.globalIndex()) && v0.dcaV0daughters() > v0tmp.dcaV0daughters()) {
            // LOGF(info, "!reject! | collision id = %d | g1 id = %d , g2 id = %d , posid1 = %d , eleid1 = %d , posid2 = %d , eleid2 = %d , pca1 = %f , pca2 = %f",
            //     collision.globalIndex(), v0.index(), v0tmp.index(), pos.globalIndex(), ele.globalIndex(), pos_tmp.globalIndex(), ele_tmp.globalIndex(), v0.dcaV0daughters(), v0tmp.dcaV0daughters());
            flag_closer = false;
            break;
          }
        } // end of v0tmp loop

        if (!flag_closer) {
          continue;
        }

        eV0Confirmation v0Status = isTrueV0(v0, pos, ele);

        fRegistry.get<TH1>(HIST("hV0Confirmation"))->Fill(v0Status);

        recalculatedVertexParameters recalculatedVertex;
        Vtx_recalculation(pos, ele, &recalculatedVertex);

        v0photons(collision.globalIndex(), v0legs.lastIndex() + 1, v0legs.lastIndex() + 2,
                  v0.x(), v0.y(), v0.z(),
                  v0.pxpos(), v0.pypos(), v0.pzpos(),
                  v0.pxneg(), v0.pyneg(), v0.pzneg(),
                  v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()), v0.dcaV0daughters()); // if v0legs is empty, lastIndex = -1.

        fillTrackTable(v0, pos);
        fillTrackTable(v0, ele);
        fillV0RecalculatedTable(v0, recalculatedVertex);
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
      fMotherSizesHisto->Fill(0.5 + static_cast<float>(lMothersIndeces.size()));
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
        // theV0.v0Id(),
        v0photons.lastIndex() + 1,
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
  void Vtx_recalculation(T lTrackPos, T lTrackNeg, recalculatedVertexParameters* recalculatedVertex)
  {
    o2::base::Propagator* prop = o2::base::Propagator::Instance(); // This singleton propagator requires some initialisation of the CCDB object.
    float bz = prop->getNominalBz();

    //*******************************************************

    // o2::track::TrackParametrizationWithError<TrackPrecision> = TrackParCov, I use the full version to have control over the data type
    o2::track::TrackParametrizationWithError<TrackPrecision> trackPosInformation = getTrackParCov(lTrackPos); // first get an object that stores Track information (positive)
    o2::track::TrackParametrizationWithError<TrackPrecision> trackNegInformation = getTrackParCov(lTrackNeg); // first get an object that stores Track information (negative)

    o2::track::TrackAuxPar helixPos(trackPosInformation, bz); // This object is a descendant of a CircleXY and stores cirlce information with respect to the magnetic field. This object uses functions and information of the o2::track::TrackParametrizationWithError<TrackPrecision> object (positive)
    o2::track::TrackAuxPar helixNeg(trackNegInformation, bz); // This object is a descendant of a CircleXY and stores cirlce information with respect to the magnetic field. This object uses functions and information of the o2::track::TrackParametrizationWithError<TrackPrecision> object (negative)

    recalculatedVertex->recalculatedConversionPoint[0] = (helixPos.xC * helixNeg.rC + helixNeg.xC * helixPos.rC) / (helixPos.rC + helixNeg.rC); // This calculates the coordinates of the conversion point as an weighted average of the two helix centers. xC and yC should be the global coordinates for the helix center as far as I understand. But you can double check the code of trackPosInformation.getCircleParamsLoc
    recalculatedVertex->recalculatedConversionPoint[1] = (helixPos.yC * helixNeg.rC + helixNeg.yC * helixPos.rC) / (helixPos.rC + helixNeg.rC); // If this calculation doesn't work check if the rotateZ function, because the "documentation" says I get global coordinates but maybe i don't.

    // I am unsure about the Z calculation but this is how it is done in AliPhysics as far as I understand
    o2::track::TrackParametrizationWithError<TrackPrecision> trackPosInformationCopy = o2::track::TrackParametrizationWithError<TrackPrecision>(trackPosInformation);
    o2::track::TrackParametrizationWithError<TrackPrecision> trackNegInformationCopy = o2::track::TrackParametrizationWithError<TrackPrecision>(trackNegInformation);

    // I think this calculation gets the closest point on the track to the conversion point
    // This alpha is a different alpha than the usual alpha and I think it is the angle between X axis and conversion point
    Double_t alphaPos = TMath::Pi() + TMath::ATan2(-(recalculatedVertex->recalculatedConversionPoint[1] - helixPos.yC), (recalculatedVertex->recalculatedConversionPoint[0] - helixPos.xC));
    Double_t alphaNeg = TMath::Pi() + TMath::ATan2(-(recalculatedVertex->recalculatedConversionPoint[1] - helixNeg.yC), (recalculatedVertex->recalculatedConversionPoint[0] - helixNeg.xC));

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
    // o2::base::PropagatorImpl<TrackPrecision>::MatCorrType::USEMatCorrLUT);
    prop->propagateToX(trackNegInformationCopy,
                       vertexNegRot.X(),
                       bz,
                       o2::base::PropagatorImpl<TrackPrecision>::MAX_SIN_PHI,
                       o2::base::PropagatorImpl<TrackPrecision>::MAX_STEP,
                       o2::base::PropagatorImpl<TrackPrecision>::MatCorrType::USEMatCorrNONE);
    // o2::base::PropagatorImpl<TrackPrecision>::MatCorrType::USEMatCorrLUT);

    // TODO: This is still off and needs to be checked...
    recalculatedVertex->recalculatedConversionPoint[2] = (trackPosInformationCopy.getZ() * helixNeg.rC + trackNegInformationCopy.getZ() * helixPos.rC) / (helixPos.rC + helixNeg.rC);
    KFPTrack kFTrackPos = createKFPTrackFromTrackParCov(trackPosInformationCopy, lTrackPos.sign(), lTrackPos.tpcNClsFound(), lTrackPos.tpcChi2NCl());
    int pdg_ePlus = -11; // e+
    KFParticle kFParticleEPlus(kFTrackPos, pdg_ePlus);

    KFPTrack kFTrackNeg = createKFPTrackFromTrackParCov(trackNegInformationCopy, lTrackNeg.sign(), lTrackNeg.tpcNClsFound(), lTrackNeg.tpcChi2NCl());
    int pdg_eMinus = 11; // e-
    KFParticle kFParticleEMinus(kFTrackNeg, pdg_eMinus);

    KFParticle gammaKF;
    gammaKF.SetConstructMethod(2);
    gammaKF.AddDaughter(kFParticleEPlus);
    gammaKF.AddDaughter(kFParticleEMinus);
    gammaKF.SetNonlinearMassConstraint(kfMassConstrain);

    if (gammaKF.GetNDF() == 0) {
      recalculatedVertex->KFParticleChi2DividedByNDF = -1.f;
    } else {
      recalculatedVertex->KFParticleChi2DividedByNDF = gammaKF.GetChi2() / gammaKF.GetNDF();
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerGammaConversions>(cfgc)};
}
