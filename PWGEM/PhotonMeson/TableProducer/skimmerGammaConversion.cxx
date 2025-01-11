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

#include <memory>
#include <map>
#include <string>
#include <vector>

// todo: remove reduantant information in GammaConversionsInfoTrue
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/gammaConvDefinitions.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

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
#include "CommonConstants/PhysicsConstants.h"

#include <TMath.h>
#include <TVector2.h>
#include "Math/Vector4D.h"

#include "Tools/KFparticle/KFUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

// using tracksAndTPCInfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::TracksCov>;
using tracksAndTPCInfo = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::TracksCovIU>;
using tracksAndTPCInfoMC = soa::Join<tracksAndTPCInfo, aod::McTrackLabels>;

struct skimmerGammaConversion {

  // configurables for CCDB access
  Configurable<std::string> ccdbPath{"ccdb-path", "GLO/GRP/GRP", "path to the ccdb object"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "path to the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<float> kfMassConstrain{"KFParticleMassConstrain", -1.f, "mass constrain for the KFParticle mother particle"};

  Configurable<int> mincrossedrows{"mincrossedrows", 10, "min. crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max. chi2/NclsTPC"};
  Configurable<float> minpt{"minpt", 0.01, "min pt for track"};
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

  Produces<aod::V0PhotonsKF> v0photonskf;
  Produces<aod::V0Legs> v0legs;
  Produces<aod::McGammasTrue> fFuncTableMcGammasFromConfirmedV0s;
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

    /// Set magnetic field for KF vertexing
    float magneticField = o2::base::Propagator::Instance()->getNominalBz();
    KFParticle::SetField(magneticField);
  }

  template <typename TTRACK, typename TKFP>
  void fillTrackTable(TTRACK const& theTrack, TKFP const& kfp)
  {
    v0legs(theTrack.collisionId(),
           theTrack.globalIndex(), theTrack.sign(),
           kfp.GetPx(), kfp.GetPy(), kfp.GetPz(), theTrack.dcaXY(), theTrack.dcaZ(),
           theTrack.tpcNClsFindable(), theTrack.tpcNClsFindableMinusFound(), theTrack.tpcNClsFindableMinusCrossedRows(), theTrack.tpcNClsShared(),
           theTrack.tpcChi2NCl(), theTrack.tpcInnerParam(), theTrack.tpcSignal(),
           theTrack.tpcNSigmaEl(), theTrack.tpcNSigmaPi(),
           theTrack.itsClusterSizes(), theTrack.itsChi2NCl(), theTrack.detectorMap(),
           theTrack.x(), theTrack.y(), theTrack.z(), theTrack.tgl());
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

  template <typename TTrack>
  bool checkV0leg(TTrack const& track)
  {
    if (track.pt() < minpt || abs(track.eta()) > maxeta) {
      return false;
    }
    if (abs(track.dcaXY()) < dcamin || dcamax < abs(track.dcaXY())) {
      return false;
    }
    if (!track.hasITS() && !track.hasTPC()) {
      return false;
    }

    if (track.hasTPC()) {
      if (track.tpcNClsCrossedRows() < mincrossedrows || track.tpcChi2NCl() > maxchi2tpc) {
        return false;
      }
      if (abs(track.tpcNSigmaEl()) > maxTPCNsigmaEl) {
        return false;
      }
    }
    return true;
  }

  template <typename TTrack, typename TCollision, typename TV0>
  void fillV0KF(TCollision const& collision, TV0 const& v0)
  {
    auto pos = v0.template posTrack_as<TTrack>(); // positive daughter
    auto ele = v0.template negTrack_as<TTrack>(); // negative daughter

    float xyz[3] = {0.f, 0.f, 0.f};
    Vtx_recalculation(o2::base::Propagator::Instance(), pos, ele, xyz);

    KFPTrack kfp_track_pos = createKFPTrackFromTrack(pos);
    KFPTrack kfp_track_ele = createKFPTrackFromTrack(ele);
    KFParticle kfp_pos(kfp_track_pos, -11);
    KFParticle kfp_ele(kfp_track_ele, 11);
    const KFParticle* GammaDaughters[2] = {&kfp_pos, &kfp_ele};

    KFParticle gammaKF;
    gammaKF.SetConstructMethod(2);
    gammaKF.Construct(GammaDaughters, 2);
    if (kfMassConstrain > -0.1) {
      gammaKF.SetNonlinearMassConstraint(kfMassConstrain);
    }
    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle KFPV(kfpVertex);

    //  LOGF(info, "recalculated vtx : x = %f , y = %f , z = %f", xyz[0], xyz[1], xyz[2]);
    //  LOGF(info, "primary vtx : x = %f , y = %f , z = %f", collision.posX(), collision.posY(), collision.posZ());

    // Transport the gamma to the recalculated decay vertex
    KFParticle gammaKF_DecayVtx = gammaKF; // with respect to (0,0,0)
    gammaKF_DecayVtx.TransportToPoint(xyz);

    //// Apply a topological constraint of the gamma to the PV. Parameters will be given at the primary vertex.
    // KFParticle gammaKF_PV = gammaKF_DecayVtx;
    // gammaKF_PV.SetProductionVertex(KFPV);

    float chi2kf = -1.f;
    if (gammaKF_DecayVtx.GetNDF() > 0) {
      chi2kf = gammaKF_DecayVtx.GetChi2() / gammaKF_DecayVtx.GetNDF();
    }

    KFParticle kfp_pos_DecayVtx = kfp_pos;
    KFParticle kfp_ele_DecayVtx = kfp_ele;
    kfp_pos_DecayVtx.TransportToPoint(xyz);
    kfp_ele_DecayVtx.TransportToPoint(xyz);

    // KFParticle kfp_pos_PV = kfp_pos_DecayVtx;
    // KFParticle kfp_ele_PV = kfp_ele_DecayVtx;
    // kfp_pos_PV.SetProductionVertex(KFPV);
    // kfp_ele_PV.SetProductionVertex(KFPV);

    // LOGF(info, "ele px = %f (original) , %f (KF at init) , %f (KF at PV) , %f (KF at SV)", ele.px(), kfp_ele.GetPx(), kfp_ele_PV.GetPx(), kfp_ele_DecayVtx.GetPx());
    // LOGF(info, "pos px = %f (original) , %f (KF at init) , %f (KF at PV) , %f (KF at SV)", pos.px(), kfp_pos.GetPx(), kfp_pos_PV.GetPx(), kfp_pos_DecayVtx.GetPx());

    // ROOT::Math::PxPyPzMVector vpos_pv(kfp_pos_PV.GetPx(), kfp_pos_PV.GetPy(), kfp_pos_PV.GetPz(), o2::constants::physics::MassElectron);
    // ROOT::Math::PxPyPzMVector vele_pv(kfp_ele_PV.GetPx(), kfp_ele_PV.GetPy(), kfp_ele_PV.GetPz(), o2::constants::physics::MassElectron);
    // ROOT::Math::PxPyPzMVector v0_pv = vpos_pv + vele_pv;

    ROOT::Math::PxPyPzMVector vpos_sv(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
    ROOT::Math::PxPyPzMVector vele_sv(kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
    ROOT::Math::PxPyPzMVector v0_sv = vpos_sv + vele_sv;
    // LOGF(info, "mee = %f (KF at PV) , %f (KF at SV)", v0_pv.M(), v0_sv.M());

    //// calculate psipair, phiv at the decay vertex
    // float phiv = getPhivPair(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz(), kfp_pos_DecayVtx.GetQ(), kfp_ele_DecayVtx.GetQ(), o2::base::Propagator::Instance()->getNominalBz());
    // float psipair = getPsiPair(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz());
    //// LOGF(info, "bz = %f , phiv = %f , psipair = %f", bz, phiv, psipair);

    float cospa_kf = cpaFromKF(gammaKF_DecayVtx, KFPV);
    float pca_kf = kfp_pos_DecayVtx.GetDistanceFromParticle(kfp_ele_DecayVtx);
    // LOGF(info, "pca = %f (DCAFitter) , %f (KF at SV)", v0.dcaV0daughters(), pca_kf);
    float alpha = v0_alpha(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz());
    float qt = v0_qt(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz());

    float v0mom = RecoDecay::sqrtSumOfSquares(gammaKF_DecayVtx.GetPx(), gammaKF_DecayVtx.GetPy(), gammaKF_DecayVtx.GetPz());
    float length = RecoDecay::sqrtSumOfSquares(gammaKF_DecayVtx.GetX() - collision.posX(), gammaKF_DecayVtx.GetY() - collision.posY(), gammaKF_DecayVtx.GetZ() - collision.posZ());
    float dca_x_v0_to_pv = (gammaKF_DecayVtx.GetX() - gammaKF_DecayVtx.GetPx() * cospa_kf * length / v0mom) - collision.posX();
    float dca_y_v0_to_pv = (gammaKF_DecayVtx.GetY() - gammaKF_DecayVtx.GetPy() * cospa_kf * length / v0mom) - collision.posY();
    float dca_z_v0_to_pv = (gammaKF_DecayVtx.GetZ() - gammaKF_DecayVtx.GetPz() * cospa_kf * length / v0mom) - collision.posZ();
    float sign_tmp = dca_y_v0_to_pv > 0 ? +1 : -1;
    float dca_xy_v0_to_pv = RecoDecay::sqrtSumOfSquares(dca_x_v0_to_pv, dca_y_v0_to_pv) * sign_tmp;

    v0photonskf(collision.globalIndex(), v0.globalIndex(), v0legs.lastIndex() + 1, v0legs.lastIndex() + 2,
                gammaKF_DecayVtx.GetX(), gammaKF_DecayVtx.GetY(), gammaKF_DecayVtx.GetZ(),
                gammaKF_DecayVtx.GetPx(), gammaKF_DecayVtx.GetPy(), gammaKF_DecayVtx.GetPz(),
                v0_sv.M(), dca_xy_v0_to_pv, dca_z_v0_to_pv,
                cospa_kf, pca_kf, alpha, qt, chi2kf);

    fillTrackTable(pos, kfp_pos_DecayVtx);
    fillTrackTable(ele, kfp_ele_DecayVtx);
  }

  // ============================ FUNCTION DEFINITIONS ====================================================

  PresliceUnsorted<aod::V0Datas> perCollision = aod::v0data::collisionId;

  void processRec(aod::Collisions const& collisions,
                  aod::BCsWithTimestamps const&,
                  aod::V0Datas const& V0s,
                  tracksAndTPCInfo const&)
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

        fillV0KF<tracksAndTPCInfo>(collision, v0);

      } // end of v0 loop
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerGammaConversion, processRec, "process reconstructed info only", true);

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  void processMc(soa::Join<aod::McCollisionLabels, aod::Collisions> const& collisions,
                 aod::McCollisions const&,
                 aod::BCsWithTimestamps const&,
                 aod::V0Datas const& theV0s,
                 tracksAndTPCInfoMC const&,
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

        if (!ele.has_mcParticle() || !pos.has_mcParticle()) {
          continue; // If no MC particle is found, skip the v0
        }

        eV0Confirmation v0Status = isTrueV0(v0, pos, ele);
        fRegistry.get<TH1>(HIST("hV0Confirmation"))->Fill(v0Status);

        fillV0KF<tracksAndTPCInfoMC>(collision, v0);
      }
    }
  }
  PROCESS_SWITCH(skimmerGammaConversion, processMc, "process reconstructed and mc info ", false);

  template <typename TV0, typename TTRACK>
  eV0Confirmation isTrueV0(TV0 const& /*theV0*/,
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
        v0photonskf.lastIndex() + 1,
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerGammaConversion>(cfgc, TaskName{"skimmer-gamma-conversion"})};
}
