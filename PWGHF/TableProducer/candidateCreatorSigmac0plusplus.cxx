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

/// \file candidateCreatorSigmac0plusplus.cxx
/// \brief Σc0,++ → Λc+(→pK-π+) π-,+ candidate builder
/// \note Λc± candidates selected from the HFLcCandidateSelector.cxx
///
/// \author Mattia Faggin <mfaggin@cern.ch>, University and INFN PADOVA

#include <string>
#include <set>
#include <vector>

#include "CCDB/BasicCCDBManager.h" // for dca recalculation
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h" // for dca recalculation
#include "DataFormatsParameters/GRPObject.h"   // for dca recalculation
#include "DetectorsBase/GeometryManager.h"     // for dca recalculation
#include "DetectorsBase/Propagator.h"          // for dca recalculation
#include "DetectorsVertexing/PVertexer.h"      // for dca recalculation
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h" // for dca recalculation

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfCandidateCreatorSigmac0plusplus {

  /// Table with Σc0,++ info
  Produces<aod::HfCandScBase> rowScCandidates;

  /// Selection of candidates Λc+
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<double> yCandLcMax{"yCandLcMax", -1., "max. candLc. Lc rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_sigmac_to_p_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cutsMassLcMax{"cutsMassLcMax", {hf_cuts_sigmac_to_p_k_pi::cuts[0], hf_cuts_sigmac_to_p_k_pi::nBinsPt, hf_cuts_sigmac_to_p_k_pi::nCutVars, hf_cuts_sigmac_to_p_k_pi::labelsPt, hf_cuts_sigmac_to_p_k_pi::labelsCutVar}, "Lc candidate selection per pT bin"};

  /// Selections on candidate soft π-,+
  Configurable<bool> applyGlobalTrkWoDcaCutsSoftPi{"applyGlobalTrkWoDcaCutsSoftPi", false, "Switch on the application of the global-track w/o dca cuts for soft pion BEFORE ALL OTHER CUSTOM CUTS"};
  Configurable<float> softPiEtaMax{"softPiEtaMax", 0.9f, "Soft pion max value for pseudorapidity (abs vale)"};
  Configurable<float> softPiChi2Max{"softPiChi2Max", 36.f, "Soft pion max value for chi2 ITS"};
  Configurable<int> softPiItsHitMap{"softPiItsHitMap", 127, "Soft pion ITS hitmap"};
  Configurable<int> softPiItsHitsMin{"softPiItsHitsMin", 1, "Minimum number of ITS layers crossed by the soft pion among those in \"softPiItsHitMap\""};
  Configurable<float> softPiDcaXYMax{"softPiDcaXYMax", 0.065, "Soft pion max dcaXY (cm)"};
  Configurable<float> softPiDcaZMax{"softPiDcaZMax", 0.065, "Soft pion max dcaZ (cm)"};

  // CCDB
  Configurable<bool> isRun2Ccdb{"isRun2Ccdb", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  HfHelper hfHelper;
  /// Cut selection object for soft π-,+
  TrackSelection softPiCuts;

  // Needed for dcaXY, dcaZ recalculation of soft pions reassigned to a new collision
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber;

  using CandidatesLc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;

  /// Filter the candidate Λc+ used for the Σc0,++ creation
  Filter filterSelectCandidateLc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc);

  /// Prepare the slicing of candidate Λc+ and pions to be used for the Σc0,++ creation
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  // Preslice<CandidatesLc> hf3ProngPerCollision = aod::track_association::collisionId;
  Preslice<CandidatesLc> hf3ProngPerCollision = aod::hf_cand::collisionId;

  HistogramRegistry histos;

  /// @brief init function, to define the soft pion selections and histograms
  /// @param
  void init(InitContext&)
  {
    auto h = histos.add<TH1>("hCounter", "", kTH1D, {{7, 0.5, 7.5}});
    h->GetXaxis()->SetBinLabel(1, "collisions");
    h->GetXaxis()->SetBinLabel(2, "soft #pi (before cuts)");
    h->GetXaxis()->SetBinLabel(3, "soft #pi (after track cuts)");
    h->GetXaxis()->SetBinLabel(4, "Lc in track loop (before cuts)");
    h->GetXaxis()->SetBinLabel(5, "Lc in track loop (after cuts)");
    h->GetXaxis()->SetBinLabel(6, "candidate #Sigma_{c} with soft #pi != Lc daughter");
    h->GetXaxis()->SetBinLabel(7, "#Sigma_{c}");

    /// process function switches
    std::array<int, 2> arrProcess = {doprocessDataTrackToCollAssoc, doprocessDataNoTrackToCollAssoc};
    int processes = std::accumulate(arrProcess.begin(), arrProcess.end(), 0);
    if (processes != 1) {
      LOG(fatal) << "Check the enabled process functions. doprocessDataTrackToCollAssoc=" << doprocessDataTrackToCollAssoc << ", doprocessDataNoTrackToCollAssoc=" << doprocessDataNoTrackToCollAssoc;
    }

    ////////////////////////////////////////
    /// set the selections for soft pion ///
    ////////////////////////////////////////

    /// apply the global-track w/o dca cuts for soft pion BEFORE ALL OTHER CUSTOM CUTS
    if (applyGlobalTrkWoDcaCutsSoftPi) {

      LOG(info) << ">>> applyGlobalTrkWoDcaCutsSoftPi==true  ==>  global-track w/o dca cuts for soft pionapplied BEFORE ALL OTHER CUSTOM CUTS <<<";

      /// same configuration as in track selection (itsMatching==1)
      softPiCuts = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, 0);

      /// remove dca cuts (applied manually after the possible track-to-collision reassociation)
      softPiCuts.SetMaxDcaXY(99999);
      softPiCuts.SetMaxDcaZ(99999);
    }

    // kinematics
    // softPiCuts.SetPtRange(0.001, 1000.); // pt
    softPiCuts.SetEtaRange(-softPiEtaMax, softPiEtaMax); // eta
    // softPiCuts.SetMaxDcaXY(softPiDcaXYMax);              // dcaXY
    // softPiCuts.SetMaxDcaZ(softPiDcaZMax);                // dcaZ

    // ITS chi2
    softPiCuts.SetMaxChi2PerClusterITS(softPiChi2Max);
    //  ITS hitmap
    std::set<uint8_t> setSoftPiItsHitMap; // = {};
    for (int idItsLayer = 0; idItsLayer < 7; idItsLayer++) {
      if (TESTBIT(softPiItsHitMap, idItsLayer)) {
        setSoftPiItsHitMap.insert(static_cast<uint8_t>(idItsLayer));
      }
    }
    LOG(info) << "### ITS hitmap for soft pion";
    LOG(info) << "    >>> setSoftPiItsHitMap.size(): " << setSoftPiItsHitMap.size();
    LOG(info) << "    >>> Custom ITS hitmap dfchecked: ";
    for (std::set<uint8_t>::iterator it = setSoftPiItsHitMap.begin(); it != setSoftPiItsHitMap.end(); it++) {
      LOG(info) << "        Layer " << static_cast<int>(*it) << " ";
    }
    LOG(info) << "############";
    softPiCuts.SetRequireITSRefit();
    softPiCuts.SetRequireHitsInITSLayers(softPiItsHitsMin, setSoftPiItsHitMap);

    /// CCDB for dcaXY, dcaZ recalculation of soft pions reassigned to another collision
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }

  /// @param trackSoftPi is the track (with dcaXY, dcaZ information)of a candidate soft-pion in the collision
  /// @param candidates are 3-prong candidates satisfying the analysis selections for Λc+ → pK-π+ (and charge conj.)
  /// @param tracks are the tracks (with dcaXY, dcaZ information) → soft-pion candidate tracks
  template <typename TRK, typename CAND>
  void makeSoftPiLcPair(TRK const& trackSoftPi, CAND const& candidatesThisColl, aod::TracksWDcaExtra const&)
  {

    /// loop over Λc+ → pK-π+ (and charge conj.) candidates
    for (const auto& candLc : candidatesThisColl) {
      histos.fill(HIST("hCounter"), 4);

      /// keep only the candidates flagged as possible Λc+ (and charge conj.) decaying into a charged pion, kaon and proton
      /// if not selected, skip it and go to the next one
      if (!(candLc.hfflag() & 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        continue;
      }
      /// keep only the candidates Λc+ (and charge conj.) within the desired rapidity
      /// if not selected, skip it and go to the next one
      if (yCandLcMax >= 0. && std::abs(hfHelper.yLc(candLc)) > yCandLcMax) {
        continue;
      }

      /// selection on the Λc+ inv. mass window we want to consider for Σc0,++ candidate creation
      auto statusSpreadMinvPKPiFromPDG = 0;
      auto statusSpreadMinvPiKPFromPDG = 0;

      // define mass window for Lc
      float mPKPiCandLcMax = -1;
      float mPiKPCandLcMax = -1;
      int pTBin = findBin(binsPt, candLc.pt());
      if (pTBin != -1) {
        mPKPiCandLcMax = cutsMassLcMax->get(pTBin, "max pKpi mass Lc");
        mPiKPCandLcMax = cutsMassLcMax->get(pTBin, "max piKp mass Lc");
      }

      if (candLc.isSelLcToPKPi() >= 1 && std::abs(hfHelper.invMassLcToPKPi(candLc) - MassLambdaCPlus) <= mPKPiCandLcMax) {
        statusSpreadMinvPKPiFromPDG = 1;
      }
      if (candLc.isSelLcToPiKP() >= 1 && std::abs(hfHelper.invMassLcToPiKP(candLc) - MassLambdaCPlus) <= mPiKPCandLcMax) {
        statusSpreadMinvPiKPFromPDG = 1;
      }
      if (statusSpreadMinvPKPiFromPDG == 0 && statusSpreadMinvPiKPFromPDG == 0) {
        /// none of the two possibilities are satisfied, therefore this candidate Lc can be skipped
        continue;
      }
      histos.fill(HIST("hCounter"), 5);

      // auto trackSoftPi = tracks.rawIteratorAt(trackId.trackId());

      //////////////////////////////////////////////////////////////////////////////////////
      ///                       Σc0,++ candidate creation                                ///
      ///                                                                                ///
      /// For each candidate Λc, let's build a Σc0,++ with the candidate soft-pion track ///
      //////////////////////////////////////////////////////////////////////////////////////

      /// Exclude the current candidate soft pion if it corresponds already to a candidate Lc prong
      int indexProng0 = candLc.template prong0_as<aod::TracksWDcaExtra>().globalIndex();
      int indexProng1 = candLc.template prong1_as<aod::TracksWDcaExtra>().globalIndex();
      int indexProng2 = candLc.template prong2_as<aod::TracksWDcaExtra>().globalIndex();
      int indexSoftPi = trackSoftPi.globalIndex();
      if (indexSoftPi == indexProng0 || indexSoftPi == indexProng1 || indexSoftPi == indexProng2) {
        continue;
      }
      histos.fill(HIST("hCounter"), 6);

      /// determine the Σc candidate charge
      int chargeLc = candLc.template prong0_as<aod::TracksWDcaExtra>().sign() + candLc.template prong1_as<aod::TracksWDcaExtra>().sign() + candLc.template prong2_as<aod::TracksWDcaExtra>().sign();
      int chargeSoftPi = trackSoftPi.sign();
      int8_t chargeSigmac = chargeLc + chargeSoftPi;
      if (std::abs(chargeSigmac) != 0 && std::abs(chargeSigmac) != 2) {
        /// this shall never happen
        LOG(fatal) << ">>> Sc candidate with charge +1 built, not possible! Charge Lc: " << chargeLc << ", charge soft pion: " << chargeSoftPi;
      }
      histos.fill(HIST("hCounter"), 7);

      /// fill the Σc0,++ candidate table
      rowScCandidates(/* general columns */
                      candLc.collisionId(),
                      /* 2-prong specific columns */
                      candLc.px(), candLc.py(), candLc.pz(),
                      trackSoftPi.px(), trackSoftPi.py(), trackSoftPi.pz(),
                      candLc.globalIndex(), trackSoftPi.globalIndex(),
                      candLc.hfflag(),
                      /* Σc0,++ specific columns */
                      chargeSigmac,
                      statusSpreadMinvPKPiFromPDG, statusSpreadMinvPiKPFromPDG);
    } /// end loop over Λc+ → pK-π+ (and charge conj.) candidates
  }   /// end makeSoftPiLcPair

  /// @brief function to loop over candidate soft pions and, for each of them, over candidate Λc+ for Σc0,++ → Λc+(→pK-π+) π- candidate reconstruction
  /// @param collision is a o2::aod::Collisions
  /// @param trackSoftPi is the track (with dcaXY, dcaZ information)of a candidate soft-pion in the collision
  /// @param tracks are the tracks (with dcaXY, dcaZ information) → soft-pion candidate tracks
  /// @param candidates are 3-prong candidates satisfying the analysis selections for Λc+ → pK-π+ (and charge conj.)
  template <bool withTimeAssoc, typename TRK>
  void createSigmaC(aod::Collisions::iterator const& collision,
                    TRK const& trackSoftPi,
                    aod::TracksWDcaExtra const& tracks,
                    CandidatesLc const& candidates,
                    aod::BCsWithTimestamps const&)
  {

    auto thisCollId = collision.globalIndex();

    histos.fill(HIST("hCounter"), 2);
    /// keep only soft-pion candidate tracks
    /// if not selected, skip it and go to the next one
    if (!softPiCuts.IsSelected(trackSoftPi)) {
      return;
    }
    if constexpr (withTimeAssoc) {
      /// dcaXY, dcaZ selections
      /// To be done separately from the others, because for reassigned tracks the dca must be recalculated
      /// TODO: to be properly adapted in case of PV refit usage
      if (trackSoftPi.collisionId() == thisCollId) {
        /// this is a track originally assigned to the current collision
        /// therefore, the dcaXY, dcaZ are those already calculated in the track-propagation workflow
        if (std::abs(trackSoftPi.dcaXY()) > softPiDcaXYMax || std::abs(trackSoftPi.dcaZ()) > softPiDcaZMax) {
          return;
        }
      } else {
        /// this is a reassigned track
        /// therefore we need to calculate the dcaXY, dcaZ with respect to this new primary vertex
        auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
        initCCDB(bc, runNumber, ccdb, isRun2Ccdb ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2Ccdb);
        auto trackParSoftPi = getTrackPar(trackSoftPi);
        o2::gpu::gpustd::array<float, 2> dcaInfo{-999., -999.};
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParSoftPi, 2.f, noMatCorr, &dcaInfo);
        if (std::abs(dcaInfo[0]) > softPiDcaXYMax || std::abs(dcaInfo[1]) > softPiDcaZMax) {
          return;
        }
      }
    } else {
      /// this is a track originally assigned to the current collision
      /// therefore, the dcaXY, dcaZ are those already calculated in the track-propagation workflow
      /// No need to consider the time-reassociated tracks, since withTimeAssoc == false here
      if (std::abs(trackSoftPi.dcaXY()) > softPiDcaXYMax || std::abs(trackSoftPi.dcaZ()) > softPiDcaZMax) {
        return;
      }
    }
    histos.fill(HIST("hCounter"), 3);

    /// loop over Λc+ → pK-π+ (and charge conj.) candidates
    if constexpr (withTimeAssoc) {
      /// need to group candidates manually
      auto candidatesThisColl = candidates.sliceBy(hf3ProngPerCollision, thisCollId);
      makeSoftPiLcPair(trackSoftPi, candidatesThisColl, tracks);
    } else {
      /// tracks and candidates already grouped by collision at the level of process function
      makeSoftPiLcPair(trackSoftPi, candidates, tracks);
    }

  } /// end createSigmaC

  /// @brief process function for Σc0,++ → Λc+(→pK-π+) π- candidate reconstruction considering also time-reassigned tracks for soft pions
  /// @param collisions are o2::aod::Collisions
  /// @param trackIndices are the indices of tracks reassigned to collisions, as obtained from the track-to-collision-associator
  /// @param tracks are the tracks (with dcaXY, dcaZ information) → soft-pion candidate tracks
  /// @param candidates are 3-prong candidates satisfying the analysis selections for Λc+ → pK-π+ (and charge conj.)
  void processDataTrackToCollAssoc(aod::Collisions const& collisions,
                                   aod::TrackAssoc const& trackIndices,
                                   aod::TracksWDcaExtra const& tracks,
                                   CandidatesLc const& candidates,
                                   aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    for (const auto& collision : collisions) {

      histos.fill(HIST("hCounter"), 1);
      // LOG(info) << "[processDataTrackToCollAssoc] Collision with globalIndex " << collision.globalIndex();

      // slice by hand the assoc. track with time per collision
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      // LOG(info) << "[processDataTrackToCollAssoc]     - number of tracks: " << trackIdsThisCollision.size();

      /// loop over tracks for soft pion
      for (const auto& trackId : trackIdsThisCollision) {
        /// slice soft pion tracks associated to the current collision
        auto trackSoftPi = trackId.track_as<aod::TracksWDcaExtra>();

        /// create SigmaC candidate with the current soft pion and Lc candidates
        createSigmaC<true>(collision, trackSoftPi, tracks, candidates, bcWithTimeStamps);

      } /// end loop over tracks for soft pion

    } // end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorSigmac0plusplus, processDataTrackToCollAssoc, "Process data using also time-reassociated tracks", true);

  /// @brief process function for Σc0,++ → Λc+(→pK-π+) π- candidate reconstruction without considering time-reassigned tracks for soft pions
  /// @param collision is a o2::aod::Collision
  /// @param tracks are the tracks (with dcaXY, dcaZ information) → soft-pion candidate tracks
  /// @param candidates are 3-prong candidates satisfying the analysis selections for Λc+ → pK-π+ (and charge conj.)
  void processDataNoTrackToCollAssoc(aod::Collision const& collision,
                                     aod::TracksWDcaExtra const& tracks,
                                     CandidatesLc const& candidates,
                                     aod::BCsWithTimestamps const& bcWithTimeStamps)
  {

    histos.fill(HIST("hCounter"), 1);
    // LOG(info) << "[processDataNoTrackToCollAssoc] Collision with globalIndex " << collision.globalIndex();
    // LOG(info) << "[processDataNoTrackToCollAssoc]     - number of tracks: " << tracks.size();
    // LOG(info) << "[processDataNoTrackToCollAssoc]     - number of Lc candidates: " << candidates.size();

    /// loop over tracks for soft pion
    /// In this case, they are already grouped by collision, using the track::CollisionId
    for (const auto& trackSoftPi : tracks) {

      /// create SigmaC candidate with the current soft pion and Lc candidates
      createSigmaC<false>(collision, trackSoftPi, tracks, candidates, bcWithTimeStamps);

    } /// end loop over tracks for soft pion
  }
  PROCESS_SWITCH(HfCandidateCreatorSigmac0plusplus, processDataNoTrackToCollAssoc, "Process data using time reassociated tracks", false);
};

/// Extends the base table with expression columns.

struct HfCandidateSigmac0plusplusMc {
  Spawns<aod::HfCandScExt> candidatesSigmac;
  Produces<aod::HfCandScMcRec> rowMCMatchScRec;
  Produces<aod::HfCandScMcGen> rowMCMatchScGen;

  using LambdacMc = soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>;
  // using LambdacMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  float zPvPosMax{1000.f};

  /// @brief init function
  void init(InitContext& initContext)
  {
    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name.compare("hf-candidate-creator-3prong") == 0) { // here we assume that the hf-candidate-creator-3prong is in the workflow
        for (const auto& option : device.options) {
          if (option.name.compare("hfEvSel.zPvPosMax") == 0) {
            zPvPosMax = option.defaultValue.get<float>();
            break;
          }
        }
        break;
      }
    }
  }

  /// @brief dummy process function, to be run on data
  /// @param
  void process(aod::Tracks const&) {}

  /// @brief process function for MC matching of Σc0,++ → Λc+(→pK-π+) π- reconstructed candidates and counting of generated ones
  /// @param candidatesSigmac reconstructed Σc0,++ candidates
  /// @param mcParticles table of generated particles
  void processMc(aod::McParticles const& mcParticles,
                 aod::TracksWMc const& tracks,
                 LambdacMc const& candsLc /*, const LambdacMcGen&*/,
                 aod::McCollisions const&)
  {

    // Match reconstructed candidates.
    candidatesSigmac->bindExternalIndices(&tracks);
    candidatesSigmac->bindExternalIndices(&candsLc);

    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t chargeSigmac = 10;

    /// Match reconstructed Σc0,++ candidates
    for (const auto& candSigmac : *candidatesSigmac) {
      indexRec = -1;
      sign = 0;
      flag = 0;
      origin = 0;
      std::vector<int> idxBhadMothers{};

      /// skip immediately the candidate Σc0,++ w/o a Λc+ matched to MC
      auto candLc = candSigmac.prongLc_as<LambdacMc>();
      if (!(std::abs(candLc.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) { /// (*)
        rowMCMatchScRec(flag, origin, -1.f, 0);
        continue;
      }

      /// matching to MC
      auto arrayDaughters = std::array{candLc.prong0_as<aod::TracksWMc>(),
                                       candLc.prong1_as<aod::TracksWMc>(),
                                       candLc.prong2_as<aod::TracksWMc>(),
                                       candSigmac.prong1_as<aod::TracksWMc>()};
      chargeSigmac = candSigmac.charge();
      if (chargeSigmac == 0) {
        /// candidate Σc0
        /// 3 levels:
        ///   1. Σc0 → Λc+ π-,+
        ///   2. Λc+ → pK-π+ direct (i) or Λc+ → resonant channel Λc± → p± K*, Λc± → Δ(1232)±± K∓ or Λc± → Λ(1520) π±  (ii)
        ///   3. in case of (ii): resonant channel to pK-π+
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kSigmaC0, std::array{+kProton, -kKPlus, +kPiPlus, -kPiPlus}, true, &sign, 3);
        if (indexRec > -1) { /// due to (*) no need to check anything for LambdaC
          flag = sign * (1 << aod::hf_cand_sigmac::DecayType::Sc0ToPKPiPi);
        }
      } else if (std::abs(chargeSigmac) == 2) {
        /// candidate Σc++
        /// 3 levels:
        ///   1. Σc0 → Λc+ π-,+
        ///   2. Λc+ → pK-π+ direct (i) or Λc+ → resonant channel Λc± → p± K*, Λc± → Δ(1232)±± K∓ or Λc± → Λ(1520) π±  (ii)
        ///   3. in case of (ii): resonant channel to pK-π+
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kSigmaCPlusPlus, std::array{+kProton, -kKPlus, +kPiPlus, +kPiPlus}, true, &sign, 3);
        if (indexRec > -1) { /// due to (*) no need to check anything for LambdaC
          flag = sign * (1 << aod::hf_cand_sigmac::DecayType::ScplusplusToPKPiPi);
        }
      }

      /// check the origin (prompt vs. non-prompt)
      if (flag != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
      }
      /// fill the table with results of reconstruction level MC matching
      if (origin == RecoDecay::OriginType::NonPrompt) {
        auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
        rowMCMatchScRec(flag, origin, bHadMother.pt(), bHadMother.pdgCode());
      } else {
        rowMCMatchScRec(flag, origin, -1.f, 0);
      }
    } /// end loop over reconstructed Σc0,++ candidates

    /// Match generated Σc0,++ candidates
    for (const auto& particle : mcParticles) {
      flag = 0;
      origin = 0;
      std::vector<int> idxBhadMothers{};

      auto mcCollision = particle.mcCollision();
      float zPv = mcCollision.posZ();
      if (zPv < -zPvPosMax || zPv > zPvPosMax) { // to avoid counting particles in collisions with Zvtx larger than the maximum, we do not match them
        rowMCMatchScGen(flag, origin, -1);
        continue;
      }

      /// 3 levels:
      ///   1. Σc0 → Λc+ π-,+
      ///   2. Λc+ → pK-π+ direct (i) or Λc+ → resonant channel Λc± → p± K*, Λc± → Δ(1232)±± K∓ or Λc± → Λ(1520) π±  (ii)
      ///   3. in case of (ii): resonant channel to pK-π+
      /// → here we check level 1. first, and then levels 2. and 3. are inherited by the Λc+ → pK-π+ MC matching in candidateCreator3Prong.cxx
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kSigmaC0, std::array{static_cast<int>(Pdg::kLambdaCPlus), static_cast<int>(kPiMinus)}, true, &sign, 1)) {
        // generated Σc0
        // for (const auto& daughter : particle.daughters_as<LambdacMcGen>()) {
        for (const auto& daughter : particle.daughters_as<aod::McParticles>()) {
          // look for Λc+ daughter decaying in pK-π+
          if (std::abs(daughter.pdgCode()) != Pdg::kLambdaCPlus)
            continue;
          // if (std::abs(daughter.flagMcMatchGen()) == (1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
          if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2)) {
            /// Λc+ daughter decaying in pK-π+ found!
            flag = sign * (1 << aod::hf_cand_sigmac::DecayType::Sc0ToPKPiPi);
            break;
          }
        }
      } else if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kSigmaCPlusPlus, std::array{static_cast<int>(Pdg::kLambdaCPlus), static_cast<int>(kPiPlus)}, true, &sign, 1)) {
        // generated Σc++
        // for (const auto& daughter : particle.daughters_as<LambdacMcGen>()) {
        for (const auto& daughter : particle.daughters_as<aod::McParticles>()) {
          // look for Λc+ daughter decaying in pK-π+
          if (std::abs(daughter.pdgCode()) != Pdg::kLambdaCPlus)
            continue;
          // if (std::abs(daughter.flagMcMatchGen()) == (1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
          if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2)) {
            /// Λc+ daughter decaying in pK-π+ found!
            flag = sign * (1 << aod::hf_cand_sigmac::DecayType::ScplusplusToPKPiPi);
            break;
          }
        }
      }

      /// check the origin (prompt vs. non-prompt)
      if (flag != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
      }
      /// fill the table with results of generation level MC matching
      if (origin == RecoDecay::OriginType::NonPrompt) {
        rowMCMatchScGen(flag, origin, idxBhadMothers[0]);
      } else {
        rowMCMatchScGen(flag, origin, -1);
      }
    } /// end loop over mcParticles
  }   /// end processMc
  PROCESS_SWITCH(HfCandidateSigmac0plusplusMc, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorSigmac0plusplus>(cfgc),
    adaptAnalysisTask<HfCandidateSigmac0plusplusMc>(cfgc)};
}
