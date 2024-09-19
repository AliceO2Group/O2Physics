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

/// \file taskFlowCharmHadrons.cxx
/// \brief Analysis task for charm hadron flow
///
/// \author S. Politanò, INFN Torino, Italy
/// \author Wu Chuntai, CUG, China

#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Qvectors.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;

enum DecayChannel { DplusToPiKPi = 0,
                    DsToKKPi,
                    DsToPiKK,
                    D0ToPiK,
                    D0ToKPi,
                    LcToPKPi,
                    LcToPiKP };

enum QvecEstimator { FV0A = 0,
                     FT0M,
                     FT0A,
                     FT0C,
                     TPCPos,
                     TPCNeg,
                     TPCTot };

struct HfTaskFlowCharmHadrons {
  Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<int> qvecDetector{"qvecDetector", 3, "Detector for Q vector estimation (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3, TPC Pos: 4, TPC Neg: 5, TPC Tot: 6)"};
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4)"};
  Configurable<int> selectionFlag{"selectionFlag", 1, "Selection Flag for hadron (e.g. 1 for skimming, 3 for topo. and kine., 7 for PID)"};
  Configurable<bool> storeMl{"storeMl", false, "Flag to store ML scores"};
  Configurable<bool> saveEpResoHisto{"saveEpResoHisto", false, "Flag to save event plane resolution histogram"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indexes of BDT scores to be stored. Two indexes max."};

  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {100, 1.78, 2.05}, ""};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {10, 0., 10.}, ""};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {10000, 0., 100.}, ""};
  ConfigurableAxis thnConfigAxisCosNPhi{"thnConfigAxisCosNPhi", {100, -1., 1.}, ""};
  ConfigurableAxis thnConfigAxisCosDeltaPhi{"thnConfigAxisCosDeltaPhi", {100, -1., 1.}, ""};
  ConfigurableAxis thnConfigAxisScalarProd{"thnConfigAxisScalarProd", {100, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisMlOne{"thnConfigAxisMlOne", {1000, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisMlTwo{"thnConfigAxisMlTwo", {1000, 0., 1.}, ""};

  using CandDsDataWMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDplusDataWMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using CandLcData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
  using CandLcDataWMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi>>;
  using CandD0DataWMl = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>>;
  using CandD0Data = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
  using CollsWithQvecs = soa::Join<aod::Collisions, aod::EvSels, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorBPoss, aod::QvectorBNegs, aod::QvectorBTots, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;

  Filter filterSelectDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;
  Filter filterSelectD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Filter filterSelectLcCandidates = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlag || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlag;

  Partition<CandDsData> selectedDsToKKPi = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag;
  Partition<CandDsData> selectedDsToPiKK = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Partition<CandDsDataWMl> selectedDsToKKPiWMl = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag;
  Partition<CandDsDataWMl> selectedDsToPiKKWMl = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Partition<CandD0Data> selectedD0ToPiK = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag;
  Partition<CandD0Data> selectedD0ToKPi = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Partition<CandD0DataWMl> selectedD0ToPiKWMl = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag;
  Partition<CandD0DataWMl> selectedD0ToKPiWMl = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Partition<CandLcData> selectedLcToPKPi = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlag;
  Partition<CandLcData> selectedLcToPiKP = aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlag;
  Partition<CandLcDataWMl> selectedLcToPKPiWMl = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlag;
  Partition<CandLcDataWMl> selectedLcToPiKPWMl = aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlag;

  SliceCache cache;
  HfHelper hfHelper;
  EventPlaneHelper epHelper;
  HfEventSelection hfEvSel; // event selection and monitoring
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "Inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality"};
    const AxisSpec thnAxisCosNPhi{thnConfigAxisCosNPhi, Form("cos(%d#varphi)", harmonic.value)};
    const AxisSpec thnAxisCosDeltaPhi{thnConfigAxisCosDeltaPhi, Form("cos(%d(#varphi - #Psi_{sub}))", harmonic.value)};
    const AxisSpec thnAxisScalarProd{thnConfigAxisScalarProd, "SP"};
    const AxisSpec thnAxisMlOne{thnConfigAxisMlOne, "Bkg score"};
    const AxisSpec thnAxisMlTwo{thnConfigAxisMlTwo, "FD score"};

    if (storeMl) {
      registry.add("hSparseFlowCharm", "THn for SP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCent, thnAxisCosNPhi, thnAxisCosDeltaPhi, thnAxisScalarProd, thnAxisMlOne, thnAxisMlTwo});
    } else {
      registry.add("hSparseFlowCharm", "THn for SP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCent, thnAxisCosNPhi, thnAxisCosDeltaPhi, thnAxisScalarProd});
    }
    registry.add("spReso/hSpResoFT0cFT0a", "hSpResoFT0cFT0a; centrality; Q_{FT0c} #bullet Q_{FT0a}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0cFV0a", "hSpResoFT0cFV0a; centrality; Q_{FT0c} #bullet Q_{FV0a}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0cTPCpos", "hSpResoFT0cTPCpos; centrality; Q_{FT0c} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0cTPCneg", "hSpResoFT0cTPCneg; centrality; Q_{FT0c} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0cTPCtot", "hSpResoFT0cTPCtot; centrality; Q_{FT0c} #bullet Q_{TPCtot}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0aFV0a", "hSpResoFT0aFV0a; centrality; Q_{FT0a} #bullet Q_{FV0a}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0aTPCpos", "hSpResoFT0aTPCpos; centrality; Q_{FT0a} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0aTPCneg", "hSpResoFT0aTPCneg; centrality; Q_{FT0a} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0aTPCtot", "hSpResoFT0aTPCtot; centrality; Q_{FT0m} #bullet Q_{TPCtot}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0mFV0a", "hSpResoFT0mFV0a; centrality; Q_{FT0m} #bullet Q_{FV0a}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0mTPCpos", "hSpResoFT0mTPCpos; centrality; Q_{FT0m} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0mTPCneg", "hSpResoFT0mTPCneg; centrality; Q_{FT0m} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0mTPCtot", "hSpResoFT0mTPCtot; centrality; Q_{FV0a} #bullet Q_{TPCtot}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFV0aTPCpos", "hSpResoFV0aTPCpos; centrality; Q_{FV0a} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFV0aTPCneg", "hSpResoFV0aTPCneg; centrality; Q_{FV0a} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFV0aTPCtot", "hSpResoFV0aTPCtot; centrality; Q_{FV0a} #bullet Q_{TPCtot}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoTPCposTPCneg", "hSpResoTPCposTPCneg; centrality; Q_{TPCpos} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});

    if (saveEpResoHisto) {
      registry.add("epReso/hEpResoFT0cFT0a", "hEpResoFT0cFT0a; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0cFV0a", "hEpResoFT0cFV0a; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0cTPCpos", "hEpResoFT0cTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0cTPCneg", "hEpResoFT0cTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0cTPCtot", "hEpResoFT0cTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0aFV0a", "hEpResoFT0aFV0a; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0aTPCpos", "hEpResoFT0aTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0aTPCneg", "hEpResoFT0aTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0aTPCtot", "hEpResoFT0aTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0mFV0a", "hEpResoFT0mFV0a; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0mTPCpos", "hEpResoFT0mTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0mTPCneg", "hEpResoFT0mTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0mTPCtot", "hEpResoFT0mTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFV0aTPCpos", "hEpResoFV0aTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFV0aTPCneg", "hEpResoFV0aTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFV0aTPCtot", "hEpResoFV0aTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoTPCposTPCneg", "hEpResoTPCposTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
    }

    hfEvSel.addHistograms(registry); // collision monitoring
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }; // end init

  /// Compute the Q vector for the candidate's tracks
  /// \param cand is the candidate
  /// \param tracksQx is the X component of the Q vector for the tracks
  /// \param tracksQy is the Y component of the Q vector for the tracks
  /// \param channel is the decay channel
  template <DecayChannel channel, typename T1>
  void getQvecDtracks(const T1& cand,
                      std::vector<float>& tracksQx,
                      std::vector<float>& tracksQy,
                      float& ampl)
  {
    // TODO: add possibility to consider different weights for the tracks, at the moment only pT is considered;
    float pXTrack0 = cand.pxProng0();
    float pYTrack0 = cand.pyProng0();
    float pTTrack0 = cand.ptProng0();
    float phiTrack0 = std::atan2(pYTrack0, pXTrack0);
    float pXTrack1 = cand.pxProng1();
    float pYTrack1 = cand.pyProng1();
    float pTTrack1 = cand.ptProng1();
    float phiTrack1 = std::atan2(pYTrack1, pXTrack1);

    tracksQx.push_back(std::cos(harmonic * phiTrack0) * pTTrack0 / ampl);
    tracksQy.push_back(std::sin(harmonic * phiTrack0) * pTTrack0 / ampl);
    tracksQx.push_back(std::cos(harmonic * phiTrack1) * pTTrack1 / ampl);
    tracksQy.push_back(std::sin(harmonic * phiTrack1) * pTTrack1 / ampl);

    if constexpr (channel != DecayChannel::D0ToPiK && channel != DecayChannel::D0ToKPi) {
      float pXTrack2 = cand.pxProng2();
      float pYTrack2 = cand.pyProng2();
      float pTTrack2 = cand.ptProng2();
      float phiTrack2 = std::atan2(pYTrack2, pXTrack2);
      tracksQx.push_back(std::cos(harmonic * phiTrack2) * pTTrack2 / ampl);
      tracksQy.push_back(std::sin(harmonic * phiTrack2) * pTTrack2 / ampl);
    }
  }

  /// Compute the delta psi in the range [0, pi/harmonic]
  /// \param psi1 is the first angle
  /// \param psi2 is the second angle
  /// \note Ported from AliAnalysisTaskSECharmHadronvn::GetDeltaPsiSubInRange
  float getDeltaPsiInRange(float psi1, float psi2)
  {
    float deltaPsi = psi1 - psi2;
    if (std::abs(deltaPsi) > constants::math::PI / harmonic) {
      if (deltaPsi > 0.)
        deltaPsi -= constants::math::TwoPI / harmonic;
      else
        deltaPsi += constants::math::TwoPI / harmonic;
    }
    return deltaPsi;
  }

  /// Fill THnSparse
  /// \param mass is the invariant mass of the candidate
  /// \param pt is the transverse momentum of the candidate
  /// \param cent is the centrality of the collision
  /// \param cosNPhi is the cosine of the n*phi angle
  /// \param cosDeltaPhi is the cosine of the n*(phi - evtPl) angle
  /// \param sp is the scalar product
  /// \param outputMl are the ML scores
  void fillThn(float& mass,
               float& pt,
               float& cent,
               float& cosNPhi,
               float& cosDeltaPhi,
               float& sp,
               std::vector<float>& outputMl)
  {
    if (storeMl) {
      registry.fill(HIST("hSparseFlowCharm"), mass, pt, cent, cosNPhi, cosDeltaPhi, sp, outputMl[0], outputMl[1]);
    } else {
      registry.fill(HIST("hSparseFlowCharm"), mass, pt, cent, cosNPhi, cosDeltaPhi, sp);
    }
  }

  /// Get the centrality
  /// \param collision is the collision with the centrality information
  float getCentrality(CollsWithQvecs::iterator const& collision)
  {
    float cent = -999.;
    switch (centEstimator) {
      case CentralityEstimator::FV0A:
        cent = collision.centFV0A();
        break;
      case CentralityEstimator::FT0M:
        cent = collision.centFT0M();
        break;
      case CentralityEstimator::FT0A:
        cent = collision.centFT0A();
        break;
      case CentralityEstimator::FT0C:
        cent = collision.centFT0C();
        break;
      default:
        LOG(warning) << "Centrality estimator not valid. Possible values are V0A, T0M, T0A, T0C. Fallback to V0A";
        cent = collision.centFV0A();
        break;
    }
    return cent;
  }

  /// Check if the collision is selected
  /// \param collision is the collision with the Q vector information
  /// \param bc is the bunch crossing with timestamp information
  /// \return true if the collision is selected, false otherwise
  template <o2::hf_centrality::CentralityEstimator centEstimator>
  bool isCollSelected(CollsWithQvecs::iterator const& collision,
                      aod::BCsWithTimestamps const&)
  {
    float centrality{-1.f};
    const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

    /// monitor the satisfied event selections
    hfEvSel.fillHistograms(collision, rejectionMask, centrality);
    return rejectionMask == 0;
  }

  /// Get the Q vector
  /// \param collision is the collision with the Q vector information
  std::vector<float> getQvec(CollsWithQvecs::iterator const& collision)
  {
    float xQVec = -999.;
    float yQVec = -999.;
    float amplQVec = -999.;
    switch (qvecDetector) {
      case QvecEstimator::FV0A:
        xQVec = collision.qvecFV0ARe();
        yQVec = collision.qvecFV0AIm();
        break;
      case QvecEstimator::FT0M:
        xQVec = collision.qvecFT0MRe();
        yQVec = collision.qvecFT0MIm();
        break;
      case QvecEstimator::FT0A:
        xQVec = collision.qvecFT0ARe();
        yQVec = collision.qvecFT0AIm();
        break;
      case QvecEstimator::FT0C:
        xQVec = collision.qvecFT0CRe();
        yQVec = collision.qvecFT0CIm();
        break;
      case QvecEstimator::TPCPos:
        xQVec = collision.qvecBPosRe();
        yQVec = collision.qvecBPosIm();
        amplQVec = collision.nTrkBPos();
        break;
      case QvecEstimator::TPCNeg:
        xQVec = collision.qvecBNegRe();
        yQVec = collision.qvecBNegIm();
        amplQVec = collision.nTrkBNeg();
        break;
      case QvecEstimator::TPCTot:
        xQVec = collision.qvecBTotRe();
        yQVec = collision.qvecBTotIm();
        amplQVec = collision.nTrkBTot();
        break;
      default:
        LOG(warning) << "Q vector estimator not valid. Please choose between FV0A, FT0M, FT0A, FT0C, TPC Pos, TPC Neg. Fallback to FV0A";
        xQVec = collision.qvecFV0ARe();
        yQVec = collision.qvecFV0AIm();
        break;
    }
    return {xQVec, yQVec, amplQVec};
  }

  /// Compute the scalar product
  /// \param collision is the collision with the Q vector information and event plane
  /// \param candidates are the selected candidates
  template <DecayChannel channel, typename T1>
  void runFlowAnalysis(CollsWithQvecs::iterator const& collision,
                       T1 const& candidates)
  {
    std::vector<float> qVecs = getQvec(collision);
    float xQVec = qVecs[0];
    float yQVec = qVecs[1];
    float amplQVec = qVecs[2];
    float evtPl = epHelper.GetEventPlane(xQVec, yQVec, harmonic);
    float cent = getCentrality(collision);
    int nProngs = 3;

    for (const auto& candidate : candidates) {
      float massCand = 0.;
      std::vector<float> outputMl = {-999., -999.};

      if constexpr (std::is_same_v<T1, CandDsData> || std::is_same_v<T1, CandDsDataWMl>) {
        switch (channel) {
          case DecayChannel::DsToKKPi:
            massCand = hfHelper.invMassDsToKKPi(candidate);
            if constexpr (std::is_same_v<T1, CandDsDataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
                outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
            }
            break;
          case DecayChannel::DsToPiKK:
            massCand = hfHelper.invMassDsToPiKK(candidate);
            if constexpr (std::is_same_v<T1, CandDsDataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
                outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
            }
            break;
          default:
            break;
        }
      } else if constexpr (std::is_same_v<T1, CandDplusData> || std::is_same_v<T1, CandDplusDataWMl>) {
        massCand = hfHelper.invMassDplusToPiKPi(candidate);
        if constexpr (std::is_same_v<T1, CandDplusDataWMl>) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
            outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
        }
      } else if constexpr (std::is_same_v<T1, CandD0Data> || std::is_same_v<T1, CandD0DataWMl>) {
        nProngs = 2;
        switch (channel) {
          case DecayChannel::D0ToPiK:
            massCand = hfHelper.invMassD0ToPiK(candidate);
            if constexpr (std::is_same_v<T1, CandD0DataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
                outputMl[iclass] = candidate.mlProbD0()[classMl->at(iclass)];
            }
            break;
          case DecayChannel::D0ToKPi:
            massCand = hfHelper.invMassD0barToKPi(candidate);
            if constexpr (std::is_same_v<T1, CandD0DataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
                outputMl[iclass] = candidate.mlProbD0bar()[classMl->at(iclass)];
            }
            break;
          default:
            break;
        }
      } else if constexpr (std::is_same_v<T1, CandLcData> || std::is_same_v<T1, CandLcDataWMl>) {
        switch (channel) {
          case DecayChannel::LcToPKPi:
            massCand = hfHelper.invMassLcToPKPi(candidate);
            if constexpr (std::is_same_v<T1, CandLcDataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
                outputMl[iclass] = candidate.mlProbLcToPKPi()[classMl->at(iclass)];
            }
            break;
          case DecayChannel::LcToPiKP:
            massCand = hfHelper.invMassLcToPiKP(candidate);
            if constexpr (std::is_same_v<T1, CandLcDataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
                outputMl[iclass] = candidate.mlProbLcToPiKP()[classMl->at(iclass)];
            }
            break;
          default:
            break;
        }
      }

      float ptCand = candidate.pt();
      float phiCand = candidate.phi();

      // If TPC is used for the SP estimation, the tracks of the hadron candidate must be removed from the TPC Q vector to avoid double counting
      if (qvecDetector == QvecEstimator::TPCNeg || qvecDetector == QvecEstimator::TPCPos) {
        float ampl = amplQVec - static_cast<float>(nProngs);
        std::vector<float> tracksQx = {};
        std::vector<float> tracksQy = {};

        getQvecDtracks<channel>(candidate, tracksQx, tracksQy, ampl);
        for (auto iTrack{0u}; iTrack < tracksQx.size(); ++iTrack) {
          xQVec -= tracksQx[iTrack];
          yQVec -= tracksQy[iTrack];
        }
      }

      float cosNPhi = std::cos(harmonic * phiCand);
      float sinNPhi = std::sin(harmonic * phiCand);
      float scalprodCand = cosNPhi * xQVec + sinNPhi * yQVec;
      float cosDeltaPhi = std::cos(harmonic * (phiCand - evtPl));

      fillThn(massCand, ptCand, cent, cosNPhi, cosDeltaPhi, scalprodCand, outputMl);
    }
  }

  // Ds with ML
  void processDsMl(CollsWithQvecs::iterator const& collision,
                   CandDsDataWMl const&)
  {
    auto candsDsToKKPiWMl = selectedDsToKKPiWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsDsToPiKKWMl = selectedDsToPiKKWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::DsToKKPi>(collision, candsDsToKKPiWMl);
    runFlowAnalysis<DecayChannel::DsToPiKK>(collision, candsDsToPiKKWMl);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDsMl, "Process Ds candidates with ML", false);

  // Ds with rectangular cuts
  void processDs(CollsWithQvecs::iterator const& collision,
                 CandDsData const&)
  {
    auto candsDsToKKPi = selectedDsToKKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsDsToPiKK = selectedDsToPiKK->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::DsToKKPi>(collision, candsDsToKKPi);
    runFlowAnalysis<DecayChannel::DsToPiKK>(collision, candsDsToPiKK);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDs, "Process Ds candidates", false);

  // Dplus with ML
  void processDplusMl(CollsWithQvecs::iterator const& collision,
                      CandDplusDataWMl const& candidatesDplus)
  {
    runFlowAnalysis<DecayChannel::DplusToPiKPi>(collision, candidatesDplus);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDplusMl, "Process Dplus candidates with ML", false);

  // Dplus with rectangular cuts
  void processDplus(CollsWithQvecs::iterator const& collision,
                    CandDplusData const& candidatesDplus)
  {
    runFlowAnalysis<DecayChannel::DplusToPiKPi>(collision, candidatesDplus);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDplus, "Process Dplus candidates", true);

  // D0 with ML
  void processD0Ml(CollsWithQvecs::iterator const& collision,
                   CandD0DataWMl const&)
  {
    auto candsD0ToPiKWMl = selectedD0ToPiKWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsD0ToKPiWMl = selectedD0ToKPiWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::D0ToPiK>(collision, candsD0ToPiKWMl);
    runFlowAnalysis<DecayChannel::D0ToKPi>(collision, candsD0ToKPiWMl);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processD0Ml, "Process D0 candidates with ML", false);

  // D0 with rectangular cuts
  void processD0(CollsWithQvecs::iterator const& collision,
                 CandD0Data const&)
  {
    auto candsD0ToPiK = selectedD0ToPiK->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsD0ToKPi = selectedD0ToKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::D0ToPiK>(collision, candsD0ToPiK);
    runFlowAnalysis<DecayChannel::D0ToKPi>(collision, candsD0ToKPi);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processD0, "Process D0 candidates", false);

  // Lc with ML
  void processLcMl(CollsWithQvecs::iterator const& collision,
                   CandLcDataWMl const&)
  {
    auto candsLcToPKPiWMl = selectedLcToPKPiWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsLcToPiKPWMl = selectedLcToPiKPWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::LcToPKPi>(collision, candsLcToPKPiWMl);
    runFlowAnalysis<DecayChannel::LcToPiKP>(collision, candsLcToPiKPWMl);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processLcMl, "Process Lc candidates with ML", false);

  // Lc with rectangular cuts
  void processLc(CollsWithQvecs::iterator const& collision,
                 CandLcData const&)
  {
    auto candsLcToPKPi = selectedLcToPKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsLcToPiKP = selectedLcToPiKP->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::LcToPKPi>(collision, candsLcToPKPi);
    runFlowAnalysis<DecayChannel::LcToPiKP>(collision, candsLcToPiKP);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processLc, "Process Lc candidates", false);

  // Resolution
  void processResolution(CollsWithQvecs::iterator const& collision,
                         aod::BCsWithTimestamps const& bcs)
  {
    if (!isCollSelected<o2::hf_centrality::CentralityEstimator::None>(collision, bcs)) {
      // no selection on the centrality is applied on purpose to allow for the resolution study in post-processing
      return;
    }

    float centrality = getCentrality(collision); // centrality not updated in the rejection mask function
    float xQVecFT0a = collision.qvecFT0ARe();
    float yQVecFT0a = collision.qvecFT0AIm();
    float xQVecFT0c = collision.qvecFT0CRe();
    float yQVecFT0c = collision.qvecFT0CIm();
    float xQVecFT0m = collision.qvecFT0MRe();
    float yQVecFT0m = collision.qvecFT0MIm();
    float xQVecFV0a = collision.qvecFV0ARe();
    float yQVecFV0a = collision.qvecFV0AIm();
    float xQVecBPos = collision.qvecBPosRe();
    float yQVecBPos = collision.qvecBPosIm();
    float xQVecBNeg = collision.qvecBNegRe();
    float yQVecBNeg = collision.qvecBNegIm();
    float xQVecBTot = collision.qvecBTotRe();
    float yQVecBTot = collision.qvecBTotIm();

    registry.fill(HIST("spReso/hSpResoFT0cFT0a"), centrality, xQVecFT0c * xQVecFT0a + yQVecFT0c * yQVecFT0a);
    registry.fill(HIST("spReso/hSpResoFT0cFV0a"), centrality, xQVecFT0c * xQVecFV0a + yQVecFT0c * yQVecFV0a);
    registry.fill(HIST("spReso/hSpResoFT0cTPCpos"), centrality, xQVecFT0c * xQVecBPos + yQVecFT0c * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFT0cTPCneg"), centrality, xQVecFT0c * xQVecBNeg + yQVecFT0c * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFT0cTPCtot"), centrality, xQVecFT0c * xQVecBTot + yQVecFT0c * yQVecBTot);
    registry.fill(HIST("spReso/hSpResoFT0aFV0a"), centrality, xQVecFT0a * xQVecFV0a + yQVecFT0a * yQVecFV0a);
    registry.fill(HIST("spReso/hSpResoFT0aTPCpos"), centrality, xQVecFT0a * xQVecBPos + yQVecFT0a * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFT0aTPCneg"), centrality, xQVecFT0a * xQVecBNeg + yQVecFT0a * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFT0aTPCtot"), centrality, xQVecFT0a * xQVecBTot + yQVecFT0a * yQVecBTot);
    registry.fill(HIST("spReso/hSpResoFT0mFV0a"), centrality, xQVecFT0m * xQVecFV0a + yQVecFT0m * yQVecFV0a);
    registry.fill(HIST("spReso/hSpResoFT0mTPCpos"), centrality, xQVecFT0m * xQVecBPos + yQVecFT0m * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFT0mTPCneg"), centrality, xQVecFT0m * xQVecBNeg + yQVecFT0m * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFT0mTPCtot"), centrality, xQVecFT0m * xQVecBTot + yQVecFT0m * yQVecBTot);
    registry.fill(HIST("spReso/hSpResoFV0aTPCpos"), centrality, xQVecFV0a * xQVecBPos + yQVecFV0a * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFV0aTPCneg"), centrality, xQVecFV0a * xQVecBNeg + yQVecFV0a * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFV0aTPCtot"), centrality, xQVecFV0a * xQVecBTot + yQVecFV0a * yQVecBTot);
    registry.fill(HIST("spReso/hSpResoTPCposTPCneg"), centrality, xQVecBPos * xQVecBNeg + yQVecBPos * yQVecBNeg);

    if (saveEpResoHisto) {
      float epFT0a = epHelper.GetEventPlane(xQVecFT0a, yQVecFT0a, harmonic);
      float epFT0c = epHelper.GetEventPlane(xQVecFT0c, yQVecFT0c, harmonic);
      float epFT0m = epHelper.GetEventPlane(xQVecFT0m, yQVecFT0m, harmonic);
      float epFV0a = epHelper.GetEventPlane(xQVecFV0a, yQVecFV0a, harmonic);
      float epBPoss = epHelper.GetEventPlane(xQVecBPos, yQVecBPos, harmonic);
      float epBNegs = epHelper.GetEventPlane(xQVecBNeg, yQVecBNeg, harmonic);
      float epBTots = epHelper.GetEventPlane(xQVecBTot, yQVecBTot, harmonic);

      registry.fill(HIST("epReso/hEpResoFT0cFT0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epFT0a)));
      registry.fill(HIST("epReso/hEpResoFT0cFV0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epFV0a)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBTots)));
      registry.fill(HIST("epReso/hEpResoFT0aFV0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epFV0a)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBTots)));
      registry.fill(HIST("epReso/hEpResoFT0mFV0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epFV0a)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBTots)));
      registry.fill(HIST("epReso/hEpResoFV0aTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFV0a, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFV0aTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFV0a, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFV0aTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFV0a, epBTots)));
      registry.fill(HIST("epReso/hEpResoTPCposTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epBPoss, epBNegs)));
    }
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processResolution, "Process resolution", false);

}; // End struct HfTaskFlowCharmHadrons

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskFlowCharmHadrons>(cfgc)};
}
