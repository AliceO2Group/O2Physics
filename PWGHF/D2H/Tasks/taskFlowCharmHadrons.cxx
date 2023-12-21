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

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Qvectors.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum DecayChannel { DplusToPiKPi = 0,
                    DsToKKPi,
                    DsToPiKK };

enum centralityEstimator { V0A = 0,
                           T0M,
                           T0A,
                           T0C };

enum qvecEstimator { FV0A = 0,
                     FT0M,
                     FT0A,
                     FT0C,
                     TPCPos,
                     TPCNeg };

struct HfTaskFlowCharmHadrons {
  Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<int> qvecDetector{"qvecDetector", 0, "Detector for Q vector estimation (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3, TPC Pos: 4, TPC Neg: 5)"};
  Configurable<int> centDetector{"centDetector", 0, "Detector for centrality estimation (V0A: 0, T0M: 1, T0A: 2, T0C: 3"};
  Configurable<int> selectionFlag{"selectionFlag", 1, "Selection Flag for hadron (e.g. 1 for skimming, 3 for topo. and kine., 7 for PID)"};
  Configurable<bool> storeMl{"storeMl", false, "Flag to store ML scores"};
  Configurable<bool> saveEpResoHisto{"saveEpResoHisto", false, "Flag to save event plane resolution histogram"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indexes of BDT scores to be stored. Two indexes max."};

  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {100, 1.78, 2.05}, ""};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {10, 0., 10.}, ""};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {10000, 0., 100.}, ""};
  ConfigurableAxis thnConfigAxisCosNPhi{"thnConfigAxisCosNPhi", {100, -1., 1.}, ""};
  ConfigurableAxis thnConfigAxisScalarProd{"thnConfigAxisScalarProd", {100, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisMlOne{"thnConfigAxisMlOne", {1000, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisMlTwo{"thnConfigAxisMlTwo", {1000, 0., 1.}, ""};

  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CollsWithQvecs = soa::Join<aod::Collisions, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorBPoss, aod::QvectorBNegs, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;

  Filter filterSelectDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;

  Partition<CandDsData> selectedDsToKKPi = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag;
  Partition<CandDsData> selectedDsToPiKK = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;

  HfHelper hfHelper;
  EventPlaneHelper epHelper;

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "Inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality"};
    const AxisSpec thnAxisCosNPhi{thnConfigAxisCosNPhi, Form("cos(%d#varphi)", harmonic.value)};
    const AxisSpec thnAxisScalarProd{thnConfigAxisScalarProd, "SP"};
    const AxisSpec thnAxisMlOne{thnConfigAxisMlOne, "Bkg score"};
    const AxisSpec thnAxisMlTwo{thnConfigAxisMlTwo, "FD score"};

    registry.add("hSparseFlowCharm", "THn for SP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCent, thnAxisCosNPhi, thnAxisScalarProd, thnAxisMlOne, thnAxisMlTwo});
    registry.add("spReso/hSpResoFT0cFT0a", "hSpResoFT0cFT0a; centrality; Q_{FT0c} #bullet Q_{FT0a}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0cFT0m", "hSpResoFT0cFT0m; centrality; Q_{FT0c} #bullet Q_{FT0m}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0cFV0m", "hSpResoFT0cFV0m; centrality; Q_{FT0c} #bullet Q_{FV0m}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0cTPCpos", "hSpResoFT0cTPCpos; centrality; Q_{FT0c} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0cTPCneg", "hSpResoFT0cTPCneg; centrality; Q_{FT0c} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0aFT0m", "hSpResoFT0aFT0m; centrality; Q_{FT0a} #bullet Q_{FT0m}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0aFV0m", "hSpResoFT0aFV0m; centrality; Q_{FT0a} #bullet Q_{FV0m}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0aTPCpos", "hSpResoFT0aTPCpos; centrality; Q_{FT0a} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0aTPCneg", "hSpResoFT0aTPCneg; centrality; Q_{FT0a} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0mFV0m", "hSpResoFT0mFV0m; centrality; Q_{FT0m} #bullet Q_{FV0m}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0mTPCpos", "hSpResoFT0mTPCpos; centrality; Q_{FT0m} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0mTPCneg", "hSpResoFT0mTPCneg; centrality; Q_{FT0m} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFV0mTPCpos", "hSpResoFV0mTPCpos; centrality; Q_{FV0m} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFV0mTPCneg", "hSpResoFV0mTPCneg; centrality; Q_{FV0m} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});

    if (saveEpResoHisto) {
      registry.add("epReso/hEpResoFT0cFT0a", "hEpResoFT0cFT0a; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0cFT0m", "hEpResoFT0cFT0m; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0cFV0m", "hEpResoFT0cFV0m; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0cTPCpos", "hEpResoFT0cTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0cTPCneg", "hEpResoFT0cTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0aFT0m", "hEpResoFT0aFT0m; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0aFV0m", "hEpResoFT0aFV0m; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0aTPCpos", "hEpResoFT0aTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0aTPCneg", "hEpResoFT0aTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0mFV0m", "hEpResoFT0mFV0m; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0mTPCpos", "hEpResoFT0mTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0mTPCneg", "hEpResoFT0mTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFV0mTPCpos", "hEpResoFV0mTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFV0mTPCneg", "hEpResoFV0mTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
    }
  }; // end init

  /// Compute the Q vector for the candidate's tracks
  /// \param cand is the candidate
  /// \param tracksQx is the X component of the Q vector for the tracks
  /// \param tracksQy is the Y component of the Q vector for the tracks
  template <typename T1>
  void getQvecDtracks(const T1& cand,
                      std::vector<float>& tracksQx,
                      std::vector<float>& tracksQy,
                      float& ampl)
  {
    // TODO: add possibility to consider different weights for the tracks, at the moment only pT is considered;
    float pXtrack0 = cand.pxProng0();
    float pYtrack0 = cand.pyProng0();
    float pTtrack0 = cand.ptProng0();
    float phiTrack0 = TMath::ATan2(pYtrack0, pXtrack0);
    float pXtrack1 = cand.pxProng1();
    float pYtrack1 = cand.pyProng1();
    float pTtrack1 = cand.ptProng1();
    float phiTrack1 = TMath::ATan2(pYtrack1, pXtrack1);
    float pXtrack2 = cand.pxProng2();
    float pYtrack2 = cand.pyProng2();
    float pTtrack2 = cand.ptProng2();
    float phiTrack2 = TMath::ATan2(pYtrack2, pXtrack2);

    tracksQx.push_back(cos(harmonic * phiTrack0) * pTtrack0 / ampl);
    tracksQy.push_back(sin(harmonic * phiTrack0) * pTtrack0 / ampl);
    tracksQx.push_back(cos(harmonic * phiTrack1) * pTtrack1 / ampl);
    tracksQy.push_back(sin(harmonic * phiTrack1) * pTtrack1 / ampl);
    tracksQx.push_back(cos(harmonic * phiTrack2) * pTtrack2 / ampl);
    tracksQy.push_back(sin(harmonic * phiTrack2) * pTtrack2 / ampl);
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
  /// \param evtPlReso is the event plane resolution
  /// \param outputMl are the ML scores
  void fillThn(float mass,
               float pt,
               float cent,
               float cosNPhi,
               float cosDeltaPhi,
               float sp,
               std::vector<float> outputMl)
  {
    registry.fill(HIST("hSparseFlowCharm"), mass, pt, cent, cosNPhi, cosDeltaPhi, sp, outputMl[0], outputMl[1]);
  }

  /// Get the centrality
  /// \param collision is the collision with the centrality information
  float getCentrality(CollsWithQvecs::iterator const& collision)
  {
    float cent = -999.;
    switch (centDetector) {
      case centralityEstimator::V0A:
        cent = collision.centFV0A();
        break;
      case centralityEstimator::T0M:
        cent = collision.centFT0M();
        break;
      case centralityEstimator::T0A:
        cent = collision.centFT0A();
        break;
      case centralityEstimator::T0C:
        cent = collision.centFT0C();
        break;
      default:
        LOG(warning) << "Centrality estimator not valid. Possible values are V0A, T0M, T0A, T0C. Fallback to V0A";
        cent = collision.centFV0A();
        break;
    }
    return cent;
  }

  /// Get the Q vector
  /// \param collision is the collision with the Q vector information
  std::vector<float> getQvec(CollsWithQvecs::iterator const& collision)
  {
    float xQVec = -999.;
    float yQVec = -999.;
    float amplQVec = -999.;
    switch (qvecDetector) {
      case qvecEstimator::FV0A:
        xQVec = collision.qvecFV0ARe();
        yQVec = collision.qvecFV0AIm();
        break;
      case qvecEstimator::FT0M:
        xQVec = collision.qvecFT0MRe();
        yQVec = collision.qvecFT0MIm();
        break;
      case qvecEstimator::FT0A:
        xQVec = collision.qvecFT0ARe();
        yQVec = collision.qvecFT0AIm();
        break;
      case qvecEstimator::FT0C:
        xQVec = collision.qvecFT0CRe();
        yQVec = collision.qvecFT0CIm();
      case qvecEstimator::TPCPos:
        xQVec = collision.qvecBPosRe();
        yQVec = collision.qvecBPosIm();
        amplQVec = collision.nTrkBPos();
        break;
      case qvecEstimator::TPCNeg:
        xQVec = collision.qvecBNegRe();
        yQVec = collision.qvecBNegIm();
        amplQVec = collision.nTrkBNeg();
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
  template <int DecayChannel, typename T1>
  void runFlowAnalysis(CollsWithQvecs::iterator const& collision,
                       T1 const& candidates)
  {
    std::vector<float> qVecs = getQvec(collision);
    float xQVec = qVecs[0];
    float yQVec = qVecs[1];
    float amplQVec = qVecs[2];
    float evtPl = epHelper.GetEventPlane(xQVec, yQVec, harmonic);
    float cent = getCentrality(collision);

    for (const auto& candidate : candidates) {
      float massCand = 0.;
      std::vector<float> outputMl = {-999., -999.};

      if constexpr (std::is_same<T1, CandDsData>::value) {
        switch (DecayChannel) {
          case DecayChannel::DsToKKPi:
            massCand = hfHelper.invMassDsToKKPi(candidate);
            if (storeMl) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
                outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
            }
            break;
          case DecayChannel::DsToPiKK:
            massCand = hfHelper.invMassDsToPiKK(candidate);
            if (storeMl) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
                outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
            }
            break;
          default:
            break;
        }
      } else if constexpr (std::is_same<T1, CandDplusData>::value) {
        massCand = hfHelper.invMassDplusToPiKPi(candidate);
        if (storeMl) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
            outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
        }
      }
      float ptCand = candidate.pt();
      float phiCand = candidate.phi();

      // If TPC is used for the SP estimation, the tracks of the hadron candidate must be removed from the TPC Q vector to avoid double counting
      if (qvecDetector == 4 || qvecDetector == 5) {
        float ampl = amplQVec - 3.;
        std::vector<float> tracksQx = {};
        std::vector<float> tracksQy = {};
        getQvecDtracks(candidate, tracksQx, tracksQy, ampl);
        for (unsigned int itrack = 0; itrack < 3; itrack++) {
          xQVec -= tracksQx[itrack];
          yQVec -= tracksQy[itrack];
        }
      }

      float cosNPhi = TMath::Cos(harmonic * phiCand);
      float sinNPhi = TMath::Sin(harmonic * phiCand);
      float scalprodCand = cosNPhi * xQVec + sinNPhi * yQVec;
      float cosDeltaPhi = TMath::Cos(harmonic * (phiCand - evtPl));

      fillThn(massCand, ptCand, cent, cosNPhi, cosDeltaPhi, scalprodCand, outputMl);
    }
  }

  // Ds
  void processDs(CollsWithQvecs::iterator const& collision,
                 CandDsData const& candidatesDs)
  {
    runFlowAnalysis<DecayChannel::DsToKKPi, Partition<CandDsData>>(collision, selectedDsToKKPi);
    runFlowAnalysis<DecayChannel::DsToPiKK, Partition<CandDsData>>(collision, selectedDsToPiKK);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDs, "Process Ds candidates", false);

  // Dplus
  void processDplus(CollsWithQvecs::iterator const& collision,
                    CandDplusData const& candidatesDplus)
  {
    runFlowAnalysis<DecayChannel::DplusToPiKPi, CandDplusData>(collision, candidatesDplus);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDplus, "Process Dplus candidates", true);

  // Resolution
  void processResolution(CollsWithQvecs::iterator const& collision)
  {
    float centrality = getCentrality(collision);
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

    registry.fill(HIST("spReso/hSpResoFT0cFT0a"), centrality, xQVecFT0c * xQVecFT0a + yQVecFT0c * yQVecFT0a);
    registry.fill(HIST("spReso/hSpResoFT0cFT0m"), centrality, xQVecFT0c * xQVecFT0m + yQVecFT0c * yQVecFT0m);
    registry.fill(HIST("spReso/hSpResoFT0cFV0m"), centrality, xQVecFT0c * xQVecFV0a + yQVecFT0c * yQVecFV0a);
    registry.fill(HIST("spReso/hSpResoFT0cTPCpos"), centrality, xQVecFT0c * xQVecBPos + yQVecFT0c * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFT0cTPCneg"), centrality, xQVecFT0c * xQVecBNeg + yQVecFT0c * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFT0aFT0m"), centrality, xQVecFT0a * xQVecFT0m + yQVecFT0a * yQVecFT0m);
    registry.fill(HIST("spReso/hSpResoFT0aFV0m"), centrality, xQVecFT0a * xQVecFV0a + yQVecFT0a * yQVecFV0a);
    registry.fill(HIST("spReso/hSpResoFT0aTPCpos"), centrality, xQVecFT0a * xQVecBPos + yQVecFT0a * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFT0aTPCneg"), centrality, xQVecFT0a * xQVecBNeg + yQVecFT0a * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFT0mFV0m"), centrality, xQVecFT0m * xQVecFV0a + yQVecFT0m * yQVecFV0a);
    registry.fill(HIST("spReso/hSpResoFT0mTPCpos"), centrality, xQVecFT0m * xQVecBPos + yQVecFT0m * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFT0mTPCneg"), centrality, xQVecFT0m * xQVecBNeg + yQVecFT0m * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFV0mTPCpos"), centrality, xQVecFV0a * xQVecBPos + yQVecFV0a * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFV0mTPCneg"), centrality, xQVecFV0a * xQVecBNeg + yQVecFV0a * yQVecBNeg);

    if (saveEpResoHisto) {
      float epFT0a = epHelper.GetEventPlane(xQVecFT0a, yQVecFT0a, harmonic);
      float epFT0c = epHelper.GetEventPlane(xQVecFT0c, yQVecFT0c, harmonic);
      float epFT0m = epHelper.GetEventPlane(xQVecFT0m, yQVecFT0m, harmonic);
      float epFV0a = epHelper.GetEventPlane(xQVecFV0a, yQVecFV0a, harmonic);
      float epBPoss = epHelper.GetEventPlane(xQVecBPos, yQVecBPos, harmonic);
      float epBNegs = epHelper.GetEventPlane(xQVecBNeg, yQVecBNeg, harmonic);

      registry.fill(HIST("epReso/hEpResoFT0cFT0a"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFT0c, epFT0a)));
      registry.fill(HIST("epReso/hEpResoFT0cFT0m"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFT0c, epFT0m)));
      registry.fill(HIST("epReso/hEpResoFT0cFV0m"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFT0c, epFV0a)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCpos"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFT0c, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCneg"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFT0c, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0aFT0m"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFT0a, epFT0m)));
      registry.fill(HIST("epReso/hEpResoFT0aFV0m"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFT0a, epFV0a)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCpos"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFT0a, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCneg"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFT0a, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0mFV0m"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFT0m, epFV0a)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCpos"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFT0m, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCneg"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFT0m, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFV0mTPCpos"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFV0a, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFV0mTPCneg"), centrality, TMath::Cos(harmonic * getDeltaPsiInRange(epFV0a, epBNegs)));
    }
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processResolution, "Process resolution", false);

}; // End struct HfTaskFlowCharmHadrons

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskFlowCharmHadrons>(cfgc)};
}
