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

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/DataModel/Qvectors.h"
#include "Common/Core/EventPlaneHelper.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum decayChannel { DplusToPiKPi = 0,
                    DsToKKPi,
                    DsToPiKK };

struct taskFlowCharmHadrons {
  Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<int> qvecDetector{"qvecDetector", 0, "Detector for Q vector estimation (FT0C: 0, FT0A: 1, FT0M: 2, FV0A: 3, TPC Pos: 4, TPC Neg: 5)"};
  Configurable<int> centDetector{"centDetector", 0, "Detector for centrality estimation (V0A: 0, T0M: 1, T0A: 2, T0C: 3"};
  Configurable<int> selectionFlag{"selectionFlag", 1, "Selection Flag for D"};
  Configurable<bool> storeMl{"storeMl", false, "Flag to store ML score"};
  Configurable<bool> saveEpResoHisto{"saveEpResoHisto", false, "Flag to save event plane resolution histogram"};
  Configurable<int> classMl{"classMl", 0, "Index of the ML class to be stored"};
  
  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {100, 1.78, 2.05}, ""};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {10, 0., 10.}, ""};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {10000, 0., 100.}, ""};
  ConfigurableAxis thnConfigAxisCosNPhi{"thnConfigAxisCosNPhi", {100, -1., 1.}, ""};
  ConfigurableAxis thnConfigAxisScalarProd{"thnConfigAxisScalarProd", {100, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisMl{"thnConfigAxisMl", {1000, 0., 1.}, ""};

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

  void init(InitContext&) {
    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "Inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality"};
    const AxisSpec thnAxisCosNPhi{thnConfigAxisCosNPhi, Form("cos(%d#varphi)", harmonic.value)};
    const AxisSpec thnAxisScalarProd{thnConfigAxisScalarProd, "SP"};
    const AxisSpec thnAxisMlScore{thnConfigAxisMl, "ML score"};

    registry.add("hSparseFlowCharm", "THn for SP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCent, thnAxisCosNPhi, thnAxisScalarProd, thnAxisMlScore});

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
  void GetQvecDtracks(const T1& cand,
                      std::vector<float>& tracksQx,
                      std::vector<float>& tracksQy,
                      float& ampl) 
  {
    // TODO: add possibility to consider different weights for the tracks, at the only pT is considered;
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
  
    tracksQx.push_back(TMath::Cos(harmonic * phiTrack0) * pTtrack0 / ampl);
    tracksQy.push_back(TMath::Sin(harmonic * phiTrack0) * pTtrack0 / ampl);
    tracksQx.push_back(TMath::Cos(harmonic * phiTrack1) * pTtrack1 / ampl);
    tracksQy.push_back(TMath::Sin(harmonic * phiTrack1) * pTtrack1 / ampl);
    tracksQx.push_back(TMath::Cos(harmonic * phiTrack2) * pTtrack2 / ampl);
    tracksQy.push_back(TMath::Sin(harmonic * phiTrack2) * pTtrack2 / ampl);
  }


  /// Compute the delta psi in the range [0, pi/harmonic]
  /// \param psi1 is the first angle
  /// \param psi2 is the second angle
  /// \note Ported from AliAnalysisTaskSECharmHadronvn::GetDeltaPsiSubInRange
  float GetDeltaPsiInRange(float psi1, float psi2)
  {
    float deltaPsi = psi1 - psi2;
    if (deltaPsi > TMath::Pi() / harmonic) {
      if (deltaPsi > 0.) deltaPsi -= 2. * TMath::Pi() / harmonic;
      else deltaPsi += 2. * TMath::Pi() / harmonic;
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
  /// \param outputMl is the ML score
  void fillThn(float mass,
               float pt,
               float cent,
               float cosNPhi,
               float cosDeltaPhi,
               float sp,
               float outputMl)
  {
    registry.fill(HIST("hSparseFlowCharm"), mass, pt, cent, cosNPhi, cosDeltaPhi, sp, outputMl);
  }


  /// Get the centrality
  /// \param collision is the collision with the centrality information
  /// \param centDetector is the detector used for the centrality estimation
  float GetCent(CollsWithQvecs::iterator const& collision, int centDetector)
  {
    float cent = -999.;
    switch (centDetector) {
      case 0: // V0M
        cent = collision.centFV0A();
        break;
      case 1: // T0M
        cent = collision.centFT0M();
        break;
      case 2: // T0A
        cent = collision.centFT0A();
        break;
      case 3: // T0C
        cent = collision.centFT0C();
        break;
      default:
        break;
    }

    return cent;
  }


  /// Compute the scalar product
  /// \param collision is the collision with the Q vector information and event plane
  /// \param candidates are the selected candidates
  template <int decayChannel, typename T1>
  void runFlowAnalysis(CollsWithQvecs::iterator const& collision,
                       T1 const& candidates)
  {
    float QvecX = -999.;
    float QvecY = -999.;
    float amplQvec = -999.;
    float evtPl = -999.;
    switch (qvecDetector.value)
    {
      case 0: // FT0c
        QvecX = collision.qvecFT0CRe();
        QvecY = collision.qvecFT0CIm();
        break;
      case 1: // FT0a
        QvecX = collision.qvecFT0ARe();
        QvecY = collision.qvecFT0AIm();
        break;
      case 2: // FT0m
        QvecX = collision.qvecFT0MRe();
        QvecY = collision.qvecFT0MIm();
        break;
      case 3: // FV0a
        QvecX = collision.qvecFV0ARe();
        QvecY = collision.qvecFV0AIm();
        break;
      case 4: // TPC Positive
        QvecX = collision.qvecBPosRe();
        QvecY = collision.qvecBPosIm();
        amplQvec = collision.nTrkBPos();
        break;
      case 5: // TPC Negative
        QvecX = collision.qvecBNegRe();
        QvecY = collision.qvecBNegIm();
        amplQvec = collision.nTrkBNeg();
        break;
      default:
        break;
    }
    evtPl = epHelper.GetEventPlane(QvecX, QvecY);

    float cent = GetCent(collision, centDetector.value);

    for (auto const& cand : candidates) {
      float massCand = 0.;
      float outputMl = -999.;

      if constexpr (std::is_same<T1, CandDsData>::value)
      {
        switch (decayChannel) {
        case decayChannel::DsToKKPi:
          massCand = hfHelper.invMassDsToKKPi(cand);
          if (storeMl) outputMl = cand.mlProbDsToPiKK()[classMl.value];
          break;
        case decayChannel::DsToPiKK:
          massCand = hfHelper.invMassDsToPiKK(cand);
          if (storeMl) outputMl = cand.mlProbDsToKKPi()[classMl.value];
          break;
        default:
          break;
        }
      }
      else if constexpr (std::is_same<T1, CandDplusData>::value)
      {
        massCand = hfHelper.invMassDplusToPiKPi(cand);
        if (storeMl) outputMl = cand.mlProbDplusToPiKPi()[classMl.value];
      }
      float ptCand = cand.pt();
      float phiCand = cand.phi();

      // TPC only
      if (qvecDetector.value == 4 || qvecDetector.value == 5) {
        float ampl = amplQvec - 3.;
        std::vector<float> tracksQx = {-999., -999., -999.};
        std::vector<float> tracksQy = {-999., -999., -999.};
        GetQvecDtracks(cand, tracksQx, tracksQy, ampl);
        for (unsigned int itrack = 0; itrack < 3; itrack++)
        {
          QvecX -= tracksQx[itrack];
          QvecY -= tracksQy[itrack];
        }
      }
      
      float cosNPhi = TMath::Cos(harmonic * phiCand);
      float sinNPhi = TMath::Sin(harmonic * phiCand);
      float scalprodCand = cosNPhi * QvecX + sinNPhi * QvecY;
      float cosDeltaPhi = TMath::Cos(harmonic * (phiCand - evtPl));

      fillThn(massCand, ptCand, cent, cosNPhi, cosDeltaPhi, scalprodCand, outputMl);
    }
  }


  // Ds
  void processDs(CollsWithQvecs::iterator const& collision,
                 CandDsData const& candidatesDs)
  {
    runFlowAnalysis<decayChannel::DsToKKPi, Partition<CandDsData>>(collision, selectedDsToKKPi);
    runFlowAnalysis<decayChannel::DsToPiKK, Partition<CandDsData>>(collision, selectedDsToPiKK);
  }
  PROCESS_SWITCH(taskFlowCharmHadrons, processDs, "Process Ds candidates", false);


  // Dplus
  void processDplus(CollsWithQvecs::iterator const& collision,
                    CandDplusData const& candidatesDplus)
  {
    runFlowAnalysis<decayChannel::DplusToPiKPi, CandDplusData>(collision, candidatesDplus);
  }
  PROCESS_SWITCH(taskFlowCharmHadrons, processDplus, "Process Dplus candidates", true);


  // Event plane
  void processEventPlaneReso(CollsWithQvecs::iterator const& collision)
  {   
      if (!saveEpResoHisto) {
        LOG(fatal) << "Event plane resolution histogram not saved. Please set saveEpResoHisto to true.";
        return;
      }

      float evCent = GetCent(collision, centDetector.value);
      float epFT0a = epHelper.GetEventPlane(collision.qvecFT0ARe(), collision.qvecFT0AIm());
      float epFT0c = epHelper.GetEventPlane(collision.qvecFT0CRe(), collision.qvecFT0CIm());
      float epFT0m = epHelper.GetEventPlane(collision.qvecFT0MRe(), collision.qvecFT0MIm());
      float epFV0a = epHelper.GetEventPlane(collision.qvecFV0ARe(), collision.qvecFV0AIm());
      float epBPoss = epHelper.GetEventPlane(collision.qvecBPosRe(), collision.qvecBPosIm());
      float epBNegs = epHelper.GetEventPlane(collision.qvecBNegRe(), collision.qvecBNegIm());
    
      registry.fill(HIST("epReso/hEpResoFT0cFT0a"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFT0c, epFT0a)));
      registry.fill(HIST("epReso/hEpResoFT0cFT0m"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFT0c, epFT0m)));
      registry.fill(HIST("epReso/hEpResoFT0cFV0m"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFT0c, epFV0a)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCpos"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFT0c, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCneg"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFT0c, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0aFT0m"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFT0a, epFT0m)));
      registry.fill(HIST("epReso/hEpResoFT0aFV0m"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFT0a, epFV0a)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCpos"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFT0a, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCneg"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFT0a, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0mFV0m"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFT0m, epFV0a)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCpos"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFT0m, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCneg"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFT0m, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFV0mTPCpos"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFV0a, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFV0mTPCneg"), evCent, TMath::Cos(harmonic * GetDeltaPsiInRange(epFV0a, epBNegs)));
  }
  PROCESS_SWITCH(taskFlowCharmHadrons, processEventPlaneReso, "Process event plane resolution", false);

}; // End struct taskFlowCharmHadrons

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<taskFlowCharmHadrons>(cfgc)};
}