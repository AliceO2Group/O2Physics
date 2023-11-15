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
#include "Common/DataModel/EvtPlanes.h"

#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum decayChannel {DplusToPiKPi = 0,
                   DsToKKPi,
                   DsToPiKK};

struct taskFlowCharmHadrons{
  Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<std::string> detector{"detector", "FT0A", "Detector name"}; // Needed only to exclude tracks in TPC case
  Configurable<float> centMin{"centMin", 30., "Minimum centrality"};
  Configurable<float> centMax{"centMax", 50., "Maximum centrality"};
  Configurable<int> selectionFlag{"selectionFlag", 7, "Selection Flag for D"};
  Configurable<std::string> charmHadron{"charmHadron", "Ds", "Charm Hadron"};
  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {100, 1.78, 2.05}, "invariant mass (GeV/#it{c}^{2})"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0., 10.}, "#it{p}_{T} (GeV/#it{c})"};

  AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
  AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
  AxisSpec thnAxisCent{100, 0., 100., "Centrality"};
  AxisSpec thnAxisScalarProd{100, 0., 10., "Scalar Product"};
  AxisSpec thnAxisCosNPhi{100, -1., 1., "cos(n#varphi)"};
  AxisSpec thnAxisSinNPhi{100, -1., 1., "sin(n#varphi)"};
  AxisSpec thnAxisSinNDeltaPhi{100, -1., 1., "sin(n#delta#varphi)"};

  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using CollsWithQvecs = soa::Filtered<soa::Join<aod::Qvectors, aod::Collisions, aod::EvtPlanes>>;

  HfHelper hfHelper;
  
  Filter filterSelectDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;
  Filter filterCentrality = aod::qvec::cent >= centMin && aod::qvec::cent < centMax;

  Partition<CandDsData> selectedDsToKKPi = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag;
  Partition<CandDsData> selectedDsToPiKK = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;

  HistogramRegistry registry{
    "registry",
    {{"hSparseFlowD", "THn for charm hadron flow", {HistType::kTHnSparseF, {{thnAxisInvMass}, {thnAxisPt}, {thnAxisCent}, {thnAxisCosNPhi},{thnAxisSinNPhi}, {thnAxisScalarProd}, {thnAxisSinNDeltaPhi}, {thnAxisScalarProd}}}}}};


  /// Compute the Q vector for the Ds candidate's tracks
  /// \param cand is the candidate
  /// \param tracksQx is the X component of the Q vector for the tracks
  /// \param tracksQy is the Y component of the Q vector for the tracks
  template <typename T1>
  void GetQvecDtracks(const T1& cand,
                      std::vector<float>& tracksQx,
                      std::vector<float>& tracksQy,
                      float& ampl) {
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
  
    tracksQx.push_back(TMath::Cos(harmonic*phiTrack0)*pTtrack0/ampl);
    tracksQy.push_back(TMath::Sin(harmonic*phiTrack0)*pTtrack0/ampl);
    tracksQx.push_back(TMath::Cos(harmonic*phiTrack1)*pTtrack1/ampl);
    tracksQy.push_back(TMath::Sin(harmonic*phiTrack1)*pTtrack1/ampl);
    tracksQx.push_back(TMath::Cos(harmonic*phiTrack2)*pTtrack2/ampl);
    tracksQy.push_back(TMath::Sin(harmonic*phiTrack2)*pTtrack2/ampl);
  }


  /// Fill THnSparse
  /// \param pt is the transverse momentum of the candidate
  /// \param mass is the invariant mass of the candidate
  /// \param cent is the centrality of the collision
  /// \param cosNPhi is the cosine of the n*phi angle
  /// \param sinNPhi is the sine of the n*phi angle
  /// \param cosDeltaPhi is the cosine of the n*(phi - evtPl) angle
  /// \param sinDeltaPhi is the sine of the n*(phi - evtPl) angle
  /// \param sp is the scalar product
  void fillThn(float pt,
               float mass,
               float cent,
               float cosNPhi,
               float sinNPhi,
               float cosDeltaPhi,
               float sinDeltaPhi,
               float sp
               )
  {
    registry.fill(HIST("hSparseFlowD"), mass, pt, cent, cosNPhi, sinNPhi, cosDeltaPhi, sinDeltaPhi, sp);
  }


  /// Compute the scalar product
  /// \param collision is the collision with the Q vector information and event plane
  /// \param candidates are the selected candidates
  template <int decayChannel, typename T1>
  void runFlowAnalysis(CollsWithQvecs::iterator const& collision,
                  const T1& candidates)
  {
    float QvecX = collision.qvecFinalRe();
    float QvecY = collision.qvecFinalIm();
    float amplQvec = collision.ampl();
    float cent = 0.; //collision.cent();  TODO: avoid ambiguity with centrality
    float evtPl = collision.evtPlFinal();

    for (auto const& cand : candidates)
    {
      float ptCand = cand.pt();
      float phiCand = cand.phi();
      float massCand = 0.;

      switch (decayChannel)
      {
        case decayChannel::DsToKKPi:
          massCand = hfHelper.invMassDsToKKPi(cand);
          break;
        case decayChannel::DsToPiKK:
          massCand = hfHelper.invMassDsToPiKK(cand);
          break;
        case decayChannel::DplusToPiKPi:
          massCand = hfHelper.invMassDplusToPiKPi(cand);
          break;
        default:
          break;
      }

      //TPC only
      if (detector.value == "TPC")
      {
        float ampl = amplQvec - 3.;
        std::vector<float> tracksQx, tracksQy;
        GetQvecDtracks(cand, tracksQx, tracksQy, ampl);
        for (unsigned int itrack = 0; itrack < 3; itrack++)
        {
          QvecX -= tracksQx[itrack];
          QvecY -= tracksQy[itrack];
        }
      }
      
      float cosNPhi = TMath::Cos(harmonic*phiCand);
      float sinNPhi = TMath::Sin(harmonic*phiCand);
      float scalprodCand = cosNPhi*QvecX + sinNPhi*QvecY;
      float cosDeltaPhi = TMath::Cos(harmonic*(phiCand - evtPl));
      float sinDeltaPhi = TMath::Sin(harmonic*(phiCand - evtPl));

      fillThn(massCand, ptCand, cent, cosNPhi, sinNPhi, cosDeltaPhi, sinDeltaPhi, scalprodCand);
    }
  }


  // Ds
  void processDs(CollsWithQvecs::iterator const& collision,
                 CandDsData const& candidatesDs)
  {
    runFlowAnalysis<decayChannel::DsToKKPi, Partition<CandDsData>>(collision, selectedDsToKKPi);
    runFlowAnalysis<decayChannel::DsToPiKK, Partition<CandDsData>>(collision, selectedDsToPiKK);

  }
  PROCESS_SWITCH(taskFlowCharmHadrons, processDs, "Process Ds candidates", true);


  // Dplus
    void processDplus(CollsWithQvecs::iterator const& collision,
                      CandDplusData const& candidatesDplus)
  {
    runFlowAnalysis<decayChannel::DplusToPiKPi, CandDplusData>(collision, candidatesDplus);
  }
  PROCESS_SWITCH(taskFlowCharmHadrons, processDplus, "Process Dplus candidates", false);

}; // End struct taskFlowCharmHadrons

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<taskFlowCharmHadrons>(cfgc)};
}