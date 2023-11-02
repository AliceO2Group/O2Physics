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

/// \file scalarProduct.cxx
/// \brief Analysis task for Scalar Product method
///
/// \author F. Chinu, Università di Torino, Italy
/// \author S. Politanò, INFN & Politecnico Torino, Italy
/// \author S. Trgolo, INFN Torino, Italy

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/DataModel/Qvectors.h"

#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
  namespace full
{
DECLARE_SOA_COLUMN(Qx, qx, float);       //! X component of the Q vector
DECLARE_SOA_COLUMN(Qy, qy, float);       //! Y component of the Q vector
DECLARE_SOA_COLUMN(Pt, pt, float);       //! Transverse momentum of prong1 (GeV/c)
DECLARE_SOA_COLUMN(M,   m, float);       //! Mass of the candidate (GeV/c^2)
DECLARE_SOA_COLUMN(Cent, cent, float);   //! Centrality of the collision
DECLARE_SOA_COLUMN(SP, sp, float);       //! Scalar product
} // namespace full
DECLARE_SOA_TABLE(HfScalarProduct, "AOD", "HFSCALARPRODUCT",
                  full::Qx,
                  full::Qy,
                  full::Pt,
                  full::Cent,
                  full::M,
                  full::SP);
} // namespace o2::aod

struct scalarProduct{
  Produces<o2::aod::HfScalarProduct> hfScalarProduct;

  Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<std::string> detector{"detector", "FT0A", "Detector name"}; // Needed only to exclude tracks in TPC case
  Configurable<float> centMin{"centMin", 30., "Minimum centrality"};
  Configurable<float> centMax{"centMax", 50., "Maximum centrality"};
  Configurable<int> selectionFlagD{"selectionFlagD", 7, "Selection Flag for D"};
  Configurable<std::string> Dmeson{"Dmeson", "Ds", "D meson"};

  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using MyQvecs = soa::Filtered<soa::Join<aod::Qvectors, aod::Collisions>>;
  
  Filter filterSelectDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagD || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagD;
  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagD;
  Filter filterCentrality = aod::qvec::cent >= centMin && aod::qvec::cent < centMax;

  // TODO: add possibility to consider different mass hypotheses
  // Partition<CandDsData> selectedDsToKKPiCand = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagD;
  // Partition<CandDsData> selectedDsToPiKKCand = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagD;


  HistogramRegistry registry{
    "registry",
    {{"hSP", "Sp;sp;entries", {HistType::kTH1F, {{100, -1.1, 1.1}}}},
     {"hSpM", "SpVsM;sp;#it{M};entries", {HistType::kTH2F, {{100, -1.1, 1.1}, {100, 1.8, 2.2}}}},
     {"hSpPt", "SpVsPt;sp;#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, -1.1, 1.1}, {100, 0., 10.}}}},
     {"hSpCent", "SpVsCent;sp;centrality;entries", {HistType::kTH2F, {{100, -1.1, 1.1}, {100, 0., 100.}}}},
     {"hSpPtCent", "SpVsPtCent;sp;#it{p}_{T} (GeV/#it{c});centrality", {HistType::kTH3F, {{100, -1.1, 1.1}, {100, 0., 10.}, {100, 0., 100.}}}}}};


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


  /// Fill histograms
  /// \param sp is the scalar product
  /// \param pt is the transverse momentum of the candidate
  /// \param cent is the centrality of the collision
  /// \param m is the mass of the candidate
  void FillHistograms(float sp,
                      float pt,
                      float cent,
                      float mass)
  {
    registry.fill(HIST("hSP"), sp);
    registry.fill(HIST("hSpM"), sp, mass);
    registry.fill(HIST("hSpPt"), sp, pt);
    registry.fill(HIST("hSpCent"), sp, cent);
    registry.fill(HIST("hSpPtCent"), sp, pt, cent);
  }


  /// Compute the invariant mass of the Ds candidate
  /// \param cand is the candidate
  template <typename T1>
  std::vector<float> invMassDs(const T1& cand) {
    std::vector<float> mCand;
    if (cand.isSelDsToKKPi() >= selectionFlagD) {
      mCand.push_back(invMassDsToKKPi(cand));
    } 
    else if (cand.isSelDsToPiKK() >= selectionFlagD) {
      mCand.push_back(invMassDsToPiKK(cand));
    } 
    return mCand;
  }

  /// Mathcing between candidates and Q vectors using the global index
  /// \param collIDCand is the global index of the candidate
  /// \param qVecs is the Q vectors
  /// \param QvecX is the X component of the Q vector
  /// \param QvecY is the Y component of the Q vector
  /// \param centQvec is the centrality of the Q vector
  bool MatchCollID(int collIDCand,
                   MyQvecs const& qVecs,
                   float& QvecX,
                   float& QvecY,
                   float& centQvec,
                   double& amplQvec)
  {
     for (auto const& qvec : qVecs) 
     {
      if (qvec.globalIndex() == collIDCand)
      {
        //GetQvectors(qvec, QvecX, QvecY);
        QvecX = qvec.qvecFinalRe();
        QvecY = qvec.qvecFinalIm();
        centQvec = qvec.cent();
        amplQvec = qvec.ampl();
        return true;
      }
     }
     return false;
  }
  

  /// Compute the scalar product
  /// \param qVecs is the Q vectors
  /// \param candidates is the candidates
  template <typename T1>
  void ComputeSPD(MyQvecs const& qVecs,
                  const T1& candidates)
  {
    float QvecX, QvecY, centQvec;
    double amplQvec;

    hfScalarProduct.reserve(candidates.size());
    for (auto const& cand : candidates)
    {
      Int_t collIDCand = cand.globalIndex();
      float ptCand = cand.pt();
      float phiCand = cand.phi();
      float mCand = 0.; // Adjust mass calculation

      if (!MatchCollID(collIDCand, qVecs, QvecX, QvecY, centQvec, amplQvec)) continue;

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
      
      float scalprod = TMath::Cos(harmonic*phiCand)*QvecX + TMath::Sin(harmonic*phiCand)*QvecY;
      FillHistograms(scalprod, ptCand, centQvec, mCand);
      hfScalarProduct(QvecX, QvecY,  ptCand, centQvec, mCand, scalprod);
    }
  }


  void processDs(MyQvecs const& qVecs,
                CandDsData const& candidatesDs)
  {
    ComputeSPD<CandDsData>(qVecs, candidatesDs);

  }
  PROCESS_SWITCH(scalarProduct, processDs, "Process Ds candidates", true);

  void processDplus(MyQvecs const& qVecs,
                    CandDplusData const& candidatesDplus)
  {
    ComputeSPD<CandDplusData>(qVecs, candidatesDplus);
  }
  PROCESS_SWITCH(scalarProduct, processDplus, "Process Dplus candidates", false);

}; // End struct ScalarProduct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<scalarProduct>(cfgc)};
}