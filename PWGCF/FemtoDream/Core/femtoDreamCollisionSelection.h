// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoDreamCollisionSelection.h
/// \brief FemtoDreamCollisionSelection - event selection within the o2femtodream framework
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMCOLLISIONSELECTION_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMCOLLISIONSELECTION_H_

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"

#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"

#include "TMath.h"

#include <iostream>
#include <string>
#include <vector>

using namespace o2::framework;

namespace o2::analysis::femtoDream
{

/// \class FemtoDreamCollisionSelection
/// \brief Small selection class to check whether a given collision fulfills the specified selections
class FemtoDreamCollisionSelection
{
 public:
  /// Destructor
  virtual ~FemtoDreamCollisionSelection() = default;

  /// Pass the selection criteria to the class
  /// \param zvtxMax Maximal value of the z-vertex
  /// \param checkTrigger whether or not to check for the trigger alias
  /// \param trig Requested trigger alias
  /// \param checkOffline whether or not to check for offline selection criteria
  void setCuts(float zvtxMax, bool checkTrigger, int trig, bool checkOffline, bool addCheckOffline, bool checkRun3, float minSphericity, float sphericityPtmin)
  {
    mCutsSet = true;
    mZvtxMax = zvtxMax;
    mCheckTrigger = checkTrigger;
    mTrigger = static_cast<triggerAliases>(trig);
    mCheckOffline = checkOffline;
    mAddCheckOffline = addCheckOffline;
    mCheckIsRun3 = checkRun3;
    mMinSphericity = minSphericity;
    mSphericityPtmin = sphericityPtmin;
  }

  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  void init(HistogramRegistry* registry)
  {
    if (!mCutsSet) {
      LOGF(error, "Event selection not set - quitting!");
    }
    mHistogramRegistry = registry;
    mHistogramRegistry->add("Event/Zvtx", "; vtx_{z} (cm); Entries", kTH1F, {{300, -12.5, 12.5}});
    mHistogramRegistry->add("Event/MultPercentile", "; Multiplicity Percentile; Entries", kTH1F, {{100, 0, 100}});
    mHistogramRegistry->add("Event/MultPercentileVSMultNTracksPV", "; Multiplicity Percentile; MultNTracks", kTH2F, {{100, 0, 100}, {200, 0, 200}});
    mHistogramRegistry->add("Event/MultNTracksPV", "; MultNTracksPV; Entries", kTH1F, {{200, 0, 200}});
    mHistogramRegistry->add("Event/MultNTracklets", "; MultNTrackslets; Entries", kTH1F, {{300, 0, 300}});
    mHistogramRegistry->add("Event/MultTPC", "; MultTPC; Entries", kTH1F, {{600, 0, 600}});
    mHistogramRegistry->add("Event/Sphericity", "; Sphericity; Entries", kTH1F, {{100, 0, 1}});
  }

  /// Print some debug information
  void printCuts()
  {
    LOG(info) << "Debug information for FemtoDreamCollisionSelection";
    LOG(info) << "Max. z-vertex: " << mZvtxMax;
    LOG(info) << "Check trigger: " << mCheckTrigger;
    LOG(info) << "Trigger: " << mTrigger;
    LOG(info) << " Check offline: " << mCheckOffline;
    LOG(info) << "Min sphericity: " << mMinSphericity;
    LOG(info) << "Min Pt (sphericity): " << mSphericityPtmin;
  }

  /// Check whether the collisions fulfills the specified selections
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return whether or not the collisions fulfills the specified selections
  template <typename C>
  bool isSelectedCollision(C const& col)
  {

    if (std::abs(col.posZ()) > mZvtxMax) {
      return false;
    }
    if (mCheckIsRun3) {
      if (mCheckOffline && !col.sel8()) {
        return false;
      }
      // all checks additional to sel8 are pending to be included in sel8, check them explicitly now and remove them once they have been added to sel8
      // kIsGoodZvtxFT0vsPV can be a dangerous cut because the default event selection value is rather tight
      // Remeber to open the cut (~4cm) with custom event selection task on hyperloop
      if (mAddCheckOffline && (!col.selection_bit(aod::evsel::kNoTimeFrameBorder) || !col.selection_bit(o2::aod::evsel::kNoITSROFrameBorder) || !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup) || !col.selection_bit(o2::aod::evsel::kIsVertexITSTPC) || !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
        return false;
      }
    } else {
      if (mCheckTrigger && !col.alias_bit(mTrigger)) {
        return false;
      }
      if (mCheckOffline && !col.sel7()) {
        return false;
      }
    }
    return true;
  }

  template <typename C, typename T, typename TC>
  bool isEmptyCollision(C const& /*col*/, T const& tracks, TC& trackCuts)
  {
    // check if there is no selected track
    for (auto const& track : tracks) {
      if (trackCuts.isSelectedMinimal(track)) {
        return false;
      }
    }
    return true;
  }

  template <typename C, typename V, typename VC, typename T>
  bool isEmptyCollision(C const& col, V const& V0s, VC& V0Cuts, T const& /*Tracks*/)
  {
    // check if there is no selected V0
    for (auto const& V0 : V0s) {
      auto postrack = V0.template posTrack_as<T>();
      auto negtrack = V0.template negTrack_as<T>();
      if (V0Cuts.isSelectedMinimal(col, V0, postrack, negtrack)) {
        return false;
      }
    }
    return true;
  }

  template <typename C, typename Casc, typename CascC, typename T>
  bool isCollisionWithoutTrkCasc(C const& col, Casc const& Cascades, CascC& CascadeCuts, T const& /*Tracks*/)
  {
    // check if there is no selected Cascade
    for (auto const& Cascade : Cascades) {
      auto postrack = Cascade.template posTrack_as<T>();
      auto negtrack = Cascade.template negTrack_as<T>();
      auto bachtrack = Cascade.template bachelor_as<T>();
      if (CascadeCuts.isSelectedMinimal(col, Cascade, postrack, negtrack, bachtrack)) {
        return false;
      }
    }
    return true;
  }

  /// Pile-up selection of PbPb collisions
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return whether or not the collisions fulfills the specified selections
  template <typename C>
  bool isPileUpCollisionPbPb(C const& col,
                             bool noSameBunchPileup, bool isGoodITSLayersAll)
  {
    if ((noSameBunchPileup && !col.selection_bit(aod::evsel::kNoSameBunchPileup)) || (isGoodITSLayersAll && !col.selection_bit(aod::evsel::kIsGoodITSLayersAll))) {
      return false;
    }

    return true;
  }

  /// occupancy selection
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return whether or not the collisions fulfills the specified selections
  template <typename C>
  bool occupancySelection(C const& col,
                          int tpcOccupancyMin, int tpcOccupancyMax)
  {
    const auto occupancy = col.trackOccupancyInTimeRange();
    if ((occupancy < tpcOccupancyMin || occupancy > tpcOccupancyMax)) {
      return false;
    }
    return true;
  }

  /// Some basic QA of the event
  /// \tparam T type of the collision
  /// \param col Collision
  template <typename T>
  void fillQA(T const& col, float cent)
  {
    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/Zvtx"), col.posZ());
      mHistogramRegistry->fill(HIST("Event/MultNTracksPV"), col.multNTracksPV());
      mHistogramRegistry->fill(HIST("Event/MultTPC"), col.multTPC());
      if (mCheckIsRun3) {
        mHistogramRegistry->fill(HIST("Event/MultPercentile"), cent);
        mHistogramRegistry->fill(HIST("Event/MultPercentileVSMultNTracksPV"), cent, col.multNTracksPV());
      } else {
        mHistogramRegistry->fill(HIST("Event/MultNTracklets"), col.multTracklets());
      }
    }
  }

  /// Initializes histograms for qn bin
  /// \param registry Histogram registry to be passed
  void initEPQA(HistogramRegistry* registry)
  {
    mHistogramQn = registry;
    mHistogramQn->add("Event/centFT0CBeforeQn", "; cent", kTH1F, {{10, 0, 100}});
    mHistogramQn->add("Event/centFT0CAfterQn", "; cent", kTH1F, {{10, 0, 100}});
    mHistogramQn->add("Event/centVsqn", "; cent; qn", kTH2F, {{10, 0, 100}, {1000, 0, 1000}});
    mHistogramQn->add("Event/centVsqnVsSpher", "; cent; qn; Sphericity", kTH3F, {{10, 0, 100}, {100, 0, 1000}, {100, 0, 1}});
    mHistogramQn->add("Event/qnBin", "; qnBin; entries", kTH1F, {{20, 0, 20}});
    mHistogramQn->add("Event/psiEP", "; #Psi_{EP} (deg); entries", kTH1F, {{100, 0, 180}});
    mHistogramQn->add("Event/epReso_FT0CTPC", "; cent; qnBin; reso_ft0c_tpc", kTH2F, {{10, 0, 100},{10,0,10}});
    mHistogramQn->add("Event/epReso_FT0ATPC", "; cent; qnBin; reso_ft0a_tpc", kTH2F, {{10, 0, 100},{10,0,10}});
    mHistogramQn->add("Event/epReso_FT0CFT0A", "; cent; qnBin; reso_ft0c_ft0a", kTH2F, {{10, 0, 100},{10,0,10}});
    mHistogramQn->add("Event/epReso_count", "; cent; qnBin; count", kTH2F, {{10, 0, 100},{10,0,10}});

    return;
  }

  /// Initializes histograms for the flow calculation
  /// \param registry Histogram registry to be passed
  void initFlow(HistogramRegistry* registry, bool doQnSeparation, int mumQnBins = 10, int centBins = 10)
  {
    if (!mCutsSet) {
      LOGF(error, "Event selection not set - quitting!");
    }

    mHistogramQn = registry;
    mHistogramQn->add("Event/hN2allQn", ";centrality; #sum Re(Q_{2,A} Q_{2,B}^{*})", kTH1F, {{centBins, 0, 100}});
    mHistogramQn->add("Event/hD2allQn", ";centrality; #sum (W_{A} W_{B})", kTH1F, {{centBins, 0, 100}});
    mHistogramQn->get<TH1>(HIST("Event/hN2allQn"))->Sumw2();
    mHistogramQn->get<TH1>(HIST("Event/hD2allQn"))->Sumw2();

    if (doQnSeparation) {
      for (int iqn(0); iqn < mumQnBins; ++iqn) {
        hN2.push_back(mHistogramQn->add(("Qn/hN2_" + std::to_string(iqn)).c_str(), ";centrality; #sum Re(Q_{2,A} Q_{2,B}^{*})", kTH1F, {{centBins, 0, 100}}));
        hD2.push_back(mHistogramQn->add(("Qn/hD2_" + std::to_string(iqn)).c_str(), ";centrality; #sum (W_{A} W_{B})", kTH1F, {{centBins, 0, 100}}));
      }
      for (int iqn(0); iqn < mumQnBins; ++iqn) {
        std::get<std::shared_ptr<TH1>>(hN2[iqn])->Sumw2();
        std::get<std::shared_ptr<TH1>>(hD2[iqn])->Sumw2();
      }
    }
    return;
  }

  /// \todo to be implemented!
  /// Compute the sphericity of an event
  /// Important here is that the filter on tracks does not interfere here!
  /// In Run 2 we used here global tracks within |eta| < 0.8
  /// \tparam T1 type of the collision
  /// \tparam T2 type of the tracks
  /// \param col Collision
  /// \param tracks All tracks
  /// \return value of the sphericity of the event
  template <typename T1, typename T2>
  float computeSphericity(T1 const& col, T2 const& tracks)
  {
    double ptTot = 0.;
    double s00 = 0.; // elements of the sphericity matrix taken form EPJC72:2124
    double s01 = 0.;
    // double s10 = 0.;
    double s11 = 0.;

    int numOfTracks = col.numContrib();
    if (numOfTracks < 3)
      return -9999.;

    for (auto const& track : tracks) {
      double pt = track.pt();
      double eta = track.eta();
      double px = track.px();
      double py = track.py();
      if (TMath::Abs(pt) < mSphericityPtmin || TMath::Abs(eta) > 0.8) {
        continue;
      }

      ptTot += pt;

      s00 += px * px / pt;
      s01 += px * py / pt;
      // s10 = s01;
      s11 += py * py / pt;
    }

    // normalize to total Pt to obtain a linear form:
    if (ptTot == 0.)
      return -9999.;
    s00 /= ptTot;
    s11 /= ptTot;
    s01 /= ptTot;

    // Calculate the trace of the sphericity matrix:
    double T = s00 + s11;
    // Calculate the determinant of the sphericity matrix:
    double D = s00 * s11 - s01 * s01; // S10 = S01

    // Calculate the eigenvalues of the sphericity matrix:
    double lambda1 = 0.5 * (T + std::sqrt(T * T - 4. * D));
    double lambda2 = 0.5 * (T - std::sqrt(T * T - 4. * D));

    if ((lambda1 + lambda2) == 0.)
      return -9999.;

    double spt = -1.;

    if (lambda2 > lambda1) {
      spt = 2. * lambda1 / (lambda1 + lambda2);
    } else {
      spt = 2. * lambda2 / (lambda1 + lambda2);
    }

    mHistogramRegistry->fill(HIST("Event/Sphericity"), spt);

    return spt;
  }

  /// Compute the qn-vector(FT0C) of an event
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return value of the qn-vector of FT0C of the event
  template <typename T>
  float computeqnVec(T const& col, int qvecMod = 0)
  {
    double qn = -999.f;
    if (qvecMod == 0){
      qn = std::sqrt(col.qvecFT0CReVec()[0] * col.qvecFT0CReVec()[0] + col.qvecFT0CImVec()[0] * col.qvecFT0CImVec()[0]) * std::sqrt(col.sumAmplFT0C());
    } else if (qvecMod == 1){
      qn = std::sqrt(col.qvecFT0AReVec()[0] * col.qvecFT0AReVec()[0] + col.qvecFT0AImVec()[0] * col.qvecFT0AImVec()[0]) * std::sqrt(col.sumAmplFT0A());
    } else {
       LOGP(error, "no selected detector of Qvec for ESE ");  
       return qn;   
    }
    return qn;
  }

  /// Compute the event plane of an event
  /// \tparam T type of the collision
  /// \param col Collision
  /// \param nmode EP in which harmonic(default 2nd harmonic)
  /// \return angle of the event plane (rad) of FT0C of the event
  template <typename T>
  float computeEP(T const& col, int nmode, int qvecMod)
  {
    double EP = -999.f;
    if (qvecMod == 0){
      EP = ((1. / nmode) * (TMath::ATan2(col.qvecFT0CImVec()[0], col.qvecFT0CReVec()[0])));
    } else if (qvecMod == 1) {
      EP = ((1. / nmode) * (TMath::ATan2(col.qvecFT0AImVec()[0], col.qvecFT0AReVec()[0])));      
    } else if (qvecMod == 2) {
      EP = ((1. / nmode) * (TMath::ATan2(col.qvecTPCallImVec()[0], col.qvecTPCallReVec()[0])));            
    } else {
       LOGP(error, "no selected detector of Qvec for EP");  
       return EP;   
    }
    
    if (EP < 0){
      EP += TMath::Pi();
    } // atan2 return in rad -pi/2-pi/2, then make it 0-pi   
    return EP;
  }

  /// Compute the event plane resolution of 3 sub-events
  /// \tparam T type of the collision
  /// \param col Collision
  /// \param nmode EP in which harmonic(default 2nd harmonic)
  template <typename T>
  void fillEPReso(T const& col, int nmode, float centrality)
  {
    const float psi_ft0c = ((1. / nmode) * (TMath::ATan2(col.qvecFT0CImVec()[0], col.qvecFT0CReVec()[0])));
    const float psi_ft0a = ((1. / nmode) * (TMath::ATan2(col.qvecFT0AImVec()[0], col.qvecFT0AReVec()[0])));
    const float psi_tpc = ((1. / nmode) * (TMath::ATan2(col.qvecTPCallImVec()[0], col.qvecTPCallReVec()[0]))); 

    mHistogramQn->fill(HIST("Event/epReso_FT0CTPC"), centrality, mQnBin + 0.f, std::cos((psi_ft0c - psi_tpc) * nmode));
    mHistogramQn->fill(HIST("Event/epReso_FT0ATPC"), centrality, mQnBin + 0.f, std::cos((psi_ft0a - psi_tpc) * nmode));
    mHistogramQn->fill(HIST("Event/epReso_FT0CFT0A"), centrality, mQnBin + 0.f, std::cos((psi_ft0c - psi_ft0a) * nmode));
    mHistogramQn->fill(HIST("Event/epReso_count"), centrality, mQnBin + 0.f);
  }

  /// \return the 1-d qn-vector separator to 2-d
  std::vector<std::vector<float>> getQnBinSeparator2D(std::vector<float> flat, const int numQnBins = 10)
  {
    size_t nBins = numQnBins + 1;

    if (flat.empty() || flat.size() % nBins != 0) {
      LOGP(error, "ConfQnBinSeparator size = {} is not divisible by {}",
           flat.size(), nBins);
      return {{-999, -999}};
    }

    size_t nCent = flat.size() / nBins;
    std::vector<std::vector<float>> res(nCent, std::vector<float>(nBins));

    for (size_t i = 0; i < nCent; ++i) {
      for (size_t j = 0; j < nBins; ++j) {
        res[i][j] = flat[i * nBins + j];
      }
    }
    return res;
  }

  /// Get the bin number of qn-vector(FT0C) of an event
  /// \param centBinWidth centrality bin width, example: per 1%, per 10% ...
  /// \return bin number of qn-vector of the event
  // add a param : bool doFillHisto ?
  int myqnBin(float centrality, float centMax, bool doFillCent, std::vector<float> qnBinSeparator, float qn, const int numQnBins, float centBinWidth = 1.f)
  {
    auto twoDSeparator = getQnBinSeparator2D(qnBinSeparator, numQnBins);
    if (twoDSeparator.empty() || twoDSeparator[0][0] == -999.) {
      LOGP(warning, "ConfQnBinSeparator not set, using default fallback!");
      return -999; // safe fallback
    }

    // if (doFillHisto)
    //   mHistogramQn->fill(HIST("Event/centFT0CBefore"), centrality);
    // add a param : bool doFillHisto ?
    int qnBin = -999;
    int mycentBin = static_cast<int>(centrality / centBinWidth);
    if (mycentBin >= static_cast<int>(centMax / centBinWidth))
      return qnBin;

    if (mycentBin > static_cast<int>(twoDSeparator.size()) - 1)
      return qnBin;

    if (doFillCent)
      mHistogramQn->fill(HIST("Event/centFT0CAfterQn"), centrality);

    for (int iqn(0); iqn < static_cast<int>(twoDSeparator[mycentBin].size()) - 1; ++iqn) {
      if (qn > twoDSeparator[mycentBin][iqn] && qn <= twoDSeparator[mycentBin][iqn + 1]) {
        qnBin = iqn;
        break;
      } else {
        continue;
      }
    }

    mQnBin = qnBin;
    return qnBin;
  }

  /// \fill event-wise informations
  template <typename T>
  void fillEPQA(T& col, float centrality, float fSpher, float qn, float psiEP, int nmode = 2)
  {
    mHistogramQn->fill(HIST("Event/centFT0CBeforeQn"), centrality);
    mHistogramQn->fill(HIST("Event/centVsqn"), centrality, qn);
    mHistogramQn->fill(HIST("Event/centVsqnVsSpher"), centrality, qn, fSpher);
    mHistogramQn->fill(HIST("Event/qnBin"), mQnBin + 0.f);
    mHistogramQn->fill(HIST("Event/psiEP"), psiEP);
    fillEPReso(col, nmode, centrality);
  }

  /// \todo to be implemented!
  /// Do cumulants for flow calculation
  /// \tparam T type of the collision
  /// \param doQnSeparation to fill flow in divied qn bins
  /// \param qnBin should be <int> in 0-9
  /// \param fEtaGap eta gap for flow cumulant
  template <typename T1, typename T2, typename TC>
  void doCumulants(T1 const& col, T2 const& tracks, TC& trackCuts, float centrality, bool doQnSeparation = false, int numQnBins = 10, float fEtaGap = 0.5f, float ptMin = 0.2f, float ptMax = 5.0f, float harmonic = 2.0f)
  {
    int numOfTracks = col.numContrib();
    if (numOfTracks < 3)
      return;

    double Q2A_re = 0., Q2A_im = 0., WA = 0.;
    double Q2B_re = 0., Q2B_im = 0., WB = 0.;

    int nA = 0, nB = 0;

    for (auto const& trk : tracks) {
      if (!trackCuts.isSelectedMinimal(trk)) {
        continue;
      }
      const double pt = trk.pt();
      const double eta = trk.eta();
      if (pt < ptMin || pt > ptMax) {
        continue;
      }

      const double w = 1.0; // TODO: NUA/NUE weight
      const double phi = trk.phi();
      const double c = w * TMath::Cos(harmonic * phi);
      const double s = w * TMath::Sin(harmonic * phi);

      if (eta > fEtaGap) {
        Q2A_re += c;
        Q2A_im += s;
        WA += w;
        nA++;
      } else if (eta < -1 * fEtaGap) {
        Q2B_re += c;
        Q2B_im += s;
        WB += w;
        nB++;
      }
    }

    // need at least 1 track on each side to form pairs; for stability, require >=2
    if (nA < 2 || nB < 2) {
      return;
    }
    const double D2_evt = WA * WB;
    if (D2_evt <= 0.) {
      return;
    }

    // N2_evt = Re(Q2A * conj(Q2B)) = Q2A_re*Q2B_re + Q2A_im*Q2B_im
    const double N2_evt = Q2A_re * Q2B_re + Q2A_im * Q2B_im;

    mHistogramQn->fill(HIST("Event/hN2allQn"), centrality, N2_evt);
    mHistogramQn->fill(HIST("Event/hD2allQn"), centrality, D2_evt);
    if (doQnSeparation && mQnBin >= 0 && mQnBin < numQnBins) {
      std::get<std::shared_ptr<TH1>>(hN2[mQnBin])->Fill(centrality, N2_evt);
      std::get<std::shared_ptr<TH1>>(hD2[mQnBin])->Fill(centrality, D2_evt);
    }
    return;
  }

 private:
  HistogramRegistry* mHistogramRegistry = nullptr; ///< For QA output
  bool mCutsSet = false;                           ///< Protection against running without cuts
  bool mCheckTrigger = false;                      ///< Check for trigger
  bool mCheckOffline = false;                      ///< Check for offline criteria (might change)
  bool mAddCheckOffline = false;                   ///< Additional check for offline criteria (added to sel8 soon)
  bool mCheckIsRun3 = false;                       ///< Check if running on Pilot Beam
  triggerAliases mTrigger = kINT7;                 ///< Trigger to check for
  float mZvtxMax = 999.f;                          ///< Maximal deviation from nominal z-vertex (cm)
  float mMinSphericity = 0.f;
  float mSphericityPtmin = 0.f;
  int mQnBin = -999;
  HistogramRegistry* mHistogramQn = nullptr; ///< For flow cumulant output
  std::vector<HistPtr> hN2;                  ///< Histograms of c22 per Qn bin
  std::vector<HistPtr> hD2;                  ///< Histograms of c22 per Qn bin
};
} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMCOLLISIONSELECTION_H_
