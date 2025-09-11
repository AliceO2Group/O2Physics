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

#include <string>
#include <iostream>
#include <vector>
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Qvectors.h"

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
                              bool noSameBunchPileup, bool isGoodZvtxFT0vsPV, 
                              bool isGoodITSLayersAll, bool noCollInRofStandard, 
                              bool noHighMultCollInPrevRof, bool noCollInTimeRangeStandard, 
                              bool /*isVertexITSTPC*/,
                              int tpcOccupancyMin, int tpcOccupancyMax)
  {
    const auto occupancy = col.trackOccupancyInTimeRange();
    if ((occupancy < tpcOccupancyMin || occupancy > tpcOccupancyMax)) {
      return false;
    }
    if ((noSameBunchPileup && !col.selection_bit(aod::evsel::kNoSameBunchPileup)) 
      || (isGoodZvtxFT0vsPV && !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) 
      || (isGoodITSLayersAll && !col.selection_bit(aod::evsel::kIsGoodITSLayersAll)) 
      || (noCollInRofStandard && !col.selection_bit(aod::evsel::kNoCollInRofStandard)) 
      || (noHighMultCollInPrevRof && !col.selection_bit(aod::evsel::kNoHighMultCollInPrevRof)) 
      || (noCollInTimeRangeStandard && !col.selection_bit(aod::evsel::kNoCollInTimeRangeStandard))
          // || (isVertexITSTPC && !col.selection_bit(aod::evsel::kIsVertexITSTPC))
    ) {
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

  /// Initializes histograms for the flow calculation
  /// \param registry Histogram registry to be passed
  void initFlow(HistogramRegistry* registry, bool doQnSeparation, int mumQnBins = 10, int binPt = 100, int binEta = 32)
  {
    if (!mCutsSet) {
      LOGF(error, "Event selection not set - quitting!");
    }
    mReQthisEvt = new TH2D("ReQthisEvt", "", binPt, 0., 5., binEta, -0.8, 0.8);
    mImQthisEvt = new TH2D("ImQthisEvt", "", binPt, 0., 5., binEta, -0.8, 0.8);
    mReQ2thisEvt = new TH2D("ReQ2thisEvt", "", binPt, 0., 5., binEta, -0.8, 0.8);
    mImQ2thisEvt = new TH2D("ImQ2thisEvt", "", binPt, 0., 5., binEta, -0.8, 0.8);
    mMQthisEvt = new TH2D("MQthisEvt", "", binPt, 0., 5., binEta, -0.8, 0.8);
    mMQWeightthisEvt = new TH2D("MQWeightthisEvt",  "", binPt, 0., 5., binEta, -0.8, 0.8);
    mHistogramQn = registry;
    mHistogramQn->add("Event/centFT0CBefore", "; cent", kTH1F, {{10, 0, 100}});
    mHistogramQn->add("Event/centFT0CAfter", "; cent", kTH1F, {{10, 0, 100}});
    mHistogramQn->add("Event/centVsqn", "; cent; qn", kTH2F, {{10, 0, 100}, {100, 0, 1000}});
    mHistogramQn->add("Event/centVsqnVsSpher", "; cent; qn; Sphericity", kTH3F, {{10, 0, 100}, {100, 0, 1000}, {100, 0, 1}});
    mHistogramQn->add("Event/qnBin", "; qnBin; entries", kTH1F, {{20, 0, 20}});
    
    mHistogramQn->add<TProfile>("Event/profileC22", "; cent; c22", kTProfile, {{10, 0, 100}});
    mHistogramQn->add<TProfile>("Event/profileC24", "; cent; c24", kTProfile, {{10, 0, 100}});
    if (doQnSeparation){      
      for (int iqn(0); iqn < mumQnBins; ++iqn) {
        mHistogramQn->add<TProfile>(("Qn/profileC22_" + std::to_string(iqn)).c_str(), "; cent; c22", kTProfile, {{10, 0, 100}});
        mHistogramQn->add(("Qn/mult_" + std::to_string(iqn)).c_str(), "; cent; c22", kTH1F, {{100, 0, 3500}});
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

  /// \todo to be implemented!
  /// Compute the qn-vector(FT0C) of an event
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return value of the qn-vector of FT0C of the event
  template <typename T>
  float computeqnVec(T const& col)
  {
    double qn = std::sqrt(col.qvecFT0CReVec()[0] * col.qvecFT0CReVec()[0] + col.qvecFT0CImVec()[0] * col.qvecFT0CImVec()[0]) * std::sqrt(col.sumAmplFT0C());
    return qn;
  }

  /// \todo to be implemented!
  /// \return the 1-d qn-vector separator to 2-d
  std::vector<std::vector<float>> getQnBinSeparator2D(std::vector<float> flat) 
  {
    constexpr size_t nBins = 11;

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

  /// \todo to be implemented!
  /// Get the bin number of qn-vector(FT0C) of an event
  /// \tparam T type of the collision
  /// \param col Collision
  /// \param centBinWidth centrality bin width, example: per 1%, per 10% ...
  /// \return bin number of qn-vector of the event
  template <typename T>
  int myqnBin(T const& col, float centMax, std::vector<float> qnBinSeparator, float fSpher, float fMult, float centBinWidth = 1.f)
  {
    auto twoDSeparator = getQnBinSeparator2D(qnBinSeparator);
    if (twoDSeparator.empty() || twoDSeparator[0][0] == -999.) {
      LOGP(warning, "ConfQnBinSeparator not set, using default fallback!");
      return -999; // safe fallback
    }

    mHistogramQn->fill(HIST("Event/centFT0CBefore"), col.centFT0C());  
    int qnBin = -999;
    float qn = computeqnVec(col);
    int mycentBin = static_cast<int>(col.centFT0C() / centBinWidth);
    if (mycentBin >= static_cast<int>(centMax / centBinWidth))
      return qnBin;

    if (mycentBin > static_cast<int>(twoDSeparator.size()) -1)
      return qnBin;

    for (int iqn(0); iqn < static_cast<int>(twoDSeparator[mycentBin].size()) - 1; ++iqn) {
      if (qn > twoDSeparator[mycentBin][iqn] && qn <= twoDSeparator[mycentBin][iqn + 1]) {
        qnBin = iqn;
        break;
      } else {
        continue;
      }
    }
    
    mHistogramQn->fill(HIST("Event/centFT0CAfter"), col.centFT0C());
    mHistogramQn->fill(HIST("Event/centVsqn"), col.centFT0C(), qn);
    mHistogramQn->fill(HIST("Event/centVsqnVsSpher"), col.centFT0C(), qn, fSpher);
    mHistogramQn->fill(HIST("Event/qnBin"), qnBin);
    if (qnBin >= 0 && qnBin < 10){
      switch (qnBin) {
        case 0:
          mHistogramQn->fill(HIST("Qn/mult_") + HIST("0"), fMult);
          break;
        case 1:
          mHistogramQn->fill(HIST("Qn/mult_") + HIST("1"), fMult);
          break;
        case 2:
          mHistogramQn->fill(HIST("Qn/mult_") + HIST("2"), fMult); 
          break;       
        case 3:
          mHistogramQn->fill(HIST("Qn/mult_") + HIST("3"), fMult);
          break;
        case 4:
          mHistogramQn->fill(HIST("Qn/mult_") + HIST("4"), fMult);
          break;
        case 5:
          mHistogramQn->fill(HIST("Qn/mult_") + HIST("5"), fMult);
          break;
        case 6:
          mHistogramQn->fill(HIST("Qn/mult_") + HIST("6"), fMult);
          break;
        case 7:
          mHistogramQn->fill(HIST("Qn/mult_") + HIST("7"), fMult);
          break;
        case 8:
          mHistogramQn->fill(HIST("Qn/mult_") + HIST("8"), fMult);
          break;
        case 9:
          mHistogramQn->fill(HIST("Qn/mult_") + HIST("9"), fMult);
          break; 
        default:
          return qnBin; // invalid qn bin       
      }
    }
    return qnBin;
  }

  /// \todo to be implemented!
  /// Fill cumulants histo for flow calculation
  /// Reset hists event-by-event
  /// \tparam T1 type of the collision
  /// \tparam T2 type of the tracks
  /// \param tracks All tracks
  template <typename T1, typename T2>
  void fillCumulants(T1 const& col, T2 const& tracks, float fHarmonic=2.f)
  { 
    int numOfTracks = col.numContrib();
    if (numOfTracks < 3)
      return ;

    mReQthisEvt->Reset();
    mImQthisEvt->Reset();
    mReQ2thisEvt->Reset(); 
    mImQ2thisEvt->Reset(); 
    mMQthisEvt->Reset();
    mMQWeightthisEvt->Reset(); 

    for (auto const& track : tracks) {
      double weight=1; // Will implement NUA&NUE correction
      double phi = track.phi();
      double pt = track.pt();
      double eta = track.eta();
      double cosnphi = weight * TMath::Cos(fHarmonic*phi);
      double sinnphi = weight * TMath::Sin(fHarmonic*phi);
      double cos2nphi = weight * TMath::Cos(2*fHarmonic*phi); 
      double sin2nphi = weight * TMath::Sin(2*fHarmonic*phi); 
      mReQthisEvt->Fill(pt, eta, cosnphi);
      mImQthisEvt->Fill(pt, eta, sinnphi);
      mReQ2thisEvt->Fill(pt,eta,cos2nphi); 
      mImQ2thisEvt->Fill(pt, eta, sin2nphi); 
      mMQthisEvt ->Fill(pt, eta);
      mMQWeightthisEvt ->Fill(pt, eta, weight); 
    }
    return;
  }

  /// \todo to be implemented!
  /// Do cumulants for flow calculation
  /// \tparam T type of the collision
  /// \param doQnSeparation to fill flow in divied qn bins
  /// \param qnBin should be <int> in 0-9
  /// \param fEtaGap eta gap for flow cumulant
  template <typename T>
  void doCumulants(T const& col, bool doQnSeparation = false, int qnBin = -999, int mumQnBinNum = 10, float fEtaGap = 0.3f, int binPt = 100, int binEta = 32)
  {
    if (mMQthisEvt->Integral(1, binPt, 1, binEta) < 2) 
      return;
  
    double allReQ  = mReQthisEvt ->Integral(1, binPt, 1, binEta);
    double allImQ  = mImQthisEvt ->Integral(1, binPt, 1, binEta);
    TComplex Q(allReQ, allImQ);
    TComplex QStar = TComplex::Conjugate(Q);
  
    double posEtaRe = mReQthisEvt->Integral(1, binPt, mReQthisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), binEta);
    double posEtaIm = mImQthisEvt->Integral(1, binPt, mImQthisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), binEta);
    if (mMQthisEvt->Integral(1, binPt,  mMQthisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), binEta) < 2) 
      return;
    float posEtaMQ = mMQWeightthisEvt->Integral(1, binPt, mMQthisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), binEta);
    TComplex posEtaQ = TComplex(posEtaRe, posEtaIm);
    TComplex posEtaQStar = TComplex::Conjugate(posEtaQ);
  
    double negEtaRe = mReQthisEvt->Integral(1, binPt, 1, mReQthisEvt->GetYaxis()->FindBin(-1*fEtaGap-1e-6));
    double negEtaIm = mImQthisEvt->Integral(1, binPt, 1, mImQthisEvt->GetYaxis()->FindBin(-1*fEtaGap-1e-6));
    if (mMQthisEvt->Integral(1, binPt, 1,  mMQthisEvt->GetYaxis()->FindBin(-1*fEtaGap-1e-6)) < 2) 
      return;
    float negEtaMQ = mMQWeightthisEvt->Integral(1, binPt, 1, mMQthisEvt->GetYaxis()->FindBin(-1*fEtaGap-1e-6));
    TComplex negEtaQ = TComplex(negEtaRe, negEtaIm);
    TComplex negEtaQStar = TComplex::Conjugate(negEtaQ);

    mHistogramQn->get<TProfile>(HIST("Event/profileC22"))->Fill(col.centFT0C(), (negEtaQ*posEtaQStar).Re()/(negEtaMQ*posEtaMQ), (negEtaMQ*posEtaMQ));
    if (doQnSeparation && qnBin >= 0 && qnBin < mumQnBinNum){
      switch (qnBin) {
        case 0:
          mHistogramQn->get<TProfile>(HIST("Qn/profileC22_") + HIST("0"))->Fill(col.centFT0C(), (negEtaQ*posEtaQStar).Re()/(negEtaMQ*posEtaMQ), (negEtaMQ*posEtaMQ));
          break;
        case 1:
          mHistogramQn->get<TProfile>(HIST("Qn/profileC22_") + HIST("1"))->Fill(col.centFT0C(), (negEtaQ*posEtaQStar).Re()/(negEtaMQ*posEtaMQ), (negEtaMQ*posEtaMQ));
          break;
        case 2:
          mHistogramQn->get<TProfile>(HIST("Qn/profileC22_") + HIST("2"))->Fill(col.centFT0C(), (negEtaQ*posEtaQStar).Re()/(negEtaMQ*posEtaMQ), (negEtaMQ*posEtaMQ));        
          break;
        case 3:
          mHistogramQn->get<TProfile>(HIST("Qn/profileC22_") + HIST("3"))->Fill(col.centFT0C(), (negEtaQ*posEtaQStar).Re()/(negEtaMQ*posEtaMQ), (negEtaMQ*posEtaMQ));
          break;
        case 4:
          mHistogramQn->get<TProfile>(HIST("Qn/profileC22_") + HIST("4"))->Fill(col.centFT0C(), (negEtaQ*posEtaQStar).Re()/(negEtaMQ*posEtaMQ), (negEtaMQ*posEtaMQ));
          break;
        case 5:
          mHistogramQn->get<TProfile>(HIST("Qn/profileC22_") + HIST("5"))->Fill(col.centFT0C(), (negEtaQ*posEtaQStar).Re()/(negEtaMQ*posEtaMQ), (negEtaMQ*posEtaMQ));
          break;
        case 6:
          mHistogramQn->get<TProfile>(HIST("Qn/profileC22_") + HIST("6"))->Fill(col.centFT0C(), (negEtaQ*posEtaQStar).Re()/(negEtaMQ*posEtaMQ), (negEtaMQ*posEtaMQ));
          break;
        case 7:
          mHistogramQn->get<TProfile>(HIST("Qn/profileC22_") + HIST("7"))->Fill(col.centFT0C(), (negEtaQ*posEtaQStar).Re()/(negEtaMQ*posEtaMQ), (negEtaMQ*posEtaMQ));
          break;
        case 8:
          mHistogramQn->get<TProfile>(HIST("Qn/profileC22_") + HIST("8"))->Fill(col.centFT0C(), (negEtaQ*posEtaQStar).Re()/(negEtaMQ*posEtaMQ), (negEtaMQ*posEtaMQ));
          break;
        case 9:
          mHistogramQn->get<TProfile>(HIST("Qn/profileC22_") + HIST("9"))->Fill(col.centFT0C(), (negEtaQ*posEtaQStar).Re()/(negEtaMQ*posEtaMQ), (negEtaMQ*posEtaMQ));
          break; 
        default:
          return; // invalid qn bin       
      }
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
  HistogramRegistry* mHistogramQn = nullptr;       ///< For flow cumulant output
  TH2D* mReQthisEvt = nullptr; ///< For flow cumulant in an event
  TH2D* mImQthisEvt = nullptr; ///< For flow cumulant in an event
  TH2D* mReQ2thisEvt = nullptr; ///< For flow cumulant in an event
  TH2D* mImQ2thisEvt = nullptr; ///< For flow cumulant in an event
  TH2D* mMQthisEvt = nullptr;   ///< For flow cumulant in an event
  TH2D* mMQWeightthisEvt = nullptr; ///< For flow cumulant in an event 
};
} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMCOLLISIONSELECTION_H_
