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

#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/TwoParticleCorrelations/Core/FilterAndAnalysisFramework.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/TwoParticleCorrelationsSkimmed.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include <DataFormatsParameters/GRPObject.h>

#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TList.h>
#include <TParameter.h>
#include <TProfile3D.h>
#include <TROOT.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis;

namespace o2::analysis::twopcorrelations
{
#define TWOPCORRLOGCOLLISIONS debug
#define TWOPCORRLOGTRACKS debug

uint64_t collisionmask = 0UL;
std::vector<uint64_t> collisionmask_opt;
uint64_t collisionmask_forced = 0UL;
uint64_t trackmask = 0UL;
std::vector<uint64_t> trackmask_opt;
uint64_t trackmask_forced = 0UL;
uint64_t pidmask = 0UL;
std::vector<uint64_t> pidmask_opt;
uint64_t pidmask_forced = 0UL;

/// \struct trackid
/// \brief Stores the track species id, its potential correction, its potential pT average, and its global index to the track table
struct trackid {
  int id;         ///< the species internal id
  float corr;     ///< correction to apply
  float ptavg;    ///< species pT average in eta phi
  uint64_t index; ///< the index to the table row
};

PWGCF::FilterAndAnalysisFramework* fFilterFramework = nullptr;

int fMultiplicityIndex = -1; //! the index to the multiplicity values array

//============================================================================================
// The two-particle configuration objects
//============================================================================================
int ptbins = 18;
float ptlow = 0.2, ptup = 2.0;
int etabins = 16;
float etalow = -0.8, etaup = 0.8;
int zvtxbins = 40;
float zvtxlow = -10.0, zvtxup = 10.0;
int phibins = 72;
float philow = 0.0;
float phiup = constants::math::TwoPI;

float phibinshift = 0.5;
float etabinwidth = (etaup - etalow) / static_cast<float>(etabins);
float phibinwidth = (phiup - philow) / static_cast<float>(phibins);
int deltaetabins = etabins * 2 - 1;
float deltaetalow = etalow - etaup, deltaetaup = etaup - etalow;
float deltaetabinwidth = (deltaetaup - deltaetalow) / static_cast<float>(deltaetabins);
int deltaphibins = phibins;
float deltaphibinwidth = constants::math::TwoPI / deltaphibins;
float deltaphilow = 0.0 - deltaphibinwidth / 2.0;
float deltaphiup = constants::math::TwoPI - deltaphibinwidth / 2.0;

bool processpairs = false;
int bankinitialcapacity = 10000;
int bankcapacitystep = 2000;
std::vector<trackid> tracksids; // the tracks id storage

PairCuts fPairCuts;              // pair suppression engine
bool fUseConversionCuts = false; // suppress resonances and conversions
bool fUseTwoTrackCut = false;    // suppress too close tracks
} // namespace o2::analysis::twopcorrelations

using namespace o2::aod::cfskim;
using namespace o2::analysis::twopcorrelations;

struct twoParticleCorrelations {

  /* the data collecting engine */
  struct DataCollectingEngine {
    //============================================================================================
    // The TwoParticleCorrelationsAnalysisTask output objects
    //============================================================================================
    /* histograms */
    TH1F* fhVertexZA;                                        //!<! the z vertex distribution for the current multiplicity/centrality class
    std::vector<TH1F*> fhN1_vsPt;                            //!<! weighted single particle distribution vs \f$p_T\f$
    std::vector<TH2F*> fhN1_vsEtaPhi;                        //!<! weighted single particle distribution vs \f$\eta,\;\phi\f$
    std::vector<TH2F*> fhSum1Pt_vsEtaPhi;                    //!<! accumulated sum of weighted \f$p_T\f$ vs \f$\eta,\;\phi\f$
    std::vector<TH3F*> fhN1_vsZEtaPhiPt;                     //!<! single particle distribution vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$
    std::vector<TH3F*> fhSum1Pt_vsZEtaPhiPt;                 //!<! accumulated sum of weighted \f$p_T\f$ vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$
    std::vector<TH3*> fhNuaNue_vsZEtaPhiPt;                  //!<! NUA+NUE correction vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$
    std::vector<TH2*> fhPtAvg_vsEtaPhi;                      //!<! average \f$p_T\f$ vs \f$\eta,\;\phi\f$
    std::vector<std::vector<TH2F*>> fhN2_vsPtPt;             //!<! weighted two particle distribution vs \f${p_T}_1, {p_T}_2\f$
    std::vector<std::vector<TH2F*>> fhN2_vsDEtaDPhi;         //!<! two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$
    std::vector<std::vector<TH2F*>> fhSum2PtPt_vsDEtaDPhi;   //!<! two-particle  \f$\sum {p_T}_1 {p_T}_2\f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$
    std::vector<std::vector<TH2F*>> fhSum2DptDpt_vsDEtaDPhi; //!<! two-particle  \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$
    std::vector<std::vector<TH2F*>> fhSupN1N1_vsDEtaDPhi;    //!<! suppressed n1n1 two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$
    std::vector<std::vector<TH2F*>> fhSupPt1Pt1_vsDEtaDPhi;  //!<! suppressed \f${p_T}_1 {p_T}_2\f$ two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$
    /* versus centrality/multiplicity  profiles */
    std::vector<TProfile*> fhN1_vsC;                        //!<! weighted single particle distribution vs event centrality/multiplicity
    std::vector<TProfile*> fhSum1Pt_vsC;                    //!<! accumulated sum of weighted \f$p_T\f$ vs event centrality/multiplicity
    std::vector<TProfile*> fhN1nw_vsC;                      //!<! un-weighted single particle distribution vs event centrality/multiplicity
    std::vector<TProfile*> fhSum1Ptnw_vsC;                  //!<! accumulated sum of un-weighted \f$p_T\f$ vs event centrality/multiplicity
    std::vector<std::vector<TProfile*>> fhN2_vsC;           //!<! weighted accumulated two particle distribution vs event centrality/multiplicity
    std::vector<std::vector<TProfile*>> fhSum2PtPt_vsC;     //!<! weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality/multiplicity
    std::vector<std::vector<TProfile*>> fhSum2DptDpt_vsC;   //!<! weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs event centrality/multiplicity
    std::vector<std::vector<TProfile*>> fhN2nw_vsC;         //!<! un-weighted accumulated two particle distribution vs event centrality/multiplicity
    std::vector<std::vector<TProfile*>> fhSum2PtPtnw_vsC;   //!<! un-weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality/multiplicity
    std::vector<std::vector<TProfile*>> fhSum2DptDptnw_vsC; //!<! un-weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ distribution vs event centrality/multiplicity

    std::vector<std::string> tname; ///< the external track names for histogram creation
    bool ccdbstored = false;

    float isCCDBstored()
    {
      return ccdbstored;
    }

    /// \brief Returns the potentially phi origin shifted phi
    /// \param phi the track azimuthal angle
    /// \return the track phi origin shifted azimuthal angle
    float GetShiftedPhi(float phi)
    {
      using namespace twopcorrelations;
      if (!(phi < phiup)) {
        return phi - constants::math::TwoPI;
      } else {
        return phi;
      }
    }

    /// \brief Returns the zero based bin index of the eta phi passed track
    /// \param t the intended track
    /// \return the zero based bin index
    ///
    /// According to the bining structure, to the track eta will correspond
    /// a zero based bin index and similarly for the track phi
    /// The returned index is the composition of both considering eta as
    /// the first index component
    /// WARNING: for performance reasons no checks are done about the consistency
    /// of track's eta and phin with the corresponding ranges so, it is suppossed
    /// the track has been accepted and it is within that ranges
    /// IF THAT IS NOT THE CASE THE ROUTINE WILL PRODUCE NONSENSE RESULTS
    template <typename TrackObject>
    int GetEtaPhiIndex(TrackObject const& t)
    {
      using namespace twopcorrelations;

      int etaix = static_cast<int>((t.eta() - etalow) / etabinwidth);
      /* consider a potential phi origin shift */
      float phi = GetShiftedPhi(t.phi());
      int phiix = static_cast<int>((phi - philow) / phibinwidth);
      return etaix * phibins + phiix;
    }

    /// \brief Returns the TH2 global index for the differential histograms
    /// \param t1 the intended track one
    /// \param t2 the intended track two
    /// \return the globl TH2 bin for delta eta delta phi
    ///
    /// WARNING: for performance reasons no checks are done about the consistency
    /// of tracks' eta and phi within the corresponding ranges so, it is suppossed
    /// the tracks have been accepted and they are within that ranges
    /// IF THAT IS NOT THE CASE THE ROUTINE WILL PRODUCE NONSENSE RESULTS
    template <typename TrackObject>
    int GetDEtaDPhiGlobalIndex(TrackObject const& t1, TrackObject const& t2)
    {
      using namespace twopcorrelations;

      /* rule: ix are always zero based while bins are always one based */
      int etaix_1 = static_cast<int>((t1.eta() - etalow) / etabinwidth);
      /* consider a potential phi origin shift */
      float phi = GetShiftedPhi(t1.phi());
      int phiix_1 = static_cast<int>((phi - philow) / phibinwidth);
      int etaix_2 = static_cast<int>((t2.eta() - etalow) / etabinwidth);
      /* consider a potential phi origin shift */
      phi = GetShiftedPhi(t2.phi());
      int phiix_2 = static_cast<int>((phi - philow) / phibinwidth);

      int deltaeta_ix = etaix_1 - etaix_2 + etabins - 1;
      int deltaphi_ix = phiix_1 - phiix_2;
      if (deltaphi_ix < 0) {
        deltaphi_ix += phibins;
      }

      return fhN2_vsDEtaDPhi[0][0]->GetBin(deltaeta_ix + 1, deltaphi_ix + 1);
    }

    void storeTrackCorrections(std::vector<TH3*> corrs)
    {
      LOGF(info, "Stored NUA&NUE corrections for %d track species", corrs.size());
      fhNuaNue_vsZEtaPhiPt = corrs;
      ccdbstored = true;
    }

    void storePtAverages(std::vector<TH2*> ptavgs)
    {
      LOGF(info, "Stored pT average for %d track species", ptavgs.size());
      fhPtAvg_vsEtaPhi = ptavgs;
      ccdbstored = true;
    }

    template <typename TrackListObject>
    void getTrackPtAvgAndCorrections(TrackListObject& tracks, std::vector<trackid>& tracksids, float zvtx)
    {
      for (auto& tid : tracksids) {
        auto const& track = tracks.iteratorAt(tid.index);
        if (!(tid.id < 0)) {
          if (fhNuaNue_vsZEtaPhiPt[tid.id] != nullptr) {
            tid.corr = fhNuaNue_vsZEtaPhiPt[tid.id]->GetBinContent(zvtx, GetEtaPhiIndex(track) + 0.5, track.pt());
          } else {
            tid.corr = 1;
          }
          if (fhPtAvg_vsEtaPhi[tid.id] != nullptr) {
            tid.ptavg = fhPtAvg_vsEtaPhi[tid.id]->GetBinContent(fhPtAvg_vsEtaPhi[tid.id]->FindBin(track.eta(), track.phi()));
          } else {
            tid.ptavg = 0;
          }
        }
      }
    }

    /// \brief fills the singles histograms in singles execution mode
    /// \param passedtracks filtered table with the accepted tracks
    /// \param tracksids id, index to the track table, correction, and average pT the passed track table
    /// \param zvtx the z vertex coordinate for singles analysis on zvtx bins
    template <typename TrackListObject>
    void processSingles(TrackListObject& passedtracks, std::vector<trackid>& tracksids, float zvtx)
    {
      for (auto& tid : tracksids) {
        auto const& track = passedtracks.iteratorAt(tid.index);
        if (!(tid.id < 0)) {
          fhN1_vsPt[tid.id]->Fill(track.pt(), tid.corr);
          fhN1_vsZEtaPhiPt[tid.id]->Fill(zvtx, GetEtaPhiIndex(track) + 0.5, track.pt(), tid.corr);
          fhSum1Pt_vsZEtaPhiPt[tid.id]->Fill(zvtx, GetEtaPhiIndex(track) + 0.5, track.pt(), track.pt() * tid.corr);
        }
      }
    }

    /// \brief fills the singles histograms in pair execution mode
    /// \param passedtracks filtered table with the accepted tracks
    /// \param tracksids id, index to the track table, correction, and average pT the passed track table
    /// \param cmul centrality - multiplicity for the collision being analyzed
    template <typename TrackListObject>
    void processTracks(TrackListObject& passedtracks, std::vector<trackid>& tracksids, float cmul)
    {
      LOGF(TWOPCORRLOGCOLLISIONS, "Processing %d tracks in a collision with cent/mult %f ", passedtracks.size(), cmul);

      /* process magnitudes */
      std::vector<double> n1(tname.size(), 0.0);       ///< weighted number of tracks for current collision
      std::vector<double> sum1Pt(tname.size(), 0.0);   ///< accumulated sum of weighted track 1 \f$p_T\f$ for current collision
      std::vector<double> n1nw(tname.size(), 0.0);     ///< not weighted number of tracks for current collision
      std::vector<double> sum1Ptnw(tname.size(), 0.0); ///< accumulated sum of not weighted track \f$p_T\f$ for current collision

      for (auto& tid : tracksids) {
        auto const& track = passedtracks.iteratorAt(tid.index);
        if (!(tid.id < 0)) {
          n1[tid.id] += tid.corr;
          sum1Pt[tid.id] += track.pt() * tid.corr;
          n1nw[tid.id] += 1;
          sum1Ptnw[tid.id] += track.pt();

          fhN1_vsEtaPhi[tid.id]->Fill(track.eta(), GetShiftedPhi(track.phi()), tid.corr);
          fhSum1Pt_vsEtaPhi[tid.id]->Fill(track.eta(), GetShiftedPhi(track.phi()), track.pt() * tid.corr);
        }
      }
      for (uint i = 0; i < n1.size(); ++i) {
        fhN1_vsC[i]->Fill(cmul, n1[i]);
        fhSum1Pt_vsC[i]->Fill(cmul, sum1Pt[i]);
        fhN1nw_vsC[i]->Fill(cmul, n1nw[i]);
        fhSum1Ptnw_vsC[i]->Fill(cmul, sum1Ptnw[i]);
      }
    }

    /// \brief fills the pair histograms in pair execution mode
    /// \param trks1 filtered table with the tracks associated to the first track in the pair
    /// \param trks2 filtered table with the tracks associated to the second track in the pair
    /// \param pix index, in the track combination histogram bank, for the passed filetered track tables
    /// \param cmul centrality - multiplicity for the collision being analyzed
    /// Be aware that at least in half of the cases traks1 and trks2 will have the same content
    template <typename TrackListObject>
    void processTrackPairs(TrackListObject& tracks, std::vector<trackid>& tracksids, float cmul, int bfield)
    {
      using namespace twopcorrelations;

      /* process pair magnitudes */
      std::vector<std::vector<double>> n2(tname.size(), std::vector<double>(tname.size(), 0.0));           ///< weighted number of track pairs for current collision
      std::vector<std::vector<double>> n2sup(tname.size(), std::vector<double>(tname.size(), 0.0));        ///< weighted number of track suppressed pairs for current collision
      std::vector<std::vector<double>> sum2PtPt(tname.size(), std::vector<double>(tname.size(), 0.0));     ///< accumulated sum of weighted pairs \f${p_T}_1 {p_T}_2\f$ for current collision
      std::vector<std::vector<double>> sum2DptDpt(tname.size(), std::vector<double>(tname.size(), 0.0));   ///< accumulated sum of weighted pairs \f$\Delta p_T \Delta p_T\f$ for current collision
      std::vector<std::vector<double>> n2nw(tname.size(), std::vector<double>(tname.size(), 0.0));         ///< not weighted number of track pairs for current collision
      std::vector<std::vector<double>> sum2PtPtnw(tname.size(), std::vector<double>(tname.size(), 0.0));   ///< accumulated sum of not weighted paisr  \f${p_T}_1 {p_T}_2\f$ for current collision
      std::vector<std::vector<double>> sum2DptDptnw(tname.size(), std::vector<double>(tname.size(), 0.0)); ///< accumulated sum of not weighted pairs \f$\Delta p_T \Delta p_T\f$ for current collision

      for (auto& tid1 : tracksids) {
        auto const& track1 = tracks.iteratorAt(tid1.index);
        if (!(tid1.id < 0)) {
          for (auto& tid2 : tracksids) {
            auto const& track2 = tracks.iteratorAt(tid2.index);
            if (!(tid2.id < 0)) {
              if (tid1.index == tid2.index) {
                /* exclude auto correlations */
                continue;
              }
              /* get the global bin for filling the differential histograms */
              int globalbin = GetDEtaDPhiGlobalIndex(track1, track2);
              if ((fUseConversionCuts && fPairCuts.conversionCuts(track1, track2)) || (fUseTwoTrackCut && fPairCuts.twoTrackCut(track1, track2, bfield))) {
                /* suppress the pair */
                fhSupN1N1_vsDEtaDPhi[tid1.id][tid2.id]->AddBinContent(globalbin, tid1.corr * tid2.corr);
                fhSupPt1Pt1_vsDEtaDPhi[tid1.id][tid2.id]->AddBinContent(globalbin, track1.pt() * track2.pt() * tid1.corr * tid2.corr);
                n2sup[tid1.id][tid2.id] += tid1.corr * tid2.corr;
              } else {
                /* count the pair */
                n2[tid1.id][tid2.id] += tid1.corr * tid2.corr;
                sum2PtPt[tid1.id][tid2.id] += track1.pt() * track2.pt() * tid1.corr * tid2.corr;
                sum2DptDpt[tid1.id][tid2.id] += (tid1.corr * track1.pt() - tid1.ptavg) * (tid2.corr * track2.pt() - tid2.ptavg);
                n2nw[tid1.id][tid2.id] += 1;
                sum2PtPtnw[tid1.id][tid2.id] += track1.pt() * track2.pt();
                sum2DptDptnw[tid1.id][tid2.id] += (track1.pt() - tid1.ptavg) * (track2.pt() - tid2.ptavg);

                fhN2_vsDEtaDPhi[tid1.id][tid2.id]->AddBinContent(globalbin, tid1.corr * tid2.corr);
                fhSum2DptDpt_vsDEtaDPhi[tid1.id][tid2.id]->AddBinContent(globalbin, (tid1.corr * track1.pt() - tid1.ptavg) * (tid2.corr * track2.pt() - tid2.ptavg));
                fhSum2PtPt_vsDEtaDPhi[tid1.id][tid2.id]->AddBinContent(globalbin, track1.pt() * track2.pt() * tid1.corr * tid2.corr);
                fhN2_vsPtPt[tid1.id][tid2.id]->Fill(track1.pt(), track2.pt(), tid1.corr * tid2.corr);
              }
            }
          }
        }
      }
      for (uint i = 0; i < tname.size(); ++i) {
        for (uint j = 0; j < tname.size(); ++j) {
          fhN2_vsC[i][j]->Fill(cmul, n2[i][j]);
          fhSum2PtPt_vsC[i][j]->Fill(cmul, sum2PtPt[i][j]);
          fhSum2DptDpt_vsC[i][j]->Fill(cmul, sum2DptDpt[i][j]);
          fhN2nw_vsC[i][j]->Fill(cmul, n2nw[i][j]);
          fhSum2PtPtnw_vsC[i][j]->Fill(cmul, sum2PtPtnw[i][j]);
          fhSum2DptDptnw_vsC[i][j]->Fill(cmul, sum2DptDptnw[i][j]);
          /* let's also update the number of entries in the differential histograms */
          fhN2_vsDEtaDPhi[i][j]->SetEntries(fhN2_vsDEtaDPhi[i][j]->GetEntries() + n2[i][j]);
          fhSum2DptDpt_vsDEtaDPhi[i][j]->SetEntries(fhSum2DptDpt_vsDEtaDPhi[i][j]->GetEntries() + n2[i][j]);
          fhSum2PtPt_vsDEtaDPhi[i][j]->SetEntries(fhSum2PtPt_vsDEtaDPhi[i][j]->GetEntries() + n2[i][j]);
          fhSupN1N1_vsDEtaDPhi[i][j]->SetEntries(fhSupN1N1_vsDEtaDPhi[i][j]->GetEntries() + n2sup[i][j]);
          fhSupPt1Pt1_vsDEtaDPhi[i][j]->SetEntries(fhSupPt1Pt1_vsDEtaDPhi[i][j]->GetEntries() + n2sup[i][j]);
        }
      }
    }

    template <typename TrackListObject>
    void processCollision(TrackListObject& Tracks, float zvtx, float centmult, int bfield)
    {
      /* TODO: the centrality should be chosen non detector dependent */
      using namespace twopcorrelations;
      /* TODO: in here we should collect the ids for the selected tracks and resize the bank if needed */

      /* now let's fill the potential correction and pT averages for the current zvtx */
      getTrackPtAvgAndCorrections(Tracks, tracksids, zvtx);

      if (!processpairs) {
        /* process single tracks */
        fhVertexZA->Fill(zvtx);
        processSingles(Tracks, tracksids, zvtx);
      } else {
        processTracks(Tracks, tracksids, centmult); /* track one */
        /* process pair magnitudes */
        processTrackPairs(Tracks, tracksids, centmult, bfield);
      }
    }

    void init(TList* fOutputList, std::vector<std::string> idnames)
    {
      LOGF(info, "Correlations processing engine::init()");
      using namespace twopcorrelations;

      /* create the histograms */
      Bool_t oldstatus = TH1::AddDirectoryStatus();
      TH1::AddDirectory(kFALSE);

      /* load the species names */
      for (auto id : idnames) {
        LOGF(info, "Adding particle species %s", id.c_str());
        tname.push_back(std::string(id.c_str()));
      }
      /* reserve the histograms place holders */
      for (uint i = 0; i < tname.size(); ++i) {
        if (!processpairs) {
          fhN1_vsPt.reserve(tname.size());
          fhN1_vsZEtaPhiPt.reserve(tname.size());
          fhSum1Pt_vsZEtaPhiPt.reserve(tname.size());
          fhNuaNue_vsZEtaPhiPt.reserve(tname.size());
          fhPtAvg_vsEtaPhi.reserve(tname.size());
        } else {
          fhN1_vsEtaPhi.reserve(tname.size());
          fhSum1Pt_vsEtaPhi.reserve(tname.size());
          fhN1_vsC.reserve(tname.size());
          fhSum1Pt_vsC.reserve(tname.size());
          fhN1nw_vsC.reserve(tname.size());
          fhSum1Ptnw_vsC.reserve(tname.size());
          fhNuaNue_vsZEtaPhiPt.reserve(tname.size());
          fhPtAvg_vsEtaPhi.reserve(tname.size());
          fhN2_vsDEtaDPhi.reserve(tname.size());
          fhSum2PtPt_vsDEtaDPhi.reserve(tname.size());
          fhSum2DptDpt_vsDEtaDPhi.reserve(tname.size());
          fhSupN1N1_vsDEtaDPhi.reserve(tname.size());
          fhSupPt1Pt1_vsDEtaDPhi.reserve(tname.size());
          fhN2_vsPtPt.reserve(tname.size());
          fhN2_vsC.reserve(tname.size());
          fhSum2PtPt_vsC.reserve(tname.size());
          fhSum2DptDpt_vsC.reserve(tname.size());
          fhN2nw_vsC.reserve(tname.size());
          fhSum2PtPtnw_vsC.reserve(tname.size());
          fhSum2DptDptnw_vsC.reserve(tname.size());
          for (uint j = 0; j < tname.size(); ++j) {
            fhN2_vsDEtaDPhi.push_back({tname.size(), nullptr});
            fhSum2PtPt_vsDEtaDPhi.push_back({tname.size(), nullptr});
            fhSum2DptDpt_vsDEtaDPhi.push_back({tname.size(), nullptr});
            fhSupN1N1_vsDEtaDPhi.push_back({tname.size(), nullptr});
            fhSupPt1Pt1_vsDEtaDPhi.push_back({tname.size(), nullptr});
            fhN2_vsPtPt.push_back({tname.size(), nullptr});
            fhN2_vsC.push_back({tname.size(), nullptr});
            fhSum2PtPt_vsC.push_back({tname.size(), nullptr});
            fhSum2DptDpt_vsC.push_back({tname.size(), nullptr});
            fhN2nw_vsC.push_back({tname.size(), nullptr});
            fhSum2PtPtnw_vsC.push_back({tname.size(), nullptr});
            fhSum2DptDptnw_vsC.push_back({tname.size(), nullptr});
          }
        }
      }

      if (!processpairs) {
        fhVertexZA = new TH1F("VertexZA", "Vertex Z; z_{vtx}", zvtxbins, zvtxlow, zvtxup);
        fOutputList->Add(fhVertexZA);
        for (uint i = 0; i < tname.size(); ++i) {
          /* histograms for each track, one and two */
          fhN1_vsPt[i] = new TH1F(TString::Format("n1_%s_vsPt", tname[i].c_str()).Data(),
                                  TString::Format("#LT n_{1} #GT;p_{t,%s} (GeV/c);#LT n_{1} #GT", tname[i].c_str()).Data(),
                                  ptbins, ptlow, ptup);
          /* we don't want the Sumw2 structure being created here */
          bool defSumw2 = TH1::GetDefaultSumw2();
          TH1::SetDefaultSumw2(false);
          fhN1_vsZEtaPhiPt[i] = new TH3F(TString::Format("n1_%s_vsZ_vsEtaPhi_vsPt", tname[i].c_str()).Data(),
                                         TString::Format("#LT n_{1} #GT;vtx_{z};#eta_{%s}#times#varphi_{%s};p_{t,%s} (GeV/c)", tname[i].c_str(), tname[i].c_str(), tname[i].c_str()).Data(),
                                         zvtxbins, zvtxlow, zvtxup, etabins * phibins, 0.0, static_cast<double>(etabins * phibins), ptbins, ptlow, ptup);
          fhSum1Pt_vsZEtaPhiPt[i] = new TH3F(TString::Format("sumPt1_%s_vsZ_vsEtaPhi_vsPt", tname[i].c_str()).Data(),
                                             TString::Format("#LT #Sigma p_{t,%s}#GT;vtx_{z};#eta_{%s}#times#varphi_{%s};p_{t,%s} (GeV/c)", tname[i].c_str(), tname[i].c_str(), tname[i].c_str(), tname[i].c_str()).Data(),
                                             zvtxbins, zvtxlow, zvtxup, etabins * phibins, 0.0, static_cast<double>(etabins * phibins), ptbins, ptlow, ptup);
          /* we return it back to previuos state */
          TH1::SetDefaultSumw2(defSumw2);

          /* the statistical uncertainties will be estimated by the subsamples method so let's get rid of the error tracking */
          fhN1_vsZEtaPhiPt[i]->SetBit(TH1::kIsNotW);
          fhN1_vsZEtaPhiPt[i]->Sumw2(false);
          fhSum1Pt_vsZEtaPhiPt[i]->SetBit(TH1::kIsNotW);
          fhSum1Pt_vsZEtaPhiPt[i]->Sumw2(false);
          fhNuaNue_vsZEtaPhiPt[i] = nullptr;
          fhPtAvg_vsEtaPhi[i] = nullptr;

          fOutputList->Add(fhN1_vsPt[i]);
          fOutputList->Add(fhN1_vsZEtaPhiPt[i]);
          fOutputList->Add(fhSum1Pt_vsZEtaPhiPt[i]);
        }
      } else {
        for (uint i = 0; i < tname.size(); ++i) {
          /* histograms for each track species */
          fhN1_vsEtaPhi[i] = new TH2F(TString::Format("n1_%s_vsEtaPhi", tname[i].c_str()).Data(),
                                      TString::Format("#LT n_{1} #GT;#eta_{%s};#varphi_{%s} (radian);#LT n_{1} #GT", tname[i].c_str(), tname[i].c_str()).Data(),
                                      etabins, etalow, etaup, phibins, philow, phiup);
          fhSum1Pt_vsEtaPhi[i] = new TH2F(TString::Format("sumPt_%s_vsEtaPhi", tname[i].c_str()).Data(),
                                          TString::Format("#LT #Sigma p_{t,%s} #GT;#eta_{%s};#varphi_{%s} (radian);#LT #Sigma p_{t,%s} #GT (GeV/c)",
                                                          tname[i].c_str(), tname[i].c_str(), tname[i].c_str(), tname[i].c_str())
                                            .Data(),
                                          etabins, etalow, etaup, phibins, philow, phiup);
          fhN1_vsC[i] = new TProfile(TString::Format("n1_%s_vsM", tname[i].c_str()).Data(),
                                     TString::Format("#LT n_{1} #GT (weighted);Centrality/Multiplicity (%%);#LT n_{1} #GT").Data(),
                                     100, 0.0, 100.0);
          fhSum1Pt_vsC[i] = new TProfile(TString::Format("sumPt_%s_vsM", tname[i].c_str()),
                                         TString::Format("#LT #Sigma p_{t,%s} #GT (weighted);Centrality/Multiplicity (%%);#LT #Sigma p_{t,%s} #GT (GeV/c)", tname[i].c_str(), tname[i].c_str()).Data(),
                                         100, 0.0, 100.0);
          fhN1nw_vsC[i] = new TProfile(TString::Format("n1Nw_%s_vsM", tname[i].c_str()).Data(),
                                       TString::Format("#LT n_{1} #GT;Centrality/Multiplicity (%%);#LT n_{1} #GT").Data(),
                                       100, 0.0, 100.0);
          fhSum1Ptnw_vsC[i] = new TProfile(TString::Format("sumPtNw_%s_vsM", tname[i].c_str()).Data(),
                                           TString::Format("#LT #Sigma p_{t,%s} #GT;Centrality/Multiplicity (%%);#LT #Sigma p_{t,%s} #GT (GeV/c)", tname[i].c_str(), tname[i].c_str()).Data(), 100, 0.0, 100.0);
          fhNuaNue_vsZEtaPhiPt[i] = nullptr;
          fhPtAvg_vsEtaPhi[i] = nullptr;
          fOutputList->Add(fhN1_vsEtaPhi[i]);
          fOutputList->Add(fhSum1Pt_vsEtaPhi[i]);
          fOutputList->Add(fhN1_vsC[i]);
          fOutputList->Add(fhSum1Pt_vsC[i]);
          fOutputList->Add(fhN1nw_vsC[i]);
          fOutputList->Add(fhSum1Ptnw_vsC[i]);
          for (uint j = 0; j < tname.size(); ++j) {
            /* histograms for each track pair combination */
            /* we don't want the Sumw2 structure being created here */
            bool defSumw2 = TH1::GetDefaultSumw2();
            TH1::SetDefaultSumw2(false);
            char pname[256];
            snprintf(pname, tname[i].length() + tname[j].length() + 1, "%s%s", tname[i].c_str(), tname[j].c_str());
            fhN2_vsDEtaDPhi[i][j] = new TH2F(TString::Format("n2_12_vsDEtaDPhi_%s", pname), TString::Format("#LT n_{2} #GT (%s);#Delta#eta;#Delta#varphi;#LT n_{2} #GT", pname),
                                             deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
            fhSum2PtPt_vsDEtaDPhi[i][j] = new TH2F(TString::Format("sumPtPt_12_vsDEtaDPhi_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s);#Delta#eta;#Delta#varphi;#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname),
                                                   deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
            fhSum2DptDpt_vsDEtaDPhi[i][j] = new TH2F(TString::Format("sumDptDpt_12_vsDEtaDPhi_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s);#Delta#eta;#Delta#varphi;#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname),
                                                     deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
            fhSupN1N1_vsDEtaDPhi[i][j] = new TH2F(TString::Format("suppn1n1_12_vsDEtaDPhi_%s", pname), TString::Format("Suppressed #LT n_{1} #GT#LT n_{1} #GT (%s);#Delta#eta;#Delta#varphi;#LT n_{1} #GT#LT n_{1} #GT", pname),
                                                  deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
            fhSupPt1Pt1_vsDEtaDPhi[i][j] = new TH2F(TString::Format("suppPtPt_12_vsDEtaDPhi_%s", pname), TString::Format("Suppressed #LT p_{t,1} #GT#LT p_{t,2} #GT (%s);#Delta#eta;#Delta#varphi;#LT p_{t,1} #GT#LT p_{t,2} #GT (GeV^{2})", pname),
                                                    deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
            /* we return it back to previuos state */
            TH1::SetDefaultSumw2(defSumw2);

            fhN2_vsPtPt[i][j] = new TH2F(TString::Format("n2_12_vsPtVsPt_%s", pname), TString::Format("#LT n_{2} #GT (%s);p_{t,1} (GeV/c);p_{t,2} (GeV/c);#LT n_{2} #GT", pname),
                                         ptbins, ptlow, ptup, ptbins, ptlow, ptup);

            fhN2_vsC[i][j] = new TProfile(TString::Format("n2_12_vsM_%s", pname), TString::Format("#LT n_{2} #GT (%s) (weighted);Centrality/Multiplicity (%%);#LT n_{2} #GT", pname), 100, 0.0, 100.0);
            fhSum2PtPt_vsC[i][j] = new TProfile(TString::Format("sumPtPt_12_vsM_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s) (weighted);Centrality/Multiplicity (%%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname), 100, 0.0, 100.0);
            fhSum2DptDpt_vsC[i][j] = new TProfile(TString::Format("sumDptDpt_12_vsM_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s) (weighted);Centrality/Multiplicity (%%);#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname), 100, 0.0, 100.0);
            fhN2nw_vsC[i][j] = new TProfile(TString::Format("n2Nw_12_vsM_%s", pname), TString::Format("#LT n_{2} #GT (%s);Centrality/Multiplicity (%%);#LT n_{2} #GT", pname), 100, 0.0, 100.0);
            fhSum2PtPtnw_vsC[i][j] = new TProfile(TString::Format("sumPtPtNw_12_vsM_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s);Centrality/Multiplicity (%%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname), 100, 0.0, 100.0);
            fhSum2DptDptnw_vsC[i][j] = new TProfile(TString::Format("sumDptDptNw_12_vsM_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s);Centrality/Multiplicity (%%);#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname), 100, 0.0, 100.0);

            /* the statistical uncertainties will be estimated by the subsamples method so let's get rid of the error tracking */
            fhN2_vsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhN2_vsDEtaDPhi[i][j]->Sumw2(false);
            fhSum2PtPt_vsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhSum2PtPt_vsDEtaDPhi[i][j]->Sumw2(false);
            fhSum2DptDpt_vsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhSum2DptDpt_vsDEtaDPhi[i][j]->Sumw2(false);
            fhSupN1N1_vsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhSupN1N1_vsDEtaDPhi[i][j]->Sumw2(false);
            fhSupPt1Pt1_vsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhSupPt1Pt1_vsDEtaDPhi[i][j]->Sumw2(false);

            fOutputList->Add(fhN2_vsDEtaDPhi[i][j]);
            fOutputList->Add(fhSum2PtPt_vsDEtaDPhi[i][j]);
            fOutputList->Add(fhSum2DptDpt_vsDEtaDPhi[i][j]);
            fOutputList->Add(fhSupN1N1_vsDEtaDPhi[i][j]);
            fOutputList->Add(fhSupPt1Pt1_vsDEtaDPhi[i][j]);
            fOutputList->Add(fhN2_vsPtPt[i][j]);
            fOutputList->Add(fhN2_vsC[i][j]);
            fOutputList->Add(fhSum2PtPt_vsC[i][j]);
            fOutputList->Add(fhSum2DptDpt_vsC[i][j]);
            fOutputList->Add(fhN2nw_vsC[i][j]);
            fOutputList->Add(fhSum2PtPtnw_vsC[i][j]);
            fOutputList->Add(fhSum2DptDptnw_vsC[i][j]);
          }
        }
      }
      TH1::AddDirectory(oldstatus);
    }
  }; // DataCollectingEngine

/* the skimming configuration */
#include "PWGCF/TwoParticleCorrelations/TableProducer/Productions/skimmingconf_20221115.h" // NOLINT

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  /* the data memebers for this task */
  /* the centrality / multiplicity limits for collecting data in this task instance */
  int ncmranges = 0;
  float* fCentMultMin = nullptr;
  float* fCentMultMax = nullptr;

  /* the data collecting engine instances */
  std::vector<DataCollectingEngine*> dataCE;

  /* the input file structure from CCDB */
  TList* ccdblst = nullptr;
  bool loadfromccdb = false;

  /* pair conversion suppression defaults */
  static constexpr float cfgPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};
  Configurable<LabeledArray<float>> cfgPairCut{"paircut", {cfgPairCutDefaults[0], 5, {"Photon", "K0", "Lambda", "Phi", "Rho"}}, "Conversion suppressions"};
  /* two tracks cut */
  Configurable<float> cfgTwoTrackCut{"twotrackcut", -1, "Two-tracks cut: -1 = off; >0 otherwise distance value (suggested: 0.02"};
  Configurable<float> cfgTwoTrackCutMinRadius{"twotrackcutminradius", 0.8f, "Two-tracks cut: radius in m from which two-tracks cut is applied"};

  Configurable<bool> cfgProcessPairs{"processpairs", false, "Process pairs: false = no, just singles, true = yes, process pairs"};
  Configurable<std::string> cfgCentSpec{"centralities", "00-05,05-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80", "Centrality/multiplicity ranges in min-max separated by commas"};

  Configurable<o2::analysis::DptDptBinningCuts> cfgBinning{"binning",
                                                           {28, -7.0, 7.0, 18, 0.2, 2.0, 16, -0.8, 0.8, 72, 0.5},
                                                           "triplets - nbins, min, max - for z_vtx, pT, eta and phi, binning plus bin fraction of phi origin shift"};
  struct : ConfigurableGroup {
    Configurable<std::string> cfgCCDBUrl{"input_ccdburl", "http://ccdb-test.cern.ch:8080", "The CCDB url for the input file"};
    Configurable<std::string> cfgCCDBPathName{"input_ccdbpath", "", "The CCDB path for the input file. Default \"\", i.e. don't load from CCDB"};
    Configurable<std::string> cfgCCDBDate{"input_ccdbdate", "20220307", "The CCDB date for the input file"};
  } cfginputfile;

  OutputObj<TList> fOutput{"DptDptCorrelationsData", OutputObjHandlingPolicy::AnalysisObject, OutputObjSourceType::OutputObjSource};

  void init(InitContext const&)
  {
    using namespace twopcorrelations;
    LOGF(info, "twoParticleCorrelations::init()");

    /* update with the configurable values */
    ptbins = cfgBinning->mPTbins;
    ptlow = cfgBinning->mPTmin;
    ptup = cfgBinning->mPTmax;
    etabins = cfgBinning->mEtabins;
    etalow = cfgBinning->mEtamin;
    etaup = cfgBinning->mEtamax;
    zvtxbins = cfgBinning->mZVtxbins;
    zvtxlow = cfgBinning->mZVtxmin;
    zvtxup = cfgBinning->mZVtxmax;
    phibins = cfgBinning->mPhibins;
    philow = 0.0f;
    phiup = constants::math::TwoPI;
    phibinshift = cfgBinning->mPhibinshift;
    processpairs = cfgProcessPairs.value;
    loadfromccdb = cfginputfile.cfgCCDBPathName->length() > 0;
    /* update the potential binning change */
    etabinwidth = (etaup - etalow) / static_cast<float>(etabins);
    phibinwidth = (phiup - philow) / static_cast<float>(phibins);

    /* the differential bining */
    deltaetabins = etabins * 2 - 1;
    deltaetalow = etalow - etaup, deltaetaup = etaup - etalow;
    deltaetabinwidth = (deltaetaup - deltaetalow) / static_cast<float>(deltaetabins);
    deltaphibins = phibins;
    deltaphibinwidth = constants::math::TwoPI / deltaphibins;
    deltaphilow = 0.0 - deltaphibinwidth / 2.0;
    deltaphiup = constants::math::TwoPI - deltaphibinwidth / 2.0;

    /* create the output directory which will own the task output */
    TList* fGlobalOutputList = new TList();
    fGlobalOutputList->SetOwner(true);
    fOutput.setObject(fGlobalOutputList);

    /* incorporate configuration parameters to the output */
    fGlobalOutputList->Add(new TParameter<Int_t>("NoBinsVertexZ", zvtxbins, 'f'));
    fGlobalOutputList->Add(new TParameter<Int_t>("NoBinsPt", ptbins, 'f'));
    fGlobalOutputList->Add(new TParameter<Int_t>("NoBinsEta", etabins, 'f'));
    fGlobalOutputList->Add(new TParameter<Int_t>("NoBinsPhi", phibins, 'f'));
    fGlobalOutputList->Add(new TParameter<Double_t>("MinVertexZ", zvtxlow, 'f'));
    fGlobalOutputList->Add(new TParameter<Double_t>("MaxVertexZ", zvtxup, 'f'));
    fGlobalOutputList->Add(new TParameter<Double_t>("MinPt", ptlow, 'f'));
    fGlobalOutputList->Add(new TParameter<Double_t>("MaxPt", ptup, 'f'));
    fGlobalOutputList->Add(new TParameter<Double_t>("MinEta", etalow, 'f'));
    fGlobalOutputList->Add(new TParameter<Double_t>("MaxEta", etaup, 'f'));
    fGlobalOutputList->Add(new TParameter<Double_t>("MinPhi", philow, 'f'));
    fGlobalOutputList->Add(new TParameter<Double_t>("MaxPhi", phiup, 'f'));
    fGlobalOutputList->Add(new TParameter<Double_t>("PhiBinShift", phibinshift, 'f'));
    fGlobalOutputList->Add(new TParameter<Bool_t>("DifferentialOutput", true, 'f'));

    /* after the parameters dump the proper phi limits are set according to the phi shift */
    phiup = phiup - phibinwidth * phibinshift;
    philow = philow - phibinwidth * phibinshift;

    /* create the data collecting engine instances according to the configured centrality/multiplicity ranges */
    {
      TObjArray* tokens = TString(cfgCentSpec.value.c_str()).Tokenize(",");
      ncmranges = tokens->GetEntries();
      fCentMultMin = new float[ncmranges];
      fCentMultMax = new float[ncmranges];
      dataCE.resize(ncmranges);

      for (int i = 0; i < ncmranges; ++i) {
        float cmmin = 0.0f;
        float cmmax = 0.0f;
        sscanf(tokens->At(i)->GetName(), "%f-%f", &cmmin, &cmmax);
        fCentMultMin[i] = cmmin;
        fCentMultMax[i] = cmmax;
        dataCE[i] = new DataCollectingEngine();

        /* create the output list for the current centrality range */
        TList* fOutputList = new TList();
        fOutputList->SetName(TString::Format("DptDptCorrelationsData-%s", tokens->At(i)->GetName()));
        fOutputList->SetOwner(true);
        /* init the data collection instance */
        /* TODO: the list of species of interest */
        dataCE[i]->init(fOutputList, std::vector<std::string>{"1", "2"});
        fGlobalOutputList->Add(fOutputList);
      }
      delete tokens;
      for (int i = 0; i < ncmranges; ++i) {
        LOGF(info, " centrality/multipliicty range: %d, low limit: %f, up limit: %f", i, fCentMultMin[i], fCentMultMax[i]);
      }
    }
    /* two-track cut and conversion suppression */
    fPairCuts.SetHistogramRegistry(nullptr); // not histogram registry for the time being, incompatible with TList when it is empty
    if (processpairs && (cfgPairCut->get("Photon") > 0 || cfgPairCut->get("K0") > 0 || cfgPairCut->get("Lambda") > 0 || cfgPairCut->get("Phi") > 0 || cfgPairCut->get("Rho") > 0)) {
      fPairCuts.SetPairCut(PairCuts::Photon, cfgPairCut->get("Photon"));
      fPairCuts.SetPairCut(PairCuts::K0, cfgPairCut->get("K0"));
      fPairCuts.SetPairCut(PairCuts::Lambda, cfgPairCut->get("Lambda"));
      fPairCuts.SetPairCut(PairCuts::Phi, cfgPairCut->get("Phi"));
      fPairCuts.SetPairCut(PairCuts::Rho, cfgPairCut->get("Rho"));
      fUseConversionCuts = true;
    }
    if (processpairs && (cfgTwoTrackCut > 0)) {
      fPairCuts.SetTwoTrackCuts(cfgTwoTrackCut, cfgTwoTrackCutMinRadius);
      fUseTwoTrackCut = true;
    }

    /* initialize access to the CCDB */
    ccdb->setURL(cfginputfile.cfgCCDBUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    /* let's initialize the skimming input machinery */
    /* collision filtering configuration */
    PWGCF::EventSelectionConfigurable eventsel(eventfilter.bfield, eventfilter.centmultsel, {}, eventfilter.zvtxsel, eventfilter.pileuprej);
    /* track filtering configuration */
    PWGCF::TrackSelectionConfigurable trksel(trackfilter.ttype, trackfilter.nclstpc, trackfilter.nxrtpc, trackfilter.nclsits, trackfilter.chi2clustpc,
                                             trackfilter.chi2clusits, trackfilter.xrofctpc, trackfilter.dcaxy, trackfilter.dcaz, trackfilter.ptrange, trackfilter.etarange);
#ifdef INCORPORATEBAYESIANPID
    PWGCF::PIDSelectionConfigurable pidsel(pidfilter.pidtpcfilter.tpcel, pidfilter.pidtpcfilter.tpcmu, pidfilter.pidtpcfilter.tpcpi, pidfilter.pidtpcfilter.tpcka, pidfilter.pidtpcfilter.tpcpr,
                                           pidfilter.pidtoffilter.tpcel, pidfilter.pidtoffilter.tpcmu, pidfilter.pidtoffilter.tpcpi, pidfilter.pidtoffilter.tpcka, pidfilter.pidtoffilter.tpcpr,
                                           pidfilter.pidbayesfilter.bayel, pidfilter.pidbayesfilter.baymu, pidfilter.pidbayesfilter.baypi, pidfilter.pidbayesfilter.bayka, pidfilter.pidbayesfilter.baypr);
#else
    PWGCF::PIDSelectionConfigurable pidsel(pidfilter.pidtpcfilter.tpcel, pidfilter.pidtpcfilter.tpcmu, pidfilter.pidtpcfilter.tpcpi, pidfilter.pidtpcfilter.tpcka, pidfilter.pidtpcfilter.tpcpr,
                                           pidfilter.pidtoffilter.tpcel, pidfilter.pidtoffilter.tpcmu, pidfilter.pidtoffilter.tpcpi, pidfilter.pidtoffilter.tpcka, pidfilter.pidtoffilter.tpcpr);
#endif
    fFilterFramework = new PWGCF::FilterAndAnalysisFramework(filterccdb.ccdburl.value, filterccdb.ccdbpath.value, filterccdb.filterdate.value);
    fFilterFramework->SetConfiguration(eventsel, trksel, pidsel, PWGCF::SelectionFilterAndAnalysis::kAnalysis);
    fFilterFramework->Init();

    collisionmask = fFilterFramework->getCollisionMask();
    collisionmask_opt = fFilterFramework->getCollisionOptMask();
    collisionmask_forced = fFilterFramework->getCollisionForcedMask();
    fMultiplicityIndex = fFilterFramework->getCollisionMultiplicityIndex();
    trackmask = fFilterFramework->getTrackMask();
    trackmask_opt = fFilterFramework->getTrackOptMask();
    trackmask_forced = fFilterFramework->getTrackForcedMask();
    pidmask = fFilterFramework->getPIDMask();
    pidmask_opt = fFilterFramework->getPIDOptMask();
    pidmask_forced = fFilterFramework->getPIDForcedMask();
    LOGF(info, "twoParticleCorrelationsFilter::init(), collision selection masks 0x%016lx, %s, and 0x%016lx and multiplicity index %d", collisionmask, fFilterFramework->printCollisionOptionalMasks().Data(), collisionmask_forced, fMultiplicityIndex);
    LOGF(info, "twoParticleCorrelationsFilter::init(), track selection masks 0x%016lx, %s, and 0x%016lx ", trackmask, fFilterFramework->printTrackOptionalMasks().Data(), trackmask_forced);
    LOGF(info, "twoParticleCorrelationsFilter::init(), PID selection masks 0x%016lx, %s, and 0x%016lx ", pidmask, fFilterFramework->printPIDOptionalMasks().Data(), pidmask_forced);
    if (collisionmask == static_cast<uint64_t>(0) || trackmask == static_cast<uint64_t>(0)) {
      LOGF(fatal, "twoParticleCorrelationsFilter::init() null masks, selecting everything!!!");
    }

    /* initial size for the tracks ids bank */
    tracksids.reserve(bankinitialcapacity);
  }

  /// \brief Get the data collecting engine index corresponding to the passed collision
  template <typename FilteredCollision>
  int getDCEindex(FilteredCollision collision)
  {
    int ixDCE = -1;
    float cm = collision.centmult()[fMultiplicityIndex];
    for (int i = 0; i < ncmranges; ++i) {
      if (cm < fCentMultMax[i]) {
        ixDCE = i;
        break;
      }
    }
    return ixDCE;
  }

  TList* getCCDBInput(const char* ccdbpath, const char* ccdbdate)
  {
    std::tm cfgtm = {};
    std::stringstream ss(ccdbdate);
    ss >> std::get_time(&cfgtm, "%Y%m%d");
    cfgtm.tm_hour = 12;
    uint64_t timestamp = std::mktime(&cfgtm) * 1000;

    TList* lst = ccdb->getForTimeStamp<TList>(ccdbpath, timestamp);
    if (lst != nullptr) {
      LOGF(info, "Correctly loaded CCDB input object");
    } else {
      LOGF(error, "CCDB input object could not be loaded");
    }
    return lst;
  }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPObject* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <typename CollisionObject, typename TracksObject>
  void processTheTask(CollisionObject const& collision, TracksObject& tracks)
  {
    using namespace twopcorrelations;
    LOGF(TWOPCORRLOGCOLLISIONS, "Received collision with mask 0x%016lx and %ld tracks", collision.selflags(), tracks.size());
    auto passOptions = [](auto options, auto mask) {
      bool all = true;
      for (auto option : options) {
        all = all && ((option & mask) != 0UL);
      }
      return all;
    };

    if ((collision.selflags() & collisionmask_forced) == collisionmask_forced && passOptions(collisionmask_opt, collision.selflags())) {
      LOGF(TWOPCORRLOGCOLLISIONS, ">> Accepted collision with mask 0x%016lx and %ld unfiltered tracks", collision.selflags(), tracks.size());
      int nAcceptedTracks = 0;
      int nRejectedTracks = 0;
      tracksids.clear();
      if (static_cast<int64_t>(tracksids.capacity()) < tracks.size()) {
        tracksids.reserve(tracksids.capacity() + bankcapacitystep);
        LOGF(warning, "Increased tracks bank size to %d", tracksids.size());
      }
      for (const auto& track : tracks) {
        int trackid = -1;
        if ((track.trackflags() & trackmask_forced) == trackmask_forced && passOptions(trackmask_opt, track.trackflags())) {
          if ((track.trackflags() & 0x1L) == 0x1L) {
            /* positive track */
            trackid = 0;
            nAcceptedTracks++;
          } else if ((track.trackflags() & 0x2L) == 0x2L) {
            /* negative track */
            trackid = 1;
            nAcceptedTracks++;
          } else {
            nRejectedTracks++;
          }
        } else {
          nRejectedTracks++;
        }
        if (!(trackid < 0)) {
          tracksids.push_back({trackid, 1.0, 0.0, static_cast<uint64_t>(track.filteredIndex())});
        }
      }
      if (nAcceptedTracks > 0) {
        int ixDCE = getDCEindex(collision);
        if (!(ixDCE < 0)) {
          if (ccdblst != nullptr && !dataCE[ixDCE]->isCCDBstored()) {
            dataCE[ixDCE]->storeTrackCorrections(std::vector<TH3*>{static_cast<TH3*>(ccdblst->FindObject(TString::Format("correction_%02d-%02d_p1", static_cast<int>(fCentMultMin[ixDCE]), static_cast<int>(fCentMultMax[ixDCE])).Data())),
                                                                   static_cast<TH3*>(ccdblst->FindObject(TString::Format("correction_%02d-%02d_m1", static_cast<int>(fCentMultMin[ixDCE]), static_cast<int>(fCentMultMax[ixDCE])).Data()))});
            dataCE[ixDCE]->storePtAverages(std::vector<TH2*>{static_cast<TH2*>(ccdblst->FindObject(TString::Format("ptavgetaphi_%02d-%02d_p", static_cast<int>(fCentMultMin[ixDCE]), static_cast<int>(fCentMultMax[ixDCE])).Data())),
                                                             static_cast<TH2*>(ccdblst->FindObject(TString::Format("ptavgetaphi_%02d-%02d_m", static_cast<int>(fCentMultMin[ixDCE]), static_cast<int>(fCentMultMax[ixDCE])).Data()))});
          }
        }
        /* magnetic field equal to zero for the time being */
        dataCE[ixDCE]->processCollision(tracks, collision.posZ(), collision.centmult()[fMultiplicityIndex], 0);
      }
      LOGF(TWOPCORRLOGCOLLISIONS, ">> Accepted collision with mask 0x%016lx, %d accepted tracks and %d rejected tracks", collision.selflags(), nAcceptedTracks, nRejectedTracks);
    } else {
    }
  }

  Filter tracksfilter = (aod::track::eta < cfgBinning->mEtamax) && (aod::track::eta > cfgBinning->mEtamin) && (aod::track::pt < cfgBinning->mPTmax) && (aod::track::pt > cfgBinning->mPTmin);

  void processRun2(soa::Join<aod::Collisions, aod::CFCollMasks>::iterator const& collision,
                   soa::Filtered<soa::Join<aod::FullTracks, aod::CFTrackMasks>>& tracks)
  {
    processTheTask(collision, tracks);
  }
  PROCESS_SWITCH(twoParticleCorrelations, processRun2, "Process Run 2, i.e. over NOT stored derived data, two particle correlations", true);

  void processRun3(aod::CFCollision const& collision, soa::Filtered<aod::CFTracks>& tracks)
  {
    processTheTask(collision, tracks);
  }
  PROCESS_SWITCH(twoParticleCorrelations, processRun3, "Process Run 3, i.e. over stored derived data, two particle correlations", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<twoParticleCorrelations>(cfgc)};
  return workflow;
}
