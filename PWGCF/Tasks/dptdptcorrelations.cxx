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

#include <CCDB/BasicCCDBManager.h>
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
#include <TVector2.h>
#include <cmath>
#include <ctime>

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/DptDptFiltered.h"
#include "PWGCF/TableProducer/dptdptfilter.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;

#define DPTDPTLOGCOLLISIONS debug
#define DPTDPTLOGTRACKS debug

namespace correlationstask
{
using namespace o2::analysis::dptdptfilter;
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
bool processmixedevents = false;
bool ptorder = false;

PairCuts fPairCuts;              // pair suppression engine
bool fUseConversionCuts = false; // suppress resonances and conversions
bool fUseTwoTrackCut = false;    // suppress too close tracks

std::vector<std::string> tname = {"O", "T"}; ///< the track names
} // namespace correlationstask

// Task for building <dpt,dpt> correlations
struct DptDptCorrelationsTask {

  /* the data collecting engine */
  template <bool smallsingles>
  struct DataCollectingEngine {
    int nspecies = 1; /* for the time being just hadrons */
    size_t nch = nspecies * 2;

    //============================================================================================
    // The DptDptCorrelationsAnalysisTask output objects
    //============================================================================================
    /* histograms */
    TH1F* fhVertexZA;                                                             //!<! the z vertex distribution for the current multiplicity/centrality class
    std::vector<TH1F*> fhN1_vsPt{nch, nullptr};                                   //!<! weighted single particle distribution vs \f$p_T\f$, for the different species
    std::vector<TH2F*> fhN1_vsEtaPhi{nch, nullptr};                               //!<! weighted single particle distribution vs \f$\eta,\;\phi\f$, for the different species
    std::vector<TH2F*> fhSum1Pt_vsEtaPhi{nch, nullptr};                           //!<! accumulated sum of weighted \f$p_T\f$ vs \f$\eta,\;\phi\f$, for the different species
    std::vector<TH3F*> fhN1_vsZEtaPhiPt{nch, nullptr};                            //!<! single particle distribution vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$, for the different species
    std::vector<TH3F*> fhSum1Pt_vsZEtaPhiPt{nch, nullptr};                        //!<! accumulated sum of weighted \f$p_T\f$ vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$, for the different species
    std::vector<TH3*> fhNuaNue_vsZEtaPhiPt{nch, nullptr};                         //!<! NUA+NUE correction vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$, for the differents species
    std::vector<TH2*> fhPtAvg_vsEtaPhi{nch, nullptr};                             //!<! average \f$p_T\f$ vs \f$\eta,\;\phi\f$, for the different species
    std::vector<std::vector<TH2F*>> fhN2_vsPtPt{nch, {nch, nullptr}};             //!<! weighted two particle distribution vs \f${p_T}_1, {p_T}_2\f$ for the different species combinations
    std::vector<std::vector<TH2F*>> fhN2_vsDEtaDPhi{nch, {nch, nullptr}};         //!<! two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$ for the different species combinations
    std::vector<std::vector<TH2F*>> fhN2cont_vsDEtaDPhi{nch, {nch, nullptr}};     //!<! two-particle distribution continuous vs \f$\Delta\eta,\;\Delta\phi\f$ for the different species combinations
    std::vector<std::vector<TH2F*>> fhSum2PtPt_vsDEtaDPhi{nch, {nch, nullptr}};   //!<! two-particle  \f$\sum {p_T}_1 {p_T}_2\f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ for the different species combinations
    std::vector<std::vector<TH2F*>> fhSum2DptDpt_vsDEtaDPhi{nch, {nch, nullptr}}; //!<! two-particle  \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ for the different species combinations
    std::vector<std::vector<TH2F*>> fhSupN1N1_vsDEtaDPhi{nch, {nch, nullptr}};    //!<! suppressed n1n1 two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$ for the different species combinations
    std::vector<std::vector<TH2F*>> fhSupPt1Pt1_vsDEtaDPhi{nch, {nch, nullptr}};  //!<! suppressed \f${p_T}_1 {p_T}_2\f$ two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$ for the different species combinations
    /* versus centrality/multiplicity  profiles */
    std::vector<TProfile*> fhN1_vsC{nch, nullptr};                               //!<! weighted single particle distribution vs event centrality/multiplicity, track 1 and 2
    std::vector<TProfile*> fhSum1Pt_vsC{nch, nullptr};                           //!<! accumulated sum of weighted \f$p_T\f$ vs event centrality/multiplicity, track 1 and 2
    std::vector<TProfile*> fhN1nw_vsC{nch, nullptr};                             //!<! un-weighted single particle distribution vs event centrality/multiplicity, track 1 and 2
    std::vector<TProfile*> fhSum1Ptnw_vsC{nch, nullptr};                         //!<! accumulated sum of un-weighted \f$p_T\f$ vs event centrality/multiplicity, track 1 and 2
    std::vector<std::vector<TProfile*>> fhN2_vsC{nch, {nch, nullptr}};           //!<! weighted accumulated two particle distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    std::vector<std::vector<TProfile*>> fhSum2PtPt_vsC{nch, {nch, nullptr}};     //!<! weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    std::vector<std::vector<TProfile*>> fhSum2DptDpt_vsC{nch, {nch, nullptr}};   //!<! weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    std::vector<std::vector<TProfile*>> fhN2nw_vsC{nch, {nch, nullptr}};         //!<! un-weighted accumulated two particle distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    std::vector<std::vector<TProfile*>> fhSum2PtPtnw_vsC{nch, {nch, nullptr}};   //!<! un-weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    std::vector<std::vector<TProfile*>> fhSum2DptDptnw_vsC{nch, {nch, nullptr}}; //!<! un-weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations

    std::vector<std::vector<std::string>> trackPairsNames = {{"OO", "OT"}, {"TO", "TT"}};
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
      using namespace correlationstask;
      using namespace o2::analysis::dptdptfilter;
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
      using namespace correlationstask;
      using namespace o2::analysis::dptdptfilter;

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
      using namespace correlationstask;
      using namespace o2::analysis::dptdptfilter;

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
      LOGF(info, "Stored NUA&NUE corrections for %d track ids", corrs.size());
      for (uint i = 0; i < corrs.size(); ++i) {
        LOGF(info, "  Stored NUA&NUE corrections %s for track id %d %s", corrs[i] != nullptr ? corrs[i]->GetName() : "nullptr", i, corrs[i] != nullptr ? "yes" : "no");
        fhNuaNue_vsZEtaPhiPt[i] = corrs[i];
        if (fhNuaNue_vsZEtaPhiPt[i] != nullptr) {
          int nbins = 0;
          double avg = 0.0;
          for (int ix = 0; ix < fhNuaNue_vsZEtaPhiPt[i]->GetNbinsX(); ++ix) {
            for (int iy = 0; iy < fhNuaNue_vsZEtaPhiPt[i]->GetNbinsY(); ++iy) {
              for (int iz = 0; iz < fhNuaNue_vsZEtaPhiPt[i]->GetNbinsZ(); ++iz) {
                nbins++;
                avg += fhNuaNue_vsZEtaPhiPt[i]->GetBinContent(ix + 1, iy + 1, iz + 1);
              }
            }
          }
          LOGF(info, "Average NUA&NUE correction for track id %d: %f", i, avg / nbins);
        }
      }
      ccdbstored = true;
    }

    void storePtAverages(std::vector<TH2*> ptavgs)
    {
      LOGF(info, "Stored pT average for %d track ids", ptavgs.size());
      for (uint i = 0; i < ptavgs.size(); ++i) {
        LOGF(info, "  Stored pT average for track id %d %s", i, ptavgs[i] != nullptr ? "yes" : "no");
        fhPtAvg_vsEtaPhi[i] = ptavgs[i];
      }
      ccdbstored = true;
    }

    template <typename TrackListObject>
    std::vector<float>* getTrackCorrections(TrackListObject const& tracks, float zvtx)
    {
      std::vector<float>* corr = new std::vector<float>(tracks.size(), 1.0f);
      int index = 0;
      for (auto& t : tracks) {
        if (fhNuaNue_vsZEtaPhiPt[t.trackacceptedid()] != nullptr) {
          (*corr)[index] = fhNuaNue_vsZEtaPhiPt[t.trackacceptedid()]->GetBinContent(fhNuaNue_vsZEtaPhiPt[t.trackacceptedid()]->FindFixBin(zvtx, GetEtaPhiIndex(t) + 0.5, t.pt()));
        }
        index++;
      }
      return corr;
    }

    template <typename TrackListObject>
    std::vector<float>* getPtAvg(TrackListObject const& tracks)
    {
      std::vector<float>* ptavg = new std::vector<float>(tracks.size(), 0.0f);
      int index = 0;
      for (auto& t : tracks) {
        if (fhPtAvg_vsEtaPhi[t.trackacceptedid()] != nullptr) {
          (*ptavg)[index] = fhPtAvg_vsEtaPhi[t.trackacceptedid()]->GetBinContent(fhPtAvg_vsEtaPhi[t.trackacceptedid()]->FindFixBin(t.eta(), t.phi()));
          index++;
        }
      }
      return ptavg;
    }

    /// \brief fills the singles histograms in singles execution mode
    /// \param passedtracks filtered table with the tracks associated to the passed index
    /// \param tix index, in the singles histogram bank, for the passed filetered track table
    template <typename TrackListObject>
    void processSingles(TrackListObject const& passedtracks, std::vector<float>* corrs, float zvtx)
    {
      int index = 0;
      for (auto& track : passedtracks) {
        float corr = (*corrs)[index];
        fhN1_vsPt[track.trackacceptedid()]->Fill(track.pt(), corr);
        if constexpr (smallsingles) {
          fhN1_vsEtaPhi[track.trackacceptedid()]->Fill(track.eta(), GetShiftedPhi(track.phi()), corr);
          fhSum1Pt_vsEtaPhi[track.trackacceptedid()]->Fill(track.eta(), GetShiftedPhi(track.phi()), track.pt() * corr);
        } else {
          fhN1_vsZEtaPhiPt[track.trackacceptedid()]->Fill(zvtx, GetEtaPhiIndex(track) + 0.5, track.pt(), corr);
          fhSum1Pt_vsZEtaPhiPt[track.trackacceptedid()]->Fill(zvtx, GetEtaPhiIndex(track) + 0.5, track.pt(), track.pt() * corr);
        }
        index++;
      }
    }

    /// \brief fills the singles histograms in pair execution mode
    /// \param passedtracks filtered table with the tracks associated to the passed index
    /// \param tix index, in the singles histogram bank, for the passed filetered track table
    /// \param cmul centrality - multiplicity for the collision being analyzed
    template <typename TrackListObject>
    void processTracks(TrackListObject const& passedtracks, std::vector<float>* corrs, float cmul)
    {
      LOGF(DPTDPTLOGCOLLISIONS, "Processing %d tracks in a collision with cent/mult %f ", passedtracks.size(), cmul);

      /* process magnitudes */
      std::vector<double> n1(nch, 0.0);       ///< weighted number of single tracks for current collision
      std::vector<double> sum1Pt(nch, 0.0);   ///< accumulated sum of weighted single track \f$p_T\f$ for current collision
      std::vector<double> n1nw(nch, 0.0);     ///< not weighted number of single tracks for current collision
      std::vector<double> sum1Ptnw(nch, 0.0); ///< accumulated sum of not weighted single \f$p_T\f$ for current collision
      int index = 0;
      for (auto& track : passedtracks) {
        float corr = (*corrs)[index];
        n1[track.trackacceptedid()] += corr;
        sum1Pt[track.trackacceptedid()] += track.pt() * corr;
        n1nw[track.trackacceptedid()] += 1;
        sum1Ptnw[track.trackacceptedid()] += track.pt();

        fhN1_vsEtaPhi[track.trackacceptedid()]->Fill(track.eta(), GetShiftedPhi(track.phi()), corr);
        fhSum1Pt_vsEtaPhi[track.trackacceptedid()]->Fill(track.eta(), GetShiftedPhi(track.phi()), track.pt() * corr);
        index++;
      }
      for (uint tid = 0; tid < nch; ++tid) {
        fhN1_vsC[tid]->Fill(cmul, n1[tid]);
        fhSum1Pt_vsC[tid]->Fill(cmul, sum1Pt[tid]);
        fhN1nw_vsC[tid]->Fill(cmul, n1nw[tid]);
        fhSum1Ptnw_vsC[tid]->Fill(cmul, sum1Ptnw[tid]);
      }
    }

    /// \brief fills the pair histograms in pair execution mode
    /// \param trks1 filtered table with the tracks associated to the first track in the pair
    /// \param trks2 filtered table with the tracks associated to the second track in the pair
    /// \param cmul centrality - multiplicity for the collision being analyzed
    /// Be aware that in most of the cases traks1 and trks2 will have the same content (exception: mixed events)
    template <bool doptorder, typename TrackOneListObject, typename TrackTwoListObject>
    void processTrackPairs(TrackOneListObject const& trks1, TrackTwoListObject const& trks2, std::vector<float>* corrs1, std::vector<float>* corrs2, std::vector<float>* ptavgs1, std::vector<float>* ptavgs2, float cmul, int bfield)
    {
      using namespace correlationstask;

      /* process pair magnitudes */
      std::vector<std::vector<double>> n2(nch, std::vector<double>(nch, 0.0));           ///< weighted number of track 1 track 2 pairs for current collision
      std::vector<std::vector<double>> n2sup(nch, std::vector<double>(nch, 0.0));        ///< weighted number of track 1 track 2 suppressed pairs for current collision
      std::vector<std::vector<double>> sum2PtPt(nch, std::vector<double>(nch, 0.0));     ///< accumulated sum of weighted track 1 track 2 \f${p_T}_1 {p_T}_2\f$ for current collision
      std::vector<std::vector<double>> sum2DptDpt(nch, std::vector<double>(nch, 0.0));   ///< accumulated sum of weighted number of track 1 tracks times weighted track 2 \f$p_T\f$ for current collision
      std::vector<std::vector<double>> n2nw(nch, std::vector<double>(nch, 0.0));         ///< not weighted number of track1 track 2 pairs for current collision
      std::vector<std::vector<double>> sum2PtPtnw(nch, std::vector<double>(nch, 0.0));   ///< accumulated sum of not weighted track 1 track 2 \f${p_T}_1 {p_T}_2\f$ for current collision
      std::vector<std::vector<double>> sum2DptDptnw(nch, std::vector<double>(nch, 0.0)); ///< accumulated sum of not weighted number of track 1 tracks times not weighted track 2 \f$p_T\f$ for current collision
      int index1 = 0;

      for (auto& track1 : trks1) {
        double ptavg_1 = (*ptavgs1)[index1];
        double corr1 = (*corrs1)[index1];
        int index2 = 0;
        for (auto& track2 : trks2) {
          /* checking the same track id condition */
          if (track1 == track2) {
            /* exclude autocorrelations */
            continue;
          }

          if constexpr (doptorder) {
            if (track2.pt() >= track1.pt()) {
              continue;
            }
          }
          /* process pair magnitudes */
          double ptavg_2 = (*ptavgs2)[index2];
          double corr2 = (*corrs2)[index2];
          double corr = corr1 * corr2;
          double dptdptnw = (track1.pt() - ptavg_1) * (track2.pt() - ptavg_2);
          double dptdptw = (corr1 * track1.pt() - ptavg_1) * (corr2 * track2.pt() - ptavg_2);

          /* get the global bin for filling the differential histograms */
          int globalbin = GetDEtaDPhiGlobalIndex(track1, track2);
          float deltaeta = track1.eta() - track2.eta();
          float deltaphi = track1.phi() - track2.phi();
          while (deltaphi >= deltaphiup) {
            deltaphi -= constants::math::TwoPI;
          }
          while (deltaphi < deltaphilow) {
            deltaphi += constants::math::TwoPI;
          }
          if ((fUseConversionCuts && fPairCuts.conversionCuts(track1, track2)) || (fUseTwoTrackCut && fPairCuts.twoTrackCut(track1, track2, bfield))) {
            /* suppress the pair */
            fhSupN1N1_vsDEtaDPhi[track1.trackacceptedid()][track2.trackacceptedid()]->AddBinContent(globalbin, corr);
            fhSupPt1Pt1_vsDEtaDPhi[track1.trackacceptedid()][track2.trackacceptedid()]->AddBinContent(globalbin, track1.pt() * track2.pt() * corr);
            n2sup[track1.trackacceptedid()][track2.trackacceptedid()] += corr;
          } else {
            /* count the pair */
            n2[track1.trackacceptedid()][track2.trackacceptedid()] += corr;
            sum2PtPt[track1.trackacceptedid()][track2.trackacceptedid()] += track1.pt() * track2.pt() * corr;
            sum2DptDpt[track1.trackacceptedid()][track2.trackacceptedid()] += dptdptw;
            n2nw[track1.trackacceptedid()][track2.trackacceptedid()] += 1;
            sum2PtPtnw[track1.trackacceptedid()][track2.trackacceptedid()] += track1.pt() * track2.pt();
            sum2DptDptnw[track1.trackacceptedid()][track2.trackacceptedid()] += dptdptnw;

            fhN2_vsDEtaDPhi[track1.trackacceptedid()][track2.trackacceptedid()]->AddBinContent(globalbin, corr);
            fhN2cont_vsDEtaDPhi[track1.trackacceptedid()][track2.trackacceptedid()]->Fill(deltaeta, deltaphi, corr);
            fhSum2DptDpt_vsDEtaDPhi[track1.trackacceptedid()][track2.trackacceptedid()]->AddBinContent(globalbin, dptdptw);
            fhSum2PtPt_vsDEtaDPhi[track1.trackacceptedid()][track2.trackacceptedid()]->AddBinContent(globalbin, track1.pt() * track2.pt() * corr);
          }
          fhN2_vsPtPt[track1.trackacceptedid()][track2.trackacceptedid()]->Fill(track1.pt(), track2.pt(), corr);
          index2++;
        }
        index1++;
      }
      for (uint pid1 = 0; pid1 < nch; ++pid1) {
        for (uint pid2 = 0; pid2 < nch; ++pid2) {
          fhN2_vsC[pid1][pid2]->Fill(cmul, n2[pid1][pid2]);
          fhSum2PtPt_vsC[pid1][pid2]->Fill(cmul, sum2PtPt[pid1][pid2]);
          fhSum2DptDpt_vsC[pid1][pid2]->Fill(cmul, sum2DptDpt[pid1][pid2]);
          fhN2nw_vsC[pid1][pid2]->Fill(cmul, n2nw[pid1][pid2]);
          fhSum2PtPtnw_vsC[pid1][pid2]->Fill(cmul, sum2PtPtnw[pid1][pid2]);
          fhSum2DptDptnw_vsC[pid1][pid2]->Fill(cmul, sum2DptDptnw[pid1][pid2]);
          /* let's also update the number of entries in the differential histograms */
          fhN2_vsDEtaDPhi[pid1][pid2]->SetEntries(fhN2_vsDEtaDPhi[pid1][pid2]->GetEntries() + n2[pid1][pid2]);
          fhSum2DptDpt_vsDEtaDPhi[pid1][pid2]->SetEntries(fhSum2DptDpt_vsDEtaDPhi[pid1][pid2]->GetEntries() + n2[pid1][pid2]);
          fhSum2PtPt_vsDEtaDPhi[pid1][pid2]->SetEntries(fhSum2PtPt_vsDEtaDPhi[pid1][pid2]->GetEntries() + n2[pid1][pid2]);
          fhSupN1N1_vsDEtaDPhi[pid1][pid2]->SetEntries(fhSupN1N1_vsDEtaDPhi[pid1][pid2]->GetEntries() + n2sup[pid1][pid2]);
          fhSupPt1Pt1_vsDEtaDPhi[pid1][pid2]->SetEntries(fhSupPt1Pt1_vsDEtaDPhi[pid1][pid2]->GetEntries() + n2sup[pid1][pid2]);
        }
      }
    }

    template <bool mixed, typename TrackOneListObject, typename TrackTwoListObject>
    void processCollision(TrackOneListObject const& Tracks1, TrackTwoListObject const& Tracks2, float zvtx, float centmult, int bfield)
    {
      using namespace correlationstask;
      std::vector<float>* corrs1;
      std::vector<float>* corrs2;
      corrs1 = getTrackCorrections(Tracks1, zvtx);
      if constexpr (mixed) {
        corrs2 = getTrackCorrections(Tracks2, zvtx);
      }

      if (!processpairs) {
        /* process single tracks */
        fhVertexZA->Fill(zvtx);
        processSingles(Tracks1, corrs1, zvtx);
        if constexpr (mixed) {
          processSingles(Tracks2, corrs2, zvtx);
        }
      } else {
        /* process track magnitudes */
        std::vector<float>* ptavgs1;
        std::vector<float>* ptavgs2;
        ptavgs1 = getPtAvg(Tracks1);
        if constexpr (mixed) {
          ptavgs2 = getPtAvg(Tracks2);
        }

        /* TODO: the centrality should be chosen non detector dependent */
        processTracks(Tracks1, corrs1, centmult);
        if constexpr (mixed) {
          processTracks(Tracks2, corrs2, centmult);
        }
        /* process pair magnitudes */
        if constexpr (mixed) {
          if (ptorder) {
            processTrackPairs<true>(Tracks1, Tracks2, corrs1, corrs2, ptavgs1, ptavgs2, centmult, bfield);
          } else {
            processTrackPairs<false>(Tracks1, Tracks2, corrs1, corrs2, ptavgs1, ptavgs2, centmult, bfield);
          }
        } else {
          if (ptorder) {
            processTrackPairs<true>(Tracks1, Tracks1, corrs1, corrs1, ptavgs1, ptavgs1, centmult, bfield);
          } else {
            processTrackPairs<false>(Tracks1, Tracks1, corrs1, corrs1, ptavgs1, ptavgs1, centmult, bfield);
          }
        }

        delete ptavgs1;
        if constexpr (mixed) {
          delete ptavgs2;
        }
      }
      delete corrs1;
      if constexpr (mixed) {
        delete corrs2;
      }
    }

    void init(TList* fOutputList)
    {
      using namespace correlationstask;
      using namespace o2::analysis::dptdptfilter;

      /* create the histograms */
      Bool_t oldstatus = TH1::AddDirectoryStatus();
      TH1::AddDirectory(kFALSE);

      if (!processpairs) {
        fhVertexZA = new TH1F("VertexZA", "Vertex Z; z_{vtx}", zvtxbins, zvtxlow, zvtxup);
        fOutputList->Add(fhVertexZA);
        for (uint i = 0; i < nch; ++i) {
          /* histograms for each track, one and two */
          fhN1_vsPt[i] = new TH1F(TString::Format("n1_%s_vsPt", tname[i].c_str()).Data(),
                                  TString::Format("#LT n_{1} #GT;p_{t,%s} (GeV/c);#LT n_{1} #GT", tname[i].c_str()).Data(),
                                  ptbins, ptlow, ptup);
          /* we don't want the Sumw2 structure being created here */
          bool defSumw2 = TH1::GetDefaultSumw2();
          if constexpr (smallsingles) {
            fhN1_vsEtaPhi[i] = new TH2F(TString::Format("n1_%s_vsEtaPhi", tname[i].c_str()).Data(),
                                        TString::Format("#LT n_{1} #GT;#eta_{%s};#varphi_{%s} (radian);#LT n_{1} #GT", tname[i].c_str(), tname[i].c_str()).Data(),
                                        etabins, etalow, etaup, phibins, philow, phiup);
            fhSum1Pt_vsEtaPhi[i] = new TH2F(TString::Format("sumPt_%s_vsEtaPhi", tname[i].c_str()).Data(),
                                            TString::Format("#LT #Sigma p_{t,%s} #GT;#eta_{%s};#varphi_{%s} (radian);#LT #Sigma p_{t,%s} #GT (GeV/c)",
                                                            tname[i].c_str(), tname[i].c_str(), tname[i].c_str(), tname[i].c_str())
                                              .Data(),
                                            etabins, etalow, etaup, phibins, philow, phiup);
          } else {
            TH1::SetDefaultSumw2(false);
            fhN1_vsZEtaPhiPt[i] = new TH3F(
              TString::Format("n1_%s_vsZ_vsEtaPhi_vsPt", tname[i].c_str()).Data(),
              TString::Format("#LT n_{1} #GT;vtx_{z};#eta_{%s}#times#varphi_{%s};p_{t,%s} (GeV/c)",
                              tname[i].c_str(),
                              tname[i].c_str(),
                              tname[i].c_str())
                .Data(),
              zvtxbins,
              zvtxlow,
              zvtxup,
              etabins * phibins,
              0.0,
              static_cast<double>(etabins * phibins),
              ptbins,
              ptlow,
              ptup);
            fhSum1Pt_vsZEtaPhiPt[i] = new TH3F(
              TString::Format("sumPt1_%s_vsZ_vsEtaPhi_vsPt", tname[i].c_str()).Data(),
              TString::Format(
                "#LT #Sigma p_{t,%s}#GT;vtx_{z};#eta_{%s}#times#varphi_{%s};p_{t,%s} (GeV/c)",
                tname[i].c_str(),
                tname[i].c_str(),
                tname[i].c_str(),
                tname[i].c_str())
                .Data(),
              zvtxbins,
              zvtxlow,
              zvtxup,
              etabins * phibins,
              0.0,
              static_cast<double>(etabins * phibins),
              ptbins,
              ptlow,
              ptup);
          }
          /* we return it back to previuos state */
          TH1::SetDefaultSumw2(defSumw2);

          /* the statistical uncertainties will be estimated by the subsamples method so let's get rid of the error tracking */
          if constexpr (smallsingles) {
            fhN1_vsEtaPhi[i]->SetBit(TH1::kIsNotW);
            fhN1_vsEtaPhi[i]->Sumw2(false);
            fhSum1Pt_vsEtaPhi[i]->SetBit(TH1::kIsNotW);
            fhSum1Pt_vsEtaPhi[i]->Sumw2(false);
          } else {
            fhN1_vsZEtaPhiPt[i]->SetBit(TH1::kIsNotW);
            fhN1_vsZEtaPhiPt[i]->Sumw2(false);
            fhSum1Pt_vsZEtaPhiPt[i]->SetBit(TH1::kIsNotW);
            fhSum1Pt_vsZEtaPhiPt[i]->Sumw2(false);
          }
          fhNuaNue_vsZEtaPhiPt[i] = nullptr;
          fhPtAvg_vsEtaPhi[i] = nullptr;

          fOutputList->Add(fhN1_vsPt[i]);
          if constexpr (smallsingles) {
            fOutputList->Add(fhN1_vsEtaPhi[i]);
            fOutputList->Add(fhSum1Pt_vsEtaPhi[i]);
          } else {
            fOutputList->Add(fhN1_vsZEtaPhiPt[i]);
            fOutputList->Add(fhSum1Pt_vsZEtaPhiPt[i]);
          }
        }
      } else {
        for (uint i = 0; i < nch; ++i) {
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
        }

        for (uint i = 0; i < nch; ++i) {
          for (uint j = 0; j < nch; ++j) {
            /* histograms for each track pair combination */
            /* we don't want the Sumw2 structure being created here */
            bool defSumw2 = TH1::GetDefaultSumw2();
            TH1::SetDefaultSumw2(false);
            const char* pname = trackPairsNames[i][j].c_str();
            fhN2_vsDEtaDPhi[i][j] = new TH2F(TString::Format("n2_12_vsDEtaDPhi_%s", pname), TString::Format("#LT n_{2} #GT (%s);#Delta#eta;#Delta#varphi;#LT n_{2} #GT", pname),
                                             deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
            fhN2cont_vsDEtaDPhi[i][j] = new TH2F(TString::Format("n2_12cont_vsDEtaDPhi_%s", pname), TString::Format("#LT n_{2} #GT (%s);#Delta#eta;#Delta#varphi;#LT n_{2} #GT", pname),
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
            fhN2cont_vsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhN2cont_vsDEtaDPhi[i][j]->Sumw2(false);
            fhSum2PtPt_vsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhSum2PtPt_vsDEtaDPhi[i][j]->Sumw2(false);
            fhSum2DptDpt_vsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhSum2DptDpt_vsDEtaDPhi[i][j]->Sumw2(false);
            fhSupN1N1_vsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhSupN1N1_vsDEtaDPhi[i][j]->Sumw2(false);
            fhSupPt1Pt1_vsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhSupPt1Pt1_vsDEtaDPhi[i][j]->Sumw2(false);

            fOutputList->Add(fhN2_vsDEtaDPhi[i][j]);
            fOutputList->Add(fhN2cont_vsDEtaDPhi[i][j]);
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

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  /* the data memebers for this task */
  /* the centrality / multiplicity limits for collecting data in this task instance */
  int ncmranges = 0;
  float* fCentMultMin = nullptr;
  float* fCentMultMax = nullptr;

  /* the data collecting engine instances */
  DataCollectingEngine<false>** dataCE;
  DataCollectingEngine<true>** dataCE_small;
  DataCollectingEngine<false>** dataCEME;

  /* the input file structure from CCDB */
  TList* ccdblst = nullptr;
  bool loadfromccdb = false;

  /* pair conversion suppression defaults */
  static constexpr float cfgPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};
  Configurable<LabeledArray<float>> cfgPairCut{"paircut", {cfgPairCutDefaults[0], 5, {"Photon", "K0", "Lambda", "Phi", "Rho"}}, "Conversion suppressions"};
  /* two tracks cut */
  Configurable<float> cfgTwoTrackCut{"twotrackcut", -1, "Two-tracks cut: -1 = off; >0 otherwise distance value (suggested: 0.02"};
  Configurable<float> cfgTwoTrackCutMinRadius{"twotrackcutminradius", 0.8f, "Two-tracks cut: radius in m from which two-tracks cut is applied"};

  Configurable<bool> cfgSmallDCE{"smalldce", true, "Use small data collecting engine for singles processing, true = yes. Default = true"};
  Configurable<bool> cfgProcessPairs{"processpairs", false, "Process pairs: false = no, just singles, true = yes, process pairs"};
  Configurable<bool> cfgProcessME{"processmixedevents", false, "Process mixed events: false = no, just same event, true = yes, also process mixed events"};
  Configurable<std::string> cfgCentSpec{"centralities", "00-05,05-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80", "Centrality/multiplicity ranges in min-max separated by commas"};

  Configurable<o2::analysis::DptDptBinningCuts> cfgBinning{"binning",
                                                           {28, -7.0, 7.0, 18, 0.2, 2.0, 16, -0.8, 0.8, 72, 0.5},
                                                           "triplets - nbins, min, max - for z_vtx, pT, eta and phi, binning plus bin fraction of phi origin shift"};
  Configurable<bool> cfgPtOrder{"ptorder", false, "enforce pT_1 < pT_2. Defalut: false"};
  struct : ConfigurableGroup {
    Configurable<std::string> cfgCCDBUrl{"input_ccdburl", "http://ccdb-test.cern.ch:8080", "The CCDB url for the input file"};
    Configurable<std::string> cfgCCDBPathName{"input_ccdbpath", "", "The CCDB path for the input file. Default \"\", i.e. don't load from CCDB"};
    Configurable<std::string> cfgCCDBDate{"input_ccdbdate", "20220307", "The CCDB date for the input file"};
  } cfginputfile;

  OutputObj<TList> fOutput{"DptDptCorrelationsData", OutputObjHandlingPolicy::AnalysisObject, OutputObjSourceType::OutputObjSource};

  void init(InitContext const&)
  {
    using namespace correlationstask;
    using namespace o2::analysis::dptdptfilter;

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
    processmixedevents = cfgProcessME.value;
    ptorder = cfgPtOrder.value;
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
    fGlobalOutputList->Add(new TParameter<int>("NoBinsPt", ptbins, 'f'));
    fGlobalOutputList->Add(new TParameter<int>("NoBinsEta", etabins, 'f'));
    fGlobalOutputList->Add(new TParameter<int>("NoBinsPhi", phibins, 'f'));
    fGlobalOutputList->Add(new TParameter<int>("NoBinsVertexZ", zvtxbins, 'f'));
    fGlobalOutputList->Add(new TParameter<float>("MinVertexZ", zvtxlow, 'f'));
    fGlobalOutputList->Add(new TParameter<float>("MaxVertexZ", zvtxup, 'f'));
    fGlobalOutputList->Add(new TParameter<float>("MinPt", ptlow, 'f'));
    fGlobalOutputList->Add(new TParameter<float>("MaxPt", ptup, 'f'));
    fGlobalOutputList->Add(new TParameter<float>("MinEta", etalow, 'f'));
    fGlobalOutputList->Add(new TParameter<float>("MaxEta", etaup, 'f'));
    fGlobalOutputList->Add(new TParameter<float>("MinPhi", philow, 'f'));
    fGlobalOutputList->Add(new TParameter<float>("MaxPhi", phiup, 'f'));
    fGlobalOutputList->Add(new TParameter<float>("PhiBinShift", phibinshift, 'f'));
    fGlobalOutputList->Add(new TParameter<float>("DifferentialOutput", true, 'f'));
    fGlobalOutputList->Add(new TParameter<bool>("SmallDCE", cfgSmallDCE.value, 'f'));

    /* after the parameters dump the proper phi limits are set according to the phi shift */
    phiup = phiup - phibinwidth * phibinshift;
    philow = philow - phibinwidth * phibinshift;

    /* create the data collecting engine instances according to the configured centrality/multiplicity ranges */
    {
      TObjArray* tokens = TString(cfgCentSpec.value.c_str()).Tokenize(",");
      ncmranges = tokens->GetEntries();
      fCentMultMin = new float[ncmranges];
      fCentMultMax = new float[ncmranges];
      dataCE = new DataCollectingEngine<false>*[ncmranges];
      if (cfgSmallDCE) {
        dataCE_small = new DataCollectingEngine<true>*[ncmranges];
      } else {
        dataCE = new DataCollectingEngine<false>*[ncmranges];
      }
      if (processmixedevents) {
        dataCEME = new DataCollectingEngine<false>*[ncmranges];
      }

      for (int i = 0; i < ncmranges; ++i) {
        auto initializeCEInstance = [&fGlobalOutputList](auto dce, auto name) {
          /* crete the output list for the passed centrality/multiplicity range */
          TList* fOutputList = new TList();
          fOutputList->SetName(name);
          fOutputList->SetOwner(true);
          /* init the data collection instance */
          dce->init(fOutputList);
          fGlobalOutputList->Add(fOutputList);
        };
        auto builSmallDCEInstance = [&initializeCEInstance](auto rg, bool me = false) {
          DataCollectingEngine<true>* dce = new DataCollectingEngine<true>();
          initializeCEInstance(dce, TString::Format("DptDptCorrelationsData%s-%s", me ? "ME" : "", rg));
          return dce;
        };
        auto buildCEInstance = [&initializeCEInstance](auto rg, bool me = false) {
          DataCollectingEngine<false>* dce = new DataCollectingEngine<false>();
          initializeCEInstance(dce, TString::Format("DptDptCorrelationsData%s-%s", me ? "ME" : "", rg));
          return dce;
        };
        float cmmin = 0.0f;
        float cmmax = 0.0f;
        sscanf(tokens->At(i)->GetName(), "%f-%f", &cmmin, &cmmax);
        fCentMultMin[i] = cmmin;
        fCentMultMax[i] = cmmax;
        if (cfgSmallDCE.value) {
          if (processpairs) {
            LOGF(fatal, "Processing pairs cannot be used with the small DCE, please configure properly!!");
          }
          dataCE_small[i] = builSmallDCEInstance(tokens->At(i)->GetName());
        } else {
          dataCE[i] = buildCEInstance(tokens->At(i)->GetName());
        }
        if (processmixedevents) {
          /* consistency check */
          if (cfgSmallDCE.value) {
            LOGF(fatal, "Mixed events cannot be used with the small DCE, please configure properly!!");
          }
          dataCEME[i] = buildCEInstance(tokens->At(i)->GetName(), true);
        }
      }
      delete tokens;
      for (int i = 0; i < ncmranges; ++i) {
        LOGF(info, " centrality/multipliicty range: %d, low limit: %f, up limit: %f", i, fCentMultMin[i], fCentMultMax[i]);
      }
    }
    /* two-track cut and conversion suppression */
    fPairCuts.SetHistogramRegistry(nullptr); // not histogram registry for the time being, incompatible with TList when it is empty
    if (processpairs && ((cfgPairCut->get("Photon") > 0) || (cfgPairCut->get("K0") > 0) || (cfgPairCut->get("Lambda") > 0) || (cfgPairCut->get("Phi") > 0) || (cfgPairCut->get("Rho") > 0))) {
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
  }

  /// \brief Get the data collecting engine index corresponding to the passed collision
  template <typename FilteredCollision>
  int getDCEindex(FilteredCollision collision)
  {
    int ixDCE = -1;
    float cm = collision.centmult();
    for (int i = 0; i < ncmranges; ++i) {
      if (cm < fCentMultMax[i]) {
        ixDCE = i;
        break;
      }
    }
    if (!(ixDCE < 0)) {
      if (cm < fCentMultMin[ixDCE]) {
        ixDCE = -1;
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
    int64_t timestamp = std::mktime(&cfgtm) * 1000;

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

  template <bool gen, typename FilterdCollision, typename FilteredTracks>
  void processSame(FilterdCollision const& collision, FilteredTracks const& tracks, uint64_t timestamp = 0)
  {
    using namespace correlationstask;

    if (ccdblst == nullptr) {
      if (loadfromccdb) {
        ccdblst = getCCDBInput(cfginputfile.cfgCCDBPathName->c_str(), cfginputfile.cfgCCDBDate->c_str());
      }
    }

    /* locate the data collecting engine for the collision centrality/multiplicity */
    int ixDCE = getDCEindex(collision);
    if (!(ixDCE < 0)) {
      if (ccdblst != nullptr && !(dataCE[ixDCE]->isCCDBstored())) {
        if constexpr (gen) {
          std::vector<TH2*> ptavgs{tname.size(), nullptr};
          for (uint isp = 0; isp < tname.size(); ++isp) {
            ptavgs[isp] = reinterpret_cast<TH2*>(ccdblst->FindObject(
              TString::Format("trueptavgetaphi_%02d-%02d_%s",
                              static_cast<int>(fCentMultMin[ixDCE]),
                              static_cast<int>(fCentMultMax[ixDCE]),
                              tname[isp].c_str())
                .Data()));
          }
          if (cfgSmallDCE.value) {
            dataCE_small[ixDCE]->storePtAverages(ptavgs);
          } else {
            dataCE[ixDCE]->storePtAverages(ptavgs);
          }
        } else {
          std::vector<TH3*> corrs{tname.size(), nullptr};
          for (uint isp = 0; isp < tname.size(); ++isp) {
            corrs[isp] = reinterpret_cast<TH3*>(ccdblst->FindObject(
              TString::Format("correction_%02d-%02d_%s",
                              static_cast<int>(fCentMultMin[ixDCE]),
                              static_cast<int>(fCentMultMax[ixDCE]),
                              tname[isp].c_str())
                .Data()));
          }
          if (cfgSmallDCE.value) {
            dataCE_small[ixDCE]->storeTrackCorrections(corrs);
          } else {
            dataCE[ixDCE]->storeTrackCorrections(corrs);
          }

          std::vector<TH2*> ptavgs{tname.size(), nullptr};
          for (uint isp = 0; isp < tname.size(); ++isp) {
            ptavgs[isp] = reinterpret_cast<TH2*>(ccdblst->FindObject(
              TString::Format("ptavgetaphi_%02d-%02d_%s",
                              static_cast<int>(fCentMultMin[ixDCE]),
                              static_cast<int>(fCentMultMax[ixDCE]),
                              tname[isp].c_str())
                .Data()));
          }
          if (cfgSmallDCE.value) {
            dataCE_small[ixDCE]->storePtAverages(ptavgs);
          } else {
            dataCE[ixDCE]->storePtAverages(ptavgs);
          }
        }
      }

      std::string generated = "";
      if constexpr (gen) {
        generated = "generated ";
      }

      LOGF(DPTDPTLOGCOLLISIONS,
           "Accepted BC id %d %scollision with cent/mult %f and %d total tracks. Assigned DCE: %d",
           collision.bcId(),
           generated.c_str(),
           collision.centmult(),
           tracks.size(),
           ixDCE);
      int bfield = 0;
      if constexpr (!gen) {
        bfield = (fUseConversionCuts || fUseTwoTrackCut) ? getMagneticField(timestamp) : 0;
      }
      if (cfgSmallDCE.value) {
        dataCE_small[ixDCE]->processCollision<false>(tracks, tracks, collision.posZ(), collision.centmult(), bfield);
      } else {
        dataCE[ixDCE]->processCollision<false>(tracks, tracks, collision.posZ(), collision.centmult(), bfield);
      }
    }
  }

  template <bool gen, typename FilterdCollision, typename FilteredTracks1, typename FilteredTracks2>
  void processMixed(FilterdCollision const& collision, FilteredTracks1 const& tracks1, FilteredTracks2 const& tracks2, uint64_t timestamp = 0)
  {
    using namespace correlationstask;

    if (ccdblst == nullptr) {
      if (loadfromccdb) {
        ccdblst = getCCDBInput(cfginputfile.cfgCCDBPathName->c_str(), cfginputfile.cfgCCDBDate->c_str());
      }
    }

    /* locate the data collecting engine for the collision centrality/multiplicity */
    int ixDCE = getDCEindex(collision);
    if (!(ixDCE < 0)) {
      if (ccdblst != nullptr && !(dataCEME[ixDCE]->isCCDBstored())) {
        if constexpr (gen) {
          std::vector<TH2*> ptavgs{tname.size(), nullptr};
          for (uint isp = 0; isp < tname.size(); ++isp) {
            ptavgs[isp] = reinterpret_cast<TH2*>(ccdblst->FindObject(
              TString::Format("trueptavgetaphi_%02d-%02d_%s",
                              static_cast<int>(fCentMultMin[ixDCE]),
                              static_cast<int>(fCentMultMax[ixDCE]),
                              tname[isp].c_str())
                .Data()));
          }
          dataCEME[ixDCE]->storePtAverages(ptavgs);
        } else {
          std::vector<TH3*> corrs{tname.size(), nullptr};
          for (uint isp = 0; isp < tname.size(); ++isp) {
            corrs[isp] = reinterpret_cast<TH3*>(ccdblst->FindObject(
              TString::Format("correction_%02d-%02d_%s",
                              static_cast<int>(fCentMultMin[ixDCE]),
                              static_cast<int>(fCentMultMax[ixDCE]),
                              tname[isp].c_str())
                .Data()));
          }
          dataCEME[ixDCE]->storeTrackCorrections(corrs);

          std::vector<TH2*> ptavgs{tname.size(), nullptr};
          for (uint isp = 0; isp < tname.size(); ++isp) {
            ptavgs[isp] = reinterpret_cast<TH2*>(ccdblst->FindObject(
              TString::Format("ptavgetaphi_%02d-%02d_%s",
                              static_cast<int>(fCentMultMin[ixDCE]),
                              static_cast<int>(fCentMultMax[ixDCE]),
                              tname[isp].c_str())
                .Data()));
          }
          dataCEME[ixDCE]->storePtAverages(ptavgs);
        }
      }

      std::string generated = "";
      if constexpr (gen) {
        generated = "generated ";
      }

      LOGF(DPTDPTLOGCOLLISIONS,
           "Accepted BC id %d mixed %scollision with cent/mult %f and %d-%d total mixed tracks. "
           "Assigned DCE: %d",
           collision.bcId(),
           generated.c_str(),
           collision.centmult(),
           tracks1.size(),
           tracks2.size(),
           ixDCE);
      int bfield = 0;
      if constexpr (!gen) {
        bfield = (fUseConversionCuts || fUseTwoTrackCut) ? getMagneticField(timestamp) : 0;
      }
      dataCEME[ixDCE]->processCollision<true>(tracks1, tracks2, collision.posZ(), collision.centmult(), bfield);
    }
  }

  Filter onlyacceptedcollisions = (aod::dptdptfilter::collisionaccepted == uint8_t(true));
  Filter onlyacceptedtracks = (aod::dptdptfilter::trackacceptedid >= int8_t(0));

  void processRecLevel(soa::Filtered<aod::DptDptCFAcceptedCollisions>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<aod::ScannedTracks>& tracks)
  {
    processSame<false>(collision, tracks, collision.bc_as<aod::BCsWithTimestamps>().timestamp());
  }
  PROCESS_SWITCH(DptDptCorrelationsTask, processRecLevel, "Process reco level correlations", false);

  void processRecLevelCheck(aod::Collisions const& collisions, aod::Tracks& tracks)
  {
    int nAssignedTracks = 0;
    int nNotAssignedTracks = 0;
    int64_t firstNotAssignedIndex = -1;
    int64_t lastNotAssignedIndex = -1;

    for (auto track : tracks) {
      if (track.has_collision()) {
        nAssignedTracks++;
      } else {
        nNotAssignedTracks++;
        if (firstNotAssignedIndex < 0) {
          firstNotAssignedIndex = track.globalIndex();
        } else {
          lastNotAssignedIndex = track.globalIndex();
        }
      }
    }
    LOGF(info, "Received %d collisions and %d tracks.", collisions.size(), tracks.size());
    LOGF(info, "  Assigned tracks %d", nAssignedTracks);
    LOGF(info, "  Not assigned tracks %d", nNotAssignedTracks);
    LOGF(info, "  First not assigned track index %d", firstNotAssignedIndex);
    LOGF(info, "  Last not assigned track index %d", lastNotAssignedIndex);
  }
  PROCESS_SWITCH(DptDptCorrelationsTask, processRecLevelCheck, "Process reco level checks", true);

  void processGenLevelCheck(aod::McCollisions const& mccollisions, aod::McParticles& particles)
  {
    int nAssignedParticles = 0;
    int nNotAssignedParticles = 0;
    int64_t firstNotAssignedIndex = -1;
    int64_t lastNotAssignedIndex = -1;

    for (auto particle : particles) {
      if (particle.has_mcCollision()) {
        nAssignedParticles++;
      } else {
        nNotAssignedParticles++;
        if (firstNotAssignedIndex < 0) {
          firstNotAssignedIndex = particle.globalIndex();
        } else {
          lastNotAssignedIndex = particle.globalIndex();
        }
      }
    }
    LOGF(info, "Received %d generated collisions and %d particles.", mccollisions.size(), particles.size());
    LOGF(info, "  Assigned tracks %d", nAssignedParticles);
    LOGF(info, "  Not assigned tracks %d", nNotAssignedParticles);
    LOGF(info, "  First not assigned track index %d", firstNotAssignedIndex);
    LOGF(info, "  Last not assigned track index %d", lastNotAssignedIndex);
  }
  PROCESS_SWITCH(DptDptCorrelationsTask, processGenLevelCheck, "Process generator level checks", true);

  void processRecLevelNotStored(
    soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
    aod::BCsWithTimestamps const&,
    soa::Filtered<soa::Join<aod::Tracks, aod::DptDptCFTracksInfo>>& tracks)
  {
    processSame<false>(collision, tracks, collision.bc_as<aod::BCsWithTimestamps>().timestamp());
  }
  PROCESS_SWITCH(DptDptCorrelationsTask,
                 processRecLevelNotStored,
                 "Process reco level correlations for not stored derived data",
                 true);

  void processGenLevel(
    soa::Filtered<aod::DptDptCFAcceptedTrueCollisions>::iterator const& collision,
    soa::Filtered<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>>& tracks)
  {
    processSame<true>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptCorrelationsTask, processGenLevel, "Process generator level correlations", false);

  void processGenLevelNotStored(
    soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>::iterator const& collision,
    soa::Filtered<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>>& particles)
  {
    processSame<true>(collision, particles);
  }
  PROCESS_SWITCH(DptDptCorrelationsTask,
                 processGenLevelNotStored,
                 "Process generator level correlations for not stored derived data",
                 false);

  std::vector<double>
    vtxBinsEdges{VARIABLE_WIDTH, -7.0f, -5.0f, -3.0f, -1.0f, 1.0f, 3.0f, 5.0f, 7.0f};
  std::vector<double> multBinsEdges{VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0, 100.1f};
  SliceCache cache;
  using BinningZVtxMultRec = ColumnBinningPolicy<aod::collision::PosZ, aod::dptdptfilter::DptDptCFCollisionCentMult>;
  BinningZVtxMultRec bindingOnVtxAndMultRec{{vtxBinsEdges, multBinsEdges}, true}; // true is for 'ignore overflows' (true by default)

  void processRecLevelMixed(soa::Filtered<aod::DptDptCFAcceptedCollisions>& collisions, aod::BCsWithTimestamps const&, soa::Filtered<aod::ScannedTracks>& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<soa::Filtered<aod::DptDptCFAcceptedCollisions>, soa::Filtered<aod::ScannedTracks>, BinningZVtxMultRec> pairreco{bindingOnVtxAndMultRec, 5, -1, collisions, tracksTuple, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

    LOGF(DPTDPTLOGCOLLISIONS, "Received %d collisions", collisions.size());
    int logcomb = 0;
    for (auto& [collision1, tracks1, collision2, tracks2] : pairreco) {
      if (logcomb < 10) {
        LOGF(DPTDPTLOGCOLLISIONS, "Received collision pair: %ld (%f, %f): %s, %ld (%f, %f): %s",
             collision1.globalIndex(), collision1.posZ(), collision1.centmult(), collision1.collisionaccepted() ? "accepted" : "not accepted",
             collision2.globalIndex(), collision2.posZ(), collision2.centmult(), collision2.collisionaccepted() ? "accepted" : "not accepted");
        logcomb++;
      }
      if (!collision1.collisionaccepted() || !collision2.collisionaccepted()) {
        LOGF(error, "Received collision pair: %ld (%f, %f): %s, %ld (%f, %f): %s",
             collision1.globalIndex(), collision1.posZ(), collision1.centmult(), collision1.collisionaccepted() ? "accepted" : "not accepted",
             collision2.globalIndex(), collision2.posZ(), collision2.centmult(), collision2.collisionaccepted() ? "accepted" : "not accepted");
      }
      processMixed<false>(collision1, tracks1, tracks2, collision1.bc_as<aod::BCsWithTimestamps>().timestamp());
    }
  }
  PROCESS_SWITCH(DptDptCorrelationsTask, processRecLevelMixed, "Process reco level mixed events correlations", false);

  void processRecLevelMixedNotStored(
    soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>& collisions,
    aod::BCsWithTimestamps const&,
    soa::Filtered<soa::Join<aod::Tracks, aod::DptDptCFTracksInfo>>& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>,
                 soa::Filtered<soa::Join<aod::Tracks, aod::DptDptCFTracksInfo>>,
                 BinningZVtxMultRec>
      pairreco{bindingOnVtxAndMultRec,
               5,
               -1,
               collisions,
               tracksTuple,
               &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

    LOGF(DPTDPTLOGCOLLISIONS, "Received %d collisions", collisions.size());
    int logcomb = 0;
    for (auto& [collision1, tracks1, collision2, tracks2] : pairreco) {
      if (logcomb < 10) {
        LOGF(DPTDPTLOGCOLLISIONS,
             "Received collision pair: %ld (%f, %f): %s, %ld (%f, %f): %s",
             collision1.globalIndex(),
             collision1.posZ(),
             collision1.centmult(),
             collision1.collisionaccepted() ? "accepted" : "not accepted",
             collision2.globalIndex(),
             collision2.posZ(),
             collision2.centmult(),
             collision2.collisionaccepted() ? "accepted" : "not accepted");
        logcomb++;
      }
      if (!collision1.collisionaccepted() || !collision2.collisionaccepted()) {
        LOGF(error,
             "Received collision pair: %ld (%f, %f): %s, %ld (%f, %f): %s",
             collision1.globalIndex(),
             collision1.posZ(),
             collision1.centmult(),
             collision1.collisionaccepted() ? "accepted" : "not accepted",
             collision2.globalIndex(),
             collision2.posZ(),
             collision2.centmult(),
             collision2.collisionaccepted() ? "accepted" : "not accepted");
      }
      processMixed<false>(collision1,
                          tracks1,
                          tracks2,
                          collision1.bc_as<aod::BCsWithTimestamps>().timestamp());
    }
  }
  PROCESS_SWITCH(DptDptCorrelationsTask,
                 processRecLevelMixedNotStored,
                 "Process reco level mixed events correlations for not stored derived data",
                 false);

  using BinningZVtxMultGen = ColumnBinningPolicy<aod::mccollision::PosZ, aod::dptdptfilter::DptDptCFCollisionCentMult>;
  BinningZVtxMultGen bindingOnVtxAndMultGen{{vtxBinsEdges, multBinsEdges}, true}; // true is for 'ignore overflows' (true by default)

  void processGenLevelMixed(soa::Filtered<aod::DptDptCFAcceptedTrueCollisions>& collisions, soa::Filtered<aod::ScannedTrueTracks>& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<soa::Filtered<aod::DptDptCFAcceptedTrueCollisions>, soa::Filtered<aod::ScannedTrueTracks>, BinningZVtxMultGen> pairgen{bindingOnVtxAndMultGen, 5, -1, collisions, tracksTuple, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

    LOGF(DPTDPTLOGCOLLISIONS, "Received %d generated collisions", collisions.size());
    int logcomb = 0;
    for (auto& [collision1, tracks1, collision2, tracks2] : pairgen) {
      if (logcomb < 10) {
        LOGF(DPTDPTLOGCOLLISIONS, "Received generated collision pair: %ld (%f, %f): %s, %ld (%f, %f): %s",
             collision1.globalIndex(), collision1.posZ(), collision1.centmult(), collision1.collisionaccepted() ? "accepted" : "not accepted",
             collision2.globalIndex(), collision2.posZ(), collision2.centmult(), collision2.collisionaccepted() ? "accepted" : "not accepted");
      }
      if (!collision1.collisionaccepted() || !collision2.collisionaccepted()) {
        LOGF(error, "Received collision pair: %ld (%f, %f): %s, %ld (%f, %f): %s",
             collision1.globalIndex(), collision1.posZ(), collision1.centmult(), collision1.collisionaccepted() ? "accepted" : "not accepted",
             collision2.globalIndex(), collision2.posZ(), collision2.centmult(), collision2.collisionaccepted() ? "accepted" : "not accepted");
      }
      processMixed<true>(collision1, tracks1, tracks2);
    }
  }
  PROCESS_SWITCH(DptDptCorrelationsTask, processGenLevelMixed, "Process generator level mixed events correlations", false);

  void processGenLevelMixedNotStored(
    soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>& collisions,
    soa::Filtered<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>>& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>,
                 soa::Filtered<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>>,
                 BinningZVtxMultGen>
      pairgen{bindingOnVtxAndMultGen,
              5,
              -1,
              collisions,
              tracksTuple,
              &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

    LOGF(DPTDPTLOGCOLLISIONS, "Received %d generated collisions", collisions.size());
    int logcomb = 0;
    for (auto& [collision1, tracks1, collision2, tracks2] : pairgen) {
      if (logcomb < 10) {
        LOGF(DPTDPTLOGCOLLISIONS,
             "Received generated collision pair: %ld (%f, %f): %s, %ld (%f, %f): %s",
             collision1.globalIndex(),
             collision1.posZ(),
             collision1.centmult(),
             collision1.collisionaccepted() ? "accepted" : "not accepted",
             collision2.globalIndex(),
             collision2.posZ(),
             collision2.centmult(),
             collision2.collisionaccepted() ? "accepted" : "not accepted");
      }
      if (!collision1.collisionaccepted() || !collision2.collisionaccepted()) {
        LOGF(error,
             "Received collision pair: %ld (%f, %f): %s, %ld (%f, %f): %s",
             collision1.globalIndex(),
             collision1.posZ(),
             collision1.centmult(),
             collision1.collisionaccepted() ? "accepted" : "not accepted",
             collision2.globalIndex(),
             collision2.posZ(),
             collision2.centmult(),
             collision2.collisionaccepted() ? "accepted" : "not accepted");
      }
      processMixed<true>(collision1, tracks1, tracks2);
    }
  }
  PROCESS_SWITCH(DptDptCorrelationsTask,
                 processGenLevelMixedNotStored,
                 "Process generator level mixed events correlations for not stored derived data",
                 false);

  /// cleans the output object when the task is not used
  void processCleaner(soa::Filtered<aod::DptDptCFAcceptedCollisions> const& colls)
  {
    LOGF(DPTDPTLOGCOLLISIONS, "Got %d new collisions", colls.size());
    fOutput->Clear();
  }
  PROCESS_SWITCH(DptDptCorrelationsTask, processCleaner, "Cleaner process for not used output", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<DptDptCorrelationsTask>(cfgc, TaskName{"DptDptCorrelationsTaskRec"}, SetDefaultProcesses{{{"processRecLevel", true}, {"processRecLevelMixed", false}, {"processCleaner", false}}}),
    adaptAnalysisTask<DptDptCorrelationsTask>(cfgc, TaskName{"DptDptCorrelationsTaskGen"}, SetDefaultProcesses{{{"processGenLevel", true}, {"processGenLevelMixed", false}, {"processCleaner", false}}})};
  return workflow;
}
