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

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/DptDptFiltered.h"
#include "PWGCF/TableProducer/dptdptfilter.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/runDataProcessing.h"
#include "DataFormatsParameters/GRPObject.h"
#include <CCDB/BasicCCDBManager.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TParameter.h>
#include <TList.h>
#include <TDirectory.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile3D.h>

#include <cmath>
#include <ctime>

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
float etabinwidth = (etaup - etalow) / float(etabins);
float phibinwidth = (phiup - philow) / float(phibins);
int deltaetabins = etabins * 2 - 1;
float deltaetalow = etalow - etaup, deltaetaup = etaup - etalow;
float deltaetabinwidth = (deltaetaup - deltaetalow) / float(deltaetabins);
int deltaphibins = phibins;
float deltaphibinwidth = constants::math::TwoPI / deltaphibins;
float deltaphilow = 0.0 - deltaphibinwidth / 2.0;
float deltaphiup = constants::math::TwoPI - deltaphibinwidth / 2.0;

bool processpairs = false;
std::string fTaskConfigurationString = "PendingToConfigure";

PairCuts fPairCuts;              // pair suppression engine
bool fUseConversionCuts = false; // suppress resonances and conversions
bool fUseTwoTrackCut = false;    // suppress too close tracks
} // namespace correlationstask

// Task for building <dpt,dpt> correlations
struct DptDptCorrelationsTask {

  /* the data collecting engine */
  struct DataCollectingEngine {
    //============================================================================================
    // The DptDptCorrelationsAnalysisTask output objects
    //============================================================================================
    /* histograms */
    TH1F* fhVertexZA;                 //!<! the z vertex distribution for the current multiplicity/centrality class
    TH1F* fhN1_vsPt[2];               //!<! weighted single particle distribution vs \f$p_T\f$, track 1 and 2
    TH2F* fhN1_vsEtaPhi[2];           //!<! weighted single particle distribution vs \f$\eta,\;\phi\f$, track 1 and 2
    TH2F* fhSum1Pt_vsEtaPhi[2];       //!<! accumulated sum of weighted \f$p_T\f$ vs \f$\eta,\;\phi\f$, track 1 and 2
    TH3F* fhN1_vsZEtaPhiPt[2];        //!<! single particle distribution vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$, track 1 and 2
    TH3F* fhSum1Pt_vsZEtaPhiPt[2];    //!<! accumulated sum of weighted \f$p_T\f$ vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$, track 1 and 2
    TH3* fhNuaNue_vsZEtaPhiPt[2];     //!<! NUA+NUE correction vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$, track 1 and 2
    TH2* fhPtAvg_vsEtaPhi[2];         //!<! average \f$p_T\f$ vs \f$\eta,\;\phi\f$, track 1 and 2
    TH2F* fhN2_vsPtPt[4];             //!<! track 1 and 2 weighted two particle distribution vs \f${p_T}_1, {p_T}_2\f$
    TH2F* fhN2_vsDEtaDPhi[4];         //!<! two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$ 1-1,1-2,2-1,2-2, combinations
    TH2F* fhSum2PtPt_vsDEtaDPhi[4];   //!<! two-particle  \f$\sum {p_T}_1 {p_T}_2\f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ 1-1,1-2,2-1,2-2, combinations
    TH2F* fhSum2DptDpt_vsDEtaDPhi[4]; //!<! two-particle  \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ 1-1,1-2,2-1,2-2, combinations
    TH2F* fhSupN1N1_vsDEtaDPhi[4];    //!<! suppressed n1n1 two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$ 1-1,1-2,2-1,2-2, combinations
    TH2F* fhSupPt1Pt1_vsDEtaDPhi[4];  //!<! suppressed \f${p_T}_1 {p_T}_2\f$ two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$ 1-1,1-2,2-1,2-2, combinations
    /* versus centrality/multiplicity  profiles */
    TProfile* fhN1_vsC[2];           //!<! weighted single particle distribution vs event centrality/multiplicity, track 1 and 2
    TProfile* fhSum1Pt_vsC[2];       //!<! accumulated sum of weighted \f$p_T\f$ vs event centrality/multiplicity, track 1 and 2
    TProfile* fhN1nw_vsC[2];         //!<! un-weighted single particle distribution vs event centrality/multiplicity, track 1 and 2
    TProfile* fhSum1Ptnw_vsC[2];     //!<! accumulated sum of un-weighted \f$p_T\f$ vs event centrality/multiplicity, track 1 and 2
    TProfile* fhN2_vsC[4];           //!<! weighted accumulated two particle distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    TProfile* fhSum2PtPt_vsC[4];     //!<! weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    TProfile* fhSum2DptDpt_vsC[4];   //!<! weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    TProfile* fhN2nw_vsC[4];         //!<! un-weighted accumulated two particle distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    TProfile* fhSum2PtPtnw_vsC[4];   //!<! un-weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    TProfile* fhSum2DptDptnw_vsC[4]; //!<! un-weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations

    /// \enum TrackPairs
    /// \brief The track combinations hadled by the class
    typedef enum {
      kOO = 0,    ///< one-one pairs
      kOT,        ///< one-two pairs
      kTO,        ///< two-one pairs
      kTT,        ///< two-two pairs
      nTrackPairs ///< the number of track pairs
    } trackpairs;

    const char* tname[2] = {"1", "2"}; ///< the external track names, one and two, for histogram creation
    const char* trackPairsNames[4] = {"OO", "OT", "TO", "TT"};
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
      if (not(phi < phiup)) {
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
    /// a zero based bin index and similarlly for the track phi
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

      int etaix = int((t.eta() - etalow) / etabinwidth);
      /* consider a potential phi origin shift */
      float phi = GetShiftedPhi(t.phi());
      int phiix = int((phi - philow) / phibinwidth);
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
      int etaix_1 = int((t1.eta() - etalow) / etabinwidth);
      /* consider a potential phi origin shift */
      float phi = GetShiftedPhi(t1.phi());
      int phiix_1 = int((phi - philow) / phibinwidth);
      int etaix_2 = int((t2.eta() - etalow) / etabinwidth);
      /* consider a potential phi origin shift */
      phi = GetShiftedPhi(t2.phi());
      int phiix_2 = int((phi - philow) / phibinwidth);

      int deltaeta_ix = etaix_1 - etaix_2 + etabins - 1;
      int deltaphi_ix = phiix_1 - phiix_2;
      if (deltaphi_ix < 0) {
        deltaphi_ix += phibins;
      }

      return fhN2_vsDEtaDPhi[kOO]->GetBin(deltaeta_ix + 1, deltaphi_ix + 1);
    }

    void storeTrackCorrections(TH3* corrs1, TH3* corrs2)
    {
      LOGF(info, "Stored NUA&NUE corrections for track one %s and track two %s", corrs1 != nullptr ? "yes" : "no", corrs2 != nullptr ? "yes" : "no");
      fhNuaNue_vsZEtaPhiPt[0] = corrs1;
      fhNuaNue_vsZEtaPhiPt[1] = corrs2;
      ccdbstored = true;
    }

    void storePtAverages(TH2* ptavgs1, TH2* ptavgs2)
    {
      LOGF(info, "Stored pT average for track one %s and track two %s", ptavgs1 != nullptr ? "yes" : "no", ptavgs2 != nullptr ? "yes" : "no");
      fhPtAvg_vsEtaPhi[0] = ptavgs1;
      fhPtAvg_vsEtaPhi[1] = ptavgs2;
      ccdbstored = true;
    }

    template <typename TrackListObject>
    std::vector<float>* getTrackCorrections(TrackListObject const& tracks, int tix, float zvtx)
    {
      std::vector<float>* corr = new std::vector<float>(tracks.size(), 1.0f);
      if (fhNuaNue_vsZEtaPhiPt[tix] != nullptr) {
        int index = 0;
        for (auto t : tracks) {
          (*corr)[index] = fhNuaNue_vsZEtaPhiPt[tix]->GetBinContent(zvtx, GetEtaPhiIndex(t) + 0.5, t.pt());
          index++;
        }
      }
      return corr;
    }

    template <typename TrackListObject>
    std::vector<float>* getPtAvg(TrackListObject const& tracks, int tix)
    {
      std::vector<float>* ptavg = new std::vector<float>(tracks.size(), 0.0f);
      if (fhPtAvg_vsEtaPhi[tix] != nullptr) {
        int index = 0;
        for (auto t : tracks) {
          (*ptavg)[index] = fhPtAvg_vsEtaPhi[tix]->GetBinContent(fhPtAvg_vsEtaPhi[tix]->FindBin(t.eta(), t.phi()));
          index++;
        }
      }
      return ptavg;
    }

    /// \brief fills the singles histograms in singles execution mode
    /// \param passedtracks filtered table with the tracks associated to the passed index
    /// \param tix index, in the singles histogram bank, for the passed filetered track table
    template <typename TrackListObject>
    void processSingles(TrackListObject const& passedtracks, std::vector<float>* corrs, int tix, float zvtx)
    {
      int index = 0;
      for (auto& track : passedtracks) {
        float corr = (*corrs)[index];
        fhN1_vsPt[tix]->Fill(track.pt(), corr);
        fhN1_vsZEtaPhiPt[tix]->Fill(zvtx, GetEtaPhiIndex(track) + 0.5, track.pt(), corr);
        fhSum1Pt_vsZEtaPhiPt[tix]->Fill(zvtx, GetEtaPhiIndex(track) + 0.5, track.pt(), track.pt() * corr);
        index++;
      }
    }

    /// \brief fills the singles histograms in pair execution mode
    /// \param passedtracks filtered table with the tracks associated to the passed index
    /// \param tix index, in the singles histogram bank, for the passed filetered track table
    /// \param cmul centrality - multiplicity for the collision being analyzed
    template <typename TrackListObject>
    void processTracks(TrackListObject const& passedtracks, std::vector<float>* corrs, int tix, float cmul)
    {
      LOGF(DPTDPTLOGCOLLISIONS, "Processing %d tracks of type %d in a collision with cent/mult %f ", passedtracks.size(), tix, cmul);

      /* process magnitudes */
      double n1 = 0;       ///< weighted number of track 1 tracks for current collision
      double sum1Pt = 0;   ///< accumulated sum of weighted track 1 \f$p_T\f$ for current collision
      double n1nw = 0;     ///< not weighted number of track 1 tracks for current collision
      double sum1Ptnw = 0; ///< accumulated sum of not weighted track 1 \f$p_T\f$ for current collision
      int index = 0;
      for (auto& track : passedtracks) {
        float corr = (*corrs)[index];
        n1 += corr;
        sum1Pt += track.pt() * corr;
        n1nw += 1;
        sum1Ptnw += track.pt();

        fhN1_vsEtaPhi[tix]->Fill(track.eta(), GetShiftedPhi(track.phi()), corr);
        fhSum1Pt_vsEtaPhi[tix]->Fill(track.eta(), GetShiftedPhi(track.phi()), track.pt() * corr);
        index++;
      }
      fhN1_vsC[tix]->Fill(cmul, n1);
      fhSum1Pt_vsC[tix]->Fill(cmul, sum1Pt);
      fhN1nw_vsC[tix]->Fill(cmul, n1nw);
      fhSum1Ptnw_vsC[tix]->Fill(cmul, sum1Ptnw);
    }

    /// \brief fills the pair histograms in pair execution mode
    /// \param trks1 filtered table with the tracks associated to the first track in the pair
    /// \param trks2 filtered table with the tracks associated to the second track in the pair
    /// \param pix index, in the track combination histogram bank, for the passed filetered track tables
    /// \param cmul centrality - multiplicity for the collision being analyzed
    /// Be aware that at least in half of the cases traks1 and trks2 will have the same content
    template <trackpairs pix, typename TrackOneListObject, typename TrackTwoListObject>
    void processTrackPairs(TrackOneListObject const& trks1, TrackTwoListObject const& trks2, std::vector<float>* corrs1, std::vector<float>* corrs2, std::vector<float>* ptavgs1, std::vector<float>* ptavgs2, float cmul, int bfield)
    {
      using namespace correlationstask;

      /* process pair magnitudes */
      double n2 = 0;           ///< weighted number of track 1 track 2 pairs for current collision
      double n2sup = 0;        ///< weighted number of track 1 track 2 suppressed pairs for current collision
      double sum2PtPt = 0;     ///< accumulated sum of weighted track 1 track 2 \f${p_T}_1 {p_T}_2\f$ for current collision
      double sum2DptDpt = 0;   ///< accumulated sum of weighted number of track 1 tracks times weighted track 2 \f$p_T\f$ for current collision
      double n2nw = 0;         ///< not weighted number of track1 track 2 pairs for current collision
      double sum2PtPtnw = 0;   ///< accumulated sum of not weighted track 1 track 2 \f${p_T}_1 {p_T}_2\f$ for current collision
      double sum2DptDptnw = 0; ///< accumulated sum of not weighted number of track 1 tracks times not weighted track 2 \f$p_T\f$ for current collision
      int index1 = 0;

      for (auto& track1 : trks1) {
        double ptavg_1 = (*ptavgs1)[index1];
        double corr1 = (*corrs1)[index1];
        int index2 = 0;
        for (auto& track2 : trks2) {
          /* checking the same track id condition */
          if constexpr (pix == kOO or pix == kTT) {
            if (track1 == track2) {
              /* exclude autocorrelations */
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
          if ((fUseConversionCuts and fPairCuts.conversionCuts(track1, track2)) or (fUseTwoTrackCut and fPairCuts.twoTrackCut(track1, track2, bfield))) {
            /* suppress the pair */
            fhSupN1N1_vsDEtaDPhi[pix]->AddBinContent(globalbin, corr);
            fhSupPt1Pt1_vsDEtaDPhi[pix]->AddBinContent(globalbin, track1.pt() * track2.pt() * corr);
            n2sup += corr;
          } else {
            /* count the pair */
            n2 += corr;
            sum2PtPt += track1.pt() * track2.pt() * corr;
            sum2DptDpt += dptdptw;
            n2nw += 1;
            sum2PtPtnw += track1.pt() * track2.pt();
            sum2DptDptnw += dptdptnw;

            fhN2_vsDEtaDPhi[pix]->AddBinContent(globalbin, corr);
            fhSum2DptDpt_vsDEtaDPhi[pix]->AddBinContent(globalbin, dptdptw);
            fhSum2PtPt_vsDEtaDPhi[pix]->AddBinContent(globalbin, track1.pt() * track2.pt() * corr);
          }
          fhN2_vsPtPt[pix]->Fill(track1.pt(), track2.pt(), corr);
        }
      }
      fhN2_vsC[pix]->Fill(cmul, n2);
      fhSum2PtPt_vsC[pix]->Fill(cmul, sum2PtPt);
      fhSum2DptDpt_vsC[pix]->Fill(cmul, sum2DptDpt);
      fhN2nw_vsC[pix]->Fill(cmul, n2nw);
      fhSum2PtPtnw_vsC[pix]->Fill(cmul, sum2PtPtnw);
      fhSum2DptDptnw_vsC[pix]->Fill(cmul, sum2DptDptnw);
      /* let's also update the number of entries in the differential histograms */
      fhN2_vsDEtaDPhi[pix]->SetEntries(fhN2_vsDEtaDPhi[pix]->GetEntries() + n2);
      fhSum2DptDpt_vsDEtaDPhi[pix]->SetEntries(fhSum2DptDpt_vsDEtaDPhi[pix]->GetEntries() + n2);
      fhSum2PtPt_vsDEtaDPhi[pix]->SetEntries(fhSum2PtPt_vsDEtaDPhi[pix]->GetEntries() + n2);
      fhSupN1N1_vsDEtaDPhi[pix]->SetEntries(fhSupN1N1_vsDEtaDPhi[pix]->GetEntries() + n2sup);
      fhSupPt1Pt1_vsDEtaDPhi[pix]->SetEntries(fhSupPt1Pt1_vsDEtaDPhi[pix]->GetEntries() + n2sup);
    }

    template <typename TrackOneListObject, typename TrackTwoListObject>
    void processCollision(TrackOneListObject const& Tracks1, TrackTwoListObject const& Tracks2, float zvtx, float centmult, int bfield)
    {
      using namespace correlationstask;
      std::vector<float>* corrs1 = getTrackCorrections(Tracks1, 0, zvtx);
      std::vector<float>* corrs2 = getTrackCorrections(Tracks2, 1, zvtx);

      if (not processpairs) {
        /* process single tracks */
        fhVertexZA->Fill(zvtx);
        processSingles(Tracks1, corrs1, 0, zvtx); /* track one */
        processSingles(Tracks2, corrs2, 1, zvtx); /* track two */
      } else {
        /* process track magnitudes */
        std::vector<float>* ptavgs1 = getPtAvg(Tracks1, 0);
        std::vector<float>* ptavgs2 = getPtAvg(Tracks2, 1);

        /* TODO: the centrality should be chosen non detector dependent */
        processTracks(Tracks1, corrs1, 0, centmult); /* track one */
        processTracks(Tracks2, corrs2, 1, centmult); /* track one */
        /* process pair magnitudes */
        processTrackPairs<kOO>(Tracks1, Tracks1, corrs1, corrs1, ptavgs1, ptavgs1, centmult, bfield);
        processTrackPairs<kOT>(Tracks1, Tracks2, corrs1, corrs2, ptavgs1, ptavgs2, centmult, bfield);
        processTrackPairs<kTO>(Tracks2, Tracks1, corrs2, corrs1, ptavgs2, ptavgs1, centmult, bfield);
        processTrackPairs<kTT>(Tracks2, Tracks2, corrs2, corrs2, ptavgs2, ptavgs2, centmult, bfield);

        delete ptavgs1;
        delete ptavgs2;
      }
      delete corrs1;
      delete corrs2;
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
        for (int i = 0; i < 2; ++i) {
          /* histograms for each track, one and two */
          fhN1_vsPt[i] = new TH1F(TString::Format("n1_%s_vsPt", tname[i]).Data(),
                                  TString::Format("#LT n_{1} #GT;p_{t,%s} (GeV/c);#LT n_{1} #GT", tname[i]).Data(),
                                  ptbins, ptlow, ptup);
          /* we don't want the Sumw2 structure being created here */
          bool defSumw2 = TH1::GetDefaultSumw2();
          TH1::SetDefaultSumw2(false);
          fhN1_vsZEtaPhiPt[i] = new TH3F(TString::Format("n1_%s_vsZ_vsEtaPhi_vsPt", tname[i]).Data(),
                                         TString::Format("#LT n_{1} #GT;vtx_{z};#eta_{%s}#times#varphi_{%s};p_{t,%s} (GeV/c)", tname[i], tname[i], tname[i]).Data(),
                                         zvtxbins, zvtxlow, zvtxup, etabins * phibins, 0.0, double(etabins * phibins), ptbins, ptlow, ptup);
          fhSum1Pt_vsZEtaPhiPt[i] = new TH3F(TString::Format("sumPt1_%s_vsZ_vsEtaPhi_vsPt", tname[i]).Data(),
                                             TString::Format("#LT #Sigma p_{t,%s}#GT;vtx_{z};#eta_{%s}#times#varphi_{%s};p_{t,%s} (GeV/c)", tname[i], tname[i], tname[i], tname[i]).Data(),
                                             zvtxbins, zvtxlow, zvtxup, etabins * phibins, 0.0, double(etabins * phibins), ptbins, ptlow, ptup);
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
        for (int i = 0; i < 2; ++i) {
          /* histograms for each track, one and two */
          fhN1_vsEtaPhi[i] = new TH2F(TString::Format("n1_%s_vsEtaPhi", tname[i]).Data(),
                                      TString::Format("#LT n_{1} #GT;#eta_{%s};#varphi_{%s} (radian);#LT n_{1} #GT", tname[i], tname[i]).Data(),
                                      etabins, etalow, etaup, phibins, philow, phiup);
          fhSum1Pt_vsEtaPhi[i] = new TH2F(TString::Format("sumPt_%s_vsEtaPhi", tname[i]).Data(),
                                          TString::Format("#LT #Sigma p_{t,%s} #GT;#eta_{%s};#varphi_{%s} (radian);#LT #Sigma p_{t,%s} #GT (GeV/c)",
                                                          tname[i], tname[i], tname[i], tname[i])
                                            .Data(),
                                          etabins, etalow, etaup, phibins, philow, phiup);
          fhN1_vsC[i] = new TProfile(TString::Format("n1_%s_vsM", tname[i]).Data(),
                                     TString::Format("#LT n_{1} #GT (weighted);Centrality/Multiplicity (%%);#LT n_{1} #GT").Data(),
                                     100, 0.0, 100.0);
          fhSum1Pt_vsC[i] = new TProfile(TString::Format("sumPt_%s_vsM", tname[i]),
                                         TString::Format("#LT #Sigma p_{t,%s} #GT (weighted);Centrality/Multiplicity (%%);#LT #Sigma p_{t,%s} #GT (GeV/c)", tname[i], tname[i]).Data(),
                                         100, 0.0, 100.0);
          fhN1nw_vsC[i] = new TProfile(TString::Format("n1Nw_%s_vsM", tname[i]).Data(),
                                       TString::Format("#LT n_{1} #GT;Centrality/Multiplicity (%%);#LT n_{1} #GT").Data(),
                                       100, 0.0, 100.0);
          fhSum1Ptnw_vsC[i] = new TProfile(TString::Format("sumPtNw_%s_vsM", tname[i]).Data(),
                                           TString::Format("#LT #Sigma p_{t,%s} #GT;Centrality/Multiplicity (%%);#LT #Sigma p_{t,%s} #GT (GeV/c)", tname[i], tname[i]).Data(), 100, 0.0, 100.0);
          fhNuaNue_vsZEtaPhiPt[i] = nullptr;
          fhPtAvg_vsEtaPhi[i] = nullptr;
          fOutputList->Add(fhN1_vsEtaPhi[i]);
          fOutputList->Add(fhSum1Pt_vsEtaPhi[i]);
          fOutputList->Add(fhN1_vsC[i]);
          fOutputList->Add(fhSum1Pt_vsC[i]);
          fOutputList->Add(fhN1nw_vsC[i]);
          fOutputList->Add(fhSum1Ptnw_vsC[i]);
        }

        for (int i = 0; i < nTrackPairs; ++i) {
          /* histograms for each track pair combination */
          /* we don't want the Sumw2 structure being created here */
          bool defSumw2 = TH1::GetDefaultSumw2();
          TH1::SetDefaultSumw2(false);
          const char* pname = trackPairsNames[i];
          fhN2_vsDEtaDPhi[i] = new TH2F(TString::Format("n2_12_vsDEtaDPhi_%s", pname), TString::Format("#LT n_{2} #GT (%s);#Delta#eta;#Delta#varphi;#LT n_{2} #GT", pname),
                                        deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
          fhSum2PtPt_vsDEtaDPhi[i] = new TH2F(TString::Format("sumPtPt_12_vsDEtaDPhi_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s);#Delta#eta;#Delta#varphi;#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname),
                                              deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
          fhSum2DptDpt_vsDEtaDPhi[i] = new TH2F(TString::Format("sumDptDpt_12_vsDEtaDPhi_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s);#Delta#eta;#Delta#varphi;#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname),
                                                deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
          fhSupN1N1_vsDEtaDPhi[i] = new TH2F(TString::Format("suppn1n1_12_vsDEtaDPhi_%s", pname), TString::Format("Suppressed #LT n_{1} #GT#LT n_{1} #GT (%s);#Delta#eta;#Delta#varphi;#LT n_{1} #GT#LT n_{1} #GT", pname),
                                             deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
          fhSupPt1Pt1_vsDEtaDPhi[i] = new TH2F(TString::Format("suppPtPt_12_vsDEtaDPhi_%s", pname), TString::Format("Suppressed #LT p_{t,1} #GT#LT p_{t,2} #GT (%s);#Delta#eta;#Delta#varphi;#LT p_{t,1} #GT#LT p_{t,2} #GT (GeV^{2})", pname),
                                               deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
          /* we return it back to previuos state */
          TH1::SetDefaultSumw2(defSumw2);

          fhN2_vsPtPt[i] = new TH2F(TString::Format("n2_12_vsPtVsPt_%s", pname), TString::Format("#LT n_{2} #GT (%s);p_{t,1} (GeV/c);p_{t,2} (GeV/c);#LT n_{2} #GT", pname),
                                    ptbins, ptlow, ptup, ptbins, ptlow, ptup);

          fhN2_vsC[i] = new TProfile(TString::Format("n2_12_vsM_%s", pname), TString::Format("#LT n_{2} #GT (%s) (weighted);Centrality/Multiplicity (%%);#LT n_{2} #GT", pname), 100, 0.0, 100.0);
          fhSum2PtPt_vsC[i] = new TProfile(TString::Format("sumPtPt_12_vsM_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s) (weighted);Centrality/Multiplicity (%%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname), 100, 0.0, 100.0);
          fhSum2DptDpt_vsC[i] = new TProfile(TString::Format("sumDptDpt_12_vsM_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s) (weighted);Centrality/Multiplicity (%%);#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname), 100, 0.0, 100.0);
          fhN2nw_vsC[i] = new TProfile(TString::Format("n2Nw_12_vsM_%s", pname), TString::Format("#LT n_{2} #GT (%s);Centrality/Multiplicity (%%);#LT n_{2} #GT", pname), 100, 0.0, 100.0);
          fhSum2PtPtnw_vsC[i] = new TProfile(TString::Format("sumPtPtNw_12_vsM_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s);Centrality/Multiplicity (%%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname), 100, 0.0, 100.0);
          fhSum2DptDptnw_vsC[i] = new TProfile(TString::Format("sumDptDptNw_12_vsM_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s);Centrality/Multiplicity (%%);#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname), 100, 0.0, 100.0);

          /* the statistical uncertainties will be estimated by the subsamples method so let's get rid of the error tracking */
          fhN2_vsDEtaDPhi[i]->SetBit(TH1::kIsNotW);
          fhN2_vsDEtaDPhi[i]->Sumw2(false);
          fhSum2PtPt_vsDEtaDPhi[i]->SetBit(TH1::kIsNotW);
          fhSum2PtPt_vsDEtaDPhi[i]->Sumw2(false);
          fhSum2DptDpt_vsDEtaDPhi[i]->SetBit(TH1::kIsNotW);
          fhSum2DptDpt_vsDEtaDPhi[i]->Sumw2(false);
          fhSupN1N1_vsDEtaDPhi[i]->SetBit(TH1::kIsNotW);
          fhSupN1N1_vsDEtaDPhi[i]->Sumw2(false);
          fhSupPt1Pt1_vsDEtaDPhi[i]->SetBit(TH1::kIsNotW);
          fhSupPt1Pt1_vsDEtaDPhi[i]->Sumw2(false);

          fOutputList->Add(fhN2_vsDEtaDPhi[i]);
          fOutputList->Add(fhSum2PtPt_vsDEtaDPhi[i]);
          fOutputList->Add(fhSum2DptDpt_vsDEtaDPhi[i]);
          fOutputList->Add(fhSupN1N1_vsDEtaDPhi[i]);
          fOutputList->Add(fhSupPt1Pt1_vsDEtaDPhi[i]);
          fOutputList->Add(fhN2_vsPtPt[i]);
          fOutputList->Add(fhN2_vsC[i]);
          fOutputList->Add(fhSum2PtPt_vsC[i]);
          fOutputList->Add(fhSum2DptDpt_vsC[i]);
          fOutputList->Add(fhN2nw_vsC[i]);
          fOutputList->Add(fhSum2PtPtnw_vsC[i]);
          fOutputList->Add(fhSum2DptDptnw_vsC[i]);
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
  DataCollectingEngine** dataCE;

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
    loadfromccdb = cfginputfile.cfgCCDBPathName->length() > 0;
    /* update the potential binning change */
    etabinwidth = (etaup - etalow) / float(etabins);
    phibinwidth = (phiup - philow) / float(phibins);

    /* the differential bining */
    deltaetabins = etabins * 2 - 1;
    deltaetalow = etalow - etaup, deltaetaup = etaup - etalow;
    deltaetabinwidth = (deltaetaup - deltaetalow) / float(deltaetabins);
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
      dataCE = new DataCollectingEngine*[ncmranges];

      for (int i = 0; i < ncmranges; ++i) {
        float cmmin = 0.0f;
        float cmmax = 0.0f;
        sscanf(tokens->At(i)->GetName(), "%f-%f", &cmmin, &cmmax);
        fCentMultMin[i] = cmmin;
        fCentMultMax[i] = cmmax;
        dataCE[i] = new DataCollectingEngine();

        /* crete the output list for the current centrality range */
        TList* fOutputList = new TList();
        fOutputList->SetName(TString::Format("DptDptCorrelationsData-%s", tokens->At(i)->GetName()));
        fOutputList->SetOwner(true);
        /* init the data collection instance */
        dataCE[i]->init(fOutputList);
        fGlobalOutputList->Add(fOutputList);
      }
      delete tokens;
      for (int i = 0; i < ncmranges; ++i) {
        LOGF(info, " centrality/multipliicty range: %d, low limit: %f, up limit: %f", i, fCentMultMin[i], fCentMultMax[i]);
      }
    }
    /* two-track cut and conversion suppression */
    fPairCuts.SetHistogramRegistry(nullptr); // not histogram registry for the time being, incompatible with TList when it is empty
    if (processpairs and (cfgPairCut->get("Photon") > 0 or cfgPairCut->get("K0") > 0 or cfgPairCut->get("Lambda") > 0 or cfgPairCut->get("Phi") > 0 or cfgPairCut->get("Rho") > 0)) {
      fPairCuts.SetPairCut(PairCuts::Photon, cfgPairCut->get("Photon"));
      fPairCuts.SetPairCut(PairCuts::K0, cfgPairCut->get("K0"));
      fPairCuts.SetPairCut(PairCuts::Lambda, cfgPairCut->get("Lambda"));
      fPairCuts.SetPairCut(PairCuts::Phi, cfgPairCut->get("Phi"));
      fPairCuts.SetPairCut(PairCuts::Rho, cfgPairCut->get("Rho"));
      fUseConversionCuts = true;
    }
    if (processpairs and (cfgTwoTrackCut > 0)) {
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
    return ixDCE;
  }

  TList* getCCDBInput(const char* ccdbpath, const char* ccdbdate)
  {
    std::tm cfgtm = {};
    std::stringstream ss(ccdbdate);
    ss >> std::get_time(&cfgtm, "%Y%m%d");
    cfgtm.tm_hour = 12;
    long timestamp = std::mktime(&cfgtm) * 1000;

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

  Filter onlyacceptedcollisions = (aod::dptdptfilter::collisionaccepted == uint8_t(true));
  Filter onlyacceptedtracks = ((aod::dptdptfilter::trackacceptedasone == uint8_t(true)) or (aod::dptdptfilter::trackacceptedastwo == uint8_t(true)));

  void processRecLevel(soa::Filtered<aod::DptDptCFAcceptedCollisions>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<aod::ScannedTracks>& tracks)
  {
    using namespace correlationstask;

    if (ccdblst == nullptr) {
      if (loadfromccdb) {
        ccdblst = getCCDBInput(cfginputfile.cfgCCDBPathName->c_str(), cfginputfile.cfgCCDBDate->c_str());
      }
    }

    /* locate the data collecting engine for the collision centrality/multiplicity */
    int ixDCE = getDCEindex(collision);
    if (not(ixDCE < 0)) {
      if (ccdblst != nullptr and not dataCE[ixDCE]->isCCDBstored()) {
        dataCE[ixDCE]->storeTrackCorrections((TH3*)ccdblst->FindObject(TString::Format("correction_%02d-%02d_p1", int(fCentMultMin[ixDCE]), int(fCentMultMax[ixDCE])).Data()),
                                             (TH3*)ccdblst->FindObject(TString::Format("correction_%02d-%02d_m1", int(fCentMultMin[ixDCE]), int(fCentMultMax[ixDCE])).Data()));
        dataCE[ixDCE]->storePtAverages((TH2*)ccdblst->FindObject(TString::Format("ptavgetaphi_%02d-%02d_p", int(fCentMultMin[ixDCE]), int(fCentMultMax[ixDCE])).Data()),
                                       (TH2*)ccdblst->FindObject(TString::Format("ptavgetaphi_%02d-%02d_m", int(fCentMultMin[ixDCE]), int(fCentMultMax[ixDCE])).Data()));
      }

      Partition<o2::aod::ScannedTracks> TracksOne = aod::dptdptfilter::trackacceptedasone == uint8_t(true);
      Partition<o2::aod::ScannedTracks> TracksTwo = aod::dptdptfilter::trackacceptedastwo == uint8_t(true);
      TracksOne.bindTable(tracks);
      TracksTwo.bindTable(tracks);

      LOGF(DPTDPTLOGCOLLISIONS, "Accepted BC id %d collision with cent/mult %f and %d total tracks. Assigned DCE: %d", collision.bcId(), collision.centmult(), tracks.size(), ixDCE);
      LOGF(DPTDPTLOGCOLLISIONS, "Accepted new collision with cent/mult %f and %d type one tracks and %d type two tracks. Assigned DCE: %d", collision.centmult(), TracksOne.size(), TracksTwo.size(), ixDCE);
      int bfield = (fUseConversionCuts or fUseTwoTrackCut) ? getMagneticField(collision.bc_as<aod::BCsWithTimestamps>().timestamp()) : 0;
      dataCE[ixDCE]->processCollision(TracksOne, TracksTwo, collision.posZ(), collision.centmult(), bfield);
    }
  }
  PROCESS_SWITCH(DptDptCorrelationsTask, processRecLevel, "Process reco level correlations", false);

  void processGenLevel(soa::Filtered<aod::DptDptCFAcceptedTrueCollisions>::iterator const& collision, soa::Filtered<aod::ScannedTrueTracks>& tracks)
  {
    using namespace correlationstask;

    if (ccdblst == nullptr) {
      if (loadfromccdb) {
        ccdblst = getCCDBInput(cfginputfile.cfgCCDBPathName->c_str(), cfginputfile.cfgCCDBDate->c_str());
      }
    }

    /* locate the data collecting engine for the collision centrality/multiplicity */
    int ixDCE = getDCEindex(collision);
    if (not(ixDCE < 0)) {
      if (ccdblst != nullptr and not dataCE[ixDCE]->isCCDBstored()) {
        dataCE[ixDCE]->storePtAverages((TH2*)ccdblst->FindObject(TString::Format("trueptavgetaphi_%02d-%02d_p", int(fCentMultMin[ixDCE]), int(fCentMultMax[ixDCE])).Data()),
                                       (TH2*)ccdblst->FindObject(TString::Format("trueptavgetaphi_%02d-%02d_m", int(fCentMultMin[ixDCE]), int(fCentMultMax[ixDCE])).Data()));
      }

      Partition<o2::aod::ScannedTrueTracks> TracksOne = aod::dptdptfilter::trackacceptedasone == uint8_t(true);
      Partition<o2::aod::ScannedTrueTracks> TracksTwo = aod::dptdptfilter::trackacceptedastwo == uint8_t(true);
      TracksOne.bindTable(tracks);
      TracksTwo.bindTable(tracks);

      LOGF(DPTDPTLOGCOLLISIONS, "Accepted BC id %d generated collision with cent/mult %f and %d total tracks. Assigned DCE: %d", collision.bcId(), collision.centmult(), tracks.size(), ixDCE);
      LOGF(DPTDPTLOGCOLLISIONS, "Accepted new generated collision with cent/mult %f and %d type one tracks and %d type two tracks. Assigned DCE: %d", collision.centmult(), TracksOne.size(), TracksTwo.size(), ixDCE);
      dataCE[ixDCE]->processCollision(TracksOne, TracksTwo, collision.posZ(), collision.centmult(), 0.0); /* TODO: this needs more thinking, suppressing pairs on filtered generator level */
    }
  }
  PROCESS_SWITCH(DptDptCorrelationsTask, processGenLevel, "Process generator level correlations", false);

  /// cleans the output object when the task is not used
  void processCleaner(aod::Collisions const& colls)
  {
    LOGF(DPTDPTLOGCOLLISIONS, "Got %d new collisions", colls.size());
    fOutput->Clear();
  }
  PROCESS_SWITCH(DptDptCorrelationsTask, processCleaner, "Cleaner process for not used output", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<DptDptCorrelationsTask>(cfgc, TaskName{"DptDptCorrelationsTaskRec"}, SetDefaultProcesses{{{"processRecLevel", true}, {"processCleaner", false}}}),
    adaptAnalysisTask<DptDptCorrelationsTask>(cfgc, TaskName{"DptDptCorrelationsTaskGen"}, SetDefaultProcesses{{{"processGenLevel", true}, {"processCleaner", false}}})};
  return workflow;
}
