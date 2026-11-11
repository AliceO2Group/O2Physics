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

/// \file twoParticleCorrelations.cxx
/// \brief Two particle correlations analysis task
/// \author victor.gonzalez.sebastian@gmail.com

#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/TwoParticleCorrelations/Core/EventSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/FilterAndAnalysisFramework.h"
#include "PWGCF/TwoParticleCorrelations/Core/PIDSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/SelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/TrackSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/TwoParticleCorrelationsSkimmed.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPObject.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TList.h>
#include <TObjArray.h>
#include <TParameter.h>

#include <sys/types.h>

#include <RtypesCore.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <sstream>
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
std::vector<uint64_t> collisionmaskOptional;
uint64_t collisionmaskForced = 0UL;
uint64_t trackmask = 0UL;
std::vector<uint64_t> trackmaskOptional;
uint64_t trackmaskForced = 0UL;
uint64_t pidmask = 0UL;
std::vector<uint64_t> pidmaskOptional;
uint64_t pidmaskForced = 0UL;

/// \struct TrackId
/// \brief Stores the track species id, its potential correction, its potential pT average, and its global index to the track table
struct TrackId {
  int id;         ///< the species internal id
  float corr;     ///< correction to apply
  float ptavg;    ///< species pT average in eta phi
  uint64_t index; ///< the index to the table row
};

pwgcf::FilterAndAnalysisFramework* fFilterFramework = nullptr;

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
std::vector<TrackId> tracksids; // the tracks id storage

PairCuts fPairCuts;              // pair suppression engine
bool fUseConversionCuts = false; // suppress resonances and conversions
bool fUseTwoTrackCut = false;    // suppress too close tracks
} // namespace o2::analysis::twopcorrelations

using namespace o2::aod::cfskim;
using namespace o2::analysis::twopcorrelations;

struct TwoParticleCorrelations {

  /* the data collecting engine */
  struct DataCollectingEngine {
    //============================================================================================
    // The TwoParticleCorrelationsAnalysisTask output objects
    //============================================================================================
    /* histograms */
    TH1F* fhVertexZA;                                       //!<! the z vertex distribution for the current multiplicity/centrality class
    std::vector<TH1F*> fhN1VsPt;                            //!<! weighted single particle distribution vs \f$p_T\f$
    std::vector<TH2F*> fhN1VsEtaPhi;                        //!<! weighted single particle distribution vs \f$\eta,\;\phi\f$
    std::vector<TH2F*> fhSum1PtVsEtaPhi;                    //!<! accumulated sum of weighted \f$p_T\f$ vs \f$\eta,\;\phi\f$
    std::vector<TH3F*> fhN1VsZEtaPhiPt;                     //!<! single particle distribution vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$
    std::vector<TH3F*> fhSum1PtVsZEtaPhiPt;                 //!<! accumulated sum of weighted \f$p_T\f$ vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$
    std::vector<TH3*> fhNuaNueVsZEtaPhiPt;                  //!<! NUA+NUE correction vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$
    std::vector<TH2*> fhPtAvgVsEtaPhi;                      //!<! average \f$p_T\f$ vs \f$\eta,\;\phi\f$
    std::vector<std::vector<TH2F*>> fhN2VsPtPt;             //!<! weighted two particle distribution vs \f${p_T}_1, {p_T}_2\f$
    std::vector<std::vector<TH2F*>> fhN2VsDEtaDPhi;         //!<! two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$
    std::vector<std::vector<TH2F*>> fhSum2PtPtVsDEtaDPhi;   //!<! two-particle  \f$\sum {p_T}_1 {p_T}_2\f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$
    std::vector<std::vector<TH2F*>> fhSum2DptDptVsDEtaDPhi; //!<! two-particle  \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$
    std::vector<std::vector<TH2F*>> fhSupN1N1VsDEtaDPhi;    //!<! suppressed n1n1 two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$
    std::vector<std::vector<TH2F*>> fhSupPt1Pt1VsDEtaDPhi;  //!<! suppressed \f${p_T}_1 {p_T}_2\f$ two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$
    /* versus centrality/multiplicity  profiles */
    std::vector<TProfile*> fhN1VsC;                        //!<! weighted single particle distribution vs event centrality/multiplicity
    std::vector<TProfile*> fhSum1PtVsC;                    //!<! accumulated sum of weighted \f$p_T\f$ vs event centrality/multiplicity
    std::vector<TProfile*> fhN1nwVsC;                      //!<! un-weighted single particle distribution vs event centrality/multiplicity
    std::vector<TProfile*> fhSum1PtnwVsC;                  //!<! accumulated sum of un-weighted \f$p_T\f$ vs event centrality/multiplicity
    std::vector<std::vector<TProfile*>> fhN2VsC;           //!<! weighted accumulated two particle distribution vs event centrality/multiplicity
    std::vector<std::vector<TProfile*>> fhSum2PtPtVsC;     //!<! weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality/multiplicity
    std::vector<std::vector<TProfile*>> fhSum2DptDptVsC;   //!<! weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs event centrality/multiplicity
    std::vector<std::vector<TProfile*>> fhN2nwVsC;         //!<! un-weighted accumulated two particle distribution vs event centrality/multiplicity
    std::vector<std::vector<TProfile*>> fhSum2PtPtnwVsC;   //!<! un-weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality/multiplicity
    std::vector<std::vector<TProfile*>> fhSum2DptDptnwVsC; //!<! un-weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ distribution vs event centrality/multiplicity

    std::vector<std::string> tname; ///< the external track names for histogram creation
    bool ccdbstored = false;

    float isCCDBstored()
    {
      return ccdbstored;
    }

    /// \brief Returns the potentially phi origin shifted phi
    /// \param phi the track azimuthal angle
    /// \return the track phi origin shifted azimuthal angle
    float getShiftedPhi(float phi)
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
    int getEtaPhiIndex(TrackObject const& t)
    {
      using namespace twopcorrelations;

      int etaix = static_cast<int>((t.eta() - etalow) / etabinwidth);
      /* consider a potential phi origin shift */
      float phi = getShiftedPhi(t.phi());
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
    int getDEtaDPhiGlobalIndex(TrackObject const& t1, TrackObject const& t2)
    {
      using namespace twopcorrelations;

      /* rule: ix are always zero based while bins are always one based */
      int etaix1 = static_cast<int>((t1.eta() - etalow) / etabinwidth);
      /* consider a potential phi origin shift */
      float phi = getShiftedPhi(t1.phi());
      int phiix1 = static_cast<int>((phi - philow) / phibinwidth);
      int etaix2 = static_cast<int>((t2.eta() - etalow) / etabinwidth);
      /* consider a potential phi origin shift */
      phi = getShiftedPhi(t2.phi());
      int phiix2 = static_cast<int>((phi - philow) / phibinwidth);

      int deltaetaix = etaix1 - etaix2 + etabins - 1;
      int deltaphiix = phiix1 - phiix2;
      if (deltaphiix < 0) {
        deltaphiix += phibins;
      }

      return fhN2VsDEtaDPhi[0][0]->GetBin(deltaetaix + 1, deltaphiix + 1);
    }

    void storeTrackCorrections(const std::vector<TH3*>& corrs)
    {
      LOGF(info, "Stored NUA&NUE corrections for %d track species", corrs.size());
      fhNuaNueVsZEtaPhiPt = corrs;
      ccdbstored = true;
    }

    void storePtAverages(const std::vector<TH2*>& ptavgs)
    {
      LOGF(info, "Stored pT average for %d track species", ptavgs.size());
      fhPtAvgVsEtaPhi = ptavgs;
      ccdbstored = true;
    }

    template <typename TrackListObject>
    void getTrackPtAvgAndCorrections(TrackListObject& tracks, std::vector<TrackId>& tracksids, float zvtx)
    {
      for (auto& tid : tracksids) {
        auto const& track = tracks.iteratorAt(tid.index);
        if (!(tid.id < 0)) {
          if (fhNuaNueVsZEtaPhiPt[tid.id] != nullptr) {
            tid.corr = fhNuaNueVsZEtaPhiPt[tid.id]->GetBinContent(zvtx, getEtaPhiIndex(track) + 0.5, track.pt());
          } else {
            tid.corr = 1;
          }
          if (fhPtAvgVsEtaPhi[tid.id] != nullptr) {
            tid.ptavg = fhPtAvgVsEtaPhi[tid.id]->GetBinContent(fhPtAvgVsEtaPhi[tid.id]->FindBin(track.eta(), track.phi()));
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
    void processSingles(TrackListObject& passedtracks, std::vector<TrackId>& tracksids, float zvtx)
    {
      for (const auto& tid : tracksids) {
        auto const& track = passedtracks.iteratorAt(tid.index);
        if (!(tid.id < 0)) {
          fhN1VsPt[tid.id]->Fill(track.pt(), tid.corr);
          fhN1VsZEtaPhiPt[tid.id]->Fill(zvtx, getEtaPhiIndex(track) + 0.5, track.pt(), tid.corr);
          fhSum1PtVsZEtaPhiPt[tid.id]->Fill(zvtx, getEtaPhiIndex(track) + 0.5, track.pt(), track.pt() * tid.corr);
        }
      }
    }

    /// \brief fills the singles histograms in pair execution mode
    /// \param passedtracks filtered table with the accepted tracks
    /// \param tracksids id, index to the track table, correction, and average pT the passed track table
    /// \param cmul centrality - multiplicity for the collision being analyzed
    template <typename TrackListObject>
    void processTracks(TrackListObject& passedtracks, std::vector<TrackId>& tracksids, float cmul)
    {
      LOGF(TWOPCORRLOGCOLLISIONS, "Processing %d tracks in a collision with cent/mult %f ", passedtracks.size(), cmul);

      /* process magnitudes */
      std::vector<double> n1(tname.size(), 0.0);       ///< weighted number of tracks for current collision
      std::vector<double> sum1Pt(tname.size(), 0.0);   ///< accumulated sum of weighted track 1 \f$p_T\f$ for current collision
      std::vector<double> n1nw(tname.size(), 0.0);     ///< not weighted number of tracks for current collision
      std::vector<double> sum1Ptnw(tname.size(), 0.0); ///< accumulated sum of not weighted track \f$p_T\f$ for current collision

      for (const auto& tid : tracksids) {
        auto const& track = passedtracks.iteratorAt(tid.index);
        if (!(tid.id < 0)) {
          n1[tid.id] += tid.corr;
          sum1Pt[tid.id] += track.pt() * tid.corr;
          n1nw[tid.id] += 1;
          sum1Ptnw[tid.id] += track.pt();

          fhN1VsEtaPhi[tid.id]->Fill(track.eta(), getShiftedPhi(track.phi()), tid.corr);
          fhSum1PtVsEtaPhi[tid.id]->Fill(track.eta(), getShiftedPhi(track.phi()), track.pt() * tid.corr);
        }
      }
      for (uint i = 0; i < n1.size(); ++i) {
        fhN1VsC[i]->Fill(cmul, n1[i]);
        fhSum1PtVsC[i]->Fill(cmul, sum1Pt[i]);
        fhN1nwVsC[i]->Fill(cmul, n1nw[i]);
        fhSum1PtnwVsC[i]->Fill(cmul, sum1Ptnw[i]);
      }
    }

    /// \brief fills the pair histograms in pair execution mode
    /// \param trks1 filtered table with the tracks associated to the first track in the pair
    /// \param trks2 filtered table with the tracks associated to the second track in the pair
    /// \param pix index, in the track combination histogram bank, for the passed filetered track tables
    /// \param cmul centrality - multiplicity for the collision being analyzed
    /// Be aware that at least in half of the cases traks1 and trks2 will have the same content
    template <typename TrackListObject>
    void processTrackPairs(TrackListObject& tracks, std::vector<TrackId>& tracksids, float cmul, int bfield)
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

      for (const auto& tid1 : tracksids) {
        auto const& track1 = tracks.iteratorAt(tid1.index);
        if (!(tid1.id < 0)) {
          for (const auto& tid2 : tracksids) {
            auto const& track2 = tracks.iteratorAt(tid2.index);
            if (!(tid2.id < 0)) {
              if (tid1.index == tid2.index) {
                /* exclude auto correlations */
                continue;
              }
              /* get the global bin for filling the differential histograms */
              int globalbin = getDEtaDPhiGlobalIndex(track1, track2);
              if ((fUseConversionCuts && fPairCuts.conversionCuts(track1, track2)) || (fUseTwoTrackCut && fPairCuts.twoTrackCut(track1, track2, bfield))) {
                /* suppress the pair */
                fhSupN1N1VsDEtaDPhi[tid1.id][tid2.id]->AddBinContent(globalbin, tid1.corr * tid2.corr);
                fhSupPt1Pt1VsDEtaDPhi[tid1.id][tid2.id]->AddBinContent(globalbin, track1.pt() * track2.pt() * tid1.corr * tid2.corr);
                n2sup[tid1.id][tid2.id] += tid1.corr * tid2.corr;
              } else {
                /* count the pair */
                n2[tid1.id][tid2.id] += tid1.corr * tid2.corr;
                sum2PtPt[tid1.id][tid2.id] += track1.pt() * track2.pt() * tid1.corr * tid2.corr;
                sum2DptDpt[tid1.id][tid2.id] += (tid1.corr * track1.pt() - tid1.ptavg) * (tid2.corr * track2.pt() - tid2.ptavg);
                n2nw[tid1.id][tid2.id] += 1;
                sum2PtPtnw[tid1.id][tid2.id] += track1.pt() * track2.pt();
                sum2DptDptnw[tid1.id][tid2.id] += (track1.pt() - tid1.ptavg) * (track2.pt() - tid2.ptavg);

                fhN2VsDEtaDPhi[tid1.id][tid2.id]->AddBinContent(globalbin, tid1.corr * tid2.corr);
                fhSum2DptDptVsDEtaDPhi[tid1.id][tid2.id]->AddBinContent(globalbin, (tid1.corr * track1.pt() - tid1.ptavg) * (tid2.corr * track2.pt() - tid2.ptavg));
                fhSum2PtPtVsDEtaDPhi[tid1.id][tid2.id]->AddBinContent(globalbin, track1.pt() * track2.pt() * tid1.corr * tid2.corr);
                fhN2VsPtPt[tid1.id][tid2.id]->Fill(track1.pt(), track2.pt(), tid1.corr * tid2.corr);
              }
            }
          }
        }
      }
      for (uint i = 0; i < tname.size(); ++i) {
        for (uint j = 0; j < tname.size(); ++j) {
          fhN2VsC[i][j]->Fill(cmul, n2[i][j]);
          fhSum2PtPtVsC[i][j]->Fill(cmul, sum2PtPt[i][j]);
          fhSum2DptDptVsC[i][j]->Fill(cmul, sum2DptDpt[i][j]);
          fhN2nwVsC[i][j]->Fill(cmul, n2nw[i][j]);
          fhSum2PtPtnwVsC[i][j]->Fill(cmul, sum2PtPtnw[i][j]);
          fhSum2DptDptnwVsC[i][j]->Fill(cmul, sum2DptDptnw[i][j]);
          /* let's also update the number of entries in the differential histograms */
          fhN2VsDEtaDPhi[i][j]->SetEntries(fhN2VsDEtaDPhi[i][j]->GetEntries() + n2[i][j]);
          fhSum2DptDptVsDEtaDPhi[i][j]->SetEntries(fhSum2DptDptVsDEtaDPhi[i][j]->GetEntries() + n2[i][j]);
          fhSum2PtPtVsDEtaDPhi[i][j]->SetEntries(fhSum2PtPtVsDEtaDPhi[i][j]->GetEntries() + n2[i][j]);
          fhSupN1N1VsDEtaDPhi[i][j]->SetEntries(fhSupN1N1VsDEtaDPhi[i][j]->GetEntries() + n2sup[i][j]);
          fhSupPt1Pt1VsDEtaDPhi[i][j]->SetEntries(fhSupPt1Pt1VsDEtaDPhi[i][j]->GetEntries() + n2sup[i][j]);
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

    void init(TList* fOutputList, const std::vector<std::string>& idnames)
    {
      LOGF(info, "Correlations processing engine::init()");
      using namespace twopcorrelations;

      /* create the histograms */
      Bool_t oldstatus = TH1::AddDirectoryStatus();
      TH1::AddDirectory(kFALSE);

      /* load the species names */
      for (const auto& id : idnames) {
        LOGF(info, "Adding particle species %s", id.c_str());
        tname.push_back(std::string(id.c_str()));
      }
      /* reserve the histograms place holders */
      for (uint i = 0; i < tname.size(); ++i) {
        if (!processpairs) {
          fhN1VsPt.reserve(tname.size());
          fhN1VsZEtaPhiPt.reserve(tname.size());
          fhSum1PtVsZEtaPhiPt.reserve(tname.size());
          fhNuaNueVsZEtaPhiPt.reserve(tname.size());
          fhPtAvgVsEtaPhi.reserve(tname.size());
        } else {
          fhN1VsEtaPhi.reserve(tname.size());
          fhSum1PtVsEtaPhi.reserve(tname.size());
          fhN1VsC.reserve(tname.size());
          fhSum1PtVsC.reserve(tname.size());
          fhN1nwVsC.reserve(tname.size());
          fhSum1PtnwVsC.reserve(tname.size());
          fhNuaNueVsZEtaPhiPt.reserve(tname.size());
          fhPtAvgVsEtaPhi.reserve(tname.size());
          fhN2VsDEtaDPhi.reserve(tname.size());
          fhSum2PtPtVsDEtaDPhi.reserve(tname.size());
          fhSum2DptDptVsDEtaDPhi.reserve(tname.size());
          fhSupN1N1VsDEtaDPhi.reserve(tname.size());
          fhSupPt1Pt1VsDEtaDPhi.reserve(tname.size());
          fhN2VsPtPt.reserve(tname.size());
          fhN2VsC.reserve(tname.size());
          fhSum2PtPtVsC.reserve(tname.size());
          fhSum2DptDptVsC.reserve(tname.size());
          fhN2nwVsC.reserve(tname.size());
          fhSum2PtPtnwVsC.reserve(tname.size());
          fhSum2DptDptnwVsC.reserve(tname.size());
          for (uint j = 0; j < tname.size(); ++j) {
            fhN2VsDEtaDPhi.push_back({tname.size(), nullptr});
            fhSum2PtPtVsDEtaDPhi.push_back({tname.size(), nullptr});
            fhSum2DptDptVsDEtaDPhi.push_back({tname.size(), nullptr});
            fhSupN1N1VsDEtaDPhi.push_back({tname.size(), nullptr});
            fhSupPt1Pt1VsDEtaDPhi.push_back({tname.size(), nullptr});
            fhN2VsPtPt.push_back({tname.size(), nullptr});
            fhN2VsC.push_back({tname.size(), nullptr});
            fhSum2PtPtVsC.push_back({tname.size(), nullptr});
            fhSum2DptDptVsC.push_back({tname.size(), nullptr});
            fhN2nwVsC.push_back({tname.size(), nullptr});
            fhSum2PtPtnwVsC.push_back({tname.size(), nullptr});
            fhSum2DptDptnwVsC.push_back({tname.size(), nullptr});
          }
        }
      }

      if (!processpairs) {
        fhVertexZA = new TH1F("VertexZA", "Vertex Z; z_{vtx}", zvtxbins, zvtxlow, zvtxup);
        fOutputList->Add(fhVertexZA);
        for (uint i = 0; i < tname.size(); ++i) {
          /* histograms for each track, one and two */
          fhN1VsPt[i] = new TH1F(TString::Format("n1_%sVsPt", tname[i].c_str()).Data(),
                                 TString::Format("#LT n_{1} #GT;p_{t,%s} (GeV/c);#LT n_{1} #GT", tname[i].c_str()).Data(),
                                 ptbins, ptlow, ptup);
          /* we don't want the Sumw2 structure being created here */
          bool defSumw2 = TH1::GetDefaultSumw2();
          TH1::SetDefaultSumw2(false);
          fhN1VsZEtaPhiPt[i] = new TH3F(TString::Format("n1_%sVsZVsEtaPhiVsPt", tname[i].c_str()).Data(),
                                        TString::Format("#LT n_{1} #GT;vtx_{z};#eta_{%s}#times#varphi_{%s};p_{t,%s} (GeV/c)", tname[i].c_str(), tname[i].c_str(), tname[i].c_str()).Data(),
                                        zvtxbins, zvtxlow, zvtxup, etabins * phibins, 0.0, static_cast<double>(etabins * phibins), ptbins, ptlow, ptup);
          fhSum1PtVsZEtaPhiPt[i] = new TH3F(TString::Format("sumPt1_%sVsZVsEtaPhiVsPt", tname[i].c_str()).Data(),
                                            TString::Format("#LT #Sigma p_{t,%s}#GT;vtx_{z};#eta_{%s}#times#varphi_{%s};p_{t,%s} (GeV/c)", tname[i].c_str(), tname[i].c_str(), tname[i].c_str(), tname[i].c_str()).Data(),
                                            zvtxbins, zvtxlow, zvtxup, etabins * phibins, 0.0, static_cast<double>(etabins * phibins), ptbins, ptlow, ptup);
          /* we return it back to previuos state */
          TH1::SetDefaultSumw2(defSumw2);

          /* the statistical uncertainties will be estimated by the subsamples method so let's get rid of the error tracking */
          fhN1VsZEtaPhiPt[i]->SetBit(TH1::kIsNotW);
          fhN1VsZEtaPhiPt[i]->Sumw2(false);
          fhSum1PtVsZEtaPhiPt[i]->SetBit(TH1::kIsNotW);
          fhSum1PtVsZEtaPhiPt[i]->Sumw2(false);
          fhNuaNueVsZEtaPhiPt[i] = nullptr;
          fhPtAvgVsEtaPhi[i] = nullptr;

          fOutputList->Add(fhN1VsPt[i]);
          fOutputList->Add(fhN1VsZEtaPhiPt[i]);
          fOutputList->Add(fhSum1PtVsZEtaPhiPt[i]);
        }
      } else {
        for (uint i = 0; i < tname.size(); ++i) {
          /* histograms for each track species */
          fhN1VsEtaPhi[i] = new TH2F(TString::Format("n1_%sVsEtaPhi", tname[i].c_str()).Data(),
                                     TString::Format("#LT n_{1} #GT;#eta_{%s};#varphi_{%s} (radian);#LT n_{1} #GT", tname[i].c_str(), tname[i].c_str()).Data(),
                                     etabins, etalow, etaup, phibins, philow, phiup);
          fhSum1PtVsEtaPhi[i] = new TH2F(TString::Format("sumPt_%sVsEtaPhi", tname[i].c_str()).Data(),
                                         TString::Format("#LT #Sigma p_{t,%s} #GT;#eta_{%s};#varphi_{%s} (radian);#LT #Sigma p_{t,%s} #GT (GeV/c)",
                                                         tname[i].c_str(), tname[i].c_str(), tname[i].c_str(), tname[i].c_str())
                                           .Data(),
                                         etabins, etalow, etaup, phibins, philow, phiup);
          fhN1VsC[i] = new TProfile(TString::Format("n1_%sVsM", tname[i].c_str()).Data(),
                                    TString::Format("#LT n_{1} #GT (weighted);Centrality/Multiplicity (%%);#LT n_{1} #GT").Data(),
                                    100, 0.0, 100.0);
          fhSum1PtVsC[i] = new TProfile(TString::Format("sumPt_%sVsM", tname[i].c_str()),
                                        TString::Format("#LT #Sigma p_{t,%s} #GT (weighted);Centrality/Multiplicity (%%);#LT #Sigma p_{t,%s} #GT (GeV/c)", tname[i].c_str(), tname[i].c_str()).Data(),
                                        100, 0.0, 100.0);
          fhN1nwVsC[i] = new TProfile(TString::Format("n1Nw_%sVsM", tname[i].c_str()).Data(),
                                      TString::Format("#LT n_{1} #GT;Centrality/Multiplicity (%%);#LT n_{1} #GT").Data(),
                                      100, 0.0, 100.0);
          fhSum1PtnwVsC[i] = new TProfile(TString::Format("sumPtNw_%sVsM", tname[i].c_str()).Data(),
                                          TString::Format("#LT #Sigma p_{t,%s} #GT;Centrality/Multiplicity (%%);#LT #Sigma p_{t,%s} #GT (GeV/c)", tname[i].c_str(), tname[i].c_str()).Data(), 100, 0.0, 100.0);
          fhNuaNueVsZEtaPhiPt[i] = nullptr;
          fhPtAvgVsEtaPhi[i] = nullptr;
          fOutputList->Add(fhN1VsEtaPhi[i]);
          fOutputList->Add(fhSum1PtVsEtaPhi[i]);
          fOutputList->Add(fhN1VsC[i]);
          fOutputList->Add(fhSum1PtVsC[i]);
          fOutputList->Add(fhN1nwVsC[i]);
          fOutputList->Add(fhSum1PtnwVsC[i]);
          for (uint j = 0; j < tname.size(); ++j) {
            /* histograms for each track pair combination */
            /* we don't want the Sumw2 structure being created here */
            bool defSumw2 = TH1::GetDefaultSumw2();
            TH1::SetDefaultSumw2(false);
            char pname[256];
            snprintf(pname, tname[i].length() + tname[j].length() + 1, "%s%s", tname[i].c_str(), tname[j].c_str());
            fhN2VsDEtaDPhi[i][j] = new TH2F(TString::Format("n2_12VsDEtaDPhi_%s", pname), TString::Format("#LT n_{2} #GT (%s);#Delta#eta;#Delta#varphi;#LT n_{2} #GT", pname),
                                            deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
            fhSum2PtPtVsDEtaDPhi[i][j] = new TH2F(TString::Format("sumPtPt_12VsDEtaDPhi_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s);#Delta#eta;#Delta#varphi;#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname),
                                                  deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
            fhSum2DptDptVsDEtaDPhi[i][j] = new TH2F(TString::Format("sumDptDpt_12VsDEtaDPhi_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s);#Delta#eta;#Delta#varphi;#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname),
                                                    deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
            fhSupN1N1VsDEtaDPhi[i][j] = new TH2F(TString::Format("suppn1n1_12VsDEtaDPhi_%s", pname), TString::Format("Suppressed #LT n_{1} #GT#LT n_{1} #GT (%s);#Delta#eta;#Delta#varphi;#LT n_{1} #GT#LT n_{1} #GT", pname),
                                                 deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
            fhSupPt1Pt1VsDEtaDPhi[i][j] = new TH2F(TString::Format("suppPtPt_12VsDEtaDPhi_%s", pname), TString::Format("Suppressed #LT p_{t,1} #GT#LT p_{t,2} #GT (%s);#Delta#eta;#Delta#varphi;#LT p_{t,1} #GT#LT p_{t,2} #GT (GeV^{2})", pname),
                                                   deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
            /* we return it back to previuos state */
            TH1::SetDefaultSumw2(defSumw2);

            fhN2VsPtPt[i][j] = new TH2F(TString::Format("n2_12VsPtVsPt_%s", pname), TString::Format("#LT n_{2} #GT (%s);p_{t,1} (GeV/c);p_{t,2} (GeV/c);#LT n_{2} #GT", pname),
                                        ptbins, ptlow, ptup, ptbins, ptlow, ptup);

            fhN2VsC[i][j] = new TProfile(TString::Format("n2_12VsM_%s", pname), TString::Format("#LT n_{2} #GT (%s) (weighted);Centrality/Multiplicity (%%);#LT n_{2} #GT", pname), 100, 0.0, 100.0);
            fhSum2PtPtVsC[i][j] = new TProfile(TString::Format("sumPtPt_12VsM_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s) (weighted);Centrality/Multiplicity (%%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname), 100, 0.0, 100.0);
            fhSum2DptDptVsC[i][j] = new TProfile(TString::Format("sumDptDpt_12VsM_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s) (weighted);Centrality/Multiplicity (%%);#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname), 100, 0.0, 100.0);
            fhN2nwVsC[i][j] = new TProfile(TString::Format("n2Nw_12VsM_%s", pname), TString::Format("#LT n_{2} #GT (%s);Centrality/Multiplicity (%%);#LT n_{2} #GT", pname), 100, 0.0, 100.0);
            fhSum2PtPtnwVsC[i][j] = new TProfile(TString::Format("sumPtPtNw_12VsM_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s);Centrality/Multiplicity (%%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname), 100, 0.0, 100.0);
            fhSum2DptDptnwVsC[i][j] = new TProfile(TString::Format("sumDptDptNw_12VsM_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s);Centrality/Multiplicity (%%);#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname), 100, 0.0, 100.0);

            /* the statistical uncertainties will be estimated by the subsamples method so let's get rid of the error tracking */
            fhN2VsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhN2VsDEtaDPhi[i][j]->Sumw2(false);
            fhSum2PtPtVsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhSum2PtPtVsDEtaDPhi[i][j]->Sumw2(false);
            fhSum2DptDptVsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhSum2DptDptVsDEtaDPhi[i][j]->Sumw2(false);
            fhSupN1N1VsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhSupN1N1VsDEtaDPhi[i][j]->Sumw2(false);
            fhSupPt1Pt1VsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
            fhSupPt1Pt1VsDEtaDPhi[i][j]->Sumw2(false);

            fOutputList->Add(fhN2VsDEtaDPhi[i][j]);
            fOutputList->Add(fhSum2PtPtVsDEtaDPhi[i][j]);
            fOutputList->Add(fhSum2DptDptVsDEtaDPhi[i][j]);
            fOutputList->Add(fhSupN1N1VsDEtaDPhi[i][j]);
            fOutputList->Add(fhSupPt1Pt1VsDEtaDPhi[i][j]);
            fOutputList->Add(fhN2VsPtPt[i][j]);
            fOutputList->Add(fhN2VsC[i][j]);
            fOutputList->Add(fhSum2PtPtVsC[i][j]);
            fOutputList->Add(fhSum2DptDptVsC[i][j]);
            fOutputList->Add(fhN2nwVsC[i][j]);
            fOutputList->Add(fhSum2PtPtnwVsC[i][j]);
            fOutputList->Add(fhSum2DptDptnwVsC[i][j]);
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
  static constexpr float CfgPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};
  Configurable<LabeledArray<float>> cfgPairCut{"cfgPairCut", {CfgPairCutDefaults[0], 5, {"Photon", "K0", "Lambda", "Phi", "Rho"}}, "Conversion suppressions"};
  /* two tracks cut */
  Configurable<float> cfgTwoTrackCut{"cfgTwoTrackCut", -1, "Two-tracks cut: -1 = off; >0 otherwise distance value (suggested: 0.02"};
  Configurable<float> cfgTwoTrackCutMinRadius{"cfgTwoTrackCutMinRadius", 0.8f, "Two-tracks cut: radius in m from which two-tracks cut is applied"};

  Configurable<bool> cfgProcessPairs{"cfgProcessPairs", false, "Process pairs: false = no, just singles, true = yes, process pairs"};
  Configurable<std::string> cfgCentSpec{"cfgCentSpec", "00-05,05-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80", "Centrality/multiplicity ranges in min-max separated by commas"};

  Configurable<o2::analysis::DptDptBinningCuts> cfgBinning{"cfgBinning",
                                                           {28, -7.0, 7.0, 18, 0.2, 2.0, 16, -0.8, 0.8, 72, 0.5},
                                                           "triplets - nbins, min, max - for z_vtx, pT, eta and phi, binning plus bin fraction of phi origin shift"};
  struct : ConfigurableGroup {
    Configurable<std::string> cfgCCDBUrl{"cfgCCDBUrl", "http://ccdb-test.cern.ch:8080", "The CCDB url for the input file"};
    Configurable<std::string> cfgCCDBPathName{"cfgCCDBPathName", "", "The CCDB path for the input file. Default \"\", i.e. don't load from CCDB"};
    Configurable<std::string> cfgCCDBDate{"cfgCCDBDate", "20220307", "The CCDB date for the input file"};
  } cfginputfile;

  OutputObj<TList> fOutput{"DptDptCorrelationsData", OutputObjHandlingPolicy::AnalysisObject, OutputObjSourceType::OutputObjSource};

  void init(InitContext const&)
  {
    using namespace twopcorrelations;
    LOGF(info, "TwoParticleCorrelations::init()");

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
    pwgcf::EventSelectionConfigurable eventsel(eventfilter.bfield, eventfilter.centmultsel, {}, eventfilter.zvtxsel, eventfilter.pileuprej);
    /* track filtering configuration */
    pwgcf::TrackSelectionConfigurable trksel(trackfilter.ttype, trackfilter.nclstpc, trackfilter.nxrtpc, trackfilter.nclsits, trackfilter.chi2clustpc,
                                             trackfilter.chi2clusits, trackfilter.xrofctpc, trackfilter.dcaxy, trackfilter.dcaz, trackfilter.ptrange, trackfilter.etarange);
#ifdef INCORPORATEBAYESIANPID
    pwgcf::PIDSelectionConfigurable pidsel(pidfilter.pidtpcfilter.tpcel, pidfilter.pidtpcfilter.tpcmu, pidfilter.pidtpcfilter.tpcpi, pidfilter.pidtpcfilter.tpcka, pidfilter.pidtpcfilter.tpcpr,
                                           pidfilter.pidtoffilter.tpcel, pidfilter.pidtoffilter.tpcmu, pidfilter.pidtoffilter.tpcpi, pidfilter.pidtoffilter.tpcka, pidfilter.pidtoffilter.tpcpr,
                                           pidfilter.pidbayesfilter.bayel, pidfilter.pidbayesfilter.baymu, pidfilter.pidbayesfilter.baypi, pidfilter.pidbayesfilter.bayka, pidfilter.pidbayesfilter.baypr);
#else
    pwgcf::PIDSelectionConfigurable pidsel(pidfilter.pidtpcfilter.tpcel, pidfilter.pidtpcfilter.tpcmu, pidfilter.pidtpcfilter.tpcpi, pidfilter.pidtpcfilter.tpcka, pidfilter.pidtpcfilter.tpcpr,
                                           pidfilter.pidtoffilter.tpcel, pidfilter.pidtoffilter.tpcmu, pidfilter.pidtoffilter.tpcpi, pidfilter.pidtoffilter.tpcka, pidfilter.pidtoffilter.tpcpr);
#endif
    fFilterFramework = new pwgcf::FilterAndAnalysisFramework(filterccdb.ccdburl.value, filterccdb.ccdbpath.value, filterccdb.filterdate.value);
    fFilterFramework->setConfiguration(eventsel, trksel, pidsel, pwgcf::SelectionFilterAndAnalysis::KAnalysis);
    fFilterFramework->init();

    collisionmask = fFilterFramework->getCollisionMask();
    collisionmaskOptional = fFilterFramework->getCollisionOptMask();
    collisionmaskForced = fFilterFramework->getCollisionForcedMask();
    fMultiplicityIndex = fFilterFramework->getCollisionMultiplicityIndex();
    trackmask = fFilterFramework->getTrackMask();
    trackmaskOptional = fFilterFramework->getTrackOptMask();
    trackmaskForced = fFilterFramework->getTrackForcedMask();
    pidmask = fFilterFramework->getPIDMask();
    pidmaskOptional = fFilterFramework->getPIDOptMask();
    pidmaskForced = fFilterFramework->getPIDForcedMask();
    LOGF(info, "twoParticleCorrelationsFilter::init(), collision selection masks 0x%016lx, %s, and 0x%016lx and multiplicity index %d", collisionmask, fFilterFramework->printCollisionOptionalMasks().Data(), collisionmaskForced, fMultiplicityIndex);
    LOGF(info, "twoParticleCorrelationsFilter::init(), track selection masks 0x%016lx, %s, and 0x%016lx ", trackmask, fFilterFramework->printTrackOptionalMasks().Data(), trackmaskForced);
    LOGF(info, "twoParticleCorrelationsFilter::init(), PID selection masks 0x%016lx, %s, and 0x%016lx ", pidmask, fFilterFramework->printPIDOptionalMasks().Data(), pidmaskForced);
    if (collisionmask == static_cast<uint64_t>(0) || trackmask == static_cast<uint64_t>(0)) {
      LOGF(fatal, "twoParticleCorrelationsFilter::init() null masks, selecting everything!!!");
    }

    /* initial size for the tracks ids bank */
    tracksids.reserve(bankinitialcapacity);
  }

  /// \brief Get the data collecting engine index corresponding to the passed collision
  template <typename FilteredCollision>
  int getDCEindex(const FilteredCollision& collision)
  {
    int ixDCE = -1;
    float cm = collision.centMult()[fMultiplicityIndex];
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
    LOGF(TWOPCORRLOGCOLLISIONS, "Received collision with mask 0x%016lx and %ld tracks", collision.selFlags(), tracks.size());
    auto passOptions = [](const auto& options, auto mask) {
      bool all = true;
      for (const auto& option : options) {
        all = all && ((option & mask) != 0UL);
      }
      return all;
    };

    if ((collision.selFlags() & collisionmaskForced) == collisionmaskForced && passOptions(collisionmaskOptional, collision.selFlags())) {
      LOGF(TWOPCORRLOGCOLLISIONS, ">> Accepted collision with mask 0x%016lx and %ld unfiltered tracks", collision.selFlags(), tracks.size());
      int nAcceptedTracks = 0;
      int nRejectedTracks = 0;
      tracksids.clear();
      if (static_cast<int64_t>(tracksids.capacity()) < tracks.size()) {
        tracksids.reserve(tracksids.capacity() + bankcapacitystep);
        LOGF(warning, "Increased tracks bank size to %d", tracksids.size());
      }
      for (const auto& track : tracks) {
        int trackid = -1;
        if ((track.trackFlags() & trackmaskForced) == trackmaskForced && passOptions(trackmaskOptional, track.trackFlags())) {
          if ((track.trackFlags() & 0x1L) == 0x1L) {
            /* positive track */
            trackid = 0;
            nAcceptedTracks++;
          } else if ((track.trackFlags() & 0x2L) == 0x2L) {
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
        dataCE[ixDCE]->processCollision(tracks, collision.posZ(), collision.centMult()[fMultiplicityIndex], 0);
      }
      LOGF(TWOPCORRLOGCOLLISIONS, ">> Accepted collision with mask 0x%016lx, %d accepted tracks and %d rejected tracks", collision.selFlags(), nAcceptedTracks, nRejectedTracks);
    } else {
    }
  }

  Filter tracksfilter = (aod::track::eta < cfgBinning->mEtamax) && (aod::track::eta > cfgBinning->mEtamin) && (aod::track::pt < cfgBinning->mPTmax) && (aod::track::pt > cfgBinning->mPTmin);

  void processRun2(soa::Join<aod::Collisions, aod::CFCollMasks>::iterator const& collision,
                   soa::Filtered<soa::Join<aod::FullTracks, aod::CFTrackMasks>> const& tracks)
  {
    processTheTask(collision, tracks);
  }
  PROCESS_SWITCH(TwoParticleCorrelations, processRun2, "Process Run 2, i.e. over NOT stored derived data, two particle correlations", true);

  void processRun3(aod::CFCollision const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    processTheTask(collision, tracks);
  }
  PROCESS_SWITCH(TwoParticleCorrelations, processRun3, "Process Run 3, i.e. over stored derived data, two particle correlations", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TwoParticleCorrelations>(cfgc)};
  return workflow;
}
