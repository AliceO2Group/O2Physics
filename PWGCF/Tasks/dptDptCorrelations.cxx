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

/// \file dptDptCorrelations.cxx
/// \brief implements two-particle correlations base data collection
/// \author victor.gonzalez.sebastian@gmail.com

#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/DptDptFiltered.h"
#include "PWGCF/TableProducer/dptDptFilter.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TableHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>

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
#include <cstdio>
#include <ctime>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;

#define DPTDPTLOGCOLLISIONS debug
#define DPTDPTLOGTRACKS debug

namespace correlationstask
{
using namespace o2::analysis::dptdptfilter;
float etabinwidth = (etaup - etalow) / static_cast<float>(etabins);
float phibinwidth = (phiup - philow) / static_cast<float>(phibins);
int deltaetabins = etabins * 2 - 1;
float deltaetalow = etalow - etaup, deltaetaup = etaup - etalow;
float deltaetabinwidth = (deltaetaup - deltaetalow) / static_cast<float>(deltaetabins);
int deltaphibins = phibins;
float deltaphibinwidth = constants::math::TwoPI / deltaphibins;
float deltaphilow = 0.0 - deltaphibinwidth / 2.0;
float deltaphiup = constants::math::TwoPI - deltaphibinwidth / 2.0;

enum HistoDimensions {
  kNONE = 0,
  k1D,
  k2D,
  k3D,
  k4D
};
HistoDimensions nNoOfDimensions = k1D; // number of dimensions for the NUA & NUE corrections
bool processpairs = false;             // process pairs analysis
bool processmixedevents = false;       // process mixed events
bool ptorder = false;                  // consider pt ordering
bool invmass = false;                  // produce the invariant mass histograms
bool corrana = false;                  // produce the correlation analysis histograms

PairCuts fPairCuts;              // pair suppression engine
bool fUseConversionCuts = false; // suppress resonances and conversions
bool fUseTwoTrackCut = false;    // suppress too close tracks

std::vector<std::string> poinames;                     ///< the species of interest names
std::vector<std::string> tnames;                       ///< the track names
std::vector<double> poimass;                           ///< the species of interest mass
std::vector<std::vector<std::string>> trackPairsNames; ///< the track pairs names
} // namespace correlationstask

// Task for building <dpt,dpt> correlations
struct DptDptCorrelations {

  /* the data collecting engine */
  template <bool smallsingles>
  struct DataCollectingEngine {
    size_t nch = correlationstask::tnames.size();

    //============================================================================================
    // The DptDptCorrelationsAnalysisTask output objects
    //============================================================================================
    /* histograms */
    TH1F* fhVertexZA;                                                            //!<! the z vertex distribution for the current multiplicity/centrality class
    std::vector<TH1F*> fhN1VsPt{nch, nullptr};                                   //!<! weighted single particle distribution vs \f$p_T\f$, for the different species
    std::vector<TH2F*> fhN1VsPtEta{nch, nullptr};                                //!<! weighted single particle distribution vs \f$p_T,\;\eta\f$, for the different species
    std::vector<TH2F*> fhN1VsEtaPhi{nch, nullptr};                               //!<! weighted single particle distribution vs \f$\eta,\;\phi\f$, for the different species
    std::vector<TH2F*> fhSum1PtVsEtaPhi{nch, nullptr};                           //!<! accumulated sum of weighted \f$p_T\f$ vs \f$\eta,\;\phi\f$, for the different species
    std::vector<TH3F*> fhN1VsZEtaPhiPt{nch, nullptr};                            //!<! single particle distribution vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$, for the different species
    std::vector<TH3F*> fhSum1PtVsZEtaPhiPt{nch, nullptr};                        //!<! accumulated sum of weighted \f$p_T\f$ vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$, for the different species
    std::vector<TH1*> fhNuaNue{nch, nullptr};                                    //!<! NUA+NUE correction for the differents species
    std::vector<TH2*> fhPtAvgVsEtaPhi{nch, nullptr};                             //!<! average \f$p_T\f$ vs \f$\eta,\;\phi\f$, for the different species
    std::vector<std::vector<TH2F*>> fhN2VsPtPt{nch, {nch, nullptr}};             //!<! weighted two particle distribution vs \f${p_T}_1, {p_T}_2\f$ for the different species combinations
    std::vector<std::vector<TH2F*>> fhN2VsDEtaDPhi{nch, {nch, nullptr}};         //!<! two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$ for the different species combinations
    std::vector<std::vector<TH2F*>> fhN2contVsDEtaDPhi{nch, {nch, nullptr}};     //!<! two-particle distribution continuous vs \f$\Delta\eta,\;\Delta\phi\f$ for the different species combinations
    std::vector<std::vector<TH2F*>> fhSum2PtPtVsDEtaDPhi{nch, {nch, nullptr}};   //!<! two-particle  \f$\sum {p_T}_1 {p_T}_2\f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ for the different species combinations
    std::vector<std::vector<TH2F*>> fhSum2DptDptVsDEtaDPhi{nch, {nch, nullptr}}; //!<! two-particle  \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ for the different species combinations
    std::vector<std::vector<TH2F*>> fhSupN1N1VsDEtaDPhi{nch, {nch, nullptr}};    //!<! suppressed n1n1 two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$ for the different species combinations
    std::vector<std::vector<TH2F*>> fhSupPt1Pt1VsDEtaDPhi{nch, {nch, nullptr}};  //!<! suppressed \f${p_T}_1 {p_T}_2\f$ two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$ for the different species combinations
    std::vector<std::vector<TH2D*>> fhInvMassDEta{nch, {nch, nullptr}};          //!<! the pair invariant mass vs delta eta
    std::vector<std::vector<TH2D*>> fhInvMassDPhi{nch, {nch, nullptr}};          //!<! the pair invariant mass vs delta phi
    /* versus centrality/multiplicity  profiles */
    std::vector<TProfile*> fhN1VsC{nch, nullptr};                               //!<! weighted single particle distribution vs event centrality/multiplicity, track 1 and 2
    std::vector<TProfile*> fhSum1PtVsC{nch, nullptr};                           //!<! accumulated sum of weighted \f$p_T\f$ vs event centrality/multiplicity, track 1 and 2
    std::vector<TProfile*> fhN1nwVsC{nch, nullptr};                             //!<! un-weighted single particle distribution vs event centrality/multiplicity, track 1 and 2
    std::vector<TProfile*> fhSum1PtnwVsC{nch, nullptr};                         //!<! accumulated sum of un-weighted \f$p_T\f$ vs event centrality/multiplicity, track 1 and 2
    std::vector<std::vector<TProfile*>> fhN2VsC{nch, {nch, nullptr}};           //!<! weighted accumulated two particle distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    std::vector<std::vector<TProfile*>> fhSum2PtPtVsC{nch, {nch, nullptr}};     //!<! weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    std::vector<std::vector<TProfile*>> fhSum2DptDptVsC{nch, {nch, nullptr}};   //!<! weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    std::vector<std::vector<TProfile*>> fhN2nwVsC{nch, {nch, nullptr}};         //!<! un-weighted accumulated two particle distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    std::vector<std::vector<TProfile*>> fhSum2PtPtnwVsC{nch, {nch, nullptr}};   //!<! un-weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations
    std::vector<std::vector<TProfile*>> fhSum2DptDptnwVsC{nch, {nch, nullptr}}; //!<! un-weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ distribution vs event centrality/multiplicity 1-1,1-2,2-1,2-2, combinations

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
      using namespace correlationstask;
      using namespace o2::analysis::dptdptfilter;
      return RecoDecay::constrainAngle(phi, philow);
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
      using namespace correlationstask;
      using namespace o2::analysis::dptdptfilter;

      int etaix = static_cast<int>((t.eta() - etalow) / etabinwidth);
      /* consider a potential phi origin shift */
      float phi = getShiftedPhi(t.phi());
      int phiix = static_cast<int>((phi - philow) / phibinwidth);
      return etaix * phibins + phiix;
    }

    /// \brief Returns the delta eta value for the differential eta
    /// \param t1 the intended track one
    /// \param t2 the intended track two
    /// \return the delta eta value for delta eta
    ///
    /// WARNING: for performance reasons no checks are done about the consistency
    /// of tracks' eta and phi within the corresponding ranges so, it is suppossed
    /// the tracks have been accepted and they are within that ranges
    /// IF THAT IS NOT THE CASE THE ROUTINE WILL PRODUCE NONSENSE RESULTS
    template <typename TrackObject>
    float getDEtaValue(TrackObject const& t1, TrackObject const& t2)
    {
      using namespace correlationstask;
      using namespace o2::analysis::dptdptfilter;

      /* rule: ix are always zero based while bins are always one based */
      int etaIx1 = static_cast<int>((t1.eta() - etalow) / etabinwidth);
      int etaIx2 = static_cast<int>((t2.eta() - etalow) / etabinwidth);

      int deltaEtaIx = etaIx1 - etaIx2 + etabins - 1;

      return deltaetalow + (deltaEtaIx + 0.5) * deltaetabinwidth;
    }

    /// \brief Returns the delta phi value for the differential phi
    /// \param t1 the intended track one
    /// \param t2 the intended track two
    /// \return the delta phi value within [-pi,pi] for delta phi
    ///
    /// WARNING: for performance reasons no checks are done about the consistency
    /// of tracks' eta and phi within the corresponding ranges so, it is suppossed
    /// the tracks have been accepted and they are within that ranges
    /// IF THAT IS NOT THE CASE THE ROUTINE WILL PRODUCE NONSENSE RESULTS
    template <typename TrackObject>
    float getDPhiValue(TrackObject const& t1, TrackObject const& t2)
    {
      using namespace correlationstask;
      using namespace o2::analysis::dptdptfilter;

      /* rule: ix are always zero based while bins are always one based */
      /* consider a potential phi origin shift */
      float phi = getShiftedPhi(t1.phi());
      int phiIx1 = static_cast<int>((phi - philow) / phibinwidth);
      /* consider a potential phi origin shift */
      phi = getShiftedPhi(t2.phi());
      int phiIx2 = static_cast<int>((phi - philow) / phibinwidth);

      int deltaPhiIx = phiIx1 - phiIx2;
      if (deltaPhiIx < 0) {
        deltaPhiIx += phibins;
      }

      float value = deltaphilow + (deltaPhiIx + 0.5) * deltaphibinwidth;

      return RecoDecay::constrainAngle(value, deltaphilow - constants::math::PI);
    }

    /// \brief Returns the TH2 global bin for the differential histograms
    /// \param t1 the intended track one
    /// \param t2 the intended track two
    /// \return the globl TH2 bin for delta eta delta phi
    ///
    /// WARNING: for performance reasons no checks are done about the consistency
    /// of tracks' eta and phi within the corresponding ranges so, it is suppossed
    /// the tracks have been accepted and they are within that ranges
    /// IF THAT IS NOT THE CASE THE ROUTINE WILL PRODUCE NONSENSE RESULTS
    template <typename TrackObject>
    int getDEtaDPhiGlobalBin(TrackObject const& t1, TrackObject const& t2)
    {
      using namespace correlationstask;
      using namespace o2::analysis::dptdptfilter;

      /* rule: ix are always zero based while bins are always one based */
      int etaIx1 = static_cast<int>((t1.eta() - etalow) / etabinwidth);
      /* consider a potential phi origin shift */
      float phi = getShiftedPhi(t1.phi());
      int phiIx1 = static_cast<int>((phi - philow) / phibinwidth);
      int etaIx2 = static_cast<int>((t2.eta() - etalow) / etabinwidth);
      /* consider a potential phi origin shift */
      phi = getShiftedPhi(t2.phi());
      int phiIx2 = static_cast<int>((phi - philow) / phibinwidth);

      int deltaEtaIx = etaIx1 - etaIx2 + etabins - 1;
      int deltaPhiIx = phiIx1 - phiIx2;
      if (deltaPhiIx < 0) {
        deltaPhiIx += phibins;
      }

      return fhN2VsDEtaDPhi[0][0]->GetBin(deltaEtaIx + 1, deltaPhiIx + 1);
    }

    /* taken from PWGCF/Core/PairCuts.h implemented by JFGO */
    template <typename TrackObject>
    double getInvMassSquared(TrackObject const& track1, double m0_1, TrackObject const& track2, double m0_2)
    {
      // calculate inv mass squared
      // same can be achieved, but with more computing time with
      /*TLorentzVector photon, p1, p2;
      p1.SetPtEtaPhiM(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), 0.510e-3);
      p2.SetPtEtaPhiM(particle->Pt(), eta[j], particle->Phi(), 0.510e-3);
      photon = p1+p2;
      photon.M()*/

      constexpr float kLARGETANTHETA = 1e10;
      constexpr float kVERYSMALLETA = 1e-10;
      float tantheta1 = kLARGETANTHETA;

      if (track1.eta() < -kVERYSMALLETA || track1.eta() > kVERYSMALLETA) {
        float expTmp = std::exp(-track1.eta());
        tantheta1 = 2.0 * expTmp / (1.0 - expTmp * expTmp);
      }

      float tantheta2 = kLARGETANTHETA;
      if (track2.eta() < -kVERYSMALLETA || track2.eta() > kVERYSMALLETA) {
        float expTmp = std::exp(-track2.eta());
        tantheta2 = 2.0 * expTmp / (1.0 - expTmp * expTmp);
      }

      float e1squ = m0_1 * m0_1 + track1.pt() * track1.pt() * (1.0 + 1.0 / tantheta1 / tantheta1);
      float e2squ = m0_2 * m0_2 + track2.pt() * track2.pt() * (1.0 + 1.0 / tantheta2 / tantheta2);

      float mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * (std::sqrt(e1squ * e2squ) - (track1.pt() * track2.pt() * (std::cos(track1.phi() - track2.phi()) + 1.0 / tantheta1 / tantheta2)));

      return mass2;
    }

    void storeTrackCorrections(std::vector<TH1*> corrs)
    {
      using namespace correlationstask;

      LOGF(info, "Storing NUA&NUE corrections for %d track ids", corrs.size());
      for (uint i = 0; i < corrs.size(); ++i) {
        if (corrs[i] != nullptr) {
          if (nNoOfDimensions != corrs[i]->GetDimension()) {
            LOGF(fatal, "  Corrections received dimensions %d for track id %d different than expected %d", corrs[i]->GetDimension(), i, static_cast<int>(nNoOfDimensions));
          } else {
            LOGF(info, "  Storing NUA&NUE corrections %s for track id %d with %d dimensions %s",
                 corrs[i] != nullptr ? corrs[i]->GetName() : "nullptr", i, static_cast<int>(nNoOfDimensions), corrs[i] != nullptr ? "yes" : "no");
          }
        }
        fhNuaNue[i] = corrs[i];
        if (fhNuaNue[i] != nullptr) {
          int nbins = 0;
          double avg = 0.0;
          for (int ix = 0; ix < fhNuaNue[i]->GetNbinsX(); ++ix) {
            if (nNoOfDimensions == k1D) {
              nbins++;
              avg += fhNuaNue[i]->GetBinContent(ix + 1);
            } else {
              for (int iy = 0; iy < fhNuaNue[i]->GetNbinsY(); ++iy) {
                if (nNoOfDimensions == k2D) {
                  nbins++;
                  avg += fhNuaNue[i]->GetBinContent(ix + 1, iy + 1);
                } else if (nNoOfDimensions == k3D || nNoOfDimensions == k4D) {
                  for (int iz = 0; iz < fhNuaNue[i]->GetNbinsZ(); ++iz) {
                    nbins++;
                    avg += fhNuaNue[i]->GetBinContent(ix + 1, iy + 1, iz + 1);
                  }
                }
              }
            }
          }
          LOGF(info, "  Average NUA&NUE correction for track id %d: %f", i, avg / nbins);
        }
      }
      ccdbstored = true;
    }

    void storePtAverages(std::vector<TH2*> ptavgs)
    {
      LOGF(info, "Stored pT average for %d track ids", ptavgs.size());
      for (uint i = 0; i < ptavgs.size(); ++i) {
        LOGF(info, "  Stored pT average for track id %d %s", i, ptavgs[i] != nullptr ? "yes" : "no");
        fhPtAvgVsEtaPhi[i] = ptavgs[i];
      }
      ccdbstored = true;
    }

    template <correlationstask::HistoDimensions nDim, typename TrackListObject>
    std::vector<float>* getTrackCorrections(TrackListObject const& tracks, float zvtx)
    {
      using namespace correlationstask;

      std::vector<float>* corr = new std::vector<float>(tracks.size(), 1.0f);
      int index = 0;
      for (const auto& t : tracks) {
        if (fhNuaNue[t.trackacceptedid()] != nullptr) {
          if constexpr (nDim == k1D) {
            (*corr)[index] = fhNuaNue[t.trackacceptedid()]->GetBinContent(fhNuaNue[t.trackacceptedid()]->FindFixBin(t.pt()));
          } else if constexpr (nDim == k2D) {
            (*corr)[index] = fhNuaNue[t.trackacceptedid()]->GetBinContent(fhNuaNue[t.trackacceptedid()]->FindFixBin(t.eta(), t.pt()));
          } else if constexpr (nDim == k3D) {
            (*corr)[index] = fhNuaNue[t.trackacceptedid()]->GetBinContent(fhNuaNue[t.trackacceptedid()]->FindFixBin(zvtx, getEtaPhiIndex(t) + 0.5, t.pt()));
          }
        }
        index++;
      }
      return corr;
    }

    template <typename TrackListObject>
    std::vector<float>* getTrackCorrections(TrackListObject const& tracks, float zvtx)
    {
      using namespace correlationstask;

      if (nNoOfDimensions == k1D) {
        return getTrackCorrections<k1D>(tracks, zvtx);
      } else if (nNoOfDimensions == k2D) {
        return getTrackCorrections<k2D>(tracks, zvtx);
      } else if (nNoOfDimensions == k3D) {
        return getTrackCorrections<k3D>(tracks, zvtx);
      }
      return getTrackCorrections<k4D>(tracks, zvtx);
    }

    template <typename TrackListObject>
    std::vector<float>* getPtAvg(TrackListObject const& tracks)
    {
      std::vector<float>* ptavg = new std::vector<float>(tracks.size(), 0.0f);
      int index = 0;
      for (auto const& t : tracks) {
        if (fhPtAvgVsEtaPhi[t.trackacceptedid()] != nullptr) {
          (*ptavg)[index] = fhPtAvgVsEtaPhi[t.trackacceptedid()]->GetBinContent(fhPtAvgVsEtaPhi[t.trackacceptedid()]->FindFixBin(t.eta(), t.phi()));
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
      for (auto const& track : passedtracks) {
        float corr = (*corrs)[index];
        fhN1VsPt[track.trackacceptedid()]->Fill(track.pt(), corr);
        if constexpr (smallsingles) {
          fhN1VsPtEta[track.trackacceptedid()]->Fill(track.eta(), track.pt(), corr);
          fhN1VsEtaPhi[track.trackacceptedid()]->Fill(track.eta(), getShiftedPhi(track.phi()), corr);
          fhSum1PtVsEtaPhi[track.trackacceptedid()]->Fill(track.eta(), getShiftedPhi(track.phi()), track.pt() * corr);
        } else {
          fhN1VsZEtaPhiPt[track.trackacceptedid()]->Fill(zvtx, getEtaPhiIndex(track) + 0.5, track.pt(), corr);
          fhSum1PtVsZEtaPhiPt[track.trackacceptedid()]->Fill(zvtx, getEtaPhiIndex(track) + 0.5, track.pt(), track.pt() * corr);
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
      for (auto const& track : passedtracks) {
        float corr = (*corrs)[index];
        n1[track.trackacceptedid()] += corr;
        sum1Pt[track.trackacceptedid()] += track.pt() * corr;
        n1nw[track.trackacceptedid()] += 1;
        sum1Ptnw[track.trackacceptedid()] += track.pt();

        fhN1VsPtEta[track.trackacceptedid()]->Fill(track.eta(), track.pt(), corr);
        fhN1VsEtaPhi[track.trackacceptedid()]->Fill(track.eta(), getShiftedPhi(track.phi()), corr);
        fhSum1PtVsEtaPhi[track.trackacceptedid()]->Fill(track.eta(), getShiftedPhi(track.phi()), track.pt() * corr);
        index++;
      }
      for (uint tid = 0; tid < nch; ++tid) {
        fhN1VsC[tid]->Fill(cmul, n1[tid]);
        fhSum1PtVsC[tid]->Fill(cmul, sum1Pt[tid]);
        fhN1nwVsC[tid]->Fill(cmul, n1nw[tid]);
        fhSum1PtnwVsC[tid]->Fill(cmul, sum1Ptnw[tid]);
      }
    }

    /// \brief fills the pair histograms in pair execution mode
    /// \param trks1 filtered table with the tracks associated to the first track in the pair
    /// \param trks2 filtered table with the tracks associated to the second track in the pair
    /// \param cmul centrality - multiplicity for the collision being analyzed
    /// Be aware that in most of the cases traks1 and trks2 will have the same content (exception: mixed events)
    template <bool doptorder, bool doinvmass, bool docorrelations, typename TrackOneListObject, typename TrackTwoListObject>
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
      int globalbin = 0;
      LOGF(debug, "Initializing globalbin to ", globalbin);

      for (auto const& track1 : trks1) {
        double ptAvg1 = (*ptavgs1)[index1];
        double corr1 = (*corrs1)[index1];
        int index2 = 0;
        for (auto const& track2 : trks2) {
          /* checking the same track id condition */
          if (track1 == track2) {
            /* exclude autocorrelations */
            index2++;
            continue;
          }

          if constexpr (doptorder) {
            if (track2.pt() >= track1.pt()) {
              index2++;
              continue;
            }
          }
          /* process pair magnitudes */
          double ptAvg2 = (*ptavgs2)[index2];
          double corr2 = (*corrs2)[index2];
          double corr = corr1 * corr2;
          double dptdptnw = (track1.pt() - ptAvg1) * (track2.pt() - ptAvg2);
          double dptdptw = (corr1 * track1.pt() - ptAvg1) * (corr2 * track2.pt() - ptAvg2);

          /* get the global bin for filling the differential histograms */
          if constexpr (docorrelations) {
            globalbin = getDEtaDPhiGlobalBin(track1, track2);
          }
          float deltaeta = track1.eta() - track2.eta();
          float deltaphi = track1.phi() - track2.phi();
          deltaphi = RecoDecay::constrainAngle(deltaphi, deltaphilow);
          if ((fUseConversionCuts && fPairCuts.conversionCuts(track1, track2)) || (fUseTwoTrackCut && fPairCuts.twoTrackCut(track1, track2, bfield))) {
            /* suppress the pair */
            if constexpr (docorrelations) {
              fhSupN1N1VsDEtaDPhi[track1.trackacceptedid()][track2.trackacceptedid()]->AddBinContent(globalbin, corr);
              fhSupPt1Pt1VsDEtaDPhi[track1.trackacceptedid()][track2.trackacceptedid()]->AddBinContent(globalbin, track1.pt() * track2.pt() * corr);
            }
            n2sup[track1.trackacceptedid()][track2.trackacceptedid()] += corr;
          } else {
            /* count the pair */
            n2[track1.trackacceptedid()][track2.trackacceptedid()] += corr;
            sum2PtPt[track1.trackacceptedid()][track2.trackacceptedid()] += track1.pt() * track2.pt() * corr;
            sum2DptDpt[track1.trackacceptedid()][track2.trackacceptedid()] += dptdptw;
            n2nw[track1.trackacceptedid()][track2.trackacceptedid()] += 1;
            sum2PtPtnw[track1.trackacceptedid()][track2.trackacceptedid()] += track1.pt() * track2.pt();
            sum2DptDptnw[track1.trackacceptedid()][track2.trackacceptedid()] += dptdptnw;

            if constexpr (docorrelations) {
              fhN2VsDEtaDPhi[track1.trackacceptedid()][track2.trackacceptedid()]->AddBinContent(globalbin, corr);
              fhN2contVsDEtaDPhi[track1.trackacceptedid()][track2.trackacceptedid()]->Fill(deltaeta, deltaphi, corr);
              fhSum2DptDptVsDEtaDPhi[track1.trackacceptedid()][track2.trackacceptedid()]->AddBinContent(globalbin, dptdptw);
              fhSum2PtPtVsDEtaDPhi[track1.trackacceptedid()][track2.trackacceptedid()]->AddBinContent(globalbin, track1.pt() * track2.pt() * corr);
              fhN2VsPtPt[track1.trackacceptedid()][track2.trackacceptedid()]->Fill(track1.pt(), track2.pt(), corr);
            }
            if constexpr (doinvmass) {
              if (!(track2.trackacceptedid() < track1.trackacceptedid())) {
                /* only 12 combinations, 21 are exactly the same */
                double invariantMass = std::sqrt(getInvMassSquared(track1, poimass[static_cast<int>(track1.trackacceptedid() / 2)], track2, poimass[static_cast<int>(track2.trackacceptedid() / 2)])) * 1000.0f;
                fhInvMassDEta[track1.trackacceptedid()][track2.trackacceptedid()]->Fill(getDEtaValue(track1, track2), invariantMass);
                fhInvMassDPhi[track1.trackacceptedid()][track2.trackacceptedid()]->Fill(getDPhiValue(track1, track2), invariantMass);
              }
            }
          }
          index2++;
        }
        index1++;
      }
      for (uint pid1 = 0; pid1 < nch; ++pid1) {
        for (uint pid2 = 0; pid2 < nch; ++pid2) {
          fhN2VsC[pid1][pid2]->Fill(cmul, n2[pid1][pid2]);
          fhSum2PtPtVsC[pid1][pid2]->Fill(cmul, sum2PtPt[pid1][pid2]);
          fhSum2DptDptVsC[pid1][pid2]->Fill(cmul, sum2DptDpt[pid1][pid2]);
          fhN2nwVsC[pid1][pid2]->Fill(cmul, n2nw[pid1][pid2]);
          fhSum2PtPtnwVsC[pid1][pid2]->Fill(cmul, sum2PtPtnw[pid1][pid2]);
          fhSum2DptDptnwVsC[pid1][pid2]->Fill(cmul, sum2DptDptnw[pid1][pid2]);
          /* let's also update the number of entries in the differential histograms */
          if constexpr (docorrelations) {
            fhN2VsDEtaDPhi[pid1][pid2]->SetEntries(fhN2VsDEtaDPhi[pid1][pid2]->GetEntries() + n2[pid1][pid2]);
            fhSum2DptDptVsDEtaDPhi[pid1][pid2]->SetEntries(fhSum2DptDptVsDEtaDPhi[pid1][pid2]->GetEntries() + n2[pid1][pid2]);
            fhSum2PtPtVsDEtaDPhi[pid1][pid2]->SetEntries(fhSum2PtPtVsDEtaDPhi[pid1][pid2]->GetEntries() + n2[pid1][pid2]);
            fhSupN1N1VsDEtaDPhi[pid1][pid2]->SetEntries(fhSupN1N1VsDEtaDPhi[pid1][pid2]->GetEntries() + n2sup[pid1][pid2]);
            fhSupPt1Pt1VsDEtaDPhi[pid1][pid2]->SetEntries(fhSupPt1Pt1VsDEtaDPhi[pid1][pid2]->GetEntries() + n2sup[pid1][pid2]);
          }
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
            /* no invariant mass analysis on a mixed event data collection */
            processTrackPairs<true, false, true>(Tracks1, Tracks2, corrs1, corrs2, ptavgs1, ptavgs2, centmult, bfield);
          } else {
            processTrackPairs<false, false, true>(Tracks1, Tracks2, corrs1, corrs2, ptavgs1, ptavgs2, centmult, bfield);
          }
        } else {
          if (ptorder) {
            if (invmass) {
              if (corrana) {
                processTrackPairs<true, true, true>(Tracks1, Tracks1, corrs1, corrs1, ptavgs1, ptavgs1, centmult, bfield);
              } else {
                processTrackPairs<true, true, false>(Tracks1, Tracks1, corrs1, corrs1, ptavgs1, ptavgs1, centmult, bfield);
              }
            } else {
              if (corrana) {
                processTrackPairs<true, false, true>(Tracks1, Tracks1, corrs1, corrs1, ptavgs1, ptavgs1, centmult, bfield);
              } else {
                processTrackPairs<true, false, false>(Tracks1, Tracks1, corrs1, corrs1, ptavgs1, ptavgs1, centmult, bfield);
              }
            }
          } else {
            if (invmass) {
              if (corrana) {
                processTrackPairs<false, true, true>(Tracks1, Tracks1, corrs1, corrs1, ptavgs1, ptavgs1, centmult, bfield);
              } else {
                processTrackPairs<false, true, false>(Tracks1, Tracks1, corrs1, corrs1, ptavgs1, ptavgs1, centmult, bfield);
              }
            } else {
              if (corrana) {
                processTrackPairs<false, false, true>(Tracks1, Tracks1, corrs1, corrs1, ptavgs1, ptavgs1, centmult, bfield);
              } else {
                processTrackPairs<false, false, false>(Tracks1, Tracks1, corrs1, corrs1, ptavgs1, ptavgs1, centmult, bfield);
              }
            }
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

    template <bool doinvmass, bool docorrelations>
    void init(TList* fOutputList)
    {
      using namespace correlationstask;
      using namespace o2::analysis::dptdptfilter;

      LOGF(info, "Do invariant mass: %s; do correlation histograms: %s", doinvmass ? "yes" : "no", docorrelations ? "yes" : "no");

      /* create the histograms */
      bool oldstatus = TH1::AddDirectoryStatus();
      TH1::AddDirectory(kFALSE);

      if (!processpairs) {
        fhVertexZA = new TH1F("VertexZA", "Vertex Z; z_{vtx}", zvtxbins, zvtxlow, zvtxup);
        fOutputList->Add(fhVertexZA);
        for (uint i = 0; i < nch; ++i) {
          /* histograms for each track, one and two */
          fhN1VsPt[i] = new TH1F(TString::Format("n1_%s_vsPt", tnames[i].c_str()).Data(),
                                 TString::Format("#LT n_{1} #GT;p_{t,%s} (GeV/c);#LT n_{1} #GT", tnames[i].c_str()).Data(),
                                 ptbins, ptlow, ptup);
          /* we don't want the Sumw2 structure being created here */
          bool defSumw2 = TH1::GetDefaultSumw2();
          if constexpr (smallsingles) {
            fhN1VsPtEta[i] = new TH2F(TString::Format("n1_%s_vsPtEta", tnames[i].c_str()).Data(),
                                      TString::Format("#LT n_{1_{%s}} #GT;#eta;#it{p}_{T}", tnames[i].c_str()).Data(),
                                      etabins, etalow, etaup, ptbins, ptlow, ptup);
            fhN1VsEtaPhi[i] = new TH2F(TString::Format("n1_%s_vsEtaPhi", tnames[i].c_str()).Data(),
                                       TString::Format("#LT n_{1} #GT;#eta_{%s};#varphi_{%s} (radian);#LT n_{1} #GT", tnames[i].c_str(), tnames[i].c_str()).Data(),
                                       etabins, etalow, etaup, phibins, philow, phiup);
            fhSum1PtVsEtaPhi[i] = new TH2F(TString::Format("sumPt_%s_vsEtaPhi", tnames[i].c_str()).Data(),
                                           TString::Format("#LT #Sigma p_{t,%s} #GT;#eta_{%s};#varphi_{%s} (radian);#LT #Sigma p_{t,%s} #GT (GeV/c)",
                                                           tnames[i].c_str(), tnames[i].c_str(), tnames[i].c_str(), tnames[i].c_str())
                                             .Data(),
                                           etabins, etalow, etaup, phibins, philow, phiup);
          } else {
            TH1::SetDefaultSumw2(false);
            fhN1VsZEtaPhiPt[i] = new TH3F(
              TString::Format("n1_%s_vsZ_vsEtaPhi_vsPt", tnames[i].c_str()).Data(),
              TString::Format("#LT n_{1} #GT;vtx_{z};#eta_{%s}#times#varphi_{%s};p_{t,%s} (GeV/c)",
                              tnames[i].c_str(),
                              tnames[i].c_str(),
                              tnames[i].c_str())
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
            fhSum1PtVsZEtaPhiPt[i] = new TH3F(
              TString::Format("sumPt1_%s_vsZ_vsEtaPhi_vsPt", tnames[i].c_str()).Data(),
              TString::Format(
                "#LT #Sigma p_{t,%s}#GT;vtx_{z};#eta_{%s}#times#varphi_{%s};p_{t,%s} (GeV/c)",
                tnames[i].c_str(),
                tnames[i].c_str(),
                tnames[i].c_str(),
                tnames[i].c_str())
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
            fhN1VsPtEta[i]->SetBit(TH1::kIsNotW);
            fhN1VsPtEta[i]->Sumw2(false);
            fhN1VsEtaPhi[i]->SetBit(TH1::kIsNotW);
            fhN1VsEtaPhi[i]->Sumw2(false);
            fhSum1PtVsEtaPhi[i]->SetBit(TH1::kIsNotW);
            fhSum1PtVsEtaPhi[i]->Sumw2(false);
          } else {
            fhN1VsZEtaPhiPt[i]->SetBit(TH1::kIsNotW);
            fhN1VsZEtaPhiPt[i]->Sumw2(false);
            fhSum1PtVsZEtaPhiPt[i]->SetBit(TH1::kIsNotW);
            fhSum1PtVsZEtaPhiPt[i]->Sumw2(false);
          }
          fhNuaNue[i] = nullptr;
          fhPtAvgVsEtaPhi[i] = nullptr;

          fOutputList->Add(fhN1VsPt[i]);
          if constexpr (smallsingles) {
            fOutputList->Add(fhN1VsPtEta[i]);
            fOutputList->Add(fhN1VsEtaPhi[i]);
            fOutputList->Add(fhSum1PtVsEtaPhi[i]);
          } else {
            fOutputList->Add(fhN1VsZEtaPhiPt[i]);
            fOutputList->Add(fhSum1PtVsZEtaPhiPt[i]);
          }
        }
      } else {
        for (uint i = 0; i < nch; ++i) {
          /* histograms for each track species */
          fhN1VsPtEta[i] = new TH2F(TString::Format("n1_%s_vsPtEta", tnames[i].c_str()).Data(),
                                    TString::Format("#LT n_{1_{%s}} #GT;#eta;#it{p}_{T}", tnames[i].c_str()).Data(),
                                    etabins, etalow, etaup, ptbins, ptlow, ptup);
          fhN1VsEtaPhi[i] = new TH2F(TString::Format("n1_%s_vsEtaPhi", tnames[i].c_str()).Data(),
                                     TString::Format("#LT n_{1} #GT;#eta_{%s};#varphi_{%s} (radian);#LT n_{1} #GT", tnames[i].c_str(), tnames[i].c_str()).Data(),
                                     etabins, etalow, etaup, phibins, philow, phiup);
          fhSum1PtVsEtaPhi[i] = new TH2F(TString::Format("sumPt_%s_vsEtaPhi", tnames[i].c_str()).Data(),
                                         TString::Format("#LT #Sigma p_{t,%s} #GT;#eta_{%s};#varphi_{%s} (radian);#LT #Sigma p_{t,%s} #GT (GeV/c)",
                                                         tnames[i].c_str(), tnames[i].c_str(), tnames[i].c_str(), tnames[i].c_str())
                                           .Data(),
                                         etabins, etalow, etaup, phibins, philow, phiup);
          fhN1VsC[i] = new TProfile(TString::Format("n1_%s_vsM", tnames[i].c_str()).Data(),
                                    TString::Format("#LT n_{1} #GT (weighted);Centrality/Multiplicity (%%);#LT n_{1} #GT").Data(),
                                    100, 0.0, 100.0);
          fhSum1PtVsC[i] = new TProfile(TString::Format("sumPt_%s_vsM", tnames[i].c_str()),
                                        TString::Format("#LT #Sigma p_{t,%s} #GT (weighted);Centrality/Multiplicity (%%);#LT #Sigma p_{t,%s} #GT (GeV/c)", tnames[i].c_str(), tnames[i].c_str()).Data(),
                                        100, 0.0, 100.0);
          fhN1nwVsC[i] = new TProfile(TString::Format("n1Nw_%s_vsM", tnames[i].c_str()).Data(),
                                      TString::Format("#LT n_{1} #GT;Centrality/Multiplicity (%%);#LT n_{1} #GT").Data(),
                                      100, 0.0, 100.0);
          fhSum1PtnwVsC[i] = new TProfile(TString::Format("sumPtNw_%s_vsM", tnames[i].c_str()).Data(),
                                          TString::Format("#LT #Sigma p_{t,%s} #GT;Centrality/Multiplicity (%%);#LT #Sigma p_{t,%s} #GT (GeV/c)", tnames[i].c_str(), tnames[i].c_str()).Data(), 100, 0.0, 100.0);
          fhNuaNue[i] = nullptr;
          fhPtAvgVsEtaPhi[i] = nullptr;
          fOutputList->Add(fhN1VsPtEta[i]);
          fOutputList->Add(fhN1VsEtaPhi[i]);
          fOutputList->Add(fhSum1PtVsEtaPhi[i]);
          fOutputList->Add(fhN1VsC[i]);
          fOutputList->Add(fhSum1PtVsC[i]);
          fOutputList->Add(fhN1nwVsC[i]);
          fOutputList->Add(fhSum1PtnwVsC[i]);
        }

        for (uint i = 0; i < nch; ++i) {
          for (uint j = 0; j < nch; ++j) {
            /* histograms for each track pair combination */
            /* we don't want the Sumw2 structure being created here */
            bool defSumw2 = TH1::GetDefaultSumw2();
            TH1::SetDefaultSumw2(false);
            const char* pname = trackPairsNames[i][j].c_str();
            if constexpr (docorrelations) {
              fhN2VsDEtaDPhi[i][j] = new TH2F(TString::Format("n2_12_vsDEtaDPhi_%s", pname), TString::Format("#LT n_{2} #GT (%s);#Delta#eta;#Delta#varphi;#LT n_{2} #GT", pname),
                                              deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
              fhN2contVsDEtaDPhi[i][j] = new TH2F(TString::Format("n2_12cont_vsDEtaDPhi_%s", pname), TString::Format("#LT n_{2} #GT (%s);#Delta#eta;#Delta#varphi;#LT n_{2} #GT", pname),
                                                  deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
              fhSum2PtPtVsDEtaDPhi[i][j] = new TH2F(TString::Format("sumPtPt_12_vsDEtaDPhi_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s);#Delta#eta;#Delta#varphi;#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname),
                                                    deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
              fhSum2DptDptVsDEtaDPhi[i][j] = new TH2F(TString::Format("sumDptDpt_12_vsDEtaDPhi_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s);#Delta#eta;#Delta#varphi;#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname),
                                                      deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
              fhSupN1N1VsDEtaDPhi[i][j] = new TH2F(TString::Format("suppn1n1_12_vsDEtaDPhi_%s", pname), TString::Format("Suppressed #LT n_{1} #GT#LT n_{1} #GT (%s);#Delta#eta;#Delta#varphi;#LT n_{1} #GT#LT n_{1} #GT", pname),
                                                   deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
              fhSupPt1Pt1VsDEtaDPhi[i][j] = new TH2F(TString::Format("suppPtPt_12_vsDEtaDPhi_%s", pname), TString::Format("Suppressed #LT p_{t,1} #GT#LT p_{t,2} #GT (%s);#Delta#eta;#Delta#varphi;#LT p_{t,1} #GT#LT p_{t,2} #GT (GeV^{2})", pname),
                                                     deltaetabins, deltaetalow, deltaetaup, deltaphibins, deltaphilow, deltaphiup);
              fhN2VsPtPt[i][j] = new TH2F(TString::Format("n2_12_vsPtVsPt_%s", pname), TString::Format("#LT n_{2} #GT (%s);p_{t,1} (GeV/c);p_{t,2} (GeV/c);#LT n_{2} #GT", pname),
                                          ptbins, ptlow, ptup, ptbins, ptlow, ptup);
            }
            if constexpr (doinvmass) {
              if (!(j < i)) {
                /* only 12 combinations, 21 are exactly the same */
                fhInvMassDEta[i][j] = new TH2D(TString::Format("n2_invMassDeta_%s", pname), TString::Format("%s invariant mass;#Delta#eta;Mass (MeV/#it{c}^{2})", pname),
                                               deltaetabins, deltaetalow, deltaetaup, 5000, 0, 5000);
                fhInvMassDPhi[i][j] = new TH2D(TString::Format("n2_invMassDphi_%s", pname), TString::Format("%s invariant mass;#Delta#varphi;Mass (MeV/#it{c}^{2})", pname),
                                               deltaphibins, deltaphilow - constants::math::PI, deltaphiup - constants::math::PI, 5000, 0, 5000);
              }
            }
            /* we return it back to previuos state */
            TH1::SetDefaultSumw2(defSumw2);

            fhN2VsC[i][j] = new TProfile(TString::Format("n2_12_vsM_%s", pname), TString::Format("#LT n_{2} #GT (%s) (weighted);Centrality/Multiplicity (%%);#LT n_{2} #GT", pname), 100, 0.0, 100.0);
            fhSum2PtPtVsC[i][j] = new TProfile(TString::Format("sumPtPt_12_vsM_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s) (weighted);Centrality/Multiplicity (%%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname), 100, 0.0, 100.0);
            fhSum2DptDptVsC[i][j] = new TProfile(TString::Format("sumDptDpt_12_vsM_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s) (weighted);Centrality/Multiplicity (%%);#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname), 100, 0.0, 100.0);
            fhN2nwVsC[i][j] = new TProfile(TString::Format("n2Nw_12_vsM_%s", pname), TString::Format("#LT n_{2} #GT (%s);Centrality/Multiplicity (%%);#LT n_{2} #GT", pname), 100, 0.0, 100.0);
            fhSum2PtPtnwVsC[i][j] = new TProfile(TString::Format("sumPtPtNw_12_vsM_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s);Centrality/Multiplicity (%%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname), 100, 0.0, 100.0);
            fhSum2DptDptnwVsC[i][j] = new TProfile(TString::Format("sumDptDptNw_12_vsM_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s);Centrality/Multiplicity (%%);#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname), 100, 0.0, 100.0);

            /* the statistical uncertainties will be estimated by the subsamples method so let's get rid of the error tracking */
            if constexpr (docorrelations) {
              fhN2VsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
              fhN2VsDEtaDPhi[i][j]->Sumw2(false);
              fhN2contVsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
              fhN2contVsDEtaDPhi[i][j]->Sumw2(false);
              fhSum2PtPtVsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
              fhSum2PtPtVsDEtaDPhi[i][j]->Sumw2(false);
              fhSum2DptDptVsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
              fhSum2DptDptVsDEtaDPhi[i][j]->Sumw2(false);
              fhSupN1N1VsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
              fhSupN1N1VsDEtaDPhi[i][j]->Sumw2(false);
              fhSupPt1Pt1VsDEtaDPhi[i][j]->SetBit(TH1::kIsNotW);
              fhSupPt1Pt1VsDEtaDPhi[i][j]->Sumw2(false);
            }
            if constexpr (doinvmass) {
              if (!(j < i)) {
                /* only 12 combinations, 21 are exactly the same */
                fhInvMassDEta[i][j]->SetBit(TH1::kIsNotW);
                fhInvMassDEta[i][j]->Sumw2(false);
                fhInvMassDPhi[i][j]->SetBit(TH1::kIsNotW);
                fhInvMassDPhi[i][j]->Sumw2(false);
              }
            }

            if constexpr (docorrelations) {
              fOutputList->Add(fhN2VsDEtaDPhi[i][j]);
              fOutputList->Add(fhN2contVsDEtaDPhi[i][j]);
              fOutputList->Add(fhSum2PtPtVsDEtaDPhi[i][j]);
              fOutputList->Add(fhSum2DptDptVsDEtaDPhi[i][j]);
              fOutputList->Add(fhSupN1N1VsDEtaDPhi[i][j]);
              fOutputList->Add(fhSupPt1Pt1VsDEtaDPhi[i][j]);
              fOutputList->Add(fhN2VsPtPt[i][j]);
            }
            if constexpr (doinvmass) {
              if (!(j < i)) {
                /* only 12 combinations, 21 are exactly the same */
                fOutputList->Add(fhInvMassDEta[i][j]);
                fOutputList->Add(fhInvMassDPhi[i][j]);
              }
            }
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

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  /* the data memebers for this task */
  /* the centrality / multiplicity limits for collecting data in this task instance */
  int ncmranges = 0;
  float* fCentMultMin = nullptr;
  float* fCentMultMax = nullptr;

  /* the data collecting engine instances */
  DataCollectingEngine<false>** dataCE;
  DataCollectingEngine<true>** dataCEsmall;
  DataCollectingEngine<false>** dataCEME;

  /* the input file structure from CCDB */
  TList* ccdblst = nullptr;
  bool loadfromccdb = false;
  std::string cfgCCDBUrl{"http://ccdb-test.cern.ch:8080"};
  std::string cfgCCDBPathNameCorrections{""};
  std::string cfgCCDBDateCorrections{"20220307"};
  std::string cfgCCDBSuffix{""};

  /* pair conversion suppression defaults */
  static constexpr float kCfgPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};
  Configurable<LabeledArray<float>> cfgPairCut{"cfgPairCut", {kCfgPairCutDefaults[0], 5, {"Photon", "K0", "Lambda", "Phi", "Rho"}}, "Conversion suppressions"};
  /* two tracks cut */
  Configurable<float> cfgTwoTrackCut{"cfgTwoTrackCut", -1, "Two-tracks cut: -1 = off; >0 otherwise distance value (suggested: 0.02"};
  Configurable<float> cfgTwoTrackCutMinRadius{"cfgTwoTrackCutMinRadius", 0.8f, "Two-tracks cut: radius in m from which two-tracks cut is applied"};

  Configurable<bool> cfgSmallDCE{"cfgSmallDCE", true, "Use small data collecting engine for singles processing, true = yes. Default = true"};
  Configurable<bool> cfgDoInvMass{"cfgDoInvMass", false, "Do the invariant mass analyis, true = yes. Default = false"};
  Configurable<bool> cfgDoCorrelations{"cfgDoCorrelations", true, "Do the correlations analysis, true = yes. Default = true"};
  Configurable<bool> cfgProcessPairs{"cfgProcessPairs", false, "Process pairs: false = no, just singles, true = yes, process pairs"};
  Configurable<bool> cfgProcessME{"cfgProcessME", false, "Process mixed events: false = no, just same event, true = yes, also process mixed events"};
  Configurable<bool> cfgPtOrder{"cfgPtOrder", false, "enforce pT_1 < pT_2. Defalut: false"};
  Configurable<int> cfgNoOfDimensions{"cfgNoOfDimensions", 1, "Number of dimensions for the NUA&NUE corrections. Default 1"};
  OutputObj<TList> fOutput{"DptDptCorrelationsData", OutputObjHandlingPolicy::AnalysisObject, OutputObjSourceType::OutputObjSource};

  void init(InitContext& initContext)
  {
    using namespace correlationstask;
    using namespace o2::analysis::dptdptfilter;

    /* create the output directory which will own the task output */
    TList* fGlobalOutputList = new TList();
    fGlobalOutputList->SetOwner(true);
    fOutput.setObject(fGlobalOutputList);

    /* check consistency and if there is something to do */
    if (doprocessCleaner) {
      if (doprocessGenLevel || doprocessGenLevelNotStored || doprocessGenLevelMixed || doprocessGenLevelMixedNotStored ||
          doprocessRecLevel || doprocessRecLevelNotStored || doprocessRecLevelMixed || doprocessRecLevelMixedNotStored) {
        LOGF(fatal, "Cleaner process is activated with other processes. Please, fix it!");
      } else {
        /* do nothing. This task will not run! */
        return;
      }
    }

    /* self configure the binning */
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mZVtxbins", zvtxbins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mZVtxmin", zvtxlow, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mZVtxmax", zvtxup, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mPTbins", ptbins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mPTmin", ptlow, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mPTmax", ptup, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mEtabins", etabins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mEtamin", etalow, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mEtamax", etaup, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mPhibins", phibins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mPhibinshift", phibinshift, false);
    philow = 0.0f;
    phiup = constants::math::TwoPI;
    processpairs = cfgProcessPairs.value;
    processmixedevents = cfgProcessME.value;
    ptorder = cfgPtOrder.value;
    invmass = cfgDoInvMass.value;
    corrana = cfgDoCorrelations.value;
    nNoOfDimensions = static_cast<HistoDimensions>(cfgNoOfDimensions.value);

    /* self configure the CCDB access to the input file */
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgCCDB.url", cfgCCDBUrl, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgCCDB.pathNameCorrections", cfgCCDBPathNameCorrections, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgCCDB.dateCorrections", cfgCCDBDateCorrections, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgCCDB.suffix", cfgCCDBSuffix, false);
    loadfromccdb = (cfgCCDBDateCorrections.length() > 0) && (cfgCCDBPathNameCorrections.length() > 0);

    /* update the potential binning change */
    etabinwidth = (etaup - etalow) / static_cast<float>(etabins);
    phibinwidth = (phiup - philow) / static_cast<float>(phibins);

    /* the differential binning */
    deltaetabins = etabins * 2 - 1;
    deltaetalow = etalow - etaup, deltaetaup = etaup - etalow;
    deltaetabinwidth = (deltaetaup - deltaetalow) / static_cast<float>(deltaetabins);
    deltaphibins = phibins;
    deltaphibinwidth = constants::math::TwoPI / deltaphibins;
    deltaphilow = 0.0 - deltaphibinwidth / 2.0;
    deltaphiup = constants::math::TwoPI - deltaphibinwidth / 2.0;

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
      /* self configure the desired species */
      o2::analysis::dptdptfilter::PIDSpeciesSelection pidselector;
      std::vector<std::string> cfgnames = {"cfgElectronPIDSelection", "cfgMuonPIDSelection", "cfgPionPIDSelection", "cfgKaonPIDSelection", "cfgProtonPIDSelection"};
      std::vector<uint8_t> spids = {0, 1, 2, 3, 4};
      for (uint i = 0; i < cfgnames.size(); ++i) {
        auto includeIt = [&pidselector, &initContext](int spid, auto name) {
          bool mUseIt = false;
          bool mExcludeIt = false;
          if (getTaskOptionValue(initContext, "dpt-dpt-filter-tracks", TString::Format("%s.mUseIt", name.c_str()).Data(), mUseIt, false) &&
              getTaskOptionValue(initContext, "dpt-dpt-filter-tracks", TString::Format("%s.mExclude", name.c_str()).Data(), mExcludeIt, false)) {
            if (mUseIt && !mExcludeIt) {
              auto cfg = new o2::analysis::TrackSelectionPIDCfg();
              cfg->mUseIt = true;
              cfg->mExclude = false;
              pidselector.addSpecies(spid, cfg);
            }
          }
        };
        includeIt(spids[i], cfgnames[i]);
      }
      uint nspecies = pidselector.getNSpecies();
      if (nspecies == 0) {
        /* unidentified analysis */
        poinames.push_back(pidselector.getHadFName());
        tnames.push_back(std::string(TString::Format("%sP", pidselector.getHadFName()).Data()));
        tnames.push_back(std::string(TString::Format("%sM", pidselector.getHadFName()).Data()));
        LOGF(info, "Incorporated species name %s to the analysis", poinames[0].c_str());
      } else {
        for (uint8_t ix = 0; ix < nspecies; ++ix) {
          poinames.push_back(std::string(pidselector.getSpeciesFName(ix)));
          tnames.push_back(std::string(TString::Format("%sP", pidselector.getSpeciesFName(ix)).Data()));
          tnames.push_back(std::string(TString::Format("%sM", pidselector.getSpeciesFName(ix)).Data()));
          poimass.push_back(pidselector.getSpeciesMass(ix));
          LOGF(info, "Incorporated species name %s with mass %f to the analysis", poinames[ix].c_str(), poimass[ix]);
        }
      }
      uint ntracknames = tnames.size();
      for (uint isp = 0; isp < ntracknames; ++isp) {
        trackPairsNames.push_back(std::vector<std::string>());
        for (uint jsp = 0; jsp < ntracknames; ++jsp) {
          trackPairsNames[isp].push_back(tnames[isp] + tnames[jsp]);
          LOGF(info, "Incorporated the pair name %s", (tnames[isp] + tnames[jsp]).c_str());
        }
      }

      /* self configure the centrality/multiplicity ranges */
      std::string centspec;
      if (getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgCentSpec", centspec, false)) {
        LOGF(info, "Got the centralities specification: %s", centspec.c_str());
        auto tokens = TString(centspec.c_str()).Tokenize(",");
        ncmranges = tokens->GetEntries();
        fCentMultMin = new float[ncmranges];
        fCentMultMax = new float[ncmranges];
        for (int i = 0; i < ncmranges; ++i) {
          float cmmin = 0.0f;
          float cmmax = 0.0f;
          sscanf(tokens->At(i)->GetName(), "%f-%f", &cmmin, &cmmax);
          fCentMultMin[i] = cmmin;
          fCentMultMax[i] = cmmax;
        }
        delete tokens;
      } else {
        LOGF(info, "No centralities specification. Setting it to: 0-100");
        ncmranges = 1;
        fCentMultMin = new float[ncmranges];
        fCentMultMax = new float[ncmranges];
        fCentMultMin[0] = 0.0f;
        fCentMultMax[0] = 100.0f;
      }
      if (cfgSmallDCE) {
        dataCEsmall = new DataCollectingEngine<true>*[ncmranges];
      } else {
        dataCE = new DataCollectingEngine<false>*[ncmranges];
      }
      if (processmixedevents) {
        dataCEME = new DataCollectingEngine<false>*[ncmranges];
      }

      for (int i = 0; i < ncmranges; ++i) {
        auto initializeCEInstance = [&fGlobalOutputList](auto dce, auto name, bool im, bool corr) {
          /* crete the output list for the passed centrality/multiplicity range */
          TList* fOutputList = new TList();
          fOutputList->SetName(name);
          fOutputList->SetOwner(true);
          /* init the data collection instance */
          if (im) {
            if (corr) {
              dce->template init<true, true>(fOutputList);
            } else {
              dce->template init<true, false>(fOutputList);
            }
          } else {
            if (corr) {
              dce->template init<false, true>(fOutputList);
            } else {
              dce->template init<false, false>(fOutputList);
            }
          }
          fGlobalOutputList->Add(fOutputList);
        };
        auto builSmallDCEInstance = [&initializeCEInstance](auto rg, bool me = false) {
          /* only for singles analysis, no sense of inv mass nor no correlations */
          DataCollectingEngine<true>* dce = new DataCollectingEngine<true>();
          initializeCEInstance(dce, TString::Format("DptDptCorrelationsData%s-%s", me ? "ME" : "", rg), false, false);
          return dce;
        };
        auto buildCEInstance = [&initializeCEInstance](auto rg, bool im, bool corr, bool me = false) {
          DataCollectingEngine<false>* dce = new DataCollectingEngine<false>();
          initializeCEInstance(dce, TString::Format("DptDptCorrelationsData%s-%s", me ? "ME" : "", rg), im, corr);
          return dce;
        };
        TString range = TString::Format("%d-%d", static_cast<int>(fCentMultMin[i]), static_cast<int>(fCentMultMax[i]));
        if (cfgSmallDCE.value) {
          if (processpairs) {
            LOGF(fatal, "Processing pairs cannot be used with the small DCE, please configure properly!!");
          }
          if (invmass) {
            LOGF(fatal, "Invariant mass cannot be used with singles in the small DCE mode, please configure properly!!");
          }
          dataCEsmall[i] = builSmallDCEInstance(range.Data());
        } else {
          if (invmass) {
            if (!processpairs) {
              LOGF(fatal, "Invariant mass cannot be used in processing singles, please configure properly!!");
            }
          }
          dataCE[i] = buildCEInstance(range.Data(), invmass, corrana);
        }
        if (processmixedevents) {
          /* consistency check */
          if (cfgSmallDCE.value) {
            LOGF(fatal, "Mixed events cannot be used with the small DCE, please configure properly!!");
          }
          if (invmass) {
            LOGF(warning, "Invariant mass will not be  used with Mixed events!!");
          }
          if (!corrana) {
            LOGF(fatal, "Mixed events makes not sense to run it without correlations, please configure properly!!");
          }
          dataCEME[i] = buildCEInstance(range.Data(), false, true, true);
        }
      }
      for (int i = 0; i < ncmranges; ++i) {
        LOGF(info, " centrality/multipliicty range: %d, low limit: %0.2f, up limit: %0.2f", i, fCentMultMin[i], fCentMultMax[i]);
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
    ccdb->setURL(cfgCCDBUrl);
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

  const char* getDimensionStr()
  {
    using namespace correlationstask;

    static constexpr std::string_view kStrDim[] = {"", "", "2D", "3D", "4D"};
    return kStrDim[nNoOfDimensions].data();
  }

  template <bool gen, typename FilterdCollision, typename FilteredTracks>
  void processSame(FilterdCollision const& collision, FilteredTracks const& tracks, uint64_t timestamp = 0)
  {
    using namespace correlationstask;

    if (ccdblst == nullptr) {
      if (loadfromccdb) {
        ccdblst = getCCDBInput(ccdb, cfgCCDBPathNameCorrections.c_str(), cfgCCDBDateCorrections.c_str(), true, cfgCCDBSuffix);
      }
    }

    /* locate the data collecting engine for the collision centrality/multiplicity */
    int ixDCE = getDCEindex(collision);
    if (!(ixDCE < 0)) {
      auto isCCDBstored = [&]() {
        if (cfgSmallDCE.value) {
          return dataCEsmall[ixDCE]->isCCDBstored();
        } else {
          return dataCE[ixDCE]->isCCDBstored();
        }
      };
      auto storePtAverages = [&](auto& ptavgs) {
        if (cfgSmallDCE.value) {
          dataCEsmall[ixDCE]->storePtAverages(ptavgs);
        } else {
          dataCE[ixDCE]->storePtAverages(ptavgs);
        }
      };
      auto storeTrackCorrections = [&](auto& corrs) {
        if (cfgSmallDCE.value) {
          dataCEsmall[ixDCE]->storeTrackCorrections(corrs);
        } else {
          dataCE[ixDCE]->storeTrackCorrections(corrs);
        }
      };
      if (ccdblst != nullptr && !(isCCDBstored())) {
        if constexpr (gen) {
          std::vector<TH2*> ptavgs{tnames.size(), nullptr};
          for (uint isp = 0; isp < tnames.size(); ++isp) {
            ptavgs[isp] = reinterpret_cast<TH2*>(ccdblst->FindObject(
              TString::Format("trueptavgetaphi_%02d-%02d_%s",
                              static_cast<int>(fCentMultMin[ixDCE]),
                              static_cast<int>(fCentMultMax[ixDCE]),
                              tnames[isp].c_str())
                .Data()));
          }
          storePtAverages(ptavgs);
        } else {
          std::vector<TH1*> corrs{tnames.size(), nullptr};
          for (uint isp = 0; isp < tnames.size(); ++isp) {
            auto hName = TString::Format("correction%s_%02d-%02d_%s", getDimensionStr(), static_cast<int>(fCentMultMin[ixDCE]), static_cast<int>(fCentMultMax[ixDCE]), tnames[isp].c_str());
            corrs[isp] = reinterpret_cast<TH1*>(ccdblst->FindObject(hName.Data()));
            if (corrs[isp] != nullptr) {
              LOGF(info, "Loaded %s", corrs[isp]->GetName());
            } else {
              LOGF(warning, "No correction histogram for species %d with name %s", isp, hName.Data());
            }
          }
          storeTrackCorrections(corrs);
          std::vector<TH2*> ptavgs{tnames.size(), nullptr};
          for (uint isp = 0; isp < tnames.size(); ++isp) {
            ptavgs[isp] = reinterpret_cast<TH2*>(ccdblst->FindObject(
              TString::Format("ptavgetaphi_%02d-%02d_%s",
                              static_cast<int>(fCentMultMin[ixDCE]),
                              static_cast<int>(fCentMultMax[ixDCE]),
                              tnames[isp].c_str())
                .Data()));
          }
          storePtAverages(ptavgs);
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
        dataCEsmall[ixDCE]->processCollision<false>(tracks, tracks, collision.posZ(), collision.centmult(), bfield);
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
        ccdblst = getCCDBInput(ccdb, cfgCCDBPathNameCorrections.c_str(), cfgCCDBDateCorrections.c_str(), true, cfgCCDBSuffix);
      }
    }

    /* locate the data collecting engine for the collision centrality/multiplicity */
    int ixDCE = getDCEindex(collision);
    if (!(ixDCE < 0)) {
      if (ccdblst != nullptr && !(dataCEME[ixDCE]->isCCDBstored())) {
        if constexpr (gen) {
          std::vector<TH2*> ptavgs{tnames.size(), nullptr};
          for (uint isp = 0; isp < tnames.size(); ++isp) {
            ptavgs[isp] = reinterpret_cast<TH2*>(ccdblst->FindObject(
              TString::Format("trueptavgetaphi_%02d-%02d_%s",
                              static_cast<int>(fCentMultMin[ixDCE]),
                              static_cast<int>(fCentMultMax[ixDCE]),
                              tnames[isp].c_str())
                .Data()));
          }
          dataCEME[ixDCE]->storePtAverages(ptavgs);
        } else {
          std::vector<TH1*> corrs{tnames.size(), nullptr};
          for (uint isp = 0; isp < tnames.size(); ++isp) {
            auto hName = TString::Format("correction%s_%02d-%02d_%s", getDimensionStr(), static_cast<int>(fCentMultMin[ixDCE]), static_cast<int>(fCentMultMax[ixDCE]), tnames[isp].c_str());
            corrs[isp] = reinterpret_cast<TH1*>(ccdblst->FindObject(hName.Data()));
            if (corrs[isp] != nullptr) {
              LOGF(info, "Loaded %s", corrs[isp]->GetName());
            } else {
              LOGF(warning, "No correction histogram for species %d with name %s", isp, hName.Data());
            }
          }
          dataCEME[ixDCE]->storeTrackCorrections(corrs);
          std::vector<TH2*> ptavgs{tnames.size(), nullptr};
          for (uint isp = 0; isp < tnames.size(); ++isp) {
            ptavgs[isp] = reinterpret_cast<TH2*>(ccdblst->FindObject(
              TString::Format("ptavgetaphi_%02d-%02d_%s",
                              static_cast<int>(fCentMultMin[ixDCE]),
                              static_cast<int>(fCentMultMax[ixDCE]),
                              tnames[isp].c_str())
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

  void processRecLevel(soa::Filtered<aod::DptDptCFAcceptedCollisions>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<aod::ScannedTracks> const& tracks)
  {
    processSame<false>(collision, tracks, collision.bc_as<aod::BCsWithTimestamps>().timestamp());
  }
  PROCESS_SWITCH(DptDptCorrelations, processRecLevel, "Process reco level correlations", false);

  void processRecLevelCheck(aod::Collisions const& collisions, aod::Tracks const& tracks)
  {
    int nAssignedTracks = 0;
    int nNotAssignedTracks = 0;
    int64_t firstNotAssignedIndex = -1;
    int64_t lastNotAssignedIndex = -1;

    for (auto const& track : tracks) {
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
  PROCESS_SWITCH(DptDptCorrelations, processRecLevelCheck, "Process reco level checks", true);

  void processGenLevelCheck(aod::McCollisions const& mccollisions, aod::McParticles const& particles)
  {
    int nAssignedParticles = 0;
    int nNotAssignedParticles = 0;
    int64_t firstNotAssignedIndex = -1;
    int64_t lastNotAssignedIndex = -1;

    for (auto const& particle : particles) {
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
  PROCESS_SWITCH(DptDptCorrelations, processGenLevelCheck, "Process generator level checks", true);

  void processRecLevelNotStored(
    soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
    aod::BCsWithTimestamps const&,
    soa::Filtered<soa::Join<aod::Tracks, aod::DptDptCFTracksInfo>> const& tracks)
  {
    processSame<false>(collision, tracks, collision.bc_as<aod::BCsWithTimestamps>().timestamp());
  }
  PROCESS_SWITCH(DptDptCorrelations, processRecLevelNotStored, "Process reco level correlations for not stored derived data", true);

  void processGenLevel(
    soa::Filtered<aod::DptDptCFAcceptedTrueCollisions>::iterator const& collision,
    soa::Filtered<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>> const& tracks)
  {
    processSame<true>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptCorrelations, processGenLevel, "Process generator level correlations", false);

  void processGenLevelNotStored(
    soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>::iterator const& collision,
    soa::Filtered<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>> const& particles)
  {
    processSame<true>(collision, particles);
  }
  PROCESS_SWITCH(DptDptCorrelations, processGenLevelNotStored, "Process generator level correlations for not stored derived data", false);

  std::vector<double>
    vtxBinsEdges{VARIABLE_WIDTH, -7.0f, -5.0f, -3.0f, -1.0f, 1.0f, 3.0f, 5.0f, 7.0f};
  std::vector<double> multBinsEdges{VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0, 100.1f};
  SliceCache cache;
  using BinningZVtxMultRec = ColumnBinningPolicy<aod::collision::PosZ, aod::dptdptfilter::DptDptCFCollisionCentMult>;
  BinningZVtxMultRec bindingOnVtxAndMultRec{{vtxBinsEdges, multBinsEdges}, true}; // true is for 'ignore overflows' (true by default)
  static constexpr int kNoOfLoggingCombinations = 10;

  void processRecLevelMixed(soa::Filtered<aod::DptDptCFAcceptedCollisions> const& collisions, aod::BCsWithTimestamps const&, soa::Filtered<aod::ScannedTracks> const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<soa::Filtered<aod::DptDptCFAcceptedCollisions>, soa::Filtered<aod::ScannedTracks>, BinningZVtxMultRec> pairreco{bindingOnVtxAndMultRec, 5, -1, collisions, tracksTuple, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

    LOGF(DPTDPTLOGCOLLISIONS, "Received %d collisions", collisions.size());
    int logcomb = 0;
    for (auto const& [collision1, tracks1, collision2, tracks2] : pairreco) {
      if (logcomb < kNoOfLoggingCombinations) {
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
  PROCESS_SWITCH(DptDptCorrelations, processRecLevelMixed, "Process reco level mixed events correlations", false);

  void processRecLevelMixedNotStored(
    soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>> const& collisions,
    aod::BCsWithTimestamps const&,
    soa::Filtered<soa::Join<aod::Tracks, aod::DptDptCFTracksInfo>> const& tracks)
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
    for (auto const& [collision1, tracks1, collision2, tracks2] : pairreco) {
      if (logcomb < kNoOfLoggingCombinations) {
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
  PROCESS_SWITCH(DptDptCorrelations, processRecLevelMixedNotStored, "Process reco level mixed events correlations for not stored derived data", false);

  using BinningZVtxMultGen = ColumnBinningPolicy<aod::mccollision::PosZ, aod::dptdptfilter::DptDptCFCollisionCentMult>;
  BinningZVtxMultGen bindingOnVtxAndMultGen{{vtxBinsEdges, multBinsEdges}, true}; // true is for 'ignore overflows' (true by default)

  void processGenLevelMixed(soa::Filtered<aod::DptDptCFAcceptedTrueCollisions> const& collisions, soa::Filtered<aod::ScannedTrueTracks> const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<soa::Filtered<aod::DptDptCFAcceptedTrueCollisions>, soa::Filtered<aod::ScannedTrueTracks>, BinningZVtxMultGen> pairgen{bindingOnVtxAndMultGen, 5, -1, collisions, tracksTuple, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

    LOGF(DPTDPTLOGCOLLISIONS, "Received %d generated collisions", collisions.size());
    int logcomb = 0;
    for (auto const& [collision1, tracks1, collision2, tracks2] : pairgen) {
      if (logcomb < kNoOfLoggingCombinations) {
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
  PROCESS_SWITCH(DptDptCorrelations, processGenLevelMixed, "Process generator level mixed events correlations", false);

  void processGenLevelMixedNotStored(soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>> const& collisions, soa::Filtered<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>> const& tracks)
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
    for (auto const& [collision1, tracks1, collision2, tracks2] : pairgen) {
      if (logcomb < kNoOfLoggingCombinations) {
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
  PROCESS_SWITCH(DptDptCorrelations, processGenLevelMixedNotStored, "Process generator level mixed events correlations for not stored derived data", false);

  /// cleans the output object when the task is not used
  void processCleaner(soa::Filtered<aod::DptDptCFAcceptedCollisions> const& colls)
  {
    LOGF(DPTDPTLOGCOLLISIONS, "Got %d new collisions", colls.size());
    fOutput->Clear();
  }
  PROCESS_SWITCH(DptDptCorrelations, processCleaner, "Cleaner process for not used output", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  o2::analysis::dptdptfilter::metadataInfo.initMetadata(cfgc);
  WorkflowSpec workflow{
    adaptAnalysisTask<DptDptCorrelations>(cfgc, TaskName{"DptDptCorrelationsRec"}, SetDefaultProcesses{{{"processRecLevel", true}, {"processRecLevelMixed", false}, {"processCleaner", false}}}),  // o2-linter: disable=name/o2-task (It is adapted multiple times)
    adaptAnalysisTask<DptDptCorrelations>(cfgc, TaskName{"DptDptCorrelationsGen"}, SetDefaultProcesses{{{"processGenLevel", false}, {"processGenLevelMixed", false}, {"processCleaner", true}}})}; // o2-linter: disable=name/o2-task (It is adapted multiple times)
  return workflow;
}
