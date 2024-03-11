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

/// \file MeanPtFlucIdentifiedMC.cxx
/// \brief Monte Carlo Study for EbyE <pt> fluctuations with moments method.
///        For charged particles and identified particles.
///
/// \author Tanu Gahlaut <tanu.gahlaut@cern.ch>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"

#include "TDatabasePDG.h"
#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

double massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
double massKa = TDatabasePDG::Instance()->GetParticle(321)->Mass();
double massPr = TDatabasePDG::Instance()->GetParticle(2212)->Mass();

struct meanPtFlucIdMC {

  Configurable<int> nPtBins{"nPtBins", 30, ""};
  Configurable<float> ptMax{"ptMax", 2.0, "maximum pT"};
  Configurable<float> ptMin{"ptMin", 0.15, "minimum pT"};
  Configurable<float> etaCut{"etaCut", 0.8, "Eta cut"};
  Configurable<float> rapCut{"rapCut", 0.5, "Rapidity Cut"};
  Configurable<float> dcaXYCut{"dcaXYCut", 0.12, "DCAxy cut"};
  Configurable<float> dcaZCut{"dcaZCut", 1.0, "DCAz cut"};
  Configurable<float> posZCut{"posZCut", 10.0, "cut for vertex Z"};
  Configurable<float> nSigCut1{"nSigCut1", 1.0, "nSigma cut (1)"};
  Configurable<float> nSigCut2{"nSigCut2", 2.0, "nSigma cut (2)"};
  Configurable<float> nSigCut3{"nSigCut3", 3.0, "nSigma cut (3)"};
  Configurable<float> nSigCut4{"nSigCut4", 4.0, "nSigma cut (4)"};
  Configurable<float> nSigCut5{"nSigCut5", 5.0, "nSigma cut (5)"};
  Configurable<float> nSigCut15{"nSigCut15", 1.5, "nSigma cut (1.5)"};
  Configurable<float> nSigCut25{"nSigCut25", 2.5, "nSigma cut (2.5)"};
  Configurable<float> piP1{"piP1", 0.65, "pion p (1)"};
  Configurable<float> piP2{"piP2", 0.70, "pion p (2)"};
  Configurable<float> piP3{"piP3", 1.40, "pion p (3)"};
  Configurable<float> piP4{"piP4", 1.70, "pion p (4)"};
  Configurable<float> kaP1{"kaP1", 0.20, "min kaon p (1)"};
  Configurable<float> kaP2{"kaP2", 0.5, "kaon p (2)"};
  Configurable<float> kaP3{"kaP3", 0.55, "kaon p (3)"};
  Configurable<float> kaP4{"kaP4", 0.60, "kaon p (4)"};
  Configurable<float> kaP5{"kaP5", 0.65, "kaon p (5)"};
  Configurable<float> kaP6{"kaP6", 1.10, "kaon p (6)"};
  Configurable<float> kaP7{"kaP7", 1.28, "kaon p (7)"};
  Configurable<float> kaP8{"kaP8", 1.50, "kaon p (8)"};
  Configurable<float> prP1{"prP1", 0.40, "min proton p (1)"};
  Configurable<float> prP2{"prP2", 0.95, "proton p (2)"};
  Configurable<float> prP3{"prP3", 1.00, "proton p (3)"};
  Configurable<float> prP4{"prP4", 1.05, "proton p (4)"};
  Configurable<float> prP5{"prP5", 1.13, "proton p (5)"};
  Configurable<float> prP6{"prP6", 1.18, "proton p (6)"};

  using MyRun3MCCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using MyRun2MCCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using TCs = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
                        aod::pidTOFFullPi, aod::pidTPCFullPi, aod::pidTOFFullPr, aod::pidTPCFullPr,
                        aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullEl, aod::pidTPCFullEl,
                        aod::pidTOFbeta, aod::McTrackLabels>;

  HistogramRegistry hist{"hist", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    const AxisSpec axisEvents{5, 0, 5, "Counts"};
    const AxisSpec axisPt{nPtBins, 0., 3., "p_{T} (GeV/c)"};
    const AxisSpec axisVtxZ{80, -20., 20., "V_{Z} (cm)"};
    const AxisSpec axisMult{100, 0, 100, "N_{ch}"};
    const AxisSpec axisMeanPt{100, 0., 3., "M(p_{T}) (GeV/c)"};
    const AxisSpec axisP{nPtBins, 0., 3., "p (GeV/c)"};
    const AxisSpec axisTPCSignal{180, 20., 200., "#frac{dE}{dx}"};
    const AxisSpec axisTOFSignal{100, 0.2, 1.2, "TOF #beta"};

    hist.add("Reco/Counts", "Counts", kTH1D, {axisEvents});
    hist.add("Reco/vtxZ", "Vertex Z ", kTH1D, {axisVtxZ});
    hist.add("Reco/Charged/h_Mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Reco/Charged/h_Pt", "p_{T} Charged Particles", kTH1D, {axisPt});
    hist.add("Reco/Charged/h_mean_pt", " <p_{T}> ", kTH1D, {axisMeanPt});
    hist.add("Reco/Charged/h_mean_Q1_Mult", " <p_{T}> vs N_{ch} ", kTProfile, {axisMult});
    hist.add("Reco/Charged/h_mean_Q1_Mult_var", " <p_{T}> vs N_{ch} ", kTProfile, {axisMult});
    hist.add("Reco/Charged/h_twopart_Mult_var", "Twopart vs N_{ch} ", kTProfile, {axisMult});
    hist.add("QA/Pion/h2_tofSignal", "TOF Signal ", kTH2D, {{axisP}, {axisTOFSignal}});
    hist.add("QA/Pion/h2_tpcSignal", "TPC Signal #frac{dE}{dx}", kTH2D, {{axisP}, {axisTPCSignal}});

    hist.addClone("Reco/Charged/", "Reco/Pion/");
    hist.addClone("Reco/Charged/", "Reco/Kaon/");
    hist.addClone("Reco/Charged/", "Reco/Proton/");

    hist.addClone("QA/Pion/", "QA/Kaon/");
    hist.addClone("QA/Pion/", "QA/Proton/");

    hist.addClone("Reco/", "Gen/");
  }

  template <typename T>
  bool selRecoColRun2(T const& col)
  {
    if (std::abs(col.posZ()) > posZCut)
      return false;

    if (!col.sel7())
      return false;

    return true;
  }

  template <typename T>
  bool selRecoColRun3(T const& col)
  {
    if (std::abs(col.posZ()) > posZCut)
      return false;

    if (!col.sel8())
      return false;

    return true;
  }

  template <typename T>
  bool selRecoTracks(T const& track)
  {
    if (track.pt() < ptMin)
      return false;

    if (track.pt() > ptMax)
      return false;

    if (track.sign() == 0)
      return false;

    if (std::abs(track.eta()) > etaCut)
      return false;

    if (std::abs(track.dcaZ()) > dcaZCut)
      return false;

    if (std::abs(track.dcaXY()) > dcaXYCut)
      return false;

    if (!track.isGlobalTrack())
      return false;

    return true;
  }

  template <typename T>
  bool selPions(T const& track)
  {
    if (std::abs(track.mcParticle().pdgCode()) == 211 && (((!track.hasTOF()) &&
                                                           ((std::abs(track.tpcNSigmaPi()) < nSigCut3 && track.p() <= piP1) || (std::abs(track.tpcNSigmaPi()) < nSigCut2 && track.p() > piP1 && track.p() <= piP2))) ||
                                                          (track.hasTOF() && std::abs(track.tpcNSigmaPi()) < nSigCut4 && std::abs(track.tofNSigmaEl()) > nSigCut1 &&
                                                           ((std::abs(track.tofNSigmaPi()) < nSigCut3 && track.p() <= piP3) || (std::abs(track.tofNSigmaPi()) < nSigCut25 && track.p() > piP3 && track.p() <= piP4) || (std::abs(track.tofNSigmaPi()) < nSigCut2 && track.p() > piP4))))) {
      if (abs(track.rapidity(massPi)) < 0.5)
        return true;
    }

    return false;
  }

  template <typename T>
  bool selKaons(T const& track)
  {
    if (std::abs(track.mcParticle().pdgCode()) == 321 && (((!track.hasTOF()) &&
                                                           ((std::abs(track.tpcNSigmaKa()) < nSigCut3 && track.pt() > kaP1 && track.p() <= kaP2) || (std::abs(track.tpcNSigmaKa()) < nSigCut25 && track.p() > kaP2 && track.p() <= kaP3) || (std::abs(track.tpcNSigmaKa()) < nSigCut2 && track.p() > kaP3 && track.p() <= kaP4) || (std::abs(track.tpcNSigmaKa()) < nSigCut15 && track.p() > kaP4 && track.p() <= kaP5))) ||
                                                          (track.hasTOF() && std::abs(track.tpcNSigmaKa()) < nSigCut4 && std::abs(track.tofNSigmaEl()) > nSigCut1 &&
                                                           ((std::abs(track.tofNSigmaKa()) < nSigCut3 && track.pt() > kaP1 && track.p() <= kaP6) || (std::abs(track.tofNSigmaKa()) < nSigCut2 && track.p() > kaP6 && track.p() <= kaP7) || (std::abs(track.tofNSigmaKa()) < nSigCut15 && track.p() > kaP7 && track.p() <= kaP8) || (std::abs(track.tofNSigmaKa()) < nSigCut1 && track.p() > kaP8))))) {
      if (abs(track.rapidity(massKa)) < 0.5)
        return true;
    }

    return false;
  }

  template <typename T>
  bool selProtons(T const& track)
  {
    if (std::abs(track.mcParticle().pdgCode()) == 2212 && (((!track.hasTOF()) &&
                                                            ((std::abs(track.tpcNSigmaPr()) < nSigCut3 && track.pt() > prP1 && track.p() <= prP2) || (std::abs(track.tpcNSigmaPr()) < nSigCut25 && track.p() > prP2 && track.p() <= prP3) ||
                                                             (std::abs(track.tpcNSigmaPr()) < nSigCut2 && track.p() > prP3 && track.p() <= prP4) || (std::abs(track.tpcNSigmaPr()) < nSigCut15 && track.p() > prP4 && track.p() <= prP5) || (std::abs(track.tpcNSigmaPr()) < nSigCut1 && track.p() > prP5 && track.p() <= prP6))) ||
                                                           (track.hasTOF() && std::abs(track.tpcNSigmaPr()) < nSigCut4 && std::abs(track.tofNSigmaEl()) > nSigCut1 && std::abs(track.tofNSigmaPr()) < nSigCut3 && track.pt() > prP1))) {
      if (abs(track.rapidity(massPr)) < 0.5)
        return true;
    }

    return false;
  }

  void moments(double pt, double* Q1, double* Q2)
  {
    *Q1 += pt;
    *Q2 += pt * pt;
  }

  void parts(double Q1, double Q2, int N, double* mean_Q1, double* twopart)
  {
    if (N > 1) {
      *mean_Q1 = Q1 / static_cast<double>(N);
      *twopart = ((Q1 * Q1) - Q2) / (static_cast<double>(N) * (static_cast<double>(N) - 1));
    }
  }

  template <typename T, typename U>
  void FillRecoHistos(T const& col, U const& tracks)
  {
    int N_Pi = 0, N_Ka = 0, N_Pr = 0;
    int Nch = 0;
    double pt_ch = 0, Q1_ch = 0, Q2_ch = 0;
    double pt_Pi = 0, Q1_Pi = 0, Q2_Pi = 0;
    double pt_Pr = 0, Q1_Pr = 0, Q2_Pr = 0;
    double pt_Ka = 0, Q1_Ka = 0, Q2_Ka = 0;
    double mean_Q1_Ch, mean_Q1_Pi, mean_Q1_Ka, mean_Q1_Pr;
    double twopart_Ch, twopart_Pi, twopart_Ka, twopart_Pr;

    for (auto& track : tracks) {
      if (selRecoTracks(track) && track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();
        if (mcParticle.isPhysicalPrimary()) {
          Nch++;
          pt_ch = track.pt();
          moments(pt_ch, &Q1_ch, &Q2_ch);
          hist.fill(HIST("Reco/Charged/h_Pt"), track.pt());

          if (selPions(track)) {
            N_Pi++;
            pt_Pi = track.pt();
            moments(pt_Pi, &Q1_Pi, &Q2_Pi);
            hist.fill(HIST("Reco/Pion/h_Pt"), track.pt());
            hist.fill(HIST("QA/Pion/h2_tpcSignal"), track.p(), track.tpcSignal());
            hist.fill(HIST("QA/Pion/h2_tofSignal"), track.p(), track.beta());
          }

          if (selKaons(track)) {
            N_Ka++;
            pt_Ka = track.pt();
            moments(pt_Ka, &Q1_Ka, &Q2_Ka);
            hist.fill(HIST("Reco/Kaon/h_Pt"), track.pt());
            hist.fill(HIST("QA/Kaon/h2_tpcSignal"), track.p(), track.tpcSignal());
            hist.fill(HIST("QA/Kaon/h2_tofSignal"), track.p(), track.beta());
          }

          if (selProtons(track)) {
            N_Pr++;
            pt_Pr = track.pt();
            moments(pt_Pr, &Q1_Pr, &Q2_Pr);
            hist.fill(HIST("Reco/Proton/h_Pt"), track.pt());
            hist.fill(HIST("QA/Proton/h2_tpcSignal"), track.p(), track.tpcSignal());
            hist.fill(HIST("QA/Proton/h2_tofSignal"), track.p(), track.beta());
          }
        }
      }
    }
    hist.fill(HIST("Reco/Counts"), 2);
    hist.fill(HIST("Reco/vtxZ"), col.posZ());

    static constexpr std::string_view dire[] = {"Reco/Charged/", "Reco/Pion/", "Reco/Kaon/", "Reco/Proton/"};

    hist.fill(HIST(dire[0]) + HIST("h_Mult"), Nch);
    hist.fill(HIST(dire[1]) + HIST("h_Mult"), N_Pi);
    hist.fill(HIST(dire[2]) + HIST("h_Mult"), N_Ka);
    hist.fill(HIST(dire[3]) + HIST("h_Mult"), N_Pr);

    parts(Q1_ch, Q2_ch, Nch, &mean_Q1_Ch, &twopart_Ch);
    if (Nch > 0 && mean_Q1_Ch != 0) {
      hist.fill(HIST(dire[0]) + HIST("h_mean_pt"), mean_Q1_Ch);
      hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult"), Nch, mean_Q1_Ch);
    }
    if (Nch > 1) {
      if (mean_Q1_Ch != 0)
        hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult_var"), Nch, mean_Q1_Ch);
      if (twopart_Ch != 0)
        hist.fill(HIST(dire[0]) + HIST("h_twopart_Mult_var"), Nch, twopart_Ch);
    }

    parts(Q1_Pi, Q2_Pi, N_Pi, &mean_Q1_Pi, &twopart_Pi);
    if (N_Pi > 0 && mean_Q1_Pi != 0) {
      hist.fill(HIST(dire[1]) + HIST("h_mean_pt"), mean_Q1_Pi);
      hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult"), Nch, mean_Q1_Pi);
    }
    if (N_Pi > 1) {
      if (mean_Q1_Pi != 0)
        hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult_var"), Nch, mean_Q1_Pi);
      if (twopart_Pi != 0)
        hist.fill(HIST(dire[1]) + HIST("h_twopart_Mult_var"), Nch, twopart_Pi);
    }

    parts(Q1_Ka, Q2_Ka, N_Ka, &mean_Q1_Ka, &twopart_Ka);
    if (N_Ka > 0 && mean_Q1_Ka != 0) {
      hist.fill(HIST(dire[2]) + HIST("h_mean_pt"), mean_Q1_Ka);
      hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult"), Nch, mean_Q1_Ka);
    }
    if (N_Ka > 1) {
      if (mean_Q1_Ka != 0)
        hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult_var"), Nch, mean_Q1_Ka);
      if (twopart_Ka != 0)
        hist.fill(HIST(dire[2]) + HIST("h_twopart_Mult_var"), Nch, twopart_Ka);
    }

    parts(Q1_Pr, Q2_Pr, N_Pr, &mean_Q1_Pr, &twopart_Pr);
    if (N_Pr > 0 && mean_Q1_Pr != 0) {
      hist.fill(HIST(dire[3]) + HIST("h_mean_pt"), mean_Q1_Pr);
      hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult"), Nch, mean_Q1_Pr);
    }
    if (N_Pr > 1) {
      if (mean_Q1_Pr != 0)
        hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult_var"), Nch, mean_Q1_Pr);
      if (twopart_Pr != 0)
        hist.fill(HIST(dire[3]) + HIST("h_twopart_Mult_var"), Nch, twopart_Pr);
    }
  }

  void process_MCRecoRun2(MyRun2MCCollisions::iterator const& col, aod::McCollisions const&, TCs const& tracks, aod::McParticles const&)
  {
    if (selRecoColRun2(col)) {
      FillRecoHistos(col, tracks);
    }
  }
  PROCESS_SWITCH(meanPtFlucIdMC, process_MCRecoRun2, "process MC Reconstructed Run-2", false);

  void process_MCRecoRun3(MyRun3MCCollisions::iterator const& col, aod::McCollisions const&, TCs const& tracks, aod::McParticles const&)
  {
    if (selRecoColRun3(col)) {
      FillRecoHistos(col, tracks);
    }
  }
  PROCESS_SWITCH(meanPtFlucIdMC, process_MCRecoRun3, "process MC Reconstructed Run-3", true);

  void process_MCGen(aod::McCollisions::iterator const& mccol, aod::McParticles const& McParticles)
  {
    int N_Pi = 0, N_Ka = 0, N_Pr = 0;
    int Nch = 0;
    double pt_ch = 0, Q1_ch = 0, Q2_ch = 0;
    double pt_Pi = 0, Q1_Pi = 0, Q2_Pi = 0;
    double pt_Pr = 0, Q1_Pr = 0, Q2_Pr = 0;
    double pt_Ka = 0, Q1_Ka = 0, Q2_Ka = 0;
    double mean_Q1_Ch, mean_Q1_Pi, mean_Q1_Ka, mean_Q1_Pr;
    double twopart_Ch, twopart_Pi, twopart_Ka, twopart_Pr;

    if (abs(mccol.posZ()) > posZCut)
      return;

    for (auto& mcParticle : McParticles) {
      if (mcParticle.isPhysicalPrimary() && mcParticle.pt() > ptMin && mcParticle.pt() < ptMax && std::abs(mcParticle.eta()) < etaCut) {
        Nch++;
        pt_ch = mcParticle.pt();
        moments(pt_ch, &Q1_ch, &Q2_ch);
        hist.fill(HIST("Gen/Charged/h_Pt"), mcParticle.pt());

        if (std::abs(mcParticle.pdgCode()) == 211 && abs(mcParticle.y()) < rapCut) {
          N_Pi++;
          pt_Pi = mcParticle.pt();
          moments(pt_Pi, &Q1_Pi, &Q2_Pi);
          hist.fill(HIST("Gen/Pion/h_Pt"), mcParticle.pt());
        }

        if (std::abs(mcParticle.pdgCode()) == 321 && abs(mcParticle.y()) < rapCut) {
          N_Ka++;
          pt_Ka = mcParticle.pt();
          moments(pt_Ka, &Q1_Ka, &Q2_Ka);
          hist.fill(HIST("Gen/Kaon/h_Pt"), mcParticle.pt());
        }

        if (std::abs(mcParticle.pdgCode()) == 2212 && abs(mcParticle.y()) < rapCut) {
          N_Pr++;
          pt_Pr = mcParticle.pt();
          moments(pt_Pr, &Q1_Pr, &Q2_Pr);
          hist.fill(HIST("Gen/Proton/h_Pt"), mcParticle.pt());
        }
      }
    }
    hist.fill(HIST("Gen/Counts"), 2);
    hist.fill(HIST("Gen/vtxZ"), mccol.posZ());
    static constexpr std::string_view dire[] = {"Gen/Charged/", "Gen/Pion/", "Gen/Kaon/", "Gen/Proton/"};

    hist.fill(HIST(dire[0]) + HIST("h_Mult"), Nch);
    hist.fill(HIST(dire[1]) + HIST("h_Mult"), N_Pi);
    hist.fill(HIST(dire[2]) + HIST("h_Mult"), N_Ka);
    hist.fill(HIST(dire[3]) + HIST("h_Mult"), N_Pr);

    parts(Q1_ch, Q2_ch, Nch, &mean_Q1_Ch, &twopart_Ch);
    if (Nch > 0 && mean_Q1_Ch != 0) {
      hist.fill(HIST(dire[0]) + HIST("h_mean_pt"), mean_Q1_Ch);
      hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult"), Nch, mean_Q1_Ch);
    }
    if (Nch > 1) {
      if (mean_Q1_Ch != 0)
        hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult_var"), Nch, mean_Q1_Ch);
      if (twopart_Ch != 0)
        hist.fill(HIST(dire[0]) + HIST("h_twopart_Mult_var"), Nch, twopart_Ch);
    }

    parts(Q1_Pi, Q2_Pi, N_Pi, &mean_Q1_Pi, &twopart_Pi);
    if (N_Pi > 0 && mean_Q1_Pi != 0) {
      hist.fill(HIST(dire[1]) + HIST("h_mean_pt"), mean_Q1_Pi);
      hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult"), Nch, mean_Q1_Pi);
    }
    if (N_Pi > 1) {
      if (mean_Q1_Pi != 0)
        hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult_var"), Nch, mean_Q1_Pi);
      if (twopart_Pi != 0)
        hist.fill(HIST(dire[1]) + HIST("h_twopart_Mult_var"), Nch, twopart_Pi);
    }

    parts(Q1_Ka, Q2_Ka, N_Ka, &mean_Q1_Ka, &twopart_Ka);
    if (N_Ka > 0 && mean_Q1_Ka != 0) {
      hist.fill(HIST(dire[2]) + HIST("h_mean_pt"), mean_Q1_Ka);
      hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult"), Nch, mean_Q1_Ka);
    }
    if (N_Ka > 1) {
      if (mean_Q1_Ka != 0)
        hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult_var"), Nch, mean_Q1_Ka);
      if (twopart_Ka != 0)
        hist.fill(HIST(dire[2]) + HIST("h_twopart_Mult_var"), Nch, twopart_Ka);
    }

    parts(Q1_Pr, Q2_Pr, N_Pr, &mean_Q1_Pr, &twopart_Pr);
    if (N_Pr > 0 && mean_Q1_Pr != 0) {
      hist.fill(HIST(dire[3]) + HIST("h_mean_pt"), mean_Q1_Pr);
      hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult"), Nch, mean_Q1_Pr);
    }
    if (N_Pr > 1) {
      if (mean_Q1_Pr != 0)
        hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult_var"), Nch, mean_Q1_Pr);
      if (twopart_Pr != 0)
        hist.fill(HIST(dire[3]) + HIST("h_twopart_Mult_var"), Nch, twopart_Pr);
    }
  }
  PROCESS_SWITCH(meanPtFlucIdMC, process_MCGen, "process MC Generated", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<meanPtFlucIdMC>(cfgc)};
}
