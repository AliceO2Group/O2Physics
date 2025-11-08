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

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <memory>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace idr2p2columns
{
DECLARE_SOA_COLUMN(BinNPIDFlag, binNpid, int8_t); // Flag tracks without proper binning as -1, and indicate type of particle 0->un-Id, 1->pion, 2->kaon, 3->proton
} // namespace idr2p2columns
DECLARE_SOA_TABLE(Flags, "AOD", "Flags", idr2p2columns::BinNPIDFlag);
} // namespace o2::aod
struct FillFlagsTable {
  Configurable<std::vector<float>> TPCnsigmacutsPi{"TPCnsigmacutsPi", {2.5}, "TPC Nsigma cuts for Pion"};
  Configurable<std::vector<float>> TPCpTrangesPi{"TPCpTrangesPi", {2}, "TPC pT ranges for Pion"};
  Configurable<std::vector<float>> TOFnsigmacutsPi{"TOFnsigmacutsPi", {2.5}, "TOF Nsigma cuts for Pion"};
  Configurable<std::vector<float>> TOFpTrangesPi{"TOFpTrangesPi", {0.6}, "TOF pT ranges for Pion"};
  Configurable<std::vector<float>> TPCnsigmacutsKa{"TPCnsigmacutsKa", {2, 1, 0.6, 2}, "TPC Nsigma cuts for Kaon"};
  Configurable<std::vector<float>> TPCpTrangesKa{"TPCpTrangesKa", {0.45, 0.55, 0.6, 2}, "TPC pT ranges for Kaon"};
  Configurable<std::vector<float>> TOFnsigmacutsKa{"TOFnsigmacutsKa", {2}, "TOF Nsigma cuts for Kaon"};
  Configurable<std::vector<float>> TOFpTrangesKa{"TOFpTrangesKa", {0.6}, "TOF pT ranges for Kaon"};
  Configurable<std::vector<float>> TPCnsigmacutsPr{"TPCnsigmacutsPr", {2.2, 1, 2.2}, "TPC Nsigma cuts for Proton"};
  Configurable<std::vector<float>> TPCpTrangesPr{"TPCpTrangesPr", {0.85, 1.1, 2}, "TPC pT ranges for Proton"};
  Configurable<std::vector<float>> TOFnsigmacutsPr{"TOFnsigmacutsPr", {2}, "TOF Nsigma cuts for Proton"};
  Configurable<std::vector<float>> TOFpTrangesPr{"TOFpTrangesPr", {1.1}, "TOF pT ranges for Proton"};

  bool PID(float trackpt, float tracknsigmatpc, float tracknsigmatof, int8_t species)
  {
    int8_t tpcindex = -1, tofindex = -1;
    auto tpcpt = (std::vector<std::vector<float>>){TPCpTrangesPi, TPCpTrangesKa, TPCpTrangesPr};
    auto tpcnsigma = (std::vector<std::vector<float>>){TPCnsigmacutsPi, TPCnsigmacutsKa, TPCnsigmacutsPr};
    auto tofpt = (std::vector<std::vector<float>>){TOFpTrangesPi, TOFpTrangesKa, TOFpTrangesPr};
    auto tofnsigma = (std::vector<std::vector<float>>){TOFnsigmacutsPi, TOFnsigmacutsKa, TOFnsigmacutsPr};
    for (std::size_t i = 0; i < tpcpt[species].size(); i++)
      if (trackpt < tpcpt[species][i]) {
        tpcindex = i;
        break;
      }
    for (std::size_t i = 0; i < tofpt[species].size(); i++)
      if (trackpt >= tofpt[species][i]) {
        tofindex = i;
        break;
      }
    if (tracknsigmatpc > tpcnsigma[species][tpcindex])
      return false;
    if ((tofindex != -1) && (tracknsigmatof > tofnsigma[species][tofindex]))
      return false;
    return true;
  }
  HistogramRegistry histos{"PID", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    const AxisSpec ptaxis{100, 0, 2, "p_T"}, nsigmaaxis{50, -6, 6, "N#sigma"}, dcaxyaxis{50, -6, 6, "DCA_{X}"}, dcazaxis{50, -6, 6, "DCA_{Z}"};
    histos.add("nsigmatpcpi", "N#sigma_{TPC} Pion", kTH2F, {ptaxis, nsigmaaxis});
    histos.add("nsigmatofpi", "N#sigma_{TOF} Pion", kTH2F, {ptaxis, nsigmaaxis});
    histos.add("nsigmatpcka", "N#sigma_{TPC} Kaon", kTH2F, {ptaxis, nsigmaaxis});
    histos.add("nsigmatofka", "N#sigma_{TOF} Kaon", kTH2F, {ptaxis, nsigmaaxis});
    histos.add("nsigmatpcpr", "N#sigma_{TPC} Proton", kTH2F, {ptaxis, nsigmaaxis});
    histos.add("nsigmatofpr", "N#sigma_{TOF} Proton", kTH2F, {ptaxis, nsigmaaxis});
    histos.add("dcaxypi", "DCA_{XY} Pion", kTH2F, {ptaxis, dcaxyaxis});
    histos.add("dcazpi", "DCA_{Z} Pion", kTH2F, {ptaxis, dcazaxis});
    histos.add("dcaxyka", "DCA_{XY} Kaon", kTH2F, {ptaxis, dcaxyaxis});
    histos.add("dcazka", "DCA_{Z} Kaon", kTH2F, {ptaxis, dcazaxis});
    histos.add("dcaxypr", "DCA_{XY} Proton", kTH2F, {ptaxis, dcaxyaxis});
    histos.add("dcazpr", "DCA_{Z} Proton", kTH2F, {ptaxis, dcazaxis});
    histos.add("ptpi", "p_T distribution Pion", kTH1D, {ptaxis});
    histos.add("ptka", "p_T distribution Kaon", kTH1D, {ptaxis});
    histos.add("ptpr", "p_T distribution Proton", kTH1D, {ptaxis});

    histos.add("recodcaxypi", "DCA_{XY} Pion", kTH2F, {ptaxis, dcaxyaxis});
    histos.add("recodcazpi", "DCA_{Z} Pion", kTH2F, {ptaxis, dcazaxis});
    histos.add("recodcaxyka", "DCA_{XY} Kaon", kTH2F, {ptaxis, dcaxyaxis});
    histos.add("recodcazka", "DCA_{Z} Kaon", kTH2F, {ptaxis, dcazaxis});
    histos.add("recodcaxypr", "DCA_{XY} Proton", kTH2F, {ptaxis, dcaxyaxis});
    histos.add("recodcazpr", "DCA_{Z} Proton", kTH2F, {ptaxis, dcazaxis});

    histos.add("genptpi", "Generated p_T distribution Pion", kTH1D, {ptaxis});
    histos.add("genptka", "Generated p_T distribution Kaon", kTH1D, {ptaxis});
    histos.add("genptpr", "Generated p_T distribution Proton", kTH1D, {ptaxis});
    histos.add("recoptpi", "Reconstructed p_T distribution Pion", kTH1D, {ptaxis});
    histos.add("recoptka", "Reconstructed p_T distribution Kaon", kTH1D, {ptaxis});
    histos.add("recoptpr", "Reconstructed p_T distribution Proton", kTH1D, {ptaxis});
    histos.add("pureidptpi", "Identifed w/o impurity p_T distribution Pion", kTH1D, {ptaxis});
    histos.add("pureidptka", "Identifed w/o impurity p_T distribution Kaon", kTH1D, {ptaxis});
    histos.add("pureidptpr", "Identifed w/o impurity p_T distribution Proton", kTH1D, {ptaxis});
  }
  Produces<aod::Flags> ftable;
  void processData(soa::Join<aod::Tracks, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::TracksExtra, aod::TracksDCA> const& tracks)
  {
    int8_t etabin, phibin, binNpid;
    for (auto track : tracks) {
      etabin = (track.eta() + 0.8) * 15; // 15= 24/1.6
      phibin = 36 * track.phi() / (2 * constants::math::PI);
      if ((etabin < 0) || (etabin >= 24) || (phibin < 0) || (phibin >= 36)) {
        binNpid = -1;
      } else {
        float nsigma_array[3] = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
        float tofnsigma_array[3] = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
        binNpid = 0;
        for (int8_t i = 0; i < 3; i++) {
          if (PID(track.pt(), fabs(nsigma_array[i]), fabs(tofnsigma_array[i]), i))
            binNpid = binNpid * 10 + i + 1;
          if (binNpid > 10) // If a track is identified as two different tracks.
          {
            if (fabs(nsigma_array[(binNpid / 10) - 1]) < fabs(nsigma_array[(binNpid % 10) - 1])) // The track is identified as the particle whose |nsigma| is the least.
              binNpid /= 10;
            else
              binNpid %= 10;
          }
        }

        switch (binNpid) {
          case 1:
            histos.fill(HIST("nsigmatpcpi"), track.pt(), nsigma_array[0]);
            histos.fill(HIST("nsigmatofpi"), track.pt(), tofnsigma_array[0]);
            histos.fill(HIST("dcaxypi"), track.pt(), track.dcaXY());
            histos.fill(HIST("dcazpi"), track.pt(), track.dcaZ());
            histos.fill(HIST("ptpi"), track.pt());
            break;
          case 2:
            histos.fill(HIST("nsigmatpcka"), track.pt(), nsigma_array[1]);
            histos.fill(HIST("nsigmatofka"), track.pt(), tofnsigma_array[1]);
            histos.fill(HIST("dcaxyka"), track.pt(), track.dcaXY());
            histos.fill(HIST("dcazka"), track.pt(), track.dcaZ());
            histos.fill(HIST("ptka"), track.pt());
            break;
          case 3:
            histos.fill(HIST("nsigmatpcpr"), track.pt(), nsigma_array[2]);
            histos.fill(HIST("nsigmatofpr"), track.pt(), tofnsigma_array[2]);
            histos.fill(HIST("dcaxypr"), track.pt(), track.dcaXY());
            histos.fill(HIST("dcazpr"), track.pt(), track.dcaZ());
            histos.fill(HIST("ptpr"), track.pt());
            break;
        }
      }
      ftable(binNpid);
    }
  }
  PROCESS_SWITCH(FillFlagsTable, processData, "Process Data", true);

  void processMC(soa::Join<aod::Tracks, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels> const& recotracks, aod::McParticles const& gentracks)
  {
    int8_t etabin, phibin, binNpid;
    for (auto track : recotracks) {
      if (track.has_mcParticle()) {
        etabin = (track.eta() + 0.8) * 15; // 15= 24/1.6
        phibin = 36 * track.phi() / (2 * constants::math::PI);
        if ((etabin < 0) || (etabin >= 24) || (phibin < 0) || (phibin >= 36)) {
          binNpid = -1;
        } else {
          float nsigma_array[3] = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
          float tofnsigma_array[3] = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
          binNpid = 0;
          for (int8_t i = 0; i < 3; i++) {
            if (PID(track.pt(), fabs(nsigma_array[i]), fabs(tofnsigma_array[i]), i))
              binNpid = binNpid * 10 + i + 1;
            if (binNpid > 10) // If a track is identified as two different tracks.
            {
              if (fabs(nsigma_array[(binNpid / 10) - 1]) < fabs(nsigma_array[(binNpid % 10) - 1])) // The track is identified as the particle whose |nsigma| is the least.
                binNpid /= 10;
              else
                binNpid %= 10;
            }
          }
          switch (binNpid) {
            case 1:
              histos.fill(HIST("nsigmatpcpi"), track.pt(), nsigma_array[0]);
              histos.fill(HIST("nsigmatofpi"), track.pt(), tofnsigma_array[0]);
              histos.fill(HIST("dcaxypi"), track.pt(), track.dcaXY());
              histos.fill(HIST("dcazpi"), track.pt(), track.dcaZ());
              histos.fill(HIST("ptpi"), track.pt());

              if (abs(track.mcParticle().pdgCode()) == 211)
                histos.fill(HIST("pureidptpi"), track.pt());
              break;
            case 2:
              histos.fill(HIST("nsigmatpcka"), track.pt(), nsigma_array[1]);
              histos.fill(HIST("nsigmatofka"), track.pt(), tofnsigma_array[1]);
              histos.fill(HIST("dcaxyka"), track.pt(), track.dcaXY());
              histos.fill(HIST("dcazka"), track.pt(), track.dcaZ());
              histos.fill(HIST("ptka"), track.pt());

              if (abs(track.mcParticle().pdgCode()) == 321)
                histos.fill(HIST("pureidptka"), track.pt());
              break;
            case 3:
              histos.fill(HIST("nsigmatpcpr"), track.pt(), nsigma_array[2]);
              histos.fill(HIST("nsigmatofpr"), track.pt(), tofnsigma_array[2]);
              histos.fill(HIST("dcaxypr"), track.pt(), track.dcaXY());
              histos.fill(HIST("dcazpr"), track.pt(), track.dcaZ());
              histos.fill(HIST("ptpr"), track.pt());

              if (abs(track.mcParticle().pdgCode()) == 2212)
                histos.fill(HIST("pureidptpr"), track.pt());
              break;
          }
          switch (abs(track.mcParticle().pdgCode())) {
            case 211:
              histos.fill(HIST("recoptpi"), track.pt());
              histos.fill(HIST("recodcaxypi"), track.pt(), track.dcaXY());
              histos.fill(HIST("recodcazpi"), track.pt(), track.dcaZ());
              break;
            case 321:
              histos.fill(HIST("recoptka"), track.pt());
              histos.fill(HIST("recodcaxyka"), track.pt(), track.dcaXY());
              histos.fill(HIST("recodcazka"), track.pt(), track.dcaZ());
              break;
            case 2212:
              histos.fill(HIST("recoptpr"), track.pt());
              histos.fill(HIST("recodcaxypr"), track.pt(), track.dcaXY());
              histos.fill(HIST("recodcazpr"), track.pt(), track.dcaZ());
              break;
          }
        }
      } else {
        binNpid = -1;
      }
      ftable(binNpid);
    }
    for (auto track : gentracks) {
      switch (abs(track.pdgCode())) {
        case 211:
          histos.fill(HIST("genptpi"), track.pt());
          break;
        case 321:
          histos.fill(HIST("genptka"), track.pt());
          break;
        case 2212:
          histos.fill(HIST("genptpr"), track.pt());
          break;
      }
    }
  }
  PROCESS_SWITCH(FillFlagsTable, processMC, "Process MC Data", false);
};

struct r2p24id {

  Configurable<float> minpT{"minpT", 0.2, "Minimum pT"};
  Configurable<float> maxpT{"maxpT", 2.0, "Maximum pT"};
  Configurable<float> trackpartition{"trackpartition", 1.0, "where(in pT) to partition"};

  Configurable<int> pid_particle1{"pid_particle1", 1, "Define particle1 type"}; // 1->Pion, 2->Kaon, 3->Proton
  Configurable<int> pid_particle2{"pid_particle2", 1, "Define particle2 type"};

  Configurable<bool> iftrackpartition{"iftrackpartition", false, "If track partition is needed"};
  Configurable<bool> ifpid{"ifpid", false, "If PID is needed"};

  struct histarray {
    std::shared_ptr<TH2> h2d_1p[2][2][2], h2d_2p[2][2][4];
    std::shared_ptr<TH1> h1d_1p[2][2];
  } hist;

  unsigned int mult1, mult2;
  uint8_t etabin1, phibin1, etabin2, phibin2;
  int8_t sign1, sign2;
  bool iftrack2;

  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  SliceCache cache;

  //-----Track&Event Selection-------------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);
  Filter ptfilter = aod::track::pt > minpT&& aod::track::pt < maxpT;
  Filter globalfilter = requireGlobalTrackInFilter();
  Filter properbinfilter = aod::idr2p2columns::binNpid != -1;

  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::Flags>>> Tracks_set1 = (ifnode(iftrackpartition.node() == true, (aod::track::pt < trackpartition), true) && ifnode(ifpid.node() == true, (aod::idr2p2columns::binNpid == pid_particle1), true));
  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::Flags>>> Tracks_set2 = (ifnode(iftrackpartition.node() == true, (aod::track::pt > trackpartition), true) && ifnode(ifpid.node() == true, (aod::idr2p2columns::binNpid == pid_particle2), true));
  //---------------------------------------------------------------------------

  void init(InitContext const&)
  {
    iftrack2 = (((int8_t)pid_particle1 != (int8_t)pid_particle2) && (static_cast<bool>(ifpid))) || (static_cast<bool>(iftrackpartition)); // denotes whether the partiton1 is different from partition2
    //-----Defining Histograms---------------------------------------------------
    const AxisSpec phi{36, 0, 2.0 * constants::math::PI, "phi"}, eta{24, -0.8, 0.8, "eta"}, etaphi1{864, 0, 864, "etaphi1"}, etaphi2{864, 0, 864, "etaphi2"};
    histos.add("h1d_n1_phi", "#phi distribution Particle", kTH1D, {phi});
    histos.add("h1d_n1_eta", "#eta distribution Particle", kTH1D, {eta});
    histos.add("h1d_n1_ptP1", "p_T for +ve_1", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM1", "p_T for -ve_1", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptP2", "p_T for +ve_2", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM2", "p_T for -ve_2", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1i_n1_multPM", "Multiplicity", kTH1I, {{200, 0, 200, "Multiplicity"}});
    histos.add("h2d_n1_etaPhiP1", "#rho_1 for +ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM1", "#rho_1 for -ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiP2", "#rho_1 for +ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM2", "#rho_1 for -ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PP", "#rho_2 for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM12", "#rho_2 for +ve_1 -ve_2 particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM21", "#rho_2 for +ve_2 -ve_1 particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2MM", "#rho_2 for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_pt_etaPhiP1", "p_T for +ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM1", "p_T for -ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiP2", "p_T for +ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM2", "p_T for -ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM12", "p_Tp_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM21", "p_Tp_T for +ve_2 -ve_1", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP", "p_Tp_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM", "p_Tp_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM12", "p_Tn for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM21", "p_Tn for +ve_2 -ve_1", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP", "p_Tn for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM", "p_Tn for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM12", "np_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM21", "np_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PP", "np_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2MM", "np_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h1d_n1_pt", "p_T", kTH1D, {{100, 0, 5, "p_T"}});
    //---------------------------------------------------------------------------
    //-----Histogram Arrays------------------------------------------------------
    hist.h2d_1p[0][0][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiM1"));
    hist.h2d_1p[0][0][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiM1"));
    hist.h2d_1p[0][1][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiP1"));
    hist.h2d_1p[0][1][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiP1"));
    hist.h2d_1p[1][0][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiM2"));
    hist.h2d_1p[1][0][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiM2"));
    hist.h2d_1p[1][1][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiP2"));
    hist.h2d_1p[1][1][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiP2"));

    hist.h2d_2p[0][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"));
    hist.h2d_2p[0][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"));
    hist.h2d_2p[0][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"));
    hist.h2d_2p[0][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"));
    if (iftrack2) {
      hist.h2d_2p[0][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM21"));
      hist.h2d_2p[0][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM21"));
      hist.h2d_2p[0][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM21"));
      hist.h2d_2p[0][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM21"));
    } else {
      hist.h2d_2p[0][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM12"));
      hist.h2d_2p[0][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM12"));
      hist.h2d_2p[0][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM12"));
      hist.h2d_2p[0][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM12"));
    }
    hist.h2d_2p[1][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM12"));
    hist.h2d_2p[1][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM12"));
    hist.h2d_2p[1][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM12"));
    hist.h2d_2p[1][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM12"));
    hist.h2d_2p[1][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"));
    hist.h2d_2p[1][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"));
    hist.h2d_2p[1][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"));
    hist.h2d_2p[1][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"));

    hist.h1d_1p[0][0] = histos.template get<TH1>(HIST("h1d_n1_ptP1"));
    hist.h1d_1p[0][1] = histos.template get<TH1>(HIST("h1d_n1_ptM1"));
    hist.h1d_1p[1][0] = histos.template get<TH1>(HIST("h1d_n1_ptP2"));
    hist.h1d_1p[1][1] = histos.template get<TH1>(HIST("h1d_n1_ptM2"));
    //---------------------------------------------------------------------------
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::Flags>> const& tracks)
  {
    auto tracks1 = Tracks_set1->sliceByCached(aod::track::collisionId, filteredCollision.globalIndex(), cache);
    auto tracks2 = Tracks_set2->sliceByCached(aod::track::collisionId, filteredCollision.globalIndex(), cache);
    mult1 = tracks1.size();
    mult2 = tracks2.size();
    if ((iftrack2 && ((mult1 < 1) || (mult2 < 1))) || ((!iftrack2) && (mult1 < 2))) // Reject Collisions without sufficient particles
      return;

    for (auto track1 : tracks) {
      histos.fill(HIST("h1d_n1_phi"), track1.phi());
      histos.fill(HIST("h1d_n1_eta"), track1.eta());
      histos.fill(HIST("h1d_n1_pt"), track1.pt());
    }
    histos.fill(HIST("h1i_n1_multPM"), mult1 + mult2);

    for (auto track1 : tracks1) {
      //---Single Particle Distribution (particle1)----------------------------------------
      sign1 = (track1.sign() + 1) / 2;
      hist.h1d_1p[0][sign1]->Fill(track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt())); // h1d_n1_pt*1
      hist.h2d_1p[0][sign1][0]->Fill(track1.eta(), track1.phi());                                // h2d_n1_etaphi*1
      hist.h2d_1p[0][sign1][1]->Fill(track1.eta(), track1.phi(), track1.pt());                   // h2d_pt_etaPhi*1
      //-----------------------------------------------------------------------
    }
    if (iftrack2) {
      for (auto track2 : tracks2) {
        //---Single Particle Distribution (particle2)----------------------------------------
        sign2 = (track2.sign() + 1) / 2;
        hist.h1d_1p[1][sign2]->Fill(track2.pt(), 1.0 / (2.0 * constants::math::PI * track2.pt())); // h1d_n1_pt*2
        hist.h2d_1p[1][sign2][0]->Fill(track2.eta(), track2.phi());                                // h2d_n1_etaphi*2
        hist.h2d_1p[1][sign2][1]->Fill(track2.eta(), track2.phi(), track2.pt());                   // h2d_pt_etaPhi*2
        //-----------------------------------------------------------------------
      }
    }
    for (auto track1 : tracks1) {
      sign1 = (track1.sign() + 1) / 2;
      etabin1 = (track1.eta() + 0.8) * 15; // 15= 24/1.6
      phibin1 = 36 * track1.phi() / (2 * constants::math::PI);

      for (auto track2 : tracks2) {

        if (track1.index() == track2.index())
          continue;

        etabin2 = (track2.eta() + 0.8) * 15; // 15= 24/1.6
        phibin2 = 36 * track2.phi() / (2 * constants::math::PI);
        sign2 = (track2.sign() + 1) / 2;

        //-----Two Particle Distribution---------------------------------------
        hist.h2d_2p[sign1][sign2][0]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);                            // h2d_n2_eta1Phi1Eta2Phi2*
        hist.h2d_2p[sign1][sign2][1]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());               // h2d_npt_eta1Phi1Eta2Phi2*
        hist.h2d_2p[sign1][sign2][2]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());               // h2d_ptn_eta1Phi1Eta2Phi2*
        hist.h2d_2p[sign1][sign2][3]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt()); // h2d_ptpt_eta1Phi1Eta2Phi2*
        //---------------------------------------------------------------------
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FillFlagsTable>(cfgc),
    adaptAnalysisTask<r2p24id>(cfgc),
  };
}
