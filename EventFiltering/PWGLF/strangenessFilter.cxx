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
//
/// \brief A filter task for strangeness triggers
//  usage:
/*
  o2-analysis-timestamp -b --aod-file AO2D.root   | \
  o2-analysis-event-selection -b | \
  o2-analysis-trackselection -b | \
  o2-analysis-trackextension -b | \
  o2-analysis-multiplicity-table -b | \
  o2-analysis-centrality-table -b | \
  o2-analysis-pid-tof -b | \
  o2-analysis-pid-tpc -b | \
  o2-analysis-weak-decay-indices -b | \
  o2-analysis-lambdakzerobuilder  --d_bz 5 -b | \
  o2-analysis-cascadebuilder  --d_bz 5 -b | \
  o2-analysis-strangeness-filter -b
*/
///
///
/// \author Chiara De Martin (chiara.de.martin@cern.ch)
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \since June 1, 2021

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/StrangenessTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

#include "../filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct strangenessFilter {

  // Recall the output table
  Produces<aod::StrangenessFilters> strgtable;

  // Define a histograms and registries
  HistogramRegistry QAHistos{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry EventsvsMultiplicity{"EventsvsMultiplicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1F> hProcessedEvents{TH1F("hProcessedEvents", "Strangeness - event filtered; Event counter; Number of events", 7, 0., 7.)};

  // Selection criteria for cascades
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> v0cospa{"v0cospa", 0.97, "V0 CosPA"}; // is it with respect to Xi decay vertex?
  Configurable<float> casccospa{"casccospa", 0.995, "V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};       // is it in sigmas?
  Configurable<float> dcacascdau{"dcacascdau", 0.8, "DCA Casc Daughters"}; // is it in sigmas?
  Configurable<float> dcamesontopv{"dcamesontopv", 0.04, "DCA Meson To PV"};
  Configurable<float> dcabaryontopv{"dcabaryontopv", 0.03, "DCA Baryon To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", 0.04, "DCA Bach To PV"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.02, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.02, "DCA Pos To PV"};
  Configurable<float> dcav0topv{"dcav0topv", 1.2, "DCA V0 To PV"};
  Configurable<float> v0radius{"v0radius", 1.2, "V0 Radius"};
  Configurable<float> v0radiusupperlimit{"v0radiusupperlimit", 34, "V0 Radius Upper Limit"};
  Configurable<float> cascradius{"cascradius", 0.6, "cascradius"};
  Configurable<float> cascradiusupperlimit{"cascradiusupperlimit", 34, "Casc Radius Upper Limit"};
  Configurable<float> rapidity{"rapidity", 2, "rapidity"};
  Configurable<float> eta{"eta", 2, "Eta"};
  Configurable<float> minpt{"minpt", 0.5, "minpt"};
  Configurable<float> etadau{"etadau", 0.8, "EtaDaughters"};
  Configurable<float> masslambdalimit{"masslambdalimit", 0.01, "masslambdalimit"}; // 0.006 Chiara
  Configurable<float> omegarej{"omegarej", 0.005, "omegarej"};
  Configurable<float> xirej{"xirej", 0.008, "xirej"}; // merge the two rejection variables into one?
  Configurable<float> ximasswindow{"ximasswindow", 0.075, "Xi Mass Window"};
  Configurable<float> omegamasswindow{"omegamasswindow", 0.075, "Omega Mass Window"}; // merge the two windows variables into one?
  Configurable<int> properlifetimefactor{"properlifetimefactor", 5, "Proper Lifetime cut"};
  Configurable<float> nsigmatpc{"nsigmatpc", 6, "N Sigmas TPC"};
  Configurable<float> nsigmatof{"nsigmatof", 5, "N Sigmas TOF (OOB condition)"};
  Configurable<bool> kint7{"kint7", 0, "Apply kINT7 event selection"};
  Configurable<bool> sel7{"sel7", 0, "Apply sel7 event selection"};
  Configurable<bool> sel8{"sel8", 0, "Apply sel8 event selection"};
  Configurable<bool> globaltrk{"globaltrk", 1, "Apply global track selection"};

  // Selections criteria for tracks
  Configurable<float> hEta{"hEta", 0.8f, "Eta range for trigger particles"};
  Configurable<float> hMinPt{"hMinPt", 1.0f, "Min pt for trigger particles"};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    AxisSpec centAxis = {centBinning, "V0M (%)"};
    AxisSpec ximassAxis = {100, 1.30f, 1.34f};
    AxisSpec omegamassAxis = {100, 1.5f, 1.8f};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};

    QAHistos.add("hCentrality", "Centrality distribution (V0M)", HistType::kTH1F, {{100, 0, 100, "V0M (%)"}});
    QAHistos.add("hVtxZAfterSel", "Vertex distribution in Z;Z (cm)", HistType::kTH1F, {{100, -20, 20}});
    QAHistos.add("hMassXiBefSel", "hMassXiBefSel", HistType::kTH1F, {ximassAxis});
    QAHistos.add("hMassXiAfterSel", "hMassXiAfterSel", HistType::kTH1F, {ximassAxis});
    QAHistos.add("hMassOmegaBefSel", "hMassOmegaBefSel", HistType::kTH1F, {omegamassAxis});
    QAHistos.add("hMassOmegaAfterSel", "hMassOmegaAfterSel", HistType::kTH1F, {omegamassAxis});
    QAHistos.add("hTOFnsigmaV0PiBefSel", "hTOFnsigmaV0PiBefSel", HistType::kTH1F, {{100, -10, +10, "TOFnsigmaV0PiBefSel"}});
    QAHistos.add("hTOFnsigmaV0PiAfterSel", "hTOFnsigmaV0PiAfterSel", HistType::kTH1F, {{100, -10, +10, "TOFnsigmaV0PiAfterSel"}});
    QAHistos.add("hTOFnsigmaPrBefSel", "hTOFnsigmaPrBefSel", HistType::kTH1F, {{100, -10, +10, "TOFnsigmaPrBefSel"}});
    QAHistos.add("hTOFnsigmaPrAfterSel", "hTOFnsigmaPrAfterSel", HistType::kTH1F, {{100, -10, +10, "TOFnsigmaPrAfterSel"}});
    QAHistos.add("hMassXiAfterSelvsPt", "hMassXiAfterSelvsPt", HistType::kTH2F, {ximassAxis, ptAxis});
    QAHistos.add("hMassOmegaAfterSelvsPt", "hMassOmegaAfterSelvsPt", HistType::kTH2F, {omegamassAxis, ptAxis});
    QAHistos.add("hTriggeredParticles", "Selected triggered particles", HistType::kTH1F, {{10, 0.5, 10.5, "Trigger counter"}});
    QAHistos.add("hPtTrigger", "hPtTrigger", HistType::kTH1F, {{300, 0, 30, "Pt of trigger particle"}});

    EventsvsMultiplicity.add("AllEventsvsMultiplicity", "Multiplicity distribution of all events", HistType::kTH1F, {centAxis});
    EventsvsMultiplicity.add("OmegaEventsvsMultiplicity", "Multiplicity distribution of events with >= 1 Omega", HistType::kTH1F, {centAxis});
    EventsvsMultiplicity.add("hXiEventsvsMultiplicity", "Multiplicity distribution of events with h + Xi", HistType::kTH1F, {centAxis});
    EventsvsMultiplicity.add("2XiEventsvsMultiplicity", "Multiplicity distribution of events with >= 2 Xi", HistType::kTH1F, {centAxis});
    EventsvsMultiplicity.add("3XiEventsvsMultiplicity", "Multiplicity distribution of events with >= 3 Xi", HistType::kTH1F, {centAxis});
    EventsvsMultiplicity.add("4XiEventsvsMultiplicity", "Multiplicity distribution of events with >= 4 Xi", HistType::kTH1F, {centAxis});
    EventsvsMultiplicity.add("SingleXiEventsvsMultiplicity", "Multiplicity distribution of events with 1 Xi (R > 24.39 cm)", HistType::kTH1F, {centAxis});

    hProcessedEvents->GetXaxis()->SetBinLabel(1, "Events processed");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "#Omega");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "high-#it{p}_{T} hadron - #Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(4, "2#Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(5, "3#Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(6, "4#Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(7, "single-#Xi");
  }

  // Filters
  Filter trackFilter = (nabs(aod::track::eta) < hEta) && (aod::track::pt > hMinPt) && (!globaltrk || requireGlobalTrackInFilter());

  // Tables
  using CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator;
  using CollisionCandidatesRun3 = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFPi, aod::pidTPCPi, aod::pidTOFPr, aod::pidTPCPr>;
  using Cascades = aod::CascDataExt;

  ////////////////////////////////////////////////////////
  ////////// Strangeness Filter - Run 2 conv /////////////
  ////////////////////////////////////////////////////////

  void processRun2(CollisionCandidates const& collision, TrackCandidates const& tracks, Cascades const& fullCasc, aod::V0sLinked const&, aod::V0Datas const& v0data, DaughterTracks& dtracks)
  {
    if (kint7 && !collision.alias()[kINT7]) {
      return;
    }
    if (sel7 && !collision.sel7()) {
      return;
    }
    if (sel8 && !collision.sel8()) {
      return;
    }

    if (TMath::Abs(collision.posZ()) > cutzvertex)
      return;

    QAHistos.fill(HIST("hVtxZAfterSel"), collision.posZ());
    QAHistos.fill(HIST("hCentrality"), collision.centRun2V0M());
    EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicity"), collision.centRun2V0M());
    hProcessedEvents->Fill(0.5);

    // Is event good? [0] = Omega, [1] = high-pT hadron + Xi, [2] = 2Xi, [3] = 3Xi, [4] = 4Xi, [5] single-Xi
    bool keepEvent[6]{false};

    // constants
    const float ctauxi = 4.91;     // from PDG
    const float ctauomega = 2.461; // from PDG

    // variables
    float xipos = -1.;
    float xiproperlifetime = -1.;
    float omegaproperlifetime = -1.;
    float xiptotmom = -1.;
    int xicounter = 0;
    int xicounterYN = 0;
    int omegacounter = 0;

    for (auto& casc : fullCasc) { // loop over cascades

      auto v0index = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0index.has_v0Data())) {
        continue; // skip those cascades for which V0 doesn't exist
      }
      auto v0 = v0index.v0Data(); // de-reference index to correct v0data in case it exists
      auto bachelor = casc.bachelor_as<DaughterTracks>();
      auto posdau = v0.posTrack_as<DaughterTracks>();
      auto negdau = v0.negTrack_as<DaughterTracks>();

      bool isXi = false;
      bool isXiYN = false;
      bool isOmega = false;

      // QA
      QAHistos.fill(HIST("hMassXiBefSel"), casc.mXi());
      QAHistos.fill(HIST("hMassOmegaBefSel"), casc.mOmega());

      // Position
      xipos = std::hypot(casc.x() - collision.posX(), casc.y() - collision.posY(), casc.z() - collision.posZ());
      // Total momentum
      xiptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      // Proper lifetime
      xiproperlifetime = RecoDecay::getMassPDG(3312) * xipos / (xiptotmom + 1e-13);
      omegaproperlifetime = RecoDecay::getMassPDG(3334) * xipos / (xiptotmom + 1e-13);

      if (casc.sign() == 1) {
        if (TMath::Abs(casc.dcapostopv()) < dcamesontopv) {
          continue;
        };
        if (TMath::Abs(casc.dcanegtopv()) < dcabaryontopv) {
          continue;
        };
        if (TMath::Abs(posdau.tpcNSigmaPi()) > nsigmatpc) {
          continue;
        };
        if (TMath::Abs(negdau.tpcNSigmaPr()) > nsigmatpc) {
          continue;
        };
        QAHistos.fill(HIST("hTOFnsigmaPrBefSel"), negdau.tofNSigmaPr());
        QAHistos.fill(HIST("hTOFnsigmaV0PiBefSel"), posdau.tofNSigmaPi());
        if (
          (TMath::Abs(posdau.tofNSigmaPi()) > nsigmatof) &&
          (TMath::Abs(negdau.tofNSigmaPr()) > nsigmatof) &&
          (TMath::Abs(bachelor.tofNSigmaPi()) > nsigmatof)) {
          continue;
        };
        QAHistos.fill(HIST("hTOFnsigmaPrAfterSel"), negdau.tofNSigmaPr());
        QAHistos.fill(HIST("hTOFnsigmaV0PiAfterSel"), posdau.tofNSigmaPi());
      } else {
        if (TMath::Abs(casc.dcanegtopv()) < dcamesontopv) {
          continue;
        };
        if (TMath::Abs(casc.dcapostopv()) < dcabaryontopv) {
          continue;
        };
        if (TMath::Abs(posdau.tpcNSigmaPr()) > nsigmatpc) {
          continue;
        };
        if (TMath::Abs(negdau.tpcNSigmaPi()) > nsigmatpc) {
          continue;
        };
        QAHistos.fill(HIST("hTOFnsigmaPrBefSel"), posdau.tofNSigmaPr());
        QAHistos.fill(HIST("hTOFnsigmaV0PiBefSel"), negdau.tofNSigmaPi());
        if (
          (TMath::Abs(posdau.tofNSigmaPr()) > nsigmatof) &&
          (TMath::Abs(negdau.tofNSigmaPi()) > nsigmatof) &&
          (TMath::Abs(bachelor.tofNSigmaPi()) > nsigmatof)) {
          continue;
        };
        QAHistos.fill(HIST("hTOFnsigmaPrAfterSel"), posdau.tofNSigmaPr());
        QAHistos.fill(HIST("hTOFnsigmaV0PiAfterSel"), negdau.tofNSigmaPi());
      }
      // this selection differes for Xi and Omegas:
      if (TMath::Abs(bachelor.tpcNSigmaPi()) > nsigmatpc) {
        continue;
      };
      if (TMath::Abs(posdau.eta()) > etadau) {
        continue;
      };
      if (TMath::Abs(negdau.eta()) > etadau) {
        continue;
      };
      if (TMath::Abs(bachelor.eta()) > etadau) {
        continue;
      };
      if (TMath::Abs(casc.dcabachtopv()) < dcabachtopv) {
        continue;
      };
      if (casc.v0radius() > v0radiusupperlimit || casc.v0radius() < v0radius) {
        continue;
      };
      if (casc.cascradius() > cascradiusupperlimit || casc.cascradius() < cascradius) {
        continue;
      };
      if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) {
        continue;
      };
      if (casc.dcaV0daughters() > dcav0dau) {
        continue;
      };
      if (casc.dcacascdaughters() > dcacascdau) {
        continue;
      };
      if (TMath::Abs(casc.mLambda() - constants::physics::MassLambda) > masslambdalimit) {
        continue;
      };
      if (TMath::Abs(casc.eta()) > eta) {
        continue;
      };

      isXi = (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospa) &&
             (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) &&
             (TMath::Abs(casc.mXi() - RecoDecay::getMassPDG(3312)) < ximasswindow) &&
             (TMath::Abs(casc.mOmega() - RecoDecay::getMassPDG(3334)) > omegarej) &&
             (xiproperlifetime < properlifetimefactor * ctauxi) &&
             (TMath::Abs(casc.yXi()) < rapidity); // add PID on bachelor
      isXiYN = (casc.cascradius() > 24.39) &&
               (TMath::Abs(casc.mXi() - RecoDecay::getMassPDG(3312)) < ximasswindow) &&
               (TMath::Abs(casc.mOmega() - RecoDecay::getMassPDG(3334)) > omegarej) &&
               (xiproperlifetime < properlifetimefactor * ctauxi) &&
               (TMath::Abs(casc.yXi()) < rapidity); // add PID on bachelor
      isOmega = (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospa) &&
                (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) &&
                (TMath::Abs(casc.mOmega() - RecoDecay::getMassPDG(3334)) < omegamasswindow) &&
                (TMath::Abs(casc.mXi() - RecoDecay::getMassPDG(3312)) > xirej) &&
                (omegaproperlifetime < properlifetimefactor * ctauomega) &&
                (TMath::Abs(casc.yOmega()) < rapidity); // add PID on bachelor

      if (isXi) {
        QAHistos.fill(HIST("hMassXiAfterSel"), casc.mXi());
        QAHistos.fill(HIST("hMassXiAfterSelvsPt"), casc.mXi(), casc.pt());
        // Count number of Xi candidates
        xicounter++;
      }
      if (isXiYN) {
        // Xis for YN interactions
        xicounterYN++;
      }
      if (isOmega) {
        QAHistos.fill(HIST("hMassOmegaAfterSel"), casc.mOmega());
        QAHistos.fill(HIST("hMassOmegaAfterSelvsPt"), casc.mOmega(), casc.pt());
        // Count number of Omega candidates
        omegacounter++;
      }
    } // end loop over cascades

    // Omega trigger definition
    if (omegacounter > 0) {
      keepEvent[0] = true;
    }

    // High-pT hadron + Xi trigger definition
    if (xicounter > 0) {
      for (auto track : tracks) { // start loop over tracks
        QAHistos.fill(HIST("hTriggeredParticles"), 1);
        QAHistos.fill(HIST("hPtTrigger"), track.pt());

        keepEvent[1] = true;
      } // end loop over tracks
    }

    // 2Xi trigger definition
    if (xicounter > 1) {
      keepEvent[2] = true;
    }

    // 3Xi trigger definition
    if (xicounter > 2) {
      keepEvent[3] = true;
    }

    // 4Xi trigger definition
    if (xicounter > 3) {
      keepEvent[4] = true;
    }

    // Single-Xi (YN) trigger definition
    if (xicounterYN > 0) {
      keepEvent[5] = true;
    }

    // Fill centrality dependent histos
    if (keepEvent[0]) {
      hProcessedEvents->Fill(1.5);
      EventsvsMultiplicity.fill(HIST("OmegaEventsvsMultiplicity"), collision.centRun2V0M());
    }
    if (keepEvent[1]) {
      hProcessedEvents->Fill(2.5);
      EventsvsMultiplicity.fill(HIST("hXiEventsvsMultiplicity"), collision.centRun2V0M());
    }
    if (keepEvent[2]) {
      hProcessedEvents->Fill(3.5);
      EventsvsMultiplicity.fill(HIST("2XiEventsvsMultiplicity"), collision.centRun2V0M());
    }
    if (keepEvent[3]) {
      hProcessedEvents->Fill(4.5);
      EventsvsMultiplicity.fill(HIST("3XiEventsvsMultiplicity"), collision.centRun2V0M());
    }
    if (keepEvent[4]) {
      hProcessedEvents->Fill(5.5);
      EventsvsMultiplicity.fill(HIST("4XiEventsvsMultiplicity"), collision.centRun2V0M());
    }
    if (keepEvent[5]) {
      hProcessedEvents->Fill(6.5);
      EventsvsMultiplicity.fill(HIST("SingleXiEventsvsMultiplicity"), collision.centRun2V0M());
    }

    // Filling the table
    strgtable(keepEvent[0], keepEvent[1], keepEvent[2], keepEvent[3], keepEvent[4], keepEvent[5]);
  }
  //
  PROCESS_SWITCH(strangenessFilter, processRun2, "Process data Run2", true);

  //////////////////////////////////////////////////////
  ////////// Strangeness Filter - Run 3 MC /////////////
  //////////////////////////////////////////////////////

  void processRun3(CollisionCandidatesRun3 const& collision, TrackCandidates const& tracks, Cascades const& fullCasc, aod::V0sLinked const&, aod::V0Datas const& v0data, DaughterTracks& dtracks)
  {
    if (sel8 && !collision.sel8()) {
      return;
    }

    if (TMath::Abs(collision.posZ()) > cutzvertex)
      return;

    QAHistos.fill(HIST("hVtxZAfterSel"), collision.posZ());
    hProcessedEvents->Fill(0.5);

    // Is event good? [0] = Omega, [1] = high-pT hadron + Xi, [2] = 2Xi, [3] = 3Xi, [4] = 4Xi, [5] single-Xi
    bool keepEvent[6]{false};

    // constants
    const float ctauxi = 4.91;     // from PDG
    const float ctauomega = 2.461; // from PDG

    // variables
    float xipos = -1.;
    float xiproperlifetime = -1.;
    float omegaproperlifetime = -1.;
    float xiptotmom = -1.;
    int xicounter = 0;
    int xicounterYN = 0;
    int omegacounter = 0;

    for (auto& casc : fullCasc) { // loop over cascades

      auto v0index = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0index.has_v0Data())) {
        continue; // skip those cascades for which V0 doesn't exist
      }
      auto v0 = v0index.v0Data(); // de-reference index to correct v0data in case it exists
      auto bachelor = casc.bachelor_as<DaughterTracks>();
      auto posdau = v0.posTrack_as<DaughterTracks>();
      auto negdau = v0.negTrack_as<DaughterTracks>();

      bool isXi = false;
      bool isXiYN = false;
      bool isOmega = false;

      // QA
      QAHistos.fill(HIST("hMassXiBefSel"), casc.mXi());
      QAHistos.fill(HIST("hMassOmegaBefSel"), casc.mOmega());

      // Position
      xipos = std::hypot(casc.x() - collision.posX(), casc.y() - collision.posY(), casc.z() - collision.posZ());
      // Total momentum
      xiptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      // Proper lifetime
      xiproperlifetime = RecoDecay::getMassPDG(3312) * xipos / (xiptotmom + 1e-13);
      omegaproperlifetime = RecoDecay::getMassPDG(3334) * xipos / (xiptotmom + 1e-13);

      if (casc.sign() == 1) {
        if (TMath::Abs(casc.dcapostopv()) < dcamesontopv) {
          continue;
        };
        if (TMath::Abs(casc.dcanegtopv()) < dcabaryontopv) {
          continue;
        };
        if (TMath::Abs(posdau.tpcNSigmaPi()) > nsigmatpc) {
          continue;
        };
        if (TMath::Abs(negdau.tpcNSigmaPr()) > nsigmatpc) {
          continue;
        };
        QAHistos.fill(HIST("hTOFnsigmaPrBefSel"), negdau.tofNSigmaPr());
        QAHistos.fill(HIST("hTOFnsigmaV0PiBefSel"), posdau.tofNSigmaPi());
        if (
          (TMath::Abs(posdau.tofNSigmaPi()) > nsigmatof) &&
          (TMath::Abs(negdau.tofNSigmaPr()) > nsigmatof) &&
          (TMath::Abs(bachelor.tofNSigmaPi()) > nsigmatof)) {
          continue;
        };
        QAHistos.fill(HIST("hTOFnsigmaPrAfterSel"), negdau.tofNSigmaPr());
        QAHistos.fill(HIST("hTOFnsigmaV0PiAfterSel"), posdau.tofNSigmaPi());
      } else {
        if (TMath::Abs(casc.dcanegtopv()) < dcamesontopv) {
          continue;
        };
        if (TMath::Abs(casc.dcapostopv()) < dcabaryontopv) {
          continue;
        };
        if (TMath::Abs(posdau.tpcNSigmaPr()) > nsigmatpc) {
          continue;
        };
        if (TMath::Abs(negdau.tpcNSigmaPi()) > nsigmatpc) {
          continue;
        };
        QAHistos.fill(HIST("hTOFnsigmaPrBefSel"), posdau.tofNSigmaPr());
        QAHistos.fill(HIST("hTOFnsigmaV0PiBefSel"), negdau.tofNSigmaPi());
        if (
          (TMath::Abs(posdau.tofNSigmaPr()) > nsigmatof) &&
          (TMath::Abs(negdau.tofNSigmaPi()) > nsigmatof) &&
          (TMath::Abs(bachelor.tofNSigmaPi()) > nsigmatof)) {
          continue;
        };
        QAHistos.fill(HIST("hTOFnsigmaPrAfterSel"), posdau.tofNSigmaPr());
        QAHistos.fill(HIST("hTOFnsigmaV0PiAfterSel"), negdau.tofNSigmaPi());
      }
      // this selection differes for Xi and Omegas:
      if (TMath::Abs(bachelor.tpcNSigmaPi()) > nsigmatpc) {
        continue;
      };
      if (TMath::Abs(posdau.eta()) > etadau) {
        continue;
      };
      if (TMath::Abs(negdau.eta()) > etadau) {
        continue;
      };
      if (TMath::Abs(bachelor.eta()) > etadau) {
        continue;
      };
      if (TMath::Abs(casc.dcabachtopv()) < dcabachtopv) {
        continue;
      };
      if (casc.v0radius() > v0radiusupperlimit || casc.v0radius() < v0radius) {
        continue;
      };
      if (casc.cascradius() > cascradiusupperlimit || casc.cascradius() < cascradius) {
        continue;
      };
      if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) {
        continue;
      };
      if (casc.dcaV0daughters() > dcav0dau) {
        continue;
      };
      if (casc.dcacascdaughters() > dcacascdau) {
        continue;
      };
      if (TMath::Abs(casc.mLambda() - constants::physics::MassLambda) > masslambdalimit) {
        continue;
      };
      if (TMath::Abs(casc.eta()) > eta) {
        continue;
      };

      isXi = (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospa) &&
             (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) &&
             (TMath::Abs(casc.mXi() - RecoDecay::getMassPDG(3312)) < ximasswindow) &&
             (TMath::Abs(casc.mOmega() - RecoDecay::getMassPDG(3334)) > omegarej) &&
             (xiproperlifetime < properlifetimefactor * ctauxi) &&
             (TMath::Abs(casc.yXi()) < rapidity); // add PID on bachelor
      isXiYN = (casc.cascradius() > 24.39) &&
               (TMath::Abs(casc.mXi() - RecoDecay::getMassPDG(3312)) < ximasswindow) &&
               (TMath::Abs(casc.mOmega() - RecoDecay::getMassPDG(3334)) > omegarej) &&
               (xiproperlifetime < properlifetimefactor * ctauxi) &&
               (TMath::Abs(casc.yXi()) < rapidity); // add PID on bachelor
      isOmega = (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospa) &&
                (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) &&
                (TMath::Abs(casc.mOmega() - RecoDecay::getMassPDG(3334)) < omegamasswindow) &&
                (TMath::Abs(casc.mXi() - RecoDecay::getMassPDG(3312)) > xirej) &&
                (omegaproperlifetime < properlifetimefactor * ctauomega) &&
                (TMath::Abs(casc.yOmega()) < rapidity); // add PID on bachelor

      if (isXi) {
        QAHistos.fill(HIST("hMassXiAfterSel"), casc.mXi());
        QAHistos.fill(HIST("hMassXiAfterSelvsPt"), casc.mXi(), casc.pt());
        // Count number of Xi candidates
        xicounter++;
      }
      if (isXiYN) {
        // Xis for YN interactions
        xicounterYN++;
      }
      if (isOmega) {
        QAHistos.fill(HIST("hMassOmegaAfterSel"), casc.mOmega());
        QAHistos.fill(HIST("hMassOmegaAfterSelvsPt"), casc.mOmega(), casc.pt());
        // Count number of Omega candidates
        omegacounter++;
      }
    } // end loop over cascades

    // Omega trigger definition
    if (omegacounter > 0) {
      keepEvent[0] = true;
    }

    // High-pT hadron + Xi trigger definition
    if (xicounter > 0) {
      for (auto track : tracks) { // start loop over tracks
        QAHistos.fill(HIST("hTriggeredParticles"), 1);
        QAHistos.fill(HIST("hPtTrigger"), track.pt());

        keepEvent[1] = true;
      } // end loop over tracks
    }

    // 2Xi trigger definition
    if (xicounter > 1) {
      keepEvent[2] = true;
    }

    // 3Xi trigger definition
    if (xicounter > 2) {
      keepEvent[3] = true;
    }

    // 4Xi trigger definition
    if (xicounter > 3) {
      keepEvent[4] = true;
    }

    // Single-Xi (YN) trigger definition
    if (xicounterYN > 0) {
      keepEvent[5] = true;
    }

    // Fill centrality dependent histos
    if (keepEvent[0]) {
      hProcessedEvents->Fill(1.5);
    }
    if (keepEvent[1]) {
      hProcessedEvents->Fill(2.5);
    }
    if (keepEvent[2]) {
      hProcessedEvents->Fill(3.5);
    }
    if (keepEvent[3]) {
      hProcessedEvents->Fill(4.5);
    }
    if (keepEvent[4]) {
      hProcessedEvents->Fill(5.5);
    }
    if (keepEvent[5]) {
      hProcessedEvents->Fill(6.5);
    }

    // Filling the table
    strgtable(keepEvent[0], keepEvent[1], keepEvent[2], keepEvent[3], keepEvent[4], keepEvent[5]);
  }
  //
  PROCESS_SWITCH(strangenessFilter, processRun3, "Process Run3", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenessFilter>(cfgc, TaskName{"lf-strangeness-filter"})};
}
