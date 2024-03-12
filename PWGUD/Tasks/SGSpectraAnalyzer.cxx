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
// \Single Gap Event Analyzer
// \author Sasha Bylinkin, alexander.bylinkin@gmail.com
// \since  April 2023

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "iostream"
#include "PWGUD/DataModel/UDTables.h"
//#include "Common/DataModel/PIDResponse.h"
//#include "PWGUD/Core/RLhelper.h"
#include <TString.h>
#include "TLorentzVector.h"
using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
#define mpion 0.1396
#define mkaon 0.4937
#define mproton 0.9383
struct SGSpectraAnalyzer {
  HistogramRegistry registry{
    "registry",
    {// Pion histograms for each eta bin and gapSide
     {"hPtPion_gap0", "Pion pT in eta [-1, 1] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_gap1", "Pion pT in eta [-1, 1] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_gap2", "Pion pT in eta [-1, 1] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hEtaPion_gap0", "Pion eta Gap 0; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hEtaPion_gap1", "Pion eta Gap 1; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hEtaPion_gap2", "Pion eta Gap 2; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hPtPion_etaBin1_gap0", "Pion pT in eta [-1, -0.5] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin2_gap0", "Pion pT in eta [-0.5, 0] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin3_gap0", "Pion pT in eta [0, 0.5] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin4_gap0", "Pion pT in eta [0.5, 1] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin1_gap1", "Pion pT in eta [-1, -0.5] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin2_gap1", "Pion pT in eta [-0.5, 0] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin3_gap1", "Pion pT in eta [0, 0.5] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin4_gap1", "Pion pT in eta [0.5, 1] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin1_gap2", "Pion pT in eta [-1, -0.5] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin2_gap2", "Pion pT in eta [-0.5, 0] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin3_gap2", "Pion pT in eta [0, 0.5] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin4_gap2", "Pion pT in eta [0.5, 1] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},

     // Kaon histograms for each eta bin and gapSide
     {"hPtKaon_gap0", "Kaon pT in eta [-1, 1] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_gap1", "Kaon pT in eta [-1, 1] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_gap2", "Kaon pT in eta [-1, 1] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hEtaKaon_gap0", "Kaon eta Gap 0; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hEtaKaon_gap1", "Kaon eta Gap 1; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hEtaKaon_gap2", "Kaon eta Gap 2; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hPtKaon_etaBin1_gap0", "Kaon pT in eta [-1, -0.5] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin2_gap0", "Kaon pT in eta [-0.5, 0] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin3_gap0", "Kaon pT in eta [0, 0.5] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin4_gap0", "Kaon pT in eta [0.5, 1] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin1_gap1", "Kaon pT in eta [-1, -0.5] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin2_gap1", "Kaon pT in eta [-0.5, 0] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin3_gap1", "Kaon pT in eta [0, 0.5] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin4_gap1", "Kaon pT in eta [0.5, 1] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin1_gap2", "Kaon pT in eta [-1, -0.5] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin2_gap2", "Kaon pT in eta [-0.5, 0] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin3_gap2", "Kaon pT in eta [0, 0.5] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin4_gap2", "Kaon pT in eta [0.5, 1] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},

     // Proton histograms for each eta bin and gapSide
     {"hPtProton_gap0", "Proton pT in eta [-1, 1] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_gap1", "Proton pT in eta [-1, 1] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_gap2", "Proton pT in eta [-1, 1] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hEtaProton_gap0", "Proton eta Gap 0; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hEtaProton_gap1", "Proton eta Gap 1; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hEtaProton_gap2", "Proton eta Gap 2; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hPtProton_etaBin1_gap0", "Proton pT in eta [-1, -0.5] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin2_gap0", "Proton pT in eta [-0.5, 0] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin3_gap0", "Proton pT in eta [0, 0.5] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin4_gap0", "Proton pT in eta [0.5, 1] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin1_gap1", "Proton pT in eta [-1, -0.5] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin2_gap1", "Proton pT in eta [-0.5, 0] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin3_gap1", "Proton pT in eta [0, 0.5] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin4_gap1", "Proton pT in eta [0.5, 1] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin1_gap2", "Proton pT in eta [-1, -0.5] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin2_gap2", "Proton pT in eta [-0.5, 0] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin3_gap2", "Proton pT in eta [0, 0.5] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin4_gap2", "Proton pT in eta [0.5, 1] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}}}};
  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  // using udtracks = soa::Join<aod::UDTracks,aod::UDTracksExtra,aod::UDTracksPID,aod::UDTracksDCA>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags>;
  // using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  // using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcs>;
  using UDCollisionFull = UDCollisionsFull::iterator;

  void process(UDCollisionFull const& collision, udtracksfull const& tracks)
  {
    TLorentzVector a;
    int gapSide = collision.gapSide();
    if (gapSide < 0 || gapSide > 2)
      return;

    for (auto& track : tracks) {
      if (!track.isPVContributor())
        break;
      float nSigmaPi = track.tpcNSigmaPi();
      float nSigmaK = track.tpcNSigmaKa();
      float nSigmaP = track.tpcNSigmaPr();

      if (std::abs(nSigmaPi) < 3.0) {
        a.SetXYZM(track.px(), track.py(), track.pz(), mpion);
        fillHistograms("Pion", a.Pt(), a.Eta(), gapSide);
      }
      if (std::abs(nSigmaK) < 3.0) {
        a.SetXYZM(track.px(), track.py(), track.pz(), mkaon);
        fillHistograms("Kaon", a.Pt(), a.Eta(), gapSide);
      }
      if (std::abs(nSigmaP) < 3.0) {
        a.SetXYZM(track.px(), track.py(), track.pz(), mproton);
        fillHistograms("Proton", a.Pt(), a.Eta(), gapSide);
      }
    }
  }

  void fillHistograms(const std::string& particleType, float pt, float eta, int gapSide)
  {
    if (particleType == "Pion") {
      if (gapSide == 0) {
        registry.get<TH1>(HIST("hPtPion_gap0"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaPion_gap0"))->Fill(eta);
      } else if (gapSide == 1) {
        registry.get<TH1>(HIST("hPtPion_gap1"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaPion_gap1"))->Fill(eta);
      } else if (gapSide == 2) {
        registry.get<TH1>(HIST("hPtPion_gap2"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaPion_gap2"))->Fill(eta);
      }
    } else if (particleType == "Kaon") {
      if (gapSide == 0) {
        registry.get<TH1>(HIST("hPtKaon_gap0"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaKaon_gap0"))->Fill(eta);
      } else if (gapSide == 1) {
        registry.get<TH1>(HIST("hPtKaon_gap1"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaKaon_gap1"))->Fill(eta);
      } else if (gapSide == 2) {
        registry.get<TH1>(HIST("hPtKaon_gap2"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaKaon_gap2"))->Fill(eta);
      }
    } else if (particleType == "Proton") {
      if (gapSide == 0) {
        registry.get<TH1>(HIST("hPtProton_gap0"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaProton_gap0"))->Fill(eta);
      } else if (gapSide == 1) {
        registry.get<TH1>(HIST("hPtProton_gap1"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaProton_gap1"))->Fill(eta);
      } else if (gapSide == 2) {
        registry.get<TH1>(HIST("hPtProton_gap2"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaProton_gap2"))->Fill(eta);
      }
    }
    std::vector<std::pair<float, float>> etaBins = {{-1, -0.5}, {-0.5, 0}, {0, 0.5}, {0.5, 1}};
    for (int i = 0; i < 4; ++i) {
      if (eta > etaBins[i].first && eta <= etaBins[i].second) {
        if (i == 0) {
          if (particleType == "Pion") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtPion_etaBin1_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtPion_etaBin1_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtPion_etaBin1_gap2"))->Fill(pt);
            }
          } else if (particleType == "Kaon") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtKaon_etaBin1_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtKaon_etaBin1_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtKaon_etaBin1_gap2"))->Fill(pt);
            }
          } else if (particleType == "Proton") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtProton_etaBin1_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtProton_etaBin1_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtProton_etaBin1_gap2"))->Fill(pt);
            }
          }
        }
        if (i == 1) {
          if (particleType == "Pion") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtPion_etaBin2_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtPion_etaBin2_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtPion_etaBin2_gap2"))->Fill(pt);
            }
          } else if (particleType == "Kaon") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtKaon_etaBin2_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtKaon_etaBin2_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtKaon_etaBin2_gap2"))->Fill(pt);
            }
          } else if (particleType == "Proton") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtProton_etaBin2_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtProton_etaBin2_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtProton_etaBin2_gap2"))->Fill(pt);
            }
          }
        }
        if (i == 2) {
          if (particleType == "Pion") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtPion_etaBin3_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtPion_etaBin3_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtPion_etaBin3_gap2"))->Fill(pt);
            }
          } else if (particleType == "Kaon") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtKaon_etaBin3_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtKaon_etaBin3_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtKaon_etaBin3_gap2"))->Fill(pt);
            }
          } else if (particleType == "Proton") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtProton_etaBin3_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtProton_etaBin3_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtProton_etaBin3_gap2"))->Fill(pt);
            }
          }
        }
        if (i == 3) {
          if (particleType == "Pion") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtPion_etaBin4_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtPion_etaBin4_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtPion_etaBin4_gap2"))->Fill(pt);
            }
          } else if (particleType == "Kaon") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtKaon_etaBin4_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtKaon_etaBin4_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtKaon_etaBin4_gap2"))->Fill(pt);
            }
          } else if (particleType == "Proton") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtProton_etaBin4_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtProton_etaBin4_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtProton_etaBin4_gap2"))->Fill(pt);
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGSpectraAnalyzer>(cfgc)};
}
