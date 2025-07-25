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

/// \file CFFilterAll.cxx
/// \brief Selection of events with triplets and pairs for femtoscopic studies
///
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "../filterTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>
#include <TMath.h>

#include "fairlogger/Logger.h"

#include <iostream>
#include <iterator>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
using FemtoFullCollision =
  soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator;

using FemtoFullTracks =
  soa::Join<aod::FullTracks, aod::TracksDCA,
            aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
            aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe>;
} // namespace o2::aod

struct CFFillterPPPhi {

  // Table for storing filter decisions
  // Leave commented for now
  /*Produces<aod::CFFilters> tags;*/

  /*event selection*/
  Configurable<bool> ConfEvtSelectZvtx{
    "ConfEvtSelectZvtx",
    true,
    "Event selection includes max. z-Vertex"};

  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};

  Configurable<bool> ConfEvtOfflineCheck{
    "ConfEvtOfflineCheck",
    false,
    "Evt sel: check for offline selection"};

  Configurable<float> ConfResoInvMassLowLimit{"ConfResoInvMassLowLimit", 1.011461, "Lower limit of the Reso invariant mass"};
  Configurable<float> ConfResoInvMassUpLimit{"ConfResoInvMassUpLimit", 1.027461, "Upper limit of the Reso invariant mass"};

  Configurable<float> Q3Max{"Q3Max", 1.f, "Max Q3"};

  /*track selection*/

  Configurable<float> ConfTrkEtaPr{"ConfTrkEtaPr", 0.85, "Et protona"};                                                             // 0.8
  Configurable<float> ConfTrkDCAxyPr{"ConfTrkDCAxyPr", 0.15, "DCAxy proton"};                                                       // 0.1
  Configurable<float> ConfTrkDCAzPr{"ConfTrkDCAzPr", 0.3, "DCAz proton"};                                                           // 0.2
  Configurable<float> ConfNClusPr{"ConfNClusPr", 70, "NClusters proton"};                                                           // 0.2
  Configurable<float> ConfNCrossedPr{"ConfNCrossedPr", 65, "NCrossedRows proton"};                                                  // 0.2
  Configurable<float> ConfTrkTPCfClsPr{"ConfTrkTPCfClsPr", 0.83, "Minimum fraction of crossed rows over findable clusters proton"}; // 0.2

  Configurable<float> ConfTrkPtPrUp{"ConfTrkPtPrUp", 6.0, "Pt_up proton"};            // 2.0
  Configurable<float> ConfTrkPtPrDown{"ConfTrkPtPrDown", 0.35, "Pt_down proton"};     // 0.5
  Configurable<float> ConfTrkPTPCPrThr{"ConfTrkPTPCPrThr", 0.8, "p_TPC,Thr proton"};  // 0.75
  Configurable<float> ConfTrkPrSigmaPID{"ConfTrkPrSigmaPID", 3.50, "n_sigma proton"}; // 3.0

  Configurable<float> ConfTrkEtaKa{"ConfTrkEtaKa", 0.85, "Eta kaon"};                                                             // 0.8
  Configurable<float> ConfTrkDCAxyKa{"ConfTrkDCAxyKa", 0.15, "DCAxy kaon"};                                                       // 0.1
  Configurable<float> ConfTrkDCAzKa{"ConfTrkDCAzKa", 0.3, "DCAz kaon"};                                                           // 0.2
  Configurable<float> ConfNClusKa{"ConfNClusKa", 70, "NClusters kaon"};                                                           // 0.2
  Configurable<float> ConfNCrossedKa{"ConfNCrossedKa", 65, "NCrossedRows kaon"};                                                  // 0.2
  Configurable<float> ConfTrkTPCfClsKa{"ConfTrkTPCfClsKa", 0.80, "Minimum fraction of crossed rows over findable clusters kaon"}; // 0.2

  Configurable<float> ConfTrkPtKaUp{"ConfTrkPtKaUp", 6.0, "Pt_up kaon"};            // 2.0
  Configurable<float> ConfTrkPtKaDown{"ConfTrkPtKaDown", 0.05, "Pt_down kaon"};     // 0.15
  Configurable<float> ConfTrkPTPCKaThr{"ConfTrkPTPCKaThr", 0.40, "p_TPC,Thr kaon"}; // 0.4
  Configurable<float> ConfTrkKaSigmaPID{"ConfTrkKaSigmaPID", 3.50, "n_sigma kaon"}; // 3.0

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {

    // histograms
    registry.add("fProcessedEvents", "CF - event filtered;;Events", HistType::kTH1F, {{3, -0.5, 2.5}});
    std::vector<std::string> eventTitles = {"all", "rejected", "ppphi"};
    for (size_t iBin = 0; iBin < eventTitles.size(); iBin++) {
      registry.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
    }

    // event cuts
    registry.add("EventCuts/fMultiplicityBefore", "Multiplicity of all processed events;Mult;Entries", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("EventCuts/fMultiplicityAfter", "Multiplicity after event cuts;Mult;Entries", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("EventCuts/fZvtxBefore", "Zvtx of all processed events;Z_{vtx};Entries", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("EventCuts/fZvtxAfter", "Zvtx after event cuts;Z_{vtx};Entries", HistType::kTH1F, {{1000, -15, 15}});

    registry.add("TrackCuts/fPtBefore", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fEtaBefore", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/fPhiBefore", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    // proton cuts
    registry.add("TrackCuts/Proton/fPProton", "Momentum of protons at PV;p (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Proton/fPTPCProton", "Momentum of protons at PV;p (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Proton/fPtProton", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Proton/fMomCorProtonDif", "Momentum correlation;p_{reco} (GeV/c); (p_{TPC} - p_{reco}) (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {600, -3, 3}}});
    registry.add("TrackCuts/Proton/fMomCorProtonRatio", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {{500, 0, 10}, {200, -1, 1}}});
    registry.add("TrackCuts/Proton/fEtaProton", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/Proton/fPhiProton", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    registry.add("TrackCuts/Proton/fDCAxyProton", "fDCAxy Proton;DCA_{XY};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/Proton/fDCAzProton", "fDCAz Proton;DCA_{Z};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/Proton/fNsigmaTPCvsPProton", "NSigmaTPC Proton;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Proton/fNsigmaTOFvsPProton", "NSigmaTOF Proton;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Proton/fNsigmaTPCTOFvsPProton", "NSigmaTPCTOF Proton;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    registry.add("TrackCuts/Proton/fTPCsClsProton", "fTPCsCls Proton;TPC Shared Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/Proton/fTPCcRowsProton", "fTPCcRows Proton;TPC Crossed Rows;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/Proton/fTrkTPCfClsProton", "fTrkTPCfCls Proton;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    registry.add("TrackCuts/Proton/fTPCnclsProton", "fTPCncls Proton;TPC Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});

    registry.add("TrackCuts/AntiProton/fPAntiProton", "Momentum of antiprotons at PV;p (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiProton/fPTPCAntiProton", "Momentum of protons at PV;p (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiProton/fPtAntiProton", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiProton/fMomCorAntiProtonDif", "Momentum correlation;p_{reco} (GeV/c); (p_{TPC} - p_{reco}) (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {600, -3, 3}}});
    registry.add("TrackCuts/AntiProton/fMomCorAntiProtonRatio", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {{500, 0, 10}, {200, -1, 1}}});
    registry.add("TrackCuts/AntiProton/fEtaAntiProton", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/AntiProton/fPhiAntiProton", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    registry.add("TrackCuts/AntiProton/fDCAxyAntiProton", "fDCAxy AntiProton;DCA_{XY};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/AntiProton/fDCAzAntiProton", "fDCAz AntiProton;DCA_{Z};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/AntiProton/fNsigmaTPCvsPAntiProton", "NSigmaTPC AntiProton;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiProton/fNsigmaTOFvsPAntiProton", "NSigmaTOF AntiProton;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiProton/fNsigmaTPCTOFvsPAntiProton", "NSigmaTPCTOF AntiProton;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    registry.add("TrackCuts/AntiProton/fTPCsClsAntiProton", "fTPCsCls AntiProton;TPC Shared Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/AntiProton/fTPCcRowsAntiProton", "fTPCcRows AntiProton;TPC Crossed Rows;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/AntiProton/fTrkTPCfClsAntiProton", "fTrkTPCfCls AntiProton;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    registry.add("TrackCuts/AntiProton/fTPCnclsAntiProton", "fTPCncls AntiProton;TPC Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});

    // kaon cuts
    registry.add("TrackCuts/Kaon/fPKaon", "Momentum of kaons at PV;p (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Kaon/fPTPCKaon", "Momentum of protons at PV;p (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Kaon/fPtKaon", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Kaon/fMomCorKaonDif", "Momentum correlation;p_{reco} (GeV/c); (p_{TPC} - p_{reco}) (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {600, -3, 3}}});
    registry.add("TrackCuts/Kaon/fMomCorKaonRatio", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {{500, 0, 10}, {200, -1, 1}}});
    registry.add("TrackCuts/Kaon/fEtaKaon", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/Kaon/fPhiKaon", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    registry.add("TrackCuts/Kaon/fDCAxyKaon", "fDCAxy Kaon;DCA_{XY};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/Kaon/fDCAzKaon", "fDCAz Kaon;DCA_{Z};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/Kaon/fNsigmaTPCvsPKaon", "NSigmaTPC Kaon;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Kaon/fNsigmaTOFvsPKaon", "NSigmaTOF Kaon;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Kaon/fNsigmaTPCTOFvsPKaon", "NSigmaTPCTOF Kaon;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    registry.add("TrackCuts/Kaon/fTPCsClsKaon", "fTPCsCls Kaon;TPC Shared Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/Kaon/fTPCcRowsKaon", "fTPCcRows Kaon;TPC Crossed Rows;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/Kaon/fTrkTPCfClsKaon", "fTrkTPCfCls Kaon;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    registry.add("TrackCuts/Kaon/fTPCnclsKaon", "fTPCncls Kaon;TPC Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});

    registry.add("TrackCuts/AntiKaon/fPAntiKaon", "Momentum of antikaons at PV;p (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiKaon/fPTPCAntiKaon", "Momentum of protons at PV;p (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiKaon/fPtAntiKaon", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiKaon/fMomCorAntiKaonDif", "Momentum correlation;p_{reco} (GeV/c); (p_{TPC} - p_{reco}) (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {600, -3, 3}}});
    registry.add("TrackCuts/AntiKaon/fMomCorAntiKaonRatio", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {{500, 0, 10}, {200, -1, 1}}});
    registry.add("TrackCuts/AntiKaon/fEtaAntiKaon", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/AntiKaon/fPhiAntiKaon", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    registry.add("TrackCuts/AntiKaon/fDCAxyAntiKaon", "fDCAxy AntiKaon;DCA_{XY};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/AntiKaon/fDCAzAntiKaon", "fDCAz AntiKaon;DCA_{Z};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/AntiKaon/fNsigmaTPCvsPAntiKaon", "NSigmaTPC AntiKaon;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiKaon/fNsigmaTOFvsPAntiKaon", "NSigmaTOF AntiKaon;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiKaon/fNsigmaTPCTOFvsPAntiKaon", "NSigmaTPCTOF AntiKaon;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    registry.add("TrackCuts/AntiKaon/fTPCsClsAntiKaon", "fTPCsCls AntiKaon;TPC Shared Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/AntiKaon/fTPCcRowsAntiKaon", "fTPCcRows AntiKaon;TPC Crossed Rows;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/AntiKaon/fTrkTPCfClsAntiKaon", "fTrkTPCfCls AntiKaon;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    registry.add("TrackCuts/AntiKaon/fTPCnclsAntiKaon", "fTPCncls AntiKaon;TPC Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});

    // phi cuts
    registry.add("TrackCuts/Phi/fPtPhiBefore", "Transverse momentum V0s;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Phi/fInvMassPhiBefore", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.8, 1.5}});

    registry.add("TrackCuts/Phi/fEtaPhiBefore", "Pseudorapidity of V0;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/Phi/fPhiPhiBefore", "Azimuthal angle of V0;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    registry.add("TrackCuts/Phi/fPtPhi", "Transverse momentum V0s;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Phi/fInvMassPhi", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.8, 1.5}});
    registry.add("TrackCuts/Phi/fEtaPhi", "Pseudorapidity of V0;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/Phi/fPhiPhi", "Azimuthal angle of V0;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    // phi daughter
    registry.add("TrackCuts/Phi/PosDaughter/Pt", "Transverse momentum Pos Daugh tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Phi/PosDaughter/Eta", "Phi Pos Daugh Eta;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/Phi/PosDaughter/Phi", "Azimuthal angle of Pos Daugh tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    registry.add("TrackCuts/Phi/NegDaughter/Pt", "Transverse momentum Neg Daugh tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Phi/NegDaughter/Eta", "Phi Neg Daugh Eta;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/Phi/NegDaughter/Phi", "Azimuthal angle of Neg Daugh tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    // triggers
    registry.add("ppphi/fMultiplicity", "Multiplicity of all triggered events;Mult;Entries", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("ppphi/fZvtx", "Zvtx of all triggered events;Z_{vtx};Entries", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("ppphi/fSE_particle", "Same Event distribution;SE;Q_{3} (GeV/c)", HistType::kTH1F, {{8000, 0, 8}});
    registry.add("ppphi/fSE_antiparticle", "Same Event distribution;SE;Q_{3} (GeV/c)", HistType::kTH1F, {{8000, 0, 8}});
    registry.add("ppphi/fProtonPtVsQ3", "Same Event distribution;SE;Q_{3} (GeV/c)", HistType::kTH1F, {{8000, 0, 8}});
    registry.add("ppphi/fPhiPtVsQ3", "Same Event distribution;SE;Q_{3} (GeV/c)", HistType::kTH1F, {{8000, 0, 8}});
    registry.add("ppphi/fAntiProtonPtVsQ3", "Same Event distribution;SE;Q_{3} (GeV/c)", HistType::kTH1F, {{8000, 0, 8}});
    registry.add("ppphi/fAntiPhiPtVsQ3", "Same Event distribution;SE;Q_{3} (GeV/c)", HistType::kTH1F, {{8000, 0, 8}});
  }

  template <typename T>
  bool isSelectedEvent(T const& col)
  {
    if (ConfEvtSelectZvtx && std::abs(col.posZ()) > ConfEvtZvtx) {
      return false;
    }
    if (ConfEvtOfflineCheck && !col.sel8()) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedTrackProton(T const& track)
  {
    bool isSelected = false;
    if (track.pt() <= ConfTrkPtPrUp.value && track.pt() >= ConfTrkPtPrDown.value && std::abs(track.eta()) <= ConfTrkEtaPr.value && std::abs(track.dcaXY()) <= ConfTrkDCAxyPr.value && std::abs(track.dcaZ()) <= ConfTrkDCAzPr.value && track.tpcNClsCrossedRows() >= ConfNCrossedPr.value && track.tpcNClsFound() >= ConfNClusPr.value && track.tpcCrossedRowsOverFindableCls() >= ConfTrkTPCfClsPr.value) {
      if (track.tpcInnerParam() < ConfTrkPTPCPrThr.value && std::abs(track.tpcNSigmaPr()) <= ConfTrkPrSigmaPID.value) {
        isSelected = true;
      }
      if (track.tpcInnerParam() >= ConfTrkPTPCPrThr.value && std::abs(std::sqrt(track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr())) <= ConfTrkPrSigmaPID.value) {
        isSelected = true;
      }
    }
    return isSelected;
  }

  template <typename T>
  bool isSelectedTrackKaon(T const& track)
  {
    bool isSelected = false;
    if (track.pt() <= ConfTrkPtKaUp.value && track.pt() >= ConfTrkPtKaDown.value && std::abs(track.eta()) <= ConfTrkEtaKa.value && std::abs(track.dcaXY()) <= ConfTrkDCAxyKa.value && std::abs(track.dcaZ()) <= ConfTrkDCAzKa.value && track.tpcNClsCrossedRows() >= ConfNCrossedKa.value && track.tpcNClsFound() >= ConfNClusKa.value && track.tpcCrossedRowsOverFindableCls() >= ConfTrkTPCfClsKa.value) {
      if (track.tpcInnerParam() < ConfTrkPTPCKaThr.value && std::abs(track.tpcNSigmaKa()) <= ConfTrkKaSigmaPID.value) {
        isSelected = true;
      }
      if (track.tpcInnerParam() >= ConfTrkPTPCKaThr.value && std::abs(std::sqrt(track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa())) <= ConfTrkKaSigmaPID.value) {
        isSelected = true;
      }
    }
    return isSelected;
  }

  float mMassProton = o2::constants::physics::MassProton;
  float mMassKaonPlus = o2::constants::physics::MassKPlus;
  float mMassKaonMinus = o2::constants::physics::MassKMinus;

  float mMassPhi = o2::constants::physics::MassPhi;

  float getkstar(const ROOT::Math::PtEtaPhiMVector part1,
                 const ROOT::Math::PtEtaPhiMVector part2)
  {
    const ROOT::Math::PtEtaPhiMVector trackSum = part1 + part2;
    const float beta = trackSum.Beta();
    const float betax =
      beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay =
      beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());
    ROOT::Math::PxPyPzMVector PartOneCMS(part1);
    ROOT::Math::PxPyPzMVector PartTwoCMS(part2);
    const ROOT::Math::Boost boostPRF =
      ROOT::Math::Boost(-betax, -betay, -betaz);
    PartOneCMS = boostPRF(PartOneCMS);
    PartTwoCMS = boostPRF(PartTwoCMS);
    const ROOT::Math::PxPyPzMVector trackRelK = PartOneCMS - PartTwoCMS;
    return 0.5 * trackRelK.P();
  }

  ROOT::Math::PxPyPzEVector getqij(const ROOT::Math::PtEtaPhiMVector parti,
                                   const ROOT::Math::PtEtaPhiMVector partj)
  {
    ROOT::Math::PxPyPzEVector vecparti(parti);
    ROOT::Math::PxPyPzEVector vecpartj(partj);
    ROOT::Math::PxPyPzEVector trackSum = vecparti + vecpartj;
    ROOT::Math::PxPyPzEVector trackDifference = vecparti - vecpartj;
    float scaling = trackDifference.Dot(trackSum) / trackSum.Dot(trackSum);
    return trackDifference - scaling * trackSum;
  }
  float getQ3(const ROOT::Math::PtEtaPhiMVector part1,
              const ROOT::Math::PtEtaPhiMVector part2,
              const ROOT::Math::PtEtaPhiMVector part3)
  {
    ROOT::Math::PxPyPzEVector q12 = getqij(part1, part2);
    ROOT::Math::PxPyPzEVector q23 = getqij(part2, part3);
    ROOT::Math::PxPyPzEVector q31 = getqij(part3, part1);
    float Q32 = q12.M2() + q23.M2() + q31.M2();
    return sqrt(-Q32);
  }

  void process(aod::FemtoFullCollision const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks)
  {
    registry.fill(HIST("fProcessedEvents"), 0);
    registry.fill(HIST("EventCuts/fMultiplicityBefore"), col.multNTracksPV());
    registry.fill(HIST("EventCuts/fZvtxBefore"), col.posZ());

    int lowQ3Triplets = 0;

    if (isSelectedEvent(col)) {

      registry.fill(HIST("EventCuts/fMultiplicityAfter"), col.multNTracksPV());
      registry.fill(HIST("EventCuts/fZvtxAfter"), col.posZ());

      std::vector<ROOT::Math::PtEtaPhiMVector> protons, antiprotons, kaons, antikaons, phi;

      // keep track of proton indices
      std::vector<int> ProtonIndex = {};
      std::vector<int> AntiProtonIndex = {};
      std::vector<int> KaonIndex = {};
      std::vector<int> AntiKaonIndex = {};

      for (auto& track : tracks) {
        registry.fill(HIST("TrackCuts/fPtBefore"), track.pt());
        registry.fill(HIST("TrackCuts/fEtaBefore"), track.eta());
        registry.fill(HIST("TrackCuts/fPhiBefore"), track.phi());

        if (isSelectedTrackProton(track)) {
          ROOT::Math::PtEtaPhiMVector temp(track.pt(), track.eta(), track.phi(), mMassProton);
          if (track.sign() > 0) {
            protons.push_back(temp);

            registry.fill(HIST("TrackCuts/Proton/fPProton"), track.p());
            registry.fill(HIST("TrackCuts/Proton/fPTPCProton"), track.tpcInnerParam());
            registry.fill(HIST("TrackCuts/Proton/fPtProton"), track.pt());
            registry.fill(HIST("TrackCuts/Proton/fMomCorProtonDif"), track.p(), track.tpcInnerParam() - track.p());
            registry.fill(HIST("TrackCuts/Proton/fMomCorProtonRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
            registry.fill(HIST("TrackCuts/Proton/fEtaProton"), track.eta());
            registry.fill(HIST("TrackCuts/Proton/fPhiProton"), track.phi());
            registry.fill(HIST("TrackCuts/Proton/fNsigmaTPCvsPProton"), track.tpcInnerParam(), track.tpcNSigmaPr());
            registry.fill(HIST("TrackCuts/Proton/fNsigmaTOFvsPProton"), track.tpcInnerParam(), track.tofNSigmaPr());
            registry.fill(HIST("TrackCuts/Proton/fNsigmaTPCTOFvsPProton"), track.tpcInnerParam(), std::abs(std::sqrt(track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr())));

            registry.fill(HIST("TrackCuts/Proton/fDCAxyProton"), track.dcaXY());
            registry.fill(HIST("TrackCuts/Proton/fDCAzProton"), track.dcaZ());
            registry.fill(HIST("TrackCuts/Proton/fTPCsClsProton"), track.tpcNClsShared());
            registry.fill(HIST("TrackCuts/Proton/fTPCcRowsProton"), track.tpcNClsCrossedRows());
            registry.fill(HIST("TrackCuts/Proton/fTrkTPCfClsProton"), track.tpcCrossedRowsOverFindableCls());
            registry.fill(HIST("TrackCuts/Proton/fTPCnclsProton"), track.tpcNClsFound());
            // ProtonIndex.push_back(track.globalIndex());
          }
          if (track.sign() < 0) {
            antiprotons.push_back(temp);

            registry.fill(HIST("TrackCuts/AntiProton/fPAntiProton"), track.p());
            registry.fill(HIST("TrackCuts/AntiProton/fPTPCAntiProton"), track.tpcInnerParam());
            registry.fill(HIST("TrackCuts/AntiProton/fPtAntiProton"), track.pt());
            registry.fill(HIST("TrackCuts/AntiProton/fMomCorAntiProtonDif"), track.p(), track.tpcInnerParam() - track.p());
            registry.fill(HIST("TrackCuts/AntiProton/fMomCorAntiProtonRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
            registry.fill(HIST("TrackCuts/AntiProton/fEtaAntiProton"), track.eta());
            registry.fill(HIST("TrackCuts/AntiProton/fPhiAntiProton"), track.phi());
            registry.fill(HIST("TrackCuts/AntiProton/fNsigmaTPCvsPAntiProton"), track.tpcInnerParam(), track.tpcNSigmaPr());
            registry.fill(HIST("TrackCuts/AntiProton/fNsigmaTOFvsPAntiProton"), track.tpcInnerParam(), track.tofNSigmaPr());
            registry.fill(HIST("TrackCuts/AntiProton/fNsigmaTPCTOFvsPAntiProton"), track.tpcInnerParam(), std::abs(std::sqrt(track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr())));

            registry.fill(HIST("TrackCuts/AntiProton/fDCAxyAntiProton"), track.dcaXY());
            registry.fill(HIST("TrackCuts/AntiProton/fDCAzAntiProton"), track.dcaZ());
            registry.fill(HIST("TrackCuts/AntiProton/fTPCsClsAntiProton"), track.tpcNClsShared());
            registry.fill(HIST("TrackCuts/AntiProton/fTPCcRowsAntiProton"), track.tpcNClsCrossedRows());
            registry.fill(HIST("TrackCuts/AntiProton/fTrkTPCfClsAntiProton"), track.tpcCrossedRowsOverFindableCls());
            registry.fill(HIST("TrackCuts/AntiProton/fTPCnclsAntiProton"), track.tpcNClsFound());
            // AntiProtonIndex.push_back(track.globalIndex());
          }
        }

        if (isSelectedTrackKaon(track)) {
          ROOT::Math::PtEtaPhiMVector temp(track.pt(), track.eta(), track.phi(), mMassKaonPlus);
          if (track.sign() > 0) {
            temp.SetM(mMassKaonPlus);
            kaons.push_back(temp);
            registry.fill(HIST("TrackCuts/Kaon/fPKaon"), track.p());
            registry.fill(HIST("TrackCuts/Kaon/fPTPCKaon"), track.tpcInnerParam());
            registry.fill(HIST("TrackCuts/Kaon/fPtKaon"), track.pt());
            registry.fill(HIST("TrackCuts/Kaon/fMomCorKaonDif"), track.p(), track.tpcInnerParam() - track.p());
            registry.fill(HIST("TrackCuts/Kaon/fMomCorKaonRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
            registry.fill(HIST("TrackCuts/Kaon/fEtaKaon"), track.eta());
            registry.fill(HIST("TrackCuts/Kaon/fPhiKaon"), track.phi());
            registry.fill(HIST("TrackCuts/Kaon/fNsigmaTPCvsPKaon"), track.tpcInnerParam(), track.tpcNSigmaKa());
            registry.fill(HIST("TrackCuts/Kaon/fNsigmaTOFvsPKaon"), track.tpcInnerParam(), track.tofNSigmaKa());
            registry.fill(HIST("TrackCuts/Kaon/fNsigmaTPCTOFvsPKaon"), track.tpcInnerParam(), std::abs(std::sqrt(track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa())));

            registry.fill(HIST("TrackCuts/Kaon/fDCAxyKaon"), track.dcaXY());
            registry.fill(HIST("TrackCuts/Kaon/fDCAzKaon"), track.dcaZ());

            registry.fill(HIST("TrackCuts/Kaon/fTPCsClsKaon"), track.tpcNClsShared());
            registry.fill(HIST("TrackCuts/Kaon/fTPCcRowsKaon"), track.tpcNClsCrossedRows());
            registry.fill(HIST("TrackCuts/Kaon/fTrkTPCfClsKaon"), track.tpcCrossedRowsOverFindableCls());
            registry.fill(HIST("TrackCuts/Kaon/fTPCnclsKaon"), track.tpcNClsFound());
            // KaonIndex.push_back(track.globalIndex());
          }
          if (track.sign() < 0) {
            temp.SetM(mMassKaonMinus);

            antikaons.push_back(temp);
            registry.fill(HIST("TrackCuts/AntiKaon/fPAntiKaon"), track.p());
            registry.fill(HIST("TrackCuts/AntiKaon/fPTPCAntiKaon"), track.tpcInnerParam());
            registry.fill(HIST("TrackCuts/AntiKaon/fPtAntiKaon"), track.pt());
            registry.fill(HIST("TrackCuts/AntiKaon/fMomCorAntiKaonDif"), track.p(), track.tpcInnerParam() - track.p());
            registry.fill(HIST("TrackCuts/AntiKaon/fMomCorAntiKaonRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
            registry.fill(HIST("TrackCuts/AntiKaon/fEtaAntiKaon"), track.eta());
            registry.fill(HIST("TrackCuts/AntiKaon/fPhiAntiKaon"), track.phi());
            registry.fill(HIST("TrackCuts/AntiKaon/fNsigmaTPCvsPAntiKaon"), track.tpcInnerParam(), track.tpcNSigmaKa());
            registry.fill(HIST("TrackCuts/AntiKaon/fNsigmaTOFvsPAntiKaon"), track.tpcInnerParam(), track.tofNSigmaKa());
            registry.fill(HIST("TrackCuts/AntiKaon/fNsigmaTPCTOFvsPAntiKaon"), track.tpcInnerParam(), std::abs(std::sqrt(track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa())));

            registry.fill(HIST("TrackCuts/AntiKaon/fDCAxyAntiKaon"), track.dcaXY());
            registry.fill(HIST("TrackCuts/AntiKaon/fDCAzAntiKaon"), track.dcaZ());
            registry.fill(HIST("TrackCuts/AntiKaon/fTPCsClsAntiKaon"), track.tpcNClsShared());
            registry.fill(HIST("TrackCuts/AntiKaon/fTPCcRowsAntiKaon"), track.tpcNClsCrossedRows());
            registry.fill(HIST("TrackCuts/AntiKaon/fTrkTPCfClsAntiKaon"), track.tpcCrossedRowsOverFindableCls());
            registry.fill(HIST("TrackCuts/AntiKaon/fTPCnclsAntiKaon"), track.tpcNClsFound());
            // AntiKaonIndex.push_back(track.globalIndex());
          }
        }

        // end track
      }

      for (const auto& postrack : kaons) {
        for (const auto& negtrack : antikaons) {

          ROOT::Math::PtEtaPhiMVector temp = postrack + negtrack;
          // temp.SetM(mMassPhi);
          registry.fill(HIST("TrackCuts/Phi/fInvMassPhiBefore"), temp.M());

          registry.fill(HIST("TrackCuts/Phi/fPtPhiBefore"), temp.pt());
          registry.fill(HIST("TrackCuts/Phi/fEtaPhiBefore"), temp.eta());
          registry.fill(HIST("TrackCuts/Phi/fPhiPhiBefore"), temp.phi());

          if ((temp.M() >= ConfResoInvMassLowLimit.value) && (temp.M() <= ConfResoInvMassUpLimit.value)) {
            // ROOT::Math::PtEtaPhiMVector temp = postrack + negtrack;
            phi.push_back(temp);

            registry.fill(HIST("TrackCuts/Phi/fPtPhi"), temp.pt());
            registry.fill(HIST("TrackCuts/Phi/fEtaPhi"), temp.eta());
            registry.fill(HIST("TrackCuts/Phi/fPhiPhi"), temp.phi());
            registry.fill(HIST("TrackCuts/Phi/fInvMassPhi"), temp.M());

            registry.fill(HIST("TrackCuts/Phi/PosDaughter/Pt"), postrack.pt());
            registry.fill(HIST("TrackCuts/Phi/PosDaughter/Eta"), postrack.eta());
            registry.fill(HIST("TrackCuts/Phi/PosDaughter/Phi"), postrack.phi());

            registry.fill(HIST("TrackCuts/Phi/NegDaughter/Pt"), negtrack.pt());
            registry.fill(HIST("TrackCuts/Phi/NegDaughter/Eta"), negtrack.eta());
            registry.fill(HIST("TrackCuts/Phi/NegDaughter/Phi"), negtrack.phi());
          }
        }
      }

      // ppphi trigger
      float Q3 = 999.f;

      for (auto iProton1 = protons.begin(); iProton1 != protons.end(); ++iProton1) {
        auto iProton2 = iProton1 + 1;
        // auto i1 = std::distance(protons.begin(), iProton1);
        for (; iProton2 != protons.end(); ++iProton2) {
          // auto i2 = std::distance(protons.begin(), iProton2);
          for (auto iPhi1 = phi.begin(); iPhi1 != phi.end(); ++iPhi1) {
            // auto i3 = std::distance(phi.begin(), iPhi1);

            Q3 = getQ3(*iProton1, *iProton2, *iPhi1);
            registry.fill(HIST("ppphi/fSE_particle"), Q3);
            registry.fill(HIST("ppphi/fProtonPtVsQ3"), Q3, (*iProton1).Pt());
            registry.fill(HIST("ppphi/fProtonPtVsQ3"), Q3, (*iProton2).Pt());
            registry.fill(HIST("ppphi/fPhiPtVsQ3"), Q3, (*iPhi1).Pt());
            if (Q3 < Q3Max.value) {
              lowQ3Triplets += 1;
            }
          }
        }
      }

      // apapphi trigger

      for (auto iAntiProton1 = antiprotons.begin(); iAntiProton1 != antiprotons.end(); ++iAntiProton1) {
        auto iAntiProton2 = iAntiProton1 + 1;
        // auto i1 = std::distance(antiprotons.begin(), iAntiProton1);
        for (; iAntiProton2 != antiprotons.end(); ++iAntiProton2) {
          // auto i2 = std::distance(antiprotons.begin(), iAntiProton2);
          for (auto iPhi1 = phi.begin(); iPhi1 != phi.end(); ++iPhi1) {
            // auto i3 = std::distance(phi.begin(), iPhi1);

            Q3 = getQ3(*iAntiProton1, *iAntiProton2, *iPhi1);
            registry.fill(HIST("ppphi/fSE_antiparticle"), Q3);
            registry.fill(HIST("ppphi/fAntiProtonPtVsQ3"), Q3, (*iAntiProton1).Pt());
            registry.fill(HIST("ppphi/fAntiProtonPtVsQ3"), Q3, (*iAntiProton2).Pt());
            registry.fill(HIST("ppphi/fAntiPhiPtVsQ3"), Q3, (*iPhi1).Pt());
            if (Q3 < Q3Max.value) {
              lowQ3Triplets += 1;
            }
          }
        }
      }

      // end event
    }

    // create tags for three body triggers
    if (lowQ3Triplets > 0) {
      registry.fill(HIST("fProcessedEvents"), 2);
      registry.fill(HIST("ppphi/fMultiplicity"), col.multNTracksPV());
      registry.fill(HIST("ppphi/fZvtx"), col.posZ());
    } else {
      registry.fill(HIST("fProcessedEvents"), 1);
    }
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<CFFillterPPPhi>(cfg)};
}
