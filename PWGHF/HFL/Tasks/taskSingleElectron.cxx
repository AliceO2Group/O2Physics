// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file taskSingleElectron.cxx
/// \brief task for electrons from heavy-flavour hadron decays
/// \author Jonghan Park (Jeonbuk National University), Seul I Jeong (Pusan National University)

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <map>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum PdgCode { kEta = 221, kOmega = 223, kEtaPrime = 331 };

enum SourceType {
  NotElec = 0,      // not electron
  DirectCharm = 1,  // electrons from prompt charm hadrons
  DirectBeauty = 2, // electrons from primary beauty hadrons
  BeautyCharm = 3,  // electrons from non-prompt charm hadrons
  DirectGamma = 4,  // electrons from direct photon
  GammaPi0 = 5,
  GammaEta = 6,
  GammaOmega = 7,
  GammaPhi = 8,
  GammaEtaPrime = 9,
  GammaRho0 = 10,
  GammaK0s = 11,
  GammaK0l = 12,
  GammaKe3 = 13,
  GammaLambda0 = 14,
  GammaSigma = 15,
  Pi0 = 16,
  Eta = 17,
  Omega = 18,
  Phi = 19,
  EtaPrime = 20,
  Rho0 = 21,
  K0s = 22,
  K0l = 23,
  Ke3 = 24,
  Lambda0 = 25,
  Sigma = 26,
  Else = 27
};

struct HfTaskSingleElectron {

  // Produces

  // Configurable
  Configurable<int> nContribMin{"nContribMin", 2, "min number of contributors"};
  Configurable<float> posZMax{"posZMax", 10., "max posZ cut"};
  Configurable<float> ptTrackMax{"ptTrackMax", 10., "max pt cut"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.5, "min pt cut"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "eta cut"};
  Configurable<int> tpcNCrossedRowMin{"tpcNCrossedRowMin", 70,
                                      "max of TPC n cluster crossed rows"};
  Configurable<float> tpcNClsFoundOverFindableMin{
      "tpcNClsFoundOverFindableMin", 0.8,
      "min # of TPC found/findable clusters"};
  Configurable<float> tpcChi2perNClMax{"tpcChi2perNClMax", 4.,
                                       "min # of tpc chi2 per clusters"};
  Configurable<int> itsIBClsMin{"itsIBClsMin", 3,
                                "min # of its clusters in IB"};
  Configurable<float> itsChi2perNClMax{"itsChi2perNClMax", 6.,
                                       "min # of tpc chi2 per clusters"};
  Configurable<float> dcaxyMax{"dcaxyMax", 1., "max of track dca in xy"};
  Configurable<float> dcazMax{"dcazMax", 2., "max of track dca in z"};
  Configurable<float> tofNSigmaMax{"tofNSigmaMax", 3., "max of tof nsigma"};
  Configurable<float> tpcNSigmaMin{"tpcNSigmaMin", -1., "min of tpc nsigma"};
  Configurable<float> tpcNSigmaMax{"tpcNSigmaMax", 3., "max of tpc nsigma"};

  Configurable<int> nBinsP{"nBinsP", 1500,
                           "number of bins of particle momentum"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};

  Configurable<int> tpcNSigmaHadronMax{"tpcNSigmaHadronMax", -3,
                                       "max of tpc hadron nsigma"};
  Configurable<int> tpcNSigmaHadronMin{"tpcNSigmaHadronMin", -5,
                                       "min of tpc hadron nsigma"};

  // SliceCache
  SliceCache cache;

  // using declarations
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using TracksEl =
      soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension,
                aod::TracksExtra, aod::TracksDCA, aod::pidTOFFullEl,
                aod::pidTPCFullEl>;
  using McTracksEl =
      soa::Join<aod::Tracks, aod::TrackExtra, aod::TracksDCA, aod::pidTOFFullEl,
                aod::pidTPCFullEl, aod::McTrackLabels>;

  // Filter
  Filter collZFilter = nabs(aod::collision::posZ) < posZMax;

  // Partition

  // ConfigurableAxis
  ConfigurableAxis axisPtEl{
      "axisPtEl",
      {VARIABLE_WIDTH, 0.5f, 0.6f, 0.7f, 0.8f,  0.9f, 1.f,   1.1f,
       1.2f,           1.3f, 1.4f, 1.5f, 1.75f, 2.0f, 2.25f, 2.5f,
       2.75f,          3.f,  3.5f, 4.0f, 5.0f,  6.0f, 8.0f,  10.0f},
      "electron pt bins"};

  // Histogram registry
  HistogramRegistry histos{
      "histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const &) {
    // AxisSpec
    const AxisSpec axisEvt{4, 0., 4., "nEvents"};
    const AxisSpec axisNCont{100, 0., 100., "nCont"};
    const AxisSpec axisPosZ{600, -30., 30., "Z_{pos}"};
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisP{nBinsP, 0., 15., "p_{T}"};
    const AxisSpec axisPt{nBinsPt, 0., 15., "p_{T}"};
    const AxisSpec axisNsig{800, -20., 20.};
    const AxisSpec axisTrackIp{4000, -0.2, 0.2, "dca"};

    // create histograms
    histos.add("nEvents", "Number of events", kTH1D, {{1, 0., 1.}});
    histos.add("VtxZ", "VtxZ; cm; entries", kTH1D, {axisPosZ});
    histos.add("etaTrack", "etaTrack; #eta; entries", kTH1D, {axisEta});
    histos.add("ptTrack",
               "#it{p}_{T} distribution of selected tracks; #it{p}_{T} "
               "(GeV/#it{c}); entries",
               kTH1D, {axisPt});

    // QA plots for trigger track selection
    histos.add("tpcNClsTrack", "tpcNClsTrack", kTH1D, {{200, 0, 200}});
    histos.add("tpcFoundFindableTrack", "", kTH1D, {{10, 0, 1}});
    histos.add("tpcChi2Track", "", kTH1D, {{100, 0, 10}});
    histos.add("itsIBClsTrack", "", kTH1D, {{10, 0, 10}});
    histos.add("itsChi2Track", "", kTH1D, {{50, 0, 50}});
    histos.add("dcaXYTrack", "", kTH1D, {{600, -3, 3}});
    histos.add("dcaZTrack", "", kTH1D, {{600, -3, 3}});

    // pid
    histos.add("tofNSigPt", "", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("tofNSigPtQA", "", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("tpcNSigP", "", kTH2D, {{axisP}, {axisNsig}});
    histos.add("tpcNSigPt", "", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("tpcNSigPAfterTofCut", "", kTH2D, {{axisP}, {axisNsig}});
    histos.add("tpcNSigPtAfterTofCut", "", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("tpcNSigPtQA", "", kTH2D, {{axisPtEl}, {axisNsig}});

    // track impact parameter
    histos.add("dcaTrack", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("dcaBeauty", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("dcaCharm", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("dcaDalitz", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("dcaConv", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("dcaHadron", "", kTH2D, {{axisPtEl}, {axisTrackIp}});

    // QA plots for MC
    histos.add("hPdgC", "", kTH1D, {{10001, -0.5, 10000.5}});
    histos.add("hPdgB", "", kTH1D, {{10001, -0.5, 10000.5}});
    histos.add("hPdgDa", "", kTH1D, {{10001, -0.5, 10000.5}});
    histos.add("hPdgCo", "", kTH1D, {{10001, -0.5, 10000.5}});
  }

  template <typename TrackType> bool trackSel(const TrackType &track) {
    if ((track.pt() > ptTrackMax) || (track.pt() < ptTrackMin)) {
      return false;
    }
    if (std::abs(track.eta()) > etaTrackMax) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < tpcNCrossedRowMin) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < tpcNClsFoundOverFindableMin) {
      return false;
    }

    if (track.tpcChi2NCl() > tpcChi2perNClMax) {
      return false;
    }

    if (!(track.itsNClsInnerBarrel() == itsIBClsMin)) {
      return false;
    }

    if (track.itsChi2NCl() > itsChi2perNClMax) {
      return false;
    }

    if (std::abs(track.dcaXY()) > dcaxyMax) {
      return false;
    }

    if (std::abs(track.dcaZ()) > dcazMax) {
      return false;
    }

    return true;
  }

  template <typename TrackType>
  int getElecSource(const TrackType &track, double &mpt, int &mpdg) {
    auto mcpart = track.mcParticle();
    if (std::abs(mcpart.pdgCode()) != kElectron) {
      return NotElec;
    }

    int motherPdg = -999;
    int grmotherPdg = -999;
    int ggrmotherPdg = -999; // mother, grand mother, grand grand mother pdg
    int motherPt = -999.;
    int grmotherPt = -999;
    int ggrmotherPt = -999.; // mother, grand mother, grand grand mother pt

    auto partMother =
        mcpart.template mothers_as<aod::McParticles>(); // first mother particle
                                                        // of electron
    auto partMotherCopy = partMother; // copy of the first mother
    auto mctrack = partMother;        // will change all the time

    motherPt = partMother.front().pt();                 // first mother pt
    motherPdg = std::abs(partMother.front().pdgCode()); // first mother pdg
    mpt = motherPt;   // copy of first mother pt
    mpdg = motherPdg; // copy of first mother pdg

    // check if electron from charm hadrons
    if ((static_cast<int>(motherPdg / 100.) % 10) == kCharm ||
        (static_cast<int>(motherPdg / 1000.) % 10) == kCharm) {

      // iterate until B hadron is found as an ancestor
      while (partMother.size()) {
        mctrack = partMother.front().template mothers_as<aod::McParticles>();
        if (mctrack.size()) {
          auto const &grmothersIdsVec = mctrack.front().mothersIds();

          if (grmothersIdsVec.empty()) {
            return DirectCharm;
          }
          grmotherPt = mctrack.front().pt();
          grmotherPdg = std::abs(mctrack.front().pdgCode());
          if ((static_cast<int>(grmotherPdg / 100.) % 10) == kBottom ||
              (static_cast<int>(grmotherPdg / 1000.) % 10) == kBottom) {
            mpt = grmotherPt;
            mpdg = grmotherPdg;
            return BeautyCharm;
          }
        }
        partMother = mctrack;
      }
    } else if ((static_cast<int>(motherPdg / 100.) % 10) == kBottom ||
               (static_cast<int>(motherPdg / 1000.) % 10) ==
                   kBottom) { // check if electron from beauty hadrons
      return DirectBeauty;
    } else if (motherPdg ==
               kGamma) { // check if electron from photon conversion
      mctrack = partMother.front().template mothers_as<aod::McParticles>();
      if (mctrack.size()) {
        auto const &grmothersIdsVec = mctrack.front().mothersIds();
        if (grmothersIdsVec.empty()) {
          return DirectGamma;
        }
        grmotherPdg = std::abs(mctrack.front().pdgCode());
        mpdg = grmotherPdg;
        mpt = mctrack.front().pt();

        partMother = mctrack;
        mctrack = partMother.front().template mothers_as<aod::McParticles>();
        if (mctrack.size()) {
          auto const &ggrmothersIdsVec = mctrack.front().mothersIds();
          if (ggrmothersIdsVec.empty()) {
            if (grmotherPdg == kPi0) {
              return GammaPi0;
            }
            if (grmotherPdg == PdgCode::kEta) {
              return GammaEta;
            }
            if (grmotherPdg == PdgCode::kOmega) {
              return GammaOmega;
            }
            if (grmotherPdg == kPhi) {
              return GammaPhi;
            }
            if (grmotherPdg == PdgCode::kEtaPrime) {
              return GammaEtaPrime;
            }
            if (grmotherPdg == kRho770_0) {
              return GammaRho0;
            }
            return Else;
          }
          ggrmotherPdg = mctrack.front().pdgCode();
          ggrmotherPt = mctrack.front().pt();
          mpdg = ggrmotherPdg;
          mpt = ggrmotherPt;
          if (grmotherPdg == kPi0) {
            if (ggrmotherPdg == kK0Short) {
              return GammaK0s;
            }
            if (ggrmotherPdg == kK0Long) {
              return GammaK0l;
            }
            if (ggrmotherPdg == kKPlus) {
              return GammaKe3;
            }
            if (ggrmotherPdg == kLambda0) {
              return GammaLambda0;
            }
            if (ggrmotherPdg == kSigmaPlus) {
              return GammaSigma;
            }
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaPi0;
          }
          if (grmotherPdg == PdgCode::kEta) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaEta;
          }
          if (grmotherPdg == PdgCode::kOmega) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaOmega;
          }
          if (grmotherPdg == kPhi) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaPhi;
          }
          if (grmotherPdg == PdgCode::kEtaPrime) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaEtaPrime;
          }
          if (grmotherPdg == kRho770_0) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaRho0;
          }
          return Else;
        }
      }
    } else { // check if electron from Dalitz decays
      mctrack = partMother.front().template mothers_as<aod::McParticles>();
      if (mctrack.size()) {
        auto const &grmothersIdsVec = mctrack.front().mothersIds();
        if (grmothersIdsVec.empty()) {
          static const std::map<int, SourceType> pdgToSource = {
              {kPi0, Pi0},
              {PdgCode::kEta, Eta},
              {PdgCode::kOmega, Omega},
              {kPhi, Phi},
              {PdgCode::kEtaPrime, EtaPrime},
              {kRho770_0, Rho0},
              {kKPlus, Ke3},
              {kK0Long, K0l}};

          auto it = pdgToSource.find(motherPdg);
          if (it != pdgToSource.end()) {
            return it->second;
          }
          return Else;
        }
        if (motherPdg == kPi0) {
          grmotherPt = mctrack.front().pt();
          grmotherPdg = mctrack.front().pdgCode();
          mpt = grmotherPt;
          mpdg = grmotherPdg;
          if (grmotherPdg == kK0Short) {
            return K0s;
          }
          if (grmotherPdg == kK0Long) {
            return K0l;
          }
          if (grmotherPdg == kKPlus) {
            return Ke3;
          }
          if (grmotherPdg == kLambda0) {
            return Lambda0;
          }
          if (grmotherPdg == kSigmaPlus) {
            return Sigma;
          }
          mpt = motherPt;
          mpdg = motherPdg;
          return Pi0;
        }
        if (motherPdg == PdgCode::kEta) {
          return Eta;
        }
        if (motherPdg == PdgCode::kOmega) {
          return Omega;
        }
        if (motherPdg == kPhi) {
          return Phi;
        }
        if (motherPdg == PdgCode::kEtaPrime) {
          return EtaPrime;
        }
        if (motherPdg == kRho770_0) {
          return Rho0;
        }
        if (motherPdg == kKPlus) {
          return Ke3;
        }
        if (motherPdg == kK0Long) {
          return K0l;
        }
        return Else;
      }
    }

    return Else;
  }

  void processData(soa::Filtered<MyCollisions>::iterator const &collision,
                   TracksEl const &tracks) {
    float const flagAnalysedEvt = 0.5;

    if (!collision.sel8()) {
      return;
    }

    if (collision.numContrib() < nContribMin) {
      return;
    }

    histos.fill(HIST("VtxZ"), collision.posZ());
    histos.fill(HIST("nEvents"), flagAnalysedEvt);

    for (const auto &track : tracks) {

      if (!trackSel(track)) {
        continue;
      }

      if (!(track.passedITSRefit() && track.passedTPCRefit())) {
        continue;
      }

      histos.fill(HIST("etaTrack"), track.eta());
      histos.fill(HIST("ptTrack"), track.pt());

      histos.fill(HIST("tpcNClsTrack"), track.tpcNClsCrossedRows());
      histos.fill(HIST("tpcFoundFindableTrack"),
                  track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("tpcChi2Track"), track.tpcChi2NCl());
      histos.fill(HIST("itsIBClsTrack"), track.itsNClsInnerBarrel());
      histos.fill(HIST("itsChi2Track"), track.itsChi2NCl());
      histos.fill(HIST("dcaXYTrack"), track.dcaXY());
      histos.fill(HIST("dcaZTrack"), track.dcaZ());

      histos.fill(HIST("tofNSigPt"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("tpcNSigP"), track.p(), track.tpcNSigmaEl());
      histos.fill(HIST("tpcNSigPt"), track.pt(), track.tpcNSigmaEl());

      if (std::abs(track.tofNSigmaEl()) > tofNSigmaMax) {
        continue;
      }
      histos.fill(HIST("tofNSigPtQA"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("tpcNSigPAfterTofCut"), track.p(), track.tpcNSigmaEl());
      histos.fill(HIST("tpcNSigPtAfterTofCut"), track.pt(),
                  track.tpcNSigmaEl());

      if (track.tpcNSigmaEl() < tpcNSigmaMin ||
          track.tpcNSigmaEl() > tpcNSigmaMax) {
        continue;
      }

      if (track.tpcNSigmaEl() < tpcNSigmaHadronMax &&
          track.tpcNSigmaEl() > tpcNSigmaHadronMin) {

        histos.fill(HIST("dcaHadron"), track.pt(), track.dcaXY());
      }

      histos.fill(HIST("tpcNSigPtQA"), track.pt(), track.tpcNSigmaEl());

      histos.fill(HIST("dcaTrack"), track.pt(), track.dcaXY());
    }
  }
  PROCESS_SWITCH(HfTaskSingleElectron, processData, "For real data", true);

  void processMc(soa::Filtered<MyCollisions>::iterator const &collision,
                 McTracksEl const &tracks, aod::McParticles const &) {
    float const flagAnalysedEvt = 0.5;

    if (!collision.sel8()) {
      return;
    }

    if (collision.numContrib() < nContribMin) {
      return;
    }

    histos.fill(HIST("VtxZ"), collision.posZ());
    histos.fill(HIST("nEvents"), flagAnalysedEvt);

    for (const auto &track : tracks) {

      if (!trackSel(track)) {
        continue;
      }

      histos.fill(HIST("etaTrack"), track.eta());
      histos.fill(HIST("ptTrack"), track.pt());

      histos.fill(HIST("tpcNClsTrack"), track.tpcNClsCrossedRows());
      histos.fill(HIST("tpcFoundFindableTrack"),
                  track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("tpcChi2Track"), track.tpcChi2NCl());
      histos.fill(HIST("itsIBClsTrack"), track.itsNClsInnerBarrel());
      histos.fill(HIST("dcaXYTrack"), track.dcaXY());
      histos.fill(HIST("dcaZTrack"), track.dcaZ());

      histos.fill(HIST("tofNSigPt"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("tpcNSigPt"), track.pt(), track.tpcNSigmaEl());

      int mpdg;   // electron source pdg code
      double mpt; // electron source pt
      int const source = getElecSource(track, mpt, mpdg);

      if (source == DirectBeauty || source == BeautyCharm) {
        histos.fill(HIST("hPdgB"), mpdg);
        histos.fill(HIST("dcaBeauty"), track.pt(), track.dcaXY());
      }

      if (source == DirectCharm) {
        histos.fill(HIST("hPdgC"), mpdg);
        histos.fill(HIST("dcaCharm"), track.pt(), track.dcaXY());
      }

      if (source >= GammaPi0 && source <= GammaSigma) {
        histos.fill(HIST("hPdgCo"), mpdg);
        histos.fill(HIST("dcaConv"), track.pt(), track.dcaXY());
      }

      if (source >= Pi0 && source <= Sigma) {
        histos.fill(HIST("hPdgDa"), mpdg);
        histos.fill(HIST("dcaDalitz"), track.pt(), track.dcaXY());
      }

      if (track.tpcNSigmaEl() < tpcNSigmaHadronMax &&
          track.tpcNSigmaEl() > tpcNSigmaHadronMin)
        histos.fill(HIST("dcaHadron"), track.pt(), track.dcaXY());

      if (std::abs(track.tofNSigmaEl()) > tofNSigmaMax) {
        continue;
      }
      histos.fill(HIST("tofNSigPtQA"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("tpcNSigPtAfterTofCut"), track.pt(),
                  track.tpcNSigmaEl());

      if (track.tpcNSigmaEl() < tpcNSigmaMin ||
          track.tpcNSigmaEl() > tpcNSigmaMax) {
        continue;
      }
      histos.fill(HIST("tpcNSigPtQA"), track.pt(), track.tpcNSigmaEl());

      histos.fill(HIST("dcaTrack"), track.pt(), track.dcaXY());
    }
  }
  PROCESS_SWITCH(HfTaskSingleElectron, processMc, "For real data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const &cfgc) {
  return WorkflowSpec{adaptAnalysisTask<HfTaskSingleElectron>(cfgc)};
}
