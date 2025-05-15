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
/// \file hyperhelium4sigmaAnalysis.cxx
/// \brief Simple check for injected hyper-helium4sigma (H4S) in MC productions
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>

#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGLF/DataModel/LFKinkDecayTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using CollisionsFull = soa::Join<aod::Collisions, aod::EvSel>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPr, aod::pidTPCFullAl, aod::pidTPCFullTr, aod::pidTPCFullPi>;
using MCLabeledTracksIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

//-------------------------------Check the decay channel of H4S-------------------------------
enum Channel {
  k2body = 0, // helium4, pion0
  k3body_p,   // triton, proton, pion0
  k3body_n,   // triton, neutron, pion+
  kNDecayChannel
};

template <class TMCTrackTo, typename TMCParticle>
Channel getDecayChannelH4S(TMCParticle const& particle)
{
  if (std::abs(particle.pdgCode()) != o2::constants::physics::Pdg::kHyperHelium4Sigma) {
    return kNDecayChannel;
  }
  bool haveAlpha = false, haveTriton = false, haveProton = false, haveNeuteron = false;
  bool haveAntiAlpha = false, haveAntiTriton = false, haveAntiProton = false, haveAntiNeuteron = false;
  bool havePionPlus = false, havePionMinus = false, havePion0 = false;
  auto daughters = particle.template daughters_as<TMCTrackTo>();
  for (const auto& mcDaughter : particle.template daughters_as<TMCTrackTo>()) {
    if (mcDaughter.pdgCode() == o2::constants::physics::Pdg::kAlpha) {
      haveAlpha = true;
    }
    if (mcDaughter.pdgCode() == -o2::constants::physics::Pdg::kAlpha) {
      haveAntiAlpha = true;
    }
    if (mcDaughter.pdgCode() == o2::constants::physics::Pdg::kTriton) {
      haveTriton = true;
    }
    if (mcDaughter.pdgCode() == -o2::constants::physics::Pdg::kTriton) {
      haveAntiTriton = true;
    }
    if (mcDaughter.pdgCode() == PDG_t::kProton) {
      haveProton = true;
    }
    if (mcDaughter.pdgCode() == -PDG_t::kProton) {
      haveAntiProton = true;
    }
    if (mcDaughter.pdgCode() == PDG_t::kNeutron) {
      haveNeuteron = true;
    }
    if (mcDaughter.pdgCode() == -PDG_t::kNeutron) {
      haveAntiNeuteron = true;
    }
    if (mcDaughter.pdgCode() == PDG_t::kPiPlus) {
      havePionPlus = true;
    }
    if (mcDaughter.pdgCode() == -PDG_t::kPiPlus) {
      havePionMinus = true;
    }
    if (mcDaughter.pdgCode() == PDG_t::kPi0) {
      havePion0 = true;
    }
  }

  if ((haveAlpha && havePion0) || (haveAntiAlpha && havePion0)) {
    return k2body;
  } else if ((haveTriton && haveProton && havePion0) || (haveAntiTriton && haveAntiProton && havePion0)) {
    return k3body_p;
  } else if ((haveTriton && haveNeuteron && havePionPlus) || (haveAntiTriton && haveAntiNeuteron && havePionMinus)) {
    return k3body_n;
  }

  return kNDecayChannel;
}
//--------------------------------------------------------------
struct Hyperhelium4sigmaAnalysis {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry registry{"registry", {}};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cutNSigmaAl{"cutNSigmaAl", 5, "NSigmaTPCAlpha"};

  void init(InitContext const&)
  {
    // Axes
    const AxisSpec vertexZAxis{100, -15., 15., "vrtx_{Z} [cm]"};
    const AxisSpec ptAxis{50, -10, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec nSigmaAxis{120, -6.f, 6.f, "n#sigma_{#alpha}"};
    const AxisSpec massAxis{100, 3.85, 4.25, "m (GeV/#it{c}^{2})"};

    registry.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    registry.add("h2MassHyperhelium4sigmaPt", "h2MassHyperhelium4sigmaPt", {HistType::kTH2F, {ptAxis, massAxis}});
    registry.add("h2NSigmaAlPt", "h2NSigmaAlPt", {HistType::kTH2F, {ptAxis, nSigmaAxis}});
  }

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               aod::KinkCands const& KinkCands, FullTracksExtIU const&)
  {
    if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
      return;
    }
    registry.fill(HIST("hVertexZRec"), collision.posZ());
    for (const auto& kinkCand : KinkCands) {
      auto dauTrack = kinkCand.trackDaug_as<FullTracksExtIU>();
      if (std::abs(dauTrack.tpcNSigmaAl()) > cutNSigmaAl) {
        continue;
      }
      float invMass = RecoDecay::m(std::array{std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()}, std::array{kinkCand.pxDaugNeut(), kinkCand.pyDaugNeut(), kinkCand.pzDaugNeut()}}, std::array{o2::constants::physics::MassAlpha, o2::constants::physics::MassPi0});
      registry.fill(HIST("h2MassHyperhelium4sigmaPt"), kinkCand.mothSign() * kinkCand.ptMoth(), invMass);
      registry.fill(HIST("h2NSigmaAlPt"), kinkCand.mothSign() * kinkCand.ptDaug(), dauTrack.tpcNSigmaAl());
    }
  }
};

//--------------------------------------------------------------
// check the performance of mcparticle
struct Hyperhelium4sigmaQa {
  // Basic checks
  HistogramRegistry registry{"registry", {}};

  ConfigurableAxis ptBins{"ptBins", {200, 0.f, 10.f}, "Binning for #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis ctBins{"ctBins", {100, 0.f, 25.f}, "Binning for c#it{t} (cm)"};
  ConfigurableAxis rigidityBins{"rigidityBins", {200, -10.f, 10.f}, "Binning for #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis nsigmaBins{"nsigmaBins", {120, -6.f, 6.f}, "Binning for n sigma"};
  ConfigurableAxis invMassBins{"invMassBins", {100, 3.85, 4.15f}, "Binning for invariant mass (GeV/#it{c}^{2})"};

  void init(InitContext&)
  {
    if (doprocessMC == true) {
      const AxisSpec ptAxis{ptBins, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec ctAxis{ctBins, "c#it{t} (cm)"};
      const AxisSpec rigidityAxis{rigidityBins, "p/z (GeV/#it{c})"};
      const AxisSpec nsigmaAxis{nsigmaBins, "TPC n#sigma"};
      const AxisSpec invMassAxis{invMassBins, "Inv Mass (GeV/#it{c}^{2})"};

      auto hCollCounter = registry.add<TH1>("hCollCounter", "hCollCounter", HistType::kTH1F, {{2, 0.0f, 2.0f}});
      registry.get<TH1>(HIST("hCollCounter"))->GetXaxis()->SetBinLabel(1, "Reconstructed Collisions");
      registry.get<TH1>(HIST("hCollCounter"))->GetXaxis()->SetBinLabel(2, "Selected");
      auto hMcCollCounter = registry.add<TH1>("hMcCollCounter", "hMcCollCounter", HistType::kTH1F, {{2, 0.0f, 2.0f}});
      registry.get<TH1>(HIST("hMcCollCounter"))->GetXaxis()->SetBinLabel(1, "MC Collisions");
      registry.get<TH1>(HIST("hMcCollCounter"))->GetXaxis()->SetBinLabel(2, "Reconstructed");

      auto hGenHyperHelium4SigmaCounter = registry.add<TH1>("hGenHyperHelium4SigmaCounter", "", HistType::kTH1F, {{10, 0.f, 10.f}});
      registry.get<TH1>(HIST("hGenHyperHelium4SigmaCounter"))->GetXaxis()->SetBinLabel(1, "H4S All");
      registry.get<TH1>(HIST("hGenHyperHelium4SigmaCounter"))->GetXaxis()->SetBinLabel(2, "Matter");
      registry.get<TH1>(HIST("hGenHyperHelium4SigmaCounter"))->GetXaxis()->SetBinLabel(3, "AntiMatter");
      registry.get<TH1>(HIST("hGenHyperHelium4SigmaCounter"))->GetXaxis()->SetBinLabel(4, "#alpha + #pi^{0}");
      registry.get<TH1>(HIST("hGenHyperHelium4SigmaCounter"))->GetXaxis()->SetBinLabel(5, "#bar{#alpha} + #pi^{0}");
      registry.get<TH1>(HIST("hGenHyperHelium4SigmaCounter"))->GetXaxis()->SetBinLabel(6, "t + p + #pi^{0}");
      registry.get<TH1>(HIST("hGenHyperHelium4SigmaCounter"))->GetXaxis()->SetBinLabel(7, "#bar{t} + #bar{p} + #pi^{0}");
      registry.get<TH1>(HIST("hGenHyperHelium4SigmaCounter"))->GetXaxis()->SetBinLabel(8, "t + n + #pi^{+}");
      registry.get<TH1>(HIST("hGenHyperHelium4SigmaCounter"))->GetXaxis()->SetBinLabel(9, "#bar{t} + #bar{n} + #pi^{+}");
      registry.get<TH1>(HIST("hGenHyperHelium4SigmaCounter"))->GetXaxis()->SetBinLabel(10, "Unexpected");

      auto hEvtSelectedHyperHelium4SigmaCounter = registry.add<TH1>("hEvtSelectedHyperHelium4SigmaCounter", "", HistType::kTH1F, {{2, 0.f, 2.f}});
      registry.get<TH1>(HIST("hEvtSelectedHyperHelium4SigmaCounter"))->GetXaxis()->SetBinLabel(1, "Generated");
      registry.get<TH1>(HIST("hEvtSelectedHyperHelium4SigmaCounter"))->GetXaxis()->SetBinLabel(2, "Survived");

      registry.add<TH1>("hGenHyperHelium4SigmaP", "", HistType::kTH1F, {ptAxis});
      registry.add<TH1>("hGenHyperHelium4SigmaPt", "", HistType::kTH1F, {ptAxis});
      registry.add<TH1>("hGenHyperHelium4SigmaCt", "", HistType::kTH1F, {ctAxis});
      registry.add<TH1>("hMcRecoInvMass", "", HistType::kTH1F, {invMassAxis});

      registry.add<TH2>("hDauHelium4TPCNSigma", "", HistType::kTH2F, {rigidityAxis, nsigmaAxis});
      registry.add<TH2>("hDauTritonTPCNSigma", "", HistType::kTH2F, {rigidityAxis, nsigmaAxis});
      registry.add<TH2>("hDauProtonTPCNSigma", "", HistType::kTH2F, {rigidityAxis, nsigmaAxis});
    }
  }

  // Configurable<bool> eventSel8Cut{"eventSel8Cut", true, "flag to enable event sel8 selection"};
  Configurable<bool> mcEventCut{"mcEventCut", true, "flag to enable mc event selection: kIsTriggerTVX and kNoTimeFrameBorder"};
  Configurable<bool> eventPosZCut{"eventPosZCut", true, "flag to enable event posZ selection"};
  Configurable<float> maxPosZ{"maxPosZ", 10.f, "max pv posZ for event selection"};

  Preslice<aod::McParticles> permcCollision = o2::aod::mcparticle::mcCollisionId;

  std::vector<int64_t> mcPartIndices;
  template <typename TTrackTable>
  void setTrackIDForMC(aod::McParticles const& particlesMC, TTrackTable const& tracks)
  {
    mcPartIndices.clear();
    mcPartIndices.resize(particlesMC.size());
    std::fill(mcPartIndices.begin(), mcPartIndices.end(), -1);
    for (const auto& track : tracks) {
      if (track.has_mcParticle()) {
        auto mcparticle = track.template mcParticle_as<aod::McParticles>();
        if (mcPartIndices[mcparticle.globalIndex()] == -1) {
          mcPartIndices[mcparticle.globalIndex()] = track.globalIndex();
        } else {
          auto candTrack = tracks.rawIteratorAt(mcPartIndices[mcparticle.globalIndex()]);
          // Use the track which has innest information (also best quality?
          if (track.x() < candTrack.x()) {
            mcPartIndices[mcparticle.globalIndex()] = track.globalIndex();
          }
        }
      }
    }
  }

  void processData(o2::aod::Collisions const&)
  {
    // dummy process function;
  }
  PROCESS_SWITCH(Hyperhelium4sigmaQa, processData, "process data", true);

  void processMC(aod::McCollisions const& mcCollisions, aod::McParticles const& particlesMC, o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels> const& collisions, MCLabeledTracksIU const& tracks)
  {
    setTrackIDForMC(particlesMC, tracks);
    std::vector<int64_t> selectedEvents(collisions.size());
    int nevts = 0;
    for (const auto& collision : collisions) {
      registry.fill(HIST("hCollCounter"), 0.5);
      if (mcEventCut && (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        continue;
      }
      if (eventPosZCut && std::abs(collision.posZ()) > maxPosZ) { // 10cm
        continue;
      }
      registry.fill(HIST("hCollCounter"), 1.5);
      selectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    selectedEvents.resize(nevts);

    for (const auto& mcCollision : mcCollisions) {
      registry.fill(HIST("hMcCollCounter"), 0.5);
      const auto evtReconstructedAndSelected = std::find(selectedEvents.begin(), selectedEvents.end(), mcCollision.globalIndex()) != selectedEvents.end();
      if (evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
        registry.fill(HIST("hMcCollCounter"), 1.5);
      } else {
        // continue;
      }

      const auto& dparticlesMC = particlesMC.sliceBy(permcCollision, mcCollision.globalIndex());

      for (const auto& mcparticle : dparticlesMC) {

        bool isMatter;
        if (mcparticle.pdgCode() == o2::constants::physics::Pdg::kHyperHelium4Sigma) {
          registry.fill(HIST("hGenHyperHelium4SigmaCounter"), 1.5);
          isMatter = true;
        } else if (mcparticle.pdgCode() == -o2::constants::physics::Pdg::kHyperHelium4Sigma) {
          registry.fill(HIST("hGenHyperHelium4SigmaCounter"), 2.5);
          isMatter = false;
        } else {
          continue;
        }

        registry.fill(HIST("hGenHyperHelium4SigmaCounter"), 0.5);
        registry.fill(HIST("hEvtSelectedHyperHelium4SigmaCounter"), 0.5);
        if (evtReconstructedAndSelected) {
          registry.fill(HIST("hEvtSelectedHyperHelium4SigmaCounter"), 1.5);
        }

        double svPos[3] = {-999, -999, -999};
        double dauHelium4Mom[3] = {-999, -999, -999};
        double dauTritonMom[3] = {-999, -999, -999};
        double dauProtonMom[3] = {-999, -999, -999};
        double dauNeuteronMom[3] = {-999, -999, -999};
        double dauChargedPionMom[3] = {-999, -999, -999};
        double dauPion0Mom[3] = {-999, -999, -999};
        auto dChannel = getDecayChannelH4S<aod::McParticles>(mcparticle);
        if (dChannel == kNDecayChannel) {
          registry.fill(HIST("hGenHyperHelium4SigmaCounter"), 9.5);
          continue;
        }
        for (const auto& mcparticleDaughter : mcparticle.daughters_as<aod::McParticles>()) {
          if (std::abs(mcparticleDaughter.pdgCode()) == o2::constants::physics::Pdg::kAlpha) {
            dauHelium4Mom[0] = mcparticleDaughter.px();
            dauHelium4Mom[1] = mcparticleDaughter.py();
            dauHelium4Mom[2] = mcparticleDaughter.pz();

            // get SV position for 2body decay
            svPos[0] = mcparticleDaughter.vx();
            svPos[1] = mcparticleDaughter.vy();
            svPos[2] = mcparticleDaughter.vz();

            if (mcPartIndices[mcparticleDaughter.globalIndex()] != -1) {
              auto track = tracks.rawIteratorAt(mcPartIndices[mcparticleDaughter.globalIndex()]);
              registry.fill(HIST("hDauHelium4TPCNSigma"), track.p() * track.sign(), track.tpcNSigmaAl());
            }
          } else if (std::abs(mcparticleDaughter.pdgCode()) == o2::constants::physics::Pdg::kTriton) {
            dauTritonMom[0] = mcparticleDaughter.px();
            dauTritonMom[1] = mcparticleDaughter.py();
            dauTritonMom[2] = mcparticleDaughter.pz();

            // get SV position for 3body decay
            svPos[0] = mcparticleDaughter.vx();
            svPos[1] = mcparticleDaughter.vy();
            svPos[2] = mcparticleDaughter.vz();

            if (mcPartIndices[mcparticleDaughter.globalIndex()] != -1) {
              auto track = tracks.rawIteratorAt(mcPartIndices[mcparticleDaughter.globalIndex()]);
              registry.fill(HIST("hDauTritonTPCNSigma"), track.p() * track.sign(), track.tpcNSigmaTr());
            }
          } else if (std::abs(mcparticleDaughter.pdgCode()) == PDG_t::kProton) {
            dauProtonMom[0] = mcparticleDaughter.px();
            dauProtonMom[1] = mcparticleDaughter.py();
            dauProtonMom[2] = mcparticleDaughter.pz();

            if (mcPartIndices[mcparticleDaughter.globalIndex()] != -1) {
              auto track = tracks.rawIteratorAt(mcPartIndices[mcparticleDaughter.globalIndex()]);
              registry.fill(HIST("hDauProtonTPCNSigma"), track.p() * track.sign(), track.tpcNSigmaPr());
            }
          } else if (std::abs(mcparticleDaughter.pdgCode()) == PDG_t::kNeutron) {
            dauNeuteronMom[0] = mcparticleDaughter.px();
            dauNeuteronMom[1] = mcparticleDaughter.py();
            dauNeuteronMom[2] = mcparticleDaughter.pz();
          } else if (std::abs(mcparticleDaughter.pdgCode()) == PDG_t::kPiPlus) {
            dauChargedPionMom[0] = mcparticleDaughter.px();
            dauChargedPionMom[1] = mcparticleDaughter.py();
            dauChargedPionMom[2] = mcparticleDaughter.pz();
          } else if (mcparticleDaughter.pdgCode() == PDG_t::kPi0) {
            dauPion0Mom[0] = mcparticleDaughter.px();
            dauPion0Mom[1] = mcparticleDaughter.py();
            dauPion0Mom[2] = mcparticleDaughter.pz();
          }
        }

        registry.fill(HIST("hGenHyperHelium4SigmaP"), mcparticle.p());
        registry.fill(HIST("hGenHyperHelium4SigmaPt"), mcparticle.pt());
        double ct = RecoDecay::sqrtSumOfSquares(svPos[0] - mcparticle.vx(), svPos[1] - mcparticle.vy(), svPos[2] - mcparticle.vz()) * o2::constants::physics::MassHyperHelium4Sigma / mcparticle.p();
        registry.fill(HIST("hGenHyperHelium4SigmaCt"), ct);

        if (dChannel == k2body) {
          if (isMatter) {
            registry.fill(HIST("hGenHyperHelium4SigmaCounter"), 3.5);
          } else {
            registry.fill(HIST("hGenHyperHelium4SigmaCounter"), 4.5);
          }
          double hyperHelium4SigmaMCMass = RecoDecay::m(std::array{std::array{dauHelium4Mom[0], dauHelium4Mom[1], dauHelium4Mom[2]}, std::array{dauPion0Mom[0], dauPion0Mom[1], dauPion0Mom[2]}}, std::array{o2::constants::physics::MassAlpha, o2::constants::physics::MassPi0});
          registry.fill(HIST("hMcRecoInvMass"), hyperHelium4SigmaMCMass);
        } else if (dChannel == k3body_p) {
          if (isMatter) {
            registry.fill(HIST("hGenHyperHelium4SigmaCounter"), 5.5);
          } else {
            registry.fill(HIST("hGenHyperHelium4SigmaCounter"), 6.5);
          }
          double hyperHelium4SigmaMCMass = RecoDecay::m(std::array{std::array{dauTritonMom[0], dauTritonMom[1], dauTritonMom[2]}, std::array{dauProtonMom[0], dauProtonMom[1], dauProtonMom[2]}, std::array{dauPion0Mom[0], dauPion0Mom[1], dauPion0Mom[2]}}, std::array{o2::constants::physics::MassTriton, o2::constants::physics::MassProton, o2::constants::physics::MassPi0});
          registry.fill(HIST("hMcRecoInvMass"), hyperHelium4SigmaMCMass);
        } else if (dChannel == k3body_n) {
          if (isMatter) {
            registry.fill(HIST("hGenHyperHelium4SigmaCounter"), 7.5);
          } else {
            registry.fill(HIST("hGenHyperHelium4SigmaCounter"), 8.5);
          }
          double hyperHelium4SigmaMCMass = RecoDecay::m(std::array{std::array{dauTritonMom[0], dauTritonMom[1], dauTritonMom[2]}, std::array{dauNeuteronMom[0], dauNeuteronMom[1], dauNeuteronMom[2]}, std::array{dauChargedPionMom[0], dauChargedPionMom[1], dauChargedPionMom[2]}}, std::array{o2::constants::physics::MassTriton, o2::constants::physics::MassNeutron, o2::constants::physics::MassPionCharged});
          registry.fill(HIST("hMcRecoInvMass"), hyperHelium4SigmaMCMass);
        }
      }
    }
  }
  PROCESS_SWITCH(Hyperhelium4sigmaQa, processMC, "do QA for MC prodcutions", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Hyperhelium4sigmaAnalysis>(cfgc),
    adaptAnalysisTask<Hyperhelium4sigmaQa>(cfgc),
  };
}
