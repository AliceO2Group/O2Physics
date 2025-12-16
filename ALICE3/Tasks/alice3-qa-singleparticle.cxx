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

///
/// \file    qa-singleparticle.cxx
/// \author  Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief   Task to monitor the single particle QA, at the particle and track level, showing the tracked and the origin of particles
///

#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TDatabasePDG.h>
#include <TMCProcess.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct Alice3SingleParticle {
  Service<o2::framework::O2DatabasePDG> pdg;
  Configurable<int> PDG{"PDG", 2212, "PDG code of the particle of interest"};
  Configurable<int> IsStable{"IsStable", 0, "Flag to check stable particles"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> ptBins{"pt-bins", 500, "Number of pT bins"};
  Configurable<float> ptMin{"pt-min", 0.f, "Lower limit in pT"};
  Configurable<float> ptMax{"pt-max", 10.f, "Upper limit in pT"};
  Configurable<int> etaBins{"eta-bins", 500, "Number of eta bins"};
  Configurable<float> etaMin{"eta-min", -3.f, "Lower limit in eta"};
  Configurable<float> etaMax{"eta-max", 3.f, "Upper limit in eta"};
  Configurable<float> yMin{"y-min", -3.f, "Lower limit in y"};
  Configurable<float> yMax{"y-max", 3.f, "Upper limit in y"};
  Configurable<int> prodBins{"prod-bins", 100, "Number of production vertex bins"};
  Configurable<float> prodMin{"prod-min", -1.f, "Lower limit in production vertex"};
  Configurable<float> prodMax{"prod-max", 1.f, "Upper limit in production vertex"};
  Configurable<int> prodBinsZ{"prod-bins-z", 100, "Number of production vertex bins"};
  Configurable<float> prodMinZ{"prod-min-z", -10.f, "Lower limit in production vertex along Z"};
  Configurable<float> prodMaxZ{"prod-max-z", 10.f, "Upper limit in production vertex along Z"};
  Configurable<float> charge{"charge", 1.f, "Particle charge to scale the reconstructed momentum"};
  Configurable<bool> doPrint{"doPrint", false, "Flag to print debug messages"};

  void init(InitContext&)
  {
    int nEnabled = 0;
    if (doprocessStandard) {
      nEnabled++;
    }
    if (doprocessParticleOnly) {
      nEnabled++;
    }
    if (doprocessNonIU) {
      nEnabled++;
    }
    if (nEnabled == 0 || nEnabled > 1) {
      LOG(fatal) << "You can't process with " << nEnabled << " process functions, pick one";
    }

    pdg->AddParticle("deuteron", "deuteron", 1.8756134, kTRUE, 0.0, 3, "Nucleus", 1000010020);
    pdg->AddAntiParticle("anti-deuteron", -1000010020);

    pdg->AddParticle("triton", "triton", 2.8089218, kTRUE, 0.0, 3, "Nucleus", 1000010030);
    pdg->AddAntiParticle("anti-triton", -1000010030);

    pdg->AddParticle("helium3", "helium3", 2.80839160743, kTRUE, 0.0, 6, "Nucleus", 1000020030);
    pdg->AddAntiParticle("anti-helium3", -1000020030);

    pdg->AddParticle("helium4", "helium4", 2.80839160743, kTRUE, 0.0, 6, "Nucleus", 1000020040);
    pdg->AddAntiParticle("anti-helium4", -1000020040);

    const TString tit = Form("%i", PDG.value);
    const AxisSpec axisVx{100, -1, 1, "Vtx_{x}"};
    const AxisSpec axisVy{100, -1, 1, "Vtx_{y}"};
    const AxisSpec axisVz{100, -20, 20, "Vtx_{z}"};
    const AxisSpec axisPDGs{100, 0.f, 100.f, "PDG code"};
    const AxisSpec axisCharge{21, -10.5f, 10.5f, "Charge (e/3)"};
    const AxisSpec axisP{ptBins, ptMax, ptMax, "#it{p} (GeV/#it{c})"};
    const AxisSpec axisPt{ptBins, ptMin, ptMax, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtMC{ptBins, ptMin, ptMax, "#it{p}_{T}^{MC} (GeV/#it{c})"};
    const AxisSpec axisPx{ptBins, -ptMax, ptMax, "#it{p}_{x} (GeV/#it{c})"};
    const AxisSpec axisPy{ptBins, -ptMax, ptMax, "#it{p}_{y} (GeV/#it{c})"};
    const AxisSpec axisPz{ptBins, -ptMax, ptMax, "#it{p}_{z} (GeV/#it{c})"};
    const AxisSpec axisEta{etaBins, etaMin, etaMax, "#it{#eta}"};
    const AxisSpec axisY{etaBins, yMin, yMax, "#it{y}"};
    const AxisSpec axisE{etaBins, 0, 1000, "E"};
    const AxisSpec axisProdx{prodBins, prodMin, prodMax, "Prod. Vertex X (cm)"};
    const AxisSpec axisPrody{prodBins, prodMin, prodMax, "Prod. Vertex Y (cm)"};
    const AxisSpec axisProdz{prodBinsZ, prodMinZ, prodMaxZ, "Prod. Vertex Z (cm)"};
    const AxisSpec axisProdRadius{prodBins, 0., 2. * prodMax, "Prod. Vertex Radius (cm)"};

    histos.add("event/VtxX", "Vertex X", kTH1D, {axisVx});
    histos.add("event/VtxY", "Vertex Y", kTH1D, {axisVy});
    histos.add("event/VtxZ", "Vertex Z", kTH1D, {axisVz});

    histos.add("particle/PDGs", "Particle PDGs", kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/PDGsPrimaries", "Particle PDGs of Primaries", kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/PDGsSecondaries", "Particle PDGs of Secondaries", kTH2D, {axisPDGs, axisCharge});
    auto h = histos.add<TH1>("particle/producedProcess", "Particle process of production " + tit, kTH1D, {{kMaxMCProcess, -0.5, -0.5 + kMaxMCProcess, "getProcess"}});
    for (int i = 0; i < kMaxMCProcess; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, TMCProcessName[i]);
    }
    histos.add("particle/Pt", "Particle Pt " + tit, kTH1D, {axisPt});
    histos.add("particle/P", "Particle P " + tit, kTH1D, {axisP});
    histos.add("particle/primariesPt", "Particle Pt (primary) " + tit, kTH1D, {axisPt});
    histos.add("particle/primariesP", "Particle P (primary) " + tit, kTH1D, {axisP});
    histos.add("particle/secondariesPt", "Particle Pt (secondary) " + tit, kTH1D, {axisPt});
    histos.add("particle/secondariesP", "Particle P (secondary) " + tit, kTH1D, {axisP});
    histos.add("particle/prodVx", "Particle Prod. Vertex X " + tit, kTH1D, {axisProdx});
    histos.add("particle/prodVy", "Particle Prod. Vertex Y " + tit, kTH1D, {axisPrody});
    histos.add("particle/prodVz", "Particle Prod. Vertex Z " + tit, kTH1D, {axisProdz});
    histos.add("particle/prodRadius", "Particle Prod. Vertex Radius " + tit, kTH1D, {axisProdRadius});
    histos.add("particle/prodVxVy", "Particle Prod. Vertex X vs Vertex Y " + tit, kTH2D, {axisProdx, axisPrody});
    histos.add("particle/prodVxVyPrm", "Particle Prod. Vertex X vs Vertex Y Primary " + tit, kTH2D, {axisProdx, axisPrody});
    histos.add("particle/prodVxVyStr", "Particle Prod. Vertex X vs Vertex Y Decay " + tit, kTH2D, {axisProdx, axisPrody});
    histos.add("particle/prodVxVyMat", "Particle Prod. Vertex X vs Vertex Y Material " + tit, kTH2D, {axisProdx, axisPrody});
    histos.add("particle/prodVxVsPt", "Particle Prod. Vertex X " + tit, kTH2D, {axisPt, axisProdx});
    histos.add("particle/prodVyVsPt", "Particle Prod. Vertex Y " + tit, kTH2D, {axisPt, axisPrody});
    histos.add("particle/prodVzVsPt", "Particle Prod. Vertex Z " + tit, kTH2D, {axisPt, axisProdz});
    histos.add("particle/prodRadiusVsPt", "Particle Prod. Vertex Radius " + tit, kTH2D, {axisPt, axisProdRadius});
    histos.add("particle/prodRadiusVsEta", "Particle Prod. Vertex Radius " + tit, kTH2D, {axisEta, axisProdRadius});
    histos.add("particle/prodRadius3DVsPt", "Particle Prod. Vertex Radius XYZ " + tit, kTH2D, {axisPt, axisProdRadius});
    histos.add("particle/Eta", "Particle Eta " + tit, kTH1D, {axisEta});
    histos.add("particle/primariesEta", "Particle Eta (primary) " + tit, kTH1D, {axisEta});
    histos.add("particle/secondariesEta", "Particle Eta (secondary) " + tit, kTH1D, {axisEta});
    histos.add("particle/Y", "Particle Y " + tit, kTH1D, {axisY});
    histos.add("particle/primariesY", "Particle Y (primary)" + tit, kTH1D, {axisY});
    histos.add("particle/secondariesY", "Particle Y (secondary)" + tit, kTH1D, {axisY});
    histos.add("particle/EvsPz", "Particle E vs Pz " + tit, kTH2D, {axisE, axisPz});
    histos.add("particle/YvzPz", "Particle Y vs Pz " + tit, kTH2D, {axisY, axisPz});
    histos.add("particle/EtavzPz", "Particle Eta vs Pz " + tit, kTH2D, {axisEta, axisPz});
    histos.add("particle/PtvzPz", "Particle Pt vs Pz " + tit, kTH2D, {axisPt, axisPz});
    histos.add("particle/PvzPz", "Particle P vs Pz " + tit, kTH2D, {axisP, axisPz});
    histos.add("particle/Px", "Particle Px " + tit, kTH1D, {axisPx});
    histos.add("particle/Py", "Particle Py " + tit, kTH1D, {axisPy});
    histos.add("particle/Pz", "Particle Pz " + tit, kTH1D, {axisPz});

    histos.add("particle/daughters/Number", "Number of Daughters " + tit, kTH1D, {{20, -0.5, 19.5}});
    histos.add("particle/daughters/PDGs", "Daughters PDGs " + tit, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/daughters/PDGsPrimaries", "Daughters PDGs Primaries of " + tit, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/daughters/PDGsSecondaries", "Daughters PDGs Secondaries of " + tit, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/daughters/prodVx", "Daughters Prod. Vertex X " + tit, kTH1D, {axisProdx});
    histos.add("particle/daughters/prodVy", "Daughters Prod. Vertex Y " + tit, kTH1D, {axisPrody});
    histos.add("particle/daughters/prodVz", "Daughters Prod. Vertex Z " + tit, kTH1D, {axisProdz});
    histos.add("particle/daughters/prodDistanceX", "Daughters Prod. distance " + tit, kTH1D, {axisProdx});
    histos.add("particle/daughters/prodDistanceY", "Daughters Prod. distance " + tit, kTH1D, {axisPrody});
    histos.add("particle/daughters/prodDistanceZ", "Daughters Prod. distance " + tit, kTH1D, {axisProdz});
    histos.add("particle/daughters/prodDistanceVsPt", "Daughters Prod. distance " + tit, kTH2D, {axisPt, axisProdRadius});
    histos.add("particle/daughters/prodRadiusVsPt", "Daughters Prod. Vertex Radius " + tit, kTH2D, {axisPt, axisProdRadius});
    histos.add("particle/daughters/prodRadius3DVsPt", "Daughters Prod. Vertex Radius XYZ " + tit, kTH2D, {axisPt, axisProdRadius});

    histos.add("particle/mothers/Number", "Number of Mothers " + tit, kTH1D, {{20, -0.5, 19.5}});
    histos.add("particle/mothers/PDGs", "Mothers PDGs " + tit, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/mothers/PDGsPrimaries", "Mothers PDGs Primaries of " + tit, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/mothers/PDGsSecondaries", "Mothers PDGs Secondaries of " + tit, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/mothers/prodVx", "Mothers Prod. Vertex X " + tit, kTH1D, {axisProdx});
    histos.add("particle/mothers/prodVy", "Mothers Prod. Vertex Y " + tit, kTH1D, {axisPrody});
    histos.add("particle/mothers/prodVz", "Mothers Prod. Vertex Z " + tit, kTH1D, {axisProdz});
    histos.add("particle/mothers/prodRadiusVsPt", "Mothers Prod. Vertex Radius " + tit, kTH2D, {axisPt, axisProdRadius});
    histos.add("particle/mothers/prodRadius3DVsPt", "Mothers Prod. Vertex Radius XYZ " + tit, kTH2D, {axisPt, axisProdRadius});

    // Go up one generation
    histos.add("particle/mothers/mothers/Number", "Number of Mothers mothers " + tit, kTH1D, {{20, -0.5, 19.5}});
    histos.add("particle/mothers/mothers/PDGs", "Mothers mothers PDGs " + tit, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/mothers/mothers/PDGsPrimaries", "Mothers mothers PDGs Primaries of " + tit, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/mothers/mothers/PDGsSecondaries", "Mothers mothers PDGs Secondaries of " + tit, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/mothers/mothers/prodVx", "Mothers mothers Prod. Vertex X " + tit, kTH1D, {axisProdx});
    histos.add("particle/mothers/mothers/prodVy", "Mothers mothers Prod. Vertex Y " + tit, kTH1D, {axisPrody});
    histos.add("particle/mothers/mothers/prodVz", "Mothers mothers Prod. Vertex Z " + tit, kTH1D, {axisProdz});
    histos.add("particle/mothers/mothers/prodRadiusVsPt", "Mothers mothers Prod. Vertex Radius " + tit, kTH2D, {axisPt, axisProdRadius});
    histos.add("particle/mothers/mothers/prodRadius3DVsPt", "Mothers mothers Prod. Vertex Radius XYZ " + tit, kTH2D, {axisPt, axisProdRadius});

    if (doprocessParticleOnly) {
      return;
    }
    histos.add("track/PDGs", "Track PDGs", kTH2D, {axisPDGs, axisCharge});
    histos.add("track/withoutParticle", "Tracks without particles", kTH1D, {axisPt});
    histos.add("track/tofPDGs", "Track wTOF PDGs", kTH2D, {axisPDGs, axisCharge});
    histos.add("track/Pt", "Track Pt " + tit, kTH1D, {axisPt});
    histos.add("track/primariesPt", "Track Pt (primary) " + tit, kTH1D, {axisPt});
    histos.add("track/secondariesPt", "Track Pt (secondary) " + tit, kTH1D, {axisPt});
    histos.add("track/PtvsMCPt", "Track Pt vs MC Pt " + tit, kTH2D, {axisPt, axisPtMC});
    histos.add("track/Eta", "Track Eta " + tit, kTH1D, {axisEta});
    histos.add("track/primariesEta", "Track Eta (primary) " + tit, kTH1D, {axisEta});
    histos.add("track/secondariesEta", "Track Eta (secondary) " + tit, kTH1D, {axisEta});
    histos.add("track/primaries", "Source for primaries " + tit, kTH2D, {axisPDGs, axisCharge});
    histos.add("track/secondaries", "Source for secondaries " + tit, kTH2D, {axisPDGs, axisCharge});
  }

  template <typename ParticleType>
  double getCharge(ParticleType const& particle)
  {
    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle) {
      return 10.f;
    }
    return pdgParticle->Charge();
  }

  template <typename ParticleType>
  const char* getPdgCodeString(ParticleType const& particle)
  {
    return Form("%i", particle.pdgCode());
  }

  void processStandard(const o2::aod::McCollisions& colls,
                       const soa::Join<o2::aod::TracksIU, o2::aod::McTrackLabels, o2::aod::TracksExtra>& tracks,
                       const aod::McParticles& mcParticles)
  {
    for (const auto& col : colls) {
      histos.fill(HIST("event/VtxX"), col.posX());
      histos.fill(HIST("event/VtxY"), col.posY());
      histos.fill(HIST("event/VtxZ"), col.posZ());
    }
    std::vector<int64_t> ParticlesOfInterest;
    for (const auto& mcParticle : mcParticles) {
      const auto& pdgString = getPdgCodeString(mcParticle);
      const auto& pdgCharge = getCharge(mcParticle);
      histos.get<TH2>(HIST("particle/PDGs"))->Fill(pdgString, pdgCharge, 1.f);
      if (mcParticle.isPhysicalPrimary()) {
        histos.get<TH2>(HIST("particle/PDGsPrimaries"))->Fill(pdgString, pdgCharge, 1.f);
      } else {
        histos.get<TH2>(HIST("particle/PDGsSecondaries"))->Fill(pdgString, pdgCharge, 1.f);
      }
      if (mcParticle.pdgCode() != PDG) {
        continue;
      }
      if (mcParticle.y() < yMin || mcParticle.y() > yMax) {
        continue;
      }
      if (mcParticle.eta() < etaMin || mcParticle.eta() > etaMax) {
        continue;
      }
      if (mcParticle.pt() < ptMin || mcParticle.pt() > ptMax) {
        continue;
      }
      if (mcParticle.vz() < prodMinZ || mcParticle.vz() > prodMaxZ) {
        continue;
      }
      histos.fill(HIST("particle/producedProcess"), mcParticle.getProcess());
      histos.fill(HIST("particle/Pt"), mcParticle.pt());
      histos.fill(HIST("particle/P"), mcParticle.p());
      histos.fill(HIST("particle/Eta"), mcParticle.eta());
      if (mcParticle.isPhysicalPrimary()) {
        histos.fill(HIST("particle/primariesPt"), mcParticle.pt());
        histos.fill(HIST("particle/primariesP"), mcParticle.p());
        histos.fill(HIST("particle/primariesEta"), mcParticle.eta());
        histos.fill(HIST("particle/primariesY"), mcParticle.y());
      } else {
        histos.fill(HIST("particle/secondariesPt"), mcParticle.pt());
        histos.fill(HIST("particle/secondariesP"), mcParticle.p());
        histos.fill(HIST("particle/secondariesEta"), mcParticle.eta());
        histos.fill(HIST("particle/secondariesY"), mcParticle.y());
      }
      histos.fill(HIST("particle/EvsPz"), mcParticle.e(), mcParticle.pz());
      histos.fill(HIST("particle/Y"), mcParticle.y());
      histos.fill(HIST("particle/YvzPz"), mcParticle.y(), mcParticle.pz());
      histos.fill(HIST("particle/EtavzPz"), mcParticle.eta(), mcParticle.pz());
      histos.fill(HIST("particle/PvzPz"), mcParticle.p(), mcParticle.pz());
      histos.fill(HIST("particle/PtvzPz"), mcParticle.pt(), mcParticle.pz());
      histos.fill(HIST("particle/Px"), mcParticle.px());
      histos.fill(HIST("particle/Py"), mcParticle.py());
      histos.fill(HIST("particle/Pz"), mcParticle.pz());
      histos.fill(HIST("particle/prodVx"), mcParticle.vx());
      histos.fill(HIST("particle/prodVy"), mcParticle.vy());
      histos.fill(HIST("particle/prodVz"), mcParticle.vz());
      if (mcParticle.has_daughters()) {
        auto daughters = mcParticle.daughters_as<aod::McParticles>();
        histos.fill(HIST("particle/daughters/Number"), daughters.size());
        for (const auto& daughter : daughters) {
          const auto& pdgStringDau = getPdgCodeString(daughter);
          const auto& pdgChargeDau = getCharge(daughter);

          histos.get<TH2>(HIST("particle/daughters/PDGs"))->Fill(pdgStringDau, pdgChargeDau, 1.f);
          if (mcParticle.isPhysicalPrimary()) {
            histos.get<TH2>(HIST("particle/daughters/PDGsPrimaries"))->Fill(pdgStringDau, pdgChargeDau, 1.f);
          } else {
            histos.get<TH2>(HIST("particle/daughters/PDGsSecondaries"))->Fill(pdgStringDau, pdgChargeDau, 1.f);
          }

          histos.fill(HIST("particle/daughters/prodVx"), daughter.vx());
          histos.fill(HIST("particle/daughters/prodVy"), daughter.vy());
          histos.fill(HIST("particle/daughters/prodVz"), daughter.vz());
          histos.fill(HIST("particle/daughters/prodDistanceX"), daughter.vx() - mcParticle.vx());
          histos.fill(HIST("particle/daughters/prodDistanceY"), daughter.vy() - mcParticle.vy());
          histos.fill(HIST("particle/daughters/prodDistanceZ"), daughter.vz() - mcParticle.vz());
          histos.fill(HIST("particle/daughters/prodDistanceVsPt"), mcParticle.pt(), std::sqrt((daughter.vx() - mcParticle.vx()) * (daughter.vx() - mcParticle.vx()) + (daughter.vy() - mcParticle.vy()) * (daughter.vy() - mcParticle.vy()) + (daughter.vz() - mcParticle.vz()) * (daughter.vz() - mcParticle.vz())));
          histos.fill(HIST("particle/daughters/prodRadiusVsPt"), mcParticle.pt(), std::sqrt(daughter.vx() * daughter.vx() + daughter.vy() * daughter.vy()));
          histos.fill(HIST("particle/daughters/prodRadius3DVsPt"), mcParticle.pt(), std::sqrt(daughter.vx() * daughter.vx() + daughter.vy() * daughter.vy() + daughter.vz() * daughter.vz()));
        }
      } else {
        histos.fill(HIST("particle/daughters/Number"), 0.f);
      }
      if (mcParticle.has_mothers()) {
        const auto& mothers = mcParticle.mothers_as<aod::McParticles>();
        for (const auto& mother : mothers) {
          const auto& pdgStringMot = getPdgCodeString(mother);
          const auto& pdgChargeMot = getCharge(mother);

          histos.get<TH2>(HIST("particle/mothers/PDGs"))->Fill(pdgStringMot, pdgChargeMot, 1.f);
          if (mcParticle.isPhysicalPrimary()) {
            histos.get<TH2>(HIST("particle/mothers/PDGsPrimaries"))->Fill(pdgStringMot, pdgChargeMot, 1.f);
          } else {
            histos.get<TH2>(HIST("particle/mothers/PDGsSecondaries"))->Fill(pdgStringMot, pdgChargeMot, 1.f);
          }

          histos.fill(HIST("particle/mothers/prodVx"), mother.vx());
          histos.fill(HIST("particle/mothers/prodVy"), mother.vy());
          histos.fill(HIST("particle/mothers/prodVz"), mother.vz());
          histos.fill(HIST("particle/mothers/prodRadiusVsPt"), mother.pt(), std::sqrt(mother.vx() * mother.vx() + mother.vy() * mother.vy()));
          histos.fill(HIST("particle/mothers/prodRadius3DVsPt"), mother.pt(), std::sqrt(mother.vx() * mother.vx() + mother.vy() * mother.vy() + mother.vz() * mother.vz()));
          if (mother.has_mothers()) {
            const auto& mothers2 = mother.mothers_as<aod::McParticles>();
            for (const auto& mother2 : mothers2) {
              const auto& pdgStringMot2 = getPdgCodeString(mother2);
              const auto& pdgChargeMot2 = getCharge(mother2);

              histos.get<TH2>(HIST("particle/mothers/mothers/PDGs"))->Fill(pdgStringMot2, pdgChargeMot2, 1.f);
              if (mcParticle.isPhysicalPrimary()) {
                histos.get<TH2>(HIST("particle/mothers/mothers/PDGsPrimaries"))->Fill(pdgStringMot2, pdgChargeMot2, 1.f);
              } else {
                histos.get<TH2>(HIST("particle/mothers/mothers/PDGsSecondaries"))->Fill(pdgStringMot2, pdgChargeMot2, 1.f);
              }

              histos.fill(HIST("particle/mothers/mothers/prodVx"), mother2.vx());
              histos.fill(HIST("particle/mothers/mothers/prodVy"), mother2.vy());
              histos.fill(HIST("particle/mothers/mothers/prodVz"), mother2.vz());
              histos.fill(HIST("particle/mothers/mothers/prodRadiusVsPt"), mother2.pt(), std::sqrt(mother2.vx() * mother2.vx() + mother2.vy() * mother2.vy()));
              histos.fill(HIST("particle/mothers/mothers/prodRadius3DVsPt"), mother2.pt(), std::sqrt(mother2.vx() * mother2.vx() + mother2.vy() * mother2.vy() + mother2.vz() * mother2.vz()));
            }
          }
        }
      }

      histos.fill(HIST("particle/prodRadius"), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
      histos.fill(HIST("particle/prodVxVy"), mcParticle.vx(), mcParticle.vy());
      if (!mcParticle.isPhysicalPrimary()) {
        if (mcParticle.getProcess() == 4) {
          histos.fill(HIST("particle/prodVxVyStr"), mcParticle.vx(), mcParticle.vy());
        } else {
          histos.fill(HIST("particle/prodVxVyMat"), mcParticle.vx(), mcParticle.vy());
        }
      } else {
        histos.fill(HIST("particle/prodVxVyPrm"), mcParticle.vx(), mcParticle.vy());
      }
      histos.fill(HIST("particle/prodVxVsPt"), mcParticle.pt(), mcParticle.vx());
      histos.fill(HIST("particle/prodVyVsPt"), mcParticle.pt(), mcParticle.vy());
      histos.fill(HIST("particle/prodVzVsPt"), mcParticle.pt(), mcParticle.vz());
      histos.fill(HIST("particle/prodRadiusVsPt"), mcParticle.pt(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
      histos.fill(HIST("particle/prodRadiusVsEta"), mcParticle.eta(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
      histos.fill(HIST("particle/prodRadius3DVsPt"), mcParticle.pt(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy() + mcParticle.vz() * mcParticle.vz()));
      ParticlesOfInterest.push_back(mcParticle.globalIndex());
    }

    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        histos.fill(HIST("track/withoutParticle"), track.pt());
        continue;
      }
      const auto mcParticle = track.mcParticle_as<aod::McParticles>();
      histos.get<TH2>(HIST("track/PDGs"))->Fill(getPdgCodeString(mcParticle), getCharge(mcParticle), 1.f);
      if (track.hasTOF()) {
        histos.get<TH2>(HIST("track/tofPDGs"))->Fill(getPdgCodeString(mcParticle), getCharge(mcParticle), 1.f);
      }
      if (IsStable.value == 0) {
        LOG(info) << mcParticle.pdgCode() << " asked for " << PDG.value;
        if (!mcParticle.has_mothers()) {
          continue;
        }
        // auto mothers = mcParticle.mothers();
        auto mothers = mcParticle.mothers_as<aod::McParticles>();
        const auto ParticleIsInteresting = std::find(ParticlesOfInterest.begin(), ParticlesOfInterest.end(), mothers[0].globalIndex()) != ParticlesOfInterest.end();
        if (!ParticleIsInteresting) {
          continue;
        }
        if (doPrint) {
          LOG(info) << "Track " << track.globalIndex() << " comes from a " << mothers[0].pdgCode() << " and is a " << mcParticle.pdgCode();
        }
      } else {
        if (mcParticle.pdgCode() != PDG) {
          continue;
        }
        histos.fill(HIST("track/Pt"), track.pt() * charge);
        histos.fill(HIST("track/PtvsMCPt"), track.pt() * charge, mcParticle.pt());
        histos.fill(HIST("track/Eta"), track.eta());
        if (mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("track/primariesPt"), track.pt() * charge);
          histos.fill(HIST("track/primariesEta"), track.eta());
        } else {
          histos.fill(HIST("track/secondariesPt"), track.pt() * charge);
          histos.fill(HIST("track/secondariesEta"), track.eta());
        }

        if (!mcParticle.has_mothers()) {
          if (doPrint) {
            LOG(info) << "Track " << track.globalIndex() << " is a " << mcParticle.pdgCode();
          }
          continue;
        }
        // auto mothers = mcParticle.mothers();
        auto mothers = mcParticle.mothers_as<aod::McParticles>();
        if (mcParticle.isPhysicalPrimary()) {
          histos.get<TH2>(HIST("track/primaries"))->Fill(getPdgCodeString(mothers[0]), getCharge(mothers[0]), 1.f);
        } else {
          histos.get<TH2>(HIST("track/secondaries"))->Fill(getPdgCodeString(mothers[0]), getCharge(mothers[0]), 1.f);
        }
        if (doPrint) {
          LOG(info) << "Track " << track.globalIndex() << " is a " << mcParticle.pdgCode() << " and comes from a " << mothers[0].pdgCode() << " and is " << (mcParticle.isPhysicalPrimary() ? "" : "not") << " a primary";
        }
      }
    }
  }
  PROCESS_SWITCH(Alice3SingleParticle, processStandard, "Process IU tracks", true);

  void processParticleOnly(const o2::aod::McCollisions& colls,
                           const aod::McParticles& mcParticles)
  {
    for (const auto& col : colls) {
      histos.fill(HIST("event/VtxX"), col.posX());
      histos.fill(HIST("event/VtxY"), col.posY());
      histos.fill(HIST("event/VtxZ"), col.posZ());
    }
    for (const auto& mcParticle : mcParticles) {
      const auto& pdgString = getPdgCodeString(mcParticle);
      const auto& pdgCharge = getCharge(mcParticle);
      histos.get<TH2>(HIST("particle/PDGs"))->Fill(pdgString, pdgCharge, 1.f);
      if (mcParticle.isPhysicalPrimary()) {
        histos.get<TH2>(HIST("particle/PDGsPrimaries"))->Fill(pdgString, pdgCharge, 1.f);
      } else {
        histos.get<TH2>(HIST("particle/PDGsSecondaries"))->Fill(pdgString, pdgCharge, 1.f);
      }
      if (mcParticle.pdgCode() != PDG) {
        continue;
      }
      if (mcParticle.y() < yMin || mcParticle.y() > yMax) {
        continue;
      }
      if (mcParticle.eta() < etaMin || mcParticle.eta() > etaMax) {
        continue;
      }
      if (mcParticle.vz() < prodMinZ || mcParticle.vz() > prodMaxZ) {
        continue;
      }
      histos.fill(HIST("particle/producedProcess"), mcParticle.getProcess());
      histos.fill(HIST("particle/Pt"), mcParticle.pt());
      histos.fill(HIST("particle/P"), mcParticle.p());
      histos.fill(HIST("particle/Eta"), mcParticle.eta());
      if (mcParticle.isPhysicalPrimary()) {
        histos.fill(HIST("particle/primariesPt"), mcParticle.pt());
        histos.fill(HIST("particle/primariesP"), mcParticle.p());
        histos.fill(HIST("particle/primariesEta"), mcParticle.eta());
        histos.fill(HIST("particle/primariesY"), mcParticle.y());
      } else {
        histos.fill(HIST("particle/secondariesPt"), mcParticle.pt());
        histos.fill(HIST("particle/secondariesP"), mcParticle.p());
        histos.fill(HIST("particle/secondariesEta"), mcParticle.eta());
        histos.fill(HIST("particle/secondariesY"), mcParticle.y());
      }
      histos.fill(HIST("particle/EvsPz"), mcParticle.e(), mcParticle.pz());
      histos.fill(HIST("particle/Y"), mcParticle.y());
      histos.fill(HIST("particle/YvzPz"), mcParticle.y(), mcParticle.pz());
      histos.fill(HIST("particle/EtavzPz"), mcParticle.eta(), mcParticle.pz());
      histos.fill(HIST("particle/PvzPz"), mcParticle.p(), mcParticle.pz());
      histos.fill(HIST("particle/PtvzPz"), mcParticle.pt(), mcParticle.pz());
      histos.fill(HIST("particle/Px"), mcParticle.px());
      histos.fill(HIST("particle/Py"), mcParticle.py());
      histos.fill(HIST("particle/Pz"), mcParticle.pz());
      histos.fill(HIST("particle/prodVx"), mcParticle.vx());
      histos.fill(HIST("particle/prodVy"), mcParticle.vy());
      histos.fill(HIST("particle/prodVz"), mcParticle.vz());
      if (mcParticle.has_daughters()) {
        auto daughters = mcParticle.daughters_as<aod::McParticles>();
        histos.fill(HIST("particle/daughters/Number"), daughters.size());
        for (const auto& daughter : daughters) {
          const auto& pdgStringDau = getPdgCodeString(daughter);
          const auto& pdgChargeDau = getCharge(daughter);

          histos.get<TH2>(HIST("particle/daughters/PDGs"))->Fill(pdgStringDau, pdgChargeDau, 1.f);
          if (mcParticle.isPhysicalPrimary()) {
            histos.get<TH2>(HIST("particle/daughters/PDGsPrimaries"))->Fill(pdgStringDau, pdgChargeDau, 1.f);
          } else {
            histos.get<TH2>(HIST("particle/daughters/PDGsSecondaries"))->Fill(pdgStringDau, pdgChargeDau, 1.f);
          }

          histos.fill(HIST("particle/daughters/prodVx"), daughter.vx());
          histos.fill(HIST("particle/daughters/prodVy"), daughter.vy());
          histos.fill(HIST("particle/daughters/prodVz"), daughter.vz());
          histos.fill(HIST("particle/daughters/prodRadiusVsPt"), mcParticle.pt(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
          histos.fill(HIST("particle/daughters/prodRadius3DVsPt"), mcParticle.pt(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy() + mcParticle.vz() * mcParticle.vz()));
        }
      } else {
        histos.fill(HIST("particle/daughters/Number"), 0.f);
      }
      if (mcParticle.has_mothers()) {
        auto mothers = mcParticle.mothers_as<aod::McParticles>();
        for (const auto& mother : mothers) {
          const auto& pdgStringMot = getPdgCodeString(mother);
          const auto& pdgChargeMot = getCharge(mother);

          histos.get<TH2>(HIST("particle/mothers/PDGs"))->Fill(pdgStringMot, pdgChargeMot, 1.f);
          if (mcParticle.isPhysicalPrimary()) {
            histos.get<TH2>(HIST("particle/mothers/PDGsPrimaries"))->Fill(pdgStringMot, pdgChargeMot, 1.f);
          } else {
            histos.get<TH2>(HIST("particle/mothers/PDGsSecondaries"))->Fill(pdgStringMot, pdgChargeMot, 1.f);
          }

          histos.fill(HIST("particle/mothers/prodVx"), mother.vx());
          histos.fill(HIST("particle/mothers/prodVy"), mother.vy());
          histos.fill(HIST("particle/mothers/prodVz"), mother.vz());
          histos.fill(HIST("particle/mothers/prodRadiusVsPt"), mother.pt(), std::sqrt(mother.vx() * mother.vx() + mother.vy() * mother.vy()));
          histos.fill(HIST("particle/mothers/prodRadius3DVsPt"), mother.pt(), std::sqrt(mother.vx() * mother.vx() + mother.vy() * mother.vy() + mother.vz() * mother.vz()));
        }
      }

      histos.fill(HIST("particle/prodRadius"), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
      histos.fill(HIST("particle/prodVxVy"), mcParticle.vx(), mcParticle.vy());
      if (!mcParticle.isPhysicalPrimary()) {
        if (mcParticle.getProcess() == 4) {
          histos.fill(HIST("particle/prodVxVyStr"), mcParticle.vx(), mcParticle.vy());
        } else {
          histos.fill(HIST("particle/prodVxVyMat"), mcParticle.vx(), mcParticle.vy());
        }
      } else {
        histos.fill(HIST("particle/prodVxVyPrm"), mcParticle.vx(), mcParticle.vy());
      }
      histos.fill(HIST("particle/prodVxVsPt"), mcParticle.pt(), mcParticle.vx());
      histos.fill(HIST("particle/prodVyVsPt"), mcParticle.pt(), mcParticle.vy());
      histos.fill(HIST("particle/prodVzVsPt"), mcParticle.pt(), mcParticle.vz());
      histos.fill(HIST("particle/prodRadiusVsPt"), mcParticle.pt(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
      histos.fill(HIST("particle/prodRadiusVsEta"), mcParticle.eta(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
      histos.fill(HIST("particle/prodRadius3DVsPt"), mcParticle.pt(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy() + mcParticle.vz() * mcParticle.vz()));
    }
  }
  PROCESS_SWITCH(Alice3SingleParticle, processParticleOnly, "Process Particle only", false);

  void processNonIU(const o2::aod::McCollisions& colls,
                    const soa::Join<o2::aod::Tracks, o2::aod::McTrackLabels, o2::aod::TracksExtra>& tracks,
                    const aod::McParticles& mcParticles)
  {
    for (const auto& col : colls) {
      histos.fill(HIST("event/VtxX"), col.posX());
      histos.fill(HIST("event/VtxY"), col.posY());
      histos.fill(HIST("event/VtxZ"), col.posZ());
    }
    std::vector<int64_t> ParticlesOfInterest;
    for (const auto& mcParticle : mcParticles) {
      const auto& pdgString = getPdgCodeString(mcParticle);
      const auto& pdgCharge = getCharge(mcParticle);
      histos.get<TH2>(HIST("particle/PDGs"))->Fill(pdgString, pdgCharge, 1.f);
      if (mcParticle.isPhysicalPrimary()) {
        histos.get<TH2>(HIST("particle/PDGsPrimaries"))->Fill(pdgString, pdgCharge, 1.f);
      } else {
        histos.get<TH2>(HIST("particle/PDGsSecondaries"))->Fill(pdgString, pdgCharge, 1.f);
      }
      if (mcParticle.pdgCode() != PDG) {
        continue;
      }
      if (mcParticle.y() < yMin || mcParticle.y() > yMax) {
        continue;
      }
      if (mcParticle.eta() < etaMin || mcParticle.eta() > etaMax) {
        continue;
      }
      if (mcParticle.pt() < ptMin || mcParticle.pt() > ptMax) {
        continue;
      }
      histos.fill(HIST("particle/Pt"), mcParticle.pt());
      histos.fill(HIST("particle/P"), mcParticle.p());
      histos.fill(HIST("particle/Eta"), mcParticle.eta());
      if (mcParticle.isPhysicalPrimary()) {
        histos.fill(HIST("particle/primariesPt"), mcParticle.pt());
        histos.fill(HIST("particle/primariesP"), mcParticle.p());
        histos.fill(HIST("particle/primariesEta"), mcParticle.eta());
        histos.fill(HIST("particle/primariesY"), mcParticle.y());
      } else {
        histos.fill(HIST("particle/secondariesPt"), mcParticle.pt());
        histos.fill(HIST("particle/secondariesP"), mcParticle.p());
        histos.fill(HIST("particle/secondariesEta"), mcParticle.eta());
        histos.fill(HIST("particle/secondariesY"), mcParticle.y());
      }
      histos.fill(HIST("particle/EvsPz"), mcParticle.e(), mcParticle.pz());
      histos.fill(HIST("particle/Y"), mcParticle.y());
      histos.fill(HIST("particle/YvzPz"), mcParticle.y(), mcParticle.pz());
      histos.fill(HIST("particle/EtavzPz"), mcParticle.eta(), mcParticle.pz());
      histos.fill(HIST("particle/PvzPz"), mcParticle.p(), mcParticle.pz());
      histos.fill(HIST("particle/PtvzPz"), mcParticle.pt(), mcParticle.pz());
      histos.fill(HIST("particle/Px"), mcParticle.px());
      histos.fill(HIST("particle/Py"), mcParticle.py());
      histos.fill(HIST("particle/Pz"), mcParticle.pz());
      histos.fill(HIST("particle/prodVx"), mcParticle.vx());
      histos.fill(HIST("particle/prodVy"), mcParticle.vy());
      histos.fill(HIST("particle/prodVz"), mcParticle.vz());
      if (mcParticle.has_daughters()) {
        auto daughters = mcParticle.daughters_as<aod::McParticles>();
        histos.fill(HIST("particle/daughters/Number"), daughters.size());
        for (const auto& daughter : daughters) {
          const auto& pdgStringDau = getPdgCodeString(daughter);
          const auto& pdgChargeDau = getCharge(daughter);

          histos.get<TH2>(HIST("particle/daughters/PDGs"))->Fill(pdgStringDau, pdgChargeDau, 1.f);
          if (mcParticle.isPhysicalPrimary()) {
            histos.get<TH2>(HIST("particle/daughters/PDGsPrimaries"))->Fill(pdgStringDau, pdgChargeDau, 1.f);
          } else {
            histos.get<TH2>(HIST("particle/daughters/PDGsSecondaries"))->Fill(pdgStringDau, pdgChargeDau, 1.f);
          }

          histos.fill(HIST("particle/daughters/prodVx"), daughter.vx());
          histos.fill(HIST("particle/daughters/prodVy"), daughter.vy());
          histos.fill(HIST("particle/daughters/prodVz"), daughter.vz());
          histos.fill(HIST("particle/daughters/prodRadiusVsPt"), mcParticle.pt(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
          histos.fill(HIST("particle/daughters/prodRadius3DVsPt"), mcParticle.pt(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy() + mcParticle.vz() * mcParticle.vz()));
        }
      } else {
        histos.fill(HIST("particle/daughters/Number"), 0.f);
      }
      if (mcParticle.has_mothers()) {
        auto mothers = mcParticle.mothers_as<aod::McParticles>();
        for (const auto& mother : mothers) {
          const auto& pdgStringMot = getPdgCodeString(mother);
          const auto& pdgChargeMot = getCharge(mother);

          histos.get<TH2>(HIST("particle/mothers/PDGs"))->Fill(pdgStringMot, pdgChargeMot, 1.f);
          if (mcParticle.isPhysicalPrimary()) {
            histos.get<TH2>(HIST("particle/mothers/PDGsPrimaries"))->Fill(pdgStringMot, pdgChargeMot, 1.f);
          } else {
            histos.get<TH2>(HIST("particle/mothers/PDGsSecondaries"))->Fill(pdgStringMot, pdgChargeMot, 1.f);
          }

          histos.fill(HIST("particle/mothers/prodVx"), mother.vx());
          histos.fill(HIST("particle/mothers/prodVy"), mother.vy());
          histos.fill(HIST("particle/mothers/prodVz"), mother.vz());
          histos.fill(HIST("particle/mothers/prodRadiusVsPt"), mother.pt(), std::sqrt(mother.vx() * mother.vx() + mother.vy() * mother.vy()));
          histos.fill(HIST("particle/mothers/prodRadius3DVsPt"), mother.pt(), std::sqrt(mother.vx() * mother.vx() + mother.vy() * mother.vy() + mother.vz() * mother.vz()));
        }
      }

      histos.fill(HIST("particle/prodRadius"), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
      histos.fill(HIST("particle/prodVxVy"), mcParticle.vx(), mcParticle.vy());
      if (!mcParticle.isPhysicalPrimary()) {
        if (mcParticle.getProcess() == 4) {
          histos.fill(HIST("particle/prodVxVyStr"), mcParticle.vx(), mcParticle.vy());
        } else {
          histos.fill(HIST("particle/prodVxVyMat"), mcParticle.vx(), mcParticle.vy());
        }
      } else {
        histos.fill(HIST("particle/prodVxVyPrm"), mcParticle.vx(), mcParticle.vy());
      }
      histos.fill(HIST("particle/prodVxVsPt"), mcParticle.pt(), mcParticle.vx());
      histos.fill(HIST("particle/prodVyVsPt"), mcParticle.pt(), mcParticle.vy());
      histos.fill(HIST("particle/prodVzVsPt"), mcParticle.pt(), mcParticle.vz());
      histos.fill(HIST("particle/prodRadiusVsPt"), mcParticle.pt(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
      histos.fill(HIST("particle/prodRadiusVsEta"), mcParticle.eta(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
      histos.fill(HIST("particle/prodRadius3DVsPt"), mcParticle.pt(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy() + mcParticle.vz() * mcParticle.vz()));
      ParticlesOfInterest.push_back(mcParticle.globalIndex());
    }

    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        histos.fill(HIST("track/withoutParticle"), track.pt());
        continue;
      }
      const auto mcParticle = track.mcParticle_as<aod::McParticles>();
      histos.get<TH2>(HIST("track/PDGs"))->Fill(getPdgCodeString(mcParticle), getCharge(mcParticle), 1.f);
      if (track.hasTOF()) {
        histos.get<TH2>(HIST("track/tofPDGs"))->Fill(getPdgCodeString(mcParticle), getCharge(mcParticle), 1.f);
      }
      if (IsStable.value == 0) {
        if (!mcParticle.has_mothers()) {
          continue;
        }
        const auto& mothers = mcParticle.mothers_as<aod::McParticles>();
        const auto& ParticleIsInteresting = std::find(ParticlesOfInterest.begin(), ParticlesOfInterest.end(), mothers[0].globalIndex()) != ParticlesOfInterest.end();
        if (!ParticleIsInteresting) {
          continue;
        }
        if (doPrint) {
          LOG(info) << "Track " << track.globalIndex() << " comes from a " << mothers[0].pdgCode() << " and is a " << mcParticle.pdgCode();
        }
      } else {
        if (mcParticle.pdgCode() != PDG) {
          continue;
        }
        histos.fill(HIST("track/Pt"), track.pt() * charge);
        histos.fill(HIST("track/PtvsMCPt"), track.pt() * charge, mcParticle.pt());
        histos.fill(HIST("track/Eta"), track.eta());
        if (mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("track/primariesPt"), track.pt() * charge);
          histos.fill(HIST("track/primariesEta"), track.eta());
        } else {
          histos.fill(HIST("track/secondariesPt"), track.pt() * charge);
          histos.fill(HIST("track/secondariesEta"), track.eta());
        }

        if (!mcParticle.has_mothers()) {
          if (doPrint) {
            LOG(info) << "Track " << track.globalIndex() << " is a " << mcParticle.pdgCode();
          }
          continue;
        }
        // auto mothers = mcParticle.mothers();
        auto mothers = mcParticle.mothers_as<aod::McParticles>();
        if (mcParticle.isPhysicalPrimary()) {
          histos.get<TH2>(HIST("track/primaries"))->Fill(getPdgCodeString(mothers[0]), getCharge(mothers[0]), 1.f);
        } else {
          histos.get<TH2>(HIST("track/secondaries"))->Fill(getPdgCodeString(mothers[0]), getCharge(mothers[0]), 1.f);
        }
        if (doPrint) {
          LOG(info) << "Track " << track.globalIndex() << " is a " << mcParticle.pdgCode() << " and comes from a " << mothers[0].pdgCode() << " and is " << (mcParticle.isPhysicalPrimary() ? "" : "not") << " a primary";
        }
      }
    }
  }
  PROCESS_SWITCH(Alice3SingleParticle, processNonIU, "Process non IU tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3SingleParticle>(cfgc)};
}
