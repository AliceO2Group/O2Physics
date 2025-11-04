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

/// \author Nicolo' Jacazio <nicolo.jacazio@cern.ch>, CERN
/// \author Alexander Kalweit <alexander.kalweit@cern.ch>, CERN

// O2 includes
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/PIDResponseTOF.h"

#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct Alice3CDeuteron {
  Configurable<float> magField{"magField", 0.5, "Magnetic field"};
  Configurable<float> minRadius{"minRadius", -100, "Minimum decay radius"};
  Configurable<float> maxRadius{"maxRadius", 100, "Maximum decay radius"};
  Configurable<float> minMomPt{"minMomPt", -100, "Minimum pT of the mother"};
  Configurable<float> minKaonPt{"minKaonPt", -100, "Minimum pT of the pion daughter"};
  Configurable<float> minPionPt{"minPionPt", -100, "Minimum pT of the kaon daughter"};
  Configurable<float> minVtxContrib{"minVtxContrib", 3, "Minimum number of contributors to the primary vertex"};
  Configurable<float> minDca{"minDca", -100, "Minimum track DCA to the primary vertex"};
  Configurable<float> minDcaDeuteron{"minDcaDeuteron", -100, "Minimum DCA of the Deuteron to the primary vertex"};
  Configurable<float> minDcaPion{"minDcaPion", -100, "Minimum DCA of the pion to the primary vertex"};
  Configurable<float> maxDca{"maxDca", 100, "Maximum track DCA to the primary vertex"};
  Configurable<float> minCpa{"minCpa", 0, "Minimum CPA"};
  Configurable<int> usePdg{"usePdg", 1, "Flag to use the PDG instead of the TOF/RICH PID"};
  Configurable<float> maxNsigmaDe{"maxNsigmaDe", 3.f, "Maximum Nsigma for deuteron"};
  Configurable<float> maxNsigmaKa{"maxNsigmaKa", 3.f, "Maximum Nsigma for kaon"};
  Configurable<float> maxNsigmaPi{"maxNsigmaPi", 3.f, "Maximum Nsigma for pion"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  o2::vertexing::DCAFitterN<3> fitter;

  void init(InitContext&)
  {

    fitter.setBz(magField);
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(1.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);

    const AxisSpec axisInvMass{1000, 2.5, 4, "Inv. Mass_{c-d}"};
    const AxisSpec axisDecayRadius{2000, 0, 0.1, "Decay radius"};
    const AxisSpec axisDecayRadiusReso{2000, -0.01, 0.01, "Decay radius resolution"};
    const AxisSpec axisPionProdRadiusXY{2000, 0, 0.01, "Pion production radius in xy"};
    const AxisSpec axisDca{5000, -0.01, 0.01, "DCA to secondary"};
    const AxisSpec axisDcaXY{5000, -0.05, 0.05, "DCA_{xy}"};
    const AxisSpec axisDcaXYProd{5000, -5e-6, 5e-6, "DCA_{xy} product"};
    const AxisSpec axisDcaZ{5000, -0.05, 0.05, "DCA_{z}"};
    const AxisSpec axisDcaZProd{5000, -5e-6, 5e-6, "DCA_{z} product"};
    const AxisSpec axisEta{100, -4, 4, "#it{#eta}"};
    const AxisSpec axisY{100, -4, 4, "#it{y}"};
    const AxisSpec axisPt{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisP{100, 0, 10, "#it{p} (GeV/#it{c})"};
    const AxisSpec axisVtxX{100, -0.1, 0.1, "Vtx_{X}"};
    const AxisSpec axisVtxY{100, -0.1, 0.1, "Vtx_{Y}"};
    const AxisSpec axisVtxZ{100, -0.1, 0.1, "Vtx_{Z}"};
    const AxisSpec axisCPA{4000, -1.1, 1.1, "CPA"};
    const AxisSpec axisNSigmaDe{1000, -10, 10, "N_{#sigma}^{TOF}(d)"};
    const AxisSpec axisNSigmaKa{1000, -10, 10, "N_{#sigma}^{TOF}(k)"};
    const AxisSpec axisNSigmaPi{1000, -10, 10, "N_{#sigma}^{TOF}(#pi)"};
    const TString tit = Form(" [%.6f, %.6f] R [%.6f, %.6f] DCA ",
                             minRadius.value, maxRadius.value,
                             minDca.value, maxDca.value);

    histos.add("event/candcuts", "cuts", kTH1D, {{10, 0, 10}});
    auto h = histos.get<TH1>(HIST("event/candcuts"));
    h->GetXaxis()->SetBinLabel(1, "magField");
    h->GetXaxis()->SetBinLabel(2, "minRadius");
    h->GetXaxis()->SetBinLabel(3, "maxRadius");
    h->GetXaxis()->SetBinLabel(4, "minMomPt");
    h->GetXaxis()->SetBinLabel(5, "minKaonPt");
    h->GetXaxis()->SetBinLabel(6, "minPionPt");
    h->GetXaxis()->SetBinLabel(7, "minVtxContrib");
    h->GetXaxis()->SetBinLabel(8, "minDca");
    h->GetXaxis()->SetBinLabel(9, "maxDca");

    h->SetBinContent(1, magField);
    h->SetBinContent(2, minRadius);
    h->SetBinContent(3, maxRadius);
    h->SetBinContent(4, minMomPt);
    h->SetBinContent(5, minKaonPt);
    h->SetBinContent(6, minPionPt);
    h->SetBinContent(7, minVtxContrib);
    h->SetBinContent(8, minDca);
    h->SetBinContent(9, maxDca);

    histos.add("cdeuteron/eta", "eta", kTH1D, {axisEta});
    histos.add("cdeuteron/y", "y", kTH1D, {axisY});
    histos.add("cdeuteron/pt", "pt", kTH1D, {axisPt});
    histos.add("cdeuteron/p", "p", kTH1D, {axisP});

    histos.add("event/vtxX", "vtxX", kTH1D, {axisVtxX});
    histos.add("event/vtxY", "vtxY", kTH1D, {axisVtxY});
    histos.add("event/vtxZ", "vtxZ", kTH1D, {axisVtxZ});
    histos.add("event/mcvtxX", "mcvtxX", kTH1D, {axisVtxX});
    histos.add("event/mcvtxY", "mcvtxY", kTH1D, {axisVtxY});
    histos.add("event/mcvtxZ", "mcvtxZ", kTH1D, {axisVtxZ});
    histos.add("event/nsigmaDe", "nsigmaDe", kTH2D, {axisPt, axisNSigmaDe});
    histos.add("event/nsigmaKa", "nsigmaKa", kTH2D, {axisPt, axisNSigmaKa});
    histos.add("event/nsigmaPi", "nsigmaPi", kTH2D, {axisPt, axisNSigmaPi});
    histos.add("event/nsigmaDecut", "nsigmaDecut", kTH2D, {axisPt, axisNSigmaDe});
    histos.add("event/nsigmaKacut", "nsigmaKacut", kTH2D, {axisPt, axisNSigmaKa});
    histos.add("event/nsigmaPicut", "nsigmaPicut", kTH2D, {axisPt, axisNSigmaPi});
    histos.add("event/track1dcaxy", "track1dcaxy Deuteron", kTH1D, {axisDcaXY});
    histos.add("event/track1dcaz", "track1dcaz Deuteron", kTH1D, {axisDcaZ});
    histos.add("event/candperdeuteron", "candperdeuteron", kTH1D, {{1000, 0, 10000}});
    histos.add("event/particlespdg", "particlespdg", kTH1D, {{100, 0, 100}});
    histos.add("event/trackspdg", "trackspdg", kTH1D, {{4, 0.5, 4.5}});
    h = histos.get<TH1>(HIST("event/trackspdg"));
    h->GetXaxis()->SetBinLabel(1, "d");
    h->GetXaxis()->SetBinLabel(2, "K");
    h->GetXaxis()->SetBinLabel(3, "#pi");
    h->GetXaxis()->SetBinLabel(4, "Rest");
    histos.add("event/multiplicity", "multiplicity", kTH1D, {{1000, 0, 10000}});

#define MakeHistos(tag)                                                                        \
  histos.add(tag "/cpa", "cpa" + tit, kTH1D, {axisCPA});                                       \
  histos.add(tag "/invmass", "invmass" + tit, kTH1D, {axisInvMass});                           \
  histos.add(tag "/invmassVsPt", "invmassVsPt" + tit, kTH2D, {axisPt, axisInvMass});           \
  histos.add(tag "/decayradius", "decayradius" + tit, kTH1D, {axisDecayRadius});               \
  histos.add(tag "/decayradiusResoX", "decayradiusResoX" + tit, kTH1D, {axisDecayRadiusReso}); \
  histos.add(tag "/decayradiusResoY", "decayradiusResoY" + tit, kTH1D, {axisDecayRadiusReso}); \
  histos.add(tag "/decayradiusResoZ", "decayradiusResoZ" + tit, kTH1D, {axisDecayRadiusReso}); \
  histos.add(tag "/decayradiusReso", "decayradiusReso" + tit, kTH1D, {axisDecayRadiusReso});   \
  histos.add(tag "/radius3xy", "radius3xy" + tit, kTH1D, {axisPionProdRadiusXY});              \
  histos.add(tag "/decaydca0", "decaydca0" + tit, kTH1D, {axisDca});                           \
  histos.add(tag "/decaydca1", "decaydca1" + tit, kTH1D, {axisDca});                           \
  histos.add(tag "/dcaxy1", "dcaxy1 Deuteron" + tit, kTH1D, {axisDcaXY});                      \
  histos.add(tag "/dcaxy2", "dcaxy2 Kaon" + tit, kTH1D, {axisDcaXY});                          \
  histos.add(tag "/dcaxy3", "dcaxy3 Pion" + tit, kTH1D, {axisDcaXY});                          \
  histos.add(tag "/dcaxy1xdcaxy2", "dcaxy1xdcaxy2" + tit, kTH1D, {axisDcaXYProd});             \
  histos.add(tag "/dcaxy3xdcaxy2", "dcaxy3xdcaxy2" + tit, kTH1D, {axisDcaXYProd});             \
  histos.add(tag "/dcaz1", "dcaz1 Deuteron" + tit, kTH1D, {axisDcaZ});                         \
  histos.add(tag "/dcaz2", "dcaz2 Kaon" + tit, kTH1D, {axisDcaZ});                             \
  histos.add(tag "/dcaz3", "dcaz3 Pion" + tit, kTH1D, {axisDcaZ});                             \
  histos.add(tag "/dcaz1xdcaz2", "dcaz1xdcaz2" + tit, kTH1D, {axisDcaZProd});                  \
  histos.add(tag "/dcaz3xdcaz2", "dcaz3xdcaz2" + tit, kTH1D, {axisDcaZProd});                  \
  histos.add(tag "/pt1", "pt1 Deuteron" + tit, kTH1D, {axisPt});                               \
  histos.add(tag "/pt2", "pt2 Kaon" + tit, kTH1D, {axisPt});                                   \
  histos.add(tag "/pt3", "pt3 Pion" + tit, kTH1D, {axisPt});                                   \
  histos.add(tag "/ptmom", "ptmom" + tit, kTH1D, {axisPt});                                    \
  histos.add(tag "/pmom", "pmom" + tit, kTH1D, {axisP});

    MakeHistos("sig");
    MakeHistos("bkg");
    MakeHistos("signocut");
    MakeHistos("bkgnocut");
    MakeHistos("sigcut");
    MakeHistos("bkgcut");

#undef MakeHistos
  }

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  SliceCache cache;
  template <bool usePID = true, typename collType, typename trackType, typename partType>
  void fillHistograms(const collType& coll, const trackType& tracks, const partType& mcParticles)
  {
    const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, coll.mcCollision().globalIndex(), cache);
    for (const auto& i : particlesInCollision) {
      histos.get<TH1>(HIST("event/particlespdg"))->Fill(Form("%i", i.pdgCode()), 1);
      if (i.pdgCode() != 12345) {
        continue;
      }
      histos.fill(HIST("cdeuteron/eta"), i.eta());
      histos.fill(HIST("cdeuteron/y"), i.y());
      histos.fill(HIST("cdeuteron/pt"), i.pt());
      histos.fill(HIST("cdeuteron/p"), i.p());
    }
    histos.fill(HIST("event/vtxX"), coll.posX());
    histos.fill(HIST("event/vtxY"), coll.posY());
    histos.fill(HIST("event/vtxZ"), coll.posZ());
    histos.fill(HIST("event/mcvtxX"), coll.mcCollision().posX());
    histos.fill(HIST("event/mcvtxY"), coll.mcCollision().posY());
    histos.fill(HIST("event/mcvtxZ"), coll.mcCollision().posZ());
    const math_utils::Point3D<float> collPos{coll.mcCollision().posX(),
                                             coll.mcCollision().posY(),
                                             coll.mcCollision().posZ()};
    // const math_utils::Point3D<float> collPos{coll.posX(),
    //                                          coll.posY(),
    //                                          coll.posZ()};
    // for (const auto& mcParticle : mcParticles) {
    //   // ParticlesOfInterest.push_back(mcParticle.globalIndex());
    // }
    int ntrks = 0;
    for (const auto& t : tracks) {
      if (!t.has_mcParticle()) {
        continue;
      }
      switch (t.template mcParticle_as<aod::McParticles>().pdgCode()) {
        case 1000010020:
          histos.fill(HIST("event/trackspdg"), 1);
          break;
        case -321:
          histos.fill(HIST("event/trackspdg"), 2);
          break;
        case 211:
          histos.fill(HIST("event/trackspdg"), 3);
          break;
        default:
          histos.fill(HIST("event/trackspdg"), 4);
      }
      ntrks++;
    }
    histos.fill(HIST("event/multiplicity"), ntrks);

    std::array<float, 2> dca1{1e10f, 1e10f};
    std::array<float, 2> dca2{1e10f, 1e10f};
    std::array<float, 2> dca3{1e10f, 1e10f};
    for (const auto& track1 : tracks) {
      if (!track1.has_mcParticle()) {
        continue;
      }
      const auto index1 = track1.globalIndex();
      int ncand = 0;
      if constexpr (usePID) {
        histos.fill(HIST("event/nsigmaDe"), track1.pt(), track1.tofNSigmaDe());
      }
      if (usePdg) {
        if (track1.template mcParticle_as<aod::McParticles>().pdgCode() != 1000010020) {
          continue;
        }
      }
      if constexpr (usePID) {
        if (abs(track1.tofNSigmaDe()) > maxNsigmaDe || track1.sign() < 0.f) {
          continue;
        }
        histos.fill(HIST("event/nsigmaDecut"), track1.pt(), track1.tofNSigmaDe());
      }
      if (!getTrackPar(track1).propagateParamToDCA(collPos,
                                                   magField * 10.f, &dca1, 100.)) {
        continue;
      }
      histos.fill(HIST("event/track1dcaxy"), dca1[0]);
      histos.fill(HIST("event/track1dcaz"), dca1[1]);

      for (const auto& track2 : tracks) {
        if (!track2.has_mcParticle()) {
          continue;
        }

        const auto index2 = track2.globalIndex();
        if (index1 == index2) {
          continue;
        }
        if constexpr (usePID) {
          histos.fill(HIST("event/nsigmaKa"), track2.pt(), track2.tofNSigmaKa());
        }
        if (usePdg) {
          if (track2.template mcParticle_as<aod::McParticles>().pdgCode() != -321) {
            continue;
          }
        }
        if constexpr (usePID) {
          if (abs(track2.tofNSigmaKa()) > maxNsigmaKa || track2.sign() > 0.f) {
            continue;
          }
          histos.fill(HIST("event/nsigmaKacut"), track2.pt(), track2.tofNSigmaKa());
        }
        if (!getTrackPar(track2).propagateParamToDCA(collPos,
                                                     magField * 10.f, &dca2, 100.)) {
          continue;
        }

        for (const auto& track3 : tracks) {
          if (!track3.has_mcParticle()) {
            continue;
          }

          const auto index3 = track3.globalIndex();
          if (index2 == index3) {
            continue;
          }
          if (index1 == index3) {
            continue;
          }
          if constexpr (usePID) {
            histos.fill(HIST("event/nsigmaPi"), track3.pt(), track3.tofNSigmaPi());
          }
          if (usePdg) {
            if (track3.template mcParticle_as<aod::McParticles>().pdgCode() != 211) {
              continue;
            }
          }
          if constexpr (usePID) {
            if (abs(track3.tofNSigmaPi()) > maxNsigmaPi || track3.sign() < 0.f) {
              continue;
            }
            histos.fill(HIST("event/nsigmaPicut"), track3.pt(), track3.tofNSigmaPi());
          }
          bool iscut = false;
          if (abs(dca1[0]) < minDca || abs(dca1[1]) < minDca) {
            iscut = true;
          }
          if (abs(dca1[0]) < minDcaDeuteron || abs(dca1[1]) < minDcaDeuteron) {
            iscut = true;
          }
          if (abs(dca1[0]) > maxDca || abs(dca1[1]) > maxDca) {
            iscut = true;
          }

          if (abs(dca2[0]) < minDca || abs(dca2[1]) < minDca) {
            iscut = true;
          }
          if (abs(dca2[0]) > maxDca || abs(dca2[1]) > maxDca) {
            iscut = true;
          }
          if (abs(dca2[0]) < minDcaPion || abs(dca2[1]) < minDcaPion) {
            iscut = true;
          }
          if (track2.pt() < minKaonPt) {
            iscut = true;
          }
          if (track3.pt() < minPionPt) {
            iscut = true;
          }
          if (!getTrackPar(track3).propagateParamToDCA(collPos,
                                                       magField * 10.f, &dca3, 100.)) {
            continue;
          }

          if (abs(dca3[0]) < minDca || abs(dca3[1]) < minDca) {
            iscut = true;
          }
          if (abs(dca3[0]) > maxDca || abs(dca3[1]) > maxDca) {
            iscut = true;
          }

          const auto mother1 = track1.template mcParticle_as<aod::McParticles>().template mothers_as<aod::McParticles>()[0];
          const auto mother2 = track2.template mcParticle_as<aod::McParticles>().template mothers_as<aod::McParticles>()[0];
          const auto mother3 = track3.template mcParticle_as<aod::McParticles>().template mothers_as<aod::McParticles>()[0];
          bool issig = true;
          if (mother1 != mother2) {
            issig = false;
          }
          if (mother1 != mother3) {
            issig = false;
          }

          auto pc1 = getTrackParCov(track1);
          auto pc2 = getTrackParCov(track2);
          auto pc3 = getTrackParCov(track3);
          if (pc1.getSigmaY2() * pc1.getSigmaZ2() - pc1.getSigmaZY() * pc1.getSigmaZY() < 0.) {
            Printf("Track 1 has issues");
            continue;
          }
          if (pc2.getSigmaY2() * pc2.getSigmaZ2() - pc2.getSigmaZY() * pc2.getSigmaZY() < 0.) {
            Printf("Track 2 has issues");
            continue;
          }
          if (pc3.getSigmaY2() * pc3.getSigmaZ2() - pc3.getSigmaZY() * pc3.getSigmaZY() < 0.) {
            Printf("Track 3 has issues");
            continue;
          }
          const int status = fitter.process(pc1, pc2, pc3);
          if (status == 0) {
            continue;
          }

          TLorentzVector v1{};
          v1.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi(), 1.8756129);

          TLorentzVector v2{};
          v2.SetPtEtaPhiM(track2.pt(), track2.eta(), track2.phi(), 0.493677);

          TLorentzVector v3{};
          v3.SetPtEtaPhiM(track3.pt(), track3.eta(), track3.phi(), 0.139570);
          v1 += v2;
          v1 += v3;
          if (v1.Pt() < minMomPt) {
            iscut = true;
          }

          // fitter.propagateTracksToVertex();
          const auto& secVtx = fitter.getPCACandidate();
          const float decay_radius = sqrt(secVtx[0] * secVtx[0] + secVtx[1] * secVtx[1] + secVtx[2] * secVtx[2]);
          if (decay_radius < minRadius) {
            iscut = true;
          }
          if (decay_radius > maxRadius) {
            iscut = true;
          }

          const float magMom = sqrt(v1.Px() * v1.Px() + v1.Py() * v1.Py() + v1.Pz() * v1.Pz());
          const float CPA = (v1.Px() * secVtx[0] + v1.Py() * secVtx[1] + v1.Pz() * secVtx[2]) / (decay_radius * magMom);
          if (abs(CPA) < minCpa) {
            iscut = true;
          }

          const float vx = mother1.vx();
          const float vy = mother1.vy();
          const float vz = mother1.vz();
          const float rmc = sqrt((secVtx[0] - vx) * (secVtx[0] - vx) + (secVtx[1] - vy) * (secVtx[1] - vy) + (secVtx[2] - vz) * (secVtx[2] - vz));
          ncand++;
          const float radius3xy = sqrt((track3.template mcParticle_as<aod::McParticles>().vx() - coll.mcCollision().posX()) * (track3.template mcParticle_as<aod::McParticles>().vx() - coll.mcCollision().posX()) +
                                       (track3.template mcParticle_as<aod::McParticles>().vy() - coll.mcCollision().posY()) * (track3.template mcParticle_as<aod::McParticles>().vy() - coll.mcCollision().posY()));

#define FillHistos(tag)                                                              \
  histos.fill(HIST(tag "/cpa"), CPA);                                                \
  histos.fill(HIST(tag "/invmass"), v1.M());                                         \
  histos.fill(HIST(tag "/invmassVsPt"), v1.Pt(), v1.M());                            \
  histos.fill(HIST(tag "/decayradius"), decay_radius);                               \
  histos.fill(HIST(tag "/decayradiusResoX"), secVtx[0] - vx);                        \
  histos.fill(HIST(tag "/decayradiusResoY"), secVtx[1] - vy);                        \
  histos.fill(HIST(tag "/decayradiusResoZ"), secVtx[2] - vz);                        \
  histos.fill(HIST(tag "/radius3xy"), radius3xy);                                    \
  histos.fill(HIST(tag "/decayradiusReso"), rmc);                                    \
  histos.fill(HIST(tag "/decaydca0"), TMath::Sqrt(fitter.getChi2AtPCACandidate(0))); \
  histos.fill(HIST(tag "/decaydca1"), TMath::Sqrt(fitter.getChi2AtPCACandidate(1))); \
  histos.fill(HIST(tag "/dcaxy1"), dca1[0]);                                         \
  histos.fill(HIST(tag "/dcaz1"), dca1[1]);                                          \
  histos.fill(HIST(tag "/dcaxy2"), dca2[0]);                                         \
  histos.fill(HIST(tag "/dcaz2"), dca2[1]);                                          \
  histos.fill(HIST(tag "/dcaxy3"), dca3[0]);                                         \
  histos.fill(HIST(tag "/dcaz3"), dca3[1]);                                          \
  histos.fill(HIST(tag "/dcaxy1xdcaxy2"), dca1[0] * dca2[0]);                        \
  histos.fill(HIST(tag "/dcaz1xdcaz2"), dca1[1] * dca2[1]);                          \
  histos.fill(HIST(tag "/dcaxy3xdcaxy2"), dca3[0] * dca2[0]);                        \
  histos.fill(HIST(tag "/dcaz3xdcaz2"), dca3[1] * dca2[1]);                          \
  histos.fill(HIST(tag "/pt1"), track1.pt());                                        \
  histos.fill(HIST(tag "/pt2"), track2.pt());                                        \
  histos.fill(HIST(tag "/pt3"), track3.pt());                                        \
  histos.fill(HIST(tag "/ptmom"), v1.Pt());                                          \
  histos.fill(HIST(tag "/pmom"), v1.P());

          if (issig) {
            FillHistos("signocut");
          } else {
            FillHistos("bkgnocut");
          }
          if (iscut) {
            if (issig) {
              FillHistos("sigcut");
            } else {
              FillHistos("bkgcut");
            }
            continue;
          }

          if (issig) {
            FillHistos("sig");
          } else {
            FillHistos("bkg");
          }
#undef FillHistos

          // fitterCasc.getTrack(1).getPxPyPzGlo(pvecbach);
        } // End loop on pions
      } // End loop on kaons
      histos.fill(HIST("event/candperdeuteron"), ncand);
    } // End loop on deuterons
  }

  void processWithPid(const soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels>::iterator& coll,
                      const o2::aod::McCollisions&,
                      const soa::Join<o2::aod::TracksIU, o2::aod::McTrackLabels, o2::aod::TracksExtra, o2::aod::TracksCov,
                                      aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullDe>& tracks,
                      const aod::McParticles& mcParticles)
  {
    fillHistograms<true>(coll, tracks, mcParticles);
  }
  PROCESS_SWITCH(Alice3CDeuteron, processWithPid, "With detector PID info", true);

  void processNoPid(const soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels>::iterator& coll,
                    const o2::aod::McCollisions&,
                    const soa::Join<o2::aod::TracksIU, o2::aod::McTrackLabels, o2::aod::TracksExtra, o2::aod::TracksCovIU>& tracks,
                    const aod::McParticles& mcParticles)
  {
    fillHistograms<false>(coll, tracks, mcParticles);
  }
  PROCESS_SWITCH(Alice3CDeuteron, processNoPid, "With detector PID info", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3CDeuteron>(cfgc)};
}
