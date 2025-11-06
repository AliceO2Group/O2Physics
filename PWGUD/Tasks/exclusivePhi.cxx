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
#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TLorentzVector.h"
#include <TString.h>

#include <iostream>
#include <vector>

using std::array;
using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// \brief Exclusive phi without PID
/// \author Simone Ragoni, Creighton
/// \author Anisa Khatun, Kansas University
/// \date 14/2/2024

struct ExclusivePhi {
  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  //_____________________________________________________________________________
  Double_t CosThetaHelicityFrame(TLorentzVector pionPositive,
                                 TLorentzVector pionNegative,
                                 TLorentzVector possibleRhoZero)
  {

    Double_t HalfSqrtSnn = 2680.;
    Double_t MassOfLead208 = 193.6823;
    Double_t MomentumBeam = TMath::Sqrt(HalfSqrtSnn * HalfSqrtSnn * 208 * 208 - MassOfLead208 * MassOfLead208);

    TLorentzVector pProjCM(0., 0., -MomentumBeam, HalfSqrtSnn * 208); // projectile
    TLorentzVector pTargCM(0., 0., MomentumBeam, HalfSqrtSnn * 208);  // target

    TVector3 beta = (-1. / possibleRhoZero.E()) * possibleRhoZero.Vect();
    TLorentzVector pPi1Dipion = pionPositive;
    TLorentzVector pPi2Dipion = pionNegative;
    TLorentzVector pProjDipion = pProjCM;
    TLorentzVector pTargDipion = pTargCM;

    pPi1Dipion.Boost(beta);
    pPi2Dipion.Boost(beta);
    pProjDipion.Boost(beta);
    pTargDipion.Boost(beta);

    TVector3 zaxis = (possibleRhoZero.Vect()).Unit();

    Double_t CosThetaHE = zaxis.Dot((pPi1Dipion.Vect()).Unit());
    return CosThetaHE;
  }
  //------------------------------------------------------------------------------------------------------
  Double_t PhiHelicityFrame(TLorentzVector muonPositive, TLorentzVector muonNegative, TLorentzVector possibleJPsi)
  {

    // Half of the energy per pair of the colliding nucleons.
    Double_t HalfSqrtSnn = 2680.;
    Double_t MassOfLead208 = 193.6823;
    Double_t MomentumBeam = TMath::Sqrt(HalfSqrtSnn * HalfSqrtSnn * 208 * 208 - MassOfLead208 * MassOfLead208);

    TLorentzVector pProjCM(0., 0., -MomentumBeam, HalfSqrtSnn * 208); // projectile
    TLorentzVector pTargCM(0., 0., MomentumBeam, HalfSqrtSnn * 208);  // target

    // Translate the dimuon parameters in the dimuon rest frame
    TVector3 beta = (-1. / possibleJPsi.E()) * possibleJPsi.Vect();
    TLorentzVector pMu1Dimu = muonPositive;
    TLorentzVector pMu2Dimu = muonNegative;
    TLorentzVector pProjDimu = pProjCM;
    TLorentzVector pTargDimu = pTargCM;
    pMu1Dimu.Boost(beta);
    pMu2Dimu.Boost(beta);
    pProjDimu.Boost(beta);
    pTargDimu.Boost(beta);

    // Axes
    TVector3 zaxis = (possibleJPsi.Vect()).Unit();
    TVector3 yaxis = ((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
    TVector3 xaxis = (yaxis.Cross(zaxis)).Unit();
    //
    // --- Calculation of the azimuthal angle (Helicity)
    //
    Double_t phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis), (pMu1Dimu.Vect()).Dot(xaxis));
    return phi;
  }
  //-----------------------------------------------------------------------------------------------------------------------
  void init(o2::framework::InitContext&)
  {
    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{12, 0., 12.}});
    TString SelectionCuts[12] = {"NoSelection", "Trackloop", "PVtracks", "|nsigmaka|<3", "|nsigmapi|>3", "|nsigmael|>3", "|nsigmamu|>3", "two tracks", "Phi-peak", "pt<0.2 GeV/c", "pt>0.2 GeV/c"};

    for (int i = 0; i < 12; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }

    registry.add("posx", "Vertex position in x", kTH1F, {{100, -0.5, 0.5}});
    registry.add("posy", "Vertex position in y", kTH1F, {{100, -0.5, 0.5}});
    registry.add("posz", "Vertex position in z", kTH1F, {{1000, -100., 100.}});

    registry.add("hTracks", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTracksITSonly", "N_{tracks ITS only}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hITSCluster", "N_{cluster}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hChi2ITSTrkSegment", "N_{cluster}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTPCCluster", "N_{cluster}", kTH1F, {{200, -0.5, 199.5}});
    registry.add("hTracksKaons", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hdEdx", "p vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon1", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon2", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon3", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon4", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon5", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon6", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon7", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon8", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon9", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});

    registry.add("hNsigEvsKa1", "NSigmaKa(t1) vs NSigmaKa (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, 0., 1000.}, {100, 0., 1000}});
    registry.add("hNsigEvsKa2", "NSigmaKa(t1) vs NSigmaKa (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, 0., 1000.}, {100, 0., 1000}});
    registry.add("hMomentum", "p_{#ka};#it{p_{trk}}, GeV/c;", kTH1F, {{100, 0., 3.}});
    registry.add("hClusterSizeAllTracks", "ClusterSizeAllTracks;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hClusterSizeMomentumCut", "ClusterSizeMomentumCut;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hClusterSizeOnlyIdentifiedKaons", "ClusterSizeOnlyIdentifiedKaons;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hClusterSizeOnlyITS", "ClusterSizeOnlyITS;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hEta1", "#eta_{#ka};#it{#eta_{trk}}, GeV/c;", kTH1F, {{100, -2., 2.}});
    registry.add("hEta2", "#eta_{#ka};#it{#eta_{trk}}, GeV/c;", kTH1F, {{100, -2., 2.}});
    registry.add("hPtPhi", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("hRapidityPhi", "Rapidity;#it{y_{KK}};", kTH1F, {{100, -2., 2.}});
    registry.add("hMassPhi", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("hMassPhiIncoherent", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});

    registry.add("hMassPtPhi", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});

    auto hSelectionCounter2 = registry.add<TH1>("hSelectionCounter2", "hSelectionCounter;;NEvents", HistType::kTH1I, {{12, 0., 12.}});
    TString SelectionCuts2[12] = {"NoSelection", "Trackloop", "PVtracks", "|nTPCCluster|<50", " track pt<0.180 GeV/c", "Kaon Band", "ITSCluster<6", "two tracks", "Phi-peak", "pt<0.2 GeV/c", "pt>0.2 GeV/c"};

    for (int i = 0; i < 12; i++) {
      hSelectionCounter2->GetXaxis()->SetBinLabel(i + 1, SelectionCuts2[i].Data());
    }

    registry.add("hMassPhiWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("hMassPtPhiWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("hMassPhiWrongMomentumWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1D, {{400, 0., 4.}});

    // Phi peak region
    registry.add("PHI/hPtPhiWithoutPID", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PHI/hMassVsPt", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("PHI/hRapidityPhiWithoutPID", "Rapidity;#it{y_{KK}};", kTH1F, {{100, -2., 2.}});
    registry.add("PHI/hCosThetaPhiWithoutPID", "CosTheta;cos(#theta);", kTH1F, {{100, -2., 2.}});
    registry.add("PHI/hPhiPhiWithoutPID", "Phi;#varphi;", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("PHI/hPtKaonVsKaon", "Pt1 vs Pt2;p_{T};p_{T};", kTH2F, {{100, 0., 3.}, {100, 0., 3.}});
    registry.add("PHI/hCostheta_Phi", "Phi vs Costheta;#it{#phi};#it{Cos#Theta};", kTH2F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}, {100, -2, 2}});
    registry.add("PHI/hMassPhiWithoutPIDPionHypothesis", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1D, {{400, 0., 4.}});
    registry.add("PHI/hPtPhiWithoutPIDPionHypothesis", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});

    // DIfferent phi topologies, 2 identified kaons, 1 identified kaon + 1 ITS track with correct selections and 4 ITS clusters
    registry.add("KaonBandPHI/hMassPtPhiIdentifiedKaons", "Raw Inv.M;#it{m_{KK}} (GeV/c^{2});Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("KaonBandPHI/hMassPhiIdentifiedKaons", "Raw Inv.M;#it{M_{KK}} (GeV/c^{2});", kTH1F, {{400, 0., 4.}});
    registry.add("KaonBandPHI/hPtPhiIdentifiedKaons", "Pt;#it{p_{t}} (GeV/c);", kTH1F, {{400, 0., 4.}});
    registry.add("KaonBandPHI/hMassPtPhiIdentifiedKaonAndITSkaon", "Raw Inv.M;#it{m_{KK}} (GeV/c^{2});Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("KaonBandPHI/hMassPhiIdentifiedKaonAndITSkaon", "Raw Inv.M;#it{m_{KK}} (GeV/c^{2});", kTH1F, {{400, 0., 4.}});
    registry.add("KaonBandPHI/hPtPhiIdentifiedKaonAndITSkaon", "Pt;#it{p_{t}} (GeV/c);", kTH1F, {{400, 0., 4.}});

    registry.add("PHI/hMassLike", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hMassUnlike", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hlikePt", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PHI/hUnlikePt", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PHI/hCoherentPhiWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hInCoherentPhiWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hCoherentMassLike", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hInCoherentMassLike", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});

    // High Mass region
    registry.add("PHIHIGH/hPtPhiWithoutPID1", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PHIHIGH/hMassVsPt1", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("PHIHIGH/hRapidityPhiWithoutPID1", "Rapidity;#it{y_{KK}};", kTH1F, {{100, -2., 2.}});
    registry.add("PHIHIGH/hCosThetaPhiWithoutPID1", "CosTheta;cos(#theta);", kTH1F, {{100, -2., 2.}});
    registry.add("PHIHIGH/hPhiPhiWithoutPID1", "Phi;#varphi;", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("PHIHIGH/hPtKaonVsKaon1", "Pt1 vs Pt2;p_{T};p_{T};", kTH2F, {{100, 0., 3.}, {100, 0., 3.}});
    registry.add("PHIHIGH/hCostheta_Phi1", "Phi vs Costheta;#it{#phi};#it{Cos#Theta};", kTH2F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}, {100, -2, 2}});
    registry.add("PHIHIGH/hMassPhiWithoutPIDPionHypothesis1", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1D, {{400, 0., 4.}});
    registry.add("PHIHIGH/hPtPhiWithoutPIDPionHypothesis1", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});

    registry.add("PHIHIGH/hMassLike1", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hMassUnlike1", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hlikePt1", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hUnlikePt1", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hCoherentPhiWithoutPID1", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hInCoherentPhiWithoutPID1", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hCoherentMassLike1", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hInCoherentMassLike1", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});

    // Low Mass region
    registry.add("PHILOW/hPtPhiWithoutPID2", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PHILOW/hMassVsPt2", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("PHILOW/hRapidityPhiWithoutPID2", "Rapidity;#it{y_{KK}};", kTH1F, {{100, -2., 2.}});
    registry.add("PHILOW/hCosThetaPhiWithoutPID2", "CosTheta;cos(#theta);", kTH1F, {{100, -2., 2.}});
    registry.add("PHILOW/hPhiPhiWithoutPID2", "Phi;#varphi;", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("PHILOW/hPtKaonVsKaon2", "Pt1 vs Pt2;p_{T};p_{T};", kTH2F, {{100, 0., 3.}, {100, 0., 3.}});
    registry.add("PHILOW/hCostheta_Phi2", "Phi vs Costheta;#it{#phi};#it{Cos#Theta};", kTH2F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}, {100, -2, 2}});
    registry.add("PHILOW/hMassPhiWithoutPIDPionHypothesis2", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1D, {{400, 0., 4.}});
    registry.add("PHILOW/hPtPhiWithoutPIDPionHypothesis2", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});

    registry.add("PHILOW/hMassLike2", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hMassUnlike2", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hlikePt2", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hUnlikePt2", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hCoherentPhiWithoutPID2", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hInCoherentPhiWithoutPID2", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hCoherentMassLike2", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hInCoherentMassLike2", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
  }

  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  // using UDCollisions = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>;
  //__________________________________________________________________________
  //  Main process
  void process(UDCollisions::iterator const& collision, udtracksfull const& tracks)
  {
    registry.fill(HIST("hSelectionCounter"), 0);

    registry.fill(HIST("posx"), collision.posX());
    registry.fill(HIST("posy"), collision.posY());
    registry.fill(HIST("posz"), collision.posZ());

    TLorentzVector phi, phiWithoutPID, phiWithKaonPID, phiWrongMomentaWithoutPID, phiWithoutPIDPionHypothesis; // lorentz vectors of tracks and the mother

    // ===================================
    // Task for phi WITH PID FROM TPC
    // ===================================
    // Here the tracks are allowed to have
    // at most 30 TPC clusters
    // TPC PID is then possible, although
    // not truly reliable
    // NB: from MC simulations coherent
    // phi SHOULD NOT leave clusters
    // in the TPC at all
    std::vector<TLorentzVector> onlyKaonTracks;
    std::vector<float> onlyKaonSigma;
    std::vector<decltype(tracks.begin())> rawKaonTracks;

    for (auto trk : tracks) {
      registry.fill(HIST("hSelectionCounter"), 1);
      if (!trk.isPVContributor()) {
        continue;
      }
      registry.fill(HIST("hSelectionCounter"), 2);

      int NFindable = trk.tpcNClsFindable();
      int NMinusFound = trk.tpcNClsFindableMinusFound();
      int NCluster = NFindable - NMinusFound;

      if (NCluster < 50) {
        continue;
      }
      registry.fill(HIST("hTPCCluster"), NCluster);
      registry.fill(HIST("hITSCluster"), trk.itsNCls());
      registry.fill(HIST("hChi2ITSTrkSegment"), trk.itsChi2NCl());

      double dEdx = trk.tpcSignal();
      registry.fill(HIST("hdEdx"), trk.tpcInnerParam() / trk.sign(), dEdx);

      TLorentzVector kaon;
      kaon.SetXYZM(trk.px(), trk.py(), trk.pz(), o2::constants::physics::MassElectron);
      // auto nSigmaKa = trk.tpcNSigmaKa();
      auto nSigmaEl = trk.tpcNSigmaEl();
      // auto nSigmaPi = trk.tpcNSigmaPi();
      // auto nSigmaMu = trk.tpcNSigmaMu();

      // if((nSigmaKa * nSigmaKa + nSigmaKa * nSigmaKa) > 9.) continue;
      /*if (fabs(nSigmaKa) > 3.)
        continue;*/
      registry.fill(HIST("hSelectionCounter"), 3);
      registry.fill(HIST("hdEdxKaon1"), trk.tpcInnerParam() / trk.sign(), dEdx);
      // if((nSigmaEl * nSigmaEl + nSigmaEl * nSigmaEl) < 9.) continue;
      // if (fabs(nSigmaPi) < 3.)continue;
      registry.fill(HIST("hSelectionCounter"), 4);
      registry.fill(HIST("hdEdxKaon2"), trk.tpcInnerParam() / trk.sign(), dEdx);
      // if((nSigmaPi * nSigmaPi + nSigmaPi * nSigmaPi) < 9.) continue;
      if (fabs(nSigmaEl) < 3.)
        continue;
      registry.fill(HIST("hSelectionCounter"), 5);
      registry.fill(HIST("hdEdxKaon3"), trk.tpcInnerParam() / trk.sign(), dEdx);
      // if((nSigmaMu * nSigmaMu + nSigmaMu * nSigmaMu) < 9.) continue;
      // if (fabs(nSigmaMu) < 3.)
      // continue;
      registry.fill(HIST("hSelectionCounter"), 6);
      registry.fill(HIST("hdEdxKaon4"), trk.tpcInnerParam() / trk.sign(), dEdx);

      onlyKaonTracks.push_back(kaon);
      onlyKaonSigma.push_back(nSigmaEl);
      rawKaonTracks.push_back(trk);

    } // trk loop

    if (onlyKaonTracks.size() == 2) {
      registry.fill(HIST("hSelectionCounter"), 7);

      for (auto kaon : onlyKaonTracks) {
        phi += kaon;
      }

      registry.fill(HIST("hdEdxKaon"), rawKaonTracks[0].tpcInnerParam() / rawKaonTracks[0].sign(), rawKaonTracks[1].tpcSignal());
      registry.fill(HIST("hNsigEvsKa1"), rawKaonTracks[0].tpcSignal(), rawKaonTracks[1].tpcSignal());
      registry.fill(HIST("hMassPtPhi"), phi.M(), phi.Pt());
      registry.fill(HIST("hRapidityPhi"), phi.Rapidity());

      if (rawKaonTracks[0].sign() != rawKaonTracks[1].sign()) {
        if ((phi.M() > 0.98) && (phi.M() < 1.05)) {
          registry.fill(HIST("hSelectionCounter"), 8);
          registry.fill(HIST("hPtPhi"), phi.Pt());

          if (phi.Pt() < 0.2) {
            registry.fill(HIST("hMassPhi"), phi.M());
            registry.fill(HIST("hSelectionCounter"), 9);
          }

          if (phi.Pt() > 0.2) {
            registry.fill(HIST("hMassPhiIncoherent"), phi.M());
            registry.fill(HIST("hSelectionCounter"), 10);
          }
        }
      }
    }

    // ===================================
    // Task for phi WITHOUT PID
    // ===================================

    // ====================================
    // Selections for events to be stored
    // ------------------------------------
    // - PV: only PV contributors
    // - only two PV contributors
    // - both PV have t.hasITS() ON
    // - only ITS tracks
    //_____________________________________
    // Create kaons WITHOUT PID
    std::vector<decltype(tracks.begin())> onlyTwoTracks;
    std::vector<decltype(tracks.begin())> onlyKaonBandPID;
    std::vector<decltype(tracks.begin())> onlyITS;
    std::vector<TLorentzVector> allTracksAreKaons;
    std::vector<TLorentzVector> allTracksArePions;
    std::vector<TLorentzVector> allTracksAreKaonsWrongMomentum;
    std::vector<TLorentzVector> allTracksAreKaonsBandPID;
    std::vector<TLorentzVector> allTracksAreITSonlyAndFourITSclusters;

    int counter = 0;
    for (auto t : tracks) {
      registry.fill(HIST("hSelectionCounter2"), 0);
      if (!t.isPVContributor()) {
        continue;
      }
      registry.fill(HIST("hSelectionCounter2"), 1);
      registry.fill(HIST("hTracks"), t.size());

      double momentum = TMath::Sqrt(t.px() * t.px() + t.py() * t.py() + t.pz() * t.pz());
      double dEdx = t.tpcSignal();

      int clusterSize[7];
      double averageClusterSize = 0.;
      for (int i = 0; i < 7; i++) { // info stored in 4 bits
        clusterSize[i] = (((1 << 4) - 1) & (t.itsClusterSizes() >> 4 * i));
        averageClusterSize += static_cast<double>(clusterSize[i]);
      }
      averageClusterSize /= 7.;
      registry.fill(HIST("hClusterSizeAllTracks"), averageClusterSize);

      int NFindable = t.tpcNClsFindable();
      int NMinusFound = t.tpcNClsFindableMinusFound();
      int NCluster = NFindable - NMinusFound;
      // registry.fill(HIST("hTPCCluster"), NCluster);

      if (NCluster > 50) {
        continue;
      }
      registry.fill(HIST("hSelectionCounter2"), 3);
      registry.fill(HIST("hdEdxKaon5"), t.tpcInnerParam() / t.sign(), dEdx);

      if (t.pt() > 0.180) {
        continue;
      }

      registry.fill(HIST("hSelectionCounter2"), 4);
      registry.fill(HIST("hMomentum"), momentum);
      registry.fill(HIST("hdEdxKaon6"), t.tpcInnerParam() / t.sign(), dEdx);
      registry.fill(HIST("hClusterSizeMomentumCut"), averageClusterSize);

      onlyTwoTracks.push_back(t);

      TLorentzVector a;
      a.SetXYZM(t.px(), t.py(), t.pz(), o2::constants::physics::MassKaonCharged);
      TLorentzVector b;
      b.SetXYZM(t.px(), t.py(), t.pz(), o2::constants::physics::MassElectron);
      TLorentzVector a2;
      a2.SetXYZM(-1. * t.px(), -1. * t.py(), -1. * t.pz(), o2::constants::physics::MassKaonCharged);

      allTracksAreKaons.push_back(a);
      allTracksArePions.push_back(b);

      if (counter < 1) {
        allTracksAreKaonsWrongMomentum.push_back(a);
      } else {
        allTracksAreKaonsWrongMomentum.push_back(a2);
      }
      counter += 1;

      bool kaonBand = false;
      if ((momentum > 0.180) && (momentum < 0.220) && (dEdx > 300)) {
        kaonBand = true;
      } else if ((momentum > 0.220) && (momentum < 0.300) && (dEdx > 180)) {
        kaonBand = true;
      } else if ((momentum > 0.300) && (momentum < 0.500) && (dEdx > 110)) {
        kaonBand = true;
      }

      if (kaonBand == true) {
        registry.fill(HIST("hSelectionCounter2"), 5);
        allTracksAreKaonsBandPID.push_back(a);
        onlyKaonBandPID.push_back(t);
        registry.fill(HIST("hdEdxKaon7"), t.tpcInnerParam() / t.sign(), dEdx);
        registry.fill(HIST("hClusterSizeOnlyIdentifiedKaons"), averageClusterSize);
      }

      if (NFindable < 1 && t.itsNCls() < 7) {
        // if((NCluster==0) && (t.itsNCls() == 4)){
        allTracksAreITSonlyAndFourITSclusters.push_back(a);
        onlyITS.push_back(t);
        registry.fill(HIST("hdEdxKaon8"), t.tpcInnerParam() / t.sign(), dEdx);
        registry.fill(HIST("hSelectionCounter2"), 6);
      }

    } // track loop

    //_____________________________________
    // Creating phis and saving all the information
    // in the case that there are ONLY 2 PV

    // if ((collision.posZ() < -10) || (collision.posZ() > 10)) {
    if (allTracksAreKaons.size() == 2) {
      registry.fill(HIST("hSelectionCounter2"), 7);
      for (auto kaon : allTracksAreKaons) {
        phiWithoutPID += kaon;
      }
      registry.fill(HIST("hTracksKaons"), allTracksAreKaons.size());
      // kaon mass hypothesis with wrong momentum for one track
      for (auto kaon : allTracksAreKaonsWrongMomentum) {
        phiWrongMomentaWithoutPID += kaon;
      }
      // pion mass hypothesis
      for (auto pion : allTracksArePions) {
        phiWithoutPIDPionHypothesis += pion;
      }

      registry.fill(HIST("hNsigEvsKa2"), onlyTwoTracks[0].tpcSignal(), onlyTwoTracks[1].tpcSignal());
      registry.fill(HIST("hEta1"), allTracksAreKaons[0].Eta());
      registry.fill(HIST("hEta2"), allTracksAreKaons[1].Eta());

      auto costhetaPhi = CosThetaHelicityFrame(allTracksAreKaons[0], allTracksAreKaons[1], phiWithoutPID);
      auto phiPhi = 1. * TMath::Pi() + PhiHelicityFrame(allTracksAreKaons[0], allTracksAreKaons[1], phiWithoutPID);

      // All invariant mass region
      registry.fill(HIST("hMassPhiWithoutPID"), phiWithoutPID.M());
      registry.fill(HIST("hMassPtPhiWithoutPID"), phiWithoutPID.M(), phiWithoutPID.Pt());
      registry.fill(HIST("hMassPhiWrongMomentumWithoutPID"), phiWrongMomentaWithoutPID.M());

      // Phi peak region
      if ((phiWithoutPID.M() > 0.98) && (phiWithoutPID.M() < 1.05)) {
        registry.fill(HIST("PHI/hPtPhiWithoutPID"), phiWithoutPID.Pt());
        registry.fill(HIST("PHI/hMassVsPt"), phiWithoutPID.M(), phiWithoutPID.Pt());
        registry.fill(HIST("PHI/hRapidityPhiWithoutPID"), phiWithoutPID.Rapidity());
        registry.fill(HIST("PHI/hPtKaonVsKaon"), allTracksAreKaons[0].Pt(), allTracksAreKaons[1].Pt());

        registry.fill(HIST("PHI/hPtPhiWithoutPIDPionHypothesis"), phiWithoutPIDPionHypothesis.Pt());
        registry.fill(HIST("PHI/hMassPhiWithoutPIDPionHypothesis"), phiWithoutPIDPionHypothesis.M());

        // unlike-sign
        if (onlyTwoTracks[0].sign() != onlyTwoTracks[1].sign()) {
          registry.fill(HIST("hSelectionCounter2"), 8);
          registry.fill(HIST("PHI/hCosThetaPhiWithoutPID"), costhetaPhi);
          registry.fill(HIST("PHI/hPhiPhiWithoutPID"), phiPhi);
          registry.fill(HIST("PHI/hCostheta_Phi"), phiPhi, costhetaPhi);
          registry.fill(HIST("PHI/hUnlikePt"), phiWithoutPID.Pt());
          registry.fill(HIST("PHI/hMassUnlike"), phiWithoutPID.M());
          if (phiWithoutPID.Pt() < 0.2) {
            registry.fill(HIST("hSelectionCounter2"), 9);
            registry.fill(HIST("PHI/hCoherentPhiWithoutPID"), phiWithoutPID.M());
          }
          if (phiWithoutPID.Pt() > 0.2) {
            registry.fill(HIST("hSelectionCounter2"), 10);
            registry.fill(HIST("PHI/hInCoherentPhiWithoutPID"), phiWithoutPID.M());
          }
        }

        // Likesign quantities
        if (onlyTwoTracks[0].sign() == onlyTwoTracks[1].sign()) {
          registry.fill(HIST("PHI/hMassLike"), phiWithoutPID.M());
          registry.fill(HIST("PHI/hlikePt"), phiWithoutPID.Pt());

          if (phiWithoutPID.Pt() < 0.2) {
            registry.fill(HIST("PHI/hCoherentMassLike"), phiWithoutPID.M());
          }
          if (phiWithoutPID.Pt() > 0.2) {
            registry.fill(HIST("PHI/hInCoherentMassLike"), phiWithoutPID.M());
          }
        }
      } // Mass cut

      // Side band above phi mass region
      if (phiWithoutPID.M() > 1.05) {

        registry.fill(HIST("PHIHIGH/hPtPhiWithoutPID1"), phiWithoutPID.Pt());
        registry.fill(HIST("PHIHIGH/hMassVsPt1"), phiWithoutPID.M(), phiWithoutPID.Pt());
        registry.fill(HIST("PHIHIGH/hRapidityPhiWithoutPID1"), phiWithoutPID.Rapidity());

        registry.fill(HIST("PHIHIGH/hPtKaonVsKaon1"), allTracksAreKaons[0].Pt(), allTracksAreKaons[1].Pt());
        registry.fill(HIST("PHIHIGH/hPtPhiWithoutPIDPionHypothesis1"), phiWithoutPIDPionHypothesis.Pt());
        registry.fill(HIST("PHIHIGH/hMassPhiWithoutPIDPionHypothesis1"), phiWithoutPIDPionHypothesis.M());

        // unlike-sign
        if (onlyTwoTracks[0].sign() != onlyTwoTracks[1].sign()) {
          registry.fill(HIST("PHIHIGH/hCosThetaPhiWithoutPID1"), costhetaPhi);
          registry.fill(HIST("PHIHIGH/hPhiPhiWithoutPID1"), phiPhi);
          registry.fill(HIST("PHIHIGH/hCostheta_Phi1"), phiPhi, costhetaPhi);

          registry.fill(HIST("PHIHIGH/hUnlikePt1"), phiWithoutPID.Pt());
          registry.fill(HIST("PHIHIGH/hMassUnlike1"), phiWithoutPID.M());
          if (phiWithoutPID.Pt() < 0.2) {
            registry.fill(HIST("PHIHIGH/hCoherentPhiWithoutPID1"), phiWithoutPID.M());
          }
          if (phiWithoutPID.Pt() > 0.2) {
            registry.fill(HIST("PHIHIGH/hInCoherentPhiWithoutPID1"), phiWithoutPID.M());
          }
        }

        // Like sign
        if (onlyTwoTracks[0].sign() == onlyTwoTracks[1].sign()) {
          // Likesign quantities
          registry.fill(HIST("PHIHIGH/hMassLike1"), phiWithoutPID.M());
          registry.fill(HIST("PHIHIGH/hlikePt1"), phiWithoutPID.Pt());

          if (phiWithoutPID.Pt() < 0.2) {
            registry.fill(HIST("PHIHIGH/hCoherentMassLike1"), phiWithoutPID.M());
          }
          if (phiWithoutPID.Pt() > 0.2) {
            registry.fill(HIST("PHIHIGH/hInCoherentMassLike1"), phiWithoutPID.M());
          }
        }
      }

      // Side band below phi mass region
      if (phiWithoutPID.M() < 0.98) {
        registry.fill(HIST("PHILOW/hPtPhiWithoutPID2"), phiWithoutPID.Pt());
        registry.fill(HIST("PHILOW/hMassVsPt2"), phiWithoutPID.M(), phiWithoutPID.Pt());
        registry.fill(HIST("PHILOW/hRapidityPhiWithoutPID2"), phiWithoutPID.Rapidity());

        registry.fill(HIST("PHILOW/hPtKaonVsKaon2"), allTracksAreKaons[0].Pt(), allTracksAreKaons[1].Pt());

        registry.fill(HIST("PHILOW/hPtPhiWithoutPIDPionHypothesis2"), phiWithoutPIDPionHypothesis.Pt());
        registry.fill(HIST("PHILOW/hMassPhiWithoutPIDPionHypothesis2"), phiWithoutPIDPionHypothesis.M());

        // unlike-sign
        registry.fill(HIST("PHILOW/hCosThetaPhiWithoutPID2"), costhetaPhi);
        registry.fill(HIST("PHILOW/hPhiPhiWithoutPID2"), phiPhi);
        registry.fill(HIST("PHILOW/hCostheta_Phi2"), phiPhi, costhetaPhi);
        registry.fill(HIST("PHILOW/hUnlikePt2"), phiWithoutPID.Pt());
        registry.fill(HIST("PHILOW/hMassUnlike2"), phiWithoutPID.M());
        if (phiWithoutPID.Pt() < 0.2) {
          registry.fill(HIST("PHILOW/hCoherentPhiWithoutPID2"), phiWithoutPID.M());
        }
        if (phiWithoutPID.Pt() > 0.2) {
          registry.fill(HIST("PHILOW/hInCoherentPhiWithoutPID2"), phiWithoutPID.M());
        }
        // like-sign
        if (onlyTwoTracks[0].sign() == onlyTwoTracks[1].sign()) {
          // Likesign quantities
          registry.fill(HIST("PHI/hMassLike2"), phiWithoutPID.M());
          registry.fill(HIST("PHI/hlikePt2"), phiWithoutPID.Pt());

          if (phiWithoutPID.Pt() < 0.2) {
            registry.fill(HIST("PHI/hCoherentMassLike2"), phiWithoutPID.M());
          }

          if (phiWithoutPID.Pt() > 0.2) {
            registry.fill(HIST("PHI/hInCoherentMassLike2"), phiWithoutPID.M());
          }
        }
      }
    } // end of two tracks only loop
    // }   // vertex cut

    if (allTracksAreKaonsBandPID.size() == 2) {

      TLorentzVector reallyPhi;
      for (auto kaon : allTracksAreKaonsBandPID) {
        reallyPhi += kaon;
      }

      registry.fill(HIST("KaonBandPHI/hMassPtPhiIdentifiedKaons"), reallyPhi.M(), reallyPhi.Pt());
      if (reallyPhi.Pt() < 0.2) {
        registry.fill(HIST("KanonBandPHI/hMassPhiIdentifiedKaons"), reallyPhi.M());
      }
      registry.fill(HIST("KaonBandPHI/hPtPhiIdentifiedKaons"), reallyPhi.Pt());
    }

    if (allTracksAreKaonsBandPID.size() == 1) {

      double momentum = TMath::Sqrt(onlyKaonBandPID[0].px() * onlyKaonBandPID[0].px() + onlyKaonBandPID[0].py() * onlyKaonBandPID[0].py() + onlyKaonBandPID[0].pz() * onlyKaonBandPID[0].pz());
      double dEdx = onlyKaonBandPID[0].tpcSignal();
      registry.fill(HIST("hdEdxKaon9"), momentum, dEdx);

      // auto ksize = allTracksAreITSonlyAndFourITSclusters.size();
      registry.fill(HIST("hTracksITSonly"), allTracksAreITSonlyAndFourITSclusters.size());

      for (std::size_t kaon = 0; kaon < allTracksAreITSonlyAndFourITSclusters.size(); kaon++) {

        int clusterSize[7];
        double averageClusterSize = 0.;
        for (int i = 0; i < 7; i++) { // info stored in 4 bits
          clusterSize[i] = (((1 << 4) - 1) & (onlyITS[kaon].itsClusterSizes() >> 4 * i));
          averageClusterSize += static_cast<double>(clusterSize[i]);
        }
        averageClusterSize /= 7.;
        registry.fill(HIST("hClusterSizeOnlyITS"), averageClusterSize);

        TLorentzVector reallyPhi;
        reallyPhi += allTracksAreKaonsBandPID[0];
        reallyPhi += allTracksAreITSonlyAndFourITSclusters[kaon];
        registry.fill(HIST("KaonBandPHI/hMassPtPhiIdentifiedKaonAndITSkaon"), reallyPhi.M(), reallyPhi.Pt());
        registry.fill(HIST("KaonBandPHI/hPtPhiIdentifiedKaonAndITSkaon"), reallyPhi.Pt());
        if (reallyPhi.Pt() < 0.2) {
          registry.fill(HIST("KaonBandPHI/hMassPhiIdentifiedKaonAndITSkaon"), reallyPhi.M());
        }
      }
    } // Kaon Band
  } // end of process

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ExclusivePhi>(cfgc)};
}
