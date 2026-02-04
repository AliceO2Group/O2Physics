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

#include "FastTracker.h"
#include "TrackUtilities.h"

#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <TAxis.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TEfficiency.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLorentzVector.h>

void drawFastTracker(float magneticField = 5.f, // in units of kGauss
                     const int nch = 100,       // number of charged particles per unit rapidity
                     const int pdg = 211)       // PDG code of the particle to track
{
  TDatabasePDG* db = TDatabasePDG::Instance();
  TParticlePDG* p = 0;
  p = db->GetParticle(pdg);
  if (!p) {
    LOG(fatal) << "Particle with PDG code " << pdg << " not found in TDatabasePDG";
    return;
  }
  const float mass = p->Mass();      // particle mass in GeV/c^2
  const float q = p->Charge() / 3.0; // charge in e

  o2::parameters::GRPMagField grpmag;
  grpmag.setFieldUniformity(true);
  grpmag.setL3Current(30000.f * (magneticField / 5.0f));
  auto field = grpmag.getNominalL3Field();
  o2::base::Propagator::initFieldFromGRP(&grpmag);

  fair::Logger::SetVerbosity(fair::Verbosity::verylow);
  o2::fastsim::FastTracker fastTracker = o2::fastsim::FastTracker();
  if (0) {
    fastTracker.SetApplyEffCorrection(false);
    Double_t x0IB = 0.001;
    Double_t x0OB = 0.01;
    Double_t xrhoIB = 2.3292e-02; // 100 mum Si
    Double_t xrhoOB = 2.3292e-01; // 1000 mum Si
    Double_t resRPhiIB = 0.00025;
    Double_t resZIB = 0.00025;
    Double_t resRPhiOB = 0.00100;
    Double_t resZOB = 0.00100;
    Double_t eff = 0.98;
    // fastTracker.AddLayer("vertex", 0.0, 250, 0, 0);                // dummy vertex for matrix calculation
    fastTracker.AddLayer("bpipe0", 0.48, 250, 0.00042, 2.772e-02); // 150 mum Be
    fastTracker.AddLayer("B00", 0.50, 250, x0IB, xrhoIB, resRPhiIB, resZIB, eff, 1);
    fastTracker.AddLayer("B01", 1.20, 250, x0IB, xrhoIB, resRPhiIB, resZIB, eff, 1);
    fastTracker.AddLayer("B02", 2.50, 250, x0IB, xrhoIB, resRPhiIB, resZIB, eff, 1);
    fastTracker.AddLayer("bpipe1", 3.7, 250, 0.0014, 9.24e-02); // 500 mum Be
    fastTracker.AddLayer("B03", 3.75, 250, x0OB, xrhoOB, resRPhiOB, resZOB, eff, 1);
    fastTracker.AddLayer("B04", 7.00, 250, x0OB, xrhoOB, resRPhiOB, resZOB, eff, 1);
    fastTracker.AddLayer("B05", 12.0, 250, x0OB, xrhoOB, resRPhiOB, resZOB, eff, 1);
    fastTracker.AddLayer("B06", 20.0, 250, x0OB, xrhoOB, resRPhiOB, resZOB, eff, 1);
    fastTracker.AddLayer("B07", 30.0, 250, x0OB, xrhoOB, resRPhiOB, resZOB, eff, 1);
    fastTracker.AddLayer("B08", 45.0, 250, x0OB, xrhoOB, resRPhiOB, resZOB, eff, 1);
    fastTracker.AddLayer("B09", 60.0, 250, x0OB, xrhoOB, resRPhiOB, resZOB, eff, 1);
    fastTracker.AddLayer("B10", 80.0, 250, x0OB, xrhoOB, resRPhiOB, resZOB, eff, 1);
    fastTracker.AddLayer("B11", 100., 250, x0OB, xrhoOB, resRPhiOB, resZOB, eff, 1);
  } else {
    std::vector<float> pixelRes{0.025, 0.025, 0.01, 0.01};
    fastTracker.AddSiliconALICE3v4(pixelRes);
  }

  fastTracker.Print();
  fastTracker.SetMagneticField(magneticField);

  TAxis ptBinning(1000, 0., 10);
  TGraph* gPt = new TGraph();
  gPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  gPt->GetYaxis()->SetTitle("Efficiency");
  TEfficiency* hEfficiency = new TEfficiency("hEfficiency", ";#it{p}_{T} (GeV/c);Efficiency", ptBinning.GetNbins(), ptBinning.GetBinLowEdge(1), ptBinning.GetBinUpEdge(ptBinning.GetNbins()));
  TH1F* hFastTrackerQA = new TH1F("hFastTrackerQA", ";#it{p}_{T} (GeV/c;Tracking code", ptBinning.GetNbins(), ptBinning.GetBinLowEdge(1), ptBinning.GetBinUpEdge(ptBinning.GetNbins()));

  TLorentzVector tlv;

  // Now we compute the efficiency for a range of pt values and plot the results
  o2::track::TrackParCov trkIn;
  o2::track::TrackParCov trkOut;
  for (int i = 1; i <= ptBinning.GetNbins(); ++i) {
    const float pt = ptBinning.GetBinCenter(i);
    tlv.SetPtEtaPhiM(pt, 0., 0., mass);
    o2::upgrade::convertTLorentzVectorToO2Track(q, tlv, {0.0, 0., 0.}, trkIn);

    bool prop = o2::base::Propagator::Instance()->propagateToR(trkIn, 0.4f, true, o2::base::Propagator::MAX_SIN_PHI, o2::base::Propagator::MAX_STEP, o2::base::Propagator::MatCorrType::USEMatCorrNONE);
    if (!prop) {
      LOG(info) << "Initial propagation to R=0 failed for pt = " << pt;
      continue;
    }

    o2::dataformats::VertexBase vertex{{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f, 0.f, 0.f, 0.f}};
    // Now we propagate to the DCA of the vertex
    if (!trkIn.propagateToDCA(vertex, magneticField)) {
      LOG(info) << "Propagation to DCA failed for pt = " << pt;
      continue;
    }

    for (int trial = 0; trial < 200; trial++) {
      hEfficiency->Fill(fastTracker.FastTrack(trkIn, trkOut, nch) > 0, pt);
    }
    int status = 4;
    status = fastTracker.FastTrack(trkIn, trkOut, nch);
    hFastTrackerQA->Fill(pt, status);
    if (status < 0) {
      LOG(debug) << " --- fatSolve: FastTrack failed with status " << status << " --- ";
      // tlv.Print();
      gPt->AddPoint(pt, 0);
      continue;
    }
    // define the efficiency
    float eff = 1.;
    for (size_t l = 1; l < fastTracker.GetNLayers(); ++l) {
      if (fastTracker.IsLayerInert(l)) {
        continue; // skip inert layers
      }
      float igoodhit = 0.f;
      igoodhit = fastTracker.GetGoodHitProb(l);
      if (igoodhit <= 0.) {
        continue;
      }
      eff *= igoodhit;
    }
    gPt->AddPoint(pt, eff);
  }

  TLatex latex;
  latex.SetTextFont(42);
  latex.SetNDC();
  latex.SetTextSize(0.037);

  new TCanvas("gPt");
  gPt->Draw("ALP");
  hEfficiency->Draw("SAME");
  latex.DrawLatex(0.1, 0.93, Form("#LT #frac{dN_{ch}}{d#eta} #GT = %i", nch));

  new TCanvas("hEfficiency");
  hEfficiency->Draw("");
  latex.DrawLatex(0.1, 0.93, Form("#LT #frac{dN_{ch}}{d#eta} #GT = %i", nch));
  new TCanvas("hFastTrackerQA");
  hFastTrackerQA->Print();
  hFastTrackerQA->Draw("");
  latex.DrawLatex(0.1, 0.93, Form("#LT #frac{dN_{ch}}{d#eta} #GT = %i", nch));
}
