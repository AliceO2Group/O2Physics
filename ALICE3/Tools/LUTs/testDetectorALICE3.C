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

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <TRandom.h>
#include <TObjArray.h>
#include "ALICE3/Core/FastTracker.h"
#include "ALICE3/Core/DetLayer.h"
#endif

void testDetectorALICE3()
{

  o2::fastsim::FastTracker* alice3 = new o2::fastsim::FastTracker();
  alice3->SetIntegrationTime(100.e-6);          // 100 ns (as in LoI)
  alice3->SetMaxRadiusOfSlowDetectors(0.00001); // no slow detectors

  alice3->SetAvgRapidity(1.5);
  alice3->SetdNdEtaCent(2000);
  // alice3->SetdNdEtaCent(0);
  alice3->SetLhcUPCscale(1.);
  alice3->SetBField(1.); // Tesla
  // alice3->SetParticleMass(0.000511); // electron (default is pion)
  // alice3->SetParticleMass(1.32); // electron (default is pion)

  // Layer properties
  Double_t x0IT = 0.001;        // 0.1%
  Double_t x0OT = 0.010;        // 1%
  Double_t resRPhiIT = 0.00025; //  2.5 mum
  Double_t resZIT = 0.00025;    //  2.5 mum
  Double_t resRPhiOT = 0.0010;  // 10 mum
  Double_t resZOT = 0.0010;     // 10 mum
  Double_t eff = 0.98;
  //
  Double_t scaleR = 1.0;
  const float length = 250;
  //
  alice3->AddLayer("vertex", 0.0, length, 0, 0);        // dummy vertex for matrix calculation
  alice3->AddLayer("bpipe0", 0.48, length, 0.00042, 0); // 150 mum Be
  alice3->AddLayer("ddd0", 0.5, length, x0IT, resRPhiIT, resZIT, eff);
  alice3->AddLayer("ddd1", 1.2, length, x0IT, resRPhiIT, resZIT, eff);
  alice3->AddLayer("ddd2", 2.5, length, x0IT, resRPhiIT, resZIT, eff);
  alice3->AddLayer("bpipe1", 5.7, length, 0.0014, 0); // 500 mum Be
  alice3->AddLayer("ddd3", 7, length, x0OT, resRPhiOT, resZOT, eff);
  alice3->AddLayer("ddd4", 10., length, x0OT, resRPhiOT, resZOT, eff);
  alice3->AddLayer("ddd5", 13., length, x0OT, resRPhiOT, resZOT, eff);
  alice3->AddLayer("ddd6", 16., length, x0OT, resRPhiOT, resZOT, eff);
  alice3->AddLayer("iTOF", 20., length, 3. * x0OT, 0.03, 0.03, 1.0);
  alice3->AddLayer("ddd7", 30. * scaleR, length, x0OT, resRPhiOT, resZOT, eff);
  alice3->AddLayer("ddd8", 45. * scaleR, length, x0OT, resRPhiOT, resZOT, eff);
  alice3->AddLayer("ddd9", 60. * scaleR, length, x0OT, resRPhiOT, resZOT, eff);
  // alice3->AddLayer("ddd8",  50.*scaleR ,  length, x0OT, resRPhiOT, resZOT,eff);
  alice3->AddLayer("ddd10", 80. * scaleR, length, x0OT, resRPhiOT, resZOT, eff);
  // alice3->AddLayer("oTOF",  85.  , 0.03, 0.15, 0.15,eff); // 5 mm / sqrt(12)
  // alice3->AddLayer("rich",  110. , 0.03, 0.1, 0.1,eff);
  //
  //

  alice3->SetMinRadTrack(20.);
  // alice3->SetAtLeastHits(5);
  // alice3->SetAtLeastCorr(0);
  // alice3->SetAtLeastFake(4);
  //
  alice3->Print();
  // alice3->SolveViaBilloir(0);
  // alice3->MakeStandardPlots(0,2,2,kTRUE);

  // return;

  /*
  // eta scan
  alice3->SetAvgRapidity(1.);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,2,kTRUE);
  alice3->SetAvgRapidity(1.7);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,3,2,kTRUE);
  alice3->SetAvgRapidity(2.5);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,6,2,kTRUE);
  */

  /*
  // Variation of VD pixel resolution
  alice3->SetResolution("ddd0",0.0003,0.0003);
  alice3->SetResolution("ddd1",0.0003,0.0003);
  alice3->SetResolution("ddd2",0.0003,0.0003);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,1,2,kTRUE);
  alice3->SetResolution("ddd0",0.0004,0.0004);
  alice3->SetResolution("ddd1",0.0004,0.0004);
  alice3->SetResolution("ddd2",0.0004,0.0004);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,2,kTRUE);
  alice3->SetResolution("ddd0",0.0005,0.0005);
  alice3->SetResolution("ddd1",0.0005,0.0005);
  alice3->SetResolution("ddd2",0.0005,0.0005);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,3,2,kTRUE);
  */

  /*
  // Test VD setup closer to LoI description
  alice3->SetRadiationLength("bpipe0",0.00071); // 250 um Be
  alice3->SetRadiationLength("ddd0",0.00032); // 30 um Si
  alice3->SetRadiationLength("ddd1",0.00032); // 30 um Si
  alice3->SetRadiationLength("ddd2",0.00032); // 30 um Si
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,2,kTRUE);
  // Test thicker VD layers
  alice3->SetRadiationLength("ddd0",0.00053); // 50 um Si
  alice3->SetRadiationLength("ddd1",0.00053); // 50 um Si
  alice3->SetRadiationLength("ddd2",0.00053); // 50 um Si
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,4,kTRUE);
  // Test thicker VD layers
  alice3->SetRadiationLength("ddd0",0.00032+0.0008); // 30 um Si + 10 Cu + 40 Ka
  alice3->SetRadiationLength("ddd1",0.00032+0.0008); // 30 um Si + 10 Cu + 40 Ka
  alice3->SetRadiationLength("ddd2",0.00032+0.0008); // 30 um Si + 10 Cu + 40 Ka
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,5,kTRUE);
  // Test thicker VD layers
  alice3->SetRadiationLength("ddd0",0.00032+0.0016); // 30 um Si + 20 Cu + 60 Ka
  alice3->SetRadiationLength("ddd1",0.00032+0.0016); // 30 um Si + 20 Cu + 60 Ka
  alice3->SetRadiationLength("ddd2",0.00032+0.0016); // 30 um Si + 20 Cu + 60 Ka
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,6,kTRUE);
  */

  /*
  // Test VD setup with Felix suggestion
  alice3->SetRadiationLength("bpipe0",0.00071); // 250 um Be
  alice3->SetRadiationLength("ddd0",0.00076); // 30 um Si + 50 um glue + carbon foam 0.03%
  alice3->SetRadiationLength("ddd1",0.00096); // 30 um Si + 50 um glue + carbon foam 0.05%
  alice3->SetRadiationLength("ddd2",0.02167); // 30 um Si + 50 um glue + carbon foam 0.05% + cold plate (1.5 mm Al2O3 2%) + 0.07% Be case
  alice3->SetRadiationLength("bpipe1",0.0023); // 800 um Be
  //TCanvas *cL = new TCanvas("cL","Layout",0,0,800,800); alice3->PlotLayout(); return;
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(0,4,2,kTRUE);
  */
  /*
  // Test VD setup with Felix suggestion (try with cold plate made of Carbon X0=20 cm)
  alice3->SetRadiationLength("bpipe0",0.00071); // 250 um Be
  alice3->SetRadiationLength("ddd0",0.00076); // 30 um Si + 50 um glue + carbon foam 0.03%
  alice3->SetRadiationLength("ddd1",0.00096); // 30 um Si + 50 um glue + carbon foam 0.05%
  alice3->SetRadiationLength("ddd2",0.00917); // 30 um Si + 50 um glue + carbon foam 0.05% + cold plate (1.5 mm carbon 0.75%) + 0.07% Be case
  alice3->SetRadiationLength("bpipe1",0.0023); // 800 um Be
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,1,kTRUE);
  */

  // Test VD setup with Felix suggestion (try with cold plate as a separate layer - gives same result)
  alice3->SetRadiationLength("bpipe0", 0.0009);      // 150 um AlBe (X0=16.6cm)
  alice3->SetRadiationLength("ddd0", 0.00076);       // 30 um Si + 50 um glue + carbon foam 0.03%
  alice3->SetRadiationLength("ddd1", 0.00096);       // 30 um Si + 50 um glue + carbon foam 0.05%
  alice3->SetRadiationLength("ddd2", 0.00167);       // 30 um Si + 50 um glue + carbon foam 0.05% + 0.07% Be case
  alice3->AddLayer("coldplate", 2.6, 250, 0.02, 0.); // (1.5 mm Al2O3 2%)
  alice3->SetRadiationLength("bpipe1", 0.0023);      // 800 um Be
  // alice3->SolveViaBilloir(0);
  // alice3->MakeStandardPlots(0, 2, 1, kTRUE);

  // Add material of petals between L0 and L1
  // alice3->AddLayer("petalcase01", 0.9, 0.0036); // 300 um AlBe / sin(30) = 600 um (0.36% X0) at 0.9 cm
  // alice3->AddLayer("petalcase01", 0.9, 0.00512, resRPhiIT, resZIT,eff); // (300 um AlBe + 30 um Si + 50 um glue + carbon foam 0.03%) / sin(30) =   (0.18% + 0.076% X0)/ sin(30) = 0.512% at 0.9 cm
  // Add material of petals between L1 and L2
  // alice3->AddLayer("petalcase12", 1.9, 0.0050); // 600 um AlBe / sin(45) = 850 um (0.50% X0) at 1.9 cm
  // alice3->SolveViaBilloir(0);
  // alice3->MakeStandardPlots(1,4,2,kTRUE);

  // Increase pipe and L0 radius by 1 mm (new square shape with mean Rin of 6mm)
  alice3->SetRadius("bpipe0", 0.58);
  alice3->SetRadius("ddd0", 0.60);
  // alice3->SolveViaBilloir(0);
  // alice3->MakeStandardPlots(1, 2, 2, kTRUE);
  /*
  // Increase Si thickess to 66um (like ITS3)
  alice3->SetRadiationLength("ddd0",0.00115); // 66 um Si + 50 um glue + carbon foam 0.03%
  alice3->SetRadiationLength("ddd1",0.00135); // 66 um Si + 50 um glue + carbon foam 0.05%
  alice3->SetRadiationLength("ddd2",0.00206); // 66 um Si + 50 um glue + carbon foam 0.05% + 0.07% Be case
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,2,4,kTRUE);
  */

  // Variation of VD pixel resolution
  alice3->SetResolution("ddd0", 0.00025, 0.0004);
  alice3->SetResolution("ddd1", 0.00025, 0.0004);
  alice3->SetResolution("ddd2", 0.00025, 0.0004);
  // alice3->SolveViaBilloir(0);
  // alice3->MakeStandardPlots(1, 4, 2, kTRUE);
  alice3->SetResolution("ddd0", 0.00025, 0.0005);
  alice3->SetResolution("ddd1", 0.00025, 0.0005);
  alice3->SetResolution("ddd2", 0.00025, 0.0005);
  // alice3->SolveViaBilloir(0);
  // alice3->MakeStandardPlots(1, 3, 2, kTRUE);
  alice3->SetResolution("ddd0", 0.00025, 0.0010);
  alice3->SetResolution("ddd1", 0.00025, 0.0010);
  alice3->SetResolution("ddd2", 0.00025, 0.0010);
  // alice3->SolveViaBilloir(0);
  // alice3->MakeStandardPlots(1, 1, 2, kTRUE);

  /*
  // Test with layer 9 at 70 cm
  alice3->SetRadius("ddd9",70.);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,2,kTRUE);
  */
  /*
  alice3->RemoveLayer("ddd9");
  alice3->SetRadius("ddd8",50.);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,2,kTRUE);


  // Test with layer 10 dead
  alice3->KillLayer("ddd10");
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,1,2,kTRUE);

  alice3->SetRadius("ddd8",45.);
  alice3->AddLayer("ddd9",  70.*scaleR ,  x0OT, resRPhiOT, resZOT,eff);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,6,2,kTRUE);
  */

  /*
  // Test VD setup with Felix suggestion (with AlBe petal)
  alice3->SetRadiationLength("bpipe0",0.0009); // 150 um AlBe (X0=16.6 cm, x/X0=0.09%)
  alice3->SetRadiationLength("ddd0",0.00076); // 30 um Si + 50 um glue + carbon foam 0.03%
  alice3->SetRadiationLength("ddd1",0.00096); // 30 um Si + 50 um glue + carbon foam 0.05%
  alice3->SetRadiationLength("ddd2",0.02169); // 30 um Si + 50 um glue + carbon foam 0.05% + cold plate (1.5 mm Al2O3 2%) + 0.09% AlBe case
  alice3->SetRadiationLength("bpipe1",0.0023); // 800 um Be
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,1,2,kTRUE);
  */
  /*
  // Test VD setup with Felix suggestion (with cold plate 0.75 mm)
  alice3->SetRadiationLength("bpipe0",0.00071); // 250 um Be
  alice3->SetRadiationLength("ddd0",0.00076); // 30 um Si + 50 um glue + carbon foam 0.03%
  alice3->SetRadiationLength("ddd1",0.00096); // 30 um Si + 50 um glue + carbon foam 0.05%
  alice3->SetRadiationLength("ddd2",0.01167); // 30 um Si + 50 um glue + carbon foam 0.05% + cold plate (0.75 mm Al2O3 1%) + 0.07% Be case
  alice3->SetRadiationLength("bpipe1",0.0023); // 800 um Be
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,1,kTRUE);
  // Test VD setup with Felix suggestion (without cold plate)
  alice3->SetRadiationLength("bpipe0",0.00071); // 250 um Be
  alice3->SetRadiationLength("ddd0",0.00076); // 30 um Si + 50 um glue + carbon foam 0.03%
  alice3->SetRadiationLength("ddd1",0.00096); // 30 um Si + 50 um glue + carbon foam 0.05%
  alice3->SetRadiationLength("ddd2",0.00167); // 30 um Si + 50 um glue + carbon foam 0.05% + cold plate (1.5 mm Al2O3 2%) + 0.07% Be case
  alice3->SetRadiationLength("bpipe1",0.0023); // 800 um Be
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,3,2,kTRUE);
  */
  /*
  // Test VD setup with Felix suggestion (outer bpipe 500 um)
  alice3->SetRadiationLength("bpipe0",0.00071); // 250 um Be
  alice3->SetRadiationLength("ddd0",0.00076); // 30 um Si + 50 um glue + carbon foam 0.03%
  alice3->SetRadiationLength("ddd1",0.00096); // 30 um Si + 50 um glue + carbon foam 0.05%
  alice3->SetRadiationLength("ddd2",0.02167); // 30 um Si + 50 um glue + carbon foam 0.05% + cold plate (1.5 mm Al2O3 2%) + 0.07% Be case
  alice3->SetRadiationLength("bpipe1",0.0014); // 500 um Be
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,1,kTRUE);
  */

  /*
  // Variation of ML pixel resolution
  alice3->SetResolution("ddd3",0.0005,0.0005);
  alice3->SetResolution("ddd4",0.0005,0.0005);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,1,2,kTRUE);
  alice3->SetResolution("ddd3",0.0015,0.0015);
  alice3->SetResolution("ddd4",0.0015,0.0015);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,2,kTRUE);
  */

  /*
  // Variation of ML thickness
  alice3->SetRadiationLength("ddd6",0.006); // 0.6%
  alice3->SetRadiationLength("ddd5",0.006); // 0.6%
  alice3->SetRadiationLength("ddd4",0.006); // 0.6%
  alice3->SetRadiationLength("ddd3",0.006); // 0.6%
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,1,2,kTRUE);
  alice3->SetRadiationLength("ddd6",0.015); // 1.5%
  alice3->SetRadiationLength("ddd5",0.015); // 1.5%
  alice3->SetRadiationLength("ddd4",0.015); // 1.5%
  alice3->SetRadiationLength("ddd3",0.015); // 1.5%
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,7,2,kTRUE);
  alice3->SetRadiationLength("ddd6",0.0005); // 0.05%
  alice3->SetRadiationLength("ddd5",0.0005); // 0.05%
  alice3->SetRadiationLength("ddd4",0.0005); // 0.05%
  alice3->SetRadiationLength("ddd3",0.0005); // 0.05%
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,6,2,kTRUE);
  alice3->SetRadiationLength("ddd6",0.002); // 0.2%
  alice3->SetRadiationLength("ddd5",0.002); // 0.2%
  alice3->SetRadiationLength("ddd4",0.002); // 0.2%
  alice3->SetRadiationLength("ddd3",0.002); // 0.2%
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,3,kTRUE);
  alice3->SetRadiationLength("ddd6",0.001); // 0.1%
  alice3->SetRadiationLength("ddd5",0.001); // 0.1%
  alice3->SetRadiationLength("ddd4",0.001); // 0.1%
  alice3->SetRadiationLength("ddd3",0.001); // 0.1%
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,8,2,kTRUE);

  alice3->AddLayer("mlbus1",  11.,  0.007*0.42);  // 600 mum Al
  alice3->AddLayer("mlbus2",  18.,  0.007*0.25);  // 600 mum Al
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,8,4,kTRUE);
  */

  /*
  // Variation of 6th layer (20 cm) thickness
  alice3->SetRadiationLength("ddd6",0.02); // 2%
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,1,kTRUE);
  alice3->SetRadiationLength("ddd6",0.03); // 3%
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,6,1,kTRUE);
  */
  /*
  // New outer pipe radius of 5.5 cm
  alice3->SetRadius("bpipe1",  5.5);
  // Variation of innermost ML layer radius
  alice3->SetRadius("ddd3",  6.0);
  alice3->SetRadius("ddd4",  9.0);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,2,kTRUE);
  // change radii
  alice3->SetRadius("ddd3",  6.0);
  alice3->SetRadius("ddd4",  10.0);
  alice3->SetRadius("ddd5",  15.0);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,3,2,kTRUE);
  // remove one layer
  alice3->SetRadius("ddd5",  12.0);
  alice3->RemoveLayer("ddd4");
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,1,2,kTRUE);
  */

  /*
  // test photon converter outside iTOF
  alice3->AddLayer("converter",  22,  0.1);  // (10%)
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,1,2,kTRUE);
  alice3->SetRadiationLength("converter", 0.2);  // (20%)
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,1,3,kTRUE);
  alice3->SetRadiationLength("converter", 0.3);  // (30%)
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,1,4,kTRUE);
  */

  /*
  // Variation of OT pixel resolution
  alice3->SetResolution("ddd5",0.0012,0.0012);
  alice3->SetResolution("ddd6",0.0012,0.0012);
  alice3->SetResolution("ddd7",0.0012,0.0012);
  alice3->SetResolution("ddd8",0.0012,0.0012);
  alice3->SetResolution("ddd9",0.0012,0.0012);
  alice3->SetResolution("ddd10",0.0012,0.0012);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,1,2,kTRUE);
  alice3->SetResolution("ddd5",0.0012,0.0024);
  alice3->SetResolution("ddd6",0.0012,0.0024);
  alice3->SetResolution("ddd7",0.0012,0.0024);
  alice3->SetResolution("ddd8",0.0012,0.0024);
  alice3->SetResolution("ddd9",0.0012,0.0024);
  alice3->SetResolution("ddd10",0.0012,0.0024);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,3,2,kTRUE);
  alice3->SetResolution("ddd5",0.0012,0.0048);
  alice3->SetResolution("ddd6",0.0012,0.0048);
  alice3->SetResolution("ddd7",0.0012,0.0048);
  alice3->SetResolution("ddd8",0.0012,0.0048);
  alice3->SetResolution("ddd9",0.0012,0.0048);
  alice3->SetResolution("ddd10",0.0012,0.0048);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,2,kTRUE);
  alice3->SetResolution("ddd5",0.0024,0.0048);
  alice3->SetResolution("ddd6",0.0024,0.0048);
  alice3->SetResolution("ddd7",0.0024,0.0048);
  alice3->SetResolution("ddd8",0.0024,0.0048);
  alice3->SetResolution("ddd9",0.0024,0.0048);
  alice3->SetResolution("ddd10",0.0024,0.0048);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,6,2,kTRUE);
  */

  /*
  // Variation of 4 outer layer radii -15%
  scaleR=0.85;
  alice3->SetRadius("ddd7",30.*scaleR);
  alice3->SetRadius("ddd8",45.*scaleR);
  alice3->SetRadius("ddd9",60.*scaleR);
  alice3->SetRadius("ddd10",80.*scaleR);
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,2,kTRUE);
  */

  /*
  // test inclusion of bRICH in tracking
  alice3->AddLayer("OTcage",  82.,  0.01);  // 1% X0
  //alice3->AddLayer("oTOF",  85.  , 0.03, 0.15, 0.15,eff); // 5 mm / sqrt(12)
  //alice3->SolveViaBilloir(0);
  //alice3->MakeStandardPlots(1,6,2,kTRUE);
  alice3->AddLayer("richRad",  90.,  0.04);  // 4% X0 (Aerogel and cooling cage)
  alice3->AddLayer("rich",  110. , 0.03, 0.1, 0.1,eff);   // 1 mm res
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,4,2,kTRUE);
  alice3->SetResolution("rich",0.05, 0.05);   // 0.5 mm res
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,3,2,kTRUE);
  alice3->SetResolution("rich",0.01, 0.01);   // 0.1 mm res
  alice3->SolveViaBilloir(0);
  alice3->MakeStandardPlots(1,1,2,kTRUE);
  */

  return;
}

// void testDetectorCurr()
// {
//   o2::fastsim::FastTracker its("ALICE", "ITS");
//   its.MakeAliceCurrent(0, 0);
//   // its.SetAtLeastCorr(4);
//   // its.SetAtLeastFake(1);
//   its.Print();
//   its.SolveViaBilloir(0);

//   its.MakeStandardPlots(0, 2, 1, kTRUE);
//   //  return;
//   its.AddTPC(0.1, 0.1);
//   its.SolveViaBilloir(0);

//   its.MakeStandardPlots(1, 1, 1, kTRUE);
// }

// void particleDependendResolution()
// {
//   // particle dependency on resolution

//   // .L DetectorK.cxx++

//   o2::fastsim::FastTracker its("ALICE", "ITS");

//   its.MakeAliceCurrent();
//   its.Print();
//   its.PlotLayout();

//   its.SolveViaBilloir(0);

//   its.SetRadius("bpipe", 2.1);
//   its.AddLayer("spd0", 2.2, 0.001, 0.0012, 0.0012);

//   TCanvas* c1 = new TCanvas("c1", "c1");

//   c1->Divide(2, 1);
//   c1->cd(1);
//   gPad->SetGridx();
//   gPad->SetGridy();
//   gPad->SetLogx(); // gPad->SetLogy();
//   c1->cd(2);
//   gPad->SetGridx();
//   gPad->SetGridy();
//   gPad->SetLogx(); // gPad->SetLogy();

//   // compare to telescope equation ?
//   //  c1->cd(1); its.GetGraphPointingResolutionTeleEqu(0,1)->Draw("AC");
//   //  c1->cd(2); its.GetGraphPointingResolutionTeleEqu(1,1)->Draw("AC");

//   its.SetParticleMass(0.140); // pion
//   its.SolveViaBilloir(0, 0);
//   c1->cd(1);
//   its.GetGraphPointingResolution(0, 1)->Draw("AC");
//   c1->cd(2);
//   its.GetGraphPointingResolution(1, 1)->Draw("AC");

//   its.SetParticleMass(0.498); // kaon
//   its.SolveViaBilloir(0, 0);
//   c1->cd(1);
//   its.GetGraphPointingResolution(0, 2)->Draw("C");
//   c1->cd(2);
//   its.GetGraphPointingResolution(1, 2)->Draw("C");

//   its.SetParticleMass(0.00051); // electron
//   its.SolveViaBilloir(0, 0);
//   c1->cd(1);
//   its.GetGraphPointingResolution(0, 3)->Draw("C");
//   c1->cd(2);
//   its.GetGraphPointingResolution(1, 3)->Draw("C");

//   its.SetParticleMass(0.938); // proton
//   its.SolveViaBilloir(0, 0);
//   c1->cd(1);
//   its.GetGraphPointingResolution(0, 4)->Draw("C");
//   c1->cd(2);
//   its.GetGraphPointingResolution(1, 4)->Draw("C");
// }
