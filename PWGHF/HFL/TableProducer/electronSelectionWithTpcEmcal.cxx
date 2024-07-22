#include "Pythia8/HeavyIons.h"
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/ColourReconnectionHooks.h"
#include "TApplication.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include <THnSparse.h>
#include <iostream>
#include <new>

using namespace std;
Double_t pi = M_PI;

int main(int argc, char *argv[]) {

  Pythia8::Pythia pythia;
  int multBins = 500;
  int fNDelPhiBins = 32;

  int seed = 0;
  int nevents = 1000;
  for (int ig = 0; ig < argc; ++ig) {
    TString str = argv[ig];
    Printf("str %d %s", ig, argv[ig]);
    if (str.Contains("--events")) {
      sscanf(argv[ig + 1], "%d", &nevents);
    }
    if (str.Contains("--seed")) {
      sscanf(argv[ig + 1], "%d", &seed);
    }
  }

  TFile *outFile =
      new TFile(Form("ex_PbPb_2760_Ev%i_Sd%i.root", nevents, seed), "RECREATE");

  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed %d", seed));
  pythia.readString("Beams:idA = 1000822080");
  pythia.readString("Beams:idB = 1000822080");
  pythia.readString("Beams:eCM = 2760.");
  pythia.readString("Beams:frameType = 1"); // the beams are colliding in their

  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10");

  pythia.readString("PartonVertex:setVertex = on");

  pythia.readString("SoftQCD:all = on");

  pythia.readString("ColourReconnection:reconnect = on");

  pythia.readString("ColourReconnection:mode = 0");

  pythia.readString("113:addChannel = 1 0.9988447 2 211 -211");
  pythia.readString("113:addChannel = 1 0.0000180 0 211 -211 211 -211 ");
  pythia.readString("113:addchannel =1 0.0001009 1 211 -211 111");
  pythia.readString(
      "223:addChannel = 1 0.0154283 2 211 -211"); // omega add channel
  pythia.readString(
      "223:addChannel = 1 0.8994773 1 211 -211 111"); // omegA channel add
  pythia.readString(
      "221: addchannel = 1  0.2274105 211 -211 111"); // eta  channel add
  pythia.readString(
      "221:addChannel =  1  0.0460021 211 -211 22"); //  eta channel add
  pythia.readString(
      "221:addchannel: = 1 0.0002680 12 211 -211 11 -11"); // eta channel add

  // for Heavy ion collision
  pythia.readString(
      "HeavyIon:SigFitDefPar =  12.32,1.63,0.19,0.0,0.0,0.0,0.0,0.0");
  pythia.readString("HeavyIon:SigFitNGen = 20");
  pythia.readString("HeavyIon:SigFitErr = "
                    "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
  pythia.init();

  Int_t nbins[5] = {32., 30., 40., 40., 150.};
  Double_t xmin[5] = {-TMath::Pi(), -1.6, 0., 0., -0.5};
  Double_t xmax[5] = {TMath::Pi(), 1.6, 20., 20., 149.5};

  THnSparseD *sparseULpipi =
      new THnSparseD("sparseULpipi", "sparse1", 5, nbins, xmin, xmax);
  THnSparseD *sparseLpipiplus = new THnSparseD(
      "sparseLpipiplus ", "sparseLpipiplus", 5, nbins, xmin, xmax);
  THnSparseD *sparseULpipioppo = new THnSparseD(
      "sparseULpipioppo", "sparseULpipioppo", 5, nbins, xmin, xmax);
  THnSparseD *sparseLpipineg =
      new THnSparseD("sparseLpipineg", "sparseLpipineg", 5, nbins, xmin, xmax);
  TH2D *pipluspinegcorr =
      new TH2D("pipluspinegcorr", "pipluspinegcorr", 32, -0.5 * TMath::Pi(),
               1.5 * TMath::Pi(), 30, -1.6, 1.6);

  THnSparseD *sparseULkk =
      new THnSparseD("sparseULkk", "sparse1", 5, nbins, xmin, xmax);
  THnSparseD *sparseLkkplus =
      new THnSparseD("sparseLkkplus ", "sparseLkkplus", 5, nbins, xmin, xmax);
  THnSparseD *sparseULkkoppo =
      new THnSparseD("sparseULkkoppo", "sparseULkkoppo", 5, nbins, xmin, xmax);
  THnSparseD *sparseLkkneg =
      new THnSparseD("sparseLkkneg", "sparseLkkneg", 5, nbins, xmin, xmax);

  THnSparseD *sparseULpp =
      new THnSparseD("sparseULpp", "sparse1", 5, nbins, xmin, xmax);
  THnSparseD *sparseLppplus =
      new THnSparseD("sparseLppplus ", "sparseLppplus", 5, nbins, xmin, xmax);
  THnSparseD *sparseULppoppo =
      new THnSparseD("sparseULppoppo", "sparseULppoppo", 5, nbins, xmin, xmax);
  THnSparseD *sparseLppneg =
      new THnSparseD("sparseLppneg", "sparseLppneg", 5, nbins, xmin, xmax);

  double *asparseULpipi = new double[5];
  double *asparseLpipiplus = new double[5];
  double *asparseULpipioppo = new double[5];
  double *asparseLpipineg = new double[5];
  double *asparseULkk = new double[5];
  double *asparseLkkplus = new double[5];
  double *asparseULkkoppo = new double[5];
  double *asparseLkkneg = new double[5];
  double *asparseULpp = new double[5];
  double *asparseLppplus = new double[5];
  double *asparseULppoppo = new double[5];
  double *asparseLppneg = new double[5];
  TH1D *multiplicity =
      new TH1D("multiplicity", "multiplicity", 150, -0.5, 149.5);
  TH1D *multiplicityV0 =
      new TH1D("multiplicityV0", "multiplicityV0", 150, -0.5, 149.5);
  TH1D *multiplicitytrk =
      new TH1D("multiplicitytrk", "multiplicitytrk", 150, -0.5, 149.5);

  TH1D *piptplustri = new TH1D("piptplustri", "piptplustri", 100, 0., 20.0);
  TH1D *piptnegtri = new TH1D("piptnegtri", "piptnegtri", 100, 0., 20.0);
  TH1D *kptplustri = new TH1D("kiptplustri", "kptplustri", 100, 0., 20.0);
  TH1D *kptnegtri = new TH1D("kptnegtri", "kptnegtri", 100, 0., 20.0);
  TH1D *proptplustri = new TH1D("proptplustri", "proptplustri", 100, 0., 20.0);
  TH1D *proptnegtri = new TH1D("proptnegtri", "proptnegtri", 100, 0., 20.0);

  for (int iEvent = 0; iEvent < nevents; ++iEvent) {

    int nCharged = 0;
    int nChargedV0 = 0;
    int nChargedntrk = 0;
    if (!pythia.next())
      continue;
    int enteries = pythia.event.size();

    Printf("Event -- %d", iEvent);
    for (int k = 0; k < pythia.event.size(); ++k) {

      if (pythia.event[k].isFinal() && pythia.event[k].isCharged() &&
          (pythia.event[k].tau() > 10 || (!pythia.event[k].canDecay()))) {
        if (TMath::Abs(pythia.event[k].eta()) < 0.8)
          ++nCharged;
        if (pythia.event[k].status() / 10 == 8) {

          if (-3.7 < pythia.event[k].eta() && pythia.event[k].eta() < -1.7)
            ++nChargedV0;
          if (pythia.event[k].eta() > 2.8 && pythia.event[k].eta() < 5.1)
            ++nChargedV0;
          if (TMath::Abs(pythia.event[k].eta()) < 0.8)
            ++nChargedntrk;
          else
            continue;
        }
      }
    }
    if (nCharged != 0 || nChargedV0 != 0 || nChargedntrk != 0) {
      multiplicity->Fill(nCharged);
      multiplicityV0->Fill(nChargedV0);
      multiplicitytrk->Fill(nChargedntrk);
    }
    int inel = 0;
    for (int j = 0; j < enteries; j++) {
      if (pythia.event[j].isCharged() && pythia.event[j].isFinal()) {
        if (pythia.event[j].pT() > 0 && pythia.event[j].pT() < 40) {
          if (abs(pythia.event[j].y()) < 0.8) {
            inel = 1;
          }
        }
      }
    }

    Double_t Pttrig = 0;
    Double_t Phitrig = 0;
    Double_t ytrig = 0;
    for (int i = 0; i < enteries; i++) {

      Pythia8::Particle &particle = pythia.event[i];

      if (inel == 1 && pythia.event[i].isCharged() &&
          pythia.event[i].isFinal()) {

        int idAbs = pythia.event[i].id();
        Pttrig = pythia.event[i].pT();
        Phitrig = pythia.event[i].phi();
        ytrig = pythia.event[i].y();
        // for pion like and unlike sign
        /////////////////// ____________correlation for charge anticharged
        /// particle
        ///__________________//////////////////////////

        if (abs(pythia.event[i].y()) <= 0.8) {
          // pion plus correlate with plus and negetive
          if (pythia.event[i].id() == 211) {

            if (pythia.event[i].pT() < 0.2 || pythia.event[i].pT() > 2.0)
              continue;

            piptplustri->Fill(Pttrig);
            Double_t diffy = -999;
            Double_t diffphi = -999;
            Double_t Ptassoc = 0;
            for (int k = 0; k < enteries; k++) {

              if (pythia.event[k].index() == pythia.event[i].index())
                continue;

              if (inel == 1 && pythia.event[k].isCharged() &&
                  pythia.event[k].isFinal()) {
                int idAbsanti = pythia.event[k].id();
                if (abs(pythia.event[k].y()) > 0.8)
                  continue;
                if (pythia.event[k].id() ==
                    -211) /// plus charge  negetive charge correlation
                {
                  if (pythia.event[k].pT() < 0.2 || pythia.event[k].pT() > 2.0)
                    continue;

                  if (pythia.event[i].pT() > pythia.event[k].pT()) {
                    diffy = pythia.event[i].y() - pythia.event[k].y();
                    diffphi = pythia.event[i].phi() - pythia.event[k].phi();
                    Ptassoc = pythia.event[k].pT();
                    if (diffphi < -TMath::Pi()) {
                      diffphi = diffphi + 2 * pi;
                    }
                    if (diffphi > TMath::Pi()) {
                      diffphi = diffphi - 2 * pi;
                    }

                    asparseULpipi[0] = diffphi;
                    asparseULpipi[1] = diffy;
                    asparseULpipi[2] = Pttrig;
                    asparseULpipi[3] = Ptassoc;
                    asparseULpipi[4] = nCharged;

                    sparseULpipi->Fill(asparseULpipi);

                    pipluspinegcorr->Fill(diffphi, diffy);
                  }
                }
                if (pythia.event[k].id() == 211) /// plus plus correlation
                {
                  if (pythia.event[k].pT() < 0.2 || pythia.event[k].pT() > 2.0)
                    continue;

                  if (pythia.event[i].pT() > pythia.event[k].pT()) {
                    diffy = pythia.event[i].y() - pythia.event[k].y();
                    diffphi = pythia.event[i].phi() - pythia.event[k].phi();

                    if (diffphi < -TMath::Pi()) {
                      diffphi = diffphi + 2 * pi;
                    }
                    if (diffphi > TMath::Pi()) {
                      diffphi = diffphi - 2 * pi;
                    }

                    asparseLpipiplus[0] = diffphi;
                    asparseLpipiplus[1] = diffy;
                    asparseLpipiplus[2] = Pttrig;
                    asparseLpipiplus[3] = pythia.event[k].pT();
                    asparseLpipiplus[4] = nCharged;

                    sparseLpipiplus->Fill(asparseLpipiplus);
                  }
                }
              }
            }
          }

          if (pythia.event[i].id() ==
              -211) { // charge negetive correlate with charge positive

            if (pythia.event[i].pT() < 0.2 || pythia.event[i].pT() > 2.0)
              continue;
            piptnegtri->Fill(Pttrig);
            Double_t diffy = -999;
            Double_t diffphi = -999;
            Double_t Ptassoc = 0;
            for (int k = 0; k < enteries; k++) {

              if (pythia.event[k].index() == pythia.event[i].index())
                continue;

              if (inel == 1 && pythia.event[k].isCharged() &&
                  pythia.event[k].isFinal()) {
                int idAbsanti = pythia.event[k].id();
                if (abs(pythia.event[k].y()) > 0.8)
                  continue;
                if (pythia.event[k].id() == -211) /// neeg neg correlation
                {

                  if (pythia.event[k].pT() < 0.2 || pythia.event[k].pT() > 2.0)
                    continue;

                  if (pythia.event[i].pT() > pythia.event[k].pT()) {
                    diffy = pythia.event[i].y() - pythia.event[k].y();
                    diffphi = pythia.event[i].phi() - pythia.event[k].phi();
                    Ptassoc = pythia.event[k].pT();
                    if (diffphi < -TMath::Pi()) {
                      diffphi = diffphi + 2 * pi;
                    }
                    if (diffphi > TMath::Pi()) {
                      diffphi = diffphi - 2 * pi;
                    }

                    asparseLpipineg[0] = diffphi;
                    asparseLpipineg[1] = diffy;
                    asparseLpipineg[2] = Pttrig;
                    asparseLpipineg[3] = Ptassoc;
                    asparseLpipineg[4] = nCharged;

                    sparseLpipineg->Fill(asparseLpipineg);
                  }
                }

                if (pythia.event[k].id() == 211) /// neg plus correlation
                {
                  if (pythia.event[k].pT() < 0.2 || pythia.event[k].pT() > 2.0)
                    continue;

                  if (pythia.event[i].pT() > pythia.event[k].pT()) {
                    diffy = pythia.event[i].y() - pythia.event[k].y();
                    diffphi = pythia.event[i].phi() - pythia.event[k].phi();
                    Ptassoc = pythia.event[k].pT();
                    if (diffphi < -TMath::Pi()) {
                      diffphi = diffphi + 2 * pi;
                    }
                    if (diffphi > TMath::Pi()) {
                      diffphi = diffphi - 2 * pi;
                    }

                    asparseULpipioppo[0] = diffphi;
                    asparseULpipioppo[1] = diffy;
                    asparseULpipioppo[2] = Pttrig;
                    asparseULpipioppo[3] = Ptassoc;
                    asparseULpipioppo[4] = nCharged;
                    sparseULpipioppo->Fill(asparseULpipioppo);
                  }
                }
              }
            }
          }
        }

        /// correlation for kaon kaon

        if (abs(pythia.event[i].y()) <= 0.7) {
          // kaon plus correlate with plus and negetive
          if (pythia.event[i].id() == 321) {

            if (pythia.event[i].pT() < 0.2 || pythia.event[i].pT() > 2.0)
              continue;

            kptplustri->Fill(Pttrig);
            Double_t diffy = -999;
            Double_t diffphi = -999;
            Double_t Ptassoc = 0;
            for (int k = 0; k < enteries; k++) {

              if (pythia.event[k].index() == pythia.event[i].index())
                continue;

              if (inel == 1 && pythia.event[k].isCharged() &&
                  pythia.event[k].isFinal()) {
                int idAbsanti = pythia.event[k].id();
                if (abs(pythia.event[k].y()) > 0.7)
                  continue;
                if (pythia.event[k].id() ==
                    -321) /// Kaon plus charge  negetive charge correlation
                {
                  if (pythia.event[k].pT() < 0.2 || pythia.event[k].pT() > 2.0)
                    continue;

                  if (pythia.event[i].pT() > pythia.event[k].pT()) {
                    diffy = pythia.event[i].y() - pythia.event[k].y();
                    diffphi = pythia.event[i].phi() - pythia.event[k].phi();
                    Ptassoc = pythia.event[k].pT();
                    if (diffphi < -TMath::Pi()) {
                      diffphi = diffphi + 2 * pi;
                    }
                    if (diffphi > TMath::Pi()) {
                      diffphi = diffphi - 2 * pi;
                    }

                    asparseULkk[0] = diffphi;
                    asparseULkk[1] = diffy;
                    asparseULkk[2] = Pttrig;
                    asparseULkk[3] = Ptassoc;
                    asparseULkk[4] = nCharged;

                    sparseULkk->Fill(asparseULkk);
                  }
                }
                if (pythia.event[k].id() ==
                    321) /// kaon plus kaon plus correlation
                {
                  if (pythia.event[k].pT() < 0.2 || pythia.event[k].pT() > 2.0)
                    continue;

                  if (pythia.event[i].pT() > pythia.event[k].pT()) {
                    diffy = pythia.event[i].y() - pythia.event[k].y();
                    diffphi = pythia.event[i].phi() - pythia.event[k].phi();

                    if (diffphi < -TMath::Pi()) {
                      diffphi = diffphi + 2 * pi;
                    }
                    if (diffphi > TMath::Pi()) {
                      diffphi = diffphi - 2 * pi;
                    }

                    asparseLkkplus[0] = diffphi;
                    asparseLkkplus[1] = diffy;
                    asparseLkkplus[2] = Pttrig;
                    asparseLkkplus[3] = pythia.event[k].pT();
                    asparseLkkplus[4] = nCharged;

                    sparseLkkplus->Fill(asparseLkkplus);
                  }
                }
              }
            }
          }

          if (pythia.event[i].id() == -321) { // kaon charge negetive correlate
                                              // with kaon charge positive

            if (pythia.event[i].pT() < 0.2 || pythia.event[i].pT() > 2.0)
              continue;
            kptnegtri->Fill(Pttrig);
            Double_t diffy = -999;
            Double_t diffphi = -999;
            Double_t Ptassoc = 0;
            for (int k = 0; k < enteries; k++) {

              if (pythia.event[k].index() == pythia.event[i].index())
                continue;

              if (inel == 1 && pythia.event[k].isCharged() &&
                  pythia.event[k].isFinal()) {
                int idAbsanti = pythia.event[k].id();
                if (abs(pythia.event[k].y()) > 0.7)
                  continue;
                if (pythia.event[k].id() ==
                    -321) /// kaon neg  kaon neg correlation
                {

                  if (pythia.event[k].pT() < 0.2 || pythia.event[k].pT() > 2.0)
                    continue;

                  if (pythia.event[i].pT() > pythia.event[k].pT()) {
                    diffy = pythia.event[i].y() - pythia.event[k].y();
                    diffphi = pythia.event[i].phi() - pythia.event[k].phi();
                    Ptassoc = pythia.event[k].pT();
                    if (diffphi < -TMath::Pi()) {
                      diffphi = diffphi + 2 * pi;
                    }
                    if (diffphi > TMath::Pi()) {
                      diffphi = diffphi - 2 * pi;
                    }

                    asparseLkkneg[0] = diffphi;
                    asparseLkkneg[1] = diffy;
                    asparseLkkneg[2] = Pttrig;
                    asparseLkkneg[3] = Ptassoc;
                    asparseLkkneg[4] = nCharged;

                    sparseLkkneg->Fill(asparseLkkneg);
                  }
                }

                if (pythia.event[k].id() == 321) /// neg plus correlation
                {
                  if (pythia.event[k].pT() < 0.2 || pythia.event[k].pT() > 2.0)
                    continue;

                  if (pythia.event[i].pT() > pythia.event[k].pT()) {
                    diffy = pythia.event[i].y() - pythia.event[k].y();
                    diffphi = pythia.event[i].phi() - pythia.event[k].phi();
                    Ptassoc = pythia.event[k].pT();
                    if (diffphi < -TMath::Pi()) {
                      diffphi = diffphi + 2 * pi;
                    }
                    if (diffphi > TMath::Pi()) {
                      diffphi = diffphi - 2 * pi;
                    }

                    asparseULkkoppo[0] = diffphi;
                    asparseULkkoppo[1] = diffy;
                    asparseULkkoppo[2] = Pttrig;
                    asparseULkkoppo[3] = Ptassoc;
                    asparseULkkoppo[4] = nCharged;
                    sparseULkkoppo->Fill(asparseULkkoppo);
                  }
                }
              }
            }
          }
        }

        if (abs(pythia.event[i].y()) <= 0.6) {
          // proton plus correlate with plus and negetive
          if (pythia.event[i].id() == 2212) {

            if (pythia.event[i].pT() < 0.5 || pythia.event[i].pT() > 2.5)
              continue;

            proptplustri->Fill(Pttrig);
            Double_t diffy = -999;
            Double_t diffphi = -999;
            Double_t Ptassoc = 0;
            for (int k = 0; k < enteries; k++) {

              if (pythia.event[k].index() == pythia.event[i].index())
                continue;

              if (inel == 1 && pythia.event[k].isCharged() &&
                  pythia.event[k].isFinal()) {
                int idAbsanti = pythia.event[k].id();
                if (abs(pythia.event[k].y()) > 0.6)
                  continue;
                if (pythia.event[k].id() ==
                    -2212) /// proton plus charge  negetive charge correlation
                {
                  if (pythia.event[k].pT() < 0.5 || pythia.event[k].pT() > 2.5)
                    continue;

                  if (pythia.event[i].pT() > pythia.event[k].pT()) {
                    diffy = pythia.event[i].y() - pythia.event[k].y();
                    diffphi = pythia.event[i].phi() - pythia.event[k].phi();
                    Ptassoc = pythia.event[k].pT();
                    if (diffphi < -TMath::Pi()) {
                      diffphi = diffphi + 2 * pi;
                    }
                    if (diffphi > TMath::Pi()) {
                      diffphi = diffphi - 2 * pi;
                    }

                    asparseULpp[0] = diffphi;
                    asparseULpp[1] = diffy;
                    asparseULpp[2] = Pttrig;
                    asparseULpp[3] = Ptassoc;
                    asparseULpp[4] = nCharged;

                    sparseULpp->Fill(asparseULpp);

                    pipluspinegcorr->Fill(diffphi, diffy);
                  }
                }
                if (pythia.event[k].id() == 2212) /// plus plus correlation
                {
                  if (pythia.event[k].pT() < 0.5 || pythia.event[k].pT() > 2.5)
                    continue;

                  if (pythia.event[i].pT() > pythia.event[k].pT()) {
                    diffy = pythia.event[i].y() - pythia.event[k].y();
                    diffphi = pythia.event[i].phi() - pythia.event[k].phi();

                    if (diffphi < -TMath::Pi()) {
                      diffphi = diffphi + 2 * pi;
                    }
                    if (diffphi > TMath::Pi()) {
                      diffphi = diffphi - 2 * pi;
                    }

                    asparseLppplus[0] = diffphi;
                    asparseLppplus[1] = diffy;
                    asparseLppplus[2] = Pttrig;
                    asparseLppplus[3] = pythia.event[k].pT();
                    asparseLppplus[4] = nCharged;

                    sparseLppplus->Fill(asparseLppplus);
                  }
                }
              }
            }
          }

          if (pythia.event[i].id() ==
              -2212) { // proton charge negetive correlate with charge positive

            if (pythia.event[i].pT() < 0.5 || pythia.event[i].pT() > 2.5)
              continue;
            proptnegtri->Fill(Pttrig);
            Double_t diffy = -999;
            Double_t diffphi = -999;
            Double_t Ptassoc = 0;
            for (int k = 0; k < enteries; k++) {

              if (pythia.event[k].index() == pythia.event[i].index())
                continue;

              if (inel == 1 && pythia.event[k].isCharged() &&
                  pythia.event[k].isFinal()) {
                int idAbsanti = pythia.event[k].id();
                if (abs(pythia.event[k].y()) > 0.6)
                  continue;
                if (pythia.event[k].id() ==
                    -2212) /// proton neg  proton neg correlation
                {

                  if (pythia.event[k].pT() < 0.5 || pythia.event[k].pT() > 2.5)
                    continue;

                  if (pythia.event[i].pT() > pythia.event[k].pT()) {
                    diffy = pythia.event[i].y() - pythia.event[k].y();
                    diffphi = pythia.event[i].phi() - pythia.event[k].phi();
                    Ptassoc = pythia.event[k].pT();
                    if (diffphi < -TMath::Pi()) {
                      diffphi = diffphi + 2 * pi;
                    }
                    if (diffphi > TMath::Pi()) {
                      diffphi = diffphi - 2 * pi;
                    }

                    asparseLppneg[0] = diffphi;
                    asparseLppneg[1] = diffy;
                    asparseLppneg[2] = Pttrig;
                    asparseLppneg[3] = Ptassoc;
                    asparseLppneg[4] = nCharged;

                    sparseLppneg->Fill(asparseLppneg);
                  }
                }

                if (pythia.event[k].id() == 2212) /// neg plus correlation
                {
                  if (pythia.event[k].pT() < 0.5 || pythia.event[k].pT() > 2.5)
                    continue;

                  if (pythia.event[i].pT() > pythia.event[k].pT()) {
                    diffy = pythia.event[i].y() - pythia.event[k].y();
                    diffphi = pythia.event[i].phi() - pythia.event[k].phi();
                    Ptassoc = pythia.event[k].pT();
                    if (diffphi < -TMath::Pi()) {
                      diffphi = diffphi + 2 * pi;
                    }
                    if (diffphi > TMath::Pi()) {
                      diffphi = diffphi - 2 * pi;
                    }

                    asparseULppoppo[0] = diffphi;
                    asparseULppoppo[1] = diffy;
                    asparseULppoppo[2] = Pttrig;
                    asparseULppoppo[3] = Ptassoc;
                    asparseULppoppo[4] = nCharged;
                    sparseULppoppo->Fill(asparseULppoppo);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  pythia.stat();

  sparseULpipi->Write();
  sparseLpipineg->Write();
  sparseLpipiplus->Write();
  sparseULpipioppo->Write();
  sparseULkk->Write();
  sparseLkkneg->Write();
  sparseLkkplus->Write();
  sparseULkkoppo->Write();
  sparseULpp->Write();
  sparseLppneg->Write();
  sparseLppplus->Write();
  sparseULppoppo->Write();

  outFile->Write();
  delete outFile;

  return 0;
}
