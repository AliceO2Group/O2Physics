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
/// @file DelphesO2LutWriter.cxx
/// @brief Porting to O2Physics of DelphesO2 code.
///        Minimal changes have been made to the original code for adaptation purposes, formatting and commented parts have been considered.
///        Relevant sources:
///                 DelphesO2/src/lutWrite.cc https://github.com/AliceO2Group/DelphesO2/blob/master/src/lutWrite.cc
/// @author: Roberto Preghenella
/// @email: preghenella@bo.infn.it
///

#include <cstdio>

#include "ALICE3/Core/DelphesO2TrackSmearer.h"
#include "ALICE3/Core/DelphesO2LutWriter.h"
#include <iostream>
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TAxis.h"
#include "TMatrixDSymEigen.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "ALICE3/Core/FastTracker.h"
#include "ALICE3/Core/TrackUtilities.h"

// #define USE_FWD_PARAM
#ifdef USE_FWD_PARAM
#include "fwdRes.C"
#endif

namespace o2::fastsim
{

void DelphesO2LutWriter::printLutWriterConfiguration()
{
  std::cout << " --- Printing configuration of LUT writer --- " << std::endl;
  std::cout << "    -> etaMaxBarrel  = " << etaMaxBarrel << std::endl;
  std::cout << "    -> usePara       = " << usePara << std::endl;
  std::cout << "    -> useDipole     = " << useDipole << std::endl;
  std::cout << "    -> useFlatDipole = " << useFlatDipole << std::endl;
}

bool DelphesO2LutWriter::fatSolve(lutEntry_t& lutEntry,
                                  float pt,
                                  float eta,
                                  const float mass,
                                  int itof,
                                  int otof,
                                  int q,
                                  const float nch)
{
  lutEntry.valid = false;

  static TLorentzVector tlv;
  tlv.SetPtEtaPhiM(pt, eta, 0., mass);
  o2::track::TrackParCov trkIn;
  o2::upgrade::convertTLorentzVectorToO2Track(q, tlv, {0., 0., 0.}, trkIn);
  o2::track::TrackParCov trkOut;
  const int status = fat.FastTrack(trkIn, trkOut, nch);
  if (status <= 0) {
    Printf(" --- fatSolve: FastTrack failed --- \n");
    tlv.Print();
    return false;
  }
  lutEntry.valid = true;
  lutEntry.itof = fat.GetGoodHitProb(itof);
  lutEntry.otof = fat.GetGoodHitProb(otof);
  for (int i = 0; i < 15; ++i)
    lutEntry.covm[i] = trkOut.getCov()[i];

  // define the efficiency
  auto totfake = 0.;
  lutEntry.eff = 1.;
  for (int i = 1; i < 20; ++i) {
    auto igoodhit = fat.GetGoodHitProb(i);
    if (igoodhit <= 0. || i == itof || i == otof)
      continue;
    lutEntry.eff *= igoodhit;
    auto pairfake = 0.;
    for (int j = i + 1; j < 20; ++j) {
      auto jgoodhit = fat.GetGoodHitProb(j);
      if (jgoodhit <= 0. || j == itof || j == otof)
        continue;
      pairfake = (1. - igoodhit) * (1. - jgoodhit);
      break;
    }
    totfake += pairfake;
  }
  lutEntry.eff2 = (1. - totfake);

  return true;
}

#ifdef USE_FWD_PARAM
bool DelphesO2LutWriter::fwdSolve(float* covm, float pt, float eta, float mass)
{
  if (fwdRes(covm, pt, eta, mass) < 0)
    return false;
  return true;
}
#else
bool DelphesO2LutWriter::fwdSolve(float*, float, float, float)
{
  return false;
}
#endif

bool DelphesO2LutWriter::fwdPara(lutEntry_t& lutEntry, float pt, float eta, float mass, float Bfield)
{
  lutEntry.valid = false;

  // parametrised forward response; interpolates between FAT at eta = 1.75 and a fixed parametrisation at eta = 4; only diagonal elements
  if (std::fabs(eta) < etaMaxBarrel || std::fabs(eta) > 4)
    return false;

  if (!fatSolve(lutEntry, pt, etaMaxBarrel, mass))
    return false;
  float covmbarrel[15] = {0};
  for (int i = 0; i < 15; ++i) {
    covmbarrel[i] = lutEntry.covm[i];
  }

  // parametrisation at eta = 4
  const double beta = 1. / std::sqrt(1 + mass * mass / pt / pt / std::cosh(eta) / std::cosh(eta));
  const float dca_pos = 2.5e-4 / std::sqrt(3); // 2.5 micron/sqrt(3)
  const float r0 = 0.5;                        // layer 0 radius [cm]
  const float r1 = 1.3;
  const float r2 = 2.5;
  const float x0layer = 0.001; // material budget (rad length) per layer
  const double sigma_alpha = 0.0136 / beta / pt * std::sqrt(x0layer * std::cosh(eta)) * (1 + 0.038 * std::log(x0layer * std::cosh(eta)));
  const double dcaxy_ms = sigma_alpha * r0 * std::sqrt(1 + r1 * r1 / (r2 - r0) / (r2 - r0));
  const double dcaxy2 = dca_pos * dca_pos + dcaxy_ms * dcaxy_ms;

  const double dcaz_ms = sigma_alpha * r0 * std::cosh(eta);
  const double dcaz2 = dca_pos * dca_pos + dcaz_ms * dcaz_ms;

  const float Leta = 2.8 / sinh(eta) - 0.01 * r0; // m
  const double relmomres_pos = 10e-6 * pt / 0.3 / Bfield / Leta / Leta * std::sqrt(720. / 15.);

  const float relmomres_barrel = std::sqrt(covmbarrel[14]) * pt;
  const float Router = 1; // m
  const float relmomres_pos_barrel = 10e-6 * pt / 0.3 / Bfield / Router / Router / std::sqrt(720. / 15.);
  const float relmomres_MS_barrel = std::sqrt(relmomres_barrel * relmomres_barrel - relmomres_pos_barrel * relmomres_pos_barrel);

  // interpolate MS contrib (rel resolution 0.4 at eta = 4)
  const float relmomres_MS_eta4 = 0.4 / beta * 0.5 / Bfield;
  const float relmomres_MS = relmomres_MS_eta4 * pow(relmomres_MS_eta4 / relmomres_MS_barrel, (std::fabs(eta) - 4.) / (4. - etaMaxBarrel));
  const float momres_tot = pt * std::sqrt(relmomres_pos * relmomres_pos + relmomres_MS * relmomres_MS); // total absolute mom reso

  // Fill cov matrix diag
  for (int i = 0; i < 15; ++i)
    lutEntry.covm[i] = 0;

  lutEntry.covm[0] = covmbarrel[0];
  if (dcaxy2 > lutEntry.covm[0])
    lutEntry.covm[0] = dcaxy2;
  lutEntry.covm[2] = covmbarrel[2];
  if (dcaz2 > lutEntry.covm[2])
    lutEntry.covm[2] = dcaz2;
  lutEntry.covm[5] = covmbarrel[5];                                // sigma^2 sin(phi)
  lutEntry.covm[9] = covmbarrel[9];                                // sigma^2 tanl
  lutEntry.covm[14] = momres_tot * momres_tot / pt / pt / pt / pt; // sigma^2 1/pt
  // Check that all numbers are numbers
  for (int i = 0; i < 15; ++i) {
    if (std::isnan(lutEntry.covm[i])) {
      Printf(" --- lutEntry.covm[%d] is NaN", i);
      return false;
    }
  }
  return true;
}

void DelphesO2LutWriter::lutWrite(const char* filename, int pdg, float field, int itof, int otof)
{

  if (useFlatDipole && useDipole) {
    Printf("Both dipole and dipole flat flags are on, please use only one of them");
    return;
  }

  // output file
  std::ofstream lutFile(filename, std::ofstream::binary);
  if (!lutFile.is_open()) {
    Printf("Did not manage to open output file!!");
    return;
  }

  // write header
  lutHeader_t lutHeader;
  // pid
  lutHeader.pdg = pdg;
  lutHeader.mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  const int q = std::abs(TDatabasePDG::Instance()->GetParticle(pdg)->Charge()) / 3;
  if (q <= 0) {
    Printf("Negative or null charge (%f) for pdg code %i. Fix the charge!", TDatabasePDG::Instance()->GetParticle(pdg)->Charge(), pdg);
    return;
  }
  lutHeader.field = field;
  // nch
  lutHeader.nchmap.log = true;
  lutHeader.nchmap.nbins = 20;
  lutHeader.nchmap.min = 0.5;
  lutHeader.nchmap.max = 3.5;
  // radius
  lutHeader.radmap.log = false;
  lutHeader.radmap.nbins = 1;
  lutHeader.radmap.min = 0.;
  lutHeader.radmap.max = 100.;
  // eta
  lutHeader.etamap.log = false;
  lutHeader.etamap.nbins = 80;
  lutHeader.etamap.min = -4.;
  lutHeader.etamap.max = 4.;
  // pt
  lutHeader.ptmap.log = true;
  lutHeader.ptmap.nbins = 200;
  lutHeader.ptmap.min = -2;
  lutHeader.ptmap.max = 2.;
  lutFile.write(reinterpret_cast<char*>(&lutHeader), sizeof(lutHeader));

  // entries
  const int nnch = lutHeader.nchmap.nbins;
  const int nrad = lutHeader.radmap.nbins;
  const int neta = lutHeader.etamap.nbins;
  const int npt = lutHeader.ptmap.nbins;
  lutEntry_t lutEntry;

  // write entries
  int nCalls = 0;
  int successfullCalls = 0;
  int failedCalls = 0;
  for (int inch = 0; inch < nnch; ++inch) {
    Printf(" --- writing nch = %d/%d", inch, nnch);
    auto nch = lutHeader.nchmap.eval(inch);
    lutEntry.nch = nch;
    fat.SetdNdEtaCent(nch);
    for (int irad = 0; irad < nrad; ++irad) {
      Printf(" --- writing irad = %d/%d", irad, nrad);
      for (int ieta = 0; ieta < neta; ++ieta) {
        nCalls++;
        Printf(" --- writing ieta = %d/%d", ieta, neta);
        auto eta = lutHeader.etamap.eval(ieta);
        lutEntry.eta = lutHeader.etamap.eval(ieta);
        for (int ipt = 0; ipt < npt; ++ipt) {
          Printf(" --- writing ipt = %d/%d", ipt, npt);
          lutEntry.pt = lutHeader.ptmap.eval(ipt);
          lutEntry.valid = true;
          if (std::fabs(eta) <= etaMaxBarrel) { // full lever arm ends at etaMaxBarrel
            Printf("Solving in the barrel");
            // printf(" --- fatSolve: pt = %f, eta = %f, mass = %f, field=%f \n", lutEntry.pt, lutEntry.eta, lutHeader.mass, lutHeader.field);
            successfullCalls++;
            if (!fatSolve(lutEntry, lutEntry.pt, lutEntry.eta, lutHeader.mass, itof, otof, q)) {
              // printf(" --- fatSolve: error \n");
              lutEntry.valid = false;
              lutEntry.eff = 0.;
              lutEntry.eff2 = 0.;
              for (int i = 0; i < 15; ++i) {
                lutEntry.covm[i] = 0.;
              }
              successfullCalls--;
              failedCalls++;
            }
          } else {
            Printf("Solving outside the barrel");
            // printf(" --- fwdSolve: pt = %f, eta = %f, mass = %f, field=%f \n", lutEntry.pt, lutEntry.eta, lutHeader.mass, lutHeader.field);
            lutEntry.eff = 1.;
            lutEntry.eff2 = 1.;
            bool retval = true;
            successfullCalls++;
            if (useFlatDipole) { // Using the parametrization at the border of the barrel
              retval = fatSolve(lutEntry, lutEntry.pt, etaMaxBarrel, lutHeader.mass, itof, otof, q);
            } else if (usePara) {
              retval = fwdPara(lutEntry, lutEntry.pt, lutEntry.eta, lutHeader.mass, field);
            } else {
              retval = fwdSolve(lutEntry.covm, lutEntry.pt, lutEntry.eta, lutHeader.mass);
            }
            if (useDipole) { // Using the parametrization at the border of the barrel only for efficiency and momentum resolution
              lutEntry_t lutEntryBarrel;
              retval = fatSolve(lutEntryBarrel, lutEntry.pt, etaMaxBarrel, lutHeader.mass, itof, otof, q);
              lutEntry.valid = lutEntryBarrel.valid;
              lutEntry.covm[14] = lutEntryBarrel.covm[14];
              lutEntry.eff = lutEntryBarrel.eff;
              lutEntry.eff2 = lutEntryBarrel.eff2;
            }
            if (!retval) {
              printf(" --- fwdSolve: error \n");
              lutEntry.valid = false;
              for (int i = 0; i < 15; ++i) {
                lutEntry.covm[i] = 0.;
              }
              successfullCalls--;
              failedCalls++;
            }
          }
          Printf("Diagonalizing");
          diagonalise(lutEntry);
          Printf("Writing");
          lutFile.write(reinterpret_cast<char*>(&lutEntry), sizeof(lutEntry_t));
        }
      }
    }
  }
  Printf(" --- finished writing LUT file %s", filename);
  Printf(" --- successfull calls: %d/%d, failed calls: %d/%d", successfullCalls, nCalls, failedCalls, nCalls);
  lutFile.close();
}

void DelphesO2LutWriter::diagonalise(lutEntry_t& lutEntry)
{
  TMatrixDSym m(5);
  for (int i = 0, k = 0; i < 5; ++i) {
    for (int j = 0; j < i + 1; ++j, ++k) {
      m(i, j) = lutEntry.covm[k];
      m(j, i) = lutEntry.covm[k];
    }
  }

  m.Print();
  TMatrixDSymEigen eigen(m);
  // eigenvalues vector
  TVectorD eigenVal = eigen.GetEigenValues();
  for (int i = 0; i < 5; ++i)
    lutEntry.eigval[i] = eigenVal[i];
  // eigenvectors matrix
  TMatrixD eigenVec = eigen.GetEigenVectors();
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      lutEntry.eigvec[i][j] = eigenVec[i][j];
  // inverse eigenvectors matrix
  eigenVec.Invert();
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      lutEntry.eiginv[i][j] = eigenVec[i][j];
}

TGraph* DelphesO2LutWriter::lutRead(const char* filename, int pdg, int what, int vs, float nch, float radius, float eta, float pt)
{
  Printf(" --- reading LUT file %s", filename);
  // vs
  static const int kNch = 0;
  static const int kEta = 1;
  static const int kPt = 2;

  // what
  static const int kEfficiency = 0;
  static const int kEfficiency2 = 1;
  static const int kEfficiencyInnerTOF = 2;
  static const int kEfficiencyOuterTOF = 3;
  static const int kPtResolution = 4;
  static const int kRPhiResolution = 5;
  static const int kZResolution = 6;

  o2::delphes::DelphesO2TrackSmearer smearer;
  smearer.loadTable(pdg, filename);
  auto lutHeader = smearer.getLUTHeader(pdg);
  map_t lutMap;
  switch (vs) {
    case kNch:
      lutMap = lutHeader->nchmap;
      break;
    case kEta:
      lutMap = lutHeader->etamap;
      break;
    case kPt:
      lutMap = lutHeader->ptmap;
      break;
  }
  auto nbins = lutMap.nbins;
  auto g = new TGraph();
  g->SetName(Form("lut_%s_%d_vs_%d_what_%d", filename, pdg, vs, what));
  g->SetTitle(Form("LUT for %s, pdg %d, vs %d, what %d", filename, pdg, vs, what));
  switch (vs) {
    case kNch:
      Printf(" --- vs = kNch");
      g->GetXaxis()->SetTitle("Nch");
      break;
    case kEta:
      Printf(" --- vs = kEta");
      g->GetXaxis()->SetTitle("#eta");
      break;
    case kPt:
      Printf(" --- vs = kPt");
      g->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      break;
    default:
      Printf(" --- error: unknown vs %d", vs);
      return nullptr;
  }
  switch (what) {
    case kEfficiency:
      Printf(" --- what = kEfficiency");
      g->GetYaxis()->SetTitle("Efficiency (%)");
      break;
    case kEfficiency2:
      Printf(" --- what = kEfficiency2");
      g->GetYaxis()->SetTitle("Efficiency2 (%)");
      break;
    case kEfficiencyInnerTOF:
      Printf(" --- what = kEfficiencyInnerTOF");
      g->GetYaxis()->SetTitle("Inner TOF Efficiency (%)");
      break;
    case kEfficiencyOuterTOF:
      Printf(" --- what = kEfficiencyOuterTOF");
      g->GetYaxis()->SetTitle("Outer TOF Efficiency (%)");
      break;
    case kPtResolution:
      Printf(" --- what = kPtResolution");
      g->GetYaxis()->SetTitle("p_{T} Resolution (%)");
      break;
    case kRPhiResolution:
      Printf(" --- what = kRPhiResolution");
      g->GetYaxis()->SetTitle("R#phi Resolution (#mum)");
      break;
    case kZResolution:
      Printf(" --- what = kZResolution");
      g->GetYaxis()->SetTitle("Z Resolution (#mum)");
      break;
    default:
      Printf(" --- error: unknown what %d", what);
      return nullptr;
  }

  bool canBeInvalid = true;
  for (int i = 0; i < nbins; ++i) {
    switch (vs) {
      case kNch:
        nch = lutMap.eval(i);
        break;
      case kEta:
        eta = lutMap.eval(i);
        break;
      case kPt:
        pt = lutMap.eval(i);
        break;
    }
    float eff = 0.;
    auto lutEntry = smearer.getLUTEntry(pdg, nch, radius, eta, pt, eff);
    if (!lutEntry->valid || lutEntry->eff == 0.) {
      if (!canBeInvalid) {
        Printf(" --- warning: it cannot be invalid");
      }
      continue;
    }
    canBeInvalid = false;

    double cen = 0.;
    switch (vs) {
      case kNch:
        cen = lutEntry->nch;
        break;
      case kEta:
        cen = lutEntry->eta;
        break;
      case kPt:
        cen = lutEntry->pt;
        break;
    }
    double val = 0.;
    switch (what) {
      case kEfficiency:
        val = lutEntry->eff * 100.; // efficiency (%)
        break;
      case kEfficiency2:
        val = lutEntry->eff2 * 100.; // efficiency (%)
        break;
      case kEfficiencyInnerTOF:
        val = lutEntry->itof * 100.; // efficiency (%)
        break;
      case kEfficiencyOuterTOF:
        val = lutEntry->otof * 100.; // efficiency (%)
        break;
      case kPtResolution:
        val = sqrt(lutEntry->covm[14]) * lutEntry->pt * 100.; // pt resolution (%)
        break;
      case kRPhiResolution:
        val = sqrt(lutEntry->covm[0]) * 1.e4; // rphi resolution (um)
        break;
      case kZResolution:
        val = sqrt(lutEntry->covm[1]) * 1.e4; // z resolution (um)
        break;
      default:
        Printf(" --- error: unknown what %d", what);
        break;
    }
    g->AddPoint(cen, val);
  }

  return g;
}
} // namespace o2::fastsim

ClassImp(o2::fastsim::DelphesO2LutWriter);
