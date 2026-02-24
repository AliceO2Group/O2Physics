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

#include "ALICE3/Core/DelphesO2LutWriter.h"

#include "ALICE3/Core/DelphesO2TrackSmearer.h"
#include "ALICE3/Core/FastTracker.h"
#include "ALICE3/Core/TrackUtilities.h"

#include "TAxis.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"

#include <cstdio>
#include <string>

// #define USE_FWD_PARAM
#ifdef USE_FWD_PARAM
#include "fwdRes.C"
#endif

namespace o2::fastsim
{

void DelphesO2LutWriter::print() const
{
  LOG(info) << " --- Printing configuration of LUT writer --- ";
  LOG(info) << "    -> etaMaxBarrel  = " << etaMaxBarrel;
  LOG(info) << "    -> usePara       = " << usePara;
  LOG(info) << "    -> useDipole     = " << useDipole;
  LOG(info) << "    -> useFlatDipole = " << useFlatDipole;
  LOG(info) << "    -> mAtLeastHits  = " << mAtLeastHits;
  LOG(info) << "    -> mAtLeastCorr  = " << mAtLeastCorr;
  LOG(info) << "    -> mAtLeastFake  = " << mAtLeastFake;
  LOG(info) << "    -> Nch Binning: = " << mNchBinning.toString();
  LOG(info) << "    -> Radius Binning: = " << mRadiusBinning.toString();
  LOG(info) << "    -> Eta Binning: = " << mEtaBinning.toString();
  LOG(info) << "    -> Pt Binning: = " << mPtBinning.toString();
  LOG(info) << " --- End of configuration --- ";
}

std::string DelphesO2LutWriter::LutBinning::toString() const
{
  std::string str = "";
  str.append(log ? "log" : "lin");
  str.append(" nbins: ");
  str.append(std::to_string(nbins));
  str.append(" min: ");
  str.append(std::to_string(min));
  str.append(" max: ");
  str.append(std::to_string(max));
  return str;
}

bool DelphesO2LutWriter::fatSolve(o2::delphes::DelphesO2TrackSmearer::lutEntry_t& lutEntry,
                                  float pt,
                                  float eta,
                                  const float mass,
                                  size_t itof,
                                  size_t otof,
                                  int q,
                                  const float nch)
{
  lutEntry.valid = false;

  static TLorentzVector tlv;
  tlv.SetPtEtaPhiM(pt, eta, 0., mass);
  o2::track::TrackParCov trkIn;
  o2::upgrade::convertTLorentzVectorToO2Track(q, tlv, {0., 0., 0.}, trkIn);
  // tlv.Print();
  // return fmt::format("X:{:+.4e} Alp:{:+.3e} Par: {:+.4e} {:+.4e} {:+.4e} {:+.4e} {:+.4e} |Q|:{:d} {:s}\n",
  //  getX(), getAlpha(), getY(), getZ(), getSnp(), getTgl(), getQ2Pt(), getAbsCharge(), getPID().getName());
  // trkIn.print();
  o2::track::TrackParCov trkOut;
  const int status = fat.FastTrack(trkIn, trkOut, nch);
  if (status <= mAtLeastHits) {
    LOGF(info, " --- fatSolve: FastTrack failed ---");
    // tlv.Print();
    return false;
  }
  LOGF(info, " --- fatSolve: FastTrack succeeded %d ---", status);
  // trkOut.print();
  lutEntry.valid = true;
  lutEntry.itof = fat.GetGoodHitProb(itof);
  lutEntry.otof = fat.GetGoodHitProb(otof);
  static constexpr int nCov = 15;
  for (int i = 0; i < nCov; ++i)
    lutEntry.covm[i] = trkOut.getCov()[i];

  // define the efficiency
  auto totfake = 0.;
  lutEntry.eff = 1.;
  for (size_t i = 1; i < fat.GetNLayers(); ++i) {
    if (fat.IsLayerInert(i))
      continue; // skip inert layers
    auto igoodhit = fat.GetGoodHitProb(i);
    if (igoodhit <= 0. || i == itof || i == otof)
      continue;
    lutEntry.eff *= igoodhit;
    auto pairfake = 0.;
    for (size_t j = i + 1; j < fat.GetNLayers(); ++j) {
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

bool DelphesO2LutWriter::fwdPara(o2::delphes::DelphesO2TrackSmearer::lutEntry_t& lutEntry, float pt, float eta, float mass, float Bfield)
{
  lutEntry.valid = false;

  // parametrised forward response; interpolates between FAT at eta = 1.75 and a fixed parametrisation at eta = 4; only diagonal elements
  static constexpr float etaLimit = 4.0f;
  if (std::fabs(eta) < etaMaxBarrel || std::fabs(eta) > etaLimit)
    return false;

  if (!fatSolve(lutEntry, pt, etaMaxBarrel, mass))
    return false;
  static constexpr int nCov = 15;
  float covmbarrel[nCov] = {0};
  for (int i = 0; i < nCov; ++i) {
    covmbarrel[i] = lutEntry.covm[i];
  }

  // parametrisation at eta = 4
  const double beta = 1. / std::sqrt(1 + mass * mass / pt / pt / std::cosh(eta) / std::cosh(eta));
  const float dcaPos = 2.5e-4 / std::sqrt(3); // 2.5 micron/sqrt(3)
  const float r0 = 0.5;                       // layer 0 radius [cm]
  const float r1 = 1.3;
  const float r2 = 2.5;
  const float x0layer = 0.001; // material budget (rad length) per layer
  const double sigmaAlpha = 0.0136 / beta / pt * std::sqrt(x0layer * std::cosh(eta)) * (1 + 0.038 * std::log(x0layer * std::cosh(eta)));
  const double dcaxyMs = sigmaAlpha * r0 * std::sqrt(1 + r1 * r1 / (r2 - r0) / (r2 - r0));
  const double dcaxy2 = dcaPos * dcaPos + dcaxyMs * dcaxyMs;

  const double dcazMs = sigmaAlpha * r0 * std::cosh(eta);
  const double dcaz2 = dcaPos * dcaPos + dcazMs * dcazMs;

  const float Leta = 2.8 / std::sinh(eta) - 0.01 * r0; // m
  const double relmomresPos = 10e-6 * pt / 0.3 / Bfield / Leta / Leta * std::sqrt(720. / 15.);

  const float relmomresBarrel = std::sqrt(covmbarrel[14]) * pt;
  const float rOuter = 1; // m
  const float relmomresPosBarrel = 10e-6 * pt / 0.3 / Bfield / rOuter / rOuter / std::sqrt(720. / 15.);
  const float relmomresMSBarrel = std::sqrt(relmomresBarrel * relmomresBarrel - relmomresPosBarrel * relmomresPosBarrel);

  // interpolate MS contrib (rel resolution 0.4 at eta = 4)
  const float relmomresMSEta4 = 0.4 / beta * 0.5 / Bfield;
  const float relmomresMS = relmomresMSEta4 * std::pow(relmomresMSEta4 / relmomresMSBarrel, (std::fabs(eta) - 4.) / (4. - etaMaxBarrel));
  const float momresTot = pt * std::sqrt(relmomresPos * relmomresPos + relmomresMS * relmomresMS); // total absolute mom reso

  // Fill cov matrix diag
  for (int i = 0; i < 15; ++i)
    lutEntry.covm[i] = 0;

  lutEntry.covm[0] = covmbarrel[0];
  if (dcaxy2 > lutEntry.covm[0])
    lutEntry.covm[0] = dcaxy2;
  lutEntry.covm[2] = covmbarrel[2];
  if (dcaz2 > lutEntry.covm[2])
    lutEntry.covm[2] = dcaz2;
  lutEntry.covm[5] = covmbarrel[5];                              // sigma^2 sin(phi)
  lutEntry.covm[9] = covmbarrel[9];                              // sigma^2 tanl
  lutEntry.covm[14] = momresTot * momresTot / pt / pt / pt / pt; // sigma^2 1/pt
  // Check that all numbers are numbers
  for (int i = 0; i < 15; ++i) {
    if (std::isnan(lutEntry.covm[i])) {
      LOGF(info, " --- lutEntry.covm[%d] is NaN", i);
      return false;
    }
  }
  return true;
}

void DelphesO2LutWriter::lutWrite(const char* filename, int pdg, float field, size_t itof, size_t otof)
{

  if (useFlatDipole && useDipole) {
    LOGF(info, "Both dipole and dipole flat flags are on, please use only one of them");
    return;
  }

  // output file
  std::ofstream lutFile(filename, std::ofstream::binary);
  if (!lutFile.is_open()) {
    LOGF(info, "Did not manage to open output file!!");
    return;
  }

  // write header
  o2::delphes::DelphesO2TrackSmearer::lutHeader_t lutHeader;
  // pid
  lutHeader.pdg = pdg;
  const TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdg);
  if (!particle) {
    LOG(fatal) << "Cannot find particle with PDG code " << pdg;
    return;
  }
  lutHeader.mass = particle->Mass();
  const int q = std::abs(particle->Charge()) / 3;
  if (q <= 0) {
    LOGF(info, "Negative or null charge (%f) for pdg code %i. Fix the charge!", particle->Charge(), pdg);
    return;
  }
  lutHeader.field = field;
  auto setMap = [](o2::delphes::DelphesO2TrackSmearer::map_t& map, LutBinning b) {
    map.log = b.log;
    map.nbins = b.nbins;
    map.min = b.min;
    map.max = b.max;
  };
  // nch
  setMap(lutHeader.nchmap, mNchBinning);
  // radius
  setMap(lutHeader.radmap, mRadiusBinning);
  // eta
  setMap(lutHeader.etamap, mEtaBinning);
  // pt
  setMap(lutHeader.ptmap, mPtBinning);

  lutFile.write(reinterpret_cast<char*>(&lutHeader), sizeof(lutHeader));

  // entries
  const int nnch = lutHeader.nchmap.nbins;
  const int nrad = lutHeader.radmap.nbins;
  const int neta = lutHeader.etamap.nbins;
  const int npt = lutHeader.ptmap.nbins;
  o2::delphes::DelphesO2TrackSmearer::lutEntry_t lutEntry;

  // write entries
  int nCalls = 0;
  int successfullCalls = 0;
  int failedCalls = 0;
  for (int inch = 0; inch < nnch; ++inch) {
    LOGF(info, " --- writing nch = %d/%d", inch, nnch);
    auto nch = lutHeader.nchmap.eval(inch);
    lutEntry.nch = nch;
    fat.SetdNdEtaCent(nch);
    for (int irad = 0; irad < nrad; ++irad) {
      LOGF(info, " --- writing irad = %d/%d", irad, nrad);
      for (int ieta = 0; ieta < neta; ++ieta) {
        LOGF(info, " --- writing ieta = %d/%d", ieta, neta);
        auto eta = lutHeader.etamap.eval(ieta);
        lutEntry.eta = lutHeader.etamap.eval(ieta);
        for (int ipt = 0; ipt < npt; ++ipt) {
          nCalls++;
          LOGF(info, " --- writing ipt = %d/%d", ipt, npt);
          lutEntry.pt = lutHeader.ptmap.eval(ipt);
          lutEntry.valid = true;
          if (std::fabs(eta) <= etaMaxBarrel) { // full lever arm ends at etaMaxBarrel
            LOGF(info, "Solving in the barrel");
            // LOGF(info, " --- fatSolve: pt = %f, eta = %f, mass = %f, field=%f", lutEntry.pt, lutEntry.eta, lutHeader.mass, lutHeader.field);
            successfullCalls++;
            if (!fatSolve(lutEntry, lutEntry.pt, lutEntry.eta, lutHeader.mass, itof, otof, q)) {
              // LOGF(info, " --- fatSolve: error");
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
            LOGF(info, "Solving outside the barrel");
            // LOGF(info, " --- fwdSolve: pt = %f, eta = %f, mass = %f, field=%f", lutEntry.pt, lutEntry.eta, lutHeader.mass, lutHeader.field);
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
              o2::delphes::DelphesO2TrackSmearer::lutEntry_t lutEntryBarrel;
              retval = fatSolve(lutEntryBarrel, lutEntry.pt, etaMaxBarrel, lutHeader.mass, itof, otof, q);
              lutEntry.valid = lutEntryBarrel.valid;
              lutEntry.covm[14] = lutEntryBarrel.covm[14];
              lutEntry.eff = lutEntryBarrel.eff;
              lutEntry.eff2 = lutEntryBarrel.eff2;
            }
            if (!retval) {
              LOGF(info, " --- fwdSolve: error");
              lutEntry.valid = false;
              for (int i = 0; i < 15; ++i) {
                lutEntry.covm[i] = 0.;
              }
              successfullCalls--;
              failedCalls++;
            }
          }
          LOGF(info, "Diagonalizing");
          diagonalise(lutEntry);
          LOGF(info, "Writing");
          lutFile.write(reinterpret_cast<char*>(&lutEntry), sizeof(o2::delphes::DelphesO2TrackSmearer::lutEntry_t));
        }
      }
    }
  }
  LOGF(info, " --- finished writing LUT file %s", filename);
  LOGF(info, " --- successfull calls: %d/%d, failed calls: %d/%d", successfullCalls, nCalls, failedCalls, nCalls);
  lutFile.close();
}

void DelphesO2LutWriter::diagonalise(o2::delphes::DelphesO2TrackSmearer::lutEntry_t& lutEntry)
{
  static constexpr int kEig = 5;
  TMatrixDSym m(kEig);
  for (int i = 0, k = 0; i < kEig; ++i) {
    for (int j = 0; j < i + 1; ++j, ++k) {
      m(i, j) = lutEntry.covm[k];
      m(j, i) = lutEntry.covm[k];
    }
  }

  // m.Print();
  TMatrixDSymEigen eigen(m);
  // eigenvalues vector
  TVectorD eigenVal = eigen.GetEigenValues();
  for (int i = 0; i < kEig; ++i)
    lutEntry.eigval[i] = eigenVal[i];
  // eigenvectors matrix
  TMatrixD eigenVec = eigen.GetEigenVectors();
  for (int i = 0; i < kEig; ++i)
    for (int j = 0; j < kEig; ++j)
      lutEntry.eigvec[i][j] = eigenVec[i][j];
  // inverse eigenvectors matrix
  eigenVec.Invert();
  for (int i = 0; i < kEig; ++i)
    for (int j = 0; j < kEig; ++j)
      lutEntry.eiginv[i][j] = eigenVec[i][j];
}

TGraph* DelphesO2LutWriter::lutRead(const char* filename, int pdg, int what, int vs, float nch, float radius, float eta, float pt)
{
  LOGF(info, " --- reading LUT file %s", filename);
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
  lutHeader->print();
  o2::delphes::DelphesO2TrackSmearer::map_t lutMap;
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
      LOGF(info, " --- vs = kNch");
      g->GetXaxis()->SetTitle("Nch");
      break;
    case kEta:
      LOGF(info, " --- vs = kEta");
      g->GetXaxis()->SetTitle("#eta");
      break;
    case kPt:
      LOGF(info, " --- vs = kPt");
      g->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      break;
    default:
      LOGF(info, " --- error: unknown vs %d", vs);
      return nullptr;
  }
  switch (what) {
    case kEfficiency:
      LOGF(info, " --- what = kEfficiency");
      g->GetYaxis()->SetTitle("Efficiency (%)");
      break;
    case kEfficiency2:
      LOGF(info, " --- what = kEfficiency2");
      g->GetYaxis()->SetTitle("Efficiency2 (%)");
      break;
    case kEfficiencyInnerTOF:
      LOGF(info, " --- what = kEfficiencyInnerTOF");
      g->GetYaxis()->SetTitle("Inner TOF Efficiency (%)");
      break;
    case kEfficiencyOuterTOF:
      LOGF(info, " --- what = kEfficiencyOuterTOF");
      g->GetYaxis()->SetTitle("Outer TOF Efficiency (%)");
      break;
    case kPtResolution:
      LOGF(info, " --- what = kPtResolution");
      g->GetYaxis()->SetTitle("p_{T} Resolution (%)");
      break;
    case kRPhiResolution:
      LOGF(info, " --- what = kRPhiResolution");
      g->GetYaxis()->SetTitle("R#phi Resolution (#mum)");
      break;
    case kZResolution:
      LOGF(info, " --- what = kZResolution");
      g->GetYaxis()->SetTitle("Z Resolution (#mum)");
      break;
    default:
      LOGF(info, " --- error: unknown what %d", what);
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
        LOGF(info, " --- warning: it cannot be invalid");
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
        val = std::sqrt(lutEntry->covm[14]) * lutEntry->pt * 100.; // pt resolution (%)
        break;
      case kRPhiResolution:
        val = std::sqrt(lutEntry->covm[0]) * 1.e4; // rphi resolution (um)
        break;
      case kZResolution:
        val = std::sqrt(lutEntry->covm[1]) * 1.e4; // z resolution (um)
        break;
      default:
        LOGF(info, " --- error: unknown what %d", what);
        break;
    }
    g->AddPoint(cen, val);
  }

  return g;
}
} // namespace o2::fastsim

ClassImp(o2::fastsim::DelphesO2LutWriter);
