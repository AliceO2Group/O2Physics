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

#include "ALICE3/Core/FlatLutWriter.h"

#include "TrackUtilities.h"

#include "ALICE3/Core/FlatTrackSmearer.h"

#include <Framework/Logger.h>
#include <ReconstructionDataFormats/Track.h>

#include <TAxis.h>
#include <TDatabasePDG.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TMatrixDSymEigen.h>
#include <TMatrixDSymfwd.h>
#include <TMatrixDfwd.h>
#include <TParticlePDG.h>
#include <TString.h>
#include <TVectorDfwd.h>

#include <Rtypes.h>

#include <cstdlib>

using namespace o2::delphes;

namespace o2::fastsim
{
void FlatLutWriter::print() const
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

std::string FlatLutWriter::LutBinning::toString() const
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

bool FlatLutWriter::fatSolve(lutEntry_t& lutEntry,
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
  tlv.SetPtEtaPhiM(pt, eta, 0.f, mass);
  o2::track::TrackParCov trkIn;
  o2::upgrade::convertTLorentzVectorToO2Track(q, tlv, {0.f, 0.f, 0.f}, trkIn);

  o2::track::TrackParCov trkOut;
  const int status = fat.FastTrack(trkIn, trkOut, nch);
  if (status <= mAtLeastHits) {
    LOGF(debug, "fatSolve: FastTrack failed with status %d (threshold %d)", status, mAtLeastHits);
    return false;
  }

  LOGF(debug, "fatSolve: FastTrack succeeded with status %d", status);

  lutEntry.valid = true;
  lutEntry.itof = fat.GetGoodHitProb(itof);
  lutEntry.otof = fat.GetGoodHitProb(otof);

  static constexpr int nCov = 15;
  for (int i = 0; i < nCov; ++i)
    lutEntry.covm[i] = trkOut.getCov()[i];

  // Define the efficiency
  auto totfake = 0.f;
  lutEntry.eff = 1.f;
  for (size_t i = 1; i < fat.GetNLayers(); ++i) {
    if (fat.IsLayerInert(i))
      continue; // skip inert layers
    auto igoodhit = fat.GetGoodHitProb(i);
    if (igoodhit <= 0.f || i == itof || i == otof)
      continue;
    lutEntry.eff *= igoodhit;
    auto pairfake = 0.f;
    for (size_t j = i + 1; j < fat.GetNLayers(); ++j) {
      auto jgoodhit = fat.GetGoodHitProb(j);
      if (jgoodhit <= 0.f || j == itof || j == otof)
        continue;
      pairfake = (1.f - igoodhit) * (1.f - jgoodhit);
      break;
    }
    totfake += pairfake;
  }
  lutEntry.eff2 = (1.f - totfake);

  return true;
}

#ifdef USE_FWD_PARAM
bool FlatLutWriter::fwdSolve(float* covm, float pt, float eta, float mass)
{
  if (fwdRes(covm, pt, eta, mass) < 0)
    return false;
  return true;
}
#else
bool FlatLutWriter::fwdSolve(float*, float, float, float)
{
  return false;
}
#endif

bool FlatLutWriter::fwdPara(lutEntry_t& lutEntry, float pt, float eta, float mass, float Bfield)
{
  lutEntry.valid = false;

  // Parametrised forward response; interpolates between FAT at eta = 1.75 and a fixed parametrisation at eta = 4
  // Only diagonal elements
  static constexpr float etaLimit = 4.0f;
  if (std::fabs(eta) < etaMaxBarrel || std::fabs(eta) > etaLimit) {
    return false;
  }

  if (!fatSolve(lutEntry, pt, etaMaxBarrel, mass)) {
    return false;
  }

  static constexpr int nCov = 15;
  float covmbarrel[nCov] = {0.f};
  for (int i = 0; i < nCov; ++i) {
    covmbarrel[i] = lutEntry.covm[i];
  }

  // Parametrisation at eta = 4
  const double beta = 1. / std::sqrt(1.0 + mass * mass / pt / pt / std::cosh(eta) / std::cosh(eta));
  const float dcaPos = 2.5e-4f / std::sqrt(3.f); // 2.5 micron/sqrt(3)
  const float r0 = 0.5f;                         // layer 0 radius [cm]
  const float r1 = 1.3f;
  const float r2 = 2.5f;
  const float x0layer = 0.001f; // material budget (rad length) per layer

  const double sigmaAlpha = 0.0136 / beta / pt * std::sqrt(x0layer * std::cosh(eta)) * (1.0 + 0.038 * std::log(x0layer * std::cosh(eta)));
  const double dcaxyMs = sigmaAlpha * r0 * std::sqrt(1.0 + r1 * r1 / (r2 - r0) / (r2 - r0));
  const double dcaxy2 = dcaPos * dcaPos + dcaxyMs * dcaxyMs;

  const double dcazMs = sigmaAlpha * r0 * std::cosh(eta);
  const double dcaz2 = dcaPos * dcaPos + dcazMs * dcazMs;

  const float Leta = 2.8f / std::sinh(eta) - 0.01f * r0; // m
  const double relmomresPos = 10e-6 * pt / 0.3 / Bfield / Leta / Leta * std::sqrt(720.0 / 15.0);

  const float relmomresBarrel = std::sqrt(covmbarrel[14]) * pt;
  const float rOuter = 1.f; // m
  const float relmomresPosBarrel = 10e-6f * pt / 0.3f / Bfield / rOuter / rOuter / std::sqrt(720.f / 15.f);
  const float relmomresMSBarrel = std::sqrt(relmomresBarrel * relmomresBarrel - relmomresPosBarrel * relmomresPosBarrel);

  // Interpolate MS contribution (rel resolution 0.4 at eta = 4)
  const float relmomresMSEta4 = 0.4f / beta * 0.5f / Bfield;
  const float relmomresMS = relmomresMSEta4 * std::pow(relmomresMSEta4 / relmomresMSBarrel, (std::fabs(eta) - 4.f) / (4.f - etaMaxBarrel));
  const float momresTot = pt * std::sqrt(relmomresPos * relmomresPos + relmomresMS * relmomresMS); // total absolute mom reso

  // Fill cov matrix diag
  for (int i = 0; i < 15; ++i) {
    lutEntry.covm[i] = 0.f;
  }

  lutEntry.covm[0] = covmbarrel[0];
  if (dcaxy2 > lutEntry.covm[0]) {
    lutEntry.covm[0] = dcaxy2;
  }
  lutEntry.covm[2] = covmbarrel[2];
  if (dcaz2 > lutEntry.covm[2]) {
    lutEntry.covm[2] = dcaz2;
  }
  lutEntry.covm[5] = covmbarrel[5];                              // sigma^2 sin(phi)
  lutEntry.covm[9] = covmbarrel[9];                              // sigma^2 tanl
  lutEntry.covm[14] = momresTot * momresTot / pt / pt / pt / pt; // sigma^2 1/pt

  // Check that all numbers are numbers
  for (int i = 0; i < 15; ++i) {
    if (std::isnan(lutEntry.covm[i])) {
      LOGF(info, "lutEntry.covm[%d] is NaN", i);
      return false;
    }
  }
  return true;
}

void FlatLutWriter::lutWrite(const char* filename, int pdg, float field, size_t itof, size_t otof)
{
  if (useFlatDipole && useDipole) {
    LOGF(error, "Both dipole and dipole flat flags are on, please use only one of them");
    return;
  }

  // Open output file for binary writing
  std::ofstream lutFile(filename, std::ofstream::binary);
  if (!lutFile.is_open()) {
    LOGF(error, "Failed to open output file: %s", filename);
    return;
  }

  // Create and write header
  lutHeader_t lutHeader;
  lutHeader.pdg = pdg;

  const TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdg);
  if (!particle) {
    LOGF(fatal, "Cannot find particle with PDG code %d", pdg);
    lutFile.close();
    return;
  }

  lutHeader.mass = particle->Mass();
  const int q = std::abs(particle->Charge()) / 3;
  if (q <= 0) {
    LOGF(error, "Negative or null charge (%f) for pdg code %d", particle->Charge(), pdg);
    lutFile.close();
    return;
  }

  lutHeader.field = field;

  // Set binning maps in header
  auto setMap = [](map_t& map, const LutBinning& b) {
    map.log = b.log;
    map.nbins = b.nbins;
    map.min = b.min;
    map.max = b.max;
  };

  setMap(lutHeader.nchmap, mNchBinning);
  setMap(lutHeader.radmap, mRadiusBinning);
  setMap(lutHeader.etamap, mEtaBinning);
  setMap(lutHeader.ptmap, mPtBinning);

  // Write header to file
  lutFile.write(reinterpret_cast<char*>(&lutHeader), sizeof(lutHeader_t));
  if (!lutFile.good()) {
    LOGF(error, "Failed to write LUT header to %s", filename);
    lutFile.close();
    return;
  }

  // Get dimensions
  const int nnch = lutHeader.nchmap.nbins;
  const int nrad = lutHeader.radmap.nbins;
  const int neta = lutHeader.etamap.nbins;
  const int npt = lutHeader.ptmap.nbins;

  LOGF(info, "Writing LUT with dimensions: nch=%d rad=%d eta=%d pt=%d (total=%zu entries)", nnch, nrad, neta, npt, static_cast<size_t>(nnch) * nrad * neta * npt);

  lutEntry_t lutEntry;
  int nCalls = 0;
  int successfulCalls = 0;
  int failedCalls = 0;

  // Write all entries sequentially
  for (int inch = 0; inch < nnch; ++inch) {
    LOGF(info, "Writing nch bin %d/%d", inch, nnch);
    auto nch = lutHeader.nchmap.eval(inch);
    lutEntry.nch = nch;
    fat.SetdNdEtaCent(nch);

    for (int irad = 0; irad < nrad; ++irad) {
      for (int ieta = 0; ieta < neta; ++ieta) {
        auto eta = lutHeader.etamap.eval(ieta);
        lutEntry.eta = eta;

        for (int ipt = 0; ipt < npt; ++ipt) {
          nCalls++;
          lutEntry.pt = lutHeader.ptmap.eval(ipt);
          lutEntry.valid = true;

          // Solve for this bin
          if (std::fabs(eta) <= etaMaxBarrel) {
            // Full lever arm region (barrel)
            LOGF(debug, "Solving barrel: pt=%f eta=%f", lutEntry.pt, lutEntry.eta);
            successfulCalls++;

            if (!fatSolve(lutEntry, lutEntry.pt, lutEntry.eta, lutHeader.mass, itof, otof, q, nch)) {
              lutEntry.valid = false;
              lutEntry.eff = 0.f;
              lutEntry.eff2 = 0.f;
              for (int i = 0; i < 15; ++i) {
                lutEntry.covm[i] = 0.f;
              }
              successfulCalls--;
              failedCalls++;
            }
          } else {
            // Forward region
            LOGF(debug, "Solving forward: pt=%f eta=%f", lutEntry.pt, lutEntry.eta);
            lutEntry.eff = 1.f;
            lutEntry.eff2 = 1.f;
            bool retval = true;
            successfulCalls++;

            if (useFlatDipole) {
              // Using the parametrization at the border of the barrel
              retval = fatSolve(lutEntry, lutEntry.pt, etaMaxBarrel, lutHeader.mass, itof, otof, q, nch);
            } else if (usePara) {
              retval = fwdPara(lutEntry, lutEntry.pt, lutEntry.eta, lutHeader.mass, field);
            } else {
              retval = fwdSolve(lutEntry.covm, lutEntry.pt, lutEntry.eta, lutHeader.mass);
            }

            if (useDipole) {
              // Using the parametrization at the border of the barrel only for efficiency and momentum resolution
              lutEntry_t lutEntryBarrel;
              retval = fatSolve(lutEntryBarrel, lutEntry.pt, etaMaxBarrel, lutHeader.mass, itof, otof, q, nch);
              lutEntry.valid = lutEntryBarrel.valid;
              lutEntry.covm[14] = lutEntryBarrel.covm[14];
              lutEntry.eff = lutEntryBarrel.eff;
              lutEntry.eff2 = lutEntryBarrel.eff2;
            }

            if (!retval) {
              LOGF(debug, "Forward solve failed");
              lutEntry.valid = false;
              for (int i = 0; i < 15; ++i) {
                lutEntry.covm[i] = 0.f;
              }
              successfulCalls--;
              failedCalls++;
            }
          }

          // Diagonalize covariance matrix
          diagonalise(lutEntry);

          // Write entry to file
          lutFile.write(reinterpret_cast<char*>(&lutEntry), sizeof(lutEntry_t));
          if (!lutFile.good()) {
            LOGF(error, "Failed to write LUT entry at index [%d,%d,%d,%d]", inch, irad, ieta, ipt);
            lutFile.close();
            return;
          }
        }
      }
    }
  }

  lutFile.close();

  LOGF(info, "Finished writing LUT file: %s", filename);
  LOGF(info, "Total calls: %d, successful: %d, failed: %d", nCalls, successfulCalls, failedCalls);
}

void FlatLutWriter::diagonalise(lutEntry_t& lutEntry)
{
  static constexpr int kEig = 5;
  TMatrixDSym m(kEig);

  for (int i = 0, k = 0; i < kEig; ++i) {
    for (int j = 0; j < i + 1; ++j, ++k) {
      m(i, j) = lutEntry.covm[k];
      m(j, i) = lutEntry.covm[k];
    }
  }

  TMatrixDSymEigen eigen(m);

  // Eigenvalues
  TVectorD eigenVal = eigen.GetEigenValues();
  for (int i = 0; i < kEig; ++i)
    lutEntry.eigval[i] = eigenVal[i];

  // Eigenvectors
  TMatrixD eigenVec = eigen.GetEigenVectors();
  for (int i = 0; i < kEig; ++i)
    for (int j = 0; j < kEig; ++j)
      lutEntry.eigvec[i][j] = eigenVec[i][j];

  // Inverse eigenvectors
  eigenVec.Invert();
  for (int i = 0; i < kEig; ++i)
    for (int j = 0; j < kEig; ++j)
      lutEntry.eiginv[i][j] = eigenVec[i][j];
}

TGraph* FlatLutWriter::lutRead(const char* filename, int pdg, int what, int vs, float nch, float radius, float eta, float pt)
{
  LOGF(info, "Reading LUT file: %s", filename);

  static const int kNch = 0;
  static const int kEta = 1;
  static const int kPt = 2;

  static const int kEfficiency = 0;
  static const int kEfficiency2 = 1;
  static const int kEfficiencyInnerTOF = 2;
  static const int kEfficiencyOuterTOF = 3;
  static const int kPtResolution = 4;
  static const int kRPhiResolution = 5;
  static const int kZResolution = 6;

  // Use TrackSmearer to load and access the LUT
  o2::delphes::TrackSmearer smearer;
  if (!smearer.loadTable(pdg, filename)) {
    LOGF(error, "Failed to load LUT from %s", filename);
    return nullptr;
  }

  const lutHeader_t* lutHeader = smearer.getLUTHeader(pdg);
  if (!lutHeader) {
    LOGF(error, "No LUT header for PDG %d", pdg);
    return nullptr;
  }

  lutHeader->print();

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
    default:
      LOGF(error, "Unknown vs: %d", vs);
      return nullptr;
  }

  auto nbins = lutMap.nbins;
  auto g = new TGraph();
  g->SetName(Form("lut_%s_%d_vs_%d_what_%d", filename, pdg, vs, what));
  g->SetTitle(Form("LUT for %s, pdg %d, vs %d, what %d", filename, pdg, vs, what));

  // Set axis labels
  switch (vs) {
    case kNch:
      LOGF(info, "Plot versus Nch");
      g->GetXaxis()->SetTitle("Nch");
      break;
    case kEta:
      LOGF(info, "Plot versus Eta");
      g->GetXaxis()->SetTitle("#eta");
      break;
    case kPt:
      LOGF(info, "Plot versus Pt");
      g->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      break;
  }

  switch (what) {
    case kEfficiency:
      LOGF(info, "Plot efficiency");
      g->GetYaxis()->SetTitle("Efficiency (%)");
      break;
    case kEfficiency2:
      LOGF(info, "Plot efficiency2");
      g->GetYaxis()->SetTitle("Efficiency2 (%)");
      break;
    case kEfficiencyInnerTOF:
      LOGF(info, "Plot inner TOF efficiency");
      g->GetYaxis()->SetTitle("Inner TOF Efficiency (%)");
      break;
    case kEfficiencyOuterTOF:
      LOGF(info, "Plot outer TOF efficiency");
      g->GetYaxis()->SetTitle("Outer TOF Efficiency (%)");
      break;
    case kPtResolution:
      LOGF(info, "Plot pt resolution");
      g->GetYaxis()->SetTitle("p_{T} Resolution (%)");
      break;
    case kRPhiResolution:
      LOGF(info, "Plot rphi resolution");
      g->GetYaxis()->SetTitle("R#phi Resolution (#mum)");
      break;
    case kZResolution:
      LOGF(info, "Plot z resolution");
      g->GetYaxis()->SetTitle("Z Resolution (#mum)");
      break;
    default:
      LOGF(error, "Unknown what: %d", what);
      delete g;
      return nullptr;
  }

  // Fill graph by iterating over one dimension
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

    float eff = 0.f;
    const lutEntry_t* lutEntry = smearer.getLUTEntry(pdg, nch, radius, eta, pt, eff);

    if (!lutEntry || !lutEntry->valid || lutEntry->eff == 0.f) {
      if (!canBeInvalid) {
        LOGF(warning, "Entry became invalid at bin %d", i);
      }
      continue;
    }
    canBeInvalid = false;

    double cen = 0.f;
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

    double val = 0.f;
    switch (what) {
      case kEfficiency:
        val = lutEntry->eff * 100.f; // efficiency (%)
        break;
      case kEfficiency2:
        val = lutEntry->eff2 * 100.f; // efficiency (%)
        break;
      case kEfficiencyInnerTOF:
        val = lutEntry->itof * 100.f; // efficiency (%)
        break;
      case kEfficiencyOuterTOF:
        val = lutEntry->otof * 100.f; // efficiency (%)
        break;
      case kPtResolution:
        val = std::sqrt(lutEntry->covm[14]) * lutEntry->pt * 100.f; // pt resolution (%)
        break;
      case kRPhiResolution:
        val = std::sqrt(lutEntry->covm[0]) * 1.e4f; // rphi resolution (um)
        break;
      case kZResolution:
        val = std::sqrt(lutEntry->covm[1]) * 1.e4f; // z resolution (um)
        break;
      default:
        LOGF(error, "Unknown what: %d", what);
        break;
    }
    g->AddPoint(cen, val);
  }

  LOGF(info, "Read %d points from LUT", g->GetN());
  return g;
}

} // namespace o2::fastsim
