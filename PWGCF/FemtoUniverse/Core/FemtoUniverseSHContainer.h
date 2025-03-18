// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoUniverseSHContainer.h
/// \brief FemtoUniverseSHContainer - Fills the Spherical Harmonics components
/// \remark This file is inherited from ~/FemtoUniverse/Core/FemtoUniverse3DContainer.h on 17/06/2024
/// \author Pritam Chakraborty, WUT Warsaw, pritam.chakraborty@pw.edu.pl8

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSESHCONTAINER_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSESHCONTAINER_H_

#include <fairlogger/Logger.h>
#include <vector>
#include <string>
#include <complex>
#include <memory>

#include "Framework/HistogramRegistry.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseMath.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseSpherHarMath.h"
#include "Math/Vector4D.h"
#include "TMath.h"
#include "TDatabasePDG.h"

namespace o2::analysis::femto_universe
{

namespace femto_universe_sh_container
{
/// Femtoscopic observable to be computed
enum Observable { kstar ///< kstar
};

/// Type of the event processind
enum EventType { same, ///< Pair from same event
                 mixed ///< Pair from mixed event
};
}; // namespace femto_universe_sh_container

/// \class femto_universe_sh_container
/// \brief Container for all histogramming related to the correlation function. The two
/// particles of the pair are passed here, and the correlation function and QA histograms
/// are filled according to the specified observable
/// \tparam eventType Type of the event (same/mixed)
/// \tparam obs Observable to be computed (k*/Q_inv/...)
template <femto_universe_sh_container::EventType eventType, femto_universe_sh_container::Observable obs>
class FemtoUniverseSHContainer
{
 public:
  /// Destructor
  virtual ~FemtoUniverseSHContainer() = default;

  /// Container for all histogramming related to the spherical harmonics of the correlation function. The two
  /// particles of the pair are passed here, and the correlation function
  /// are filled according to the specified observable
  /// \tparam T type of the configurable for the axis configuration
  /// \param registry Histogram registry to be passed
  /// \param kstarbins k* binning for the histograms
  template <typename T>
  void init(HistogramRegistry* registry, T& kstarbins, int /*maxl*/)
  {
    kStarBins = kstarbins;
    std::string femtoObs1D;

    std::vector<int> fels(kMaxJM);
    std::vector<int> fems(kMaxJM);
    std::vector<int> felsi(kMaxJM);
    std::vector<int> femsi(kMaxJM);

    // Fill in els and ems table
    int el = 0;
    int em = 0;
    int il = 0;
    do {

      fels[il] = el;
      fems[il] = em;
      felsi[il] = static_cast<int>(el);
      femsi[il] = static_cast<int>(em);
      em++;
      il++;
      if (em > el) {
        el++;
        em = -el;
      }
    } while (el <= kMaxL);

    kHistogramRegistry = registry;
    femtoObs1D = "#it{q} (GeV/#it{c})";

    framework::AxisSpec femtoObsAxis1D = {kstarbins, femtoObs1D.c_str()};
    std::string folderName = static_cast<std::string>(kFolderSuffix[kEventType]) + static_cast<std::string>(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kRecon]);
    std::string suffix;

    for (int ihist = 0; ihist < kMaxJM; ihist++) {
      if (femsi[ihist] < 0) {
        suffix = "Ylm" + std::to_string(felsi[ihist]) + std::to_string(felsi[ihist] - femsi[ihist]);
      } else {
        suffix = "Ylm" + std::to_string(felsi[ihist]) + std::to_string(femsi[ihist]);
      }

      if (kFolderSuffix[kEventType] == kFolderSuffix[0]) {
        fnumsreal[ihist] = kHistogramRegistry->add<TH1>(("NumRe" + suffix).c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
        fnumsimag[ihist] = kHistogramRegistry->add<TH1>(("NumIm" + suffix).c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
      } else {
        fdensreal[ihist] = kHistogramRegistry->add<TH1>(("DenRe" + suffix).c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
        fdensimag[ihist] = kHistogramRegistry->add<TH1>(("DenIm" + suffix).c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
      }
    }

    if (kFolderSuffix[kEventType] == kFolderSuffix[0]) {
      std::string bufnameNum = "CovNum";
      fcovnum = kHistogramRegistry->add<TH3>((bufnameNum).c_str(), "; x; y; z", kTH3D, {{kstarbins}, {(kMaxJM * 2), -0.5, ((static_cast<float>(kMaxJM) * 2.0 - 0.5))}, {(kMaxJM * 2), -0.5, ((static_cast<float>(kMaxJM) * 2.0 - 0.5))}});
    } else if (kFolderSuffix[kEventType] == kFolderSuffix[1]) {
      std::string bufnameDen = "CovDen";
      fcovden = kHistogramRegistry->add<TH3>((bufnameDen).c_str(), "; x; y; z", kTH3D, {{kstarbins}, {(kMaxJM * 2), -0.5, ((static_cast<float>(kMaxJM) * 2.0 - 0.5))}, {(kMaxJM * 2), -0.5, ((static_cast<float>(kMaxJM) * 2.0 - 0.5))}});
    }

    fbinctn = new TH1D(TString("BinCountNum"), "Bin Occupation (Numerator)", static_cast<int>(kStarBins[0]), kStarBins[1], kStarBins[2]);
    fbinctd = new TH1D(TString("BinCountDen"), "Bin Occupation (Denominator)", static_cast<int>(kStarBins[0]), kStarBins[1], kStarBins[2]);
  }

  /// Set the PDG codes of the two particles involved
  /// \param pdg1 PDG code of particle one
  /// \param pdg2 PDG code of particle two
  void setPDGCodes(const int pdg1, const int pdg2)
  {
    kMassOne = TDatabasePDG::Instance()->GetParticle(pdg1)->Mass();
    kMassTwo = TDatabasePDG::Instance()->GetParticle(pdg2)->Mass();
    kPDGOne = pdg1;
    kPDGTwo = pdg2;
  }

  /// To compute the bin value for cavariance matrix
  /// \param qbin value of the qth k* bin
  /// \param ilmzero
  /// \param zeroimag
  /// \param ilmprim
  /// \param primimag
  int getBin(int qbin, int ilmzero, int zeroimag, int ilmprim, int primimag)
  {
    return qbin * kMaxJM * kMaxJM * 4 + (ilmprim * 2 + primimag) * kMaxJM * 2 + ilmzero * 2 + zeroimag;
  }

  /// Templated function to compute the necessary observables and fill the histograms for respective Spherical Harmonic
  /// \tparam T type of the femtouniverseparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param ChosenEventType same or mixed event
  /// \param maxl Maximum valie of L component of the spherical harmonics
  template <bool isMC, typename T>
  void addEventPair(T const& part1, T const& part2, uint8_t ChosenEventType, int /*maxl*/, bool isiden)
  {
    std::vector<std::complex<double>> fYlmBuffer(kMaxJM);
    std::vector<double> f3d;
    f3d = FemtoUniverseMath::newpairfunc(part1, kMassOne, part2, kMassTwo, isiden);

    const float kv = f3d[0];
    const float qout = f3d[1];
    const float qside = f3d[2];
    const float qlong = f3d[3];

    int nqbin = fbinctn->GetXaxis()->FindFixBin(kv) - 1;

    FemtoUniverseSpherHarMath kYlm;
    kYlm.doYlmUpToL(kMaxL, qout, qside, qlong, fYlmBuffer.data());

    if (ChosenEventType == femto_universe_sh_container::EventType::same) {
      for (int ihist = 0; ihist < kMaxJM; ihist++) {
        fnumsreal[ihist]->Fill(kv, real(fYlmBuffer[ihist]));
        fnumsimag[ihist]->Fill(kv, -imag(fYlmBuffer[ihist]));
        fbinctn->Fill(kv, 1.0);
      }

      if (nqbin < fbinctn->GetNbinsX()) {
        for (int ilmzero = 0; ilmzero < kMaxJM; ilmzero++) {
          for (int ilmprim = 0; ilmprim < kMaxJM; ilmprim++) {
            fcovmnum[getBin(nqbin, ilmzero, 0, ilmprim, 0)] += (real(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]));
            fcovmnum[getBin(nqbin, ilmzero, 0, ilmprim, 1)] += (real(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]));
            fcovmnum[getBin(nqbin, ilmzero, 1, ilmprim, 0)] += (-imag(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]));
            fcovmnum[getBin(nqbin, ilmzero, 1, ilmprim, 1)] += (-imag(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]));
          }
        }
      }
    } else if (ChosenEventType == femto_universe_sh_container::EventType::mixed) {
      for (int ihist = 0; ihist < kMaxJM; ihist++) {
        fdensreal[ihist]->Fill(kv, real(fYlmBuffer[ihist]));
        fdensimag[ihist]->Fill(kv, -imag(fYlmBuffer[ihist]));
      }
      if (nqbin < fbinctn->GetNbinsX()) {
        for (int ilmzero = 0; ilmzero < kMaxJM; ilmzero++) {
          for (int ilmprim = 0; ilmprim < kMaxJM; ilmprim++) {
            fcovmden[getBin(nqbin, ilmzero, 0, ilmprim, 0)] += (real(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]));
            fcovmden[getBin(nqbin, ilmzero, 0, ilmprim, 1)] += (real(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]));
            fcovmden[getBin(nqbin, ilmzero, 1, ilmprim, 0)] += (-imag(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]));
            fcovmden[getBin(nqbin, ilmzero, 1, ilmprim, 1)] += (-imag(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]));
          }
        }
      }
    }
  }

  /// Function to fill covariance matrix in 3D histograms
  /// \param ChosenEventType same or mixed event
  /// \param MaxJM Maximum value of J
  void packCov(uint8_t ChosenEventType, int /*MaxJM*/)
  {
    if (ChosenEventType == femto_universe_sh_container::EventType::same) {
      for (int ibin = 1; ibin <= fcovnum->GetNbinsX(); ibin++) {
        for (int ilmz = 0; ilmz < kMaxJM * 2; ilmz++) {
          for (int ilmp = 0; ilmp < kMaxJM * 2; ilmp++) {
            auto bin = getBin(ibin - 1, ilmz / 2, ilmz % 2, ilmp / 2, ilmp % 2);
            auto value = fcovmnum[bin];
            fcovnum->SetBinContent(ibin, ilmz + 1, ilmp + 1, value);
          }
        }
      }
    } else if (ChosenEventType == femto_universe_sh_container::EventType::mixed) {
      for (int ibin = 1; ibin <= fcovden->GetNbinsX(); ibin++) {
        for (int ilmz = 0; ilmz < kMaxJM * 2; ilmz++) {
          for (int ilmp = 0; ilmp < kMaxJM * 2; ilmp++) {
            auto bin = getBin(ibin - 1, ilmz / 2, ilmz % 2, ilmp / 2, ilmp % 2);
            auto value = fcovmden[bin];
            fcovden->SetBinContent(ibin, ilmz + 1, ilmp + 1, value);
          }
        }
      }
    }
  }

 private:
  std::array<std::shared_ptr<TH1>, 10> fnumsreal{};
  std::array<std::shared_ptr<TH1>, 10> fnumsimag{};
  std::array<std::shared_ptr<TH1>, 10> fdensreal{};
  std::array<std::shared_ptr<TH1>, 10> fdensimag{};

  std::shared_ptr<TH3> fcovnum{};
  std::shared_ptr<TH3> fcovden{};

  TH1D* fbinctn;
  TH1D* fbinctd;

  static constexpr int kMaxL = 1;
  static constexpr int kMaxJM = (kMaxL + 1) * (kMaxL + 1);

  std::array<float, (kMaxJM * kMaxJM * 4 * 100)> fcovmnum{}; ///< Covariance matrix for the numerator
  std::array<float, (kMaxJM * kMaxJM * 4 * 100)> fcovmden{}; ///< Covariance matrix for the numerator

 protected:
  HistogramRegistry* kHistogramRegistry = nullptr;                                  ///< For QA output
  static constexpr std::string_view kFolderSuffix[2] = {"SameEvent", "MixedEvent"}; ///< Folder naming for the output according to kEventType
  static constexpr int kEventType = eventType;                                      ///< Type of the event (same/mixed, according to FEMTOUNIVERSESHCONTAINER::EventType)
  float kMassOne = 0.f;                                                             ///< PDG mass of particle 1
  float kMassTwo = 0.f;                                                             ///< PDG mass of particle 2
  int kPDGOne = 0;                                                                  ///< PDG code of particle 1
  int kPDGTwo = 0;                                                                  ///< PDG code of particle 2
  std::vector<double> kStarBins;
};

} // namespace o2::analysis::femto_universe

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSESHCONTAINER_H_
