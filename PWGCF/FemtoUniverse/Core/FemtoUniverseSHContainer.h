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

/// \file FemtoUniverseSHContainer..h
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

using namespace o2::framework;

namespace o2::analysis::femtoUniverse
{

namespace femtoUniverseSHContainer
{
/// Femtoscopic observable to be computed
enum Observable { kstar ///< kstar
};

/// Type of the event processind
enum EventType { same, ///< Pair from same event
                 mixed ///< Pair from mixed event
};
}; // namespace femtoUniverseSHContainer

/// \class femtoUniverseSHContainer
/// \brief Container for all histogramming related to the correlation function. The two
/// particles of the pair are passed here, and the correlation function and QA histograms
/// are filled according to the specified observable
/// \tparam eventType Type of the event (same/mixed)
/// \tparam obs Observable to be computed (k*/Q_inv/...)
template <femtoUniverseSHContainer::EventType eventType, femtoUniverseSHContainer::Observable obs>
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
    KStarBins = kstarbins;
    std::string femtoObs1D;

    std::vector<int> fels(fMaxJM);
    std::vector<int> fems(fMaxJM);
    std::vector<int> felsi(fMaxJM);
    std::vector<int> femsi(fMaxJM);

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
    } while (el <= fMaxL);

    mHistogramRegistry = registry;
    femtoObs1D = "#it{q} (GeV/#it{c})";

    framework::AxisSpec femtoObsAxis1D = {kstarbins, femtoObs1D.c_str()};
    std::string folderName = static_cast<std::string>(mFolderSuffix[mEventType]) + static_cast<std::string>(o2::aod::femtouniverseMCparticle::MCTypeName[o2::aod::femtouniverseMCparticle::MCType::kRecon]);
    std::string suffix;

    for (int ihist = 0; ihist < fMaxJM; ihist++) {
      if (femsi[ihist] < 0) {
        suffix = "Ylm" + std::to_string(felsi[ihist]) + std::to_string(felsi[ihist] - femsi[ihist]);
      } else {
        suffix = "Ylm" + std::to_string(felsi[ihist]) + std::to_string(femsi[ihist]);
      }

      if (mFolderSuffix[mEventType] == mFolderSuffix[0]) {
        fnumsreal[ihist] = mHistogramRegistry->add<TH1>(("NumRe" + suffix).c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
        fnumsimag[ihist] = mHistogramRegistry->add<TH1>(("NumIm" + suffix).c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
      } else {
        fdensreal[ihist] = mHistogramRegistry->add<TH1>(("DenRe" + suffix).c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
        fdensimag[ihist] = mHistogramRegistry->add<TH1>(("DenIm" + suffix).c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
      }
    }

    if (mFolderSuffix[mEventType] == mFolderSuffix[0]) {
      std::string bufnameNum = "CovNum";
      fcovnum = mHistogramRegistry->add<TH3>((bufnameNum).c_str(), "; x; y; z", kTH3D, {{kstarbins}, {(fMaxJM * 2), -0.5, ((static_cast<float>(fMaxJM) * 2.0 - 0.5))}, {(fMaxJM * 2), -0.5, ((static_cast<float>(fMaxJM) * 2.0 - 0.5))}});
    } else if (mFolderSuffix[mEventType] == mFolderSuffix[1]) {
      std::string bufnameDen = "CovDen";
      fcovden = mHistogramRegistry->add<TH3>((bufnameDen).c_str(), "; x; y; z", kTH3D, {{kstarbins}, {(fMaxJM * 2), -0.5, ((static_cast<float>(fMaxJM) * 2.0 - 0.5))}, {(fMaxJM * 2), -0.5, ((static_cast<float>(fMaxJM) * 2.0 - 0.5))}});
    }

    fbinctn = new TH1D(TString("BinCountNum"), "Bin Occupation (Numerator)", static_cast<int>(KStarBins[0]), KStarBins[1], KStarBins[2]);
    fbinctd = new TH1D(TString("BinCountDen"), "Bin Occupation (Denominator)", static_cast<int>(KStarBins[0]), KStarBins[1], KStarBins[2]);
  }

  /// Set the PDG codes of the two particles involved
  /// \param pdg1 PDG code of particle one
  /// \param pdg2 PDG code of particle two
  void setPDGCodes(const int pdg1, const int pdg2)
  {
    mMassOne = TDatabasePDG::Instance()->GetParticle(pdg1)->Mass();
    mMassTwo = TDatabasePDG::Instance()->GetParticle(pdg2)->Mass();
    mPDGOne = pdg1;
    mPDGTwo = pdg2;
  }

  /// To compute the bin value for cavariance matrix
  /// \param qbin value of the qth k* bin
  /// \param ilmzero
  /// \param zeroimag
  /// \param ilmprim
  /// \param primimag
  int GetBin(int qbin, int ilmzero, int zeroimag, int ilmprim, int primimag)
  {
    return qbin * fMaxJM * fMaxJM * 4 + (ilmprim * 2 + primimag) * fMaxJM * 2 + ilmzero * 2 + zeroimag;
  }

  /// Templated function to compute the necessary observables and fill the histograms for respective Spherical Harmonic
  /// \tparam T type of the femtouniverseparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param ChosenEventType same or mixed event
  /// \param maxl Maximum valie of L component of the spherical harmonics
  template <bool isMC, typename T>
  void AddEventPair(T const& part1, T const& part2, uint8_t ChosenEventType, int /*maxl*/)
  {
    // int fMaxL = 2;
    // int fMaxJM = (2+1)*(2+1);
    std::vector<std::complex<double>> fYlmBuffer(fMaxJM);
    std::vector<double> f3d;
    f3d = FemtoUniverseMath::getpairmom3d(part1, mMassOne, part2, mMassTwo, true, true);

    // const float qstar = f3d[0];
    const float qout = f3d[1];
    const float qside = f3d[2];
    const float qlong = f3d[3];

    double kv = sqrt(qout * qout + qside * qside + qlong * qlong);
    int nqbin = fbinctn->GetXaxis()->FindFixBin(kv) - 1;

    FemtoUniverseSpherHarMath Ylm;
    Ylm.YlmUpToL(fMaxL, qout, qside, qlong, fYlmBuffer.data());

    if (ChosenEventType == femtoUniverseSHContainer::EventType::same) {
      for (int ihist = 0; ihist < fMaxJM; ihist++) {
        fnumsreal[ihist]->Fill(kv, real(fYlmBuffer[ihist]));
        fnumsimag[ihist]->Fill(kv, -imag(fYlmBuffer[ihist]));
        fbinctn->Fill(kv, 1.0);
      }

      if (nqbin < fbinctn->GetNbinsX()) {
        for (int ilmzero = 0; ilmzero < fMaxJM; ilmzero++) {
          for (int ilmprim = 0; ilmprim < fMaxJM; ilmprim++) {
            fcovmnum[GetBin(nqbin, ilmzero, 0, ilmprim, 0)] += (real(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]));
            fcovmnum[GetBin(nqbin, ilmzero, 0, ilmprim, 1)] += (real(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]));
            fcovmnum[GetBin(nqbin, ilmzero, 1, ilmprim, 0)] += (-imag(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]));
            fcovmnum[GetBin(nqbin, ilmzero, 1, ilmprim, 1)] += (-imag(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]));
          }
        }
      }
    } else if (ChosenEventType == femtoUniverseSHContainer::EventType::mixed) {
      for (int ihist = 0; ihist < fMaxJM; ihist++) {
        fdensreal[ihist]->Fill(kv, real(fYlmBuffer[ihist]));
        fdensimag[ihist]->Fill(kv, -imag(fYlmBuffer[ihist]));
      }
      if (nqbin < fbinctn->GetNbinsX()) {
        for (int ilmzero = 0; ilmzero < fMaxJM; ilmzero++) {
          for (int ilmprim = 0; ilmprim < fMaxJM; ilmprim++) {
            fcovmden[GetBin(nqbin, ilmzero, 0, ilmprim, 0)] += (real(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]));
            fcovmden[GetBin(nqbin, ilmzero, 0, ilmprim, 1)] += (real(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]));
            fcovmden[GetBin(nqbin, ilmzero, 1, ilmprim, 0)] += (-imag(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]));
            fcovmden[GetBin(nqbin, ilmzero, 1, ilmprim, 1)] += (-imag(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]));
          }
        }
      }
    }
  }

  /// Function to fill covariance matrix in 3D histograms
  /// \param ChosenEventType same or mixed event
  /// \param MaxJM Maximum value of J
  void PackCov(uint8_t ChosenEventType, int /*MaxJM*/)
  {
    if (ChosenEventType == femtoUniverseSHContainer::EventType::same) {
      for (int ibin = 1; ibin <= fcovnum->GetNbinsX(); ibin++) {
        for (int ilmz = 0; ilmz < fMaxJM * 2; ilmz++) {
          for (int ilmp = 0; ilmp < fMaxJM * 2; ilmp++) {
            auto bin = GetBin(ibin - 1, ilmz / 2, ilmz % 2, ilmp / 2, ilmp % 2);
            auto value = fcovmnum[bin];
            fcovnum->SetBinContent(ibin, ilmz + 1, ilmp + 1, value);
          }
        }
      }
    } else if (ChosenEventType == femtoUniverseSHContainer::EventType::mixed) {
      for (int ibin = 1; ibin <= fcovden->GetNbinsX(); ibin++) {
        for (int ilmz = 0; ilmz < fMaxJM * 2; ilmz++) {
          for (int ilmp = 0; ilmp < fMaxJM * 2; ilmp++) {
            auto bin = GetBin(ibin - 1, ilmz / 2, ilmz % 2, ilmp / 2, ilmp % 2);
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

  static constexpr int fMaxL = 1;
  static constexpr int fMaxJM = (fMaxL + 1) * (fMaxL + 1);

  std::array<float, (fMaxJM * fMaxJM * 4 * 100)> fcovmnum{}; ///< Covariance matrix for the numerator
  std::array<float, (fMaxJM * fMaxJM * 4 * 100)> fcovmden{}; ///< Covariance matrix for the numerator

 protected:
  HistogramRegistry* mHistogramRegistry = nullptr;                                  ///< For QA output
  static constexpr std::string_view mFolderSuffix[2] = {"SameEvent", "MixedEvent"}; ///< Folder naming for the output according to mEventType
  static constexpr int mEventType = eventType;                                      ///< Type of the event (same/mixed, according to FEMTOUNIVERSESHCONTAINER::EventType)
  float mMassOne = 0.f;                                                             ///< PDG mass of particle 1
  float mMassTwo = 0.f;                                                             ///< PDG mass of particle 2
  int mPDGOne = 0;                                                                  ///< PDG code of particle 1
  int mPDGTwo = 0;                                                                  ///< PDG code of particle 2
  std::vector<double> KStarBins;
};

} // namespace o2::analysis::femtoUniverse

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSESHCONTAINER_H_
