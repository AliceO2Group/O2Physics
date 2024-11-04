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

/// \file FemtoUniversePairSHCentMultKt.h
/// \brief FemtoUniversePairSHCentMultKt - Spherical Harmonics in mult, kT bins
/// \author Pritam Chakraborty, WUT, pritam.chakraborty@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRSHCENTMULTKT_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRSHCENTMULTKT_H_

#include <vector>
#include <string>
#include <complex>
#include <memory>
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;

namespace o2::analysis::femtoUniverse
{

/// \class FemtoUniversePairSHCentMultKt
/// \brief Container for all histogramming related to the spherical harmonics of
/// the correlation function. The two particles of the pair are passed here, and
/// the correlation function are filled according to the specified observable
/// \tparam eventType Type of the event (same/mixed)
/// \tparam obs Observable to be computed (k*/Q_inv/...)
template <femtoUniverseSHContainer::EventType eventType,
          femtoUniverseSHContainer::Observable obs>
class PairSHCentMultKt
{
 public:
  virtual ~PairSHCentMultKt() = default;
  /// @brief
  /// \tparam t1
  /// \param registry Histogram registry to be passed
  /// \param kstarbins Binning of k*
  /// \param centmultbins Number of multiplicity bins
  /// \param ktbins Number of kT bins
  /// \param maxl Maximum valie of L component of the spherical harmonics
  template <typename t1>
  void init(HistogramRegistry* registry, t1& kstarbins, t1& centmultbins,
            t1& ktbins, int /*maxl*/)
  {
    PairSHCentMultKtRegistry = registry;
    AxisSpec kstarAxis = {kstarbins, "#it{k*} (GeV/#it{c})"};
    KStarBins = kstarbins;

    CentMultBins = centmultbins;
    KtBins = ktbins;
    KtBins.erase(KtBins.begin());
    CentMultBins.erase(CentMultBins.begin());

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

    femtoObs1D = "#it{q} (GeV/#it{c})";

    framework::AxisSpec femtoObsAxis1D = {kstarbins, femtoObs1D.c_str()};

    for (int i = 0; i < static_cast<int>(CentMultBins.size() - 1); i++) {
      int lowBin = static_cast<int>((CentMultBins[i]));
      int highBin = static_cast<int>((CentMultBins[i + 1]));

      std::string HistTitle =
        "mult_" + std::to_string(lowBin) + "-" + std::to_string(highBin);
      std::string HistSuffix1 =
        std::to_string(static_cast<int>(CentMultBins[i]));
      std::string HistSuffix2 =
        std::to_string(static_cast<int>(CentMultBins[i + 1]));
      std::string HistFolderMult = "mult_" + HistSuffix1 + "_" + HistSuffix2;

      for (int j = 0; j < static_cast<int>(KtBins.size() - 1); j++) {
        int ktlowBin = static_cast<int>(KtBins[j]);
        int kthighBin = static_cast<int>(KtBins[j + 1]);

        std::string HistTitlekT =
          "kT_" + std::to_string(ktlowBin) + "-" + std::to_string(kthighBin);
        std::string HistSuffixkT1 =
          std::to_string(static_cast<int>(KtBins[j] * 100.0));
        std::string HistSuffixkT2 =
          std::to_string(static_cast<int>(KtBins[j + 1] * 100.0));
        std::string HistFolderkT = "kT_" + HistSuffixkT1 + "_" + HistSuffixkT2;

        std::string suffix;
        fbinctn[i][j] = new TH1D(
          TString("BinCountNum"), "Bin Occupation (Numerator)",
          static_cast<int>(KStarBins[0]), KStarBins[1], KStarBins[2]);
        fbinctd[i][j] = new TH1D(
          TString("BinCountDen"), "Bin Occupation (Denominator)",
          static_cast<int>(KStarBins[0]), KStarBins[1], KStarBins[2]);

        for (int ihist = 0; ihist < fMaxJM; ihist++) {
          if (femsi[ihist] < 0) {
            suffix = "Ylm" + std::to_string(felsi[ihist]) +
                     std::to_string(felsi[ihist] - femsi[ihist]);
          } else {
            suffix = "Ylm" + std::to_string(felsi[ihist]) +
                     std::to_string(femsi[ihist]);
          }
          if (mFolderSuffix[mEventType] == mFolderSuffix[0]) {
            fnumsreal[i][j][ihist] = PairSHCentMultKtRegistry->add<TH1>(
              (HistFolderMult + "/" + HistFolderkT + "/" + "NumRe" + suffix)
                .c_str(),
              ("; " + femtoObs1D + "; Entries").c_str(), kTH1D,
              {femtoObsAxis1D});
            fnumsimag[i][j][ihist] = PairSHCentMultKtRegistry->add<TH1>(
              (HistFolderMult + "/" + HistFolderkT + "/" + "NumIm" + suffix)
                .c_str(),
              ("; " + femtoObs1D + "; Entries").c_str(), kTH1D,
              {femtoObsAxis1D});
          } else {
            fdensreal[i][j][ihist] = PairSHCentMultKtRegistry->add<TH1>(
              (HistFolderMult + "/" + HistFolderkT + "/" + "DenRe" + suffix)
                .c_str(),
              ("; " + femtoObs1D + "; Entries").c_str(), kTH1D,
              {femtoObsAxis1D});
            fdensimag[i][j][ihist] = PairSHCentMultKtRegistry->add<TH1>(
              (HistFolderMult + "/" + HistFolderkT + "/" + "DenIm" + suffix)
                .c_str(),
              ("; " + femtoObs1D + "; Entries").c_str(), kTH1D,
              {femtoObsAxis1D});
          }
        }

        if (mFolderSuffix[mEventType] == mFolderSuffix[0]) {
          std::string bufnameNum = "CovNum";
          fcovnum[i][j] = PairSHCentMultKtRegistry->add<TH3>(
            (HistFolderMult + "/" + HistFolderkT + "/" + bufnameNum).c_str(),
            "; x; y; z", kTH3D,
            {{kstarbins},
             {(fMaxJM * 2), -0.5, ((static_cast<float>(fMaxJM) * 2.0 - 0.5))},
             {(fMaxJM * 2), -0.5,
              ((static_cast<float>(fMaxJM) * 2.0 - 0.5))}});
        } else if (mFolderSuffix[mEventType] == mFolderSuffix[1]) {
          std::string bufnameDen = "CovDen";
          fcovden[i][j] = PairSHCentMultKtRegistry->add<TH3>(
            (HistFolderMult + "/" + HistFolderkT + "/" + bufnameDen).c_str(),
            "; x; y; z", kTH3D,
            {{kstarbins},
             {(fMaxJM * 2), -0.5, ((static_cast<float>(fMaxJM) * 2.0 - 0.5))},
             {(fMaxJM * 2), -0.5,
              ((static_cast<float>(fMaxJM) * 2.0 - 0.5))}});
        }
      }
    }
  }

  /// Templated function to access different multiplicity directory and call
  /// fill_kT_NumDen \param part1 particle 1 \param part2 particle 2 \param
  /// ChosenEventType Same or Mixed evet type \param maxl Maximum valie of L
  /// component of the spherical harmonics \param multval Multiplicity value
  /// \param ktval kT value
  template <typename T>
  void fill_mult_NumDen(T const& part1, T const& part2, uint8_t ChosenEventType,
                        int maxl, int multval, float ktval, bool isiden)
  {
    int multbinval;
    int absmultval = multval;

    if ((absmultval >= CentMultBins[0]) && (absmultval < CentMultBins[1])) {
      multbinval = 0;
    } else if (absmultval < CentMultBins[2]) {
      multbinval = 1;
    } else if (absmultval < CentMultBins[3]) {
      multbinval = 2;
    } else if (absmultval < CentMultBins[4]) {
      multbinval = 3;
    } else {
      return;
    }
    // std::cout<<"multbinval "<<multbinval<<std::endl;
    fill_kT_NumDen(part1, part2, ChosenEventType, maxl, multbinval, ktval, isiden);
  }

  /// Templated function to access different kT directory and call AddEventPair
  /// \param part1 particle 1
  /// \param part2 particle 2
  /// \param ChosenEventType Same or Mixed evet type
  /// \param maxl Maximum valie of L component of the spherical harmonics
  /// \param multval Multiplicity value
  /// \param ktval kT value
  template <typename T>
  void fill_kT_NumDen(T const& part1, T const& part2, uint8_t ChosenEventType,
                      int maxl, int multval, float ktval, bool isiden)
  {
    int ktbinval = -1;
    if ((ktval >= KtBins[0]) && (ktval < KtBins[1])) {
      ktbinval = 0;
    } else if (ktval < KtBins[2]) {
      ktbinval = 1;
    } else if (ktval < KtBins[3]) {
      ktbinval = 2;
    } else if (ktval < KtBins[4]) {
      ktbinval = 3;
    } else {
      return;
    }
    AddEventPair(part1, part2, ChosenEventType, maxl, multval, ktbinval, isiden);
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
    return qbin * fMaxJM * fMaxJM * 4 + (ilmprim * 2 + primimag) * fMaxJM * 2 +
           ilmzero * 2 + zeroimag;
  }

  /// Templated function to compute the necessary observables and fill the
  /// histograms for respective Spherical Harmonic \tparam T type of the
  /// femtouniverseparticle \param part1 Particle one \param part2 Particle two
  /// \param ChosenEventType Same or Mixed evet type
  /// \param maxl Maximum valie of L component of the spherical harmonics
  /// \param multval Multiplicity value
  /// \param ktval kT value
  template <typename T>
  void AddEventPair(T const& part1, T const& part2, uint8_t ChosenEventType,
                    int /*maxl*/, int multval, int ktval, bool isiden)
  {
    int fMultBin = multval;
    int fKtBin = ktval;
    std::vector<std::complex<double>> fYlmBuffer(fMaxJM);
    std::vector<double> f3d;
    f3d = FemtoUniverseMath::newpairfunc(part1, mMassOne, part2, mMassTwo,
                                         isiden);

    const float qout = f3d[1];
    const float qside = f3d[2];
    const float qlong = f3d[3];

    double kv = sqrt(qout * qout + qside * qside + qlong * qlong);
    int nqbin = fbinctn[0][0]->GetXaxis()->FindFixBin(kv) - 1;

    FemtoUniverseSpherHarMath Ylm;
    Ylm.YlmUpToL(fMaxL, qout, qside, qlong, fYlmBuffer.data());

    if (ChosenEventType == femtoUniverseSHContainer::EventType::same) {
      for (int ihist = 0; ihist < fMaxJM; ihist++) {
        fnumsreal[fMultBin][fKtBin][ihist]->Fill(kv, real(fYlmBuffer[ihist]));
        fnumsimag[fMultBin][fKtBin][ihist]->Fill(kv, -imag(fYlmBuffer[ihist]));
        fbinctn[fMultBin][fKtBin]->Fill(kv, 1.0);
      }
      if (nqbin < fbinctn[0][0]->GetNbinsX()) {
        for (int ilmzero = 0; ilmzero < fMaxJM; ilmzero++) {
          for (int ilmprim = 0; ilmprim < fMaxJM; ilmprim++) {
            fcovmnum[fMultBin][fKtBin][GetBin(nqbin, ilmzero, 0, ilmprim, 0)] +=
              (real(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]));
            fcovmnum[fMultBin][fKtBin][GetBin(nqbin, ilmzero, 0, ilmprim, 1)] +=
              (real(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]));
            fcovmnum[fMultBin][fKtBin][GetBin(nqbin, ilmzero, 1, ilmprim, 0)] +=
              (-imag(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]));
            fcovmnum[fMultBin][fKtBin][GetBin(nqbin, ilmzero, 1, ilmprim, 1)] +=
              (-imag(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]));
          }
        }
      }
    } else if (ChosenEventType == femtoUniverseSHContainer::EventType::mixed) {
      for (int ihist = 0; ihist < fMaxJM; ihist++) {
        fdensreal[fMultBin][fKtBin][ihist]->Fill(kv, real(fYlmBuffer[ihist]));
        fdensimag[fMultBin][fKtBin][ihist]->Fill(kv, -imag(fYlmBuffer[ihist]));
        fbinctd[fMultBin][fKtBin]->Fill(kv, 1.0);
      }
      if (nqbin < fbinctd[0][0]->GetNbinsX()) {
        for (int ilmzero = 0; ilmzero < fMaxJM; ilmzero++) {
          for (int ilmprim = 0; ilmprim < fMaxJM; ilmprim++) {
            fcovmden[fMultBin][fKtBin][GetBin(nqbin, ilmzero, 0, ilmprim, 0)] +=
              (real(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]));
            fcovmden[fMultBin][fKtBin][GetBin(nqbin, ilmzero, 0, ilmprim, 1)] +=
              (real(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]));
            fcovmden[fMultBin][fKtBin][GetBin(nqbin, ilmzero, 1, ilmprim, 0)] +=
              (-imag(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]));
            fcovmden[fMultBin][fKtBin][GetBin(nqbin, ilmzero, 1, ilmprim, 1)] +=
              (-imag(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]));
          }
        }
      }
    }
  }

  /// Function to fill covariance matrix in 3D histograms
  /// \param ChosenEventType same or mixed event
  /// \param MaxJM Maximum value of J
  /// \param multval Multiplicity value
  /// \param ktval kT value
  void PackCov(uint8_t ChosenEventType, int /*MaxJM*/, int multval, int ktval)
  {
    int fMultBin = multval;
    int fKtBin = ktval;
    if (ChosenEventType == femtoUniverseSHContainer::EventType::same) {
      for (int ibin = 1; ibin <= fcovnum[0][0]->GetNbinsX(); ibin++) {
        for (int ilmz = 0; ilmz < fMaxJM * 2; ilmz++) {
          for (int ilmp = 0; ilmp < fMaxJM * 2; ilmp++) {
            auto bin = GetBin(ibin - 1, ilmz / 2, ilmz % 2, ilmp / 2, ilmp % 2);
            auto value = fcovmnum[fMultBin][fKtBin][bin];
            fcovnum[fMultBin][fKtBin]->SetBinContent(ibin, ilmz + 1, ilmp + 1,
                                                     value);
          }
        }
      }
    } else if (ChosenEventType == femtoUniverseSHContainer::EventType::mixed) {
      for (int ibin = 1; ibin <= fcovden[0][0]->GetNbinsX(); ibin++) {
        for (int ilmz = 0; ilmz < fMaxJM * 2; ilmz++) {
          for (int ilmp = 0; ilmp < fMaxJM * 2; ilmp++) {
            auto bin = GetBin(ibin - 1, ilmz / 2, ilmz % 2, ilmp / 2, ilmp % 2);
            auto value = fcovmden[fMultBin][fKtBin][bin];
            fcovden[fMultBin][fKtBin]->SetBinContent(ibin, ilmz + 1, ilmp + 1,
                                                     value);
          }
        }
      }
    }
  }

  /// Function to acces each multiplicity and kT directory and call PackCov
  /// \param ChosenEventType same or mixed event
  /// \param MaxJM Maximum value of J
  void fill_mult_kT_Cov(uint8_t ChosenEventType, int MaxJM)
  {
    for (int multbinvalcov = 0;
         multbinvalcov < static_cast<int>(CentMultBins.size() - 2);
         multbinvalcov++) {
      for (int ktbinvalcov = 0;
           ktbinvalcov < static_cast<int>(KtBins.size() - 1); ktbinvalcov++) {
        PackCov(ChosenEventType, MaxJM, multbinvalcov, ktbinvalcov);
      }
    }
  }

 private:
  std::array<std::array<std::array<std::shared_ptr<TH1>, 10>, 4>, 4>
    fnumsreal{};
  std::array<std::array<std::array<std::shared_ptr<TH1>, 10>, 4>, 4>
    fnumsimag{};
  std::array<std::array<std::array<std::shared_ptr<TH1>, 10>, 4>, 4>
    fdensreal{};
  std::array<std::array<std::array<std::shared_ptr<TH1>, 10>, 4>, 4>
    fdensimag{};

  TH1D* fbinctn[10][10];
  TH1D* fbinctd[10][10];

  static constexpr int fMaxL = 2;
  static constexpr int fMaxJM = (fMaxL + 1) * (fMaxL + 1);

  std::array<std::array<std::array<float, (fMaxJM * fMaxJM * 4 * 100)>, 4>, 4>
    fcovmnum{}; ///< Covariance matrix for the numerator
  std::array<std::array<std::array<float, (fMaxJM * fMaxJM * 4 * 100)>, 4>, 4>
    fcovmden{}; ///< Covariance matrix for the numerator

  std::array<std::array<std::shared_ptr<TH3>, 4>, 4> fcovnum{};
  std::array<std::array<std::shared_ptr<TH3>, 4>, 4> fcovden{};

 protected:
  HistogramRegistry* PairSHCentMultKtRegistry = nullptr;
  static constexpr std::string_view mFolderSuffix[2] = {
    "SameEvent",
    "MixedEvent"}; ///< Folder naming for the output according to mEventType
  static constexpr int mEventType =
    eventType;          ///< Type of the event (same/mixed, according to
                        ///< FEMTOUNIVERSESHCONTAINER::EventType)
  float mMassOne = 0.f; ///< PDG mass of particle 1
  float mMassTwo = 0.f; ///< PDG mass of particle 2
  int mPDGOne = 0;      ///< PDG code of particle 1
  int mPDGTwo = 0;      ///< PDG code of particle 2
  std::vector<double> CentMultBins;
  std::vector<double> KtBins;
  std::vector<double> KStarBins;
  bool UseKt = false;
  bool Use3D = false;
};

} // namespace o2::analysis::femtoUniverse

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRSHCENTMULTKT_H_
