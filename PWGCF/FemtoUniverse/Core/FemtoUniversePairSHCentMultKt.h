// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseSHContainer.h"

#include "Framework/HistogramRegistry.h"

#include <complex>
#include <memory>
#include <string>
#include <vector>

// using namespace o2::constants::physics;

namespace o2::analysis::femto_universe
{

/// \class FemtoUniversePairSHCentMultKt
/// \brief Container for all histogramming related to the spherical harmonics of
/// the correlation function. The two particles of the pair are passed here, and
/// the correlation function are filled according to the specified observable
/// \tparam eventType Type of the event (same/mixed)
/// \tparam obs Observable to be computed (k*/Q_inv/...)
template <femto_universe_sh_container::EventType eventType,
          femto_universe_sh_container::Observable obs>
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
            t1& ktbins, bool isqinvfill, int /*maxl*/)
  {
    pairSHCentMultKtRegistry = registry;
    AxisSpec kstarAxis = {kstarbins, "#it{k*} (GeV/#it{c})"};
    kStarBins = kstarbins;

    centMultBins = centmultbins;
    ktBins = ktbins;
    ktBins.erase(ktBins.begin());
    centMultBins.erase(centMultBins.begin());

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

    femtoObs1D = "#it{q} (GeV/#it{c})";

    framework::AxisSpec femtoObsAxis1D = {kstarbins, femtoObs1D.c_str()};

    for (int i = 0; i < static_cast<int>(centMultBins.size() - 1); i++) {
      int lowBin = static_cast<int>((centMultBins[i]));
      int highBin = static_cast<int>((centMultBins[i + 1]));

      std::string histTitle = "mult_" + std::to_string(lowBin) + "-" + std::to_string(highBin);
      std::string histSuffix1 = std::to_string(static_cast<int>(centMultBins[i]));
      std::string histSuffix2 = std::to_string(static_cast<int>(centMultBins[i + 1]));
      std::string histFolderMult = "mult_" + histSuffix1 + "_" + histSuffix2;

      for (int j = 0; j < static_cast<int>(ktBins.size() - 1); j++) {
        int ktlowBin = static_cast<int>(ktBins[j]);
        int kthighBin = static_cast<int>(ktBins[j + 1]);

        std::string histTitlekT = "kT_" + std::to_string(ktlowBin) + "-" + std::to_string(kthighBin);
        std::string histSuffixkT1 = std::to_string(static_cast<int>(ktBins[j] * 100.0));
        std::string histSuffixkT2 = std::to_string(static_cast<int>(ktBins[j + 1] * 100.0));
        std::string histFolderkT = "kT_" + histSuffixkT1 + "_" + histSuffixkT2;

        std::string suffix;
        fbinctn[i][j] = new TH1D(TString("BinCountNum"), "Bin Occupation (Numerator)", static_cast<int>(kStarBins[0]), kStarBins[1], kStarBins[2]);
        fbinctd[i][j] = new TH1D(TString("BinCountDen"), "Bin Occupation (Denominator)", static_cast<int>(kStarBins[0]), kStarBins[1], kStarBins[2]);

        for (int ihist = 0; ihist < kMaxJM; ihist++) {
          if (femsi[ihist] < 0) {
            suffix = "Ylm" + std::to_string(felsi[ihist]) +
                     std::to_string(felsi[ihist] - femsi[ihist]);
          } else {
            suffix = "Ylm" + std::to_string(felsi[ihist]) + std::to_string(femsi[ihist]);
          }
          // std::cout<<"ihist "<<ihist<<" "<<suffix.c_str()<<std::endl;
          if (FolderSuffix[EventType] == FolderSuffix[0]) {
            fnumsreal[i][j][ihist] = pairSHCentMultKtRegistry->add<TH1>(
              (histFolderMult + "/" + histFolderkT + "/" + "NumRe" + suffix).c_str(),
              ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
            fnumsimag[i][j][ihist] = pairSHCentMultKtRegistry->add<TH1>(
              (histFolderMult + "/" + histFolderkT + "/" + "NumIm" + suffix).c_str(),
              ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
          } else {
            fdensreal[i][j][ihist] = pairSHCentMultKtRegistry->add<TH1>(
              (histFolderMult + "/" + histFolderkT + "/" + "DenRe" + suffix).c_str(),
              ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
            fdensimag[i][j][ihist] = pairSHCentMultKtRegistry->add<TH1>(
              (histFolderMult + "/" + histFolderkT + "/" + "DenIm" + suffix).c_str(),
              ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
          }
        }

        if (FolderSuffix[EventType] == FolderSuffix[0]) {
          std::string bufnameNum = "CovNum";
          fcovnum[i][j] = pairSHCentMultKtRegistry->add<TH3>((histFolderMult + "/" + histFolderkT + "/" + bufnameNum).c_str(), "; x; y; z", kTH3D,
                                                             {{kstarbins},
                                                              {(kMaxJM * 2), -0.5, ((static_cast<float>(kMaxJM) * 2.0 - 0.5))},
                                                              {(kMaxJM * 2), -0.5,
                                                               ((static_cast<float>(kMaxJM) * 2.0 - 0.5))}});
          fcovnum[i][j]->Sumw2();
        } else if (FolderSuffix[EventType] == FolderSuffix[1]) {
          std::string bufnameDen = "CovDen";
          fcovden[i][j] = pairSHCentMultKtRegistry->add<TH3>((histFolderMult + "/" + histFolderkT + "/" + bufnameDen).c_str(), "; x; y; z", kTH3D,
                                                             {{kstarbins},
                                                              {(kMaxJM * 2), -0.5, ((static_cast<float>(kMaxJM) * 2.0 - 0.5))},
                                                              {(kMaxJM * 2), -0.5,
                                                               ((static_cast<float>(kMaxJM) * 2.0 - 0.5))}});
        }
        if (isqinvfill) {
          if (FolderSuffix[EventType] == FolderSuffix[0]) {
            std::string bufnameNum = "h1DNum";
            fnums1D[i][j] = pairSHCentMultKtRegistry->add<TH1>(
              (histFolderMult + "/" + histFolderkT + "/" + bufnameNum).c_str(),
              ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
            fnums1D[i][j]->Sumw2();
          } else if (FolderSuffix[EventType] == FolderSuffix[1]) {
            std::string bufnameNum = "h1DDen";
            fdens1D[i][j] = pairSHCentMultKtRegistry->add<TH1>(
              (histFolderMult + "/" + histFolderkT + "/" + bufnameNum).c_str(),
              ("; " + femtoObs1D + "; Entries").c_str(), kTH1D, {femtoObsAxis1D});
            fdens1D[i][j]->Sumw2();
          }
        }
      }
    }
  }

  /// Templated function to access different multiplicity directory and call
  /// fillkTNumDen \param part1 particle 1 \param part2 particle 2 \param
  /// ChosenEventType Same or Mixed evet type \param maxl Maximum valie of L
  /// component of the spherical harmonics \param multval Multiplicity value
  /// \param ktval kT value
  template <typename T>
  void fillMultNumDen(T const& part1, T const& part2, uint8_t ChosenEventType,
                      int maxl, int multval, float ktval, bool isiden, bool isqinvfill)
  {
    int multbinval;
    int absmultval = multval;

    if ((absmultval >= centMultBins[0]) && (absmultval < centMultBins[1])) {
      multbinval = 0;
    } else if ((absmultval >= centMultBins[1]) && (absmultval < centMultBins[2])) {
      multbinval = 1;
    } else if ((absmultval >= centMultBins[2]) && (absmultval < centMultBins[3])) {
      multbinval = 2;
    } else if ((absmultval >= centMultBins[3]) && (absmultval < centMultBins[4])) {
      multbinval = 3;
    } else {
      return;
    }
    // std::cout<<"multbinval "<<multbinval<<std::endl;
    fillkTNumDen(part1, part2, ChosenEventType, maxl, multbinval, ktval, isiden, isqinvfill);
  }

  /// Templated function to access different kT directory and call addEventPair
  /// \param part1 particle 1
  /// \param part2 particle 2
  /// \param ChosenEventType Same or Mixed evet type
  /// \param maxl Maximum valie of L component of the spherical harmonics
  /// \param multval Multiplicity value
  /// \param ktval kT value
  template <typename T>
  void fillkTNumDen(T const& part1, T const& part2, uint8_t ChosenEventType,
                    int maxl, int multval, float ktval, bool isiden, bool isqinvfill)
  {
    int ktbinval = -1;
    if (ktval >= ktBins[0] && ktval < ktBins[1]) {
      ktbinval = 0;
    } else if (ktval >= ktBins[1] && ktval < ktBins[2]) {
      ktbinval = 1;
    } else if (ktval >= ktBins[2] && ktval < ktBins[3]) {
      ktbinval = 2;
    } else if (ktval >= ktBins[3] && ktval < ktBins[4]) {
      ktbinval = 3;
    } else if (ktval >= ktBins[4] && ktval < ktBins[5]) {
      ktbinval = 4;
    } else if (ktval >= ktBins[5] && ktval < ktBins[6]) {
      ktbinval = 5;
    } else if (ktval >= ktBins[6] && ktval < ktBins[7]) {
      ktbinval = 6;
    } else {
      return;
    }
    addEventPair(part1, part2, ChosenEventType, maxl, multval, ktbinval, isiden, isqinvfill);
  }

  /// Set the PDG codes of the two particles involved
  /// \param pdg1 PDG code of particle one
  /// \param pdg2 PDG code of particle two
  void setPionPairMass()
  {
    mMassOne = o2::constants::physics::MassPiPlus; // FIXME: Get from the PDG service of the common header
    mMassTwo = o2::constants::physics::MassPiPlus; // FIXME: Get from the PDG service of the common header
  }

  /// To compute the bin value for cavariance matrix
  /// \param qbin value of the qth k* bin
  /// \param ilmzero
  /// \param zeroimag
  /// \param ilmprim
  /// \param primimag
  int getBin(int qbin, int ilmzero, int zeroimag, int ilmprim, int primimag)
  {
    return qbin * kMaxJM * kMaxJM * 4 + (ilmprim * 2 + primimag) * kMaxJM * 2 +
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
  void addEventPair(T const& part1, T const& part2, uint8_t ChosenEventType,
                    int /*maxl*/, int multval, int ktval, bool isiden, bool isqinvfill)
  {
    int fMultBin = multval;
    int fKtBin = ktval;
    std::vector<std::complex<double>> fYlmBuffer(kMaxJM);
    std::vector<double> f3d;
    setPionPairMass();
    f3d = FemtoUniverseMath::newpairfunc(part1, mMassOne, part2, mMassTwo,
                                         isiden);

    const float qout = f3d[1];
    const float qside = f3d[2];
    const float qlong = f3d[3];

    double kv = std::sqrt(qout * qout + qside * qside + qlong * qlong);

    // int nqbin = fbinctn[0][0]->GetXaxis()->FindFixBin(kv);
    // int nqbinnotfix = fbinctn[0][0]->GetXaxis()->FindBin(kv);

    FemtoUniverseSpherHarMath kYlm;
    kYlm.doYlmUpToL(kMaxL, qout, qside, qlong, fYlmBuffer.data());

    if (ChosenEventType == femto_universe_sh_container::EventType::same) {
      for (int ihist = 0; ihist < kMaxJM; ihist++) {
        fnumsreal[fMultBin][fKtBin][ihist]->Fill(kv, real(fYlmBuffer[ihist]));
        fnumsimag[fMultBin][fKtBin][ihist]->Fill(kv, -imag(fYlmBuffer[ihist]));
        fbinctn[fMultBin][fKtBin]->Fill(kv, 1.0);
      }
      if (isqinvfill) {
        fnums1D[fMultBin][fKtBin]->Fill(f3d[0]);
      }
      for (int ilmzero = 0; ilmzero < kMaxJM * 2; ilmzero++) {
        for (int ilmprim = 0; ilmprim < kMaxJM * 2; ilmprim++) {
          if ((ilmzero % 2) == 0 && (ilmprim % 2) == 0) {
            fcovnum[fMultBin][fKtBin]->Fill(kv, ilmzero, ilmprim, (real(fYlmBuffer[ilmzero / 2]) * real(fYlmBuffer[ilmprim / 2])));
          } else if ((ilmzero % 2) == 0 && (ilmprim % 2) == 1) {
            fcovnum[fMultBin][fKtBin]->Fill(kv, ilmzero, ilmprim, (real(fYlmBuffer[ilmzero / 2]) * -imag(fYlmBuffer[ilmprim / 2])));
          } else if ((ilmzero % 2) == 1 && (ilmprim % 2) == 0) {
            fcovnum[fMultBin][fKtBin]->Fill(kv, ilmzero, ilmprim, (-imag(fYlmBuffer[ilmzero / 2]) * real(fYlmBuffer[ilmprim / 2])));
          } else if ((ilmzero % 2) == 1 && (ilmprim % 2) == 1) {
            fcovnum[fMultBin][fKtBin]->Fill(kv, ilmzero, ilmprim, (-imag(fYlmBuffer[ilmzero / 2]) * -imag(fYlmBuffer[ilmprim / 2])));
          }
        }
      }
    } else if (ChosenEventType == femto_universe_sh_container::EventType::mixed) {
      for (int ihist = 0; ihist < kMaxJM; ihist++) {
        fdensreal[fMultBin][fKtBin][ihist]->Fill(kv, real(fYlmBuffer[ihist]));
        fdensimag[fMultBin][fKtBin][ihist]->Fill(kv, -imag(fYlmBuffer[ihist]));
        fbinctd[fMultBin][fKtBin]->Fill(kv, 1.0);
      }
      if (isqinvfill) {
        fdens1D[fMultBin][fKtBin]->Fill(f3d[0]);
      }
      for (int ilmzero = 0; ilmzero < kMaxJM * 2; ilmzero++) {
        for (int ilmprim = 0; ilmprim < kMaxJM * 2; ilmprim++) {
          if ((ilmzero % 2) == 0 && (ilmprim % 2) == 0) {
            fcovden[fMultBin][fKtBin]->Fill(kv, ilmzero, ilmprim, (real(fYlmBuffer[ilmzero / 2]) * real(fYlmBuffer[ilmprim / 2])));
          } else if ((ilmzero % 2) == 0 && (ilmprim % 2) == 1) {
            fcovden[fMultBin][fKtBin]->Fill(kv, ilmzero, ilmprim, (real(fYlmBuffer[ilmzero / 2]) * -imag(fYlmBuffer[ilmprim / 2])));
          } else if ((ilmzero % 2) == 1 && (ilmprim % 2) == 0) {
            fcovden[fMultBin][fKtBin]->Fill(kv, ilmzero, ilmprim, (-imag(fYlmBuffer[ilmzero / 2]) * real(fYlmBuffer[ilmprim / 2])));
          } else if ((ilmzero % 2) == 1 && (ilmprim % 2) == 1) {
            fcovden[fMultBin][fKtBin]->Fill(kv, ilmzero, ilmprim, (-imag(fYlmBuffer[ilmzero / 2]) * -imag(fYlmBuffer[ilmprim / 2])));
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
  void packCov(uint8_t ChosenEventType, int /*MaxJM*/, int multval, int ktval)
  {
    int fMultBin = multval;
    int fKtBin = ktval;
    if (ChosenEventType == femto_universe_sh_container::EventType::same) {
      for (int ibin = 1; ibin <= fcovnum[0][0]->GetNbinsX(); ibin++) {
        for (int ilmz = 0; ilmz < kMaxJM * 2; ilmz++) {
          for (int ilmp = 0; ilmp < kMaxJM * 2; ilmp++) {
            auto bin = getBin(ibin - 1, ilmz / 2, ilmz % 2, ilmp / 2, ilmp % 2);
            auto value = fcovmnum[fMultBin][fKtBin][bin];
            fcovnum[fMultBin][fKtBin]->SetBinContent(ibin, ilmz + 1, ilmp + 1, value);
          }
        }
      }
    } else if (ChosenEventType == femto_universe_sh_container::EventType::mixed) {
      for (int ibin = 1; ibin <= fcovden[0][0]->GetNbinsX(); ibin++) {
        for (int ilmz = 0; ilmz < kMaxJM * 2; ilmz++) {
          for (int ilmp = 0; ilmp < kMaxJM * 2; ilmp++) {
            auto bin = getBin(ibin - 1, ilmz / 2, ilmz % 2, ilmp / 2, ilmp % 2);
            auto value = fcovmden[fMultBin][fKtBin][bin];
            fcovden[fMultBin][fKtBin]->SetBinContent(ibin, ilmz + 1, ilmp + 1, value);
          }
        }
      }
    }
  }

  /// Function to acces each multiplicity and kT directory and call packCov
  /// \param ChosenEventType same or mixed event
  /// \param MaxJM Maximum value of J
  void fillMultkTCov(uint8_t ChosenEventType, int MaxJM)
  {
    for (int multbinvalcov = 0;
         multbinvalcov < static_cast<int>(centMultBins.size() - 1);
         multbinvalcov++) {
      for (int ktbinvalcov = 0;
           ktbinvalcov < static_cast<int>(ktBins.size() - 1); ktbinvalcov++) {
        packCov(ChosenEventType, MaxJM, multbinvalcov, ktbinvalcov);
      }
    }
  }

 private:
  std::array<std::array<std::array<std::shared_ptr<TH1>, 10>, 7>, 4> fnumsreal{};
  std::array<std::array<std::array<std::shared_ptr<TH1>, 10>, 7>, 4> fnumsimag{};
  std::array<std::array<std::array<std::shared_ptr<TH1>, 10>, 7>, 4> fdensreal{};
  std::array<std::array<std::array<std::shared_ptr<TH1>, 10>, 7>, 4> fdensimag{};
  std::array<std::array<std::shared_ptr<TH1>, 7>, 4> fnums1D{};
  std::array<std::array<std::shared_ptr<TH1>, 7>, 4> fdens1D{};

  TH1D* fbinctn[10][10];
  TH1D* fbinctd[10][10];

  static constexpr int kMaxL = 2;
  static constexpr int kMaxJM = (kMaxL + 1) * (kMaxL + 1);

  std::array<std::array<std::array<float, (kMaxJM * kMaxJM * 4 * 100)>, 7>, 4> fcovmnum{}; ///< Covariance matrix for the numerator
  std::array<std::array<std::array<float, (kMaxJM * kMaxJM * 4 * 100)>, 7>, 4> fcovmden{}; ///< Covariance matrix for the numerator

  std::array<std::array<std::shared_ptr<TH3>, 7>, 4> fcovnum{};
  std::array<std::array<std::shared_ptr<TH3>, 7>, 4> fcovden{};

 protected:
  HistogramRegistry* pairSHCentMultKtRegistry = nullptr;
  static constexpr std::string_view FolderSuffix[2] = {"SameEvent", "MixedEvent"}; ///< Folder naming for the output according to EventType
  static constexpr int EventType = eventType;                                      ///< Type of the event (same/mixed, according to FEMTOUNIVERSESHCONTAINER::EventType)
  float mMassOne = 0.f;                                                            ///< PDG mass of particle 1
  float mMassTwo = 0.f;                                                            ///< PDG mass of particle 2
  int mPDGOne = 0;                                                                 ///< PDG code of particle 1
  int mPDGTwo = 0;                                                                 ///< PDG code of particle 2
  std::vector<double> centMultBins;
  std::vector<double> ktBins;
  std::vector<double> kStarBins;
  bool useKt = false;
  bool use3D = false;
};

} // namespace o2::analysis::femto_universe

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRSHCENTMULTKT_H_
