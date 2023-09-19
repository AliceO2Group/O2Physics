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
// O2 includes

/// \file HFFilterHelpers.cxx
/// \brief Header file with definition of variables, methods, and tables used in the HFFilter.cxx task
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Marcel Lesch <marcel.lesch@tum.de>, TUM
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
/// \author Biao Zhang <biao.zhang@cern.ch>, CCNU
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University

#include "HFFilterHelpers.cxx"

namespace o2::aod
{
namespace hffilters
{

/// Update the TPC PID baesd on the spline of particles
/// \param track is a track parameter
/// \param mMassPar is the mass of particles
/// \param paraBetheBloch  vector for the parameters of BetheBloch formula
/// \return updated nsigma value for TPC PID
template <typename T>
double HfFilterHelper::getTPCSplineCalib(const T& track, const float mMassPar, const std::vector<double> paraBetheBloch)
{
  auto bgScaling = 1 / mMassPar;
  double expBethe = tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() * bgScaling), paraBetheBloch[0], paraBetheBloch[1], paraBetheBloch[2], paraBetheBloch[3], paraBetheBloch[4]);
  double expSigma = expBethe * paraBetheBloch[5];
  return static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
}

/// Single-track cuts for bachelor track of beauty candidates
/// \param track is a track parameter
/// \param trackPar is a track parameter
/// \param dca is the 2d array with dcaXY and dcaZ of the track
/// \param pTMinSoftPion min pT for soft pions
/// \param pTMinBeautyBachelor min pT for beauty bachelor pions
/// \param pTBinsTrack pT bins for dca cuts
/// \param cutsSingleTrackBeauty cuts for all tracks
/// \return a flag that encodes the selection for soft pions BIT(kSoftPion), tracks for beauty BIT(kForBeauty), or soft pions for beauty BIT(kSoftPionForBeauty)
template <typename T, typename T1, typename T2, typename T3, typename T4>
int8_t HfFilterHelper::isSelectedTrackForSoftPionOrBeauty(const T track, const T1& trackPar, const T2& dca, const float& pTMinSoftPion, const float& pTMinBeautyBachelor, const T3& pTBinsTrack, const T4& cutsSingleTrackBeauty)
{

  int8_t retValue{BIT(kSoftPion) | BIT(kForBeauty) | BIT(kSoftPionForBeauty)};

  if (!track.isGlobalTrackWoDCA()) {
    return kRejected;
  }

  auto pT = trackPar.getPt();
  auto pTBinTrack = findBin(pTBinsTrack, pT);
  if (pTBinTrack == -1) {
    return kRejected;
  }

  if (pT < pTMinSoftPion) { // soft pion should be less stringent than usual tracks
    return kRejected;
  }

  if (std::fabs(trackPar.getEta()) > 0.8) {
    return kRejected;
  }

  if (std::fabs(dca[1]) > 2.f) {
    return kRejected;
  }

  // below only regular beauty tracks, not required for soft pions
  if (pT < pTMinBeautyBachelor) {
    CLRBIT(retValue, kForBeauty);
  }

  if (std::fabs(dca[0]) < cutsSingleTrackBeauty.get(pTBinTrack, 0u)) { // minimum DCAxy
    CLRBIT(retValue, kForBeauty);
    CLRBIT(retValue, kSoftPionForBeauty);
  }
  if (std::fabs(dca[0]) > cutsSingleTrackBeauty.get(pTBinTrack, 1u)) { // maximum DCAxy
    CLRBIT(retValue, kForBeauty);
    CLRBIT(retValue, kSoftPionForBeauty);
  }

  return retValue;
}

/// Basic selection of proton candidates
/// \param track is a track
/// \param trackPar is a track parameter
/// \param femtoMinProtonPt min pT for proton candidates
/// \param femtoMaxNsigmaProton max Nsigma for proton candidates
/// \param setTPCCalib flag to activate TPC PID postcalibrations
/// \param hMapProton map of nSigma mean and sigma calibrations for proton
/// \param hSplineProton spline of proton and anti-proton calibrations
/// \param activateQA flag to activate the filling of QA histos
/// \param hProtonTPCPID histo with NsigmaTPC vs. p
/// \param hProtonTOFPID histo with NsigmaTOF vs. p
/// \return true if track passes all cuts
template <typename T1, typename T2, typename H2, typename H3>
bool HfFilterHelper::isSelectedProton4Femto(const T1& track, const T2& trackPar, const float& femtoMinProtonPt, const float& femtoMaxNsigmaProton, const int setTPCCalib, H3 hMapProton, const std::array<std::vector<double>, 2>& hSplineProton, const int& activateQA, H2 hProtonTPCPID, H2 hProtonTOFPID)
{
  if (trackPar.getPt() < femtoMinProtonPt) {
    return false;
  }

  if (std::fabs(trackPar.getEta()) > 0.8) {
    return false;
  }

  if (!track.isGlobalTrack()) {
    return false; // use only global tracks
  }

  float NSigmaTPC = track.tpcNSigmaPr();
  float NSigmaTOF = track.tofNSigmaPr();

  if (setTPCCalib == 1) {
    NSigmaTPC = getTPCPostCalib(hMapProton, track, kPr);
  } else if (setTPCCalib == 2) {
    if (track.sign() > 0) {
      NSigmaTPC = getTPCSplineCalib(track, massProton, hSplineProton[0]);
    } else {
      NSigmaTPC = getTPCSplineCalib(track, massProton, hSplineProton[1]);
    }
  }

  float NSigma = std::sqrt(NSigmaTPC * NSigmaTPC + NSigmaTOF * NSigmaTOF);

  if (NSigma > femtoMaxNsigmaProton) {
    return false;
  }

  if (activateQA > 1) {
    hProtonTPCPID->Fill(track.p(), NSigmaTPC);
    hProtonTOFPID->Fill(track.p(), NSigmaTOF);
  }

  return true;
}

/// Basic selection of proton candidates for Lc
/// \param track is a track
/// \param nsigmaTPCProtonLc max NsigmaTPC for proton candidates
/// \param nsigmaTOFProtonLc max NsigmaTOF for proton candidates
/// \param setTPCCalib flag to activate TPC PID postcalibrations
/// \param hMapProton map of nSigma mean and sigma calibrations for proton
/// \param hSplineProton spline of proton and anti-proton calibrations
/// \return true if track passes all cuts
template <typename T, typename H3>
bool HfFilterHelper::isSelectedProton4CharmBaryons(const T& track, const float& nsigmaTPCProtonLc, const float& nsigmaTOFProtonLc, const int setTPCCalib, H3 hMapProton, const std::array<std::vector<double>, 2>& hSplineProton)
{
  float NSigmaTPC = track.tpcNSigmaPr();
  float NSigmaTOF = track.tofNSigmaPr();

  if (setTPCCalib == 1) {
    NSigmaTPC = getTPCPostCalib(hMapProton, track, kPr);
  } else if (setTPCCalib == 2) {
    if (track.sign() > 0) {
      NSigmaTPC = getTPCSplineCalib(track, massProton, hSplineProton[0]);
    } else {
      NSigmaTPC = getTPCSplineCalib(track, massProton, hSplineProton[1]);
    }
  }

  if (std::fabs(NSigmaTPC) > nsigmaTPCProtonLc) {
    return false;
  }
  if (track.hasTOF() && std::fabs(NSigmaTOF) > nsigmaTOFProtonLc) {
    return false;
  }

  return true;
}

/// Basic selection of kaon candidates for charm candidates
/// \param track is a track
/// \param nsigmaTPCKaon3Prong max NsigmaTPC for kaon candidates
/// \param nsigmaTOFKaon3Prong max NsigmaTOF for kaon candidates
/// \param setTPCCalib flag to activate TPC PID postcalibrations
/// \param hMapKaon map of nSigma mean and sigma calibrations for kaon
/// \param hSplineKaon spline of kaon and anti-kaon calibrations
/// \return true if track passes all cuts
template <typename T, typename H3>
bool HfFilterHelper::isSelectedKaon4Charm3Prong(const T& track, const float& nsigmaTPCKaon3Prong, const float& nsigmaTOFKaon3Prong, const int setTPCCalib, H3 hMapKaon, const std::array<std::vector<double>, 2>& hSplineKaon)
{
  float NSigmaTPC = track.tpcNSigmaKa();
  float NSigmaTOF = track.tofNSigmaKa();

  if (setTPCCalib == 1) {
    NSigmaTPC = getTPCPostCalib(hMapKaon, track, kKa); // use pion correction map for kaon for the moment
  } else if (setTPCCalib == 2) {
    if (track.sign() > 0) {
      NSigmaTPC = getTPCSplineCalib(track, massK, hSplineKaon[0]);
    } else {
      NSigmaTPC = getTPCSplineCalib(track, massK, hSplineKaon[1]);
    }
  }

  if (std::fabs(NSigmaTPC) > nsigmaTPCKaon3Prong) {
    return false;
  }
  if (track.hasTOF() && std::fabs(NSigmaTOF) > nsigmaTOFKaon3Prong) {
    return false;
  }

  return true;
}

/// Basic additional selection of D+ candidates
/// \param trackOppositeCharge is the opposite charge track
/// \param nsigmaTPCKaon3Prong max NsigmaTPC for kaon candidates
/// \param nsigmaTOFKaon3Prong max NsigmaTOF for kaon candidates
/// \param setTPCCalib flag to activate TPC PID postcalibrations
/// \param hMapKaon map of nSigma mean and sigma calibrations for kaon
/// \param hSplineKaon spline of kaon and anti-kaon calibrations
/// \return BIT(0) for Kpipi
template <typename T, typename H3>
int8_t HfFilterHelper::isDplusPreselected(const T& trackOppositeCharge, const float& nsigmaTPCKaon3Prong, const float& nsigmaTOFKaon3Prong, const int setTPCCalib, H3 hMapKaon, std::array<std::vector<double>, 2>& hSplineKaon)
{
  int8_t retValue = 0;

  // check PID of opposite charge track
  if (!isSelectedKaon4Charm3Prong(trackOppositeCharge, nsigmaTPCKaon3Prong, nsigmaTOFKaon3Prong, setTPCCalib, hMapKaon, hSplineKaon)) {
    return retValue;
  }

  retValue |= BIT(0);
  return retValue;
}

/// Basic additional selection of Ds candidates
/// \param pTrackSameChargeFirst is the first same-charge track momentum
/// \param pTrackSameChargeSecond is the second same-charge track momentum
/// \param pTrackOppositeCharge is the opposite charge track momentum
/// \param trackOppositeCharge is the opposite charge track
/// \param nsigmaTPCKaon3Prong max NsigmaTPC for kaon candidates
/// \param nsigmaTOFKaon3Prong max NsigmaTOF for kaon candidates
/// \param setTPCCalib flag to activate TPC PID postcalibrations
/// \param hMapKaon map of nSigma mean and sigma calibrations for kaon
/// \param hSplineKaon spline of kaon and anti-kaon calibrations
/// \return BIT(0) for KKpi, BIT(1) for piKK
template <typename P, typename T, typename H3>
int8_t HfFilterHelper::isDsPreselected(const P& pTrackSameChargeFirst, const P& pTrackSameChargeSecond, const P& pTrackOppositeCharge, const T& trackOppositeCharge, const float& nsigmaTPCKaon3Prong, const float& nsigmaTOFKaon3Prong, const int setTPCCalib, H3 hMapKaon, const std::array<std::vector<double>, 2>& hSplineKaon)
{
  int8_t retValue = 0;

  // check PID of opposite charge track
  if (!isSelectedKaon4Charm3Prong(trackOppositeCharge, nsigmaTPCKaon3Prong, nsigmaTOFKaon3Prong, setTPCCalib, hMapKaon, hSplineKaon)) {
    return retValue;
  }

  // check delta-mass for phi resonance
  auto invMassKKFirst = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge}, std::array{massK, massK});
  auto invMassKKSecond = RecoDecay::m(std::array{pTrackSameChargeSecond, pTrackOppositeCharge}, std::array{massK, massK});

  if (std::fabs(invMassKKFirst - massPhi) < 0.02) {
    retValue |= BIT(0);
  }
  if (std::fabs(invMassKKSecond - massPhi) < 0.02) {
    retValue |= BIT(1);
  }

  return retValue;
}

/// Basic additional selection of Lc->pKpi and Xic->pKpi candidates
/// \param trackSameChargeFirst is the first same-charge track
/// \param trackSameChargeSecond is the second same-charge track
/// \param trackOppositeCharge is the opposite charge track
/// \param nsigmaTPCProtonLc max NsigmaTPC for proton candidates
/// \param nsigmaTOFProtonLc max NsigmaTOF for proton candidates
/// \param nsigmaTPCKaon3Prong max NsigmaTPC for kaon candidates
/// \param nsigmaTOFKaon3Prong max NsigmaTOF for kaon candidates
/// \param setTPCCalib flag to activate TPC PID postcalibrations
/// \param hMapProton map of nSigma mean and sigma calibrations for proton
/// \param hSplineProton spline of proton and anti-proton calibrations
/// \param hMapKaon map of nSigma mean and sigma calibrations for kaon
/// \param hSplineKaon spline of kaon and anti-kaon calibrations
/// \return BIT(0) for pKpi, BIT(1) for piKp
template <typename T, typename H3>
int8_t HfFilterHelper::isCharmBaryonPreselected(const T& trackSameChargeFirst, const T& trackSameChargeSecond, const T& trackOppositeCharge, const float& nsigmaTPCProtonLc, const float& nsigmaTOFProtonLc, const float& nsigmaTPCKaon3Prong, const float& nsigmaTOFKaon3Prong, const int setTPCCalib, H3 hMapProton, const std::array<std::vector<double>, 2>& hSplineProton, H3 hMapKaon, const std::array<std::vector<double>, 2>& hSplineKaon)
{
  int8_t retValue = 0;
  // check PID of opposite charge track
  if (!isSelectedKaon4Charm3Prong(trackOppositeCharge, nsigmaTPCKaon3Prong, nsigmaTOFKaon3Prong, setTPCCalib, hMapKaon, hSplineKaon)) {
    return retValue;
  }

  if (isSelectedProton4CharmBaryons(trackSameChargeFirst, nsigmaTPCProtonLc, nsigmaTOFProtonLc, setTPCCalib, hMapProton, hSplineProton)) {
    retValue |= BIT(0);
  }
  if (isSelectedProton4CharmBaryons(trackSameChargeSecond, nsigmaTPCProtonLc, nsigmaTOFProtonLc, setTPCCalib, hMapProton, hSplineProton)) {
    retValue |= BIT(1);
  }

  return retValue;
}

/// Basic additional selection of D0 candidates
/// \param trackPos is the positive track
/// \param trackNeg is the negative track
/// \param nsigmaTPCPionKaonDzero max NsigmaTPC for pion/kaon candidates
/// \param nsigmaTOFPionKaonDzero max NsigmaTOF for pion/kaon candidates
/// \param setTPCCalib flag to activate TPC PID postcalibrations
/// \param hMapPion map of nSigma mean and sigma calibrations for pion
/// \param hSplinePion spline of pion and anti-pion calibrations
/// \param hSplineKaon spline of kaon and anti-kaon calibrations
/// \return BIT(0) for D0, BIT(1) for D0bar
template <typename T, typename H3>
int8_t HfFilterHelper::isDzeroPreselected(const T& trackPos, const T& trackNeg, const float& nsigmaTPCPionKaonDzero, const float& nsigmaTOFPionKaonDzero, const int& setTPCCalib, H3 hMapPion, const std::array<std::vector<double>, 2>& hSplinePion, const std::array<std::vector<double>, 2>& hSplineKaon)
{
  int8_t retValue = 0;

  float NSigmaPiTPCPos = trackPos.tpcNSigmaPi();
  float NSigmaPiTOFPos = trackPos.tofNSigmaPi();
  float NSigmaKaTPCPos = trackPos.tpcNSigmaKa();
  float NSigmaKaTOFPos = trackPos.tofNSigmaKa();

  float NSigmaPiTPCNeg = trackNeg.tpcNSigmaPi();
  float NSigmaPiTOFNeg = trackNeg.tofNSigmaPi();
  float NSigmaKaTPCNeg = trackNeg.tpcNSigmaKa();
  float NSigmaKaTOFNeg = trackNeg.tofNSigmaKa();

  if (setTPCCalib == 1) {
    NSigmaPiTPCPos = getTPCPostCalib(hMapPion, trackPos, kPi);
    NSigmaPiTPCNeg = getTPCPostCalib(hMapPion, trackNeg, kPi);
    NSigmaKaTPCPos = getTPCPostCalib(hMapPion, trackPos, kKa); // use pion correction map for kaon for the moment
    NSigmaKaTPCNeg = getTPCPostCalib(hMapPion, trackNeg, kKa); // use pion correction map for kaon for the moment
  } else if (setTPCCalib == 2) {
    NSigmaPiTPCPos = getTPCSplineCalib(trackPos, massPi, hSplinePion[0]);
    NSigmaPiTPCNeg = getTPCSplineCalib(trackNeg, massPi, hSplinePion[1]);
    NSigmaKaTPCPos = getTPCSplineCalib(trackPos, massK, hSplineKaon[0]);
    NSigmaKaTPCNeg = getTPCSplineCalib(trackNeg, massK, hSplineKaon[1]);
  }

  if ((std::fabs(NSigmaPiTPCPos) <= nsigmaTPCPionKaonDzero && (!trackPos.hasTOF() || std::fabs(NSigmaPiTOFPos) <= nsigmaTOFPionKaonDzero)) && (std::fabs(NSigmaKaTPCNeg) <= nsigmaTPCPionKaonDzero && (!trackNeg.hasTOF() || std::fabs(NSigmaKaTOFNeg) <= nsigmaTOFPionKaonDzero))) {
    retValue |= BIT(0);
  }
  if ((std::fabs(NSigmaPiTPCNeg) <= nsigmaTPCPionKaonDzero && (!trackNeg.hasTOF() || std::fabs(NSigmaPiTOFNeg) <= nsigmaTOFPionKaonDzero)) && (std::fabs(NSigmaKaTPCPos) <= nsigmaTPCPionKaonDzero && (!trackPos.hasTOF() || std::fabs(NSigmaKaTOFPos) <= nsigmaTOFPionKaonDzero))) {
    retValue |= BIT(1);
  }

  return retValue;
}

/// Mass selection of D0 candidates to build Bplus candidates
/// \param pTrackPos is the positive track momentum
/// \param pTrackNeg is the negative track momentum
/// \param ptD is the pt of the D0 meson candidate
/// \param isSelected is the flag containing the selection tag for the D0 candidate
/// \param deltaMassCharmHadronForBeauty is the maximum delta mass value (for candidates with pT < 10 GeV/c)
/// \param activateQA flag to activate the filling of QA histos
/// \param hMassVsPt histo with invariant mass vs pt
/// \return 1 for D0, 2 for D0bar, 3 for both
template <typename T, typename H2>
int8_t HfFilterHelper::isSelectedD0InMassRange(const T& pTrackPos, const T& pTrackNeg, const float& ptD, const float& phiD, int8_t isSelected, const float& deltaMassCharmHadronForBeauty, const int& activateQA, H2 hMassVsPt)
{
  int8_t retValue = 0;
  if (TESTBIT(isSelected, 0)) {
    auto invMassD0 = RecoDecay::m(std::array{pTrackPos, pTrackNeg}, std::array{massPi, massK});
    if (activateQA) {
      hMassVsPt->Fill(ptD, invMassD0);
    }
    if (std::fabs(invMassD0 - massD0) < deltaMassCharmHadronForBeauty || ptD > 10) {
      retValue |= BIT(0);
    }
  }
  if (TESTBIT(isSelected, 1)) {
    auto invMassD0bar = RecoDecay::m(std::array{pTrackPos, pTrackNeg}, std::array{massK, massPi});
    if (activateQA) {
      hMassVsPt->Fill(ptD, invMassD0bar);
    }
    if (std::fabs(invMassD0bar - massD0) < deltaMassCharmHadronForBeauty || ptD > 10) {
      retValue |= BIT(1);
    }
  }

  return retValue;
}

/// Mass selection of D+ candidates to build B0 candidates
/// \param pTrackSameChargeFirst is the first same-charge track momentum
/// \param pTrackSameChargeFirst is the second same-charge track momentum
/// \param pTrackSameChargeFirst is the opposite charge track momentum
/// \param ptD is the pt of the D+ meson candidate
/// \param deltaMassCharmHadronForBeauty is the maximum delta mass value (for candidates with pT < 10 GeV/c)
/// \param activateQA flag to activate the filling of QA histos
/// \param hMassVsPt histo with invariant mass vs pt
/// \return BIT(0) (==1) for D+, 0 otherwise
template <typename T, typename H2>
int8_t HfFilterHelper::isSelectedDplusInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD, const float& phiD, const float& deltaMassCharmHadronForBeauty, const int& activateQA, H2 hMassVsPt)
{
  auto invMassDplus = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackSameChargeSecond, pTrackOppositeCharge}, std::array{massPi, massPi, massK});
  if (activateQA) {
    hMassVsPt->Fill(ptD, invMassDplus);
  }

  if (std::fabs(invMassDplus - massDPlus) > deltaMassCharmHadronForBeauty && ptD > 0) {
    return 0;
  }

  return BIT(0);
}

/// Mass selection of of Ds candidates to build Bs candidates
/// \param pTrackSameChargeFirst is the first same-charge track momentum
/// \param pTrackSameChargeFirst is the second same-charge track momentum
/// \param pTrackSameChargeFirst is the opposite charge track momentum
/// \param ptD is the pt of the Ds meson candidate
/// \param isSelected is the flag containing the selection tag for the Ds candidate
/// \param deltaMassCharmHadronForBeauty is the maximum delta mass value (for candidates with pT < 10 GeV/c)
/// \param activateQA flag to activate the filling of QA histos
/// \param hMassVsPt histo with invariant mass vs pt
/// \return BIT(0) for KKpi, BIT(1) for piKK, BIT(2) for phipi, BIT(3) for piphi
template <typename T, typename H2>
int8_t HfFilterHelper::isSelectedDsInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD, const float& phiD, int8_t isSelected, const float& deltaMassCharmHadronForBeauty, const int& activateQA, H2 hMassVsPt)
{
  int8_t retValue = 0;
  if (TESTBIT(isSelected, 0)) {
    auto invMassDsToKKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massK, massK, massPi});
    if (activateQA) {
      hMassVsPt->Fill(ptD, invMassDsToKKPi);
    }
    if (std::fabs(invMassDsToKKPi - massDs) < deltaMassCharmHadronForBeauty || ptD > 10) {
      retValue |= BIT(0);
    }
  }
  if (TESTBIT(isSelected, 1)) {
    auto invMassDsToPiKK = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massK, massK});
    if (activateQA) {
      hMassVsPt->Fill(ptD, invMassDsToPiKK);
    }
    if (std::fabs(invMassDsToPiKK - massDs) < deltaMassCharmHadronForBeauty || ptD > 10) {
      retValue |= BIT(1);
    }
  }

  return retValue;
}

/// Mass selection of Lc candidates to build Lb candidates
/// \param pTrackSameChargeFirst is the first same-charge track momentum
/// \param pTrackSameChargeSecond is the second same-charge track momentum
/// \param pTrackOppositeCharge is the opposite charge track momentum
/// \param ptLc is the pt of the D0 meson candidate
/// \param isSelected is the flag containing the selection tag for the D0 candidate
/// \param deltaMassCharmHadronForBeauty is the maximum delta mass value (for candidates with pT < 10 GeV/c)
/// \param activateQA flag to activate the filling of QA histos
/// \param hMassVsPt histo with invariant mass vs pt
/// \return BIT(0) for pKpi with mass cut, BIT(1) for piKp with mass cut
template <typename T, typename H2>
int8_t HfFilterHelper::isSelectedLcInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptLc, const float& phiLc, const int8_t isSelected, const float& deltaMassCharmHadronForBeauty, const int& activateQA, H2 hMassVsPt)
{
  int8_t retValue = 0;
  if (TESTBIT(isSelected, 0)) {
    auto invMassLcToPKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massProton, massK, massPi});
    if (activateQA) {
      hMassVsPt->Fill(ptLc, invMassLcToPKPi);
    }
    if (std::fabs(invMassLcToPKPi - massLc) < deltaMassCharmHadronForBeauty || ptLc > 10) {
      retValue |= BIT(0);
    }
  }
  if (TESTBIT(isSelected, 1)) {
    auto invMassLcToPiKP = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massK, massProton});
    if (activateQA) {
      hMassVsPt->Fill(ptLc, invMassLcToPiKP);
    }
    if (std::fabs(invMassLcToPiKP - massLc) < deltaMassCharmHadronForBeauty || ptLc > 10) {
      retValue |= BIT(1);
    }
  }

  return retValue;
}

/// Mass selection of Xic candidates to build Lb candidates
/// \param pTrackSameChargeFirst is the first same-charge track momentum
/// \param pTrackSameChargeSecond is the second same-charge track momentum
/// \param pTrackOppositeCharge is the opposite charge track momentum
/// \param ptXic is the pt of the D0 meson candidate
/// \param isSelected is the flag containing the selection tag for the D0 candidate
/// \param deltaMassCharmHadronForBeauty is the maximum delta mass value (for candidates with pT < 10 GeV/c)
/// \param activateQA flag to activate the filling of QA histos
/// \param hMassVsPt histo with invariant mass vs pt
/// \return BIT(0) for pKpi with mass cut, BIT(1) for piKp with mass cut
template <typename T, typename H2>
int8_t HfFilterHelper::isSelectedXicInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptXic, const float& phiXic, const int8_t isSelected, const float& deltaMassCharmHadronForBeauty, const int& activateQA, H2 hMassVsPt)
{
  int8_t retValue = 0;
  if (TESTBIT(isSelected, 0)) {
    auto invMassXicToPKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massProton, massK, massPi});
    if (activateQA) {
      hMassVsPt->Fill(ptXic, invMassXicToPKPi);
    }
    if (std::fabs(invMassXicToPKPi - massXic) < deltaMassCharmHadronForBeauty || ptXic > 10) {
      retValue |= BIT(0);
    }
  }
  if (TESTBIT(isSelected, 1)) {
    auto invMassXicToPiKP = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massK, massProton});
    if (activateQA) {
      hMassVsPt->Fill(ptXic, invMassXicToPiKP);
    }
    if (std::fabs(invMassXicToPiKP - massXic) < deltaMassCharmHadronForBeauty || ptXic > 10) {
      retValue |= BIT(1);
    }
  }

  return retValue;
}

/// Basic selection of V0 candidates
/// \param v0 is the v0 candidate
/// \param dauTracks is a 2-element array with positive and negative V0 daughter tracks
/// \param collision is the current collision
/// \param minGammaCosinePa is the minimum required cosp of the gamma
/// \param minV0CosinePa is the minimum required cosp of K0S/Lambda
/// \param minV0Radius is the minimum required K0S/Lambda radius
/// \param maxNsigmaPrForLambda is the maximum allowed nSigma TPC/TOF for protons in Lambda decays (applied only if PID info available)
/// \param deltaMassK0s is the maximum allowed delta mass for K0S
/// \param deltaMassLambda is the maximum allowed delta mass for Lambda
/// \param setTPCCalib flag to activate TPC PID postcalibrations
/// \param hMapProton map of nSigma mean and sigma calibrations for proton
/// \param hSplineProton spline of proton and anti-proton calibrations
/// \param activateQA flag to fill QA histos
/// \param hV0Selected is the pointer to the QA histo for selected gammas
/// \param hArmPod is the pointer to an array of QA histo AP plot before selection
/// \return an integer passes all cuts
template <typename V0, typename Coll, typename T, typename H2, typename H3>
int8_t HfFilterHelper::isSelectedV0(const V0& v0, const array<T, 2>& dauTracks, const Coll& collision, const float& minGammaCosinePa, const float& minV0CosinePa, const float& minV0Radius, const float& maxNsigmaPrForLambda, const float& deltaMassK0s, const float& deltaMassLambda, const int& setTPCCalib, H3 hMapProton, const std::array<std::vector<double>, 2>& hSplineProton, const int& activateQA, H2 hV0Selected, std::array<H2, 4>& hArmPod)
{
  int8_t isSelected{BIT(kPhoton) | BIT(kK0S) | BIT(kLambda) | BIT(kAntiLambda)};

  if (activateQA > 1) {
    for (int iV0{kPhoton}; iV0 < kNV0; ++iV0) {
      hV0Selected->Fill(0., iV0);
    }
  }

  // eta of daughters
  if (std::fabs(dauTracks[0].eta()) > 1. || std::fabs(dauTracks[1].eta()) > 1.) { // cut all V0 daughters with |eta| > 1.
    if (activateQA > 1) {
      for (int iV0{kPhoton}; iV0 < kNV0; ++iV0) {
        hV0Selected->Fill(1., iV0);
      }
    }
    return kRejected;
  }

  // V0 radius
  if (v0.v0radius() < 0. || v0.v0radius() > 180.) {
    CLRBIT(isSelected, kPhoton);
    if (activateQA > 1) {
      hV0Selected->Fill(2., kPhoton);
    }
  }
  if (v0.v0radius() < minV0Radius) {
    for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(2., iV0);
      }
    }
  }

  auto v0CosinePa = v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
  // cosine of pointing angle
  if (TESTBIT(isSelected, kPhoton) && v0CosinePa < minGammaCosinePa) {
    CLRBIT(isSelected, kPhoton);
    if (activateQA > 1) {
      hV0Selected->Fill(3., kPhoton);
    }
  }
  for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
    if (TESTBIT(isSelected, iV0) && v0CosinePa < minV0CosinePa) {
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(3., iV0);
      }
    }
  }

  // armenteros-podolanski / mass
  if (TESTBIT(isSelected, kPhoton) && (std::pow(v0.alpha() / 0.95, 2) + std::pow(v0.qtarm() / 0.05, 2)) >= 1) {
    CLRBIT(isSelected, kPhoton);
    if (activateQA > 1) {
      hV0Selected->Fill(4., kPhoton);
    }
  }
  if (TESTBIT(isSelected, kK0S) && std::fabs(v0.mK0Short() - massK0S) > deltaMassK0s) {
    CLRBIT(isSelected, kK0S);
    if (activateQA > 1) {
      hV0Selected->Fill(4., kK0S);
    }
  }
  if (TESTBIT(isSelected, kLambda) && std::fabs(v0.mLambda() - massLambda) > deltaMassLambda) {
    CLRBIT(isSelected, kLambda);
    if (activateQA > 1) {
      hV0Selected->Fill(4., kLambda);
    }
  }
  if (TESTBIT(isSelected, kAntiLambda) && std::fabs(v0.mAntiLambda() - massLambda) > deltaMassLambda) {
    CLRBIT(isSelected, kAntiLambda);
    if (activateQA > 1) {
      hV0Selected->Fill(4., kAntiLambda);
    }
  }

  // DCA V0 and V0 daughters
  for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
    if (TESTBIT(isSelected, iV0) && v0.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > 0.1f) { // we want only primary V0s
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(5., iV0);
      }
    }
    if (TESTBIT(isSelected, iV0) && (v0.dcaV0daughters() > 1.f || std::fabs(v0.dcapostopv()) < 0.05f || std::fabs(v0.dcanegtopv()) < 0.05f)) {
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(6., iV0);
      }
    }
  }

  // psi pair (photon only)
  if (TESTBIT(isSelected, kPhoton) && std::fabs(v0.psipair()) > 0.1) {
    CLRBIT(isSelected, kPhoton);
    if (activateQA > 1) {
      hV0Selected->Fill(7., kPhoton);
    }
  }

  // PID (Lambda/AntiLambda only)
  float nSigmaPrTpc[2] = {dauTracks[0].tpcNSigmaPr(), dauTracks[1].tpcNSigmaPr()};
  float nSigmaPrTof[2] = {dauTracks[0].tofNSigmaPr(), dauTracks[1].tofNSigmaPr()};
  if (setTPCCalib == 1) {
    for (int iDau{0}; iDau < 2; ++iDau) {
      nSigmaPrTpc[iDau] = getTPCPostCalib(hMapProton, dauTracks[iDau], kPr);
    }
  } else if (setTPCCalib == 2) {
    for (int iDau{0}; iDau < 2; ++iDau) {
      nSigmaPrTpc[iDau] = getTPCSplineCalib(dauTracks[iDau], massProton, hSplineProton[iDau]);
    }
  }

  if (TESTBIT(isSelected, kLambda) && ((dauTracks[0].hasTPC() && std::fabs(nSigmaPrTpc[0]) > maxNsigmaPrForLambda) || (dauTracks[0].hasTOF() && std::fabs(nSigmaPrTof[0]) > maxNsigmaPrForLambda))) {
    CLRBIT(isSelected, kLambda);
    if (activateQA > 1) {
      hV0Selected->Fill(8., kLambda);
    }
  }
  if (TESTBIT(isSelected, kAntiLambda) && ((dauTracks[1].hasTPC() && std::fabs(nSigmaPrTpc[1]) > maxNsigmaPrForLambda) || (dauTracks[1].hasTOF() && std::fabs(nSigmaPrTof[1]) > maxNsigmaPrForLambda))) {
    CLRBIT(isSelected, kAntiLambda);
    if (activateQA > 1) {
      hV0Selected->Fill(8., kAntiLambda);
    }
  }

  if (activateQA) {
    for (int iV0{kPhoton}; iV0 < kNV0; ++iV0) {
      if (TESTBIT(isSelected, iV0)) {
        hArmPod[iV0]->Fill(v0.alpha(), v0.qtarm());
        if (activateQA > 1) {
          hV0Selected->Fill(9., iV0);
        }
      }
    }
  }

  return isSelected;
}

/// Basic selection of cascade candidates
/// \param casc is the cascade candidate
/// \param v0 is the cascade daughter
/// \param dauTracks is a 3-element array with bachelor, positive and negative V0 daughter tracks
/// \param collision is the collision
/// \param minPtXiBachelor is the minimum required pT for the cascade bachelor
/// \param deltaMassXi is the maximum delta mass for the Xi
/// \param deltaMassLambda is the maximum delta mass for the Lambda daughter
/// \param cosPAXi is the minimum value of cosPA for the cascade
/// \param cosPALambda is the minimum value of cosPA for the Lambda daughter
/// \param DCAxyXi is the maximum DCAxy of the Xi to the PV
/// \param maxNsigma is the maximum number of sigma to accept a given PID hypothesis
/// \param setTPCCalib flag to activate TPC PID postcalibrations
/// \param hMapProton map of nSigma mean and sigma calibrations for proton
/// \param hMapPion map of nSigma mean and sigma calibrations for pion
/// \param hSplineProton spline of proton and anti-proton calibrations
/// \param hSplinePion spline of pion and anti-pion calibrations
/// \return true if cascade passes all cuts
template <typename Casc, typename V0, typename T, typename Coll, typename H3>
bool HfFilterHelper::isSelectedCascade(const Casc& casc, const V0& v0, const array<T, 3>& dauTracks, const Coll& collision, const float& minPtXiBachelor, const float& deltaMassXi, const float& deltaMassLambda, const float& cosPAXi, const float& cosPALambda, const float& DCAxyXi, const float& maxNsigma, const int& setTPCCalib, H3 hMapPion, H3 hMapProton, const std::array<std::vector<double>, 2>& hSplinePion, const std::array<std::vector<double>, 2>& hSplineProton)
{
  // eta of daughters
  if (std::fabs(dauTracks[0].eta()) > 1. || std::fabs(dauTracks[1].eta()) > 1. || std::fabs(dauTracks[2].eta()) > 1.) { // cut all V0 daughters with |eta| > 1.
    return false;
  }

  // V0 radius
  if (v0.v0radius() < 1.2) {
    return false;
  }

  // cascade radius
  if (casc.cascradius() < 0.6) {
    return false;
  }

  // V0 cosp
  if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cosPALambda) {
    return false;
  }

  // cascade cosp
  if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cosPAXi) {
    return false;
  }

  // cascade DCAxy to PV
  if (std::fabs(casc.dcaXYCascToPV()) > DCAxyXi) {
    return false;
  }

  // Xi bachelor min pT
  if (dauTracks[0].pt() < minPtXiBachelor) {
    return false;
  }

  // dau dca
  if (std::fabs(casc.dcaV0daughters()) > 1.f || std::fabs(casc.dcacascdaughters()) > 1.f) {
    return false;
  }

  // cascade mass
  if (std::fabs(casc.mXi() - massXi) > deltaMassXi) {
    return false;
  }

  // V0 mass
  if (std::fabs(casc.mLambda() - massLambda) > deltaMassLambda) {
    return false;
  }

  // PID
  float nSigmaPrTpc[3] = {-999., dauTracks[1].tpcNSigmaPr(), dauTracks[2].tpcNSigmaPr()};
  float nSigmaPrTof[3] = {-999., dauTracks[1].tofNSigmaPr(), dauTracks[2].tofNSigmaPr()};
  float nSigmaPiTpc[3] = {dauTracks[0].tpcNSigmaPi(), dauTracks[1].tpcNSigmaPi(), dauTracks[2].tpcNSigmaPi()};
  float nSigmaPiTof[3] = {dauTracks[0].tofNSigmaPi(), dauTracks[1].tofNSigmaPi(), dauTracks[2].tofNSigmaPi()};
  if (setTPCCalib == 1) {
    for (int iDau{0}; iDau < 3; ++iDau) {
      nSigmaPiTpc[iDau] = getTPCPostCalib(hMapPion, dauTracks[iDau], kPi);
      if (iDau == 0) {
        continue;
      }
      nSigmaPrTpc[iDau] = getTPCPostCalib(hMapProton, dauTracks[iDau], kPr);
    }
  } else if (setTPCCalib == 2) {
    for (int iDau{0}; iDau < 3; ++iDau) {
      nSigmaPiTpc[iDau] = getTPCSplineCalib(dauTracks[iDau], massPi, (dauTracks[iDau].sign() > 0) ? hSplinePion[0] : hSplinePion[1]);
      if (iDau == 0) {
        continue;
      }
      nSigmaPrTpc[iDau] = getTPCSplineCalib(dauTracks[iDau], massProton, (dauTracks[iDau].sign() > 0) ? hSplineProton[0] : hSplineProton[1]);
    }
  }

  // PID to V0 tracks
  if (dauTracks[0].sign() < 0) { // Xi-
    if ((dauTracks[1].hasTPC() && std::fabs(nSigmaPrTpc[1]) > maxNsigma) && (dauTracks[1].hasTOF() && std::fabs(nSigmaPrTof[1]) > maxNsigma)) {
      return false;
    }
    if ((dauTracks[2].hasTPC() && std::fabs(nSigmaPiTpc[2]) > maxNsigma) && (dauTracks[2].hasTOF() && std::fabs(nSigmaPiTof[2]) > maxNsigma)) {
      return false;
    }
  } else if (dauTracks[0].sign() > 0) { // Xi+
    if ((dauTracks[2].hasTPC() && std::fabs(nSigmaPrTpc[2]) > maxNsigma) && (dauTracks[2].hasTOF() && std::fabs(nSigmaPrTof[2]) > maxNsigma)) {
      return false;
    }
    if ((dauTracks[1].hasTPC() && std::fabs(nSigmaPiTpc[1]) > maxNsigma) && (dauTracks[1].hasTOF() && std::fabs(nSigmaPiTof[1]) > maxNsigma)) {
      return false;
    }
  }

  // bachelor PID
  if ((dauTracks[0].hasTPC() && std::fabs(nSigmaPiTpc[0]) > maxNsigma) && (dauTracks[0].hasTOF() && std::fabs(nSigmaPiTof[0]) > maxNsigma)) {
    return false;
  }

  // additional track cuts
  for (auto const& dauTrack : dauTracks) {
    //  TPC clusters selections
    if (dauTrack.tpcNClsFound() < 70) { // TODO: put me as a configurable please
      return false;
    }
    if (dauTrack.tpcNClsCrossedRows() < 70) {
      return false;
    }
    if (dauTrack.tpcCrossedRowsOverFindableCls() < 0.8) {
      return false;
    }
  }

  return true;
}

/// Single-track cuts for bachelor track of charm baryon candidates
/// \param track is a track
/// \param dca is the 2d array with dcaXY and dcaZ of the track
/// \param minPt is the minimum pT
/// \param pTBinsTrack pt bins for DCA cuts
/// \param cutsSingleTrack cuts for all tracks
/// \param maxNsigmaTPC is the maximum nSigma TPC for pions and kaons
/// \param maxNsigmaTOF is the maximum nSigma TOF for pions and kaons
/// \param setTPCCalib flag to activate TPC PID postcalibrations
/// \param hMapPion map of nSigma mean and sigma calibrations for pion
/// \param hSplinePion spline of pion and anti-pion calibrations
/// \param hSplineKaon spline of kaon and anti-kaon calibrations
/// \return 0 if rejected, or a bitmap that contains the information whether it is selected as pion and/or kaon
template <typename T, typename T2, typename T3, typename T4, typename H3>
int8_t HfFilterHelper::isSelectedBachelorForCharmBaryon(const T& track, const T2& dca, const float& minPt, const T3& pTBinsTrack, const T4& cutsSingleTrack, const float& maxNsigmaTPC, const float& maxNsigmaTOF, const int& setTPCCalib, H3 hMapPion, const std::array<std::vector<double>, 2>& hSplinePion, const std::array<std::vector<double>, 2>& hSplineKaon)
{
  int8_t retValue{BIT(kPionForCharmBaryon) | BIT(kKaonForCharmBaryon)};

  if (!track.isGlobalTrackWoDCA()) {
    return kRejected;
  }

  if (track.pt() < minPt) {
    return kRejected;
  }

  auto pTBinTrack = findBin(pTBinsTrack, track.pt());
  if (pTBinTrack == -1) {
    return kRejected;
  }

  if (std::fabs(dca[0]) < cutsSingleTrack.get(pTBinTrack, 0u)) {
    return kRejected; // minimum DCAxy
  }
  if (std::fabs(dca[0]) > cutsSingleTrack.get(pTBinTrack, 1u)) {
    return kRejected; // maximum DCAxy
  }

  if (std::fabs(dca[1]) > 2.f) {
    return kRejected; // maximum DCAz
  }

  if (track.tpcNClsFound() < 70) {
    return kRejected;
  }

  if (track.itsNCls() < 3) {
    return kRejected;
  }

  float nSigmaPiTpc = track.tpcNSigmaPi();
  float nSigmaKaTpc = track.tpcNSigmaKa();
  float nSigmaPiTof = track.tofNSigmaPi();
  float nSigmaKaTof = track.tofNSigmaKa();
  if (setTPCCalib == 1) {
    nSigmaPiTpc = getTPCPostCalib(hMapPion, track, kPi);
    nSigmaKaTpc = getTPCPostCalib(hMapPion, track, kKa);
  } else if (setTPCCalib == 2) {
    nSigmaPiTpc = getTPCSplineCalib(track, massPi, (track.sign() > 0) ? hSplinePion[0] : hSplinePion[1]);
    nSigmaKaTpc = getTPCSplineCalib(track, massK, (track.sign() > 0) ? hSplineKaon[0] : hSplineKaon[1]);
  }

  if ((track.hasTPC() && std::fabs(nSigmaPiTpc) > maxNsigmaTPC) && (track.hasTOF() && std::fabs(nSigmaPiTof) > maxNsigmaTOF)) {
    CLRBIT(retValue, kPionForCharmBaryon);
  }
  if ((track.hasTPC() && std::fabs(nSigmaKaTpc) > maxNsigmaTPC) && (track.hasTOF() && std::fabs(nSigmaKaTof) > maxNsigmaTOF)) {
    CLRBIT(retValue, kKaonForCharmBaryon);
  }

  return retValue;
}

/// BDT selections
/// \param scores is a 3-element array with BDT out scores
/// \param thresholdBDTScores is the LabelledArray containing the BDT cut values
/// \return 0 if rejected, otherwise bitmap with BIT(RecoDecay::OriginType::Prompt) and/or BIT(RecoDecay::OriginType::NonPrompt) on
template <typename T, typename U>
int8_t HfFilterHelper::isBDTSelected(const T& scores, const U& thresholdBDTScores)
{
  int8_t retValue = 0;
  if (scores.size() < 3) {
    return retValue;
  }

  if (scores[0] > thresholdBDTScores.get(0u, 0u)) {
    return retValue;
  }
  if (scores[1] > thresholdBDTScores.get(0u, 1u)) {
    retValue |= BIT(RecoDecay::OriginType::Prompt);
  }
  if (scores[2] > thresholdBDTScores.get(0u, 2u)) {
    retValue |= BIT(RecoDecay::OriginType::NonPrompt);
  }

  return retValue;
}

/// Computation of the relative momentum between particle pairs
/// \param pTrack is the track momentum array
/// \param ProtonMass is the mass of a proton
/// \param CharmCandMomentum is the three momentum of a charm candidate
/// \param CharmMass is the mass of the charm hadron
/// \return relative momentum of pair
template <typename T>
T HfFilterHelper::computeRelativeMomentum(const std::array<T, 3>& pTrack, const std::array<T, 3>& CharmCandMomentum, const T& CharmMass)
{
  ROOT::Math::PxPyPzMVector part1(pTrack[0], pTrack[1], pTrack[2], massProton);
  ROOT::Math::PxPyPzMVector part2(CharmCandMomentum[0], CharmCandMomentum[1], CharmCandMomentum[2], CharmMass);

  ROOT::Math::PxPyPzMVector trackSum = part1 + part2;
  ROOT::Math::Boost boostv12{trackSum.BoostToCM()};
  ROOT::Math::PxPyPzMVector part1CM = boostv12(part1);
  ROOT::Math::PxPyPzMVector part2CM = boostv12(part2);
  ROOT::Math::PxPyPzMVector trackRelK = part1CM - part2CM;

  T kStar = 0.5 * trackRelK.P();
  return kStar;
} // float computeRelativeMomentum(const T& track, const std::array<float, 3>& CharmCandMomentum, const float& CharmMass)

/// Computation of the number of candidates in an event that do not share daughter tracks
/// \return 0 or 1 in case of less than 2 independent candidates in a single event, 2 otherwise
template <typename T>
int HfFilterHelper::computeNumberOfCandidates(std::vector<std::vector<T>> indices)
{
  if (indices.size() < 2) {
    return indices.size();
  }

  std::vector<int> numIndependentCand{};
  for (auto iCand{0u}; iCand < indices.size(); ++iCand) {
    int nIndependent = 0;
    for (auto iCandSecond{0u}; iCandSecond < indices.size(); ++iCandSecond) {
      if (iCand == iCandSecond) {
        continue;
      } else {
        bool hasOverlap = false;
        for (auto idxFirst{0u}; idxFirst < indices[iCand].size(); ++idxFirst) {
          for (auto idxSecond{0u}; idxSecond < indices[iCandSecond].size(); ++idxSecond) {
            if (indices[iCand][idxFirst] == indices[iCandSecond][idxSecond]) {
              hasOverlap = true;
              break;
            }
          }
        }
        if (!hasOverlap) {
          nIndependent++;
        }
      }
    }
    numIndependentCand.push_back(nIndependent);
  }
  std::sort(numIndependentCand.begin(), numIndependentCand.end());

  if (numIndependentCand.back() == 0) {
    return numIndependentCand.back();
  }

  return 2;
}

/// ML helper methods

/// Iinitialisation of ONNX session
/// \param onnxFile is the onnx file name
/// \param partName is the particle name
/// \param env is the ONNX environment
/// \param sessionOpt is the ONNX session options
/// \param inputShapes is the input shape
/// \param dataType is the data type (1=float, 11=double)
/// \param loadModelsFromCCDB is the flag to decide whether the ONNX file is read from CCDB or not
/// \param ccdbApi is the CCDB API
/// \param mlModelPathCCDB is the model path in CCDB
/// \param timestampCCDB is the CCDB timestamp
/// \return the pointer to the ONNX Ort::Experimental::Session
Ort::Experimental::Session* HfFilterHelper::initONNXSession(std::string& onnxFile, std::string partName, Ort::Env& env, Ort::SessionOptions& sessionOpt, std::vector<std::vector<int64_t>>& inputShapes, int& dataType, bool loadModelsFromCCDB, o2::ccdb::CcdbApi& ccdbApi, std::string mlModelPathCCDB, int64_t timestampCCDB)
{
  // hard coded, we do not let the user change this
  sessionOpt.SetIntraOpNumThreads(1);
  sessionOpt.SetInterOpNumThreads(1);
  Ort::Experimental::Session* session = nullptr;

  std::map<std::string, std::string> metadata;
  bool retrieveSuccess = true;
  if (loadModelsFromCCDB) {
    retrieveSuccess = ccdbApi.retrieveBlob(mlModelPathCCDB + partName, ".", metadata, timestampCCDB, false, onnxFile);
  }
  if (retrieveSuccess) {
    session = new Ort::Experimental::Session{env, onnxFile, sessionOpt};
    inputShapes = session->GetInputShapes();
    if (inputShapes[0][0] < 0) {
      LOGF(warning, Form("Model for %s with negative input shape likely because converted with hummingbird, setting it to 1.", partName.data()));
      inputShapes[0][0] = 1;
    }

    Ort::TypeInfo typeInfo = session->GetInputTypeInfo(0);
    auto tensorInfo = typeInfo.GetTensorTypeAndShapeInfo();
    dataType = tensorInfo.GetElementType();
  } else {
    LOG(fatal) << "Error encountered while fetching/loading the ML model from CCDB! Maybe the ML model doesn't exist yet for this runnumber/timestamp?";
  }

  return session;
}

/// Iinitialisation of ONNX session
/// \param inputFeatures is the vector with input features
/// \param session is the ONNX Ort::Experimental::Session
/// \param inputShapes is the input shape
/// \return the array with the three output scores
template <typename T>
std::array<T, 3> HfFilterHelper::predictONNX(std::vector<T>& inputFeatures, std::shared_ptr<Ort::Experimental::Session>& session, std::vector<std::vector<int64_t>>& inputShapes)
{
  std::array<T, 3> scores{-1., 2., 2.};
  std::vector<Ort::Value> inputTensor{};
  inputTensor.push_back(Ort::Experimental::Value::CreateTensor<T>(inputFeatures.data(), inputFeatures.size(), inputShapes[0]));

  // double-check the dimensions of the input tensor
  if (inputTensor[0].GetTensorTypeAndShapeInfo().GetShape()[0] > 0) { // vectorial models can have negative shape if the shape is unknown
    assert(inputTensor[0].IsTensor() && inputTensor[0].GetTensorTypeAndShapeInfo().GetShape() == inputShapes[0]);
  }
  try {
    auto outputTensor = session->Run(session->GetInputNames(), inputTensor, session->GetOutputNames());
    assert(outputTensor.size() == session->GetOutputNames().size() && outputTensor[1].IsTensor());
    auto typeInfo = outputTensor[1].GetTensorTypeAndShapeInfo();
    assert(typeInfo.GetElementCount() == 3); // we need multiclass
    scores[0] = outputTensor[1].GetTensorMutableData<T>()[0];
    scores[1] = outputTensor[1].GetTensorMutableData<T>()[1];
    scores[2] = outputTensor[1].GetTensorMutableData<T>()[2];
  } catch (const Ort::Exception& exception) {
    LOG(error) << "Error running model inference: " << exception.what();
  }

  return scores;
}

/// PID postcalibrations

/// load the TPC spline from the CCDB
/// \param ccdbApi is Api for CCDB
/// \param bunchCrossing is the timestamp of bunchcrossing for the run number
/// \param ccdbPath  is the path on CCDB
/// \return a vector include parameters for BetheBloch formula
std::vector<double> HfFilterHelper::setValuesBB(o2::ccdb::CcdbApi& ccdbApi, aod::BCsWithTimestamps::iterator const& bunchCrossing, const std::string ccdbPath)
{
  map<string, string> metadata;
  auto hSpline = ccdbApi.retrieveFromTFileAny<TH1F>(ccdbPath, metadata, bunchCrossing.timestamp());

  if (!hSpline) {
    LOG(fatal) << "File from CCDB in path " << ccdbPath << " was not found for run " << bunchCrossing.runNumber();
  }

  TAxis* axis = hSpline->GetXaxis();
  std::vector<double> parsBB{static_cast<double>(hSpline->GetBinContent(axis->FindBin("bb1"))),
                             static_cast<double>(hSpline->GetBinContent(axis->FindBin("bb2"))),
                             static_cast<double>(hSpline->GetBinContent(axis->FindBin("bb3"))),
                             static_cast<double>(hSpline->GetBinContent(axis->FindBin("bb4"))),
                             static_cast<double>(hSpline->GetBinContent(axis->FindBin("bb5"))),
                             static_cast<double>(hSpline->GetBinContent(axis->FindBin("Resolution")))};
  return parsBB;
}

/// compute TPC postcalibrated nsigma based on calibration histograms from CCDB
/// \param hCalibMean calibration histograms of mean from CCDB
/// \param hCalibSigma calibration histograms of sigma from CCDB
/// \param track is the track
/// \param pidSpecies is the PID species
/// \return the corrected Nsigma value for the PID species
template <typename T, typename H3>
float HfFilterHelper::getTPCPostCalib(const array<H3, 2>& hCalibMap, const T& track, const int pidSpecies)
{
  auto tpcNCls = track.tpcNClsFound();
  auto tpcPin = track.tpcInnerParam();
  auto eta = track.eta();
  auto tpcNSigma = 0.;

  if (pidSpecies == kKa) {
    tpcNSigma = track.tpcNSigmaKa();
  } else if (pidSpecies == kPi) {
    tpcNSigma = track.tpcNSigmaPi();
  } else if (pidSpecies == kPr) {
    tpcNSigma = track.tpcNSigmaPr();
  } else {
    LOG(fatal) << "Wrong PID Species be selected, please check!";
  }
  auto binTPCNCls = hCalibMap[0]->GetXaxis()->FindBin(tpcNCls);
  binTPCNCls = (binTPCNCls == 0 ? 1 : binTPCNCls);
  binTPCNCls = std::min(hCalibMap[0]->GetXaxis()->GetNbins(), binTPCNCls);
  auto binPin = hCalibMap[0]->GetYaxis()->FindBin(tpcPin);
  binPin = (binPin == 0 ? 1 : binPin);
  binPin = std::min(hCalibMap[0]->GetYaxis()->GetNbins(), binPin);
  auto binEta = hCalibMap[0]->GetZaxis()->FindBin(eta);
  binEta = (binEta == 0 ? 1 : binEta);
  binEta = std::min(hCalibMap[0]->GetZaxis()->GetNbins(), binEta);

  auto mean = hCalibMap[0]->GetBinContent(binTPCNCls, binPin, binEta);
  auto width = hCalibMap[1]->GetBinContent(binTPCNCls, binPin, binEta);

  return (tpcNSigma - mean) / width;
}
}
}