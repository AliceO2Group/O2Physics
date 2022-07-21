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

/// \file tpcSkimsTableCreator.cxx
/// \brief Task to produce table with clean selections for TPC PID calibration
///
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>
/// \author Christian Sonnabend <christian.sonnabend@cern.ch>
/// \author Jeremy Wilkinson <jeremy.wilkinson@cern.ch>
/// O2
#include <CCDB/BasicCCDBManager.h>
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
/// O2Physics
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/StrangenessTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"

#include "tpcSkimsTableCreator.h"
/// ROOT
#include "TRandom3.h"
#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::dataformats;

struct TreeWriterTpcV0 {

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TrackSelection>;
  using Coll = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;

  /// Tables to be produced
  Produces<o2::aod::SkimmedTPCV0Tree> rowTPCTree;

  /// Configurables
  Configurable<float> cutPAV0{"cutPAV0", 0., "Cut on the cos(pointing angle) of the decay"};
  Configurable<float> nSigmaTOFdautrack{"nSigmaTOFdautrack", 5., "n-sigma TOF cut on the daughter tracks. Set 0 to switch it off."};
  Configurable<float> nClNorm{"nClNorm", 152., "Number of cluster normalization. Run 2: 159, Run 3 152"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  /// Configurables downsampling
  Configurable<double> dwnSmplFactor_Pi{"dwnSmplFactor_Pi", 1., "downsampling factor for pions, default fraction to keep is 1."};
  Configurable<double> dwnSmplFactor_Pr{"dwnSmplFactor_Pr", 1., "downsampling factor for protons, default fraction to keep is 1."};
  Configurable<double> dwnSmplFactor_El{"dwnSmplFactor_El", 1., "downsampling factor for electrons, default fraction to keep is 1."};
  Configurable<float> sqrtSNN{"sqrt_s_NN", 0., "sqrt(s_NN), used for downsampling with the Tsallis distribution"};
  Configurable<float> downsamplingTsalisPions{"downsamplingTsalisPions", -1., "Downsampling factor to reduce the number of pions"};
  Configurable<float> downsamplingTsalisProtons{"downsamplingTsalisProtons", -1., "Downsampling factor to reduce the number of protons"};
  Configurable<float> downsamplingTsalisElectrons{"downsamplingTsalisElectrons", -1., "Downsampling factor to reduce the number of electrons"};
  /// Configurables kaon
  Configurable<float> invariantMassCutK0Short{"invariantMassCutK0Short", 0.5, "Mass cut for K0short"};
  Configurable<float> cutQTK0min{"cutQTK0min", 0.1075, "Minimum qt for K0short"};
  Configurable<float> cutQTK0max{"cutQTK0max", 0.215, "Minimum qt for K0short"};
  Configurable<float> cutAPK0pinner{"cutAPK0pinner", 0.199, "inner momentum for K0short"};
  Configurable<float> cutAPK0pouter{"cutAPK0pouter", 0.215, "outer momentum for K0short"};
  Configurable<float> cutAPK0alinner{"cutAPK0alinner", 0.8, "inner alpha for K0short"};
  Configurable<float> cutAPK0alouter{"cutAPK0alouter", 0.9, "outer alpha for K0short"};
  /// Configurables lambda/antilambda
  Configurable<float> invariantMassCutLambda{"invariantMassCutLambda", 0.5, "Mass cut for Lambda and Anti-Lambda"};
  Configurable<float> cutQTmin{"cutQTmin", 0.03, "the minimum value of qt for Lambda"};
  Configurable<float> cutAlphaminL{"cutAlphaminL", 0.35, "minimum alpha for Lambda"};
  Configurable<float> cutAlphamaxL{"cutAlphamaxL", 0.7, "maximum alpha for Lambda"};
  Configurable<float> cutAlphaminAL{"cutAlphaminAL", -0.7, "minimum alpha for Anti-Lambda"};
  Configurable<float> cutAlphamaxAL{"cutAlphamaxAL", -0.35, "maximum alpha for Anti-Lambda"};
  Configurable<float> cutAPL0{"cutAPL0", 0.10, "cutAPL par0"};
  Configurable<float> cutAPL1{"cutAPL1", 0.69, "cutAPL par1"};
  Configurable<float> cutAPL2{"cutAPL2", 0.5, "cutAPL par2"};
  /// Configurables gamma
  Configurable<float> gammaAsymmetryMax{"gammaAsymmetryMax", 0.95, "maximum photon asymetry."};
  Configurable<float> gammaPsiPairMax{"gammaPsiPairMax", 0.1, "maximum psi angle of the track pair. "};
  Configurable<float> gammaQtPtMultiplicator{"gammaQtPtMultiplicator", 0.11, "Multiply pt of V0s by this value to get the 2nd denominator in the armenteros cut. The products maximum value is fV0QtMax."};
  Configurable<float> gammaQtMax{"gammaQtMax", 0.040, "the maximum value of the product, that is the maximum qt"};
  Configurable<float> gammaV0RadiusMin{"gammaV0RadiusMin", 5., "the minimum V0 radius of the gamma decay"};
  Configurable<float> gammaV0RadiusMax{"gammaV0RadiusMax", 180., "the maximum V0 radius of the gamma decay"};
  Configurable<float> cutAlphaG1{"cutAlphaG1", 0.35, "maximum alpha gamma decay cut 1"};
  Configurable<float> cutAlphaG2min{"cutAlphaG2min", 0.6, "minimum alpha gamma decay cut 2"};
  Configurable<float> cutAlphaG2max{"cutAlphaG2max", 0.8, "maximum alpha gamma decay cut 2"};
  Configurable<float> cutQTG{"cutQTG", 0.04, "maximum qt gamma decay"};

  /// Kaon selection
  template <typename C, typename V0>
  bool selectionKaon(C const& collision, V0 const& v0)
  {
    // initialise dynamic variables
    const float alpha = v0.alpha();
    const float qt = v0.qtarm();
    const float cosPA = v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
    const float q_K = {cutAPK0pinner * sqrt(abs(1 - ((alpha * alpha) / (cutAPK0alinner * cutAPK0alinner))))};
    /// Armenteros-Podolanski cut
    if ((qt < cutQTK0min) || (qt > cutQTK0max) || (qt < q_K)) {
      return false;
    }
    /// Cut on cosine pointing angle
    if (cosPA < cutPAV0) {
      return false;
    }
    /// Invariant mass cut
    if (abs(v0.mK0Short() - 0.497611) > invariantMassCutK0Short) {
      return false;
    }

    return true;
  };

  /// Pion daughter selection
  template <typename T>
  bool selectionPion(T const& track)
  {
    /// TOF cut daughter tracks
    if (nSigmaTOFdautrack != 0. && abs(track.tofNSigmaPi()) > nSigmaTOFdautrack) {
      return false;
    }
    /// Pion downsampling
    if (downsamplingTsalisPions > 0. && downsampleTsalisCharged(track.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion]) == 0) {
      return false;
    }
    return true;
  }

  /// Lambda selection
  template <typename C, typename V0>
  bool selectionLambda(C const& collision, V0 const& v0)
  {
    // initialise dynamic variables

    const float alpha = v0.alpha();
    const float qt = v0.qtarm();
    const float cosPA = v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ());

    const float q_L = cutAPL0 * sqrt(abs(1 - (((alpha - cutAPL1) * (alpha - cutAPL1)) / (cutAPL2 * cutAPL2))));
    /// Armenteros-Podolanski cut
    if ((alpha < cutAlphaminL) || (alpha > cutAlphamaxL) || (qt < cutQTmin) || (qt > q_L)) {
      return false;
    }
    /// Cut on cosine pointing angle
    if (cosPA < cutPAV0) {
      return false;
    }
    /// Invariant mass cut
    if (abs(v0.mLambda() - 1.115683) > invariantMassCutLambda) {
      return false;
    }
    return true;
  };

  /// Proton daughter selection
  template <typename T>
  bool selectionProton(T const& track)
  {
    /// TOF cut daughter tracks
    if (nSigmaTOFdautrack != 0. && abs(track.tofNSigmaPr()) > nSigmaTOFdautrack) {
      return false;
    }
    /// Proton downsampling
    if (downsamplingTsalisProtons > 0. && downsampleTsalisCharged(track.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton]) == 0) {
      return false;
    }
    return true;
  }

  /// Antilambda selection
  template <typename C, typename V0>
  bool selectionAntiLambda(C const& collision, V0 const& v0)
  {
    // initialise dynamic variables
    const float alpha = v0.alpha();
    const float qt = v0.qtarm();
    const float cosPA = v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ());

    const float q_L = cutAPL0 * sqrt(abs(1 - (((alpha + cutAPL1) * (alpha + cutAPL1)) / (cutAPL2 * cutAPL2))));
    /// Armenteros-Podolanski cut
    if ((alpha < cutAlphaminAL) || (alpha > cutAlphamaxAL) || (qt < cutQTmin) || (qt > q_L)) {
      return false;
    }
    /// Cut on cosine pointing angle
    if (cosPA < cutPAV0) {
      return false;
    }
    /// Invariant mass cut
    if (abs(v0.mAntiLambda() - 1.115683) > invariantMassCutLambda) {
      return false;
    }
    return true;
  };

  /// Gamma selection
  template <typename C, typename V0>
  bool selectionGamma(C const& collision, V0 const& v0)
  {
    // initialise dynamic variables
    const float alpha = v0.alpha();
    const float qt = v0.qtarm();
    const float cosPA = v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
    const float pT = v0.pt();
    float lQtMaxPtDep = gammaQtPtMultiplicator * pT;
    if (lQtMaxPtDep > gammaQtMax) {
      lQtMaxPtDep = gammaQtMax;
    }

    /// Armenteros-Podolanski cut
    if ((((qt > gammaQtMax) || (TMath::Abs(alpha) > cutAlphaG1)) || ((qt > cutQTG) || (TMath::Abs(alpha) < cutAlphaG2min) || (TMath::Abs(alpha) > cutAlphaG2max)))) {
      return false;
    }
    if ((TMath::Power(alpha / gammaAsymmetryMax, 2) + TMath::Power(qt / lQtMaxPtDep, 2)) < 1) {
      return false;
    }
    if (v0.v0radius() < gammaV0RadiusMin || v0.v0radius() > gammaV0RadiusMax) {
      return false;
    }
    /// Cut on cosine pointing angle
    if (cosPA < cutPAV0) {
      return false;
    }
    if ((TMath::Abs(v0.psipair()) < gammaPsiPairMax)) {
      return false;
    }
    return true;
  };

  /// Electron daughter selection
  template <typename T>
  bool selectionElectron(T const& track)
  {
    /// Electron downsampling
    if (downsamplingTsalisElectrons > 0. && downsampleTsalisCharged(track.pt(), downsamplingTsalisElectrons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Electron]) == 0) {
      return false;
    }
    return true;
  }

  /// Funktion to fill skimmed tables
  template <typename T, typename C, typename V0>
  void fillSkimmedV0Table(V0 const& v0, T const& track, C const& collision, const float nSigmaTPC, const float nSigmaTOF, const float dEdxExp, const o2::track::PID::ID id, double dwnSmplFactor)
  {

    const double ncl = track.tpcNClsFound();
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const int multTPC = collision.multTPC();

    const float alpha = v0.alpha();
    const float qt = v0.qtarm();
    const float cosPA = v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
    const float pT = v0.pt();

    double pseudoRndm = track.pt() * 1000. - (long)(track.pt() * 1000);
    if (pseudoRndm < dwnSmplFactor) {
      rowTPCTree(track.tpcSignal(),
                 1. / dEdxExp,
                 track.tpcInnerParam(),
                 track.tgl(),
                 track.signed1Pt(),
                 track.eta(),
                 track.phi(),
                 track.y(),
                 mass,
                 bg,
                 multTPC / 11000.,
                 std::sqrt(nClNorm / ncl),
                 id,
                 nSigmaTPC,
                 nSigmaTOF,
                 alpha,
                 qt,
                 cosPA,
                 pT);
    }
  };

  double tsalisCharged(double pt, double mass, double sqrts)
  {
    const double a = 6.81, b = 59.24;
    const double c = 0.082, d = 0.151;
    double mt = sqrt(mass * mass + pt * pt);
    double n = a + b / sqrts;
    double T = c + d / sqrts;
    double p0 = n * T;
    double result = pow((1. + mt / p0), -n);
    return result;
  };

  /// Random downsampling trigger function using Tsalis/Hagedorn spectra fit (sqrt(s) = 62.4 GeV to 13 TeV)
  /// as in https://iopscience.iop.org/article/10.1088/2399-6528/aab00f/pdf
  TRandom3* fRndm = new TRandom3(0);
  int downsampleTsalisCharged(double pt, double factor1Pt, double sqrts, double mass)
  {
    double prob = tsalisCharged(pt, mass, sqrts) * pt;
    double probNorm = tsalisCharged(1., mass, sqrts);
    int triggerMask = 0;
    if ((fRndm->Rndm() * ((prob / probNorm) * pt * pt)) < factor1Pt)
      triggerMask = 1;
    return triggerMask;
  };

  /// Event selection
  template <typename CollisionType, typename TrackType>
  bool isEventSelected(const CollisionType& collision, const TrackType& tracks)
  {
    if (applyEvSel == 1) {
      if (!collision.sel7()) {
        return false;
      }
    } else if (applyEvSel == 2) {
      if (!collision.sel8()) {
        return false;
      }
    }
    return true;
  };

  void init(o2::framework::InitContext& initContext)
  {
  }

  void process(Coll::iterator const& collision, Trks const& tracks, aod::V0Datas const& v0s)
  {
    /// Check event slection
    if (!isEventSelected(collision, tracks)) {
      return;
    }

    rowTPCTree.reserve(tracks.size());

    /// Loop over v0 candidates
    for (auto v0 : v0s) {
      auto posTrack = v0.posTrack_as<Trks>();
      auto negTrack = v0.negTrack_as<Trks>();
      /// Fill table for Kaons
      if (selectionKaon(collision, v0)) {
        if (selectionPion(posTrack)) {
          fillSkimmedV0Table(v0, posTrack, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignal()), o2::track::PID::Pion, dwnSmplFactor_Pi);
        }
        if (selectionPion(negTrack)) {
          fillSkimmedV0Table(v0, negTrack, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignal()), o2::track::PID::Pion, dwnSmplFactor_Pi);
        }
      }
      /// Fill table for Lambdas
      if (selectionLambda(collision, v0)) {
        if (selectionProton(posTrack)) {
          fillSkimmedV0Table(v0, posTrack, collision, posTrack.tpcNSigmaPr(), posTrack.tofNSigmaPr(), posTrack.tpcExpSignalPr(posTrack.tpcSignal()), o2::track::PID::Proton, dwnSmplFactor_Pr);
        }
        if (selectionPion(negTrack)) {
          fillSkimmedV0Table(v0, negTrack, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignal()), o2::track::PID::Pion, dwnSmplFactor_Pi);
        }
      }
      /// Fill table for Antilambdas
      if (selectionAntiLambda(collision, v0)) {
        if (selectionPion(posTrack)) {
          fillSkimmedV0Table(v0, posTrack, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignal()), o2::track::PID::Pion, dwnSmplFactor_Pi);
        }
        if (selectionProton(negTrack)) {
          fillSkimmedV0Table(v0, negTrack, collision, negTrack.tpcNSigmaPr(), negTrack.tofNSigmaPr(), negTrack.tpcExpSignalPr(negTrack.tpcSignal()), o2::track::PID::Proton, dwnSmplFactor_Pr);
        }
      }
      /// Fill table for Gammas
      if (selectionGamma(collision, v0)) {
        if (selectionElectron(posTrack)) {
          fillSkimmedV0Table(v0, posTrack, collision, posTrack.tpcNSigmaEl(), posTrack.tofNSigmaEl(), posTrack.tpcExpSignalEl(posTrack.tpcSignal()), o2::track::PID::Electron, dwnSmplFactor_El);
        }
        if (selectionElectron(negTrack)) {
          fillSkimmedV0Table(v0, negTrack, collision, negTrack.tpcNSigmaEl(), negTrack.tofNSigmaEl(), negTrack.tpcExpSignalEl(negTrack.tpcSignal()), o2::track::PID::Electron, dwnSmplFactor_El);
        }
      }
    } /// Loop V0 candidates
  }   /// process
};    /// struct TreeWriterTpcV0

struct TreeWriterTPCTOF {
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TrackSelection>;
  using Coll = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;

  /// Tables to be produced
  Produces<o2::aod::SkimmedTPCTOFTree> rowTPCTOFTree;

  /// Configurables
  Configurable<float> nClNorm{"nClNorm", 152., "Number of cluster normalization. Run 2: 159, Run 3 152"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  /// Proton
  Configurable<float> maxMomTPCOnlyPr{"maxMomTPCOnlyPr", 0.6, "Maximum momentum for TPC only cut proton"};
  Configurable<float> nSigmaTPCOnlyPr{"nSigmaTPCOnlyPr", 4., "number of sigma for TPC only cut proton"};
  Configurable<float> nSigmaTPC_TPCTOF_Pr{"nSigmaTPC_TPCTOF_Pr", 4., "number of sigma for TPC cut for TPC and TOF combined proton"};
  Configurable<float> nSigmaTOF_TPCTOF_Pr{"nSigmaTOF_TPCTOF_Pr", 3., "number of sigma for TOF cut for TPC and TOF combined proton"};
  Configurable<double> dwnSmplFactor_Pr{"dwnSmplFactor_Pr", 1., "downsampling factor for protons, default fraction to keep is 1."};
  /// Kaon
  Configurable<float> maxMomTPCOnlyKa{"maxMomTPCOnlyKa", 0.3, "Maximum momentum for TPC only cut kaon"};
  Configurable<float> nSigmaTPCOnlyKa{"nSigmaTPCOnlyKa", 4., "number of sigma for TPC only cut kaon"};
  Configurable<float> nSigmaTPC_TPCTOF_Ka{"nSigmaTPC_TPCTOF_Ka", 4., "number of sigma for TPC cut for TPC and TOF combined kaon"};
  Configurable<float> nSigmaTOF_TPCTOF_Ka{"nSigmaTOF_TPCTOF_Ka", 3., "number of sigma for TOF cut for TPC and TOF combined kaon"};
  Configurable<double> dwnSmplFactor_Ka{"dwnSmplFactor_Ka", 1., "downsampling factor for kaons, default fraction to keep is 1."};
  /// Pion
  Configurable<float> maxMomTPCOnlyPi{"maxMomTPCOnlyPi", 0.5, "Maximum momentum for TPC only cut pion"};
  Configurable<float> nSigmaTPCOnlyPi{"nSigmaTPCOnlyPi", 4., "number of sigma for TPC only cut pion"};
  Configurable<float> nSigmaTPC_TPCTOF_Pi{"nSigmaTPC_TPCTOF_Pi", 4., "number of sigma for TPC cut for TPC and TOF combined pion"};
  Configurable<float> nSigmaTOF_TPCTOF_Pi{"nSigmaTOF_TPCTOF_Pi", 4., "number of sigma for TOF cut for TPC and TOF combined pion"};
  Configurable<double> dwnSmplFactor_Pi{"dwnSmplFactor_Pi", 1., "downsampling factor for pions, default fraction to keep is 1."};
  /// pT dependent downsampling
  Configurable<float> sqrtSNN{"sqrt_s_NN", 0., "sqrt(s_NN), used for downsampling with the Tsallis distribution"};
  Configurable<float> downsamplingTsalisProtons{"downsamplingTsalisProtons", -1., "Downsampling factor to reduce the number of protons"};
  Configurable<float> downsamplingTsalisKaons{"downsamplingTsalisKaons", -1., "Downsampling factor to reduce the number of kaons"};
  Configurable<float> downsamplingTsalisPions{"downsamplingTsalisPions", -1., "Downsampling factor to reduce the number of pions"};

  double tsalisCharged(double pt, double mass, double sqrts)
  {
    const double a = 6.81, b = 59.24;
    const double c = 0.082, d = 0.151;
    double mt = sqrt(mass * mass + pt * pt);
    double n = a + b / sqrts;
    double T = c + d / sqrts;
    double p0 = n * T;
    double result = pow((1. + mt / p0), -n);
    return result;
  };

  /// Random downsampling trigger function using Tsalis/Hagedorn spectra fit (sqrt(s) = 62.4 GeV to 13 TeV)
  /// as in https://iopscience.iop.org/article/10.1088/2399-6528/aab00f/pdf
  TRandom3* fRndm = new TRandom3(0);
  bool downsampleTsalisCharged(double pt, float factor1Pt, double sqrts, double mass)
  {
    if (factor1Pt < 0.) {
      return true;
    }
    double prob = tsalisCharged(pt, mass, sqrts) * pt;
    double probNorm = tsalisCharged(1., mass, sqrts);
    if ((fRndm->Rndm() * ((prob / probNorm) * pt * pt)) > factor1Pt) {
      return false;
    } else {
      return true;
    }
  };

  /// Function to fill trees
  template <typename T, typename C>
  void fillSkimmedTPCTOFTable(T const& track, C const& collision, const float nSigmaTPC, const float nSigmaTOF, const float dEdxExp, const o2::track::PID::ID id, double dwnSmplFactor)
  {

    const double ncl = track.tpcNClsFound();
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const int multTPC = collision.multTPC();

    double pseudoRndm = track.pt() * 1000. - (long)(track.pt() * 1000);
    if (pseudoRndm < dwnSmplFactor) {
      rowTPCTOFTree(track.tpcSignal(),
                    1. / dEdxExp,
                    track.tpcInnerParam(),
                    track.tgl(),
                    track.signed1Pt(),
                    track.eta(),
                    track.phi(),
                    track.y(),
                    mass,
                    bg,
                    multTPC / 11000.,
                    std::sqrt(nClNorm / ncl),
                    id,
                    nSigmaTPC,
                    nSigmaTOF);
    }
  };

  /// Event selection
  template <typename CollisionType, typename TrackType>
  bool isEventSelected(const CollisionType& collision, const TrackType& tracks)
  {
    if (applyEvSel == 1) {
      if (!collision.sel7()) {
        return false;
      }
    } else if (applyEvSel == 2) {
      if (!collision.sel8()) {
        return false;
      }
    }
    return true;
  };

  /// Track selection
  template <typename CollisionType, typename TrackType>
  bool isTrackSelected(const CollisionType& collision, const TrackType& track)
  {
    if (!track.isGlobalTrack()) { // Skipping non global tracks
      return false;
    }
    if (!track.hasITS()) { // Skipping tracks without ITS
      return false;
    }
    if (!track.hasTPC()) { // Skipping tracks without TPC
      return false;
    }
    return true;
  };

  void init(o2::framework::InitContext& initContext)
  {
  }
  void process(Coll::iterator const& collision, Trks const& tracks)
  {
    /// Check event selection
    if (!isEventSelected(collision, tracks)) {
      return;
    }
    rowTPCTOFTree.reserve(tracks.size());
    for (auto const& trk : tracks) {

      /// Check track selection
      if (!isTrackSelected(collision, trk)) {
        continue;
      }
      /// Fill tree for protons
      if (trk.tpcInnerParam() < maxMomTPCOnlyPr && abs(trk.tpcNSigmaPr()) < nSigmaTPCOnlyPr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignal()), o2::track::PID::Proton, dwnSmplFactor_Pr);
      } else if (trk.tpcInnerParam() > maxMomTPCOnlyPr && abs(trk.tofNSigmaPr()) < nSigmaTOF_TPCTOF_Pr && abs(trk.tpcNSigmaPr()) < nSigmaTPC_TPCTOF_Pr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignal()), o2::track::PID::Proton, dwnSmplFactor_Pr);
      }
      /// Fill tree for kaons
      if (trk.tpcInnerParam() < maxMomTPCOnlyKa && abs(trk.tpcNSigmaKa()) < nSigmaTPCOnlyKa && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignal()), o2::track::PID::Kaon, dwnSmplFactor_Ka);
      } else if (trk.tpcInnerParam() > maxMomTPCOnlyKa && abs(trk.tofNSigmaKa()) < nSigmaTOF_TPCTOF_Ka && abs(trk.tpcNSigmaKa()) < nSigmaTPC_TPCTOF_Ka && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignal()), o2::track::PID::Kaon, dwnSmplFactor_Ka);
      }
      /// Fill tree pions
      if (trk.tpcInnerParam() < maxMomTPCOnlyPi && abs(trk.tpcNSigmaPi()) < nSigmaTPCOnlyPi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignal()), o2::track::PID::Pion, dwnSmplFactor_Pi);
      } else if (trk.tpcInnerParam() > maxMomTPCOnlyPi && abs(trk.tofNSigmaPi()) < nSigmaTOF_TPCTOF_Pi && abs(trk.tpcNSigmaPi()) < nSigmaTPC_TPCTOF_Pi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignal()), o2::track::PID::Pion, dwnSmplFactor_Pi);
      }
    } /// Loop tracks
  }   /// process
};    /// struct TreeWriterTPCTOF

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<TreeWriterTPCTOF>(cfgc)};
  workflow.push_back(adaptAnalysisTask<TreeWriterTpcV0>(cfgc));
  return workflow;
}