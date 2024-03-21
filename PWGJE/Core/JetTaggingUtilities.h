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

/// \file JetTaggingUtilities.h
/// \brief Jet tagging related utilities
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Hanseo Park <hanseo.park@cern.ch>

#ifndef PWGJE_CORE_JETTAGGINGUTILITIES_H_
#define PWGJE_CORE_JETTAGGINGUTILITIES_H_

#include <cmath>
#include <limits>
#include <numeric>
#include <tuple>
#include <vector>

#include "TF1.h"
#include "Framework/Logger.h"
#include "Common/Core/RecoDecay.h"
#include "PWGJE/Core/JetUtilities.h"

enum JetTaggingSpecies {
  none = 0,
  charm = 1,
  beauty = 2,
  lightflavour = 3,
  lightquark = 4,
  gluon = 5
};

namespace jettaggingutilities
{
const int cmTomum = 10000; // using cm -> #mum for impact parameter (dca)

/**
 * returns the globalIndex of the earliest mother of a particle in the shower. returns -1 if a suitable mother is not found
 *
 * @param particle MCParticle whose mother is to be found
 */
template <typename T>
int getOriginalMotherIndex(const typename T::iterator& particle)
{

  // if (!particle) {
  //   return -1.0;
  // }
  auto mother = particle;

  while (mother.has_mothers()) {

    mother = mother.template mothers_first_as<T>();

    int motherStatusCode = std::abs(mother.getGenStatusCode());

    if (motherStatusCode == 23 || motherStatusCode == 33 || motherStatusCode == 43 || motherStatusCode == 63) {
      return mother.globalIndex();
    }
  }
  return -1.0;
}

/**
 * returns the globalIndex of the earliest HF mother of a particle in the shower. returns -1 if a suitable mother is not found. Should be used only on already identified HF particles
 *
 * @param hfparticle MCParticle whose mother is to be found
 */

template <typename T>
int getOriginalHFMotherIndex(const typename T::iterator& hfparticle)
{

  // if (!hfparticle) {
  //   return -1.0;
  // }
  auto mother = hfparticle;

  while (mother.has_mothers()) {

    mother = mother.template mothers_first_as<T>();

    int motherStatusCode = std::abs(mother.getGenStatusCode());

    if (motherStatusCode == 23 || motherStatusCode == 33 || motherStatusCode == 43 || motherStatusCode == 63 || (motherStatusCode == 51 && mother.template mothers_first_as<T>().pdgCode() == 21)) {
      return mother.globalIndex();
    }
  }
  return -1.0;
}

/**
 * checks if atrack in a reco level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The first track originating from an HF shower can be extracted by reference
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param hftrack track passed as reference which is then replaced by the first track that originated from an HF shower
 */
template <typename T, typename U, typename V>
int jetTrackFromHFShower(T const& jet, U const& tracks, V const& particles, typename U::iterator& hftrack)
{

  bool hasMcParticle = false;
  int origin = -1;
  for (auto& track : jet.template tracks_as<U>()) {
    if (!track.has_mcParticle()) {
      continue;
    }
    hasMcParticle = true;
    auto const& particle = track.template mcParticle_as<V>();
    origin = RecoDecay::getCharmHadronOrigin(particles, particle, true);
    if (origin == 1 || origin == 2) { // 1=charm , 2=beauty
      hftrack = track;
      if (origin == 1) {
        return JetTaggingSpecies::charm;
      }
      if (origin == 2) {
        return JetTaggingSpecies::beauty;
      }
    }
  }
  if (hasMcParticle) {
    return JetTaggingSpecies::lightflavour;
  } else {
    return JetTaggingSpecies::none;
  }
}

/**
 * checks if a particle in a generator level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The first particle originating from an HF shower can be extracted by reference
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param hfparticle particle passed as reference which is then replaced by the first track that originated from an HF shower
 */
template <typename T, typename U>
int jetParticleFromHFShower(T const& jet, U const& particles, typename U::iterator& hfparticle)
{

  int origin = -1;
  for (auto& particle : jet.template tracks_as<U>()) {
    origin = RecoDecay::getCharmHadronOrigin(particles, particle, true);
    if (origin == 1 || origin == 2) { // 1=charm , 2=beauty
      hfparticle = particle;
      if (origin == 1) {
        return JetTaggingSpecies::charm;
      }
      if (origin == 2) {
        return JetTaggingSpecies::beauty;
      }
    }
  }
  return JetTaggingSpecies::lightflavour;
}

/**
 * returns if a reco level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The requirement is that the jet contains a particle from an HF shower and that the original HF quark is within dRMax of the jet axis in eta-phi
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param dRMax maximum distance in eta-phi of initiating heavy-flavour quark from the jet axis
 */

template <typename T, typename U, typename V>
int mcdJetFromHFShower(T const& jet, U const& tracks, V const& particles, float dRMax = 0.25)
{

  typename U::iterator hftrack;
  int origin = jetTrackFromHFShower(jet, tracks, particles, hftrack);
  if (origin == JetTaggingSpecies::charm || origin == JetTaggingSpecies::beauty) {
    if (!hftrack.has_mcParticle()) {
      return JetTaggingSpecies::none;
    }
    auto const& hfparticle = hftrack.template mcParticle_as<V>();

    int originalHFMotherIndex = getOriginalHFMotherIndex<V>(hfparticle);
    if (originalHFMotherIndex > -1.0) {

      if (jetutilities::deltaR(jet, particles.iteratorAt(originalHFMotherIndex)) < dRMax) {

        return origin;

      } else {
        return JetTaggingSpecies::none;
      }

    } else {
      return JetTaggingSpecies::none;
    }

  } else {

    return JetTaggingSpecies::lightflavour;
  }
}

/**
 * checks if a generator level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The requirement is that the jet contains a particle from an HF shower and that the original HF quark is within dRMax of the jet axis in eta-phi
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param dRMax maximum distance in eta-phi of initiating heavy-flavour quark from the jet axis
 */

template <typename T, typename U, typename V>
int mcpJetFromHFShower(T const& jet, U const& particles, float dRMax = 0.25)
{

  typename U::iterator hfparticle;
  int origin = jetParticleFromHFShower(jet, particles, hfparticle);
  if (origin == JetTaggingSpecies::charm || origin == JetTaggingSpecies::beauty) {

    int originalHFMotherIndex = getOriginalHFMotherIndex<U>(hfparticle);
    if (originalHFMotherIndex > -1.0) {

      if (jetutilities::deltaR(jet, particles.iteratorAt(originalHFMotherIndex)) < dRMax) {

        return origin;

      } else {
        return JetTaggingSpecies::none;
      }

    } else {
      return JetTaggingSpecies::none;
    }

  } else {

    return JetTaggingSpecies::lightflavour;
  }
}

/**
 * returns the pdg code of the original scattered parton closest to the jet axis, with the restriction that the parton and jet axis must be within dRMax in eta-phi
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param dRMax maximum distance in eta-phi of initiating heavy-flavour quark from the jet axis
 */

template <typename T, typename U>
int jetOrigin(T const& jet, U const& particles, float dRMax = 0.25)
{
  bool firstPartonFound = false;
  typename U::iterator parton1;
  typename U::iterator parton2;
  for (auto& particle : particles) {
    if (std::abs(particle.getGenStatusCode() == 23)) {
      if (!firstPartonFound) {
        parton1 = particle;
        firstPartonFound = true;
      } else {
        parton2 = particle;
      }
    }
  }

  float dR1 = jetutilities::deltaR(jet, parton1);
  float dR2 = jetutilities::deltaR(jet, parton2);

  if (dR1 <= dR2 && dR1 < dRMax) {

    return parton1.pdgCode();
  }
  if (dR2 <= dR1 && dR2 < dRMax) {

    return parton2.pdgCode();
  }

  return 0;
}

/**
 * return geometric sign which is calculated scalar product between jet axis with DCA (track propagated to PV )
 * positive and negative value are expected from primary vertex
 * positive value is expected from secondary vertex
 *
 * @param collision which is needed external table of collision due to postion X and Y
 * @param jet
 * @param track which is needed each DCA_X and Y which is measured in jettaggerhfExtension.cxx
 */
template <typename T, typename U, typename V>
int getGeoSign(T const& collision, U const& jet, V const& track)
{
  auto trackPar = getTrackPar(track);
  auto xyz = trackPar.getXYZGlo();
  auto dcaX = xyz.X() - collision.posX();
  auto dcaY = xyz.Y() - collision.posY();
  auto sign = TMath::Sign(1, dcaX * jet.px() + dcaY * jet.py() + track.dcaZ() * jet.pz());
  if (sign < -1 || sign > 1)
    LOGF(info, Form("Sign is %d", sign));
  return sign;
}

void calculateDcaXYZ(float& dcaXYZ, float& sigmaDcaXYZ2, float dcaXY, float dcaZ, float cYY, float cZY, float cZZ, float sigmaDcaXY2, float sigmaDcaZ2)
{
  dcaXYZ = std::sqrt(dcaXY * dcaXY + dcaZ * dcaZ);
  float dFdxy = 2 * dcaXY / dcaXYZ;
  float dFdz = 2 * dcaZ / dcaXYZ;
  sigmaDcaXYZ2 = std::abs(cYY * dFdxy * dFdxy + cZZ * dFdz * dFdz + 2 * cZY * dFdxy * dFdz);
}

/**
 * Generates and configures a resolution function for the impact parameter distribution of jets,
 * tailored to the specified jet flavour. The resolution function is a composite of exponential
 * and Gaussian components, designed to model the distribution's characteristics accurately.
 *
 * @param jetflavour An identifier for the jet flavour (e.g., 1 for c-jet, 2 for b-jet, 3 for light-flavour jet),
 *                   used to select the appropriate parameter set for the resolution function.
 */
// TODO: The fitting function obtained with fewer events locally. We will change to parameters with higher statistics in the future
template <typename T = float>
TF1* getResolutionFunction(T const& jetflavour)
{
  TF1* fResoFunc = nullptr;
  fResoFunc = new TF1("fResoFunc", "expo(0)+expo(2)+expo(4)+gaus(6)", -40, 0);
  switch (static_cast<int>(jetflavour)) {
    case 1: // c-jet
      fResoFunc->SetParameters(3.56397, 0.0977181, 4.34329, 0.497469, 1.29482, 0.494429, 7350.86, 0.27023, 1.0822);
      break;
    case 2: // b-jet
      fResoFunc->SetParameters(5.28562, 0.754341, 2.67886, 0.0556572, 4.26838, 173.447, 1484.69, 0.0359118, -0.923933);
      break;
    case 3: // lf-jet (light-flavour jet)
      fResoFunc->SetParameters(7.69109, 1.02418, 1.91596, 0.0589323, -5.80351, 1.06489, 6282.07, -0.0558066, 0.888488);
      break;
    default:
      // Handle unexpected jet flavour, possibly with a warning or default function
      fResoFunc->SetParameters(8.07756, 0.89488, 3.73556, 0.0628505, -0.0302361, 0.865316, 14170.4, 0.0148564, 0.915118);
      break;
  }

  return fResoFunc;
}

/**
 * Calculates the probability of a given track being associated with a jet, based on the geometric
 * sign and the resolution function of the jet's impact parameter significance. This probability
 * helps in distinguishing between tracks likely originating from the primary vertex and those from
 * secondary vertices, aiding in jet flavor tagging.
 *
 * @param fResoFuncjet The resolution function for the jet, used to model the distribution of impact
 *                     parameter significances for tracks associated with the jet.
 * @param collision The collision event data, providing context for calculating geometric signs.
 * @param jet The specific jet being analyzed.
 * @param track The track for which the probability is being calculated.
 * @param minSignImpXYSig The minimum significance of the impact parameter in the XY plane, used as
 *                        the lower limit for integration of the resolution function. Defaults to -10.
 * @return The calculated probability of the track being associated with the jet, based on its
 *         impact parameter significance.
 */
template <typename T, typename U, typename V, typename W>
float getTrackProbability(T const& fResoFuncjet, U const& collision, V const& jet, W const& track, const float& minSignImpXYSig = -10)
{
  float probTrack = 0.;
  auto varSignImpXYSig = getGeoSign(collision, jet, track) * TMath::Abs(track.dcaXY()) / TMath::Sqrt(track.sigmaDcaXY2());
  probTrack = fResoFuncjet->Integral(minSignImpXYSig, -1 * TMath::Abs(varSignImpXYSig)) / fResoFuncjet->Integral(minSignImpXYSig, 0);

  return probTrack;
}

/**
 * Computes the jet probability (JP) for a given jet, considering only tracks with a positive geometric
 * sign. JP is calculated using the product of individual track probabilities and the sum of logarithmic
 * terms derived from these probabilities, providing a measure for the likelihood of the jet being
 * associated with a particular flavor based on its constituent tracks' impact parameters.
 *
 * @param fResoFuncjet The resolution function for the jet, applied to each track within the jet to
 *                     assess its probability based on the impact parameter significance.
 * @param collision The collision event data, necessary for geometric sign calculations.
 * @param jet The jet for which the probability is being calculated.
 * @param tracks The collection of tracks associated with the jet.
 * @return The jet probability (JP), indicating the likelihood of the jet's association with a
 *         specific flavor. Returns -1 if the jet contains fewer than two tracks with a positive
 *         geometric sign.
 */
template <typename T, typename U, typename V, typename W, typename X>
float getJetProbability(T const& fResoFuncjet, U const& collision, V const& jet, W const& jtracks, X const& tracks)
{
  std::vector<float> jetTracksPt;
  float trackjetProb = 1.;

  for (auto& jtrack : jet.template tracks_as<W>()) {
    auto track = jtrack.template track_as<X>();

    float probTrack = getTrackProbability(fResoFuncjet, collision, jet, track);

    auto geoSign = getGeoSign(collision, jet, track);
    if (geoSign > 0) { // only take positive sign track for JP calculation
      trackjetProb *= TMath::Abs(probTrack);
      jetTracksPt.push_back(track.pt());
    }
  }

  float JP = -1.;
  if (jetTracksPt.size() < 2)
    return -1;

  float sumjetProb = 0.;
  for (int i = 0; i < jetTracksPt.size(); i++) {
    sumjetProb += (TMath::Power(-1 * TMath::Log(trackjetProb), i) / TMath::Factorial(i));
  }

  JP = trackjetProb * sumjetProb;
  return JP;
}

}; // namespace jettaggingutilities

#endif // PWGJE_CORE_JETTAGGINGUTILITIES_H_
