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
#include <algorithm>
#include <functional>
#include <memory>

#include "TF1.h"
#include "Framework/Logger.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
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

//________________________________________________________________________
bool isBHadron(int pc)
{
  std::vector<int> bPdG = {511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 20523, 20513, 20523, 515, 525, 531, 10531, 533, 10533,
                           20533, 535, 541, 10541, 543, 10543, 20543, 545, 551, 10551, 100551, 110551, 200551, 210551, 553, 10553, 20553,
                           30553, 100553, 110553, 120553, 130553, 200553, 210553, 220553, 300553, 9000533, 9010553, 555, 10555, 20555,
                           100555, 110555, 120555, 200555, 557, 100557, 5122, 5112, 5212, 5222, 5114, 5214, 5224, 5132, 5232, 5312, 5322,
                           5314, 5324, 5332, 5334, 5142, 5242, 5412, 5422, 5414, 5424, 5342, 5432, 5434, 5442, 5444, 5512, 5522, 5514, 5524,
                           5532, 5534, 5542, 5544, 5554};

  return (std::find(bPdG.begin(), bPdG.end(), std::abs(pc)) != bPdG.end());
}
//________________________________________________________________________
bool isCHadron(int pc)
{
  std::vector<int> bPdG = {411, 421, 10411, 10421, 413, 423, 10413, 10423, 20431, 20423, 415, 425, 431, 10431, 433, 10433, 20433, 435, 441,
                           10441, 100441, 443, 10443, 20443, 100443, 30443, 9000443, 9010443, 9020443, 445, 100445, 4122, 4222, 4212, 4112,
                           4224, 4214, 4114, 4232, 4132, 4322, 4312, 4324, 4314, 4332, 4334, 4412, 4422, 4414, 4424, 4432, 4434, 4444};

  return (std::find(bPdG.begin(), bPdG.end(), std::abs(pc)) != bPdG.end());
}

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
int jetTrackFromHFShower(T const& jet, U const& /*tracks*/, V const& particles, typename U::iterator& hftrack, bool searchUpToQuark)
{

  bool hasMcParticle = false;
  int origin = -1;
  for (auto& track : jet.template tracks_as<U>()) {
    if (!track.has_mcParticle()) {
      continue;
    }
    hasMcParticle = true;
    auto const& particle = track.template mcParticle_as<V>();
    origin = RecoDecay::getCharmHadronOrigin(particles, particle, searchUpToQuark);
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
int jetParticleFromHFShower(T const& jet, U const& particles, typename U::iterator& hfparticle, bool searchUpToQuark)
{

  int origin = -1;
  for (const auto& particle : jet.template tracks_as<U>()) {
    origin = RecoDecay::getCharmHadronOrigin(particles, particle, searchUpToQuark);
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
int mcdJetFromHFShower(T const& jet, U const& tracks, V const& particles, float dRMax = 0.25, bool searchUpToQuark = false)
{

  typename U::iterator hftrack;
  int origin = jetTrackFromHFShower(jet, tracks, particles, hftrack, searchUpToQuark);
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

template <typename T, typename U>
int mcpJetFromHFShower(T const& jet, U const& particles, float dRMax = 0.25, bool searchUpToQuark = false)
{

  typename U::iterator hfparticle;
  int origin = jetParticleFromHFShower(jet, particles, hfparticle, searchUpToQuark);
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
 * return the jet flavor: 0 for lf-jet, 1 for c-jet, 2 for b-jet
 *
 * @param AnyJet the jet that we need to study its flavor
 * @param AllMCParticles a vector of all the mc particles stack
 */
template <typename AnyJet, typename AllMCParticles>
int16_t getJetFlavor(AnyJet const& jet, AllMCParticles const& mcparticles)
{
  bool charmQuark = false;
  for (auto& mcpart : mcparticles) {
    int pdgcode = mcpart.pdgCode();
    if (TMath::Abs(pdgcode) == 21 || (TMath::Abs(pdgcode) >= 1 && TMath::Abs(pdgcode) <= 5)) {
      double dR = jetutilities::deltaR(jet, mcpart);

      if (dR < jet.r() / 100.f) {
        if (TMath::Abs(pdgcode) == 5) {
          return JetTaggingSpecies::beauty; // Beauty jet
        } else if (TMath::Abs(pdgcode) == 4) {
          charmQuark = true;
        }
      }
    }
  }

  if (charmQuark) {
    return JetTaggingSpecies::charm; // Charm jet
  }

  return JetTaggingSpecies::lightflavour; // Light flavor jet
}

/**
 * return the jet flavor if it finds a HF hadron inside the jet: 0 for lf-jet, 1 for c-jet, 2 for b-jet
 *
 * @param AnyJet the jet that we need to study its flavor
 * @param AllMCParticles a vector of all the mc particles stack
 */
template <typename AnyJet, typename AllMCParticles>
int16_t getJetFlavorHadron(AnyJet const& jet, AllMCParticles const& mcparticles)
{
  bool charmHadron = false;

  for (auto& mcpart : mcparticles) {
    int pdgcode = mcpart.pdgCode();
    if (isBHadron(pdgcode) || isCHadron(pdgcode)) {
      double dR = jetutilities::deltaR(jet, mcpart);

      if (dR < jet.r() / 100.f) {
        if (isBHadron(pdgcode)) {
          return JetTaggingSpecies::beauty; // Beauty jet
        } else if (isCHadron(pdgcode)) {
          charmHadron = true;
        }
      }
    }
  }

  if (charmHadron) {
    return JetTaggingSpecies::charm; // Charm jet
  }

  return JetTaggingSpecies::lightflavour; // Light flavor jet
}

/**
 * return acceptance of track about DCA xy and z due to cut for QualityTracks
 */
template <typename T>
bool trackAcceptanceWithDca(T const& track, float trackDcaXYMax, float trackDcaZMax)
{
  if (std::abs(track.dcaXY()) > trackDcaXYMax)
    return false;
  if (std::abs(track.dcaZ()) > trackDcaZMax)
    return false;
  return true;
}

/**
 * retrun acceptance of prong about chi2 and error of decay length due to cut for high quality secondary vertex
 */
template <typename T>
bool prongAcceptance(T const& prong, float prongChi2PCAMin, float prongChi2PCAMax, float prongsigmaLxyMax, float prongIPxyMin, float prongIPxyMax, bool doXYZ)
{
  if (prong.chi2PCA() < prongChi2PCAMin)
    return false;
  if (prong.chi2PCA() > prongChi2PCAMax)
    return false;
  if (!doXYZ) {
    if (prong.errorDecayLengthXY() > prongsigmaLxyMax)
      return false;
    if (std::abs(prong.impactParameterXY()) < prongIPxyMin)
      return false;
    if (std::abs(prong.impactParameterXY()) > prongIPxyMax)
      return false;
  } else {
    if (prong.errorDecayLength() > prongsigmaLxyMax)
      return false;
    // TODO
    if (std::abs(prong.impactParameterXY()) < prongIPxyMin)
      return false;
    if (std::abs(prong.impactParameterXY()) > prongIPxyMax)
      return false;
  }
  return true;
}

/**
 * return geometric sign which is calculated scalar product between jet axis with DCA (track propagated to PV )
 * positive and negative value are expected from primary vertex
 * positive value is expected from secondary vertex
 *
 * @param jet
 * @param jtrack which is needed aod::JTrackExtras
 */
template <typename T, typename U>
int getGeoSign(T const& jet, U const& jtrack)
{
  auto sign = TMath::Sign(1, jtrack.dcaX() * jet.px() + jtrack.dcaY() * jet.py() + jtrack.dcaZ() * jet.pz());
  if (sign < -1 || sign > 1)
    LOGF(info, Form("Sign is %d", sign));
  return sign;
}

/**
 * Orders the tracks associated with a jet based on signed impact parameter significance and stores them
 * in a vector in descending order.
 */
template <typename T, typename U, typename Vec = std::vector<float>>
void orderForIPJetTracks(T const& jet, U const& /*jtracks*/, float const& trackDcaXYMax, float const& trackDcaZMax, Vec& vecSignImpSig)
{
  for (auto& jtrack : jet.template tracks_as<U>()) {
    if (!trackAcceptanceWithDca(jtrack, trackDcaXYMax, trackDcaZMax))
      continue;
    auto geoSign = getGeoSign(jet, jtrack);
    auto varSignImpXYSig = geoSign * std::abs(jtrack.dcaXY()) / jtrack.sigmadcaXY();
    vecSignImpSig.push_back(varSignImpXYSig);
  }
  std::sort(vecSignImpSig.begin(), vecSignImpSig.end(), std::greater<float>());
}

/**
 * Checks if a jet is greater than the given tagging working point based on the signed impact parameter significances
 */
template <typename T, typename U>
bool isGreaterThanTaggingPoint(T const& jet, U const& jtracks, float const& trackDcaXYMax, float const& trackDcaZMax, float const& taggingPoint = 1.0, int const& cnt = 1)
{
  if (cnt == 0) {
    return true; // untagged
  }
  std::vector<float> vecSignImpSig;
  orderForIPJetTracks(jet, jtracks, trackDcaXYMax, trackDcaZMax, vecSignImpSig);
  if (vecSignImpSig.size() > static_cast<std::vector<float>::size_type>(cnt) - 1) {
    for (int i = 0; i < cnt; i++) {
      if (vecSignImpSig[i] < taggingPoint) { // tagger point set
        return false;
      }
    }
  } else {
    return false;
  }
  return true;
}

/**
 * Creates and sets the parameters of a resolution function (TF1) for constituents of jet based on signed impact parameter significnace in plane XY
 * This function is typically used to set up a resolution function for jet tagging purposes.
 *
 * @param vecParams A vector containing the parameters for the resolution function.
 * @return A unique pointer to the TF1 resolution function with the parameters set.
 */
template <typename T>
std::unique_ptr<TF1> setResolutionFunction(T const& vecParams)
{
  std::unique_ptr<TF1> fResoFunc(new TF1("fResoFunc", "gaus(0)+expo(3)+expo(5)+expo(7)", -40, 0));
  for (typename T::size_type i = 0; i < vecParams.size(); i++) {
    fResoFunc->SetParameter(i, vecParams[i]);
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
 * @param track The track for which the probability is being calculated.
 * @param minSignImpXYSig The minimum significance of the impact parameter in the XY plane, used as
 *                        the lower limit for integration of the resolution function. Defaults to -40.
 * @return The calculated probability of the track being associated with the jet, based on its
 *         impact parameter significance.
 */
template <typename T, typename U>
float getTrackProbability(T const& fResoFuncjet, U const& track, const float& minSignImpXYSig = -40)
{
  float probTrack = 0;
  auto varSignImpXYSig = TMath::Abs(track.dcaXY()) / track.sigmadcaXY();
  if (-varSignImpXYSig < minSignImpXYSig)
    varSignImpXYSig = -minSignImpXYSig - 0.01; // To avoid overflow for integral
  probTrack = fResoFuncjet->Integral(minSignImpXYSig, -varSignImpXYSig) / fResoFuncjet->Integral(minSignImpXYSig, 0);

  return probTrack;
}

/**
 * Computes the jet probability (JP) for a given jet, considering only tracks with a positive geometric
 * sign. JP is calculated using the product of individual track probabilities and the sum of logarithmic
 * terms derived from these probabilities, providing a measure for the likelihood of the jet being
 * associated with a particular flavor based on its constituent tracks' impact parameters.
 *
 * @param fResoFuncjet: The resolution function for the jet, applied to each track within the jet to
 *                     assess its probability based on the impact parameter significance.
 * @param collision: The collision event data, necessary for geometric sign calculations.
 * @param jet: The jet for which the probability is being calculated.
 * @param jtracks: Tracks in jets
 * @param cnt: ordering number of impact parameter cnt=0: untagged, cnt=1: first, cnt=2: seconde, cnt=3: third.
 * @param tagPoint: tagging working point which is selected by condiered efficiency and puriy
 * @param minSignImpXYSig: To avoid over fitting
 * @return The jet probability (JP), indicating the likelihood of the jet's association with a
 *         specific flavor. Returns -1 if the jet contains fewer than two tracks with a positive
 *         geometric sign.
 */
template <typename T, typename U, typename V>
float getJetProbability(T const& fResoFuncjet, U const& jet, V const& jtracks, float const& trackDcaXYMax, float const& trackDcaZMax, const int& cnt, const float& tagPoint = 1.0, const float& minSignImpXYSig = -10)
{
  if (!(isGreaterThanTaggingPoint(jet, jtracks, trackDcaXYMax, trackDcaZMax, tagPoint, cnt)))
    return -1;

  std::vector<float> jetTracksPt;
  float trackjetProb = 1.;

  for (auto& jtrack : jet.template tracks_as<V>()) {
    if (!trackAcceptanceWithDca(jtrack, trackDcaXYMax, trackDcaZMax))
      continue;

    float probTrack = getTrackProbability(fResoFuncjet, jtrack, minSignImpXYSig);

    auto geoSign = getGeoSign(jet, jtrack);
    if (geoSign > 0) { // only take positive sign track for JP calculation
      trackjetProb *= probTrack;
      jetTracksPt.push_back(jtrack.pt());
    }
  }

  float JP = -1.;
  if (jetTracksPt.size() < 2)
    return -1;

  float sumjetProb = 0.;
  for (std::vector<float>::size_type i = 0; i < jetTracksPt.size(); i++) {
    sumjetProb += (TMath::Power(-1 * TMath::Log(trackjetProb), static_cast<int>(i)) / TMath::Factorial(i));
  }

  JP = trackjetProb * sumjetProb;
  return JP;
}

// For secaondy vertex method utilites
class bjetCandSV
{
 public:
  bjetCandSV() = default;

  bjetCandSV(float xpv, float ypv, float zpv, float xsv, float ysv, float zsv,
             float pxVal, float pyVal, float pzVal, float eVal, float mVal, float chi2Val,
             float errDecayLength, float errDecayLengthXY,
             float rSecVertex, float ptVal, float pVal,
             std::array<float, 3> pVec, float etaVal, float phiVal,
             float yVal, float decayLen, float decayLenXY,
             float decayLenNorm, float decayLenXYNorm,
             float cpaVal, float impParXY)
    : m_xPVertex(xpv), m_yPVertex(ypv), m_zPVertex(zpv), m_xSecondaryVertex(xsv), m_ySecondaryVertex(ysv), m_zSecondaryVertex(zsv), m_px(pxVal), m_py(pyVal), m_pz(pzVal), m_e(eVal), m_m(mVal), m_chi2PCA(chi2Val), m_errorDecayLength(errDecayLength), m_errorDecayLengthXY(errDecayLengthXY), m_rSecondaryVertex(rSecVertex), m_pt(ptVal), m_p(pVal), m_pVector(pVec), m_eta(etaVal), m_phi(phiVal), m_y(yVal), m_decayLength(decayLen), m_decayLengthXY(decayLenXY), m_decayLengthNormalised(decayLenNorm), m_decayLengthXYNormalised(decayLenXYNorm), m_cpa(cpaVal), m_impactParameterXY(impParXY)
  {
  }

  float xPVertex() const { return m_xPVertex; }
  float yPVertex() const { return m_yPVertex; }
  float zPVertex() const { return m_zPVertex; }

  float xSecondaryVertex() const { return m_xSecondaryVertex; }
  float ySecondaryVertex() const { return m_ySecondaryVertex; }
  float zSecondaryVertex() const { return m_zSecondaryVertex; }

  float px() const { return m_px; }
  float py() const { return m_py; }
  float pz() const { return m_pz; }
  float e() const { return m_e; }
  float m() const { return m_m; }
  float chi2PCA() const { return m_chi2PCA; }

  float errorDecayLength() const { return m_errorDecayLength; }
  float errorDecayLengthXY() const { return m_errorDecayLengthXY; }

  float rSecondaryVertex() const { return m_rSecondaryVertex; }
  float pt() const { return m_pt; }
  float p() const { return m_p; }

  std::array<float, 3> pVector() const { return m_pVector; }

  float eta() const { return m_eta; }
  float phi() const { return m_phi; }
  float y() const { return m_y; }

  float decayLength() const { return m_decayLength; }
  float decayLengthXY() const { return m_decayLengthXY; }
  float decayLengthNormalised() const { return m_decayLengthNormalised; }
  float decayLengthXYNormalised() const { return m_decayLengthXYNormalised; }

  float cpa() const { return m_cpa; }
  float impactParameterXY() const { return m_impactParameterXY; }

 private:
  float m_xPVertex, m_yPVertex, m_zPVertex;
  float m_xSecondaryVertex, m_ySecondaryVertex, m_zSecondaryVertex;
  float m_px, m_py, m_pz, m_e, m_m, m_chi2PCA;
  float m_errorDecayLength, m_errorDecayLengthXY;
  float m_rSecondaryVertex, m_pt, m_p;
  std::array<float, 3> m_pVector;
  float m_eta, m_phi, m_y;
  float m_decayLength, m_decayLengthXY, m_decayLengthNormalised, m_decayLengthXYNormalised;
  float m_cpa, m_impactParameterXY;
};

template <typename ProngType, typename JetType>
bjetCandSV jetFromProngMaxDecayLength(const JetType& jet, float const& prongChi2PCAMin, float const& prongChi2PCAMax, float const& prongsigmaLxyMax, float const& prongIPxyMin, float const& prongIPxyMax, const bool& doXYZ = false)
{
  float xPVertex = 0.0f;
  float yPVertex = 0.0f;
  float zPVertex = 0.0f;
  float xSecondaryVertex = 0.0f;
  float ySecondaryVertex = 0.0f;
  float zSecondaryVertex = 0.0f;
  float px = 0.0f;
  float py = 0.0f;
  float pz = 0.0f;
  float e = 0.0f;
  float m = 0.0f;
  float chi2PCA = 0.0f;
  float errorDecayLength = 0.0f;
  float errorDecayLengthXY = 0.0f;

  float rSecondaryVertex = 0.0f;
  float pt = 0.0f;
  float p = 0.0f;
  std::array<float, 3> pVector = {0.0f, 0.0f, 0.0f};
  float eta = 0.0f;
  float phi = 0.0f;
  float y = 0.0f;
  float decayLength = 0.0f;
  float decayLengthXY = 0.0f;
  float decayLengthNormalised = 0.0f;
  float decayLengthXYNormalised = 0.0f;
  float cpa = 0.0f;
  float impactParameterXY = 0.0f;

  float maxSxy = -1.0f;

  for (const auto& prong : jet.template secondaryVertices_as<ProngType>()) {
    float Sxy = -1.;
    if (!doXYZ) {
      Sxy = prong.decayLengthXY() / prong.errorDecayLengthXY();
    } else {
      Sxy = prong.decayLength() / prong.errorDecayLength();
    }

    if (maxSxy < Sxy) {
      maxSxy = Sxy;

      xPVertex = prong.xPVertex();
      yPVertex = prong.yPVertex();
      zPVertex = prong.zPVertex();
      xSecondaryVertex = prong.xSecondaryVertex();
      ySecondaryVertex = prong.ySecondaryVertex();
      zSecondaryVertex = prong.zSecondaryVertex();
      px = prong.px();
      py = prong.py();
      pz = prong.pz();
      e = prong.e();
      m = prong.m();
      chi2PCA = prong.chi2PCA();
      errorDecayLength = prong.errorDecayLength();
      errorDecayLengthXY = prong.errorDecayLengthXY();
      rSecondaryVertex = prong.rSecondaryVertex();
      pt = prong.pt();
      p = prong.p();
      pVector = prong.pVector();
      eta = prong.eta();
      phi = prong.phi();
      y = prong.y();
      decayLength = prong.decayLength();
      decayLengthXY = prong.decayLengthXY();
      decayLengthNormalised = prong.decayLengthNormalised();
      decayLengthXYNormalised = prong.decayLengthXYNormalised();
      cpa = prong.cpa();
      impactParameterXY = prong.impactParameterXY();
    }
  }

  return bjetCandSV(
    xPVertex, yPVertex, zPVertex,
    xSecondaryVertex, ySecondaryVertex, zSecondaryVertex,
    px, py, pz, e, m, chi2PCA,
    errorDecayLength, errorDecayLengthXY,
    rSecondaryVertex, pt, p,
    pVector, eta, phi,
    y, decayLength, decayLengthXY,
    decayLengthNormalised, decayLengthXYNormalised,
    cpa, impactParameterXY);
}

template <typename T, typename U>
bool isTaggedJetSV(T const jet, U const& /*prongs*/, float const& prongChi2PCAMin, float const& prongChi2PCAMax, float const& prongsigmaLxyMax, float const& prongIPxyMin, float const& prongIPxyMax, float const& doXYZ = false, float const& tagPointForSV = 15.)
{
  auto bjetCand = jetFromProngMaxDecayLength<U>(jet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, doXYZ);
  if (!prongAcceptance(bjetCand, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, doXYZ))
    return false;
  if (!doXYZ) {
    auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
    if (maxSxy < tagPointForSV)
      return false;
  } else {
    auto maxSxyz = bjetCand.decayLength() / bjetCand.errorDecayLength();
    if (maxSxyz < tagPointForSV)
      return false;
  }
  return true;
}

}; // namespace jettaggingutilities

#endif // PWGJE_CORE_JETTAGGINGUTILITIES_H_
