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
/// \file   TrackUtilities.h
/// \author Nicolò Jacazio, Universita del Piemonte Orientale (IT)
/// \brief  Set of utilities for the ALICE3 track handling
/// \since  May 21, 2025
///

#ifndef ALICE3_CORE_TRACKUTILITIES_H_
#define ALICE3_CORE_TRACKUTILITIES_H_

#include <ReconstructionDataFormats/Track.h>

#include <TLorentzVector.h>

#include <array>
#include <cmath>
#include <span>
#include <vector>

namespace o2::upgrade
{

class OTFParticle
{
 public:
  OTFParticle() = default;

  template <typename TParticle>
  explicit OTFParticle(const TParticle& particle)
  {
    mPdgCode = particle.pdgCode();
    mGlobalIndex = particle.globalIndex();
    mCollisionId = particle.mcCollisionId();
    mPx = particle.px();
    mPy = particle.py();
    mPz = particle.pz();
    mE = particle.e();
    mVx = particle.vx();
    mVy = particle.vy();
    mVz = particle.vz();
    mIsFromMcParticles = true;
    if (particle.has_mothers()) {
      mIndicesMother = {particle.mothersIds().front(), particle.mothersIds().back()};
    }
  }

  // Setters
  void setIsAlive(const bool isAlive) { mIsAlive = isAlive; }
  void setIsPrimary(const bool isPrimary) { mIsPrimary = isPrimary; }
  void setCollisionId(const int collisionId) { mCollisionId = collisionId; }
  void setPDG(const int pdg) { mPdgCode = pdg; }
  void setIndicesMother(const int start, const int stop) { mIndicesMother = {start, stop}; }
  void setIndicesDaughter(const int start, const int stop) { mIndicesDaughter = {start, stop}; }
  void setVxVyVz(const float vx, const float vy, const float vz)
  {
    mVx = vx;
    mVy = vy;
    mVz = vz;
  }
  void setPxPyPzE(const float px, const float py, const float pz, const float e)
  {
    mPx = px;
    mPy = py;
    mPz = pz;
    mE = e;
  }

  // Getters
  int pdgCode() const { return mPdgCode; }
  int globalIndex() const { return mGlobalIndex; }
  int collisionId() const { return mCollisionId; }
  bool isAlive() const { return mIsAlive; }
  bool isPrimary() const { return mIsPrimary; }
  bool isFromMcParticles() const { return mIsFromMcParticles; }
  float weight() const
  {
    static constexpr float Weight = 1.f;
    return Weight;
  }
  uint8_t flags() const
  {
    static constexpr uint8_t Flags = 1;
    return Flags;
  }
  float vt() const
  {
    static constexpr float Vt = 1.f;
    return Vt;
  }
  int statusCode() const
  {
    static constexpr int StatusCode = 1;
    return StatusCode;
  }
  float vx() const { return mVx; }
  float vy() const { return mVy; }
  float vz() const { return mVz; }
  float px() const { return mPx; }
  float py() const { return mPy; }
  float pz() const { return mPz; }
  float e() const { return mE; }
  float radius() const { return std::hypot(mVx, mVy); }
  float r() const { return radius(); }
  float pt() const { return std::hypot(mPx, mPy); }
  float p() const { return std::hypot(mPx, mPy, mPz); }
  float phi() const { return o2::constants::math::PI + std::atan2(-1.0f * py(), -1.0f * px()); }
  float eta() const
  {
    // As https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1943
    static constexpr float Tolerance = 1e-7f;
    if ((p() - mPz) < Tolerance) {
      return (mPz < 0.0f) ? -100.0f : 100.0f;
    } else {
      return 0.5f * std::log((p() + mPz) / (p() - mPz));
    }
  }

  float y() const
  {
    // As https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1922
    static constexpr float Tolerance = 1e-7f;
    if ((e() - mPz) < Tolerance) {
      return (mPz < 0.0f) ? -100.0f : 100.0f;
    } else {
      return 0.5f * std::log((mE + mPz) / (mE - mPz));
    }
  }

  bool hasDaughters() const { return (mIndicesDaughter[0] > 0); }
  bool hasMothers() const { return (mIndicesMother[0] > 0); }
  int getMotherIndexStart() const { return mIndicesMother[0]; }
  int getMotherIndexStop() const { return mIndicesMother[1]; }
  int getDaughterIndexStart() const { return mIndicesDaughter[0]; }
  int getDaughterIndexStop() const { return mIndicesDaughter[1]; }
  std::array<int, 2> getMothers() const { return mIndicesMother; }
  std::array<int, 2> getDaughters() const { return mIndicesDaughter; }
  std::span<const int> getMotherSpan() const { return hasMothers() ? std::span<const int>(mIndicesMother.data(), 2) : std::span<const int>(); }

  bool hasNaN() const
  {
    return std::isnan(mPx) || std::isnan(mPy) || std::isnan(mPz) || std::isnan(mE) ||
           std::isnan(mVx) || std::isnan(mVy) || std::isnan(mVz);
  }

  bool hasIndex() const
  {
    return (mGlobalIndex != -1);
  }

 private:
  int mPdgCode{}, mGlobalIndex{-1};
  int mCollisionId{};
  float mE{};
  float mVx{}, mVy{}, mVz{};
  float mPx{}, mPy{}, mPz{};
  bool mIsAlive{}, mIsFromMcParticles{false};
  bool mIsPrimary{};

  std::array<int, 2> mIndicesMother{-1, -1}, mIndicesDaughter{-1, -1};
};

/// Function to convert a TLorentzVector into a perfect Track
/// \param charge particle charge (integer)
/// \param particle the particle to convert (TLorentzVector)
/// \param productionVertex where the particle was produced
/// \param o2track the address of the resulting TrackParCov
void convertTLorentzVectorToO2Track(const int charge,
                                    const TLorentzVector particle,
                                    const std::vector<double> productionVertex,
                                    o2::track::TrackParCov& o2track);

/// Function to convert a TLorentzVector into a perfect Track
/// \param pdgCode particle pdg
/// \param particle the particle to convert (TLorentzVector)
/// \param productionVertex where the particle was produced
/// \param o2track the address of the resulting TrackParCov
/// \param pdg the pdg service
template <typename PdgService>
void convertTLorentzVectorToO2Track(int pdgCode,
                                    TLorentzVector particle,
                                    std::vector<double> productionVertex,
                                    o2::track::TrackParCov& o2track,
                                    const PdgService& pdg)
{
  const auto pdgInfo = pdg->GetParticle(pdgCode);
  int charge = 0;
  if (pdgInfo != nullptr) {
    charge = pdgInfo->Charge() / 3;
  }
  convertTLorentzVectorToO2Track(charge, particle, productionVertex, o2track);
}

/// Function to convert a OTFParticle into a perfect Track
/// \param particle the particle to convert (OTFParticle)
/// \param o2track the address of the resulting TrackParCov
/// \param pdg the pdg service
template <typename PdgService>
void convertOTFParticleToO2Track(const OTFParticle& particle,
                                 o2::track::TrackParCov& o2track,
                                 const PdgService& pdg)
{
  static TLorentzVector tlv;
  tlv.SetPxPyPzE(particle.px(), particle.py(), particle.pz(), particle.e());
  convertTLorentzVectorToO2Track(particle.pdgCode(), tlv, {particle.vx(), particle.vy(), particle.vz()}, o2track, pdg);
}

/// Function to convert a McParticle into a perfect Track
/// \param particle the particle to convert (mcParticle)
/// \param o2track the address of the resulting TrackParCov
/// \param pdg the pdg service
template <typename McParticleType, typename PdgService>
void convertMCParticleToO2Track(McParticleType& particle,
                                o2::track::TrackParCov& o2track,
                                const PdgService& pdg)
{
  static TLorentzVector tlv;
  tlv.SetPxPyPzE(particle.px(), particle.py(), particle.pz(), particle.e());
  convertTLorentzVectorToO2Track(particle.pdgCode(), tlv, {particle.vx(), particle.vy(), particle.vz()}, o2track, pdg);
}

/// Function to convert a McParticle into a perfect Track
/// \param particle the particle to convert (mcParticle)
/// \param o2track the address of the resulting TrackParCov
/// \param pdg the pdg service
template <typename McParticleType, typename PdgService>
o2::track::TrackParCov convertMCParticleToO2Track(McParticleType& particle,
                                                  const PdgService& pdg)
{
  o2::track::TrackParCov o2track;
  convertMCParticleToO2Track(particle, o2track, pdg);
  return o2track;
}

} // namespace o2::upgrade

#endif // ALICE3_CORE_TRACKUTILITIES_H_
