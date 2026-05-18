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
/// \file OTFParticle.h
/// \brief Basic class to hold information regarding a mc particle to be used in fast simulation
/// \author Jesper Karlsson Gumprecht <jesper.gumprecht@cern.ch>
///

#ifndef ALICE3_CORE_OTFPARTICLE_H_
#define ALICE3_CORE_OTFPARTICLE_H_

#include <CommonConstants/MathConstants.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <span>

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
    mVt = particle.vt();
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
  void setProductionTime(const float vt) { mVt = vt; }
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
    return Flags; // todo
  }
  int statusCode() const
  {
    static constexpr int StatusCode = 1;
    return StatusCode; // todo
  }
  float vx() const { return mVx; }
  float vy() const { return mVy; }
  float vz() const { return mVz; }
  float vt() const { return mVt; }
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
  int getMotherIndexStart() const { return mIndicesMother[0]; }
  int getMotherIndexStop() const { return mIndicesMother[1]; }
  int getDaughterIndexStart() const { return mIndicesDaughter[0]; }
  int getDaughterIndexStop() const { return mIndicesDaughter[1]; }
  std::array<int, 2> getMothers() const { return mIndicesMother; }
  std::array<int, 2> getDaughters() const { return mIndicesDaughter; }
  std::span<const int> getMotherSpan() const { return hasMothers() ? std::span<const int>(mIndicesMother.data(), 2) : std::span<const int>(); }

  // Checks
  bool hasDaughters() const { return (mIndicesDaughter[0] > 0); }
  bool hasMothers() const { return (mIndicesMother[0] > 0); }
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
  float mVx{}, mVy{}, mVz{}, mVt{};
  float mPx{}, mPy{}, mPz{}, mE{};
  bool mIsAlive{}, mIsFromMcParticles{false};
  bool mIsPrimary{};
  std::array<int, 2> mIndicesMother{-1, -1}, mIndicesDaughter{-1, -1};
};

} // namespace o2::upgrade
#endif // ALICE3_CORE_OTFPARTICLE_H_
