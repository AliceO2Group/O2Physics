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
#include <bitset>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>

namespace o2::upgrade
{

enum class DecayerBits { ProducedByDecayer = 0,
                         IsPrimary,
                         IsAlive };

class OTFParticle
{
 public:
  OTFParticle() = default;

  template <typename TParticle>
  explicit OTFParticle(const TParticle& particle) : mPdgCode(particle.pdgCode()),
                                                    mGlobalIndex(particle.globalIndex()),
                                                    mCollisionId(particle.mcCollisionId()),
                                                    mVx(particle.vx()),
                                                    mVy(particle.vy()),
                                                    mVz(particle.vz()),
                                                    mVt(particle.vt()),
                                                    mPx(particle.px()),
                                                    mPy(particle.py()),
                                                    mPz(particle.pz()),
                                                    mE(particle.e()),
                                                    mStatusCode(particle.statusCode()),
                                                    mFlag(particle.flags())
  {
    setBitOff(DecayerBits::ProducedByDecayer);
    if (particle.has_mothers()) {
      mIndicesMother = {particle.mothersIds().front(), particle.mothersIds().back()};
    }
    if (particle.has_daughters()) {
      mIndicesDaughter = {particle.daughtersIds().front(), particle.daughtersIds().back()};
    }
    if constexpr (requires { particle.decayerBits(); }) {
      mBits = particle.decayerBits();
    } else {
      // If we are here, we created particle in the standard workflow -- without secondaries
      // Then we should set all particles as physical primaries accordingly
      setBitOn(DecayerBits::IsPrimary);
    }
  }

  // Setters
  void setCollisionId(const int collisionId) { mCollisionId = collisionId; }
  void setPDG(const int pdg) { mPdgCode = pdg; }
  void setIndicesMother(const int start, const int stop) { mIndicesMother = {start, stop}; }
  void setIndicesDaughter(const int start, const int stop) { mIndicesDaughter = {start, stop}; }
  void setProductionTime(const float vt) { mVt = vt; }
  void setFlags(uint8_t flag) { mFlag = flag; }
  void setDecayRadius(const float decayRadius) { mDecayRadius = decayRadius; }
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
  void setIndexOffset(const std::size_t offset)
  {
    static constexpr int NotFound = -1;
    mIndicesMother[0] = (mIndicesMother[0] >= 0) ? mIndicesMother[0] + static_cast<int>(offset) : NotFound;
    mIndicesMother[1] = (mIndicesMother[1] >= 0) ? mIndicesMother[1] + static_cast<int>(offset) : NotFound;
    mIndicesDaughter[0] = (mIndicesDaughter[0] >= 0) ? mIndicesDaughter[0] + static_cast<int>(offset) : NotFound;
    mIndicesDaughter[1] = (mIndicesDaughter[1] >= 0) ? mIndicesDaughter[1] + static_cast<int>(offset) : NotFound;
  }

  // Getters
  [[nodiscard]] int pdgCode() const { return mPdgCode; }
  [[nodiscard]] int globalIndex() const { return mGlobalIndex; }
  [[nodiscard]] int collisionId() const { return mCollisionId; }
  [[nodiscard]] bool isAlive() const { return checkBit(DecayerBits::IsAlive); }
  [[nodiscard]] bool isPrimary() const { return checkBit(DecayerBits::IsPrimary); }
  [[nodiscard]] bool isFromMcParticles() const { return !checkBit(DecayerBits::ProducedByDecayer); }
  [[nodiscard]] float weight() const
  {
    static constexpr float Weight = 1.f;
    return Weight;
  }
  [[nodiscard]] uint8_t flags() const { return mFlag; }
  [[nodiscard]] int statusCode() const { return mStatusCode; }
  [[nodiscard]] float vx() const { return mVx; }
  [[nodiscard]] float vy() const { return mVy; }
  [[nodiscard]] float vz() const { return mVz; }
  [[nodiscard]] float vt() const { return mVt; }
  [[nodiscard]] float px() const { return mPx; }
  [[nodiscard]] float py() const { return mPy; }
  [[nodiscard]] float pz() const { return mPz; }
  [[nodiscard]] float e() const { return mE; }
  [[nodiscard]] float radius() const { return std::hypot(mVx, mVy); }
  [[nodiscard]] float decayRadius() const { return mDecayRadius; }
  [[nodiscard]] float r() const { return radius(); }
  [[nodiscard]] float pt() const { return std::hypot(mPx, mPy); }
  [[nodiscard]] float p() const { return std::hypot(mPx, mPy, mPz); }
  [[nodiscard]] float phi() const { return o2::constants::math::PI + std::atan2(-1.0f * py(), -1.0f * px()); }
  [[nodiscard]] float eta() const
  {
    // Conditionally defined to avoid FPEs
    // As https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1959
    static constexpr float Tolerance = 1e-7f;
    if ((p() - mPz) < Tolerance) {
      return (mPz < 0.0f) ? -100.0f : 100.0f;
    }
    return 0.5f * std::log((p() + mPz) / (p() - mPz));
  }
  [[nodiscard]] float y() const
  {
    // Conditionally defined to avoid FPEs
    // As https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1980
    static constexpr float Tolerance = 1e-7f;
    if ((e() - mPz) < Tolerance) {
      return (mPz < 0.0f) ? -100.0f : 100.0f;
    }
    return 0.5f * std::log((mE + mPz) / (mE - mPz));
  }
  [[nodiscard]] int getMotherIndexStart() const { return mIndicesMother[0]; }
  [[nodiscard]] int getMotherIndexStop() const { return mIndicesMother[1]; }
  [[nodiscard]] int getDaughterIndexStart() const { return mIndicesDaughter[0]; }
  [[nodiscard]] int getDaughterIndexStop() const { return mIndicesDaughter[1]; }
  [[nodiscard]] const std::array<int, 2>& getMothers() const { return mIndicesMother; }
  [[nodiscard]] const std::array<int, 2>& getDaughters() const { return mIndicesDaughter; }
  [[nodiscard]] std::span<const int> getMotherSpan() const { return hasMothers() ? std::span<const int>(mIndicesMother.data(), 2) : std::span<const int>(); }

  // Checks
  [[nodiscard]] bool hasDaughters() const { return (mIndicesDaughter[0] >= 0); }
  [[nodiscard]] bool hasMothers() const { return (mIndicesMother[0] >= 0); }
  [[nodiscard]] bool hasNaN() const
  {
    return std::isnan(mPx) || std::isnan(mPy) || std::isnan(mPz) || std::isnan(mE) ||
           std::isnan(mVx) || std::isnan(mVy) || std::isnan(mVz);
  }
  [[nodiscard]] bool hasIndex() const
  {
    return (mGlobalIndex != -1);
  }

  // Bits
  [[nodiscard]] bool checkBit(DecayerBits bit) const { return mBits.test(static_cast<size_t>(bit)); }
  void setBit(DecayerBits bit, bool value = true) { mBits.set(static_cast<size_t>(bit), value); }
  void setBitOn(DecayerBits bit) { mBits.set(static_cast<size_t>(bit), true); }
  void setBitOff(DecayerBits bit) { mBits.set(static_cast<size_t>(bit), false); }

  [[nodiscard]] const std::bitset<8>& getBits() const { return mBits; }
  [[nodiscard]] uint8_t getBitsValue() const { return static_cast<uint8_t>(mBits.to_ulong()); }
  void setBits(std::bitset<8> bits) { mBits = bits; }

 private:
  int mPdgCode{}, mGlobalIndex{-1};
  int mCollisionId{-1};
  float mVx{}, mVy{}, mVz{}, mVt{};
  float mPx{}, mPy{}, mPz{}, mE{};
  float mDecayRadius{-1};

  int mStatusCode{};
  uint8_t mFlag{};
  std::bitset<8> mBits;
  std::array<int, 2> mIndicesMother{-1, -1}, mIndicesDaughter{-1, -1};
};

} // namespace o2::upgrade
#endif // ALICE3_CORE_OTFPARTICLE_H_
