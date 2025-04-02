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
/// \file HNMUtilities.h
///
/// \brief This code provides helper functions for the reconstruction of heavy neutral mesons (omega and eta meson) via their three pion decay
///
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt
///

#ifndef PWGEM_PHOTONMESON_UTILS_HNMUTILITIES_H_
#define PWGEM_PHOTONMESON_UTILS_HNMUTILITIES_H_

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <vector>
#include "TVector3.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/EventSelection.h"
#include "EMCALBase/Geometry.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/filterTables.h"

using namespace o2::aod::pwgem::photonmeson;

// -------> Struct to store photons from EMC clusters or V0s
namespace o2::aod::pwgem::photonmeson::hnmutilities
{
struct Photon {
  Photon(float px, float py, float pz, bool isFromConversion) : px(px), py(py), pz(pz), pt(std::sqrt(px * px + py * py)), isFromConversion(isFromConversion)
  {
    photon.SetPxPyPzE(px, py, pz, std::sqrt(px * px + py * py + pz * pz));
  }

  static Photon fromPxPyPz(float px, float py, float pz)
  {
    Photon photon(px, py, pz, true);
    return photon;
  }

  static Photon fromEtaPhiEnergy(float eta, float phi, float energy)
  {
    float theta = 2 * std::atan2(std::exp(-eta), 1);
    float px = energy * std::sin(theta) * std::cos(phi);
    float py = energy * std::sin(theta) * std::sin(phi);
    float pz = energy * std::cos(theta);
    Photon photon(px, py, pz, false);
    return photon;
  }

  ROOT::Math::PxPyPzEVector photon;
  float px, py, pz, pt;
  bool isFromConversion;
};

// -------> Struct to store gamma gamma pairs (pi0 or eta meson candidates)
struct GammaGammaPair {
  GammaGammaPair(Photon p1, Photon p2) : p1(p1), p2(p2)
  {
    vGG = p1.photon + p2.photon;
  }
  Photon p1, p2;
  ROOT::Math::PxPyPzEVector vGG;

  bool isPi0 = false;
  bool isEta = false;
  ushort reconstructionType;
  void setReconstructionType(ushort type) { reconstructionType = type; }

  float m() const { return vGG.M(); }
  float pT() const { return vGG.Pt(); }
};

// -------> Enum to specify how the heavy neutral meson mass should be corrected based on the PDG mass of its light neutral meson decay daughter
enum MassCorrectionType {
  kNoHNMMassCorrection = 0,
  kSubDeltaPi0 = 1,
  kSubLambda = 2
};

struct HeavyNeutralMeson {
  HeavyNeutralMeson(GammaGammaPair* gg, float eTracks, float pxTracks, float pyTracks, float pzTracks) : gg(gg)
  {
    vHeavyNeutralMeson.SetPxPyPzE(gg->vGG.Px() + pxTracks, gg->vGG.Py() + pyTracks, gg->vGG.Pz() + pzTracks, gg->vGG.e() + eTracks);
  }

  GammaGammaPair* gg = nullptr;
  ROOT::Math::PxPyPzEVector vHeavyNeutralMeson;

  float m(int massCorrectionType) const
  {
    float massHNM = vHeavyNeutralMeson.M();
    switch (massCorrectionType) {
      case kNoHNMMassCorrection: // No mass correction
        break;
      case kSubDeltaPi0: // Subtract the mass difference of the reconstructed light neutral meson mass to the PDG mass
        massHNM -= gg->m();
        massHNM += (gg->isPi0 ? constants::physics::MassPi0 : 0.547862);
        break;
      case kSubLambda: // Subtract an opening angle dependent mass correction (see https://arxiv.org/abs/2502.19956 for details)
        LOGF(warning, "SubLambdaMassCorrection not yet implemented!");
        break;
      default:
        LOGF(fatal, "Unknown mass correction type %d", massCorrectionType);
    }
    return massHNM;
  }
  float pT() const { return vHeavyNeutralMeson.Pt(); }
};

/// \brief Reconstruct light neutral mesons from EMC clusters and V0s and fill them into the vGGs vector
template <typename C, typename V>
void reconstructGGs(C clusters, V v0s, std::vector<GammaGammaPair>& vGGs)
{
  std::vector<Photon> vPhotons;
  for (const auto& cluster : clusters)
    vPhotons.push_back(Photon::fromEtaPhiEnergy(cluster.eta(), cluster.phi(), cluster.e()));

  for (const auto& v0 : v0s)
    vPhotons.push_back(Photon::fromPxPyPz(v0.px(), v0.py(), v0.pz()));

  vGGs.clear();
  // loop over all photon combinations and build meson candidates
  for (unsigned int ig1 = 0; ig1 < vPhotons.size(); ++ig1) {
    for (unsigned int ig2 = ig1 + 1; ig2 < vPhotons.size(); ++ig2) {
      GammaGammaPair lightMeson(vPhotons[ig1], vPhotons[ig2]); // build lightMeson from photons
      if (vPhotons[ig1].isFromConversion && vPhotons[ig2].isFromConversion)
        lightMeson.setReconstructionType(photonpair::kPCMPCM);
      else if (!vPhotons[ig1].isFromConversion && !vPhotons[ig2].isFromConversion)
        lightMeson.setReconstructionType(photonpair::kEMCEMC);
      else
        lightMeson.setReconstructionType(photonpair::kPCMEMC);

      vGGs.push_back(lightMeson);
    }
  }
}

/// \brief Reconstruct heavy neutral mesons from tracks and GG candidates and fill them into the vHNMs vector
template <typename Track>
void reconstructHeavyNeutralMesons(Track const& tracks, std::vector<GammaGammaPair>& vGGs, std::vector<HeavyNeutralMeson>& vHNMs)
{
  vHNMs.clear();
  for (const auto& posTrack : tracks) {
    if (!posTrack.isGlobalTrack() || posTrack.sign() < 0)
      continue;
    for (const auto& negTrack : tracks) {
      if (!negTrack.isGlobalTrack() || negTrack.sign() > 0)
        continue;
      for (size_t iGG = 0; iGG < vGGs.size(); iGG++) {
        HeavyNeutralMeson heavyNeutralMeson(&vGGs.at(iGG), posTrack.energy(constants::physics::MassPiPlus) + negTrack.energy(constants::physics::MassPiMinus), posTrack.px() + negTrack.px(), posTrack.py() + negTrack.py(), posTrack.pz() + negTrack.pz());
        vHNMs.push_back(heavyNeutralMeson);
      }
    }
  }
}

} // namespace o2::aod::pwgem::photonmeson::hnmutilities

#endif // PWGEM_PHOTONMESON_UTILS_HNMUTILITIES_H_
