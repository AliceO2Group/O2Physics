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
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/filterTables.h"

namespace o2::aod::pwgem::photonmeson::utils::hnmutilities
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

namespace ReconstructionType
{
enum ReconstructionType {
  kPCM = 0,
  kPCMEMC = 1,
  kEMC = 2
};
} // namespace ReconstructionType

enum MassCorrectionType {
  kNoHNMMassCorrection = 0,
  kSubDeltaPi0 = 1,
  kSubLambda = 2
};

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
      case kSubDeltaPi0:
        massHNM -= gg->m();
        massHNM += (gg->isPi0 ? constants::physics::MassPi0 : 0.547862);
        break;
      case kSubLambda:
        LOGF(warning, "SubLambdaMassCorrection not yet implemented!");
        break;
      default:
        LOGF(fatal, "Unknown mass correction type %d", massCorrectionType);
    }
    return massHNM;
  }
  float pT() const { return vHeavyNeutralMeson.Pt(); }
};

/// \brief Process light neutral meson candidates (pi0 or eta), calculate invariant mass and pT and fill histograms
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
        lightMeson.setReconstructionType(ReconstructionType::kPCM);
      else if (!vPhotons[ig1].isFromConversion && !vPhotons[ig2].isFromConversion)
        lightMeson.setReconstructionType(ReconstructionType::kEMC);
      else
        lightMeson.setReconstructionType(ReconstructionType::kPCMEMC);

      vGGs.push_back(lightMeson);
    }
  }
}

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
      for (unsigned long iGG = 0; iGG < vGGs.size(); iGG++) {
        HeavyNeutralMeson heavyNeutralMeson(&vGGs.at(iGG), posTrack.energy(constants::physics::MassPiPlus) + negTrack.energy(constants::physics::MassPiMinus), posTrack.px() + negTrack.px(), posTrack.py() + negTrack.py(), posTrack.pz() + negTrack.pz());
        vHNMs.push_back(heavyNeutralMeson);
      }
    }
  }
}

} // namespace o2::aod::pwgem::photonmeson::utils::hnmutilities

#endif // PWGEM_PHOTONMESON_UTILS_HNMUTILITIES_H_
