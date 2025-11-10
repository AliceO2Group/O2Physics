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

/// \file LFHypernucleiKfTables.h
/// \brief Slim hypernuclei kf tables
/// \author Janik Ditzel <jditzel@cern.ch> and Michael Hartung <mhartung@cern.ch>

#ifndef PWGLF_DATAMODEL_LFHYPERNUCLEIKFTABLES_H_
#define PWGLF_DATAMODEL_LFHYPERNUCLEIKFTABLES_H_

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cmath>

namespace o2::aod
{
namespace hykfmccoll
{
DECLARE_SOA_COLUMN(PassedEvSel, passedEvSel, bool); //!
}
DECLARE_SOA_TABLE(HypKfMcColls, "AOD", "HYPKFMCCOLL",
                  o2::soa::Index<>,
                  hykfmccoll::PassedEvSel,
                  mccollision::PosX,
                  mccollision::PosY,
                  mccollision::PosZ);
using HypKfMcColl = HypKfMcColls::iterator;

namespace hykfmc
{
DECLARE_SOA_INDEX_COLUMN(HypKfMcColl, hypKfMcColl);
DECLARE_SOA_COLUMN(Species, species, int8_t);                   //!
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool); //!
DECLARE_SOA_COLUMN(Svx, svx, float);                            //!
DECLARE_SOA_COLUMN(Svy, svy, float);                            //!
DECLARE_SOA_COLUMN(Svz, svz, float);                            //!
DECLARE_SOA_COLUMN(Occupancy, occupancy, int);                  //!
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);                  //!
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) { return RecoDecay::pt(std::array{px, py}); });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, [](float E, float pz) { return 0.5 * std::log((E + pz) / (E - pz)); });
DECLARE_SOA_DYNAMIC_COLUMN(Mass, mass, [](float E, float px, float py, float pz) { return std::sqrt(E * E - px * px - py * py - pz * pz); });
DECLARE_SOA_DYNAMIC_COLUMN(IsMatter, isMatter, [](int pdgCode) { return pdgCode > 0; });
} // namespace hykfmc

DECLARE_SOA_TABLE(HypKfMcParts, "AOD", "HYPKFMCPART",
                  o2::soa::Index<>,
                  hykfmc::HypKfMcCollId,
                  hykfmc::Species,
                  mcparticle::PdgCode,
                  hykfmc::IsPhysicalPrimary,
                  mcparticle::Px,
                  mcparticle::Py,
                  mcparticle::Pz,
                  mcparticle::E,
                  hykfmc::Svx,
                  hykfmc::Svy,
                  hykfmc::Svz,
                  hykfmc::Pt<mcparticle::Px, mcparticle::Py>,
                  hykfmc::Y<mcparticle::E, mcparticle::Pz>,
                  hykfmc::Mass<mcparticle::E, mcparticle::Px, mcparticle::Py, mcparticle::Pz>,
                  hykfmc::IsMatter<mcparticle::PdgCode>);
using HypKfMcPart = HypKfMcParts::iterator;

DECLARE_SOA_TABLE(HypKfColls, "AOD", "HYPKFCOLL",
                  o2::soa::Index<>,
                  hykfmccoll::PassedEvSel,
                  hykfmc::HypKfMcCollId,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  cent::CentFT0A,
                  cent::CentFT0C,
                  cent::CentFT0M,
                  hykfmc::Occupancy,
                  hykfmc::RunNumber);
using HypKfColl = HypKfColls::iterator;

namespace hykftrk
{
DECLARE_SOA_INDEX_COLUMN(HypKfColl, hypKfColl);
DECLARE_SOA_COLUMN(Rigidity, rigidity, float);              //!
DECLARE_SOA_COLUMN(TpcNcluster, tpcNcluster, float);        //!
DECLARE_SOA_COLUMN(TpcNsigma, tpcNsigma, float);            //!
DECLARE_SOA_COLUMN(TpcNsigmaNhp, tpcNsigmaNhp, float);      //!
DECLARE_SOA_COLUMN(TpcNsigmaNlp, tpcNsigmaNlp, float);      //!
DECLARE_SOA_COLUMN(TofMass, tofMass, float);                //!
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool); //!
DECLARE_SOA_COLUMN(SubMass, subMass, float);                //!
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float pt, float phi) { return (double)pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float pt, float phi) { return (double)pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float pt, float eta) { return (double)pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float pt, float eta) { return pt * std::cosh(eta); }); //
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, [](float pt, float eta, float mass) { return std::log((RecoDecay::sqrtSumOfSquares(mass, pt * std::cosh(eta)) + pt * std::sinh(eta)) / RecoDecay::sqrtSumOfSquares(mass, pt)); });
DECLARE_SOA_DYNAMIC_COLUMN(Lambda, lambda, [](float eta) { return 1. / std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(ItsNcluster, itsNcluster, [](uint32_t itsClusterSizes) {
  uint8_t n = 0;
  for (uint8_t i = 0; i < 0x08; i++) {
    if (itsClusterSizes >> (4 * i) & 15)
      n++;
  }
  return n;
});
DECLARE_SOA_DYNAMIC_COLUMN(ItsFirstLayer, itsFirstLayer, [](uint32_t itsClusterSizes) {
  for (int i = 0; i < 0x08; i++) {
    if (itsClusterSizes >> (4 * i) & 15)
      return i;
  }
  return -999;
});
DECLARE_SOA_DYNAMIC_COLUMN(ItsMeanClsSize, itsMeanClsSize, [](uint32_t itsClusterSizes) {
  int sum = 0, n = 0;
  for (int i = 0; i < 0x08; i++) {
    sum += (itsClusterSizes >> (4 * i) & 15);
    if (itsClusterSizes >> (4 * i) & 15)
      n++;
  }
  return n > 0 ? static_cast<float>(sum) / n : 0.f;
});
} // namespace hykftrk

DECLARE_SOA_TABLE(HypKfTracks, "AOD", "HYPKFTRACK",
                  o2::soa::Index<>,
                  hykfmc::Species,
                  track::Pt,
                  track::Eta,
                  track::Phi,
                  track::DcaXY,
                  track::DcaZ,
                  hykftrk::TpcNcluster,
                  track::TPCChi2NCl,
                  track::ITSClusterSizes,
                  track::ITSChi2NCl,
                  hykftrk::Rigidity,
                  track::TPCSignal,
                  hykftrk::TpcNsigma,
                  hykftrk::TpcNsigmaNhp,
                  hykftrk::TpcNsigmaNlp,
                  hykftrk::TofMass,
                  hykftrk::IsPVContributor,
                  hykftrk::Px<track::Pt, track::Phi>,
                  hykftrk::Py<track::Pt, track::Phi>,
                  hykftrk::Pz<track::Pt, track::Eta>,
                  hykftrk::P<track::Pt, track::Eta>,
                  hykftrk::Lambda<track::Eta>,
                  hykftrk::ItsNcluster<track::ITSClusterSizes>,
                  hykftrk::ItsFirstLayer<track::ITSClusterSizes>,
                  hykftrk::ItsMeanClsSize<track::ITSClusterSizes>);
using HypKfTrack = HypKfTracks::iterator;

DECLARE_SOA_TABLE(HypKfSubDs, "AOD", "HYPKFSUBD",
                  o2::soa::Index<>,
                  hykftrk::SubMass);
using HypKfSubD = HypKfSubDs::iterator;

DECLARE_SOA_TABLE(HypKfDaughtAdds, "AOD", "HYPKFDAUGHTADD",
                  o2::soa::Index<>,
                  track::X,
                  track::Y,
                  track::Z,
                  mcparticle::Px,
                  mcparticle::Py,
                  mcparticle::Pz);
using HypKfDaughtAdd = HypKfDaughtAdds::iterator;

namespace hykfhyp
{
DECLARE_SOA_INDEX_COLUMN(HypKfColl, hypKfColl);
DECLARE_SOA_INDEX_COLUMN(HypKfMcPart, hypKfMcPart);
DECLARE_SOA_ARRAY_INDEX_COLUMN(HypKfDaughtAdd, hypKfDaughtAdd);
DECLARE_SOA_ARRAY_INDEX_COLUMN(HypKfTrack, hypKfTrack);
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(HypDaughter, hypDaughter, int, "HypKfHypNucs");
DECLARE_SOA_ARRAY_INDEX_COLUMN(HypKfSubD, hypKfSubD);
DECLARE_SOA_COLUMN(Primary, primary, bool);        //!
DECLARE_SOA_COLUMN(Mass, mass, float);             //!
DECLARE_SOA_COLUMN(Px, px, float);                 //!
DECLARE_SOA_COLUMN(Py, py, float);                 //!
DECLARE_SOA_COLUMN(Pz, pz, float);                 //!
DECLARE_SOA_COLUMN(DcaToPvXY, dcaToPvXY, float);   //!
DECLARE_SOA_COLUMN(DcaToPvZ, dcaToPvZ, float);     //!
DECLARE_SOA_COLUMN(DcaToVtxXY, dcaToVtxXY, float); //!
DECLARE_SOA_COLUMN(DcaToVtxZ, dcaToVtxZ, float);   //!
DECLARE_SOA_COLUMN(Chi2, chi2, float);             //!
DECLARE_SOA_COLUMN(DevToPvXY, devToPvXY, float);   //!
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) { return RecoDecay::phi(std::array{px, py}); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float px, float py, float pz) { return RecoDecay::p(px, py, pz); }); //
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, [](float px, float py, float pz, float mass) { return RecoDecay::y(std::array{px, py, pz}, mass); });
DECLARE_SOA_DYNAMIC_COLUMN(McTrue, mcTrue, [](int hypKfMcPartId) { return hypKfMcPartId >= 0; });
DECLARE_SOA_DYNAMIC_COLUMN(IsMatter, isMatter, [](int8_t species) { return species > 0; });
DECLARE_SOA_DYNAMIC_COLUMN(Cascade, cascade, [](int hypDaughter) { return hypDaughter >= 0; });
} // namespace hykfhyp

DECLARE_SOA_TABLE(HypKfHypNucs, "AOD", "HYPKFHYPNUC",
                  o2::soa::Index<>,
                  hykfhyp::HypKfMcPartId,
                  hykfhyp::HypKfCollId,
                  hykfhyp::HypKfTrackIds,
                  hykfhyp::HypKfDaughtAddIds,
                  hykfhyp::HypDaughterId,
                  hykfhyp::HypKfSubDIds,
                  hykfmc::Species,
                  hykfhyp::Primary,
                  hykfhyp::Mass,
                  hykfhyp::Px,
                  hykfhyp::Py,
                  hykfhyp::Pz,
                  hykfhyp::DcaToPvXY,
                  hykfhyp::DcaToPvZ,
                  hykfhyp::DevToPvXY,
                  hykfhyp::DcaToVtxXY,
                  hykfhyp::DcaToVtxZ,
                  hykfhyp::Chi2,
                  hykfmc::Svx,
                  hykfmc::Svy,
                  hykfmc::Svz,
                  hykfhyp::Y<hykfhyp::Px, hykfhyp::Py, hykfhyp::Pz, hykfhyp::Mass>,
                  hykfhyp::Pt<hykfhyp::Px, hykfhyp::Py>,
                  hykfhyp::Eta<hykfhyp::Px, hykfhyp::Py, hykfhyp::Pz>,
                  hykfhyp::Phi<hykfhyp::Px, hykfhyp::Py>,
                  hykfhyp::P<hykfhyp::Px, hykfhyp::Py, hykfhyp::Pz>,
                  hykfhyp::McTrue<hykfhyp::HypKfMcPartId>,
                  hykfhyp::IsMatter<hykfmc::Species>,
                  hykfhyp::Cascade<hykfhyp::HypDaughterId>);
using HypKfHypNuc = HypKfHypNucs::iterator;
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFHYPERNUCLEIKFTABLES_H_
