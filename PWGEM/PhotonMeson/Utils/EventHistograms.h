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

/// \header file for histograms
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_EVENTHISTOGRAMS_H_
#define PWGEM_PHOTONMESON_UTILS_EVENTHISTOGRAMS_H_
using namespace o2::framework;

namespace o2::aod::pwgem::photonmeson::utils::eventhistogram
{
void addEventHistograms(HistogramRegistry* fRegistry, bool doFlow)
{
  // event info
  auto hCollisionCounter = fRegistry->add<TH1>("Event/before/hCollisionCounter", "collision counter;;Number of events", kTH1F, {{10, 0.5, 10.5}}, false);
  hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
  hCollisionCounter->GetXaxis()->SetBinLabel(2, "No TF border");
  hCollisionCounter->GetXaxis()->SetBinLabel(3, "No ITS ROF border");
  hCollisionCounter->GetXaxis()->SetBinLabel(4, "No Same Bunch Pileup");
  hCollisionCounter->GetXaxis()->SetBinLabel(5, "Is Vertex ITSTPC");
  hCollisionCounter->GetXaxis()->SetBinLabel(6, "Is Good Zvtx FT0vsPV");
  hCollisionCounter->GetXaxis()->SetBinLabel(7, "FT0AND");
  hCollisionCounter->GetXaxis()->SetBinLabel(8, "sel8");
  hCollisionCounter->GetXaxis()->SetBinLabel(9, "|Z_{vtx}| < 10 cm");
  hCollisionCounter->GetXaxis()->SetBinLabel(10, "accepted");

  fRegistry->add("Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
  fRegistry->add("Event/before/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
  fRegistry->add("Event/before/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
  fRegistry->add("Event/before/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", kTH2F, {{200, 0, 200000}, {60, 0, 60000}}, false);
  fRegistry->add("Event/before/hCentFT0A", "hCentFT0A;centrality FT0A (%)", kTH1F, {{110, 0, 110}}, false);
  fRegistry->add("Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1F, {{110, 0, 110}}, false);
  fRegistry->add("Event/before/hCentFT0M", "hCentFT0M;centrality FT0M (%)", kTH1F, {{110, 0, 110}}, false);
  fRegistry->add("Event/before/hCentFT0CvsMultNTracksPV", "hCentFT0CvsMultNTracksPV;centrality FT0C (%);N_{track} to PV", kTH2F, {{110, 0, 110}, {500, 0, 5000}}, false);
  fRegistry->add("Event/before/hMultFT0CvsMultNTracksPV", "hMultFT0CvsMultNTracksPV;mult. FT0C;N_{track} to PV", kTH2F, {{60, 0, 60000}, {500, 0, 5000}}, false);
  fRegistry->add("Event/before/hMultFT0CvsOccupancy", "hMultFT0CvsOccupancy;mult. FT0C;N_{track} in time range", kTH2F, {{60, 0, 60000}, {2000, 0, 20000}}, false);

  if (doFlow) {
    fRegistry->add("Event/before/hQ2xFT0M_CentFT0C", "hQ2xFT0M_CentFT0C;centrality FT0C (%);Q_{2,x}^{FT0M}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFT0M_CentFT0C", "hQ2yFT0M_CentFT0C;centrality FT0C (%);Q_{2,y}^{FT0M}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xFT0A_CentFT0C", "hQ2xFT0A_CentFT0C;centrality FT0C (%);Q_{2,x}^{FT0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFT0A_CentFT0C", "hQ2yFT0A_CentFT0C;centrality FT0C (%);Q_{2,y}^{FT0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xFT0C_CentFT0C", "hQ2xFT0C_CentFT0C;centrality FT0C (%);Q_{2,x}^{FT0C}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFT0C_CentFT0C", "hQ2yFT0C_CentFT0C;centrality FT0C (%);Q_{2,y}^{FT0C}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xFV0A_CentFT0C", "hQ2xFV0A_CentFT0C;centrality FT0C (%);Q_{2,x}^{FV0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFV0A_CentFT0C", "hQ2yFV0A_CentFT0C;centrality FT0C (%);Q_{2,y}^{FV0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xBPos_CentFT0C", "hQ2xBPos_CentFT0C;centrality FT0C (%);Q_{2,x}^{BPos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yBPos_CentFT0C", "hQ2yBPos_CentFT0C;centrality FT0C (%);Q_{2,y}^{BPos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xBNeg_CentFT0C", "hQ2xBNeg_CentFT0C;centrality FT0C (%);Q_{2,x}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yBNeg_CentFT0C", "hQ2yBNeg_CentFT0C;centrality FT0C (%);Q_{2,y}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);

    fRegistry->add("Event/before/hEP2FT0M_CentFT0C", "2nd harmonics event plane FT0M;centrality FT0C (%);#Psi_{2}^{FT0M} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2FT0A_CentFT0C", "2nd harmonics event plane FT0A;centrality FT0C (%);#Psi_{2}^{FT0A} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2FT0C_CentFT0C", "2nd harmonics event plane FT0C;centrality FT0C (%);#Psi_{2}^{FT0C} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2FV0A_CentFT0C", "2nd harmonics event plane FV0A;centrality FT0C (%);#Psi_{2}^{FV0A} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2BPos_CentFT0C", "2nd harmonics event plane BPos;centrality FT0C (%);#Psi_{2}^{BPos} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2BNeg_CentFT0C", "2nd harmonics event plane BNeg;centrality FT0C (%);#Psi_{2}^{BNeg} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);

    fRegistry->add("Event/before/hQ2FT0MQ2BPos_CentFT0C", "Q_{2}^{FT0M} #upoint Q_{2}^{BPos};centrality FT0C (%);Q_{2}^{FT0M} #upoint Q_{2}^{BPos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2FT0MQ2BNeg_CentFT0C", "Q_{2}^{FT0M} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{FT0M} #upoint Q_{2}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2BPosQ2BNeg_CentFT0C", "Q_{2}^{BPos} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{BPos} #upoint Q_{2}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false); // this is common for FT0M, FT0A, FT0C, FV0A resolution.

    fRegistry->add("Event/before/hQ2FT0CQ2BPos_CentFT0C", "Q_{2}^{FT0C} #upoint Q_{2}^{BPos};centrality FT0C (%);Q_{2}^{FT0C} #upoint Q_{2}^{BPos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2FT0CQ2BNeg_CentFT0C", "Q_{2}^{FT0C} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{FT0C} #upoint Q_{2}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);

    fRegistry->add("Event/before/hQ2FT0AQ2BPos_CentFT0C", "Q_{2}^{FT0A} #upoint Q_{2}^{BPos};centrality FT0C (%);Q_{2}^{FT0A} #upoint Q_{2}^{BPos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2FT0AQ2BNeg_CentFT0C", "Q_{2}^{FT0A} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{FT0A} #upoint Q_{2}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);

    fRegistry->add("Event/before/hQ2FV0AQ2BPos_CentFT0C", "Q_{2}^{FV0A} #upoint Q_{2}^{BPos};centrality FT0C (%);Q_{2}^{FV0A} #upoint Q_{2}^{BPos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2FV0AQ2BNeg_CentFT0C", "Q_{2}^{FV0A} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{FV0A} #upoint Q_{2}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
  }
  fRegistry->addClone("Event/before/", "Event/after/");
}

template <const int ev_id, typename TCollision>
void fillEventInfo(HistogramRegistry* fRegistry, TCollision const& collision, const bool doFlow, const float weight = 1.f)
{
  static constexpr std::string_view event_types[2] = {"before/", "after/"};
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 1.0);
  if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 2.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 3.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 4.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 5.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 6.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 7.0);
  }
  if (collision.sel8()) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 8.0);
  }
  if (abs(collision.posZ()) < 10.0) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 9.0);
  }
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hZvtx"), collision.posZ());

  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPV"), collision.multNTracksPV());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPVeta1"), collision.multNTracksPVeta1());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0"), collision.multFT0A(), collision.multFT0C());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0A"), collision.centFT0A());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0C"), collision.centFT0C());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0M"), collision.centFT0M());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0CvsMultNTracksPV"), collision.centFT0C(), collision.multNTracksPV());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0CvsMultNTracksPV"), collision.multFT0C(), collision.multNTracksPV());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0CvsOccupancy"), collision.multFT0C(), collision.trackOccupancyInTimeRange());

  if (doFlow) {
    std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
    std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
    std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
    std::array<float, 2> q2fv0a = {collision.q2xfv0a(), collision.q2yfv0a()};
    std::array<float, 2> q2bpos = {collision.q2xbpos(), collision.q2ybpos()};
    std::array<float, 2> q2bneg = {collision.q2xbneg(), collision.q2ybneg()};

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFT0M_CentFT0C"), collision.centFT0C(), collision.q2xft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFT0M_CentFT0C"), collision.centFT0C(), collision.q2yft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFT0A_CentFT0C"), collision.centFT0C(), collision.q2xft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFT0A_CentFT0C"), collision.centFT0C(), collision.q2yft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFT0C_CentFT0C"), collision.centFT0C(), collision.q2xft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFT0C_CentFT0C"), collision.centFT0C(), collision.q2yft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFV0A_CentFT0C"), collision.centFT0C(), collision.q2xfv0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFV0A_CentFT0C"), collision.centFT0C(), collision.q2yfv0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xBPos_CentFT0C"), collision.centFT0C(), collision.q2xbpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yBPos_CentFT0C"), collision.centFT0C(), collision.q2ybpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xBNeg_CentFT0C"), collision.centFT0C(), collision.q2xbneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yBNeg_CentFT0C"), collision.centFT0C(), collision.q2ybneg());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FT0M_CentFT0C"), collision.centFT0C(), collision.ep2ft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FT0A_CentFT0C"), collision.centFT0C(), collision.ep2ft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FT0C_CentFT0C"), collision.centFT0C(), collision.ep2ft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FV0A_CentFT0C"), collision.centFT0C(), collision.ep2fv0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2BPos_CentFT0C"), collision.centFT0C(), collision.ep2bpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2BNeg_CentFT0C"), collision.centFT0C(), collision.ep2bneg());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FT0MQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0m, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FT0MQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0m, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FT0AQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0a, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FT0AQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0a, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FT0CQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0c, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FT0CQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0c, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FV0AQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2fv0a, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FV0AQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2fv0a, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2BPosQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2bpos, q2bneg));
  }
}

} // namespace o2::aod::pwgem::photonmeson::utils::eventhistogram

#endif // PWGEM_PHOTONMESON_UTILS_EVENTHISTOGRAMS_H_
