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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <cmath>
#include <array>
#include <utility>

#include <TF1.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

#define NPARTICLES 8

const std::vector<double> ptBins = {0.0, 0.05, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                                    0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4,
                                    1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                    4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 20.0};
const double mass[NPARTICLES] = {0.000510999, 0.139570, 0.493677, 0.938272, 1.87561, 2.8083916, 2.8089211, 3.7273794};
enum PType : uint8_t {
  kEl,
  kPi,
  kKa,
  kPr,
  kDe,
  kHe,
  kTr,
  kAl,
  kNull
};

// function to convert eta to y
double eta2y(double pt, double m, double eta)
{
  double mt = sqrt(m * m + pt * pt);
  return asinh(pt / mt * sinh(eta));
}

int getparticleint(int pdgcode)
{
  if (pdgcode == 11) {
    return kEl;
  } else if (pdgcode == 211) {
    return kPi;
  } else if (pdgcode == 321) {
    return kKa;
  } else if (pdgcode == 2212) {
    return kPr;
  } else if (pdgcode == 1000010020) {
    return kDe;
  } else if (pdgcode == 1000020030) {
    return kHe;
  } else if (pdgcode == 1000010030) {
    return kTr;
  } else if (pdgcode == 1000020040) {
    return kAl;
  } else {
    return kNull;
  }
}

// No track selection --> only event selection here
struct TrackCheckTaskEvSel {

  Configurable<float> cfgCutVZ{"cfgCutVZ", 10.f, "option to configure z-vertex cut"};
  Configurable<bool> cfgIsRun3{"cfgIsRun3", true, "option to set Run 3 format"};

  HistogramRegistry histograms{"histograms"};

  void init(InitContext&)
  {
    AxisSpec dcaAxis = {800, -4., 4.};
    AxisSpec yAxis = {40, -2., 2.};

    histograms.add("GenAftEvSel/hGenTrkPrimAftEvSel", "Gen Prim tracks AftEvSel (charged); #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_el", "Gen Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_pi", "Gen Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_ka", "Gen Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_pr", "Gen Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_de", "Gen Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_he", "Gen Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_tr", "Gen Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_al", "Gen Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH2F, {{ptBins}, yAxis}});

    histograms.add("RecAftEvSel/hRecTrkPrimAftEvSel", "Reco Prim tracks AftEvSel (charged); #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_el", "Reco Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_pi", "Reco Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_ka", "Reco Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_pr", "Reco Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_de", "Reco Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_he", "Reco Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_tr", "Reco Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_al", "Reco Prim tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
  }

  // Filters
  Filter collfilter = nabs(aod::collision::posZ) < cfgCutVZ;
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& col,
               soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels>& tracks, aod::McParticles const& mcParticles)
  {

    // event selection
    if (cfgIsRun3) {
      if (!col.sel8()) {
        return;
      }
    } else {
      if (!col.sel7()) {
        return;
      }
    }

    // Loop on tracks
    for (auto& track : tracks) {

      if (!track.has_mcParticle()) {
        continue;
      }
      const auto particle = track.mcParticle_as<aod::McParticles>();
      int pdgcode = fabs(particle.pdgCode());

      if (!particle.isPhysicalPrimary()) {
        continue;
      }

      // calculate rapidity
      int pint = getparticleint(pdgcode);

      if (pint == kNull) {
        continue;
      }

      double y_gen = eta2y(particle.pt(), mass[pint], particle.eta());
      double y_rec = eta2y(track.pt(), mass[pint], track.eta());

      // generated
      if (pdgcode == 11) {
        histograms.fill(HIST("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_el"), particle.pt(), y_gen);
      } else if (pdgcode == 211) {
        histograms.fill(HIST("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_pi"), particle.pt(), y_gen);
      } else if (pdgcode == 321) {
        histograms.fill(HIST("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_ka"), particle.pt(), y_gen);
      } else if (pdgcode == 2212) {
        histograms.fill(HIST("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_pr"), particle.pt(), y_gen);
      } else if (pdgcode == 1000010020) {
        histograms.fill(HIST("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_de"), particle.pt(), y_gen);
      } else if (pdgcode == 1000020030) {
        histograms.fill(HIST("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_he"), particle.pt(), y_gen);
      } else if (pdgcode == 1000010030) {
        histograms.fill(HIST("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_tr"), particle.pt(), y_gen);
      } else if (pdgcode == 1000020040) {
        histograms.fill(HIST("GenAftEvSel/hGenTrkPrimAftEvSel_truepid_al"), particle.pt(), y_gen);
      }

      histograms.fill(HIST("GenAftEvSel/hGenTrkPrimAftEvSel"), particle.pt(), y_gen);

      // reconstructed
      if (pdgcode == 11) {
        histograms.fill(HIST("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_el"), track.pt(), y_rec);
      } else if (pdgcode == 211) {
        histograms.fill(HIST("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_pi"), track.pt(), y_rec);
      } else if (pdgcode == 321) {
        histograms.fill(HIST("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_ka"), track.pt(), y_rec);
      } else if (pdgcode == 2212) {
        histograms.fill(HIST("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_pr"), track.pt(), y_rec);
      } else if (pdgcode == 1000010020) {
        histograms.fill(HIST("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_de"), track.pt(), y_rec);
      } else if (pdgcode == 1000020030) {
        histograms.fill(HIST("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_he"), track.pt(), y_rec);
      } else if (pdgcode == 1000010030) {
        histograms.fill(HIST("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_tr"), track.pt(), y_rec);
      } else if (pdgcode == 1000020040) {
        histograms.fill(HIST("RecAftEvSel/hRecTrkPrimAftEvSel_truepid_al"), track.pt(), y_rec);
      }

      histograms.fill(HIST("RecAftEvSel/hRecTrkPrimAftEvSel"), track.pt(), y_rec);
    }
  }
};

// event selection + track selection here
struct TrackCheckTaskEvSelTrackSel {

  Configurable<float> cfgCutVZ{"cfgCutVZ", 10.f, "option to configure z-vertex cut"};
  Configurable<bool> cfgIsRun3{"cfgIsRun3", true, "option to set Run 3 format"};

  HistogramRegistry histograms{"histograms"};

  void init(InitContext&)
  {
    AxisSpec dcaAxis = {800, -4., 4.};
    AxisSpec yAxis = {40, -2., 2.};

    histograms.add("GenAftTrkSel/hGenTrkPrimAftTrkSel", "Gen Prim tracks AftTrkSel (charged); #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_el", "Gen Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_pi", "Gen Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_ka", "Gen Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_pr", "Gen Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_de", "Gen Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_he", "Gen Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_tr", "Gen Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});
    histograms.add("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_al", "Gen Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y", {kTH2F, {{ptBins}, yAxis}});

    histograms.add("RecAftTrkSel/hRecTrkPrimAftTrkSel", "Reco Prim tracks AftTrkSel (charged); #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_el", "Reco Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_pi", "Reco Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_ka", "Reco Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_pr", "Reco Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_de", "Reco Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_he", "Reco Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_tr", "Reco Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_al", "Reco Prim tracks AftTrkSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});

    histograms.add("RecAftTrkSel/hRecTrkSecAftTrkSel", "Reco Sec tracks AftEvSel (charged); #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_el", "Reco Sec tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_pi", "Reco Sec tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_ka", "Reco Sec tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_pr", "Reco Sec tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_de", "Reco Sec tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_he", "Reco Sec tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_tr", "Reco Sec tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
    histograms.add("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_al", "Reco Sec tracks AftEvSel; #it{p}_{T} (GeV/#it{c}); y; DCA_{xy} (cm)", {kTH3F, {{ptBins}, yAxis, dcaAxis}});
  }

  // Filters
  Filter collfilter = nabs(aod::collision::posZ) < cfgCutVZ;
  Filter trackfilter = requireGlobalTrackInFilter();
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& col,
               soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                       aod::TrackSelection, aod::McTrackLabels>>& tracks,
               aod::McParticles& mcParticles)
  {

    // event selection
    if (cfgIsRun3) {
      if (!col.sel8()) {
        return;
      }
    } else {
      if (!col.sel7()) {
        return;
      }
    }

    // Loop on tracks
    for (auto& track : tracks) {

      if (!track.has_mcParticle()) {
        continue;
      }
      const auto particle = track.mcParticle_as<aod::McParticles>();
      int pdgcode = fabs(particle.pdgCode());

      int pint = getparticleint(pdgcode);
      if (pint == kNull) {
        continue;
      }

      double y_gen = eta2y(particle.pt(), mass[pint], particle.eta());
      double y_rec = eta2y(track.pt(), mass[pint], track.eta());

      if (particle.isPhysicalPrimary()) {
        if (pdgcode == 11) {
          histograms.fill(HIST("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_el"), particle.pt(), y_gen);
        } else if (pdgcode == 211) {
          histograms.fill(HIST("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_pi"), particle.pt(), y_gen);
        } else if (pdgcode == 321) {
          histograms.fill(HIST("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_ka"), particle.pt(), y_gen);
        } else if (pdgcode == 2212) {
          histograms.fill(HIST("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_pr"), particle.pt(), y_gen);
        } else if (pdgcode == 1000010020) {
          histograms.fill(HIST("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_de"), particle.pt(), y_gen);
        } else if (pdgcode == 1000020030) {
          histograms.fill(HIST("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_he"), particle.pt(), y_gen);
        } else if (pdgcode == 1000010030) {
          histograms.fill(HIST("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_tr"), particle.pt(), y_gen);
        } else if (pdgcode == 1000020040) {
          histograms.fill(HIST("GenAftTrkSel/hGenTrkPrimAftTrkSel_truepid_al"), particle.pt(), y_gen);
        }

        histograms.fill(HIST("GenAftTrkSel/hGenTrkPrimAftTrkSel"), particle.pt(), y_gen);

        // reconstructed
        if (pdgcode == 11) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_el"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 211) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_pi"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 321) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_ka"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 2212) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_pr"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 1000010020) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_de"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 1000020030) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_he"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 1000010030) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_tr"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 1000020040) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkPrimAftTrkSel_truepid_al"), track.pt(), y_rec, track.dcaXY());
        }

        histograms.fill(HIST("RecAftTrkSel/hRecTrkPrimAftTrkSel"), track.pt(), y_rec, track.dcaXY());
      } else {
        // reconstructed
        if (pdgcode == 11) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_el"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 211) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_pi"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 321) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_ka"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 2212) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_pr"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 1000010020) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_de"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 1000020030) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_he"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 1000010030) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_tr"), track.pt(), y_rec, track.dcaXY());
        } else if (pdgcode == 1000020040) {
          histograms.fill(HIST("RecAftTrkSel/hRecTrkSecAftTrkSel_truepid_al"), track.pt(), y_rec, track.dcaXY());
        }
        histograms.fill(HIST("RecAftTrkSel/hRecTrkSecAftTrkSel"), track.pt(), y_rec, track.dcaXY());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TrackCheckTaskEvSel>(cfgc, TaskName{"track-histos-evsel"}),
    adaptAnalysisTask<TrackCheckTaskEvSelTrackSel>(cfgc, TaskName{"track-histos-evsel-trksel"})};
}
