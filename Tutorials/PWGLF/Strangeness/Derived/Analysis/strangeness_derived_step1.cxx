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
/// \brief Step4 of the Strangeness tutorial
/// \author Romain Schotter
/// based on the original codes from:
/// \author Nepeivoda Roman (roman.nepeivoda@cern.ch)
/// \author Chiara De Martin (chiara.de.martin@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/EventSelection.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 0
// Starting point: loop over all cascades and fill invariant mass histogram
// STEP 1
// Apply selections on topological variables of Cascades

struct strangeness_derived_tutorial {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rXi{"xi", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rOmega{"omega", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurable parameters for cascade selection
  Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.98, "Casc CosPA"};
  Configurable<double> cascadesetting_v0cospa{"cascadesetting_v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "DCA cascade daughters"};
  Configurable<float> cascadesetting_dcav0dau{"cascadesetting_dcav0dau", 1.0, "DCA v0 daughters"};
  Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.06, "DCA bachelor to PV"};
  Configurable<float> cascadesetting_dcapostopv{"cascadesetting_dcapostopv", 0.06, "DCA positive to PV"};
  Configurable<float> cascadesetting_dcanegtopv{"cascadesetting_dcanegtopv", 0.06, "DCA negative to PV"};
  Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "minimum V0 DCA to PV"};
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.5, "cascradius"};
  Configurable<float> cascadesetting_v0radius{"cascadesetting_v0radius", 1.2, "v0radius"};
  Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "v0 mass window"};
  Configurable<float> cascadesetting_competingmassrej{"cascadesetting_competingmassrej", 0.008, "Competing mass rejection"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec XiMassAxis = {100, 1.28f, 1.36f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec OmegaMassAxis = {100, 1.63f, 1.7f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    // Xi/Omega reconstruction
    rXi.add("hMassXi", "hMassXi", {HistType::kTH1F, {XiMassAxis}});
    rXi.add("hMassXiSelected", "hMassXiSelected", {HistType::kTH1F, {XiMassAxis}});

    rOmega.add("hMassOmega", "hMassOmega", {HistType::kTH1F, {OmegaMassAxis}});
    rOmega.add("hMassOmegaSelected", "hMassOmegaSelected", {HistType::kTH1F, {OmegaMassAxis}});

    // Xi/Omega topological cuts
    rXi.add("hCascDCAV0Daughters", "hCascDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}});
    rXi.add("hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}});

    rOmega.add("hCascDCAV0Daughters", "hCascDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}});
    rOmega.add("hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  // Filters on Cascades
  // Cannot filter on dynamic columns
  Filter preFilterCascades = (aod::cascdata::dcaV0daughters < cascadesetting_dcav0dau &&
                              nabs(aod::cascdata::dcapostopv) > cascadesetting_dcapostopv &&
                              nabs(aod::cascdata::dcanegtopv) > cascadesetting_dcanegtopv &&
                              nabs(aod::cascdata::dcabachtopv) > cascadesetting_dcabachtopv &&
                              aod::cascdata::dcacascdaughters < cascadesetting_dcacascdau);

  void process(soa::Filtered<soa::Join<aod::StraCollisions, aod::StraEvSels>>::iterator const& collision,
               soa::Filtered<soa::Join<aod::CascCores, aod::CascExtras>> const& Cascades)
  {
    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());

    // Cascades
    for (const auto& casc : Cascades) {
      rXi.fill(HIST("hMassXi"), casc.mXi());
      rOmega.fill(HIST("hMassOmega"), casc.mOmega());

      // Cut on dynamic columns
      if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadesetting_cospa)
        continue;
      if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadesetting_v0cospa)
        continue;
      if (std::abs(casc.mLambda() - o2::constants::physics::MassLambda) > cascadesetting_v0masswindow)
        continue;
      if (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) < cascadesetting_mindcav0topv)
        continue;
      if (casc.cascradius() < cascadesetting_cascradius)
        continue;
      if (casc.v0radius() < cascadesetting_v0radius)
        continue;

      // Fill histograms! (if possible)
      rXi.fill(HIST("hMassXiSelected"), casc.mXi());

      rXi.fill(HIST("hCascDCAV0Daughters"), casc.dcaV0daughters());
      rXi.fill(HIST("hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));

      if (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > cascadesetting_competingmassrej) { // competing mass rejection, only in case of Omega
        rOmega.fill(HIST("hMassOmegaSelected"), casc.mOmega());

        rOmega.fill(HIST("hCascDCAV0Daughters"), casc.dcaV0daughters());
        rOmega.fill(HIST("hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangeness_derived_tutorial>(cfgc)};
}
