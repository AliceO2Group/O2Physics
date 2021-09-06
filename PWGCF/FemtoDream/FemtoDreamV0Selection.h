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

/// \file FemtoDreamV0Selection.h
/// \brief Definition of the FemtoDreamV0Selection
/// \author Valentina Mantovani Sarti, TU München valentina.mantovani-sarti@tum.de
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#ifndef ANALYSIS_TASKS_PWGCF_FEMTODREAM_FEMTODREAMV0SELECTION_H_
#define ANALYSIS_TASKS_PWGCF_FEMTODREAM_FEMTODREAMV0SELECTION_H_

#include "FemtoDreamObjectSelection.h"
#include "FemtoDreamTrackSelection.h"
#include "FemtoDreamSelection.h"

#include "ReconstructionDataFormats/PID.h"
#include "Common/Core/RecoDecay.h"
#include "Framework/HistogramRegistry.h"

using namespace o2::framework;

namespace o2::analysis::femtoDream
{
namespace femtoDreamV0Selection
{
/// The different selections this task is capable of doing
enum V0Sel { kpTV0Min, //!Min. p_T (GeV/c)
             kpTV0Max, //!Max. p_T (GeV/c)
             kDCAV0DaughMax,
             kCPAV0Min,
             kTranRadV0Min,
             kTranRadV0Max,
             kDecVtxMax };
enum ChildTrackType { kPosTrack,
                      kNegTrack };
} // namespace femtoDreamV0Selection

/// \class FemtoDreamV0Selection
/// \brief Cut class to contain and execute all cuts applied to V0s
class FemtoDreamV0Selection : public FemtoDreamObjectSelection<float, femtoDreamV0Selection::V0Sel>
{
 public:
  FemtoDreamV0Selection() : nPtV0MinSel(0),
                            nPtV0MaxSel(0),
                            nDCAV0DaughMax(0),
                            nCPAV0Min(0),
                            nTranRadV0Min(0),
                            nTranRadV0Max(0),
                            nDecVtxMax(0),
                            pTV0Min(9999999.),
                            pTV0Max(-9999999.),
                            DCAV0DaughMax(-9999999.),
                            CPAV0Min(9999999.),
                            TranRadV0Min(9999999.),
                            TranRadV0Max(-9999999.),
                            DecVtxMax(-9999999.){};
  /// Initializes histograms for the task
  template <o2::aod::femtodreamparticle::ParticleType part, o2::aod::femtodreamparticle::ParticleType daugh, typename cutContainerType>
  void init(HistogramRegistry* registry);

  template <typename C, typename V, typename T>
  bool isSelectedMinimal(C const& col, V const& v0, T const& posTrack, T const& negTrack);

  /// \todo for the moment the PID of the tracks is factored out into a separate field, hence 5 values in total
  template <typename cutContainerType, typename C, typename V, typename T>
  std::array<cutContainerType, 5> getCutContainer(C const& col, V const& v0, T const& posTrack, T const& negTrack);

  template <o2::aod::femtodreamparticle::ParticleType part, o2::aod::femtodreamparticle::ParticleType daugh, typename C, typename V, typename T>
  void fillQA(C const& col, V const& v0, T const& posTrack, T const& negTrack);

  template <typename T1, typename T2>
  void setChildCuts(femtoDreamV0Selection::ChildTrackType child, T1 selVal, T2 selVar, femtoDreamSelection::SelectionType selType)
  {
    if (child == femtoDreamV0Selection::kPosTrack) {
      PosDaughTrack.setSelection(selVal, selVar, selType);
    } else if (child == femtoDreamV0Selection::kNegTrack) {
      NegDaughTrack.setSelection(selVal, selVar, selType);
    }
  }
  template <typename T>
  void setChildPIDSpecies(femtoDreamV0Selection::ChildTrackType child, T& pids)
  {
    if (child == femtoDreamV0Selection::kPosTrack) {
      PosDaughTrack.setPIDSpecies(pids);
    } else if (child == femtoDreamV0Selection::kNegTrack) {
      NegDaughTrack.setPIDSpecies(pids);
    }
  }

 private:
  int nPtV0MinSel;
  int nPtV0MaxSel;
  int nDCAV0DaughMax;
  int nCPAV0Min;
  int nTranRadV0Min;
  int nTranRadV0Max;
  int nDecVtxMax;
  float pTV0Min;
  float pTV0Max;
  float DCAV0DaughMax;
  float CPAV0Min;
  float TranRadV0Min;
  float TranRadV0Max;
  float DecVtxMax;
  FemtoDreamTrackSelection PosDaughTrack;
  FemtoDreamTrackSelection NegDaughTrack;

}; // namespace femtoDream

template <o2::aod::femtodreamparticle::ParticleType part, o2::aod::femtodreamparticle::ParticleType daugh, typename cutContainerType>
void FemtoDreamV0Selection::init(HistogramRegistry* registry)
{
  if (registry) {
    mHistogramRegistry = registry;
    fillSelectionHistogram<part>();
    fillSelectionHistogram<daugh>();

    /// \todo this should be an automatic check in the parent class, and the return type should be templated
    int nSelections = getNSelections();
    if (8 * sizeof(cutContainerType) < nSelections) {
      LOG(FATAL) << "FemtoDreamV0Cuts: Number of selections to large for your container - quitting!";
    }
    std::string folderName = static_cast<std::string>(o2::aod::femtodreamparticle::ParticleTypeName[part]);
    /// \todo initialize histograms for children tracks of v0s
    mHistogramRegistry->add((folderName + "/pThist").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{1000, 0, 10}});
    mHistogramRegistry->add((folderName + "/etahist").c_str(), "; #eta; Entries", kTH1F, {{1000, -1, 1}});
    mHistogramRegistry->add((folderName + "/phihist").c_str(), "; #phi; Entries", kTH1F, {{1000, 0, 2. * M_PI}});
    mHistogramRegistry->add((folderName + "/dcaDauToVtx").c_str(), "; DCADaug_{Vtx} (cm); Entries", kTH1F, {{1000, 0, 10}});
    mHistogramRegistry->add((folderName + "/transRadius").c_str(), "; #it{r}_{xy} (cm); Entries", kTH1F, {{1500, 0, 150}});
    mHistogramRegistry->add((folderName + "/decayVtxXPV").c_str(), "; #it{iVtx}_{x} (cm); Entries", kTH1F, {{2000, 0, 200}});
    mHistogramRegistry->add((folderName + "/decayVtxYPV").c_str(), "; #it{iVtx}_{y} (cm)); Entries", kTH1F, {{2000, 0, 200}});
    mHistogramRegistry->add((folderName + "/decayVtxZPV").c_str(), "; #it{iVtx}_{z} (cm); Entries", kTH1F, {{2000, 0, 200}});
    mHistogramRegistry->add((folderName + "/cpa").c_str(), "; #it{cos(#alpha)}; Entries", kTH1F, {{1000, 0.9, 1.}});
    mHistogramRegistry->add((folderName + "/cpapTBins").c_str(), "; #it{p}_{T} (GeV/#it{c}); #it{cos(#alpha)}", kTH2F, {{8, 0.3, 4.3}, {1000, 0.9, 1.}});

    PosDaughTrack.init<aod::femtodreamparticle::ParticleType::kV0Child, aod::femtodreamparticle::cutContainerType>(mHistogramRegistry, "Pos");
    NegDaughTrack.init<aod::femtodreamparticle::ParticleType::kV0Child, aod::femtodreamparticle::cutContainerType>(mHistogramRegistry, "Neg");
  }
  /// check whether the most open cuts are fulfilled - most of this should have already be done by the filters
  nPtV0MinSel = getNSelections(femtoDreamV0Selection::kpTV0Min);
  nPtV0MaxSel = getNSelections(femtoDreamV0Selection::kpTV0Max);
  nDCAV0DaughMax = getNSelections(femtoDreamV0Selection::kDCAV0DaughMax);
  nCPAV0Min = getNSelections(femtoDreamV0Selection::kCPAV0Min);
  nTranRadV0Min = getNSelections(femtoDreamV0Selection::kTranRadV0Min);
  nTranRadV0Max = getNSelections(femtoDreamV0Selection::kTranRadV0Max);
  nDecVtxMax = getNSelections(femtoDreamV0Selection::kDecVtxMax);

  pTV0Min = getMinimalSelection(femtoDreamV0Selection::kpTV0Min, femtoDreamSelection::kLowerLimit);
  pTV0Max = getMinimalSelection(femtoDreamV0Selection::kpTV0Max, femtoDreamSelection::kUpperLimit);
  DCAV0DaughMax = getMinimalSelection(femtoDreamV0Selection::kDCAV0DaughMax, femtoDreamSelection::kUpperLimit);
  CPAV0Min = getMinimalSelection(femtoDreamV0Selection::kCPAV0Min, femtoDreamSelection::kLowerLimit);
  TranRadV0Min = getMinimalSelection(femtoDreamV0Selection::kTranRadV0Min, femtoDreamSelection::kLowerLimit);
  TranRadV0Max = getMinimalSelection(femtoDreamV0Selection::kTranRadV0Max, femtoDreamSelection::kUpperLimit);
  DecVtxMax = getMinimalSelection(femtoDreamV0Selection::kDecVtxMax, femtoDreamSelection::kAbsUpperLimit);
}

template <typename C, typename V, typename T>
bool FemtoDreamV0Selection::isSelectedMinimal(C const& col, V const& v0, T const& posTrack, T const& negTrack)
{
  const auto signPos = posTrack.sign();
  const auto signNeg = negTrack.sign();
  if (signPos < 0 || signNeg > 0) {
    printf("-Something wrong in isSelectedMinimal--\n");
    printf("ERROR - Wrong sign for V0 daughters\n");
  }
  const float pT = v0.pt();
  const std::vector<float> decVtx = {v0.x(), v0.y(), v0.z()};
  const float tranRad = v0.v0radius();
  const float dcaDaughv0 = v0.dcaV0daughters();
  const float cpav0 = v0.v0cosPA(col.posX(), col.posY(), col.posZ());

  if (nPtV0MinSel > 0 && pT < pTV0Min) {
    return false;
  }
  if (nPtV0MaxSel > 0 && pT > pTV0Max) {
    return false;
  }
  if (nDCAV0DaughMax > 0 && dcaDaughv0 > DCAV0DaughMax) {
    return false;
  }
  if (nCPAV0Min > 0 && cpav0 < CPAV0Min) {
    return false;
  }
  if (nTranRadV0Min > 0 && tranRad < TranRadV0Min) {
    return false;
  }
  if (nTranRadV0Max > 0 && tranRad > TranRadV0Max) {
    return false;
  }
  for (int i = 0; i < decVtx.size(); i++) {
    if (nDecVtxMax > 0 && decVtx.at(i) > DecVtxMax) {
      return false;
    }
  }
  const auto dcaXYpos = posTrack.dcaXY();
  const auto dcaZpos = posTrack.dcaZ();
  const auto dcapos = std::sqrt(pow(dcaXYpos, 2.) + pow(dcaZpos, 2.));
  if (!PosDaughTrack.isSelectedMinimal(posTrack)) {
    return false;
  }
  if (!NegDaughTrack.isSelectedMinimal(negTrack)) {
    return false;
  }
  return true;
}

/// the CosPA of V0 needs as argument the posXYZ of collisions vertex so we need to pass the collsion as well
template <typename cutContainerType, typename C, typename V, typename T>
std::array<cutContainerType, 5> FemtoDreamV0Selection::getCutContainer(C const& col, V const& v0, T const& posTrack, T const& negTrack)
{
  auto outputPosTrack = PosDaughTrack.getCutContainer<cutContainerType>(posTrack);
  auto outputNegTrack = NegDaughTrack.getCutContainer<cutContainerType>(negTrack);
  cutContainerType output = 0;
  size_t counter = 0;

  const auto pT = v0.pt();
  const auto tranRad = v0.v0radius();
  const auto dcaDaughv0 = v0.dcaV0daughters();
  const auto cpav0 = v0.v0cosPA(col.posX(), col.posY(), col.posZ());
  const std::vector<float> decVtx = {v0.x(), v0.y(), v0.z()};

  float observable;
  for (auto& sel : mSelections) {
    const auto selVariable = sel.getSelectionVariable();
    if (selVariable == femtoDreamV0Selection::kDecVtxMax) {
      for (size_t i = 0; i < decVtx.size(); ++i) {
        auto decVtxValue = decVtx.at(i);
        sel.checkSelectionSetBit(decVtxValue, output, counter);
      }
    } else {
      switch (selVariable) {
        case (femtoDreamV0Selection::kpTV0Min):
        case (femtoDreamV0Selection::kpTV0Max):
          observable = pT;
          break;
        case (femtoDreamV0Selection::kDCAV0DaughMax):
          observable = dcaDaughv0;
          break;
        case (femtoDreamV0Selection::kCPAV0Min):
          observable = cpav0;
          break;
        case (femtoDreamV0Selection::kTranRadV0Min):
        case (femtoDreamV0Selection::kTranRadV0Max):
          observable = tranRad;
          break;
        case (femtoDreamV0Selection::kDecVtxMax):
          break;
      }
      sel.checkSelectionSetBit(observable, output, counter);
    }
  }
  return {output, outputPosTrack.at(0), outputPosTrack.at(1), outputNegTrack.at(0), outputNegTrack.at(1)};
}

template <o2::aod::femtodreamparticle::ParticleType part, o2::aod::femtodreamparticle::ParticleType daugh, typename C, typename V, typename T>
void FemtoDreamV0Selection::fillQA(C const& col, V const& v0, T const& posTrack, T const& negTrack)
{
  if (mHistogramRegistry) {
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/pThist"), v0.pt());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/etahist"), v0.eta());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/phihist"), v0.phi());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/dcaDauToVtx"), v0.dcaV0daughters());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/transRadius"), v0.v0radius());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/decayVtxXPV"), v0.x());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/decayVtxYPV"), v0.y());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/decayVtxZPV"), v0.z());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/cpa"), v0.v0cosPA(col.posX(), col.posY(), col.posZ()));
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/cpapTBins"), v0.pt(), v0.v0cosPA(col.posX(), col.posY(), col.posZ()));
  }
  /// \TODO to add possibility to write to Daughters folders
  // PosDaughTrack.fillQA<aod::femtodreamparticle::ParticleType::kV0Child>(posTrack, "Pos");
  // NegDaughTrack.fillQA<aod::femtodreamparticle::ParticleType::kV0Child>(negTrack, "Neg");
}

} // namespace o2::analysis::femtoDream

#endif /* ANALYSIS_TASKS_PWGCF_FEMTODREAM_FEMTODREAMV0SELECTION_H_ */
