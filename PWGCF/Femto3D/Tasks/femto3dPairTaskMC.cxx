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
//
/// \brief Femto3D pair mixing task
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 31 May 2023

#include <vector>
#include <TParameter.h>
#include <TH1F.h>

#include "PWGCF/Femto3D/Core/femto3dPairTask.h"
#include "PWGCF/Femto3D/DataModel/singletrackselector.h"
#include "TLorentzVector.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Expressions.h"
#include "Framework/StaticFor.h"
#include "MathUtils/Utils.h"
#include "Common/DataModel/Multiplicity.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct FemtoCorrelationsMC {
  // using allinfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPr, aod::TOFSignal, aod::TracksDCA, aod::pidTOFFullPr, aod::pidTOFbeta, aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullDe, aod::pidTPCFullDe>; // aod::pidTPCPr
  /// Construct a registry object with direct declaration
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> _min_P{"min_P", 0.0, "lower mometum limit"};
  Configurable<float> _max_P{"max_P", 100.0, "upper mometum limit"};
  Configurable<float> _eta{"eta", 100.0, "abs eta value limit"};
  Configurable<float> _dcaXY{"dcaXY", 10.0, "abs dcaXY value limit"};
  Configurable<float> _dcaZ{"dcaZ", 10.0, "abs dcaZ value limit"};
  Configurable<int16_t> _tpcNClsFound{"minTpcNClsFound", 0, "minimum allowed number of TPC clasters"};
  Configurable<float> _tpcChi2NCl{"tpcChi2NCl", 100.0, "upper limit for chi2 value of a fit over TPC clasters"};
  Configurable<float> _tpcCrossedRowsOverFindableCls{"tpcCrossedRowsOverFindableCls", 0, "lower limit of TPC CrossedRows/FindableCls value"};
  Configurable<int> _tpcNClsShared{"maxTpcNClsShared", 100, "maximum allowed number of TPC shared clasters"};
  Configurable<int> _itsNCls{"minItsNCls", 0, "minimum allowed number of ITS clasters"};
  Configurable<float> _itsChi2NCl{"itsChi2NCl", 100.0, "upper limit for chi2 value of a fit over ITS clasters"};
  Configurable<float> _vertexZ{"VertexZ", 10.0, "abs vertexZ value limit"};

  Configurable<int> _sign_1{"sign_1", 1, "sign of the first particle in a pair"};
  Configurable<int> _particlePDG_1{"particlePDG_1", 2212, "PDG code of the first particle in a pair to perform PID for (only pion, kaon, proton and deurton are supported now)"};
  Configurable<std::vector<float>> _tpcNSigma_1{"tpcNSigma_1", std::vector<float>{-3.0f, 3.0f}, "first particle PID: Nsigma range in TPC before the TOF is used"};
  Configurable<float> _PIDtrshld_1{"PIDtrshld_1", 10.0, "first particle PID: value of momentum from which the PID is done with TOF (before that only TPC is used)"};
  Configurable<std::vector<float>> _tofNSigma_1{"tofNSigma_1", std::vector<float>{-3.0f, 3.0f}, "first particle PID: Nsigma range in TOF"};
  Configurable<float> _tpcNSigmaResidual_1{"tpcNSigmaResidual_1", 5, "first particle PID: residual TPC Nsigma cut (abs. value) to use with the TOF"};

  Configurable<int> _sign_2{"sign_2", 1, "sign of the second particle in a pair"};
  Configurable<int> _particlePDG_2{"particlePDG_2", 2212, "PDG code of the second particle in a pair to perform PID for (only pion, kaon, proton and deurton are supported now)"};
  Configurable<std::vector<float>> _tpcNSigma_2{"tpcNSigma_2", std::vector<float>{-3.0f, 3.0f}, "second particle PID: Nsigma range in TPC before the TOF is used"};
  Configurable<float> _PIDtrshld_2{"PIDtrshld_2", 10.0, "second particle PID: value of momentum from which the PID is done with TOF (before that only TPC is used)"};
  Configurable<std::vector<float>> _tofNSigma_2{"tofNSigma_2", std::vector<float>{-3.0f, 3.0f}, "second particle PID: Nsigma range in TOF"};
  Configurable<float> _tpcNSigmaResidual_2{"tpcNSigmaResidual_2", 5, "second particle PID: residual TPC Nsigma cut (abs. value) to use with the TOF"};

  Configurable<int> _particlePDGtoReject{"particlePDGtoRejectFromSecond", 0, "applied only if the particles are non-identical and only to the second particle in the pair!!!"};
  Configurable<std::vector<float>> _rejectWithinNsigmaTOF{"rejectWithinNsigmaTOF", std::vector<float>{-0.0f, 0.0f}, "TOF rejection Nsigma range for the particle specified with PDG to be rejected"};

  Configurable<float> _radiusTPC{"radiusTPC", 1.2, "TPC radius to calculate phi_star for"};

  ConfigurableAxis CFkStarBinning{"CFkStarBinning", {500, 0.005, 5.005}, "k* binning of the res. matrix (Nbins, lowlimit, uplimit)"};

  Configurable<int> _vertexNbinsToMix{"vertexNbinsToMix", 10, "Number of vertexZ bins for the mixing"};
  Configurable<int> _multNsubBins{"multSubBins", 10, "number of sub-bins to perform the mixing within"};

  bool IsIdentical;

  std::pair<int, std::vector<float>> TPCcuts_1;
  std::pair<int, std::vector<float>> TOFcuts_1;

  std::pair<int, std::vector<float>> TPCcuts_2;
  std::pair<int, std::vector<float>> TOFcuts_2;

  using FilteredCollisions = aod::SingleCollSels;
  using FilteredTracks = soa::Join<aod::SingleTrackSels, aod::SingleTrkMCs>;

  typedef std::shared_ptr<soa::Filtered<FilteredTracks>::iterator> trkType;
  typedef std::shared_ptr<soa::Filtered<FilteredCollisions>::iterator> colType;

  std::map<int64_t, std::vector<trkType>> selectedtracks_1;
  std::map<int64_t, std::vector<trkType>> selectedtracks_2;
  std::map<std::pair<int, float>, std::vector<colType>> mixbins;

  std::unique_ptr<o2::aod::singletrackselector::FemtoPair<trkType>> Pair = std::make_unique<o2::aod::singletrackselector::FemtoPair<trkType>>();

  Filter pFilter = o2::aod::singletrackselector::p > _min_P&& o2::aod::singletrackselector::p < _max_P;
  Filter etaFilter = nabs(o2::aod::singletrackselector::eta) < _eta;

  Filter tpcTrkFilter = o2::aod::singletrackselector::tpcNClsFound >= _tpcNClsFound &&
                        o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedTpcChi2NCl) < _tpcChi2NCl &&
                        o2::aod::singletrackselector::unPack<singletrackselector::binning::rowsOverFindable>(o2::aod::singletrackselector::storedTpcCrossedRowsOverFindableCls) > _tpcCrossedRowsOverFindableCls;

  Filter itsTrkFilter = o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedItsChi2NCl) < _itsChi2NCl;

  Filter vertexFilter = nabs(o2::aod::singletrackselector::posZ) < _vertexZ;

  void init(o2::framework::InitContext&)
  {

    IsIdentical = (_sign_1 * _particlePDG_1 == _sign_2 * _particlePDG_2);

    Pair->SetIdentical(IsIdentical);
    Pair->SetPDG1(_particlePDG_1);
    Pair->SetPDG2(_particlePDG_2);

    TPCcuts_1 = std::make_pair(_particlePDG_1, _tpcNSigma_1);
    TOFcuts_1 = std::make_pair(_particlePDG_1, _tofNSigma_1);
    TPCcuts_2 = std::make_pair(_particlePDG_2, _tpcNSigma_2);
    TOFcuts_2 = std::make_pair(_particlePDG_2, _tofNSigma_2);

    registry.add("FirstParticle/dcaxyz_vs_pt_primary", "dcaxyz_vs_pt_primary", kTH3F, {{100, 0., 5., "pt"}, {250, -1., 1., "DCA_XY(pt) primary"}, {250, -1., 1., "DCA_Z(pt) primary"}});
    registry.add("FirstParticle/dcaxyz_vs_pt_weakdecay", "dcaxyz_vs_pt_weakdecay", kTH3F, {{100, 0., 5., "pt"}, {250, -1., 1., "DCA_XY(pt) weakdecay"}, {250, -1., 1., "DCA_Z(pt) weakdecay"}});
    registry.add("FirstParticle/dcaxyz_vs_pt_material", "dcaxyz_vs_pt_material", kTH3F, {{100, 0., 5., "pt"}, {250, -1., 1., "DCA_XY(pt) material"}, {250, -1., 1., "DCA_Z(pt) material"}});

    registry.add("FirstParticle/pSpectraEl", "pSpectraEl", kTH1F, {{100, 0., 5., "p"}});
    registry.add("FirstParticle/pSpectraMu", "pSpectraMu", kTH1F, {{100, 0., 5., "p"}});
    registry.add("FirstParticle/pSpectraPi", "pSpectraPi", kTH1F, {{100, 0., 5., "p"}});
    registry.add("FirstParticle/pSpectraKa", "pSpectraKa", kTH1F, {{100, 0., 5., "p"}});
    registry.add("FirstParticle/pSpectraPr", "pSpectraPr", kTH1F, {{100, 0., 5., "p"}});
    registry.add("FirstParticle/pSpectraDe", "pSpectraDe", kTH1F, {{100, 0., 5., "p"}});
    registry.add("FirstParticle/pSpectraAll", "pSpectrAll", kTH1F, {{100, 0., 5., "p"}});

    if (!IsIdentical) {
      registry.add("SecondParticle/dcaxyz_vs_pt_primary", "dcaxyz_vs_pt_primary", kTH3F, {{100, 0., 5., "pt"}, {200, -1., 1., "DCA_XY(pt) primary"}, {200, -1., 1., "DCA_Z(pt) primary"}});
      registry.add("SecondParticle/dcaxyz_vs_pt_weakdecay", "dcaxyz_vs_pt_weakdecay", kTH3F, {{100, 0., 5., "pt"}, {200, -1., 1., "DCA_XY(pt) weakdecay"}, {200, -1., 1., "DCA_Z(pt) weakdecay"}});
      registry.add("SecondParticle/dcaxyz_vs_pt_material", "dcaxyz_vs_pt_material", kTH3F, {{100, 0., 5., "pt"}, {200, -1., 1., "DCA_XY(pt) material"}, {200, -1., 1., "DCA_Z(pt) material"}});

      registry.add("SecondParticle/pSpectraEl", "pSpectraEl", kTH1F, {{100, 0., 5., "p"}});
      registry.add("SecondParticle/pSpectraMu", "pSpectraMu", kTH1F, {{100, 0., 5., "p"}});
      registry.add("SecondParticle/pSpectraPi", "pSpectraPi", kTH1F, {{100, 0., 5., "p"}});
      registry.add("SecondParticle/pSpectraKa", "pSpectraKa", kTH1F, {{100, 0., 5., "p"}});
      registry.add("SecondParticle/pSpectraPr", "pSpectraPr", kTH1F, {{100, 0., 5., "p"}});
      registry.add("SecondParticle/pSpectraDe", "pSpectraDe", kTH1F, {{100, 0., 5., "p"}});
      registry.add("SecondParticle/pSpectraAll", "pSpectrAll", kTH1F, {{100, 0., 5., "p"}});
    }

    registry.add("ResolutionMatrix", "ResolutionMatrix_rec(gen)", kTH2F, {{CFkStarBinning, "k*_gen (GeV/c)"}, {CFkStarBinning, "k*_rec (GeV/c)"}});
    registry.add("DoubleTrackEffects", "DoubleTrackEffects_deta(dphi*)", kTH2F, {{200, -M_PI, M_PI, "dphi*"}, {200, -0.5, 0.5, "deta"}});
  }

  template <typename Type>
  void fillEtaPhi(Type const& tracks)
  {                                              // template for particles from the same collision identical
    for (int ii = 0; ii < tracks.size(); ii++) { // nested loop for all the combinations
      for (int iii = ii + 1; iii < tracks.size(); iii++) {

        Pair->SetPair(tracks[ii], tracks[iii]);

        registry.fill(HIST("DoubleTrackEffects"), Pair->GetPhiStarDiff(_radiusTPC), Pair->GetEtaDiff());
        Pair->ResetPair();
      }
    }
  }

  template <typename Type>
  void fillEtaPhi(Type const& tracks1, Type const& tracks2)
  { // template for particles from the same collision non-identical
    for (auto ii : tracks1) {
      for (auto iii : tracks2) {

        Pair->SetPair(ii, iii);

        registry.fill(HIST("DoubleTrackEffects"), Pair->GetPhiStarDiff(_radiusTPC), Pair->GetEtaDiff());
        Pair->ResetPair();
      }
    }
  }

  template <typename Type>
  void fillResMatrix(Type const& tracks1, Type const& tracks2)
  { // template for ME
    for (auto ii : tracks1) {
      for (auto iii : tracks2) {
        Pair->SetPair(ii, iii);

        TLorentzVector first4momentumGen;
        first4momentumGen.SetPtEtaPhiM(ii->pt_MC(), ii->eta_MC(), ii->phi_MC(), particle_mass(_particlePDG_1));
        TLorentzVector second4momentumGen;
        second4momentumGen.SetPtEtaPhiM(iii->pt_MC(), iii->eta_MC(), iii->phi_MC(), particle_mass(_particlePDG_2));

        registry.fill(HIST("ResolutionMatrix"), o2::aod::singletrackselector::GetKstarFrom4vectors(first4momentumGen, second4momentumGen, IsIdentical), Pair->GetKstar());
        Pair->ResetPair();
      }
    }
  }

  void process(soa::Filtered<FilteredCollisions> const& collisions, soa::Filtered<FilteredTracks> const& tracks)
  {
    if (_particlePDG_1 == 0 || _particlePDG_2 == 0)
      LOGF(fatal, "One of passed PDG is 0!!!");

    int trackPDG, trackOrigin;

    for (auto track : tracks) {
      if (abs(track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().posZ()) > _vertexZ)
        continue;
      if (track.tpcNClsShared() > _tpcNClsShared || track.itsNCls() < _itsNCls)
        continue;

      if (track.sign() == _sign_1 && (track.p() < _PIDtrshld_1 ? o2::aod::singletrackselector::TPCselection(track, TPCcuts_1) : o2::aod::singletrackselector::TOFselection(track, TOFcuts_1, _tpcNSigmaResidual_1))) {
        trackOrigin = track.origin();
        switch (trackOrigin) {
          case 0:
            registry.fill(HIST("FirstParticle/dcaxyz_vs_pt_primary"), track.pt(), track.dcaXY(), track.dcaZ());
            break;
          case 1:
            registry.fill(HIST("FirstParticle/dcaxyz_vs_pt_weakdecay"), track.pt(), track.dcaXY(), track.dcaZ());
            break;
          case 2:
            registry.fill(HIST("FirstParticle/dcaxyz_vs_pt_material"), track.pt(), track.dcaXY(), track.dcaZ());
            break;
        }
        if (abs(track.dcaXY()) > _dcaXY || abs(track.dcaZ()) > _dcaZ)
          continue;

        selectedtracks_1[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track)); // filling the map: eventID <-> selected particles1

        trackPDG = abs(track.pdgCode());

        registry.fill(HIST("FirstParticle/pSpectraAll"), track.p_MC());

        switch (trackPDG) {
          case 11:
            registry.fill(HIST("FirstParticle/pSpectraEl"), track.p_MC());
            break;
          case 13:
            registry.fill(HIST("FirstParticle/pSpectraMu"), track.p_MC());
            break;
          case 211:
            registry.fill(HIST("FirstParticle/pSpectraPi"), track.p_MC());
            break;
          case 321:
            registry.fill(HIST("FirstParticle/pSpectraKa"), track.p_MC());
            break;
          case 2212:
            registry.fill(HIST("FirstParticle/pSpectraPr"), track.p_MC());
            break;
          case 1000010020:
            registry.fill(HIST("FirstParticle/pSpectraDe"), track.p_MC());
            break;
        }
      }

      if (IsIdentical) {
        continue;
      } else if (track.sign() != _sign_2 && !TOFselection(track, std::make_pair(_particlePDGtoReject, _rejectWithinNsigmaTOF)) && (track.p() < _PIDtrshld_2 ? o2::aod::singletrackselector::TPCselection(track, TPCcuts_2) : o2::aod::singletrackselector::TOFselection(track, TOFcuts_2, _tpcNSigmaResidual_2))) { // filling the map: eventID <-> selected particles2 if (see condition above ^)
        trackOrigin = track.origin();
        switch (trackOrigin) {
          case 0:
            registry.fill(HIST("SecondParticle/dcaxyz_vs_pt_primary"), track.pt(), track.dcaXY(), track.dcaZ());
            break;
          case 1:
            registry.fill(HIST("SecondParticle/dcaxyz_vs_pt_weakdecay"), track.pt(), track.dcaXY(), track.dcaZ());
            break;
          case 2:
            registry.fill(HIST("SecondParticle/dcaxyz_vs_pt_material"), track.pt(), track.dcaXY(), track.dcaZ());
            break;
        }
        if (abs(track.dcaXY()) > _dcaXY || abs(track.dcaZ()) > _dcaZ)
          continue;

        selectedtracks_2[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track)); // filling the map: eventID <-> selected particles2

        trackPDG = abs(track.pdgCode());

        registry.fill(HIST("SecondParticle/pSpectraAll"), track.p_MC());

        switch (trackPDG) {
          case 11:
            registry.fill(HIST("SecondParticle/pSpectraEl"), track.p_MC());
            break;
          case 13:
            registry.fill(HIST("SecondParticle/pSpectraMu"), track.p_MC());
            break;
          case 211:
            registry.fill(HIST("SecondParticle/pSpectraPi"), track.p_MC());
            break;
          case 321:
            registry.fill(HIST("SecondParticle/pSpectraKa"), track.p_MC());
            break;
          case 2212:
            registry.fill(HIST("SecondParticle/pSpectraPr"), track.p_MC());
            break;
          case 1000010020:
            registry.fill(HIST("SecondParticle/pSpectraDe"), track.p_MC());
            break;
        }
      }
    }

    for (auto collision : collisions) {
      if (selectedtracks_1.find(collision.globalIndex()) == selectedtracks_1.end()) {
        if (IsIdentical)
          continue;
        else if (selectedtracks_2.find(collision.globalIndex()) == selectedtracks_2.end())
          continue;
      }
      int vertexBinToMix = std::floor((collision.posZ() + _vertexZ) / (2 * _vertexZ / _vertexNbinsToMix));
      int centBinToMix = std::floor(collision.multPerc() / (100.0 / _multNsubBins));

      mixbins[std::pair<int, int>{vertexBinToMix, centBinToMix}].push_back(std::make_shared<decltype(collision)>(collision));
    }

    //====================================== filling deta(dphi*) & res. matrix starts here ======================================

    if (IsIdentical) { //====================================== identical ======================================

      for (auto i = mixbins.begin(); i != mixbins.end(); i++) { // iterating over all vertex&mult bins

        for (int indx1 = 0; indx1 < (i->second).size(); indx1++) { // iterating over all selected collisions with selected tracks

          auto col1 = (i->second)[indx1];

          Pair->SetMagField1(col1->magField());
          Pair->SetMagField2(col1->magField());

          fillEtaPhi(selectedtracks_1[col1->index()]); // filling deta(dphi*) -- SE identical

          for (int indx2 = indx1 + 1; indx2 < (i->second).size(); indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = (i->second)[indx2];

            Pair->SetMagField2(col2->magField());
            fillResMatrix(selectedtracks_1[col1->index()], selectedtracks_1[col2->index()]); // filling res. matrix -- ME identical
          }
        }
      }

    } else { //====================================== non-identical ======================================

      for (auto i = mixbins.begin(); i != mixbins.end(); i++) { // iterating over all vertex&mult bins

        for (int indx1 = 0; indx1 < (i->second).size(); indx1++) { // iterating over all selected collisions with selected tracks1

          auto col1 = (i->second)[indx1];

          Pair->SetMagField1(col1->magField());
          Pair->SetMagField2(col1->magField());

          fillEtaPhi(selectedtracks_1[col1->index()], selectedtracks_2[col1->index()]); // filling deta(dphi*) -- SE non-identical

          for (int indx2 = indx1 + 1; indx2 < (i->second).size(); indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = (i->second)[indx2];

            Pair->SetMagField2(col2->magField());
            fillResMatrix(selectedtracks_1[col1->index()], selectedtracks_2[col2->index()]); // filling res. matrix -- ME non-identical
          }
        }
      }

    } //====================================== end of mixing non-identical ======================================

    // clearing up
    for (auto i = selectedtracks_1.begin(); i != selectedtracks_1.end(); i++)
      (i->second).clear();
    selectedtracks_1.clear();

    if (!IsIdentical) {
      for (auto i = selectedtracks_2.begin(); i != selectedtracks_2.end(); i++)
        (i->second).clear();
      selectedtracks_2.clear();
    }

    for (auto i = mixbins.begin(); i != mixbins.end(); i++)
      (i->second).clear();
    mixbins.clear();
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FemtoCorrelationsMC>(cfgc)};
}
