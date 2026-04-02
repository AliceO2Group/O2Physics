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

/// \file particleOriginAnalysis.cxx
/// \brief Classifies filtered physical primary identified particles as prompt
///        (from string fragmentation) or from resonance decay by inspecting
///        the MC mother chain. Produces per-species prompt/decay fractions
///        and mother identity distributions vs pT and centrality.
///        Designed for generator-level MC studies of balance function
///        origin decomposition across collision systems.
///        Consumes the filtered tables produced by dptDptFilter.
/// \author victor.gonzalez.sebastian@gmail.com

#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include "PWGCF/DataModel/DptDptFiltered.h"
#include "PWGCF/TableProducer/dptDptFilter.h"

#include "Common/Core/TableHelper.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include <TH2F.h>
#include <TH3F.h>
#include <TPDGCode.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;

#define FORMATSTRING(theformat, theparams...) TString::Format(theformat, theparams).Data()

namespace particleorigintask
{
using namespace o2::analysis::dptdptfilter;

// ============================================================================
// Compact encoding of the mother resonance identity
// ============================================================================
enum MotherSpecies {
  kMothPrompt = 0,  ///< not from decay (prompt from hadronization)
  kMothRho0,        ///< rho0(770)        - PDG 113
  kMothRhoCharged,  ///< rho+/-(770)      - PDG 213
  kMothOmega782,    ///< omega(782)       - PDG 223
  kMothEta,         ///< eta              - PDG 221
  kMothEtaPrime,    ///< eta'(958)        - PDG 331
  kMothKStar892Ch,  ///< K*(892)+/-       - PDG 323
  kMothKStar892_0,  ///< K*(892)0         - PDG 313
  kMothPhi1020,     ///< phi(1020)        - PDG 333
  kMothDelta,       ///< Delta(1232) all  - PDG 1114,2114,2214,2224
  kMothSigmaStar,   ///< Sigma*(1385)     - PDG 3114,3214,3224
  kMothLambda1520,  ///< Lambda(1520)     - PDG 3124
  kMothF0_980,      ///< f0(980)          - PDG 9010221
  kMothOtherMeson,  ///< other meson mothers
  kMothOtherBaryon, ///< other baryon mothers
  kNMotherSpecies
};

static const char* motherLabel[kNMotherSpecies] = {
  "prompt", "#rho^{0}", "#rho^{#pm}", "#omega", "#eta", "#eta'",
  "K*^{#pm}", "K*^{0}", "#phi",
  "#Delta", "#Sigma*", "#Lambda(1520)",
  "f_{0}(980)", "other meson", "other baryon"};

// PDG codes used from TPDGCode.h (ROOT):
//   kRho770_0 (113), kRho770Plus (213), kLambda1520 (3124),
//   kPiPlus (211), kKPlus (321), kProton (2212)
// PDG codes used from CommonConstants/PhysicsConstants.h (O2):
//   o2::constants::physics::Pdg::kEta (221), kOmega (223), kEtaPrime (331),
//   kK0Star892 (313), kKPlusStar892 (323), kPhi (333)
// PDG codes NOT in either header - defined here:
static constexpr int KPdgDeltaMinusMinus = 1114; // o2-linter: disable=pdg/explicit-code(not existing)
static constexpr int KPdgDelta0 = 2114;          // o2-linter: disable=pdg/explicit-code(not existing)
static constexpr int KPdgDeltaPlus = 2214;       // o2-linter: disable=pdg/explicit-code(not existing)
static constexpr int KPdgDeltaPlusPlus = 2224;   // o2-linter: disable=pdg/explicit-code(not existing)
static constexpr int KPdgSigmaStarMinus = 3114;  // o2-linter: disable=pdg/explicit-code(not existing)
static constexpr int KPdgSigmaStar0 = 3214;      // o2-linter: disable=pdg/explicit-code(not existing)
static constexpr int KPdgSigmaStarPlus = 3224;   // o2-linter: disable=pdg/explicit-code(not existing)
static constexpr int KPdgF0_980 = 9010221;       // o2-linter: disable=pdg/explicit-code(not existing),name/function-variable(clashes with f0),name/constexpr-constant(clashes with f0)

/// PDG codes below this threshold are mesons; at or above are baryons
static constexpr int KPdgBaryonThreshold = 1000; // o2-linter: disable=pdg/explicit-code(not a PDG code)

/// \brief Encode absolute PDG code of a mother into compact bin index
inline int encodeMotherPDG(int absPdg)
{
  using namespace o2::constants::physics;
  switch (absPdg) {
    case kRho770_0:
      return kMothRho0;
    case kRho770Plus: /* kRho770Minus is just -213, we use abs */
      return kMothRhoCharged;
    case Pdg::kOmega:
      return kMothOmega782;
    case Pdg::kEta:
      return kMothEta;
    case Pdg::kEtaPrime:
      return kMothEtaPrime;
    case Pdg::kKPlusStar892:
      return kMothKStar892Ch;
    case Pdg::kK0Star892:
      return kMothKStar892_0;
    case Pdg::kPhi:
      return kMothPhi1020;
    case KPdgDeltaMinusMinus:
    case KPdgDelta0:
    case KPdgDeltaPlus:
    case KPdgDeltaPlusPlus:
      return kMothDelta;
    case KPdgSigmaStarMinus:
    case KPdgSigmaStar0:
    case KPdgSigmaStarPlus:
      return kMothSigmaStar;
    case kLambda1520:
      return kMothLambda1520;
    case KPdgF0_980:
      return kMothF0_980;
    default:
      return (absPdg < KPdgBaryonThreshold) ? kMothOtherMeson : kMothOtherBaryon;
  }
}

/// PDG codes above this threshold correspond to hadrons (mesons and baryons).
/// Below are quarks (1-6), leptons (11-16), gauge bosons (21-25), and
/// special/internal generator codes.
static constexpr int KPdgHadronThreshold = 100; // o2-linter: disable=pdg/explicit-code(not a PDG code)

// ============================================================================
// Classification utilities
// ============================================================================

/// \brief Check if a mother is a beam or initial-state particle
///        using its status code rather than just its PDG code
inline bool isBeamParticle(auto const& mother)
{
  int status = mother.getGenStatusCode();
  // Pythia status codes 11-19: beam particles and event-setup entries
  // HepMC status 4: beam particles
  static constexpr int BEAMSTATUSCODEMIN = 11;
  static constexpr int BEAMSTATUSCODEMAX = 19;
  static constexpr int HEPMCBEAMSTATUSCODE = 4;

  return (status >= BEAMSTATUSCODEMIN && status <= BEAMSTATUSCODEMAX) || mother.getHepMCStatusCode() == HEPMCBEAMSTATUSCODE;
}

/// \brief Classify a physical primary by its immediate mother
/// \return {isFromDecay, motherPdgCode}
template <typename McParticleObject, typename McParticlesTable>
inline std::pair<bool, int> classifyImmediate(McParticleObject const& particle,
                                              McParticlesTable const& /*mcParticles*/)
{
  if (!particle.has_mothers()) {
    return {false, 0};
  }
  auto mother = particle.template mothers_first_as<McParticlesTable>();
  int absMomPdg = std::abs(mother.pdgCode());
  if (absMomPdg > KPdgHadronThreshold && !isBeamParticle(mother)) {
    return {true, mother.pdgCode()};
  }
  return {false, 0};
}

/// \brief Trace back to the earliest hadronic ancestor
/// \return {isFromDecay, originalAncestorPdgCode}
template <typename McParticleObject, typename McParticlesTable>
inline std::pair<bool, int> classifyAncestor(McParticleObject const& particle,
                                             McParticlesTable const& mcParticles)
{
  /* reuse classifyImmediate for the first step */
  auto [isFromDecay, immediatePdg] = classifyImmediate(particle, mcParticles);
  if (!isFromDecay) {
    return {false, 0};
  }

  /* the immediate mother is hadronic; now walk up to find the original resonance */
  auto current = particle.template mothers_first_as<McParticlesTable>();
  int originalPdg = current.pdgCode();
  constexpr int KMaxDepth = 20;
  for (int depth = 0; current.has_mothers() && depth < KMaxDepth; ++depth) {
    auto grandparent = current.template mothers_first_as<McParticlesTable>();
    if (std::abs(grandparent.pdgCode()) <= KPdgHadronThreshold || isBeamParticle(grandparent)) {
      break;
    }
    originalPdg = grandparent.pdgCode();
    current = grandparent;
  }
  return {true, originalPdg};
}

} // namespace particleorigintask

// ============================================================================
// The analysis task
// ============================================================================
struct ParticleOriginAnalysis {

  /* the histogram registry */
  HistogramRegistry registry{"ParticleOriginData", {}, OutputObjHandlingPolicy::AnalysisObject};

  /* the species and track names from configuration */
  std::vector<std::string> poinames; ///< the species of interest names
  std::vector<std::string> tnames;   ///< the track names

  /* the number of track ids from configuration */
  size_t nch = 0;

  /* histogram pointers for direct access */
  /* per track id histograms: indexed by trackacceptedid */
  static constexpr int KNMo = particleorigintask::kNMotherSpecies;

  std::vector<std::shared_ptr<TH1>> fhPromptVsPt;         ///< prompt counts vs pT, per track id
  std::vector<std::shared_ptr<TH1>> fhDecayVsPt;          ///< from-decay counts vs pT, per track id
  std::vector<std::shared_ptr<TH2>> fhPromptVsCentVsPt;   ///< prompt counts vs (cent, pT), per track id
  std::vector<std::shared_ptr<TH2>> fhDecayVsCentVsPt;    ///< from-decay counts vs (cent, pT), per track id
  std::vector<std::shared_ptr<TH3>> fhMotherVsPtVsCent;   ///< immediate mother (encoded) vs pT vs cent, per track id
  std::vector<std::shared_ptr<TH3>> fhAncestorVsPtVsCent; ///< earliest ancestor (encoded) vs pT vs cent, per track id
  std::vector<std::shared_ptr<TH1>> fhMotherPDG;          ///< immediate mother |PDG| (fine, 1D safety net), per track id
  std::vector<std::shared_ptr<TH1>> fhAncestorPDG;        ///< earliest ancestor |PDG| (fine, 1D safety net), per track id

  void init(InitContext& initContext)
  {
    using namespace particleorigintask;
    using namespace o2::analysis::dptdptfilter;

    /* self configure the binning from the filter task */
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mPTbins", ptbins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mPTmin", ptlow, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mPTmax", ptup, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mEtabins", etabins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mEtamin", etalow, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mEtamax", etaup, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mZVtxbins", zvtxbins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mZVtxmin", zvtxlow, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgBinning.mZVtxmax", zvtxup, false);

    /* self configure the desired species */
    o2::analysis::dptdptfilter::PIDSpeciesSelection pidselector;
    std::vector<std::string> cfgnames = {"cfgElectronPIDSelection", "cfgMuonPIDSelection", "cfgPionPIDSelection", "cfgKaonPIDSelection", "cfgProtonPIDSelection"};
    std::vector<uint8_t> spids = {0, 1, 2, 3, 4};
    for (uint i = 0; i < cfgnames.size(); ++i) {
      auto includeIt = [&pidselector, &initContext](int spid, auto name) {
        bool mUseIt = false;
        bool mExcludeIt = false;
        if (getTaskOptionValue(initContext, "dpt-dpt-filter-tracks", TString::Format("%s.mUseIt", name.c_str()).Data(), mUseIt, false) &&
            getTaskOptionValue(initContext, "dpt-dpt-filter-tracks", TString::Format("%s.mExclude", name.c_str()).Data(), mExcludeIt, false)) {
          if (mUseIt && !mExcludeIt) {
            auto cfg = new o2::analysis::TrackSelectionPIDCfg();
            cfg->mUseIt = true;
            cfg->mExclude = false;
            pidselector.addSpecies(spid, cfg);
          }
        }
      };
      includeIt(spids[i], cfgnames[i]);
    }
    uint nspecies = pidselector.getNSpecies();
    if (nspecies == 0) {
      /* unidentified analysis */
      poinames.push_back(pidselector.getHadFName());
      tnames.push_back(std::string(TString::Format("%sP", pidselector.getHadFName()).Data()));
      tnames.push_back(std::string(TString::Format("%sM", pidselector.getHadFName()).Data()));
      LOGF(info, "Incorporated species name %s to the analysis", poinames[0].c_str());
    } else {
      for (uint8_t ix = 0; ix < nspecies; ++ix) {
        poinames.push_back(std::string(pidselector.getSpeciesFName(ix)));
        tnames.push_back(std::string(TString::Format("%sP", pidselector.getSpeciesFName(ix)).Data()));
        tnames.push_back(std::string(TString::Format("%sM", pidselector.getSpeciesFName(ix)).Data()));
        LOGF(info, "Incorporated species name %s to the analysis", poinames[ix].c_str());
      }
    }
    nch = tnames.size();

    /* self configure centrality/multiplicity ranges from the filter task */
    int nCentBins = 1;
    std::vector<double> centBinEdges = {0.0f, 100.0f};
    std::string centspec;
    if (getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgCentSpec", centspec, false)) {
      LOGF(info, "Got the centralities specification: %s", centspec.c_str());
      auto tokens = TString(centspec.c_str()).Tokenize(",");
      nCentBins = tokens->GetEntries();
      centBinEdges.clear();
      for (int i = 0; i < nCentBins; ++i) {
        double cmmin = 0.0;
        double cmmax = 0.0;
        sscanf(tokens->At(i)->GetName(), "%lf-%lf", &cmmin, &cmmax);
        if (i == 0) {
          centBinEdges.push_back(cmmin);
        }
        centBinEdges.push_back(cmmax);
      }
      delete tokens;
    } else {
      LOGF(info, "No centralities specification. Setting it to: 0-100");
    }

    /* build the centrality axis with variable bin edges */
    const AxisSpec centAxis{centBinEdges, "centrality (%)"};
    const AxisSpec ptAxis{ptbins, ptlow, ptup, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec motherAxis{KNMo, -0.5f, KNMo - 0.5f, "mother species"};
    const AxisSpec pdgAxis{100, 0.5f, 100.5f, "species"};
    const AxisSpec zvtxAxis{zvtxbins, zvtxlow, zvtxup, "#it{z}_{vtx}"};

    /* resize the histogram vectors */
    fhPromptVsPt.resize(nch, nullptr);
    fhDecayVsPt.resize(nch, nullptr);
    fhPromptVsCentVsPt.resize(nch, nullptr);
    fhDecayVsCentVsPt.resize(nch, nullptr);
    fhMotherVsPtVsCent.resize(nch, nullptr);
    fhAncestorVsPtVsCent.resize(nch, nullptr);
    fhMotherPDG.resize(nch, nullptr);
    fhAncestorPDG.resize(nch, nullptr);

    /* per track id histograms */
    for (uint i = 0; i < nch; ++i) {
      const char* tname = tnames[i].c_str();

      fhPromptVsPt[i] = registry.add<TH1>(
        FORMATSTRING("PromptVsPt_%s", tname),
        FORMATSTRING("Prompt %s;#it{p}_{T} (GeV/#it{c});counts", tname),
        kTH1F, {ptAxis});

      fhDecayVsPt[i] = registry.add<TH1>(
        FORMATSTRING("DecayVsPt_%s", tname),
        FORMATSTRING("From decay %s;#it{p}_{T} (GeV/#it{c});counts", tname),
        kTH1F, {ptAxis});

      fhPromptVsCentVsPt[i] = registry.add<TH2>(
        FORMATSTRING("PromptVsCentVsPt_%s", tname),
        FORMATSTRING("Prompt %s;centrality (%%);#it{p}_{T} (GeV/#it{c})", tname),
        kTH2F, {centAxis, ptAxis});

      fhDecayVsCentVsPt[i] = registry.add<TH2>(
        FORMATSTRING("DecayVsCentVsPt_%s", tname),
        FORMATSTRING("From decay %s;centrality (%%);#it{p}_{T} (GeV/#it{c})", tname),
        kTH2F, {centAxis, ptAxis});

      fhMotherVsPtVsCent[i] = registry.add<TH3>(
        FORMATSTRING("MotherVsPtVsCent_%s", tname),
        FORMATSTRING("Immediate mother of %s;mother;#it{p}_{T} (GeV/#it{c});centrality (%%)", tname),
        kTH3F, {motherAxis, ptAxis, centAxis});

      fhAncestorVsPtVsCent[i] = registry.add<TH3>(
        FORMATSTRING("AncestorVsPtVsCent_%s", tname),
        FORMATSTRING("Earliest ancestor of %s;ancestor;#it{p}_{T} (GeV/#it{c});centrality (%%)", tname),
        kTH3F, {motherAxis, ptAxis, centAxis});

      /* label the encoded mother/ancestor axis */
      for (int im = 0; im < KNMo; ++im) {
        fhMotherVsPtVsCent[i]->GetXaxis()->SetBinLabel(im + 1, motherLabel[im]);
        fhAncestorVsPtVsCent[i]->GetXaxis()->SetBinLabel(im + 1, motherLabel[im]);
      }

      fhMotherPDG[i] = registry.add<TH1>(
        FORMATSTRING("MotherPDG_%s", tname),
        FORMATSTRING("Immediate mother PDG of %s from decay;PDG code;counts", tname),
        kTH1F, {pdgAxis});

      fhAncestorPDG[i] = registry.add<TH1>(
        FORMATSTRING("AncestorPDG_%s", tname),
        FORMATSTRING("Earliest ancestor PDG of %s from decay;PDG code;counts", tname),
        kTH1F, {pdgAxis});
    }
  }

  /// \brief Fill the origin classification histograms for one accepted particle
  /// \param tid the trackacceptedid from the filter
  /// \param centmult the centrality/multiplicity percentile
  /// \param pt the particle transverse momentum
  /// \param isFromDecay true if immediate mother is a hadron
  /// \param immediatePdg the PDG code of the immediate mother
  /// \param isFromDecayFull true if earliest ancestor is a hadron
  /// \param ancestorPdg the PDG code of the earliest hadronic ancestor
  void fillOrigin(int tid, float centmult, float pt,
                  bool isFromDecay, int immediatePdg,
                  bool isFromDecayFull, int ancestorPdg)
  {
    using namespace particleorigintask;

    if (isFromDecay) {
      fhDecayVsPt[tid]->Fill(pt);
      fhDecayVsCentVsPt[tid]->Fill(centmult, pt);
      int encodedMother = encodeMotherPDG(std::abs(immediatePdg));
      fhMotherVsPtVsCent[tid]->Fill(static_cast<float>(encodedMother), pt, centmult);
      fhMotherPDG[tid]->Fill(TString::Format("%d", immediatePdg).Data(), 1.0);
    } else {
      fhPromptVsPt[tid]->Fill(pt);
      fhPromptVsCentVsPt[tid]->Fill(centmult, pt);
      fhMotherVsPtVsCent[tid]->Fill(static_cast<float>(kMothPrompt), pt, centmult);
    }
    if (isFromDecayFull) {
      int encodedAncestor = encodeMotherPDG(std::abs(ancestorPdg));
      fhAncestorVsPtVsCent[tid]->Fill(static_cast<float>(encodedAncestor), pt, centmult);
      fhAncestorPDG[tid]->Fill(TString::Format("%d", ancestorPdg).Data(), 1.0);
    } else {
      fhAncestorVsPtVsCent[tid]->Fill(static_cast<float>(kMothPrompt), pt, centmult);
    }
  }

  Filter onlyacceptedcollisions = (aod::dptdptfilter::collisionaccepted == uint8_t(true));

  /// \brief Process: one accepted generator-level collision at a time
  /// with its associated MC particles joined with filter information.
  /// Follows the processGenLevelNotStored / processSame pattern.
  void process(soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>::iterator const& collision,
               soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo> const& particles)
  {
    using namespace particleorigintask;
    using namespace o2::analysis::dptdptfilter;

    float centmult = collision.centmult();

    for (auto const& particle : particles) {
      int8_t tid = particle.trackacceptedid();
      if (tid < 0) {
        /* particle not accepted */
        continue;
      }
      float pt = particle.pt();

      /* classify origin using the McParticles mother chain */
      auto [isFromDecay, immediatePdg] = classifyImmediate(particle, particles);
      auto [isFromDecayFull, ancestorPdg] = classifyAncestor(particle, particles);

      /* fill histograms indexed by trackacceptedid */
      fillOrigin(tid, centmult, pt, isFromDecay, immediatePdg,
                 isFromDecayFull, ancestorPdg);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ParticleOriginAnalysis>(cfgc)};
}
