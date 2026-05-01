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
/// \file onTheFlyDecayer.cxx
/// \brief pre-processing for on-the-fly analysis
/// \author Jesper Karlsson Gumprecht <jesper.gumprecht@cern.ch>
///

#include "ALICE3/Core/Decayer.h"
#include "ALICE3/Core/TrackUtilities.h"
#include "ALICE3/DataModel/OTFMCParticle.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <sys/types.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <gsl/span>
#include <map>
#include <ostream>
#include <span>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

static constexpr int NumDecays = 7;
static constexpr int NumParameters = 1;
static constexpr int DefaultParameters[NumDecays][NumParameters]{{1}, {1}, {1}, {1}, {1}, {1}, {1}};
static constexpr float Tolerance = 1e-7f;
static const std::vector<std::string> ParameterNames{"enable"};
static const std::vector<std::string> ParticleNames{"K0s",
                                                    "Lambda",
                                                    "Anti-Lambda",
                                                    "Xi",
                                                    "Anti-Xi",
                                                    "Omega",
                                                    "Anti-Omega"};

static const std::vector<int> pdgCodes{PDG_t::kK0Short,
                                       PDG_t::kLambda0,
                                       PDG_t::kLambda0Bar,
                                       PDG_t::kXiMinus,
                                       PDG_t::kXiPlusBar,
                                       PDG_t::kOmegaMinus,
                                       PDG_t::kOmegaPlusBar};

enum class V0 { kNeg,
                kPos,
                kSize };
enum class AntiV0 { kPos,
                    kNeg,
                    kSize };
enum class Cascade { kBach,
                     kV0,
                     kNeg,
                     kPos,
                     kSize };
enum class AntiCascade { kBach,
                         kV0,
                         kPos,
                         kNeg,
                         kSize };

struct OnTheFlyDecayer {
  Produces<aod::McPartWithDaus> tableMcParticlesWithDau;

  o2::upgrade::Decayer decayer;
  Service<o2::framework::O2DatabasePDG> pdgDB;
  std::map<int, std::vector<o2::upgrade::OTFParticle>> mDecayDaughters;

  Configurable<int> seed{"seed", 0, "Set seed for particle decayer"};
  Configurable<float> magneticField{"magneticField", 20., "Magnetic field (kG)"};
  Configurable<float> maxEta{"maxEta", 2.5, "Only decay particles that appear within selected eta range"};
  Configurable<LabeledArray<int>> enabledDecays{"enabledDecays",
                                                {DefaultParameters[0], NumDecays, NumParameters, ParticleNames, ParameterNames},
                                                "Enable option for particle to be decayed: 0 - no, 1 - yes"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis axisRadius{"axisRadius", {200, 0, 100}, "Radial axis"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};

  ConfigurableAxis axisLogRadius{"axisLogRadius", {VARIABLE_WIDTH, 0.0f, 0.01f, 0.0104713f, 0.0109648f, 0.0114815f, 0.0120226f, 0.0125893f, 0.0131826f, 0.0138038f, 0.0144544f, 0.0151356f, 0.0158489f, 0.0165959f, 0.017378f, 0.018197f, 0.0190546f, 0.0199526f, 0.020893f, 0.0218776f, 0.0229087f, 0.0239883f, 0.0251189f, 0.0263027f, 0.0275423f, 0.0288403f, 0.0301995f, 0.0316228f, 0.0331131f, 0.0346737f, 0.0363078f, 0.0380189f, 0.0398107f, 0.0416869f, 0.0436516f, 0.0457088f, 0.047863f, 0.0501187f, 0.0524807f, 0.0549541f, 0.057544f, 0.060256f, 0.0630957f, 0.0660693f, 0.0691831f, 0.0724436f, 0.0758578f, 0.0794328f, 0.0831764f, 0.0870964f, 0.0912011f, 0.0954993f, 0.1f, 0.104713f, 0.109648f, 0.114815f, 0.120226f, 0.125893f, 0.131826f, 0.138038f, 0.144544f, 0.151356f, 0.158489f, 0.165959f, 0.17378f, 0.18197f, 0.190546f, 0.199526f, 0.20893f, 0.218776f, 0.229087f, 0.239883f, 0.251189f, 0.263027f, 0.275423f, 0.288403f, 0.301995f, 0.316228f, 0.331131f, 0.346737f, 0.363078f, 0.380189f, 0.398107f, 0.416869f, 0.436516f, 0.457088f, 0.47863f, 0.501187f, 0.524807f, 0.549541f, 0.57544f, 0.60256f, 0.630957f, 0.660693f, 0.691831f, 0.724436f, 0.758578f, 0.794328f, 0.831764f, 0.870964f, 0.912011f, 0.954993f, 1.0f, 1.04713f, 1.09648f, 1.14815f, 1.20226f, 1.25893f, 1.31826f, 1.38038f, 1.44544f, 1.51356f, 1.58489f, 1.65959f, 1.7378f, 1.8197f, 1.90546f, 1.99526f, 2.0893f, 2.18776f, 2.29087f, 2.39883f, 2.51189f, 2.63027f, 2.75423f, 2.88403f, 3.01995f, 3.16228f, 3.31131f, 3.46737f, 3.63078f, 3.80189f, 3.98107f, 4.16869f, 4.36516f, 4.57088f, 4.7863f, 5.01187f, 5.24807f, 5.49541f, 5.7544f, 6.0256f, 6.30957f, 6.60693f, 6.91831f, 7.24436f, 7.58578f, 7.94328f, 8.31764f, 8.70964f, 9.12011f, 9.54993f, 10.0f, 10.4713f, 10.9648f, 11.4815f, 12.0226f, 12.5893f, 13.1826f, 13.8038f, 14.4544f, 15.1356f, 15.8489f, 16.5959f, 17.378f, 18.197f, 19.0546f, 19.9526f, 20.893f, 21.8776f, 22.9087f, 23.9883f, 25.1189f, 26.3027f, 27.5423f, 28.8403f, 30.1995f, 31.6228f, 33.1131f, 34.6737f, 36.3078f, 38.0189f, 39.8107f, 41.6869f, 43.6516f, 45.7088f, 47.863f, 50.1187f, 52.4807f, 54.9541f, 57.544f, 60.256f, 63.0957f, 66.0693f, 69.1831f, 72.4436f, 75.8578f, 79.4328f, 83.1764f, 87.0964f, 91.2011f, 95.4993f, 100.0f}, "Radial axis"};

  ConfigurableAxis axisLogPt{"axisLogPt", {VARIABLE_WIDTH, 0.0f, 0.001f, 0.00104713f, 0.00109648f, 0.00114815f, 0.00120226f, 0.00125893f, 0.00131826f, 0.00138038f, 0.00144544f, 0.00151356f, 0.00158489f, 0.00165959f, 0.0017378f, 0.0018197f, 0.00190546f, 0.00199526f, 0.0020893f, 0.00218776f, 0.00229087f, 0.00239883f, 0.00251189f, 0.00263027f, 0.00275423f, 0.00288403f, 0.00301995f, 0.00316228f, 0.00331131f, 0.00346737f, 0.00363078f, 0.00380189f, 0.00398107f, 0.00416869f, 0.00436516f, 0.00457088f, 0.0047863f, 0.00501187f, 0.00524807f, 0.00549541f, 0.0057544f, 0.0060256f, 0.00630957f, 0.00660693f, 0.00691831f, 0.00724436f, 0.00758578f, 0.00794328f, 0.00831764f, 0.00870964f, 0.00912011f, 0.00954993f, 0.01f, 0.0104713f, 0.0109648f, 0.0114815f, 0.0120226f, 0.0125893f, 0.0131826f, 0.0138038f, 0.0144544f, 0.0151356f, 0.0158489f, 0.0165959f, 0.017378f, 0.018197f, 0.0190546f, 0.0199526f, 0.020893f, 0.0218776f, 0.0229087f, 0.0239883f, 0.0251189f, 0.0263027f, 0.0275423f, 0.0288403f, 0.0301995f, 0.0316228f, 0.0331131f, 0.0346737f, 0.0363078f, 0.0380189f, 0.0398107f, 0.0416869f, 0.0436516f, 0.0457088f, 0.047863f, 0.0501187f, 0.0524807f, 0.0549541f, 0.057544f, 0.060256f, 0.0630957f, 0.0660693f, 0.0691831f, 0.0724436f, 0.0758578f, 0.0794328f, 0.0831764f, 0.0870964f, 0.0912011f, 0.0954993f, 0.1f, 0.104713f, 0.109648f, 0.114815f, 0.120226f, 0.125893f, 0.131826f, 0.138038f, 0.144544f, 0.151356f, 0.158489f, 0.165959f, 0.17378f, 0.18197f, 0.190546f, 0.199526f, 0.20893f, 0.218776f, 0.229087f, 0.239883f, 0.251189f, 0.263027f, 0.275423f, 0.288403f, 0.301995f, 0.316228f, 0.331131f, 0.346737f, 0.363078f, 0.380189f, 0.398107f, 0.416869f, 0.436516f, 0.457088f, 0.47863f, 0.501187f, 0.524807f, 0.549541f, 0.57544f, 0.60256f, 0.630957f, 0.660693f, 0.691831f, 0.724436f, 0.758578f, 0.794328f, 0.831764f, 0.870964f, 0.912011f, 0.954993f, 1.0f, 1.04713f, 1.09648f, 1.14815f, 1.20226f, 1.25893f, 1.31826f, 1.38038f, 1.44544f, 1.51356f, 1.58489f, 1.65959f, 1.7378f, 1.8197f, 1.90546f, 1.99526f, 2.0893f, 2.18776f, 2.29087f, 2.39883f, 2.51189f, 2.63027f, 2.75423f, 2.88403f, 3.01995f, 3.16228f, 3.31131f, 3.46737f, 3.63078f, 3.80189f, 3.98107f, 4.16869f, 4.36516f, 4.57088f, 4.7863f, 5.01187f, 5.24807f, 5.49541f, 5.7544f, 6.0256f, 6.30957f, 6.60693f, 6.91831f, 7.24436f, 7.58578f, 7.94328f, 8.31764f, 8.70964f, 9.12011f, 9.54993f, 10.0f}, "pt axis for QA histograms"};

  struct McParticleAlice3 {
    McParticleAlice3() = default;
    ~McParticleAlice3() = default;
    McParticleAlice3(const McParticleAlice3& src) = default;
    McParticleAlice3(int collisionId,
                     int pdgCode,
                     int statusCode,
                     int flags,
                     int mother0,
                     int mother1,
                     int daughter0,
                     int daughter1,
                     float weight,
                     float px, float py, float pz, float e,
                     float vx, float vy, float vz, float vt,
                     float phi, float eta, float pt, float p, float y,
                     bool isAlive, bool isPrimary) : collisionId(collisionId),
                                                     pdgCode(pdgCode),
                                                     statusCode(statusCode),
                                                     flags(flags),
                                                     mothersIds{mother0, mother1},
                                                     daughtersIdSlice{daughter0, daughter1},
                                                     weight(weight),
                                                     px(px),
                                                     py(py),
                                                     pz(pz),
                                                     e(e),
                                                     vx(vx),
                                                     vy(vy),
                                                     vz(vz),
                                                     vt(vt),
                                                     phi(phi),
                                                     eta(eta),
                                                     pt(pt),
                                                     p(p),
                                                     y(y),
                                                     isAlive(isAlive),
                                                     isPrimary(isPrimary) {}

    bool hasNaN() const
    {
      return std::isnan(px) || std::isnan(py) || std::isnan(pz) || std::isnan(e) ||
             std::isnan(vx) || std::isnan(vy) || std::isnan(vz) || std::isnan(vt) ||
             std::isnan(phi) || std::isnan(eta) || std::isnan(pt) || std::isnan(p) ||
             std::isnan(y) || std::isnan(weight);
    }

    int collisionId;
    int pdgCode;
    int statusCode;
    int flags;
    int mothersIds[2];
    int daughtersIdSlice[2];
    float weight;
    float px, py, pz, e;
    float vx, vy, vz, vt;
    float phi, eta, pt, p, y;
    bool isAlive;
    bool isPrimary;
  };

  template <typename TEnumerate>
  int idx(TEnumerate enumerate)
  {
    return static_cast<int>(enumerate);
  }

  bool checkDecayChannel(const int pdgMother, const std::vector<o2::upgrade::OTFParticle>& daus)
  {
    switch (pdgMother) {
      case PDG_t::kK0Short:
        if (daus.size() != static_cast<size_t>(V0::kSize)) {
          return false;
        }
        return (daus[idx(V0::kNeg)].pdgCode() == PDG_t::kPiMinus &&
                daus[idx(V0::kPos)].pdgCode() == PDG_t::kPiPlus);
      case PDG_t::kLambda0:
        if (daus.size() != static_cast<size_t>(V0::kSize)) {
          return false;
        }
        return (daus[idx(V0::kNeg)].pdgCode() == PDG_t::kPiMinus &&
                daus[idx(V0::kPos)].pdgCode() == PDG_t::kProton);
      case PDG_t::kLambda0Bar:
        if (daus.size() != static_cast<size_t>(AntiV0::kSize)) {
          return false;
        }
        return (daus[idx(AntiV0::kNeg)].pdgCode() == PDG_t::kProtonBar &&
                daus[idx(AntiV0::kPos)].pdgCode() == PDG_t::kPiPlus);
      case PDG_t::kXiMinus:
        if (daus.size() != static_cast<size_t>(Cascade::kSize)) {
          return false;
        }
        return (daus[idx(Cascade::kBach)].pdgCode() == PDG_t::kPiMinus &&
                daus[idx(Cascade::kV0)].pdgCode() == PDG_t::kLambda0 &&
                daus[idx(Cascade::kNeg)].pdgCode() == PDG_t::kPiMinus &&
                daus[idx(Cascade::kPos)].pdgCode() == PDG_t::kProton);
      case PDG_t::kXiPlusBar:
        if (daus.size() != static_cast<size_t>(AntiCascade::kSize)) {
          return false;
        }
        return (daus[idx(AntiCascade::kBach)].pdgCode() == PDG_t::kPiPlus &&
                daus[idx(AntiCascade::kV0)].pdgCode() == PDG_t::kLambda0Bar &&
                daus[idx(AntiCascade::kPos)].pdgCode() == PDG_t::kPiPlus &&
                daus[idx(AntiCascade::kNeg)].pdgCode() == PDG_t::kProtonBar);
      case PDG_t::kOmegaMinus:
        if (daus.size() != static_cast<size_t>(Cascade::kSize)) {
          return false;
        }
        return (daus[idx(Cascade::kBach)].pdgCode() == PDG_t::kKMinus &&
                daus[idx(Cascade::kV0)].pdgCode() == PDG_t::kLambda0 &&
                daus[idx(Cascade::kNeg)].pdgCode() == PDG_t::kPiMinus &&
                daus[idx(Cascade::kPos)].pdgCode() == PDG_t::kProton);
      case PDG_t::kOmegaPlusBar:
        if (daus.size() != static_cast<size_t>(AntiCascade::kSize)) {
          return false;
        }
        return (daus[idx(AntiCascade::kBach)].pdgCode() == PDG_t::kKPlus &&
                daus[idx(AntiCascade::kV0)].pdgCode() == PDG_t::kLambda0Bar &&
                daus[idx(AntiCascade::kPos)].pdgCode() == PDG_t::kPiPlus &&
                daus[idx(AntiCascade::kNeg)].pdgCode() == PDG_t::kProtonBar);
      default:
        return false;
    }
  }

  std::vector<int> mEnabledDecays;
  void init(o2::framework::InitContext&)
  {
    LOG(info) << "Initializing on-the-fly-decayer.";
    LOG(info) << "Using seed: " << seed;
    LOG(info) << "Using magnetic field: " << magneticField;
    decayer.setSeed(seed);
    decayer.setBField(magneticField);
    for (int i = 0; i < NumDecays; ++i) {
      if (enabledDecays->get(ParticleNames[i].c_str(), "enable")) {
        LOG(info) << "Decay enabled: " << ParticleNames[i].c_str();
        mEnabledDecays.push_back(pdgCodes[i]);
      }
    }

    auto hNaNBookkeeping = histos.add<TH1>("hNaNBookkeeping", "hNaNBookkeeping", kTH1D, {{2, -0.5, 1.5}});
    hNaNBookkeeping->GetXaxis()->SetBinLabel(1, "OK");
    hNaNBookkeeping->GetXaxis()->SetBinLabel(2, "NaN");

    histos.add("K0S/hGenK0S", "hGenK0S;pT (GeV)", kTH1D, {axisPt});
    histos.add("K0S/hPosDauDecayRadius", "hPosDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("K0S/hNegDauDecayRadius", "hNegDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});

    histos.add("Lambda/hGenLambda", "hGenLambda;pT (GeV)", kTH1D, {axisPt});
    histos.add("Lambda/hPosDauDecayRadius", "hPosDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("Lambda/hNegDauDecayRadius", "hNegDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});

    histos.add("AntiLambda/hGenAntiLambda", "hGenAntiLambda;pT (GeV)", kTH1D, {axisPt});
    histos.add("AntiLambda/hPosDauDecayRadius", "hPosDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("AntiLambda/hNegDauDecayRadius", "hNegDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});

    histos.add("Xi/hGenXi", "hGenXi;pT (GeV)", kTH1D, {axisPt});
    histos.add("Xi/hV0DauDecayRadius", "hV0DauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("Xi/hBachDauDecayRadius", "hBachDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("Xi/hPosDauDecayRadius", "hPosDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("Xi/hNegDauDecayRadius", "hNegDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});

    histos.add("AntiXi/hGenAntiXi", "hGenAntiXi;pT (GeV)", kTH1D, {axisPt});
    histos.add("AntiXi/hV0DauDecayRadius", "hV0DauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("AntiXi/hBachDauDecayRadius", "hBachDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("AntiXi/hPosDauDecayRadius", "hPosDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("AntiXi/hNegDauDecayRadius", "hNegDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});

    histos.add("Omega/hGenOmega", "hGenOmega;pT (GeV)", kTH1D, {axisPt});
    histos.add("Omega/hV0DauDecayRadius", "hV0DauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("Omega/hBachDauDecayRadius", "hBachDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("Omega/hPosDauDecayRadius", "hPosDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("Omega/hNegDauDecayRadius", "hNegDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});

    histos.add("AntiOmega/hGenAntiOmega", "hGenAntiOmega;pT (GeV)", kTH1D, {axisPt});
    histos.add("AntiOmega/hV0DauDecayRadius", "hV0DauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("AntiOmega/hBachDauDecayRadius", "hBachDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("AntiOmega/hPosDauDecayRadius", "hPosDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});
    histos.add("AntiOmega/hNegDauDecayRadius", "hNegDauDecayRadius;Radius decay vtx 2D (cm); Daughter pT (GeV)", kTH2D, {axisLogRadius, axisLogPt});

    histos.add("Secondaries/hGenEl", "hGenEl;pT (GeV)", kTH1D, {axisPt});
    histos.add("Secondaries/hGenMu", "hGenMu;pT (GeV)", kTH1D, {axisPt});
    histos.add("Secondaries/hGenPi", "hGenPi;pT (GeV)", kTH1D, {axisPt});
    histos.add("Secondaries/hGenKa", "hGenKa;pT (GeV)", kTH1D, {axisPt});
    histos.add("Secondaries/hGenPr", "hGenPr;pT (GeV)", kTH1D, {axisPt});
  }

  bool canDecay(const int pdgCode)
  {
    return std::find(mEnabledDecays.begin(), mEnabledDecays.end(), pdgCode) != mEnabledDecays.end();
  }

  std::vector<McParticleAlice3> mcParticlesAlice3;
  void process(aod::McCollision const&, aod::McParticles const& mcParticles)
  {
    mDecayDaughters.clear();
    mcParticlesAlice3.clear();
    u_int64_t nStoredDaughters = 0;
    for (int index{0}; index < static_cast<int>(mcParticles.size()); ++index) {
      const auto& particle = mcParticles.rawIteratorAt(index);
      std::vector<o2::upgrade::OTFParticle> decayDaughters, decayStack;
      if (canDecay(particle.pdgCode()) && std::abs(particle.eta()) < maxEta) {
        o2::upgrade::OTFParticle mother(particle);
        decayStack = decayer.decayParticle(pdgDB, mother);
        while (!decayStack.empty()) {
          o2::upgrade::OTFParticle otfParticle = decayStack.back();
          decayStack.pop_back();

          const bool stable = !canDecay(otfParticle.pdgCode());
          otfParticle.setIsAlive(stable);
          decayDaughters.push_back(otfParticle);

          if (stable) {
            continue;
          }

          std::vector<o2::upgrade::OTFParticle> daughters = decayer.decayParticle(pdgDB, otfParticle);
          for (const o2::upgrade::OTFParticle& dau : daughters) {
            decayStack.push_back(dau);
          }
        }

        if (decayDaughters.empty()) {
          LOG(error) << "Attempted to decay " << particle.pdgCode() << " but resulting vector of daugthers were empty";
          continue;
        }

        switch (particle.pdgCode()) { // Do QA on specific channels
          case PDG_t::kK0Short:
            if (!checkDecayChannel(PDG_t::kK0Short, decayDaughters)) {
              break;
            }
            histos.fill(HIST("K0S/hGenK0S"), particle.pt());
            histos.fill(HIST("K0S/hPosDauDecayRadius"), decayDaughters[idx(V0::kPos)].radius(), decayDaughters[idx(V0::kPos)].pt());
            histos.fill(HIST("K0S/hNegDauDecayRadius"), decayDaughters[idx(V0::kNeg)].radius(), decayDaughters[idx(V0::kNeg)].pt());
            break;

          case PDG_t::kLambda0:
            if (!checkDecayChannel(PDG_t::kLambda0, decayDaughters)) {
              break;
            }
            histos.fill(HIST("Lambda/hGenLambda"), particle.pt());
            histos.fill(HIST("Lambda/hPosDauDecayRadius"), decayDaughters[idx(V0::kPos)].radius(), decayDaughters[idx(V0::kPos)].pt());
            histos.fill(HIST("Lambda/hNegDauDecayRadius"), decayDaughters[idx(V0::kNeg)].radius(), decayDaughters[idx(V0::kNeg)].pt());
            break;

          case PDG_t::kLambda0Bar:
            if (!checkDecayChannel(PDG_t::kLambda0Bar, decayDaughters)) {
              break;
            }
            histos.fill(HIST("AntiLambda/hGenAntiLambda"), particle.pt());
            histos.fill(HIST("AntiLambda/hPosDauDecayRadius"), decayDaughters[idx(AntiV0::kPos)].radius(), decayDaughters[idx(AntiV0::kPos)].pt());
            histos.fill(HIST("AntiLambda/hNegDauDecayRadius"), decayDaughters[idx(AntiV0::kNeg)].radius(), decayDaughters[idx(AntiV0::kNeg)].pt());
            break;

          case PDG_t::kXiMinus:
            if (!checkDecayChannel(PDG_t::kXiMinus, decayDaughters)) {
              break;
            }
            histos.fill(HIST("Xi/hGenXi"), particle.pt());
            histos.fill(HIST("Xi/hBachDauDecayRadius"), decayDaughters[idx(Cascade::kBach)].radius(), decayDaughters[idx(Cascade::kBach)].pt());
            histos.fill(HIST("Xi/hV0DauDecayRadius"), decayDaughters[idx(Cascade::kV0)].radius(), decayDaughters[idx(Cascade::kV0)].pt());
            histos.fill(HIST("Xi/hPosDauDecayRadius"), decayDaughters[idx(Cascade::kPos)].radius(), decayDaughters[idx(Cascade::kPos)].pt());
            histos.fill(HIST("Xi/hNegDauDecayRadius"), decayDaughters[idx(Cascade::kNeg)].radius(), decayDaughters[idx(Cascade::kNeg)].pt());
            break;

          case PDG_t::kXiPlusBar:
            if (!checkDecayChannel(PDG_t::kXiPlusBar, decayDaughters)) {
              break;
            }
            histos.fill(HIST("AntiXi/hGenAntiXi"), particle.pt());
            histos.fill(HIST("AntiXi/hBachDauDecayRadius"), decayDaughters[idx(AntiCascade::kBach)].radius(), decayDaughters[idx(AntiCascade::kBach)].pt());
            histos.fill(HIST("AntiXi/hV0DauDecayRadius"), decayDaughters[idx(AntiCascade::kV0)].radius(), decayDaughters[idx(AntiCascade::kV0)].pt());
            histos.fill(HIST("AntiXi/hPosDauDecayRadius"), decayDaughters[idx(AntiCascade::kPos)].radius(), decayDaughters[idx(AntiCascade::kPos)].pt());
            histos.fill(HIST("AntiXi/hNegDauDecayRadius"), decayDaughters[idx(AntiCascade::kNeg)].radius(), decayDaughters[idx(AntiCascade::kNeg)].pt());
            break;

          case PDG_t::kOmegaMinus:
            if (!checkDecayChannel(PDG_t::kOmegaMinus, decayDaughters)) {
              break;
            }
            histos.fill(HIST("Omega/hGenOmega"), particle.pt());
            histos.fill(HIST("Omega/hBachDauDecayRadius"), decayDaughters[idx(Cascade::kBach)].radius(), decayDaughters[idx(Cascade::kBach)].pt());
            histos.fill(HIST("Omega/hV0DauDecayRadius"), decayDaughters[idx(Cascade::kV0)].radius(), decayDaughters[idx(Cascade::kV0)].pt());
            histos.fill(HIST("Omega/hPosDauDecayRadius"), decayDaughters[idx(Cascade::kPos)].radius(), decayDaughters[idx(Cascade::kPos)].pt());
            histos.fill(HIST("Omega/hNegDauDecayRadius"), decayDaughters[idx(Cascade::kNeg)].radius(), decayDaughters[idx(Cascade::kNeg)].pt());
            break;

          case PDG_t::kOmegaPlusBar:
            if (!checkDecayChannel(PDG_t::kOmegaPlusBar, decayDaughters)) {
              break;
            }
            histos.fill(HIST("AntiOmega/hGenAntiOmega"), particle.pt());
            histos.fill(HIST("AntiOmega/hBachDauDecayRadius"), decayDaughters[idx(AntiCascade::kBach)].radius(), decayDaughters[idx(AntiCascade::kBach)].pt());
            histos.fill(HIST("AntiOmega/hV0DauDecayRadius"), decayDaughters[idx(AntiCascade::kV0)].radius(), decayDaughters[idx(AntiCascade::kV0)].pt());
            histos.fill(HIST("AntiOmega/hPosDauDecayRadius"), decayDaughters[idx(AntiCascade::kPos)].radius(), decayDaughters[idx(AntiCascade::kPos)].pt());
            histos.fill(HIST("AntiOmega/hNegDauDecayRadius"), decayDaughters[idx(AntiCascade::kNeg)].radius(), decayDaughters[idx(AntiCascade::kNeg)].pt());
            break;

          default:
            break;
        }
      }

      int daughtersIdSlice[2];
      if (canDecay(particle.pdgCode())) {
        daughtersIdSlice[0] = idx(mcParticles.size() + nStoredDaughters);
        daughtersIdSlice[1] = idx(mcParticles.size() + nStoredDaughters + decayDaughters.size());
      } else {
        daughtersIdSlice[0] = idx(particle.daughtersIds()[0]);
        daughtersIdSlice[1] = idx(particle.daughtersIds()[1]);
      }

      mDecayDaughters.emplace(index, decayDaughters);
      nStoredDaughters += decayDaughters.size();

      const float phi = o2::constants::math::PI + std::atan2(-1.0f * particle.py(), -1.0f * particle.px());
      float eta; // As https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1922
      const float pt = std::sqrt(particle.px() * particle.px() + particle.py() * particle.py());
      const float p = std::sqrt(particle.px() * particle.px() + particle.py() * particle.py() + particle.pz() * particle.pz());
      float y; // As https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1943

      if ((p - particle.pz()) < Tolerance) {
        eta = (particle.pz() < 0.0f) ? -100.0f : 100.0f;
      } else {
        eta = 0.5f * std::log((p + particle.pz()) / (p - particle.pz()));
      }

      if ((particle.e() - particle.pz()) < Tolerance) {
        y = (particle.pz() < 0.0f) ? -100.0f : 100.0f;
      } else {
        y = 0.5f * std::log((particle.e() + particle.pz()) / (particle.e() - particle.pz()));
      }

      // TODO: Particle status code
      // TODO: Expression columns
      auto mothers = particle.mothersIds();
      int mother0 = mothers.size() > 0 ? mothers[0] : -1;
      int mother1 = mothers.size() > 1 ? mothers[1] : mother0;
      mcParticlesAlice3.push_back(McParticleAlice3{particle.mcCollisionId(), particle.pdgCode(), particle.statusCode(),
                                                   particle.flags(), mother0, mother1,
                                                   daughtersIdSlice[0], daughtersIdSlice[1], particle.weight(),
                                                   particle.px(), particle.py(), particle.pz(), particle.e(),
                                                   particle.vx(), particle.vy(), particle.vz(), particle.vt(),
                                                   phi, eta, pt, p, y, !canDecay(particle.pdgCode()), true});
    }

    int daughtersIdSlice[2] = {-1, -1};
    for (const auto& [index, decayDaughters] : mDecayDaughters) {
      for (const auto& dau : decayDaughters) {
        if (index >= mcParticles.size()) {
          LOG(error) << "--- Index " << index << " out of bounds for mcParticles table of size " << mcParticles.size() << std::endl;
          continue;
        }

        const float phi = o2::constants::math::PI + std::atan2(-1.0f * dau.py(), -1.0f * dau.px());
        float eta; // Conditional as https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1922
        const float pt = std::sqrt(dau.px() * dau.px() + dau.py() * dau.py());
        const float p = std::sqrt(dau.px() * dau.px() + dau.py() * dau.py() + dau.pz() * dau.pz());
        float y; // Conditional as https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1943

        if ((p - dau.pz()) < Tolerance) {
          eta = (dau.pz() < 0.0f) ? -100.0f : 100.0f;
        } else {
          eta = 0.5f * std::log((p + dau.pz()) / (p - dau.pz()));
        }

        if ((dau.e() - dau.pz()) < Tolerance) {
          y = (dau.pz() < 0.0f) ? -100.0f : 100.0f;
        } else {
          y = 0.5f * std::log((dau.e() + dau.pz()) / (dau.e() - dau.pz()));
        }

        switch (dau.pdgCode()) {
          case PDG_t::kElectron:
            histos.fill(HIST("Secondaries/hGenEl"), pt);
            break;

          case PDG_t::kMuonMinus:
            histos.fill(HIST("Secondaries/hGenMu"), pt);
            break;

          case PDG_t::kPiPlus:
            histos.fill(HIST("Secondaries/hGenPi"), pt);
            break;

          case PDG_t::kKPlus:
            histos.fill(HIST("Secondaries/hGenKa"), pt);
            break;

          case PDG_t::kProton:
            histos.fill(HIST("Secondaries/hGenPr"), pt);
            break;

          default:
            break;
        }

        // TODO: Particle status code
        // TODO: Expression columns
        // TODO: vt
        auto mother = mcParticles.rawIteratorAt(index);
        mcParticlesAlice3.push_back(McParticleAlice3{mother.mcCollisionId(), dau.pdgCode(), 1,
                                                     -1, index, index, daughtersIdSlice[0], daughtersIdSlice[1], mother.weight(),
                                                     dau.px(), dau.py(), dau.pz(), dau.e(),
                                                     dau.vx(), dau.vy(), dau.vz(), mother.vt(),
                                                     phi, eta, pt, p, y, dau.isAlive(), false});
      }
    }

    for (const auto& particle : mcParticlesAlice3) {
      if (particle.hasNaN()) {
        histos.fill(HIST("hNaNBookkeeping"), 1);
        continue;
      }

      histos.fill(HIST("hNaNBookkeeping"), 0);
      std::span<const int> motherSpan(particle.mothersIds, 2);
      tableMcParticlesWithDau(particle.collisionId, particle.pdgCode, particle.statusCode,
                              particle.flags, motherSpan, particle.daughtersIdSlice, particle.weight,
                              particle.px, particle.py, particle.pz, particle.e,
                              particle.vx, particle.vy, particle.vz, particle.vt,
                              particle.phi, particle.eta, particle.pt, particle.p, particle.y,
                              particle.isAlive, particle.isPrimary);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyDecayer>(cfgc)};
}
