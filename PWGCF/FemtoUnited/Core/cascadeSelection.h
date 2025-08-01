// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file cascadeSelection.h
/// \brief cascade selection
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_CASCADESELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_CASCADESELECTION_H_

#include "PWGCF/FemtoUnited/Core/baseSelection.h"
#include "PWGCF/FemtoUnited/Core/dataTypes.h"
#include "PWGCF/FemtoUnited/Core/modes.h"

#include "Framework/Configurable.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace o2::analysis::femtounited
{
namespace cascadeselection
{

struct ConfCascadeFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CascadeFilters");
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 99.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massXiMin{"massXiMin", 1.2f, "Minimum Xi mass"};
  o2::framework::Configurable<float> massXiMax{"massXiMax", 1.4f, "Maximum Xi mass"};
  o2::framework::Configurable<float> rejectMassXiMin{"rejectMassXiMin", 1.317f, "Reject Minimum Xi mass for Omega hypothesis"};
  o2::framework::Configurable<float> rejectMassXiMax{"rejectMassXiMax", 1.325f, "Rejection Maximum Xi mass for Omega hypothesis"};
  o2::framework::Configurable<float> massOmegaMin{"massOmegaMin", 1.5f, "Minimum Omega mass"};
  o2::framework::Configurable<float> massOmegaMax{"massOmegaMax", 1.9f, "Maximum Omega mass"};
  o2::framework::Configurable<float> rejectMassOmegaMin{"rejectMassOmegaMin", 1.668f, "Reject minimum Omega mass for Xi hypothesis"};
  o2::framework::Configurable<float> rejectMassOmegaMax{"rejectMassOmegaMax", 1.676f, "Reject maximum Omega mass for Xi hypothesis"};
  o2::framework::Configurable<float> massLambdaMin{"massLambdaMin", 1.0f, "Minimum Lambda mass"};
  o2::framework::Configurable<float> massLambdaMax{"massLambdaMax", 1.2f, "Maximum Lambda mass"};
};

#define CASCADE_DEFAULT_BITS                                                                                                                               \
  o2::framework::Configurable<std::vector<float>> cascadeCpaMin{"cascadeCpaMin", {0.95f}, "Minimum cosine of pointing angle"};                             \
  o2::framework::Configurable<std::vector<float>> cascadeTransRadMin{"cascadeTransRadMin", {0.9f}, "Minimum transverse radius (cm)"};                      \
  o2::framework::Configurable<std::vector<float>> cascadeDcaDauMax{"cascadeDcaDauMax", {0.25f}, "Maximum DCA between the daughters at decay vertex (cm)"}; \
  o2::framework::Configurable<std::vector<float>> lambdaCpaMin{"lambdaCpaMin", {0.78f}, "Minimum cosine of pointing angle"};                               \
  o2::framework::Configurable<std::vector<float>> lambdaTransRadMin{"lambdaTransRadMin", {0.9f}, "Minimum transverse radius (cm)"};                        \
  o2::framework::Configurable<std::vector<float>> lambdaDcaDauMax{"lambdaDcaDauMax", {0.5f}, "Maximum DCA between the daughters at decay vertex (cm)"};    \
  o2::framework::Configurable<std::vector<float>> lambdaDcaToPvMin{"lambdaDcaToPvMin", {0.3f}, "Minimum DCA between the lambda and primary vertex"};       \
  o2::framework::Configurable<std::vector<float>> dauAbsEtaMax{"dauAbsEtaMax", {0.8f}, "Minimum DCA of the daughters from primary vertex (cm)"};           \
  o2::framework::Configurable<std::vector<float>> dauDcaMin{"dauDcaMin", {0.05f}, "Minimum DCA of the daughters from primary vertex (cm)"};                \
  o2::framework::Configurable<std::vector<float>> dauTpcClustersMin{"dauTpcClustersMin", {80.f}, "Minimum number of TPC clusters for daughter tracks"};    \
  o2::framework::Configurable<std::vector<float>> posDauTpc{"posDauTpc", {5.f}, "Maximum |nsimga_Pion/Proton| TPC for positive daughter tracks"};          \
  o2::framework::Configurable<std::vector<float>> negDauTpc{"negDauTpc", {5.f}, "Maximum |nsimga_Pion/Proton| TPC for negative daughter tracks"};

struct ConfXiBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("XiBits");
  CASCADE_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> bachelorTpcPion{"bachelorTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for bachelor tracks"};
};

struct ConfOmegaBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("OmegaBits");
  CASCADE_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> bachelorTpcKaon{"bachelorTpcKaon", {5.f}, "Maximum |nsimga_Kaon| TPC for bachelor tracks"};
};

#undef CASCADE_DEFAULT_BITS

#define CASCADE_DEFAULT_SELECTION(defaultMassMin, defaultMassMax, defaultPdgCode)                              \
  o2::framework::Configurable<int> pdgCode{"pdgCode", defaultPdgCode, "Track PDG code"};                       \
  o2::framework::Configurable<int> sign{"sign", 1, "Sign of the Lambda (+1 for Lambda and -1 for Antilambda"}; \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                        \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                                      \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                                   \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                                    \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};                                     \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};        \
  o2::framework::Configurable<float> massMin{"massMin", defaultMassMin, "Minimum invariant mass for Cascade"}; \
  o2::framework::Configurable<float> massMax{"massMax", defaultMassMax, "Maximum invariant mass for Cascade"}; \
  o2::framework::Configurable<o2::aod::femtodatatypes::CascadeMaskType> mask{"mask", 0, "Bitmask for cascade selection"};

struct ConfXiSelection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("XiSelection");
  CASCADE_DEFAULT_SELECTION(1.22, 1.42, 3212)
};

struct ConfOmegaSelection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("OmegaSelection");
  CASCADE_DEFAULT_SELECTION(1.57, 1.77, 3334)
};

/// The different selections this task is capable of doing
enum CascadeSels {
  // selections for cascades
  kCascadeCpaMin,      ///< Min. CPA (cosine pointing angle)
  kCascadeDcaDaughMax, ///< Max. DCA of the daughers at decay vertex
  kCascadeTransRadMin, ///< max. transverse radius

  // selection for lambda daughter
  kLambdaCpaMin,      ///< Min. DCA of the lambda daughers at primary vertex
  kLambdaDcaDauMax,   ///< TPC PID for daughters (Pion/Proton)
  kLambdaTransRadMin, ///< Min. number of TPC clusters of daughter
  kLambdaDcaToPvMin,  ///< Min. DCA to primary vertex of daughter lambda

  // selection for bachelor/daugthers
  kDauAbsEtaMax, ///< Min. DCA of the daughers/bachelor at primary vertex
  kDauTpcClsMin, ///< Min. number of TPC clusters of daughters/bachelor
  kDauDcaMin,    ///< TPC Pion PID for negative daughter

  // PID selection for cascade bachelor
  kBachelorTpcPion, ///< TPC Pion PID for bachelor
  kBachelorTpcKaon, ///< TPC Kaon PID for bachelor
                    ///
  // PID selection for lambda daughers
  kPosDauTpc, ///< TPC PID for positive daughter
  kNegDauTpc, ///< TPC PID for negative daughter

  kCascadeSelsMax
};

const std::string omegaSelsName = std::string("Omega Selection Object");
const std::string xiSelsName = std::string("Xi Selection Object");
const std::unordered_map<CascadeSels, std::string> CascadeSelsNames = {
  {kCascadeCpaMin, "Cascade CPA Min"},
  {kCascadeDcaDaughMax, "Cascade DCA Daughters Max"},
  {kCascadeTransRadMin, "Cascade Transverse Radius Min"},

  {kLambdaCpaMin, "Lambda CPA Min"},
  {kLambdaDcaDauMax, "Lambda DCA Daughter Max"},
  {kLambdaTransRadMin, "Lambda Transverse Radius Min"},
  {kLambdaDcaToPvMin, "Lambda DCA to PV Min"},

  {kDauAbsEtaMax, "Daughter Abs Eta Max"},
  {kDauTpcClsMin, "Daughter TPC Clusters Min"},
  {kDauDcaMin, "Daughter DCA Min"},

  {kBachelorTpcPion, "Bachelor TPC Pion PID"},
  {kBachelorTpcKaon, "Bachelor TPC Kaon PID"},

  {kPosDauTpc, "Positive Daughter TPC PID"},
  {kNegDauTpc, "Negative Daughter TPC PID"},

  {kCascadeSelsMax, "Cascade Selections Max"}};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
template <modes::Cascade cascadeType>
class CascadeSelection : public BaseSelection<float, o2::aod::femtodatatypes::CascadeMaskType, kCascadeSelsMax>
{
 public:
  CascadeSelection() {}
  virtual ~CascadeSelection() = default;

  template <typename T1, typename T2>
  void configure(T1 const& config, T2 const& filter)
  {
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      mXiMassLowerLimit = filter.massXiMin.value;
      mXiMassUpperLimit = filter.massXiMax.value;
      mOmegaMassLowerLimit = filter.rejectMassOmegaMin.value;
      mOmegaMassUpperLimit = filter.rejectMassOmegaMax.value;
      this->addSelection(config.bachelorTpcPion.value, kBachelorTpcPion, limits::kAbsUpperLimit, true, true);
    }
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kOmega)) {
      mOmegaMassLowerLimit = filter.massOmegaMin.value;
      mOmegaMassUpperLimit = filter.massOmegaMax.value;
      mXiMassLowerLimit = filter.rejectMassXiMin.value;
      mXiMassUpperLimit = filter.rejectMassXiMax.value;
      this->addSelection(config.bachelorTpcKaon.value, kBachelorTpcKaon, limits::kAbsUpperLimit, true, true);
    }

    this->addSelection(config.posDauTpc.value, kPosDauTpc, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.negDauTpc.value, kNegDauTpc, limits::kAbsUpperLimit, true, true);

    this->addSelection(config.cascadeCpaMin.value, kCascadeCpaMin, limits::kLowerLimit, true, true);
    this->addSelection(config.cascadeTransRadMin.value, kCascadeTransRadMin, limits::kLowerLimit, true, true);
    this->addSelection(config.cascadeDcaDauMax.value, kCascadeDcaDaughMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.lambdaCpaMin.value, kLambdaCpaMin, limits::kLowerLimit, true, true);
    this->addSelection(config.lambdaTransRadMin.value, kLambdaTransRadMin, limits::kLowerLimit, true, true);
    this->addSelection(config.lambdaDcaDauMax.value, kLambdaDcaDauMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.lambdaDcaToPvMin.value, kLambdaDcaToPvMin, limits::kLowerLimit, true, true);
    this->addSelection(config.dauAbsEtaMax.value, kDauAbsEtaMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.dauDcaMin.value, kDauDcaMin, limits::kAbsLowerLimit, true, true);
    this->addSelection(config.dauTpcClustersMin.value, kDauTpcClsMin, limits::kLowerLimit, true, true);
  };

  template <typename T1, typename T2, typename T3>
  void applySelections(T1 const& cascade, T2 const& /*tracks*/, T3 const& col)
  {
    this->reset();
    // cascade selections
    this->evaluateObservable(kCascadeCpaMin, cascade.casccosPA(col.posX(), col.posY(), col.posZ()));
    this->evaluateObservable(kCascadeDcaDaughMax, cascade.dcacascdaughters());
    this->evaluateObservable(kCascadeTransRadMin, cascade.cascradius());

    // lambda selection
    this->evaluateObservable(kLambdaCpaMin, cascade.v0cosPA(col.posX(), col.posY(), col.posZ()));
    this->evaluateObservable(kLambdaDcaDauMax, cascade.dcaV0daughters());
    this->evaluateObservable(kLambdaTransRadMin, cascade.v0radius());
    this->evaluateObservable(kLambdaDcaToPvMin, cascade.dcav0topv(col.posX(), col.posY(), col.posZ()));

    auto bachelor = cascade.template bachelor_as<T2>();
    auto posDaughter = cascade.template posTrack_as<T2>();
    auto negDaughter = cascade.template negTrack_as<T2>();

    // daughter selections
    std::array<float, 3> etaDaughters = {std::fabs(bachelor.eta()), std::fabs(posDaughter.eta()), std::fabs(negDaughter.eta())};
    this->evaluateObservable(kDauAbsEtaMax, *std::max_element(etaDaughters.begin(), etaDaughters.end()));

    std::array<float, 3> dcaDaughters = {std::hypot(bachelor.dcaXY(), bachelor.dcaZ()),
                                         std::hypot(posDaughter.dcaXY(), posDaughter.dcaZ()),
                                         std::hypot(negDaughter.dcaXY(), negDaughter.dcaZ())};
    this->evaluateObservable(kDauDcaMin, *std::min_element(dcaDaughters.begin(), dcaDaughters.end()));

    std::array<float, 3> clustersDaughters = {1.f * bachelor.tpcNClsFound(), 1.f * posDaughter.tpcNClsFound(), 1.f * negDaughter.tpcNClsFound()};
    this->evaluateObservable(kDauTpcClsMin, *std::min_element(clustersDaughters.begin(), clustersDaughters.end()));

    // bachelor pid selection
    // check both pion and kaon PID for xi and omega
    this->evaluateObservable(kBachelorTpcPion, bachelor.tpcNSigmaPi());
    this->evaluateObservable(kBachelorTpcKaon, bachelor.tpcNSigmaKa());

    // depending on the charge, we check lambda and antilambda hypothesis
    if (cascade.sign() < 0) {
      this->evaluateObservable(kPosDauTpc, posDaughter.tpcNSigmaPr());
      this->evaluateObservable(kNegDauTpc, negDaughter.tpcNSigmaPi());
    } else if (cascade.sign() > 0) {
      this->evaluateObservable(kPosDauTpc, posDaughter.tpcNSigmaPi());
      this->evaluateObservable(kNegDauTpc, negDaughter.tpcNSigmaPr());
    } else {
      LOG(warn) << "Encountered Cascade candidate with 0 charge";
    }

    this->assembleBitmask();
  };

  template <typename T>
  bool checkHypothesis(T const& cascadeCandidate)
  {
    // no need to check PID of the bachelor/daughters here, they are set as minimal cuts
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      return (mXiMassLowerLimit < cascadeCandidate.mXi() && mXiMassUpperLimit > cascadeCandidate.mXi()) &&           // inside xi mass window
             (cascadeCandidate.mOmega() < mOmegaMassLowerLimit || cascadeCandidate.mOmega() > mOmegaMassUpperLimit); // outside omega mass window
    }
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      return (mOmegaMassLowerLimit < cascadeCandidate.mOmega() && mOmegaMassUpperLimit > cascadeCandidate.mOmega()) && // inside omega mass window
             (cascadeCandidate.mXi() < mXiMassLowerLimit || cascadeCandidate.mXi() > mXiMassUpperLimit);               // outside xi mass window
    }
    return false;
  }

 protected:
  float mXiMassLowerLimit = 0.f;
  float mXiMassUpperLimit = 999.f;

  float mOmegaMassLowerLimit = 0.f;
  float mOmegaMassUpperLimit = 999.f;
};
} // namespace cascadeselection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_CASCADESELECTION_H_
