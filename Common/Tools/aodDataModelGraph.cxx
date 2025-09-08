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

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>

#include <fmt/printf.h>

#include <cstdio>
#include <string>
#include <utility>

using namespace o2::framework;
using namespace o2::aod;
using namespace o2::soa;

static int count = 0;
static int width = 10;
static int height = 10;

void inline graphSize()
{
  fmt::printf(
    R"(
size="%d,%d";
)",
    width, height);
}

struct Style {
  const char* color;
  const char* background;
  const char* fontcolor;
  const char* headerfontcolor;
  const char* headerbgcolor;
  const char* methodcolor;
  const char* methodbgcolor;
  const char* indexcolor;
  const char* indexbgcolor;
};

static Style styles[] = {
  {"black", "gray80", "black", "black", "gray70", "black", "gray60", "black", "gray50"},
  {"/reds9/2", "/reds9/4", "white", "white", "/reds9/7", "black", "/reds9/6", "/reds9/1", "/reds9/5"},
  {"/greens9/2", "/greens9/4", "white", "white", "/greens9/7", "black", "/greens9/6", "/greens9/1", "/greens9/5"},
  {"/blues9/2", "/blues9/4", "white", "white", "/blues9/7", "black", "/blues9/6", "/blues9/1", "/blues9/5"},
};

Style const& getDefaultStyle()
{
  return styles[0];
}

enum StyleType : int {
  DEFAULT = 0,
  RED = 1,
  GREEN = 2,
  BLUE = 3,
};

static std::array<std::pair<std::string, StyleType>, 14> tableStyles = {
  std::make_pair("HfTrackIndexProng", StyleType::BLUE),
  std::make_pair("HfCandProng", StyleType::BLUE),
  std::make_pair("pidResp", StyleType::GREEN),
  std::make_pair("Mults", StyleType::GREEN),
  std::make_pair("CentRun2V0Ms", StyleType::GREEN),
  std::make_pair("Timestamps", StyleType::GREEN),
  std::make_pair("Jet", StyleType::BLUE),
  std::make_pair("Mc", StyleType::RED),
  std::make_pair("V0Datas", StyleType::GREEN),
  std::make_pair("CascData", StyleType::GREEN),
  std::make_pair("TrackSelection", StyleType::GREEN),
  std::make_pair("TracksDCA", StyleType::GREEN),
  std::make_pair("Transient", StyleType::GREEN),
  std::make_pair("Extension", StyleType::GREEN),
};

template <typename T>
Style getStyleFor()
{
  auto label = MetadataTrait<T>::metadata::tableLabel();
  auto entry = std::find_if(tableStyles.begin(), tableStyles.end(),
                            [&](auto&& x) {
                              if (std::string(label).find(x.first) != std::string::npos) {
                                return true;
                              }
                              return false;
                            });
  if (entry != tableStyles.end()) {
    auto value = *entry;
    return styles[value.second];
  }
  return styles[StyleType::DEFAULT];
}

void inline nodeEmpty()
{
  fmt::printf(
    R"(node[shape=none,height=0,width=0,label=""])");
}

void inline nodeNormal()
{
  fmt::printf(
    R"(node[shape=plain,style=filled,fillcolor=gray95])");
}

void inline graphHeader(char const* type, char const* name)
{
  fmt::printf(R"(%s %s {
edge[dir=back, arrowtail=empty]
)",
              type, name);
  nodeNormal();
}

void inline graphFooter()
{
  fmt::printf("}\n");
}

template <typename T>
void displayEntity();

template <typename... Ts>
void displayOriginals(pack<Ts...>)
{
  graphHeader("subgraph", fmt::format("cluster_{}", count++).c_str());
  fmt::printf("label = %s;\n", MetadataTrait<pack_element_t<1, pack<Ts...>>>::metadata::tableLabel());
  (..., displayEntity<Ts>());
  graphFooter();
}

template <typename C>
void printColumn(char const* fg, char const* bg)
{
  fmt::printf("<TR><TD color='%s' bgcolor='%s'>%s</TD></TR>", fg, bg, C::columnLabel());
}

template <typename... C>
int displayPersistentColumns(pack<C...>, const char* fg, const char* bg)
{
  int n = 0;
  ([&]() {
    if constexpr (o2::soa::is_persistent_v<C> && !o2::soa::is_index_column_v<C>) {
      printColumn<C>(fg, bg);
      ++n;
    }
  }(),
   ...);
  if (n > 0) {
    fmt::printf("%s", "\n");
  }
  return n;
}

template <typename... C>
int displayDynamicColumns(pack<C...>, const char* fg, const char* bg)
{
  int n = 0;
  ([&]() {
    if constexpr (!o2::soa::is_persistent_v<C>) {
      printColumn<C>(fg, bg);
      ++n;
    }
  }(),
   ...);
  if (n > 0) {
    fmt::printf("%s", "\n");
  }
  return n;
}

template <typename... C>
void displayIndexColumns(pack<C...>, char const* fg, char const* bg)
{
  ([&]() {
    if constexpr (o2::soa::is_index_column_v<C>) {
      printColumn<C>(fg, bg);
    }
  }(),
   ...);
  fmt::printf("%s", "\n");
}

template <typename C, typename T>
void printIndex()
{
  if constexpr (!is_type_with_originals_v<typename C::binding_t>) {
    auto a = MetadataTrait<typename C::binding_t>::metadata::tableLabel();
    auto b = MetadataTrait<T>::metadata::tableLabel();
    fmt::printf("%s -> %s []\n", a, b);
  } else {
    using main_original = pack_element_t<1, typename C::binding_t::originals>;
    auto a = MetadataTrait<main_original>::metadata::tableLabel();
    auto b = MetadataTrait<T>::metadata::tableLabel();
    fmt::printf("%s -> %s []\n", a, b);
  }
}

template <typename T, typename... C>
void dumpIndex(pack<C...>)
{
  ([&]() {
    if constexpr (o2::soa::is_index_column_v<C>) {
      printIndex<C, T>();
    }
  }(),
   ...);
  fmt::printf("%s", "\n");
}

template <typename T>
void displayTable()
{
  auto style = getStyleFor<T>();
  auto label = MetadataTrait<T>::metadata::tableLabel();
  fmt::printf(R"(%s[color="%s" cellpadding="0" fillcolor="%s" fontcolor="%s" label = <
<TABLE cellpadding='2' cellspacing='0' cellborder='0' ><TH cellpadding='0' bgcolor="black"><TD bgcolor="%s"><font color="%s">%s</font></TD></TH>)",
              label, style.color, style.background, style.fontcolor, style.headerbgcolor, style.headerfontcolor, label);
  if (displayPersistentColumns(typename T::table_t::columns{}, style.color, style.background) > 0) {
    fmt::printf("%s", "HR");
  }
  if (displayDynamicColumns(typename T::table_t::columns{}, style.methodcolor, style.methodbgcolor) > 0) {
    fmt::printf("%s", "HR");
  }
  displayIndexColumns(typename T::table_t::columns{}, style.indexcolor, style.indexbgcolor);
  fmt::printf("%s", "</TABLE>\n>]\n");
  dumpIndex<T>(typename T::table_t::columns{});
}

template <typename T>
void displayEntity()
{
  if constexpr (is_soa_join_t<T>::value) {
    displayOriginals(typename T::originals{});
  } else {
    displayTable<T>();
  }
}

template <typename... T>
void displayEntities()
{
  graphHeader("subgraph", fmt::format("cluster_{}", count++).c_str());
  (..., displayEntity<T>());
  graphFooter();
}

int main(int, char**)
{
  graphHeader("digraph", "hierarchy");
  graphSize();
  fmt::printf(R"(compound = true;
)");

  displayEntity<BCs>();
  /// rank trick to avoid BCs moving
  nodeEmpty();
  fmt::printf(R"({rank = same; BCs -> root[style=invis];};)");
  nodeNormal();

  displayEntity<Zdcs>();
  displayEntity<FT0s>();
  displayEntity<FV0As>();
  displayEntity<FDDs>();
  displayEntity<HMPIDs>();

  displayEntities<Collisions, CentRun2V0Ms, Mults, Timestamps>();
  displayEntity<McCollisions>();
  displayEntity<McCollisionLabels>();

  displayEntity<Calos>();
  displayEntity<CaloTriggers>();
  displayEntity<McCaloLabels>();

  displayEntity<FV0Cs>();
  displayEntity<Run2BCInfos>();

  displayEntities<Tracks, TracksCov, TracksExtra, TracksDCA, TrackSelection,
                  pidTOFFullEl, pidTOFFullMu, pidTOFFullPi,
                  pidTOFFullKa, pidTOFFullPr, pidTOFFullDe,
                  pidTOFFullTr, pidTOFFullHe, pidTOFFullAl,
                  pidTPCFullEl, pidTPCFullMu, pidTPCFullPi,
                  pidTPCFullKa, pidTPCFullPr, pidTPCFullDe,
                  pidTPCFullTr, pidTPCFullHe, pidTPCFullAl>();
  displayEntity<AmbiguousTracks>();
  displayEntity<AmbiguousMFTTracks>();

  displayEntity<McParticles>();
  displayEntity<McTrackLabels>();

  displayEntity<ChargedJets>();
  displayEntity<ChargedJetConstituents>();

  displayEntities<V0s, V0Datas>();

  displayEntities<Cascades, CascDataFull>();

  displayEntities<MFTTracks, FwdTracks, FwdTracksCov>();

  displayEntities<Hf2Prongs, HfCand2Prong>();
  displayEntities<Hf3Prongs, HfCand3Prong>();

  graphFooter();
  return 0;
}
