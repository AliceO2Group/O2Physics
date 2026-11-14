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

/// \file checkTables.C
/// \brief Checks which strangeness tables (and which versions) are present in a derived AO2D
///        and reports which steps of the V0 tutorial (o2-analysistutorial-lf-strangeness-derived-v0-*) can run.
///
/// Usage:  root -l -b -q 'checkTables.C("data/AO2D.root")' or something like root -l -b -q "checkTables.C(\"$(cat input_data.txt)\")"

#include <TCollection.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TString.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

struct RequiredTable {
  std::string base;  // tree name without version, e.g. "o2straevsels"
  int version{0};    // version required by the current O2Physics data model
  int firstStep{-1}; // first step that uses it (-1 = skeleton)
  std::string hint;  // converter to use if the version in the file is older
};

int findVersion(TDirectory* dir, const std::string& base)
{
  int found = -1;
  TIter next(dir->GetListOfKeys());
  while (TKey* key = dynamic_cast<TKey*>(next())) {
    TString name = key->GetName();
    name.ToLower();
    if (name == base.c_str()) {
      found = std::max(found, 0);
    } else if (name.BeginsWith((base + "_").c_str())) {
      TString v = name(base.size() + 1, name.Length());
      if (v.IsDigit()) {
        found = std::max(found, v.Atoi());
      }
    }
  }
  return found;
}

void checkTables(const char* fname = "data/AO2D.root")
{
  const std::vector<RequiredTable> required = {
    {.base = "o2stracollision", .version = 0, .firstStep = -1, .hint = ""},
    {.base = "o2straevsels", .version = 6, .firstStep = -1, .hint = "o2-analysis-lf-straevselsconverterN (chain up to straevselsconverter6)"},
    {.base = "o2v0core", .version = 0, .firstStep = 0, .hint = ""},
    {.base = "o2v0extra", .version = 0, .firstStep = 2, .hint = ""},
    {.base = "o2dautrackextra", .version = 3, .firstStep = 2, .hint = "o2-analysis-lf-stradautracksextraconverterN (up to stradautracksextraconverter3)"},
    {.base = "o2dautracktpcpid", .version = 1, .firstStep = 2, .hint = "o2-analysis-lf-stradautrackstpcpidconverter"},
    // step 3: tables read as input by the helper task strangenesstofpid,
    // not by the analysis itself: the tutorial only uses V0s
    {.base = "o2strastamps", .version = 1, .firstStep = 3, .hint = "o2-analysis-lf-strastampsconverter"},
    {.base = "o2straevtimes", .version = 1, .firstStep = 3, .hint = "o2-analysis-lf-stradautrackstofpidconverter3"},
    {.base = "o2dautracktofpid", .version = 3, .firstStep = 3, .hint = "o2-analysis-lf-stradautrackstofpidconverterN (up to stradautrackstofpidconverter4)"},
    {.base = "o2v0collref", .version = 0, .firstStep = 3, .hint = ""},
    {.base = "o2casccore", .version = 0, .firstStep = 3, .hint = ""},
    {.base = "o2cascextra", .version = 0, .firstStep = 3, .hint = ""},
    {.base = "o2casccollref", .version = 0, .firstStep = 3, .hint = ""},
    // step 4: MC
    {.base = "o2v0coremclabel", .version = 0, .firstStep = 4, .hint = ""},
    {.base = "o2v0mccore", .version = 2, .firstStep = 4, .hint = "o2-analysis-lf-v0coresconverter / v0coresconverter2"},
    {.base = "o2v0mccollref", .version = 0, .firstStep = 4, .hint = ""},
    {.base = "o2stracolllabels", .version = 0, .firstStep = 4, .hint = ""},
    {.base = "o2stramccollision", .version = 2, .firstStep = 4, .hint = "o2-analysis-lf-stramccollisionconverter / stramccollisionconverter2"},
  };

  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    std::cout << "Cannot open " << fname << std::endl;
    return;
  }
  TDirectory* df = nullptr;
  TIter next(f->GetListOfKeys());
  while (TKey* key = dynamic_cast<TKey*>(next())) {
    if (TString(key->GetName()).BeginsWith("DF_")) {
      df = f->GetDirectory(key->GetName());
      break;
    }
  }
  if (!df) {
    std::cout << "No DF_* folder in " << fname << std::endl;
    return;
  }

  std::cout << "\nTrees in " << df->GetName() << ":\n";
  TIter nextTree(df->GetListOfKeys());
  while (TKey* key = dynamic_cast<TKey*>(nextTree())) {
    std::cout << "  " << key->GetName() << "\n";
  }

  int firstBrokenStep = 99;    // first step that cannot run, even with converters
  int firstConverterStep = 99; // first step that needs converters
  std::cout << "\nTables required by the tutorial:\n";
  for (const auto& t : required) {
    int v = findVersion(df, t.base);
    TString status = "ok";
    if (v < 0) {
      status = "MISSING";
      firstBrokenStep = std::min(firstBrokenStep, t.firstStep);
    } else if (v < t.version) {
      status = Form("v%d < v%d required", v, t.version);
      firstConverterStep = std::min(firstConverterStep, t.firstStep);
    } else if (v > t.version) {
      status = Form("v%d > v%d: O2Physics is older than the data, update it", v, t.version);
      firstBrokenStep = std::min(firstBrokenStep, t.firstStep);
    }
    TString step = t.firstStep < 0 ? "skeleton" : Form("step%d", t.firstStep);
    std::cout << Form("  %-20s (from %-8s) %s", t.base.c_str(), step.Data(), status.Data());
    if (status != "ok" && !t.hint.empty() && v >= 0) {
      std::cout << "\n        -> " << t.hint;
    }
    std::cout << "\n";
  }

  std::cout << "\nSummary:\n";
  for (int s = -1; s <= 4; s++) {
    TString name = s < 0 ? "skeleton" : Form("step%d", s);
    TString res = s >= firstBrokenStep      ? "CANNOT run with this file"
                  : s >= firstConverterStep ? "ok, but converters are needed"
                                            : "ok";
    std::cout << Form("  %-8s %s\n", name.Data(), res.Data());
  }
  std::cout << std::endl;
  f->Close();
}
