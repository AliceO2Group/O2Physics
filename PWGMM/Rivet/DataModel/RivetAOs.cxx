// Copyright 2023-2099 CERN and copyright holders of ALICE O2.
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
#include "RivetAOs.h"
#include <YODA/IO.h>
#include <TSystem.h>
#include <TClass.h>
#include <vector>
#include <string>
#include <Rivet/AnalysisHandler.hh>
#include <Rivet/Tools/RivetPaths.hh>
// Sigh - linter needs additional headers already included

namespace o2
{
namespace rivet
{
//----------------------------------------------------------------
RivetAOs::RivetAOs()
  : mData("")
{
  mLoadPaths.SetOwner(true);
}
//----------------------------------------------------------------
RivetAOs::RivetAOs(const RivetAOs& o)
  : mData(o.mData),
    mEquivalent(o.mEquivalent)
{
  mLoadPaths.SetOwner(true);
}
//----------------------------------------------------------------
void RivetAOs::AddLoadPath(const char* path)
{
  mLoadPaths.Add(new TObjString(path));
}
//----------------------------------------------------------------
void RivetAOs::ReadIn(std::istream& in)
{
  mData.String() = std::string((std::istreambuf_iterator<char>(in)),
                               std::istreambuf_iterator<char>());
}
//----------------------------------------------------------------
void RivetAOs::WriteOut(std::ostream& out) const
{
  out << mData.GetString() << std::flush;
}
//----------------------------------------------------------------
void RivetAOs::SaveAs(const char* filename, Option_t*) const
{
  std::ofstream out(filename);
  WriteOut(out);
}
//----------------------------------------------------------------
void RivetAOs::Print(Option_t*) const
{
  std::cout << "=== Rivet Analysis Objects:\n";
  WriteOut(std::cout);
  std::cout << "=== " << std::endl;
}
#if 0
//----------------------------------------------------------------
void RivetAOs::merge(o2::mergers::MergeInterface* const other)
{
  Info("merge", "O2 merge interface member function called");
  RivetAOs* o = dynamic_cast<RivetAOs*>(other);
  if (!o) { // Sigh - annoying linter prefers `!` over `not`
    throw std::runtime_error("Cannot merge with object not "
                             "of class RivetAOs");
  }
  std::vector<std::string> tmps;
  // Sigh - annoying linter prefers `!` over `not`
  if (!mData.GetString().IsNull()) {
    tmps.push_back(WriteTmp());
  }
  // Sigh - annoying linter prefers `!` over `not`
  if (!o->Data().IsNull()) {
    tmps.push_back(o->WriteTmp());
  }

  for (auto f : tmps) {
    Info("Merge", "  %s", f.c_str());
    if (gSystem->AccessPathName(f.c_str())) {
      Warning("Merge", "Does not exist!");
    }
  }

  for (auto o : mLoadPaths) {
    Rivet::addAnalysisLibPath(o->GetName());
    Rivet::addAnalysisDataPath(o->GetName());
  }

  Rivet::AnalysisHandler ah("");
  if (MergeYoda(ah, tmps, mEquivalent, true) <= 0)
    throw std::runtime_error("Failed to merge YODA outputs");
}
#endif
//----------------------------------------------------------------
Long64_t RivetAOs::Merge(TCollection* collection, Option_t* options)
{
  // Info("Merge", "ROOT merge interface member function called");
  std::vector<std::string> tmps;
  bool equiv = TString(options).Contains("equiv");

  try {
    if (!mData.GetString().IsNull()) {
      tmps.push_back(WriteTmp());
    }

    for (auto* o : *collection) {
      if (!o->IsA()->InheritsFrom(RivetAOs::Class())) {
        Error("Merge", "Cannot merge Output with object of class %s",
              o->ClassName());
        return -1;
      }

      RivetAOs* ro = dynamic_cast<RivetAOs*>(o);
      tmps.push_back(ro->WriteTmp());

      mLoadPaths.Merge(&ro->mLoadPaths);
    }
  } catch (...) {
    Error("Merge", "Failed to write some temporary files");
    return -1;
  }
  Info("Merge", "Will merge content of ");
  for (auto f : tmps) {
    Info("Merge", "  %s", f.c_str());
    if (gSystem->AccessPathName(f.c_str())) {
      Warning("Merge", "Does not exist!");
    }
  }

  for (auto o : mLoadPaths) {
    Rivet::addAnalysisLibPath(o->GetName());
    Rivet::addAnalysisDataPath(o->GetName());
  }

  Rivet::AnalysisHandler ah("");
  return MergeYoda(ah, tmps, equiv, true);
}
//----------------------------------------------------------------
Long64_t RivetAOs::Merge(TCollection* collection)
{
  return Merge(collection, mEquivalent ? "equiv" : "");
}
//----------------------------------------------------------------
Long_t RivetAOs::MergeYoda(Rivet::AnalysisHandler& ah,
                           const std::vector<std::string>& fnms,
                           bool equiv,
                           bool clean)
{
  Info("MergeYoda", "Merging %zu files w/YODA AOs with equiv=%d",
       fnms.size(), equiv);
  if (fnms.size() <= 0) {
    Warning("MergeYoda", "No files to merge");
    return 0;
  }
  try {
#if RIVET_VERSION_CODE >= 30104
    std::vector<std::string> dels;
    std::vector<std::string> adds;
    std::vector<std::string> match;
    std::vector<std::string> umatch;
    ah.mergeYodas(fnms, dels, adds, match, umatch, equiv);
#else
    std::vector<std::string> dels;
    std::vector<std::string> adds;
    ah.mergeYodas(fnms, dels, adds, equiv);
#endif
  } catch (std::exception& e) {
    Fatal("Merge", "Failed to merge data: %s", e.what());
    return -1;
  }

  std::stringstream out;
  const std::vector<YODA::AnalysisObjectPtr> output = ah.getYodaAOs(true);
  try {
    YODA::write(out, output.begin(), output.end(), "yoda");
  } catch (std::exception& e) {
    Error("Merge", "failed to write out data: %s", e.what());
    return -1;
  } catch (...) { //< YODA::WriteError&
    Fatal("Merge", "Failed to write merged objects");
    return -1;
  }
  ReadIn(out);

  if (clean) {
    for (auto fn : fnms) {
      // Continue if the file isn't there
      if (gSystem->AccessPathName(fn.c_str()))
        continue;
      Info("MergeYoda", "Removing temporary file %s", fn.c_str());
      gSystem->Unlink(fn.c_str());
    }
  }
  return fnms.size();
}
//----------------------------------------------------------------
std::string RivetAOs::WriteTmp() const
{
  TString base("rivetoutput");
  FILE* fp = gSystem->TempFileName(base, ".");

  int wn = fwrite(mData.GetString().Data(), 1,
                  mData.GetString().Length(), fp);
  if (wn != mData.GetString().Length())
    throw std::runtime_error("Bad write");

  fclose(fp);

  // Rename file so that it has ending ".yoda" - otherwise YODA
  // cannot figure out what to do - sigh!
  std::string ret(base.Data());
  ret += ".yoda";
  gSystem->Rename(base.Data(), ret.c_str());

  return ret;
}
} // namespace rivet
} // namespace o2
//
// EOF
//
