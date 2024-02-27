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
#ifndef PWGMM_RIVET_DATAMODEL_RIVETAOS_H_
#define PWGMM_RIVET_DATAMODEL_RIVETAOS_H_
#include <TObject.h>
#include <TList.h>
#include <TObjString.h>
#include <vector>
#include <string>
// #include <Mergers/MergeInterface.h>

// Forward declaration
namespace Rivet
{
class AnalysisHandler;
}

namespace o2
{
namespace rivet
{
/**  Container of output of Rivet Analyses.  Objects of this class
 *   contain a single string containing the YODA formatted @c
 *   Rivet::AnalysisObject.  This string can be written to file for
 *   further processing.  Objects of this class can be merged (as
 *   required by the @c TSelector protocol) using the facilities of @c
 *   Rivet::AnalysisHandler.  This allows us to parallelize Rivet jobs
 *   (e.g., on the Grid), merge the outputs, and then perform the @c
 *   Rivet::Analysis::finalize step in the normal terminate
 *   (`AliAnalysisTask::Terminate`) step.
 *
 *   @ingroup o2-rivet
 *
 *   @todo Probably need to derive from
 *   o2::mergers::MergeInterface to be mergable.
 */
class RivetAOs : public TObject //, public o2::mergers::MergeInterface
{
 public:
  /** Default constructor */
  RivetAOs();
  /** Copy constructor */
  RivetAOs(const RivetAOs& o);
  /** Destructor */
  virtual ~RivetAOs() {}

  /** Whether merging assumes equivalent output */
  void SetEquivalent(bool equiv) { mEquivalent = equiv; }
  /** Whether merging assumes equivalent output */
  bool GetEquivalent() const { return mEquivalent; }
  /** Add a load path */
  void ClearLoads() { mLoadPaths.Delete(); }
  /** Add a load path */
  void AddLoadPath(const char* path);
  /** Set used load paths */
  const TList& GetLoadPaths() const { return mLoadPaths; }
  /** Set used load paths */
  TList& GetLoadPaths() { return mLoadPaths; }

  /**
   * Read content of passed stream into our data member.
   *
   * Note, the chosen implementation isn't the most efficient but is
   * pretty straight forward (see
   * https://insanecoding.blogspot.com/2011/11/how-to-read-in-file-in-c.html).
   *
   * However, since this member function isn't called that often and
   * since the data we read in isn't humongous, we prefer the simpler
   * approach.
   *
   * @param in Input stream to read from
   */
  void ReadIn(std::istream& in);
  /** Writes our data to an output stream.
   *
   *  @param out Output stream
   */
  void WriteOut(std::ostream& out) const;
  /** Save data to a named file
   *
   *  @param filename Name of file to write to
   */
  void SaveAs(const char* filename = "Rivet.yoda",
              Option_t* = "") const override; // *MENU*
  /** Print the content of this object
   *
   *  @param options Not used
   */
  void Print(Option_t* options = "") const override; // *MENU*
  /*  Implementation of o2::mergers::MergeInterface
   *
   *  Note that this merges one object at a time, which isn't
   *  terribly efficient.  Throws exceptions in case of errors
   */
  // void merge(o2::mergers::MergeInterface* const other) override;
  /**
   * Merge objects of this class.  This is done by writing temporary
   * files of each object and then using @c
   * Rivet::AnalysisHandler::mergeYodas on those temporary files.
   *
   * Note, @c Rivet::AnalysisHandler::mergeYodas actually calls
   * finalize on all analyses, so they better be well-behaved for this
   * to work.
   *
   * @param collection Collection of objects to merge
   * @param options    If it contains "equiv" assume equivalent merge
   *
   * @return Number of merged objects or negative in case of errors
   */
  Long64_t Merge(TCollection* collection, Option_t* options);
  /**
   * Merge objects of this class.  This is done by writing temporary
   * files of each object and then using @c
   * Rivet::AnalysisHandler::mergeYodas on those temporary files.
   *
   * @param collection Collection of objects to merge
   *
   * @return Number of merged objects or negative in case of errors
   */
  Long64_t Merge(TCollection* collection);
  /**
   * Get the data stored
   */
  const TString& Data() const { return mData.GetString(); }

 protected:
  /**
   * Merge data in yoda files @a fnms.  The result of merging (and @c
   * Rivet::Analysis::finalize) is stored in this object.
   *
   * @param ah Rivet::AnalysisHandler to use
   * @param fnms Vector of filenames (strings)
   * @param equiv Merge assuming equivalence of data
   *
   * @return Number of files merged or negative in case of error
   */
  Long_t MergeYoda(Rivet::AnalysisHandler& ah,
                   const std::vector<std::string>& fnms,
                   bool equiv,
                   bool clean = false);
  /**
   * Write content to a tempoary file and return the file name written
   * to.  This is public because it will be called from @c RivetTask.
   *
   * @return Filename of the temporary file
   */
  std::string WriteTmp() const;

  TObjString mData; // The data
  TList mLoadPaths; // Load paths used
  bool mEquivalent = true;

  ClassDefOverride(RivetAOs, 1);
};
} // namespace rivet
} // namespace o2
#endif // PWGMM_RIVET_DATAMODEL_RIVETAOS_H_
//
// EOF
//
