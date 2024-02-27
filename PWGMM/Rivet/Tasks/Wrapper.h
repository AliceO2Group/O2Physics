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
#ifndef PWGMM_RIVET_TASKS_WRAPPER_H_
#define PWGMM_RIVET_TASKS_WRAPPER_H_
#include "RivetAOs.h"
// Linter wants a bunch of header files already included from elsewere
// - really silly.
#include <YODA/IO.h>
#include <Framework/AnalysisManagers.h>
#include <HepMC3/GenEvent.h>
#include <algorithm>
#include <list>
#include <string>
#include <memory>
#include <vector>
#include <map>
#include <filesystem>
#include <regex>
#include <Rivet/AnalysisHandler.hh>
#include <Rivet/Tools/RivetPaths.hh>

namespace o2
{
namespace rivet
{
/** Import here so we do not need to use fully qualified name */
template <typename T>
using Configurable = o2::framework::Configurable<T>;

/** Wrapper around Rivet */
struct Wrapper {
  /** configurations */
  struct : o2::framework::ConfigurableGroup {
    /** @{
        @name Settings */
    Configurable<double> crossSection{"rivet-cross-section", -1,
                                      "Override EG cross-section"};
    Configurable<bool> mergeEquiv{"rivet-merge-equivalent", true,
                                  "Whether merging assume "
                                  "equivalent output"};
    Configurable<bool> ignoreBeams{"rivet-ignore-beams", false,
                                   "Ignore beam requirements"};
    Configurable<bool> pwd{"rivet-pwd", false,
                           "Add current directory load path"};
    Configurable<bool> finalize{"rivet-finalize", true,
                                "Run Rivet::Analysis::finalize on output"};
    Configurable<std::string> dump{"rivet-dump", "",
                                   "Dump YODA objects to disk"};
    Configurable<std::string> anas{"rivet-analysis", "",
                                   "Comma separated list of analyses. "
                                   "An analysis name may be followed by "
                                   "a list of options as "
                                   "[:key=value]*"};
    Configurable<std::string> paths{"rivet-load-paths", "",
                                    "Colon separated list of load paths. "
                                    "Rivet will look for plugins (analyses) "
                                    "and data files in these directories"};
    Configurable<std::string> pres{"rivet-pre-loads", "",
                                   "Comma separated list of preloads. "
                                   "Rivet will load these YODA files "
                                   "before running the analyses"};
    Configurable<std::string> srcs{"rivet-sources", "",
                                   "Comma separated list of sources. "
                                   "Each source will be compiled into a "
                                   "Rivet plugin before execution"};
    Configurable<std::string> flags{"rivet-flags", "",
                                    "Extra compiler flags when compiling "
                                    "plugins"};
    Configurable<std::string> log{"rivet-log", "WARNING",
                                  "Set Rivet logging level "
                                  "(TRACE,DEBUG,INFO,WARN,ERROR,CRITICAL)"};
  } configs;

  /** Type of a HepMC3 event */
  using Event = HepMC3::GenEvent;
  /** Type of Rivet analysis handler */
  using Handler = Rivet::AnalysisHandler;
  /** Type of pointer to Rivet analysis handler */
  using HandlerPtr = std::shared_ptr<Handler>;
  /** Pointer to output object */
  using RivetAOsPtr = std::shared_ptr<RivetAOs>;
  /** Our analysis handler */
  HandlerPtr mHandler;
  /** Type of a list of strings */
  using StringList = std::list<std::string>;
  /** List of analyses to run, including options */
  StringList mAnalyses;
  /** List of load paths */
  StringList mLoadPaths;
  /** List of sources to compile  */
  StringList mSources;
  /** List of data files to preload  */
  StringList mPreLoads;
  /** Pointer to output object */
  RivetAOsPtr mOutput;
  /** Run finalize on output */
  bool mFinalize = false;

  /** Append elements in string @a s, separated by character @a
   * sep, to the list @a l
   *
   * @param l   List to append to
   * @param s   String of elements to append to @a l
   * @param sep Separator of elements in @a s
   */
  void appendToList(StringList& l,
                    Configurable<std::string>& s,
                    char sep)
  {
    std::stringstream str(s.value);
    std::string part;
    while (std::getline(str, part, sep))
      l.push_back(part);
  }
  /** Format a list @a l with separator @a sep
   *
   * @param l   List to format
   * @param sep Separator used in format
   * @return String-formatted list
   */
  std::string formatList(const std::string& heading,
                         const StringList& l) const
  {
    if (l.size() <= 0) {
      return heading;
    }

    std::stringstream s;
    s << heading;
    std::string sep = "\n" + std::string(heading.length(), ' ');

    std::copy(l.begin(), l.end(),
              std::ostream_iterator<std::string>(s, sep.c_str()));
    auto ret = s.str();
    ret.erase(ret.length() - sep.length());

    return ret;
  }
  /** Initialize this.  This sets
   *
   * - the analyses to run
   * - the loads paths to use
   * - the pre-loaded data to read
   * - the sources  to compile
   *
   * and builds the sources into plug-ins
   *
   * @param output Shared pointer to the output object
   */
  void init(RivetAOsPtr output)
  {
    appendToList(mAnalyses, configs.anas, ',');
    appendToList(mLoadPaths, configs.paths, ':');
    appendToList(mPreLoads, configs.pres, ',');
    appendToList(mSources, configs.srcs, ',');
    if (configs.pwd.value) {
      mLoadPaths.push_front(std::filesystem::current_path());
    }

    std::string logLvl = configs.log;

    LOG(info) << "=== o2::rivet::Wrapper ===\n"
              << formatList("Analyses:   ", mAnalyses) << "\n"
              << formatList("Load paths: ", mLoadPaths) << "\n"
              << "PWD:        " << configs.pwd.value << "\n"
              << formatList("Pre-loads:  ", mPreLoads) << "\n"
              << formatList("Sources:    ", mSources) << "\n"
              << "Logging:    " << logLvl + "\n";

    mOutput = output;
    mFinalize = configs.finalize;

    setLogLevel(logLvl);
    initLoadPaths();
    buildPlugins();
  }
  /** Utility function to expand environment variables in paths,
      etc. Note, environment variables _must_ by given like `${NAME}`
      - i.e., the initial dollar sign `$` _and_ start and end braces
      `{}' _must_ be given. */
  static void expandEnvVars(std::string& src)
  {
    static std::regex re("\\$\\{([^}]*)\\}");
    std::smatch match;

    while (std::regex_search(src, match, re)) {
      auto eval = std::getenv(match[1].str().c_str());
      std::string val = (eval ? eval : "");
      src.replace(match[0].first, match[0].second, val);
    }
  }
  /** Add paths, for binary code and data, to the Rivet search paths
   */
  void initLoadPaths()
  {
    // Add O2Physics paths ($O2PHYSICS_ROOT/lib and
    // $O2PHYSICS_ROOT/share/O2Physics/Rivet) to load paths.
    auto o2p = std::getenv("O2PHYSICS_ROOT");
    if (!o2p) {
      // Possibly try O2 installation too
      o2p = std::getenv("O2_ROOT");
    }
    if (o2p) {
      std::string so2p(o2p);
      LOG(info) << "Adding O2Physics paths at " << so2p;
      mLoadPaths.push_back(so2p + "/lib");
      mLoadPaths.push_back(so2p + "/share/O2Physics/Rivet");
    }
    // Alternatively, one could use the location of the application
    // itself to set the appropriate paths - see
    //
    // https://github.com/gpakosz/whereami/blob/master/src/whereami.c
    for (auto p : mLoadPaths) {
      expandEnvVars(p);
      LOG(info) << "Adding " << p << " to load path";
      Rivet::addAnalysisLibPath(p);
      Rivet::addAnalysisDataPath(p);
      mOutput->AddLoadPath(p.c_str());
    }
  }
  /** Buid analyses plugins */
  void buildPlugins()
  {
    for (auto s : mSources) {
      expandEnvVars(s);
      LOG(info) << "Compiling " << s;

      // Use filesystem interface
      std::filesystem::path p = s;

      // Build command
      std::stringstream c;
      c << "rivet-build Rivet" << p.stem().string() << ".so "
        << s << " " << configs.flags.value;
      LOG(info) << "Compilation command: \"" << c.str() << "\"";
      int ret = std::system(c.str().c_str());
      if (ret != 0) {
        throw std::runtime_error("Failed to combile " + s);
      }
    }
  }
  /** Find numeric log-level correspondinf to string log level.
   *
   * @param lvl String value
   * @return log level correspondig to lvl or -1
   */
  int findLogLevel(std::string const& lvl) const
  {
    const std::map<std::string, int> lvls = {
      {"TRACE", Rivet::Log::TRACE},
      {"DEBUG", Rivet::Log::DEBUG},
      {"INFO", Rivet::Log::INFO},
      {"WARN", Rivet::Log::WARN},
      {"WARNING", Rivet::Log::WARN},
      {"ERROR", Rivet::Log::ERROR},
      {"CRITICAL", Rivet::Log::CRITICAL},
      {"FATAL", Rivet::Log::CRITICAL}};

    std::string logLvl = lvl;
    std::transform(logLvl.begin(), logLvl.end(), logLvl.begin(),
                   [](char c) { return std::toupper(c); });

    auto iter = lvls.find(lvl);
    if (iter == lvls.end()) {
      return -1;
    }
    return iter->second;
  }

  /** Set log level of rivet */
  void setLogLevel(std::string const& lvl)
  {
    int logVal = findLogLevel(lvl);
    if (logVal < 0) {
      return; // No change
    }
    Rivet::Log::setLevel("Rivet.AnalysisHandler", logVal);
  }
  /** Initialize the Rivet hander
   *
   * This registers the analyses to run and loads the pre-loaded
   * data.
   */
  void initHandler()
  {
    LOG(info) << "Creating analysis handler";

    mHandler = std::make_shared<Handler>();
    mHandler->setIgnoreBeams(configs.ignoreBeams.value);
    if (configs.crossSection.value > 0) {
      mHandler->setCrossSection(configs.crossSection.value, 0, true);
    }

    for (auto a : mAnalyses) {
      LOG(info) << "Adding analysis " << a;
      mHandler->addAnalysis(a);
    }

    for (auto d : mPreLoads) {
      expandEnvVars(d);
      LOG(info) << "Adding preloaded data " << d;
      mHandler->readData(d);
    }
  }
  /** Process an event.  The @c HepMC3::GenEvent passed is passed
   * on to Rivet, and Rivet then runs the analyses @c event
   * methods on that event.
   *
   * @param event A @c HepMC3::GenEvent object to analyse.
   */
  void process(const Event& event)
  {
    if (!mHandler) { // I prefer to use `not` instead of `!` but
                     // linter apparently likes to use the less-clear
                     // operator - sigh!
      initHandler();
      mHandler->init(event);
    }
    if (mAnalyses.size() <= 0) {
      LOG(warning) << "No analysis run" << std::endl;
    }
    mHandler->analyze(event);
  }
  /** At end of a run.
   *
   * This will get the YODA analysis objects (AOs) and possibly
   * run finalize on them for the registerd analyses.  It can
   * also output the AOs to disk, if so instructed.
   */
  void postRun()
  {
    if (!mHandler) { // I prefer to use `not` instead of `!` but
                     // linter apparently likes to use the less-clear
                     // operator - sigh!
      LOG(warning) << "No handler, nothing to do";
      return;
    }

    if (mOutput) {
      // Get analysis objects.  True argument means also get raw
      // objects too.
      auto aos = mHandler->getYodaAOs(true);

      // Create a string stream and write Yoda's to that and then
      // stream those to the output data object.
      std::stringstream out;
      YODA::write(out, aos.begin(), aos.end(), "yoda");

      mOutput->ReadIn(out);
      mOutput->SetEquivalent(configs.mergeEquiv.value);
    }

    if (configs.dump->empty()) {
      return;
    }

    // Write output to disk
    LOG(info) << "Writing to \"" << configs.dump.value << "\""
              << (mFinalize ? " including finalizsed objects" : "");
    if (mFinalize) {
      mHandler->finalize();
    }
    auto aos = mHandler->getYodaAOs(true);
    YODA::write(configs.dump.value, aos.begin(), aos.end());
  }
};
} // namespace rivet
// Create specialisation for output manager
namespace framework
{
/** This specialisation of @c o2::framework::OutputManager ensures
 * that we can call the post-processing routine of @c
 * o2::rivet::Wrapper and thus ensure that the YODA string is
 * written to the output object.
 *
 * The O2 framework (via @c o2::framework::adoptAnalysisTask<T>)
 * inspects the members of the pass class (@c T) and creates @c
 * o2::framework::OutputManager callbacks for every member.  The
 * default template for this does nothing.
 *
 * Thus, to delegate a call to a member of the analysis task (of
 * class @c T), we can specialise the @c
 * o2::framework::OutputManager template on the member type.  We
 * will then effectively have call-backs for
 *
 * - @c appendOutput - when the task is constructed
 * - @c prepare - when a new set of data is recieved
 * - @c finalize - when a set of data has been processed
 * - @c postRun - when the run is over
 *
 * Concretely, we use the @c postRun set the YODA string on the
 * output object
 * */

template <>
struct OutputManager<o2::rivet::Wrapper> {
  using Target = o2::rivet::Wrapper;
  static bool appendOutput(std::vector<OutputSpec>&,
                           Target&,
                           uint32_t)
  {
    // LOG(info) << "\n"
    //    << "**************************************\n"
    //    << "Append output for o2::rivet::Wrapper\n"
    //    << "**************************************";
    return true;
  }
  static bool prepare(ProcessingContext&, Target&)
  {
    // LOG(info) << "\n"
    //    << "*****************************************\n"
    //    << "Prepare for data for o2::rivet::Wrapper\n"
    //    << "*****************************************";
    return true;
  }
  static bool postRun(EndOfStreamContext&, Target& t)
  {
    LOG(info) << "\n"
              << "****************************************\n"
              << "Post run output for o2::rivet::Wrapper\n"
              << "****************************************";
    t.postRun();
    return true;
  }
  static bool finalize(ProcessingContext&, Target& t)
  {
    // LOG(info) << "\n"
    //    << "**************************************\n"
    //    << "Finalize data for o2::rivet::Wrapper\n"
    //    << "**************************************";
    return true;
  }
};

/** Specialisation to pull in configurables from the wrapper
 *
 * Ideally, the wrapper should simply derive from
 * ConfigurableGroup and all should flow automatically, but that
 * doesn't work for some reason.
 */
template <>
struct OptionManager<o2::rivet::Wrapper> {
  using Target = o2::rivet::Wrapper;
  static bool
    appendOption(std::vector<o2::framework::ConfigParamSpec>& options,
                 Target& target)
  {
    // LOG(info) << "Appending wrapper options";
    OptionManager<ConfigurableGroup>::appendOption(options, target.configs);
    return true;
  }
  static bool
    prepare(o2::framework::InitContext& ic, Target& target)
  {
    // LOG(info) << "Preparing wrapepr options";
    OptionManager<ConfigurableGroup>::prepare(ic, target.configs);
    return true;
  }
};
} // namespace framework

} // namespace o2

#endif // PWGMM_RIVET_TASKS_WRAPPER_H_
