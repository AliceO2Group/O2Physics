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

#include <regex>
#include <TObjArray.h>

#include <fairlogger/Logger.h>
#include "PWGCF/TwoParticleCorrelations/Core/SkimmingConfigurableCuts.h"

using namespace o2;
using namespace o2::analysis::PWGCF;

/// \brief These are the implemented basic bricks
/// If more are implemented the list must be expanded and the
/// corresponding brick construction implemented
template <typename TValueToFilter>
const char* CutBrick<TValueToFilter>::mgImplementedbricks[] = {"lim", "fnlim", "th", "fnth", "rg", "fnrg", "xrg", "fnxrg", "mrg", "cwv"};

///////////////////////////////////////////////////////////////////////////////////////
template <typename TValueToFilter>
CutBrick<TValueToFilter>* CutBrick<TValueToFilter>::constructBrick(const char* name, const char* regex, const std::set<std::string>& allowed)
{
  LOGF(info, "Construct brick %s from RE %s", name, regex);
  CutBrick<TValueToFilter>* thebrick = nullptr;

  bool found = false;
  for (auto& bname : allowed) {
    if (TString(regex).BeginsWith(bname)) {
      found = true;
      break;
    }
  }
  if (not found) {
    LOGF(fatal, "CutBrick<TValueToFilter>* CutBrick<TValueToFilter>::constructBrick", "Wrong RE: %s, trying to construct a not allowed basic cut brick", regex);
  }
  TString brickregex = TString(name) + "{" + TString(regex) + "}";

  if (TString(regex).BeginsWith("lim")) {
    thebrick = new CutBrickLimit<TValueToFilter>(brickregex);
  } else if (TString(regex).BeginsWith("fnlim")) {
    thebrick = new CutBrickFnLimit<TValueToFilter>(brickregex);
  } else if (TString(regex).BeginsWith("th")) {
    thebrick = new CutBrickThreshold<TValueToFilter>(brickregex);
  } else if (TString(regex).BeginsWith("fnth")) {
    thebrick = new CutBrickFnThreshold<TValueToFilter>(brickregex);
  } else if (TString(regex).BeginsWith("rg")) {
    thebrick = new CutBrickRange<TValueToFilter>(brickregex);
  } else if (TString(regex).BeginsWith("fnrg")) {
    thebrick = new CutBrickFnRange<TValueToFilter>(brickregex);
  } else if (TString(regex).BeginsWith("xrg")) {
    thebrick = new CutBrickExtToRange<TValueToFilter>(brickregex);
  } else if (TString(regex).BeginsWith("fnxrg")) {
    thebrick = new CutBrickFnExtToRange<TValueToFilter>(brickregex);
  } else if (TString(regex).BeginsWith("mrg")) {
    thebrick = new CutBrickSelectorMultipleRanges<TValueToFilter>(brickregex);
  } else if (TString(regex).BeginsWith("cwv")) {
    thebrick = new CutWithVariations<TValueToFilter>(brickregex);
  } else {
    LOGF(fatal, "CutBrick<TValueToFilter>* CutBrick<TValueToFilter>::constructBrick", "Wrong RE: %s, trying to construct an unknown basic cut brick", regex);
  }
  return thebrick;
}

/// Default constructor
template <typename TValueToFilter>
CutBrick<TValueToFilter>::CutBrick()
  : TNamed(),
    mState(kPASSIVE),
    mMode(kUNSELECTED)
{
}

/// Named constructor
/// \param name The name of the brick
template <typename TValueToFilter>
CutBrick<TValueToFilter>::CutBrick(const char* name, const char* title)
  : TNamed(name, title),
    mState(kPASSIVE),
    mMode(kUNSELECTED)
{
}

templateClassImp(CutBrick);

///////////////////////////////////////////////////////////////////////////////////////
/// Default constructor
template <typename TValueToFilter>
CutBrickLimit<TValueToFilter>::CutBrickLimit()
  : CutBrick<TValueToFilter>(),
    mLimit(TValueToFilter(0))
{
}

/// Named constructor
/// \param name The name of the brick
/// \param value The limit value
template <typename TValueToFilter>
CutBrickLimit<TValueToFilter>::CutBrickLimit(const char* name, const TValueToFilter& value)
  : CutBrick<TValueToFilter>(name, TString::Format("%s{lim{%f}}", name, float(value))),
    mLimit(value)
{
}

/// \brief Cut string constructor
/// \param cutstr The cuts string
template <typename TValueToFilter>
CutBrickLimit<TValueToFilter>::CutBrickLimit(const TString& cutstr)
  : CutBrick<TValueToFilter>(),
    mLimit(TValueToFilter(0))
{
  ConstructCutFromString(cutstr);
}

/// \brief Construct the cut from a cut string
/// \param cutstr The cut string
/// The cut string should have the structure
///    name{lim{val}}
/// If the cut string is correctly parsed the cut is correctly built
/// if not a fatal exception is raised
template <typename TValueToFilter>
void CutBrickLimit<TValueToFilter>::ConstructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());

  std::regex cutregex("^(\\w+)\\{lim\\{((?:-?[\\d]+\\.?[\\d]*)|(?:-?[\\d]*\\.?[\\d]+))}}$", std::regex_constants::ECMAScript | std::regex_constants::icase);
  std::string in(cutstr.Data());
  std::smatch m;

  std::regex_search(in, m, cutregex);
  if (m.empty() or (m.size() < 3)) {
    LOGF(fatal, "CutBrickLimit<TValueToFilter>::ConstructCutFromString", "Wrong RE: %s, use pT{lim{2.0}} for instance", cutstr.Data());
  } else {
    this->SetName(m[1].str().c_str());
    this->SetTitle(cutstr.Data());
    mLimit = TValueToFilter(std::stod(m[2]));
  }
}

/// \brief Returns wether the cut brick is incorporated in the selection chain
/// \return true if the cut brick is incorporated
template <typename TValueToFilter>
std::vector<bool> CutBrickLimit<TValueToFilter>::IsArmed()
{
  std::vector<bool> res;
  if (this->mMode == this->kSELECTED) {
    res.push_back(true);
  } else {
    res.push_back(false);
  }
  return res;
}

templateClassImp(CutBrickLimit);
template class o2::analysis::PWGCF::CutBrickLimit<int>;
template class o2::analysis::PWGCF::CutBrickLimit<float>;

///////////////////////////////////////////////////////////////////////////////////////
/// Default constructor
template <typename TValueToFilter>
CutBrickFnLimit<TValueToFilter>::CutBrickFnLimit()
  : CutBrickLimit<TValueToFilter>(),
    mFunction{}
{
}

/// Named constructor
/// \param name The name of the brick
/// \param fn The function which will provide limit value
template <typename TValueToFilter>
CutBrickFnLimit<TValueToFilter>::CutBrickFnLimit(const char* name, const TF1& fn)
  : CutBrickLimit<TValueToFilter>(),
    mFunction(fn)
{
  this->SetName(name);
}

/// \brief Cut string constructor
/// \param cutstr The cuts string
template <typename TValueToFilter>
CutBrickFnLimit<TValueToFilter>::CutBrickFnLimit(const TString& cutstr)
  : CutBrickLimit<TValueToFilter>(),
    mFunction{}
{
  ConstructCutFromString(cutstr);
}

/// \brief Construct the cut from a cut string
/// \param cutstr The cut string
/// The cut string should have the structure
///    name{fnlim{function}}
/// If the cut string is correctly parsed the cut is correctly built
/// if not a fatal exception is raised
template <typename TValueToFilter>
void CutBrickFnLimit<TValueToFilter>::ConstructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());

  std::regex cutregex("^(\\w+)\\{fnlim\\{(\\w+)=([\\w\\-\\*\\/\\+\\()\\.]+)}}$", std::regex_constants::ECMAScript | std::regex_constants::icase);
  std::string in(cutstr.Data());
  std::smatch m;

  std::regex_search(in, m, cutregex);
  if (m.empty() or (m.size() < 4)) {
    LOGF(fatal, "CutBrickFnLimit<TValueToFilter>::ConstructCutFromString", "Wrong RE: %s, use pT{fnlim{myfn=2.0*sin(x)/x}} for instance", cutstr.Data());
  } else {
    this->SetName(m[2].str().c_str());
    this->SetTitle(cutstr.Data());
    mFunction = TF1(m[1].str().c_str(), m[3].str().c_str(), 0, 1, TF1::EAddToList::kNo);
    if (not mFunction.IsValid()) {
      LOGF(fatal, "CutBrickFnLimit<TValueToFilter>::ConstructCutFromString", "Wrong function expression: %s, use pT{fnlim{myfn=2.0*sin(x)/x}} for instance", cutstr.Data());
    }
  }
}

templateClassImp(CutBrickFnLimit);
template class o2::analysis::PWGCF::CutBrickFnLimit<int>;
template class o2::analysis::PWGCF::CutBrickFnLimit<float>;

///////////////////////////////////////////////////////////////////////////////////////
/// Default constructor
template <typename TValueToFilter>
CutBrickThreshold<TValueToFilter>::CutBrickThreshold()
  : CutBrick<TValueToFilter>(),
    mThreshold(0)
{
}

/// Named constructor
/// \param name The name of the brick
/// \param value The threshold value
template <typename TValueToFilter>
CutBrickThreshold<TValueToFilter>::CutBrickThreshold(const char* name, const TValueToFilter& value)
  : CutBrick<TValueToFilter>(name, TString::Format("%s{th{%f}}", name, float(value))),
    mThreshold(value)
{
}

/// \brief Cut string constructor
/// \param cutstr The cuts string
template <typename TValueToFilter>
CutBrickThreshold<TValueToFilter>::CutBrickThreshold(const TString& cutstr)
  : CutBrick<TValueToFilter>(),
    mThreshold(0)
{
  ConstructCutFromString(cutstr);
}

/// \brief Construct the cut from a cut string
/// \param cutstr The cut string
/// The cut string should have the structure
///    name{th{val}}
/// If the cut string is correctly parsed the cut is correctly built
/// if not a fatal exception is raised
template <typename TValueToFilter>
void CutBrickThreshold<TValueToFilter>::ConstructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());

  std::regex cutregex("^(\\w+)\\{th\\{((?:-?[\\d]+\\.?[\\d]*)|(?:-?[\\d]*\\.?[\\d]+))}}$", std::regex_constants::ECMAScript | std::regex_constants::icase);
  std::string in(cutstr.Data());
  std::smatch m;

  std::regex_search(in, m, cutregex);
  if (m.empty() or (m.size() < 3)) {
    LOGF(fatal, "CutBrickThreshold<TValueToFilter>::ConstructCutFromString", "Wrong RE: %s, use pT{th{0.2}} for instance", cutstr.Data());
  } else {
    this->SetName(m[1].str().c_str());
    this->SetTitle(cutstr.Data());
    mThreshold = TValueToFilter(std::stod(m[2]));
  }
}

/// \brief Returns wether the cut brick is incorporated in the selection chain
/// \return true if the cut brick is incorporated
template <typename TValueToFilter>
std::vector<bool> CutBrickThreshold<TValueToFilter>::IsArmed()
{
  std::vector<bool> res;
  if (this->mMode == this->kSELECTED) {
    res.push_back(true);
  } else {
    res.push_back(false);
  }
  return res;
}

templateClassImp(CutBrickThreshold);
template class o2::analysis::PWGCF::CutBrickThreshold<int>;
template class o2::analysis::PWGCF::CutBrickThreshold<float>;

///////////////////////////////////////////////////////////////////////////////////////
/// Default constructor
template <typename TValueToFilter>
CutBrickFnThreshold<TValueToFilter>::CutBrickFnThreshold()
  : CutBrickThreshold<TValueToFilter>(),
    mFunction{}
{
}

/// Named constructor
/// \param name The name of the brick
/// \param fn The function which will provide the threshold value
template <typename TValueToFilter>
CutBrickFnThreshold<TValueToFilter>::CutBrickFnThreshold(const char* name, const TF1& fn)
  : CutBrickThreshold<TValueToFilter>(),
    mFunction(fn)
{
  this->SetName(name);
}

/// \brief Cut string constructor
/// \param cutstr The cuts string
template <typename TValueToFilter>
CutBrickFnThreshold<TValueToFilter>::CutBrickFnThreshold(const TString& cutstr)
  : CutBrickThreshold<TValueToFilter>(),
    mFunction{}
{
  ConstructCutFromString(cutstr);
}

/// \brief Construct the cut from a cut string
/// \param cutstr The cut string
/// The cut string should have the structure
///    name{fnth{function}}
/// If the cut string is correctly parsed the cut is correctly built
/// if not a fatal exception is raised
template <typename TValueToFilter>
void CutBrickFnThreshold<TValueToFilter>::ConstructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());

  std::regex cutregex("^(\\w+)\\{fnth\\{(\\w+)=([\\w\\-\\*\\/\\+\\()\\.]+)}}$", std::regex_constants::ECMAScript | std::regex_constants::icase);
  std::string in(cutstr.Data());
  std::smatch m;

  std::regex_search(in, m, cutregex);
  if (m.empty() or (m.size() < 4)) {
    LOGF(fatal, "CutBrickFnThreshold<TValueToFilter>::ConstructCutFromString", "Wrong RE: %s, use pT{fnth{myfn=2.0*sin(x)/x}} for instance", cutstr.Data());
  } else {
    this->SetName(m[2].str().c_str());
    this->SetTitle(cutstr.Data());
    mFunction = TF1(m[1].str().c_str(), m[3].str().c_str(), 0, 1, TF1::EAddToList::kNo);
    if (not mFunction.IsValid()) {
      LOGF(fatal, "CutBrickFnThreshold<TValueToFilter>::ConstructCutFromString", "Wrong function expression: %s, use pT{fnth{myfn=2.0*sin(x)/x}} for instance", cutstr.Data());
    }
  }
}

templateClassImp(CutBrickFnThreshold);
template class o2::analysis::PWGCF::CutBrickFnThreshold<int>;
template class o2::analysis::PWGCF::CutBrickFnThreshold<float>;

///////////////////////////////////////////////////////////////////////////////////////
/// Default constructor
template <typename TValueToFilter>
CutBrickRange<TValueToFilter>::CutBrickRange()
  : CutBrick<TValueToFilter>(),
    mLow(0),
    mUp(0)
{
}

/// Named constructor
/// \param name The name of the brick
/// \param low The low value for the cut range
/// \param high The high value for the cut range
template <typename TValueToFilter>
CutBrickRange<TValueToFilter>::CutBrickRange(const char* name, const TValueToFilter& low, const TValueToFilter& high)
  : CutBrick<TValueToFilter>(name, TString::Format("%s{rg{%f,%f}}", name, float(low), float(high))),
    mLow(low),
    mUp(high)
{
}

/// \brief Cut string constructor
/// \param cutstr The cuts string
template <typename TValueToFilter>
CutBrickRange<TValueToFilter>::CutBrickRange(const TString& cutstr)
  : CutBrick<TValueToFilter>(),
    mLow(0),
    mUp(0)
{
  ConstructCutFromString(cutstr);
}

/// \brief Construct the cut from a cut string
/// \param cutstr The cut string
/// The cut string should have the structure
///    name{rg{low,high}}
/// If the cut string is correctly parsed the cut is correctly built
/// if not a fatal exception is raised
template <typename TValueToFilter>
void CutBrickRange<TValueToFilter>::ConstructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());

  std::regex cutregex("^(\\w+)\\{rg\\{((?:-?[\\d]+\\.?[\\d]*)|(?:-?[\\d]*\\.?[\\d]+)),((?:-?[\\d]+\\.?[\\d]*)|(?:-?[\\d]*\\.?[\\d]+))}}$", std::regex_constants::ECMAScript | std::regex_constants::icase);
  std::string in(cutstr.Data());
  std::smatch m;

  std::regex_search(in, m, cutregex);
  if (m.empty() or (m.size() < 4)) {
    LOGF(fatal, "CutBrickRange<TValueToFilter>::ConstructCutFromString", "Wrong RE: %s, use pT{rg{0.2,2.0}} for instance", cutstr.Data());
  } else {
    this->SetName(m[1].str().c_str());
    this->SetTitle(cutstr.Data());
    mLow = TValueToFilter(std::stod(m[2]));
    mUp = TValueToFilter(std::stod(m[3]));
  }
}

/// \brief Returns wether the cut brick is incorporated in the selection chain
/// \return true if the cut brick is incorporated
template <typename TValueToFilter>
std::vector<bool> CutBrickRange<TValueToFilter>::IsArmed()
{
  std::vector<bool> res;
  if (this->mMode == this->kSELECTED) {
    res.push_back(true);
  } else {
    res.push_back(false);
  }
  return res;
}

templateClassImp(CutBrickRange);
template class o2::analysis::PWGCF::CutBrickRange<int>;
template class o2::analysis::PWGCF::CutBrickRange<float>;

///////////////////////////////////////////////////////////////////////////////////////
/// Default constructor
template <typename TValueToFilter>
CutBrickFnRange<TValueToFilter>::CutBrickFnRange()
  : CutBrickRange<TValueToFilter>(),
    mLowFunction{},
    mUpFunction{}
{
}

/// Named constructor
/// \param name The name of the brick
/// \param lowfn The function which will provide the low value for the cut range
/// \param upfn The function which will provide the upper value for the cut range
template <typename TValueToFilter>
CutBrickFnRange<TValueToFilter>::CutBrickFnRange(const char* name, const TF1& lowfn, const TF1& upfn)
  : CutBrickRange<TValueToFilter>(),
    mLowFunction{lowfn},
    mUpFunction{upfn}
{
  this->SetName(name);
}

/// \brief Cut string constructor
/// \param cutstr The cuts string
template <typename TValueToFilter>
CutBrickFnRange<TValueToFilter>::CutBrickFnRange(const TString& cutstr)
  : CutBrickRange<TValueToFilter>(),
    mLowFunction{},
    mUpFunction{}
{
  ConstructCutFromString(cutstr);
}

/// \brief Construct the cut from a cut string
/// \param cutstr The cut string
/// The cut string should have the structure
///    name{fnrg{lowfn,highfn}}
/// If the cut string is correctly parsed the cut is correctly built
/// if not a fatal exception is raised
template <typename TValueToFilter>
void CutBrickFnRange<TValueToFilter>::ConstructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());

  std::regex cutregex("^(\\w+)\\{fnrg\\{(\\w+)=([\\w\\-\\*\\/\\+\\()\\.]+),([\\w\\-\\*\\/\\+\\()\\.]+)}}$", std::regex_constants::ECMAScript | std::regex_constants::icase);
  std::string in(cutstr.Data());
  std::smatch m;

  std::regex_search(in, m, cutregex);
  if (m.empty() or (m.size() < 5)) {
    LOGF(fatal, "CutBrickFnRange<TValueToFilter>::ConstructCutFromString", "Wrong RE: %s, use pT{fnrg{myfn=2.0*sin(x)/x}} for instance", cutstr.Data());
  } else {
    this->SetName(m[2].str().c_str());
    this->SetTitle(cutstr.Data());
    mLowFunction = TF1(TString::Format("%s_low", m[1].str().c_str()), m[3].str().c_str(), 0, 1, TF1::EAddToList::kNo);
    mUpFunction = TF1(TString::Format("%s_up", m[1].str().c_str()), m[4].str().c_str(), 0, 1, TF1::EAddToList::kNo);
    if (not mLowFunction.IsValid() or not mUpFunction.IsValid()) {
      LOGF(fatal, "CutBrickFnRange<TValueToFilter>::ConstructCutFromString", "Wrong function expression: %s, use pT{fnrg{myfn=2.0*sin(x)/x}} for instance", cutstr.Data());
    }
  }
}

templateClassImp(CutBrickFnRange);
template class o2::analysis::PWGCF::CutBrickFnRange<int>;
template class o2::analysis::PWGCF::CutBrickFnRange<float>;

///////////////////////////////////////////////////////////////////////////////////////
/// Default constructor
template <typename TValueToFilter>
CutBrickExtToRange<TValueToFilter>::CutBrickExtToRange()
  : CutBrick<TValueToFilter>(),
    mLow(0),
    mUp(0)
{
}

/// Named constructor
/// \param name The name of the brick
/// \param low The low value for the cut excluded range
/// \param high The high value for the cut excluded range
template <typename TValueToFilter>
CutBrickExtToRange<TValueToFilter>::CutBrickExtToRange(const char* name, const TValueToFilter& low, const TValueToFilter& high)
  : CutBrick<TValueToFilter>(name, TString::Format("%s{xrg{%f,%f}}", name, float(low), float(high))),
    mLow(low),
    mUp(high)
{
}

/// \brief Cut string constructor
/// \param cutstr The cuts string
template <typename TValueToFilter>
CutBrickExtToRange<TValueToFilter>::CutBrickExtToRange(const TString& cutstr)
  : CutBrick<TValueToFilter>(),
    mLow(0),
    mUp(0)
{
  ConstructCutFromString(cutstr);
}

/// \brief Construct the cut from a cut string
/// \param cutstr The cut string
/// The cut string should have the structure
///    name{xrg{low,high}}
/// If the cut string is correctly parsed the cut is correctly built
/// if not a fatal exception is raised
template <typename TValueToFilter>
void CutBrickExtToRange<TValueToFilter>::ConstructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());

  std::regex cutregex("^(\\w+)\\{xrg\\{((?:-?[\\d]+\\.?[\\d]*)|(?:-?[\\d]*\\.?[\\d]+)),((?:-?[\\d]+\\.?[\\d]*)|(?:-?[\\d]*\\.?[\\d]+))}}$", std::regex_constants::ECMAScript | std::regex_constants::icase);
  std::string in(cutstr.Data());
  std::smatch m;

  std::regex_search(in, m, cutregex);
  if (m.empty() or (m.size() < 4)) {
    LOGF(fatal, "CutBrickExtToRange<TValueToFilter>::ConstructCutFromString", "Wrong RE: %s, use minv{xrg{0.02,0.04}} for instance", cutstr.Data());
  } else {
    this->SetName(m[1].str().c_str());
    this->SetTitle(cutstr.Data());
    mLow = TValueToFilter(std::stod(m[2]));
    mUp = TValueToFilter(std::stod(m[3]));
  }
}

/// \brief Returns wether the cut brick is incorporated in the selection chain
/// \return true if the cut brick is incorporated
template <typename TValueToFilter>
std::vector<bool> CutBrickExtToRange<TValueToFilter>::IsArmed()
{
  std::vector<bool> res;
  if (this->mMode == this->kSELECTED) {
    res.push_back(true);
  } else {
    res.push_back(false);
  }
  return res;
}

templateClassImp(CutBrickExtToRange);
template class o2::analysis::PWGCF::CutBrickExtToRange<int>;
template class o2::analysis::PWGCF::CutBrickExtToRange<float>;

///////////////////////////////////////////////////////////////////////////////////////
/// Default constructor
template <typename TValueToFilter>
CutBrickFnExtToRange<TValueToFilter>::CutBrickFnExtToRange()
  : CutBrickExtToRange<TValueToFilter>(),
    mLowFunction{},
    mUpFunction{}
{
}

/// Named constructor
/// \param name The name of the brick
/// \param lowfn The function which will provide the low value for the cut excluded range
/// \param upfn The function which will provide the upper value for the cut excluded range
template <typename TValueToFilter>
CutBrickFnExtToRange<TValueToFilter>::CutBrickFnExtToRange(const char*, const TF1& lowfn, const TF1& upfn)
  : CutBrickExtToRange<TValueToFilter>(),
    mLowFunction{lowfn},
    mUpFunction{upfn}
{
}

/// \brief Cut string constructor
/// \param cutstr The cuts string
template <typename TValueToFilter>
CutBrickFnExtToRange<TValueToFilter>::CutBrickFnExtToRange(const TString& cutstr)
  : CutBrickExtToRange<TValueToFilter>(),
    mLowFunction{},
    mUpFunction{}
{
  ConstructCutFromString(cutstr);
}

/// \brief Construct the cut from a cut string
/// \param cutstr The cut string
/// The cut string should have the structure
///    name{fnxrg{low,high}}
/// If the cut string is correctly parsed the cut is correctly built
/// if not a fatal exception is raised
template <typename TValueToFilter>
void CutBrickFnExtToRange<TValueToFilter>::ConstructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());

  std::regex cutregex("^(\\w+)\\{fnxrg\\{(\\w+)=([\\w\\-\\*\\/\\+\\()\\.]+),([\\w\\-\\*\\/\\+\\()\\.]+)}}$", std::regex_constants::ECMAScript | std::regex_constants::icase);
  std::string in(cutstr.Data());
  std::smatch m;

  std::regex_search(in, m, cutregex);
  if (m.empty() or (m.size() < 5)) {
    LOGF(fatal, "CutBrickFnExtToRange<TValueToFilter>::ConstructCutFromString", "Wrong RE: %s, use pT{fnxrg{myfn=2.0*sin(x)/x}} for instance", cutstr.Data());
  } else {
    this->SetName(m[2].str().c_str());
    this->SetTitle(cutstr.Data());
    mLowFunction = TF1(TString::Format("%s_low", m[1].str().c_str()), m[3].str().c_str(), 0, 1, TF1::EAddToList::kNo);
    mUpFunction = TF1(TString::Format("%s_up", m[1].str().c_str()), m[4].str().c_str(), 0, 1, TF1::EAddToList::kNo);
    if (not mLowFunction.IsValid() or not mUpFunction.IsValid()) {
      LOGF(fatal, "CutBrickFnExtToRange<TValueToFilter>::ConstructCutFromString", "Wrong function expression: %s, use pT{fnxrg{myfn=2.0*sin(x)/x}} for instance", cutstr.Data());
    }
  }
}

templateClassImp(CutBrickFnExtToRange);
template class o2::analysis::PWGCF::CutBrickFnExtToRange<int>;
template class o2::analysis::PWGCF::CutBrickFnExtToRange<float>;

///////////////////////////////////////////////////////////////////////////////////////
/// Default constructor
template <typename TValueToFilter>
CutBrickSelectorMultipleRanges<TValueToFilter>::CutBrickSelectorMultipleRanges()
  : CutBrick<TValueToFilter>()
{
}

/// Named constructor
/// \param name The name of the brick
/// \param edges Vector with the ranges edges
template <typename TValueToFilter>
CutBrickSelectorMultipleRanges<TValueToFilter>::CutBrickSelectorMultipleRanges(const char* name, const std::vector<TValueToFilter>& edges)
  : CutBrick<TValueToFilter>(name, name)
{
  TString title = name;
  title += "{";
  bool first = true;
  for (auto edge : edges) {
    mEdges.push_back(edge);
    if (first) {
      title += TString::Format("%.2f", float(edge));
    } else {
      title += TString::Format(",%.2f", float(edge));
    }
  }
  title += "}";
  this->SetTitle(title.Data());
  for (unsigned int i = 1; i < mEdges.size(); ++i) {
    mActive.push_back(false);
  }
}

/// \brief Cut string constructor
/// \param cutstr The cuts string
template <typename TValueToFilter>
CutBrickSelectorMultipleRanges<TValueToFilter>::CutBrickSelectorMultipleRanges(const TString& cutstr)
  : CutBrick<TValueToFilter>()
{
  ConstructCutFromString(cutstr);
}

/// \brief Construct the cut from a cut string
/// \param cutstr The cut string
/// The cut string should have the structure
///    name{mrg{edge,edge, ...,edge}}
/// If the cut string is correctly parsed the cut is correctly built
/// if not a fatal exception is raised
template <typename TValueToFilter>
void CutBrickSelectorMultipleRanges<TValueToFilter>::ConstructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());

  std::regex cutregex("^(\\w+)\\{mrg\\{([a-zA-Z]\\w+)((?:,((?:-?[\\d]+\\.?[\\d]*)|(?:-?[\\d]*\\.?[\\d]+))){2,})}}$", std::regex_constants::ECMAScript | std::regex_constants::icase);
  std::string in(cutstr.Data());
  std::smatch m;

  bool res = std::regex_search(in, m, cutregex);
  if (not res or m.empty() or (m.size() < 5)) {
    LOGF(fatal, "CutBrickSelectorMultipleRanges<TValueToFilter>::ConstructCutFromString", "Wrong RE: %s, use centmult{mrg{V0M,0,5,10,20,30,40,50,60,70,80}} for instance", cutstr.Data());
  } else {
    this->SetName(m[2].str().c_str());
    this->SetTitle(cutstr.Data());
    /* list of edges */
    TObjArray* tokens = TString(m[3]).Tokenize(",");
    for (int i = 0; i < tokens->GetEntries(); ++i) {
      mEdges.push_back(TValueToFilter(TString(tokens->At(i)->GetName()).Atof()));
    }
    for (unsigned int i = 1; i < mEdges.size(); ++i) {
      mActive.push_back(false);
    }
    delete tokens;
  }
}

/// \brief Returns wether the cut brick is incorporated in the selection chain
/// \return true if the cut brick is incorporated
template <typename TValueToFilter>
std::vector<bool> CutBrickSelectorMultipleRanges<TValueToFilter>::IsArmed()
{
  std::vector<bool> res;
  res.reserve(Length());
  if (this->mMode == this->kSELECTED) {
    for (unsigned int i = 0; i < mActive.size(); ++i) {
      res.push_back(true);
    }
  } else {
    for (unsigned int i = 0; i < mActive.size(); ++i) {
      res.push_back(false);
    }
  }
  return res;
}

/// \brief Filter the passed value to update the brick status accordingly
/// \param value The value to filter
/// \return true if the value passed the cut false otherwise
template <typename TValueToFilter>
std::vector<bool> CutBrickSelectorMultipleRanges<TValueToFilter>::Filter(const TValueToFilter& value)
{
  if ((mEdges.front() <= value) and (value < mEdges.back())) {
    this->mState = this->kACTIVE;
    unsigned int last = mActive.size();
    for (unsigned int i = 0; i < mActive.size(); ++i) {
      if (value < mEdges[i + 1]) {
        mActive[i] = true;
        last = i;
        break;
      } else {
        mActive[i] = false;
      }
    }
    for (unsigned int i = last + 1; i < mActive.size(); ++i) {
      mActive[i] = false;
    }
  } else {
    this->mState = this->kPASSIVE;
    for (unsigned int i = 0; i < mActive.size(); ++i) {
      mActive[i] = false;
    }
  }
  return std::vector<bool>(mActive);
}

templateClassImp(CutBrickSelectorMultipleRanges);
template class o2::analysis::PWGCF::CutBrickSelectorMultipleRanges<int>;
template class o2::analysis::PWGCF::CutBrickSelectorMultipleRanges<float>;

///////////////////////////////////////////////////////////////////////////////////////
/// Default constructor
template <typename TValueToFilter>
CutWithVariations<TValueToFilter>::CutWithVariations()
  : CutBrick<TValueToFilter>(),
    mAllowSeveralDefaults(false),
    mDefaultBricks{},
    mVariationBricks{}
{
}

/// Named constructor
/// \param name The name of the brick
/// \param cutstr The string associated with the cut
/// \param severaldefaults The cut should allow multiple defaults values or not
template <typename TValueToFilter>
CutWithVariations<TValueToFilter>::CutWithVariations(const char* name, const char* cutstr, bool severaldefaults)
  : CutBrick<TValueToFilter>(name, cutstr),
    mAllowSeveralDefaults(severaldefaults),
    mDefaultBricks{},
    mVariationBricks{}
{
}

/// \brief Cut string constructor
/// \param cutstr The cuts string
template <typename TValueToFilter>
CutWithVariations<TValueToFilter>::CutWithVariations(const TString& cutstr)
  : CutBrick<TValueToFilter>(),
    mAllowSeveralDefaults(false),
    mDefaultBricks{},
    mVariationBricks{}
{
  ConstructCutFromString(cutstr);
}

/// \brief Construct the cut from a cut string
/// \param cutstr The cuts string
/// The cut string should have the structure
///    name{def,def,..,def[;alt,alt,...,alt]}
/// where each of the def and alt are basic cut bricks
/// If the cut string is correctly parsed the cut is correctly built
/// if not a fatal exception is raised
template <typename TValueToFilter>
void CutWithVariations<TValueToFilter>::ConstructCutFromString(const TString& cutstr)
{
  /* let's catch the first level */
  LOGF(info, "Cut with variations: cut string: %s", cutstr.Data());

  std::regex cutregex("^(\\w+)\\{cwv\\{([\\w\\d.,:{}=\\-\\+\\*\\/]+)}}$", std::regex_constants::ECMAScript | std::regex_constants::icase);
  std::string in(cutstr.Data());
  std::smatch m;

  bool res = std::regex_search(in, m, cutregex);
  if (not res or m.empty() or (m.size() < 3)) {
    LOGF(fatal, "CutWithVariations<TValueToFilter>::ConstructCutFromString", "Wrong RE: %s, use pT{cwv{rg{0.2,10.0}}} for instance", cutstr.Data());
  }
  this->SetName(m[1].str().c_str());
  this->SetTitle(cutstr.Data());

  /* let's split default and variations */
  TObjArray* lev1toks = TString(m[2]).Tokenize(":");
  if (lev1toks->GetEntries() > 2) {
    LOGF(fatal, "CutWithVariations<TValueToFilter>::ConstructCutFromString", "Wrong RE: %s, use pT{cwv{rg{0.2,10.0}}} for instance", cutstr.Data());
  }
  bool atleastonearmed = false;
  auto addCuts = [&](TList& cutlist, std::string cuttxt, bool reqflag) {
    std::smatch m;
    while (cuttxt.length() > 0) {
      std::set<std::string> allowed = {"lim", "th", "rg", "xrg", "mrg", "fnlim", "fnth", "fnrg", "fnxrg"};
      std::regex cutregex("(\\w+\\{[\\w.,=\\-\\+\\*\\/]+}(?:-(?:no|yes))*)");
      bool res = regex_search(cuttxt, m, cutregex);
      if (not res or m.empty() or m.size() != 2) {
        LOGF(fatal, "CutWithVariations<TValueToFilter>::ConstructCutFromString", "Cut with variations malformed RE %s", cuttxt.c_str());
      }
      TString brickre = m[1].str().c_str();
      bool isarmed = brickre.EndsWith("-yes");
      if (isarmed) {
        brickre.Remove(brickre.Index("-yes"), strlen("-yes"));
      } else if (brickre.EndsWith("-no")) {
        brickre.Remove(brickre.Index("-no"), strlen("-no"));
      } else if (reqflag) {
        LOGF(fatal, "CutWithVariations<TValueToFilter>::ConstructCutFromString", "Wrong RE: %s, alternatives not correctly flagged", cuttxt.c_str());
      }
      CutBrick<TValueToFilter>* brick = CutBrick<TValueToFilter>::constructBrick(this->GetName(), brickre, allowed);
      brick->Arm(isarmed);
      atleastonearmed = atleastonearmed || isarmed;
      cutlist.Add(brick);
      cuttxt = m.suffix();
    }
  };
  /* let's handle the default values */
  {
    LOGF(info, "Cut with variations, extracting defaults: cut string: %s", lev1toks->At(0)->GetName());
    addCuts(mDefaultBricks, lev1toks->At(0)->GetName(), true);
    if (mDefaultBricks.GetEntries() > 1) {
      /* TODO: several default options for track type and for track pid selection */
      LOGF(fatal, "CutWithVariations<TValueToFilter>::ConstructCutFromString", "Wrong RE: %s, several defaults only for trktype or trkpid pending of implementation", cutstr.Data());
    }
  }
  /* let's now handle the variations if any */
  if (lev1toks->GetEntries() > 1) {
    LOGF(info, "Cut with variations, extracting variations: cut string: %s", lev1toks->At(1)->GetName());
    addCuts(mVariationBricks, lev1toks->At(1)->GetName(), true);
  }
  this->Arm(atleastonearmed);
  delete lev1toks;
}

/// \brief Stores the brick with a default value for the cut
/// \param brick pointer to the brick to incorporate
/// \returns true if the brick was successfully added
/// If several default values are allowed it is only required
/// that the name of the new default brick were unique
/// If only one default value is allowed it is required that
/// no previous default were stored
/// If any of the above conditions fails the brick is not
/// added and false is returned
template <typename TValueToFilter>
bool CutWithVariations<TValueToFilter>::AddDefaultBrick(CutBrick<TValueToFilter>* brick)
{
  if (mAllowSeveralDefaults) {
    if (mDefaultBricks.FindObject(brick->GetName())) {
      return false;
    } else {
      mDefaultBricks.Add(brick);
      return true;
    }
  } else {
    if (mDefaultBricks.GetEntries() > 0) {
      return false;
    } else {
      mDefaultBricks.Add(brick);
      return true;
    }
  }
}

/// \brief Stores the brick with the variation to the default value for the cut
/// \param brick pointer to the brick to incorporate
/// \returns true if the brick was successfully added
/// It is required that the brick name were unique in the list of variation values brick
template <typename TValueToFilter>
bool CutWithVariations<TValueToFilter>::AddVariationBrick(CutBrick<TValueToFilter>* brick)
{
  if (mVariationBricks.FindObject(brick->GetName())) {
    return false;
  } else {
    mVariationBricks.Add(brick);
    return true;
  }
}

/// \brief Returns wether the cut brick is incorporated in the selection chain
/// \return true if the cut brick is incorporated
template <typename TValueToFilter>
std::vector<bool> CutWithVariations<TValueToFilter>::IsArmed()
{
  std::vector<bool> res;
  res.reserve(Length());
  if (this->mMode == this->kSELECTED) {
    int nArmedDefaults = 0;
    int nArmedVariants = 0;
    for (int i = 0; i < mDefaultBricks.GetEntries(); ++i) {
      std::vector<bool> tmp = ((CutBrick<TValueToFilter>*)mDefaultBricks.At(i))->IsArmed();
      res.insert(res.end(), tmp.begin(), tmp.end());
      for (bool iarmed : tmp) {
        if (iarmed) {
          nArmedDefaults++;
        }
      }
    }
    for (int i = 0; i < mVariationBricks.GetEntries(); ++i) {
      std::vector<bool> tmp = ((CutBrick<TValueToFilter>*)mVariationBricks.At(i))->IsArmed();
      res.insert(res.end(), tmp.begin(), tmp.end());
      for (bool iarmed : tmp) {
        if (iarmed) {
          nArmedVariants++;
        }
      }
    }
    if (nArmedDefaults > 1 or nArmedVariants > 1 or (nArmedDefaults + nArmedVariants) > 1) {
      LOGF(fatal, "CutWithVariations<TValueToFilter>::IsArmed(%s), More than one alternative selected. Default armed %d, variants armed %d", this->GetName(), nArmedDefaults, nArmedVariants);
    }
  } else {
    for (int i = 0; i < Length(); ++i) {
      res.push_back(false);
    }
  }
  return res;
}

/// Filters the passed value
/// The bricks on the default values list and in the variation
/// values list will change to active or passive accordingly to the passed value
/// \param value The value to filter
/// \returns true if the value activated any of the bricks

template <typename TValueToFilter>
std::vector<bool> CutWithVariations<TValueToFilter>::Filter(const TValueToFilter& value)
{
  std::vector<bool> res;
  res.reserve(Length());
  for (int i = 0; i < mDefaultBricks.GetEntries(); ++i) {
    std::vector<bool> tmp = ((CutBrick<TValueToFilter>*)mDefaultBricks.At(i))->Filter(value);
    res.insert(res.end(), tmp.begin(), tmp.end());
  }
  for (int i = 0; i < mVariationBricks.GetEntries(); ++i) {
    std::vector<bool> tmp = ((CutBrick<TValueToFilter>*)mVariationBricks.At(i))->Filter(value);
    res.insert(res.end(), tmp.begin(), tmp.end());
  }
  return res;
}

/// Return the length needed to code the cut
/// The length is in brick units. The actual length is implementation dependent
/// \returns Cut length in units of bricks
template <typename TValueToFilter>
int CutWithVariations<TValueToFilter>::Length()
{
  /* TODO: should a single default cut without variations return zero length? */
  int length = 0;
  for (int i = 0; i < mDefaultBricks.GetEntries(); ++i) {
    length += ((CutBrick<TValueToFilter>*)mDefaultBricks.At(i))->Length();
  }
  for (int i = 0; i < mVariationBricks.GetEntries(); ++i) {
    length += ((CutBrick<TValueToFilter>*)mVariationBricks.At(i))->Length();
  }
  return length;
}

/// Virtual function. Return the index of the armed brick within this brick
/// \returns The index of the armed brick within this brick. Default -1
template <typename TValueToFilter>
int CutWithVariations<TValueToFilter>::getArmedIndex()
{
  int index = -1;
  if (this->mMode == this->kSELECTED) {
    auto checkBrickList = [&index](auto& brklst) {
      bool found = false;
      for (int i = 0; not found and i < brklst.GetEntries(); ++i) {
        index++;
        std::vector<bool> tmp = ((CutBrick<TValueToFilter>*)brklst.At(i))->IsArmed();
        for (bool iarmed : tmp) {
          if (iarmed) {
            found = true;
            break;
          }
        }
      }
      return found;
    };
    if (not checkBrickList(mDefaultBricks)) {
      if (not checkBrickList(mVariationBricks)) {
        LOGF(fatal, "CutWithVariations::getArmedIndex(). There should be at least one brick armed");
        return -1;
      }
    }
  }
  return index;
}

templateClassImp(o2::analysis::CutWithVariations);
template class o2::analysis::PWGCF::CutWithVariations<float>;
template class o2::analysis::PWGCF::CutWithVariations<int>;

/// Default constructor
SpecialCutBrick::SpecialCutBrick()
  : TNamed(),
    mState(kPASSIVE),
    mMode(kUNSELECTED)
{
}

/// Named constructor
/// \param name The name of the brick
SpecialCutBrick::SpecialCutBrick(const char* name, const char* title)
  : TNamed(name, title),
    mState(kPASSIVE),
    mMode(kUNSELECTED)
{
}

ClassImp(SpecialCutBrick);

/// \brief Constructor from regular expression
TrackSelectionBrick::TrackSelectionBrick(const TString& regex) : SpecialCutBrick()
{
  bool armed = false;
  TString name = regex;
  if (regex.EndsWith("-yes")) {
    name.Remove(regex.Index("-", strlen("-yes")));
    armed = true;
  } else if (regex.EndsWith("-no")) {
    name.Remove(regex.Index("-", strlen("-no")));
    armed = false;
  }
  SetName(name);
  SetTitle(regex);
  if (name.EqualTo("FB1LHC2010")) {
    constructFB1LHC2010();
  } else if (name.EqualTo("FB1")) {
    constructFB1LHC2011();
  } else if (name.EqualTo("FB32LHC2010")) {
    constructFB32LHC2010();
  } else if (name.EqualTo("FB32")) {
    constructFB32LHC2011();
  } else if (name.EqualTo("FB64LHC2010")) {
    constructFB64LHC2010();
  } else if (name.EqualTo("FB64")) {
    constructFB64LHC2011();
  } else {
    LOGF(fatal, "TrackSelectionBrick::TrackSelectionBrick", "Wrong RE: %s, trying to construct an unknown track type selector", regex.Data());
  }
  this->Arm(armed);
}

// Default TPC only track selection according to LHC2010
void TrackSelectionBrick::constructFB1LHC2010()
{
  SetTrackType(o2::aod::track::Run2Track);
  SetRequireGoldenChi2(true);
  SetMinNClustersTPC(50);
  SetMaxChi2PerClusterTPC(4.f);
  SetMaxDcaXY(2.4f);
  SetMaxDcaZ(3.2f);
  /* TODO: 2D DCA cut */
}

// Default track selection requiring one hit in the SPD DCAxy according to LHC2010
void TrackSelectionBrick::constructFB32LHC2010()
{
  SetTrackType(o2::aod::track::Run2Track);
  SetRequireITSRefit(true);
  SetRequireTPCRefit(true);
  SetRequireGoldenChi2(true);
  SetMinNCrossedRowsTPC(70);
  SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
  SetMaxChi2PerClusterTPC(4.f);
  SetRequireHitsInITSLayers(1, {0, 1}); // one hit in any SPD layer
  SetMaxChi2PerClusterITS(36.f);
  SetMaxDcaXYPtDep([](float pt) { return 0.0182f + 0.0350f / pow(pt, 1.01f); });
  SetMaxDcaZ(2.f);
}

// Default track selection requiring no hit in the SPD and one in the innermost DCAxy according to LHC2010
// SDD -> complementary tracks to global selection
void TrackSelectionBrick::constructFB64LHC2010()
{
  constructFB32LHC2010();
  ResetITSRequirements();
  SetRequireNoHitsInITSLayers({0, 1}); // no hit in SPD layers
  SetRequireHitsInITSLayers(1, {2});   // one hit in first SDD layer
}

// Default TPC only track selection according to LHC2011
void TrackSelectionBrick::constructFB1LHC2011()
{
  /* the same as for LHC2010 */
  constructFB1LHC2010();
}

// Default track selection requiring one hit in the SPD DCAxy according to LHC2011
void TrackSelectionBrick::constructFB32LHC2011()
{
  SetTrackType(o2::aod::track::Run2Track);
  SetRequireITSRefit(true);
  SetRequireTPCRefit(true);
  SetRequireGoldenChi2(true);
  SetMinNCrossedRowsTPC(70);
  SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
  SetMaxChi2PerClusterTPC(4.f);
  SetRequireHitsInITSLayers(1, {0, 1}); // one hit in any SPD layer
  SetMaxChi2PerClusterITS(36.f);
  SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / pow(pt, 1.1f); });
  SetMaxDcaZ(2.f);
}

// Default track selection requiring no hit in the SPD and one in the innermost DCAxy according to LHC2011
// SDD -> complementary tracks to global selection
void TrackSelectionBrick::constructFB64LHC2011()
{
  constructFB32LHC2011();
  ResetITSRequirements();
  SetRequireNoHitsInITSLayers({0, 1}); // no hit in SPD layers
  SetRequireHitsInITSLayers(1, {2});   // one hit in first SDD layer
}

bool TrackSelectionBrick::FulfillsITSHitRequirements(uint8_t itsClusterMap)
{
  constexpr uint8_t bit = 1;
  for (auto& itsRequirement : mRequiredITSHits) {
    auto hits = std::count_if(itsRequirement.second.begin(), itsRequirement.second.end(), [&](auto&& requiredLayer) { return itsClusterMap & (bit << requiredLayer); });
    if ((itsRequirement.first == -1) && (hits > 0)) {
      return false; // no hits were required in specified layers
    } else if (hits < itsRequirement.first) {
      return false; // not enough hits found in specified layers
    }
  }
  return true;
}

const std::string TrackSelectionBrick::mCutNames[static_cast<int>(TrackSelectionBrick::TrackCuts::kNCuts)] =
  {
    "TrackType",
    "TPCNCls",
    "TPCCrossedRowsOverNCls",
    "TPCRefit",
    "ITSNCls",
    "ITSChi2NDF",
    "ITSRefit",
    "ITSHits",
    "GoldenChi2",
    "DCAxy",
    "DCAz"};

/// \brief Returns wether the cut brick is incorporated in the selection chain
/// \return true if the cut brick is incorporated
std::vector<bool> TrackSelectionBrick::IsArmed()
{
  std::vector<bool> res;
  if (this->mMode == this->kSELECTED) {
    res.push_back(true);
  } else {
    res.push_back(false);
  }
  return res;
}

ClassImp(TrackSelectionBrick);
