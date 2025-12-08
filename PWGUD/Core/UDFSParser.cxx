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

#include "UDFSParser.h"

#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/BunchFilling.h"
#include "Framework/Logger.h"

#include <string>
#include <vector>

// -----------------------------------------------------------------------------
UDFSParser::UDFSParser(const char* filename)
{
  fisActive = readFS(filename);
  if (!fisActive) {
    LOGF(error, "<UDFSParser> Error while reading %s", filename);
  }
}

// -----------------------------------------------------------------------------
// read csv file with filling scheme definition
// compute [empty, A-side, C-side, BB] BCs according to ALICE counting
bool UDFSParser::readFS(const char* filename)
{
  // open file
  std::ifstream infile(filename);
  if (!infile) {
    LOGF(info, "<UDFSParser> file '%s' is not properly opened.", filename);
    return false;
  }

  // read file line by line
  int bucket;
  bool inb1 = false;
  bool inb2 = false;
  std::string line;
  while (std::getline(infile, line)) {
    // look for 'HEAD ON COLLISIONS FOR B1'
    if (line.find("FOR B1") != std::string::npos) {
      inb1 = true;
      inb2 = false;
    } else if (line.find("FOR B2") != std::string::npos) {
      inb1 = false;
      inb2 = true;
    }
    if (inb1 || inb2) {
      // 5 comma separated items
      auto toks = tokenize(line);
      if (toks.size() == 5) {
        if (!isNumber(toks[0])) {
          continue;
        }
        bucket = stoi(toks[0]) / 10;
        if (inb1) {
          // fill fP2BCsA and fP2BCsBB
          fP2BCsA.push_back(o2::constants::lhc::LHCBunch2P2BC(bucket, o2::constants::lhc::BeamDirection(0)));
          if (isNumber(toks[2])) {
            fP2BCsBB.push_back(o2::constants::lhc::LHCBunch2P2BC(bucket, o2::constants::lhc::BeamDirection(0)));
          }
        } else {
          // fill fP2BCsC
          fP2BCsC.push_back(o2::constants::lhc::LHCBunch2P2BC(bucket, o2::constants::lhc::BeamDirection(1)));
        }
      }
    }
  }

  // fill fP2BCsE
  for (int bcnum = 0; bcnum < o2::constants::lhc::LHCMaxBunches; bcnum++) {
    if (!isP2BCA(bcnum) && !isP2BCC(bcnum) && !isP2BCBB(bcnum)) {
      fP2BCsE.push_back(bcnum);
    }
  }
  return true;
}

// -----------------------------------------------------------------------------
std::string UDFSParser::patternString(int ibeam)
{

  // A or C beam?
  std::vector<int> P2BCs2use;
  switch (ibeam) {
    case 0:
      P2BCs2use = fP2BCsA;
      break;
    case 1:
      P2BCs2use = fP2BCsC;
      break;
    default:
      LOGF(info, "Value of parameter ibeam must be 0 or 1!");
      return std::string("Not defined");
  }

  // convert P2BCs to buckets
  std::vector<int> buckets;
  for (auto bcnum : P2BCs2use) {
    auto bucket = o2::constants::lhc::P2BC2LHCBunch(bcnum, o2::constants::lhc::BeamDirection(ibeam)) * 10;
    LOGF(debug, "bcnum %d, bucket %d", bcnum, bucket);
    buckets.push_back(bucket);
  }

  return o2::BunchFilling::buckets2PatternString(buckets, o2::constants::lhc::BeamDirection(ibeam));
}

// -----------------------------------------------------------------------------
void UDFSParser::Print()
{
  LOGF(info, "Number of colliding BCs: %d", fP2BCsBB.size());
  LOGF(info, "  P2BC, BucketA, BucketC");
  for (auto bcnum : fP2BCsBB) {
    LOGF(info, "  %d, %d, %d", bcnum,
         o2::constants::lhc::P2BC2LHCBunch(bcnum, o2::constants::lhc::BeamDirection(0)),
         o2::constants::lhc::P2BC2LHCBunch(bcnum, o2::constants::lhc::BeamDirection(1)));
  }
}

// -----------------------------------------------------------------------------
bool UDFSParser::isNumber(std::string s)
{
  return !s.empty() && std::find_if(s.begin(),
                                    s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}

// -----------------------------------------------------------------------------
std::string UDFSParser::trim(std::string str,
                             std::string whitespace = " \t")
{
  const auto strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos)
    return ""; // no content

  const auto strEnd = str.find_last_not_of(whitespace);
  const auto strRange = strEnd - strBegin + 1;

  return str.substr(strBegin, strRange);
}

// -----------------------------------------------------------------------------
std::vector<std::string> UDFSParser::tokenize(std::string& str,
                                              std::string separator)
{
  // vector with tokens
  char* token;
  std::vector<std::string> tokens;

  // tokenize
  char* p = str.data();
  while ((token = strtok_r(p, separator.data(), &p))) {
    tokens.push_back(trim(std::string(token)));
  }
  return tokens;
}

// -----------------------------------------------------------------------------
bool UDFSParser::isInVector(int num, std::vector<int> vec)
{
  return !fisActive || find(vec.begin(), vec.end(), num) != vec.end();
}

// -----------------------------------------------------------------------------
bool UDFSParser::isP2BCE(int bcnum)
{
  return isInVector(bcnum, fP2BCsE);
}

// -----------------------------------------------------------------------------
bool UDFSParser::isP2BCA(int bcnum)
{
  return isInVector(bcnum, fP2BCsA);
}

// -----------------------------------------------------------------------------
bool UDFSParser::isP2BCC(int bcnum)
{
  return isInVector(bcnum, fP2BCsC);
}

// -----------------------------------------------------------------------------
bool UDFSParser::isP2BCBB(int bcnum)
{
  return isInVector(bcnum, fP2BCsBB);
}

// -----------------------------------------------------------------------------
int UDFSParser::P2BCtype(int bcnum)
{
  if (isP2BCBB(bcnum)) {
    return 3;
  } else if (isP2BCA(bcnum)) {
    return 1;
  } else if (isP2BCC(bcnum)) {
    return 2;
  } else if (isP2BCE(bcnum)) {
    return 0;
  } else {
    return -1;
  }
}

// -----------------------------------------------------------------------------
