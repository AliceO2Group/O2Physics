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

#ifndef PWGUD_CORE_UDFSPARSER_H_
#define PWGUD_CORE_UDFSPARSER_H_

// #include <gandiva/projector.h>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// fillingScheme can be downloaded from
//
class UDFSParser
{
 public:
  // constructor/destructor
  UDFSParser() {}
  explicit UDFSParser(const char* filename);
  ~UDFSParser() {}

  // read filleingScheme from file
  bool readFS(const char* filename);

  // check type of BC
  // 0: Empty
  // 1: A-side
  // 2: C-side
  // 3: BB
  //-1: noDef
  int P2BCtype(int bcnum);

  // check P2BC
  bool isP2BCE(int nBC);
  bool isP2BCA(int nBC);
  bool isP2BCC(int nBC);
  bool isP2BCBB(int nBC);

  // get pattern string
  std::string patternString(int ibeam);

  void Print();

 private:
  // has been loaded with filling scheme
  bool fisActive = false;

  // vectors with BCs of different type
  std::vector<int> fP2BCsE;  // empty
  std::vector<int> fP2BCsA;  // A-side
  std::vector<int> fP2BCsC;  // C-side
  std::vector<int> fP2BCsBB; // BB

  // helper functions for string parsing
  bool isNumber(std::string s);
  std::string trim(std::string str, std::string whitespace);
  std::vector<std::string> tokenize(std::string& str, std::string separator = ",");
  bool isInVector(int num, std::vector<int> vec);

  // ClassDefNV(UDFSParser, 1);
};

#endif // PWGUD_CORE_UDFSPARSER_H_
