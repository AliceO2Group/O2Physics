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

#ifndef PWGMM_MULT_CORE_INCLUDE_FUNCTIONS_H_
#define PWGMM_MULT_CORE_INCLUDE_FUNCTIONS_H_

namespace pwgmm::mult
{
template <typename MCC>
concept has_hepmc_xs = requires(MCC::iterator const& mcc) {
  mcc.xsectGen();
};

template <typename MCC>
concept has_hepmc_pid = requires(MCC::iterator const& mcc) {
  mcc.processId();
};

template <typename MCC>
concept has_hepmc_pdf = requires(MCC::iterator const& mcc) {
  mcc.pdfId1();
  mcc.pdfId2();
};

template <typename MCC>
concept has_hepmc_hi = requires(MCC::iterator const& mcc) {
  mcc.ncollHard();
  mcc.ncoll();
};

template <typename C>
concept has_FT0C = requires(C::iterator const& c) {
  c.centFT0C();
};

template <typename IC>
concept iterator_with_FT0C = requires(IC const& c) {
  c.centFT0C();
};

template <typename C>
concept has_centFT0CVariant1 = requires(C::iterator const& c) {
  c.centFT0CVariant1();
};

template <typename IC>
concept iterator_with_centFT0CVariant1 = requires(IC const& c) {
  c.centFT0CVariant1();
};

template <typename C>
concept has_FT0M = requires(C::iterator const& c) {
  c.centFT0M();
};

template <typename IC>
concept iterator_with_FT0M = requires(IC const& c) {
  c.centFT0M();
};

template <typename C>
concept has_centNGlobal = requires(C::iterator const& c) {
  c.centNGlobal();
};

template <typename IC>
concept iterator_with_centNGlobal = requires(IC const& c) {
  c.centNGlobal();
};

template <typename C>
concept has_centMFT = requires(C::iterator const& c) {
  c.centMFT();
};

template <typename IC>
concept iterator_with_centMFT = requires(IC const& c) {
  c.centMFT();
};

template <typename C>
concept has_genFT0C = requires(C::iterator const& c) {
  c.gencentFT0C();
};

template <typename C>
concept has_genFT0M = requires(C::iterator const& c) {
  c.gencentFT0M();
};

template <typename C>
concept iterator_with_genFT0C = requires(C const& c) {
  c.gencentFT0C();
};

template <typename C>
concept iterator_with_genFT0M = requires(C const& c) {
  c.gencentFT0M();
};

template <typename C>
concept has_reco_cent = has_FT0C<C> || has_centFT0CVariant1<C> || has_FT0M<C> || has_centNGlobal<C> || has_centMFT<C>;

template <typename C>
concept has_gen_cent = has_genFT0C<C> && has_genFT0M<C>;

template <typename MCC>
concept has_Centrality = requires(MCC::iterator const& mcc) {
  mcc.centrality();
};

template <typename MCC>
concept iterator_with_Centrality = requires(MCC const& mcc) {
  mcc.centrality();
};

template <typename C>
  requires(!(iterator_with_FT0C<C> || iterator_with_centFT0CVariant1<C> || iterator_with_FT0M<C> || iterator_with_centNGlobal<C> || iterator_with_centMFT<C>))
static float getRecoCent(C const&)
{
  return -1;
}

template <iterator_with_FT0C C>
static float getRecoCent(C const& collision)
{
  return collision.centFT0C();
}

template <iterator_with_centFT0CVariant1 C>
static float getRecoCent(C const& collision)
{
  return collision.centFT0CVariant1();
}

template <iterator_with_FT0M C>
static float getRecoCent(C const& collision)
{
  return collision.centFT0M();
}

template <iterator_with_centNGlobal C>
static float getRecoCent(C const& collision)
{
  return collision.centNGlobal();
}

template <iterator_with_centMFT C>
static float getRecoCent(C const& collision)
{
  return collision.centMFT();
}

template <typename MCC>
  requires(!iterator_with_genFT0C<MCC>)
static float getGenCentFT0C(MCC const&)
{
  return -1;
}

template <typename MCC>
  requires(!iterator_with_genFT0M<MCC>)
static float getGenCentFT0M(MCC const&)
{
  return -1;
}

template <iterator_with_genFT0C MCC>
static float getGenCentFT0C(MCC const& mccollision)
{
  return mccollision.gencentFT0C();
}

template <iterator_with_genFT0M MCC>
static float getGenCentFT0M(MCC const& mccollision)
{
  return mccollision.gencentFT0M();
}

template <typename MCC>
  requires(!iterator_with_Centrality<MCC>)
static float getSimCent(MCC const&)
{
  return -1;
}

template <iterator_with_Centrality MCC>
static float getSimCent(MCC const& mccollision)
{
  return mccollision.centrality();
}
} // namespace pwgmm::mult

#endif // PWGMM_MULT_CORE_INCLUDE_FUNCTIONS_H_
