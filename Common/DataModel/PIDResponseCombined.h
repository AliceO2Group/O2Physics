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

///
/// \file   PIDResponseCombined.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Set of tables, tasks and utilities to provide the interface between
///         the analysis data model and the combined PID response
///

#ifndef COMMON_DATAMODEL_PIDRESPONSECOMBINED_H_
#define COMMON_DATAMODEL_PIDRESPONSECOMBINED_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <ReconstructionDataFormats/PID.h>

#include <cstdint>

namespace o2::aod
{

namespace pidbayes
{
typedef int8_t binned_prob_t;
// Bayesian probabilities with reduced size
DECLARE_SOA_COLUMN(BayesEl, bayesEl, binned_prob_t);                //! Bayesian probability for electron expressed in %
DECLARE_SOA_COLUMN(BayesMu, bayesMu, binned_prob_t);                //! Bayesian probability for muon expressed in %
DECLARE_SOA_COLUMN(BayesPi, bayesPi, binned_prob_t);                //! Bayesian probability for pion expressed in %
DECLARE_SOA_COLUMN(BayesKa, bayesKa, binned_prob_t);                //! Bayesian probability for kaon expressed in %
DECLARE_SOA_COLUMN(BayesPr, bayesPr, binned_prob_t);                //! Bayesian probability for proton expressed in %
DECLARE_SOA_COLUMN(BayesDe, bayesDe, binned_prob_t);                //! Bayesian probability for deuteron expressed in %
DECLARE_SOA_COLUMN(BayesTr, bayesTr, binned_prob_t);                //! Bayesian probability for triton expressed in %
DECLARE_SOA_COLUMN(BayesHe, bayesHe, binned_prob_t);                //! Bayesian probability for helium3 expressed in %
DECLARE_SOA_COLUMN(BayesAl, bayesAl, binned_prob_t);                //! Bayesian probability for alpha expressed in %
DECLARE_SOA_COLUMN(BayesProb, bayesProb, binned_prob_t);            //! Bayesian probability of the most probable ID
DECLARE_SOA_COLUMN(BayesID, bayesID, o2::track::pid_constants::ID); //! Most probable ID

} // namespace pidbayes

// Table for each particle hypothesis
DECLARE_SOA_TABLE(pidBayesEl, "AOD", "pidBayesEl", //! Binned (in percentage) Bayesian probability of having a Electron
                  pidbayes::BayesEl);
DECLARE_SOA_TABLE(pidBayesMu, "AOD", "pidBayesMu", //! Binned (in percentage) Bayesian probability of having a Muon
                  pidbayes::BayesMu);
DECLARE_SOA_TABLE(pidBayesPi, "AOD", "pidBayesPi", //! Binned (in percentage) Bayesian probability of having a Pion
                  pidbayes::BayesPi);
DECLARE_SOA_TABLE(pidBayesKa, "AOD", "pidBayesKa", //! Binned (in percentage) Bayesian probability of having a Kaon
                  pidbayes::BayesKa);
DECLARE_SOA_TABLE(pidBayesPr, "AOD", "pidBayesPr", //! Binned (in percentage) Bayesian probability of having a Proton
                  pidbayes::BayesPr);
DECLARE_SOA_TABLE(pidBayesDe, "AOD", "pidBayesDe", //! Binned (in percentage) Bayesian probability of having a Deuteron
                  pidbayes::BayesDe);
DECLARE_SOA_TABLE(pidBayesTr, "AOD", "pidBayesTr", //! Binned (in percentage) Bayesian probability of having a Triton
                  pidbayes::BayesTr);
DECLARE_SOA_TABLE(pidBayesHe, "AOD", "pidBayesHe", //! Binned (in percentage) Bayesian probability of having a Helium3
                  pidbayes::BayesHe);
DECLARE_SOA_TABLE(pidBayesAl, "AOD", "pidBayesAl", //! Binned (in percentage) Bayesian probability of having a Alpha
                  pidbayes::BayesAl);

// Table for the most probable particle
DECLARE_SOA_TABLE(pidBayes, "AOD", "pidBayes", pidbayes::BayesProb, pidbayes::BayesID); //! Index of the most probable ID and its bayesian probability

} // namespace o2::aod

#endif // COMMON_DATAMODEL_PIDRESPONSECOMBINED_H_
