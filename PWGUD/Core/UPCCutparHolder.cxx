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

#include "UPCCutparHolder.h"

// setters
void UPCCutparHolder::setUseFwdCuts(bool useFwdCuts) { fUseFwdCuts = useFwdCuts; }
void UPCCutparHolder::setTrackType(int trackType) { fTrackType = trackType; }
void UPCCutparHolder::setFwdPtLow(float fwdPtLow) { fFwdPtLow = fwdPtLow; }
void UPCCutparHolder::setFwdPtHigh(float fwdPtHigh) { fFwdPtHigh = fwdPtHigh; }
void UPCCutparHolder::setFwdEtaLow(float fwdEtaLow) { fFwdEtaLow = fwdEtaLow; }
void UPCCutparHolder::setFwdEtaHigh(float fwdEtaHigh) { fFwdEtaHigh = fwdEtaHigh; }
void UPCCutparHolder::setMuonRAtAbsorberEndLow(float muonRAtAbsorberEndLow) { fMuonRAtAbsorberEndLow = muonRAtAbsorberEndLow; }
void UPCCutparHolder::setMuonRAtAbsorberEndHigh(float muonRAtAbsorberEndHigh) { fMuonRAtAbsorberEndHigh = muonRAtAbsorberEndHigh; }
void UPCCutparHolder::setMuonPDcaHighFirst(float muonPDcaHighFirst) { fMuonPDcaHighFirst = muonPDcaHighFirst; }
void UPCCutparHolder::setMuonPDcaHighSecond(float muonPDcaHighSecond) { fMuonPDcaHighSecond = muonPDcaHighSecond; }
void UPCCutparHolder::setFwdChi2Low(float fwdChi2Low) { fFwdChi2Low = fwdChi2Low; }
void UPCCutparHolder::setFwdChi2High(float fwdChi2High) { fFwdChi2High = fwdChi2High; }
void UPCCutparHolder::setUseBarCuts(bool useBarCuts) { fUseBarCuts = useBarCuts; }
void UPCCutparHolder::setBarPtLow(float barPtLow) { fBarPtLow = barPtLow; }
void UPCCutparHolder::setBarPtHigh(float barPtHigh) { fBarPtHigh = barPtHigh; }
void UPCCutparHolder::setBarEtaLow(float barEtaLow) { fBarEtaLow = barEtaLow; }
void UPCCutparHolder::setBarEtaHigh(float barEtaHigh) { fBarEtaHigh = barEtaHigh; }
void UPCCutparHolder::setITSNClusLow(int ITSNClusLow) { fITSNClusLow = ITSNClusLow; }
void UPCCutparHolder::setITSNClusHigh(int ITSNClusHigh) { fITSNClusHigh = ITSNClusHigh; }
void UPCCutparHolder::setITSChi2Low(float ITSChi2Low) { fITSChi2Low = ITSChi2Low; }
void UPCCutparHolder::setITSChi2High(float ITSChi2High) { fITSChi2High = ITSChi2High; }
void UPCCutparHolder::setTPCNClsLow(int TPCNClsLow) { fTPCNClsLow = TPCNClsLow; }
void UPCCutparHolder::setTPCNClsHigh(int TPCNClsHigh) { fTPCNClsHigh = TPCNClsHigh; }
void UPCCutparHolder::setTPCChi2Low(float TPCChi2Low) { fTPCChi2Low = TPCChi2Low; }
void UPCCutparHolder::setTPCChi2High(float TPCChi2High) { fTPCChi2High = TPCChi2High; }
void UPCCutparHolder::setCheckMaxDcaXY(bool checkMaxDcaXY) { fCheckMaxDcaXY = checkMaxDcaXY; }
void UPCCutparHolder::setDcaZLow(float dcaZLow) { fDcaZLow = dcaZLow; }
void UPCCutparHolder::setDcaZHigh(float dcaZHigh) { fDcaZHigh = dcaZHigh; }
void UPCCutparHolder::setRequireTOF(bool requireTOF) { fRequireTOF = requireTOF; }
void UPCCutparHolder::setRequireITSTPC(bool requireITSTPC) { fRequireITSTPC = requireITSTPC; }
void UPCCutparHolder::setProduceITSITS(bool produceITSITS) { fProduceITSITS = produceITSITS; }
void UPCCutparHolder::setMaxNContrib(int maxNContrib) { fMaxNContrib = maxNContrib; }
void UPCCutparHolder::setAmbigSwitch(int ambigSwitch) { fAmbigSwitch = ambigSwitch; }

// getters
bool UPCCutparHolder::getUseFwdCuts() const { return fUseFwdCuts; }
int UPCCutparHolder::getTrackType() const { return fTrackType; }
float UPCCutparHolder::getFwdPtLow() const { return fFwdPtLow; }
float UPCCutparHolder::getFwdPtHigh() const { return fFwdPtHigh; }
float UPCCutparHolder::getFwdEtaLow() const { return fFwdEtaLow; }
float UPCCutparHolder::getFwdEtaHigh() const { return fFwdEtaHigh; }
float UPCCutparHolder::getMuonRAtAbsorberEndLow() const { return fMuonRAtAbsorberEndLow; }
float UPCCutparHolder::getMuonRAtAbsorberEndHigh() const { return fMuonRAtAbsorberEndHigh; }
float UPCCutparHolder::getMuonPDcaHighFirst() const { return fMuonPDcaHighFirst; }
float UPCCutparHolder::getMuonPDcaHighSecond() const { return fMuonPDcaHighSecond; }
float UPCCutparHolder::getFwdChi2Low() const { return fFwdChi2Low; }
float UPCCutparHolder::getFwdChi2High() const { return fFwdChi2High; }
bool UPCCutparHolder::getUseBarCuts() const { return fUseBarCuts; }
float UPCCutparHolder::getBarPtLow() const { return fBarPtLow; }
float UPCCutparHolder::getBarPtHigh() const { return fBarPtHigh; }
float UPCCutparHolder::getBarEtaLow() const { return fBarEtaLow; }
float UPCCutparHolder::getBarEtaHigh() const { return fBarEtaHigh; }
int UPCCutparHolder::getITSNClusLow() const { return fITSNClusLow; }
int UPCCutparHolder::getITSNClusHigh() const { return fITSNClusHigh; }
float UPCCutparHolder::getITSChi2Low() const { return fITSChi2Low; }
float UPCCutparHolder::getITSChi2High() const { return fITSChi2High; }
int UPCCutparHolder::getTPCNClsLow() const { return fTPCNClsLow; }
int UPCCutparHolder::getTPCNClsHigh() const { return fTPCNClsHigh; }
float UPCCutparHolder::getTPCChi2Low() const { return fTPCChi2Low; }
float UPCCutparHolder::getTPCChi2High() const { return fTPCChi2High; }
bool UPCCutparHolder::getCheckMaxDcaXY() const { return fCheckMaxDcaXY; }
float UPCCutparHolder::getDcaZLow() const { return fDcaZLow; }
float UPCCutparHolder::getDcaZHigh() const { return fDcaZHigh; }
bool UPCCutparHolder::getRequireTOF() const { return fRequireTOF; }
bool UPCCutparHolder::getRequireITSTPC() const { return fRequireITSTPC; }
bool UPCCutparHolder::getProduceITSITS() const { return fProduceITSITS; }
int UPCCutparHolder::getMaxNContrib() const { return fMaxNContrib; }
int UPCCutparHolder::getAmbigSwitch() const { return fAmbigSwitch; }
