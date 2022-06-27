"""
\author Luca Micheletti <luca.micheletti@cern.ch>, CERN
"""
import yaml
import sys
import argparse
from array import array
from os import path
import ROOT
from ROOT import TFile, TF1, TH1F, TTree
from ROOT import gRandom
from DQFitter import DQFitter

def GenerateTutorialSample():
    '''
    This method create the sample for the tutorial
    '''
    nEvents = 100000
    SigOverBkg = 0.1
    fOut = TFile("tutorial.root", "RECREATE")

    funcMassBkg = TF1("funcMassBkg", "expo", 0., 5.)
    funcMassBkg.SetParameter(0, 0.00)
    funcMassBkg.SetParameter(1, -0.5)

    funcMassSig = TF1("funcMassSig", "gaus", 0., 5.)
    funcMassSig.SetParameter(0, 1.0)
    funcMassSig.SetParameter(1, 3.1)
    funcMassSig.SetParameter(2, 0.07)

    histMass = TH1F("histMass", "histMass", 100, 0., 5.)
    histMass.FillRandom("funcMassBkg", int(nEvents - (nEvents * SigOverBkg)))
    histMass.FillRandom("funcMassSig", int(nEvents * SigOverBkg))
    histMass.Write()

    m = array('f', [0.])
    tree = TTree("data", "data")
    tree.Branch("m", m, "m/F")

    for iEvent in range(0, nEvents):
        seed = gRandom.Rndm()
        if seed > SigOverBkg:
            m[0] = funcMassBkg.GetRandom()
        else:
            m[0] = funcMassSig.GetRandom()
        tree.Fill()
    tree.Write()

    fOut.Close()

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--gen_tutorial", help="generate tutorial sample", action="store_true")
    parser.add_argument("--run_fit", help="run the multi trial", action="store_true")
    args = parser.parse_args()

    print('Loading task configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading task configuration: Done!')

    if args.gen_tutorial:
        GenerateTutorialSample()

    if args.run_fit:
        dqFitter = DQFitter(inputCfg["input"]["input_file_name"], inputCfg["input"]["input_name"])
        dqFitter.SetFitConfig(inputCfg["input"]["pdf_dictionary"])
        dqFitter.MultiTrial()
main()