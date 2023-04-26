#!/usr/bin/env python3

import argparse
import os
import json
from os.path import join, abspath, exists
import re
import subprocess
import ROOT as root


def clean(dir,name,printtext):
    if printtext:
        print("==> Delete simulation artefacts in the output directory\n")
    for f in os.listdir(dir):
        pattern_del=name+"_|o2simtopology|MCStepLogger"
        #pattern_keep="Kine.root|serverlog|mergerlog|workerlog"
        pattern_keep="Kine.root|serverlog"
        if (re.search(pattern_del, f)) and (not re.search(pattern_keep, f) ):
            os.remove(abspath(join(dir,f)))


def read_sim_config(args):
    IniName = "GeneratorEMCocktail.ini"
    generatorName = "${O2DPG_ROOT}/MC/config/PWGEM/external/generator/GeneratorEMCocktailV2.C"
    generatorFunc = "GenerateEMCocktail"

    # matches GeneratorParamEMlibV2.h
    GeneratorParamEMlibV2_CollisionSystem={"kpp900GeV":0x000, "kpp2760GeV":0x64, "kpp7TeV":0xC8, "kpPb":0x12C, "kPbPb":0x190}
    GeneratorParamEMlibV2_Centrality={"kpp":0x0, "k0005":0x1, "k0510":0x2, "k1020":0x3, "k2030":0x4, "k3040":0x5, "k4050":0x6, "k5060":0x7, "k0010":0x8, "k2040":0x9, "k4060":0xA, "k6080":0xB, "k0020":0xC, "k0040":0xD, "k2080":0xE, "k4080":0xF, "k2050":0x10, "kCentralities":0x11}

    allKeys = ["collisionSystem", "centrality", "decayMode","selectedMothers", "paramFile", "paramFileDir", "numberOfParticles", "minPt", "maxPt", "pythiaErrorTolerance", "externalDecayer", "decayLongLived", "dynamicalPtRange", "useYWeights", "paramV2FileDir", "toFixEP", "yGenRange", "useLMeeDecaytable", "weightingMode"]
    if not exists(args.sim_config_file):
        print("ERROR: File "+args.sim_config_file+" not found")
        return 1

    with open(args.sim_config_file, "r") as f:
        simConfigs = json.load(f)

    print("==> Reading simualtion configuration from "+abspath(args.sim_config_file))

    missing_key = False
    for key in allKeys:
        if not key in simConfigs:
            print("ERROR: Missing key \'"+key+"\'")
            missing_key = True
    if missing_key:
        return 1

    collisionSystem      = GeneratorParamEMlibV2_CollisionSystem[simConfigs["collisionSystem"]]
    centrality            = GeneratorParamEMlibV2_Centrality[simConfigs["centrality"]]
    decayMode             = simConfigs["decayMode"]
    selectedMothers       = simConfigs["selectedMothers"]
    paramFile             = simConfigs["paramFile"]
    paramFileDir          = simConfigs["paramFileDir"]
    numberOfParticles     = simConfigs["numberOfParticles"]
    minPt                 = simConfigs["minPt"]
    maxPt                 = simConfigs["maxPt"]
    pythiaErrorTolerance  = simConfigs["pythiaErrorTolerance"]
    externalDecayer       = simConfigs["externalDecayer"]
    decayLongLived        = simConfigs["decayLongLived"]
    dynamicalPtRange      = simConfigs["dynamicalPtRange"]
    useYWeights           = simConfigs["useYWeights"]
    paramV2FileDir        = simConfigs["paramV2FileDir"]
    toFixEP               = simConfigs["toFixEP"]
    yGenRange             = simConfigs["yGenRange"]
    useLMeeDecaytable     = simConfigs["useLMeeDecaytable"]
    weightingMode           = simConfigs["weightingMode"]

    print(f"collisionSystem\t\t{collisionSystem}\ncentrality\t\t{centrality}\ndecayMode\t\t{decayMode}\nselectedMothers\t\t{selectedMothers}\nparamFile\t\t{paramFile}\nparamFileDir\t\t{paramFileDir}\nnumberOfParticles\t{numberOfParticles}\nminPt\t\t\t{minPt}\nmaxPt\t\t\t{maxPt}\npythiaErrorTolerance\t{pythiaErrorTolerance}")
    print(f"externalDecayer\t\t{externalDecayer}\ndecayLongLived\t\t{decayLongLived}\ndynamicalPtRange\t{dynamicalPtRange}\nuseYWeights\t\t{useYWeights}\nparamV2FileDir\t{paramV2FileDir}\ntoFixEP\t\t\t{toFixEP}\nyGenRange\t\t{yGenRange}")
    print(f"useLMeeDecaytable\t{useLMeeDecaytable}\nweightingMode\t\t{weightingMode}")
    ini_file = abspath(join(args.output,IniName))
    print("==> Writing simulation configuration to "+ini_file)
    os.makedirs(os.path.dirname(ini_file), exist_ok=True)
    with open(f"{ini_file}",'w') as f:
        f.write(f"[GeneratorExternal]\nfileName = {generatorName}\nfuncName={generatorFunc}(")
        f.write(f"{collisionSystem},{centrality},{decayMode},{selectedMothers},\"{paramFile}\",\"{paramFileDir}\",{numberOfParticles},{minPt},{maxPt},{pythiaErrorTolerance},{externalDecayer},{decayLongLived},{dynamicalPtRange},{useYWeights},\"{paramV2FileDir}\",{toFixEP},{yGenRange},\"{useLMeeDecaytable}\",{weightingMode}")
        f.write(")")
    return ini_file

def run_only_simulation(args):
    sim_log_name = "simlog"
    name="o2sim"

    output_dir = abspath(args.output)
    output_name = abspath(join(args.output,name))
    ini_file = abspath(args.sim_ini_file)
    log_file = abspath(join(args.output,sim_log_name))

    current_dir = os.getcwd()
    if not exists(output_dir):
        os.makedirs(output_dir)
    os.chdir(output_dir)


    print("\n==> Running cocktail simulation")
    command_to_run = ['o2-sim','-n',args.nEvents,'-g','external','-o',output_name,'--configFile',ini_file,'--noGeant']
    print("running with: ",end='')
    print(*command_to_run)
    print("redirect output to "+log_file)
    print("\nrunning...")
    f = open(log_file, "w")
    simproc = subprocess.Popen(['o2-sim','-n',args.nEvents,'-g','external','-o',output_name,'--configFile',ini_file,'--noGeant'], stdout=f, stderr=subprocess.STDOUT)
    simproc.communicate()
    f.close()
    print("...done\n")

    os.chdir(current_dir)

    if args.clean:
        clean(args.output,name,True)

    if not simproc.returncode==0:
        print("ERROR: o2sim finished with error. See "+abspath(join(args.output,name+"_serverlog"))+" and "+log_file+"for details\n")
        return 1
    else:
        print("Cocktail simulation finished successful. Result in "+output_name+"_Kine.root\n")
        return 0


def run_only_analysis_task(args):
    sim_log_name = "readerlog"
    ana_log_name = "analog"
    name="kineReader"

    print("\n==> Running cocktail analysis task")

    output_dir = abspath(args.output)
    output_name = join(output_dir,name)
    current_dir = os.getcwd()

    kin_file = abspath(args.input)
    json_file = abspath(args.ana_config_file)
    sim_log_file = join(output_dir,sim_log_name)
    ana_log_file = join(output_dir,ana_log_name)

    if not exists(output_dir):
        os.makedirs(output_dir)
    os.chdir(output_dir)

    f = root.TFile(kin_file)
    tree = f.Get("o2sim")
    nEvents = tree.GetEntries()
    print(f"Found {nEvents} events in {kin_file}")

    sim_command_to_run = ['o2-sim','-g','extkinO2','-n',str(nEvents),'-o',output_name,'--extKinFile',kin_file,'--noGeant','--forwardKine','--noDiscOutput']
    print("running reader with: ",end='')
    print(*sim_command_to_run)
    print("redirect reader output to "+sim_log_file)

    fs = open(sim_log_file, "w")
    simproc = subprocess.Popen(sim_command_to_run, stdout=fs, stderr=subprocess.STDOUT)

    proxy_command_to_run = ['o2-sim-mctracks-proxy','--nevents',str(nEvents),'--o2sim-pid',str(simproc.pid)]
    ana_command_to_run = ['o2-analysis-em-lmee-lf-cocktail', '--configuration','json://'+json_file,'-b']
    print("running analysis with: ",end='')
    print(*proxy_command_to_run,end='')
    print(" | ",end='')
    print(*ana_command_to_run)
    print("redirect analysis output to "+ana_log_file)

    fa = open(ana_log_file,"w")
    proxyproc = subprocess.Popen(proxy_command_to_run, stdout=subprocess.PIPE)
    anaproc = subprocess.Popen(ana_command_to_run, stdin=proxyproc.stdout, stdout=fa, stderr=subprocess.STDOUT)

    print("\nrunning...\n")
    simproc.communicate()
    if not simproc.returncode==0:
        print("ERROR: Reader finished with error. See "+abspath(join(args.output,name+"_serverlog"))+" and "+sim_log_file+" for details\n")
        proxyproc.kill();
        anaproc.kill();
    proxyproc.communicate()
    anaproc.communicate()
    fs.close()
    fa.close()
    print("...done\n")


    if not anaproc.returncode==0:
        print("ERROR: Analysis task finished with error. See "+abspath(join(args.output,"analog"))+" for details\n")
    else:
        print("Analysis task finished successful. Histograms in "+output_dir+"/AnalysisResults.root\n")

    os.chdir(current_dir)

    clean(args.output,name,False)

    if (not simproc.returncode) and (not anaproc.returncode):
        return 0
    else:
        return 1

def run_full(args):
    sim_log_name = "simlog"
    ana_log_name = "analog"
    name="o2sim"

    print("\n==> Running cocktail simulation + analysis task")

    output_dir = abspath(args.output)
    output_name = join(output_dir,name)
    current_dir = os.getcwd()

    ini_file = abspath(args.sim_ini_file)
    json_file = abspath(args.ana_config_file)
    sim_log_file = join(output_dir,sim_log_name)
    ana_log_file = join(output_dir,ana_log_name)

    if not exists(output_dir):
        os.makedirs(output_dir)
    os.chdir(output_dir)

    sim_command_to_run = ['o2-sim','-g','external','-n',args.nEvents,'-o',output_name,'--configFile',ini_file,'--noGeant','--forwardKine']
    if not args.save_kine:
        sim_command_to_run.append('--noDiscOutput')
    print("running simulation with: ",end='')
    print(*sim_command_to_run)
    print("redirect sim output to "+sim_log_file)

    fs = open(sim_log_file, "w")
    simproc = subprocess.Popen(sim_command_to_run, stdout=fs, stderr=subprocess.STDOUT)

    proxy_command_to_run = ['o2-sim-mctracks-proxy','--nevents',args.nEvents,'--o2sim-pid',str(simproc.pid)]
    ana_command_to_run = ['o2-analysis-em-lmee-lf-cocktail','--configuration','json://'+json_file,'-b']
    print("running analysis with: ",end='')
    print(*proxy_command_to_run,end='')
    print(" | ",end='')
    print(*ana_command_to_run)
    print("redirect analysis output to "+ana_log_file)

    fa = open(ana_log_file,"w")
    proxyproc = subprocess.Popen(proxy_command_to_run, stdout=subprocess.PIPE)
    anaproc = subprocess.Popen(ana_command_to_run, stdin=proxyproc.stdout, stdout=fa, stderr=subprocess.STDOUT)

    print("\nrunning...\n")
    simproc.communicate()
    if not simproc.returncode==0:
        print("ERROR: o2sim finished with error. See "+abspath(join(args.output,name+"_serverlog"))+" and "+sim_log_file+" for details\n")
        proxyproc.kill();
        anaproc.kill();
    else:
        print("Cocktail simulation finished successful.",end='')
        if args.save_kine:
            print(" Output written to "+output_name+"_Kine.root.",end='')
        print(" Waiting for analysis task to finish...\n")
    proxyproc.communicate()
    anaproc.communicate()
    fs.close()
    fa.close()
    print("...done\n")


    if not anaproc.returncode==0:
        print("ERROR: Analysis task finished with error. See "+abspath(join(args.output,"analog"))+" for details\n")
    else:
        print("Analysis task finished successful. Histograms in "+output_dir+"/AnalysisResults.root\n")

    os.chdir(current_dir)

    if args.clean:
        clean(args.output,name,True)

    if (not simproc.returncode) and (not anaproc.returncode):
        return 0
    else:
        return 1

def main():
    parser = argparse.ArgumentParser()
    common_sim_parser = argparse.ArgumentParser(add_help=False)
    sim_input = common_sim_parser.add_mutually_exclusive_group(required=True)
    common_sim_parser.add_argument("--nEvents","-n",help="Number of events to be generated",required=True)
    sim_input.add_argument("--sim-ini-file",dest="sim_ini_file",help="INI file for the generator")
    sim_input.add_argument("--sim-config-file",dest="sim_config_file",help="JSON config file for the generator (will be converted to INI file)")
    common_sim_parser.add_argument("--output","-o",default=os.getcwd(),help="Output directory")
    common_sim_parser.add_argument("--clean", dest="clean", action="store_true",help="Delete unwanted files after finishing the job")
    common_ana_parser = argparse.ArgumentParser(add_help=False)
    common_ana_parser.add_argument("--ana-config-file",dest="ana_config_file",help="Config file for the analysis task",required=True)
    sub_parsers = parser.add_subparsers(dest="command",required=True)
    sim_parser = sub_parsers.add_parser("sim-only",parents=[common_sim_parser])
    ana_parser = sub_parsers.add_parser("ana-only",parents=[common_ana_parser])
    ana_parser.add_argument("--input","-i",help="Analysis input file",required=True)
    ana_parser.add_argument("--output","-o",default=os.getcwd(),help="Output directory")
    full_parser = sub_parsers.add_parser("full",parents=[common_sim_parser,common_ana_parser])
    full_parser.add_argument("--save-kine",dest="save_kine",action="store_true",help="Write the generator output to o2sim_Kine.root")
    args = parser.parse_args()

    print("\n######################")
    print("## LMee LF Cocktail ##")
    print("######################")

    print("\n==> Running in mode \'"+args.command+"\'\n")

    if ( ((args.command=="full") or (args.command=="sim-only")) and args.sim_config_file):
        args.sim_ini_file=read_sim_config(args)
        if args.sim_ini_file == 1:
            return 1

    if (args.command=="full"):
        return run_full(args)
    if (args.command=="sim-only"):
        return run_only_simulation(args)
    if (args.command=="ana-only"):
        return run_only_analysis_task(args)

main()