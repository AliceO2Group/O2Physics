#!/bin/awk -f

# Fix include style.
# NB: Run before sorting.

{
    if ($1 ~ /^#include/) {
        h = substr($2, 2, length($2) - 2)
        if ( h ~ /^(PWG[A-Z]{2}|Common|ALICE3|DPG|EventFiltering|Tools|Tutorials)\/.*\.h/ ) { $2 = "\""h"\"" } # O2Physics
        else if ( h ~ /^(Algorithm|CCDB|Common[A-Z]|DataFormats|DCAFitter|Detectors|EMCAL|Field|Framework|FT0|FV0|GlobalTracking|GPU|ITS|MathUtils|MFT|MCH|MID|PHOS|PID|ReconstructionDataFormats|SimulationDataFormat|TOF|TPC|ZDC).*\/.*\.h/ ) { $2 = "<"h">" } # O2
        else if ( h ~ /^(T[A-Z]|Math\/|Roo[A-Z])[[:alnum:]\/]+\.h/ ) { $2 = "<"h">" } # ROOT
        else if ( h ~ /^KF[A-Z][[:alnum:]]+\.h/ ) { $2 = "<"h">" } # KFParticle
        else if ( h ~ /^(fastjet\/|onnxruntime)/ ) { $2 = "<"h">" } # FastJet, ONNX
        else if ( h ~ /^.*DataModel\// ) { $2 = "\""h"\"" } # incomplete path to DataModel
        else if ( h ~ /^([[:alnum:]_]+\/)+[[:alnum:]_]+\.h/ ) { $2 = "<"h">" } # other third-party
        else if ( $2 ~ /^".*\./ ) { } # other local-looking file
        else if ( h ~ /^[[:lower:]_]+\.h/ ) { $2 = "<"h">" } # C system
        else if ( h ~ /^[[:lower:]_\/]+/ ) { $2 = "<"h">" } # C++ system
    }
    print
}
