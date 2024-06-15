#!/bin/bash

# \brief Bash script to run the azimuthal correlation analysis
# \usage ./AnalysisExecution.sh

# Create the output directory
OUTDIR="Output_AnalysisExecution"
if [ ! -d "$OUTDIR" ]; then
    mkdir Output_AnalysisExecution
fi

# Correlation extraction
root -b -l <<EOF &> Output_AnalysisExecution/stdoutExtractCorrel.log
Printf("[INFO] EXTRACT CORRELATIONS");
.L DhCorrelationExtraction.cxx+
.x ExtractOutputCorrel.C("config_CorrAnalysis_DsToKKPi.json")
.q
EOF

# Correlation fit
root -b -l <<EOF &> Output_AnalysisExecution/stdoutFitCorrel.log
Printf("[INFO] FIT CORRELATIONS");
.L DhCorrelationFitter.cxx+
.x FitCorrel.C("config_CorrAnalysis_DsToKKPi.json")
.q
EOF
