# ALICE O2Physics PID Feature Extractor

A comprehensive C++ task for the O2Physics framework that extracts particle identification (PID) features from ALICE AO2D files. This tool processes track data and generates high-quality features suitable for machine learning-based particle identification.

## Overview

The **PIDFeatureExtractor** combines information from multiple ALICE detectors (TPC and TOF) to create a rich feature set for distinguishing between different particle types (pions, kaons, protons, and electrons). The extracted features are saved in both ROOT TTree and CSV formats for easy access and analysis.

### Key Features

- **Multi-Detector Integration**: Combines TPC (dE/dx) and TOF (time-of-flight) PID information
- **CCDB Integration**: Automatically fetches unavailable features from the ALICE Conditions Database (CCDB)
- **Bayesian Probability Computation**: Calculates combined PID probabilities using Gaussian likelihood in n-sigma space
- **Flexible Output**: Exports to both ROOT TTree and CSV formats simultaneously
- **Quality Control**: Includes QC histograms for track kinematics and detector response
- **MC Truth Matching**: Includes Monte Carlo truth information for simulated data
- **Configurable Selection**: User-adjustable kinematic cuts (pT, η range) via JSON configuration
- **Track Quality Metrics**: Stores DCA and TPC fit quality information

## Requirements

### O2Physics Environment

- **O2Physics framework**: Latest version with PID response tables and CCDB access
- **ROOT**: Version 6.x or later
- **CMake**: 3.x or later
- **C++ Standard**: C++17 or later
- **CCDB Access**: Network connection to ALICE CCDB for fetching PID calibrations
- **bash**: For running the `run.sh` execution script

### Required Data Tables

The task expects the following input tables from AO2D files. Some tables may be fetched from CCDB if not present in the file:

| Table | Source | Purpose | Fallback |
|-------|--------|---------|----------|
| `aod::Tracks` | AO2D | Base track properties (momentum, angles) | Required |
| `aod::TracksExtra` | AO2D | Extended track information | Required |
| `aod::TracksDCA` | AO2D | Impact parameters (DCA) | Required |
| `aod::pidTPCPi/Ka/Pr/El` | AO2D/CCDB | TPC n-sigma values for each particle species | CCDB |
| `aod::pidTOFPi/Ka/Pr/El` | AO2D/CCDB | TOF n-sigma values for each particle species | CCDB |
| `aod::pidTOFmass` | AO2D/CCDB | TOF reconstructed mass | CCDB |
| `aod::pidTOFbeta` | AO2D/CCDB | TOF beta (v/c) measurement | CCDB |
| `aod::McTrackLabels` | AO2D | MC truth matching (optional, for simulated data) | Optional |
| `aod::McParticles` | AO2D | MC particle information | Optional |

**Note:** If PID tables are not available in the AO2D file, the framework automatically retrieves PID calibrations from CCDB using the collision timestamp to access the correct calibration period.

## Extracted Features

### Kinematic Variables (11 features)

| Variable | Type | Range | Unit | Description |
|----------|------|-------|------|-------------|
| `event_id` | int | - | - | Unique collision event identifier |
| `track_id` | int | - | - | Track index within event |
| `px`, `py`, `pz` | float | - | GeV/c | Cartesian momentum components |
| `pt` | float | 0.1-20 | GeV/c | Transverse momentum |
| `p` | float | - | GeV/c | Total momentum |
| `eta` | float | -1.5 to 1.5 | - | Pseudorapidity |
| `phi` | float | -π to π | rad | Azimuthal angle |
| `theta` | float | 0 to π | rad | Polar angle |
| `charge` | int | ±1 | - | Track charge |
| `track_type` | int | 0-2 | - | Track classification |

### TPC Detector Features (7 features)

| Variable | Type | Range | Unit | Description | Source |
|----------|------|-------|------|-------------|--------|
| `tpc_signal` | float | 0-300 | - | Specific ionization (dE/dx) | AO2D |
| `tpc_nsigma_pi` | float | - | σ | n-sigma deviation from pion | AO2D/CCDB |
| `tpc_nsigma_ka` | float | - | σ | n-sigma deviation from kaon | AO2D/CCDB |
| `tpc_nsigma_pr` | float | - | σ | n-sigma deviation from proton | AO2D/CCDB |
| `tpc_nsigma_el` | float | - | σ | n-sigma deviation from electron | AO2D/CCDB |
| `tpc_nclusters` | int | 0-160 | - | Number of TPC clusters | AO2D |
| `tpc_chi2` | float | - | - | TPC track fit chi-square/ndf | AO2D |

**TPC Features Source:** n-sigma values are computed from `tpc_signal` and PID calibrations (from AO2D or CCDB). If not in AO2D, calibration data is fetched from CCDB using the collision timestamp.

### TOF Detector Features (6 features)

| Variable | Type | Range | Unit | Description | Source |
|----------|------|-------|------|-------------|--------|
| `tof_beta` | float | 0-1.2 | - | Velocity over speed of light | AO2D/CCDB |
| `tof_mass` | float | -0.2-2.0 | GeV/c² | Reconstructed mass | AO2D/CCDB |
| `tof_nsigma_pi` | float | - | σ | n-sigma deviation from pion | AO2D/CCDB |
| `tof_nsigma_ka` | float | - | σ | n-sigma deviation from kaon | AO2D/CCDB |
| `tof_nsigma_pr` | float | - | σ | n-sigma deviation from proton | AO2D/CCDB |
| `tof_nsigma_el` | float | - | σ | n-sigma deviation from electron | AO2D/CCDB |

**TOF Features Source:** If not available in AO2D file, the framework fetches calibration and response parameters from CCDB. Beta and mass can be recomputed from raw TOF information and length measurement using CCDB calibrations.

### Bayesian PID Features (4 features)

| Variable | Type | Range | Unit | Description |
|----------|------|-------|------|-------------|
| `bayes_prob_pi` | float | 0-1 | - | Probability of being pion |
| `bayes_prob_ka` | float | 0-1 | - | Probability of being kaon |
| `bayes_prob_pr` | float | 0-1 | - | Probability of being proton |
| `bayes_prob_el` | float | 0-1 | - | Probability of being electron |

**Note**: Bayesian probabilities sum to 1.0 and are computed using Gaussian likelihoods in n-sigma space (from either AO2D or CCDB-derived values) with configurable priors.

### Track Quality Features (2 features)

| Variable | Type | Unit | Description |
|----------|------|------|-------------|
| `dca_xy` | float | cm | Distance of closest approach in xy-plane |
| `dca_z` | float | cm | Distance of closest approach along beam |

### Detector Availability Flags (2 features)

| Variable | Type | Description |
|----------|------|-------------|
| `has_tpc` | bool | Track has valid TPC information |
| `has_tof` | bool | Track has valid TOF information |

### Monte Carlo Truth (4 features, simulated data only)

| Variable | Type | Description |
|----------|------|-------------|
| `mc_pdg` | int | PDG code of true particle |
| `mc_px`, `mc_py`, `mc_pz` | float | True momentum components |

**Total: 39 features per track**

## Installation

### 1. Clone the Repository

```bash
cd ~/O2Physics  # or your O2Physics installation directory
git clone <repository-url> pid-extractor
cd pid-extractor
```

### 2. Verify Directory Structure

Ensure your repository has the following structure:

```
pid-extractor/
├── CMakeLists.txt
├── PIDFeatureExtractor.cxx
├── myConfigExtractor.json
├── run.sh
└── README.md
```

### 3. Set Executable Permissions

Make the `run.sh` script executable:

```bash
chmod +x run.sh
```

## Configuration

### Configuration File: `myConfigExtractor.json`

All task parameters are configured through the **`myConfigExtractor.json`** file located in the task directory. This JSON file specifies all runtime options for the PID feature extractor.

#### Configuration File Format

```json
{
  "output_path": "pid_features",
  "export_csv": true,
  "export_root": true,
  "eta_min": -1.5,
  "eta_max": 1.5,
  "pt_min": 0.1,
  "pt_max": 20.0,
  "ccdb_url": "http://alice-ccdb.cern.ch"
}
```

#### Configuration Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_path` | string | `pid_features` | Base path for output files (without extension) |
| `export_csv` | boolean | `true` | Enable CSV export of features |
| `export_root` | boolean | `true` | Enable ROOT file export of features |
| `eta_min` | float | `-1.5` | Minimum pseudorapidity cut for track selection |
| `eta_max` | float | `1.5` | Maximum pseudorapidity cut for track selection |
| `pt_min` | float | `0.1` | Minimum transverse momentum cut (GeV/c) |
| `pt_max` | float | `20.0` | Maximum transverse momentum cut (GeV/c) |
| `ccdb_url` | string | `http://alice-ccdb.cern.ch` | CCDB server URL for fetching PID calibrations |

#### Example Configurations

**Example 1: High-pT kaon extraction**

```json
{
  "output_path": "kaons_highpt",
  "export_csv": true,
  "export_root": true,
  "eta_min": -0.8,
  "eta_max": 0.8,
  "pt_min": 3.0,
  "pt_max": 20.0,
  "ccdb_url": "http://alice-ccdb.cern.ch"
}
```

**Example 2: ROOT-only output for large datasets**

```json
{
  "output_path": "batch_output",
  "export_csv": false,
  "export_root": true,
  "eta_min": -1.5,
  "eta_max": 1.5,
  "pt_min": 0.1,
  "pt_max": 20.0,
  "ccdb_url": "http://alice-ccdb.cern.ch"
}
```

**Example 3: Using alternative CCDB**

```json
{
  "output_path": "ccdb_test",
  "export_csv": true,
  "export_root": true,
  "eta_min": -1.5,
  "eta_max": 1.5,
  "pt_min": 0.1,
  "pt_max": 20.0,
  "ccdb_url": "http://alice-ccdb-test.cern.ch"
}
```

## Usage

### Quick Start

The task runs within the O2Physics framework using the provided **`run.sh`** execution script. All configuration is read from **`myConfigExtractor.json`**.

#### Basic Execution

```bash
./run.sh
```

This command:
1. Reads configuration from `myConfigExtractor.json`
2. Initializes the O2Physics environment
3. Launches the PID feature extractor
4. Processes AO2D data according to configured parameters
5. Generates ROOT and/or CSV output files

#### With AO2D File Input

```bash
./run.sh --aod-file AO2D.root
```

#### With Multiple Input Files

```bash
./run.sh --aod-file file1.root file2.root file3.root
```

#### Troubleshooting Execution

If you encounter permission errors:

```bash
bash ./run.sh
```

### Running in O2Physics Environment

The `run.sh` script should be executed within an active O2Physics environment. To ensure proper setup:

```bash
# Source O2Physics environment (if not already sourced)
source ~/O2Physics/setup.sh

# Make script executable
chmod +x run.sh

# Run the task
./run.sh
```

### Modifying Configuration

#### Method 1: Edit Configuration File

Modify `myConfigExtractor.json` before running:

```bash
# Edit the configuration file
nano myConfigExtractor.json

# Run with new configuration
./run.sh
```

#### Method 2: Environment Variables

You can override configuration parameters via environment variables (if supported by `run.sh`):

```bash
export OUTPUT_PATH="custom_output"
export PT_MIN=2.0
export PT_MAX=10.0
./run.sh
```

### Advanced Usage Examples

**Example 1: Process test data with verbose output**

```bash
# Edit myConfigExtractor.json for test parameters
./run.sh --aod-file test_data.root --verbose
```

**Example 2: High-pT analysis**

1. Update `myConfigExtractor.json`:
```json
{
  "output_path": "high_pt_pions",
  "pt_min": 5.0,
  "pt_max": 20.0,
  "export_csv": true,
  "export_root": true
}
```

2. Run the script:
```bash
./run.sh --aod-file physics_data.root
```

**Example 3: Batch processing multiple files**

```bash
# Configure for batch output
cat > myConfigExtractor.json << EOF
{
  "output_path": "batch_results",
  "export_csv": false,
  "export_root": true,
  "eta_min": -0.9,
  "eta_max": 0.9,
  "pt_min": 0.5,
  "pt_max": 10.0
}
EOF

# Run on multiple files
./run.sh --aod-file data_run1.root data_run2.root data_run3.root
```

## Data Flow and CCDB Integration

### Processing Pipeline

```
AO2D File
    ↓
Load Configuration from myConfigExtractor.json
    ↓
Initialize O2Physics Workflow
    ↓
Load Tracks + Available PID Tables
    ↓
Missing Features? → Query CCDB (using timestamp)
    ↓
Compute Bayesian PID (using AO2D or CCDB n-sigma values)
    ↓
Apply Track Selection (eta, pT cuts from config)
    ↓
Fill QC Histograms
    ↓
Write to ROOT TTree + CSV (based on config)
```

### CCDB Timestamp-Based Access

The task uses the collision timestamp from each event to query CCDB for the correct calibration period:

1. **Event Timestamp**: Read from collision data
2. **CCDB Query**: Fetch PID calibrations valid for that timestamp
3. **Feature Computation**: Use CCDB calibrations if AO2D values unavailable
4. **Caching**: Calibrations are cached per run to minimize CCDB queries

### Calibration Objects Retrieved from CCDB

If not in AO2D, the task fetches:

- **TPC PID Response**: Bethe-Bloch curve parameters and n-sigma calculation
- **TOF PID Response**: TOF expected times and mass resolution
- **Track Propagation**: Path length and TOF time calibrations
- **Detector Resolution**: TPC and TOF resolution parameters for n-sigma computation

## Output Files

### ROOT Output (`<outputPath>.root`)

Contains a TTree named `pid_features` with one entry per track. The tree includes all 39 features as separate branches. All feature values are either from AO2D or computed using CCDB calibrations.

**Access in ROOT:**

```cpp
TFile file("pid_features.root");
TTree* tree = (TTree*)file.Get("pid_features");

// Plot TPC dE/dx vs pT (TPC signal always from AO2D)
tree->Draw("tpc_signal:pt", "has_tpc==1", "scatter");

// Select pions with high confidence (probabilities computed from CCDB if needed)
tree->Draw("eta", "bayes_prob_pi > 0.9", "hist");
```

### CSV Output (`<outputPath>.csv`)

Plain text comma-separated values format with header row. All 39 features are included. Features are either from AO2D or computed using CCDB calibrations. The filename is determined by the `output_path` parameter in `myConfigExtractor.json`.

**Example CSV structure:**

```csv
event_id,track_id,px,py,pz,pt,p,eta,phi,theta,charge,track_type,...
0,0,1.234,-0.567,2.345,1.456,2.678,-0.123,1.456,2.345,1,0,...
0,1,0.987,0.654,1.234,1.123,1.567,0.456,2.123,0.987,-1,0,...
```

### Quality Control Histograms

The task also produces QC histograms in the ROOT file:

- `QC/nTracks`: Total number of processed tracks
- `QC/pt`: pT distribution
- `QC/eta`: Pseudorapidity distribution
- `QC/tpc_dEdx_vs_pt`: TPC dE/dx vs pT (2D)
- `QC/tof_beta_vs_p`: TOF beta vs momentum (2D)
- `QC/mass_vs_p`: TOF mass vs momentum (2D)

## Machine Learning Integration

### Data Loading in Python

```python
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import xgboost as xgb

# Load CSV file (features are from AO2D or CCDB-derived)
# Output filename based on output_path from myConfigExtractor.json
df = pd.read_csv("pid_features.csv")

# Filter out invalid data (missing detector info)
df_valid = df[(df['has_tpc'] == True) & (df['has_tof'] == True)]

# Prepare features (exclude MC truth for real data)
feature_cols = [col for col in df.columns 
                if col not in ['event_id', 'track_id', 'mc_pdg', 'mc_px', 'mc_py', 'mc_pz']]

X = df_valid[feature_cols].values
y = df_valid['mc_pdg'].values  # For simulated data

# Standardize features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Train model
model = xgb.XGBClassifier(n_estimators=100, max_depth=6)
model.fit(X_scaled, y)
```

### Loading ROOT File in Python

```python
import uproot
import pandas as pd

# Open ROOT file (output filename based on output_path from config)
file = uproot.open("pid_features.root")
tree = file["pid_features"]

# Convert to pandas DataFrame
df = tree.arrays(library="pd")

# Access specific branches
bayes_probs = df[['bayes_prob_pi', 'bayes_prob_ka', 'bayes_prob_pr', 'bayes_prob_el']]
```

## Algorithm Details

### Bayesian PID Calculation

The task computes Bayesian probabilities combining TPC and TOF information. n-sigma values are either from AO2D or computed using CCDB calibrations:

**Likelihood Calculation:**

For each particle hypothesis (π, K, p, e):

```
L_i = exp(-0.5 * (ns_TPC_i² + ns_TOF_i²))
```

Where `ns_TPC_i` and `ns_TOF_i` are n-sigma deviations from expected values (from AO2D or CCDB).

**Prior Probabilities:**

Default priors (configurable in source code):
- Pions: 1.0
- Kaons: 0.2
- Protons: 0.1
- Electrons: 0.05

**Final Probabilities:**

```
P(i|TPC,TOF) = (L_i * Prior_i) / Σ(L_j * Prior_j)
```

This ensures the four probabilities sum to exactly 1.0.

### Kinematic Calculations

- **Transverse momentum**: \(p_T = \sqrt{p_x^2 + p_y^2}\)
- **Total momentum**: \(p = \sqrt{p_x^2 + p_y^2 + p_z^2}\)
- **Pseudorapidity**: \(\eta = -\ln(\tan(\theta/2))\)
- **Polar angle**: \(\theta = 2 \arctan(e^{-\eta})\)

### Invalid Data Handling

- Missing TPC: All TPC variables set to `-999`
- Missing TOF: All TOF variables set to `-999`
- Missing CCDB connection: Task logs warning and uses available AO2D values
- No MC match: `mc_pdg` set to `0`, momentum components set to `0`
- Invalid TOF n-sigma (NaN): Treated as 0 contribution in Bayesian calculation

## Troubleshooting

### Issue: `./run.sh: command not found` or permission denied

**Causes:**
- Script is not executable
- Running from wrong directory
- Shell incompatibility

**Solution:**

```bash
# Make executable
chmod +x run.sh

# Run explicitly with bash
bash ./run.sh

# Or run from correct directory
cd ~/pid-extractor
./run.sh
```

### Issue: Configuration file not found

**Error Message:** `myConfigExtractor.json not found`

**Causes:**
- File doesn't exist in current directory
- Running from wrong location

**Solution:**

```bash
# Verify file exists
ls -la myConfigExtractor.json

# Ensure you're in the correct directory
cd ~/pid-extractor
./run.sh
```

### Issue: CCDB Connection Error

**Error Message:** `Cannot connect to CCDB at <url>`

**Causes:**
- Network connectivity issues
- CCDB server down
- Incorrect CCDB URL in configuration

**Solution:**

1. Verify CCDB connectivity:
```bash
curl http://alice-ccdb.cern.ch/
```

2. Update `myConfigExtractor.json` with alternative CCDB:
```json
{
  "ccdb_url": "http://alice-ccdb-test.cern.ch"
}
```

3. Run again:
```bash
./run.sh
```

### Issue: Invalid JSON configuration

**Error Message:** `JSON parse error` or similar

**Causes:**
- Syntax error in `myConfigExtractor.json`
- Missing quotes or commas
- Invalid data types

**Solution:**

1. Validate JSON syntax:
```bash
python3 -m json.tool myConfigExtractor.json
```

2. Fix any reported errors and retry:
```bash
./run.sh
```

### Issue: Features appear as -999 or invalid

**Cause:** PID tables not in AO2D and CCDB fetch failed

**Solution:**
- Check network/CCDB access
- Verify collision timestamp is within valid range
- Ensure AO2D file contains basic track information

### Issue: No tracks extracted

**Causes:**
- Track selection cuts too restrictive (check `eta_min`, `eta_max`, `pt_min`, `pt_max` in config)
- Input data doesn't contain required track tables

**Solution:**

Edit `myConfigExtractor.json` to relax selection criteria:

```json
{
  "eta_min": -1.5,
  "eta_max": 1.5,
  "pt_min": 0.1,
  "pt_max": 20.0
}
```

Then run:
```bash
./run.sh
```

### Issue: Many -999 values in TPC/TOF output

**Cause:** Tracks lack TPC or TOF information (detector not active for those tracks)

**Solution:**
- Use track selection filters to require both detectors:
```python
df_valid = df[(df['has_tpc'] == True) & (df['has_tof'] == True)]
```

### Issue: CSV file is very large

**Cause:** Processing too many tracks or both output formats enabled

**Solutions:**

Modify `myConfigExtractor.json`:

```json
{
  "export_csv": false,
  "export_root": true
}
```

Or apply stricter cuts:

```json
{
  "pt_min": 2.0,
  "pt_max": 10.0,
  "eta_min": -0.8,
  "eta_max": 0.8
}
```

Then run:
```bash
./run.sh
```

### Issue: ROOT TTree has no entries

**Cause:** Track selection filtered out all tracks

**Solution:**
- Verify input data quality
- Reduce selection criteria in `myConfigExtractor.json`
- Check that tracks pass kinematic cuts
- Verify CCDB calibrations loaded successfully

## Performance Notes

- **Processing Speed**: ~10,000 tracks/second (depends on system and CCDB latency)
- **Memory Usage**: ~500 MB for typical dataset
- **File Size**: ~50-100 bytes per track in CSV format
- **CCDB Overhead**: Minimal after initial calibration load per run (~1-2 seconds)

### Optimisation Tips

1. **Use ROOT format for large datasets** (set `export_csv: false` in config)
2. **Process in batches** for better memory management
3. **Apply kinematic cuts** in `myConfigExtractor.json` to reduce output size
4. **Disable MC truth** for real data to save space
5. **Reuse CCDB connections** when processing multiple files from same run

## Contributing

### Adding New Features

To add a new feature:

1. Declare member variable in the struct
2. Create TTree branch in `init()`
3. Add CSV header column
4. Fill variable in `process()`
5. Update documentation

Example:

```cpp
// 1. Member variable
float my_new_feature;

// 2. In init()
featureTree->Branch("my_new_feature", &my_new_feature);

// 3. In process()
my_new_feature = t.someNewMethod();

// 4. Update CSV
csvFile << my_new_feature << ",";
```

### Modifying Configuration Options

To add a new configurable parameter:

1. Add parameter to `myConfigExtractor.json`
2. Create `Configurable<T>` variable in task struct
3. Use parameter value in processing logic
4. Document in this README

### Modifying Bayesian Priors

Edit the priors array in the `process()` function of the task:

```cpp
float priors[4] = {1.f, 0.2f, 0.1f, 0.05f};  // π, K, p, e
```

### Fetching Additional CCDB Objects

If you need additional calibrations:

1. Query CCDB for the object path
2. Load in the `init()` function using ccdbApi
3. Cache calibration data for the processing period
4. Use in `process()` function for feature computation

## Repository Structure

```
pid-extractor/
├── README.md                          # This file - complete documentation
├── PIDFeatureExtractor.cxx            # Main task implementation
├── myConfigExtractor.json             # Configuration file (edit this for options)
├── run.sh                             # Execution script (run with: ./run.sh)
├── CMakeLists.txt                     # Build configuration
└── LICENSE                            # License file
```

## References

- [ALICE O2Physics Framework](https://github.com/AliceO2Group/O2Physics)
- [ALICE Conditions Database (CCDB)](https://alice-ccdb.cern.ch/)
- [ALICE Detector Performance](https://arxiv.org/abs/1910.14400)
- [PID with Machine Learning in ALICE](https://arxiv.org/abs/2204.13255)
- [O2Physics CCDB Integration Guide](https://github.com/AliceO2Group/O2Physics/wiki/CCDB)

## Author

High-Energy Physics Collaboration - ALICE Experiment

## License

MIT License (or appropriate license for your institution)

## Support

For issues, questions, or suggestions:
- Open an issue on GitHub
- Contact the ALICE physics working group
- Check existing documentation
- Consult CCDB documentation for calibration issues

## Quick Reference

### To run the task:
```bash
./run.sh
```

### To modify parameters:
1. Edit `myConfigExtractor.json`
2. Run `./run.sh`

### To see what configuration is active:
```bash
cat myConfigExtractor.json
```

### To verify output:
```bash
ls -lh pid_features.root pid_features.csv
```

---

**Last Updated:** 2025-11-13  
**Task Version:** 1.0.0  
**O2Physics Compatibility:** Latest  
**CCDB Support:** Integrated  
**Run Method:** `./run.sh` with `myConfigExtractor.json`
