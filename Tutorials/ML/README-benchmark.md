# ONNX inference benchmark tutorial

`benchmarkOnnxModel.cxx` registers `o2-analysis-benchmark-onnx-model`.
It runs once with deterministic synthetic data, writes one CSV row, and exits.
It accepts one dense floating-point input of shape `[batch, features]` and one
or more outputs. Static dimensions must match the options; symbolic dimensions
are filled from the options. Output types are left to ONNX Runtime.

## Build and AMD prerequisites

Build O2Physics normally after loading its O2 environment. This target links
`O2::ML`, `O2::GPUCommon`, `ONNXRuntime::ONNXRuntime`, and the platform dynamic
loader. No HIP headers are required for building the task. The O2 installation
must contain `ML/3rdparty/GPUORTFloat16.h` and its GPUCommon headers.

For GPU runs, ONNX Runtime must contain **MIGraphXExecutionProvider**, and its
MIGraphX/ROCm shared libraries, including `libamdhip64.so`, must be discoverable
at runtime. The task checks provider availability and disables CPU fallback:
unsupported graphs fail rather than producing a misleading GPU benchmark.
No changes to O2's OrtInterface are needed for this task; the current checked-out
OrtInterface does not attach an AMD execution provider, so this example uses
ORT directly.

Check the actual hydra driver, ROCm, MIGraphX and ORT combination first. MI100
(gfx908) and MI50 (gfx906) do not have identical software support. MI50 is absent
from current supported ROCm hardware lists, so it needs a compatible older
stack or a separately validated build. The old ROCm execution provider was
removed from ONNX Runtime starting with 1.23; this task uses MIGraphX.

References:
- https://onnxruntime.ai/docs/execution-providers/MIGraphX-ExecutionProvider.html
- https://onnxruntime.ai/docs/execution-providers/ROCm-ExecutionProvider.html
- https://rocm.docs.amd.com/projects/install-on-linux/en/latest/reference/system-requirements.html

## Run the three comparisons

Use equivalent FP32 and FP16 exports of the same network. Replace the example
feature count with the actual model input width. Use separate output files;
each invocation overwrites its output CSV.

```bash
o2-analysis-benchmark-onnx-model -b --model-path model-fp32.onnx \
  --provider cpu --fp16 false --batch-size 1024 --features 5 \
  --cpu-threads 1 --warmup 20 --iterations 200 --output-csv cpu-fp32.csv

o2-analysis-benchmark-onnx-model -b --model-path model-fp32.onnx \
  --provider migraphx --gpu-id 0 --fp16 false --batch-size 1024 --features 5 \
  --warmup 20 --iterations 200 --output-csv gpu-fp32.csv

o2-analysis-benchmark-onnx-model -b --model-path model-fp16.onnx \
  --provider migraphx --gpu-id 0 --fp16 true --batch-size 1024 --features 5 \
  --warmup 20 --iterations 200 --output-csv gpu-fp16.csv
```

`gpu-id` is the HIP logical ID after any `HIP_VISIBLE_DEVICES` or
`ROCR_VISIBLE_DEVICES` masking. Run one benchmark process at a time on an idle
GPU. Repeat with the same batch sizes and CPU thread count, and compare several
runs. Small batches can be dominated by host overhead and transfers.

`--fp16 true` selects O2's `o2::OrtDataType::Float16_t` input storage and checks
that the model input is FLOAT16. It does not convert model weights or prove that
every internal operation is FP16. Export/inspect the model accordingly. FP16
models retaining FP32 input interfaces need an adapted input path; this small
tutorial intentionally rejects them. CPU runs require FP32 input.

## What the CSV measures

- `samples_per_second = batch_size * iterations / seconds`.
- Timing covers synchronous `Session::Run`, input/output transfers, and output
  allocation/destruction. It excludes session creation, MIGraphX compilation,
  warmup, input generation, FP16 conversion, telemetry, and CSV I/O. This is
  host-to-host inference throughput, not GPU kernel-only throughput.
- GPU memory fields are **device-wide used bytes**, sampled through
  `hipMemGetInfo` before session creation, after warmup, and after timed runs.
  They include other processes and runtime allocations; they are neither
  model-only memory nor a transient peak. The baseline follows HIP context
  initialization. CPU runs leave these fields empty. GPU telemetry failures
  fail the run rather than reporting zero usage.
- Precision labels describe the requested model/input variant, not measured
  arithmetic precision. Strict CPU fallback rejection can reject graphs that
  would otherwise run with some CPU operators.

This benchmark does not check prediction accuracy. Compare real model outputs
between the exports separately before interpreting performance gains.

## Adding the AO2D input later

Planned input on hydra:

```text
/lustre/alice/users/marin/R3A/TPCTrees/TPCTrees_LHC24as_250224/AO2D_merge_LHC24as_updatedBB_250224.root
```

The file has not been inspected here. Do not pass `--aod-file` to this standalone
synthetic task yet. First establish the AO2D tables, feature columns and order,
training normalization, selections, and missing-value policy. Then add an
analysis `process(...)` input with the appropriate table definitions, build
contiguous batches from its rows, and reuse the inference path. Keep session
creation in `init`, aggregate timing/sample counts across batches, and write
results at end of stream. Handle the final partial batch explicitly (dynamic
batch, or pad and count only valid samples). Keep data reading/preprocessing
outside inference timing or report a separate full-pipeline metric.
