Migrate the specified file (or all files mentioned in the conversation) from `Service<o2::ccdb::BasicCCDBManager>` to the declarative CCDB table approach.

## Background

The old approach uses `Service<o2::ccdb::BasicCCDBManager>` and calls `ccdb->getForTimeStamp<T>(path, timestamp)` at runtime. The new approach declares CCDB columns and timestamped tables using macros, so the framework fetches objects automatically and exposes them as columns on BC rows.

**New API summary:**

```cpp
// In namespace o2::aod (or a sub-namespace):
DECLARE_SOA_CCDB_COLUMN(StructName, getterName, ConcreteType, "CCDB/Object/Path");

// ... or, when the object needs fixing up after deserialisation, the _FULL form. The 0 is
// the run-dependence mode; the finaliser is the trailing argument (see the two sections
// "Objects needing post-deserialisation fixup" and "Run-dependent objects"):
DECLARE_SOA_CCDB_COLUMN_FULL(StructName, "fStructName", getterName, ConcreteType, "CCDB/Object/Path", 0,
                             [](ConcreteType* o) { return fixUp(o); });

DECLARE_SOA_TIMESTAMPED_TABLE(TableName, aod::Timestamps, o2::aod::timestamp::Timestamp, 1, "TABLEDESC",
                              ns::StructName, ns::OtherColumn);

// ... or, when the object is constant across something coarser than a timestamp, the
// uniform form (see "Uniformity: how often the object can change"):
DECLARE_SOA_UNIFORM_TABLE(TableName, aod::Timestamps, o2::aod::timestamp::Timestamp,
                          aod::BCs, o2::aod::bc::RunNumber, 1, "TABLEDESC",
                          ns::StructName);

// In the task — basic usage:
using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::TableName>;
void process(MyBCs const& bcs) {
  for (auto const& bc : bcs) {
    auto const& obj = bc.getterName();   // reference to cached deserialized object; treat as immutable
  }
}
```

**Configurable CCDB paths** (`ConfigurableCCDBPath<Column>`):

`ConfigurableCCDBPath<Column>` is a `Configurable<std::string>` named `"ccdb:" + Column::mLabel` (where `mLabel = "f" + StructName`), defaulting to `Column::query`. Declare it **only** when a task consumes a *shared* column and needs a different default path than the shared declaration carries. It does not make the path settable from the task's config block — see Known gaps. No O2Physics task currently declares one; the model migration explains why it deliberately does not (`propagationServiceV2.cxx:70-73`).

Why it is a default and not an override: `ArrowSupport.cxx:641-666` copies the `ConfigParamSpec` — `opt.defaultValue`, not a runtime value — from the task onto `internal-dpl-aod-ccdb` at topology-build time, before any config JSON is read. The task's own `ccdb:fXxx` option is auto-registered (`AnalysisTask.h:575`) but nothing reads it for path resolution.

It is purely declarative: declaring it is sufficient, and the accessor stays `bc.getterName()`. Do **not** log `.value` as "the path in use" — it reflects the task device's own option, which the fetcher never reads, so it prints exactly the value that is being ignored.

Required headers (add if missing): `<Framework/ASoA.h>`, `<Framework/AnalysisDataModel.h>`, `<Framework/Configurable.h>`
Headers to remove (if no longer needed): `<CCDB/BasicCCDBManager.h>`

## What to do

Read the target file(s) and perform the following migration. Do NOT do a complete migration if the patterns are ambiguous or out of scope — instead note what was skipped and why.

### Step 1 — Inventory

Find every `ccdb->getForTimeStamp<T>(path, ts)` call (and variants like `fCCDB->getForTimeStamp`, `mCcdb->getForTimeStamp`). For each call record:
- The concrete C++ type `T`
- The CCDB path string (may be a `Configurable` variable — record the default value and the Configurable's name)
- The timestamp source (BC timestamp, computed value, etc.)
- Where the result is used

Also record the other call shapes, which each need a different decision: `getForRun` (see Known gaps — not the same query), `getSpecific<T>(path, ts, metadata)` (only `runNumber` metadata is expressible, via run-dependence mode 1; anything else — skip), and `get<T>(path)` relying on a previously-set manager timestamp (treat as global/init-time unless that timestamp came from a BC).

**Deduplicate**: for the same (type, path) pair, declare only one CCDB column. Multiple call sites collapse into multiple uses of the same getter.

### Step 2 — Identify scope

Determine whether each fetch is:
- **Per-BC/per-collision** (called inside `process()` with a timestamp from a BC) — these can be migrated
- **Per-run** (called once when `runNumber` changes, guarded by `mRunNumber == ...`) — migratable, but declare the table with `DECLARE_SOA_UNIFORM_TABLE(..., aod::BCs, o2::aod::bc::RunNumber, ...)` so the fetcher queries once per run (see "Uniformity"). The timestamp default queries once per distinct BC timestamp, which is the opposite of what the `mRunNumber` guard was doing.
- **Optional** (the task calls `ccdb->setFatalWhenNull(false)` and branches on a null result — common for MC without a calibration) — a column cannot express absence; a missing object aborts the job. Leave these as-is and note it.
- **Global/init-time** (called in `init()` with a fixed timestamp, not keyed to a BC) — these **cannot** be migrated to CCDB tables; leave them as-is and note this

Skip the migration for any global/init-time fetches. Skip the whole file if all fetches are global.

### Step 3 — Declare CCDB columns and table

**First check whether a shared column already exists**, and join it instead of redeclaring: `GLO/Config/GRPMagField`, `GLO/Calib/MeanVertex`, `GLO/Config/GRPECS`, `GLO/Config/GRPLHCIF` are in `Common/DataModel/GloCCDBObjects.h` (`aod::GloCCDBObjects`, `aod::GrpLHCIFCCDBObjects`); `GLO/Param/MatLUT` is in `aod::GeomCCDBObjects` in the same header; `TPC/Calib/VDriftTgl` in `Common/DataModel/TpcCCDBObjects.h`; TrackTuner calibrations in `Common/DataModel/TrackTunerCCDBObjects.h`. Declare a per-task table only for objects with no shared column.

Otherwise, in the `o2::aod` namespace (or a private sub-namespace inside the file, before the task struct), declare:

```cpp
namespace o2::aod
{
namespace myccdbtask   // use a short, unique snake_case name derived from the task name
{
DECLARE_SOA_CCDB_COLUMN(StructName, getterName, fully::qualified::ConcreteType, "CCDB/Path"); //!
// one per unique (type, path) pair
} // namespace myccdbtask

DECLARE_SOA_TIMESTAMPED_TABLE(MyTaskCCDBObjects, aod::Timestamps, o2::aod::timestamp::Timestamp, 1, "MYTASKCCDB", //!
                              myccdbtask::StructName /*, ... */);
} // namespace o2::aod
```

Before writing the declaration, settle three things per column — each has its own section
below, and getting them wrong is silent rather than loud:

1. **Does the object need fixing up after deserialisation?** If so use `DECLARE_SOA_CCDB_COLUMN_FULL`
   with a finaliser — see "Objects needing post-deserialisation fixup".
2. **How often can the object change?** Timestamp (the default) or run — see "Uniformity: how
   often the object can change". Choose from the object's validity, not from how the old code
   happened to fetch it.
3. **Is the path the same for every run?** If it varies by period, declare the mapping in the
   query string instead of porting the run-range `if/else` — see "Paths that vary by run".

Rules for naming:
- `StructName` / `getterName`: derive from the type name, e.g. `GRPMagField` / `grpMagField`, `MeanVertex` / `meanVertex`
- **The label `f<StructName>` is a workflow-global key**, not a per-table one: it names the path option on the fetcher, and two columns sharing a `StructName` share one path — first-seen wins, no warning (`AnalysisCCDBHelpers.cxx:91-100`, `ArrowSupport.cxx:665`). Harmless while the paths agree, silent divergence otherwise. If your path differs from an existing column of that name, pick a distinct `StructName`.
- Table name: `<TaskStruct>CCDBObjects`, e.g. `SkimmerDalitzEECCDBObjects`
- `_Desc_` string: short ALL-CAPS string unique within the binary (≤ 16 chars to fit the AOD descriptor), e.g. `"DALZECC"`, `"TOFCALIB"`
- Namespace: lowercase snake-case derived from the task name (avoid collisions with other CCDB column namespaces in the file)
- Use the **default value** of any `Configurable` path as the compile-time path in the `DECLARE_SOA_CCDB_COLUMN` macro; if the path has no obvious default, leave a `// TODO: verify path` comment

### Step 4 — Update the task struct

1. **Remove** `Service<o2::ccdb::BasicCCDBManager> ccdb;` (and any variant field name)
2. **Remove** `int mRunNumber;` (or similar run-caching variables) **only if** their sole purpose was to guard CCDB re-fetches
3. **Remove** `ccdb->setURL(...)`, `ccdb->setCaching(...)`, `ccdb->setLocalObjectValidityChecking()`, `ccdb->setCreatedNotAfter(...)`, `ccdb->setFatalWhenNull(...)` from `init()`
4. **Remove** the entire `initCCDB()`/`initMagField()` helper method if it only did CCDB fetching; otherwise remove just the CCDB lines from it
5. **Handle path Configurables** — for each `Configurable<std::string>` that held a CCDB path:
   - If its default equals the compile-time path you put in the column (the normal case): **remove** it. `AnalysisTask.h:575` auto-registers `ccdb:fXxx` on every subscribing task with default `Column::query`, so a `ConfigurableCCDBPath` would only restate that default. Declare `ConfigurableCCDBPath<ns::ColumnName>` only when consuming a *shared* column whose declared path differs from this task's old default, and comment that it sets the fetcher's default, not a runtime value.
   - If the path was never intended to be user-facing (e.g. internal fixed paths): **remove** it outright; the compile-time path in `DECLARE_SOA_CCDB_COLUMN` is sufficient.
   - Always remove Configurables that were only used for CCDB manager setup and not for paths: `ccdb-url`, `ccdb-no-later-than`, `skipGRPOquery`, `d_bz_input` (if only used to bypass CCDB), etc.
6. **Remove** cached pointer member variables (e.g. `GRPMagField* grpmag = nullptr`) if they were only populated by CCDB fetches that are now replaced

### Step 5 — Update process() signatures

Define one alias near the top of the task or just below the table declaration:
```cpp
using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::MyTaskCCDBObjects>;
```

Then for each `process()` that used to call `getForTimeStamp`:

- If `process()` already takes `aod::BCsWithTimestamps const&` directly: change it to `MyBCs const&`.
- If `process()` accesses BCs via `collision.bc_as<aod::BCsWithTimestamps>()`: add `MyBCs const&` to the process signature (so the framework knows to provide it) and replace the `bc_as<>` type with `MyBCs`.
- If `process()` does not currently mention BCs but called `ccdb->getForTimeStamp(path, collision.bc_as<...>().timestamp())`: add `MyBCs const&` to the signature and obtain the BC via `collision.bc_as<MyBCs>()`.
- Replace every `ccdb->getForTimeStamp<T>(path, ts)` call with `bc.getterName()`. The returned reference is to a cached deserialized object; treat it as immutable.
- Null-pointer checks (`if (!grpmag)`) on the result become unnecessary — the framework guarantees the object is present (or the task fails early). Remove them.
- If a helper template like `initCCDB(collision)` was called per-collision, inline its remaining (non-CCDB) work or drop it.
- When configuring from the first BC, guard the empty table first: `if (bcs.size() == 0) { return; }` (`StandardCCDBLoader.h:69`). Empty BC tables occur, and `bcs.begin()` / `iteratorAt(0)` on one is undefined.

### Step 6 — Fix includes

- Remove `#include <CCDB/BasicCCDBManager.h>` if no other code in the file still uses `BasicCCDBManager`
- Ensure `#include <Framework/ASoA.h>` is present (may already be included transitively)
- Keep all type headers (e.g. `<DataFormatsParameters/GRPMagField.h>`) since they are still needed for the concrete type

### Step 7 — Final review

After making changes:
- Check that every remaining use of `ccdb` / `fCCDB` / `mCcdb` has been handled
- Check that `mRunNumber` (or similar) is fully removed if unused
- Check that any leftover `Configurable<std::string>` for a path is either replaced by `ConfigurableCCDBPath<>` or removed
- Search for stale references to removed Configurables (e.g. `grpmagPath.value` lingering in log messages — switch to `grpMagFieldPath.value`)
- If `init()` is now empty, it can be removed
- Note any patterns that were intentionally skipped

## Important limitations — tell the user if any apply

- **Configurable paths**: after migration a path is a compile-time constant plus one option on the fetcher device. It is **not** settable from the task's config block or from a Hyperloop wagon parameter — read the Known gaps entry before telling the user a path is configurable.
- **`getRunDuration()` calls**: these use `BasicCCDBManager` statically and are unrelated to per-BC fetching — do not touch them.
- **`ctpRateFetcher` / other helpers**: out of scope.
- **Multiple tasks in one file**: tasks can share a single CCDB table declaration if they need the same objects; otherwise each task gets its own with a unique `_Desc_`.
- **Non-BC timestamps**: if the timestamp comes from something other than a BC, the migration is non-trivial — flag it instead of forcing it. This is the single most common blocker in practice. `Common/Tools/EventSelectionModule.h:244` computes `ts = sorTimestamp / 2 + eorTimestamp / 2` (mid-run, from `getRunDuration` / `AggregatedRunInfo`) and fetches `EventSelectionParams`, `ITS/Config/AlpideParam`, `TriggerAliases` and `ITS/Calib/TimeDeadMap` at it. A BC-keyed column fetches at each BC's own timestamp instead, so migrating these silently changes which object version is served whenever an object is revised mid-run. They need a run-uniform table before they can move — and note that one queries at the run's first BC, not mid-run, so it is equivalent only for objects with a single version per run.
- **Global/init-time fetches** (e.g. `efficiencyGlobal.cxx` style): not migratable — the timestamped-table mechanism requires a row in a BC-keyed table.
- **Magnetic-field side effects**: tasks that compute `d_bz` from a fetched `GRPMagField` and seed a propagator can keep that logic, just sourcing the object from `bc.grpMagField()` instead of `ccdb->getForTimeStamp(...)`.

## Lessons learned (established in-tree, with references)

### Why this migration matters beyond tidiness

The old per-task path Configurable is device-name-scoped, which made it a silent-divergence trap: `propagationServiceV2` shared `ccdb.lutPath` with `propagationService` (`Common/Tools/StandardCCDBLoader.h:45`, default `GLO/Param/MatLUT`), but every config setting `GLO/Param/MatLUTInner` did so under a `propagation-service` block (8 of them, all under `Tutorials/`), so V2 silently used the full LUT — different material corrections, no warning.

Migration makes the path one option on the fetcher device and turns two *tasks* declaring conflicting defaults into a warning (`ArrowSupport.cxx:641-666`). That is the real gain. It does not fix the scoping itself: **runtime** overrides remain device-scoped, just aimed at `internal-dpl-aod-ccdb` instead of the task — see Known gaps before telling a user a migrated path is "now configurable".

### Objects needing post-deserialisation fixup

Some objects are not usable straight out of the ROOT streamer. `MatLayerCylSet` is a `FlatObject`: its internal pointers are unfixed and its voxel lookup unbuilt until `MatLayerCylSet::rectifyPtrFromFile()` runs. Use the `_FULL` form, which carries the finaliser (the plain `DECLARE_SOA_CCDB_COLUMN` passes an identity one):

```cpp
DECLARE_SOA_CCDB_COLUMN_FULL(MatLUT, "fMatLUT", matLUT, o2::base::MatLayerCylSet, "GLO/Param/MatLUT", 0, //!
                             [](o2::base::MatLayerCylSet* lut) { return o2::base::MatLayerCylSet::rectifyPtrFromFile(lut); });
```

The argument before the finaliser is the run-dependence mode (`0` = query by timestamp, the default the plain macro passes; see "Run-dependent objects"). The finaliser must be the **last** macro argument (commas in a lambda body are absorbed by `__VA_ARGS__`), has signature `T* (*)(T*)`, and runs on the receiving device once per (re)deserialisation, before the object is ever handed out. Ownership contract: whatever it returns is what the column cache later `delete`s, so a finaliser returning a *different* instance must dispose of the one it was given.

Do **not** put this fixup in the task. There is no `finaliseCCDB` hook on the analysis path (`adaptAnalysisTask` wires only `EndOfStream`, `AnalysisTask.h:610-619`; grep confirms zero uses of `finaliseCCDB` in O2Physics).

### Uniformity: how often the object can change

Every CCDB table declares a *uniformity column*: rows sharing its value resolve to the same
object, so the fetcher queries once per distinct value instead of once per row.
`DECLARE_SOA_TIMESTAMPED_TABLE` defaults it to the timestamp column, which is the
pre-existing behaviour — every distinct timestamp may yield a different object.

Pick it from the object's real validity, and only then:

| Object changes ... | Uniformity | Declare with |
| --- | --- | --- |
| within a run (calibrations, drift velocity) | timestamp (default) | `DECLARE_SOA_TIMESTAMPED_TABLE` |
| per run or per period (geometry, material, per-period calibrations) | `aod::BCs` / `aod::bc::RunNumber` | `DECLARE_SOA_UNIFORM_TABLE` |

Worked examples in the tree: `aod::TpcCalibCCDBObjects` keeps the timestamp default because
the TPC drift velocity genuinely varies within a run; `aod::GeomCCDBObjects` and
`aod::TrackTunerCCDBObjects` are run-uniform.

The run number lives on `aod::BCs` and the timestamp on `aod::Timestamps`; the macro takes both,
and the fetcher fatals if they are not row-aligned (anything joinable with the BCs is fine). The
uniformity column must be an int32/int64/uint64 (`AnalysisCCDBHelpers.cxx:246-260`).

### Paths that vary by run: declare a mapping, not code

A column's path may be a plain path, or a mapping from uniformity value to path:

```
"520259-529691=…/pp2023/pass4/vsPhi;559348-559387=…/ppRef/polarity_positive;fallback"
```

Ranges are inclusive; either bound may be omitted (`-hi=path`, `lo-=path`); entries are
separated by `;`; an entry without `=` is an explicit fallback. **A value matching no range
is fatal**, deliberately — silently substituting another period's calibration is the failure
mode this whole mechanism exists to prevent. A string with no `=` is a plain path, so
existing columns are unaffected.

**The key matched against the ranges is the table's uniformity value**, so a run-range mapping
requires `DECLARE_SOA_UNIFORM_TABLE` on `aod::bc::RunNumber`. On a timestamp-uniform table the
ranges match nothing: fatal if there is no fallback entry, and — worse — every row silently
takes the fallback if there is one.

The mapping is *data* in the schema metadata, so the run ranges stop being compiled in: the
whole mapping is replaceable through the `ccdb:fXxx` option, subject to the Known gaps scoping.

This replaces hand-written run-range tables. `TrackTuner::getPathInputFileAutomaticFromCCDB()`
is the model case: ~50 lines of `else if (lo <= runNumber && runNumber <= hi)` became the
declaration in `Common/DataModel/TrackTunerCCDBObjects.h`. When porting one, **derive the
mapping mechanically and diff it against the source** — first-match-wins must reproduce the
`if/else` order, which matters whenever ranges overlap (in TrackTuner, one PbPb range sits
inside a pp range and must stay *after* it).

### Serving migrated and un-migrated callers from one module

Shared modules must keep working for tasks that have not migrated. Detect the capability
rather than adding a configuration flag:

```cpp
auto const& bc = collision.template bc_as<TBCs>();
if constexpr (requires { bc.vdriftTgl(); }) {
  mVDriftMgr.update(bc.vdriftTgl());          // column path
} else {
  mVDriftMgr.update(bc.timestamp());          // legacy CCDB query
}
```

The discarded branch is not instantiated, so an un-migrated caller compiles exactly as before
and a migrated one never references the CCDB manager. `strangenessBuilderModule::updateVDrift`
uses this. Where a whole function parameter falls away, add an overload of different arity
that forwards (see "Shared module signatures") and put a `static_assert` with a readable
message on the ccdb-free one, so calling it with an unjoined BC table names the missing table
instead of failing somewhere inside the template.

### Two path settings must never both be live

After migration the column is the single source of truth for a path. If the task still has an
old `Configurable<std::string>` for the same object, **fail loudly when both are set** rather
than silently preferring one — that divergence is exactly the bug this migration exists to
kill. `TrackPropagationModule::init` fatals when `trackTuner.pathInputFile` is non-empty while
the calibrations come from columns, naming the option to use instead (`ccdb:fTrackTunerDca`).

Caveat: this test only works for Configurables whose default is empty. One with a non-empty
default cannot be distinguished from an unset one, so that hole stays open until the framework
can report whether an option was explicitly set.

### Grouping columns into tables

One table per **family of objects used together with similar validity intervals** — not one per consuming task. Geometry and material description (`GLO/Param/MatLUT`, and later `GLO/Config/GeometryAligned`, `GLO/Config/Geometry`, `<DET>/Calib/Align`; see `GRPGeomRequest` in `O2/Detectors/Base/src/GRPGeomHelper.cxx:44-60`) is one family with essentially static validity. The GRP family changes per run, and `GRPMagField` is requested per timeframe in O2 (`GRPGeomHelper.cxx:72`). Splitting on that boundary keeps a task from fetching a multi-hundred-MB LUT it never asked for.

**Several timestamped tables can be joined onto the same BCs.** `soa::Join<aod::BCsWithTimestamps, aod::GloCCDBObjects, aod::GeomCCDBObjects>` works: the duplicated `aod::Timestamps` is deduplicated when `originals` is merged (`ASoA.h:131-140`), giving 4 originals, and every accessor resolves. Do not invent per-use-case tables to work around a limitation that does not exist.

### Global state is not a lookup

Migrating removes CCDB *queries*, not side effects. Two things stay:

- `Propagator::initFieldFromGRP()` rebuilds or rescales a `MagneticField`, attaches it to `TGeoGlobalMagField::Instance()` and locks it (`O2/Detectors/Base/src/Propagator.cxx:107-149`). Keep it guarded on run change.
- `Propagator::Instance()->setMatLUT()` is a pointer store, so it is cheaper to redo unconditionally every timeframe — and doing so picks up a relocated column buffer for free instead of dangling.

Everything else (mean vertex, run number) should become a direct read at the point of use, with no cached member and no `initCCDB()` helper. A cached pointer plus a "did the buffer move?" check is strictly worse than reading the column fresh.

`Propagator` cannot itself become a column value: private constructor, deleted copy/move, singleton `Instance()` (`Propagator.h:157-201`).

### Shared module signatures

If a shared module takes a `StandardCCDBLoader`, change it to take the values it actually uses (`int runNumber`, `MeanVertexObject const*`) and keep a thin forwarding overload for un-migrated callers, so V1 tasks stay byte-identical. `TrackPropagationModule::fillTrackTables` does this — the two overloads differ in arity, so overload resolution is unambiguous.

### What the migration does and does not buy

The fetcher downloads once into a shm cache and the column stores `(handle, segment, size)` (`AnalysisCCDBHelpers.cxx:320-330`). What is shared is the **serialised blob**; each consumer still streams its own heap copy in the column getter. So expect fewer downloads, one configuration point and cross-device consistency — but not a per-device RSS reduction. For a `FlatObject` like the LUT, real memory sharing needs a zero-copy path (`FlatObject::setActualBufferAddress`) that does not exist yet.

Measured on `propagationServiceV2` (`HF_LHC23_pass4_Thin_small_2P3PDstar`, `daily-20260922-0000-1`), migrated wagon 61527 against un-migrated baseline wagon 61526:

| | baseline 766059 | migrated 766456 |
| --- | --- | --- |
| propagation device `cpuUsedAbsolute` | 3,972,078 | 2,866,181 (**−28 %**) |
| propagation device peak PSS | 160 MB | 168 MB |
| train peak PSS | 1340 MB | 1366 MB |
| CCDB queries from the task device | 3.03 MB fetched | 0 — all via the fetcher |

Quote the CPU number when justifying a migration; do not promise memory.

### Known gaps in the mechanism

- **Runtime path overrides are scoped to the fetcher device, so most existing config mechanisms cannot reach them.** `ccdb:fXxx` is registered as an option on `internal-dpl-aod-ccdb`, and DPL resolves options strictly per device: `DataProcessingDevice.cxx:2615` builds the retriever keyed on the device `name`, and `ConfigurationOptionsRetriever.cxx:47` does `getRecursive(mMainKey)` with no global or wildcard section. A value placed in any other device's block is unreachable. Hyperloop writes wagon parameters into the wagon's own device block, and the config JSONs in the tree are written the same way, so **the natural way to set a migrated path silently does nothing**. The four levers that do work today: pass `--ccdb:fXxx` on the workflow command line (DPL forwards CLI options to every device declaring them, which is what `propagationServiceV2.cxx:73` tells users to do); hand-write an `internal-dpl-aod-ccdb` block into the config JSON; edit the compile-time path in the column declaration; or declare `ConfigurableCCDBPath` in a task to change that default — but note that ArrowSupport takes the *first* device in workflow order carrying the option, and every subscriber carries it with the compile-time default, so a `ConfigurableCCDBPath` in a task that is not the first subscriber loses to that default with only a warning. Making this usable from Hyperloop needs a Hyperloop-side change (emit `ccdb:*` into an `internal-dpl-aod-ccdb` block); making it usable generally needs the fetcher to consult the declaring task's runtime value, which the build-time propagation in `ArrowSupport.cxx:641-666` cannot do.

  Observed end-to-end (Hyperloop test 766459): the wagon carried `ccdb:fMatLUT = GLO/Param/MatLUTInner` and `configuration.json` shipped it, inside the `propagation-service-v2` block; `dpl-config.json` then held the key **twice**, `MatLUTInner` under `propagation-service-v2` and `MatLUT` under `internal-dpl-aod-ccdb`. The job fetched the full LUT.

  **Verify a path override by payload, never by config** — reading it back from the wagon or from `configuration.json` only proves it was *stored*. Three checks prove it was *used*: the fetcher's `[internal-dpl-aod-ccdb] ccdb reads <url>` log lines, the object UUID in that URL, and `reader.ccdb-cache-fetched-bytes` in `performanceMetrics_processed.json`. Grepping for the path *name* is not one of them.

- **Run-dependent objects need an explicit mode, and no in-tree column sets one.** The run-dependence argument of `DECLARE_SOA_CCDB_COLUMN_FULL` (`0` = by timestamp; `1` additionally sends the run as `runNumber` metadata; `2` queries by run instead of timestamp — `CCDBFetcherHelper.cxx:189-198`) *is* honoured: the fetcher reads it per field, takes the run from the uniformity value, and fatals if the table is not uniform in `fRunNumber` (`AnalysisCCDBHelpers.cxx:288-303`). What is missing is use: every in-tree column is mode 0. In particular `GLO/Config/GRPECS` is requested run-dependently in O2 (`GRPGeomHelper.cxx:66`) but `ccdbGlo::GRPECSObject` is mode 0 inside the timestamp-uniform `aod::GloCCDBObjects` — do not rely on it for run-keyed lookups until it moves to a run-uniform table with mode 1.
- **`getForRun` is not the same query.** `BasicCCDBManager::getForRun` resolves the run duration and queries at *mid-run* (`BasicCCDBManager.h:364-381`); a column queries at each BC's timestamp. Identical for objects with one version per run, divergent otherwise.
- **Row cardinality, not query count.** The uniformity column already collapses the *queries* to one per distinct value, but the table still carries one row per BC per column — a `FixedSizeList<int64,3>`, 24 B, rebuilt every timeframe. Collapsing the rows too needs a non-extension table plus lookup by value at the consumer, which does not exist yet. So a run-uniform table costs the same arrow memory as before; what it saves is the fetching.
- **Multi-run dataframes.** Skimmed datasets can span runs. Every existing consumer configures from `bcs.begin()` and applies it to the whole DF (`propagationServiceV2.cxx`, `StandardCCDBLoader.h:70-77`, `strangenessBuilderModule.h:884`), which is wrong for such a DF. Migrating preserves this bug unless it is fixed deliberately — do not claim the migration fixes it.

### Practical gotchas

- `DECLARE_SOA_CCDB_COLUMN` expands to a non-template member using `TClass::GetClass` and `TBufferFile`, and `ASoA.h` neither includes nor forward-declares them. The header that expands the macro must have `<TClass.h>` and `<TBufferFile.h>` in scope — include them explicitly (`// IWYU pragma: keep`) rather than relying on `AnalysisDataModel.h` to drag them in transitively.
- The getter holds one deserialised object per column *type*, replaced (and the previous one `delete`d) as soon as a row with a different payload is read (`ASoA.h:2413-2432`). A reference from `bc1.x()` dangles after `bc2.x()` when the two rows resolve to different objects — another run in a multi-run DF, or a mid-run revision. Read at the point of use; never hold the reference across rows or timeframes. This is also why a pointer stored into global state (`setMatLUT`) must be refreshed every timeframe.
- A failed fetch is fatal, not silent: if `extractCCDBPayload` returns null the getter aborts naming the type, the path and the `ccdb:` option to check. A mistyped path therefore stops the job rather than dereferencing null.
- Do not add a `sources` member to a table's metadata struct. It makes the struct satisfy both `soa::with_sources` and `soa::with_sources_generator`, and `getInputMetadata` becomes ambiguous.
- Device options are matched by device *name*. Never look a task's own option up by a hardcoded name (`device.name == "propagation-service"` silently matched nothing in `propagation-service-v2`); take the running device from `initContext.services().get<DeviceSpec const>()`. Spell the type out rather than using `auto`, or the pre-existing `option.defaultValue.get<bool>()` becomes a dependent name and needs `template`.
- Verify with the *control*: when changing a shared header, compile an un-migrated consumer too. A new error appearing in both is yours; the same errors in both means you changed nothing for them.

$ARGUMENTS
