Migrate the specified file (or all files mentioned in the conversation) from `Service<o2::ccdb::BasicCCDBManager>` to the declarative CCDB table approach.

## Background

The old approach uses `Service<o2::ccdb::BasicCCDBManager>` and calls `ccdb->getForTimeStamp<T>(path, timestamp)` at runtime. The new approach declares CCDB columns and timestamped tables using macros, so the framework fetches objects automatically and exposes them as columns on BC rows.

**New API summary:**

```cpp
// In namespace o2::aod (or a sub-namespace):
DECLARE_SOA_CCDB_COLUMN(StructName, getterName, ConcreteType, "CCDB/Object/Path");

// ... or, when the object needs fixing up after deserialisation, the _FULL form, whose
// trailing argument is the finaliser (see "Objects needing post-deserialisation fixup"):
DECLARE_SOA_CCDB_COLUMN_FULL(StructName, "fStructName", getterName, ConcreteType, "CCDB/Object/Path",
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

If the original task used a `Configurable<std::string>` to supply the CCDB path, the path can remain user-overridable after migration using `ConfigurableCCDBPath<Column>`. This is a typed `Configurable<std::string>` whose option name is automatically set to `"ccdb:" + Column::mLabel` (where `mLabel = "f" + StructName`), defaulting to the compile-time path in the column declaration. The framework reads this option name when resolving CCDB URLs, so users can still redirect the path via JSON config.

The Configurable is purely declarative for the path-override mechanism: declaring it is sufficient — you do **not** pass `.value` to a getter or fetcher. The accessor remains `bc.getterName()`. The `.value` member is still available if the task wants to log the resolved path.

```cpp
struct MyTask {
  // Replaces: Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "..."};
  ConfigurableCCDBPath<ns::GRPMagField> grpMagFieldPath;   // option name = "ccdb:fGRPMagField"

  void process(MyBCs const& bcs) {
    auto const& grpmag = bcs.iteratorAt(0).grpMagField();   // path override is honoured automatically
    LOGP(info, "Using GRPMagField path: {}", grpMagFieldPath.value);
  }
};
```

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

**Deduplicate**: for the same (type, path) pair, declare only one CCDB column. Multiple call sites collapse into multiple uses of the same getter.

### Step 2 — Identify scope

Determine whether each fetch is:
- **Per-BC/per-collision** (called inside `process()` with a timestamp from a BC) — these can be migrated
- **Per-run** (called once when `runNumber` changes, guarded by `mRunNumber == ...`) — these can be migrated; the framework caches per unique timestamp automatically
- **Global/init-time** (called in `init()` with a fixed timestamp, not keyed to a BC) — these **cannot** be migrated to CCDB tables; leave them as-is and note this

Skip the migration for any global/init-time fetches. Skip the whole file if all fetches are global.

### Step 3 — Declare CCDB columns and table

In the `o2::aod` namespace (or a private sub-namespace inside the file, before the task struct), declare:

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
   - If the path was used as the sole argument to `getForTimeStamp` and the user may want to override it at runtime: **replace** it with `ConfigurableCCDBPath<ns::ColumnName>` (e.g. `ConfigurableCCDBPath<ns::GRPMagField> grpMagFieldPath;`). The member name should match the getter for clarity. Keep a comment explaining what path it controls.
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

- **Configurable paths**: CCDB column paths are compile-time constants in the macro. Add `ConfigurableCCDBPath<Column>` to allow runtime override; its default is `Column::query` so it always agrees with the macro by construction.
- **`getRunDuration()` calls**: these use `BasicCCDBManager` statically and are unrelated to per-BC fetching — do not touch them.
- **`ctpRateFetcher` / other helpers**: out of scope.
- **Multiple tasks in one file**: tasks can share a single CCDB table declaration if they need the same objects; otherwise each task gets its own with a unique `_Desc_`.
- **Non-BC timestamps**: if the timestamp comes from something other than a BC, the migration is non-trivial — flag it instead of forcing it. This is the single most common blocker in practice. `Common/Tools/EventSelectionModule.h:243` computes `ts = sorTimestamp / 2 + eorTimestamp / 2` (mid-run, from `getRunDuration` / `AggregatedRunInfo`) and fetches `EventSelectionParams`, `ITS/Config/AlpideParam`, `TriggerAliases` and `ITS/Calib/TimeDeadMap` at it. A BC-keyed column fetches at each BC's own timestamp instead, so migrating these silently changes which object version is served whenever an object is revised mid-run. They need a run-keyed table before they can move.
- **Global/init-time fetches** (e.g. `efficiencyGlobal.cxx` style): not migratable — the timestamped-table mechanism requires a row in a BC-keyed table.
- **Magnetic-field side effects**: tasks that compute `d_bz` from a fetched `GRPMagField` and seed a propagator can keep that logic, just sourcing the object from `bc.grpMagField()` instead of `ccdb->getForTimeStamp(...)`.

## Lessons learned (established in-tree, with references)

### Why this migration matters beyond tidiness

The per-task path Configurable is a silent-divergence trap. `propagationService` and `propagationServiceV2` share the identical `ccdb.lutPath` Configurable (`Common/Tools/StandardCCDBLoader.h:45`, default `GLO/Param/MatLUT`), but config JSONs key overrides by *device name*. Every config in the tree carries a `propagation-service` block setting `GLO/Param/MatLUTInner` and no `propagation-service-v2` block, so V2 silently fell back to the full LUT — different material corrections, no warning. After migration the path is one option on the fetcher device, and two tasks disagreeing produces a warning (`ArrowSupport.cxx:641-666`) instead of silence.

### Objects needing post-deserialisation fixup

Some objects are not usable straight out of the ROOT streamer. `MatLayerCylSet` is a `FlatObject`: its internal pointers are unfixed and its voxel lookup unbuilt until `MatLayerCylSet::rectifyPtrFromFile()` runs. Use the `_FULL` form, which carries the finaliser (the plain `DECLARE_SOA_CCDB_COLUMN` passes an identity one):

```cpp
DECLARE_SOA_CCDB_COLUMN_FULL(MatLUT, "fMatLUT", matLUT, o2::base::MatLayerCylSet, "GLO/Param/MatLUT", //!
                             [](o2::base::MatLayerCylSet* lut) { return o2::base::MatLayerCylSet::rectifyPtrFromFile(lut); });
```

The finaliser must be the **last** macro argument (commas in a lambda body are absorbed by `__VA_ARGS__`), has signature `T* (*)(T*)`, and runs on the receiving device once per (re)deserialisation, before the object is ever handed out. Ownership contract: whatever it returns is what the column cache later `delete`s, so a finaliser returning a *different* instance must dispose of the one it was given.

Do **not** put this fixup in the task. There is no `finaliseCCDB` hook on the analysis path (`adaptAnalysisTask` wires only `EndOfStream`, `AnalysisTask.h:610-619`; grep confirms zero uses of `finaliseCCDB` in O2Physics), and even if there were, an opt-in hook means a task that forgets it gets a silently broken object.

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

Two consequences worth knowing before choosing:

- The uniformity column may live in a **different table** from the timestamp — the run number
  is on `aod::BCs`, the timestamp on `aod::Timestamps`. Both are handed to the fetcher
  automatically (the table's `generateSources()` merges their originals) and read positionally.
- Positional reading is only sound if the two sources are **row-aligned**. ASoA encodes no
  type-level relation between tables that merely have equal row counts, so this cannot be a
  `static_assert`; the fetcher compares the two column lengths and fatals on a mismatch.
  Anything joinable with the BCs is fine.

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

The mapping is *data*, carried in the schema metadata. That matters: the CCDB fetcher is a
separate device and must not depend on code from the task that declared the column, so a
resolver lambda would not do. It also means the run ranges stop being compiled in — the whole
mapping is replaceable at runtime through the `ccdb:fXxx` option.

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

**Several timestamped tables can be joined onto the same BCs.** `soa::Join<aod::BCsWithTimestamps, aod::GloCCDBObjects, aod::GeomCCDBObjects>` works: the duplicated `aod::Timestamps` is deduplicated when `originals` is merged (`ASoA.h:172-186`), giving 4 originals, and every accessor resolves. Do not invent per-use-case tables to work around a limitation that does not exist.

### Global state is not a lookup

Migrating removes CCDB *queries*, not side effects. Two things stay:

- `Propagator::initFieldFromGRP()` rebuilds or rescales a `MagneticField`, attaches it to `TGeoGlobalMagField::Instance()` and locks it (`O2/Detectors/Base/src/Propagator.cxx:107-149`). Keep it guarded on run change.
- `Propagator::Instance()->setMatLUT()` is a pointer store, so it is cheaper to redo unconditionally every timeframe — and doing so picks up a relocated column buffer for free instead of dangling.

Everything else (mean vertex, run number) should become a direct read at the point of use, with no cached member and no `initCCDB()` helper. A cached pointer plus a "did the buffer move?" check is strictly worse than reading the column fresh.

`Propagator` cannot itself become a column value: private constructor, deleted copy/move, singleton `Instance()` (`Propagator.h:157-201`).

### Shared module signatures

If a shared module takes a `StandardCCDBLoader`, change it to take the values it actually uses (`int runNumber`, `MeanVertexObject const*`) and keep a thin forwarding overload for un-migrated callers, so V1 tasks stay byte-identical. `TrackPropagationModule::fillTrackTables` does this — the two overloads differ in arity, so overload resolution is unambiguous.

### What the migration does and does not buy

The fetcher downloads once into a shm cache and the column stores `(handle, segment, size)` (`AnalysisCCDBHelpers.cxx:213-222`). What is shared is the **serialised blob**; each consumer still streams its own heap copy in the column getter. So expect fewer downloads, one configuration point and cross-device consistency — but not a per-device RSS reduction. For a `FlatObject` like the LUT, real memory sharing needs a zero-copy path (`FlatObject::setActualBufferAddress`) that does not exist yet.

### Known gaps in the mechanism

- **Run-dependent objects are not served correctly.** The analysis fetcher still hardcodes `.runNumber = 1, .runDependent = 0` for every column, even though `CCDBFetcherHelper.cxx:189-195` implements the run-dependent query paths. `GLO/Config/GRPECS` is marked "Run dependent !!!" in O2 and already has a column — verify before relying on it. Now that a run-uniform table gives the fetcher a run number per row, wiring this through is small and worth doing.
- **`getForRun` is not the same query.** `BasicCCDBManager::getForRun` resolves the run duration and queries at *mid-run* (`BasicCCDBManager.h:364-374`); a column queries at each BC's timestamp. Identical for objects with one version per run, divergent otherwise.
- **Row cardinality, not query count.** The uniformity column already collapses the *queries* to one per distinct value, but the table still carries one row per BC per column — a `FixedSizeList<int64,3>`, 24 B, rebuilt every timeframe. Collapsing the rows too needs a non-extension table plus lookup by value at the consumer, which does not exist yet. So a run-uniform table costs the same arrow memory as before; what it saves is the fetching.
- **Multi-run dataframes.** Skimmed datasets can span runs. Every existing consumer configures from `bcs.begin()` and applies it to the whole DF (`propagationServiceV2.cxx`, `StandardCCDBLoader.h:70-77`, `strangenessBuilderModule.h:850`), which is wrong for such a DF. Migrating preserves this bug unless it is fixed deliberately — do not claim the migration fixes it.

### Practical gotchas

- `DECLARE_SOA_CCDB_COLUMN` expands to code using `TClass` and `TBufferFile`, but `ASoA.h` only sees them forward-declared. A translation unit that includes the column header without otherwise pulling in `<TClass.h>` and `<TBufferFile.h>` fails to compile. Include them if needed.
- A failed fetch is fatal, not silent: if `extractCCDBPayload` returns null the getter aborts naming the type, the path and the `ccdb:` option to check. A mistyped path therefore stops the job rather than dereferencing null.
- Do not add a `sources` member to a table's metadata struct. It makes the struct satisfy both `soa::with_sources` and `soa::with_sources_generator`, and `getInputMetadata` becomes ambiguous.
- Device options are matched by device *name*. Never look a task's own option up by a hardcoded name (`device.name == "propagation-service"` silently matched nothing in `propagation-service-v2`); take the running device from `initContext.services().get<DeviceSpec const>()`. Spell the type out rather than using `auto`, or the pre-existing `option.defaultValue.get<bool>()` becomes a dependent name and needs `template`.
- Verify with the *control*: when changing a shared header, compile an un-migrated consumer too. A new error appearing in both is yours; the same errors in both means you changed nothing for them.

$ARGUMENTS
