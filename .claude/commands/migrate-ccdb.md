Migrate the specified file (or all files mentioned in the conversation) from `Service<o2::ccdb::BasicCCDBManager>` to the declarative CCDB table approach.

## Background

The old approach uses `Service<o2::ccdb::BasicCCDBManager>` and calls `ccdb->getForTimeStamp<T>(path, timestamp)` at runtime. The new approach declares CCDB columns and timestamped tables using macros, so the framework fetches objects automatically and exposes them as columns on BC rows.

**New API summary:**

```cpp
// In namespace o2::aod (or a sub-namespace):
DECLARE_SOA_CCDB_COLUMN(StructName, getterName, ConcreteType, "CCDB/Object/Path");

DECLARE_SOA_TIMESTAMPED_TABLE(TableName, aod::Timestamps, o2::aod::timestamp::Timestamp, 1, "TABLEDESC",
                              ns::StructName, ns::OtherColumn);

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
- **Non-BC timestamps**: if the timestamp comes from something other than a BC (e.g. computed manually), the migration is non-trivial — flag it instead of forcing it.
- **Global/init-time fetches** (e.g. `efficiencyGlobal.cxx` style): not migratable — the timestamped-table mechanism requires a row in a BC-keyed table.
- **Magnetic-field side effects**: tasks that compute `d_bz` from a fetched `GRPMagField` and seed a propagator can keep that logic, just sourcing the object from `bc.grpMagField()` instead of `ccdb->getForTimeStamp(...)`.

$ARGUMENTS
