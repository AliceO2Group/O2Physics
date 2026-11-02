# Strangeness table converters

Converter tasks migrate old versions of the strangeness derived-data tables to
newer ones, so that analyses written against the current data model can still
run over older derived data. Each task reads one or more legacy tables and
produces the corresponding newer version, filling any column that did not exist
in the old table with a dummy value.

Add the workflow you need to your job with `o2-analysis-lf-<name>`; only the
converters whose *input* tables are actually present in your derived data need
to be enabled.

The table below is generated automatically from the sources.

| Converter | Reads | Produces | Description |
|---|---|---|---|
| `stracentconverter.cxx` | `StraCents_000` | `StraCents_001` | Converts Stra Cents from 000 to 001 |
| `stracentconverter2.cxx` | `StraCents_001` | `StraCents_002` | Converts Stra Cents from 001 to 002 |
| `stradautracksconverter.cxx` | `V0Cores` + `V0Extras` + `V0TOFs`, `CascCores` + `CascExtras` + `CascTOFs`, `DauTrackExtras` *(=DauTrackExtras_003)* | `DauTrackTOFPIDs_000` |  |
| `stradautracksextraconverter.cxx` | `DauTrackExtras_000` | `DauTrackExtras_001` | Converts V0 version 001 to 002 |
| `stradautracksextraconverter2.cxx` | `DauTrackExtras_001` | `DauTrackExtras_002` | Converts daughter TracksExtra from 1 to 2 |
| `stradautracksextraconverter3.cxx` | `DauTrackExtras_002` | `DauTrackExtras_003` | Converts daughter TracksExtra from 2 to 3 |
| `stradautrackstofpidconverter.cxx` | `V0Cores` + `V0Extras` + `V0TOFs`, `CascCores` + `CascExtras` + `CascTOFs`, `DauTrackExtras` *(=DauTrackExtras_003)* | `DauTrackTOFPIDs` *(=DauTrackTOFPIDs_003)* |  |
| `stradautrackstofpidconverter2.cxx` | `StraCollisions`, `DauTrackExtras` *(=DauTrackExtras_003)* + `DauTrackTOFPIDs_000`, `V0CollRefs` + `V0Cores` + `V0Extras` | `DauTrackTOFPIDs_001`<br>`StraEvTimes_000` | converts DauTrackTOFPIDs_000 to _001 |
| `stradautrackstofpidconverter3.cxx` | `DauTrackTOFPIDs_001`, `StraEvTimes_000` | `DauTrackTOFPIDs_002`<br>`StraEvTimes_001` | converts DauTrackTOFPIDs_001 to _002 |
| `stradautrackstofpidconverter4.cxx` | `DauTrackTOFPIDs_002`, `StraCollisions` + `StraStamps` *(=StraStamps_001)* | `DauTrackTOFPIDs_003` | converts DauTrackTOFPIDs_002 to _003 |
| `stradautrackstpcpidconverter.cxx` | `DauTrackTPCPIDs_000` | `DauTrackTPCPIDs_001` | converts DauTrackTOFPIDs_000 to _001 |
| `straevselextrasconverter.cxx` | **processAll**: `StraEvSels_005` + `StraEvSelExtras_000`<br>**processStraEvSelsOnly**: `StraEvSels_005` *(off by default)* | `StraEvSelExtras_001` |  |
| `straevselextrasconverter2.cxx` | `StraEvSels` *(=StraEvSels_006)* | `StraEvSelExtras` *(=StraEvSelExtras_001)* | Produce dummy StraEvSelExtras for analysis subscribing to StraEvSelExtras but not saved in the strangeness derived data (typically when running over pp strangeness derived data) |
| `straevselsconverter.cxx` | `StraEvSels_000` + `StraRawCents_004` | `StraEvSels_001` | Converts Stra Event selections from 000 to 001 |
| `straevselsconverter2.cxx` | `StraEvSels_001` | `StraEvSels_002` | Converts Stra Event selections from 000 to 001 |
| `straevselsconverter2rawcents.cxx` | `StraEvSels_001` | `StraRawCents_004` | Converts V0 version 001 to 002 |
| `straevselsconverter2rawcents2.cxx` | `StraEvSels_002` | `StraRawCents_004` | Converts V0 version 001 to 002 |
| `straevselsconverter2rawcents3.cxx` | `StraEvSels_003` | `StraRawCents_004` | Converts V0 version 001 to 002 |
| `straevselsconverter3.cxx` | `StraEvSels_002` | `StraEvSels_003` | Converts Stra Event selections from 000 to 001 |
| `straevselsconverter4.cxx` | `StraEvSels_003` | `StraEvSels_004` | Converts Stra Event selections from 000 to 001 |
| `straevselsconverter5.cxx` | `StraEvSels_004` + `StraStamps` *(=StraStamps_001)* | `StraEvSels_005` | Converts Stra Event selections from 004 to 005 |
| `straevselsconverter6.cxx` | `StraEvSels_005` | `StraEvSels_006` | Converts Stra Event selections from 005 to 006 |
| `stramccollisionconverter.cxx` | `StraMCCollisions_000` | `StraMCCollisions_001` | Converts V0 version 001 to 002 |
| `stramccollisionconverter2.cxx` | `StraMCCollisions_001` | `StraMCCollisions_002` | Converts V0 version 001 to 002 |
| `stramccollmultconverter.cxx` | `StraMCCollMults_000` | `StraMCCollMults_001` | Converts V0 version 001 to 002 |
| `stramccollmultconverter2.cxx` | `StraMCCollMults_001` | `StraMCCollMults_002` | Converts V0 version 001 to 002 |
| `strarawcentsconverter.cxx` | **process000to001**: `StraRawCents_000` *(off by default)*<br>**process002to003**: `StraRawCents_002` *(off by default)* | `StraRawCents_001`<br>`StraRawCents_003` | Converts V0 version 001 to 002 |
| `strarawcentsconverter2v4.cxx` | `StraRawCents_003` | `StraRawCents_004` | Converts V0 version 001 to 002 |
| `strastampsconverter.cxx` | `StraStamps_000` | `StraStamps_001` | Converts Stra Stamps from 000 to 001 |
| `v0coresconverter.cxx` | `V0MCCores_000` | `V0MCCores_001` | Converts V0 version 001 to 002 |
| `v0coresconverter2.cxx` | `V0MCCores_001` | `V0MCCores_002` | Converts V0 version 001 to 002 |
| `v0mlscoresconverter.cxx` | `V0Cores` | `V0GammaMLScores`<br>`V0LambdaMLScores`<br>`V0AntiLambdaMLScores`<br>`V0K0ShortMLScores` | Converts V0 version 001 to 002 |
| `zdcneutronsconverter.cxx` | `StraEvSels` *(=StraEvSels_006)* | `ZDCNeutrons`<br>`ZDCNMCCollRefs` |  |

Notes:

* Tables joined with `+` are subscribed together in a single argument
  (`soa::Join<...>`); a comma separates independent arguments.
* `Name *(=Name_00N)*` means the task subscribes to (or produces) the
  unversioned alias, which currently resolves to version `N` in the data model.
* Tasks with more than one `process` function list each one separately;
  a switch marked *off by default* must be enabled explicitly.
