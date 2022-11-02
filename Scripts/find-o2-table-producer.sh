#!/bin/bash

# Find the workflow that produces a given table.

# Check that we are inside the O2 or the O2Physics directory.
[[ "$PWD/" != *"/O2"*"/"* ]] && { echo "You must be inside the O2 or the O2Physics directory."; exit 1; }
[ ! "$1" ] && { echo "Provide a table name."; exit 1; }
# Find files that produce the table.
table="$1"
echo "Table $table is produced in:"
files=$(grep -r -i --include="*.cxx" -E "<$table>|<aod::$table>|<o2::aod::$table>" | grep -E 'Produces<|Spawns<|Builds<' | cut -d: -f1 | sort -u)
for f in $files; do
  # Extract the workflow name from the CMakeLists.txt in the same directory.
  wf=$(grep -B 1 "$(basename "$f")" "$(dirname "$f")/CMakeLists.txt" | head -n 1 | cut -d\( -f2)
  echo "$wf in $f"
done
