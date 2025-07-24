#!/usr/bin/env python3

# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

"""!
@brief  Find dependencies required to produce a given table or to run a given workflow.

This scripts finds out which workflows can produce a given table and which tables are needed to run a given workflow.
These dependencies are reported in the terminal in the form of a nested tree
and optionally visualised graphically in a figure.
- Depending on the specified maximum number of levels of the workflow tree, one can request any dependency depth,
  from the direct dependencies to the full workflow topology leading to the provided table(s) and/or workflow(s).
- By default, all possible ways of producing a table and all possible ways of running a workflow are considered.
  This can be customised by providing a list of tables and workflows that should be excluded from the search.
- Since the "Couldn't get TTree" error message in O2 reports the lowercase name of the missing table,
  table names are treated as case-insensitive by default
  so that one can just copy-paste the name from the error message.

@author Vít Kučera <vit.kucera@cern.ch>, Inha University
@date   2022-12-10
"""

import argparse
import glob
import hashlib
import json
import os
import subprocess as sp  # nosec B404
import sys


def eprint(*args, **kwargs):
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def msg_err(message: str):
    """Print an error message."""
    eprint("\x1b[1;31mError: %s\x1b[0m" % message)


def msg_fatal(message: str):
    """Print an error message and exit."""
    msg_err(message)
    sys.exit(1)


def msg_warn(message: str):
    """Print a warning message."""
    eprint("\x1b[1;36mWarning:\x1b[0m %s" % message)


def load_workflows_from_json():
    """Load all workflows from JSON files."""
    try:
        dir_o2p = os.environ["O2PHYSICS_ROOT"]
    except KeyError:
        msg_fatal("O2Physics environment is not loaded.")
    dir_json = f"{dir_o2p}/share/dpl"
    db_wf = {}
    for file_json in glob.glob(f"{dir_json}/*.json"):
        # print(file_json)
        # Get the workflow name from the JSON file name
        workflow = os.path.basename(file_json).split(".")[0]
        try:
            with open(file_json, "r", encoding="utf8", errors="ignore") as j:
                specs_wf = json.load(j)
                db_wf[workflow] = specs_wf["workflow"]
        except FileNotFoundError:
            msg_fatal("JSON file not found.")
        # print(workflow)
    return db_wf


def format_table_name(description: str, subspec: int):
    """Format table description name, including potential versions."""
    if not subspec:
        return description
    return f"{description}_{subspec:03d}"


def get_devices(specs_wf: dict):
    """Get the list of devices of a given workflow loaded from a JSON file."""
    return [d["name"] for d in specs_wf]


def get_inputs(specs_wf: dict, device=""):
    """Get the list of input tables of a given workflow loaded from a JSON file.
    If a device name is provided, only inputs of that device are considered."""
    list_inputs = []
    for dev in specs_wf:
        if device and dev["name"] != device:
            continue
        list_inputs += [
            format_table_name(i["description"], i["subspec"]) for i in dev["inputs"] if i["origin"] == "AOD"
        ]
    return list(dict.fromkeys(list_inputs))  # Remove duplicities


def get_outputs(specs_wf: dict, device=""):
    """Get the list of output tables of a given workflow loaded from a JSON file.
    If a device name is provided, only outputs of that device are considered."""
    list_outputs = []
    # Loop over devices
    for dev in specs_wf:
        if device and dev["name"] != device:
            continue
        list_outputs += [
            format_table_name(i["description"], i["subspec"]) for i in dev["outputs"] if i["origin"] == "AOD"
        ]
    return list(dict.fromkeys(list_outputs))  # Remove duplicities


def print_workflows(dic_wf_all: dict, list_wf=None):
    """Print properties of a given workflow or workflows in the simplified dictionary.
    If no workflow name is provided, print all."""

    def print_wf(dic_wf: dict, wf: str):
        """Print a single workflow"""
        print(wf)
        # print(dic_wf)
        # Loop over devices
        for dev, dic_dev in dic_wf.items():
            print(f"  device: {dev}")
            print(f"    inputs:  {dic_dev['inputs']}")
            print(f"    outputs: {dic_dev['outputs']}")

    if list_wf:
        for wf in list_wf:
            print_wf(dic_wf_all[wf], wf)
    else:
        for wf, dic_wf in dic_wf_all.items():
            print_wf(dic_wf, wf)


def get_table_producers(table: str, dic_wf_all: dict, case_sensitive=False, reverse=False):
    """Find all workflows that have this table as output."""
    list_producers = []
    if not case_sensitive:
        table = table.lower()
    # Loop over workflows
    for wf, dic_wf in dic_wf_all.items():
        # Loop over devices
        for dev in dic_wf:
            outputs = [o if case_sensitive else o.lower() for o in dic_wf[dev]["inputs" if reverse else "outputs"]]
            if table in outputs:
                list_producers.append(wf)
    return list(dict.fromkeys(list_producers))  # Remove duplicities


def get_workflow_outputs(wf: str, dic_wf_all: dict):
    """Get list of workflow outputs from the simplified dictionary."""
    list_outputs = []
    # Loop over devices
    for dev in dic_wf_all[wf]:
        list_outputs += dic_wf_all[wf][dev]["outputs"]
    return list(dict.fromkeys(list_outputs))  # Remove duplicities


def get_workflow_inputs(wf: str, dic_wf_all: dict):
    """Get list of workflow inputs from the simplified dictionary."""
    list_inputs = []
    # list_outputs = get_workflow_outputs(wf, dic_wf_all)
    # Loop over devices
    for dev in dic_wf_all[wf]:
        list_inputs += dic_wf_all[wf][dev]["inputs"]
    list_inputs = list(dict.fromkeys(list_inputs))  # Remove duplicities
    # l = [d for d in l if d not in list_outputs] # avoid circular dependence (can be legit)
    return list_inputs


def get_tree_for_workflow(
    wf: str, dic_wf_all: dict, dic_wf_tree=None, case_sensitive=False, level=0, levels_max=0, reverse=False
):
    """Get the dependency tree of tables and workflows needed to run this workflow."""
    # print(level, levels_max)
    if dic_wf_tree is None:
        dic_wf_tree = {}
    if wf not in dic_wf_all:
        msg_fatal(f"Workflow {wf} not found")
    if wf not in dic_wf_tree:
        dic_wf_tree[wf] = dic_wf_all[wf]
    if reverse:
        inputs = get_workflow_outputs(wf, dic_wf_all)
        symbol_direction = "->"
    else:
        inputs = get_workflow_inputs(wf, dic_wf_all)
        symbol_direction = "<-"
    if inputs:
        print(f"{level * '    '}{wf} {symbol_direction} {inputs}")
        if levels_max < 0 or level < levels_max:
            for tab in inputs:
                producers = get_table_producers(tab, dic_wf_all, case_sensitive, reverse)
                if producers:
                    print(f"{level * '    ' + '  '}{tab} {symbol_direction} {producers}")
                    for p in producers:
                        if p not in dic_wf_tree:  # avoid infinite recursion
                            get_tree_for_workflow(
                                p, dic_wf_all, dic_wf_tree, case_sensitive, level + 1, levels_max, reverse
                            )
    return dic_wf_tree


def get_tree_for_table(tab: str, dic_wf_all: dict, dic_wf_tree=None, case_sensitive=False, levels_max=0, reverse=False):
    """Get the dependency tree of tables and workflows needed to produce this table."""
    if dic_wf_tree is None:
        dic_wf_tree = {}
    producers = get_table_producers(tab, dic_wf_all, case_sensitive, reverse)
    symbol_direction = "<-"
    if reverse:
        symbol_direction = "->"
    if producers:
        print(f"{tab} {symbol_direction} {producers}")
        if levels_max == 0:  # Add producers in the dependency dictionary.
            for p in producers:
                if p not in dic_wf_tree:
                    dic_wf_tree[p] = dic_wf_all[p]
        else:  # Search for more dependencies if needed.
            print("\nWorkflow dependency tree:\n")
            for p in producers:
                get_tree_for_workflow(p, dic_wf_all, dic_wf_tree, case_sensitive, 0, levels_max, reverse)
    else:
        print(f'No {"consumers" if reverse else "producers"} found')
    return dic_wf_tree


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Find dependencies required to produce a given table or to run a given workflow."
    )
    parser.add_argument(
        "-t", dest="table", type=str, nargs="+", help="table(s) for normal (backward) search (i.e. find producers)"
    )
    parser.add_argument(
        "-w", dest="workflow", type=str, nargs="+", help="workflow(s) for normal (backward) search (i.e. find inputs)"
    )
    parser.add_argument(
        "-T", dest="table_rev", type=str, nargs="+", help="table(s) for reverse (forward) search (i.e. find consumers)"
    )
    parser.add_argument(
        "-W",
        dest="workflow_rev",
        type=str,
        nargs="+",
        help="workflow(s) for reverse (forward) search (i.e. find outputs)",
    )
    parser.add_argument(
        "-c",
        dest="case",
        action="store_true",
        help="be case-sensitive with table names",
    )
    parser.add_argument(
        "-g",
        dest="suffix",
        type=str,
        choices=["pdf", "svg", "png"],
        help="make a topology graph in a given format",
    )
    parser.add_argument(
        "-x",
        dest="exclude",
        type=str,
        nargs="+",
        help="tables and workflows to exclude",
    )
    parser.add_argument(
        "-l",
        dest="levels",
        type=int,
        default=0,
        help="maximum number of workflow tree levels (default = 0, include all if < 0)",
    )
    args = parser.parse_args()
    if not (args.table or args.workflow or args.table_rev or args.workflow_rev):
        parser.error("Provide table(s) and/or workflow(s)")
    tables = args.table
    workflows = args.workflow
    tables_rev = args.table_rev
    workflows_rev = args.workflow_rev
    case_sensitive = args.case
    graph_suffix = args.suffix
    list_exclude = args.exclude
    n_levels = args.levels

    # Load all workflows from JSON files
    dic_wf_all_full = load_workflows_from_json()
    # Extract only needed info and make a simplified dictionary
    dic_wf_all_simple = {}
    for wf, dic_wf in dic_wf_all_full.items():
        # Skip excluded workflows
        if list_exclude and wf in list_exclude:
            continue
        dic_wf_all_simple[wf] = {}
        list_dev = get_devices(dic_wf)
        for dev in list_dev:
            dic_wf_all_simple[wf][dev] = {}
            list_inputs = get_inputs(dic_wf, dev)
            list_outputs = get_outputs(dic_wf, dev)
            # Skip excluded tables
            if list_exclude:
                list_inputs = [i for i in list_inputs if i not in list_exclude]
                list_outputs = [o for o in list_outputs if o not in list_exclude]
            dic_wf_all_simple[wf][dev]["inputs"] = list_inputs
            dic_wf_all_simple[wf][dev]["outputs"] = list_outputs
    # print_workflows(dic_wf_all_simple)
    # return

    # Dictionary with dependencies
    dic_deps = {}

    # Find table dependencies
    for t, reverse in zip((tables, tables_rev), (False, True)):
        if t:
            for table in t:
                print(f"\nTable: {table}\n")
                if not table:
                    msg_fatal("Bad table")
                # producers = get_table_producers(table, dic_wf_all_simple, case_sensitive)
                # if not producers:
                #     print("No producers found")
                #     return
                # print(producers)
                # print_workflows(dic_wf_all_simple, producers)
                get_tree_for_table(table, dic_wf_all_simple, dic_deps, case_sensitive, n_levels, reverse)

    # Find workflow dependencies
    for w, reverse in zip((workflows, workflows_rev), (False, True)):
        if w:
            for workflow in w:
                print(f"\nWorkflow: {workflow}\n")
                if not workflow:
                    msg_fatal("Bad workflow")
                # print_workflows(dic_wf_all_simple, [workflow])
                get_tree_for_workflow(workflow, dic_wf_all_simple, dic_deps, case_sensitive, 0, n_levels, reverse)

    # Print the tree dictionary with dependencies
    # print("\nTree\n")
    # print(dic_deps)

    # Produce topology graph.
    if graph_suffix and dic_deps:
        names_all = []
        for names in (tables, tables_rev, workflows, workflows_rev):
            if names:
                names_all += names
        names_all = list(dict.fromkeys(names_all))  # Remove duplicities
        basename = "_".join(names_all)
        # Set a short file name when the full name would be longer than 255 characters.
        if len(basename) > 251:
            basename = "o2_dependencies_" + hashlib.sha1(basename.encode(), usedforsecurity=False).hexdigest()
        ext_graph = graph_suffix
        path_file_dot = basename + ".gv"
        path_file_graph = basename + "." + ext_graph
        print(f"\nMaking dot file in: {path_file_dot}")
        dot = "digraph {\n"
        dot += "  node [shape=box, fontname=Courier, fontsize=20]\n"
        dot += "  ranksep=2 // vertical node separation\n"
        # dot += "  edge [dir=back] // inverted arrow direction\n"
        # dot += "  rankdir=BT // bottom to top drawing\n"
        # lines with dependencies
        dot_deps = ""
        # subgraph for tables
        dot_tables = "  subgraph tables {\n"
        dot_tables += "    node [fillcolor=lightgrey,style=filled]\n"
        list_tables = []
        # subgraph for workflows
        dot_workflows = "  subgraph workflows {\n"
        dot_workflows += '    node [fillcolor=papayawhip,style="filled,rounded"]\n'
        for wf in dic_deps:
            # Hyphens are not allowed in node names.
            node_wf = wf.replace("-", "_")
            # Remove the workflow prefix.
            # label_wf = wf.replace("o2-analysis-", "")
            label_wf = wf
            # Replace hyphens with line breaks to save horizontal space.
            # label_wf = label_wf.replace("-", "\\n")
            dot_workflows += f'    {node_wf} [label="{label_wf}"]\n'
            inputs = get_workflow_inputs(wf, dic_deps)
            outputs = get_workflow_outputs(wf, dic_deps)
            list_tables += inputs + outputs
            nodes_in = " ".join(inputs)
            nodes_out = " ".join(outputs)
            dot_deps += f"  {{{nodes_in}}} -> {node_wf} -> {{{nodes_out}}}\n"
        list_tables = list(dict.fromkeys(list_tables))  # Remove duplicities
        for table in list_tables:
            dot_tables += f"    {table}\n"
        dot_tables += "  }\n"
        dot_workflows += "  }\n"
        dot += dot_workflows + dot_tables + dot_deps
        dot += "}\n"
        try:
            with open(path_file_dot, "w", encoding="utf-8") as file_dot:
                file_dot.write(dot)
        except IOError:
            msg_fatal(f"Failed to open file {path_file_dot}")
        cmd = f"dot -T{ext_graph} {path_file_dot} -o {path_file_graph}"
        try:
            print(f"Making graph in: {path_file_graph}")
            sp.run(cmd, shell=True, check=True)  # nosec B602
        except sp.CalledProcessError:
            msg_fatal(f"Failed to execute: {cmd}")


if __name__ == "__main__":
    main()
