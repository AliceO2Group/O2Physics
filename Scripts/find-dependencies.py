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

"""
Find dependencies required to produce a given table or to run a given workflow.
Author: Vít Kučera on 2022-12-10
"""

import json
import argparse
import sys
import os
import glob
import subprocess as sp  # nosec B404


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


def join_strings(obj):
    """Return strings concatenated into one."""
    if isinstance(obj, str):
        return obj
    elif isinstance(obj, list):
        return " ".join(obj)
    else:
        msg_fatal("Cannot convert %s into a string" % type(obj))


def load_workflows():
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
        workflow = os.path.basename(file_json).split('.')[0]
        try:
            with open(file_json, 'r') as j:
                specs_wf = json.load(j)
                db_wf[workflow] = specs_wf["workflow"]
        except FileNotFoundError:
            msg_fatal("JSON file not found.")
        # print(workflow)
    return db_wf


def print_workflows(dic_wf_all : dict, list_wf=None):
    """Print properties of a given workflow or workflows in the simplified dictionary."""
    def print_wf(dic_wf : dict, wf : str):
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


def get_devices(specs_wf : dict):
    """Get the list of devices of a given workflow loaded from a JSON files"""
    return [d["name"] for d in specs_wf]


def get_inputs(specs_wf : dict, device=""):
    """Get the list of input tables of a given workflow loaded from a JSON files.
    If a device names is provided, only inputs of that device are considered."""
    l = []
    for dev in specs_wf:
        if device and dev["name"] != device:
            continue
        l += [i['binding'] for i in dev["inputs"] if i["origin"] == "AOD"]
    return list(dict.fromkeys(l)) # Remove duplicities


def get_outputs(specs_wf : dict, device=""):
    """Get the list of output tables of a given workflow loaded from a JSON files.
    If a device names is provided, only outputs of that device are considered."""
    l = []
    # Loop over devices
    for dev in specs_wf:
        if device and dev["name"] != device:
            continue
        l += [i['binding'] for i in dev["outputs"] if i["origin"] == "AOD"]
    return list(dict.fromkeys(l)) # Remove duplicities


def get_table_producers(table : str, dic_wf_all : dict, case_sensitive=False):
    """Find all workflows that have this table as output."""
    l = []
    if not case_sensitive:
        table = table.lower()
    # Loop over workflows
    for wf, dic_wf in dic_wf_all.items():
        # Loop over devices
        for dev in dic_wf:
            outputs = [o if case_sensitive else o.lower() for o in dic_wf[dev]["outputs"]]
            if table in outputs:
                l.append(wf)
    return l


def get_workflow_outputs(wf : str, dic_wf_all : dict):
    """Get list of workflow outputs from the simplified dictionary"""
    l = []
    # Loop over devices
    for dev in dic_wf_all[wf]:
        l += dic_wf_all[wf][dev]["outputs"]
    return list(dict.fromkeys(l)) # Remove duplicities


def get_workflow_inputs(wf : str, dic_wf_all : dict):
    """Get list of workflow inputs from the simplified dictionary"""
    l = []
    # list_outputs = get_workflow_outputs(wf, dic_wf_all)
    # Loop over devices
    for dev in dic_wf_all[wf]:
        l += dic_wf_all[wf][dev]["inputs"]
    l = list(dict.fromkeys(l)) # Remove duplicities
    # l = [d for d in l if d not in list_outputs] # avoid circular dependence (can be legit)
    return l


def get_tree_for_workflow(wf, dic_wf_all, dic_wf_tree=None, case_sensitive=False, level=0, levels_max=0):
    """Get the dependency tree of workflows needed to run this workflow"""
    # print(level, levels_max)
    if dic_wf_tree is None:
        dic_wf_tree = {}
    if wf not in dic_wf_all:
        msg_fatal(f"Workflow {wf} not found")
    if wf not in dic_wf_tree:
        dic_wf_tree[wf] = dic_wf_all[wf]
    inputs = get_workflow_inputs(wf, dic_wf_all)
    if inputs:
        print(f"{level * '    '}{wf} <- {inputs}")
        if levels_max <= -1 or level < levels_max:
            for tab in inputs:
                producers = get_table_producers(tab, dic_wf_all, case_sensitive)
                if producers:
                    print(f"{(level) * '    ' + '  '}{tab} <- {producers}")
                    for p in producers:
                        if p not in dic_wf_tree: # avoid infinite recursion
                            get_tree_for_workflow(p, dic_wf_all, dic_wf_tree, case_sensitive, level + 1, levels_max)
    return dic_wf_tree


def get_tree_for_table(tab, dic_wf_all, dic_wf_tree=None, case_sensitive=False, levels_max=0):
    """Get the dependency tree of workflows needed to produce this table"""
    if dic_wf_tree is None:
        dic_wf_tree = {}
    producers = get_table_producers(tab, dic_wf_all, case_sensitive)
    if producers:
        print(f"{tab} <- {producers}\n")
        print("Workflow dependency tree:\n")
        for p in producers:
            get_tree_for_workflow(p, dic_wf_all, dic_wf_tree, case_sensitive, 0, levels_max)
    else:
        print("No producers found")
    return dic_wf_tree


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Find dependencies required to produce a given table or to run a given workflow."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-t", dest="table", help="table"
    )
    group.add_argument(
        "-w", dest="workflow", help="workflow"
    )
    parser.add_argument(
        "-c", dest="case", action="store_true", help="be case-sensitive with table names"
    )
    parser.add_argument(
        "-g", dest="suffix", type=str, choices=["pdf", "svg", "png"], help="make a topology graph if suffix provided"
    )
    parser.add_argument(
        "-x", dest="exclude", type=str, nargs='+', help="tables and workflows to exclude"
    )
    parser.add_argument(
        "-l", dest="levels", type=int, default=0, help="maximum number of workflow tree levels"
    )
    args = parser.parse_args()
    table = args.table
    workflow = args.workflow
    case_sensitive = args.case
    graph_suffix = args.suffix
    list_exclude = args.exclude
    n_levels = args.levels

    # Load all workflows from JSON files
    dic_wf_all_full = load_workflows()
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

    # Find table producers
    if table:
        print(f"Table: {table}\n")
        # producers = get_table_producers(table, dic_wf_all_simple, case_sensitive)
        # if not producers:
        #     print("No producers found")
        #     return
        # print(producers)
        # print_workflows(dic_wf_all_simple, producers)
        dic_deps = get_tree_for_table(table, dic_wf_all_simple, None, case_sensitive, n_levels)

    # Find workflow dependencies
    if workflow:
        print(f"Workflow: {workflow}\n")
        # if workflow not in dic_wf_all_simple:
        #     msg_fatal(f"Workflow {workflow} not found")
        # print_workflows(dic_wf_all_simple, [workflow])
        dic_deps = get_tree_for_workflow(workflow, dic_wf_all_simple, None, case_sensitive, 0, n_levels)

    # Print the tree dictionary with dependencies
    # print("\nTree")
    # print(dic_deps)

    # Produce topology graph.
    if graph_suffix and dic_deps:
        basename = workflow if workflow else table
        ext_graph = graph_suffix
        path_file_dot = basename + ".gv"
        path_file_graph = basename + "." + ext_graph
        print(f"\nMaking dot file in: {path_file_dot}")
        dot = "digraph {\n"
        dot += "  ranksep=2 // vertical node separation\n"
        dot += '  node [shape=box, fontname=Courier, fontsize=20]\n'
        # dot += "  edge [dir=back] // inverted arrow direction\n"
        # dot += "  rankdir=BT // bottom to top drawing\n"
        # lines with dependencies
        dot_deps = ""
        # subgraph for tables
        dot_tables = '  subgraph tables {\n'
        # dot_tables += '    label="tables"\n'
        dot_tables += '    node [fillcolor=lightgrey,style=filled]\n'
        list_tables = []
        # subgraph for workflows
        dot_workflows = '  subgraph workflows {\n'
        # dot_wf += '    label="workflows"\n'
        dot_workflows += '    node [fillcolor=papayawhip,style="filled,rounded"]\n'
        for wf in dic_deps:
            # Hyphens are not allowed in node names.
            node_wf = wf.replace("-", "_")
            # Remove the workflow prefix.
            # label_wf = wf.replace("o2-analysis-", "")
            label_wf = wf
            # Replace hyphens with line breaks to save horizontal space.
            # label_wf = label_wf.replace("-", "\\n")
            dot_workflows += '    %s [label="%s"]\n' % (node_wf, label_wf)
            inputs = get_workflow_inputs(wf, dic_deps)
            outputs = get_workflow_outputs(wf, dic_deps)
            list_tables += inputs + outputs
            nodes_in = join_strings(inputs).replace("-", "_")
            nodes_out = join_strings(outputs).replace("-", "_")
            dot_deps += f"  {{{nodes_in}}} -> {node_wf} -> {{{nodes_out}}}\n"
        list_tables = list(dict.fromkeys(list_tables)) # Remove duplicities
        for table in list_tables:
            dot_tables += f"    {table}\n"
        dot_tables += "  }\n"
        dot_workflows += "  }\n"
        dot += dot_workflows + dot_tables + dot_deps
        dot += "}\n"
        try:
            with open(path_file_dot, "w") as file_dot:
                file_dot.write(dot)
        except IOError:
            msg_fatal("Failed to open file " + path_file_dot)
        cmd = f"dot -T{ext_graph} {path_file_dot} -o {path_file_graph}"
        try:
            print(f"Making graph in: {path_file_graph}")
            sp.run(cmd, shell=True, check=True)  # nosec B602
        except sp.CalledProcessError:
            msg_fatal(f"executing: {cmd}")


main()