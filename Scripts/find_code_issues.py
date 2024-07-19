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
@brief  Find issues in code
@author Vít Kučera <vit.kucera@cern.ch>, Inha University
@date   2024-07-14
"""

from abc import ABC #, abstractmethod
import argparse
import os
import sys


class TestSpec(ABC):
    """Prototype of the test class"""
    name = "Test template" # name of the test
    message = "Test failed" # error message
    suffixes = [] # suffixes of files to test
    per_line = True # test lines separately one by one

    def file_matches(self, path: str) -> bool:
        """Test whether the path matches the pattern for files to test"""
        return path.endswith(tuple(self.suffixes)) if self.suffixes else True

    # @abstractmethod
    def test_line(self, line: str) -> bool:
        """Test a line."""
        raise NotImplementedError()

    # @abstractmethod
    def test_file(self, path : str, content) -> bool:
        """Test a file in a way that cannot be done line by line."""
        raise NotImplementedError()

    def run(self, path : str, content, github: bool = False) -> bool:
        """Test file content."""
        # print(content)
        passed = True
        if not self.file_matches(path):
            return passed
        # print(f"Running test {self.name} for {path} with {len(content)} lines")
        if self.per_line:
            for i, line in enumerate(content):
                line = line.strip()
                # print(i + 1, line)
                if not self.test_line(line):
                    passed = False
                    print(f"{path}:{i + 1}: {self.message} [{self.name}]")
                    if github:
                        print(f"::error file={path},line={i + 1},title=[{self.name}]::{self.message}")
        else:
            passed = self.test_file(path, content)
            if not passed:
                print(f"{path}: {self.message} [{self.name}]")
                if github:
                    print(f"::error file={path},title=[{self.name}]::{self.message}")
        return passed


##########################
# Implementations of tests
##########################

# Bad code

class TestIOStream(TestSpec):
    name = "Include iostream"
    message = "Including iostream is not allowed. Use O2 logging instead."
    suffixes = [".h", ".cxx"]

    def test_line(self, line: str) -> bool:
        if line.startswith("//"):
            return True
        return not line.strip().startswith("#include <iostream>")


class TestUsingStd(TestSpec):
    name = "Importing std names"
    message = "Importing std names is not allowed."
    suffixes = [".h", ".cxx"]

    def test_line(self, line: str) -> bool:
        if line.startswith("//"):
            return True
        return not line.strip().startswith("using std::")


class TestUsingDirectives(TestSpec):
    name = "Using directives"
    message = "Using directives are not allowed in headers."
    suffixes = [".h"]

    def test_line(self, line: str) -> bool:
        if line.startswith("//"):
            return True
        return not line.strip().startswith("using namespace")


class TestStdPrefix(TestSpec):
    name = "std prefix"
    message = "Use std:: prefix for STD methods."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if line.startswith("//"):
            return True
        return not " abs(" in line


class TestROOT(TestSpec):
    """Test use of ROOT entities"""
    name = "ROOT entities"
    message = "Consider replacing ROOT entities with STD C++ or O2 entities."
    suffixes = [".h", ".cxx"]
    keywords = ["TMath", "Double_t", "Float_t", "Int_t", "Bool_t"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and not "Macros/" in path

    def test_line(self, line: str) -> bool:
        if line.startswith("//"):
            return True
        for k in self.keywords:
            if k in line:
                return False
        return True


class TestPI(TestSpec):
    """Test use of PI"""
    name = "Use of PI"
    message = "Consider using the PI constant (and its multiples) defined in o2::constants::math."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if line.startswith("//"):
            return True
        return not "M_PI" in line


class TestLogging(TestSpec):
    """Test logging"""
    name = "Use of logging"
    message = "Consider using O2 logging with LOG, LOGF, LOGP."
    suffixes = [".h", ".cxx"]
    keywords = ["Printf(", "printf(", "cout <", "cin >"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and not "Macros/" in path

    def test_line(self, line: str) -> bool:
        if line.startswith("//"):
            return True
        for k in self.keywords:
            if k in line:
                return False
        return True


# Naming conventions


class TestNameFunction(TestSpec):
    """Test function names.
    Likely to report false positives.
    Can accidentally spot names of variable too but it is fine because same conventions apply.
    """
    name = "function names"
    message = "Use camelCase for names of functions and variables."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        line = line.strip()
        if line.startswith("//"):
            return True
        # Look for "Type Name(..."
        # regexp ?
        words = line.split("(")
        if len(words) < 2:
            return True
        words = words[0].split()
        if len(words) != 2:
            return True
        if not words[1].isalnum():
            return True
        if words[0] in ("return", ":", "#define", "#if", "new", "virtual", "case"):
            return True
        if words[0].endswith(".template"):
            return True
        if words[0][-1] == ",": # multiple template arguments
            return True
        # Extract function name.
        function_name = words[1]
        if "::" in function_name:
            function_name = function_name.split("::")[1]
        # The actual test comes here.
        # if not (function_name[0].islower() and not "_" in function_name and not "-" in function_name):
        #     print(line, words)
        #     return False
        # return True
        return function_name[0].islower() and not "_" in function_name and not "-" in function_name


class TestNameMacro(TestSpec):
    """Test macro names"""
    name = "macro names"
    message = "Use capital letters and underscores in macro names."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        line = line.strip()
        if line.startswith("//"):
            return True
        if not line.startswith("#define "):
            return True
        # Extract macro name.
        macro_name = line.split()[1]
        if "(" in macro_name:
            macro_name = macro_name.split("(")[0]
        # The actual test comes here.
        return macro_name.isupper() and not "-" in macro_name


class TestNameNamespace(TestSpec):
    """Test names of namespaces."""
    name = "namespace names"
    message = "Use lowercase letters and underscores in names of namespaces."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        line = line.strip()
        if line.startswith("//"):
            return True
        if not line.startswith("namespace "):
            return True
        # Extract namespace name.
        namespace_name = line.split()[1]
        # if not namespace_name.islower():
        #     print(namespace_name)
        #     return False
        # return True
        # The actual test comes here.
        return namespace_name.islower() and not "-" in namespace_name


class TestNameUpperCamelCase(TestSpec):
    """Test CamelCase of a name."""
    keyword = "key"
    name = f"{keyword} CamelCase"
    message = f"Use CamelCase for names of {keyword}."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        line = line.strip()
        if line.startswith("//"):
            return True
        if not line.startswith(f"{self.keyword} "):
            return True
        # Extract object name.
        words = line.split()
        if not words[1].isalnum(): # "struct : ...", "enum { ..."
            return True
        object_name = words[1]
        if words[1] == "class" and len(words) > 2: # enum class
            object_name = words[2]
        # if not object_name[0].isupper():
        #     print(object_name)
        #     return False
        # return True
        # The actual test comes here.
        return object_name[0].isupper() and not "_" in object_name and not "-" in object_name


class TestNameEnum(TestNameUpperCamelCase):
    """Test names of enumerators.
    Issue with 'enum class'"""
    keyword = "enum"
    name = f"{keyword} names"
    message = f"Use CamelCase for names of enumerators and their values."


class TestNameClass(TestNameUpperCamelCase):
    """Test names of classes."""
    keyword = "class"
    name = f"{keyword} names"
    message = f"Use CamelCase for names of classes."


class TestNameStruct(TestNameUpperCamelCase):
    """Test names of structs."""
    keyword = "struct"
    name = f"{keyword} names"
    message = f"Use CamelCase for names of structs."


class TestNameFileCpp(TestSpec):
    """Test names of C++ files."""
    name = "C++ file names"
    message = "Use camelCase or CamelCase for names of C++ files. See the O2 naming conventions for details."
    suffixes = [".h", ".cxx", ".C"]
    per_line = False

    def test_file(self, path : str, content) -> bool:
        file_name = os.path.basename(path)
        return not "_" in file_name and not "-" in file_name


class TestNameFilePython(TestSpec):
    """Test names of Python files."""
    name = "Python file names"
    message = "Use lowercase letters and underscores in names of Python files."
    suffixes = [".py", ".ipynb"]
    per_line = False

    def test_file(self, path : str, content) -> bool:
        file_name = os.path.basename(path)
        return file_name.islower() and not "-" in file_name


class TestNameWorkflow(TestSpec):
    """Test names of O2 workflows."""
    name = "O2 workflow names"
    message = "Name of a workflow must consist of hyphenated lowercase words and must match the name of the workflow file."
    suffixes = ["CMakeLists.txt"]
    per_line = False

    def test_file(self, path : str, content) -> bool:
        passed = True
        workflow_name = ""
        for i, line in enumerate(content):
            if not line.startswith("o2physics_add_dpl_workflow("):
                continue
            # Extract workflow name.
            workflow_name = line.strip().split("(")[1]
            if not workflow_name.islower() or "_" in workflow_name:
                passed = False
                print(f"{path}:{i + 1}: Invalid workflow name: {workflow_name}.")
                # if github:
                #     print(f"::error file={path},line={i + 1},title=[{self.name}]::{self.message}")
                continue
            # Extract workflow file name.
            next_line = content[i + 1].strip()
            words = next_line.split()
            if words[0] != "SOURCES":
                passed = False
                print(f"{path}:{i + 2}: Did not find sources for workflow: {workflow_name}.")
                continue
            workflow_file_name = os.path.basename(words[1])
            expected_workflow_file_name = "".join([w.title() if w[0].isnumeric() else w.capitalize() for w in workflow_name.split("-")]) + ".cxx"
            expected_workflow_file_name = f"{expected_workflow_file_name[0].lower()}{expected_workflow_file_name[1:]}"
            if expected_workflow_file_name != workflow_file_name:
                passed = False
                print(f"{path}:{i + 1}: Workflow name {workflow_name} does not match the workflow file name {workflow_file_name} (expected {expected_workflow_file_name}).")
        return passed


# PWG-specific


class TestHfConstAuto(TestSpec):
    name = "PWGHF: const auto"
    message = "Use \"const auto\" instead of \"auto const\""
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and "PWGHF/" in path

    def test_line(self, line: str) -> bool:
        return not "auto const" in line


class TestHfNameStruct(TestSpec):
    name = "PWGHF: struct name"
    message = "Names of PWGHF structs must start with \"Hf\"."
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and "PWGHF/" in path

    def test_line(self, line: str) -> bool:
        line = line.strip()
        if line.startswith("//"):
            return True
        if not line.startswith("struct "):
            return True
        # Extract struct name.
        words = line.split()
        if not words[1].isalnum(): # "struct : ...", "enum { ..."
            return True
        struct_name = words[1]
        if words[1] == "class" and len(words) > 2: # enum class
            struct_name = words[2]
        # The actual test comes here.
        return struct_name.startswith("Hf")


class TestHfStructMembers(TestSpec):
    """Test order of struct members."""
    name = "PWGHF: struct members"
    message = "Order struct members."
    suffixes = [".cxx"]
    per_line = False
    member_order = ["Spawns<", "Builds<", "Produces<", "Configurable<", "HfHelper hfHelper", "Service<", "using ", "Filter ", "Preslice<", "Partition<", "ConfigurableAxis ", "AxisSpec ", "HistogramRegistry ", "OutputObj<", "void init(", "void process"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and "PWGHF/" in path

    def test_file(self, path : str, content) -> bool:
        passed = True
        dic_struct = {}
        struct_name = ""
        for i, line in enumerate(content):
            if line.startswith("//"):
                continue
            if line.startswith("struct "):
                struct_name = line.strip().split()[1]
                dic_struct[struct_name] = {}
                continue
            if not struct_name:
                continue
            for member in self.member_order:
                if line.startswith(f"  {member}"):
                    if member not in dic_struct[struct_name]:
                        dic_struct[struct_name][member] = []
                    dic_struct[struct_name][member].append(i + 1)
                    break
        # print(dic_struct)
        last_line_last_member = 0 # number of the last line of the previous member category
        index_last_member = 0 # index of the previous member category in member_order
        for struct_name in dic_struct:
            for i_m, member in enumerate(self.member_order):
                if member not in dic_struct[struct_name]:
                    continue
                first_line = min(dic_struct[struct_name][member])
                last_line = max(dic_struct[struct_name][member])
                if first_line < last_line_last_member:
                    passed = False
                    print(f"{path}:{first_line}: {struct_name}: {member.strip()} appears too early (before end of {self.member_order[index_last_member].strip()}).")
                last_line_last_member = last_line
                index_last_member = i_m
        return passed


class TestHfNameFileWorkflow(TestSpec):
    """PWGHF: Test names of workflow files."""
    name = "PWGHF: workflow file names"
    message = "Name of a workflow file must match the name of the main task in it. See the PWGHF O2 naming conventions for details."
    suffixes = [".cxx"]
    per_line = False

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and "PWGHF/" in path and not "Macros/" in path

    def test_file(self, path : str, content) -> bool:
        file_name = os.path.basename(path).rstrip(".cxx")
        if file_name[0].isupper():
            return False
        base_struct_name = f"Hf{file_name[0].upper()}{file_name[1:]}" # expected base of struct names
        # print(f"For file {file_name} expecting to find {base_struct_name}.")
        struct_names = []
        for line in content:
            if not line.startswith("struct "):
                continue
            # Extract struct name.
            words = line.split()
            if not words[1].isalnum(): # "struct : ..."
                continue
            struct_name = words[1]
            struct_names.append(struct_name)
        # print(f"Found structs: {struct_names}.")
        for struct in struct_names:
            if struct.startswith(base_struct_name):
                return True
        return False


# End of test implementations


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Run tests on files to find code issues."
    )
    parser.add_argument("paths", type=str, nargs="+", help="File path(s)")
    parser.add_argument(
        "-g",
        dest="github",
        action="store_true",
        help="Print messages as GitHub annotations",
    )
    args = parser.parse_args()

    tests = []

    # Bad code
    enable_bad_practice = True
    enable_bad_practice = False
    if enable_bad_practice:
        tests.append(TestIOStream())
        tests.append(TestUsingDirectives())
        tests.append(TestUsingStd())
        tests.append(TestStdPrefix())
        tests.append(TestROOT())
        tests.append(TestPI())
        tests.append(TestLogging())

    # Naming conventions
    enable_naming = True
    enable_naming = False
    if enable_naming:
        tests.append(TestNameFunction())
        tests.append(TestNameMacro())
        tests.append(TestNameNamespace())
        tests.append(TestNameEnum())
        tests.append(TestNameClass())
        tests.append(TestNameStruct())
        tests.append(TestNameFileCpp())
        tests.append(TestNameFilePython())
        tests.append(TestNameWorkflow())

    # PWGHF
    enable_pwghf = True
    enable_pwghf = False
    if enable_pwghf:
        tests.append(TestHfConstAuto())
        tests.append(TestHfStructMembers())
        tests.append(TestHfNameStruct())
        tests.append(TestHfNameFileWorkflow())

    test_names = [t.name for t in tests]
    suffixes = tuple(set([s for test in tests for s in test.suffixes])) # all suffixes from all enabled tests
    passed = True # global status

    print("Enabled tests:", test_names)

    print(f"Testing {len(args.paths)} files.")
    for path in args.paths:
        # Skip not tested files.
        if not path.endswith(suffixes):
            continue
        try:
            with open(path, "r") as file:
                content = file.readlines()
                for test in tests:
                    if not test.run(path, content, args.github):
                        passed = False
                    # if not passed:
                    #     print(f"File {path} failed the test {test.name}.")
                    # print(f"Test {'passed' if passed else 'failed'}.")
        except IOError:
            print(f"Failed to open file {path}")
            sys.exit(1)

    if passed:
        print("All tests passed.")
    else:
        print("Issues have been found.")
        sys.exit(1)


if __name__ == "__main__":
    main()
