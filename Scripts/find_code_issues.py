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
import re
import sys


github_mode = False # GitHub mode

def is_camel_case(name: str) -> bool:
    """forExample or ForExample"""
    return not "_" in name and not "-" in name


def is_upper_camel_case(name: str) -> bool:
    """ForExample"""
    return name[0].isupper() and is_camel_case(name)


def is_lower_camel_case(name: str) -> bool:
    """forExample"""
    return name[0].islower() and is_camel_case(name)


def is_kebab_case(name: str) -> bool:
    """for-example"""
    return name.islower() and not "_" in name


def is_snake_case(name: str) -> bool:
    """for_example"""
    return name.islower() and not "-" in name


def is_screaming_snake_case(name: str) -> bool:
    """FOR_EXAMPLE"""
    return name.isupper() and not "-" in name


def print_error(path: str, line: int, title: str, message: str) -> str:
    """Format and print error message."""
    # return # Use to suppress error message when counting speed.
    str_line = "" if line is None else f"{line}:"
    print(f"{path}:{str_line} {message} [{title}]") # terminal format
    if github_mode:
        str_line = "" if line is None else f",line={line}"
        print(f"::error file={path}{str_line},title=[{title}]::{message}") # GitHub action format


def is_comment_cpp(line: str) -> bool:
    """Test whether a line is a C++ comment."""
    return line.startswith(("//", "/*"))


class TestSpec(ABC):
    """Prototype of a test class"""
    name = "test template" # short name of the test
    message = "Test failed" # error message
    suffixes = [] # suffixes of files to test
    per_line = True # Test lines separately one by one.

    def file_matches(self, path: str) -> bool:
        """Test whether the path matches the pattern for files to test."""
        return path.endswith(tuple(self.suffixes)) if self.suffixes else True

    # @abstractmethod
    def test_line(self, line: str) -> bool:
        """Test a line."""
        raise NotImplementedError()

    # @abstractmethod
    def test_file(self, path : str, content) -> bool:
        """Test a file in a way that cannot be done line by line."""
        raise NotImplementedError()

    def run(self, path : str, content) -> bool:
        """Run the test."""
        # print(content)
        passed = True
        if not self.file_matches(path):
            return passed
        # print(f"Running test {self.name} for {path} with {len(content)} lines")
        if self.per_line:
            for i, line in enumerate(content):
                line = line.strip()
                if not line:
                    continue
                # print(i + 1, line)
                if not self.test_line(line):
                    passed = False
                    print_error(path, i + 1, self.name, self.message)
        else:
            passed = self.test_file(path, content)
            if not passed:
                print_error(path, None, self.name, self.message)
        return passed


##########################
# Implementations of tests
##########################

# Bad practice

class TestIOStream(TestSpec):
    """Detect included iostream."""
    name = "include iostream"
    message = "Including iostream is not allowed. Use O2 logging instead."
    suffixes = [".h", ".cxx"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        return not line.strip().startswith("#include <iostream>")


class TestUsingStd(TestSpec):
    """Detect importing names from the std namespace."""
    name = "import std names"
    message = "Importing names from the std namespace is not allowed."
    suffixes = [".h", ".cxx"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        return not line.strip().startswith("using std::")


class TestUsingDirectives(TestSpec):
    """Detect using directives in headers."""
    name = "using directives"
    message = "Using directives are not allowed in headers."
    suffixes = [".h"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        return not line.strip().startswith("using namespace")


class TestStdPrefix(TestSpec):
    """Detect missing std:: prefix for common names from the std namespace."""
    name = "std prefix"
    message = "Use std:: prefix for names from the std namespace."
    suffixes = [".h", ".cxx", ".C"]
    prefix_bad = "[^\w:.]"
    patterns = ["vector<", "array[<\{\(]", "f?abs\(", "sqrt\(", "pow\(", "min\(", "max\(", "log\(", "exp\(", "sin\(", "cos\(", "tan\(", "atan\(", "atan2\("]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        for pattern in self.patterns:
            if re.search(f"{self.prefix_bad}{pattern}", line):
                return False
            # occurrences = re.findall(pattern, line)
            # if occurrences:
            #     print(occurrences)
            #     return False
        return True


class TestROOT(TestSpec):
    """Detect use of unnecessary ROOT entities."""
    name = "ROOT entities"
    message = "Consider replacing ROOT entities with STD C++ or O2 entities."
    suffixes = [".h", ".cxx"]
    keywords = ["TMath", "Double_t", "Float_t", "Int_t", "Bool_t"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and not "Macros/" in path

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        for k in self.keywords:
            if k in line:
                return False
        return True


class TestPI(TestSpec):
    """Detect use of external PI."""
    name = "external PI"
    message = "Consider using the PI constant (and its multiples) defined in o2::constants::math."
    suffixes = [".h", ".cxx"]
    keywords = ["M_PI", "TMath::Pi", "TMath::TwoPi"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and not "Macros/" in path

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        for k in self.keywords:
            if k in line:
                return False
        return True


class TestLogging(TestSpec):
    """Detect non-O2 logging."""
    name = "logging"
    message = "Consider using O2 logging (LOG, LOGF, LOGP)."
    suffixes = [".h", ".cxx"]
    keywords = ["Printf(", "printf(", "cout <", "cin >"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and not "Macros/" in path

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        for k in self.keywords:
            if k in line:
                return False
        return True


class TestConstRefInForLoop(TestSpec):
    """Test const refs in range-based for loops."""
    name = "const ref in for loop"
    message = "Use constant references for non-modified iterators in range-based for loops."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not line.startswith("for (") or " : " not in line:
            return True
        line = line[:line.index(" : ")] # keep only the iterator part
        return True if re.search("\([\w]* ?const ?[\w]*&", line) else False


class TestConstRefInSubscription(TestSpec):
    """Test const refs in process function subscriptions."""
    name = "const ref in process"
    message = "Use constant references for table subscriptions in process functions."
    suffixes = [".cxx"]
    # suffixes = [".h"]
    per_line = False

    def test_file(self, path : str, content) -> bool:
        passed = True
        n_parens_opened = 0 # number of opened parentheses
        arguments = "" # process function arguments
        line_process = 0
        for i, line in enumerate(content):
            line = line.strip()
            if is_comment_cpp(line):
                continue
            if re.search("^void process[\w]*\(", line):
                line_process = (i + 1)
                i_closing = line.rfind(")")
                i_start = line.find("(") + 1
                i_end = i_closing if i_closing != -1 else len(line)
                arguments = line[i_start:i_end] # get arguments between parentheses
                n_parens_opened = line.count("(") - line.count(")")
            elif n_parens_opened > 0:
                i_closing = line.rfind(")")
                i_start = 0
                i_end = i_closing if i_closing != -1 else len(line)
                arguments += " " + line[i_start:i_end] # get arguments between parentheses
                n_parens_opened += line.count("(") - line.count(")")
            if line_process > 0 and n_parens_opened == 0:
                # process arguments
                # sanitise template arguments
                template_args = re.findall("<[\w:, ]*>", arguments)
                if template_args:
                    for arg in template_args:
                        if ", " in arg:
                            arguments = arguments.replace(arg, arg.replace(", ", ":"))
                words = arguments.split(", ")
                # test
                for arg in words:
                    if not re.search("[\w<>:]* ?const ?[\w<>:]*&", arg):
                        passed = False
                        print_error(path, i + 1, self.name, f"Argument {arg} is not const&.")
                line_process = 0
        return passed


# Naming conventions
# Reference: https://rawgit.com/AliceO2Group/CodingGuidelines/master/naming_formatting.html


class TestNameFunctionVariable(TestSpec):
    """Test names of functions and of most variables.
    Might report false positives.
    Does not detect multiple variable declarations, i.e. "type name1, name2;"
    Does not detect function arguments on the same line as the function declaration.
    Does not check capitalisation for constexpr because of special rules for constants. See TestNameConstant.
    """
    name = "function/variable names"
    message = "Use lowerCamelCase for names of functions and variables."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        # Look for declarations of functions and variables.

        # Strip away irrelevant remainders of the line after the object name.
        # For functions, stripping after "(" is enough but this way we also identify many declarations of variables.
        for keyword in ("(", "{", ";", " = ", "//", "/*"):
            if keyword in line:
                line = line[:line.find(keyword)]

        # Check the words.
        words = line.split()

        # number of words
        if len(words) < 2:
            return True

        # First word starts with a letter.
        if not words[0][0].isalpha():
            return True

        # Reject false positives with same structure.
        if words[0] in ("return", "if", "else", "new", "delete", "delete[]", "case", "typename", "using", "typedef", "enum", "namespace", "struct", "class"):
            return True
        if len(words) > 2 and words[1] in ("typename", "class", "struct"):
            return True

        # Identify the position of the name for cases "name[n + m]".
        funval_name = words[-1] # expecting the name in the last word
        if funval_name.endswith("]") and "[" not in funval_name: # it's an array and we do not have the name before "[" here
            opens_brackets = ["[" in w for w in words]
            index_name = opens_brackets.index(True) # the name is in the first element with "["
            funval_name = words[index_name]
            words = words[:(index_name + 1)] # Strip away words after the name.
            if len(words) < 2: # Check the adjusted number of words.
                return True

        # All words before the name start with an alphanumeric character (underscores not allowed).
        # Rejects expressions, e.g. * = += << }, but accepts numbers in array declarations.
        if not all(w[0].isalnum() for w in words[:-1]):
            return True

        # Extract function/variable name.
        if "[" in funval_name: # Remove brackets for arrays.
            funval_name = funval_name[:funval_name.find("[")]
        if "::" in funval_name: # Remove the class prefix for methods.
            funval_name = funval_name.split("::")[-1]

        # Check the name candidate.

        # Names of variables and functions are identifiers.
        if not funval_name.isidentifier(): # should be same as ^[\w]+$
            return True

        # print(f"{line} -> {funval_name}")
        # return True
        # The actual test comes here.
        if "constexpr" in words[:-1]:
            return is_camel_case(funval_name)
        return is_lower_camel_case(funval_name)


class TestNameMacro(TestSpec):
    """Test macro names."""
    name = "macro names"
    message = "Use SCREAMING_SNAKE_CASE for macro names."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not line.startswith("#define "):
            return True
        # Extract macro name.
        macro_name = line.split()[1]
        if "(" in macro_name:
            macro_name = macro_name.split("(")[0]
        # The actual test comes here.
        return is_screaming_snake_case(macro_name)


class TestNameConstant(TestSpec):
    """Test constexpr constant names."""
    name = "constexpr constant names"
    message = "Use UpperCamelCase for constexpr constant names. Names of special constants may be prefixed with \"k\"."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        words = line.split()
        if not "constexpr" in words or not "=" in words:
            return True
        # Extract constant name.
        words = words[:words.index("=")] # keep only words before "="
        constant_name = words[-1] # last word before "="
        if constant_name.endswith("]") and "[" not in constant_name: # it's an array and we do not have the name before "[" here
            opens_brackets = ["[" in w for w in words]
            constant_name = words[opens_brackets.index(True)] # the name is in the first element with "["
        # The actual test comes here.
        if constant_name.startswith("k"): # exception for special constants
            constant_name = constant_name[1:] # test the name without "k"
        return is_upper_camel_case(constant_name)


class TestNameNamespace(TestSpec):
    """Test names of namespaces."""
    name = "namespace names"
    message = "Use snake_case for names of namespaces."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not line.startswith("namespace "):
            return True
        # Extract namespace name.
        namespace_name = line.split()[1]
        # The actual test comes here.
        return is_snake_case(namespace_name)


class TestNameUpperCamelCase(TestSpec):
    """Base class for a test of UpperCamelCase names."""
    keyword = "key"
    name = f"{keyword} UpperCamelCase"
    message = f"Use UpperCamelCase for names of {keyword}."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not line.startswith(f"{self.keyword} "):
            return True
        # Extract object name.
        words = line.split()
        if not words[1].isalnum(): # "struct : ...", "enum { ..."
            return True
        object_name = words[1]
        if object_name == "class" and len(words) > 2: # enum class ...
            object_name = words[2]
        # The actual test comes here.
        return is_upper_camel_case(object_name)


class TestNameEnum(TestNameUpperCamelCase):
    """Test names of enumerators."""
    keyword = "enum"
    name = f"{keyword} names"
    message = f"Use UpperCamelCase for names of enumerators and their values."


class TestNameClass(TestNameUpperCamelCase):
    """Test names of classes."""
    keyword = "class"
    name = f"{keyword} names"
    message = f"Use UpperCamelCase for names of classes."


class TestNameStruct(TestNameUpperCamelCase):
    """Test names of structs."""
    keyword = "struct"
    name = f"{keyword} names"
    message = f"Use UpperCamelCase for names of structs."


class TestNameFileCpp(TestSpec):
    """Test names of C++ files."""
    name = "C++ file names"
    message = "Use lowerCamelCase or UpperCamelCase for names of C++ files. See the O2 naming conventions for details."
    suffixes = [".h", ".cxx", ".C"]
    per_line = False

    def test_file(self, path : str, content) -> bool:
        file_name = os.path.basename(path)
        return is_camel_case(file_name)


class TestNameFilePython(TestSpec):
    """Test names of Python files."""
    name = "Python file names"
    message = "Use snake_case for names of Python files."
    suffixes = [".py", ".ipynb"]
    per_line = False

    def test_file(self, path : str, content) -> bool:
        file_name = os.path.basename(path)
        return is_snake_case(file_name)


class TestNameWorkflow(TestSpec):
    """Test names of O2 workflows."""
    name = "O2 workflow names"
    message = "Use kebab-case for names of workflows and match the name of the workflow file."
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
            if not is_kebab_case(workflow_name):
                passed = False
                print_error(path, i + 1, self.name, f"Invalid workflow name: {workflow_name}.")
                continue
            # Extract workflow file name.
            next_line = content[i + 1].strip()
            words = next_line.split()
            if words[0] != "SOURCES":
                passed = False
                print_error(path, i + 2, self.name, f"Did not find sources for workflow: {workflow_name}.")
                continue
            workflow_file_name = os.path.basename(words[1]) # the actual file name
            # Generate the file name matching the workflow name.
            expected_workflow_file_name = "".join([w.title() if w[0].isnumeric() else w.capitalize() for w in workflow_name.split("-")]) + ".cxx"
            expected_workflow_file_name = f"{expected_workflow_file_name[0].lower()}{expected_workflow_file_name[1:]}" # start with lowercase letter
            # Compare the actual and expected file names.
            if expected_workflow_file_name != workflow_file_name:
                passed = False
                print_error(path, i + 1, self.name, f"Workflow name {workflow_name} does not match the workflow file name {workflow_file_name} (expected {expected_workflow_file_name}).")
        return passed


# PWG-specific


class TestHfConstAuto(TestSpec):
    """PWGHF: Detect swapped const auto."""
    name = "PWGHF: const auto"
    message = "Use \"const auto\" instead of \"auto const\"."
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and "PWGHF/" in path

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        return not "auto const" in line


class TestHfNameStructClass(TestSpec):
    """PWGHF: Test names of structs and classes."""
    name = "PWGHF: struct/class names"
    message = "Names of PWGHF structs and classes must start with \"Hf\"."
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and "PWGHF/" in path

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not line.startswith(("struct ", "class ")):
            return True
        # Extract struct/class name.
        words = line.split()
        if not words[1].isalnum(): # "struct : ..."
            return True
        struct_name = words[1]
        # The actual test comes here.
        return struct_name.startswith("Hf")


class TestHfStructMembers(TestSpec):
    """PWGHF: Test order of struct members.
    Caveat: Does not see Configurables in ConfigurableGroup."""
    name = "PWGHF: struct member order"
    message = "Declare struct members in the conventional order. See the PWGHF coding guidelines."
    suffixes = [".cxx"]
    per_line = False
    member_order = ["Spawns<", "Builds<", "Produces<", "Configurable<", "HfHelper ", "SliceCache ", "Service<", "using ", "Filter ", "Preslice<", "Partition<", "ConfigurableAxis ", "AxisSpec ", "HistogramRegistry ", "OutputObj<", "void init(", "void process"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and "PWGHF/" in path

    def test_file(self, path : str, content) -> bool:
        passed = True
        dic_struct = {}
        struct_name = ""
        for i, line in enumerate(content):
            if line.strip().startswith("//"):
                continue
            if line.startswith("struct "): # expecting no indentation
                struct_name = line.strip().split()[1]
                dic_struct[struct_name] = {}
                continue
            if not struct_name:
                continue
            # Save line numbers of members of the current struct for each category.
            for member in self.member_order:
                if line.startswith(f"  {member}"): # expecting single-level indentation for direct members
                    if member not in dic_struct[struct_name]:
                        dic_struct[struct_name][member] = []
                    dic_struct[struct_name][member].append(i + 1) # save line number
                    break
        # print(dic_struct)
        # Detect members declared in a wrong order.
        last_line_last_member = 0 # line number of the last member of the previous member category
        index_last_member = 0 # index of the previous member category in the member_order list
        for struct_name in dic_struct:
            for i_m, member in enumerate(self.member_order):
                if member not in dic_struct[struct_name]:
                    continue
                first_line = min(dic_struct[struct_name][member]) # line number of the first member of this category
                last_line = max(dic_struct[struct_name][member]) # line number of the last member of this category
                if first_line < last_line_last_member: # The current category starts before the end of the previous category.
                    passed = False
                    print_error(path, first_line, self.name, f"{struct_name}: {member.strip()} appears too early (before end of {self.member_order[index_last_member].strip()}).")
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
        struct_names = [] # actual struct names in the file
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
        for struct_name in struct_names:
            if struct_name.startswith(base_struct_name):
                return True
        return False


class TestHfNameConfigurable(TestSpec):
    """PWGHF: Test names of configurables."""
    name = "PWGHF: Configurable names"
    message = "Use camelCase for Configurable names and use the same name for the struct member as for the JSON string."
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and "PWGHF/" in path and not "Macros/" in path

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not line.startswith("Configurable"):
            return True
        # Extract Configurable name.
        words = line.split()
        if words[2] == "=": # expecting Configurable... nameCpp = {"nameJson",
            name_cpp = words[1] # nameCpp
            name_json = words[3][1:] # expecting "nameJson",
        else:
            names = words[1].split("{") # expecting Configurable... nameCpp{"nameJson",
            name_cpp = names[0] # nameCpp
            name_json = names[1] # expecting "nameJson",
        if name_json[0] != "\"": # JSON name is not a literal string.
            return True
        name_json = name_json.strip("\",") # expecting nameJson
        # The actual test comes here.
        return is_lower_camel_case(name_cpp) and name_cpp == name_json


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
        help="Print messages also as GitHub annotations",
    )
    args = parser.parse_args()
    if args.github:
        global github_mode
        github_mode = True

    tests = [] # list of activated tests

    # Bad practice
    enable_bad_practice = True
    if enable_bad_practice:
        # tests.append(TestIOStream())
        # tests.append(TestUsingStd())
        # tests.append(TestUsingDirectives())
        tests.append(TestStdPrefix())
        # tests.append(TestROOT())
        # tests.append(TestPI())
        # tests.append(TestLogging())
        # tests.append(TestConstRefInForLoop())
        # tests.append(TestConstRefInSubscription())

    # Naming conventions
    enable_naming = False
    if enable_naming:
        tests.append(TestNameFunctionVariable())
        tests.append(TestNameMacro())
        tests.append(TestNameConstant())
        tests.append(TestNameNamespace())
        tests.append(TestNameEnum())
        tests.append(TestNameClass())
        tests.append(TestNameStruct())
        tests.append(TestNameFileCpp())
        tests.append(TestNameFilePython())
        tests.append(TestNameWorkflow())

    # PWGHF
    enable_pwghf = False
    if enable_pwghf:
        tests.append(TestHfConstAuto())
        tests.append(TestHfNameStructClass())
        tests.append(TestHfStructMembers())
        tests.append(TestHfNameFileWorkflow())
        tests.append(TestHfNameConfigurable())

    test_names = [t.name for t in tests] # short names of activated tests
    suffixes = tuple(set([s for test in tests for s in test.suffixes])) # all suffixes from all enabled tests
    passed = True # global result of all tests

    # Report overview before running.
    print(f"Testing {len(args.paths)} files.")
    print("Enabled tests:", test_names)
    print("Suffixes of tested files:", suffixes)

    # Test files.
    for path in args.paths:
        # Skip not tested files.
        if not path.endswith(suffixes):
            continue
        try:
            with open(path, "r") as file:
                content = file.readlines()
                for test in tests:
                    if not test.run(path, content):
                        passed = False
                    # if not passed:
                    #     print(f"File {path} failed the test {test.name}.")
                    # print(f"Test {'passed' if passed else 'failed'}.")
        except IOError:
            print(f"Failed to open file {path}.")
            sys.exit(1)

    # Report result.
    if passed:
        print("All tests passed.")
    else:
        print("Issues have been found.")
        sys.exit(1)


if __name__ == "__main__":
    main()
