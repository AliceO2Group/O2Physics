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
@brief  O2 linter (Find issues in O2 code)
@author Vít Kučera <vit.kucera@cern.ch>, Inha University
@date   2024-07-14
"""

from abc import ABC
import argparse
import os
import re
import sys
from typing import Union

github_mode = False # GitHub mode

def is_camel_case(name: str) -> bool:
    """forExample or ForExample"""
    return not "_" in name and not "-" in name and not " " in name


def is_upper_camel_case(name: str) -> bool:
    """ForExample"""
    return name[0].isupper() and is_camel_case(name)


def is_lower_camel_case(name: str) -> bool:
    """forExample"""
    return name[0].islower() and is_camel_case(name)


def is_kebab_case(name: str) -> bool:
    """for-example"""
    return name.islower() and not "_" in name and not " " in name


def is_snake_case(name: str) -> bool:
    """for_example"""
    return name.islower() and not "-" in name and not " " in name


def is_screaming_snake_case(name: str) -> bool:
    """FOR_EXAMPLE"""
    return name.isupper() and not "-" in name and not " " in name


def kebab_case_to_camel_case_u(line: str) -> str:
    """Convert kebab-case string to UpperCamelCase string."""
    return "".join([w.title() if w[0].isnumeric() else w.capitalize() for w in line.split("-")])


def kebab_case_to_camel_case_l(line: str) -> str:
    """Convert kebab-case string to lowerCamelCase string."""
    new_line = kebab_case_to_camel_case_u(line)
    return f"{new_line[0].lower()}{new_line[1:]}" # start with lowercase letter


def print_error(path: str, line: Union[int, None], title: str, message: str):
    """Format and print error message."""
    # return # Use to suppress error message when counting speed.
    str_line = "" if line is None else f"{line}:"
    print(f"{path}:{str_line} {message} [{title}]") # terminal format
    if github_mode:
        str_line = "" if line is None else f",line={line}"
        print(f"::error file={path}{str_line},title=[{title}]::{message}") # GitHub annotation format


def is_comment_cpp(line: str) -> bool:
    """Test whether a line is a C++ comment."""
    return line.strip().startswith(("//", "/*"))


def remove_comment_cpp(line: str) -> str:
    """Remove C++ comments from the end of a line."""
    for keyword in ("//", "/*"):
        if keyword in line:
            line = line[:line.index(keyword)]
    return line.strip()


class TestSpec(ABC):
    """Prototype of a test class"""
    name = "test-template" # short name of the test
    message = "Test failed" # error message
    suffixes: "list[str]" = [] # suffixes of files to test
    per_line = True # Test lines separately one by one.

    def file_matches(self, path: str) -> bool:
        """Test whether the path matches the pattern for files to test."""
        return path.endswith(tuple(self.suffixes)) if self.suffixes else True

    def is_disabled(self, line: str, prefix_comment="//") -> bool:
        """Detect whether the test is explicitly disabled."""
        prefix_disable = "o2-linter: disable="
        for prefix in [prefix_comment, prefix_disable]:
            if not prefix in line:
                return False
            line = line[(line.index(prefix) + len(prefix)):] # Strip away part before prefix.
        if self.name in line:
            return True
        return False

    def test_line(self, line: str) -> bool:
        """Test a line."""
        raise NotImplementedError()

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
                if self.is_disabled(line):
                    continue
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
    name = "include-iostream"
    message = "Including iostream is discouraged. Use O2 logging instead."
    suffixes = [".h", ".cxx"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        return not line.startswith("#include <iostream>")


class TestUsingStd(TestSpec):
    """Detect importing names from the std namespace."""
    name = "import-std-name"
    message = "Importing names from the std namespace is not allowed in headers."
    suffixes = [".h"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        return not line.startswith("using std::")


class TestUsingDirectives(TestSpec):
    """Detect using directives in headers."""
    name = "using-directive"
    message = "Using directives are not allowed in headers."
    suffixes = [".h"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        return not line.startswith("using namespace")


class TestStdPrefix(TestSpec):
    """Detect missing std:: prefix for common names from the std namespace."""
    name = "std-prefix"
    message = "Use std:: prefix for names from the std namespace."
    suffixes = [".h", ".cxx", ".C"]
    prefix_bad = r"[^\w:\.\"]"
    patterns = [r"vector<", r"array[<\{\(]", r"f?abs\(", r"sqrt\(", r"pow\(", r"min\(", r"max\(", r"log\(", r"exp\(", r"sin\(", r"cos\(", r"tan\(", r"atan\(", r"atan2\("]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        for pattern in self.patterns:
            iterators = re.finditer(fr"{self.prefix_bad}{pattern}", line)
            matches = [(it.start(), it.group()) for it in iterators]
            if not matches:
                continue
            if not "\"" in line: # Found a match which cannot be inside a string.
                return False
            # Ignore matches inside strings.
            for match in matches:
                n_quotes_before = line.count("\"", 0, match[0]) # Count quotation marks before the match.
                if n_quotes_before % 2: # If odd, we are inside a string and we should ignore this match.
                    continue
                # We are not inside a string and this match is valid.
                return False
        return True


class TestROOT(TestSpec):
    """Detect use of unnecessary ROOT entities."""
    name = "root-entity"
    message = "Consider replacing ROOT entities with STD C++ or O2 entities."
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and not "Macros/" in path

    def test_line(self, line: str) -> bool:
        pattern = r"TMath|(U?(Int|Char|Short)|Double(32)?|Float(16)?|U?Long(64)?|Bool)_t"
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        return re.search(pattern, line) is None


class TestPI(TestSpec):
    """Detect use of external PI."""
    name = "external-pi"
    message = "Consider using the PI constant (and its multiples) defined in o2::constants::math."
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and not "Macros/" in path

    def test_line(self, line: str) -> bool:
        pattern = r"M_PI|TMath::(Two)?Pi"
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        return re.search(pattern, line) is None


class TestPdgDatabase(TestSpec):
    """Detect use of TDatabasePDG."""
    name = "pdg/database"
    message = "Direct use of TDatabasePDG is not allowed. Use o2::constants::physics::Mass... or Service<o2::framework::O2DatabasePDG>."
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and not "Macros/" in path

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        return not "TDatabasePDG" in line


class TestPdgCode(TestSpec):
    """Detect use of hard-coded PDG codes."""
    name = "pdg/explicit-code"
    message = "Avoid using hard-coded PDG codes. Use named values from PDG_t or o2::constants::physics::Pdg instead."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        if re.search(r"->(GetParticle|Mass)\([+-]?[0-9]+\)", line):
            return False
        match = re.search(r"[Pp][Dd][Gg][\w]* ={1,2} [+-]?([0-9]+);", line)
        if match:
            code = match.group(1)
            if code not in ("0", "1", "999"):
                return False
        return True


class TestPdgMass(TestSpec):
    """Detect unnecessary call of Mass() for a known PDG code."""
    name = "pdg/known-mass"
    message = "Consider using o2::constants::physics::Mass... instead of calling a database method for a known PDG code."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        pattern_pdg_code = r"[+-]?(k[A-Z][a-zA-Z0-9]*|[0-9]+)"
        if re.search(fr"->GetParticle\({pattern_pdg_code}\)->Mass\(\)", line):
            return False
        if re.search(fr"->Mass\({pattern_pdg_code}\)", line):
            return False
        return True


class TestLogging(TestSpec):
    """Detect non-O2 logging."""
    name = "logging"
    message = "Consider using O2 logging (LOG, LOGF, LOGP)."
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and not "Macros/" in path

    def test_line(self, line: str) -> bool:
        pattern = r"^([Pp]rintf\(|(std::)?cout <)"
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        return re.search(pattern, line) is None


class TestConstRefInForLoop(TestSpec):
    """Test const refs in range-based for loops."""
    name = "const-ref-in-for-loop"
    message = "Use constant references for non-modified iterators in range-based for loops."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        if not re.match(r"for \(.* :", line):
            return True
        line = line[:line.index(" :")] # keep only the iterator part
        return re.search(r"(\w const|const \w+)& ", line) is not None


class TestConstRefInSubscription(TestSpec):
    """Test const refs in process function subscriptions."""
    name = "const-ref-in-process"
    message = "Use constant references for table subscriptions in process functions."
    suffixes = [".cxx"]
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
            if self.is_disabled(line):
                continue
            if "//" in line: # Remove comment. (Ignore /* to avoid truncating at /*parameter*/.)
                line = line[:line.index("//")]
            if re.search(r"^void process[\w]*\(", line):
                line_process = (i + 1)
                i_closing = line.rfind(")")
                i_start = line.index("(") + 1
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
                template_args = re.findall(r"<[\w:, ]*>", arguments)
                if template_args:
                    for arg in template_args:
                        if ", " in arg:
                            arguments = arguments.replace(arg, arg.replace(", ", ":"))
                words = arguments.split(", ")
                # test
                for arg in words:
                    if not re.search(r"([\w>] const|const [\w<>:]+)&", arg):
                        passed = False
                        print_error(path, i + 1, self.name, f"Argument {arg} is not const&.")
                line_process = 0
        return passed


# Documentation
# Reference: https://rawgit.com/AliceO2Group/CodingGuidelines/master/comments_guidelines.html


class TestDocumentationFile(TestSpec):
    """Test mandatory documentation of C++ files."""
    name = "doc/file"
    message = "Provide mandatory file documentation."
    suffixes = [".h", ".cxx", ".C"]
    per_line = False

    def test_file(self, path : str, content) -> bool:
        passed = False
        doc_items = []
        doc_items.append({"keyword": "file", "pattern": fr"{os.path.basename(path)}$", "found": False})
        doc_items.append({"keyword": "brief", "pattern": r"\w.* \w", "found": False}) # at least two words
        doc_items.append({"keyword": "author", "pattern": r"[\w]+", "found": False})
        doc_prefix = "///"
        n_lines_copyright = 11
        last_doc_line = n_lines_copyright

        for i, line in enumerate(content):
            if i < n_lines_copyright: # Skip copyright lines.
                continue
            if line.strip() and not line.startswith(doc_prefix): # Stop at the first non-empty non-doc line.
                break
            if line.startswith(doc_prefix):
                last_doc_line = (i + 1)
            for item in doc_items:
                if re.search(fr"^{doc_prefix} [\\@]{item['keyword']} +{item['pattern']}", line):
                    item["found"] = True
                    # print_error(path, i + 1, self.name, f"Found \{item['keyword']}.")
                    break
            if all(item["found"] for item in doc_items): # All items have been found.
                passed = True
                break
        if not passed:
            for item in doc_items:
                if not item["found"]:
                    print_error(path, last_doc_line, self.name, f"Documentation for \\{item['keyword']} is missing, incorrect or misplaced.")
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
    name = "name/function-variable"
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
                line = line[:line.index(keyword)]

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
            funval_name = funval_name[:funval_name.index("[")]
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
    name = "name/macro"
    message = "Use SCREAMING_SNAKE_CASE for names of macros."
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
    name = "name/constexpr-constant"
    message = "Use UpperCamelCase for names of constexpr constants. Names of special constants may be prefixed with \"k\"."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        words = line.split()
        if not "constexpr" in words or not "=" in words:
            return True
        # Extract constant name.
        words = words[:words.index("=")] # keep only words before "="
        constant_name = words[-1] # last word before "="
        if constant_name.endswith("]") and "[" not in constant_name: # it's an array and we do not have the name before "[" here
            opens_brackets = ["[" in w for w in words]
            constant_name = words[opens_brackets.index(True)] # the name is in the first element with "["
        if "[" in constant_name: # Remove brackets for arrays.
            constant_name = constant_name[:constant_name.index("[")]
        if "::" in constant_name: # Remove the class prefix for methods.
            constant_name = constant_name.split("::")[-1]
        # The actual test comes here.
        if constant_name.startswith("k"): # exception for special constants
            constant_name = constant_name[1:] # test the name without "k"
        return is_upper_camel_case(constant_name)


class TestNameColumn(TestSpec):
    """Test names of O2 columns."""
    name = "name/o2-column"
    message = "Use UpperCamelCase for names of O2 columns and matching lowerCamelCase names for their getters."
    suffixes = [".h", ".cxx"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not (match := re.match(r"DECLARE(_[A-Z]+)*_COLUMN(_[A-Z]+)*\(", line)):
            return True
        # Extract names of the column type and getter.
        line = remove_comment_cpp(line)
        line = line[len(match.group()):].strip() # Extract part after "(".
        if not (match := re.match(r"([^,]+), ([^,\) ]+)", line)):
            print(f"Failed to extract column type and getter from \"{line}\".")
            return False
        column_type_name = match.group(1)
        column_getter_name = match.group(2)
        # print(f"Got \"{column_type_name}\" \"{column_getter_name}\"")
        # return True
        if column_type_name[0] == "_": # probably a macro variable
            return True
        # The actual test comes here.
        if not is_upper_camel_case(column_type_name):
            return False
        if not is_lower_camel_case(column_getter_name):
            return False
        return f"{column_type_name[0].lower()}{column_type_name[1:]}" == column_getter_name


class TestNameTable(TestSpec):
    """Test names of O2 tables."""
    name = "name/o2-table"
    message = "Use UpperCamelCase for names of O2 tables."
    suffixes = [".h", ".cxx"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not (match := re.match(r"DECLARE(_[A-Z]+)*_TABLES?(_[A-Z]+)*\(", line)):
            return True
        # Extract names of the column type and getter.
        line = remove_comment_cpp(line)
        line = line[len(match.group()):].strip() # Extract part after "(".
        if not (match := re.match(r"([^,\) ]+)", line)):
            print(f"Failed to extract table type from \"{line}\".")
            return False
        table_type_name = match.group(1)
        # print(f"Got \"{table_type_name}\"")
        # return True
        if table_type_name[0] == "_": # probably a macro variable
            return True
        # The actual test comes here.
        return is_upper_camel_case(table_type_name)


class TestNameNamespace(TestSpec):
    """Test names of namespaces."""
    name = "name/namespace"
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
    name = f"name/{keyword}"
    message = f"Use UpperCamelCase for names of {keyword}."
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not line.startswith(f"{self.keyword} "):
            return True
        line = remove_comment_cpp(line)
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
    name = f"name/enum"
    message = f"Use UpperCamelCase for names of enumerators and their values."


class TestNameClass(TestNameUpperCamelCase):
    """Test names of classes."""
    keyword = "class"
    name = f"name/class"
    message = f"Use UpperCamelCase for names of classes."


class TestNameStruct(TestNameUpperCamelCase):
    """Test names of structs."""
    keyword = "struct"
    name = f"name/struct"
    message = f"Use UpperCamelCase for names of structs."


class TestNameFileCpp(TestSpec):
    """Test names of C++ files."""
    name = "name/file-cpp"
    message = "Use lowerCamelCase or UpperCamelCase for names of C++ files. See the O2 naming conventions for details."
    suffixes = [".h", ".cxx", ".C"]
    per_line = False

    def test_file(self, path : str, content) -> bool:
        file_name = os.path.basename(path)
        # workflow file
        if re.search(r"(TableProducer|Tasks)/.*\.cxx", path):
            return is_lower_camel_case(file_name)
        # data model file
        if "DataModel/" in path:
            return is_upper_camel_case(file_name)
        # utility file
        if re.search(r"[Uu]til(ity|ities|s)?", file_name):
            return is_lower_camel_case(file_name)
        return is_camel_case(file_name)


class TestNameFilePython(TestSpec):
    """Test names of Python files."""
    name = "name/file-python"
    message = "Use snake_case for names of Python files."
    suffixes = [".py", ".ipynb"]
    per_line = False

    def test_file(self, path : str, content) -> bool:
        file_name = os.path.basename(path)
        return is_snake_case(file_name)


class TestNameWorkflow(TestSpec):
    """Test names of O2 workflows."""
    name = "name/o2-workflow"
    message = "Use kebab-case for names of workflows and match the name of the workflow file."
    suffixes = ["CMakeLists.txt"]
    per_line = False

    def test_file(self, path : str, content) -> bool:
        passed = True
        workflow_name = ""
        for i, line in enumerate(content):
            if not line.startswith("o2physics_add_dpl_workflow("):
                continue
            if self.is_disabled(line, "#"):
                continue
            # Extract workflow name.
            workflow_name = line.strip().split("(")[1].split()[0]
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
            expected_workflow_file_name = kebab_case_to_camel_case_l(workflow_name) + ".cxx"
            # Compare the actual and expected file names.
            if expected_workflow_file_name != workflow_file_name:
                passed = False
                print_error(path, i + 1, self.name, f"Workflow name {workflow_name} does not match its file name {workflow_file_name}. (Matches {expected_workflow_file_name}.)")
        return passed


class TestNameDevice(TestSpec):
    """Test names of devices."""
    name = "name/o2-device"
    message = "Specify device name only when it cannot be derived from the struct name. Only append to the default name."
    suffixes = [".cxx"]
    per_line = False

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path)

    def test_file(self, path : str, content) -> bool:
        is_inside_define = False
        is_inside_adapt = False
        struct_name = ""
        struct_templated = False # Is the struct templated?
        n_parens_opened = 0 # number of opened parentheses
        passed = True
        for i, line in enumerate(content):
            if not line.strip():
                continue
            if self.is_disabled(line, "#"):
                continue
            if is_comment_cpp(line):
                continue
            line = remove_comment_cpp(line)
            # Wait for defineDataProcessing.
            if not is_inside_define:
                if not re.match(r"(o2::)?(framework::)?WorkflowSpec defineDataProcessing\(", line):
                    continue
                # print(f"{i + 1}: Entering define.")
                is_inside_define = True
            # Return at the end of defineDataProcessing.
            if is_inside_define and line[0] == "}":
                # print(f"{i + 1}: Exiting define.")
                break
            # Wait for adaptAnalysisTask.
            if not is_inside_adapt:
                if (index := line.find("adaptAnalysisTask<")) == -1:
                    continue
                # print(f"{i + 1}: Entering adapt.")
                is_inside_adapt = True
                line = line[(index + len("adaptAnalysisTask<")):]
                # Extract struct name.
                if not (match := re.match(r"([^>]+)", line)):
                    print_error(path, i + 1, self.name, f"Failed to extract struct name from \"{line}\".")
                    return False
                struct_name = match.group(1)
                if (index := struct_name.find("<")) > -1:
                    struct_templated = True
                    # print(f"{i + 1}: Got templated struct name {struct_name}")
                    struct_name = struct_name[:index]
                # print(f"{i + 1}: Got struct name {struct_name}")
                line = line[(line.index(struct_name) + len(struct_name)):]
            if is_inside_adapt:
                n_parens_opened += line.count("(") - line.count(")")
                # print(f"{i + 1}: {n_parens_opened} opened parens")
                if n_parens_opened <= 0:
                    # print(f"{i + 1}: Exiting adapt.")
                    is_inside_adapt = False
                # Find device name.
                if not "TaskName{" in line:
                    continue
                passed = False
                # Extract device name.
                if not (match := re.search(r"TaskName\{\"([^\}]+)\"\}", line)):
                    print_error(path, i + 1, self.name, f"Failed to extract device name from \"{line}\".")
                    return False
                device_name = match.group(1)
                # print(f"{i + 1}: Got struct \"{struct_name}\" with device name \"{device_name}\".")
                # Test device name.
                expected_struct_name = kebab_case_to_camel_case_u(device_name)
                if struct_name == expected_struct_name:
                    print_error(path, i + 1, self.name, f"Specified device name {device_name} matches exactly the struct name {struct_name}. TaskName seems unnecessary.")
                # if template, check that device name is extension of struct name
                elif expected_struct_name.startswith(struct_name):
                    print_error(path, i + 1, self.name, f"Specified device name {device_name} is an extension of the struct name {struct_name}. Is the struct templated? {struct_templated}")
                else:
                    print_error(path, i + 1, self.name, f"Specified device name {device_name} does not match the struct name {struct_name}. (Matching {expected_struct_name})")
        return passed


# PWG-HF


class TestHfConstAuto(TestSpec):
    """PWGHF: Detect swapped const auto."""
    name = "pwghf/const-auto"
    message = "Use \"const auto\" instead of \"auto const\"."
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and "PWGHF/" in path

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        return not "auto const" in line


class TestHfNameStructClass(TestSpec):
    """PWGHF: Test names of structs and classes."""
    name = "pwghf/name/struct-class"
    message = "Names of PWGHF structs and classes must start with \"Hf\"."
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and "PWGHF/" in path

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not line.startswith(("struct ", "class ")):
            return True
        line = remove_comment_cpp(line)
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
    name = "pwghf/struct-member-order"
    message = "Declare struct members in the conventional order. See the PWGHF coding guidelines."
    suffixes = [".cxx"]
    per_line = False
    member_order = ["Spawns<", "Builds<", "Produces<", "Configurable<", "HfHelper ", "SliceCache ", "Service<", "using ", "Filter ", "Preslice<", "Partition<", "ConfigurableAxis ", "AxisSpec ", "HistogramRegistry ", "OutputObj<", "void init(", "void process"]

    def file_matches(self, path: str) -> bool:
        return TestSpec.file_matches(self, path) and "PWGHF/" in path

    def test_file(self, path : str, content) -> bool:
        passed = True
        dic_struct:dict[str, dict] = {}
        struct_name = ""
        for i, line in enumerate(content):
            if is_comment_cpp(line):
                continue
            if line.startswith("struct "): # expecting no indentation
                line = remove_comment_cpp(line)
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
    name = "pwghf/name/workflow-file"
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
    name = "pwghf/name/configurable"
    message = "Use lowerCamelCase for names of configurables and use the same name for the struct member as for the JSON string."
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
        tests.append(TestIOStream())
        tests.append(TestUsingStd())
        tests.append(TestUsingDirectives())
        tests.append(TestStdPrefix())
        tests.append(TestROOT())
        tests.append(TestPI())
        tests.append(TestPdgDatabase())
        tests.append(TestPdgCode())
        tests.append(TestPdgMass())
        tests.append(TestLogging())
        tests.append(TestConstRefInForLoop())
        tests.append(TestConstRefInSubscription())

    # Documentation
    enable_documentation = True
    if enable_documentation:
        tests.append(TestDocumentationFile())

    # Naming conventions
    enable_naming = True
    if enable_naming:
        tests.append(TestNameFunctionVariable())
        tests.append(TestNameMacro())
        tests.append(TestNameConstant())
        tests.append(TestNameColumn())
        tests.append(TestNameTable())
        tests.append(TestNameNamespace())
        tests.append(TestNameEnum())
        tests.append(TestNameClass())
        tests.append(TestNameStruct())
        tests.append(TestNameFileCpp())
        tests.append(TestNameFilePython())
        tests.append(TestNameWorkflow())
        tests.append(TestNameDevice())

    # PWG-HF
    enable_pwghf = True
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
    # print(args.paths)
    print("Enabled tests:", test_names)
    print("Suffixes of tested files:", suffixes)
    # print(f"Github annotations: {github_mode}.")

    # Test files.
    for path in args.paths:
        # print(f"Processing path \"{path}\".")
        # Skip not tested files.
        if not path.endswith(suffixes):
            # print(f"Skipping path \"{path}\".")
            continue
        try:
            with open(path, "r") as file:
                content = file.readlines()
                for test in tests:
                    result = test.run(path, content)
                    if not result:
                        passed = False
                    # print(f"File \"{path}\" {'passed' if result else 'failed'} the test {test.name}.")
        except IOError:
            print(f"Failed to open file \"{path}\".")
            sys.exit(1)

    # Report global result.
    if passed:
        print("All tests passed.")
    else:
        print("Issues have been found.")
        sys.exit(1)


if __name__ == "__main__":
    main()
