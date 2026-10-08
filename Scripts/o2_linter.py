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
@brief  O2 linter (Find O2-specific issues in O2 code)
@author Vít Kučera <vit.kucera@cern.ch>, Inha University
@date   2024-07-14
"""

import argparse
import os
import re
import sys
from enum import Enum
from pathlib import Path
from typing import Union

github_mode = False  # GitHub mode
prefix_disable = "o2-linter: disable="  # prefix for disabling tests
file_config = "o2linter_config"  # name of the configuration file (applied per directory)
# If this file exists in the path of the tested file,
# failures of tests listed in this file will not make the linter fail.


# issue severity levels
class Severity(Enum):
    WARNING = 1
    ERROR = 2
    DEFAULT = ERROR


# strings for error messages
message_levels = {
    Severity.WARNING: "warning",
    Severity.ERROR: "error",
}


class Reference(Enum):
    ISO_CPP = 1
    LLVM = 2
    GOOGLE = 3
    PY_ZEN = 4
    PY_PEP8 = 5
    O2 = 6
    PWG_HF = 7
    LINTER_1 = 8
    LINTER_2 = 9


references_list: "list[tuple[Reference, str, str]]" = [
    (Reference.ISO_CPP, "C++ Core Guidelines", "https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines"),
    (Reference.LLVM, "LLVM Coding Standards", "https://llvm.org/docs/CodingStandards.html"),
    (Reference.GOOGLE, "Google C++ Style Guide", "https://google.github.io/styleguide/cppguide.html"),
    (Reference.PY_ZEN, "The Zen of Python", "https://peps.python.org/pep-0020/"),
    (Reference.PY_PEP8, "Style Guide for Python Code", "https://peps.python.org/pep-0008/"),
    (Reference.O2, "ALICE O2 Coding Guidelines", "https://github.com/AliceO2Group/CodingGuidelines"),
    (
        Reference.PWG_HF,
        "PWG-HF guidelines",
        "https://aliceo2group.github.io/analysis-framework/docs/advanced-specifics/pwghf.html#contribute",
    ),
    (
        Reference.LINTER_1,
        "Proposal of the O2 linter",
        "https://indico.cern.ch/event/1482467/#29-development-of-the-o2-linte",
    ),
    (
        Reference.LINTER_2,
        "Update of the O2 linter",
        "https://indico.cern.ch/event/1513748/#29-o2-linter-development",
    ),
]

references: "dict[Reference, dict]" = {name: {"title": title, "url": url} for name, title, url in references_list}


def is_camel_case(name: str) -> bool:
    """forExample or ForExample"""
    return "_" not in name and "-" not in name and " " not in name


def is_upper_camel_case(name: str) -> bool:
    """ForExample"""
    if not name:
        return False
    return name[0].isupper() and is_camel_case(name)


def is_lower_camel_case(name: str) -> bool:
    """forExample"""
    if not name:
        return False
    return name[0].islower() and is_camel_case(name)


def is_kebab_case(name: str) -> bool:
    """for-example"""
    if not name:
        return False
    return name.islower() and "_" not in name and " " not in name


def is_snake_case(name: str) -> bool:
    """for_example"""
    if not name:
        return False
    return name.islower() and "-" not in name and " " not in name


def is_screaming_snake_case(name: str) -> bool:
    """FOR_EXAMPLE"""
    if not name:
        return False
    return name.isupper() and "-" not in name and " " not in name


def kebab_case_to_camel_case_u(line: str) -> str:
    """Convert kebab-case string to UpperCamelCase string."""
    return "".join([w.title() if w[0].isnumeric() else w.capitalize() for w in line.split("-")])


def kebab_case_to_camel_case_l(line: str) -> str:
    """Convert kebab-case string to lowerCamelCase string."""
    new_line = kebab_case_to_camel_case_u(line)
    return f"{new_line[0].lower()}{new_line[1:]}"  # start with lowercase letter


def camel_case_to_kebab_case(line: str) -> str:
    """Convert CamelCase string to kebab-case string.
    As done in O2/Framework/Foundation/include/Framework/TypeIdHelpers.h:type_to_task_name
    """
    if not line.strip():
        return line
    new_line = []
    for i, c in enumerate(line):
        if i > 0 and c.isupper() and line[i - 1] != "-":
            new_line.append("-")
        new_line.append(c.lower())
    return "".join(new_line)


def is_comment_cpp(line: str) -> bool:
    """Test whether a line is a C++ comment."""
    return line.strip().startswith(("//", "/*"))


def remove_comment_cpp(line: str) -> str:
    """Remove C++ comments from the end of a line."""
    for keyword in ("//", "/*"):
        if keyword in line:
            line = line[: line.index(keyword)]
    return line.strip()


def block_ranges(line: str, char_open: str, char_close: str) -> "list[list[int]]":
    """Get list of index ranges of longest blocks opened with char_open and closed with char_close."""
    # print(f"Looking for {char_open}{char_close} blocks in \"{line}\".")
    # print(line)
    list_ranges: list[list[int]] = []
    if not all((line, len(char_open) == 1, len(char_close) == 1)):
        return list_ranges

    def direction(char: str) -> int:
        if char == char_open:
            return 1
        if char == char_close:
            return -1
        return 0

    list_levels = []  # list of block levels (net number of opened blocks)
    level_sum = 0  # current block level (sum of previous directions)
    for char in line:
        list_levels.append(level_sum := level_sum + direction(char))
    level_min = min(list_levels)  # minimum level (!= 0 if line has opened blocks)
    # Look for openings (level_min + 1) and closings (level_min).
    index_start = -1
    is_opened = False
    # print(list_levels)
    for i, level in enumerate(list_levels):
        if not is_opened and level > level_min:
            is_opened = True
            index_start = i
            # print(f"Opening at {i}")
        elif is_opened and (level == level_min or i == len(list_levels) - 1):
            is_opened = False
            list_ranges.append([index_start, i])
            # print(f"Closing at {i}")
    # print(f"block_ranges: Found {len(list_ranges)} blocks: {list_ranges}.")
    if is_opened:
        print("block_ranges: Block left opened.")
    return list_ranges


def get_tolerated_tests(path: str) -> "list[str]":
    """Get the list of tolerated tests.

    Looks for the configuration file.
    Starts in the test file directory and iterates through parents.
    """
    tests: list[str] = []
    for directory in Path(path).resolve().parents:
        path_tests = directory / file_config
        if path_tests.is_file():
            with path_tests.open() as content:
                tests = [line.strip() for line in content.readlines() if line.strip()]
                print(f"{path}:1: info: Tolerating tests from {path_tests}. {tests}")
            break
    return tests


class TestSpec:
    """Prototype of a test class"""

    name: str = "test-template"  # short name of the test
    message: str = "Test failed"  # error message
    rationale: str = "Rationale missing"  # brief rationale explanation
    references: "list[Reference]" = []  # list of references relevant for this rule
    severity_default: Severity = Severity.DEFAULT
    severity_current: Severity = Severity.DEFAULT
    suffixes: "list[str]" = []  # suffixes of files to test
    per_line: bool = True  # Test lines separately one by one.
    tolerated: bool = False  # flag for tolerating issues
    n_issues: int = 0  # issue counter
    n_disabled: int = 0  # counter of disabled issues
    n_tolerated: int = 0  # counter of tolerated issues

    def file_matches(self, path: str) -> bool:
        """Test whether the path matches the pattern for files to test."""
        return path.endswith(tuple(self.suffixes)) if self.suffixes else True

    def is_disabled(self, line: str, prefix_comment="//") -> bool:
        """Detect whether the test is explicitly disabled."""
        for prefix in [prefix_comment, prefix_disable]:
            if prefix not in line:
                return False
            line = line[(line.index(prefix) + len(prefix)) :]  # Strip away part before prefix.
        if self.name in line:
            self.n_disabled += 1
            # Look for a comment with a reason for disabling.
            if re.search(r" \(.{3,}\)", line):
                return True
        return False

    def print_error(self, path: str, line: Union[int, None], message: str):
        """Format and print error message."""
        # return # Use to suppress error messages.
        line = line or 1
        # terminal format
        print(f"{path}:{line}: {message_levels[self.severity_current]}: {message} [{self.name}]")
        if github_mode and not self.tolerated:  # Annotate only not tolerated issues.
            # GitHub annotation format
            print(f"::{message_levels[self.severity_current]} file={path},line={line},title=[{self.name}]::{message}")

    def test_line(self, line: str) -> bool:
        """Test a line."""
        raise NotImplementedError()

    def test_file(self, path: str, content) -> bool:
        """Test a file in a way that cannot be done line by line."""
        raise NotImplementedError()

    def run(self, path: str, content: "list[str]") -> bool:
        """Run the test."""
        # print(content)
        passed = True
        self.severity_current = Severity.WARNING if self.tolerated else self.severity_default
        if not self.file_matches(path):
            return passed
        # print(f"Running test {self.name} for {path} with {len(content)} lines")
        if self.per_line:
            for i, line in enumerate(content):
                if not isinstance(self, TestUsingDirective):  # Keep the indentation if needed.
                    line = line.strip()
                if not line:
                    continue
                # print(i + 1, line)
                if self.is_disabled(line):
                    continue
                if not self.test_line(line):
                    passed = False
                    if self.tolerated:
                        self.n_tolerated += 1
                    else:
                        self.n_issues += 1
                    self.print_error(path, i + 1, self.message)
        else:
            passed = self.test_file(path, content)
            if not passed:
                if self.tolerated:
                    self.n_tolerated += 1
                else:
                    self.n_issues += 1
                self.print_error(path, None, self.message)
        return passed or self.tolerated


##########################
# Implementations of tests
##########################

# Bad practice


class TestIoStream(TestSpec):
    """Detect included iostream."""

    name = "include-iostream"
    message = "Do not include iostream. Use O2 logging instead."
    rationale = "Performance. Avoid injection of static constructors. Consistent logging."
    references = [Reference.LLVM, Reference.LINTER_1]
    suffixes = [".h", ".cxx"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        return not line.startswith("#include <iostream>")


class TestUsingStd(TestSpec):
    """Detect importing names from the std namespace."""

    name = "import-std-name"
    message = "Do not import names from the std namespace in headers."
    rationale = "Code safety. Avoid namespace pollution with common names."
    references = [Reference.LINTER_1]
    suffixes = [".h"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        return not line.startswith("using std::")


class TestUsingDirective(TestSpec):
    """Detect using directives in headers."""

    name = "using-directive"
    message = "Do not put using directives at global scope in headers."
    rationale = "Code safety. Avoid namespace pollution."
    references = [Reference.O2, Reference.ISO_CPP, Reference.LLVM, Reference.GOOGLE, Reference.LINTER_1]
    suffixes = [".h"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        return not line.startswith("using namespace")


class TestStdPrefix(TestSpec):
    """Detect missing std:: prefix for common names from the std namespace."""

    name = "std-prefix"
    message = "Use std:: prefix for names from the std namespace."
    rationale = "Code clarity, safety and portability. Avoid ambiguity (e.g. abs)."
    references = [Reference.LLVM, Reference.LINTER_1]
    suffixes = [".h", ".cxx", ".C"]
    prefix_bad = r"[^\w:\.\"]"
    patterns = [
        r"vector<",
        r"array[<\{\(]",
        r"f?abs\(",
        r"sqrt\(",
        r"pow\(",
        r"min\(",
        r"max\(",
        r"log(2|10)?\(",
        r"exp\(",
        r"a?(sin|cos|tan)h?\(",
        r"atan2\(",
        r"erfc?\(",
        r"hypot\(",
    ]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        for pattern in self.patterns:
            iterators = re.finditer(rf"{self.prefix_bad}{pattern}", line)
            matches = [(it.start(), it.group()) for it in iterators]
            if not matches:
                continue
            if '"' not in line:  # Found a match which cannot be inside a string.
                return False
            # Ignore matches inside strings.
            for match in matches:
                n_quotes_before = line.count('"', 0, match[0])  # Count quotation marks before the match.
                if not n_quotes_before % 2:  # If even, we are not inside a string and this match is valid.
                    return False
        return True


class TestRootEntity(TestSpec):
    """Detect unnecessary use of ROOT entities."""

    name = "root/entity"
    message = "Replace ROOT entities with equivalents from standard C++ or from O2."
    rationale = "Code simplicity and maintainability. O2 is not a ROOT code."
    references = [Reference.ISO_CPP, Reference.LINTER_1, Reference.PY_ZEN]
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return super().file_matches(path) and "Macros/" not in path

    def test_line(self, line: str) -> bool:
        pattern = (
            r"TMath::(Abs|Sqrt|Power|Min|Max|Log(2|10)?|Exp|A?(Sin|Cos|Tan)H?|ATan2|Erfc?|Hypot)\(|"
            r"(U?(Int|Char|Short)|Double(32)?|Float(16)?|U?Long(64)?|Bool)_t"
        )
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        return re.search(pattern, line) is None


class TestRootLorentzVector(TestSpec):
    """Detect use of TLorentzVector."""

    name = "root/lorentz-vector"
    message = (
        "Do not use the TLorentzVector legacy class. "
        "Use std::array with RecoDecay methods or the ROOT::Math::LorentzVector template instead."
    )
    rationale = "Performance. Use up-to-date tools."
    references = [Reference.LINTER_2]
    suffixes = [".h", ".cxx"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        return "TLorentzVector" not in line


class TestPi(TestSpec):
    """Detect use of external pi."""

    name = "external-pi"
    message = "Use the PI constant (and its multiples and fractions) defined in o2::constants::math."
    rationale = "Code maintainability."
    references = [Reference.LINTER_1, Reference.PY_ZEN]
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return super().file_matches(path) and "Macros/" not in path

    def test_line(self, line: str) -> bool:
        pattern = r"[^\w]M_PI|TMath::(Two)?Pi"
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        return re.search(pattern, line) is None


class TestTwoPiAddSubtract(TestSpec):
    """Detect adding/subtracting of 2 pi."""

    name = "two-pi-add-subtract"
    message = "Use RecoDecay::constrainAngle to restrict angle to a given range."
    rationale = "Code maintainability and safety. Use existing tools."
    references = [Reference.ISO_CPP, Reference.LINTER_1, Reference.PY_ZEN]
    suffixes = [".h", ".cxx"]

    def test_line(self, line: str) -> bool:
        pattern_two_pi = (
            r"(2(\.0*f?)? \* (M_PI|TMath::Pi\(\)|(((o2::)?constants::)?math::)?PI)|"
            r"(((o2::)?constants::)?math::)?TwoPI|TMath::TwoPi\(\))"
        )
        pattern = rf"[\+-]=? {pattern_two_pi}"
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        return re.search(pattern, line) is None


class TestPiMultipleFraction(TestSpec):
    """Detect multiples/fractions of pi for existing equivalent constants."""

    name = "pi-multiple-fraction"
    message = "Use multiples/fractions of PI defined in o2::constants::math."
    rationale = "Code maintainability."
    references = [Reference.LINTER_1, Reference.PY_ZEN]
    suffixes = [".h", ".cxx"]

    def test_line(self, line: str) -> bool:
        pattern_pi = r"(M_PI|TMath::(Two)?Pi\(\)|(((o2::)?constants::)?math::)?(Two)?PI)"
        pattern_multiple = r"(2(\.0*f?)?|0\.2?5f?) \* "  # * 2, 0.25, 0.5
        pattern_fraction = r" / ((2|3|4)([ ,;\)]|\.0*f?))"  # / 2, 3, 4
        pattern = rf"{pattern_multiple}{pattern_pi}[^\w]|{pattern_pi}{pattern_fraction}"
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        return re.search(pattern, line) is None


class TestPdgDatabase(TestSpec):
    """Detect use of TDatabasePDG."""

    name = "pdg/database"
    message = (
        "Do not use TDatabasePDG directly. "
        "Use o2::constants::physics::Mass... or Service<o2::framework::O2DatabasePDG> instead."
    )
    rationale = "Performance."
    references = [Reference.LINTER_1]
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return super().file_matches(path) and "Macros/" not in path

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        return "TDatabasePDG" not in line


class TestPdgExplicitCode(TestSpec):
    """Detect hard-coded PDG codes."""

    name = "pdg/explicit-code"
    message = "Avoid hard-coded PDG codes. Use named values from PDG_t or o2::constants::physics::Pdg instead."
    rationale = "Code comprehensibility, readability, maintainability and safety."
    references = [Reference.O2, Reference.ISO_CPP, Reference.LINTER_1]
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        if re.search(r"->(GetParticle|Mass)\([+-]?[0-9]+\)", line):
            return False
        match = re.search(r"[Pp][Dd][Gg].* !?={1,2} [+-]?([0-9]+)", line)
        if match:
            code = match.group(1)
            if code not in ("0", "1", "999"):
                return False
        return True


class TestPdgExplicitMass(TestSpec):
    """Detect hard-coded particle masses."""

    name = "pdg/explicit-mass"
    message = "Avoid hard-coded particle masses. Use o2::constants::physics::Mass... instead."
    rationale = "Code comprehensibility, readability, maintainability and safety."
    references = [Reference.O2, Reference.ISO_CPP, Reference.LINTER_2]
    suffixes = [".h", ".cxx"]
    masses: "list[str]" = []  # list of mass values to detect

    def __init__(self) -> None:
        super().__init__()

        # List of masses of commonly used particles
        self.masses.append(r"0\.000511")  # e
        self.masses.append(r"0\.105")  # μ
        self.masses.append(r"0\.139")  # π+
        self.masses.append(r"0\.493")  # K+
        self.masses.append(r"0\.497")  # K0
        self.masses.append(r"0\.938")  # p
        self.masses.append(r"0\.939")  # n
        self.masses.append(r"1\.115")  # Λ
        self.masses.append(r"1\.864")  # D0
        self.masses.append(r"2\.286")  # Λc
        self.masses.append(r"3\.096")  # J/ψ

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        iterators = re.finditer(rf"(^|\D)({'|'.join(self.masses)})", line)
        matches = [(it.start(), it.group(2)) for it in iterators]
        if not matches:
            return True
        if '"' not in line:  # Found a match which cannot be inside a string.
            return False
        # Ignore matches inside strings.
        for match in matches:
            n_quotes_before = line.count('"', 0, match[0])  # Count quotation marks before the match.
            if not n_quotes_before % 2:  # If even, we are not inside a string and this match is valid.
                return False
        return True


class TestPdgKnownMass(TestSpec):
    """Detect unnecessary call of Mass() for a known PDG code."""

    name = "pdg/known-mass"
    message = "Use o2::constants::physics::Mass... instead of calling a database method for a known PDG code."
    rationale = "Performance."
    references = [Reference.LINTER_1]
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        pattern_pdg_code = r"[+-]?(k[A-Z][a-zA-Z0-9]*|[0-9]+)"
        if re.search(rf"->GetParticle\({pattern_pdg_code}\)->Mass\(\)", line):
            return False
        return not re.search(rf"->Mass\({pattern_pdg_code}\)", line)


class TestLogging(TestSpec):
    """Detect non-O2 logging."""

    name = "logging"
    message = "Use O2 logging (LOG, LOGF, LOGP)."
    rationale = "Logs easy to read and process."
    references = [Reference.LINTER_1]
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return super().file_matches(path) and "Macros/" not in path

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
    rationale = "Performance, code comprehensibility and safety."
    references = [Reference.O2, Reference.ISO_CPP, Reference.LLVM]
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        if not re.match(r"for \(.* :", line):
            return True
        line = line[: line.index(" :")]  # keep only the iterator part
        return re.search(r"(\w const|const \w+)& ", line) is not None


class TestConstRefInSubscription(TestSpec):
    """Test const refs in process function subscriptions.
    Test only top-level process functions (called process() or has PROCESS_SWITCH).
    """

    name = "const-ref-in-process"
    message = "Use constant references for table subscriptions in process functions."
    rationale = "Performance, code comprehensibility and safety."
    references = [Reference.O2, Reference.ISO_CPP, Reference.LINTER_1]
    suffixes = [".cxx"]
    per_line = False

    def test_file(self, path: str, content) -> bool:
        passed = True
        n_parens_opened = 0  # number of opened parentheses
        arguments = ""  # process function arguments
        line_process = 0  # line number of the process function
        # Find names of all top-level process functions.
        names_functions = ["process"]  # names of allowed process functions to test
        for i, line in enumerate(content):
            line = line.strip()
            if is_comment_cpp(line):
                continue
            if not line.startswith("PROCESS_SWITCH"):
                continue
            words = line.split()
            if len(words) < 2:
                passed = False
                self.print_error(
                    path,
                    i + 1,
                    "Failed to get the process function name. Keep it on the same line as the switch.",
                )
                continue
            names_functions.append(words[1][:-1])  # Remove the trailing comma.
            # self.print_error(path, i + 1, f"Got process function name {words[1][:-1]}.")
        # Test process functions.
        for i, line in enumerate(content):
            line = line.strip()
            if is_comment_cpp(line):
                continue
            if self.is_disabled(line):
                continue
            if "//" in line:  # Remove comment. (Ignore /* to avoid truncating at /*parameter*/.)
                line = line[: line.index("//")]
            if (match := re.match(r"void (process[\w]*)\(", line)) and match.group(1) in names_functions:
                line_process = i + 1
                i_closing = line.rfind(")")
                i_start = line.index("(") + 1
                i_end = i_closing if i_closing != -1 else len(line)
                arguments = line[i_start:i_end]  # get arguments between parentheses
                n_parens_opened = line.count("(") - line.count(")")
            elif n_parens_opened > 0:
                i_closing = line.rfind(")")
                i_start = 0
                i_end = i_closing if i_closing != -1 else len(line)
                arguments += " " + line[i_start:i_end]  # get arguments between parentheses
                n_parens_opened += line.count("(") - line.count(")")
            if line_process > 0 and n_parens_opened == 0:
                # Process arguments.
                # Sanitise arguments with spaces between <>.
                for start, end in block_ranges(arguments, "<", ">"):
                    arg = arguments[start : (end + 1)]
                    # print(f"Found argument \"{arg}\" in [{start}, {end}]")
                    if ", " in arg:
                        arguments = arguments.replace(arg, arg.replace(", ", "__"))
                # Extract arguments.
                words = arguments.split(", ")
                # Test each argument.
                for arg in words:
                    if not re.search(r"([\w>] const|const [\w<>:]+)&", arg):
                        passed = False
                        self.print_error(path, i + 1, f"Argument {arg} is not const&.")
                line_process = 0
        return passed


class TestWorkflowOptions(TestSpec):
    """Detect usage of workflow options in defineDataProcessing. (Not supported on AliHyperloop.)"""

    name = "o2-workflow-options"
    message = (
        "Do not use workflow options to customise workflow topology composition in defineDataProcessing. "
        "Use process function switches or metadata instead."
    )
    rationale = "Not supported on AliHyperloop."
    references = [Reference.LINTER_1]
    suffixes = [".cxx"]
    per_line = False

    def test_file(self, path: str, content) -> bool:
        is_inside_define = False  # Are we inside defineDataProcessing?
        for _i, line in enumerate(content):  # pylint: disable=unused-variable
            if not line.strip():
                continue
            if self.is_disabled(line):
                continue
            if is_comment_cpp(line):
                continue
            line = remove_comment_cpp(line)
            # Wait for defineDataProcessing.
            if not is_inside_define:
                if not re.match(r"((o2::)?framework::)?WorkflowSpec defineDataProcessing\(", line):
                    continue
                # print(f"{i + 1}: Entering define.")
                is_inside_define = True
            # Return at the end of defineDataProcessing.
            if is_inside_define and line[0] == "}":
                # print(f"{i + 1}: Exiting define.")
                break
            # Detect options.
            if ".options()" in line:
                return False
        return True


class TestMagicNumber(TestSpec):
    """Detect magic numbers."""

    name = "magic-number"
    message = "Avoid magic numbers in expressions. Assign the value to a clearly named variable or constant."
    rationale = "Code comprehensibility, maintainability and safety."
    references = [Reference.O2, Reference.ISO_CPP, Reference.LINTER_2]
    suffixes = [".h", ".cxx", ".C"]
    pattern_compare = r"([<>]=?|[!=]=)"
    pattern_number = r"[\+-]?([\d\.]+(e[\+-]?\d+)?)f?"

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        iterators = re.finditer(
            rf" {self.pattern_compare} {self.pattern_number}|\W{self.pattern_number} {self.pattern_compare} ", line
        )
        matches = [(it.start(), it.group(2), it.group(4)) for it in iterators]
        if not matches:
            return True
        # Ignore matches inside strings.
        for match in matches:
            n_quotes_before = line.count('"', 0, match[0])  # Count quotation marks before the match.
            if n_quotes_before % 2:  # If odd, we are inside a string and we should ignore this match.
                continue
            # We are not inside a string and this match is valid.
            for match_n in (match[1], match[2]):
                # Accept only 0 or 1 (int or float).
                if (match_n is not None) and (re.match(r"[01](\.0?)?$", match_n) is None):
                    return False
        return True


# Documentation
# Reference: https://rawgit.com/AliceO2Group/CodingGuidelines/master/comments_guidelines.html


class TestDocumentationFile(TestSpec):
    """Test mandatory documentation of C++ files."""

    name = "doc/file"
    message = "Provide mandatory file documentation."
    rationale = "Code comprehensibility. Collaboration."
    references = [Reference.O2, Reference.LINTER_1]
    suffixes = [".h", ".cxx", ".C"]
    per_line = False

    def test_file(self, path: str, content) -> bool:
        passed = False
        doc_items = []
        doc_items.append({"keyword": "file", "pattern": rf"{os.path.basename(path)}$", "found": False})
        doc_items.append({"keyword": "brief", "pattern": r"\w.* \w", "found": False})  # at least two words
        doc_items.append({"keyword": "author", "pattern": r"[\w]+", "found": False})
        doc_prefix = "///"
        n_lines_copyright = 11
        last_doc_line = n_lines_copyright

        for i, line in enumerate(content):
            if i < n_lines_copyright:  # Skip copyright lines.
                continue
            if line.strip() and not line.startswith(doc_prefix):  # Stop at the first non-empty non-doc line.
                break
            if line.startswith(doc_prefix):
                last_doc_line = i + 1
            for item in doc_items:
                if re.search(rf"^{doc_prefix} [\\@]{item['keyword']} +{item['pattern']}", line):
                    item["found"] = True
                    # self.print_error(path, i + 1, f"Found \{item['keyword']}.")
                    break
            if all(item["found"] for item in doc_items):  # All items have been found.
                passed = True
                break
        if not passed:
            for item in doc_items:
                if not item["found"]:
                    self.print_error(
                        path,
                        last_doc_line,
                        f"Documentation for \\{item['keyword']} is missing, incorrect or misplaced.",
                    )
        return passed


# Naming conventions
# Reference: https://rawgit.com/AliceO2Group/CodingGuidelines/master/naming_formatting.html

rationale_names = "Code readability, comprehensibility and searchability."
references_names = [Reference.O2, Reference.LINTER_1]


class TestNameFunctionVariable(TestSpec):
    """Test names of functions and of most variables.
    Might report false positives.
    Does not detect multiple variable declarations, i.e. "type name1, name2;"
    Does not detect function arguments on the same line as the function declaration.
    Does not detect multi-line declarations.
    Does not check capitalisation for constexpr because of special rules for constants. See TestNameConstant.
    """

    name = "name/function-variable"
    message = "Use lowerCamelCase for names of functions and variables."
    rationale = rationale_names
    references = references_names
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        # Look for declarations of functions and variables.

        # Strip away irrelevant remainders of the line after the object name.
        # For functions, stripping after "(" is enough but this way we also identify many declarations of variables.
        for keyword in ("(", "{", ";", " = ", "//", "/*"):
            if keyword in line:
                line = line[: line.index(keyword)]

        # Check the words.
        words = line.split()

        # number of words
        if len(words) < 2:
            return True

        # First word starts with a letter.
        if not words[0][0].isalpha():
            return True

        # Reject false positives with same structure.
        if words[0] in (
            "return",
            "if",
            "else",
            "new",
            "delete",
            "delete[]",
            "case",
            "typename",
            "using",
            "typedef",
            "enum",
            "namespace",
            "struct",
            "class",
            "explicit",
            "concept",
            "throw",
        ):
            return True
        if len(words) > 2 and words[1] in ("typename", "class", "struct"):
            return True

        # Identify the position of the name for cases "name[n + m]".
        funval_name = words[-1]  # expecting the name in the last word
        if (
            funval_name.endswith("]") and "[" not in funval_name
        ):  # it's an array and we do not have the name before "[" here
            opens_brackets = ["[" in w for w in words]
            if not any(opens_brackets):  # The opening "[" is not on this line. We have to give up.
                return True
            index_name = opens_brackets.index(True)  # the name is in the first element with "["
            funval_name = words[index_name]
            words = words[: (index_name + 1)]  # Strip away words after the name.
            if len(words) < 2:  # Check the adjusted number of words.
                return True

        # All words before the name start with an alphanumeric character (underscores not allowed).
        # Rejects expressions, e.g. * = += << }, but accepts numbers in array declarations.
        if not all(w[0].isalnum() for w in words[:-1]):
            return True

        # Extract function/variable name.
        if "[" in funval_name:  # Remove brackets for arrays.
            funval_name = funval_name[: funval_name.index("[")]
        if "::" in funval_name:  # Remove the class prefix for methods.
            funval_name = funval_name.split("::")[-1]

        # Check the name candidate.

        # Names of variables and functions are identifiers.
        if not funval_name.isidentifier():  # should be same as ^[\w]+$
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
    message = "Use SCREAMING_SNAKE_CASE for names of macros. Leading and double underscores are not allowed."
    rationale = rationale_names
    references = references_names
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
        if macro_name.startswith("_"):
            return False
        if "__" in macro_name:
            return False
        return is_screaming_snake_case(macro_name)


class TestNameConstant(TestSpec):
    """Test constexpr constant names."""

    name = "name/constexpr-constant"
    message = "Use UpperCamelCase for names of constexpr constants."
    rationale = rationale_names
    references = references_names
    suffixes = [".h", ".cxx", ".C"]

    def __init__(self) -> None:
        super().__init__()
        keyword = r"(.+ )"  # e.g. "static "
        type_val = r"([\w:<>+\-*\/, ]+ )"  # value type e.g. "std::array<Type, n + 1> "
        prefix = r"(\w+::)"  # prefix with namespace or class, e.g. "MyClass::"
        name_val = r"(\w+)"  # name of the constant
        array = r"(\[.*\])"  # array declaration: "[...]"
        assignment = r"( =|\(\d|{)"  # value assignment, e.g. " = 2", " = expression", "(2)", "{2}", "{{...}}"
        self.pattern = re.compile(rf"{keyword}?constexpr {type_val}?{prefix}*{name_val}{array}?{assignment}")

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        line = remove_comment_cpp(line)
        if not (match := self.pattern.match(line)):
            return True
        constant_name = match.group(4)
        # The actual test comes here.
        return is_upper_camel_case(constant_name)


class TestNameColumn(TestSpec):
    """Test names of O2 columns."""

    name = "name/o2-column"
    message = "Use UpperCamelCase for names of O2 columns and matching lowerCamelCase names for their getters."
    rationale = rationale_names
    references = references_names
    suffixes = [".h", ".cxx"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not (match := re.match(r"DECLARE(_[A-Z]+)*_COLUMN(_[A-Z]+)*\(", line)):
            return True
        # Extract names of the column type and getter.
        line = remove_comment_cpp(line)
        line = line[len(match.group()) :].strip()  # Extract part after "(".
        if not (match := re.match(r"([^,]+), ([^,\) ]+)", line)):
            print(f'Failed to extract column type and getter from "{line}".')
            return False
        column_type_name = match.group(1)
        column_getter_name = match.group(2)
        # print(f"Got \"{column_type_name}\" \"{column_getter_name}\"")
        # return True
        if column_type_name[0] == "_":  # probably a macro variable
            return True
        if "#" in column_type_name:  # Remove "#" for strings in macros.
            column_type_name = column_type_name[: column_type_name.index("#")]
        if "#" in column_getter_name:  # Remove "#" for strings in macros.
            column_getter_name = column_getter_name[: column_getter_name.index("#")]
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
    rationale = rationale_names
    references = references_names
    suffixes = [".h", ".cxx"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not (match := re.match(r"DECLARE(_[A-Z]+)*_TABLES?(_[A-Z]+)*\(", line)):
            return True
        # Extract names of the table type.
        line = remove_comment_cpp(line)
        line = line[len(match.group()) :].strip()  # Extract part after "(".
        if not (match := re.match(r"([^,\) ]+)", line)):
            print(f'Failed to extract table type from "{line}".')
            return False
        table_type_name = match.group(1)
        # print(f"Got \"{table_type_name}\"")
        # return True
        # Check for a version suffix.
        if match := re.match(r"(.*)_([0-9]{3})", line):
            table_type_name = match.group(1)
            # table_version = match.group(2)
            # print(f"Got versioned table \"{table_type_name}\", version {table_version}")
        if table_type_name[0] == "_":  # probably a macro variable
            return True
        if "#" in table_type_name:  # Remove "#" for strings in macros.
            table_type_name = table_type_name[: table_type_name.index("#")]
        # The actual test comes here.
        return is_upper_camel_case(table_type_name)


class TestNameNamespace(TestSpec):
    """Test names of namespaces."""

    name = "name/namespace"
    message = "Use snake_case for names of namespaces. Double underscores are not allowed."
    rationale = rationale_names
    references = references_names
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not line.startswith("namespace "):
            return True
        # Extract namespace name.
        namespace_name = line.split()[1]
        if namespace_name == "{":  # ignore anonymous namespaces
            return True
        # The actual test comes here.
        if "__" in namespace_name:
            return False
        return is_snake_case(namespace_name)


class TestNameType(TestSpec):
    """Test names of defined types."""

    name = "name/type"
    message = "Use UpperCamelCase for names of defined types (including concepts)."
    rationale = rationale_names
    references = references_names
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not (match := re.match(r"(using|concept) (\w+) = ", line)):
            return True
        # Extract type name.
        type_name = match.group(2)
        # The actual test comes here.
        return is_upper_camel_case(type_name)


class TestNameUpperCamelCase(TestSpec):
    """Base class for a test of UpperCamelCase names."""

    keyword = "key"
    name = f"name/{keyword}"
    message = f"Use UpperCamelCase for names of {keyword}."
    rationale = rationale_names
    references = references_names
    suffixes = [".h", ".cxx", ".C"]

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not (match := re.match(rf"{self.keyword}( (class|struct))? (\w+)", line)):
            return True
        # Extract object name.
        object_name = match.group(3)
        # The actual test comes here.
        return is_upper_camel_case(object_name)


class TestNameEnum(TestNameUpperCamelCase):
    """Test names of enumerators."""

    keyword = "enum"
    name = "name/enum"
    message = "Use UpperCamelCase for names of enumerators and their values."


class TestNameClass(TestNameUpperCamelCase):
    """Test names of classes."""

    keyword = "class"
    name = "name/class"
    message = "Use UpperCamelCase for names of classes."


class TestNameStruct(TestNameUpperCamelCase):
    """Test names of structs."""

    keyword = "struct"
    name = "name/struct"
    message = "Use UpperCamelCase for names of structs."


class TestNameFileCpp(TestSpec):
    """Test names of C++ files."""

    name = "name/file-cpp"
    message = "Use lowerCamelCase or UpperCamelCase for names of C++ files. See the O2 naming conventions for details."
    rationale = rationale_names
    references = references_names
    suffixes = [".h", ".cxx", ".C"]
    per_line = False

    def test_file(self, path: str, content) -> bool:
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
    rationale = rationale_names
    references = [Reference.LINTER_1, Reference.PY_PEP8]
    suffixes = [".py", ".ipynb"]
    per_line = False

    def test_file(self, path: str, content) -> bool:
        file_name = os.path.basename(path)
        return is_snake_case(file_name)


class TestNameWorkflow(TestSpec):
    """Test names of O2 workflows."""

    name = "name/o2-workflow"
    message = "Use kebab-case for names of workflows and match the name of the workflow file."
    rationale = f"{rationale_names} Correspondence workflow ↔ file."
    references = references_names
    suffixes = ["CMakeLists.txt"]
    per_line = False

    def test_file(self, path: str, content) -> bool:
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
                self.print_error(path, i + 1, f"Invalid workflow name: {workflow_name}.")
                continue
            # Extract workflow file name.
            next_line = content[i + 1].strip()
            words = next_line.split()
            if words[0] != "SOURCES":
                passed = False
                self.print_error(path, i + 2, f"Did not find sources for workflow: {workflow_name}.")
                continue
            workflow_file_name = os.path.basename(words[1])  # the actual file name
            # Generate the file name matching the workflow name.
            expected_workflow_file_name = kebab_case_to_camel_case_l(workflow_name) + ".cxx"
            # Compare the actual and expected file names.
            if expected_workflow_file_name != workflow_file_name:
                passed = False
                self.print_error(
                    path,
                    i + 1,
                    f"Workflow name {workflow_name} does not match its file name {workflow_file_name}. "
                    f"(Matches {expected_workflow_file_name}.)",
                )
        return passed


class TestNameTask(TestSpec):
    """Test explicit task names.
    Detect usage of TaskName, check whether it is needed and whether the task name matches the struct name."""

    name = "name/o2-task"
    message = "Specify task name only when it cannot be derived from the struct name. Only append to the default name."
    rationale = f"{rationale_names} Correspondence struct ↔ device."
    references = [Reference.LINTER_1]
    suffixes = [".cxx"]
    per_line = False

    def test_file(self, path: str, content) -> bool:
        is_inside_define = False  # Are we inside defineDataProcessing?
        is_inside_adapt = False  # Are we inside adaptAnalysisTask?
        struct_name = ""
        struct_templated = False  # Is the struct templated?
        n_parens_opened = 0  # number of opened parentheses
        passed = True
        for i, line in enumerate(content):
            if not line.strip():
                continue
            if self.is_disabled(line):
                continue
            if is_comment_cpp(line):
                continue
            line = remove_comment_cpp(line)
            # Wait for defineDataProcessing.
            if not is_inside_define:
                if not re.match(r"((o2::)?framework::)?WorkflowSpec defineDataProcessing\(", line):
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
                line = line[(index + len("adaptAnalysisTask<")) :]
                # Extract struct name.
                if not (match := re.match(r"([^>]+)", line)):
                    self.print_error(path, i + 1, f'Failed to extract struct name from "{line}".')
                    return False
                struct_name = match.group(1)
                if (index := struct_name.find("<")) > -1:
                    struct_templated = True
                    # print(f"{i + 1}: Got templated struct name {struct_name}")
                    struct_name = struct_name[:index]
                # print(f"{i + 1}: Got struct name {struct_name}")
                line = line[(line.index(struct_name) + len(struct_name)) :]
            if is_inside_adapt:
                n_parens_opened += line.count("(") - line.count(")")
                # print(f"{i + 1}: {n_parens_opened} opened parens")
                if n_parens_opened <= 0:
                    # print(f"{i + 1}: Exiting adapt.")
                    is_inside_adapt = False
                # Find explicit task name.
                if "TaskName{" not in line:
                    continue
                passed = False
                # Extract explicit task name.
                if not (match := re.search(r"TaskName\{\"([^\}]+)\"\}", line)):
                    self.print_error(path, i + 1, f'Failed to extract explicit task name from "{line}".')
                    return False
                task_name = match.group(1)
                # print(f"{i + 1}: Got struct \"{struct_name}\" with task name \"{task_name}\".")
                # Test explicit task name.
                device_name_from_struct_name = camel_case_to_kebab_case(
                    struct_name
                )  # default device name, in absence of TaskName
                device_name_from_task_name = camel_case_to_kebab_case(
                    task_name
                )  # actual device name, generated from TaskName
                struct_name_from_device_name = kebab_case_to_camel_case_u(
                    device_name_from_task_name
                )  # struct name matching the TaskName
                if not is_kebab_case(device_name_from_task_name):
                    self.print_error(
                        path,
                        i + 1,
                        f"Specified task name {task_name} produces an invalid device name "
                        f"{device_name_from_task_name}.",
                    )
                    passed = False
                elif device_name_from_struct_name == device_name_from_task_name:
                    # If the task name results in the same device name as the struct name would,
                    # TaskName is redundant and should be removed.
                    self.print_error(
                        path,
                        i + 1,
                        f"Specified task name {task_name} and the struct name {struct_name} produce "
                        f"the same device name {device_name_from_struct_name}. TaskName is redundant.",
                    )
                    passed = False
                elif device_name_from_struct_name.replace("-", "") == device_name_from_task_name.replace("-", ""):
                    # If the device names generated from the task name and from the struct name differ in hyphenation,
                    # capitalisation of the struct name should be fixed and TaskName should be removed.
                    # (special cases: alice3-, -2prong)
                    self.print_error(
                        path,
                        i + 1,
                        f"Device names {device_name_from_task_name} and {device_name_from_struct_name} generated "
                        f"from the specified task name {task_name} and from the struct name {struct_name}, "
                        f"respectively, differ in hyphenation. Consider fixing capitalisation of the struct name "
                        f"to {struct_name_from_device_name} and removing TaskName.",
                    )
                    passed = False
                elif device_name_from_task_name.startswith(device_name_from_struct_name):
                    # If the device name generated from the task name is an extension of the device name generated
                    # from the struct name, accept it if the struct is templated. If the struct is not templated,
                    # extension is acceptable if adaptAnalysisTask is called multiple times for the same struct.
                    if not struct_templated:
                        self.print_error(
                            path,
                            i + 1,
                            f"Device name {device_name_from_task_name} from the specified task name "
                            f"{task_name} is an extension of the device name {device_name_from_struct_name} "
                            f"from the struct name {struct_name} but the struct is not templated. "
                            "Is it adapted multiple times?",
                        )
                        passed = False
                    # else:
                    #     self.print_error(path, i + 1, f"Device name {device_name_from_task_name} from
                    # the specified task name {task_name} is an extension of the device name
                    # {device_name_from_struct_name} from the struct name {struct_name} and the struct is templated.
                    # All good")
                else:
                    # Other cases should be rejected.
                    self.print_error(
                        path,
                        i + 1,
                        f"Specified task name {task_name} produces device name {device_name_from_task_name} "
                        f"which does not match the device name {device_name_from_struct_name} from "
                        f"the struct name {struct_name}. (Matching struct name {struct_name_from_device_name})",
                    )
                    passed = False
        return passed


class TestNameFileWorkflow(TestSpec):
    """Test names of workflow files."""

    name = "name/workflow-file"
    message = (
        "Name of a workflow file must match the name of the main struct in it (without the PWG prefix). "
        '(Class implementation files should be in "Core" directories.)'
    )
    rationale = f"{rationale_names} Correspondence file ↔ struct."
    references = [Reference.LINTER_1]
    suffixes = [".cxx"]
    per_line = False

    def file_matches(self, path: str) -> bool:
        return super().file_matches(path) and "/Core/" not in path

    def test_file(self, path: str, content) -> bool:
        file_name = os.path.basename(path)[:-4]  # file name without suffix
        base_struct_name = f"{file_name[0].upper()}{file_name[1:]}"  # expected base of struct names
        if match := re.search("PWG([A-Z]{2})/", path):
            name_pwg = match.group(1)
            prefix_pwg = name_pwg.capitalize()
            if name_pwg in ("HF"):
                base_struct_name = rf"{prefix_pwg}{base_struct_name}"  # mandatory PWG prefix
            else:
                base_struct_name = rf"({prefix_pwg})?{base_struct_name}"  # optional PWG prefix
        # print(f"For file {file_name} expecting to find {base_struct_name}.")
        struct_names = []  # actual struct names in the file
        for line in content:
            if self.is_disabled(line):
                return True
            if not line.startswith("struct "):
                continue
            # Extract struct name.
            words = line.split()
            if not words[1].isidentifier():  # "struct : ..."
                continue
            struct_name = words[1]
            struct_names.append(struct_name)
        # print(f"Found structs: {struct_names}.")
        return any(re.match(base_struct_name, struct_name) for struct_name in struct_names)


class TestNameConfigurable(TestSpec):
    """Test names of configurables."""

    name = "name/configurable"
    message = (
        "Use lowerCamelCase for names of configurables and use the same name "
        "for the struct member as for the JSON string. (Declare the type and names on the same line.)"
    )
    rationale = f"{rationale_names} Correspondence C++ code ↔ JSON."
    references = [Reference.O2, Reference.LINTER_1]
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return super().file_matches(path) and "Macros/" not in path

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not (match := re.match(r"((o2::)?framework::)?Configurable(\w+|<.+>) (\w+)( = )?{([^,{]+),", line)):
            return not re.match(r"((o2::)?framework::)?Configurable", line)
        # Extract Configurable name.
        name_cpp = match.group(4)  # nameCpp
        name_json = match.group(6)  # expecting "nameJson"
        if name_json[0] != '"':  # JSON name is not a literal string.
            return True
        name_json = name_json[1:-1]  # Strip away quotation marks.
        # The actual test comes here.
        return is_lower_camel_case(name_cpp) and name_cpp == name_json


# PWG-HF

references_hf = [Reference.LINTER_1, Reference.PWG_HF]


class TestHfNameStructClass(TestSpec):
    """PWGHF: Test names of structs and classes."""

    name = "pwghf/name/struct-class"
    message = 'Names of PWGHF structs and classes must start with "Hf".'
    rationale = f"{rationale_names} Correspondence device ↔ workflow."
    references = references_hf
    suffixes = [".h", ".cxx"]

    def file_matches(self, path: str) -> bool:
        return super().file_matches(path) and "PWGHF/" in path and "Macros/" not in path

    def test_line(self, line: str) -> bool:
        if is_comment_cpp(line):
            return True
        if not line.startswith(("struct ", "class ")):
            return True
        line = remove_comment_cpp(line)
        # Extract struct/class name.
        words = line.split()
        if not words[1].isalnum():  # "struct : ..."
            return True
        struct_name = words[1]
        # The actual test comes here.
        return struct_name.startswith("Hf")


class TestHfNameFileTask(TestSpec):
    """PWGHF: Test names of task workflow files."""

    name = "pwghf/name/task-file"
    message = 'Name of a PWGHF task workflow file must start with "task".'
    rationale = rationale_names
    references = references_hf
    suffixes = [".cxx"]
    per_line = False

    def file_matches(self, path: str) -> bool:
        return super().file_matches(path) and "PWGHF/" in path and "Macros/" not in path

    def test_file(self, path: str, content) -> bool:
        file_name = os.path.basename(path)
        return not ("/Tasks/" in path and not file_name.startswith("task"))


class TestHfStructMembers(TestSpec):
    """PWGHF: Test order of struct members.
    Caveat: Does not see Configurables in ConfigurableGroup."""

    name = "pwghf/struct-member-order"
    message = "Declare struct members in the conventional order. See the PWGHF coding guidelines."
    rationale = rationale_names
    references = references_hf
    suffixes = [".cxx"]
    per_line = False
    member_order = [
        "Spawns<",
        "Builds<",
        "Produces<",
        "Configurable<",
        "HfHelper ",
        "SliceCache ",
        "Service<",
        "using ",
        "Filter ",
        "Preslice<",
        "PresliceUnsorted<",
        "Partition<",
        "ConfigurableAxis ",
        "AxisSpec ",
        "HistogramRegistry ",
        "OutputObj<",
        "void init(",
        "void process",
    ]

    def file_matches(self, path: str) -> bool:
        return super().file_matches(path) and "PWGHF/" in path

    def test_file(self, path: str, content) -> bool:
        passed = True
        dic_structs: dict[str, dict] = {}
        struct_name = ""
        for i, line in enumerate(content):
            if is_comment_cpp(line):
                continue
            if line.startswith("struct "):  # expecting no indentation
                line = remove_comment_cpp(line)
                struct_name = line.strip().split()[1]
                dic_structs[struct_name] = {}
                continue
            if not struct_name:
                continue
            # Save line numbers of members of the current struct for each category.
            for member in self.member_order:
                if line.startswith(f"  {member}"):  # expecting single-level indentation for direct members
                    if member not in dic_structs[struct_name]:
                        dic_structs[struct_name][member] = []
                    dic_structs[struct_name][member].append(i + 1)  # save line number
                    break
        # print(dic_struct)
        # Detect members declared in a wrong order.
        last_line_last_member = 0  # line number of the last member of the previous member category
        index_last_member = 0  # index of the previous member category in the member_order list
        for struct_name, dic_struct in dic_structs.items():
            for i_m, member in enumerate(self.member_order):
                if member not in dic_struct:
                    continue
                first_line = min(dic_struct[member])  # line number of the first member of this category
                last_line = max(dic_struct[member])  # line number of the last member of this category
                if (
                    first_line < last_line_last_member
                ):  # The current category starts before the end of the previous category.
                    passed = False
                    self.print_error(
                        path,
                        first_line,
                        f"{struct_name}: {member.strip()} appears too early "
                        f"(before end of {self.member_order[index_last_member].strip()}).",
                    )
                last_line_last_member = last_line
                index_last_member = i_m
        return passed


# End of test implementations


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="O2 linter (Find O2-specific issues in O2 code)")
    parser.add_argument("paths", type=str, nargs="+", help="File path(s)")
    parser.add_argument(
        "-g",
        dest="github",
        action="store_true",
        help="Print messages also as GitHub annotations",
    )
    args = parser.parse_args()
    if args.github:
        global github_mode  # pylint: disable=global-statement  # noqa: PLW0603
        github_mode = True

    tests: list[TestSpec] = []  # list of activated tests

    # Bad practice
    enable_bad_practice = True
    if enable_bad_practice:
        tests.append(TestIoStream())
        tests.append(TestUsingStd())
        tests.append(TestUsingDirective())
        tests.append(TestStdPrefix())
        tests.append(TestRootEntity())
        tests.append(TestRootLorentzVector())
        tests.append(TestPi())
        tests.append(TestTwoPiAddSubtract())
        tests.append(TestPiMultipleFraction())
        tests.append(TestPdgDatabase())
        tests.append(TestPdgExplicitCode())
        tests.append(TestPdgExplicitMass())
        tests.append(TestPdgKnownMass())
        tests.append(TestLogging())
        tests.append(TestConstRefInForLoop())
        tests.append(TestConstRefInSubscription())
        tests.append(TestWorkflowOptions())
        tests.append(TestMagicNumber())

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
        tests.append(TestNameType())
        tests.append(TestNameEnum())
        tests.append(TestNameClass())
        tests.append(TestNameStruct())
        tests.append(TestNameFileCpp())
        tests.append(TestNameFilePython())
        tests.append(TestNameWorkflow())
        tests.append(TestNameTask())
        tests.append(TestNameFileWorkflow())
        tests.append(TestNameConfigurable())

    # PWG-HF
    enable_pwghf = True
    if enable_pwghf:
        tests.append(TestHfNameStructClass())
        tests.append(TestHfNameFileTask())
        tests.append(TestHfStructMembers())

    test_names = [t.name for t in tests]  # short names of activated tests
    suffixes = tuple({s for test in tests for s in test.suffixes})  # all suffixes from all enabled tests
    passed = True  # global result of all tests
    n_files_bad = dict.fromkeys(test_names, 0)  # counter of files with issues

    # Report overview before running.
    print(f"Testing {len(args.paths)} files.")
    # print(args.paths)
    print("Enabled tests:", test_names)
    print("Suffixes of tested files:", sorted(suffixes))
    # print(f"Github annotations: {github_mode}.")

    # Test files.
    for path in args.paths:
        # print(f"Processing path \"{path}\".")
        # Skip not tested files.
        if not path.endswith(suffixes):
            # print(f"Skipping path \"{path}\".")
            continue
        try:
            with open(path, encoding="utf-8") as file:
                tolerated_tests = get_tolerated_tests(path)
                content = file.readlines()
                for test in tests:
                    test.tolerated = test.name in tolerated_tests
                    result = test.run(path, content)
                    if not result:
                        n_files_bad[test.name] += 1
                        passed = False
                    # print(f"File \"{path}\" {'passed' if result else 'failed'} the test {test.name}.")
        except OSError:
            print(f'Failed to open file "{path}".')
            sys.exit(1)

    # Report results for tests that failed or were disabled or were tolerated.
    n_issues, n_disabled, n_tolerated = 0, 0, 0  # global counters
    if not passed or any(n > 0 for n in (test.n_disabled + test.n_tolerated for test in tests)):
        print("\nResults for failed, tolerated and disabled tests")
        len_max = max(len(name) for name in test_names)
        print(f"test{' ' * (len_max - len('test'))}\tissues\ttolerated\tdisabled\tbad files\trationale")
        print("-" * len_max)
        ref_names = []
        for test in tests:
            if any(n > 0 for n in (test.n_issues, test.n_disabled, test.n_tolerated, n_files_bad[test.name])):
                ref_ids = sorted(ref.value for ref in test.references)
                ref_names += test.references
                print(
                    f"{test.name}{' ' * (len_max - len(test.name))}\t{test.n_issues}\t{test.n_tolerated}"
                    f"\t\t{test.n_disabled}\t\t{n_files_bad[test.name]}\t\t{test.rationale} {ref_ids}"
                )
            n_issues += test.n_issues
            n_disabled += test.n_disabled
            n_tolerated += test.n_tolerated
        print("-" * len_max)
        # Print the totals.
        name_total = "total"
        print(f"{name_total}{' ' * (len_max - len(name_total))}\t{n_issues}\t{n_tolerated}\t\t{n_disabled}")
        # Print list of references for listed tests.
        print("\nReferences")
        ref_names = list(dict.fromkeys(ref_names))
        for ref_name, data in references.items():
            if ref_name in ref_names:
                print(f"[{ref_name.value}]\t{data['title']}. <{data['url']}>.")

    # Report global result.
    title_result = "O2 linter result"
    if passed:
        msg_result = "All tests passed."
        if github_mode:
            print(f"\n::notice title={title_result}::{msg_result}")
        else:
            print(f"\n{title_result}: {msg_result}")
    else:
        msg_result = "Issues have been found."
        msg_disable = (
            f'Exceptionally, you can disable a test for a line by adding a comment with "{prefix_disable}"'
            " followed by the name of the test and parentheses with a reason for the exception."
        )
        msg_tolerate = f'To tolerate certain issues in a directory, add a line with the test name in "{file_config}".'
        if github_mode:
            print(f"\n::error title={title_result}::{msg_result}")
            print(f"::notice::{msg_disable}")
            print(f"::notice::{msg_tolerate}")
        else:
            print(f"\n{title_result}: {msg_result}")
            print(msg_disable)
            print(msg_tolerate)

    # Make results available to the GitHub actions.
    if github_mode:
        try:
            with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as fh:
                print(f"n_issues={n_issues}", file=fh)
                print(f"n_disabled={n_disabled}", file=fh)
                print(f"n_tolerated={n_tolerated}", file=fh)
        except KeyError:
            print("Skipping writing in GITHUB_OUTPUT.")

    # Print tips.
    print(
        "\nTip: You can run the O2 linter locally from the O2Physics directory with: python3 Scripts/o2_linter.py <files>"
    )

    if not passed:
        sys.exit(1)


if __name__ == "__main__":
    main()
