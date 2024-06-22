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

include_guard()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-error")

# Enabled warnings supported by Clang and GCC
set(O2PHYSICS_WARNINGS_COMMON "NSObject-attribute;address;address-of-packed-member;attribute-warning;attributes;bool-operation;builtin-macro-redefined;char-subscripts;class-conversion;comment;conversion-null;cpp;dangling-else;delete-incomplete;delete-non-virtual-dtor;deprecated;deprecated-copy;deprecated-declarations;div-by-zero;empty-body;endif-labels;enum-compare;exceptions;format-extra-args;format-zero-length;frame-address;free-nonheap-object;ignored-attributes;ignored-qualifiers;inaccessible-base;infinite-recursion;int-in-bool-context;int-to-pointer-cast;invalid-offsetof;logical-not-parentheses;main;memset-transposed-args;misleading-indentation;mismatched-new-delete;missing-field-initializers;narrowing;noexcept-type;nonnull;odr;parentheses;pessimizing-move;pointer-compare;pragmas;psabi;range-loop-construct;redundant-move;reorder;return-type;sequence-point;shift-count-negative;shift-count-overflow;sizeof-array-argument;sizeof-array-div;sizeof-pointer-div;sizeof-pointer-memaccess;string-compare;switch;switch-bool;tautological-compare;trigraphs;uninitialized;unknown-pragmas;unused;unused-but-set-parameter;unused-but-set-variable;unused-function;unused-label;unused-parameter;unused-result;unused-value;unused-variable;varargs;vexing-parse")
# TODO: sign-compare

# Enabled warnings supported by Clang only
set(O2PHYSICS_WARNINGS_CLANG "#pragma-messages;#warnings;CFString-literal;CL4;IndependentClass-attribute;absolute-value;abstract-final-class;address-of-temporary;aix-compat;align-mismatch;alloca-with-align-alignof;always-inline-coroutine;ambiguous-delete;ambiguous-ellipsis;ambiguous-macro;ambiguous-member-template;ambiguous-reversed-operator;analyzer-incompatible-plugin;anon-enum-enum-conversion;anonymous-pack-parens;arc;arc-bridge-casts-disallowed-in-nonarc;arc-non-pod-memaccess;arc-performSelector-leaks;arc-retain-cycles;arc-unsafe-retained-assign;argument-outside-range;argument-undefined-behaviour;array-bounds;array-parameter;asm;asm-operand-widths;assume;atimport-in-framework-header;atomic-access;atomic-alignment;atomic-memory-ordering;atomic-property-with-user-defined-accessor;attribute-packed-for-bitfield;auto-disable-vptr-sanitizer;auto-storage-class;auto-var-id;availability;avr-rtlib-linking-quirks;backend-plugin;backslash-newline-escape;binding-in-condition;bitfield-constant-conversion;bitfield-width;bitwise-conditional-parentheses;bitwise-instead-of-logical;bitwise-op-parentheses;block-capture-autoreleasing;bool-conversion;bool-conversions;braced-scalar-init;branch-protection;bridge-cast;builtin-assume-aligned-alignment;builtin-memcpy-chk-size;builtin-requires-header;call-to-pure-virtual-from-ctor-dtor;called-once-parameter;cast-calling-convention;cast-of-sel-type;cast-qual-unrelated;clang-cl-pch;class-varargs;cmse-union-leak;comments;compare-distinct-pointer-types;compound-token-split;compound-token-split-by-macro;concepts-ts-compat;conditional-type-mismatch;config-macros;constant-conversion;constant-evaluated;constant-logical-operand;constexpr-not-const;conversion;coroutine;coroutine-missing-unhandled-exception;ctu;cuda-compat;cxx-attribute-extension;dangling;dangling-field;dangling-gsl;dangling-initializer-list;darwin-sdk-settings;dealloc-in-category;debug-compression-unavailable;defaulted-function-deleted;delegating-ctor-cycles;delete-abstract-non-virtual-dtor;delete-non-abstract-non-virtual-dtor;deprecate-lax-vec-conv-all;deprecated-altivec-src-compat;deprecated-anon-enum-enum-conversion;deprecated-array-compare;deprecated-attributes;deprecated-builtins;deprecated-comma-subscript;deprecated-copy-with-user-provided-copy;deprecated-coroutine;deprecated-enum-compare;deprecated-enum-compare-conditional;deprecated-enum-enum-conversion;deprecated-enum-float-conversion;deprecated-experimental-coroutine;deprecated-increment-bool;deprecated-non-prototype;deprecated-objc-isa-usage;deprecated-objc-pointer-introspection;deprecated-objc-pointer-introspection-performSelector;deprecated-pragma;deprecated-register;deprecated-static-analyzer-flag;deprecated-type;deprecated-volatile;deprecated-writable-strings;distributed-object-modifiers;division-by-zero;dll-attribute-on-redeclaration;dllexport-explicit-instantiation-decl;dllimport-static-field-def;dtor-name;dtor-typedef;duplicate-decl-specifier;duplicate-protocol;dynamic-class-memaccess;dynamic-exception-spec;elaborated-enum-base;elaborated-enum-class;empty-decomposition;empty-init-stmt;encode-type;enum-compare-conditional;enum-compare-switch;enum-conversion;enum-enum-conversion;enum-float-conversion;enum-too-large;excess-initializers;expansion-to-defined;explicit-initialize-call;export-unnamed;export-using-directive;extern-c-compat;extern-initializer;extra-qualification;extra-tokens;final-dtor-non-final-class;final-macro;fixed-point-overflow;flag-enum;for-loop-analysis;format;format-insufficient-args;format-invalid-specifier;format-security;format=2;fortify-source;frame-larger-than;frame-larger-than=;framework-include-private-from-public;function-def-in-objc-container;function-multiversion;fuse-ld-path;gcc-compat;global-isel;gnu;gnu-alignof-expression;gnu-array-member-paren-init;gnu-designator;gnu-folding-constant;gnu-inline-cpp-without-extern;gnu-null-pointer-arithmetic;gnu-static-float-init;gnu-string-literal-operator-template;gnu-variable-sized-type-not-at-end;gpu-maybe-wrong-side;header-guard;hip-only;hlsl-extensions;ignored-availability-without-sdk-settings;ignored-optimization-argument;ignored-pragma-intrinsic;ignored-pragmas;ignored-reference-qualifiers;implicit;implicit-const-int-float-conversion;implicit-conversion-floating-point-to-bool;implicit-exception-spec-mismatch;implicit-fixed-point-conversion;implicit-float-conversion;implicit-function-declaration;implicit-int;implicit-int-float-conversion;implicitly-unsigned-literal;include-next-absolute-path;include-next-outside-header;incompatible-exception-spec;incompatible-function-pointer-types;incompatible-library-redeclaration;incompatible-ms-struct;incompatible-pointer-types;incompatible-pointer-types-discards-qualifiers;incompatible-property-type;incompatible-sysroot;incomplete-framework-module-declaration;incomplete-implementation;incomplete-module;incomplete-setjmp-declaration;incomplete-umbrella;inconsistent-dllimport;inconsistent-missing-override;increment-bool;initializer-overrides;injected-class-name;inline-asm;inline-namespace-reopened-noninline;inline-new-delete;instantiation-after-specialization;int-conversion;int-conversions;int-to-void-pointer-cast;integer-overflow;interrupt-service-routine;invalid-command-line-argument;invalid-constexpr;invalid-iboutlet;invalid-initializer-from-system-header;invalid-ios-deployment-target;invalid-no-builtin-names;invalid-noreturn;invalid-or-nonexistent-directory;invalid-partial-specialization;invalid-pp-token;invalid-source-encoding;invalid-token-paste;jump-seh-finally;keyword-compat;knr-promoted-parameter;large-by-value-copy;linker-warnings;literal-conversion;literal-range;local-type-template-args;logical-op-parentheses;macro-redefined;main-return-type;malformed-warning-check;many-braces-around-scalar-init;max-unsigned-zero;memsize-comparison;microsoft;microsoft-abstract;microsoft-anon-tag;microsoft-cast;microsoft-const-init;microsoft-default-arg-redefinition;microsoft-drectve-section;microsoft-enum-forward-reference;microsoft-exception-spec;microsoft-exists;microsoft-explicit-constructor-call;microsoft-extra-qualification;microsoft-goto;microsoft-inaccessible-base;microsoft-include;microsoft-mutable-reference;microsoft-pure-definition;microsoft-sealed;microsoft-static-assert;microsoft-template;microsoft-template-shadow;microsoft-union-member-reference;microsoft-unqualified-friend;microsoft-using-decl;microsoft-void-pseudo-dtor;misexpect;mismatched-parameter-types;mismatched-return-types;mismatched-tags;missing-braces;missing-constinit;missing-declarations;missing-exception-spec;missing-method-return-type;missing-noescape;missing-prototype-for-cc;missing-selector-name;missing-sysroot;misspelled-assumption;module-conflict;module-file-config-mismatch;module-file-extension;module-import-in-extern-c;modules-ambiguous-internal-linkage;modules-import-nested-redundant;move;msvc-include;msvc-not-found;multichar;multiple-move-vbase;new-returns-null;noderef;non-c-typedef-for-linkage;non-gcc;non-literal-null-conversion;non-pod-varargs;non-power-of-two-alignment;nonportable-include-path;nonportable-vector-initialization;nontrivial-memaccess;nsconsumed-mismatch;nsreturns-mismatch;null-arithmetic;null-character;null-conversion;null-dereference;null-pointer-arithmetic;null-pointer-subtraction;nullability;nullability-completeness;nullability-completeness-on-arrays;nullability-declspec;nullability-inferred-on-nested-type;objc-autosynthesis-property-ivar-name-match;objc-bool-constant-conversion;objc-boxing;objc-circular-container;objc-cocoa-api;objc-designated-initializers;objc-dictionary-duplicate-keys;objc-flexible-array;objc-forward-class-redefinition;objc-literal-compare;objc-literal-conversion;objc-macro-redefinition;objc-method-access;objc-missing-super-calls;objc-multiple-method-names;objc-noncopy-retain-block-property;objc-nonunified-exceptions;objc-property-implementation;objc-property-implicit-mismatch;objc-property-matches-cocoa-ownership-rule;objc-property-no-attribute;objc-property-synthesis;objc-protocol-method-implementation;objc-protocol-property-synthesis;objc-protocol-qualifiers;objc-readonly-with-setter-property;objc-redundant-api-use;objc-redundant-literal-use;objc-root-class;objc-signed-char-bool;objc-signed-char-bool-implicit-float-conversion;objc-string-compare;objc-string-concatenation;objc-unsafe-perform-selector;opencl-unsupported-rgba;openmp;openmp-51-extensions;openmp-clauses;openmp-loop-form;openmp-mapping;openmp-target;option-ignored;ordered-compare-function-pointers;out-of-line-declaration;out-of-scope-function;overloaded-shift-op-parentheses;overloaded-virtual;override-init;override-module;overriding-t-option;parentheses-equality;partial-availability;pass-failed;pch-date-time;pedantic-macros;pointer-arith;pointer-bool-conversion;pointer-integer-compare;pointer-sign;pointer-to-enum-cast;pointer-to-int-cast;pointer-type-mismatch;potentially-direct-selector;potentially-evaluated-expression;pragma-clang-attribute;pragma-once-outside-header;pragma-pack;pragma-system-header-outside-header;predefined-identifier-outside-function;private-extern;private-header;private-module;profile-instr-out-of-date;profile-instr-unprofiled;property-access-dot-syntax;property-attribute-mismatch;protocol;protocol-property-synthesis-ambiguity;qualified-void-return-type;readonly-iboutlet-property;receiver-expr;receiver-forward-class;redeclared-class-member;redundant-consteval-if;register;reinterpret-base-class;reorder-ctor;reorder-init-list;requires-super-attribute;reserved-user-defined-literal;restrict-expansion;return-stack-address;return-type-c-linkage;rewrite-not-bool;rtti;sarif-format-unstable;section;self-assign;self-assign-field;self-assign-overloaded;self-move;semicolon-before-method-body;sentinel;serialized-diagnostics;shadow;shadow-all;shadow-ivar;shift-negative-value;shift-op-parentheses;shift-overflow;signed-unsigned-wchar;sizeof-array-decay;slash-u-filename;slh-asm-goto;sometimes-uninitialized;source-mgr;source-uses-openmp;stack-exhausted;stack-protector;static-float-init;static-in-inline;static-inline-explicit-instantiation;static-local-in-inline;static-self-init;stdlibcxx-not-found;strict-potentially-direct-selector;strict-prototypes;string-concatenation;string-plus-char;string-plus-int;strlcpy-strlcat-size;strncat-size;suspicious-bzero;suspicious-memaccess;swift-name-attribute;sync-fetch-and-nand-semantics-changed;target-clones-mixed-specifiers;tautological-bitwise-compare;tautological-constant-compare;tautological-constant-out-of-range-compare;tautological-objc-bool-compare;tautological-overlap-compare;tautological-pointer-compare;tautological-undefined-compare;tcb-enforcement;tentative-definition-incomplete-type;type-safety;typedef-redefinition;typename-missing;unable-to-open-stats-file;unaligned-qualifier-implicit-cast;unavailable-declarations;undefined-bool-conversion;undefined-inline;undefined-internal;undefined-var-template;underaligned-exception-object;unevaluated-expression;unguarded-availability;unguarded-availability-new;unicode;unicode-homoglyph;unicode-whitespace;unicode-zero-width;uninitialized-const-reference;unknown-argument;unknown-assumption;unknown-attributes;unknown-cuda-version;unknown-directives;unknown-escape-sequence;unknown-sanitizers;unknown-warning-option;unnamed-type-template-args;unneeded-internal-declaration;unqualified-std-cast-call;unreachable-code;unreachable-code-aggressive;unreachable-code-generic-assoc;unsequenced;unsupported-abi;unsupported-abs;unsupported-availability-guard;unsupported-cb;unsupported-floating-point-opt;unsupported-friend;unsupported-gpopt;unsupported-nan;unsupported-target-opt;unsupported-visibility;unusable-partial-specialization;unused-command-line-argument;unused-comparison;unused-const-variable;unused-getter-return-value;unused-lambda-capture;unused-local-typedef;unused-private-field;unused-property-ivar;unused-volatile-lvalue;user-defined-literals;user-defined-warnings;variadic-macros;vec-elem-size;visibility;void-pointer-to-enum-cast;void-pointer-to-int-cast;void-ptr-dereference;wasm-exception-spec;writable-strings;write-strings;xor-used-as-pow")

# Enabled warnings supported by GCC only
set(O2PHYSICS_WARNINGS_GCC "aggressive-loop-optimizations;analyzer-double-fclose;analyzer-double-free;analyzer-exposure-through-output-file;analyzer-file-leak;analyzer-free-of-non-heap;analyzer-malloc-leak;analyzer-mismatching-deallocation;analyzer-null-argument;analyzer-null-dereference;analyzer-possible-null-argument;analyzer-possible-null-dereference;analyzer-shift-count-negative;analyzer-shift-count-overflow;analyzer-stale-setjmp-buffer;analyzer-tainted-allocation-size;analyzer-tainted-array-index;analyzer-tainted-divisor;analyzer-tainted-offset;analyzer-tainted-size;analyzer-unsafe-call-within-signal-handler;analyzer-use-after-free;analyzer-use-of-pointer-in-stale-stack-frame;analyzer-use-of-uninitialized-value;analyzer-write-to-const;analyzer-write-to-string-literal;array-bounds=1;array-compare;array-parameter=2;bool-compare;builtin-declaration-mismatch;cannot-profile;cast-function-type;catch-value=1;class-memaccess;clobbered;coverage-invalid-line-number;coverage-mismatch;dangling-pointer=2;format-contains-nul;format-diag;format-overflow=1;format-truncation=1;format=1;if-not-aligned;implicit-fallthrough=3;inherited-variadic-ctor;init-list-lifetime;init-self;interference-size;invalid-memory-model;literal-suffix;lto-type-mismatch;maybe-uninitialized;memset-elt-size;mismatched-dealloc;missing-attributes;missing-profile;missing-requires;missing-template-keyword;multistatement-macros;non-template-friend;nonnull-compare;openmp-simd;overflow;packed-bitfield-compat;packed-not-aligned;pmf-conversions;prio-ctor-dtor;restrict;return-local-addr;scalar-storage-order;sized-deallocation;strict-aliasing=3;strict-overflow=1;stringop-overflow=2;stringop-overread;stringop-truncation;subobject-linkage;switch-outside-range;switch-unreachable;sync-nand;terminate;tsan;type-limits;unused-local-typedefs;use-after-free=2;virtual-move-assign;vla-parameter;volatile-register-var;zero-length-bounds")

# Function to build a list of warning flags from their names
function(o2_build_warning_flags)
  cmake_parse_arguments(PARSE_ARGV 0 A "" "PREFIX;OUTPUTVARNAME" "WARNINGS")
  if(A_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Unexpected unparsed arguments: ${A_UNPARSED_ARGUMENTS}")
  endif()
  list(TRANSFORM A_WARNINGS STRIP)
  list(TRANSFORM A_WARNINGS PREPEND ${A_PREFIX})
  string(JOIN " " OUTPUT ${A_WARNINGS})
  set(${A_OUTPUTVARNAME} ${OUTPUT} PARENT_SCOPE)
endfunction()

message(STATUS "O2PHYSICS_WARNINGS_AS_ERRORS: ${O2PHYSICS_WARNINGS_AS_ERRORS}")

# Treat warnings as errors.
if(O2PHYSICS_WARNINGS_AS_ERRORS)
  # Add warnings for all platforms.
  o2_build_warning_flags(PREFIX "-Werror="
                         OUTPUTVARNAME O2PHYSICS_CXX_WARNINGS_COMMON
                         WARNINGS ${O2PHYSICS_WARNINGS_COMMON})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${O2PHYSICS_CXX_WARNINGS_COMMON}")
  if(APPLE)
    # Add warnings for macOS only.
    o2_build_warning_flags(PREFIX "-Werror="
                           OUTPUTVARNAME O2PHYSICS_CXX_WARNINGS_APPLE
                           WARNINGS ${O2PHYSICS_WARNINGS_CLANG})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${O2PHYSICS_CXX_WARNINGS_APPLE}")
  elseif(UNIX)
    # Add warnings for Linux only.
    o2_build_warning_flags(PREFIX "-Werror="
                           OUTPUTVARNAME O2PHYSICS_CXX_WARNINGS_UNIX
                           WARNINGS ${O2PHYSICS_WARNINGS_GCC})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${O2PHYSICS_CXX_WARNINGS_UNIX}")
  endif()
endif()

IF (ENABLE_TIMETRACE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftime-trace")
ENDIF()

set(CMAKE_CXX_FLAGS_COVERAGE "-g -O2 -fprofile-arcs -ftest-coverage")
set(CMAKE_C_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_COVERAGE}")
set(CMAKE_Fortran_FLAGS_COVERAGE "-g -O2 -fprofile-arcs -ftest-coverage")
set(CMAKE_LINK_FLAGS_COVERAGE "--coverage -fprofile-arcs  -fPIC")

MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_COVERAGE
    CMAKE_C_FLAGS_COVERAGE
    CMAKE_Fortran_FLAGS_COVERAGE
    CMAKE_LINK_FLAGS_COVERAGE)

#Check the compiler and set the compile and link flags
IF (NOT CMAKE_BUILD_TYPE)
  Message(STATUS "Set BuildType to DEBUG")
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
ENDIF (NOT CMAKE_BUILD_TYPE)

IF(ENABLE_CASSERT) #For the CI, we want to have <cassert> assertions enabled
    set(CMAKE_CXX_FLAGS_RELEASE "-O2")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
ELSE()
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
    if (CMAKE_BUILD_TYPE STREQUAL "RELEASE" OR CMAKE_BUILD_TYPE STREQUAL "RELWITHDEBINFO")
      set(FAIR_MIN_SEVERITY "info")
    endif()
ENDIF()
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g")
# make sure Debug build not optimized (does not seem to work without CACHE + FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "Debug mode build flags" FORCE)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "Debug mode build flags" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0" CACHE STRING "Debug mode build flags" FORCE)

if(APPLE)
elseif(UNIX)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined") # avoid undefined in our libs
endif()

message(STATUS "Using build type: ${CMAKE_BUILD_TYPE} - CXXFLAGS: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
