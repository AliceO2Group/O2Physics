include_guard()

find_package(ONNXRuntime::ONNXRuntime CONFIG)

if (NOT ONNXRuntime::ONNXRuntime_FOUND)
  message(WARNING "No config found for ONNXRuntime::ONNXRuntime. Trying pkgconfig version")
  find_package(PkgConfig)
  pkg_check_modules(ONNXRuntime REQUIRED IMPORTED_TARGET GLOBAL libonnxruntime)

  if (TARGET PkgConfig::ONNXRuntime)
     message(WARNING "Creating alias ONNXRuntime::ONNXRuntime")
     add_library(ONNXRuntime::ONNXRuntime ALIAS PkgConfig::ONNXRuntime)
  endif()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(ONNXRuntime::ONNXRuntime
                                    REQUIRED_VARS ONNXRuntime_LINK_LIBRARIES)
endif()
