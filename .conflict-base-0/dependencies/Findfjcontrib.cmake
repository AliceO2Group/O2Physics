
if (TARGET FastJet::Contrib)
   # nothing to do, target is already there
   set(fjcontrib_FOUND TRUE)
   return()
endif()

find_library(fjcontrib_LIBRARY_SHARED NAMES libfastjetcontribfragile.dylib libfastjetcontribfragile.so PATHS ${fjcontrib_ROOT}/lib)

if (fjcontrib_LIBRARY_SHARED)
      add_library(fjcontrib SHARED IMPORTED)
      get_filename_component(libdir ${fjcontrib_LIBRARY_SHARED} DIRECTORY)
      get_filename_component(incdir ${libdir}/../include ABSOLUTE)
      set_target_properties(fjcontrib PROPERTIES IMPORTED_LOCATION
                            ${fjcontrib_LIBRARY_SHARED}
                            INTERFACE_LINK_DIRECTORIES ${libdir}
                            INTERFACE_INCLUDE_DIRECTORIES ${incdir})
      set_target_properties(fjcontrib PROPERTIES IMPORTED_GLOBAL TRUE)
      add_library(FastJet::Contrib ALIAS fjcontrib)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(fjcontrib REQUIRED_VARS fjcontrib_LIBRARY_SHARED)

mark_as_advanced(fjcontrib_LIBRARY_SHARED)

