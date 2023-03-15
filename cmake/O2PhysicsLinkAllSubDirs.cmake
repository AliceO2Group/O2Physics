function(o2physics_link_all_subdirs mainDir currentDir installDir)
  # Find the relative directory
  string(REPLACE "${mainDir}" "" relDir "${currentDir}")
  # Splite the relative directory to a list of subdirectories
  string(REPLACE "/" ";" listDirs "${relDir}")
  foreach(dir ${listDirs})
    # Link each subdirectory inside include
    # This helps to use #include "A/B/C/x.h" from the flat installation of header files
    message("Linking ${dir} in ${installDir}")
    execute_process(
      COMMAND ln -sfh . ${dir}
      WORKING_DIRECTORY ${installDir}
      )
  endforeach()
endfunction()
