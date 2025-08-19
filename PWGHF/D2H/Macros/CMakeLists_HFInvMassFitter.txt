cmake_minimum_required(VERSION 3.10)

project(HFInvMassFitter)

set(HFFITTER_RAPIDJSON_INCLUDE_DIRS "" CACHE STRING "Location of rapidjson include directories")

find_package(ROOT REQUIRED COMPONENTS RooFit RooFitCore)

include_directories(${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR} ${HFFITTER_RAPIDJSON_INCLUDE_DIRS})

set(SOURCES
  HFInvMassFitter.cxx
  runMassFitter.C
)

add_executable(runMassFitter ${SOURCES} "HFInvMassFitter.h")

ROOT_GENERATE_DICTIONARY(G__HFInvMassFitter
  HFInvMassFitter.h LINKDEF HFInvMassFitterLinkDef.h
  MODULE runMassFitter
)

target_link_libraries(runMassFitter PRIVATE ${ROOT_LIBRARIES} ROOT::EG ROOT::RooFit ROOT::RooFitCore)
