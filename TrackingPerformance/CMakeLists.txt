#[[
Copyright (c) 2020-2024 Key4hep-Project.

This file is part of Key4hep.
See https://key4hep.github.io/key4hep-doc/ for further info.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
]]

set(PackageName TrackingPerformance)

project(${PackageName})

include(ExternalData)
list(APPEND ExternalData_URL_TEMPLATES
  "https://key4hep.web.cern.ch:443/testFiles/k4RecTracker/%(hash)"
)

file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cc
  ${PROJECT_SOURCE_DIR}/src/*.cpp
  ${PROJECT_SOURCE_DIR}/src/components/*.cpp
)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

find_package(ROOT REQUIRED COMPONENTS Core RIO Tree MathCore MathMore Graf Graf3d Hist)

gaudi_add_module(${PackageName}
  SOURCES ${sources}
  LINK
    Gaudi::GaudiKernel
    EDM4HEP::edm4hep
    k4FWCore::k4FWCore
    ROOT::MathCore
    ROOT::MathMore
    ROOT::Physics
    ROOT::Graf
    ROOT::Graf3d
    ROOT::Hist
)

target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
 )

set_target_properties(${PackageName} PROPERTIES
  PUBLIC_HEADER "${headers}"
)

file(GLOB test_python_scripts
  ${PROJECT_SOURCE_DIR}/test/*.py
)

set(test_shell_scripts
  ${PROJECT_SOURCE_DIR}/test/testTrackingValidation.sh
)

set(test_scripts
  ${test_python_scripts}
  ${test_shell_scripts}
)

file(COPY ${test_scripts} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test)

install(TARGETS ${PackageName}
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}" COMPONENT dev
)
message(STATUS "CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")
install(FILES ${test_scripts} DESTINATION test)

SET(test_name "testTrackingValidation")

ExternalData_Add_Test(testTrackingValidation
  NAME ${test_name}
  COMMAND sh +x TrackingPerformance/test/testTrackingValidation.sh
          DATA{${CMAKE_SOURCE_DIR}/TrackingPerformance/test/inputFiles/SimpleGatrIDEAv3o1.onnx}
)

set_test_env(${test_name})
set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
ExternalData_Add_Target(${test_name})

