set(proj cpprestsdk)

if(MITK_USE_${proj})
  set(${proj}_DEPENDS ${proj} Boost zlib)

  if(DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR})
    message(FATAL_ERROR "${proj}_DIR variable is defined but corresponds to non-existing directory!")
  endif()

  if(NOT DEFINED ${proj}_DIR)
    ExternalProject_Add(${proj}
      GIT_REPOSITORY https://github.com/Microsoft/cpprestsdk.git
      GIT_TAG v2.10.10
      SOURCE_SUBDIR Release
      CMAKE_ARGS ${ep_common_args}
      CMAKE_CACHE_ARGS ${ep_common_cache_args}
        -DBUILD_SAMPLES:BOOL=OFF
        -DBUILD_TESTS:BOOL=OFF
      CMAKE_CACHE_DEFAULT_ARGS ${ep_common_cache_default_args}
    )

    ExternalProject_Add_Step(${proj} git-submodule
      COMMAND ${GIT_EXECUTABLE} submodule update --init -- Release/libs/websocketpp
      WORKING_DIRECTORY <SOURCE_DIR>
      DEPENDEES patch
      DEPENDERS configure
    )

    set(${proj}_DIR ${ep_prefix})
  else()
    mitkMacroEmptyExternalProject(${proj} "${proj_DEPENDENCIES}")
  endif()
endif()
