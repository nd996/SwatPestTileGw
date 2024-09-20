file(GLOB sources *.f*)
list(FILTER sources EXCLUDE REGEX ".*modparm\\.f$")

add_executable(swat_681 main.f ${sources})
add_custom_command(TARGET swat_681 
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:swat_681> ../../)
