file(GLOB SRCS *.cpp)
file(GLOB HDRS *.h)

include(GoogleTest)
add_executable(testAlgorithms ${SRCS} ${HDRS})

target_link_libraries(testAlgorithms
    Analysis
    gtest
   # libgmock
)
set_property(TARGET testAlgorithms PROPERTY FOLDER "Tests")
gtest_discover_tests(testAlgorithms
 EXTRA_ARGS --gtest_color=yes)

install(TARGETS testAlgorithms DESTINATION bin)