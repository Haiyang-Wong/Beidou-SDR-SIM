cmake_minimum_required(VERSION 3.15)

project(BDS_SDR_SIM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Threads REQUIRED)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS TRUE CACHE INTERNAL "")
find_package(UHD QUIET)

if (MSVC)
    add_compile_options(/utf-8)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wno-trigraphs)
endif()

set(GNSS_CORE_SOURCES
    src/bds_sdr.cpp
    src/bds_sig.cpp
    src/channel.cpp
    src/geodesy.cpp
    src/gnss_time.cpp
    src/inav_msg.cpp
    src/rinex.cpp
)

add_library(gnss_core
    ${GNSS_CORE_SOURCES}
)

target_include_directories(gnss_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(gnss_core
    PUBLIC
        Threads::Threads
)

if (UHD_FOUND)
    target_compile_definitions(gnss_core PUBLIC BDS_SDR_WITH_UHD)
    target_include_directories(gnss_core PUBLIC ${UHD_INCLUDE_DIRS})
    target_link_libraries(gnss_core PUBLIC ${UHD_LIBRARIES})
    message(VERBOSE "UHD found: USRP realtime streaming enabled")
else()
    message(VERBOSE "UHD not found: build keeps file output; --usrp will report an error at runtime")
endif()

add_executable(${PROJECT_NAME}
    src/main.cpp
)

target_link_libraries(${PROJECT_NAME}
    PRIVATE gnss_core
)
