# 指定最低可用的CMake版本, 3.15对vscode的CMake插件支持较好
cmake_minimum_required(VERSION 3.15)

# 定义项目名称BDS_SDR_SIM, 并指定使用C++语言
project(BDS_SDR_SIM LANGUAGES CXX)

# 设置C++标准为C++17
set(CMAKE_CXX_STANDARD 17)
# 强制使用所指定的C++标准, 如果编译器不支持C++17则报错
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# 禁止使用编译器特有的扩展，如GNU扩展，强制使用标准C++
set(CMAKE_CXX_EXTENSIONS OFF)

file(GLOB GNSS_CORE_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

# 构建名称为gnss_core的target，类型为静态库，其包含src目录下的所有.cpp文件
add_library(gnss_core
    ${GNSS_CORE_SOURCES}
)

# 为gnss_core库指定include搜索路径
target_include_directories(gnss_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# if (MSVC)
#     add_compile_options(/utf-8)
# endif()


if (MSVC)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()


# 定义名称为BDS_SDR_SIM的可执行文件target，其入口源文件为src/main.cpp
add_executable(${PROJECT_NAME}
    src/main.cpp
)

# 链接核心库（依赖关系声明）
target_link_libraries(${PROJECT_NAME}
    PRIVATE gnss_core
)
