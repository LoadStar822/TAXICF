cmake_minimum_required(VERSION 3.10)
project(TAXICF VERSION 0.1.0)
add_definitions(-DTAXICF_VERSION="${PROJECT_VERSION}")
set(PROJECT_ROOT_DIR "${CMAKE_SOURCE_DIR}")
add_definitions(-DPROJECT_ROOT_DIR="${PROJECT_ROOT_DIR}")


set(CMAKE_BUILD_TYPE Release)

set( CMAKE_CXX_STANDARD 20 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

if( NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" )
    message( FATAL_ERROR
        "Compiler id '${CMAKE_CXX_COMPILER_ID}' is not supported, please \
        check the documentation." )
endif()



# 设置编译选项
include(CheckCXXCompilerFlag)

# 检测是否支持 AVX2 指令集
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)

if(COMPILER_SUPPORTS_AVX2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -mavx2 -flto")
    message(STATUS "AVX2 instructions are supported and enabled.")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -flto")
    message(WARNING "AVX2 instructions are not supported. Proceeding without AVX2 optimization.")
endif()

# 设置SeqAn3库的路径
set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/libs/seqan3/share/cmake/seqan3")
set(SEQAN3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/libs/seqan3/include")

find_package( BZip2 REQUIRED )
find_package( ZLIB REQUIRED )
find_package(Threads REQUIRED)
find_package (seqan3 REQUIRED)

# OpenMP（可选）
find_package(OpenMP)


# 添加可执行文件
add_executable(TAXICF
    TAXICF.cpp
    src/build/TAXICFBuild.cpp
    src/classify/TAXICFClassify.cpp
    libs/xxhash.c
)


# 链接库
target_link_libraries(TAXICF
    PRIVATE
    BZip2::BZip2
    ZLIB::ZLIB
    seqan3::seqan3
    Threads::Threads
)

# 设置包含目录
target_include_directories(TAXICF
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/build
    ${PROJECT_SOURCE_DIR}/src/classify
    ${PROJECT_SOURCE_DIR}/src/utils
    ${PROJECT_SOURCE_DIR}/libs
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/cereal/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/sdsl-lite/include
)


if (OpenMP_CXX_FOUND)
    target_link_libraries(TAXICF PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP found: enabling parallel build.")
else()
    message(WARNING "OpenMP not found: building without parallelism.")
endif()


# ICF 单元测试可执行文件
add_executable(icf_filter_tests
    tests/icf/interleaved_cuckoo_filter_test.cpp
    libs/xxhash.c
)

target_link_libraries(icf_filter_tests
    PRIVATE
    BZip2::BZip2
    ZLIB::ZLIB
    seqan3::seqan3
    Threads::Threads
)

if (OpenMP_CXX_FOUND)
    target_link_libraries(icf_filter_tests PRIVATE OpenMP::OpenMP_CXX)
endif()

target_include_directories(icf_filter_tests
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/build
    ${PROJECT_SOURCE_DIR}/src/classify
    ${PROJECT_SOURCE_DIR}/src/utils
    ${PROJECT_SOURCE_DIR}/libs
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/cereal/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/sdsl-lite/include
)

add_executable(syncmer_tests
    tests/utils/syncmer_test.cpp
)

target_link_libraries(syncmer_tests
    PRIVATE
    BZip2::BZip2
    ZLIB::ZLIB
    seqan3::seqan3
    Threads::Threads
)

if (OpenMP_CXX_FOUND)
    target_link_libraries(syncmer_tests PRIVATE OpenMP::OpenMP_CXX)
endif()

target_include_directories(syncmer_tests
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/utils
    ${PROJECT_SOURCE_DIR}/libs
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/cereal/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/sdsl-lite/include
)

# 忽略 register 关键字警告
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wno-register)
    add_compile_options(-Wno-sign-compare) # 忽略 signed/unsigned 比较警告
endif()


# 安装规则
install(TARGETS TAXICF
    DESTINATION bin
)
