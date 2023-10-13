configure_file(config.h.in config.h)

add_library(tiifhe SHARED tiifhe.h bgv.h bgv.c
	${CMAKE_CURRENT_BINARY_DIR}/config.h
	random.h random.c util.h util.c)
target_link_libraries(tiifhe crypto ssl)
target_include_directories(tiifhe PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(tiifhe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(TIIFHE_MSG STREQUAL "mpz")
	target_sources(tiifhe PRIVATE msg_mpz.h msg_mpz.c)
	target_link_libraries(tiifhe gmp)
else()
	message(FATAL_ERROR "missing message implementation (possible values: mpz)")
endif()

if(TIIFHE_RING STREQUAL "hexl")
	target_sources(tiifhe PRIVATE ring_hexl.h ring_hexl.cpp)
	target_link_libraries(tiifhe hexl)
elseif(TIIFHE_RING STREQUAL "i64")
	target_sources(tiifhe PRIVATE ring_i64.h ring_i64.c)
else()
	message(FATAL_ERROR "missing ring implementation (possible values: hexl, i64)")
endif()
