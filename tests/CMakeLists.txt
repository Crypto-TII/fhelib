set(SRC ${PROJECT_SOURCE_DIR}/src)
set(CFG ${CMAKE_CURRENT_BINARY_DIR}/../src)
add_compile_options(-fsanitize=address)
add_link_options(-fsanitize=address)

add_executable(random random.c
	${SRC}/util.h ${SRC}/util.c)
target_link_libraries(random crypto ssl)
target_include_directories(random PRIVATE ${SRC})
target_include_directories(random PRIVATE ${CFG})

add_executable(ring_i64 ring_i64.c
	${CFG}/config.h
	${SRC}/random.h ${SRC}/random.c
	${SRC}/util.h ${SRC}/util.c)
target_link_libraries(ring_i64 crypto gmp ssl)
target_include_directories(ring_i64 PRIVATE ${SRC})
target_include_directories(ring_i64 PRIVATE ${CFG})

add_executable(msg_mpz msg_mpz.c
	${CFG}/config.h
	${SRC}/util.h ${SRC}/util.c)
target_link_libraries(msg_mpz gmp)
target_include_directories(msg_mpz PRIVATE ${SRC})
target_include_directories(msg_mpz PRIVATE ${CFG})

add_executable(bgv bgv.c)
target_link_libraries(bgv tiifhe)
