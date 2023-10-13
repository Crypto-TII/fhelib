# default parameters
if(NOT TIIFHE_N)
	set(TIIFHE_N "(1 << 15)")
endif()

if(NOT TIIFHE_T)
	set(TIIFHE_T "\"340282366920938463463374607431773454337\"")
endif()

if (NOT TIIFHE_Q)
	set(TIIFHE_Q "\\
	288230376154267649, \\
	288230376155185153, \\
	288230376155250689, \\
	288230376156758017, \\
	288230376157413377, \\
	288230376158396417, \\
	288230376160755713, \\
	288230376161280001, \\
	288230376161673217, \\
	288230376161738753, \\
	288230376162459649, \\
	288230376164294657")
endif()

if (NOT TIIFHE_P)
	set(TIIFHE_P "\\
	288230376166850561, \\
	288230376167047169, \\
	288230376168030209, \\
	288230376168947713")
endif()

if (NOT TIIFHE_OMEGA)
	set(TIIFHE_OMEGA 3)
endif()

# set other variables
function(add result x y)
	math(EXPR value "${x} + ${y}")
	set(${result} ${value} PARENT_SCOPE)
endfunction()

function(ceil result x y)
	math(EXPR value "((${x} - 1) / ${y}) + 1")
	set(${result} ${value} PARENT_SCOPE)
endfunction()

function(len result input)
	string(REPLACE "," ";" list "${input}")
	list(LENGTH list length)
	set(${result} ${length} PARENT_SCOPE)
endfunction()

len(TIIFHE_QLEN ${TIIFHE_Q})
len(TIIFHE_PLEN ${TIIFHE_P})
add(TIIFHE_QPLEN ${TIIFHE_QLEN} ${TIIFHE_PLEN})

ceil(TIIFHE_KMAX ${TIIFHE_QLEN} ${TIIFHE_OMEGA})
if (TIIFHE_PLEN GREATER TIIFHE_KMAX)
	set(TIIFHE_KMAX ${TIIFHE_PLEN})
endif()
