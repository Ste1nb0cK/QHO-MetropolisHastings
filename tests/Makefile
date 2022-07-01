SHELL = /bin/sh
#specify .o names
OBJ_DIR = ./build_tests
#ma1n_code.o is not included here to make compiling the tests easier
OBJS =  action_change.o fill.o  Metropolis_Hasting_Parallel.o Metropolis_Hasting_Serial.o transform_u_to_x.o \
 x_sq_Metropolis_Hasting_Serial.o
TESTSDIR = ./tests
#specify to the compiler where to find the headers
INCLUDES = -I "../include"
#compiler specifications
CXX = mpic++
CXXFLAGS = -Wall -O3
SANITIZERS = # -fsanitize=leak -fsanitize=address -fsanitize=leak
DEBUGFLAG = -g
PROFILEFLAG= -pg
#File list used for each test
FILESTESTSIZE = x_sq_Metropolis_Hasting_Serial.cpp action_change.cpp transform_u_to_x.cpp
TESTFLAG = $(shell pkg-config --cflags catch2)

#specify profile report name

test.x: build_tests/test.o
	${CXX} ${CXXFLAGS} ${PROFILEFLAG} ${SANITIZERS} ${INCLUDES} -o $@ \
	$(shell find ${OBJ_DIR} -name '*.o')



build_tests/test.o: test.cpp
	@echo "#Building test.o"
	${CXX} ${CXXFLAGS} ${PROFILEFLAG} ${SANITIZERS} ${INCLUDES} ${TESTFLAG} \
	-c test.cpp -o $@

clean:
	-rm -f ${OBJ_DIR}/*.o foo* *.out *.x
