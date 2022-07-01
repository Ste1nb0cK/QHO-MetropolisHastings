SHELL = /bin/sh
#specify .o names
SRC_DIR = ./src
OBJ_DIR = ./build
TEST_OBJ_DIR = ./tests/build_tests
SRCS = $(shell find ${SRC_DIR} -name '*.cpp')
OBJS = action_change.o fill.o main.o Metropolis_Hasting_Parallel.o Metropolis_Hasting_Serial.o transform_u_to_x.o \
x_sq_Metropolis_Hasting_Parallel.o x_sq_Metropolis_Hasting_Serial.o

TESTSDIR = ./tests
#specify to the compiler where to find the headers
INCLUDES = -I "./include"
#compiler specifications
CXX = mpic++
CXXFLAGS = -Wall
OPTIMIZATIONFLAG = -O3
SANITIZERS = -fsanitize=leak -fsanitize=address -fsanitize=leak
DEBUGFLAG = -g
PROFILEFLAG= -pg
#File list used for each test
FILESTESTSIZE = clouster_matrix.cpp percolation.cpp index_matrix.cpp
TESTFLAG = $(shell pkg-config --cflags catch2)
#specify profile report name
PROF_REPORT = prof_report.txt

#specify where to find the .cpp files
vpath %.cpp ./src
vpath %.o ./build
foo: ${OBJS}
	${CXX} ${CXXFLAGS} ${PROFILEFLAG} ${SANITIZERS} ${INCLUDES} -o $@ \
	$(shell find ${OBJ_DIR} -name '*.o')

clean:
	-rm -f ${OBJ_DIR}/*.o foo* *.out ${PROF_REPORT}
	-cd tests; make clean

#implicit rule for making .o from .cpp
.cpp.o:
	@ -mkdir build
	${CXX} ${CXXFLAGS} ${PROFILEFLAG} ${SANITIZERS} ${INCLUDES} -c $< -o  ${OBJ_DIR}/$@

#format .cpp files using clang format
format:
	clang-format -i --dry-run src/*.cpp
	clang-format -i --dry-run include/*.hpp
#debug the code
debug:
	${CXX} ${DEBUGFLAG}  ${INCLUDES} ${SRCS} -o foo_debug
	chmod +x init_debug.sh
	./init_debug.sh
report:
	pdflatex main.tex
simu1:
	./foo 4 0.6 1

profile:
	./foo 8 0.5 2
	gprof foo gmon.out>${PROF_REPORT}
	cat ${PROF_REPORT}

test:
	cp -v -u ./build/*.o tests/build_tests/
	rm tests/build_tests/main1_code.o
	cd tests; make; ./test.x
