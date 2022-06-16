







AN_setter=${HOME}/AN_common/AN-common-config
io_setter=${HOME}/root_macros/io_lib/iolib-config
ccflg=`${FASTJET3}/fastjet-config --cxxflags` `root-config --cflags` `${io_setter} -I` `${AN_setter} -I`  -I${ROOUNFOLD}/src
ldflg=`${FASTJET3}/fastjet-config --libs` `root-config --glibs` `${io_setter} -L` `${AN_setter} -L`  -L${ROOUNFOLD} -lRooUnfold

LIB_FASTJET=`${FASTJET3}/fastjet-config --cxxflags --libs`
LIB_ROOT=`root-config --cflags --glibs`
LIB_TRI= ${LIB_ROOT} ${LIB_FASTJET} `${io_setter} -L` `${AN_setter} -L` -L${ROOUNFOLD} -lRooUnfold


# compilation option
CC=g++
CFLAGS=-std=c++11 -O3 -Wno-deprecated
CFLAGS_CHECK=-std=c++11 -O0 -Wno-deprecated -g

bin/main: obj/events.o \
          obj/main.o   \
          obj/MemTimeProgression.o \
          obj/test_loop.o \
          obj/rooResF.o \
		  obj/friend_fnc.o \
          obj/rooResF_check.o \
          obj/big_rooResF.o \
          obj/cut_rooResF.o \
          obj/rat_rooResF.o \
          obj/crazy.o \
          obj/list_ids.o \
          obj/large_ones.o \
          obj/crazy_limit.o \
          obj/sane_bins.o \
          obj/bins70.o \
          obj/jet_draw.o \
          obj/five_rooResF.o \
          obj/thesis_emb.o \
          obj/eta_match.o \
          obj/sys_err.o
	${CC} ${CFLAGS} -o $@ $^ ${LIB_TRI} 

obj/events.o: src/events.cxx src/events.h src/lists.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/main.o: src/main.cxx src/events.h src/MemTimeProgression.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/MemTimeProgression.o: src/MemTimeProgression.cxx src/MemTimeProgression.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/test_loop.o: src/test_loop.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/rooResF.o: src/rooResF.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/friend_fnc.o: src/friend_fnc.cxx src/friend_fnc.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/rooResF_check.o: src/rooResF_check.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/big_rooResF.o: src/big_rooResF.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/cut_rooResF.o: src/cut_rooResF.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/rat_rooResF.o: src/rat_rooResF.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/crazy.o: src/crazy.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/list_ids.o: src/list_ids.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/large_ones.o: src/large_ones.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/crazy_limit.o: src/crazy_limit.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/sane_bins.o: src/sane_bins.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/bins70.o: src/bins70.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/jet_draw.o: src/jet_draw.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/five_rooResF.o: src/five_rooResF.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/thesis_emb.o: src/thesis_emb.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/eta_match.o: src/eta_match.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/sys_err.o: src/sys_err.cxx src/events.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@
