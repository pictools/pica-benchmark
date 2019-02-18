CXX=icpc
CXXFLAGS_COMMON = -O3 -DNDEBUG -g -qopenmp -qopt-zmm-usage=high
CXXFLAGS= $(CXXFLAGS_COMMON) -xCOMMON-AVX512
CXXFLAGS_RCP=$(CXXFLAGS_COMMON) -xCORE-AVX512 -no-prec-div
INCLUDE_FLAGS=-I ../pica/include -I ../src -I ../thirdParty
OUT_DIR=../bin

all: 01-pusher-baseline 02-pusher-vectorized 03-pusher-vectorized-soa 01-full-baseline 02-full-supercells 03-full-supercells-aos 04-full-supercells-vec

01-pusher-baseline: benchmark-utility
	$(CXX) -o $(OUT_DIR)/01-pusher-baseline $(CXXFLAGS) $(INCLUDE_FLAGS) ../src/benchmark/pusher/01-baseline/Main.cpp benchmark-utility.a

02-pusher-vectorized: benchmark-utility
	$(CXX) -o $(OUT_DIR)/02-pusher-vectorized $(CXXFLAGS) $(INCLUDE_FLAGS) ../src/benchmark/pusher/02-vectorized/Main.cpp benchmark-utility.a

03-pusher-vectorized-soa: benchmark-utility
	$(CXX) -o $(OUT_DIR)/03-pusher-vectorized $(CXXFLAGS) $(INCLUDE_FLAGS) ../src/benchmark/pusher/03-vectorized-soa/Main.cpp benchmark-utility.a

04-pusher-vectorized-soa-rcp: benchmark-utility
	$(CXX) -o $(OUT_DIR)/04-pusher-vectorized-soa-rcp $(CXXFLAGS_RCP) $(INCLUDE_FLAGS) ../src/benchmark/pusher/04-vectorized-soa-rcp/Main.cpp benchmark-utility.a

01-full-baseline:
	$(CXX) -o $(OUT_DIR)/01-full-baseline $(CXXFLAGS) $(INCLUDE_FLAGS) ../src/benchmark/full/01-baseline/Main.cpp benchmark-utility.a

02-full-supercells:
	$(CXX) -o $(OUT_DIR)/02-full-supercells $(CXXFLAGS) $(INCLUDE_FLAGS) ../src/benchmark/full/02-supercells/Main.cpp benchmark-utility.a

03-full-supercells-aos:
	$(CXX) -o $(OUT_DIR)/03-full-supercells-aos $(CXXFLAGS) $(INCLUDE_FLAGS) ../src/benchmark/full/03-supercells-AoS/Main.cpp benchmark-utility.a

04-full-supercells-vec:
	$(CXX) -o $(OUT_DIR)/04-full-supercells-vec $(CXXFLAGS) $(INCLUDE_FLAGS) ../src/benchmark/full/04-supercells-vec/Main.cpp benchmark-utility.a

benchmark-utility: prepare
	$(CXX) -c $(CXXFLAGS) $(INCLUDE_FLAGS) ../src/utility/Output.cpp ../src/utility/Parser.cpp ../src/utility/Timer.cpp ../src/utility/Random.cpp
	ar rc benchmark-utility.a Output.o Parser.o Timer.o Random.o

prepare:
	mkdir -p $(OUT_DIR)

clean:
	rm -r $(OUT_DIR)
	rm *.o
