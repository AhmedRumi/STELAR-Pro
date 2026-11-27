CXX ?= g++
NVCC ?= /usr/local/cuda/bin/nvcc
CUDA_ARCH ?= all-major
BUILD ?= build
CPPFLAGS += -Iinclude -DSTELAR_SOURCE_ROOT='"$(abspath ..)"'
CXXFLAGS ?= -O3 -g -std=c++17 -ffp-contract=off -Wall -Wextra -pthread
CORE = src/core.cpp src/s4.cpp src/s6.cpp
OBJECTS = $(patsubst src/%.cpp,$(BUILD)/%.o,$(CORE))
.PHONY: all cpu gpu test test-gpu clean
all: cpu
cpu: bin/stelar-pro-cpp-cpu
gpu: bin/stelar-pro-cpp
$(BUILD) bin:
	mkdir -p $@
$(BUILD)/%.o: src/%.cpp include/stelar/core.hpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
$(BUILD)/cuda.o: cuda/backend.cu cuda/glibc_math_compat.hpp include/stelar/core.hpp | $(BUILD)
	$(NVCC) -std=c++17 -O3 -lineinfo -arch=$(CUDA_ARCH) -Xcompiler=-include,$(abspath cuda/glibc_math_compat.hpp) $(CPPFLAGS) -c $< -o $@
bin/stelar-pro-cpp-cpu: $(OBJECTS) $(BUILD)/main.o $(BUILD)/cuda_stub.o | bin
	$(CXX) $(CXXFLAGS) $^ -o $@
bin/stelar-pro-cpp: $(OBJECTS) $(BUILD)/main.o $(BUILD)/cuda.o | bin
	$(NVCC) -std=c++17 -O3 -arch=$(CUDA_ARCH) -Xcompiler=-pthread $^ -o $@
$(BUILD)/test_s4: test/test_s4.cpp $(OBJECTS) $(BUILD)/cuda_stub.o | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@
$(BUILD)/test_core: test/test_core.cpp $(OBJECTS) $(BUILD)/cuda_stub.o | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@
$(BUILD)/test_s6: test/test_s6.cpp src/s6.cpp $(filter-out $(BUILD)/s6.o,$(OBJECTS)) $(BUILD)/cuda_stub.o | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $(filter-out src/s6.cpp,$^) -o $@
$(BUILD)/test_cuda: test/test_cuda.cpp $(OBJECTS) $(BUILD)/cuda.o | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $(BUILD)/test_cuda.o
	$(NVCC) -std=c++17 -Xcompiler=-pthread $(BUILD)/test_cuda.o $(OBJECTS) $(BUILD)/cuda.o -o $@
test: cpu $(BUILD)/test_s4 $(BUILD)/test_core $(BUILD)/test_s6
	$(BUILD)/test_core
	$(BUILD)/test_s4
	$(BUILD)/test_s6
	python3 test/parity.py --native bin/stelar-pro-cpp-cpu --output $(BUILD)/parity
test-gpu: gpu $(BUILD)/test_cuda
	$(BUILD)/test_cuda
clean:
	$(RM) $(BUILD)/*.o $(BUILD)/test_s4 bin/stelar-pro-cpp bin/stelar-pro-cpp-cpu
