# NVCC = nvcc
# NVCC_FLAGS = -O3 -Xcompiler -fPIC -shared

NVCC = nvcc
NVCC_FLAGS = -O3 -Xcompiler -fPIC -shared \
  -gencode arch=compute_70,code=sm_70 \
  -gencode arch=compute_75,code=sm_75 \
  -gencode arch=compute_80,code=sm_80 \
  -gencode arch=compute_86,code=sm_86 \
  -gencode arch=compute_89,code=sm_89 \
  -gencode arch=compute_90,code=sm_90 \
  -gencode arch=compute_90,code=compute_90


TARGET = libweight_calc.so
SRC = weight_calc.cu

$(TARGET): $(SRC)
	$(NVCC) $(NVCC_FLAGS) -o $@ $<

clean:
	rm -f $(TARGET) 