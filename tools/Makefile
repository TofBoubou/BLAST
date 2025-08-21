CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall
MPP_ROOT = ../libs/mutationpp
EIGEN_ROOT = ../libs/eigen
INCLUDES = -I$(MPP_ROOT)/install/include -I$(EIGEN_ROOT)
LDFLAGS = -L$(MPP_ROOT)/install/lib -lmutation++ -Wl,-rpath,$(MPP_ROOT)/install/lib

TARGET = temp_enthalpy_converter

all: $(TARGET)

$(TARGET): temp_enthalpy_converter.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: all clean