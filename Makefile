# Simple Makefile for BLAST project

CXX ?= g++
CXXFLAGS ?= -std=c++23 -O2 -Wall -Wextra -Wfatal-errors -Iinclude -Ilibs/eigen -Ilibs/yaml-cpp/include -Ilibs/mutationpp/install/include
LDFLAGS ?= -Llibs/yaml-cpp/build -Llibs/mutationpp/install/lib -lyaml-cpp -lmutation++

SOURCES := $(shell find src -name '*.cpp')
OBJECTS := $(patsubst src/%.cpp,build/%.o,$(SOURCES))
TARGET ?= blast

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

build/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf build $(TARGET)