


CXXFLAGS= -std=c++17 -pedantic -Wall -fopenmp
CXXOPTS= -O0 -g

INC_PATH= -I${HOME}/include
LIB_PATH= -L${HOME}/lib
LIB= -lm -lblas -llapack -lfftw3

SRC=src
BIN=bin
BUILD=build


TESTS=tests
TEST_BIN=test_bin


DEPS := $(wildcard $(BUILD)/*.d)
ifneq ($(DEPS),)
include $(DEPS)
endif

# Builds objects from source
$(BUILD)/%.o: $(SRC)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXOPTS) $(INC_PATH) -MMD -MP $< -c -o $@


$(TEST_BIN)/%: $(TESTS)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXOPTS) $(INC_PATH) $(LDFLAGS) $(LIB_PATH) $< -o $@ $(LIB)

# Generic linker for executables
%: $(BUILD)/%.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIB_PATH) -o $(BIN)/$@ $< $(LIB)

# Disallow secondary compilation
%: %.cpp

.SECONDARY:

.PHONY: clean


clean:
	rm $(BUILD)/**

