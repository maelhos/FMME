
.PHONY: default
default: ku

cpp_poly_source_files := $(shell find src -name *.cpp)
cpp_poly_object_files := $(patsubst src/%.cpp, build/%.o, $(cpp_poly_source_files))
INCLUDE_PATH := -I include

CPP_INCLUDE_PATHS := -I pml/ntl-extras/include -I ntl/include
ntl_extras_object_files := $(shell find pml/ntl-extras/lib -name *.o)

OUT := ku

OMP_LIB := -fopenmp
LIBS := -lgmp -lm -lssl -lcrypto -pthread
OPTS := -march=native -ffast-math -std=c++20 -L ntl/src/ntl.a -O3 -g
CC := g++

$(cpp_poly_object_files): build/%.o : src/%.cpp
	@mkdir -p $(dir $@) && \
	$(CC) $(OMP_LIB) $(OPTS) -c -Wall $(INCLUDE_PATH) $(CPP_INCLUDE_PATHS) $(patsubst build/%.o, src/%.cpp, $@) -o $@ $(LIBS)

.PHONY: ku
ku: $(cpp_poly_object_files) 
	@$(CC) $(OMP_LIB) $(OPTS) -Wall $(CPP_INCLUDE_PATHS) $(INCLUDE_PATH) -o $(OUT) $(cpp_poly_object_files) $(ntl_extras_object_files) ntl/src/ntl.a $(LIBS)

clean:
	@rm -rfv build
	@mkdir build

	@echo "CLEAN"