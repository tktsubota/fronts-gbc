# Variables
BUILD_DIR = build
SRC_DIR = .
CXX_SOURCES = $(shell find $(SRC_DIR) -name '*.cpp' -or -name '*.cxx')
BUILD_TYPE ?= Release
# Armadillo install prefix (see setup.sh)
PREFIX ?= $(HOME)/opt/armadillo

# Build target: generate build files and compile
build: $(BUILD_DIR)/Makefile
	$(MAKE) -C $(BUILD_DIR)

# Generate build system if not present or changes are detected
$(BUILD_DIR)/Makefile: $(CXX_SOURCES) CMakeLists.txt
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake ../$(SRC_DIR) -DCMAKE_PREFIX_PATH=$(PREFIX) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE)

# Clean build directory
clean:
	rm -rf $(BUILD_DIR)

rebuild: clean build

.PHONY: all build clean rebuild