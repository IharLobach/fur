# ---- User-tweakable ---------------------------------------------------------
PYTHON ?= python3
CXX     ?= clang++
# ---------------------------------------------------------------------------

# Derived values from the active Python (inside your conda env)
PY_VER_NUM   := $(shell $(PYTHON) -c 'import sys; print(f"{sys.version_info[0]}{sys.version_info[1]}")')
PY_INCLUDE   := $(shell $(PYTHON) -c 'import sysconfig; print(sysconfig.get_paths()["include"])')
NUMPY_INC    := $(shell $(PYTHON) -c 'import numpy; print(numpy.get_include())')
EXT_SUFFIX   := $(shell $(PYTHON) -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')
CONDA_PREFIX := $(shell $(PYTHON) -c 'import os; print(os.environ.get("CONDA_PREFIX",""))')

# Boost.Python library names from conda-forge (match your Python version)
BOOST_PY_LIB    := boost_python$(PY_VER_NUM)
BOOST_NUMPY_LIB := boost_numpy$(PY_VER_NUM)

# Common flags (Apple Silicon / macOS)
CXXFLAGS := -O3 -std=c++17 -fPIC \
            -I$(PY_INCLUDE) -I$(NUMPY_INC) -I$(CONDA_PREFIX)/include
 # Let Python decide how to link on macOS (often uses -undefined dynamic_lookup)
PY_LDFLAGS := $(shell $(PYTHON) -c 'import sysconfig; print(sysconfig.get_config_var("LDFLAGS") or "")')

LDFLAGS  := -shared $(PY_LDFLAGS) -L$(CONDA_PREFIX)/lib
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  # Ensure Python symbols are resolved at runtime (no hard link to libpython)
  LDFLAGS += -Wl,-undefined,dynamic_lookup
endif
LDLIBS   := -l$(BOOST_PY_LIB) -l$(BOOST_NUMPY_LIB)

SRC      := coherent_modes_cpp.cpp
TARGET   := coherent_modes_cpp$(EXT_SUFFIX)
OBJ      := $(SRC:.cpp=.o)

all: $(TARGET)

$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

clean:
	rm -f $(OBJ) $(TARGET)

print-vars:
	@echo PY_INCLUDE=$(PY_INCLUDE)
	@echo NUMPY_INC=$(NUMPY_INC)
	@echo EXT_SUFFIX=$(EXT_SUFFIX)
	@echo CONDA_PREFIX=$(CONDA_PREFIX)
	@echo BOOST_PY_LIB=$(BOOST_PY_LIB)
	@echo BOOST_NUMPY_LIB=$(BOOST_NUMPY_LIB)
