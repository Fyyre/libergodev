CC      ?= cc
AR      ?= ar
PYTHON  ?= python3
LATEX   ?= pdflatex

CFLAGS  := -O3 -Wall -Wextra -Werror -pedantic -std=c99 -fPIC -Iinclude
LDFLAGS := -lm

BUILD_DIR := build
LIB_SO    := $(BUILD_DIR)/libergodev.so
LIB_A     := $(BUILD_DIR)/libergodev.a
TEST_BIN  := $(BUILD_DIR)/test_ergodev
PDF_OUT   := docs/symbolic_ergodeviation.pdf

.PHONY: all check test-c test-py pdf clean

all: $(LIB_SO) $(LIB_A) $(TEST_BIN)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Build object file for library (excludes main() in ergodev.c via -DERGODEV_LIB_ONLY if desired,
# or links directly; here we compile the library object without main using a macro guard)
$(BUILD_DIR)/ergodev.o: src/ergodev.c include/ergodev.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -DERGODEV_NO_MAIN -c src/ergodev.c -o $@

$(LIB_SO): $(BUILD_DIR)/ergodev.o
	$(CC) -shared $^ $(LDFLAGS) -o $@

$(LIB_A): $(BUILD_DIR)/ergodev.o
	$(AR) rcs $@ $^

$(TEST_BIN): src/ergodev.c include/ergodev.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) src/ergodev.c $(LDFLAGS) -o $@

test-c: $(TEST_BIN)
	@echo "=== Running C99 O(1) Streaming Engine Self-Test ==="
	./$(TEST_BIN)

test-py:
	@echo ""
	@echo "=== Running Python Exact Rational & Monte Carlo Verifier ==="
	$(PYTHON) scratch/verify_es.py

check: test-c test-py

pdf:
	cd docs && $(LATEX) -interaction=nonstopmode symbolic_ergodeviation.tex && $(LATEX) -interaction=nonstopmode symbolic_ergodeviation.tex

clean:
	rm -rf $(BUILD_DIR) docs/*.aux docs/*.log docs/*.out $(PDF_OUT)