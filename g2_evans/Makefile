# Evans validado (FLINT >= 3, que inclui arb/acb).  Debian/Ubuntu: apt install libflint-dev
# Se o FLINT estiver em outro prefixo: make FLINT_PREFIX=/caminho
FLINT_PREFIX ?= /usr
CC ?= gcc
CFLAGS ?= -O2 -march=native -Wall -Wno-misleading-indentation
CPPFLAGS += -I$(FLINT_PREFIX)/include
LDFLAGS += -L$(FLINT_PREFIX)/lib -Wl,-rpath,$(FLINT_PREFIX)/lib
LDLIBS = -lflint -lmpfr -lgmp -lm

evans: src/evans.c src/tm.inc src/lohner.inc
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ src/evans.c $(LDFLAGS) $(LDLIBS)

test: evans
	./evans point 1 1 0.5 1 100 -1.5 2 2
	./evans box 1 1 0.5 0.5 0.5 5 5 6 0.25 16 3
	python3 tests/check_reference.py
	python3 tests/test_lamstar_beta.py
	python3 tests/test_idempotente.py
	python3 tests/test_ancoragem.py

clean:
	rm -f evans
