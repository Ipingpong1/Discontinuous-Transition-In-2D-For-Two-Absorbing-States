CC ?= gcc
CFLAGS ?= -O3 -march=native -fopenmp -Wall -Wextra -std=c11

dnsim: dnsim.c
	$(CC) $(CFLAGS) -o $@ $< -lm

test: dnsim
	./dnsim --selftest

clean:
	rm -f dnsim
.PHONY: test clean
