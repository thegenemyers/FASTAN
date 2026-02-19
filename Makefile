DEST_DIR = ~/bin

CFLAGS = -g -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

CC = gcc

ALL = FasTAN 

all: $(ALL)

GDB.c: gene_core.c gene_core.h
GDB.h: gene_core.h

FasTAN: FasTAN.c alncode.c alncode.h ANO.c ANO.h GDB.c GDB.h ONElib.c ONElib.h align.c align.h
	$(CC) $(CFLAGS) -o FasTAN FasTAN.c alncode.c align.c ANO.c GDB.c gene_core.c ONElib.c -lm -lz

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f FasTAN.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf FasTAN.tar.gz LICENSE README.md Makefile *.h *.c
