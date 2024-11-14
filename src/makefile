
gcc = gcc

flags = -g -O3

libs = -lm -lgsl -lgslcblas -lfftw3

arc = bbh.c  alloc.c  ic.c  outscalar.c  outtensor.c  outvector.c \
      del.c  lie.c  inverse.c  connections.c  ricci.c  metric.c  \
      cn_scheme.c  diverge.c 

all:
	@$(gcc) $(flags) -o execut $(arc) $(libs)
