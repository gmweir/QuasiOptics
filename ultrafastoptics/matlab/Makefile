SHELL=/bin/sh
MEX=mex
OBJS=stackgdgradmex.mexglx stackgdfastmex.mexglx stackgdgradfastmex.mexglx

%.mexglx: %.c
	$(MEX) $<

all: $(OBJS)

clean:
	rm $(OBJS)
