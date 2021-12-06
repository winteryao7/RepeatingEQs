SUBS = sacio.o fft.o Complex.o
PROGRAM = wfcc sac_wfcc

all: $(PROGRAM)

$(PROGRAM): %:%.o $(SUBS)
	gcc -o ./$@ $@.o $(SUBS) -lm

bak:
	tar -zcvf wfcc_sac.tar.gz *.c *.h Makefile

clean:
	rm *.o
