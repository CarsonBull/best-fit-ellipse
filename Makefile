TARGETS = efit eparam
CFLAGS += -std=c99 -Wall -D_GNU_SOURCE
LDLIBS += -lm

.PHONY: all clean

all: $(TARGETS)

clean:
	$(RM) -r $(TARGETS) *.o *.dSYM *~

efit: efit.o dsm.o
eparam: eparam.o dsm.o
