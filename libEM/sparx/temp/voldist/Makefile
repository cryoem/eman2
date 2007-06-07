CC     = mpic++
CFLAGS = -O3

.cpp.o:
	$(CC) $(CFLAGS) -c $<

OBJS = runcartrec.o utilcomm2d.o

runcartrec: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o runcartrec

clean:
	rm -f *.o runcartrec
