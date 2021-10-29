# Makefile
# ISV/FSV

# set these for the test program
DATAPOINTS=2000
BLOCKSIZE=50
ITERS=10000

CC=gcc

# Only one of the CFLAGS should be uncommented
# uncomment the next one during testing
#CFLAGS=-ansi -Wall -pedantic -g -DISV_CHK
# uncomment the following for fast code for real runs
CFLAGS=-ansi -Wall -pedantic -O2

LDFLAGS=-lgsl -lgslcblas

EXEC=test_isv gendata_isv isvprog gendata_fsv fsvprog foreprog


# running codes

isv.tab: isvprog isvdata.dat
	nice ./isvprog $(ITERS) $(DATAPOINTS) $(BLOCKSIZE) > isv.tab
	echo "Done."

isvdata.dat: gendata_isv
	gendata_isv $(DATAPOINTS) > isvdata.dat

test: test_isv
	./test_isv

fsv.tab: fsvprog fsvdata.dat
	nice ./fsvprog $(ITERS) $(DATAPOINTS) $(BLOCKSIZE) > fsv.tab
	echo "Done."

forecast.tab: foreprog fsvdata.dat
	nice ./foreprog $(ITERS) $(DATAPOINTS) $(BLOCKSIZE) > forecast.tab

fsvdata.dat: gendata_fsv
	gendata_fsv $(DATAPOINTS) > fsvdata.dat


# compiling programs

isvprog: isvprog.o isv.o isv.h
	$(CC) $(CFLAGS) isvprog.o isv.o -o isvprog $(LDFLAGS)

gendata_isv: gendata_isv.o
	$(CC) $(CFLAGS) gendata_isv.o -o gendata_isv $(LDFLAGS)

test_isv: test_isv.o isv.o isv.h
	$(CC) $(CFLAGS) test_isv.o isv.o -o test_isv $(LDFLAGS)

fsvprog: fsvprog.o isv.o fsv.o isv.h fsv.h
	$(CC) $(CFLAGS) fsvprog.o isv.o fsv.o -o fsvprog $(LDFLAGS)

foreprog: foreprog.o forecast_fsv.o isv.o fsv.o isv.h fsv.h forecast_fsv.h
	$(CC) $(CFLAGS) foreprog.o forecast_fsv.o isv.o fsv.o -o foreprog $(LDFLAGS)

gendata_fsv: gendata_fsv.o fsv.h
	$(CC) $(CFLAGS) gendata_fsv.o -o gendata_fsv $(LDFLAGS)


# individual objects (*.o) are built by _implicit_ rules...


# misc utilities
clean:
	rm -f core *.o *~ *.tab *.out *.ind *.dat $(EXEC)

edit:
	gnuclient README.txt Makefile *.h *.c &

print:
	a2ps README.txt Makefile *.h *.c

# eof

