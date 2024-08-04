#OBJECTS = chkcmdline.o readMAT.o readpdb.o sch.o Gcounter.o TOZ.o KAITEN.o VEC_kei.o kaiten.o cd_func.o SF_counter.o sch_tri.o Mcounter.o MVandSCORE.o lsfit.o det_sf.o enephobic_ca.o high_reso_doc.o

CC = gcc
#COM =  readpdb.o  cd_func.o chkcmdline.o rmsd.o
COM =  chkcmdline.o mrc.o readpdb.o mrcfft.o
OBJECTS = $(COM)
HEADER = struct.h func.h mrc.h
PG = VESPER_Power

OPT = -O3 -msse4.2 -msse -msse2 -msse3 -mmmx -ftree-vectorize -ftree-vectorizer-verbose=3 -march=native -funroll-loops -fopenmp

$(PG) : main.c $(OBJECTS) $(HEADER) $(COM)
	$(CC) $(CFLAGS) $(OPT) -o $(PG) main.c $(OBJECTS) -lm -lfftw3_threads -lfftw3f


chkcmdline.o : chkcmdline.c $(HEADER)
	$(CC) -O3 -c chkcmdline.c

readpdb.o : readpdb.c $(HEADER)
	$(CC) $(OPT) -c readpdb.c -lm

mrc.o : mrc.c $(HEADER)
	$(CC) $(OPT) -c mrc.c -lm

mrcfft.o : mrcfft.c $(HEADER) mrcfft.h
	$(CC) $(CFLAGS) $(OPT) -c mrcfft.c -lm -lfftw3_threads -lfftw3f

clean:
	find . -name "*.o" -type f -delete
