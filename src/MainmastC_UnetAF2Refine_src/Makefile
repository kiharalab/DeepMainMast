#OBJECTS = chkcmdline.o readMAT.o readpdb.o sch.o Gcounter.o TOZ.o KAITEN.o VEC_kei.o kaiten.o cd_func.o SF_counter.o sch_tri.o Mcounter.o MVandSCORE.o lsfit.o det_sf.o enephobic_ca.o high_reso_doc.o

CC = gcc
#COM =  readpdb.o  cd_func.o chkcmdline.o rmsd.o
COM =  chkcmdline.o mrc.o readpdb.o thread.o dp.o rmsd.o model.o
OBJECTS = $(COM)
HEADER = struct.h func.h mrc.h thread.h
PG = MainmastC_UnetAF2Refine

OPT = -O3 -msse4.2 -msse -msse2 -msse3 -mmmx -ftree-vectorize -ftree-vectorizer-verbose=3 -ffast-math -march=native -funroll-loops -fopenmp

$(PG) : main.c $(OBJECTS) $(HEADER) $(COM)
	$(CC) -o $(PG) main.c $(OBJECTS) -lm $(OPT) 


chkcmdline.o : chkcmdline.c $(HEADER)
	$(CC) -O3 -c chkcmdline.c

readpdb.o : readpdb.c $(HEADER)
	$(CC) -O3 -c readpdb.c -lm $(OPT)

mrc.o : mrc.c $(HEADER)
	$(CC) -O3 -c mrc.c -lm $(OPT)

thread.o : thread.c $(HEADER)
	$(CC) -O3 -c thread.c -lm $(OPT)

dp.o : dp.c $(HEADER)
	$(CC) -O3 -c dp.c -lm $(OPT)

rmsd.o : rmsd.c rmsd.h
	$(CC) -O3 -c rmsd.c -lm $(OPT)

model.o : model.c model.h thread.h
	$(CC) -O3 -c model.c -lm $(OPT)

clean:
	find . -name "*.o" -type f -delete
