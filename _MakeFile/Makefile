BAMTOOLS=../Bamtools/include
BAMTOOLS_LD=../Bamtools/lib
#RootPath = ../ScaffoldGapFillerTest/
#NeedleManPath = ../Depends/needleman_wunsch/
#NormDependsPath = ../Depends/
#INCLUDES:= ../Depends ../Depends/needleman_wunsch ../Bamtools/include/  # 头文件目录

CC = g++
CFLAGS = -O3 -Wall -static -I$(BAMTOOLS) -L$(BAMTOOLS_LD) -Wl,-rpath,$(BAMTOOLS_LD)

OBJ = ./obj/

all: ScaffoldFillerTest

needleman_wunsch.o: needleman_wunsch.c needleman_wunsch.h
	$(CC) $(CFLAGS) -c needleman_wunsch.c

clsalgorithm.o : clsalgorithm.cpp clsalgorithm.h  needleman_wunsch.h
	$(CC) $(CFLAGS) -c clsalgorithm.cpp

clscorestructure.o: clscorestructure.cpp clscorestructure.h
	$(CC) $(CFLAGS) -c clscorestructure.cpp

clsscaffoldfiller.o: clsscaffoldfiller.cpp clsscaffoldfiller.h
	$(CC) $(CFLAGS) -c clsscaffoldfiller.cpp

alignment_scoring.o: alignment_scoring.c alignment_scoring.h
	$(CC) $(CFLAGS) -c alignment_scoring.c

alignment.o: alignment.c alignment.h
	$(CC) $(CFLAGS) -c alignment.c
	
bam_parse.o : bam_parse.cpp bam_parse.h
	$(CC) $(CFLAGS) -c bam_parse.cpp

local_alignment.o: local_alignment.cpp local_alignment.h
	$(CC) $(CFLAGS) -c local_alignment.cpp

fai_parser.o : fai_parser.cpp fai_parser.h
	$(CC) $(CFLAGS) -c fai_parser.cpp

KmerUtils.o :KmerUtils.cpp KmerUtils.h
	$(CC) $(CFLAGS) -c KmerUtils.cpp

smith-waterman.o : smith-waterman.cpp smith-waterman.h
	$(CC) $(CFLAGS) -c smith-waterman.cpp

fasta_parser.o : fasta_parser.cpp fasta_parser.h
	$(CC) $(CFLAGS) -c fasta_parser.cpp

main.o: main.cpp 
	$(CC) $(CFLAGS) -c main.cpp

ScaffoldFillerTest: clsalgorithm.o clscorestructure.o clsscaffoldfiller.o alignment.o alignment_scoring.o needleman_wunsch.o bam_parse.o local_alignment.o \ fai_parser.o KmerUtils.o smith-waterman.o fasta_parser.o main.o 
	$(CC) $(CFLAGS) -o ScaffoldFiller clsalgorithm.o clscorestructure.o clsscaffoldfiller.o alignment.o alignment_scoring.o needleman_wunsch.o bam_parse.o local_alignment.o fai_parser.o KmerUtils.o smith-waterman.o fasta_parser.o main.o \
	-lbamtools -lz


RM-F := rm -f
.PHONY : all deps objs clean veryclean rebuild info
clean :	
	@$(RM-F) *.o
	@$(RM-F) *.d

veryclean: clean
	@$(RM-F) ScaffoldFiller


