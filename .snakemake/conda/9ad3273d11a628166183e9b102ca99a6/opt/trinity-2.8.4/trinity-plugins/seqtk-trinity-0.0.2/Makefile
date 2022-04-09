
all:seqtk-trinity

seqtk-trinity:seqtk.c khash.h kseq.h
		gcc $(CFLAGS) seqtk.c -o $@ -lz -lm

clean:
		rm -fr gmon.out *.o ext/*.o a.out seqtk trimadap *~ *.a *.dSYM session*
