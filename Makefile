

#FFLAGS = -O2 -e
FFLAGS = -g -e
LOADSPEC2ZRT = spec2zrt.o ffts.o 
#
#
#SACLIB = /usr/local/SAC/lib
SACLIB = /usr/local/sac2000/lib
#
#SACLIBNAME = sac
SACLIBNAME = sac2000
#
#
#spec2zrt:	$(LOADSPEC2ZRT)
#	f77 $(FFLAGS) $(LOADSPEC2ZRT) -o spec2zrt /disk19/sac/lib/sac2000.a
spec2zrt:	$(LOADSPEC2ZRT)
	f77 $(FFLAGS) $(LOADSPEC2ZRT) -o spec2zrt -L$(SACLIB) -lsac2000
spec2zrt.o: spec2zrt.f
	f77 -c -e -u  $(FFLAGS) spec2zrt.f

mmarktime:	mmarktime.o
	f77 $(FFLAGS) -o mmarktime mmarktime.o -L$(SACLIB) -l$(SACLIBNAME)

