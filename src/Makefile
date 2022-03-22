F90 = ftn
F77 = ftn
OPTIMIZE =  -O4 -w 
FFTWOMP = -lfftw3f_omp
MODFLAG = -module 
OMPLIB = -qopenmp 

CC = cc
C++ = CC
CCOPTIMIZE = -w 
CFITSIO_LIBDIR=/global/homes/m/${USER}/cfitsio/lib/

CCOPTIONS = $(CDEFS) $(CCOPTIMIZE)

# PKS2MAP MODULE AND OBJECT FILES
pmdir  = ./pks2map
p2m_mods = \
         $(pmdir)/halomodels.o\
         $(pmdir)/mpivars.o\
         $(pmdir)/textlib.o\
         $(pmdir)/random.o\
         $(pmdir)/healpixvars.o\
         $(pmdir)/flatskyvars.o\
         $(pmdir)/fitsvars.o\
         $(pmdir)/cosmology.o\
         $(pmdir)/diffuse_profiles.o\
         $(pmdir)/line_profile.o\
         $(pmdir)/integrate_profiles.o\
         $(pmdir)/maptable.o\
         $(pmdir)/haloproject.o\
         $(pmdir)/pksc.o


p2m_objs = $(pmdir)/pks2map.o 	 

mmt_mods = \
         $(pmdir)/textlib.o\
         $(pmdir)/cosmology.o\
         $(pmdir)/halomodels.o\
         $(pmdir)/diffuse_profiles.o\
         $(pmdir)/integrate_profiles.o\
         $(pmdir)/maptable.o

mmt_objs = $(pmdir)/make_maptable.o

hpxdir = ./healpix
hpx_mods = \
            $(hpxdir)/healpix_types.o\
            $(hpxdir)/cgetEnvironment.o\
            $(hpxdir)/extension.o\
            $(hpxdir)/long_intrinsic.o\
            $(hpxdir)/misc_utils.o \
            $(hpxdir)/num_rec.o\
            $(hpxdir)/bit_manipulation.o\
            $(hpxdir)/indmed.o\
            $(hpxdir)/statistics.o\
            $(hpxdir)/pix_tools.o\
            $(hpxdir)/fitstools.o\
            $(hpxdir)/head_fits.o

bindir = ../bin

OPTIONS       =  $(OPTIMIZE) $(MODFLAG)$(hpxdir) $(OMPLIB)

EXEC_t = make_maptable
OBJS_t = $(mmt_mods) $(mmt_objs)

EXEC_p = pks2map
OBJS_p = $(hpx_mods) $(p2m_mods) $(p2m_objs)

.SUFFIXES: .o .f .f90 .F90 .c .C

$(pmdir)/%.o: $(pmdir)/%.f90 
	$(F90) $(OPTIONS) -c $< -o $@
$(pmdir)/%.o: $(pmdir)/%.C 
	$(C++) $(OPTIONS) -c $< -o $@

$(hpxdir)/%.o: $(hpxdir)/%.F90
	$(F90) $(OPTIONS) -c $< -o $@
$(hpxdir)/%.o: $(hpxdir)/%.f90
	$(F90) $(OPTIONS) -c $< -o $@
$(hpxdir)/%.o: $(hpxdir)/%.c
	$(CC)  $(CCOPTIONS) -c $< -o $@


$(EXEC_t): $(OBJS_t) 
	$(F90) $(OPTIONS) $(OBJS_t) -lm -o  $(bindir)/$(EXEC_t)  
$(EXEC_p): $(OBJS_p) 
	$(F90) $(OPTIONS) $(OBJS_p) -L$(CFITSIO_LIBDIR) -lm -lcfitsio -o\
	$(bindir)/$(EXEC_p)  

EXEC = $(bindir)/$(EXEC_t) $(bindir)/$(EXEC_p)
OBJS = $(OBJS_t) $(OBJS_p)

clean:
	rm -f $(EXEC) $(OBJS) $(pmdir)/*.mod $(pmdir)/*.o $(hpxdir)/*.mod $(hpxdir)/*.o 

