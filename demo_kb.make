PROJECT=int2         # for historical reasons we always like to 
                     # call the executable int2, but you could set it to 
                     # something more descriptive

###HOST=linux
###HOST=macosx
HOST=linux-openmp
###HOST=macosx-openmp

include make.inc



.PHONY: $(PROJECT) clean list

.f.$(OBJSUF):
	$(FC) $(FFLAGS) $<

.f.$(MODSUF):
	$(FC) $(FFLAGS) $<

.SUFFIXES: $(MODSUF) .$(OBJSUF) .f .c

# SOURCE FILE LIST
#
vpath %.f .:..

FMODS = 

FSRCS = nufft2d_demof90_kb.f nufft2df90_gs.f nufft2df90_kb.f i0eva.f dirft2d_omp.f dfft.f next235.f
##FSRCS = nufft2d_demof90_kb.f nufft2df90.f nufft2df90_kb.f i0eva.f dirft2d_omp.f dfft.f next235.f

ifeq ($(WITH_SECOND),1) 
FSRCS += second-r8.f
endif

#
# object files list
MODS    =  $(FMODS:.f=.$(MODSUF)) 
OBJS    =  $(FMODS:.f=.$(OBJSUF)) $(FSRCS:.f=.$(OBJSUF)) 
#

$(PROJECT):   $(MODS)   $(OBJS)
	rm -f $(PROJECT)
	$(FLINK) $(OBJS)
	./$(PROJECT)
#
clean: 
	rm -f $(OBJS)
# 
list: $(FSRCS)
	echo $^
#
pack: $(FSRCS)
	cat $^ > _tmp_.pack.f
#

