#---------------------------------------------------------------------
include ../../util/Make.inc.local
#---------------------------------------------------------------------
BASE_LIB = $(LIBDIR)/lib_base.a
SXBASE_LIB = $(BASE_LIB)
MODEL_LIB = $(LIBDIR)/lib_land.a
SXMODEL_LIB = $(MODEL_LIB)

RM      = rm -f
STRIP   = strip
SIZE    = size

#-
MODS1 = mo_land_constants.f90 \
		mo_land_memory.f90 \
		mo_land_parameters.f90

OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)

#-
all:
	@$(M_K) libbase
	@$(M_K) m_all
	@echo land is OK

land:
	@$(M_K) m_all
	@echo land is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

libbase:
	(cd ../base; $(M_K) -f Makefile)

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	@$(RM) $*.o

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(mo_land_constants.o): \
  $(BASE_LIB)(mo_kind.o)

$(MODEL_LIB)(mo_land_parameters.o): \
  $(BASE_LIB)(mo_kind.o)

$(MODEL_LIB)(mo_land_memory.o): \
  $(BASE_LIB)(mo_kind.o)
