#---------------------------------------------------------------------
include ../../util/Make.inc.local
#---------------------------------------------------------------------
SHARED_LIB = $(LIBDIR)/lib_shared.a
SXSHARED_LIB = $(SHARED_LIB)
BASE_LIB = $(LIBDIR)/lib_base.a
SXBASE_LIB = $(BASE_LIB)
ATM_LIB = $(LIBDIR)/lib_atmland.a
SXATM_LIB = $(ATM_LIB)
MODEL_LIB = $(LIBDIR)/lib_land.a
SXMODEL_LIB = $(MODEL_LIB)

RM      = rm -f
STRIP   = strip
SIZE    = size

#-
MODS1 = mo_land_constants.f90 \
		mo_land_memory.f90 \
		mo_land_parameters.f90 \
        mo_land_util.f90

OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)

#-
all:
	@$(M_K) libatm
	@$(M_K) libbase
	@$(M_K) libshared
	@$(M_K) m_all
	@echo land is OK

land:
	@$(M_K) m_all
	@echo land is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

libatm:
	(cd ../atmland; $(M_K) -f Makefile atmland)

libbase:
	(cd ../base; $(M_K) -f Makefile)

libshared:
	(cd ../shared; $(M_K) -f Makefile shared)

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	@$(RM) $*.o

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(mo_land_constants.o): \
  $(SHARED_LIB)(mo_kind.o)

$(MODEL_LIB)(mo_land_parameters.o): \
  $(SHARED_LIB)(mo_kind.o)

$(MODEL_LIB)(mo_land_memory.o): \
  $(SHARED_LIB)(mo_kind.o)

$(MODEL_LIB)(mo_land_util.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_memory.o) \
  $(BASE_LIB)(mo_base_util_functions.o) \
  $(BASE_LIB)(mo_base_constants.o) \
  $(ATM_LIB)(mo_atmland_memory.o) \
  $(MODEL_LIB)(mo_land_parameters.o)

