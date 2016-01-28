#---------------------------------------------------------------------
include ../../util/Make.inc.local
#---------------------------------------------------------------------
BASE_LIB = $(LIBDIR)/lib_base.a
SXBASE_LIB = $(BASE_LIB)
MODEL_LIB = $(LIBDIR)/lib_atmland.a
SXMODEL_LIB = $(MODEL_LIB)

#-
RM      = rm -f
STRIP   = strip
SIZE    = size

#-
MODS1 = mo_atmland_memory.f90 \
        mo_atmland_constants.f90 \
        mo_atmland_forcing.f90 \
        mo_atmland_util.f90

OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)

#-
all:
	@$(M_K) libbase
	@$(M_K) m_all
	@echo atmland is OK

atmland:
	@$(M_K) m_all
	@echo atmland is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

libbase:
	(cd ../base; $(M_K) -f Makefile)

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	@$(RM) $*.o

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(mo_atmland_memory.o): \
  $(BASE_LIB)(mo_kind.o)

$(MODEL_LIB)(mo_atmland_constants.o): \
  $(BASE_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_constants.o)

$(MODEL_LIB)(mo_atmland_forcing.o): \
  $(BASE_LIB)(mo_kind.o) \
  $(MODEL_LIB)(mo_atmland_constants.o)

$(MODEL_LIB)(mo_atmland_util.o): \
  $(BASE_LIB)(mo_kind.o) \
  $(MODEL_LIB)(mo_atmland_constants.o)
