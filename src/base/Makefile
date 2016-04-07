#---------------------------------------------------------------------
include ../../util/Make.inc.local
#---------------------------------------------------------------------
SHARED_LIB = $(LIBDIR)/lib_shared.a
SXSHARED_LIB = $(SHARED_LIB)
MODEL_LIB = $(LIBDIR)/lib_base.a
SXMODEL_LIB = $(MODEL_LIB)

#-
RM      = rm -f
STRIP   = strip
SIZE    = size

#-
MODS1 = mo_base_constants.f90 \
        mo_base_memory.f90 \
        mo_base_util_functions.f90 \
        mo_jsb_config_class.f90

OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)

#-
all:
	@$(M_K) libshared
	@$(M_K) m_all
	@echo base is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

libshared:
	(cd ../shared; $(M_K) -f Makefile shared)

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	@$(RM) $*.o

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(mo_base_constants.o): \
  $(SHARED_LIB)(mo_kind.o)

$(MODEL_LIB)(mo_base_memory.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(MODEL_LIB)(mo_base_constants.o)

$(MODEL_LIB)(mo_base_util_functions.o): \
  $(SHARED_LIB)(mo_kind.o)

$(MODEL_LIB)(mo_jsb_config_class.o): \
  $(SHARED_LIB)(mo_io_units.o)

