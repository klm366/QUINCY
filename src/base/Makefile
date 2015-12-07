#---------------------------------------------------------------------
include ../../util/Make.inc.local
#---------------------------------------------------------------------
MODEL_LIB = $(LIBDIR)/lib_base.a
SXMODEL_LIB = $(MODEL_LIB)

#-
RM      = rm -f
STRIP   = strip
SIZE    = size

#-
MODS1 = mo_kind.f90 \
        mo_base_constants.f90 \
        mo_base_memory.f90

OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)

#-
all:
	@$(M_K) m_all
	@echo base is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	@$(RM) $*.o

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(mo_base_constants.o): \
  $(MODEL_LIB)(mo_kind.o) 

$(MODEL_LIB)(mo_base_memory.o): \
  $(MODEL_LIB)(mo_kind.o) \
  $(MODEL_LIB)(mo_base_constants.o) 

