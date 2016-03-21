#---------------------------------------------------------------------
include ../../util/Make.inc.local
#---------------------------------------------------------------------
BASE_LIB = $(LIBDIR)/lib_base.a
SXBASE_LIB = $(BASE_LIB)
MODEL_LIB = $(LIBDIR)/lib_shared.a
SXMODEL_LIB = $(MODEL_LIB)

#-
RM      = rm -f
STRIP   = strip
SIZE    = size

#-
MODS1 = mo_namelist.f90 \
        mo_exception.f90 \
        mo_filename.f90 \
        mo_io_units.f90 \
        mo_util_string.f90

OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)

#-
all:
	@$(M_K) libbase
	@$(M_K) m_all
	@echo shared is OK

shared:
	@$(M_K) m_all
	@echo shared is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

libbase:
	(cd ../base; $(M_K) -f Makefile)

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	@$(RM) $*.o

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(mo_namelist.o): \
  $(MODEL_LIB)(mo_util_string.o) \
  $(MODEL_LIB)(mo_filename.o) \
  $(MODEL_LIB)(mo_exception.o)

$(MODEL_LIB)(mo_exception.o): \
  $(MODEL_LIB)(mo_io_units.o)

$(MODEL_LIB)(mo_filename.o): \
  $(MODEL_LIB)(mo_exception.o)

$(MODEL_LIB)(mo_util_string.o): \
  $(BASE_LIB)(mo_kind.o)


