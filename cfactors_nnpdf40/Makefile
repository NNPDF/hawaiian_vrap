all: rmxeq mkxeq
.PHONY: all

# main programs to be compiled 
MAIN = runE866 runE866ratio runE605 runE906ratio runCfactors

############################## do not change #############################

SHELL=/bin/bash

GCC=g++

LDFLAGS += -O3 -std=c++11
LDFLAGS += $(shell lhapdf-config --ldflags)

GLINKER=$(GCC) 
PGMSF= $(MAIN)

-include $(addsuffix .d,$(PGMS))

# rule to compile source programs
$(addsuffix .o,$(PGMSF)): %.o: %.cc Makefile
	$(GCC) $< -c $(LDFLAGS) 

# rule to link object files
$(MAIN): %: %.o  Makefile
	$(GLINKER) $<  $(LDFLAGS) -o $@

# produce executables
mkxeq: $(MAIN)

# remove old executables
rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"		

# clean directory 

clean: rmxeq
	@ -rm -rf *.o *~ $(MAIN)
.PHONY: clean

############################################################################
