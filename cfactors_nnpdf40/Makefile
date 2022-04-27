all: rmxeq mkxeq
.PHONY: all

# main programs to be compiled 
EXEGGUTORS = runE866 runE866ratio runE605 runE906ratio
CFACGEN = runCfactors
MAIN = $(EXEGGUTORS) $(CFACGEN)

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

# run'em all
run: mkxeq
	@echo "Running all cfactors, this could take long..."
	@mkdir -p output
	@echo "Running NLO..."
	@$(foreach var,$(EXEGGUTORS),./$(var) 1;)
	@echo "Running NNLO..."
	@$(foreach var,$(EXEGGUTORS),./$(var) 1 1;)
	@echo "Computing cfactors..."
	@./$(CFACGEN)
	@echo "... done"


# remove old executables
rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"		

# clean directory 

clean: rmxeq
	@ -rm -rf *.o *~ $(MAIN)
.PHONY: clean run

############################################################################
