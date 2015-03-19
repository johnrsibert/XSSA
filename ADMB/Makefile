# identify some extra file name suffixes
.SUFFIXES:

# tell make not to delete these intermediate targets
.PRECIOUS: %.c %.cpp %.o

# make some special PHONY targets
.PHONY: all help rules clean dox

ifndef ADMB_HOME
  ADMB_HOME=/usr/local/admb
endif

ADMB_XML=/home/jsibert/admb/trunk/tests/xml
VPATH=$(TPL_DIR):$(SRC_DIR):$(ADMB_XML)


# make the "safe" version by default
# to make "the optimized" version, type  `make OPT=TRUE
ifeq ($(OPT),TRUE)
  FLAGS = -f -g
else
  FLAGS = -g
  OPT = FALSE
endif

 
BIN_DIR = ../run

xssams: xssams.tpl
	@echo ------- Translate $<
	admb  -r $(FLAGS) $<
	@echo -------- move $@
	mv -f $@ $(BIN_DIR)

noKxssams: noKxssams.tpl
	@echo ------- Translate $<
	admb  -r $(FLAGS) $<
	@echo -------- move $@
	mv -f $@ $(BIN_DIR)

#####################################################################

rules:
	@echo OPT = $(OPT)
	@echo CURDIR = $(CURDIR)
	@echo PWD = $(PWD)
	@echo OSTYPE = $(OSTYPE)
	@echo OS     = $(OS)
	@echo TERM = $(TERM)
	@echo BIN_DIR = $(BIN_DIR)
	@echo ADMB_HOME = $(ADMB_HOME)
	@echo ADMB_VERSION = $(ADMB_VERSION)
	@echo FLAGS = $(FLAGS)

clean:
	@rm -vf *.o
	@rm -vf *.obj
	@rm -vf *.cpp
	@rm -vf *.htp
	@rm -fv xssams
	@rm -fv noKxssams
	@rm -fv make.out
	



