DISTRIB=$(PWD)
include $(DISTRIB)/config.mk

CMD_OPTIMIZED_COMPILE=make $@ OPTIMIZE_FLAGS="-O3 -DNDEBUG" DEBUG_FLAGS="-g"


all:
	+make -f Makefile.version.mk DEBUG= FLAGS="-O3 -DNDEBUG"
	+make -f Makefile.version.mk DEBUG=Debug FLAGS="-g"

dep:
	make -f Makefile.version.mk $@

### all ### provide help by default
help:
	@echo ""
	@echo "Makefile for MNS"
	@echo "written by Yvan Mokwinski"
	@echo ""
	@echo ""
	@echo "make all"
	@echo "     build optimized and debug programs"
	@echo ""
	@echo "make build_doc"
	@echo "     build code documentation using Doxygen"
	@echo "     the  main page of the generated documentation is "
	@echo "     doc/html/html/index.html"
	@echo ""
	@echo "make build_report"
	@echo "     build the report for the practicum"
	@echo ""
	@echo "make clean"
	@echo "     clean object files"
	@echo ""
	@echo "make realclean"
	@echo "     clean object files and Doxygen documentation"
	@echo ""
	@echo ""

Doc:
	cd doc;doxygen Doxyfile

clean_doc:
	rm -rf doc/html

build_report:
	cd report;make $@

clean_report:
	cd report;make $@

cleandistrib:
	rm -rf $(PLATFORM)_$(CC)
	rm -rf $(PLATFORM)_$(CC)Debug
	rm -rf bin/*
	rm -rf lib/*
clean:cleandistrib


realclean:cleandistrib clean_doc clean_report


#
# TESTING
#
test:
	+make -f Makefile.version.mk $@

