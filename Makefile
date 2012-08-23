.PHONY: all
all: Adaptor_trimmer Guess_fastq_format Quality_trimmer
targets=Adaptor_trimmer Guess_fastq_format Quality_trimmer
CXX := g++

INC = -I SeqAn1.3 -I Boost1.50

DEBUG_FLAGS = -O0 -g3 
RELEASE_FLAGS = -O3 
NOASSERT_FLAGS = -DNDEBUG

CXXFLAGS = -O3 $(INC)

FIND :=$(shell which find)

LINKS = $(shell $(FIND) libs -type f)
GENERAL_LIST = $(shell $(FIND) data) \
			   Makefile LICENSE VERSION AUTHORS README.rst INSTALL.rst \
			   copy_seqan_dependence.py \
			   $(shell $(FIND) . -name "*.sh")
SRC_PKG_LIST = $(wildcard *.h) \
               $(wildcard *.cc) \
               $(shell $(FIND) SeqAn1.3 Boost1.50) \
               $(GENERAL_LIST)
ALL_FILES_SEQAN = $(shell $(FIND) SeqAn1.3 -name "*.h")
ALL_FILES_BOOST = $(shell $(FIND) Boost1.50 -name "*.hpp")
ALL_FILES = $(ALL_FILES_SEQAN) $(ALL_FILES_BOOST)
VERSION = $(shell cat VERSION)
Adaptor_trimmer: Fasta_reader.h seq.h $(ALL_FILES)
	$(CXX) $(CXXFLAGS) $(LINKS) Adaptor_trimmer.cc -o $@
Guess_fastq_format: tools.h $(ALL_FILES)
	$(CXX) $(CXXFLAGS) $(LINKS) Guess_fastq_format.cc -o $@
Quality_trimmer: tools.h $(ALL_FILES) 
	$(CXX) $(CXXFLAGS) $(LINKS) Quality_trimmer.cc -o $@

Adaptor_trimmer.tgz: $(SRC_PKG_LIST)
	rm -rf Adaptor_trimmer.tgz
	rm -rf Adaptor_trimmer-*
	chmod a+x *.sh
	rm -rf tmp
	mkdir tmp
	mkdir tmp/Adaptor_trimmer-$(VERSION)
	tar -czvf tmp.tgz $(SRC_PKG_LIST)
	mv tmp.tgz tmp/Adaptor_trimmer-$(VERSION)
	cd tmp/Adaptor_trimmer-$(VERSION) ; tar -xzvf tmp.tgz ; rm -f tmp.tgz
	cd tmp ; tar -czvf $@ Adaptor_trimmer-$(VERSION)
	cp tmp/$@ .
	rm -rf tmp
.PHONY: clean
clean:
	rm -rf Adaptor_trimmer-*
	rm -rf Adaptor_trimmer.tgz
	rm -rf *.swp
	rm -rf $(targets)
