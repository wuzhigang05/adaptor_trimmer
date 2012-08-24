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
# default link with mac libraries
LINKS = $(shell $(FIND) libs/mac -type f)

ifneq (,$(findstring Linux,$(shell uname)))
LINKS = $(shell $(FIND) libs/linux -type f)
endif
SRC_DIR = src

GENERAL_LIST = $(shell $(FIND) data) \
			   $(shell $(FIND) libs) \
			   $(shell $(FIND) scripts) \
			   $(wildcard *.md) \
			   $(wildcard *.rst) \
			   Makefile LICENSE VERSION AUTHORS 
SRC_PKG_LIST = $(shell $(FIND) src) \
               $(shell $(FIND) SeqAn1.3 Boost1.50) \
               $(GENERAL_LIST)
ALL_FILES_SEQAN = $(shell $(FIND) SeqAn1.3 -name "*.h")
ALL_FILES_BOOST = $(shell $(FIND) Boost1.50 -name "*.hpp")
ALL_FILES = $(ALL_FILES_SEQAN) $(ALL_FILES_BOOST)
VERSION = $(shell cat VERSION)
Adaptor_trimmer: $(SRC_DIR)/Adaptor_trimmer.cc $(SRC_DIR)/Fasta_reader.h $(SRC_DIR)/seq.h $(ALL_FILES) 
	@echo compiling: $@
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LINKS) 
	@echo compiling: $@ done
	mkdir bin
	mv $@ bin
Guess_fastq_format: $(SRC_DIR)/Guess_fastq_format.cc $(SRC_DIR)/tools.h $(ALL_FILES)
	@echo compiling: $@
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LINKS) 
	@echo compiling: $@ done
	mv $@ bin
Quality_trimmer: $(SRC_DIR)/Quality_trimmer.cc $(SRC_DIR)/tools.h $(ALL_FILES) 
	@echo compiling: $@
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LINKS)
	@echo compiling: $@ done
	mv $@ bin

Adaptor_trimmer.tgz: $(SRC_PKG_LIST)
	rm -rf Adaptor_trimmer.tgz
	rm -rf Adaptor_trimmer-*
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
	rm -rf bin
