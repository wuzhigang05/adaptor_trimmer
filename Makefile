PROG := Adaptor_trimmer
.PHONY: all
all: $(PROG) Guess_fastq_format Quality_trimmer
targets=$(PROG) Guess_fastq_format Quality_trimmer
CXX := g++

INC = -I SeqAn1.3 -I Boost1.50

DEBUG_FLAGS = -O0 -g3 
RELEASE_FLAGS = -O3 
NOASSERT_FLAGS = -DNDEBUG

CXXFLAGS = $(RELEASE_FLAGS) $(INC)
DEBUG_CXXFLAGS = -ggdb $(INC)
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
$(PROG): $(SRC_DIR)/$(PROG).cc $(SRC_DIR)/Fasta_reader.h $(SRC_DIR)/seq.h $(ALL_FILES) 
	@echo compiling: $@
	$(CXX) $(CXXFLAGS) -o $@ $< $(LINKS) 
	@echo compiling: $@ done
	rm -rf bin
	mkdir bin
	mv $@ bin
debug: $(SRC_DIR)/Adaptor_trimmer.cc $(SRC_DIR)/Fasta_reader.h $(SRC_DIR)/seq.h $(ALL_FILES) 
	$(CXX) $(DEBUG_CXXFLAGS)  -o Adaptor_trimmer $(SRC_DIR)/Adaptor_trimmer.cc $(LINKS) 
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

$(PROG).tgz: $(SRC_PKG_LIST)
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
	rm -f Alignment_log*
	rm -f with_5_adaptor no_5_adaptor
	rm -rf Adaptor_trimmer.dSYM
	rm -rf Quality_trimmer.dSYM
	rm -rf Guess_fastq_format.dSYM

.PHONY: test
PREFIX := ../
test:
	dir = test
	$(shell rm -rf $(dir))
	$(shell mkdir $(dir)) 
	$(shell cd $(dir))
#	dir=$(shell `pwd`)
	@echo $(dir)
	$(shell cp $(PREFIX)bin/Adaptor_trimmer .)
	@echo "Test Adaptor_trimmer dynamic programming mode (take input from STDIN)"
	cat $(PREFIX)data/adaptor_test_data.fastq $(PREFIX)data/adaptor_test_data.fastq | ./Adaptor_trimmer -I -o with_5_adaptor -n no_5_adaptor  -5 IamasINGLEADAPT -3 IAMARiGHTADAPTOR -f fastq -l 0 -r 0
	@echo "Test Adaptor_trimmer dynamic programming mode (take input from STDIN) done ..."
#	@echo Test Adaptor_trimmer dynamic programming mode (take input from file)
#	./Adaptor_trimmer -I -o with_5_adaptor -n no_5_adaptor -i data/adaptor_test_data.fastq data/adaptor_test_data.fastq  -5 IamasINGLEADAPT -3 IAMARiGHTADAPTOR -f fastq
#	@echo Test Adaptor_trimmer dynamic programming mode (take input from file) done ...
#	@echo Test Adaptor_trimmer IUPAC mode
#	cat data/AS10.fastq | ./Adaptor_trimmer  -I -5 GYGCASCAGKCGMGAAW -o with_5_adaptor -n no_5_adaptor -U -f fastq
#	cat data/AS10.fastq | ./Adaptor_trimmer  -I -5 G[CT]GCA[CG]CAG[GT]CG[CA]GAA[AT] -o with_5_adaptor -n no_5_adaptor -U -f fastq
#	@echo Test Adaptor_trimmer IUPAC mode done ...
#	@echo Test Adaptor_trimmer using leading and tailing bases mode
#	./Adaptor_trimmer  data/adaptor_test_data.fastq -H 12 -t 4 -o with_5_adaptor -f fastq
#	@echo Test Adaptor_trimmer using leading and tailing bases mode done
#	@echo Test Guess_fastq_format
#	./Guess_fastq_format data/FS2.fastq
#	@echo Test Guess_fastq_format done ...
#	@echo Test Quality_trimmer
#	./Quality_trimmer -f "fastq-sanger" data/AS10.fastq -c 20 -l 100  -s seqs_lessthan20.fastq >seqs_nolessthan20.fastq
#	@echo Test Quality_trimmer done ...
#	@echo Test demultiplexing mode
#	cat data/multiplexing.fastq | ./Adaptor_trimmer -I -5 ACGCCGCAgagtttgatcntggctcag -5 ACGTGTTAgagtttgatcntggctcag -f fastq -o with_5_adaptor
#	@echo Test demultiplexing mode done ...
