Introduction
------------

Adaptor trimming is an essential task before doing any serious analysis. 
Most of current available adaptor trimming tools, which do not allow mismatches and indels. 
However, sequencing error is unavoidable in most, if not all, of sequencing platforms. 
Using these adatpor trimming tools will leave adaptors carrying sequencing 
error untrimmed. This undoubtedly will cause artifically and non-biologically meaningful bias, 
which is undesirable. Additionally, current adatpor trimming  tools seems also having trouble to
handle sequences having ambigous IUPAC nucleotides(one position have several alternative 
nucleotides). They are very widely used in studies aims to investigate microbial 
community diversity. If ambiguous IUPAC nucleotide only occur in single position, one can 
circumvent the limit of no mismatch and indel exerted by current adaptor trimming tools by 
running the adaptor trimming tools couple times with different combinations of adatpors. 
However, if ambiguous IUPAC nucleotides occurs in multiple positions, like the universal IUPAC
primer used to amply Archea 16S contains 5 IUPAC nucleotide position, it would be extremely 
curbersome to trim (at least 2^^5 combination). We need a versatile adatpor trimming tool that 
is able to handling sequencing error and ambiguous IUPAC nucleotides to cope with the rapid 
advancement of sequecning technology. Adatpor_trimmer is a trimming tool impletemented in C++ 
that handles these two problems efficicently. It has three modes: 1) the fixed adaptor mode, 
which is based on dynamic programming, handles sequence contains mismatches and indels 
efficiently; 2) the IUPAC mode, which is based on regular expression search engine, processes 
ambiguous IUPAC adaptor efficiently; the other mode is simply can cut any specified number of 
heading and tailing bases. The tool as a whole is very easy to use because we have carefully 
chosen the default value for all optional options but also gives the user the flexibility of 
customizing the tool to tail their need by making the options available. For example, 
It supports both fasta and fastq formats and you can use Adaptor_trimmer as a UNIX command 
tool by piping. Most importantly, adaptor_trimmer is fast and competes most avaible tools.

Main feature
============

  1. allowing mismatch and indel;
  2. support IUPAC adaptor trimming;
  3. support case-insensitive
  4. support cut leading and tailing bases
  5. support simple regular expression adaptor trimming

Example Usage
=============

**Test Adaptor_trimmer fixed primer mode.**

  * take input from stdin arbitrary number of fastq using pipe and cut both 5' and 3' adaptors 
    and force exact match both for 5' adaptor (via -l option) and 3' adaptor (via -r option).
    write sequences with adaptor found to with_5_adaptor and sequences with no adaptor found 
    to no_5_adaptor, both of which will be write to STDOUT by default.

        cat data/adaptor_test_data.fastq data/adaptor_test_data.fastq | Adaptor_trimmer -I -o with_5_adaptor -n no_5_adaptor  -5 IamasINGLEADAPT -3 IAMARiGHTADAPTOR -f fastq -l 0 -r 0

  * same as above but take input from arbitrary number of files using default editing distance, 
    which is 20% of respective adaptor length. This is the parameter most of user should use. 
    Don't worry about we are setting the eidt distance too high, if there are several alignments 
    between adaptor and sequence meeting the requirement, program will always report the best 
    alignment, which means if there is exact match. The program will only report the exact 
    match position.

        Adaptor_trimmer -I -o with_5_adaptor -n no_5_adaptor -i data/adaptor_test_data.fastq data/adaptor_test_data.fastq  -5 IamasINGLEADAPT 
-3 IAMARiGHTADAPTOR -f fastq

** Fold long FASTA/Q lines and remove FASTA/Q comments: **

  * a
  * b

        seqtk seq -Cl60 in.fa > out.fa

* Convert multi-line FASTQ to 4-line FASTQ:

        seqtk seq -l0 in.fq > out.fq

* Reverse complement FASTA/Q:

        seqtk seq -r in.fq > out.fq

* Extract sequences with names in file `name.lst`, one sequence name per line:

        seqtk subseq in.fq name.lst > out.fq

* Extract sequences in regions contained in file `reg.bed`:

        seqtk subseq in.fa reg.bed > out.fa

* Mask regions in `reg.bed` to lowercases:

        seqtk seq -M reg.bed in.fa > out.fa

* Subsample 10000 read pairs from two large paired FASTQ files (remember to use the same random seed to keep pairing):

        seqtk sample -s100 read1.fq 10000 > sub1.fq
        seqtk sample -s100 read2.fq 10000 > sub2.fq

* Trim low-quality bases from both ends using the Phred algorithm:

        seqtk trimfq in.fq > out.fq

* Trim 5bp from the left end of each read and 10bp from the right end:

        seqtk trimfq -b 5 -e 10 in.fa > out.fa
