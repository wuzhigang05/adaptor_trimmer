Introduction
=============
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
==============
  1. allowing mismatch and indel;
  2. support IUPAC adaptor trimming;
  3. support case-insensitive
  4. support cut leading and tailing bases
  5. support simple regular expression adaptor trimming

Example Usage
==============
1. test Adaptor_trimmer fixed primer mode
  1. take input from stdin arbitrary number of fastq using pipe and cut both 5' and 3' adaptors 
     and force exact match both for 5' adaptor (via -l option) and 3' adaptor (via -r option).
     write sequences with adaptor found to with_5_adaptor and sequences with no adaptor found 
     to no_5_adaptor, both of which will be write to STDOUT by default. ``
     $ cat data/adaptor_test_data.fastq data/adaptor_test_data.fastq | Adaptor_trimmer -I -o with_5_adaptor -n no_5_adaptor  -5 IamasINGLEADAPT -3 IAMARiGHTADAPTOR -f fastq -l 0 -r 0``
  2. same as above but take input from arbitrary number of files using default editing distance, 
     which is 20% of respective adaptor length. This is the parameter most of user should use. 
     Don't worry about we are setting the eidt distance too high, if there are several alignments 
     between adaptor and sequence meeting the requirement, program will always report the best 
     alignment, which means if there is exact match. The program will only report the exact 
     match position. $ Adaptor_trimmer -I -o with_5_adaptor -n no_5_adaptor -i data/adaptor_test_data.fastq data/adaptor_test_data.fastq  -5 IamasINGLEADAPT -3 IAMARiGHTADAPTOR -f fastq
2. test Adaptor_trimmer IUPAC mode
  1. cut the 5' adaptor with case-insensitive (-I) and IUPAC (-U) mode ON. $ cat data/BS1.fastq | Adaptor_trimmer  -I -5 GAGTTTGATCNTGGCTCAG  -o with_5_adaptor -n no_5_adaptor -U -f fastq
  2. Adaptor_trimmer also support simple regular expression only allowing use of square brackets 
     [] to denote alternative nucleotides. $  cat data/BS1.fastq | Adaptor_trimmer  -I -5 GAGTTTGATC[ACGT]TGGCTCAG  -o with_5_adaptor -n no_5_adaptor -U -f fastq
3. Performance test: comparison of dynamic programming mode and regular expression mode
  1. IUPAC mode, which internally using regular expression
   $ time Adaptor_trimmer data/FS11.fastq -I -5 TGGAGGGCAAGTCTGGTG  -o with_5_adaptor -n no_5_adaptor  -f fastq -U
    real 0m0.723s
  2. dynamic programming mode
   $ time Adaptor_trimmer data/FS11.fastq -I -5 TGGAGGGCAAGTCTGGTG  -o with_5_adaptor -n no_5_adaptor  -f fastq -l 0
    real 0m6.700s
  CONCLUSION: if you just wanna using the exact match, then I recommend you using the IUPAC mode 
            because it's 10 times faster than the dynamic program way.
4. test Adaptor_trimmer using leading and tailing bases mode:
  1. read input from file
   $ Adaptor_trimmer  data/adaptor_test_data.fastq -H 12 -t 4 -o with_5_adaptor -f fastq
  2. read input from STDIN (pipe) 
   $ cat data/adaptor_test_data.fastq | Adaptor_trimmer -H 12 -t 4 -o with_5_adaptor -f fastq
