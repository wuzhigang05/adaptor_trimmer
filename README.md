**Introduction**
================

Adaptor trimming is an essential task before doing any serious analysis. 
Most of current available adaptor trimming tools, which do not allow mismatches and indels. 
However, sequencing error is unavoidable in most, if not all, of sequencing platforms. 
Using these adaptor trimming tools will leave adaptors carrying sequencing 
error untrimmed. This undoubtedly will cause artifically and non-biologically meaningful bias, 
which is undesirable. Additionally, current adaptor trimming  tools seems also having trouble to
handle sequences having ambigous IUPAC nucleotides(one position have several alternative 
nucleotides). They are very widely used in studies aims to investigate microbial 
community diversity. If ambiguous IUPAC nucleotide only occur in single position, one can 
circumvent the limit of not allowing mismatch and indel exerted by current adaptor trimming tools 
by running the adaptor trimming tools couple times with different combinations of adaptors. 
However, if ambiguous IUPAC nucleotides occurs in multiple positions, such as the universal IUPAC
primer used to amply Archea 16S contains 5 IUPAC nucleotide position 
(GYGCASCAGKCGMGAAW which is equivalent to G[CT]GCA[CG]CAG[GT]CG[CA]GAA[AT]). it would be extremely 
cumbersome to trim (at least 2^^5 combination). We need a versatile adaptor trimming tool that 
is able to handling sequencing error and ambiguous IUPAC nucleotides to cope with the rapid 
advancement of sequecning technology. adaptor_trimmer is a trimming tool impletemented in C++ 
that handles these two problems efficicently. It has three modes: 1) the dynamic programmming mode, 
which is based on dynamic programming, handles sequence contains mismatches and indels 
efficiently; 2) the IUPAC mode, which is based on regular expression search engine, processes 
ambiguous IUPAC adaptor efficiently; in the other mode one can simply cut any specified number 
of heading and tailing bases. Generally, the tool is very straightfoward to use. We have 
carefully chosen the default value for all optional options so that it requires minimal effort 
from the user. But also we tried to give the user the flexibility of customizing the tool to 
tail their need by making these options available. For example, it takes both FASTA and FASTQ 
formats and you can use Adaptor_trimmer as a typical UNIX command tool which can take I/O 
by piping. Most importantly, adaptor_trimmer is fast and competes most avaible tools.

**Main feature**
----------------

  1. allowing mismatch and indel;
  2. supporting IUPAC adaptor trimming;
  3. supporting case-insensitive
  4. supporting cut leading and tailing bases
  5. supporting simple regular expression adaptor trimming

**Example usage**
-----------------

** Test Adaptor_trimmer dynamic programming mode. **

  * Here we instruct the program to take input from arbitrary number of FASTQs using pipe, 
    cut both 5' and 3' adaptors and force it to use exact match for both 5' adaptor (via -l 
    option) and 3' adaptor (via -r option). We also instruct the program to write sequences 
    with adaptor being found to with_5_adaptor and sequences with no adaptor being found to 
    no_5_adaptor, both of which will be write to STDOUT by default. For illustration purpose, 
    here I am catting from two same FASTQs, you can of course cat from different FASTQs in 
    reality.

        cat data/adaptor_test_data.fastq data/adaptor_test_data.fastq | ./Adaptor_trimmer -I -o with_5_adaptor -n no_5_adaptor  -5 IamasINGLEADAPT -3 IAMARiGHTADAPTOR -f fastq -l 0 -r 0

  * Same as above but taking input from an arbitrary number of files (-i) and trimming adaptors 
    using default paramters (editing distance), which is 20% of respective adaptor length. 
    **This is the parameter most of user should use.** Don't worry about we are setting the 
    eidt distance too high here, if there are several alignments between adaptor and sequence 
    meeting the requirement, program will always **ONLY** report the best alignment. 

        ./Adaptor_trimmer -I -o with_5_adaptor -n no_5_adaptor -i data/adaptor_test_data.fastq data/adaptor_test_data.fastq  -5 IamasINGLEADAPT 
-3 IAMARiGHTADAPTOR -f fastq

** Test Adaptor_trimmer IUPAC mode. **

  * cut the 5' adaptor with case-insensitive (-I) and IUPAC (-U) mode ON.
    
        cat data/AS10.fastq | ./Adaptor_trimmer  -I -5 GYGCASCAGKCGMGAAW -o with_5_adaptor -n no_5_adaptor -U -f fastq
 
  * Adaptor_trimmer also support simple regular expression only allowing use of square brackets 
    [] to denote alternative nucleotides. 

        cat data/AS10.fastq | ./Adaptor_trimmer  -I -5 G[CT]GCA[CG]CAG[GT]CG[CA]GAA[AT] -o with_5_adaptor -n no_5_adaptor -U -f fastq  
 
** Performance test: comparison of dynamic programming mode and regular expression mode. **

  * IUPAC mode, which internally using regular expression 
        
        time ./Adaptor_trimmer data/FS2.fastq -I -5 TGGAGGGCAAGTCTGGTG  -o with_5_adaptor -n no_5_adaptor  -f fastq -U
        
        real 0m0.216s

  * dynamic programming mode 
        
        time ./Adaptor_trimmer data/FS2.fastq -I -5 TGGAGGGCAAGTCTGGTG  -o with_5_adaptor -n no_5_adaptor  -f fastq -l 0
        
        real 0m3.347s

  * CONCLUSION: if you want to use the exact match anyway, then I recommend you using the IUPAC mode because it's 10 times 
    faster than the dynamic program way.

** Test Adaptor_trimmer using leading and tailing bases mode. **

  * Below I instruct Adaptor_trimmer to take input from file (note I didn't use specify -i here, the program will automatically take data/adaptor_test_data.fastq as input file), cut off 12 leading bases (-H) and tailing 4 bases (-t) and write the trimmed reads to with_5_adaptor(-o).
        
        ./Adaptor_trimmer  data/adaptor_test_data.fastq -H 12 -t 4 -o with_5_adaptor -f fastq

  * Below I instruct Adaptor_trimmer to take  input from STDIN (pipe) and do the same thing as above.
        
        cat data/adaptor_test_data.fastq | ./Adaptor_trimmer -H 12 -t 4 -o with_5_adaptor -f fastq

**Companion tools**
-----------------

After using Adaptor_trimmer, we already get our clean sequence free of adaptor(s) contamination. Very often in NGS data 
analysis, we want to extract the high quality region from these sequences. To fit this need, we provide two accompanying 
tools along with Adaptor_trimmer. They are **Guess_fastq_format and Quality_trimmer. **


**Guess_fastq_format**

  Guess_fastq_format as its indicates is able to tell the format of FASTQ in terms of fastq-sanger, fastq-solex, 
  fastq-illumina. If your sequence is generated from Roche 454 then its highly possible that your sequence is in 
  fastq-sanger format. You can refer to Nucleic Acids Res. 2010 April; 38(6): 1767–1771. for more information 
  about FASTQ format. Below is a typical usage of Guess_fastq_format:
    
    ./Guess_fastq_format data/FS2.fastq 


**Quality_trimmer**
  
  Given a sequence quality cutoff Quality_trimmer is able to extract the sub-maximum array of sequence positions out.
  In addition to that, one also need to specify the format of the fastq being examined. As mentioned above, you can 
  get this information by running Guess_fastq_format program. Below is a typical usage of Quality_trimmer. I instructed 
  the program to extract the maximum contiguous region with quality score >= 20 (-c) from data/AS10.fastq file which is in 
  fastq-sanger format (-f) and write the extracted result which meet the length requirement of 100 (-l) to 
  seqs_nolessthan20.fastq and the sequences failed length filtering to seqs_lessthan20.fastq (-s).

    ./Quality_trimmer -f "fastq-sanger" data/AS10.fastq -c 20 -l 100  -s seqs_lessthan20.fastq >seqs_nolessthan20.fastq



**Using three tools as a whole**

    ./Guess_fastq_format data/AS10.fastq

    cat data/AS10.fastq | ./Adaptor_trimmer  -I -5 GYGCASCAGKCGMGAAW  -n no_5_adaptor -U -f fastq |  ./Quality_trimmer -f "fastq-sanger" -c 20 -l 100  -s seqs_lessthan20.fastq >seqs_nolessthan20.fastq

**Other possible usage cases**

 * Demultiplexing

  One can use Adaptor_trimmer to demultiplexing. Below is an example usage:

    for i in ATCTGACTGG ATCTAGCTGC ....; do ./Adaptor_trimmer -i multiplexed1.fastq multiplexed2.fastq -5 $i -o $i.fastq -l 1 -n no_barcode_found.fastq -f fastq; done
  
  Code explanation: You can specify barcode as many as you would in the ecllipse position. Above I set the edit distance to be 1. 
  You can set your own according to your situation. By the way, all the example showing above regarding Adaptor_trimmer are 
  dealing with fastq files, if you want to handle fasta file, just leave off the "-f fastq" option.

**Suggestions and Bugs**

Please send your comments and bug reports to zhigang.wu@email.ucr.edu


