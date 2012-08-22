/*
 * =====================================================================================
 *
 *       Filename:  Quality_trimmer.cc
 *
 *    Description:  Given a format and a quality cutoff, the program will return a longest
 *                  contiguous segment of sequence, in which the minimum quality score for
 *                  each nucleotide is >= specified quality cutoff.
 *
 *        Version:  1.0
 *        Created:  08/20/2012 10:33:23
 *       Revision:  none
 *       Compiler:  gcc
 *           Test:  Comparison of my program and DynamicTrim.pl yields similar result.
 *           Quality_trimmer -f "fastq-sanger" <test.fastq >hello.fastq
 *           perl ~/self_shell/DynamicTrim.pl -h 19 test.fastq
 *                  diff hello.fastq test.fastq.trimmed
 *
 *
 *         Author:  Zhigang Wu (UCR), zhigang.wu@email.ucr.edu
 *   Organization:  
 *
 * =====================================================================================
 */
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <algorithm>
#include <map>
#include <numeric>
#include <typeinfo>
#include "tools.h"
#include <err.h>
#define all(c) c.begin(), c.end()
#define P(v) cerr << #v << v << endl;
#define tr(container, it) for (typeof(container.begin()) it = container.begin(); it != container.end(); it++)

#include "tools.h"
#define PROG_NAME "Quality_trimmer"
using namespace std;
using namespace seqan;
namespace po = boost::program_options;
const boost::regex sanger_expr ("[!\"#$%&'()*+,-./0-9:]");
// its very unlikely for fastq-sanger have a score over 40 (ASCII code 'I').
// Based on this, we constructed below solexa_illumina_expr, if any qual string
// match below regular expression, then we will determine it as sollexa or
// illumina format. However we are not insterested in further differentiating them 
// because the offset for both of them is same (64).
const boost::regex solexa_illumina_expr ("[J-Z\\\[\\]^_`a-z]");
boost::match_results<std::string::const_iterator> what;
bool guess_format_sanger(const string & qual)
{
  return boost::regex_search(qual.begin(), qual.end(), what, sanger_expr);
}
bool guess_format_solexa_illumina(const string & qual)
{
  return boost::regex_search(qual.begin(), qual.end(), what, solexa_illumina_expr);
}

/* below is like a generator function in python, which will read one record
 * per time*/
int readRecord_stdin(CharString & id, CharString & query, CharString & qual, 
                  seqan::RecordReader<std::istream, seqan::SinglePass<> > & reader,
                  const char * format)
{
  if (strcmp(format, "fastq") == 0)
  {
    if(readRecord(id, query, qual, reader, Fastq()) == 0)
      return 1;                                 /* indicates read record success */
    return 0;
  }
  else if (strcmp(format, "fasta") == 0)
  {
    if (readRecord(id, query, reader, Fasta()) == 0)
    { 
      qual = ""; /* set qual default value to empty for fasta format */
      return 1;
    }
    return 0;
  }
  else
    errx(1, "bad format specified here in file: ", __FILE__, " Line: ", __LINE__);
}
/* below is like a generator function in python, which will read one record
 * per time*/
int readRecord_Stream(CharString & id, CharString & query, CharString & qual, 
                  seqan::SequenceStream & seqStream,
                  const char * format)
{
  if (strcmp(format, "fastq") == 0)
  {
    if(readRecord(id, query, qual, seqStream) == 0)
      return 1;                                 /* indicates read record success */
    return 0;
  }
  else if (strcmp(format, "fasta") == 0)
  {
    if (readRecord(id, query, seqStream) == 0 )
    { 
      qual = ""; /* set qual default value to empty for fasta format */
      return 1;
    }
    return 0;
  }
  else
    errx(1, "bad format specified here in file: ", __FILE__, " Line: ", __LINE__);
}

pair<vector<int>::iterator, vector<int>::iterator> find_maximum_subarray( vector<int> & c )
{

  typedef vector<int>::iterator Tci;
  Tci start_previous = c.begin();               /* used to hold the to be returned start iterator */
  Tci start = c.begin();
  Tci end;                                      /* used to hold the to be returned end iterator */
  Tci tmp = c.begin()+1;
  int max_sofar = 0;
  int accum_sum_previous = 0;
  int accum_sum_current = 0;
  for(; tmp != c.end(); ++tmp)
  {
    accum_sum_current = accumulate(start, tmp, 0);
    accum_sum_previous = accumulate(start, --tmp, 0);
    tmp++;
    if (accum_sum_current < accum_sum_previous) /* case encoutner negative */
      start = tmp;                              
    if (accum_sum_current > max_sofar)
    {
      start_previous = start;
      end = tmp;
    }
    max_sofar = max(accum_sum_current, max_sofar);
  }
  return make_pair<Tci, Tci>(start_previous, end);
}

int char2int( char c)
{
  return int(c);
}
vector<int> get_quality_array(const CharString & qual)
{
  string str_qual = string(toCString(qual));
  vector<int> v;
  v.reserve(str_qual.length());
  transform(str_qual.begin(), str_qual.end(), back_inserter(v), char2int);
  return v;
}

//ofstream tie_file2ostream(string filename)
//{
//  ofstream out;
//  out.open(filename.c_str());
//  if (!out.is_open())
//    errx(1, "Error in open file %s" , __FILE__, __LINE__);
//  return out;
//}
int main (int argc, char * argv[])
{
  po::options_description desc(
"\n" PROG_NAME " Part of NGS C++ Toolkit by Zhigang Wu zhigang.wu@email.ucr.edu.\n"
"Given a fastq file and a quality cutoff, the program will gurantee to return a longest\n"
"contiguous segment of sequence, in which the minimum quality score for each nucleotide \n"
"is >= specified quality cutoff. If you know your fastq format clearly in which format. \n"
"That's please specify it. If you don't know, using Guess_fastq_fromat to guess it first.\n" 
"The format matters becuase fastq-sanger format use an offset of 33 while all others \n"
"(fastq-illumina, fastq-Solexa) are using an offset of 64. If you are using Roche 454 \n"
"sequencing, it's highly possible that the fastq is in fastq-sanger format. Please see \n"
"Nucleic Acids Res. 2010 April; 38(6): 1767–1771. for more information about fastq format.\n"
   );
  desc.add_options()
    ("help,h", "print this help page")
    ("input,i", po::value<vector<string> >()->default_value(vector<string>(1, "stdin"), 
                                                            "STDIN"), 
     "input file(s), has to be in fastq format")
    ("format,f", po::value<string>(), 
     "the format of input file, valid formats include fastq-sanger and fastq-illumina "
     "fastq-solexa")
    ("length,l", po::value<int>()->default_value(0, "0"), 
     "length >= specifed value will be output")
    ("output-long,o", po::value<string>(), "output file. default [STDOUT]")
    ("output-short,s", po::value<string>(), "output file below the length cutoff. "
     "default [STDOUT]")
    ("cutoff,c", po::value<int>()->default_value(20, "20"), 
     "the quality cutoff. Each position of the trimmed reads will have a quality score "
     ">= cutoff specified here");

    
  po::positional_options_description p;
  p.add("input", -1);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
      options(desc).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    cerr << desc << endl;
    exit(1);
  }
//  check whether there is any valid input specified
  if (vm["input"].as<vector<string> >().size() ==1 && 
      vm["input"].as<vector<string> >()[0] == string("stdin") 
      && isatty(fileno(stdin)))
  {
    cerr << "You have to specify at least one input file OR direct STDIN to this program" 
      << endl;
    cerr << "For directint STDIN you can either use '|' or <infile" << endl;
    cerr << desc << endl;
    exit(1);
  }
// check whether format specified is valid or not
  if (vm.count("format"))
  {
    string format = vm["format"].as<string>();
    if (!(format == string("fastq-sanger") || format == string("fastq-illumina") || 
          format == string("fastq-solexa")))
      errx(1, "in-valid format specified");
  }
  else
    errx(1, "you have to specify the format for the input fastq file");


// setup output stream                                                
  ofstream out;
  ofstream out_short;
  if (vm.count("output_long"))
  {
    out.open(vm["output_long"].as<string>().c_str());
    if (!out.is_open())
    {
      errx(1, "Error in open file %s" , __FILE__, __LINE__);
    }
  }
  if (vm.count("output_short"))
  {
    out_short.open(vm["output_short"].as<string>().c_str());
    if (!out_short.is_open())
    {
      errx(1, "Error in open file %s" , __FILE__, __LINE__);
    }
  }
  std::ostream & OS_long = vm.count("output_long") ? out : std::cout;
  std::ostream & OS_short = vm.count("output_short") ? out_short: std::cout;

// setup the offset
  int offset = (vm["format"].as<string>() == string("fastq-sanger") ) ? 33 : 64;

// start processing
  vector<int> vi;
  vector<int> vi_offset_negated;
  pair<vector<int>::iterator, vector<int>::iterator> tmp;
  int clip_start;
  int clip_end;
  int length;
  CharString id, query, qual;
  vector<int> stat;
  /* from stdin */
  if (vm["input"].as<vector<string> >().size() == 1 &&
      vm["input"].as<vector<string> >()[0] == string("stdin") )
  {
    if (!isatty(fileno(stdin)))
    {
      seqan::RecordReader<std::istream, seqan::SinglePass<> > reader(std::cin);
      while(readRecord_stdin(id, query, qual, reader, "fastq"))
      {
        if (qual == "")
          errx(1, "Input file is not a valid fastq file");
        vi = get_quality_array(qual);
        tr(vi, it)                              /* for each value negate the offset */
          *it = *it - offset;
        tr(vi, it)                              /* prepare the 1 and -1 array */
          *it = *it >= vm["cutoff"].as<int>() ? 1 : -1;
        tmp = find_maximum_subarray(vi);
        clip_start = tmp.first-vi.begin();
        clip_end = tmp.second-vi.begin();
        length = clip_end - clip_start;
        qual = infix(qual, clip_start, clip_end);
        query = infix(query, clip_start, clip_end);
        stat.push_back(length);
        if (length >= vm["length"].as<int>())
          writeRecord(OS_long, id, query, qual, Fastq());
        else
          writeRecord(OS_short, id, query, qual, Fastq());
      }
    }
  }
  else
  {
    string filename;
    for (int i = 0; i < vm["input"].as<vector<string> >().size(); ++i)
    {
      filename = vm["input"].as<vector<string> >()[i]; 
      SequenceStream seqStream(filename.c_str());
      while(readRecord_Stream(id, query, qual, seqStream, "fastq") )
      {
        if (qual == "")
          errx(1, "Input file is not a valid fastq file");
        vi = get_quality_array(qual);
        tr(vi, it)                              /* for each value negate the offset */
          *it = *it - offset;
        tr(vi, it)                              /* prepare the 1 and -1 array */
          *it = *it >= vm["cutoff"].as<int>() ? 1 : -1;
        tmp = find_maximum_subarray(vi);
        clip_start = tmp.first-vi.begin();
        clip_end = tmp.second-vi.begin();
        length = clip_end - clip_start;
        qual = infix(qual, clip_start, clip_end);
        query = infix(query, clip_start, clip_end);
        stat.push_back(length);
        if (length >= vm["length"].as<int>())
          writeRecord(OS_long, id, query, qual, Fastq());
        else
          writeRecord(OS_short, id, query, qual, Fastq());
      }
    }
  }
  cerr << "mean length is:   " << mean<int>(stat) << endl;
  cerr << "median length is: " << median<int>(stat) << endl;
  return 0;

}


