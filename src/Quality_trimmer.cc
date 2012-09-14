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
#include <fstream>
#include <map>
#include <numeric>
#include <typeinfo>
#include "tools.h"
#include <err.h>
#include "readRecord.h"
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

template <typename TStream>
string guess_format(TStream & stream)
{
  Seqrecord seq;
  std::string format = "";
  while(ReadRecord(stream, seq, FASTQ()))
  {
    if (seq.Qual == "")
      errx(1, "Input file is not a valid fastq file: %s has empty qual str", seq.Seq.c_str());
    if (guess_format_sanger(string(seq.Qual.c_str())))
    {
      format = "fastq_sanger";
      break;
    }
    if (guess_format_solexa_illumina(string(seq.Qual.c_str())))
    {
      format = "fastq_solexa_illumina";
      break;
    }
  }

  if(format != string(""))
    cerr << "format is: " << format << endl;
  else
  {
    cerr << "We are unable to determine it's format. " << endl
      << "Its highly possible the input sequences is in fastq-sanger"
      << " format with quality score very high range from 26-40." 
      << "we will go ahead use this fastq-sanger format"<< endl;
    format =  "fastq_sanger";
  }
  return format;
}

template <typename TStream>
void process(TStream & stream, const po::variables_map & vm, std::string & format, 
              ostream & OS_long, ostream & OS_short, vector<int> & stat)
{
// setup the offset
  int offset = 0;
  if (format == string("fastq_sanger"))
     offset = 33;
  else
     offset = 64;
  vector<int> vi;
  vector<int> vi_offset_negated;
  pair<vector<int>::iterator, vector<int>::iterator> tmp;
  int clip_start;
  int clip_end;
  int length;

  Seqrecord seq;

  while(ReadRecord(stream, seq, FASTQ()))
  {
    if (seq.Qual == "")
      errx(1, "Input file is not a valid fastq file: %s has empty qual str", seq.ID.c_str());
    vi = get_quality_array(seq.Qual);
    tr(vi, it)                              /* for each value negate the offset */
      *it = *it - offset;
    tr(vi, it)                              /* prepare the 1 and -1 array */
      *it = *it >= vm["cutoff"].as<int>() ? 1 : -1;
    tmp = find_maximum_subarray(vi);
    clip_start = tmp.first-vi.begin();
    clip_end = tmp.second-vi.begin();
    length = clip_end - clip_start;
    seq.Qual = seq.Qual.substr(clip_start, clip_end - clip_start);
    seq.Seq = seq.Seq.substr(clip_start, clip_end - clip_start);
    stat.push_back(length);
    if (length >= vm["length"].as<int>())
      WriteRecord(OS_long, seq, FASTQ());
    else
      WriteRecord(OS_short, seq, FASTQ());
  }
}


int main (int argc, char * argv[])
{
  po::options_description desc(
"\n" PROG_NAME " Part of NGS C++ Toolkit by Zhigang Wu zhigang.wu@email.ucr.edu.\n"
"Given a fastq file and a quality cutoff, the program will gurantee to return a longest\n"
"contiguous segment of sequence, in which the minimum quality score for each nucleotide \n"
"is >= specified quality cutoff. You have to specify the fastq format, which you can get\n"
"by running Guess_fastq_fromat. The format matters becuase fastq-sanger format use an\n"
"offset of 33 while all others (fastq-illumina, fastq-Solexa) are using an offset of 64.\n"
"If you are using Roche 454 sequencing, it's highly possible that the fastq is in \n"
"fastq-sanger format. Please see Nucleic Acids Res. 2010 April; 38(6): 1767–1771.\n"
"for more information about fastq format.\n"
   );
  desc.add_options()
    ("help,h", "print this help page")
    ("input,i", po::value<vector<string> >()->default_value(vector<string>(1, "stdin"), 
                                                            "STDIN"), 
     "input file(s), has to be in fastq format")
    ("format,f", po::value<string>(), 
     "the format of input file, valid formats include fastq_sanger and fastq_illumina "
     "fastq_solexa")
    ("length,l", po::value<int>()->default_value(0, "0"), 
     "length >= specifed value will be output")
    ("output_long,o", po::value<string>(), "output file. default [STDOUT]")
    ("output_short,s", po::value<string>(), "output file below the length cutoff. "
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
//  if (vm.count("format"))
//  {
//    string format = vm["format"].as<string>();
//    if (!(format == string("fastq-sanger") || format == string("fastq-illumina") || 
//          format == string("fastq-solexa")))
//      errx(1, "in-valid format specified");
//  }
//  else
//    errx(1, "you have to specify the format for the input fastq file");


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

  std::vector<int> stat;
// start processing
  /* from stdin */
  if (vm["input"].as<vector<string> >().size() == 1 &&
      vm["input"].as<vector<string> >()[0] == string("stdin") )
  {
    if (!isatty(fileno(stdin)))
      if (vm.count("format"))
      {
        std::string format = vm["format"].as<std::string>();
        process(std::cin, vm, format, OS_long, OS_short, stat);
      }
      else
        errx(1, "You have to specify format if you take input from STDIN");
    else
      errx(1, "The program is hanging and waiting for input from STDIN");
  }
  else
  {
    string filename;
    for (int i = 0; i < vm["input"].as<vector<string> >().size(); ++i)
    {
      filename = vm["input"].as<vector<string> >()[i]; 
      ifstream stream;
      stream.open(filename.c_str());
      if (!stream.good())
        errx(1, "cannot open file %s for read", filename.c_str());

      std::string format;
      if (!vm.count("format"))
      {
        format = guess_format(stream);
        stream.close();
        stream.open(filename.c_str());
      }
      else
        format = vm["format"].as<std::string>();
      if (!stream.good())
        errx(1, "cannot open file %s for read", filename.c_str());
      process(stream, vm, format, OS_long, OS_short, stat);
    }
  }
  cerr << "mean length is:   " << mean<int>(stat) << endl;
  cerr << "median length is: " << median<int>(stat) << endl;
  return 0;

}


