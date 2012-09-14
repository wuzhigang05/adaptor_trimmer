/*
 * =====================================================================================
 *
 *       Filename:  Guess_fastq_format.cc 
 *
 *    Description:  Give a fastq file, guess whether it's in fastq-sanger or 
 *                  fastq-solex and fastq-illumina format.
 *
 *        Version:  1.0
 *        Created:  08/20/2012 10:33:23
 *       Revision:  none
 *       Compiler:  gcc
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
#include <err.h>
#include "readRecord.h"
//#include "tools.h"
#define PROG_NAME "Guess_fastq_format"
using namespace std;
using namespace seqan;
namespace po = boost::program_options;
const boost::regex sanger_expr ("[!\"#$%&'()*+,-./0123456789:]");
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
    format =  "fastq-sanger";
  }
  return format;
}

int main (int argc, char * argv[])
{
  po::options_description desc(
"\n" PROG_NAME " Part of NGS C++ Toolkit by Zhigang Wu zhigang.wu@email.ucr.edu.\n"
"Give a fastq file, guess whether it's in fastq-sanger or fastq-solex and fastq-illumina\n"
"format. The main thing program tries to determine is if its a fastq-sanger format or not,\n"
"becase all others (fastq-illumina, fastq-Solexa) are using offset of 64. However, \n"
"fastq-sanger is using 33 as offset. If you are using Roche 454 sequencing, it's highly \n"
"possible that the fastq is in fastq-sanger format. \n"
"Please see Nucleic Acids Res. 2010 April; 38(6): 1767–1771. for more about fastq format.\n"
   );
  desc.add_options()
    ("help,h", "print this help page")
    ("input,i", po::value<vector<string> >()->default_value(vector<string>(1, "stdin"), 
                                                            "STDIN"), 
     "input file(s), has to be in fastq format")
    ("format,f", po::value<string>(), 
     "the format of input file, valid formats include fastq-sanger and fastq-illumina "
     "fastq-solexa")
//    ("sample,s", po::value<int>()->default_value(1000, "1000"), 
//     "Number of sequence records"
//     "used for guess format")
    ;
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
    cerr << "For directing STDIN you can either use '|' or <infile" << endl;
    cerr << desc << endl;
    exit(1);
  }
// check whether format specified is valid or not
  string format = "";
  int count = 0;
  if (vm["input"].as<vector<string> >().size() == 1 &&
      vm["input"].as<vector<string> >()[0] == string("stdin") )
  {
    if (!isatty(fileno(stdin)))
      guess_format(std::cin);
    else
      errx(1, "The program is hanging and waiting for input from STDIN");
  }
  else
  {
    string filename;
    for (int i = 0; i < vm["input"].as<vector<string> >().size(); ++i)
    {
      filename = vm["input"].as<vector<string> >()[i]; 
      ifstream stream(filename.c_str());
      if (!stream.good())
        errx(1, "cannot open file %s for read", filename.c_str());
      guess_format(stream);
    }
  }
  return 0;

}

