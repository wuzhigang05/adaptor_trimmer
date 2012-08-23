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
#include "tools.h"
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
    ("sample,s", po::value<int>()->default_value(1000, "1000"), 
     "Number of sequence records"
     "used for guess format");
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
  CharString id, query, qual;
  string format = "";
  int count = 0;
  if (vm["input"].as<vector<string> >().size() == 1 &&
      vm["input"].as<vector<string> >()[0] == string("stdin") )
  {
    if (!isatty(fileno(stdin)))
    {
      seqan::RecordReader<std::istream, seqan::SinglePass<> > reader(std::cin);
      while(readRecord_stdin(id, query, qual, reader, "fastq") && 
          count++ < vm["sample"].as<int>() )
      {
        if (qual == "")
          errx(1, "Input file is not a valid fastq file");
        if (guess_format_sanger(string(toCString(qual))))
        {
          format = "fastq_sanger";
          break;
        }
        if (guess_format_solexa_illumina(string(toCString(qual))))
        {
          format = "fastq_solexa_illumina";
          break;
        }
      }
      if(format != string(""))
        cerr << "format is: " << format << endl;
      else
      {
        cerr << "We are unable to determine it's format based on "
          << vm["sample"].as<int>() << " sequences." << endl
          << "Its highly possible the input sequences is in fastq-sanger"
          << " format with quality score very high range from 26-40." 
          << "we will go ahead use this fastq-sanger format"<< endl;
        format =  "fastq-sanger";
      }
    }
    else
    {
      cerr << "The program is hanging and waiting for input from STDIN" << endl;
      exit(1);
    }
  }
  else
  {
    string filename;
    for (int i = 0; i < vm["input"].as<vector<string> >().size(); ++i)
    {
      filename = vm["input"].as<vector<string> >()[i]; 
      SequenceStream seqStream(filename.c_str());
      while(readRecord_Stream(id, query, qual, seqStream, "fastq") && 
          count++ < vm["sample"].as<int>() )
      {
        if (qual == "")
          errx(1, "Input file is not a valid fastq file");
        if (guess_format_sanger(string(toCString(qual))))
        {
          format = "fastq_sanger";
          break;
        }
        if (guess_format_solexa_illumina(string(toCString(qual))))
        {
          format = "fastq_solexa_illumina";
          break;
        }
      }
      if(format != string(""))
        cerr << "format is: " << format << endl;
      else
      {
        cerr << "We are unable to determine it's format based on "
          << vm["sample"].as<int>() << " sequences." << endl
          << "Its highly possible the input sequences is in fastq-sanger"
          << " format with quality score very high ranging from 26 to 40." 
          << "We will go ahead use fastq-sanger as its format"<< endl;
        format =  "fastq-sanger";
      }
    }
  }
  return 0;

}

