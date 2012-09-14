/*
 * =====================================================================================
 *
 *       Filename:  readRecord.cpp
 *
 *    Description:  common interface for reading seqs in both fastq and fasta format
 *
 *        Version:  1.0
 *        Created:  09/11/2012 14:22:59
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zhigang Wu (UCR), zhigang.wu@email.ucr.edu
 *   Organization:  UC Riverside
 *
 * =====================================================================================
 */

#ifndef READRECROD_H
#define READRECROD_H
#include <iostream> 
#include <string>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <err.h>
#include "seq.h"
class FASTA {};
class FASTQ {};

template <typename TStream>
// return >0 for successfully read a seq, 0 for failure.
//int ReadRecord(TStream & stream, TString &id, TString & seq, TString & qual)
int ReadRecord(TStream & stream, Seqrecord & seqobj, const FASTA &)
{
  char c;
  std::string line;
  seqobj.Seq.clear();
  if (stream.eof())
    return 0;
  if ((c = stream.get()) != '>')
  {
    std::cerr << "ERROR: first character of line is not '>'" << std::endl;  
    return 0;
  }
  else
    getline(stream, seqobj.ID, '\n');

  while(!stream.eof())
  {
    c = stream.peek();
    if (c == '>' || stream.eof())
      return 1;
    else
    {
      getline(stream, line, '\n');
      seqobj.Seq += line;
    }
  }
}

template <typename TStream>
int ReadRecord(TStream & stream, Seqrecord & seqobj, const FASTQ &)
{
  char c;
  std::string line;
  seqobj.Seq.clear(); // reinitialize the Seq
  seqobj.Qual.clear(); // reinitialize the Qual

  if (stream.eof())
    return 0;
  if ((c = stream.get()) != '@')
  {
    std::cerr << "ERROR: first character of line is not '@'" << std::endl;  
    return 0;
  }
  else
    getline(stream, seqobj.ID, '\n');


/* read in the seq, it doestn't matter whether its line wrapped or not
 */
  while(!stream.eof())
  {
    c = stream.peek();
    if (c == '+')
    {
      getline(stream, line, '\n');              /* discard this line */
      break;
    }
    else
    {
      getline(stream, line, '\n');   
      seqobj.Seq += line;
    }
  }

/* read in qual portion, the difficulty is to differentiate wether the @ at the beginning
 * of is a quality str or a new start of seq
 */
  while(!stream.eof())
  {
    c = stream.peek();
// this line is new fastq start rather than quality str
    if (c == '@' && seqobj.Seq.length() == seqobj.Qual.length() || stream.eof())
      return 1;
    else
    {
      getline(stream, line, '\n');
      seqobj.Qual += line;
    }
  }
}


// Wrapper for two ReadRecord in readRecord.h
template <typename TStream>
int ReadRecord(TStream &stream, Seqrecord & seq,
                const char * format)
{
  if (strcmp(format, "fastq") == 0)
    return ReadRecord(stream, seq, FASTQ());   /* 1 indicates read record success; 0 for fail */
  else if (strcmp(format, "fasta") == 0)
  {
    seq.Qual = ""; /* set qual default value to empty for fasta format */
    return ReadRecord(stream, seq, FASTA());   /* 1 indicates read record success; 0 for fail */
  }
  else
    errx(1, "bad format specified here in file: ", __FILE__, " Line: ", __LINE__);
}

/* write FASTA without line wrapping
 */
template <typename TStream>
void WriteRecord (TStream &stream, const Seqrecord & seqobj, const FASTA &)
{
        stream << '>' << seqobj.ID << std::endl
          << seqobj.Seq << std::endl;
}

/* write FASTQ without line wrapping
 */
template <typename TStream>
void WriteRecord (TStream &stream, const Seqrecord & seqobj, const FASTQ &)
{
    stream << '@' << seqobj.ID << std::endl
    << seqobj.Seq << std::endl
    << '+' << std::endl 
    << seqobj.Qual << std::endl; 
}

/* Wrapper for FASTA/Q without line wrapping
 */
template <typename TStream>
void WriteRecord (TStream &stream, const Seqrecord & seqobj, const char * format)
{
  if(strcmp(format, "fasta") == 0)
    WriteRecord(stream, seqobj, FASTA());
  else if(strcmp(format, "fastq") == 0)
    WriteRecord(stream, seqobj, FASTQ());
  else
    errx(1, "bad format specified here in file: ", __FILE__, " Line: ", __LINE__);
}

//int main (int argc, char** argv)
//{
//  assert (argc ==3);
//  std::ifstream in(argv[1]);
//  Seqrecord seqobj;
//  time_t start = clock();
//  double sec;
//  int rec = 0;
//  if (strcmp(argv[2], "fasta") == 0 )
//  {
//    if (in.good())
//    {
//      while(ReadRecord(in, seqobj, FASTA()))
//      {
//        rec++;
//        std::cout << '>' << seqobj.ID << std::endl
//          << seqobj.Seq << std::endl;
//      }
//    }
//    else
//    {
//      while(ReadRecord(std::cin, seqobj, FASTA()))
//      {
//        rec++;
//        std::cout << '>' << seqobj.ID << std::endl
//          << seqobj.Seq << std::endl;
//      }
//    }
//  }
//  else if (strcmp(argv[2], "fastq") == 0)
//  {
//    if (in.good())
//    {
//      while(ReadRecord(in, seqobj, FASTQ()))
//      {
//        rec++;
//        std::cout << '@' << seqobj.ID << std::endl
//          << seqobj.Seq << std::endl 
//          << '+' << std::endl
//          <<seqobj.Qual << std::endl;
//      }
//    }
//    else
//    {
//      while(ReadRecord(std::cin, seqobj, FASTQ()))
//      {
//        rec++;
//        std::cout << '@' << seqobj.ID << std::endl
//          << seqobj.Seq << std::endl
//          << "+" << std::endl 
//          << seqobj.Qual << std::endl;
//      }
//    }
//  }
//
////  sec =  (float) (time(NULL) - start);
//  sec =  (double) (clock() - start) / CLOCKS_PER_SEC ;
//  std::cerr << "Saw " << rec << " lines in " << std::setprecision(6) << sec << " seconds." << std::endl;
//}
#endif // READRECROD_H
