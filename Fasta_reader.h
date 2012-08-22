/*
 * =====================================================================================
 *
 *       Filename:  Fasta_reader.h
 *
 *    Description:  This is a file that can read in fasta either from STDIN or
 *                  a file.
 *
 *        Version:  1.0
 *        Created:  08/09/2012 15:42:59
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zhigang Wu (UCR), zhigang.wu@email.ucr.edu
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef   __FASTA_READER_H__
#define   __FASTA_READER_H_
#include <string>
#include "seq.h"
#include <err.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
//#include <seqan/basic.h>
//#include <seqan/file.h>
#include <string.h>
/* 
 * read from file pointer, one record per time
 */
//class Seqrecord
//{
//public:
//  std::string ID;
//  std::string Seq;
//  Seqrecord(std::string id, std::string seq):ID(id), Seq(seq) {}
//  Seqrecord():ID(""), Seq("") {}
//
//};

int Fasta_reader(FILE * fp, Seqrecord & seqobj)
{

  const int bufsize = 1024;
  char buf [bufsize];
  char next_char;
/* 
 * do not do any error checking here, you are supposed to check the 
 */
  next_char = fgetc(fp);
  ungetc(next_char, fp);
/*   
 *   if first char is not > then return false
 */
  if (next_char != '>')
  {
//    std::cerr << __FILE__ 
//      << " First character is not '>'" 
//      <<__LINE__ << std::endl;
//    std::cerr << "Any lines following this line will not be procssed" << std::endl;
    return 0;
  }
  fgets(buf, bufsize, fp);
  char * name = strtok(buf+1, " \t\n"); // skip the '>'
  seqobj.ID = name; 
  std::string seq;
  while((next_char = fgetc(fp)) != EOF)
  {
    if (next_char == '>')
    { 
      ungetc(next_char, fp); //its very important to put is back
      break;
    }

    seq += next_char;
  }
/* 
 * replace the white space
 */
 std::string::size_type pos;
  while( (pos = seq.find("\n")) != std::string::npos)
    seq.erase(pos, 1);
  seqobj.Seq = seq;
  return 1;
      
}


/* 
 * read from file 
 */
int Fastq_reader(FILE *fp, Seqrecord & seqobj)
{
  const int bufsize = 1024;
  char buf [bufsize];
  char next_char;
/* 
 * do not do any error checking here, you are supposed to check before here
 */
  next_char = fgetc(fp);
  ungetc(next_char, fp);
/*   
 *   if first char is not > then return false
 */
  if (next_char != '@')
  {
//    errx(0, __FILE__ " First character is not @" __LINE__ );
//    std::cerr  << " First character is not @" << std::endl;
//    std::cerr << "Any lines following this line won't be procssed" << std::endl;
    return 0;
  }
  fgets(buf, bufsize, fp);
  char * name = strtok(buf+1, " \t\n"); // skip the '>'
  seqobj.ID = name; 
  while(fgets(buf, bufsize, fp) != NULL)
  {
    if (buf[strlen(buf)-1] == '\n')             /* remove the \n char */
      buf[strlen(buf)-1] = '\0';
    seqobj.Seq = buf;                           /* get the seq */
    if (fgets(buf, bufsize, fp) != NULL)        /* get the + sign */
    {
      if (buf[0] != '+')
        errx(0, __FILE__ " First character is not '+' " \
            "Using other tool to format it first.", __LINE__ );
    }

    else
    {
      std::cerr << "In complete sequenc record found, "
        << seqobj.ID << " without + sign" __FILE__ << __LINE__ << std::endl;
      std::cerr << "Any lines following this line won't be procssed" << std::endl;
      return 0;
    }

    if (fgets(buf, bufsize, fp) != NULL)        /* get the qual string */
    {
      if (buf[strlen(buf)-1] == '\n')             /* remove the \n char */
        buf[strlen(buf)-1] = '\0';
      seqobj.Qual = buf;
    }
    else
    {
      std::cerr << "In complete sequenc record found, "
        << seqobj.ID << " without quality string" __FILE__ << __LINE__ << std::endl;
      std::cerr << "Any lines following this line won't be procssed" << std::endl;
      return 0;
    }
    return 1;                                   /* read seqrecord success, singal */
  }
}

int Fasta_reader(FILE *fp, Seqrecord & seqobj, 
                const char * format)
{
  if (strcmp(format, "fasta") == 0)
    Fasta_reader(fp, seqobj);
  else if (strcmp(format, "fastq") == 0)
    Fastq_reader(fp, seqobj);
  else
    errx(0, __FILE__, "invalid format specified %d", __LINE__);
}
void writeSeqrecord(std::ostream & os, 
                    const Seqrecord & seqobj, 
                    const char * format = "fasta")
{
  if (strcmp(format, "fasta") == 0)
  {
    os << ">" << seqobj.ID << std::endl
              << seqobj.Seq << std::endl;
  }
  else if(strcmp(format, "fastq") == 0)
  {
    os << "@" << seqobj.ID << std::endl
              << seqobj.Seq << std::endl
              << "+" << std::endl
              << seqobj.Qual << std::endl;
//    seqan::writeRecord(os, seqobj.ID, seqobj.Seq, seqobj.Qual, seqan::Fastq());
  }
  else
  {
    std::cerr << "invalid format. valid formats: fasta, fastq" << std::endl;
    exit(1);
  }
}
#endif

