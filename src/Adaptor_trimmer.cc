/*
 * =====================================================================================
 *
 *       Filename:  Adaptor_trimmer.cc
 *
 *    Description:  This programs aims to provide a flexiable tool for doing a very basic
 *                  task in NGS sequence analysis, which is removing the adaptor sequences.
 *
 *        Version:  1.0
 *        Created:  08/11/2012 13:38:49
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Zhigang Wu (UCR), zhigang.wu@email.ucr.edu
 *   Organization:  Department of Botany and Plant Sciences, UC Riverside 
 *
 * =====================================================================================
 */
#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <fstream>
#include <string>
#include <iomanip>
#include <seqan/basic.h>
#include <boost/program_options.hpp>
#include <math.h>
#include <iomanip>
#include <err.h>
#include "Fasta_reader.h"
#include "readRecord.h"
#include <algorithm>
#include <iterator>
#include <vector>
#include <map>
#include <sys/stat.h>
#include <stdio.h>
#include <boost/regex.hpp>

using namespace std;
using namespace seqan;
namespace po = boost::program_options;

typedef CharString TSequence;
typedef Align<TSequence, ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;
typedef Iterator<TRow>::Type TIterator;

typedef enum {UNDEFINED_5, FOUND_5_ADAPT, NO_5_ADAPT} adapt_5;
typedef enum {UNDEFINED_3, FOUND_3_ADAPT, NO_3_ADAPT } adapt_3;
static int have_3_adaptor_only_counter = 0;
static int have_5_adaptor_only_counter = 0;
static int have_5_and_3_adaptor_counter = 0;
static int no_5_and_3_adaptor_counter = 0;
static int sum = 0;                            /* the sum of above four, which should = total */
static int total = 0;                          /* total number of sequence processed */
static std::string format;
/* 
 * return 1 if found adaptor
 * return 0 if no adaptor found
 */
int judge_adaptor(adapt_5 & d5, adapt_3 & d3)
{
  if (d5 != UNDEFINED_5 && d3 != UNDEFINED_3)
  {
    switch(d5)
    {
      case NO_5_ADAPT:
        switch(d3)
        {
          case NO_3_ADAPT:
            ++no_5_and_3_adaptor_counter;
            return 0;
          case FOUND_3_ADAPT:
            ++have_3_adaptor_only_counter;
            return 1;
        }
      case FOUND_5_ADAPT:
        switch(d3)
        {
          case NO_3_ADAPT:
            ++have_5_adaptor_only_counter;
            return 1;
          case FOUND_3_ADAPT:
            ++have_5_and_3_adaptor_counter;
            return 1;
        }
    }
  }
  else if (d5 != UNDEFINED_5 )                  /* case d3 == UNDEFINED_3 */
  {
    switch(d5)
    {
      case NO_5_ADAPT:
        ++no_5_and_3_adaptor_counter;
        return 0;
      case FOUND_5_ADAPT:
        ++have_5_adaptor_only_counter;
        return 1;
    }
  }
  else if (d3 != UNDEFINED_3)                   /* case d5 == UNDEFINED_5 */
  {
    switch(d3)
    {
      case NO_3_ADAPT:
        ++no_5_and_3_adaptor_counter;
        return 0;
      case FOUND_3_ADAPT:
        ++have_3_adaptor_only_counter;
        return 1;
    }
   
  }

}


/* count how many gaps in the specified row of pairwise alignment
 * if i = 0 for the up row, 1 for below row
 */
int num_of_gaps(TAlign & Ali, int i = 0)
{
  TIterator start = iter(row(Ali, i), 0);
  TIterator end   = iter(row(Ali, i), length(row(Ali, i)));
  unsigned gaps = 0;
  while( start != end)
  {
    if (isGap(start))
      ++gaps; 
    ++start;
  }
  return gaps;
}

vector<int> get_regex_match_StartEnd(const string & str, 
                                     const string & pattern,
                                     bool CaseInSensitive)
{
  boost::regex e;
  if (CaseInSensitive)
    e = boost::regex(pattern, boost::regex::basic|boost::regex::icase);
  else
    e = boost::regex(pattern, boost::regex::basic);
  std::string::const_iterator start, end, initial;
  start = str.begin();
  initial = start;
  end = str.end();
  vector<int> v;
  boost::match_results<std::string::const_iterator> what;
  if (boost::regex_search(start, end, what, e))
  {
    v.push_back(what[0].first - initial); 
    v.push_back(what[0].second - initial);
  }
  return v;
}

//int regEx_match(CharString & id, CharString & query, CharString & qual, const CharString & db, 
int regEx_match(Seqrecord & seq, const string & db, 
                int mismatch, const string & direction, bool CaseInSensitive, 
                ostream & os_align)
{
  vector <int> start_end_pos = get_regex_match_StartEnd(seq.Seq, 
                                                        db, 
                                                        CaseInSensitive );
  if( start_end_pos.size() == 0)
    return 0;
  else
  {
    int start_pos = start_end_pos[0];
    int end_pos = start_end_pos[1];
//    std::string match_str = infix(seq.Seq, start_pos, end_pos);
    std::string match_str = seq.Seq.substr(start_pos, end_pos - start_pos);
    if (direction == string("five")) 
    {
      os_align << ">" << seq.ID << "_five" << endl;
      os_align << "Distance from Seq head to 5' adaptor end is: " << end_pos << endl;
      os_align << match_str << endl;

//      query = suffix(query, end_pos);
      seq.Seq = seq.Seq.substr(end_pos);
      if (seq.Qual !=  "")        /* if seq.Qual == "" means this is fasta sequence */
//        seq.Qual = suffix(seq.Qual, end_pos);
        seq.Qual = seq.Qual.substr(end_pos);
      return 1;
    }
    else if (direction == string("three"))
    {
      os_align << ">" << seq.ID << "_three" << endl;
      os_align << "Distance from Seq tail to 3' adaptor start is: " 
        << length(seq.Seq) - start_pos << endl;
      os_align << match_str << endl;

//      query = prefix(query, start_pos); 
      seq.Seq = seq.Seq.substr(0, start_pos); 
      if (seq.Qual != "")        /* if seq.Qual == "" means this is fasta sequence */
//        qual = prefix(qual, start_pos);
        seq.Qual = seq.Qual.substr(0, start_pos); 
      return 1;
    }
    else
      errx(1, "bad format specified here in file: ", __FILE__, " Line: ", __LINE__);
  }

}

//int do_alignment(CharString & id, CharString & query, CharString & qual, const CharString & db, 
int do_alignment(Seqrecord & seq, const string & db, 
                  int mismatch, const string & direction, bool CaseInSensitive, 
                  ostream & os_align)
/*   note if direction == five, then db has to be the five prime adaptor sequence,
 *   if direction  == three, then db has to be the three primer adaptor sequence
 */
{
  TAlign Ali;
  resize(rows(Ali), 2);
  /*  make a local copy of query and db if CaseInSensitive is set.
   *  This seems cubersome. But We are doing this in order to overcome
   *  the side effect of toUpper function, which modifies the string in
   *  place rather than return a copy of original string*/
//  string query_uppercase = query;
  string query_uppercase = seq.Seq;
  string db_uppercase = db;

  if (CaseInSensitive)
  {
    toUpper(query_uppercase);
    toUpper(db_uppercase);
  }
/*   
 *   make sure db are not empyt sequence do it outside of this function
 */
  assignSource(row(Ali, 0), query_uppercase);
  assignSource(row(Ali, 1), db_uppercase);
  Score <int> scoring(1, -1, -1, -1);
  int target_score = length(db_uppercase) - 2 * mismatch; 
  LocalAlignmentEnumerator<Score<int>, Unbanded> enumerator(scoring, target_score);
  if(nextLocalAlignment(Ali, enumerator))
  {
    unsigned start_pos = clippedBeginPosition(row(Ali, 0));
    unsigned d_end =clippedEndPosition(row(Ali, 1)); 
    int d_gaps = num_of_gaps(Ali, 1);
    int gaps = num_of_gaps(Ali, 0);
//    unsigned end_pos = clippedEndPosition(row(Ali, 0)) - gaps + length(db) - d_end + d_gaps; /* inclusive */
    unsigned end_pos = clippedEndPosition(row(Ali, 0)) - gaps + db.length() - d_end + d_gaps; /* inclusive */
    //                                                          ------------------------------
    //                                                          if the best alignment does not extend to the end of adaptor then add these up
    //                              below is the case
    //                        startpos = 12, end_pos = 23, d_start = 0, d_end = 11
    //            0          10        20
    //            TCAGATCTTCTATGGAGGGCAAGCNGNCCGAGCGAGCCGCGGTAATCCCAGCTCCAATAGCG
    //     query(read):       TGGAGGGCAAG
    //                        |||||||||||
    //     db(adaptor):       TGGAGGGCAAG 
    //     Full Adaptor:      TGGAGGGCAAGTCTGGTG                                               
    if (direction == string("five")) 
    {
      os_align << ">" << seq.ID << "_five" << endl;
      os_align << ">" << seq.ID << "_five" << endl;
      os_align << "Best alignment score is: " << getScore(enumerator) << endl;
      os_align << "Distance from Seq head to 5' adaptor end is: " << end_pos << endl;
      os_align << Ali << endl;

//      seq.Seq = suffix(seq.Seq, end_pos);
      seq.Seq = seq.Seq.substr(end_pos);
      if (seq.Qual !=  "")        /* if seq.Qual == "" means this is fasta sequence */
//        seq.Qual = suffix(seq.Qual, end_pos);
        seq.Qual = seq.Qual.substr(end_pos);
      return 1;
    }
    else if (direction == string("three"))
    {
      os_align << ">" << seq.ID << "_three" << endl;
      os_align << "Best alignment score is: " << getScore(enumerator) << endl;
      os_align << "Distance from Seq tail to 3' adaptor start is: " 
//        << length(query) - start_pos << endl;
        << seq.Seq.length() - start_pos << endl;
      os_align << Ali << endl;

//      seq.Seq = prefix(seq.Seq, start_pos); 
      seq.Seq = seq.Seq.substr(0, start_pos); 
      if (seq.Qual != "")        /* if seq.Qual == "" means this is fasta sequence */
//        seq.Qual = prefix(seq.Qual, start_pos);
        seq.Qual = seq.Qual.substr(0, start_pos); 
      return 1;
    }
    else
      errx(1, "bad format specified here in file: ", __FILE__, " Line: ", __LINE__);
  }
  else
    return 0;                                   /* no adaptor found */
}

/* 
 * wrapper for calling do_alignment, this function could have been incorporated
 * into do_alignment function, but do_alignment function already been so long.
 * I do not want to have a very long function. Briefly, if detected adaptor,
 * d5 will be set up correspondingly.
 */
//void examine_5(CharString & id, CharString & query, CharString & qual, 
//    const CharString & adaptor, int num_mismatches,  const string & direction, 
void examine_5(Seqrecord & seqobj, 
    const string & adaptor, int num_mismatches,  const string & direction, 
    bool CaseInSensitive, ostream & os_align, adapt_5 & d5, bool IUPAC)
{
  if (!IUPAC)
//    if (do_alignment(id, query, qual, adaptor, num_mismatches, direction, 
    if (do_alignment(seqobj, adaptor, num_mismatches, direction, 
          CaseInSensitive, os_align))
      d5 = FOUND_5_ADAPT;
    else
      d5 = NO_5_ADAPT;
  else
    if (regEx_match(seqobj, adaptor, num_mismatches, direction, 
          CaseInSensitive, os_align))
      d5 = FOUND_5_ADAPT;
    else
      d5 = NO_5_ADAPT;

}

template <typename TIter>
void print(std::ostream & os, TIter  begin, TIter  end, const char * message = "", const char * sep = "\n")
{
  if (strcmp(message, "") != 0)
    os << message << ": " << sep;
  typedef typename std::iterator_traits<TIter>::value_type T;
  copy(begin, end, ostream_iterator<T>(os, sep)); 
  os << std::endl;
}

/* similar to above function, this function only examines whether the query has
 * significant similarity with the 3' adaptor
 */
//void examine_3(CharString & id, CharString & query, CharString & qual, 
void examine_3(Seqrecord & seqobj, 
    const string & adaptor, int num_mismatches,  const string & direction, 
    bool CaseInSensitive, ostream & os_align, adapt_3 & d3, bool IUPAC)
{
  if(!IUPAC)
//    if (do_alignment(id, query, qual, adaptor, num_mismatches, direction, 
    if (do_alignment(seqobj, adaptor, num_mismatches, direction, 
          CaseInSensitive, os_align))
      d3 = FOUND_3_ADAPT;
    else
      d3 = NO_3_ADAPT;
  else
//    if (regEx_match(id, query, qual, adaptor, num_mismatches, direction, 
    if (regEx_match(seqobj, adaptor, num_mismatches, direction, 
          CaseInSensitive, os_align))
      d3 = FOUND_3_ADAPT;
    else
      d3 = NO_3_ADAPT;
}
  

void report_paramter_setting (const po::variables_map & vm, 
                              int width = 55, ostream & out = cerr)
{
  out << "Here is what you set:" << endl;
  out << setw(width) << left << "Input file:" 
       << setw(width) << left;
  copy(vm["input"].as<vector<string> >().begin(), vm["input"].as<vector<string> >().end(),
      ostream_iterator<string>(out, " "));
  out << endl;

  out << setw(width) << left << "Format is:" << setw(width) << left << format << std::endl;
  if (vm.count("five"))
  {
    out << setw(width) << left << "5' adaptor:" << setw(width) << left;
    print(out,
          vm["five"].as<std::vector<std::string> >().begin(),
          vm["five"].as<std::vector<std::string> >().end(),
          "", " ");
  }
  if (vm.count("three"))
    out << setw(width) << left << "3' adaptor:" 
         << setw(width) << left << vm["three"].as<string>()<< endl;
  if (vm.count("out_with_adaptor"))
    out << setw(width) << left << "output file for seqs with adaptor:" 
         << setw(width) << left << vm["out_with_adaptor"].as<string>()<< endl;
  if (vm.count("out_no_adaptor"))
    out << setw(width) << left << "output file for seqs with no adaptor found:" 
         << setw(width) << left << vm["out_no_adaptor"].as<string>()<< endl;
  if (vm.count("out_align") && ! (vm.count("header") || vm.count("tail")))
    out << setw(width) << left << "out file for alignment between adaptor(s) and seq:" 
         << setw(width) << left << vm["out_align"].as<string>()<< endl;
  if (vm.count("case-insensitive"))
     out << setw(width) << left << "case-insensitive is:" 
         << setw(width) << left << "ON"<< endl;
  else
      out << setw(width) << left << "case-insensitive is:" 
         << setw(width) << left << "OFF"<< endl;
  if (vm.count("IUPAC"))
     out << setw(width) << left << "IUPAC is:" 
         << setw(width) << left << "ON"<< endl;
  else
      out << setw(width) << left << "IUPAC is:" 
         << setw(width) << left << "OFF"<< endl;
  if (vm.count("head"))
     out << setw(width) << left << "head is:" 
         << setw(width) << left << vm["head"].as<int>() << endl;
  if (vm.count("tail"))
     out << setw(width) << left << "tail is:" 
         << setw(width) << left << vm["tail"].as<int>()<< endl;
}

/* using seqan interface write seq to ostream
 */
//void write_seq2stream(ostream & os, const CharString & id,
//                    const CharString & seq, const CharString & qual)
//{ 
//  if (seq != "")              /* if seq is empty means this is a solely adaptor sequence */
//  {                                             /* Here we skipped printing these empty sequences */
//    if( qual != "" )
//      writeRecord(os, id, seq, qual, Fastq());
//    else
//      writeRecord(os, id, seq, Fasta());
//  }
//}

void get_time_stamp(char * buffer)
{
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  strftime (buffer,80,"D%dH%HM%MS%S.txt",timeinfo);
}
/* get the allowed number of mismatches and indels for 5'
 */
//int get_allowed_mismatch_indel(const po::variables_map & vm, const CharString & adaptor,
int get_allowed_mismatch_indel(const po::variables_map & vm, const string & adaptor,
                               const char * direction)
{
  int allowed_gap_mismatches;
//  int adaptor_len = length(adaptor);
  int adaptor_len = adaptor.length();
  float pct = vm["percent"].as<float>();
  if (strcmp(direction, "five") == 0)
    if (vm.count("five-mismatch"))
    {
      if ( vm["five-mismatch"].as<int>() >= adaptor_len / 2 ) /* throw warning exceed cutoff*/
        cerr << "Warning: Your allowed number of mismatches and indels for 5' is too high!!"
          << endl << "The result may contains false positive, sequence without adaptors will"
          << " be cut at any positions." << endl << " Set this value less than" 
          << adaptor_len / 2 << " through -l option" << endl;
      allowed_gap_mismatches = vm["five-mismatch"].as<int>();
    }
    else
      allowed_gap_mismatches = int(ceil(adaptor_len * pct)); /* by default 0.2 shortest adapt */
  else if (strcmp(direction, "three") == 0)
    if (vm.count("three-mismatch"))
    {
      if ( vm["three-mismatch"].as<int>() >= adaptor_len / 2 ) /* throw warning exceed cutoff*/
        cerr << "Warning: Your allowed number of mismatches and indels for 3' is too high!!"
          << endl << "The result may contains false positive, sequence without adaptors will"
          << " be cut at any positions." << endl << " Set this value less than" 
          << adaptor_len / 2 << " through -l option" << endl;
      allowed_gap_mismatches = vm["three-mismatch"].as<int>();
    }
    else
      allowed_gap_mismatches = int(ceil(adaptor_len * pct)); /* by default 20% shortest adapt */
  else
    errx(1, "bad format specified here in file: ", __FILE__, " Line: ", __LINE__);

  return allowed_gap_mismatches;

}

map<std::string, std::string> IUPAC_map()
{
  map<std::string, std::string> m;
  m["R"] = "[AG]";
  m["Y"] = "[CT]";
  m["S"] = "[GC]";
  m["W"] = "[AT]";
  m["K"] = "[GT]";
  m["M"] = "[AC]";
  m["B"] = "[CGT]";
  m["D"] = "[AGT]";
  m["H"] = "[ACT]";
  m["V"] = "[ACG]";
  m["N"] = "[ACGT]";
  return m;
}

//CharString translate_IUPAC(const string & db)
string translate_IUPAC(const string & db)
{
  string str_db = db;
  map<std::string, std::string> m = IUPAC_map();
  typedef map<std::string, std::string> Mtype;
  Mtype::iterator it = m.begin();
  for(; it != m.end(); ++it)
    str_db = boost::regex_replace(str_db, boost::regex(it->first), it->second);
//  string new_db = str_db;
  return str_db;
}

void Judge_Write(ostream & OS_with_adaptor, ostream & OS_no_adaptor, Seqrecord & seq, 
                  adapt_5 & d5, adapt_3 & d3, const char * format)
{
    if (judge_adaptor(d5,d3))
//        write_seq2stream(OS_with_adaptor, id, query, qual);
      WriteRecord(OS_with_adaptor, seq, format);
    else /* case whole the seq length < specified to be cut leading length*/
//        write_seq2stream(OS_no_adaptor, id, query, qual);
      WriteRecord(OS_no_adaptor, seq, format);
}

template <typename TStream>
int cut_leading_tailing_stream( TStream &stream, 
                      const po::variables_map & vm, 
                      ostream & OS_with_adaptor,
                      ostream & OS_no_adaptor,
                      adapt_5 & d5, adapt_3 & d3,
                      const char * format)
{
  if (vm.count("head") && vm.count("tail"))
  {
    int start_pos = vm["head"].as<int>();
    int end_pos = vm["tail"].as<int>();
    Seqrecord seq;
    while(ReadRecord(stream, seq, format))
    {
      ++total;
//      id = seq.ID; query = seq.Seq; qual = seq.Qual;  
      if ( start_pos < length(seq.Seq) - 1 && /* check at least < seq length */
          length(seq.Seq) - end_pos - start_pos > vm["length-cutoff"].as<int>())
        /* length(query)-endpos = tobe cut position */
      {
//        query = infix(query, start_pos, length(query) - end_pos);
        seq.Seq = seq.Seq.substr(start_pos, length(seq.Seq) - end_pos - start_pos);
        if (seq.Qual != "")
//          qual  = infix(qual, start_pos, length(query) - end_pos);
          seq.Qual = seq.Qual.substr(start_pos, length(seq.Qual) - end_pos - start_pos);
        d5 = FOUND_5_ADAPT;
        d3 = FOUND_3_ADAPT;
      }
      Judge_Write(OS_with_adaptor, OS_no_adaptor, seq, d5, d3, format);
    }
  }
  else if (vm.count("head") )
  {
    int start_pos = vm["head"].as<int>();
    Seqrecord seq;
//    CharString id, query, qual;
    while(ReadRecord(stream, seq, format))
    {
      ++total;
//      id = seq.ID; query = seq.Seq; qual = seq.Qual;  
      if ( start_pos < length(seq.Seq) - 1) /* check at least < seq length */
      {
//        query = suffix(query, start_pos);
        seq.Seq = seq.Seq.substr(start_pos);
        if (seq.Qual != "")
//          qual  = suffix(qual, start_pos);
          seq.Qual = seq.Qual.substr(start_pos);
        d5 = FOUND_5_ADAPT;
      }
      Judge_Write(OS_with_adaptor, OS_no_adaptor, seq, d5, d3, format);
    }
  }
  else if (vm.count("tail") )
  {
    int end_pos = vm["tail"].as<int>();
    Seqrecord seq;
//    CharString id, query, qual;
    while(ReadRecord(stream, seq, format))
    {
      ++total;
//      id = seq.ID; query = seq.Seq; qual = seq.Qual; 
      if ( length(seq.Seq) - end_pos > vm["length-cutoff"].as<int>() ) 
      {
//        query = prefix(query, end_pos);
        seq.Seq = seq.Seq.substr(0, end_pos); 
        if (seq.Qual != "")
//          qual  = prefix(qual, end_pos);
          seq.Qual = seq.Qual.substr(0, end_pos); 
        d3 = FOUND_3_ADAPT;
      }
      Judge_Write(OS_with_adaptor, OS_no_adaptor, seq, d5, d3, format);
    }
  }
}

bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
        return true;
    return false;
} 

void RemovePreExistedFiles(const po::variables_map & vm, const char * format)
{
  for (int i = 0; i < vm["five"].as<std::vector<std::string> >().size(); ++i)
  {
    string five_adaptor = vm["five"].as<std::vector<std::string> >()[i];
    if (vm.count("out_with_adaptor"))
    {
      char filename[200];
      int n = sprintf(filename, "%s%s.%s", five_adaptor.c_str(), vm["out_with_adaptor"].as<string>().c_str(), format);
      if (fileExists(std::string (filename)))
      {
        std::cerr << filename << " exists in your current direcotry" << std::endl;
        if(remove(filename) != 0)
          perror( "Error deleting file");
        else
          std::cerr << "delete " << filename << " successfully" << std::endl;
      }
    }
  }
}

std::map<std::string, std::ostream*> GetMapOfFile2Ofstream (const po::variables_map & vm, const char * format)
{
  std::map<std::string, std::ostream*> m; // m used to contain the std::map of filename to ofstream pointer
  for (int i = 0; i < vm["five"].as<std::vector<std::string> >().size(); ++i)
  {
    std::string five_adaptor = vm["five"].as<std::vector<std::string> >()[i];
    if (vm.count("out_with_adaptor"))
    {
      char filename[200];
      int n = sprintf(filename, "%s%s.%s", five_adaptor.c_str(), vm["out_with_adaptor"].as<string>().c_str(), format);
      if (!m.count(five_adaptor))
      {
        m[five_adaptor] = new ofstream(filename, fstream::app);
      }
      else
        errx(1, "Error: %s already in map!\n", five_adaptor.c_str());
    }
    else
      m[five_adaptor] = &std::cout;
  }
  return m;
}

template <typename TStream>
void read_process_5_3_fastaq(TStream & stream, 
    const std::string & five_adaptor, const std::string & three_adaptor, 
    int five_allowed_mismatch_indel, int three_allowed_mismatch_indel, 
    bool case_insensitive, ostream & OS_alignment, ostream & OS_with_adaptor,
    ostream & OS_no_adaptor, adapt_3 & d3, adapt_5 & d5, bool IUPAC, const char * format)
{
  Seqrecord seq;
//  CharString id, query, qual;
  while(ReadRecord(stream, seq, format))
  {
    ++total;
//    id = seq.ID; query = seq.Seq; qual = seq.Qual;  
    examine_5(seq, five_adaptor, five_allowed_mismatch_indel, "five", 
        case_insensitive, OS_alignment, d5, IUPAC);
    examine_3(seq, three_adaptor, three_allowed_mismatch_indel, "three", 
        case_insensitive, OS_alignment, d3, IUPAC);
    Judge_Write(OS_with_adaptor, OS_no_adaptor, seq, d5, d3, format);
  }
}

template <typename TStream>
void read_process_5_fastaq(TStream & stream, 
//    const CharString & five_adaptor, const CharString & three_adaptor, 
    const std::string & five_adaptor, const std::string & three_adaptor, 
    int five_allowed_mismatch_indel, int three_allowed_mismatch_indel, 
    bool case_insensitive, ostream & OS_alignment, ostream & OS_with_adaptor,
    ostream & OS_no_adaptor, adapt_3 & d3, adapt_5 & d5, bool IUPAC, const char * format)
{
  Seqrecord seq;
//  CharString id, query, qual;
  while(ReadRecord(stream, seq, format))
  {
    ++total;
//    id = seq.ID; query = seq.Seq; qual = seq.Qual;  
//    examine_5(id, query, qual, five_adaptor, five_allowed_mismatch_indel, "five", 
    examine_5(seq,five_adaptor, five_allowed_mismatch_indel, "five", 
        case_insensitive, OS_alignment, d5, IUPAC);
    Judge_Write(OS_with_adaptor, OS_no_adaptor, seq, d5, d3, format);
  }
}


template <typename TStream>
void read_process_3_fastaq(TStream & stream, 
//    const CharString & five_adaptor, const CharString & three_adaptor, 
    const std::string & five_adaptor, const std::string & three_adaptor, 
    int five_allowed_mismatch_indel, int three_allowed_mismatch_indel, 
    bool case_insensitive, ostream & OS_alignment, ostream & OS_with_adaptor,
    ostream & OS_no_adaptor, adapt_3 & d3, adapt_5 & d5, bool IUPAC, const char * format)
{
  Seqrecord seq;
  CharString id, query, qual;
  while(ReadRecord(stream, seq, format))
  {
    ++total;
    id = seq.ID; query = seq.Seq; qual = seq.Qual;  
    examine_3(seq,three_adaptor, three_allowed_mismatch_indel, "three", 
        case_insensitive, OS_alignment, d3, IUPAC);
    Judge_Write(OS_with_adaptor, OS_no_adaptor, seq, d5, d3, format);
  }
}

/* This is the main interface that process sequences from STDIN
 */
template <typename TStream>
void TrimmingSeq_from_stream( TStream & stream,
                             const po::variables_map & vm, 
                             const bool & case_insensitive,
                             adapt_5 & d5, adapt_3 & d3, 
                             ostream & OS_alignment, 
                             bool inspect_5, 
                             bool inspect_3,
                             const char * format)
{
  bool IUPAC = vm.count("IUPAC") ? true : false; 
  ofstream out1;
  ofstream out;

  if (vm.count("out_no_adaptor") )
  {
    std::string filename = vm["out_no_adaptor"].as<string>();
    out1.open(filename.c_str());
    if (!out1.is_open())
      errx(1, "Error in open file %s", filename.c_str());
  }
  std::ostream & OS_no_adaptor = vm.count("out_no_adaptor") ? out1 : std::cout;

  if (vm["five"].as<std::vector<std::string> >().size() <= 1) /* using default OS_with_adaptor */
  {
    std::string five_adaptor = vm.count("five") ? vm["five"].as<std::vector<std::string> >()[0] : "";
    std::string three_adaptor = vm.count("three") ? vm["three"].as<string>() : "";

    if (vm.count("out_with_adaptor"))
    {
      out.open(vm["out_with_adaptor"].as<string>().c_str());
      if (!out.is_open())
        errx(1, "Error in open file %s" , __FILE__, __LINE__);
    }
    std::ostream & OS_with_adaptor = vm.count("out_with_adaptor") ? out : std::cout;

    if (IUPAC)                    /* if IUPAC is true translate the IUPAC letter first */
    {                                             /* after translation will be like */
      five_adaptor = translate_IUPAC(five_adaptor); /*   in   = GWGTTTGAAG */
      three_adaptor = translate_IUPAC(three_adaptor); /* out = G[AT]GTTTGAAG */
    }
    int five_allowed_mismatch_indel = get_allowed_mismatch_indel(vm, five_adaptor, "five");
    int three_allowed_mismatch_indel = get_allowed_mismatch_indel(vm, three_adaptor, "three");
//    CharString id, query, qual;
//    Seqrecord seq;
    if (inspect_5 && inspect_3)
    {
      read_process_5_3_fastaq(stream, five_adaptor, three_adaptor, 
                              five_allowed_mismatch_indel, 
                              three_allowed_mismatch_indel, 
                              case_insensitive, OS_alignment, 
                              OS_with_adaptor, OS_no_adaptor, d3, d5,IUPAC, format);
    }
    else if (inspect_5)
    {
      read_process_5_fastaq(stream, five_adaptor, three_adaptor, 
                              five_allowed_mismatch_indel, 
                              three_allowed_mismatch_indel, 
                              case_insensitive, OS_alignment, 
                              OS_with_adaptor, OS_no_adaptor, d3, d5,IUPAC, format);
    }
    else if (inspect_3)
    {
      read_process_3_fastaq(stream, five_adaptor, three_adaptor, 
                              five_allowed_mismatch_indel, 
                              three_allowed_mismatch_indel, 
                              case_insensitive, OS_alignment, 
                              OS_with_adaptor, OS_no_adaptor, d3, d5,IUPAC, format);
    }
    else
      errx(1, "one of inspect_3 and inspect_5 has to be true");
  }
  else
  {
    std::string three_adaptor = vm.count("three") ? vm["three"].as<string>() : "";
  // m used to contain the std::map of filename to ofstream pointer
    std::map<std::string, std::ostream*> m = GetMapOfFile2Ofstream(vm, format); 
    if (inspect_5 && inspect_3)
    {
      Seqrecord seq;
//      CharString id, query, qual;
      while(ReadRecord(stream, seq, format))
      {
        ++total;
//        id = seq.ID; query = seq.Seq; qual = seq.Qual;  
        std::string five_adaptor;
        for (int i = 0; i < vm["five"].as<std::vector<std::string> >().size(); ++i)
        {
          five_adaptor = vm["five"].as<std::vector<std::string> >()[i];
          if (IUPAC)                    /* if IUPAC is true translate the IUPAC letter first */
            five_adaptor = translate_IUPAC(five_adaptor); /*   in   = GWGTTTGAAG */
          int five_allowed_mismatch_indel = get_allowed_mismatch_indel(vm, five_adaptor, "five");
          examine_5(seq,five_adaptor, five_allowed_mismatch_indel, "five", 
              case_insensitive, OS_alignment, d5, IUPAC);
          if (d5 == FOUND_5_ADAPT) // stop if already found the adaptor
            break;
        }
        if (IUPAC)                              /* need to check 3' adaptor */
          three_adaptor = translate_IUPAC(three_adaptor); /* out = G[AT]GTTTGAAG */
        int three_allowed_mismatch_indel = get_allowed_mismatch_indel(vm, three_adaptor, "three");
        examine_3(seq,three_adaptor, three_allowed_mismatch_indel, "three", 
            case_insensitive, OS_alignment, d3, IUPAC);

        std::ostream & OS_with_adaptor = *m[five_adaptor];     /* overwrites default OS_with_adaptor */
        Judge_Write(OS_with_adaptor, OS_no_adaptor, seq, d5, d3, format);
      }
    }
    else if (inspect_5)
    {
      Seqrecord seq;
//      CharString id, query, qual;
      while(ReadRecord(stream, seq, format))
      {
        ++total;
//        id = seq.ID; query = seq.Seq; qual = seq.Qual;  
        std::string five_adaptor; 
        for (int i = 0; i < vm["five"].as<std::vector<std::string> >().size(); ++i)
        {
          five_adaptor = vm["five"].as<std::vector<std::string> >()[i];
          if (IUPAC)                    /* if IUPAC is true translate the IUPAC letter first */
            five_adaptor = translate_IUPAC(five_adaptor); /*   in   = GWGTTTGAAG */
          int five_allowed_mismatch_indel = get_allowed_mismatch_indel(vm, five_adaptor, "five");
          examine_5(seq,five_adaptor, five_allowed_mismatch_indel, "five", 
              case_insensitive, OS_alignment, d5, IUPAC);
          if (d5 == FOUND_5_ADAPT) // stop if already found the adaptor
            break;
        }
        std::ostream & OS_with_adaptor = *m[five_adaptor];     /* overwrites default OS_with_adaptor */
        Judge_Write(OS_with_adaptor, OS_no_adaptor, seq, d5, d3, format);
      }
    }
    else if (inspect_3)
    {
      if (vm.count("out_with_adaptor"))
      {
        out.open(vm["out_with_adaptor"].as<string>().c_str());
        if (!out.is_open())
          errx(1, "Error in open file %s" , __FILE__, __LINE__);
      }
      std::ostream & OS_with_adaptor = vm.count("out_with_adaptor") ? out : std::cout;

      Seqrecord seq;
//      CharString id, query, qual;
      while(ReadRecord(stream, seq, format))
      {
        ++total;
//        id = seq.ID; query = seq.Seq; qual = seq.Qual;  
        if (IUPAC)                              /* need to check 3' adaptor */
          three_adaptor = translate_IUPAC(three_adaptor); /* out = G[AT]GTTTGAAG */
        int three_allowed_mismatch_indel = get_allowed_mismatch_indel(vm, three_adaptor, "three");
        examine_3(seq,three_adaptor, three_allowed_mismatch_indel, "three", 
            case_insensitive, OS_alignment, d3, IUPAC);
        Judge_Write(OS_with_adaptor, OS_no_adaptor, seq, d5, d3, format);
      }
    }
    if (vm.count("out_with_adaptor"))           /* only need to delete when -o is supplied */
    {
      for(typeof(m.begin()) it = m.begin(); it != m.end(); ++it) /* clean up the map */
        delete it->second;
    }
  }

}

template <typename TStream>
string  get_Fastaq_format(TStream & stream)
{
  std::string format;
  if (stream.peek() == '>')
    format = "fasta";
  else if (stream.peek() == '@')
    format = "fastq";
  else
    errx(1, "unable to determine format");
  return format;
}

std::string get_format_app(const po::variables_map & vm)
{
  std::string format;
  if (vm["input"].as<vector<string> >().size() == 1 &&
      vm["input"].as<vector<string> >()[0] == string("stdin") )
  {
    if (!isatty(fileno(stdin)))
      format = get_Fastaq_format(std::cin);
  }
  else
  {
    string filename;
    filename = vm["input"].as<vector<string> >()[0]; 
    ifstream stream(filename.c_str());
    if (!stream.good())
      errx(1, "cannot open file %s for read", filename.c_str());
    format = get_Fastaq_format(stream);
    stream.close();
  }
  return format;
}

int main (int argc, char * argv[])
{
  bool case_insensitive;
  bool IUPAC;
  char buf[80];
  get_time_stamp(buf);
  string alignment_log = string("Alignment_log_") + string(buf);
  po::options_description desc(
      "\nPart of NGS C++ Toolkit by Zhigang Wu zhigang.wu@email.ucr.edu.\n\n"
      );
  desc.add_options()
    ("help,h", "Print this help page. True if present. [Boolean]")
//    ("five,5", po::value<string>(), 
    ("five,5", po::value<vector<string> >(), 
     "Five prime adaptor sequence. Any sequence follows this 5' adaptor will "
     "be retained and any sequence precedes 5' adaptor (including 5' adaptor) "
     "will be trimmed off. [str]")
    ("three,3", po::value<string>(), 
     "Three prime adaptor sequence. Any sequence precedes 3' adaptor will be retained "
     "and any sequence follows 3' adaptor (including 3' adaptor) will be "
     "trimmed off. [str]")
    ("input,i", po::value<vector<string> >()->default_value(vector<string>(1, "stdin"), 
                                                            "STDIN"), 
     "Input file, has to be in fastq or fasta format. [str]")
//    ("format,f", po::value<string>()->default_value("fasta", "fasta"), 
//     "the format of input file, valid formats include fastq and fasta [str]")
    ("five-mismatch,l", po::value<int>(), 
     "Allowed number of mismatches or gaps between the read sequence "
     "and the 5' adaptor sequence. By setting this value to 1, we are "
     "allowing two mismatches between adaptor and read sequence or 1 deletion "
     "in the read sequence or 2 insertions in the read sequence. So, if you set "
     "this value to 2, you multiply by 2 to get a sense of how many of mismatches &"
     "indels allowed in the alignment. Since the program gurantees only report the "
     "the best alignment, so setting this value the larger the better. Doing so is "
     "fine for those sequence with adaptor sequences. However, for those reads without "
     "adaptor squences will be wrongly trimmed. If you set this value, e.g. half of "
     "of the adaptor length, for those without adaptors an unpredictable position will "
     "be reported. In other words, this way will cause false positive. So, by "
     "default this value is set to 20% percent of the of 5' adaptor length. In "
     "circumstances where you want to do exact match, set this value to 0. [int]")
    ("three-mismatch,r", po::value<int>(), 
     "Generally same as -l option but for 3' adaptor sequence. "
     "By default this value is set to 20% percent of 3' adaptor "
     "sequence if you specified a 3' adaptor. [int]")
    ("percent,p", po::value<float>()->default_value(0.15, "0.15"), 
     "Percent of adaptor length used for maximum tolerable five-mismatch and "
     "three-mismatch. For example, if we set this number to 0.15, and the length of 5' "
     "adaptor is 20 and the length of 3' is 10, then 3 (0.15*20) and 2 ceil((0.15*10)) will "
     "be used as the value for the five-mismatch and three-mismatch, respectively. [float]")
    ("out_with_adaptor,o", po::value<string>(), 
     "Clean sequence with adaptor trimmed will be write to this file. default: STDOUT [str]")
    ("out_no_adaptor,n", po::value<string>(), 
     "Sequence without adaptor being found will be write to this file. default: STDOUT [str]")
    ("out_align,a", po::value<string>()->default_value(alignment_log, 
                                                       "Alignment_log_....txt"), 
     "If an adaptor is found within a read sequence, the alignment between adaptor(s) "
     "and this sequence and the distance between adaptor(s) and ends will be write "
     "to this file. The file is generated automatically you don't need to specify anything.")
    ("case-insensitive,I", po::value(&case_insensitive)->zero_tokens(), 
     "Toggle to switch case insensitive default is ON. True if present. [Boolean]")
    ("IUPAC,U", po::value(&IUPAC)->zero_tokens(), 
     "Toggle to start IUPAC match, default is OFF. Note if you set IUPAC mode, no mismatch "
     "indel will be allowed. This option is uncompatible with -l and -r option. "
     "True if present. [Boolean]")
    ("head,H", po::value<int>(), 
     "Cut the leading n base from the input seq. If this option is set, it will ignore"
     "-l -r -U, -I options which in this case make no sense. The trimmed seq will be print to"
     "-o, or STDOUT, those witout being trimmed due to  will be print to -n or STDOUT. [int]")
    ("tail,t", po::value<int>(), "Cut the tailing n base from the input seq. Others are same"
     "as -H option. [int]")
    ("length-cutoff", po::value<int>()->default_value(0, "0"), 
     "suppress trimmed sequence length less than this value from output. [int]")
    ;
  po::positional_options_description p;
//  p.add("five", -1);
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
  if(!vm.count("five") && !vm.count("three"))
  {
    if (!(vm.count("head") || vm.count("tail")))
    {
      cerr << "You have to at least supply one adaptor sequence via -5 and (or) -3 \n"
        << "Or at least supply a number of leading and tailing bases -H and (or) -t "
        << "to be trimmed"
        << endl;
      cerr << desc << endl;
      exit(1);
    }
  }
  format = get_format_app(vm);
  report_paramter_setting(vm);
/* 
 * start perform local alignment
 */
  
  adapt_5 d5;                                   /* d5 = UNDEFINED_5 by default */
  adapt_3 d3;                                   /* d3 = UNDEFINED_3 by default */

  ofstream out;
  ofstream out1;

  if (vm.count("out_with_adaptor"))
  {
    out.open(vm["out_with_adaptor"].as<string>().c_str());
    if (!out.is_open())
      errx(1, "Error in open file %s" , __FILE__, __LINE__);
  }
  std::ostream & OS_with_adaptor = vm.count("out_with_adaptor") ? out : std::cout;
  if (vm.count("out_no_adaptor") )
  {
    out1.open(vm["out_no_adaptor"].as<string>().c_str());
    if (!out1.is_open())
      errx(1, "Error in open file %s" , __FILE__, __LINE__);
  }
  std::ostream & OS_no_adaptor = vm.count("out_no_adaptor") ? out1 : std::cout;

  if(vm.count("head") || vm.count("tail"))/* if any of this two set, ignore all other options*/
  {
    if (vm["input"].as<vector<string> >().size() == 1 &&
        vm["input"].as<vector<string> >()[0] == string("stdin") )
    {
      if (!isatty(fileno(stdin)))
      {
        cut_leading_tailing_stream(std::cin, vm, OS_with_adaptor, OS_no_adaptor,d5, d3, 
                                  format.c_str());
        return 0;
      }
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
        cut_leading_tailing_stream(stream, vm, OS_with_adaptor, OS_no_adaptor,d5, d3, 
                                    format.c_str());
      }
      return 0;
    }

  }
// I put this below vm.count(header) and v.count(tail)
// because if any of them is true, we will not write alignment report
  ofstream out_alignment(vm["out_align"].as<string>().c_str()); 
  std::ostream & OS_alignment = out_alignment;  

  bool inspect_5 = vm.count("five") ? true : false;
  bool inspect_3 = vm.count("three") ? true : false;  
  if (inspect_5 && inspect_3)
  {
    // delete the pre-existing files in current directory

    RemovePreExistedFiles(vm, format.c_str());

    if (vm["input"].as<vector<string> >().size() == 1 &&
        vm["input"].as<vector<string> >()[0] == string("stdin") )
    {
      if (!isatty(fileno(stdin)))
        TrimmingSeq_from_stream(std::cin, vm, case_insensitive, d5, d3, OS_alignment,
                               inspect_5, inspect_3,
                                  format.c_str());
      else
        errx(1, "The program is hanging and waiting for input from STDIN");
    }
    else                                        /* read from file */
    {

      string filename;
      for (int i = 0; i < vm["input"].as<vector<string> >().size(); ++i)
      {
        filename = vm["input"].as<vector<string> >()[i]; 
        ifstream stream(filename.c_str());
        if (!stream.good())
          errx(1, "cannot open file %s for read", filename.c_str());
        TrimmingSeq_from_stream(stream, vm, case_insensitive, d5, d3, OS_alignment,
                              inspect_5, inspect_3,
                                  format.c_str());
      }
    }
  }
  else if (inspect_5)
  {
    bool inspect_5 = true;
    bool inspect_3 = false;    /* don't need to look at 3' end because 3' adaptor specified */
    // delete the pre-existing files in current directory
    RemovePreExistedFiles(vm, vm["format"].as<string>().c_str());

   if (vm["input"].as<vector<string> >().size() == 1 &&
      vm["input"].as<vector<string> >()[0] == string("stdin") )
    {
      if (!isatty(fileno(stdin)))
        TrimmingSeq_from_stream(std::cin, vm, case_insensitive, d5, d3, OS_alignment,
                               inspect_5, inspect_3,
                                  format.c_str());
      else
        errx(1, "The program is hanging and waiting for input from STDIN");
    }
    else                                        /* read from file */
    {
      string filename;
      for (int i = 0; i < vm["input"].as<vector<string> >().size(); ++i)
      {
        filename = vm["input"].as<vector<string> >()[i]; 
        ifstream stream(filename.c_str());
        if (!stream.good())
          errx(1, "cannot open file %s for read", filename.c_str());
        TrimmingSeq_from_stream(stream,vm, case_insensitive, d5, d3, OS_alignment,
                              inspect_5, inspect_3,
                                  format.c_str());
      }
     
    }

  }
  else if (inspect_3)
  {
    bool inspect_5 = false;   /* don't need to look at 5' end because 5' adaptor specified  */
    bool inspect_3 = true;                   
    if (vm["input"].as<vector<string> >().size() == 1 &&
      vm["input"].as<vector<string> >()[0] == string("stdin") )
    {
      if (!isatty(fileno(stdin)))
        TrimmingSeq_from_stream(std::cin, vm, case_insensitive, d5, d3, OS_alignment,
                               inspect_5, inspect_3,
                                  format.c_str());
      else
        errx(1, "The program is hanging and waiting for input from STDIN");
    }
    else                                        /* read from file */
    {
      string filename;
      for (int i = 0; i < vm["input"].as<vector<string> >().size(); ++i)
      {
        filename = vm["input"].as<vector<string> >()[i]; 
        ifstream stream(filename.c_str());
        if (!stream.good())
          errx(1, "cannot open file %s for read", filename.c_str());
        TrimmingSeq_from_stream(stream,vm, case_insensitive, d5, d3, OS_alignment,
                              inspect_5, inspect_3,
                                  format.c_str());
      }
    }
  }
  cerr << setw(20) << left << "Number of Seq only has 5' adaptor: " 
    << setw(10) << right << have_5_adaptor_only_counter << endl;
  cerr << setw(20) << left << "Number of Seq only has 3' adaptor: " 
    << setw(10) << right << have_3_adaptor_only_counter << endl;
  cerr << setw(20) << left << "Number of Seq has both adaptors: " 
    << setw(10) << right << have_5_and_3_adaptor_counter << endl;
  cerr << setw(20) << left << "Number of Seq without adaptor found: " 
    << setw(10) << right << no_5_and_3_adaptor_counter << endl;
  cerr << setw(20) << left << "Total Seq processed: "
    << setw(10) << right << total << endl;
  cerr << setw(20) << left << "Total Seq sum: " 
    << setw(10) << right << have_5_and_3_adaptor_counter + \
    have_5_adaptor_only_counter + have_3_adaptor_only_counter + \
    no_5_and_3_adaptor_counter
    << endl;
/* total should equal to sum
 */
  return 0;
}

