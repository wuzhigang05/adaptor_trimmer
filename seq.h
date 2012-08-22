/*
 * =====================================================================================
 *
 *       Filename:  Seq.cc
 *
 *    Description:  This is a Sequence Class interface
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
#ifndef   __SEQ_H__
#define   __SEQ_H__
#include <string>

/*-----------------------------------------------------------------------------
 *  actually, I can just declare the Seq as a Structure
 *-----------------------------------------------------------------------------*/

//template <typename Tid = std::string, 
//         typename Tseq = std::string>
//class Seqrecord
//{
//public:
//  Tid ID;
//  Tseq Seq;
//  Seqrecord(Tid id, Tseq seq):ID(id), Seq(seq) {}
//  Seqrecord():ID(Tid()), Seq(Tseq()) {}
//
//};
class Seqrecord
{
public:
  std::string ID;
  std::string Seq;
  std::string Qual;
  Seqrecord(std::string id    , 
            std::string seq   = "",
            std::string qual  = ""):
                ID(id), Seq(seq), Qual(qual) {}
  Seqrecord():ID(""), Seq(""), Qual("") {}

};
#endif

