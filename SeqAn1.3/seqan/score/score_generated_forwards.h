// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================


#ifndef SEQAN_HEADER_SCORE_GENERATED_FORWARDS_H 
#define SEQAN_HEADER_SCORE_GENERATED_FORWARDS_H 

//////////////////////////////////////////////////////////////////////////////
// NOTE: This file is automatically generated by build_forwards.py
//       Do not edit this file manually!
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// CLASSES
//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN {

//____________________________________________________________________________
// Blosum30_

struct Blosum30_;       	// "core/include/seqan/score/score_matrix_data.h"(65)

//____________________________________________________________________________
// Blosum45_

struct Blosum45_;       	// "core/include/seqan/score/score_matrix_data.h"(133)

//____________________________________________________________________________
// Blosum62_

struct Blosum62_;       	// "core/include/seqan/score/score_matrix_data.h"(200)

//____________________________________________________________________________
// Blosum80_

struct Blosum80_;       	// "core/include/seqan/score/score_matrix_data.h"(267)

//____________________________________________________________________________
// Pam120_

struct Pam120_;       	// "core/include/seqan/score/score_matrix_data.h"(404)

//____________________________________________________________________________
// Pam200_

struct Pam200_;       	// "core/include/seqan/score/score_matrix_data.h"(473)

//____________________________________________________________________________
// Pam250_

struct Pam250_;       	// "core/include/seqan/score/score_matrix_data.h"(542)

//____________________________________________________________________________
// Pam40_

struct Pam40_;       	// "core/include/seqan/score/score_matrix_data.h"(334)

//____________________________________________________________________________
// Score

template <typename TValue , typename TSpec > class Score;       	// "core/include/seqan/score/score_base.h"(54)

//____________________________________________________________________________
// ScoreMatrix

template <typename TSequenceValue , typename TSpec > struct ScoreMatrix;       	// "core/include/seqan/score/score_matrix.h"(49)

//____________________________________________________________________________
// ScoringMatrixData_

template <typename TValue, typename TSequenceValue, typename TSpec> struct ScoringMatrixData_;       	// "core/include/seqan/score/score_matrix.h"(45)

//____________________________________________________________________________
// Vtml200_

struct Vtml200_;       	// "core/include/seqan/score/score_matrix_data.h"(611)

} //namespace SEQAN_NAMESPACE_MAIN
//////////////////////////////////////////////////////////////////////////////

namespace seqan {

//____________________________________________________________________________
// TagScoreMatrixFile_

struct TagScoreMatrixFile_;       	// "core/include/seqan/score/score_matrix_io.h"(65)

} //namespace seqan


//////////////////////////////////////////////////////////////////////////////
// TYPEDEFS
//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN {

//____________________________________________________________________________
// Blosum30

typedef Score<int, ScoreMatrix<AminoAcid, Blosum30_> > Blosum30;       	// "core/include/seqan/score/score_matrix_data.h"(77)

//____________________________________________________________________________
// Blosum45

typedef Score<int, ScoreMatrix<AminoAcid, Blosum45_> > Blosum45;       	// "core/include/seqan/score/score_matrix_data.h"(145)

//____________________________________________________________________________
// Blosum62

typedef Score<int, ScoreMatrix<AminoAcid, Blosum62_> > Blosum62;       	// "core/include/seqan/score/score_matrix_data.h"(212)

//____________________________________________________________________________
// Blosum80

typedef Score<int, ScoreMatrix<AminoAcid, Blosum80_> > Blosum80;       	// "core/include/seqan/score/score_matrix_data.h"(279)

//____________________________________________________________________________
// EditDistanceScore

typedef Score<int, EditDistance> EditDistanceScore;       	// "core/include/seqan/score/score_edit.h"(75)

//____________________________________________________________________________
// Pam120

typedef Score<int, ScoreMatrix<AminoAcid, Pam120_> > Pam120;       	// "core/include/seqan/score/score_matrix_data.h"(416)

//____________________________________________________________________________
// Pam200

typedef Score<int, ScoreMatrix<AminoAcid, Pam200_> > Pam200;       	// "core/include/seqan/score/score_matrix_data.h"(485)

//____________________________________________________________________________
// Pam250

typedef Score<int, ScoreMatrix<AminoAcid, Pam250_> > Pam250;       	// "core/include/seqan/score/score_matrix_data.h"(554)

//____________________________________________________________________________
// Pam40

typedef Score<int, ScoreMatrix<AminoAcid, Pam40_> > Pam40;       	// "core/include/seqan/score/score_matrix_data.h"(346)

//____________________________________________________________________________
// SimpleScore

typedef Score<int, Simple> SimpleScore;       	// "core/include/seqan/score/score_simple.h"(118)

//____________________________________________________________________________
// Vtml200

typedef Score<int, ScoreMatrix<AminoAcid, Vtml200_> > Vtml200;       	// "core/include/seqan/score/score_matrix_data.h"(623)

} //namespace SEQAN_NAMESPACE_MAIN
//////////////////////////////////////////////////////////////////////////////

namespace seqan {

//____________________________________________________________________________
// ScoreMatrixFile

typedef Tag<TagScoreMatrixFile_> ScoreMatrixFile;       	// "core/include/seqan/score/score_matrix_io.h"(66)

} //namespace seqan


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN {

//____________________________________________________________________________
// score

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2> inline TValue score(Score<TValue, TSpec> const & me, TPos1 pos1, TPos2 pos2, TSeq1 const &seq1, TSeq2 const &seq2);       	// "core/include/seqan/score/score_base.h"(206)
template <typename TValue, typename TSequenceValue, typename TSpec, typename TVal1, typename TVal2> inline TValue score(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc, TVal1 val1, TVal2 val2);       	// "core/include/seqan/score/score_matrix.h"(140)
template <typename TValue, typename TSpec, typename TVal1, typename TVal2> inline TValue score(Score<TValue, TSpec> const & me, TVal1 left, TVal2 right);       	// "core/include/seqan/score/score_simple.h"(353)

//____________________________________________________________________________
// scoreGap

template <typename TValue, typename TSpec> inline TValue scoreGap(Score<TValue, TSpec> const & me);       	// "core/include/seqan/score/score_simple.h"(309)

//____________________________________________________________________________
// scoreGapExtend

template <typename TValue> inline TValue scoreGapExtend(Score<TValue, EditDistance> &);       	// "core/include/seqan/score/score_edit.h"(111)
template <typename TValue> inline TValue scoreGapExtend(Score<TValue, EditDistance> const &);       	// "core/include/seqan/score/score_edit.h"(119)
template <typename TValue, typename TSpec> inline TValue scoreGapExtend(Score<TValue, TSpec> const & me);       	// "core/include/seqan/score/score_simple.h"(220)

//____________________________________________________________________________
// scoreGapExtendHorizontal

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2> inline TValue scoreGapExtendHorizontal( Score<TValue, TSpec> const & me, TPos1, TPos2, TSeq1 const &, TSeq2 const &);       	// "core/include/seqan/score/score_base.h"(120)

//____________________________________________________________________________
// scoreGapExtendVertical

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2> inline TValue scoreGapExtendVertical( Score<TValue, TSpec> const & me, TPos1, TPos2, TSeq1 const &, TSeq2 const &);       	// "core/include/seqan/score/score_base.h"(140)

//____________________________________________________________________________
// scoreGapHorizontal

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2> inline TValue scoreGapHorizontal( Score<TValue, TSpec> const & me, TPos1, TPos2, TSeq1 const &, TSeq2 const &);       	// "core/include/seqan/score/score_base.h"(160)

//____________________________________________________________________________
// scoreGapOpen

template <typename TValue> inline TValue scoreGapOpen(Score<TValue, EditDistance> &);       	// "core/include/seqan/score/score_edit.h"(127)
template <typename TValue> inline TValue scoreGapOpen(Score<TValue, EditDistance> const &);       	// "core/include/seqan/score/score_edit.h"(135)
template <typename TValue, typename TSpec> inline TValue scoreGapOpen(Score<TValue, TSpec> const & me);       	// "core/include/seqan/score/score_simple.h"(263)

//____________________________________________________________________________
// scoreGapOpenHorizontal

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2> inline TValue scoreGapOpenHorizontal( Score<TValue, TSpec> const & me, TPos1, TPos2, TSeq1 const &, TSeq2 const &);       	// "core/include/seqan/score/score_base.h"(79)

//____________________________________________________________________________
// scoreGapOpenVertical

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2> inline TValue scoreGapOpenVertical( Score<TValue, TSpec> const & me, TPos1, TPos2, TSeq1 const &, TSeq2 const &);       	// "core/include/seqan/score/score_base.h"(100)

//____________________________________________________________________________
// scoreGapVertical

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2> inline TValue scoreGapVertical( Score<TValue, TSpec> const & me, TPos1, TPos2, TSeq1 const &, TSeq2 const &);       	// "core/include/seqan/score/score_base.h"(180)

//____________________________________________________________________________
// scoreMatch

template <typename TValue> inline TValue scoreMatch(Score<TValue, EditDistance> &);       	// "core/include/seqan/score/score_edit.h"(79)
template <typename TValue> inline TValue scoreMatch(Score<TValue, EditDistance> const &);       	// "core/include/seqan/score/score_edit.h"(87)
template <typename TValue, typename TSpec> inline TValue scoreMatch(Score<TValue, TSpec> const & me);       	// "core/include/seqan/score/score_simple.h"(136)

//____________________________________________________________________________
// scoreMismatch

template <typename TValue> inline TValue scoreMismatch(Score<TValue, EditDistance> &);       	// "core/include/seqan/score/score_edit.h"(95)
template <typename TValue> inline TValue scoreMismatch(Score<TValue, EditDistance> const &);       	// "core/include/seqan/score/score_edit.h"(103)
template <typename TValue, typename TSpec> inline TValue scoreMismatch(Score<TValue, TSpec> const & me);       	// "core/include/seqan/score/score_simple.h"(177)

//____________________________________________________________________________
// setDefaultScoreMatrix

template <typename TValue, typename TSequenceValue, typename TSpec, typename TTag> inline void setDefaultScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, TTag);       	// "core/include/seqan/score/score_matrix.h"(190)
template <typename TValue, typename TSequenceValue, typename TSpec> inline void setDefaultScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, Default);       	// "core/include/seqan/score/score_matrix.h"(207)

//____________________________________________________________________________
// setScore

template <typename TValue, typename TSequenceValue, typename TSpec, typename TVal1, typename TVal2, typename T> inline void setScore(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, TVal1 val1, TVal2 val2, T score);       	// "core/include/seqan/score/score_matrix.h"(165)

//____________________________________________________________________________
// setScoreGap

template <typename TValue, typename TSpec> inline void setScoreGap(Score<TValue, TSpec> & me, TValue const & value);       	// "core/include/seqan/score/score_simple.h"(330)

//____________________________________________________________________________
// setScoreGapExtend

template <typename TValue, typename TSpec> inline void setScoreGapExtend(Score<TValue, TSpec> & me, TValue const & value);       	// "core/include/seqan/score/score_simple.h"(241)

//____________________________________________________________________________
// setScoreGapOpen

template <typename TValue, typename TSpec> inline void setScoreGapOpen(Score<TValue, TSpec> & me, TValue const & value);       	// "core/include/seqan/score/score_simple.h"(284)

//____________________________________________________________________________
// setScoreMatch

template <typename TValue, typename TSpec> inline void setScoreMatch(Score<TValue, TSpec> & me, TValue const & value);       	// "core/include/seqan/score/score_simple.h"(156)

//____________________________________________________________________________
// setScoreMismatch

template <typename TValue, typename TSpec> inline void setScoreMismatch(Score<TValue, TSpec> & me, TValue const & value);       	// "core/include/seqan/score/score_simple.h"(198)

} //namespace SEQAN_NAMESPACE_MAIN
//////////////////////////////////////////////////////////////////////////////

namespace seqan {

//____________________________________________________________________________
// _sprintfValue

inline void _sprintfValue(char * buf, unsigned val);       	// "core/include/seqan/score/score_matrix_io.h"(346)
inline void _sprintfValue(char * buf, int val);       	// "core/include/seqan/score/score_matrix_io.h"(356)
inline void _sprintfValue(char * buf, float val);       	// "core/include/seqan/score/score_matrix_io.h"(366)
inline void _sprintfValue(char * buf, double val);       	// "core/include/seqan/score/score_matrix_io.h"(377)

//____________________________________________________________________________
// _sscanfValue

inline void _sscanfValue(const char * buf, unsigned int & val);       	// "core/include/seqan/score/score_matrix_io.h"(99)
inline void _sscanfValue(const char * buf, int & val);       	// "core/include/seqan/score/score_matrix_io.h"(116)
inline void _sscanfValue(const char * buf, float & val);       	// "core/include/seqan/score/score_matrix_io.h"(133)
inline void _sscanfValue(const char * buf, double & val);       	// "core/include/seqan/score/score_matrix_io.h"(150)

//____________________________________________________________________________
// loadScoreMatrix

template <typename TValue, typename TSequenceValue, typename TSpec> inline void loadScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, CharString const & filename);       	// "core/include/seqan/score/score_matrix_io.h"(296)
template <typename TValue, typename TSequenceValue, typename TSpec, typename TMeta> inline void loadScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, CharString const & filename, TMeta & meta);       	// "core/include/seqan/score/score_matrix_io.h"(314)

//____________________________________________________________________________
// read

template <typename TStream, typename TValue, typename TSequenceValue, typename TSpec> int read(RecordReader<TStream, SinglePass<> > & reader, Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, ScoreMatrixFile const & );       	// "core/include/seqan/score/score_matrix_io.h"(188)
template <typename TStream, typename TValue, typename TSequenceValue, typename TSpec> inline void read(RecordReader<TStream, SinglePass<> > & fl, Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc);       	// "core/include/seqan/score/score_matrix_io.h"(274)

//____________________________________________________________________________
// readMeta

template <typename TStream, typename TMeta> int readMeta(RecordReader<TStream, SinglePass<> > & reader, TMeta & meta, ScoreMatrixFile const & );       	// "core/include/seqan/score/score_matrix_io.h"(166)

//____________________________________________________________________________
// write

template <typename TStream, typename TValue, typename TSequenceValue, typename TSpec> inline int write(TStream & stream, Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc, CharString const & meta);       	// "core/include/seqan/score/score_matrix_io.h"(404)
template <typename TStream, typename TValue, typename TSequenceValue, typename TSpec> inline int write(TStream & stream, Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc);       	// "core/include/seqan/score/score_matrix_io.h"(494)

} //namespace seqan

#endif

