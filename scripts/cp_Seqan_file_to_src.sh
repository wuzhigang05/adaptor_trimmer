rm -rf SeqAn1.3
mkdir SeqAn1.3
mkdir SeqAn1.3/seqan

for i in `get_what_Seqan_file_used.sh seqan`; do cp -r ~/software/Cplusplus_libs/seqan-trunk/core/include/$i* SeqAn1.3/seqan/; done
python copy_seqan_dependence.py /Users/kaitang/software/Cplusplus_libs/seqan-trunk/extras/include/seqan SeqAn1.3/seqan/
# need to be further tested below script
# python copy_dependence.py seqan/sequence.h seqan/basic.h seqan/seq_io.h ~/software/Cplusplus_libs/seqan-trunk/core/include/ Test/


