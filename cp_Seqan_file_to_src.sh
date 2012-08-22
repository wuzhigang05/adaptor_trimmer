rm -rf SeqAn1.3
mkdir SeqAn1.3
mkdir SeqAn1.3/seqan
for i in `get_what_Seqan_file_used.sh seqan`; do cp -r ~/software/Cplusplus_libs/seqan-trunk/core/include/$i* SeqAn1.3/seqan/; done
cp ~/software/Cplusplus_libs/seqan-trunk/LICENSE SeqAn1.3
cp ~/software/Cplusplus_libs/seqan-trunk/README SeqAn1.3
cp -r /Users/kaitang/software/Cplusplus_libs/seqan-trunk/core/include/seqan/platform* SeqAn1.3/seqan/
cp -r /Users/kaitang/software/Cplusplus_libs/seqan-trunk/core/include/seqan/misc* SeqAn1.3/seqan/
