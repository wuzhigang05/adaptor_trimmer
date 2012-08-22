pattern=$1
grep "#include <$pattern" -h  *.h *.cc | grep '//' -v | sort | uniq | sed 's_[^ ]*[ ]*<\([^ ]*\)>_\1_g' | tr '.' '\t' | cut -f1

