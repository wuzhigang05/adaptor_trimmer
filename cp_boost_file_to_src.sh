rm -rf Boost1.50
mkdir Boost1.50
bcp boost/regex.hpp boost/program_options Boost1.50/ --boost=/Users/kaitang/software/Cplusplus_libs/boost/boost_1_50_0/
# bcp omits to copy the program_options.hpp, so below I copy it manually
cp /Users/kaitang/software/Cplusplus_libs/boost/boost_1_50_0/boost/program_options.hpp Boost1.50/boost/
# copy the license
cp ~/software/Cplusplus_libs/boost/boost_1_50_0/LICENSE_1_0.txt Boost1.50
