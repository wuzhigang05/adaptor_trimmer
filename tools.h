/* This file contains many functions that I wrote for convenience purpose
 * Copyright by Zhigang Wu
 * Department of Botany and Plant Sciences
 * University of California
 * Riverside, CA, 92521
 */

#ifndef __TOOL_H__
#define __TOOL_H__
#include <iostream>
#include <algorithm>
#include <err.h>
#include <string>
#include <vector>
#include <numeric>
template<typename T>
double mean (std::vector<T> & v)
{
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  double mean = sum / v.size();
  return mean;
//  double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
//  double stdev = std::sqrt(sq_sum / v.size() - mean * mean);
}

template<typename T>
double median (std::vector<T> & v)
{
  std::sort(v.begin(),v.end()); 
  if (v.size() % 2 == 0)
  {
    return double (v[v.size()/2 - 1] + v[v.size()/2])/2;
  }
  else
  {
    return  (double) v[v.size()/2];
  }
}


template<typename T>
std::ostream & print_v(std::ostream & out, const std::vector<T> & m, const char* separator = "\n")
{
  copy(m.begin(), m.end(), std::ostream_iterator<T>(out, separator));
  return out;
}

int require(std::string & arg, const char * value)
{
  if (arg != std::string(value))
  {
    errx (1, "ERROR: ", arg.c_str(), " != ", value);
  }
}
// below copy_if is like the filter function 
template<class InputIterator, class OutputIterator, class UnaryPredicate>
OutputIterator copy_if(InputIterator first, InputIterator last, 
                    OutputIterator d_first, UnaryPredicate pred)
{
    while (first != last) {
        if(pred(*first))
            *d_first++ = *first;
         first++;
    }
    return d_first;
}
/* 
 * readline from FILE, newline being removed
 */
int readline(FILE * IN, char * buf, int bufsize)
{
  if(fgets(buf, bufsize, IN) != NULL)
  {
    if (buf[strlen(buf)-1] == '\n')
      buf[strlen(buf)-1] = '\0';
    return 1;
  }
  return 0;
}

/* below is an example of using copy_if */
//std::vector<int> even_odd(std::vector<int>& v, 
//                      const char * even="yes")
//{
//  std::vector<int> tmp;
//  if (even == "yes")
//  {
//    copy_if(v.begin(), v.end(),
//        back_inserter(tmp), isEven);
//  }else{
//    copy_if(v.begin(), v.end(),
//        back_inserter(tmp), isOdd);
//  }
//  return tmp;
//}


/* foreach function is like the map function of */

/* transform function is also each more versatile than for_each function */

#endif __TOOL_H__
