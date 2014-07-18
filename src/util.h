#ifndef UTIL_H
#define UTIL_H
#include <utility>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <cassert>


//
// Print the elements of a map nicely
//
template <class T, class Q>
void
PrintMap(
  std::ostream & ostr,
  const std::map<T,Q> & m,
  const std::string & key_value_sep,
  const std::string & pair_sep)
{
  for(typename std::map<T,Q>::const_iterator M = m.begin();
      M != m.end();
      ++M)
  {
    ostr << M->first << key_value_sep << M->second << pair_sep;
  }
}


std::string TrimRight(const std::string & str);
std::string TrimLeft(const std::string & str);
std::string Trim(const std::string & str);


//==========================================================
// Error messages
//==========================================================
void DIE(const std::string & msg);
void WARN(const std::string & msg);
void DIE_IF(bool bad, const std::string & msg);


//
// Write a number then backspace
//
void WriteStatusNumber(std::ostream & out, unsigned n);

//
// Return the same string UPPERCASED
//
std::string Upcase(const std::string &);

//
// Split string str into fields separated by sep
// field numbers start at 0
//
int SplitString(const std::string &, char, std::vector<std::string> &);

std::string SetAsString(const std::set<std::string> &,  const std::string & = " ");

std::string VectorAsString(const std::vector<std::string> &, const std::string & = " ");

#endif
