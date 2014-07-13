#ifndef UTIL_H
#define UTIL_H
#include <utility>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <cassert>

using namespace std;


//
// Print the elements of a map nicely
//
template <class T, class Q>
void
PrintMap(
  ostream & ostr,
  const map<T,Q> & m,
  const string & key_value_sep,
  const string & pair_sep)
{
  for(typename map<T,Q>::const_iterator M = m.begin();
      M != m.end();
      ++M)
  {
    ostr << M->first << key_value_sep << M->second << pair_sep;
  }
}



inline 
string 
TrimRight(const string & str)
{
    string tmp = str;
    return tmp.erase(tmp.find_last_not_of(" ") + 1);
}

inline
string
TrimLeft(const string & str)
{
    string tmp = str;
    return tmp.erase(0, tmp.find_first_not_of(" "));
}

inline
string
Trim(const string & str)
{
    string tmp = str;
    return TrimLeft(TrimRight(str));
}

//==========================================================
// Error messages
//==========================================================

//
// Die with the message
//
inline
void
DIE(const string & msg)
{
  cerr << "fatal: " << msg << endl;
  exit(3);
}

inline
void
WARN(const string & msg)
{
    cerr << "warning: " << msg << endl;
}

//
// Die, with a message, if bad is true
//
inline
void
DIE_IF(bool bad, const string & msg)
{
  if(bad) DIE(msg);
}

//
// Die, with a message, if bad is true
//
inline
void
ERROR_IF(bool bad, int line_number, const string & msg)
{
  if(bad) 
  {
    cerr << line_number << ": "; 
    DIE(msg);
  }
}

//
// Write a number then backspace
//
void
WriteStatusNumber(
  ostream & out,
  unsigned n);

//
// Return the same string UPPERCASED
//
string
Upcase(const string &);

//
// Split string str into fields separated by sep
// field numbers start at 0
//
int
SplitString(
  const string &, 
  char,
  vector<string> &);

string
SetAsString(const set<string> &,  const string & = " ");

string
VectorAsString( const vector<string> &, const string & = " ");

#endif
