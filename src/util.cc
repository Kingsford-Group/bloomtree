#include "util.h"
#include <sstream>
#include <fenv.h>
#include <signal.h>


//
// Write out a number and the backspace over it
//
void
WriteStatusNumber(
  ostream & out,
  unsigned n)
{
  ostringstream s;
  s<<n;
  out << s.str();
  for(unsigned i = 0; i < s.str().length(); i++) out << "\b";
}

//
// Return the same string UPPERCASED
//
string
Upcase(const string & s)
{
  string out;
  for(unsigned i = 0; i < s.length(); i++) out += toupper(s[i]);
  return out;
}

string
SetAsString(const set<string> & s, const string & delim)
{
    string tmp = "";
    for(set<string>::iterator I = s.begin();
        I != s.end();
        ++I)
    {
        tmp += ((I == s.begin())?"":delim) + *I;
    }
    return tmp;
}

//
// Split string str into fields separated by sep
// field numbers start at 0
//
int
SplitString(
  const string & str1, 
  char sep,
  vector<string> & fields)
{
  string str = str1;
  fields.clear();
  string::size_type pos;
  while((pos=str.find(sep)) != string::npos)
  {
    fields.push_back(str.substr(0,pos));
    str.erase(0,pos+1);  
  }
  fields.push_back(str);
  return fields.size();
}



string
VectorAsString(
    const vector<string> & in,
    const string & delim
    )
{
    ostringstream oss;
    for(vector<string>::const_iterator I = in.begin();
        I != in.end();
        ++I)
    {
        if(I != in.begin()) oss << delim;
        oss << *I;
    }
    return oss.str();
}
