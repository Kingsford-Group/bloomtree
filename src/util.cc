#include "util.h"
#include <sstream>
#include <fenv.h>
#include <signal.h>

std::string TrimRight(const std::string & str) {
    std::string tmp = str;
    return tmp.erase(tmp.find_last_not_of(" ") + 1);
}

std::string TrimLeft(const std::string & str) {
    std::string tmp = str;
    return tmp.erase(0, tmp.find_first_not_of(" "));
}

std::string Trim(const std::string & str) {
    std::string tmp = str;
    return TrimLeft(TrimRight(str));
}

//
// Die with the message
//
void DIE(const std::string & msg) {
    std::cerr << "fatal: " << msg << std::endl;
    exit(3);
}

void WARN(const std::string & msg) {
    std::cerr << "warning: " << msg << std::endl;
}

//
// Die, with a message, if bad is true
//
void DIE_IF(bool bad, const std::string & msg) {
    if(bad) DIE(msg);
}

//
// Write out a number and the backspace over it
//
void WriteStatusNumber(std::ostream & out, unsigned n) {
    std::ostringstream s;
    s<<n;
    out << s.str();
    for(unsigned i = 0; i < s.str().length(); i++) out << "\b";
}

//
// Return the same string UPPERCASED
//
std::string Upcase(const std::string & s) {
    std::string out;
    for(unsigned i = 0; i < s.length(); i++) out += toupper(s[i]);
    return out;
}

std::string SetAsString(const std::set<std::string> & s, const std::string & delim) {
    std::string tmp = "";
    for(std::set<std::string>::iterator I = s.begin();
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
int SplitString(const std::string & str1, char sep, std::vector<std::string> & fields) {
    std::string str = str1;
    fields.clear();
    std::string::size_type pos;
    while((pos=str.find(sep)) != std::string::npos) {
        fields.push_back(str.substr(0,pos));
        str.erase(0,pos+1);  
    }
    fields.push_back(str);
    return fields.size();
}


std::string VectorAsString(const std::vector<std::string> & in, const std::string & delim) {
    std::ostringstream oss;
    for(std::vector<std::string>::const_iterator I = in.begin();
        I != in.end();
        ++I)
    {
        if(I != in.begin()) oss << delim;
        oss << *I;
    }
    return oss.str();
}
