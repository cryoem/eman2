#if 0
#ifndef __Tokenizer__
#define  __Tokenizer__

#include <string>
using std::string;

class Tokenizer
{  public:
     Tokenizer(const string& str, 
               const string& delim = WHITE,
               bool want_delim = false)
     : target(str), position(0), delimiter(delim), 
       delimToken(want_delim) {  }
     bool moreToken() const;
     string nextToken();
     string setDelimiter(const string& delim)
     {  delimiter = delim; }
     int tokenCount() const;

   private:
     string next();
     string nextAll();
     static const string WHITE;
     int position;
     const string& target;
     string delimiter;
     bool delimToken;
};
#endif //  __Tokenizer__
#endif
