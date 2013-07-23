#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>

class Exception
{
public:
    Exception(std::string _message = "Standart exception");
    std::string getMessage();
private:
    std::string message;
};


class NoLinesException: public Exception
{
public:
    NoLinesException(std::string _message): Exception(_message)
    {}

};

#endif // EXCEPTION_H
