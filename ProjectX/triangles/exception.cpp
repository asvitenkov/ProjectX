
#include "exception.h"

Exception::Exception(std::string _message): message(_message)
{
}

std::string Exception::getMessage()
{
    return message;
}
