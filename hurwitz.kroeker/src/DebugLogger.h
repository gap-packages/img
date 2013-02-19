#pragma once


#include <iostream>
/// @note idea from http://bytes.com/topic/c/answers/127843-null-output-stream
struct nullstream:
    std::ostream {
        struct nullbuf: std::streambuf 
            {
                int overflow(int c) { return traits_type::not_eof(c); }
            } m_sbuf;
        nullstream() : std::ios(&m_sbuf), std::ostream(&m_sbuf)
     {}
};

// todo: man kann moeglicherweise #I bei allen ausgaben davorsetzen, wenn man den StreamOperator spezialisiert wie z.B. nullstream
class DebugLogger
{
public:
    static int level_g ;

    static nullstream  ns_g;

    static  void log(std::string message)
    {
        if (level_g >0)
            std::cerr << message;
    
    }

    static std::ostream & logStream()
    {
        if (level_g==0)
            return ns_g;
        return std::cerr;
    }

    static void setLevel(int level)
    {
        level_g=level;
    }

};

