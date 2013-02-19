//
// C++ Interface: OutputMode
//
// Description: 
//
//
// Author: J. Kröker <kroeker@mathpc26>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//

namespace RationalMapSearch
{

/// @note: typesafe enum needs c++0x compiling option!

     /*
        enum class OutputMode
        {
            GAPOutput,
            M2Output,
            defaultOutput
        };
    */

        /// is an enum emulation; not completely error-prone but acceptable to provide backward compatibility instead of enums.
    class OutputMode
    {
  
    
        private :

            int enumVal_m;

            OutputMode(int val)
            {
                enumVal_m=val;
            }

         protected:
            
            static OutputMode getDefaultOutputMode()
            {
                return OutputMode(1);
            }

            static OutputMode getGAPOutputMode()
            {
                return OutputMode(2);
            }
            static OutputMode getM2OutputOutputMode()
            {
                return OutputMode(3);
            }
       public :
        static const OutputMode GAPOutput;
        static const OutputMode M2Output;
        static const OutputMode defaultOutput;
  
        bool operator==(const OutputMode & mode)  const     
            {
                return mode.enumVal_m==enumVal_m;
            }
            bool operator!=(const OutputMode & mode)  const   
            {
                return mode.enumVal_m!=enumVal_m;
            }
          
    };
}
