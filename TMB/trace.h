#ifndef __TRACE__
#define __TRACE__

#undef HERE
/**
\def HERE
Indicates the line number of the file which the statement appears.
*/
#define HERE std::cout  << "reached " << __LINE__ << " in " << __FILE__ << "\n";

#undef HALT
/**
\def HALT
Prints the file name and line number and exits with exict code = 1.
*/
#define HALT std::cout <<"\nBailing out in file"<<__FILE__<<" at line " <<__LINE__<< std::endl; exit(1);

#undef TRACE
/**
\def TRACE
Prints the value of its argument, file name and line number.
*/
#define TRACE(object) std::cout << "line " << __LINE__ << ", file " << __FILE__ << ", " << #object " = " << object << "\n";

#undef TTRACE
/**
\def TTRACE
Prints the value of two arguments (note the double 'T'), file name and line number.
*/
#define TTRACE(o1,o2) std::cout << "line " << __LINE__ << ", file " << __FILE__ << ", " << #o1 " = " << o1<< ", " << #o2 " = " << o2 << "\n";

#undef ASSERT
/** 
\def ASSERT
It the argument is logically false, prints the file name, line number and value of argument and exits with exict code = 1.
*/
#define ASSERT(object) if (!object) { std::cout << "ASSERT: line = " << __LINE__ << " file = " << __FILE__ << " " << #object << " = " << object << " (false)\n"; exit(1); }
#endif
