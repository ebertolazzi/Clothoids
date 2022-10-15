/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2013                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Utils.hh"

#if defined(UTILS_OS_WINDOWS)
  // MINGW is contained in WINDOWS 
  #include "SystemUtils_win32.cxx"
#elif defined(UTILS_OS_OSX)
  #include "SystemUtils_osx.cxx"
#elif defined(UTILS_OS_LINUX)
  #include "SystemUtils_linux.cxx"
#else
  #error "unsupported OS!"
#endif

namespace Utils {

  #if defined(UTILS_OS_OSX) || defined(UTILS_OS_LINUX)

  std::string
  get_date() {
    char   buffer[20];
    time_t rawtime;
    time( &rawtime );
    struct tm timeinfo;
    localtime_r( &rawtime, &timeinfo );
    strftime( buffer, 20, "%F", &timeinfo );
    return std::string{buffer};
  }

  std::string
  get_day_time() {
    char   buffer[20];
    time_t rawtime;
    time( &rawtime );
    struct tm timeinfo;
    localtime_r( &rawtime, &timeinfo );
    strftime( buffer, 20, "%T", &timeinfo );
    return std::string{buffer};
  }

  std::string
  get_day_time_and_date() {
    char   buffer[20];
    time_t rawtime;
    time( &rawtime );
    struct tm timeinfo;
    localtime_r( &rawtime, &timeinfo );
    strftime( buffer, 20, "%T %F", &timeinfo );
    return std::string{buffer};
  }

  std::string
  get_log_date_time() {
    char   buffer[100];
    time_t rawtime;
    time( &rawtime );
    struct tm timeinfo;
    localtime_r( &rawtime, &timeinfo );
    strftime( buffer, 100, "date_%Y-%m-%d_time_%H-%M-%S", &timeinfo );
    return std::string{buffer};
  }

  #endif

}
