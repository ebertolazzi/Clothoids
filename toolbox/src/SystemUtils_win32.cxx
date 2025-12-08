/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#pragma once

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <direct.h>
#include <fstream>

#ifdef UTILS_OS_MINGW
#include <sys/stat.h>
#endif

namespace Utils
{

  using std::map;
  using std::string;
  using std::vector;

  static char const digits[] = "0123456789ABCDEF";

  // Will be constructed exactly once ...
  static class WinSockStub
  {
  public:
    WinSockStub( void )
    {
      WORD    wVersionRequested;
      WSADATA wsaData;
      // Using MAKEWORD macro, Winsock version request 2.2
      wVersionRequested = MAKEWORD( 2, 2 );
      /* int wsaerr =*/WSAStartup( wVersionRequested, &wsaData );
      return;
    }
  } __winsock_stub;

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  bool
  get_environment( string_view ename, string & res )
  {
    char  buffer[1024];
    DWORD var_size = GetEnvironmentVariable( ename.data(), buffer, 1024 );
    res            = string{ buffer };
    return var_size != 0;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  void
  set_environment( string_view ename, string_view newval, bool /* overwrite */ )
  {
    SetEnvironmentVariable( ename.data(), newval.data() );
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  // Fetches the MAC address
  void
  get_MAC_address( map<string, string> & addr )
  {
    IP_ADAPTER_INFO AdapterInfo[16];                   // Allocate information for up to 16 NICs
    DWORD           dwBufLen = sizeof( AdapterInfo );  // Save the memory size of buffer

    DWORD dwStatus = GetAdaptersInfo( AdapterInfo,  // [out] buffer to receive data
                                      &dwBufLen );  // [in] size of receive data buffer

    assert( dwStatus == ERROR_SUCCESS );  // Verify return value is valid, no buffer overflow

    PIP_ADAPTER_INFO pAdapterInfo = AdapterInfo;  // Contains pointer to current adapter info
    addr.clear();

    do
    {
      if ( pAdapterInfo->Type == MIB_IF_TYPE_ETHERNET || pAdapterInfo->Type == IF_TYPE_IEEE80211 )
      {
        unsigned char const * MACData{ pAdapterInfo->Address };
        string                str;
        for ( int k{ 0 }; k < 6; ++k )
        {
          str += digits[MACData[k] / 16];
          str += digits[MACData[k] % 16];
          if ( k < 5 ) str += ':';
        }
        // pAdapterInfo->Description
        addr[pAdapterInfo->AdapterName] = str;
      }
      pAdapterInfo = pAdapterInfo->Next;  // Progress through linked list
    } while ( pAdapterInfo );  // Terminate if last adapter
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  string
  get_host_name()
  {
    char szHostName[1024];
    bool ok = gethostname( szHostName, 1024 ) == 0;
    if ( ok ) return string( szHostName );
    return string( "" );
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  // Fetches the MAC address
  void
  get_IP_address( vector<string> & addr )
  {
    IP_ADAPTER_INFO AdapterInfo[16];
    DWORD           dwBufLen = sizeof( AdapterInfo );
    DWORD           dwStatus = GetAdaptersInfo( AdapterInfo, &dwBufLen );
    assert( dwStatus == ERROR_SUCCESS );  // Verify return value is valid, no buffer overflow

    PIP_ADAPTER_INFO pAdapterInfo = AdapterInfo;  // Contains pointer to current adapter info
    addr.clear();

    do
    {
      addr.emplace_back( pAdapterInfo->IpAddressList.IpAddress.String );
      pAdapterInfo = pAdapterInfo->Next;  // Progress through linked list
    } while ( pAdapterInfo );
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  string
  get_date()
  {
    SYSTEMTIME st;
    GetSystemTime( &st );
    return fmt::format( "{:04}-{:02}-{:02}", st.wYear, st.wMonth, st.wDay );
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  string
  get_day_time()
  {
    SYSTEMTIME st;
    GetSystemTime( &st );
    return fmt::format( "{:02}-{:02}-{:02}", st.wHour, st.wMinute, st.wSecond );
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  string
  get_log_date_time()
  {
    SYSTEMTIME st;
    GetSystemTime( &st );
    return fmt::format( "date_{:04}-{:02}-{:02}_time_{:02}-{:02}-{:02}", st.wYear, st.wMonth, st.wDay, st.wHour,
                        st.wMinute, st.wSecond );
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  string
  get_day_time_and_date()
  {
    SYSTEMTIME st;
    GetSystemTime( &st );
    return fmt::format( "{:02}-{:02}-{:02} {:04}-{:02}-{:02}", st.wHour, st.wMinute, st.wSecond, st.wYear, st.wMonth,
                        st.wDay );
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  string
  get_user_name()
  {
#ifdef UTILS_OS_MINGW
    char const * USER = getenv( "USER" );
    UTILS_ASSERT( USER != nullptr, "get_user_name, undefined enviroment `USER`" );
    return string{ USER };
#else
    char buffer[1024];
    GetEnvironmentVariable( "USER", buffer, 1024 );
    if ( buffer[0] == '\0' ) GetEnvironmentVariable( "USERNAME", buffer, 1024 );
    return string( buffer );
#endif
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  string
  get_home_directory()
  {
#ifdef UTILS_OS_MINGW
    char const * HOME = getenv( "HOME" );
    UTILS_ASSERT( HOME != nullptr, "get_home_directory, undefined enviroment `HOME`" );
    return string{ HOME };
#else
    char buffer[1024];
    GetEnvironmentVariable( "HOMEDRIVE", buffer, DWORD( 1024 ) );
    string res = buffer;
    GetEnvironmentVariable( "HOMEPATH", buffer, DWORD( 1024 ) );
    res += buffer;
    return res;
#endif
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  std::string
  get_executable_path_name()
  {
    char  pathName[1024];
    DWORD pathNameCapacity = 1024;
    GetModuleFileNameA( nullptr, pathName, pathNameCapacity );
    return std::string( pathName );
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */

  bool
  check_if_file_exists( string_view fname )
  {
#ifdef UTILS_OS_MINGW
    struct stat buffer;
    if ( stat( fnam.data(), &buffer ) == 0 ) return S_ISREG( buffer.st_mode );
    return false;
#else
    DWORD ftyp = GetFileAttributesA( fname.data() );
    if ( ftyp == INVALID_FILE_ATTRIBUTES ) return false;  // something is wrong with your path!
    if ( ftyp & FILE_ATTRIBUTE_DIRECTORY ) return false;  // this is a directory!
    return true;                                          // this is not a directory!
#endif
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  bool
  check_if_dir_exists( string_view dirname )
  {
#ifdef UTILS_OS_MINGW
    struct stat buffer;
    if ( stat( dirname.data(), &buffer ) == 0 ) return S_ISDIR( buffer.st_mode );
    return false;
#else
    DWORD ftyp = GetFileAttributesA( dirname.data() );
    if ( ftyp == INVALID_FILE_ATTRIBUTES ) return false;  // something is wrong with your path!
    if ( ftyp & FILE_ATTRIBUTE_DIRECTORY ) return true;   // this is a directory!
    return false;                                         // this is not a directory!
#endif
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  bool
  make_directory( string_view dirname, unsigned /* mode */ )
  {
    bool ok = check_if_dir_exists( dirname );
    if ( ok ) return false;
    CreateDirectoryA( dirname.data(), NULL );
    return true;
  }

}  // namespace Utils

#endif
