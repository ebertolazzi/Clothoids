#ifndef DOXYGEN_SHOULD_SKIP_THIS

//#include <stdio.h>
#include <string.h>
#include <net/if_dl.h>
#include <net/if.h>
#include <ifaddrs.h>
#include <netdb.h>
#include <cstdio>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <fstream>
#include <dirent.h>

#include <mach-o/dyld.h>

#include <string>
#include <vector>
#include <map>

namespace Utils {

  using std::string;
  using std::vector;
  using std::map;

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  string
  get_host_name() {
    char szHostName[1024];
    bool ok = gethostname( szHostName, 1024 ) == 0;
    if ( ok ) return string{szHostName};
    return string{""};
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  void
  get_MAC_address( map<string,string> & addr ) {
    struct ifaddrs *ifap, *ifaptr;

    int ier = getifaddrs(&ifap);
    UTILS_ASSERT(
      ier >= 0,
      "No network interface found for licence managing, getifaddrs return = {}\n", ier
    );

    addr.clear();
    for ( ifaptr = ifap;
          ifaptr != nullptr;
          ifaptr = ifaptr->ifa_next ) {
      if ( strncmp( "en", ifaptr->ifa_name, 2 ) == 0 &&
           ifaptr->ifa_addr->sa_family == AF_LINK ) {
        char saddr[100];
        struct sockaddr_dl * tmp = reinterpret_cast<struct sockaddr_dl *>((ifaptr)->ifa_addr);
        unsigned char * ptr = reinterpret_cast<unsigned char *>(LLADDR(tmp));
        snprintf(
          saddr, 100, "%02x:%02x:%02x:%02x:%02x:%02x",
          ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5]
        );
        addr[ifaptr->ifa_name] = saddr;
      }
    }
    freeifaddrs(ifap);
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  // Fetches the IP address
  void
  get_IP_address( vector<string> & addr ) {
    addr.clear ();
    char szHostName[128];
    if( gethostname(szHostName, 128) == 0 ) {
      // Get host adresses
      struct hostent * pHost;
      pHost = gethostbyname(szHostName);

      for( int i = 0; pHost!= nullptr && pHost->h_addr_list[i]!= nullptr; i++ ) {
        unsigned char * h = reinterpret_cast<unsigned char*>(pHost->h_addr_list[i]);
        string str{""};
        for( int j = 0; j < pHost->h_length; j++ ) {
          if( j > 0 ) str += '.';
          str += fmt::format("{}",h[j]);
        }
        addr.push_back( str );
        // str now contains one local IP address - do whatever you want to do with it (probably add it to a list)
      }
    }
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  string
  get_user_name() {
    char const * USER = getenv("USER");
    UTILS_ASSERT( USER != nullptr, "get_user_name, undefined enviroment `USER`" );
    return string{USER};
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  string
  get_home_directory() {
    char const * HOME = getenv("HOME");
    UTILS_ASSERT( HOME != nullptr, "get_home_directory, undefined enviroment `HOME`" );
    return string{HOME};
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  std::string
  get_executable_path_name() {
    uint32_t pathNameSize = 0;
    _NSGetExecutablePath( nullptr, &pathNameSize );
    std::vector<char> res(pathNameSize+1);
    std::fill(res.begin(),res.end(),'\0');
    if ( _NSGetExecutablePath( res.data(), &pathNameSize) != 0 ) return "NOT FOUND";
    return std::string{res.data()};
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  bool
  check_if_file_exists( char const * fname ) {
    struct stat buffer;
    return (stat (fname, &buffer) == 0);
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  bool
  check_if_dir_exists( char const * dirname ) {
    DIR *pDir = opendir(dirname);
    bool bExists = false;
    if ( pDir != nullptr ) {
      bExists = true;
      (void) closedir (pDir);
    }
    return bExists;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  bool
  make_directory( char const * dirname, unsigned mode ) {
    bool ok = check_if_dir_exists( dirname );
    if ( !ok ) ok = mkdir( dirname, mode_t(mode) ) == 0;
    return ok;
  }

}

#endif
