#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <stdio.h>
#include <string.h>

#include <net/if.h>
#include <ifaddrs.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <unistd.h>

#include <sys/ioctl.h>
#include <sys/socket.h>
#include <sys/stat.h>

#include <fstream>
#include <dirent.h>

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
    if ( ok ) return string(szHostName);
    return string{""};
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  void
  get_MAC_address( map<string,string> & mac_addr ) {
    char   buf[8192] = {0};
    struct ifconf ifc = {0};
    struct ifreq *ifr = nullptr;
    char   ip[INET6_ADDRSTRLEN] = {0};

    // Get a socket handle.
    int sck = socket(PF_INET, SOCK_DGRAM, 0);
    UTILS_ASSERT0( sck >= 0, "get_MAC_address call of socket failed\n" );

    // Query available interfaces.
    ifc.ifc_len = sizeof(buf);
    ifc.ifc_buf = buf;
    UTILS_ASSERT0(
      ioctl(sck, SIOCGIFCONF, &ifc) >= 0,
      "get_MAC_address call of ioctl failed\n"
    );

    // Iterate through the list of interfaces.
    ifr = ifc.ifc_req;
    int nInterfaces = ifc.ifc_len / sizeof(struct ifreq);
    UTILS_ASSERT0( nInterfaces > 0, "get_MAC_address no interface found\n" );

    for( int i = 0; i < nInterfaces; ++i ) {

      struct ifreq * item = &ifr[i];
      if ( strncmp( "eth", item->ifr_name, 3 ) == 0 ||
           strncmp( "en",  item->ifr_name, 2 ) == 0 ) {
        struct sockaddr * addr = &(item->ifr_addr);
        char   macp[19];

        /* Get the IP address*/
        UTILS_ASSERT0(
          ioctl(sck, SIOCGIFADDR, item) >= 0,
          "get_MAC_address call of ioctl failed\n"
        );

        UTILS_ASSERT0(
          inet_ntop(AF_INET, &(((struct sockaddr_in *)addr)->sin_addr), ip, sizeof ip) != nullptr,
          "get_MAC_address call of inet_ntop failed\n"
        );

        // Get the MAC address
        UTILS_ASSERT0(
          ioctl(sck, SIOCGIFHWADDR, item) >= 0,
          "get_MAC_address call of ioctl failed\n"
        );

        snprintf(
          macp, 19, "%02x:%02x:%02x:%02x:%02x:%02x",
          static_cast<unsigned char>( item->ifr_hwaddr.sa_data[0] ),
          static_cast<unsigned char>( item->ifr_hwaddr.sa_data[1] ),
          static_cast<unsigned char>( item->ifr_hwaddr.sa_data[2] ),
          static_cast<unsigned char>( item->ifr_hwaddr.sa_data[3] ),
          static_cast<unsigned char>( item->ifr_hwaddr.sa_data[4] ),
          static_cast<unsigned char>( item->ifr_hwaddr.sa_data[5] )
        );

        mac_addr[item->ifr_name] = macp;
      }
    }
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
        unsigned char * h = reinterpret_cast<unsigned char*>( pHost->h_addr_list[i] );
        string str = "";
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
    char     pathName[1024];
    uint32_t pathNameCapacity = 1024;
    uint32_t pathNameSize     = uint32_t( readlink("/proc/self/exe", pathName, pathNameCapacity - 1) );
    pathName[pathNameSize]    = '\0';
    return std::string{pathName};
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
