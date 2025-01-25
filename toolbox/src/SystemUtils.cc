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

#if defined(__llvm__) || defined(__clang__)
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wduplicate-enum"
#endif

#include "Utils.hh"
#include "Utils_fmt.hh"

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

  /*!
   * \addtogroup OS
   * @{
   */

  //! Fetches the value of an environment variable.
  /*!
   * This function retrieves the value of the environment variable with the name
   * specified by `ename` and stores it in the reference `res`.
   *
   * \param ename Name of the environment variable to retrieve.
   * \param res Reference to a string where the result will be stored.
   * \return True if the environment variable was found and its value retrieved,
   *         false otherwise.
   */
   bool get_environment(char const ename[], string &res);

   //! Sets the value of an environment variable.
   /*!
    * This function sets the value of the environment variable specified by
    * `ename` to `newval`. If the variable already exists, it will be overwritten
    * if the `overwrite` flag is true.
    *
    * \param ename Name of the environment variable to set.
    * \param newval The new value to set for the environment variable.
    * \param overwrite Flag to indicate whether the environment variable should be
    *                  overwritten if it already exists.
    */
    void set_environment(char const ename[], char const newval[], bool overwrite);

    //! Retrieves the hostname of the system.
    /*!
     * \return A string containing the system's hostname.
     */
    string get_host_name();

    //! Retrieves the MAC addresses of network interfaces.
    /*!
     * This function retrieves the MAC addresses for all available network
     * interfaces on the system and stores them in the provided map, with
     * interface names as the keys and MAC addresses as the values.
     *
     * \param mac_addr A reference to a map where the MAC addresses will be stored.
     */
    void get_MAC_address(map<string,string> &mac_addr);

    //! Fetches the IP addresses of the system.
    /*!
     * This function retrieves all IP addresses associated with the current system's
     * network interfaces and stores them in the provided vector `addr`.
     *
     * \param addr A reference to a vector where the IP addresses will be stored.
     */
    void get_IP_address(vector<string> &addr);

    //! Retrieves the username of the current user.
    /*!
     * \return A string containing the username of the current user.
     */
    string get_user_name();

    //! Retrieves the home directory of the current user.
    /*!
     * \return A string containing the home directory path of the current user.
     */
    string get_home_directory();

    //! Retrieves the full path to the current executable.
    /*!
     * \return A string containing the full path to the executable.
     */
    std::string get_executable_path_name();

    //! Checks if a file exists.
    /*!
     * \param fname The path to the file to check.
     * \return True if the file exists and is a regular file, false otherwise.
     */
    bool check_if_file_exists(char const *fname);

    //! Checks if a directory exists.
    /*!
     * \param dirname The path to the directory to check.
     * \return True if the directory exists and is valid, false otherwise.
     */
    bool check_if_dir_exists(char const *dirname);

    //! Creates a directory if it does not exist.
    /*!
     * This function creates a directory with the specified mode if it does not
     * already exist.
     *
     * \param dirname The path to the directory to create.
     * \param mode The permissions mode to set for the new directory.
     * \return True if the directory was created or already exists, false otherwise.
     */
    bool make_directory(char const *dirname, unsigned mode);

  /*! @} */  // End of OS group

  #if defined(UTILS_OS_OSX) || defined(UTILS_OS_LINUX)

  string
  get_date() {
    char   buffer[20];
    time_t rawtime;
    time( &rawtime );
    struct tm timeinfo;
    localtime_r( &rawtime, &timeinfo );
    strftime( buffer, 20, "%F", &timeinfo );
    return string{buffer};
  }

  string
  get_day_time() {
    char   buffer[20];
    time_t rawtime;
    time( &rawtime );
    struct tm timeinfo;
    localtime_r( &rawtime, &timeinfo );
    strftime( buffer, 20, "%T", &timeinfo );
    return string{buffer};
  }

  string
  get_day_time_and_date() {
    char   buffer[20];
    time_t rawtime;
    time( &rawtime );
    struct tm timeinfo;
    localtime_r( &rawtime, &timeinfo );
    strftime( buffer, 20, "%T %F", &timeinfo );
    return string{buffer};
  }

  string
  get_log_date_time() {
    char   buffer[100];
    time_t rawtime;
    time( &rawtime );
    struct tm timeinfo;
    localtime_r( &rawtime, &timeinfo );
    strftime( buffer, 100, "date_%Y-%m-%d_time_%H-%M-%S", &timeinfo );
    return string{buffer};
  }

  #endif

}
