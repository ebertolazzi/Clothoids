#
# Check for a configuration file on a upper directory.
# This permits to use a unique configuration file for
# large projects.
# On a local project use the default in this file.
#
if File.exists?(File.expand_path('../Rakefile_configure.rb', File.dirname(__FILE__))) then
  # found in the root of the local project
  require_relative '../Rakefile_configure.rb'
elsif File.exists?(File.expand_path('../../Rakefile_configure.rb', File.dirname(__FILE__))) then
  # found in the upper project
  require_relative '../../Rakefile_configure.rb'
elsif File.exists?(File.expand_path('../../cmake_utils/Rakefile_configure.rb', File.dirname(__FILE__))) then
  # found in the upper project under cmake_utils
  require_relative '../../cmake_utils/Rakefile_configure.rb'
else
  #-------------------------
  COMPILE_DEBUG      = false
  COMPILE_DYNAMIC    = false
  COMPILE_EXECUTABLE = true
  # compiler on WINDOWS
  # possible values are: :vs, :nmake, :mingw, :msys
  WIN_COMPILER = :nmake # The default is to use Visual Studio with nmake
  #-------------------------
end
