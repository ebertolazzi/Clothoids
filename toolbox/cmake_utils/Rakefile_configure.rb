#
# Check for a configuration file on a upper directory.
# This permits to use a unique configuration file for
# large projects.
# On a local project use the default in this file.
#
if File.exist?(File.expand_path('../Rakefile_configure.rb', File.dirname(__FILE__))) then
  # found in the root of the local project
  require_relative '../Rakefile_configure.rb'
elsif File.exist?(File.expand_path('../../Rakefile_configure.rb', File.dirname(__FILE__))) then
  # found in the upper project
  require_relative '../../Rakefile_configure.rb'
elsif File.exist?(File.expand_path('../../cmake_utils/Rakefile_configure.rb', File.dirname(__FILE__))) then
  # found in the upper project under cmake_utils
  require_relative '../../cmake_utils/Rakefile_configure.rb'
else
  #-------------------------
  COMPILE_DEBUG      = false
  COMPILE_DYNAMIC    = false
  COMPILE_EXECUTABLE = true
  #-------------------------
end
