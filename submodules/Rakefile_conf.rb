# check for a conf file on a upper directory
if File.exists?(File.expand_path('../Rakefile_conf.rb', File.dirname(__FILE__))) then
  require_relative '../Rakefile_conf.rb'
else
  case RUBY_PLATFORM
  when /mingw|mswin|aarch64-linux-gnu/ # dont paralellize on raspberry
    PARALLEL = ''
    QUIET    = ''
  else
    require 'etc'
    cmakeversion = %x( cmake --version ).scan(/\d+\.\d+\.\d+/).last
    mm = cmakeversion.split('.');
    if mm[0].to_i > 3 || (mm[0].to_i == 3 && mm[1].to_i >= 12) then
      PARALLEL = "--parallel #{Etc.nprocessors} "
      QUIET    = '-- --quiet '
    else
      PARALLEL = ''
      QUIET    = ''
    end
  end
end

#-------------------------
COMPILE_DEBUG      = false unless defined?(COMPILE_DEBUG)
COMPILE_DYNAMIC    = false unless defined?(COMPILE_DYNAMIC)
COMPILE_EXECUTABLE = false unless defined?(COMPILE_EXECUTABLE)
# compiler on WINDOWS
USE_NMAKE = true  unless defined?(USE_NMAKE)
USE_MINGW = false unless defined?(USE_MINGW)
USE_MSYS  = false unless defined?(USE_MSYS)
RUN_CPACK = false unless defined?(RUN_CPACK)
#-------------------------
