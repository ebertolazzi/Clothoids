# check for a conf file on a upper directory
if Dir.exist?('../Rakefile_conf.rb') then
  require_relative '../Rakefile_conf.rb'
else
  COMPILE_DEBUG      = false
  COMPILE_DYNAMIC    = false
  COMPILE_EXECUTABLE = false
end
