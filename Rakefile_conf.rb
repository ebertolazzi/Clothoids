# check for a conf file on a upper directory
if File.exists?('../Rakefile_conf.rb') then
  require_relative '../Rakefile_conf.rb'
else
  COMPILE_DEBUG      = false
  COMPILE_DYNAMIC    = false
  COMPILE_EXECUTABLE = false

  cmakeversion = %x( cmake --version ).scan(/\d+\.\d+/).last
  mm = cmakeversion.split('.');
  if mm[0].to_i > 3 || (mm[0].to_i == 3 && mm[1].to_i >= 12) then
    PARALLEL = '--parallel 8 '
    QUIET    = '-- --quiet '
  else
    PARALLEL = ''
    QUIET    = ''
  end
end
