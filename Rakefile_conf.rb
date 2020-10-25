# check for a conf file on a upper directory
if File.exists?(File.expand_path('../Rakefile_conf.rb', File.dirname(__FILE__))) then
  require_relative '../Rakefile_conf.rb'
else
  COMPILE_DEBUG      = false
  COMPILE_DYNAMIC    = false
  COMPILE_EXECUTABLE = true

  case RUBY_PLATFORM
  when /mingw|mswin/
    PARALLEL = ''
    QUIET    = ''
  else
    cmakeversion = %x( cmake --version ).scan(/\d+\.\d+\.\d+/).last
    mm = cmakeversion.split('.');
    if mm[0].to_i > 3 || (mm[0].to_i == 3 && mm[1].to_i >= 12) then
      PARALLEL = '--parallel 8 '
      QUIET    = '-- --quiet '
    else
      PARALLEL = ''
      QUIET    = ''
    end
  end
end
