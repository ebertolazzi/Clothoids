#
#
#

%w(colorize rake fileutils).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

require "rake/clean"

CLEAN.include   ["./**/*.o", "./**/*.obj", "./bin/**/example*", "./build"]
CLOBBER.include []
CLEAN.exclude('**/[cC][oO][rR][eE]')

verbose(false)

cmakeversion = %x( cmake --version ).scan(/\d+\.\d+/).last
if cmakeversion >= "3.12" then
  PARALLEL = '--parallel 8 '
else
  PARALLEL = ''
end

task :default => [:build]

LIB_NAME="GenericContainer"

desc "run tests"
task :run do
	sh "./bin/example1"
	sh "./bin/example2"
	sh "./bin/example3"
	sh "./bin/example4"
	sh "./bin/example5"
	sh "./bin/example6"
	sh "./bin/example7"
	sh "./bin/example8"
	sh "./bin/example9"
	sh "./bin/example10"
	sh "./bin/example11"
end

desc "run tests"
task :run_win do
	sh "./bin/Release/example1"
	sh "./bin/Release/example2"
	sh "./bin/Release/example3"
	sh "./bin/Release/example4"
	sh "./bin/Release/example5"
	sh "./bin/Release/example6"
	sh "./bin/Release/example7"
	sh "./bin/Release/example8"
	sh "./bin/Release/example9"
	sh "./bin/Release/example10"
	sh "./bin/Release/example11"
end

desc "compile for UNIX/OSX"
task :build, [:executable] do |t, args|
  args.with_defaults(:executable => "no" )

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  puts "\n\nPrepare #{LIB_NAME} project".green

  if args.executable.to_s == "yes" then
    sh 'cmake -DBUILD_EXECUTABLE=true -DCMAKE_INSTALL_PREFIX:PATH=../lib ' + PARALLEL + '..'
  else
    sh 'cmake -DCMAKE_INSTALL_PREFIX:PATH=../lib ' + PARALLEL + '..'
  end

  puts "\n\nBuild #{LIB_NAME} Debug".green
  sh 'cmake --build . --config Debug --target install ' + PARALLEL
  FileUtils.cp "../lib/lib#{LIB_NAME}.a", "../lib/lib#{LIB_NAME}_debug.a"

  puts "\n\nBuild #{LIB_NAME} Release".green
  sh 'cmake --build . --config Release --target install ' + PARALLEL
  FileUtils.cd '..'

end

desc "compile for Visual Studio [default year=2017 bits=x64]"
task :build_win, [:year, :bits, :executable] do |t, args|
  args.with_defaults( :year => "2017", :bits => "x64", :executable => "no" )

  puts "\n\nPrepare #{LIB_NAME} project".green

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  tmp = " -DBITS=#{args.bits} -DYEAR=#{args.year} "

  if args.executable.to_s == "yes" then
    tmp += ' -DBUILD_EXECUTABLE=true'
  end

  tmp += ' -DCMAKE_INSTALL_PREFIX:PATH=..\lib ..'

  win32_64 = ''
  case args.bits
  when /x64/
    win32_64 = ' Win64'
  end

  case args.year
  when "2010"
    sh 'cmake -G "Visual Studio 10 2010' + win32_64 +'" ' + tmp
  when "2012"
    sh 'cmake -G "Visual Studio 11 2012' + win32_64 +'" ' + tmp
  when "2013"
    sh 'cmake -G "Visual Studio 12 2013' + win32_64 +'" ' + tmp
  when "2015"
    sh 'cmake -G "Visual Studio 14 2015' + win32_64 +'" ' + tmp
  when "2017"
    sh 'cmake -G "Visual Studio 15 2017' + win32_64 +'" ' + tmp
  else
    puts "Visual Studio year ``#{year}'' not supported!\n";
  end

  FileUtils.mkdir_p "../lib"
  sh 'cmake --build . --config Release --target install ' + PARALLEL

  libname = "#{LIB_NAME}_vs#{args.year}_#{args.bits}"

  puts "\n\nBuild #{LIB_NAME} Debug".green
  sh 'cmake --build . --config Debug --target install ' + PARALLEL
  FileUtils.cp "Debug/#{LIB_NAME}.lib", "../lib/#{libname}_debug.lib"

  puts "\n\nBuild #{LIB_NAME} Release".green
  FileUtils.cp "Release/#{LIB_NAME}.lib", "../lib/#{libname}.lib"

  FileUtils.cd '..'

end
