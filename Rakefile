#
#
#

%w( colorize rake fileutils ).each do |gem|
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

LIB_NAME="Clothoids"

task :default => [:build]

desc "run tests"
task :run do
  sh "./bin/testBiarc"
  sh "./bin/testDistance"
  sh "./bin/testG2"
  sh "./bin/testG2plot"
  sh "./bin/testG2stat"
  sh "./bin/testG2stat2arc"
  sh "./bin/testG2statCLC"
  sh "./bin/testIntersect"
  sh "./bin/testPolyline"
  sh "./bin/testTriangle2D"
end

desc "run tests"
task :run_win do
  sh "./bin/Release/testBiarc"
  sh "./bin/Release/testDistance"
  sh "./bin/Release/testG2"
  sh "./bin/Release/testG2plot"
  sh "./bin/Release/testG2stat"
  sh "./bin/Release/testG2stat2arc"
  sh "./bin/Release/testG2statCLC"
  sh "./bin/Release/testIntersect"
  sh "./bin/Release/testPolyline"
  sh "./bin/Release/testTriangle2D"
end


desc "compile for UNIX/OSX [default GC='./GC',executable='no']"
task :build, [:gc_dir,:executable] do |t, args|

  args.with_defaults( :gc_dir => './GC', :executable => "no" )

  if args.gc_dir == './GC' then
    puts "\n\nBuild submodule GenericContainer".green
    FileUtils.rm_rf "GC"
    sh "git clone -b develop --depth 1 https://github.com/ebertolazzi/GenericContainer.git GC"
    FileUtils.cd "GC"
    FileUtils.cp "../CMakeLists-cflags.txt", "CMakeLists-cflags.txt"
    sh "rake build"
    FileUtils.cd ".."
  else
    puts "\n\nUse GenericContainer at #{args.gc_dir}".green
  end

  FileUtils.rm_rf   "lib"
  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  cmd = 'cmake -DCMAKE_INSTALL_PREFIX:PATH=../lib'
  if args.executable.to_s == "yes" then
    cmd += ' -DBUILD_EXECUTABLE=true'
  end
  cmd += " .."

  puts "\n\nPrepare #{LIB_NAME} project".green

  sh cmd

  puts "\n\nCompile #{LIB_NAME} Debug".green
  sh 'cmake --build . --config Debug --target install'
  FileUtils.cp "../lib/lib#{LIB_NAME}.a", "../lib/lib#{LIB_NAME}_debug.a"  

  puts "\n\nCompile #{LIB_NAME} Release".green
  sh 'cmake --build . --config Release --target install'
  FileUtils.cd '..'

end

desc "compile for Visual Studio [default year=2017 bits=x64 gc_dir=./GC executable=no]"
task :build_win, [:year, :bits, :gc_dir, :executable ] do |t, args|
  args.with_defaults( :year       => "2017",
                      :bits       => "x64",
                      :gc_dir     => './GC',
                      :executable => 'no', )

  if args.gc_dir == './GC' then
    puts "\n\nBuild submodule GenericContainer".green
    FileUtils.rm_rf "GC"
    sh "git clone -b develop --depth 1 https://github.com/ebertolazzi/GenericContainer.git GC"
    FileUtils.cd "GC"
    FileUtils.cp "../CMakeLists-cflags.txt", "CMakeLists-cflags.txt"
    sh "rake build"
    FileUtils.cd ".."
  else
    puts "\n\nUse GenericContainer at #{args.gc_dir}".green
  end

  puts "\n\nPrepare #{LIB_NAME} project".green
  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  tmp = " -DBITS=#{args.bits} -DYEAR=#{args.year} " +
        ' -DCMAKE_INSTALL_PREFIX:PATH=..\lib'

  if args.executable.to_s == "yes" then
    tmp += ' -DBUILD_EXECUTABLE=true'
  end
  tmp += " .."

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
    puts "Visual Studio year #{year} not supported!\n";
  end

  FileUtils.mkdir_p "../lib"

  puts "\n\nCompile #{LIB_NAME} Debug".green
  sh 'cmake --build . --config Debug --target install'
  FileUtils.cp "Debug/#{LIB_NAME}.lib", "../lib/#{LIB_NAME}_vs#{args.year}_#{args.bits}_debug.lib"

  puts "\n\nCompile #{LIB_NAME} Release".green
  sh 'cmake --build . --config Release  --target install'
  FileUtils.cp "Release/#{LIB_NAME}.lib", "../lib/#{LIB_NAME}_vs#{args.year}_#{args.bits}.lib"  

  FileUtils.cd '..'

end
