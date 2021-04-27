%w(colorize rake fileutils).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

require_relative "./Rakefile_common.rb"

task :default => [:build]

TESTS = [
  "testBiarc",
  "testDistance",
  "testG2",
  "testG2plot",
  "testG2stat",
  "testG2stat2arc",
  "testG2statCLC",
  "testIntersect",
  "testPolyline",
  "testTriangle2D"
]

"run tests on linux/osx"
task :run do
  TESTS.each do |cmd|
    sh "./bin/#{cmd}"
  end
end

desc "run tests (Release) on windows"
task :run_win do
  TESTS.each do |cmd|
    sh "bin\\Release\\#{cmd}.exe"
  end
end

desc "run tests (Debug) on windows"
task :run_win_debug do
  TESTS.each do |cmd|
    sh "bin\\Debug\\#{cmd}.exe"
  end
end

desc "compile for Visual Studio [default year=2017, bits=x64]"
task :build_win, [:year, :bits] do |t, args|

  args.with_defaults( :year => "2017", :bits => "x64" )

  Rake::Task[:win_3rd].invoke(args.year,args.bits,args.lapack)

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmake_cmd = win_vs(args.bits,args.year)
  if COMPILE_EXECUTABLE then
    cmake_cmd += ' -DBUILD_EXECUTABLE:VAR=true '
  else
    cmake_cmd += ' -DBUILD_EXECUTABLE:VAR=false '
  end
  if COMPILE_DYNAMIC then
    cmake_cmd += ' -DBUILD_SHARED:VAR=true '
  else
    cmake_cmd += ' -DBUILD_SHARED:VAR=false '
  end

  FileUtils.mkdir_p "../lib/lib"
  FileUtils.mkdir_p "../lib/bin"
  FileUtils.mkdir_p "../lib/bin/"+args.bits
  FileUtils.mkdir_p "../lib/dll"
  FileUtils.mkdir_p "../lib/include"

  if COMPILE_DEBUG then
    sh cmake_cmd + ' -DCMAKE_BUILD_TYPE:VAR=Debug --loglevel=WARNING ..'
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh cmake_cmd + ' -DCMAKE_BUILD_TYPE:VAR=Release --loglevel=WARNING ..'
    sh 'cmake  --build . --config Release  --target install '+PARALLEL+QUIET
  end
  FileUtils.cd '..'
end

desc "compile for OSX"
task :build, [:os] do |t, args|

  args.with_defaults( :os => "osx" )

  case args.os
  when 'osx'
    Rake::Task[:osx_3rd].invoke()
  when 'linux'
    Rake::Task[:linux_3rd].invoke()
  end

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmake_cmd = "cmake "

  if COMPILE_EXECUTABLE then
    cmake_cmd += '-DBUILD_EXECUTABLE:VAR=true '
  else
    cmake_cmd += '-DBUILD_EXECUTABLE:VAR=false '
  end
  if COMPILE_DYNAMIC then
    cmake_cmd += '-DBUILD_SHARED:VAR=true '
  else
    cmake_cmd += '-DBUILD_SHARED:VAR=false '
  end

  if COMPILE_DEBUG then
    sh cmake_cmd + '-DCMAKE_BUILD_TYPE:VAR=Debug .. ' #--loglevel=WARNING ..'
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh cmake_cmd + '-DCMAKE_BUILD_TYPE:VAR=Release .. ' #--loglevel=WARNING ..'
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end

  FileUtils.cd '..'
end

desc "compile for LINUX"
task :build_linux do
  Rake::Task[:build].invoke("linux")
end

desc "compile for OSX"
task :build_osx do
  Rake::Task[:build].invoke("osx")
end

desc 'install third parties for osx'
task :osx_3rd do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
  FileUtils.cp 'CMakeLists-cflags.txt', 'submodules/Utils/CMakeLists-cflags.txt'
  FileUtils.cp 'CMakeLists-cflags.txt', 'submodules/quarticRootsFlocke/CMakeLists-cflags.txt'
  FileUtils.cp 'CMakeLists-cflags.txt', 'submodules/GenericContainer/CMakeLists-cflags.txt'
  FileUtils.cd 'submodules'
  puts "\n\nSUBMODULES (for CLOTHOIDS)\n\n".green
  sh "rake build_osx"
  FileUtils.cd '..'
end

desc 'install third parties for linux'
task :linux_3rd do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
  FileUtils.cp 'CMakeLists-cflags.txt', 'submodules/Utils/CMakeLists-cflags.txt'
  FileUtils.cp 'CMakeLists-cflags.txt', 'submodules/quarticRootsFlocke/CMakeLists-cflags.txt'
  FileUtils.cp 'CMakeLists-cflags.txt', 'submodules/GenericContainer/CMakeLists-cflags.txt'
  FileUtils.cd 'submodules'
  puts "\n\nSUBMODULES (for CLOTHOIDS)\n\n".green
  sh "rake build_linux"
  FileUtils.cd '..'
end

desc "compile for Visual Studio [default year=2017, bits=x64]"
task :win_3rd, [:year, :bits] do |t, args|
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
  FileUtils.cp 'CMakeLists-cflags.txt', 'submodules/Utils/CMakeLists-cflags.txt'
  FileUtils.cp 'CMakeLists-cflags.txt', 'submodules/quarticRootsFlocke/CMakeLists-cflags.txt'
  FileUtils.cp 'CMakeLists-cflags.txt', 'submodules/GenericContainer/CMakeLists-cflags.txt'
  args.with_defaults( :year => "2017", :bits => "x64" )
  FileUtils.cd 'submodules'
  puts "\n\nSUBMODULES (for CLOTHOIDS)\n\n".green
  sh "rake build_win[#{args.year},#{args.bits}]"
  FileUtils.cd '..'
end

task :clean_osx do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
end

task :clean_linux do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
end

task :clean_win do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
end
