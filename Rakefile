%w(colorize rake fileutils).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

case RUBY_PLATFORM
when /darwin/
  OS = :mac
when /linux/
  OS = :linux
when /cygwin|mswin|mingw|bccwin|wince|emx/
  OS = :win
end

require_relative "./Rakefile_common.rb"

file_base = File.expand_path(File.dirname(__FILE__)).to_s+'/lib'

cmd_cmake_build = ""
if COMPILE_EXECUTABLE then
  cmd_cmake_build += ' -DEB_ENABLE_TESTS:VAR=ON '
else
  cmd_cmake_build += ' -DEB_ENABLE_TESTS:VAR=OFF '
end
if COMPILE_DYNAMIC then
  cmd_cmake_build += ' -DEB_BUILD_SHARED:VAR=ON '
else
  cmd_cmake_build += ' -DEB_BUILD_SHARED:VAR=OFF '
end
if COMPILE_DEBUG then
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Debug --loglevel=STATUS '
else
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Release --loglevel=STATUS '
end
cmd_cmake_build += " -DEB_INSTALL_LOCAL=ON "

FileUtils.cp './cmake/CMakeLists-cflags.txt', 'submodules/Utils/cmake/CMakeLists-cflags.txt'
FileUtils.cp './cmake/CMakeLists-cflags.txt', 'submodules/quarticRootsFlocke/cmake/CMakeLists-cflags.txt'

desc "default task --> build"
task :default => :build

desc "run tests"
task :test do
  FileUtils.cd "build"
  sh 'ctest --output-on-failure'
  FileUtils.cd '..'
end

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
    sh "./bin/#{cmd}" if File.exist?( "./bin/#{cmd}" )
  end
end

desc "run tests (Release) on windows"
task :run_win do
  TESTS.each do |cmd|
    sh "bin\\Release\\#{cmd}.exe" if File.exist?( "bin\\Release\\#{cmd}.exe" )
  end
end

desc "run tests (Debug) on windows"
task :run_win_debug do
  TESTS.each do |cmd|
    sh "bin\\Debug\\#{cmd}.exe"
  end
end

desc "build lib"
task :build do
  puts "UTILS build".green
  case OS
  when :mac
    Rake::Task[:build_osx].invoke
  when :linux
    Rake::Task[:build_linux].invoke
  when :win
    Rake::Task[:build_win].invoke
  end
end

desc "compile for Visual Studio [default year=2017, bits=x64]"
task :build_win, [:year, :bits] do |t, args|

  args.with_defaults( :year => "2017", :bits => "x64" )

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = win_vs(args.bits,args.year) + cmd_cmake_build

  puts "run CMAKE for CLOTHOIDS".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMAKE for CLOTHOIDS".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end

  if RUN_CPACK then
    puts "run CPACK for ROOTS".yellow
    sh 'cpack -C CPackConfig.cmake'
    sh 'cpack -C CPackSourceConfig.cmake'
  end

  FileUtils.cd '..'
end

desc "compile for OSX"
task :build_osx do

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = "cmake " + cmd_cmake_build

  puts "run CMAKE for CLOTHOIDS".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMAKE for CLOTHOIDS".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake  --build . --config Release  --target install '+PARALLEL+QUIET
  end

  if RUN_CPACK then
    puts "run CPACK for ROOTS".yellow
    sh 'cpack -C CPackConfig.cmake'
    sh 'cpack -C CPackSourceConfig.cmake'
  end

  FileUtils.cd '..'
end

desc "compile for LINUX"
task :build_linux => :build_osx

task :clean_osx do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
end

task :clean_linux => :clean_osx
task :clean_win => :clean_osx
