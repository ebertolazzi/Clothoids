require_relative "./cmake_utils/Rakefile_common.rb"

CLEAN.include   ["./**/*.o", "./**/*.obj", "./bin/**/example*", "./build"]
CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }
CLOBBER.include []

desc "default task --> build"
task :default => :build

desc "compile for Visual Studio [default year=2017, bits=x64]"
task :build_win, [:year, :bits] do |t, args|

  args.with_defaults( :year => "2017", :bits => "x64" )

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  cmd_cmake = cmake_generation_command(args.bits,args.year) + cmd_cmake_build()

  puts "run CMAKE for CLOTHOIDS".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMAKE for CLOTHOIDS".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end

  FileUtils.cd '..'
end

desc "compile for OSX/LINUX/MINGW"
task :build_osx_linux_mingw do

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

  FileUtils.cd '..'
end

task :clean_osx_linux_mingw do
  FileUtils.rm_rf 'build'
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
end

desc "compile for OSX"
task :build_osx => :build_osx_linux_mingw

desc "compile for LINUX"
task :build_linux => :build_osx_linux_mingw

desc "compile for MINGW"
task :build_mingw => :build_osx_linux_mingw

task :clean_osx   => :clean_osx_linux_mingw
task :clean_linux => :clean_osx_linux_mingw
task :clean_mingw => :clean_osx_linux_mingw
task :clean_win   => :clean_osx_linux_mingw

desc 'pack for OSX/LINUX/MINGW/WINDOWS'
task :cpack do
  FileUtils.cd "build"
  puts "run CPACK for CLOTHOIDS".yellow
  sh 'cpack -C CPackConfig.cmake'
  sh 'cpack -C CPackSourceConfig.cmake'
  FileUtils.cd ".."
end
