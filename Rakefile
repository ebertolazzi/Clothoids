if File.exist?(File.expand_path('./cmake_utils/Rakefile_common.rb', File.dirname(__FILE__))) then
  require_relative "./cmake_utils/Rakefile_common.rb"
else
  require_relative "../Rakefile_common.rb"
end

CLEAN.include   ["./**/*.o", "./**/*.obj", "./bin/**/example*", "./build"]
CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }
CLOBBER.include []

desc "default task --> build"
task :default => :build

task :mingw_pacman do
  sh 'pacman -S development'
  sh 'pacman -S mingw-w64-x86_64-toolchain'
  sh 'pacman -S mingw-w64-x86_64-cmake'
  sh 'pacman -S mingw-w64-x86_64-ninja'
end

desc "compile for Visual Studio"
task :build_win do
  # check architecture
  case `where cl.exe`.chop
  when /(x64|amd64)\\cl\.exe/
    VS_ARCH = 'x64'
  when /(bin|x86|amd32)\\cl\.exe/
    VS_ARCH = 'x86'
  else
    raise RuntimeError, "Cannot determine architecture for Visual Studio".red
  end

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  puts "run CMAKE for CLOTHOIDS".yellow
  sh "cmake -G Ninja -DBITS:VAR=#{VS_ARCH} " + cmd_cmake_build() + ' ..'

  puts "compile with CMAKE for CLOTHOIDS".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL
  end

  FileUtils.cd '..'
end

desc "compile for OSX/LINUX/MINGW"
task :build_osx_linux_mingw do

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  puts "run CMAKE for CLOTHOIDS".yellow
  sh "cmake -G Ninja " + cmd_cmake_build + ' ..'

  puts "compile with CMAKE for CLOTHOIDS".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  else
    sh 'cmake  --build . --config Release  --target install '+PARALLEL
  end

  FileUtils.cd '..'
end

task :clean_osx_linux_mingw do
  FileUtils.rm_rf 'build'
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
end

task :build_osx   => :build_osx_linux_mingw do end
task :build_linux => :build_osx_linux_mingw do end
task :build_mingw => :build_osx_linux_mingw do end

task :clean_osx   => :clean_osx_linux_mingw do end
task :clean_linux => :clean_osx_linux_mingw do end
task :clean_mingw => :clean_osx_linux_mingw do end
task :clean_win   => :clean_osx_linux_mingw do end

desc 'pack for OSX/LINUX/MINGW/WINDOWS'
task :cpack do
  FileUtils.cd "build"
  puts "run CPACK for CLOTHOIDS".yellow
  sh 'cpack -C CPackConfig.cmake'
  sh 'cpack -C CPackSourceConfig.cmake'
  FileUtils.cd ".."
end
