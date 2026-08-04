begin
  require 'colorize'
rescue LoadError
  class String
    def green  = self
    def yellow = self
    def magenta = self
  end
end

require 'fileutils'
require 'rake/clean'
require 'shellwords'


CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }

# Configurazione esterna (opzionale)
if File.exist?(File.expand_path('../Rakefile_configure.rb', File.dirname(__FILE__)))
  require_relative '../Rakefile_configure.rb'
elsif File.exist?(File.expand_path('../../Rakefile_configure.rb', File.dirname(__FILE__)))
  require_relative '../../Rakefile_configure.rb'
else
  COMPILE_DEBUG      = false
  COMPILE_DYNAMIC    = false
  COMPILE_EXECUTABLE = true
end

# ----------------------------------------------------------------------------
# Rilevamento OS
# ----------------------------------------------------------------------------
case RUBY_PLATFORM
when /darwin/
  OS = :mac
when /linux|cygwin/
  OS = :linux
when /msys/
  OS = :mingw
else
  OS = :win
end

def build_type
  COMPILE_DEBUG ? 'Debug' : 'Release'
end

def install_prefix
  File.expand_path('lib', __dir__)
end

def project_root
  File.expand_path(__dir__)
end

def build_dir
  File.join(project_root, 'build')
end

def cmake_generator
  'Ninja'
end

def cmake_configure_command(enable_tests: false)
  args = [
    "-G", cmake_generator,
    "-B", build_dir,
    "-DCMAKE_BUILD_TYPE=#{build_type}",
    "-DCMAKE_INSTALL_PREFIX=#{install_prefix}",
    "-DCMAKE_INSTALL_LIBDIR=lib",
    "-DCMAKE_INSTALL_INCLUDEDIR=include",
    "-DCMAKE_INSTALL_BINDIR=bin",
    "-DBUILD_SHARED_LIBS=#{COMPILE_DYNAMIC ? 'ON' : 'OFF'}",
    "-DBUILD_TESTING=#{enable_tests ? 'ON' : 'OFF'}",
    "-DCLOTHOIDS_INSTALL=ON",
    "-DCLOTHOIDS_BUILD_BENCHMARKS=OFF",
    "-DCLOTHOIDS_STRICT_WARNINGS=OFF",
    "-DCLOTHOIDS_POPULATE_TOOLBOX=OFF",
    "-DCLOTHOIDS_ALLOW_NETWORK_FETCH=OFF",
    "-DUTILS_UPDATE_3RDPARTY=OFF",
    project_root
  ]
  args.join(' ')
end

def cmake_build_command(target = nil)
  cmd = ["--build", build_dir, "--config", build_type]
  cmd += ["--target", target] if target
  cmd += ["--parallel"]
  cmd.join(' ')
end

def cleanup_duplicate_root_libraries
  lib_root = File.join(project_root, 'lib')
  return unless Dir.exist?(lib_root)

  FileUtils.rm_rf File.join(lib_root, 'Users')

  patterns = [
    'libClothoids*.a',
    'libClothoids*.dylib',
    'libClothoids*.so',
    'libClothoids*.dll'
  ]

  patterns.each do |pattern|
    Dir.glob(File.join(lib_root, pattern)).each do |path|
      next if File.dirname(path) == File.join(lib_root, 'lib')
      puts "Removing duplicate root library #{path}".yellow if respond_to?(:yellow)
      FileUtils.rm_f path
    end
  end
end

# ----------------------------------------------------------------------------
# Task principali
# ----------------------------------------------------------------------------
desc "default task --> build"
task :default => :build

desc "git clean reset"
task :git_clean do
  sh "git reset --hard"
  sh "git clean -d -x -f"
end

desc "Configure CMake (without building)"
task :configure do
  puts "Configuring CMake in #{build_dir} without tests...".green
  sh "cmake " + cmake_configure_command(enable_tests: false)
end

desc "Configure CMake with tests enabled"
task :configure_tests do
  puts "Configuring CMake in #{build_dir} with tests enabled...".green
  sh "cmake " + cmake_configure_command(enable_tests: true)
end

desc "Build and install the library (does NOT compile tests)"
task :build => :configure do
  puts "Compiling and installing library...".green
  sh "cmake " + cmake_build_command('install')
  cleanup_duplicate_root_libraries
end

desc "Compile all tests"
task :compile_tests => :configure_tests do
  puts "Compiling test executables...".green
  sh "cmake " + cmake_build_command('Clothoids_all_tests')
end

desc "Compile all targets with tests enabled"
task :compile_all => :configure_tests do
  puts "Compiling all targets with tests enabled...".green
  sh "cmake " + cmake_build_command
end

desc "Run all tests (compiles tests if needed)"
task :test => :compile_tests do
  puts "Running tests...".green
  sh "ctest --test-dir #{build_dir} --build-config #{build_type} --output-on-failure"
end

desc "Build and run all tests (alias for test)"
task :run => :test do
  # same as test
end

desc "Clean build artifacts (keeps lib/ and lib3rd/)"
task :clean_build do
  puts "Removing build directory only...".green
  FileUtils.rm_rf build_dir
end

desc "Full clean (removes build, lib, lib3rd, bin, object files)"
task :clean do
  puts "Full cleaning...".green
  FileUtils.rm_rf build_dir
  FileUtils.rm_rf File.join(project_root, 'lib')
  FileUtils.rm_rf File.join(project_root, 'lib3rd')
  FileUtils.rm_rf File.join(project_root, 'bin')
  FileUtils.rm_f Dir.glob(File.join(project_root, '**', '*.o'))
  FileUtils.rm_f Dir.glob(File.join(project_root, '**', '*.obj'))
end

desc "Package using CPack"
task :cpack do
  puts "Creating packages...".green
  Dir.chdir(build_dir) do
    sh "cpack -C #{build_type} CPackConfig.cmake"
    sh "cpack -C #{build_type} CPackSourceConfig.cmake"
  end
end

CLOBBER.include []
