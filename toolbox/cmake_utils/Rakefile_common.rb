#
#
#
%w(colorize fileutils pathname rubygems/package net/http zip zlib uri openssl).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

require 'rake/clean'
# avoid to remove file "core" (in Eigen inclusion)
CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }

require_relative "./Rakefile_configure.rb"

#    ___  ____
#   / _ \/ ___|
#  | | | \___ \
#  | |_| |___) |
#   \___/|____/
#
case RUBY_PLATFORM
when /darwin/
  OS = :mac
when /linux/
  OS = :linux
when /cygwin|mswin|mingw|bccwin|wince|emx/
  # detect windows shell not in msys2
  OS = :win
when /msys/
  # msys2 envirorment to compile with MINGW
  OS = :mingw
end
#
#    ____ __  __    _    _  _______
#   / ___|  \/  |  / \  | |/ / ____|
#  | |   | |\/| | / _ \ | ' /|  _|
#  | |___| |  | |/ ___ \| . \| |___
#   \____|_|  |_/_/   \_\_|\_\_____|
#
case OS
when :linux,:mac
  require 'etc'
  cmakeversion = %x( cmake --version ).scan(/\d+\.\d+\.\d+/).last
  mm = cmakeversion.split('.');
  if mm[0].to_i > 3 || (mm[0].to_i == 3 && mm[1].to_i >= 12) then
    PARALLEL = "--parallel #{Etc.nprocessors} "
    QUIET    = '-- --quiet '
  else
    PARALLEL = ''
    QUIET    = ''
  end
else
  PARALLEL = ''
  QUIET    = ''
end

def cmd_cmake_build()
  res = ""
  if COMPILE_EXECUTABLE then
    res += ' -DUTILS_ENABLE_TESTS:VAR=ON '
  else
    res += ' -DUTILS_ENABLE_TESTS:VAR=OFF '
  end
  if COMPILE_DYNAMIC then
    res += ' -DUTILS_BUILD_SHARED:VAR=ON '
  else
    res += ' -DUTILS_BUILD_SHARED:VAR=OFF '
  end
  if COMPILE_DEBUG then
    res += ' -DCMAKE_BUILD_TYPE:VAR=Debug --loglevel=STATUS '
  else
    res += ' -DCMAKE_BUILD_TYPE:VAR=Release --loglevel=STATUS '
  end
  return res
end

desc "default task --> build"
task :default => :build

desc "git submodule reset"
task :git_submodules do
  #sh "git clean -d -x -f"
  sh "git reset --hard"
  sh "git submodule update --init --recursive"
  sh "git submodule sync --recursive"
  sh "git submodule foreach --recursive git reset --hard"
  sh "git submodule foreach --recursive git clean -d -x -f"
  # estrae i sottomoduli alla corretta versione!
  sh "git submodule update --checkout --recursive"
end

desc "git clean reset"
task :git_clean do
  sh "git reset --hard"
  sh "git clean -d -x -f"
end


#   ____  _   _ _   _
#  |  _ \| | | | \ | |
#  | |_) | | | |  \| |
#  |  _ <| |_| | |\  |
#  |_| \_\\___/|_| \_|
#
desc "run tests"
task :run do
  puts "Run tests".green
  case OS
  when :mac,:linux,:mingw
    Dir.glob('./bin/*.exe').each do |exe|
      next unless File.exist?(exe)
      puts "execute #{exe}".yellow
      sh exe
    end
  when :win
    Dir.glob('./bin/*.exe').each do |exe|
      next unless File.exist?(exe)
      puts "execute #{exe}".yellow
      system(exe)
    end
  end
end

desc "run tests"
task :test do
  FileUtils.cd "build"
  sh 'ctest --output-on-failure'
  FileUtils.cd '..'
end

#   ____  _   _ ___ _     ____
#  | __ )| | | |_ _| |   |  _ \
#  |  _ \| | | || || |   | | | |
#  | |_) | |_| || || |___| |_| |
#  |____/ \___/|___|_____|____/
#
desc "build"
task :build do
  case OS
  when :mac
    puts "Build (osx)".green
    Rake::Task[:build_osx].invoke
  when :linux
    puts "Build (linux)".green
    Rake::Task[:build_linux].invoke
  when :win
    puts "Build (windows)".green
    Rake::Task[:build_win].invoke
  when :mingw
    puts "Build (mingw)".green
    Rake::Task[:build_mingw].invoke
  end
end

desc "clean"
task :clean do
  case OS
  when :mac
    puts "Clean (osx)".green
    Rake::Task[:clean_osx].invoke
  when :linux
    puts "Clean (linux)".green
    Rake::Task[:clean_linux].invoke
  when :win
    puts "Clean (windows)".green
    Rake::Task[:clean_win].invoke
  when :mingw
    puts "Clean (mingw)".green
    Rake::Task[:clean_mingw].invoke
  end
end

desc "default task --> build"
task :default => :build

#    ____ ___  __  __ ____ ___ _     _____ ____
#   / ___/ _ \|  \/  |  _ \_ _| |   | ____|  _ \
#  | |  | | | | |\/| | |_) | || |   |  _| | |_) |
#  | |__| |_| | |  | |  __/| || |___| |___|  _ <
#   \____\___/|_|  |_|_|  |___|_____|_____|_| \_\
#
def cmake_generation_command( bits, year )

  tmp = " -DBITS:VAR=#{bits} "

  case WIN_COMPILER
  when :nmake
    tmp = 'cmake -G "NMake Makefiles" ' + tmp
  when :mingw
    tmp = 'cmake -G "MinGW Makefiles" ' + tmp
  when :msys
    tmp = 'cmake -G "MSYS Makefiles" ' + tmp
  else
    win32_64  = ''
    win32_64_ = '-A Win32'
    case bits
    when /x64/
      win32_64  = ' Win64'
      win32_64_ = ' -A x64'
    end

    case year
    when "2010"
      tmp = 'cmake -G "Visual Studio 10 2010' + win32_64 +'" ' + tmp
    when "2012"
      tmp = 'cmake -G "Visual Studio 11 2012' + win32_64 +'" ' + tmp
    when "2013"
      tmp = 'cmake -G "Visual Studio 12 2013' + win32_64 +'" ' + tmp
    when "2015"
      tmp = 'cmake -G "Visual Studio 14 2015' + win32_64 +'" ' + tmp
    when "2017"
      tmp = 'cmake -G "Visual Studio 15 2017' + win32_64 +'" ' + tmp
    when "2019"
      tmp = 'cmake -G "Visual Studio 16 2019"' + win32_64_ + tmp
    else
      puts "Visual Studio year #{year} not supported!\n";
      return ""
    end
  end
  return tmp
end

#   _   _ _____ _______        _____  ____  _  __
#  | \ | | ____|_   _\ \      / / _ \|  _ \| |/ /
#  |  \| |  _|   | |  \ \ /\ / / | | | |_) | ' /
#  | |\  | |___  | |   \ V  V /| |_| |  _ <| . \
#  |_| \_|_____| |_|    \_/\_/  \___/|_| \_\_|\_\
#

#
# https://stackoverflow.com/questions/6934185/ruby-net-http-following-redirects/6934503
#
def url_resolve(
  uri_str,
  agent        = 'curl/7.43.0',
  max_attempts = 10,
  timeout      = 10
)
  attempts = 0
  cookie   = nil

  until attempts >= max_attempts
    attempts += 1

    url  = URI.parse(uri_str)
    http = Net::HTTP.new(url.host, url.port)
    http.open_timeout = timeout
    http.read_timeout = timeout
    path = url.path
    path = '/' if path == ''
    path += '?' + url.query unless url.query.nil?

    params = { 'User-Agent' => agent, 'Accept' => '*/*' }
    params['Cookie'] = cookie unless cookie.nil?
    request = Net::HTTP::Get.new(path, params)

    if url.instance_of?(URI::HTTPS)
      http.use_ssl = true
      http.verify_mode = OpenSSL::SSL::VERIFY_NONE
    end
    response = http.request(request)

    case response
    when Net::HTTPSuccess then
      break
    when Net::HTTPRedirection then
      location = response['Location']
      cookie   = response['Set-Cookie']
      new_uri  = URI.parse(location)
      uri_str  = if new_uri.relative? then url + location else new_uri.to_s end
    else
      raise 'Unexpected response: ' + response.inspect
    end
  end
  raise 'Too many http redirects' if attempts == max_attempts
  uri_str
  # response.body
end

def url_download( url_address, filename )
  if File.exist?(filename)
    puts "file: `#{filename}` already downloaded"
  else
    puts "downloading: #{filename}..."
    uri_str = url_resolve(url_address)
    uri     = URI( uri_str )
    File.binwrite( filename, Net::HTTP.get(uri) )
    puts "done"
  end
end

#
# https://stackoverflow.com/questions/856891/unzip-zip-tar-tag-gz-files-with-ruby
#
def extract_tgz( tar_gz_archive, destination = '.' )
  tar_longlink = '././@LongLink'
  Gem::Package::TarReader.new( Zlib::GzipReader.open tar_gz_archive ) do |tar|
    dest = nil
    tar.each do |entry|
      if entry.full_name == tar_longlink
        dest = File.join destination, entry.read.strip
        next
      end
      dest ||= File.join destination, entry.full_name
      if entry.directory?
        File.delete dest if File.file? dest
        FileUtils.mkdir_p dest, :mode => entry.header.mode, :verbose => false
      elsif entry.file?
        FileUtils.rm_rf dest if File.directory? dest
        FileUtils.mkdir_p File.dirname(dest), :mode => 0777, :verbose => false
        File.open dest, "wb" do |f| f.print entry.read end
        FileUtils.chmod entry.header.mode, dest, :verbose => false
      elsif entry.header.typeflag == '' #file?
        File.open dest, "wb" do |f| f.print entry.read end
        FileUtils.chmod entry.header.mode, dest, :verbose => false
      elsif entry.header.typeflag == '2' #Symlink!
        File.symlink entry.header.linkname, dest
      elsif entry.header.typeflag == 'g' && entry.full_name == "pax_global_header"
        puts "Skip entry: #{entry.full_name} type: #{entry.header.typeflag}."
      else
        puts "Unkown tar entry: #{entry.full_name} type: #{entry.header.typeflag}."
      end
      dest = nil
    end
  end
end

#
# https://stackoverflow.com/questions/19754883/how-to-unzip-a-zip-file-containing-folders-and-files-in-rails-while-keeping-the
#
def extract_zip( filename, destination_path='.' )
  if Zip.constants.include? :File
    zzfile = Zip::File
  else
    zzfile = Zip::ZipFile
  end
  zzfile.open(filename) do |zip_file|
    zip_file.each do |f|
      f_path=File.join(destination_path, f.name)
      FileUtils.mkdir_p(File.dirname(f_path))
      zip_file.extract(f, f_path) unless File.exist?(f_path)
    end
  end
end
