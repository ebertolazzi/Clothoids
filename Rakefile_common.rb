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

require_relative "./Rakefile_conf.rb"

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

def win_vs( bits, year )

  tmp = " -DBITS:VAR=#{bits} -DYEAR:VAR=#{year} "

  if USE_NMAKE then
    tmp = 'cmake -G "NMake Makefiles" ' + tmp
  elsif USE_MINGW then
    tmp = 'cmake -G "MinGW Makefiles" ' + tmp
  elsif USE_MSYS then
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
