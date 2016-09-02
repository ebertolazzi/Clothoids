require 'fileutils'

tmp = `gcc -dumpversion`
if tmp <= "4.9.0" then
  File.open("tmp.txt", "w") do |fout|
    File.open("src/GenericContainerConfig.hh","r") do |fin|
      fin.each do |line|
        line.gsub!(/^\#define GENERIC_CONTAINER_USE_REGEX/,"//#define GENERIC_CONTAINER_USE_REGEX")
        fout.puts line
        puts line
      end
    end
  end
  FileUtils.mv("tmp.txt","src/GenericContainerConfig.hh")
end
