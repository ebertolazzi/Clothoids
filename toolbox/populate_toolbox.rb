#!/usr/bin/ruby -w

require 'fileutils'

FileUtils.rm_rf   "src"
#FileUtils.mkdir_p "src"
#FileUtils.mkdir_p "src/Clothoids"
#FileUtils.mkdir_p "src/Utils/fmt"
#FileUtils.mkdir_p "src/Utils/zstream"
#FileUtils.mkdir_p "src/Utils/mingw-std-threads"

FileUtils.rm_rf   "bin"
FileUtils.mkdir_p "bin"

FileUtils.cp_r "../src/.", "./src";
FileUtils.cp_r "../submodules/quarticRootsFlocke/src/.", "./src";
FileUtils.cp_r "../submodules/Utils/src/.",              "./src";

lst = Dir["../doc/*"]
lst.each do |filename|
  FileUtils.cp filename, "./doc/" + File.basename(filename);
end

FileUtils.cp "../license.txt", "license.txt"
