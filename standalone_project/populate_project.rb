#!/usr/bin/ruby -w

require 'fileutils'

FileUtils.rm_rf   "src"
FileUtils.mkdir_p "src"
FileUtils.mkdir_p "src/Clothoids"
FileUtils.mkdir_p "src/Utils/fmt"

FileUtils.rm_rf   "bin"
FileUtils.mkdir_p "bin"

lst = Dir["../src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../src/Clothoids/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Clothoids/" + File.basename(filename);
end

lst = Dir["../submodules/quarticRootsFlocke/src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../submodules/quarticRootsFlocke/src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end

lst = Dir["../submodules/UtilsLite/src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../submodules/UtilsLite/src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../submodules/UtilsLite/src/Utils/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/" + File.basename(filename);
end
lst = Dir["../submodules/UtilsLite/src/Utils/*.c*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/" + File.basename(filename);
end
lst = Dir["../submodules/UtilsLite/src/Utils/fmt/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/fmt/" + File.basename(filename);
end
lst = Dir["../submodules/UtilsLite/src/Utils/fmt/*.c*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/fmt/" + File.basename(filename);
end
lst = Dir["../submodules/UtilsLite/src/Utils/zstream/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/zstream/" + File.basename(filename);
end

lst = Dir["../doc/*"]
lst.each do |filename|
  FileUtils.cp filename, "./doc/" + File.basename(filename);
end

FileUtils.cp "../license.txt", "license.txt"
