#!/usr/bin/ruby -w

require 'fileutils'

FileUtils.rm_rf   "src"
FileUtils.mkdir_p "src"
FileUtils.mkdir_p "src/Clothoids"
FileUtils.mkdir_p "src/Utils/fmt"
FileUtils.mkdir_p "src/Utils/zstream"
FileUtils.mkdir_p "src/Utils/mingw-std-threads"

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

lst = Dir["../submodules/Utils/src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/Utils/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/Utils/*.c*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/Utils/fmt/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/fmt/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/Utils/fmt/*.c*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/fmt/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/Utils/zstream/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/zstream/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/Utils/mingw-std-threads/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/mingw-std-threads/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/lib3rd/include/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end

FileUtils.cp "../submodules/Utils/src/Utils/mingw-std-threads/LICENSE",
             "./src/Utils/mingw-std-threads/LICENSE";
FileUtils.cp "../submodules/Utils/src/Utils/mingw-std-threads/README.md",
             "./src/Utils/mingw-std-threads/README.md";

lst = Dir["../doc/*"]
lst.each do |filename|
  FileUtils.cp filename, "./doc/" + File.basename(filename);
end

FileUtils.cp "../license.txt", "license.txt"
