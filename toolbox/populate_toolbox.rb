#!/usr/bin/ruby -w

require 'fileutils'

FileUtils.rm_rf   "src"
FileUtils.mkdir_p "src"

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

lst = Dir["../GC/src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../GC/src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end

lst = Dir["../submodules/quarticRootsFlocke/src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../submodules/quarticRootsFlocke/src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end

FileUtils.rm_rf   "lib"
FileUtils.mkdir_p "lib"
lst = Dir["../matlab/*.m"]
lst.each do |filename|
  FileUtils.cp filename, "./lib/" + File.basename(filename);
end

FileUtils.rm_rf   "src_mex"
FileUtils.mkdir_p "src_mex"
lst = Dir["../src_mex/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src_mex/" + File.basename(filename);
end
lst = Dir["../src_mex/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src_mex/" + File.basename(filename);
end

FileUtils.rm_rf   "tests"
FileUtils.mkdir_p "tests"
lst = Dir["../tests-matlab/*.m*"]
lst.each do |filename|
  FileUtils.cp filename, "./tests/" + File.basename(filename);
end

lst = Dir["../doc/*"]
lst.each do |filename|
  FileUtils.cp filename, "./doc/" + File.basename(filename);
end

FileUtils.cp "../license.txt", "license.txt"
