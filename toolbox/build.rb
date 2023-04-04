require_relative 'populate_toolbox.rb'
#require 'FileUtils'

FileUtils.rm_rf "build"
FileUtils.mkdir "build"
FileUtils.cd "build" do
  system("cmake -G Ninja ..");
  system("ninja");
end
