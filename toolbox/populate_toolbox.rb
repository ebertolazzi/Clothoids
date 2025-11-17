#!/usr/bin/ruby -w

require 'fileutils'

FileUtils.rm_rf "src"
Dir.glob("bin/*.mex*").each { |file| File.delete(file)}

FileUtils.cp_r  "../src/.",                                 "./src";
FileUtils.cp_r  "../submodules/quarticRootsFlocke/src/.",   "./src";
FileUtils.cp_r  "../submodules/UtilsLite/src/.",            "./src";
FileUtils.cp_r  "../submodules/GenericContainer/src/.",     "./src";
FileUtils.cp_r  "../submodules/GenericContainer/include/.", "./src";

# elimino dipendenze da Eigen
# FileUtils.rm_rf "./src/Eigen";
FileUtils.rm_rf "./src/Utils_Poly.cc";
FileUtils.rm_rf "./src/Utils_GG2D.cc";
FileUtils.rm_rf "./src/Utils_NelderMead.cc";
FileUtils.rm_rf "./src/Utils_HJPatternSearch.cc";
FileUtils.rm_rf "./src/GenericContainerSupport.cc";
FileUtils.rm_rf "./src/GenericContainerInterface_C.cc";
FileUtils.rm_rf "./src/GenericContainerSerialize.cc";
FileUtils.rm_rf "./src/GenericContainerTables.cc";

#lst = Dir["../doc/*"]
#lst.each do |filename|
#  FileUtils.cp filename, "./doc/" + File.basename(filename);
#end

FileUtils.cp "../license.txt", "license.txt"
