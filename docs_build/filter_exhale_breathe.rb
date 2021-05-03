require 'fileutils'

dir     = ARGV[0]
project = ARGV[1]

puts "dir     = #{dir}"
puts "project = #{project}"

Dir.glob(dir+"/**/*.rst").each do |f|
  puts "filter: #{f}"
  out = "";
  File.open(f,"r") do |file|
    file.each_line do |line|
      line.gsub!(/(\.\. *doxygen.*::.*)/) { |m| m+"\n   :project: #{project}" };
      out += line;
    end
  end
  File.open(f,"w") do |file|
    file.write(out)
  end
end
