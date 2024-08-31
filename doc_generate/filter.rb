require 'fileutils'

index='./doxygen2/html/genindex.html'
system("sed -i .bak 's/(C\+\+ class)/(MATLAB class)/g' #{index}");
system("sed -i .bak 's/(C\+\+ function)/(MATLAB function)/g' #{index}");
#system("sed -i .bak 's/(C\+\+ class)/(class)/g' #{index}");
#system("sed -i .bak 's/(C\+\+ function)/(function)/g' #{index}");
#FileUtils.rm "#{index}.bak"

Dir.glob("./doxygen2/html/*.html").each do |f|
  puts "filter: #{f}"
  out = "";
  File.open(f,"r") do |file|
    file.each_line do |line|
      line.gsub!(/<span class="pre">in<\/span> *<em><span class="pre">(self|ignoredArg)<\/span><\/em>,*/,'');
      line.gsub!(/ *<span class="pre">in<\/span> */,'');
      line.gsub!(/<em><span class="pre">(self|ignoredArg)<\/span><\/em>,? */,'');
      line.gsub!(/<span class="pre">(self|ignoredArg)\)<\/span>/,'<span class="sig-paren">)</span>');
      line.gsub!(
        /<span class="pre">(in<\/span> *)/,
        '<span class="sig-paren">()</span>'
      )
      line.gsub!(
        /<span class="pre">\(in<\/span> *<span class="pre">(self|ignoredArg)\)<\/span>/,
        '<span class="sig-paren">()</span>'
      )
      out += line;
    end
  end
  File.open(f,"w") do |file|
    file.write(out)
  end
end
