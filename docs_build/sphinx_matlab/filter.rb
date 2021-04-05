require 'fileutils'

index='./_build/html/genindex.html'
system("sed -i .bak 's/(C\+\+ class)/(MATLAB class)/g' #{index}");
system("sed -i .bak 's/(C\+\+ function)/(MATLAB function)/g' #{index}");
FileUtils.rm "#{index}.bak"
# 
# to_be_removed  = '<span class="pre">in<\/span> *<em><span class="pre">self<\/span><\/em>,*'
# to_be_removed0 = ' *<span class="pre">in<\/span> *'
# to_be_removed1 = '<em><span class="pre">self<\/span><\/em>,* *'
# to_be_removed2 = '<span class="pre">self)<\/span>'
# 
# to_be_removed3 = '<span class="pre">(in<\/span> *)'
# to_be_removed4 = '<span class="pre">(in<\/span> *<span class="pre">self)<\/span>'
# subs3 = '<span class="sig-paren">()<\/span>'

Dir.glob("./_build/html/api-matlab/*.html").each do |f|
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
