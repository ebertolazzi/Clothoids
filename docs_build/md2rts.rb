#
system("pandoc #{ARGV[0]}.md --from markdown --to rst -s -o #{ARGV[0]}.rst");
