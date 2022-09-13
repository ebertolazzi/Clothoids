puts "Setup submodules"
system('git submodule init')
system('git submodule update')
system('git submodule sync')
system('git submodule foreach --recursive git submodule init')
system('git submodule foreach --recursive git submodule update')
system('git submodule foreach --recursive git submodule sync')

puts ARGV

if ARGV.size() > 0 && ARGV[0] == "--last" then
  puts "\nUpdate submodules to last version"
  system('git submodule foreach --recursive git pull')
end
