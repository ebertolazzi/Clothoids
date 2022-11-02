require 'yaml'

puts "Setup submodules"
system('git submodule init')
system('git submodule update')
system('git submodule sync')
system('git submodule foreach --recursive git submodule init')
system('git submodule foreach --recursive git submodule update')
system('git submodule foreach --recursive git submodule sync')

branches = YAML.load_file("./sub_branches.yaml")
branches.each do |dir, branch|
  puts "\n\nChecking out branch #{branch} in #{dir}"
  if branch.include? ":" then
    res = branch.split(':');
    system("(cd #{dir} && git reset --hard && git fetch --tags && git checkout '#{res[1]}' )")
  else
    system("(cd #{dir} && git reset --hard && git checkout #{branch} && git pull && git reset --hard)")
  end
end

puts ARGV

if ARGV.size() > 0 && ARGV[0] == "--last" then
  puts "\nUpdate submodules to last version"
  system('git submodule foreach --recursive git pull')
end

system('ruby submodules/Utils/setup.rb')
system('ruby submodules/quarticRootsFlocke/setup.rb')
