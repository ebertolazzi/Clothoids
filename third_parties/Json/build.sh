# os_mac, os_linux
URL="https://github.com/Tencent/rapidjson.git"
if [ -d rapidjson ];
then
  echo "rapidjson already downloaded"
else
  git clone --depth 1 $URL
fi

cp -r rapidjson/include ../../lib3rd
