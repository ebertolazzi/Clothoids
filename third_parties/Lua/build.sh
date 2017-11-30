DIR=lua-5.3.4
FILE="$DIR.tar.gz"
URL="http://www.lua.org/ftp/$DIR.tar.gz"
if [ -f $FILE ];
then
  echo "$FILE already downloaded"
else
  curl $URL > $FILE
fi
rm -rf lua
tar -zxvf $FILE
mv $DIR lua

cp -f CMakeLists.txt lua/CMakeLists.txt
mkdir -p lua/build
cd lua/build ; cmake .. ; make ; make install ; cd ../..
