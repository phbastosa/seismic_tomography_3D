utils="../../src/essentials/utils.cpp"
flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $utils utilsTest.cpp $flags -o utilsTest.exe

./utilsTest.exe

rm *.o *.exe

python3 verifyUtilsTest.py