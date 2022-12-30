utils="../../essentials/utils.cpp"
flags="-fopenmp -fast -acc -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $utils utilsTest.cpp $flags -o utilsTest.exe

./utilsTest.exe

rm *.o *.exe

python3 verifyUtilsTest.py