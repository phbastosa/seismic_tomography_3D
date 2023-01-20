utils="../../src/essentials/utils.cpp"
model="../../src/essentials/model.cpp"
flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $utils $model modelTest.cpp $flags -o modelTest.exe

./modelTest.exe

rm *.o *.exe

python3 verifyModelTest.py

rm outputs/*.bin