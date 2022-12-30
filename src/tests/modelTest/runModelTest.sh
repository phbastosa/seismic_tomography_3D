utils="../../essentials/utils.cpp"
model="../../essentials/model.cpp"
flags="-fopenmp -fast -acc -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $utils $model modelTest.cpp $flags -o modelTest.exe

./modelTest.exe

rm *.o *.exe

python3 verifyModelTest.py

rm outputs/*.bin