### Usage and tests for utils class

## utilsTest.cpp
    - Min Max function tested 
    - Read and write binary files tested
    - Find a parameter inside a text file texted for all types
    - Sparse matrix struct tested 
    - Simple linear system tested using least square conjugate gradient solver
    - Sparse derivative matrix generator tested
    - Smoothing filters tested
    - Trilinear interpolation tested using a known 3D scalar function

<h3 align="center"> Smoothing examples implemented </h3>

![Screenshot from 2022-12-10 18-02-39](https://user-images.githubusercontent.com/44127778/206875580-6efe79d9-3df2-426a-8ed1-848fcf2a1a8a.png)


<h3 align="center"> Function used on trilinear interpolation </h3>

![Screenshot from 2022-12-10 18-14-32](https://user-images.githubusercontent.com/44127778/206875535-9046df18-a821-4085-bc0e-8a982a08bb6a.png)

## Terminal outputs

    Min integer between 5 and 3 is: 3
    Min integer between 3 and 5 is: 3

    Max integer between 2 and 6 is: 6
    Max integer between 6 and 2 is: 6

    Min float between 5.3 and 3.5 is: 3.5
    Min float between 3.5 and 5.3 is: 3.5

    Max float between 2.6 and 6.2 is: 6.2
    Max float between 6.2 and 2.6 is: 6.2

    Min float between 5.3, 3.3 and 3.5 is: 3.3
    Max float between 3.5, 3.3 and 5.3 is: 5.3

    Min float between 5.5, 5.3, 3.3 and 3.5 is: 3.3

    Testing pick up parameter from file

    parameter1 = fileNameTest
    parameter2 = 1512
    parameter3 = 10.457
    parameter4 = 30.4, 50.2, 23.8
    parameter5 = 1

    Testing sparse matrix

    Objective: Solve the linear system below
    A =         x =     b =
        2 1 3       1        6
        5 4 2       1       11
        4 2 1       1        7
    Target functions:
        - Struct sparseMatrix
        - sparse_lscg()

    Solution: [1, 1, 1]

    A row = 0, A col = 0, A value = 2
    A row = 0, A col = 1, A value = 1
    A row = 0, A col = 2, A value = 3
    A row = 1, A col = 0, A value = 5
    A row = 1, A col = 1, A value = 4
    A row = 1, A col = 2, A value = 2
    A row = 2, A col = 0, A value = 4
    A row = 2, A col = 1, A value = 2
    A row = 2, A col = 2, A value = 1

    x = 1, 1, 1

    Testing derivative matrix generator

    L 0th order = identity
    row = 0, col = 0, value = 1
    row = 1, col = 1, value = 1
    row = 2, col = 2, value = 1
    row = 3, col = 3, value = 1
    row = 4, col = 4, value = 1

    L 1st order
    row = 0, col = 0, value = -1
    row = 0, col = 1, value = 1
    row = 1, col = 1, value = -1
    row = 1, col = 2, value = 1
    row = 2, col = 2, value = -1
    row = 2, col = 3, value = 1
    row = 3, col = 3, value = -1
    row = 3, col = 4, value = 1

    L 2nd order
    row = 0, col = 0, value = 1
    row = 0, col = 1, value = -2
    row = 0, col = 2, value = 1
    row = 1, col = 1, value = 1
    row = 1, col = 2, value = -2
    row = 1, col = 3, value = 1
    row = 2, col = 2, value = 1
    row = 2, col = 3, value = -2
    row = 2, col = 4, value = 1

    Testing smoothing filters

    Smoothed volumes written

    Testing trilinear interpolation

    F(z,x,y) = 1500 + 90*z + 30*x - 30*y

    Point: (z = 0.5, 0.7, 0.2)

    knowing: F(0,0,0), F(0,0,1), F(0,1,0), F(1,0,0)
             F(1,1,0), F(1,0,1), F(1,1,0), F(1,1,1)

    Exact solution = 1560
    Interpolated solution = 1560