## Test and usage of tomography class

### GenerateTomoModels.py
    - It generates a true model to apply eikonal and get observed data
    - And generate a initial model reference to get a reconstructed models

### tomographyTest.cpp 
    - Simple model reconstruction
    - 7 iterations of least squares conjugate gradient  
    - First order Tikhonov with 1000 of regularization parameter     

### verifyTomographyTest.py 
    - Plot of reconstructed models per iteration
    - Plot of convergency curve

<h3 align="center"> Model reconstruction </h3>

![testFigure](https://user-images.githubusercontent.com/44127778/207900555-650c599b-6764-41cf-b3dc-43632d36a7b1.png)


<h3 align="center"> Objective function decreasing </h3>

![convergency](https://user-images.githubusercontent.com/44127778/207900646-255e1040-0225-4746-8a2e-8a570ec5524a.png)