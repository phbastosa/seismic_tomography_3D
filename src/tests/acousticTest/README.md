## Test and usage of tomography class

### buildModel.py
    - It generates a layered model to apply wave equation kernel

### acousticTest.cpp 
    - Full carpet OBN modeling with reciprocity 
    - Single node positioned on seafloor and 576 shots working as receivers    
    - Source Ricker with 30 Hz of max frequency, dt = 1 ms and 1 s of modeling
    - Cerjan ABC to damp all directions with low memory usase scheme
    - 8E2T difference finite operators 

### parametersTest.txt
    - To set all parameters needed in modeling

### verifyAcousticTest.py 
    - Plot of source wavelet
    - Plot of Cerjan ABC dampers
    - Plot of layered model
    - Plot of seismogram  

### Figures

![Screenshot from 2022-12-26 20-15-36](https://user-images.githubusercontent.com/44127778/209588373-1b2627ec-70a8-48a7-ac13-7eb9627069b8.png)

![Screenshot from 2022-12-26 20-15-48](https://user-images.githubusercontent.com/44127778/209588378-d5ccda5a-f50c-469f-9df3-21ad6d8cfb0d.png)

![Screenshot from 2022-12-26 20-19-52](https://user-images.githubusercontent.com/44127778/209588385-144bc151-175a-4b54-92bf-82774be6aa38.png)

![Screenshot from 2022-12-26 20-20-20](https://user-images.githubusercontent.com/44127778/209588394-10bf26f2-0da5-483f-ba77-18feccea1c65.png)