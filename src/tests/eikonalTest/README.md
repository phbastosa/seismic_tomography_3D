# Usage and tests for eikonal class

## eikonal.cpp
    - 3D refraction travel times accuracy study in high contrasted media
    
    - Experiment shown as expanded abstract and oral presentation in SimBGf 2022

    - Methods used:
        Podvin & Lecomte (1991)
        Jeong & Whitaker (2008)
        Noble, Gesret and Belayouni (2014)

    - Comparison with the analitical equation for refracted waves

## Paper can be found on SBGf [website](https://sbgf.org.br/mysbgf/eventos/expanded_abstracts/IX_SimBGf/session/M%C3%A9todos%20Geof%C3%ADsicos%20e%20Geof%C3%ADsica%20Computacional/3D%20refraction%20travel%20times%20accuracy%20study%20in%20high%20contrasted%20media.pdf).  


<h3 align="center"> Scheme implemented to verify accuracy of travel times </h3>

![modelGeometry](https://user-images.githubusercontent.com/44127778/206926084-dc082f40-5883-4e3c-9ccb-256134a5c31e.png)

<h3 align="center"> Comparison results of methods analized </h3>

![shot5](https://user-images.githubusercontent.com/44127778/206926089-ef1889de-97bc-4448-bb6c-4c8d172c87e3.png)

### Eikonal equation runtime test

## Asus x5570zd: Ryzen 5 3500 and GTX 1050 4GB mobile

    Model dimensions:
    Samples in x: 881 -> 22000 m
    Samples in y: 881 -> 22000 m
    Samples in z: 45 -> 1100 m

    Central shot applied in position:
    x = 11000 m
    y = 11000 m
    elevation = 0 m

    Circular geometry applied in configuration:
    x center = 11000 m
    y center = 11000 m
    elevation = 0 m
    offset = 10000 m

    ----------------- Run time ------------------

    Podvin & Lecomte (1991) time = 25.1027 s.
    Jeong & Whitaker (2008) time = 20.4786 s.
    Noble, Gesret and Belayouni (2014) time = 39.8723 s.