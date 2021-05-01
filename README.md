# Non-Gaussian Surface Simulations

**Author:** Sebastian Korsak

Special thanks to Max Pierini and his repository: https://github.com/maxdevblock/j_johnson_M, and my prof. V. Costantoudis who gave me the code in MATLAB.

The method we followed is described on the paper: *Numerical Simulation of 3D Rough Surfaces and Analysis
of Interfacial Contact Characteristics*, by Guoqing Yang, Baotong Li, Yang Wang and Jun Hong. Here we have some of the steps as are writen in this paper.

# How to use it

Just import the library,

```python
import non_gaussian_surfaces as seb
```

and type,

```python
seb.SAimage_fft_2(500,3,1,3,20,20,0.9,800.0,corr=True)
```

# Results

![3D Non Gaussian Rough Surface](https://imgur.com/jrfH5Fn)

![Binary Image Extracted from the 3D Gaussian Surface](https://imgur.com/nKt1YX2)

![Correlation Functions](https://imgur.com/7IA6yya)
