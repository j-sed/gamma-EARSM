# gamma-EARSM
A three-equation model based around explicit algebraic Reynolds stress modeling is proposed.The model is built on an existing and tested model, namely, EARSM by Hellsten. As a third equation, it incorporates an intermittency equation for the prediction of laminar-turbulent transitions in a boundary layer from $\gamma-SST$ by Menter. The aim of this new formulation is to take the best from these two models, concretely the weak anisotropy and an ability to predict boundary layer transitions. The model was calibrated on a set of flat plate flows under various conditions and its performance was tested on simple cases of external and internal aerodynamics.
## Code
The code is developed from published code of @furstj, namely his EARSM and gamma-SST turbulence models. Please, check the code out at https://github.com/furstj/myTurbulenceModels 
## How to use
* Compile by executing wmake in root directory
* To use in simulation add one of these to controlDict file
* * libs ("libmyCompressibleTurbulenceModels.so");
* * libs ("libmyInompressibleTurbulenceModels.so");

