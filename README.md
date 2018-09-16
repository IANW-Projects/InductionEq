# Induction equation

This is a suit for nummerically solving the nonlinear magnetic induction equation with OpenCL. It is intended as a tool to investigate different properties of the linear and nonlinear induction equation for user specified testcases. All numerical aspects are encapsulated in OpenCL kernels. The OpenCL host side is abstracted with the help of MatCL, an OpenCL interface for MathWorks Matlab. This provides the user with and intuitive and easy way of handling and processing input and output data without any intricate knowledge of OpenCL and allows for interactive development. Support for Julia is currently under development and will be provided in the near future. Usage of the OpenCL kernels is not limited to Matlab or Julia, on the contrary they can be used with any kind of host code or application that supports OpenCL. 

This project is the base of current mathematical and physics research. Hence, special emphasis is placed on the computational mathematics, meaning the form of discretization, methods of different order for spacial discretization and time integration and admissable boundary conditions for the linear and nonlinear induction equation. Results will be published in the near future. 

Exemplary testcases for the linear and nonlinear induction equation are provided in the folder examples. They inlcude a 3D convergence study for the linear induction equation, inspired by the study of Koley et al. (2012) and two testcases for the nonlinear equation. A more detailed describtion of how to use the script is provided in the files itself. 

This is still very much work in progress. If you have any questions or want to contribute feel free to contact us.

## Prerequisites & Setup

To run the examples the following must be installed:

 - OpenCL Driver (CPU, GPU) 
 - OpenCL SDK 
 - Mathworks Matlab
 - MatCL [Available on GitHub](https://github.com/philipheinisch/MatCL)
 
 For ease of use you can add `MatCL` to the search path of Matlab.
 
 ## License

This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.
