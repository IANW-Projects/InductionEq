%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

% Use this script to test divergence cleaning with the Conjugate Gradient
% (CG) method

clear
close all
clc

addpath('../matlab')

%prepare_vars() initializes containers (maps) in which all variables (label, value) are
%stored. Variables are grouped together by purpose:
%I_Mesh: Contains all variables concerning the discrete mesh, e.g.
%stepsize, box size
%I_IEq: Contains all variables relevant for the induction equation, e.g.
%discretization form, switch for additional artificial dissipation
%I_DC: Contains all variables concerning divergence cleaning, e.g.
%discretization, error threshold
%I_Tech: Contains all variables of technical nature, e.g. the device,
%computation precision
%I_RunOps: Specifies the parameters of a computation, e.g. which variables
%will be saved, the testcase 

%Variables in capital letters are program specific defines which are set by
%openCL compile settings. If the value is a string in captial letters
%beginning with USE it acts as a switch, e.g. to switch between different
%discretizations. In case the key is written in capital letters the program
%define will be set to the corresponding value of the key.

divergence_cleaning.prepare_vars();

global I_Mesh I_DC I_Tech I_RunOps I_Results

% Mesh related parameters. The minimum number of nodes per direction is
% constrained by the number of boundary nodes, e.g. the 4th order method 
% has 2*4 boundary nodes which means the minimum number of nodes amounts to 8. 
N = uint16(40); I_Mesh('NODES_X') = N; I_Mesh('NODES_Y') = N; I_Mesh('NODES_Z') = uint16(12);
I_Mesh('XMIN') = -1.0; I_Mesh('XMAX') = 1.0;
I_Mesh('YMIN') = -1.0; I_Mesh('YMAX') = 1.0;
I_Mesh('ZMIN') = -1.0; I_Mesh('ZMAX') = 1.0;


% Divergence cleaning related variables
% Choose how the laplace operator will be discretized
I_DC('divergence_cleaning_form') = 'USE_LAPLACE_WIDE_STENCIL_DIRICHLET';
%USE_LAPLACE_WIDE_STENCIL_LNS USE_LAPLACE_WIDE_STENCIL_DIRICHLET USE_LAPLACE_NARROW_STENCIL_DIRICHLET
%Specify a threshold for the divergence norm
I_DC('absolute_error_threshold') = 1e-3;
%The divergence cleaner will exit after max_iterations even if the error
%threshold is not reached
I_DC('max_iterations') = 100;

% Technical details
% Specify the OpenCL device that will be used for computation
I_Tech('device') = 1;

%Switch between single floating point and double floating point precision
I_Tech('REAL') = 'double';
I_Tech('REAL4') = sprintf('%s4',I_Tech('REAL'));

%Compiler based optimizations
if strcmp(I_Tech('REAL'),'float')
    I_Tech('optimizations') = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only -cl-single-precision-constant';
else
    I_Tech('optimizations') = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only';
end

% Option relevant for run
% Defines the order 
I_RunOps('order') = 4; % 2, 4, 6
I_RunOps('operator_form') = 'classical'; % 'classical' or 'extended' operators
%I_RunOps('testcase') = 'rotation_2D'; %TODO: Implement testcases for
%divergence cleaning

% Optional plotting parameters (2D plots). 
% Choose the cross section with
% 'x', 'y', 'z'
% If you want to plot multiple cross sections use
% 'xy', 'xz', 'yz', 'xyz'
I_RunOps('plot_field_b') = 0;
I_RunOps('plot_divergence') = 0;
I_RunOps('plot_phi') = 0;
%If set to 1 the magnetic field will be saved to I_Results('field_b')
I_RunOps('save_fields') = 0;

%Initialize the magnetic field and additional fields needed for divergence
%cleaning.Also calculates and sets additional variables, e.g. stepsize,
%local work group size, etc...
[field_b, DC_fields] = divergence_cleaning.initialize();

%Takes the initialized fields reduces the divergence
field_b = divergence_cleaning.clean_div_CG(field_b, DC_fields);

fprintf('Divergence Norm:  %.15e (final), %.15e (initial)\n', ...
        I_Results('divergence_norm_final'), I_Results('divergence_norm_initial'));
fprintf('  Residual Norm:  %.15e\n', I_Results('residual_norm'));
fprintf('after %5d iterations (%f seconds)\n\n', ...
        I_Results('iterations'), I_Results('runtime'));
