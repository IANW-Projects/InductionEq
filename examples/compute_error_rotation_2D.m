%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

%Use this script to specify the general conditions and the testcase you
%want to investigate. It calls functions to initializes all variables and fields relevant for
%the testcase and advance the numerical solution in time. Depending on the
%properties that shall be investigated additional functionalities can be
%enabled, e.g. divergence cleaning.

%clc;
clear;
close all;

addpath('../matlab')

%prepare_vars() initializes containers (maps) in which all variables (label, value) are
%stored. Variables are grouped together by purpose:
%I_Mesh: Contains all variables concerning the discrete mesh, e.g.
%stepsize, box size
%I_TI: Contains all variables concerning time integration, e.g. time
%integrator, timestep
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
%beginning with USE it acts as a switch, e.g. to enable artificial dissipation
%or to switch between different discretizations. In case the key is written
%in capital letters the program define will be set to the corresponding
%value of the key.
induction_eq.prepare_vars();


global I_Mesh I_TI I_IEq I_DC I_Tech I_RunOps I_Results

% Mesh related parameters. The minimum number of nodes per direction is
% constrained by the number of boundary nodes, e.g. the 4th order method 
% has 2*4 boundary nodes which means the minimum number of nodes amounts to 8. 
N = uint32(40);
I_Mesh('NODES_X') = N; I_Mesh('NODES_Y') = N; I_Mesh('NODES_Z') = uint32(12);
I_Mesh('XMIN') = -1.0; I_Mesh('XMAX') = 1.0;
I_Mesh('YMIN') = -1.0; I_Mesh('YMAX') = 1.0;
I_Mesh('ZMIN') = -1.0; I_Mesh('ZMAX') = 1.0;

% Time integration related variables
I_TI('cfl') = 0.95; %Define the Courant–Friedrichs–Lewy condition
I_TI('final_time') = 2*pi; 
%Chose the time integrator. Below is a list of up to date available
%options:
% SSPRK33, SSPRK104,
% KennedyCarpenterLewis2R54C, CalvoFrancoRandez2R64,
% CarpenterKennedy2N54, ToulorgeDesmet2N84F
I_TI('time_integrator') = 'CarpenterKennedy2N54';


% Induction equation related variables
%Specify how the three part of the linear induction equation shall be
%computed. 
I_IEq('form_uibj') = 'USE_UIBJ_PRODUCT'; % PRODUCT, SPLIT, CENTRAL
I_IEq('form_source') = 'USE_SOURCE_CENTRAL'; % CENTRAL, SPLIT, ZERO
I_IEq('form_ujbi') = 'USE_UJBI_CENTRAL'; % SPLIT, CENTRAL, PRODUCT
%Enable or disable Hall-Term
I_IEq('hall_term') = 'NONE'; % NONE, USE_HALL
%Enable or disable artificial dissipation
I_IEq('dissipation') = 'NONE'; %NONE, USE_ARTIFICIAL_DISSIPATION
%Specify what kind of artifical dissipation shall be used
%USE_ADAPTIVE_DISSIPATION, USE_FIRST_ORDER_DISSIPATION, USE_HIGH_ORDER_DISSIPATION
I_IEq('dissipation_form') = 'USE_ADAPTIVE_DISSIPATION'; 
%Additional parameters needed for adaptive dissipation. For typical values
%see Svärd et al. (2009).
I_IEq('MP2MIN') = 1;
I_IEq('MP2MAX') = 10;
I_IEq('CMIN') = 1;
I_IEq('CMAX') = 10;


% Divergence cleaning related variables
%Enable divergence cleaning (true) or disable divergence cleaning (false)
I_DC('divergence_cleaning') = false;
%Choose how the laplace operator will be discretized
I_DC('divergence_cleaning_form') = 'USE_LAPLACE_WIDE_STENCIL_LNS';
%USE_LAPLACE_WIDE_STENCIL_LNS USE_LAPLACE_WIDE_STENCIL_DIRICHLET USE_LAPLACE_NARROW_STENCIL_DIRICHLET
%Specify a threshold for the divergence norm
I_DC('absolute_error_threshold') = 1e-3;
%The divergence cleaner will exit after max_iterations even if the error
%threshold is not reached
I_DC('max_iterations') = 50;


% Technical details
% Specify the OpenCL device that will be used for computation
I_Tech('device') = 1;

%Switch between single floating point and double floating point precision
I_Tech('REAL') = 'double'; % float
I_Tech('REAL4') = sprintf('%s4',I_Tech('REAL')); %Vector datatype

%Compiler based optimizations
if strcmp(I_Tech('REAL'),'float')
    I_Tech('optimizations') = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only -cl-single-precision-constant';
else
    I_Tech('optimizations') = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only';
end

% Options relevant for a run
% Defines the order 
I_RunOps('order') = 4; % 2, 4, 6
% Specify the testcase. The name of the testcase has to be equal to the
% name of a header file which contains a function describing the initial state.
% Example testcases are:
% rotation_2D, rotation_3D, alfven_periodic_2D, hall_travelling_wave, hall_periodic
I_RunOps('testcase') = 'rotation_2D';
I_RunOps('variable_u') = false; % must be set to true if a variable velocity is used
I_RunOps('periodic') = 'NONE'; % 'NONE', 'USE_PERIODIC'; must be set to 'USE_PERIODIC'
                                       % if periodic boundary conditions should be used
% Optional plotting parameters (2D plots). 
% Choose the cross section with
% 'x', 'y', 'z'
% If you want to plot multiple cross sections use
% 'xy', 'xz', 'yz', 'xyz'
I_RunOps('plot_numerical_solution') = 'xyz';
I_RunOps('plot_analytical_solution') = 'xyz';
I_RunOps('plot_difference') = 'xyz';
I_RunOps('plot_divergence') = 'xyz';
%If set to 1 the magnetic field will be saved to I_Results('field_b')
I_RunOps('save_fields') = false;

%Initialize the magnetic field, the velocity field and the density field
%according to the specified testcase. Also calculates and sets additional
%variables, e.g. stepsize, timestep, local work group size, etc...
[field_b_init, field_u_init, field_rho_init] = induction_eq.initialize();

fprintf('Testcase: %s \nOrder: %d \nTime integrator: %s\nDT: %.16e   N_STEPS: %5d   FINAL_TIME: %.16e\nDX: %.16e   NODES_X: %5d\nDY: %.16e   NODES_Y: %5d\nDZ: %.16e   NODES_Z: %5d\nDissipation: %s \nDivergence cleaning: %d \nREAL: %s\n',...
        I_RunOps('testcase'), I_RunOps('order'), I_TI('time_integrator'), I_TI('DT'), I_TI('num_steps'), I_TI('final_time'), I_Mesh('DX'), I_Mesh('NODES_X'), I_Mesh('DY'), I_Mesh('NODES_Y'), I_Mesh('DZ'), I_Mesh('NODES_Z'), I_IEq('dissipation'), I_DC('divergence_cleaning'), I_Tech('REAL'));
%%
% Takes the initialized fields and advances the solution in time 
induction_eq.compute_numerical_solution(field_b_init, field_u_init, field_rho_init);

fprintf('Divergence Norm: %.15e    Total Energy: %.15e\nRelative Error: %.15f %%\n\n', ...
        I_Results('divergence_norm'), I_Results('energy'), 100*I_Results('rel_err'))
