%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

function [field_b] = clean_div(field_b, DC_fields)


global I_DC I_Tech I_Mesh

field_b_div = DC_fields.divB;
field_phi = DC_fields.phi;
field_laplace = DC_fields.laplace;
field_r = DC_fields.r;
field_u = DC_fields.u;
output = DC_fields.output;

field_b_div(:) = 0;
field_phi(:) = 0;
field_laplace(:) = 0;
field_r(:) = 0;
field_u(:) = 0;
output(:) = 0;

iterations = 0;
err_phi = Inf;

%Variable to log kernel runtime
iterations = 0;
kernel_time = 0;

g_range_projection = [I_Mesh('NODES_X') , I_Mesh('NODES_Y'), I_Mesh('NODES_Z')];



%% Init CG-Cycle

% Compute divergence of B === right-hand side
kernel_time = kernel_time + cl_run_kernel(I_Tech('device'), 'calc_div_divcleaning', I_DC('g_range'), I_DC('l_range'), field_b, field_b_div, 0);

% Compute initial residual
field_r(:) = field_b_div(:);
output(:) = 0;
kernel_time = kernel_time + cl_run_kernel(I_Tech('device'), 'dot_product', I_Tech('g_range'), I_Tech('l_range'), field_r, field_r, output, 0);
residual_norm2 = sum(output);
old_residual_norm2 = 1;
I_Results('divergence_norm_initial') = sqrt(residual_norm2);

% Set tolerance
tolerance = I_DC('absolute_error_threshold')*I_DC('absolute_error_threshold');

%%
%Start CG_cycle
while (residual_norm2 > tolerance) && (iterations < I_DC('max_iterations'))

    % update search direction
    beta_ = residual_norm2 / old_residual_norm2;
    kernel_time = kernel_time + cl_run_kernel(I_Tech('device'), 'calc_axpy', I_DC('g_range'), I_DC('l_range'), field_u, field_r, beta_, 0);

    % compute Laplace
    kernel_time = kernel_time + cl_run_kernel(I_Tech('device'), 'calc_laplace_divcleaning', I_DC('g_range'), I_DC('l_range'), field_u, field_laplace, 0);

    output(:) = 0;
    kernel_time = kernel_time + cl_run_kernel(I_Tech('device'), 'dot_product', I_Tech('g_range'), I_Tech('l_range'), field_u, field_laplace, output, 0);
    alpha_ = residual_norm2 / sum(output);

    % improve solution phi and residual
    kernel_time = kernel_time + cl_run_kernel(I_Tech('device'), 'calc_xpay', I_DC('g_range'), I_DC('l_range'), field_phi, field_u, alpha_, 0);
    kernel_time = kernel_time + cl_run_kernel(I_Tech('device'), 'calc_xpay', I_DC('g_range'), I_DC('l_range'), field_r, field_laplace, -alpha_, 0);

    old_residual_norm2 = residual_norm2;
    output(:) = 0;
    kernel_time = kernel_time + cl_run_kernel(I_Tech('device'), 'dot_product', I_Tech('g_range'), I_Tech('l_range'), field_r, field_r, output, 0);
    residual_norm2 = sum(output);

    % update counters
    iterations = iterations + 1;

end

%% Project magnetic field
kernel_time = kernel_time + cl_run_kernel(I_Tech('device'), 'projector_divcleaning', g_range_projection, I_DC('l_range'), field_b, field_phi, 0);

end
