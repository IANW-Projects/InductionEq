%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

function [field_b] = clean_div_CG(field_b, DC_fields)


global I_DC I_Tech I_Mesh I_Results I_RunOps

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
tic
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
elapsed_time = toc;
I_Results('residual_norm') = sqrt(residual_norm2);
I_Results('iterations') = iterations;

%% Project magnetic field
kernel_time = kernel_time + cl_run_kernel(I_Tech('device'), 'projector_divcleaning', g_range_projection, I_DC('l_range'), field_b, field_phi, 0);

I_Results('runtime') = elapsed_time;
I_Results('kernel_runtime') = kernel_time;

cl_run_kernel(I_Tech('device'), 'calc_div', g_range_projection, I_DC('l_range'), field_b, field_b_div, 0);
output(:) = 0;
cl_run_kernel(I_Tech('device'), 'norm2_S', I_Tech('g_range'), I_Tech('l_range'), field_b_div, output, 0);
I_Results('divergence_norm_final') = sqrt(sum(output));

% save fields
if I_RunOps('save_fields') == 1
    I_Results('field_b') = field_b;
    I_Results('field_phi') = field_phi;
    I_Results('field_b_div') = field_b_div;
end

%Optional plots

if ismember(lower(char(I_RunOps('plot_field_b'))),{'x','y','z','xy', 'xz', 'yz', 'xyz'})
    plot_2D(field_b, I_RunOps('plot_field_b'),...
        I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z'), 'B-Field after cleaning');
end

if ismember(lower(char(I_RunOps('plot_divergence'))),{'x','y','z','xy', 'xz', 'yz', 'xyz'})
    plot_2D(field_b_div, I_RunOps('plot_divergence'),...
        I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z'), 'Divergence after cleaning');
end

if ismember(lower(char(I_RunOps('plot_phi'))),{'x','y','z','xy', 'xz', 'yz', 'xyz'})
    plot_2D(field_phi, I_RunOps('plot_phi'),...
        I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z'), 'Phi');
end

end
