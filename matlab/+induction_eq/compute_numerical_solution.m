%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

function [] = compute_numerical_solution(field_b, field_u, field_rho)

% Takes the initialized fields and advances the solution in time

global I_Mesh I_TI I_IEq I_DC I_Tech I_RunOps I_Results

% setup fields, settings etc.
[kernel_path_list, settings, IEq_fields, time_integrator_num_fields, DC_fields, RK_Step] = induction_eq.compute_numerical_solution_setup(field_b, field_u, field_rho);

% Initialise logging variables
kernel_runtime = 0;
RK_block_size = 15000;

% Compile kernel
[~] = cl_run_kernel(I_Tech('device'), kernel_path_list, settings);

% Switch case for time integrator
tic
switch time_integrator_num_fields
    case 2
        field_b = IEq_fields.field_b;
        field_b2 = IEq_fields.field_b2;
        field_curlB_rho = IEq_fields.field_curlB_rho;
        field_u = IEq_fields.field_u;
        field_rho = IEq_fields.field_rho;
        current_time = IEq_fields.current_time;

        if I_DC('divergence_cleaning')
            for i = 1:I_TI('num_steps')
                t = cl_run_kernel(I_Tech('device'), RK_Step, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_curlB_rho, ...
                                  field_u, field_rho, current_time, [0 0 0 0 0]);
                kernel_runtime = kernel_runtime + t;

                [field_b] = induction_eq.clean_div(field_b, DC_fields);
            end
        else % I_DC('divergence_cleaning') == false
            num_steps_run = I_TI('num_steps');
            while num_steps_run > RK_block_size
                kernel_list = repmat(RK_Step, 1, RK_block_size);

                t = cl_run_kernel(I_Tech('device'), kernel_list, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_curlB_rho, ...
                                  field_u, field_rho, current_time, [0 0 0 0 0]);
                kernel_runtime =  kernel_runtime + t;
                num_steps_run = num_steps_run - RK_block_size;
            end
            if num_steps_run > 0
                kernel_list = repmat(RK_Step, 1, num_steps_run);

                t = cl_run_kernel(I_Tech('device'), kernel_list, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_curlB_rho, ...
                                  field_u, field_rho, current_time, [0 0 0 0 0]);
                kernel_runtime = kernel_runtime + t;
            end
        end

    case 3
        field_b = IEq_fields.field_b;
        field_b2 = IEq_fields.field_b2;
        field_b3 = IEq_fields.field_b3;
        field_curlB_rho = IEq_fields.field_curlB_rho;
        field_u = IEq_fields.field_u;
        field_rho = IEq_fields.field_rho;
        current_time = IEq_fields.current_time;

        if I_DC('divergence_cleaning')
            for i = 1:I_TI('num_steps')
                t = cl_run_kernel(I_Tech('device'), RK_Step, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_b3, field_curlB_rho, ...
                                  field_u, field_rho, current_time, [0 0 0 0 0]);
                kernel_runtime = kernel_runtime + t;

                [field_b] = induction_eq.clean_div(field_b, DC_fields);
            end
        else % I_DC('divergence_cleaning') == false
            num_steps_run = I_TI('num_steps');
            while num_steps_run > RK_block_size
                kernel_list = repmat(RK_Step, 1, RK_block_size);

                t = cl_run_kernel(I_Tech('device'), kernel_list, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_b3, field_curlB_rho, ...
                                  field_u, field_rho, current_time, [0 0 0 0 0]);
                kernel_runtime =  kernel_runtime + t;
                num_steps_run = num_steps_run - RK_block_size;
            end
            if num_steps_run > 0
                kernel_list = repmat(RK_Step, 1, num_steps_run);

                t = cl_run_kernel(I_Tech('device'), kernel_list, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_b3, field_curlB_rho, ...
                                  field_u, field_rho, current_time, [0 0 0 0 0]);
                kernel_runtime = kernel_runtime + t;
            end
        end

    otherwise
      error('Wrong value of time_integrator_num_fields = %d', time_integrator_num_fields)
end
runtime = toc;
fprintf('Elapsed time is %f seconds.\n', runtime);

% Save total runtime and kernel runtime
I_Results('runtime') = runtime;
I_Results('kernel_runtime') = kernel_runtime;

% save fields
if I_RunOps('save_fields')
    I_Results('field_b') = field_b;
end

field_b_ana = IEq_fields.field_b2;
field_divB = IEq_fields.field_divB;
norm2_output = IEq_fields.norm2_output;

%Calculate analytical solution and error
current_time(1) = I_TI('final_time');
cl_run_kernel(I_Tech('device'), 'analytical_b', I_IEq('g_range'), I_IEq('l_range'), field_b_ana, current_time,0);

norm2_output(:) = 0;
cl_run_kernel(I_Tech('device'), 'norm2_diff', I_Tech('g_range'), I_Tech('l_range'), field_b, field_b_ana, norm2_output, 0);
I_Results('abs_err') = sqrt(sum(norm2_output));

norm2_output(:) = 0;
cl_run_kernel(I_Tech('device'), 'norm2', I_Tech('g_range'), I_Tech('l_range'), field_b_ana, norm2_output, 0);
I_Results('rel_err') = I_Results('abs_err') / sqrt(sum(norm2_output));

if I_RunOps('save_fields')
  I_Results('field_b_ana') = field_b_ana;
end

%Calculate divergence
cl_run_kernel(I_Tech('device'), 'calc_div', I_IEq('g_range'), I_IEq('l_range'), field_b, field_divB, 0);
norm2_output(:) = 0;
cl_run_kernel(I_Tech('device'), 'norm2_S', I_Tech('g_range'), I_Tech('l_range'), field_divB, norm2_output, 0);
I_Results('divergence_norm') = sqrt(sum(norm2_output));

%Calculate energy
norm2_output(:) = 0;
cl_run_kernel(I_Tech('device'), 'norm2', I_Tech('g_range'), I_Tech('l_range'), field_b, norm2_output, 0);
I_Results('energy') = sum(norm2_output);

%Optional plots
if I_RunOps('plot_numerical_solution') == 1
    plot_2D(field_b, 'Z', I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z'), 'Numerical Solution');
end

if I_RunOps('plot_analytical_solution') == 1
    plot_2D(field_b_ana, 'Z', I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z'), 'Analytical Solution');
end

if I_RunOps('plot_difference') == 1
    field_diff = field_b(:,:) - field_b_ana(:,:);
    plot_2D(field_diff, 'Z', I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z'), 'Diff. of Analytical and Num. Solution');
end

if I_RunOps('plot_divergence') == 1
    plot_2D(field_divB, 'Z', I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z'), 'Divergence of Num. Solution');
end

end
