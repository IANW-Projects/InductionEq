%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

function [] = compute_numerical_solution(field_b, field_u, field_rho)

% Takes the initialized fields and advances the solution in time

global I_Mesh I_TI I_IEq I_DC I_Tech I_RunOps I_Results

% setup fields, settings etc.
[kernel_path_list, settings, IEq_fields, time_integrator_num_fields, DC_fields, RK_Step] = induction_eq.compute_numerical_solution_setup(field_b, field_u, field_rho);

% Initialise logging variables
kernel_runtime = 0;
RK_block_size = 15000;
if I_RunOps('save_integrals_over_time')
  energy = IEq_fields.energy;
  L2error_B = IEq_fields.L2error_B;
  L2error_divB = IEq_fields.L2error_divB;
end
norm2_output = IEq_fields.norm2_output;

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
        field_divB = IEq_fields.field_divB;

        if I_DC('divergence_cleaning') && ~I_RunOps('save_integrals_over_time')
            for step = 1:I_TI('num_steps')
                t = cl_run_kernel(I_Tech('device'), RK_Step, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_curlB_rho, ...
                                  field_u, field_rho, current_time, 0);
                kernel_runtime = kernel_runtime + t;

                [field_b] = induction_eq.clean_div(field_b, DC_fields);
            end
        elseif I_DC('divergence_cleaning') && I_RunOps('save_integrals_over_time')
            for step = 1:I_TI('num_steps')
                norm2_output(:) = 0;
                cl_run_kernel(I_Tech('device'), 'norm2', I_Tech('g_range'), I_Tech('l_range'), field_b, norm2_output, 0);
                energy(step) = sum(norm2_output);

                cl_run_kernel(I_Tech('device'), 'calc_div', I_IEq('g_range'), I_IEq('l_range'), field_b, field_divB, 0);
                norm2_output(:) = 0;
                cl_run_kernel(I_Tech('device'), 'norm2_S', I_Tech('g_range'), I_Tech('l_range'), field_divB, norm2_output, 0);
                L2error_divB(step) = sqrt(sum(norm2_output));

                cl_run_kernel(I_Tech('device'), 'analytical_b', I_IEq('g_range'), I_IEq('l_range'), field_b2, current_time, 0);
                norm2_output(:) = 0;
                cl_run_kernel(I_Tech('device'), 'norm2_diff', I_Tech('g_range'), I_Tech('l_range'), field_b, field_b2, norm2_output, 0);
                L2error_B(step) = sqrt(sum(norm2_output));

                t = cl_run_kernel(I_Tech('device'), RK_Step, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_curlB_rho, ...
                                  field_u, field_rho, current_time, 0);
                kernel_runtime = kernel_runtime + t;

                [field_b] = induction_eq.clean_div(field_b, DC_fields);
            end
        elseif I_RunOps('save_integrals_over_time')
            for step = 1:I_TI('num_steps')
                norm2_output(:) = 0;
                cl_run_kernel(I_Tech('device'), 'norm2', I_Tech('g_range'), I_Tech('l_range'), field_b, norm2_output, 0);
                energy(step) = sum(norm2_output);

                cl_run_kernel(I_Tech('device'), 'calc_div', I_IEq('g_range'), I_IEq('l_range'), field_b, field_divB, 0);
                norm2_output(:) = 0;
                cl_run_kernel(I_Tech('device'), 'norm2_S', I_Tech('g_range'), I_Tech('l_range'), field_divB, norm2_output, 0);
                L2error_divB(step) = sqrt(sum(norm2_output));

                cl_run_kernel(I_Tech('device'), 'analytical_b', I_IEq('g_range'), I_IEq('l_range'), field_b2, current_time, 0);
                norm2_output(:) = 0;
                cl_run_kernel(I_Tech('device'), 'norm2_diff', I_Tech('g_range'), I_Tech('l_range'), field_b, field_b2, norm2_output, 0);
                L2error_B(step) = sqrt(sum(norm2_output));

                t = cl_run_kernel(I_Tech('device'), RK_Step, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_curlB_rho, ...
                                  field_u, field_rho, current_time, 0);
                kernel_runtime = kernel_runtime + t;
            end
        else % I_DC('divergence_cleaning') == false == I_RunOps('save_integrals_over_time')
            num_steps_run = I_TI('num_steps');
            while num_steps_run > RK_block_size
                kernel_list = repmat(RK_Step, 1, RK_block_size);

                t = cl_run_kernel(I_Tech('device'), kernel_list, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_curlB_rho, ...
                                  field_u, field_rho, current_time, 0);
                kernel_runtime =  kernel_runtime + t;
                num_steps_run = num_steps_run - RK_block_size;
            end
            if num_steps_run > 0
                kernel_list = repmat(RK_Step, 1, num_steps_run);

                t = cl_run_kernel(I_Tech('device'), kernel_list, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_curlB_rho, ...
                                  field_u, field_rho, current_time, 0);
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
        field_divB = IEq_fields.field_divB;

        if I_DC('divergence_cleaning') && ~I_RunOps('save_integrals_over_time')
            for step = 1:I_TI('num_steps')
                t = cl_run_kernel(I_Tech('device'), RK_Step, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_b3, field_curlB_rho, ...
                                  field_u, field_rho, current_time, 0);
                kernel_runtime = kernel_runtime + t;

                [field_b] = induction_eq.clean_div(field_b, DC_fields);
            end
        elseif I_DC('divergence_cleaning') && I_RunOps('save_integrals_over_time')
            for step = 1:I_TI('num_steps')
                norm2_output(:) = 0;
                cl_run_kernel(I_Tech('device'), 'norm2', I_Tech('g_range'), I_Tech('l_range'), field_b, norm2_output, 0);
                energy(step) = sum(norm2_output);

                cl_run_kernel(I_Tech('device'), 'calc_div', I_IEq('g_range'), I_IEq('l_range'), field_b, field_divB, 0);
                norm2_output(:) = 0;
                cl_run_kernel(I_Tech('device'), 'norm2_S', I_Tech('g_range'), I_Tech('l_range'), field_divB, norm2_output, 0);
                L2error_divB(step) = sqrt(sum(norm2_output));

                cl_run_kernel(I_Tech('device'), 'analytical_b', I_IEq('g_range'), I_IEq('l_range'), field_b2, current_time, 0);
                norm2_output(:) = 0;
                cl_run_kernel(I_Tech('device'), 'norm2_diff', I_Tech('g_range'), I_Tech('l_range'), field_b, field_b2, norm2_output, 0);
                L2error_B(step) = sqrt(sum(norm2_output));

                t = cl_run_kernel(I_Tech('device'), RK_Step, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_b3, field_curlB_rho, ...
                                  field_u, field_rho, current_time, 0);
                kernel_runtime = kernel_runtime + t;

                [field_b] = induction_eq.clean_div(field_b, DC_fields);
            end
        elseif I_RunOps('save_integrals_over_time')
            for step = 1:I_TI('num_steps')
                norm2_output(:) = 0;
                cl_run_kernel(I_Tech('device'), 'norm2', I_Tech('g_range'), I_Tech('l_range'), field_b, norm2_output, 0);
                energy(step) = sum(norm2_output);

                cl_run_kernel(I_Tech('device'), 'calc_div', I_IEq('g_range'), I_IEq('l_range'), field_b, field_divB, 0);
                norm2_output(:) = 0;
                cl_run_kernel(I_Tech('device'), 'norm2_S', I_Tech('g_range'), I_Tech('l_range'), field_divB, norm2_output, 0);
                L2error_divB(step) = sqrt(sum(norm2_output));

                cl_run_kernel(I_Tech('device'), 'analytical_b', I_IEq('g_range'), I_IEq('l_range'), field_b2, current_time, 0);
                norm2_output(:) = 0;
                cl_run_kernel(I_Tech('device'), 'norm2_diff', I_Tech('g_range'), I_Tech('l_range'), field_b, field_b2, norm2_output, 0);
                L2error_B(step) = sqrt(sum(norm2_output));

                t = cl_run_kernel(I_Tech('device'), RK_Step, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_b3, field_curlB_rho, ...
                                  field_u, field_rho, current_time, 0);
                kernel_runtime = kernel_runtime + t;
            end
        else % I_DC('divergence_cleaning') == false == I_RunOps('save_integrals_over_time')
            num_steps_run = I_TI('num_steps');
            while num_steps_run > RK_block_size
                kernel_list = repmat(RK_Step, 1, RK_block_size);

                t = cl_run_kernel(I_Tech('device'), kernel_list, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_b3, field_curlB_rho, ...
                                  field_u, field_rho, current_time, 0);
                kernel_runtime =  kernel_runtime + t;
                num_steps_run = num_steps_run - RK_block_size;
            end
            if num_steps_run > 0
                kernel_list = repmat(RK_Step, 1, num_steps_run);

                t = cl_run_kernel(I_Tech('device'), kernel_list, I_IEq('g_range'), I_IEq('l_range'), ...
                                  field_b, field_b2, field_b3, field_curlB_rho, ...
                                  field_u, field_rho, current_time, 0);
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

%Calculate analytical solution and error
current_time(1) = I_TI('final_time');
cl_run_kernel(I_Tech('device'), 'analytical_b', I_IEq('g_range'), I_IEq('l_range'), field_b_ana, current_time, 0);

norm2_output(:) = 0;
cl_run_kernel(I_Tech('device'), 'norm2_diff', I_Tech('g_range'), I_Tech('l_range'), field_b, field_b_ana, norm2_output, 0);
I_Results('abs_err') = sqrt(sum(norm2_output));
if I_RunOps('save_integrals_over_time')
    L2error_B(I_TI('num_steps')+1) = I_Results('abs_err');
    I_Results('L2error_B_over_time') = L2error_B;
    I_Results('time') = linspace(0, I_TI('final_time'), length(L2error_B));
end

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
if I_RunOps('save_integrals_over_time')
    L2error_divB(I_TI('num_steps')+1) = I_Results('divergence_norm');
    I_Results('L2error_divB_over_time') = L2error_divB;
end

%Calculate energy
norm2_output(:) = 0;
cl_run_kernel(I_Tech('device'), 'norm2', I_Tech('g_range'), I_Tech('l_range'), field_b, norm2_output, 0);
I_Results('energy') = sum(norm2_output);
if I_RunOps('save_integrals_over_time')
    energy(I_TI('num_steps')+1) = I_Results('energy');
    I_Results('energy_over_time') = energy;
end


%Optional plots
if ismember(lower(char(I_RunOps('plot_numerical_solution'))),{'x','y','z','xy', 'xz', 'yz', 'xyz'})
    plot_2D(field_b, I_RunOps('plot_numerical_solution'),...
        I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z'), 'Numerical Solution');
end

if ismember(lower(char(I_RunOps('plot_analytical_solution'))),{'x','y','z','xy', 'xz', 'yz', 'xyz'})
    plot_2D(field_b_ana, I_RunOps('plot_analytical_solution'),...
        I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z'), 'Analytical Solution');
end

if ismember(lower(char(I_RunOps('plot_difference'))),{'x','y','z','xy', 'xz', 'yz', 'xyz'})
    field_diff = field_b(:,:) - field_b_ana(:,:);
    plot_2D(field_diff, I_RunOps('plot_difference'),...
        I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z'), 'Diff. of Analytical and Num. Solution');
end

if ismember(lower(char(I_RunOps('plot_divergence'))),{'x','y','z','xy', 'xz', 'yz', 'xyz'})
    plot_2D(field_divB, I_RunOps('plot_divergence'),...
        I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z'), 'Divergence of Num. Solution');
end


end
