import math
import numpy as np
#import pyopencl as cl
import cl

# from +InductionEq/preparse_vars.m
def preparse_vars():
	global I_Mesh, I_TI, I_IEq, I_DC, I_Tech, I_RunOps, I_Results
	I_Mesh = {'NODES_X':0, 'NODES_Y':0, 'NODES_Z':0, 'DX':0, 'DY':0, 'DZ':0, 'XMIN':0, 'XMAX':0, 'YMAX':0, 'YMIN':0, 'ZMIN':0, 'ZMAX':0}
	I_TI = {'cfl':0, 'final_time':0, 'time_integrator':'\0', 'DT':0, 'num_steps':0}
	I_IEq = {'form_uibj':0, 'form_source':0, 'form_ujbi':0, 'dissipation':'NONE', 'dissipation_form':'\0', 'HO_DISSIPATION_FACTOR':1, 'hall_term':'NONE', 'g_range':0, 'l_range':0, 'MP2MIN':0, 'MP2MAX':0, 'CMIN':0, 'CMAX':0}
	I_DC = {'use':0, 'form':'\0', 'BNODES':0, 'error_threshold':0, 'g_range':0, 'l_range':0,'max_iterations':0}
	I_Tech = {'device':0, 'REAL':'\0','num_nodes_pad':0, 'num_groups':0,'W_SIZE':0, 'g_range':0, 'l_range':0, 'optimizations':'\0', 'context' : None}
	I_RunOps = {'order':0, 'operator_form':'classical', 'testcase':'\0', 'variable_u':False, 'periodic':'NONE', 'plot_energy':0, 'plot_numerical_solution':0,'plot_analytical_solution':0, 'plot_difference':0, 'plot_divergence':0, 'save_fields':False, 'save_integrals_over_time':False}
	I_Results = {'abs_err':0, 'rel_err':0, 'field_b':0, 'field_b_ana':0, 'energy':0, 'divergence_norm':0, 'runtime':0, 'kernel_runtime':0, 'energy_over_time':0, 'L2error_B_over_time':0, 'L2error_divB_over_time':0, 'time':0 }

# from generate_settings.m

def generate_settings(map, keys):
    # Atomatically generates a formatted settings string that contains all
    # compiler optimizations and OpenCL defines of map with the given keys where 
    # keys is a cell array, e.g. settings_mesh = generate_settings(I_Mesh, {'DX'; 'DY'; 'DZ'}

	settings = '';
	for key in keys:
		if key in map:
			if not(type(map[key]) is float or type(map[key]) is int):
				if key==  'optimizations':
					settings = settings + map[key]
				elif 'USE' in map[key]:
					settings = settings + ' -D' +  map[key] + '=1'
				else:
					settings = settings + ' -D' + key + '=' + str(map[key])
			else:
				settings = settings + ' -D' + key + "=" + str(map[key])
		else:
			print('Unknown identifier ' +  keys[i] + '\n')
	return settings




# from +InductionEq/initialize.m
def initialize():
	global I_Mesh, I_TI, I_IEq, I_DC, I_Tech, I_RunOps


# *NEW* generate opencl context with PyOpenCl
	cl.create_ctx()

# a python float is a 64bit float => float means double
	if I_RunOps['periodic'] == 'USE_PERIODIC':
		DX = float(I_Mesh['XMAX'] - I_Mesh['XMIN']) / float(I_Mesh['NODES_X']);
		DY = float(I_Mesh['YMAX'] - I_Mesh['YMIN']) / float(I_Mesh['NODES_Y']);
		DZ = float(I_Mesh['ZMAX'] - I_Mesh['ZMIN']) / float(I_Mesh['NODES_Z']);
	else:
		DX = float(I_Mesh['XMAX'] - I_Mesh['XMIN']) / (float(I_Mesh['NODES_X'])-1);
		DY = float(I_Mesh['YMAX'] - I_Mesh['YMIN']) / (float(I_Mesh['NODES_Y'])-1);
		DZ = float(I_Mesh['ZMAX'] - I_Mesh['ZMIN']) / (float(I_Mesh['NODES_Z'])-1);


	I_Mesh['DX'] = DX; I_Mesh['DY'] = DY; I_Mesh['DZ'] = DZ;
	num_nodes = I_Mesh['NODES_X']*I_Mesh['NODES_Y']*I_Mesh['NODES_Z'];

# later ...
#    type = dev_type(I_Tech('device'));
#    if (strcmp(type{1},'CPU'))
#        group_size = 16;
#    else
#        if (num_nodes > lw_size(I_Tech('device')))
#            group_size = lw_size(I_Tech('device'));
#        else
#          group_size = 2^floor(log(num_nodes) / log(2));
#        end
#    end

# FOR MY GPU
	group_size = 1024

	group_size = float(group_size);
	num_nodes_pad = math.ceil(float(num_nodes)/group_size)*group_size
	num_groups = math.ceil(num_nodes_pad/group_size)

	I_Tech['num_nodes_pad'] = int(num_nodes_pad)

	I_Tech['num_groups'] = int(num_groups);
	I_Tech['W_SIZE'] = int(group_size);


# Initialise fields: magnetic, veclocity, and charge density
	if I_Tech['REAL'] == 'float':
		field_u_init = np.zeros((4, I_Tech['num_nodes_pad']), dtype=np.float32)
		field_b_init = np.zeros((4, I_Tech['num_nodes_pad']), dtype=np.float32)
		field_rho_init = np.zeros((1, I_Tech['num_nodes_pad']), dtype=np.float32)+1
	else:
		field_u_init = np.zeros((4, I_Tech['num_nodes_pad']), dtype=np.float64)
		field_b_init = np.zeros((4, I_Tech['num_nodes_pad']), dtype=np.float64)
		field_rho_init = np.zeros((1, I_Tech['num_nodes_pad']), dtype=np.float64)+1


# Generate settings needed for computation of the initial state for the
# specified testcase.
	settings_tech = generate_settings(I_Tech, ['REAL', 'REAL4', 'optimizations']);
	settings_mesh = generate_settings(I_Mesh, ['DX', 'DY', 'DZ', 'NODES_X', 'NODES_Y', 'NODES_Z','XMIN', 'XMAX', 'YMIN', 'YMAX','ZMIN', 'ZMAX'])
	settings = settings_tech + settings_mesh


	kernel_path_list =  []
	kernel_path_list.append('../include/'+I_RunOps['testcase']+'.h')
	kernel_path_list.append('../include/utils.h')
	kernel_path_list.append('../kernel/kernel_init.cl')

	cl.compile_kernels(kernel_path_list, settings)

    # Set global and local range for calculation of the inital state
    # The global and local range can either be a vector of type uint32 or
    # double. Double is supported due to the fact that by default Matlab
    # stores all variables as double.
	I_IEq['g_range'] = np.array([I_Mesh['NODES_X'], I_Mesh ['NODES_Y'], I_Mesh['NODES_Z']], dtype=np.int32)
	I_IEq['l_range'] = np.array([0], dtype=np.int32)

    # Compute initial state for magnetic field, velocity field and density
    # field
	cl.run_kernel('init_b', I_IEq['g_range'], I_IEq['l_range'], field_b_init)
	cl.run_kernel('init_u', I_IEq['g_range'], I_IEq['l_range'], field_u_init)
	cl.run_kernel('init_rho', I_IEq['g_range'], I_IEq['l_range'], field_rho_init)

    # Compute time step from CFL number, min step size and maximal velocity
	u_max = max(np.sqrt(np.power(field_u_init[1,:], 2) + np.power(field_u_init[2,:],2) + np.power(field_u_init[3,:],2)))
	dt = I_TI['cfl'] * min([DX, DY, DZ]) / float(u_max)
	num_steps = math.ceil(I_TI['final_time']/dt)
	dt = I_TI['final_time'] / num_steps

	I_TI['DT'] = dt
	I_TI['num_steps'] = num_steps

    # Global and local range for dot product and norm
	I_Tech['g_range'] = np.array([I_Tech['num_nodes_pad'], 1, 1], dtype=np.uint32)
	I_Tech['l_range'] = np.array([I_Tech['W_SIZE'], 1, 1], dtype=np.uint32)

    # Depending on the discretization of the laplace operator for
    # divergence cleaning set the number of boundary nodes
	if I_DC['divergence_cleaning_form'] == 'USE_LAPLACE_WIDE_STENCIL_LNS':
        	I_DC['BNODES'] = 0;
	elif (I_DC['divergence_cleaning_form'] == 'USE_LAPLACE_WIDE_STENCIL_DIRICHLET') or (I_DC['divergence_cleaning_form'] == 'USE_LAPLACE_NARROW_STENCIL_DIRICHLET'):
		I_DC['BNODES'] = 1;
	else:
		print('ERROR: Option I_DC(''divergence_cleaning_form'') = '+str(I_DC['divergence_cleaning_form']) + ' unknown.')

    # Set global and local range for divergence cleaning
	I_DC['g_range'] = np.array([I_Mesh['NODES_X']-2*I_DC['BNODES'], I_Mesh['NODES_Y']-2*I_DC['BNODES'], I_Mesh['NODES_Z']-2*I_DC['BNODES']], dtype=np.int32)
	I_DC['l_range'] = np.array([0], dtype=np.int32)

	return field_b_init, field_u_init, field_rho_init

def compute_numerical_solution_setup(field_b, field_u, field_rho):
	global I_Mesh, I_TI, I_IEq, I_DC, I_Tech, I_RunOps

	# Memory allocation for fields
	if I_Tech['REAL'] == 'float':
	# Initialize auxillary fields needed for time stepping
		field_b2 = np.zeros((4, I_Tech['num_nodes_pad']), dtype=np.float32)
		field_b3 = np.zeros((4, I_Tech['num_nodes_pad']), dtype=np.float32)
		field_divB = np.zeros((1, I_Tech['num_nodes_pad']), dtype=np.float32)
		field_curlB_rho = np.zeros((4, I_Tech['num_nodes_pad']), dtype=np.float32)
		current_time = np.zeros(2, dtype=np.float32) # contains two elements: time at the beginning of a step and of a stage
		norm2_output = np.zeros((1, I_Tech['num_groups']), dtype=np.float32)
	# Allocate memory for field used for divergence cleaning
		if I_DC['divergence_cleaning']:
        #Field for divergence of b-field and phi
			field_divB = np.zeros((1, I_Tech['num_nodes_pad']), dtype=np.float32)
			field_phi = np.zeros((1, I_Tech['num_nodes_pad']), dtype=np.float32)
        #Fields needed for CG
			field_laplace = np.zeros((1, I_Tech['num_nodes_pad']) , dtype=np.float32);
			field_r = np.zeros((1, I_Tech['num_nodes_pad']), dtype=np.float32);
			field_u_div_cleaning = np.zeros((1, I_Tech['num_nodes_pad']), dtype = np.float32);
    
		energy = np.zeros((1, I_TI['num_steps']+1), dtype=np.float32);
		L2error_B = np.zeros((1, I_TI['num_steps']+1), dtype=np.float32)
		L2error_divB = np.zeros((1, I_TI['num_steps']+1), dtype=np.float32);
	else:
       # Initialize auxillary fields needed for time stepping
                field_b2 = np.zeros((4, I_Tech['num_nodes_pad']), dtype=np.float64)
                field_b3 = np.zeros((4, I_Tech['num_nodes_pad']), dtype=np.float64)
                field_divB = np.zeros((1, I_Tech['num_nodes_pad']), dtype=np.float64)
                field_curlB_rho = np.zeros((4, I_Tech['num_nodes_pad']), dtype=np.float64)
                current_time = np.zeros(2, dtype=np.float64) # contains two elements: time at the beginning of a step and of a stage
                norm2_output = np.zeros((1, I_Tech['num_groups']), dtype=np.float64)
          # Allocate memory for field used for divergence cleaning
                if I_DC['divergence_cleaning']:
        #Field for divergence of b-field and phi
                        field_divB = np.zeros((1, I_Tech['num_nodes_pad']), dtype=np.float64)
                        field_phi = np.zeros((1, I_Tech['num_nodes_pad']), dtype=np.float64)
        #Fields needed for CG
                        field_laplace = np.zeros((1, I_Tech['num_nodes_pad']) , dtype=np.float64)
                        field_r = np.zeros((1, I_Tech['num_nodes_pad']), dtype=np.float64);
                        field_u_div_cleaning = np.zeros((1, I_Tech['num_nodes_pad']), dtype =np.float64)
    
                energy = np.zeros((1, I_TI['num_steps']+1), dtype=np.float64);
                L2error_B = np.zeros((1, I_TI['num_steps']+1), dtype=np.float64)
                L2error_divB = np.zeros((1, I_TI['num_steps']+1), dtype=np.float64);

	kernel_path_list = []

# Include header file containing the coefficients of the respective order
	if I_RunOps['order'] == 2:
		if I_RunOps['operator_form'] == 'extended':
			kernel_path_list.append('../include/2ndOrderExtended.h')
		else: # classical operators
			kernel_path_list.append('../include/2ndOrder.h')

	elif I_RunOps['order'] == 4:
		if I_RunOps['operator_form'] ==  'extended':
			kernel_path_list.append('../include/4thOrderExtended.h')
		else: # classical operators
			kernel_path_list.append('../include/4thOrder.h')

	elif I_RunOps['order'] == 6:
		if I_RunOps['operator_form'] == 'extended':
			kernel_path_list.append('../include/6thOrderExtended.h')
		else: # classical operators
			kernel_path_list.append('../include/6thOrder.h')

	else:
    		print('Specify order \n')

	kernel_path_list.append('../include/utils.h')
	kernel_path_list.append('../include/' + I_RunOps['testcase'] + '.h')
	kernel_path_list.append('../kernel/kernel_init.cl')

	kernel_path_list.append('../kernel/SBP_operator.cl')

	kernel_path_list.append('../kernel/kernel_operator.cl')

	# Include files for divergence cleaning if desired
	if I_DC['divergence_cleaning']:
		pass
    #Prepare temporary fields for divergence cleaning
#		DC_fields = struct('divB',field_divB,'phi',field_phi,'laplace', field_laplace,'r', field_r,'u', field_u_div_cleaning,'output', norm2_output);
#    kernel_path_list = [kernel_path_list, {'../kernel/kernel_div_cleaning.cl'}];
#    kernel_path_list = [kernel_path_list, {'../kernel/kernel_CG.cl'}];
#    settings_dc = generate_settings(I_DC, {'BNODES', 'divergence_cleaning_form'});
	else:
		DC_fields = 0
		settings_dc =''

# Include files for artificial dissipation if desired
	if I_IEq['dissipation'] ==  'USE_ARTIFICIAL_DISSIPATION':
		pass
#		kernel_path_list = [kernel_path_list, {'../include/artificial_dissipation.h'}];
#		kernel_path_list = [kernel_path_list, {'../kernel/artificial_dissipation.cl'}];
#		settings_dissipation = generate_settings(I_IEq, {'dissipation', 'dissipation_form', 'HO_DISSIPATION_FACTOR'});

	if I_IEq['dissipation_form'] == 'USE_ADAPTIVE_DISSIPATION':
		pass
#		settings_dissipation = strcat(settings_dissipation, generate_settings(I_IEq, {'MP2MIN'; 'MP2MAX'; 'CMIN'; 'CMAX'}));
	else:
    		settings_dissipation = ''

	kernel_path_list.append('../kernel/induction_Eq_volume.cl')
	kernel_path_list.append('../kernel/induction_Eq_surface.cl')
	kernel_path_list.append('../kernel/kernel_dot_product.cl')
	kernel_path_list.append('../kernel/kernel_time_integrator.cl')

# Specify compile settings
	settings_tech = generate_settings(I_Tech, ['REAL', 'REAL4', 'W_SIZE', 'optimizations'])

	settings_mesh = generate_settings(I_Mesh,  ['DX', 'DY', 'DZ','NODES_X', 'NODES_Y', 'NODES_Z', 'XMIN', 'XMAX', 'YMIN', 'YMAX','ZMIN', 'ZMAX'])

	settings_induction_eq = generate_settings(I_IEq, ['form_uibj', 'form_source', 'form_ujbi', 'hall_term'])

	settings_time_integration = generate_settings(I_TI, ['DT'])

	settings_runops = generate_settings(I_RunOps, ['periodic'])


	settings = settings_induction_eq + settings_time_integration + settings_tech + settings_mesh + settings_dissipation + settings_dc + settings_runops


	if I_TI['time_integrator'] ==  'SSPRK33':
        # 'SSPRK33_1', 'SSPRK33_2a', 'SSPRK33_2b', 'SSPRK33_3a', 'SSPRK33_3b'};
        # calc_curlB_rho = {'calc_curlB_rho_1_3args', 'calc_curlB_rho_2_3args', 'calc_curlB_rho_3_3args'};
		RK_Step = [];
		if I_IEq['hall_term'] ==  'USE_HALL':
			RK_Step.append('calc_curlB_rho_1_3args')
		if I_RunOps['variable_u']:
			RK_Step.append('calc_u_3_args')
		RK_Step.append('SSPRK33_1')

		RK_Step.append('SSPRK33_2a')
		if I_IEq['hall_term']== 'USE_HALL':
			RK_Step.append('calc_curlB_rho_2_3args')
		if I_RunOps['variable_u']:
			RK_Step.append('calc_u_3_args')
		RK_Step.append('SSPRK33_2b')

		RK_Step.append('SSPRK33_3a')
		if I_IEq['hall_term'] ==  'USE_HALL':
			RK_Step.append('calc_curlB_rho_3_3args')
		if I_RunOps['variable_u']:
			RK_Step.append('calc_u_3_args')
		RK_Step.append('SSPRK33_3b')

		RK_Step.append('calc_time_3_args')

		time_integrator_num_fields = 3;

# switch time_integrator_num_fields
	if time_integrator_num_fields == 2:
		IEq_fields = {'field_b':field_b, \
		'field_b2': field_b2, \
		'field_curlB_rho': field_curlB_rho, \
		'field_u': field_u, \
		'field_rho':field_rho, \
		'current_time':current_time, \
		'field_divB': field_divB, \
		'norm2_output':norm2_output, \
		'energy': energy, \
		'L2error_B': L2error_B, \
		'L2error_divB': L2error_divB}
	elif time_integrator_num_fields == 3:
		IEq_fields = {'field_b': field_b, \
		'field_b2': field_b2, \
		'field_b3': field_b3, \
		'field_curlB_rho': field_curlB_rho, \
		'field_u': field_u, \
		'field_rho': field_rho, \
		'current_time': current_time, \
		'field_divB': field_divB, \
                'norm2_output': norm2_output,\
                'energy': energy,\
                'L2error_B': L2error_B, \
                'L2error_divB': L2error_divB}
	else:
        	print('Wrong value of time_integrator_num_fields =' + str( time_integrator_num_fields))


	return kernel_path_list, settings, IEq_fields, time_integrator_num_fields, DC_fields, RK_Step

def compute_numerical_solution(field_b, field_u, field_rho):
	global I_mesh, I_TI, I_IEq, I_Tech, I_RunOps, I_Results

# setup fields, settings etc.
	(kernel_path_list, settings, IEq_fields, time_integrator_num_fields, DC_fields, RK_Step) = compute_numerical_solution_setup(field_b, field_u, field_rho)

# Initialise logging variables
	kernel_runtime = 0;
	RK_block_size = 15000;
	if I_RunOps['save_integrals_over_time']:
		energy = IEq_fields['energy']
		L2error_B = IEq_fields['L2error_B']
		L2error_divB = IEq_fields['L2error_divB']

	norm2_output = IEq_fields['norm2_output']

# Compile kernel
	cl.compile_kernels(kernel_path_list, settings);

# Switch case for time integrator

	if time_integrator_num_fields == 2:
		print("time_integrator_num_fields = 2 not ready yet \n")
	elif time_integrator_num_fields == 3:
		field_b = IEq_fields['field_b']
		field_b2 = IEq_fields['field_b2']
		field_b3 = IEq_fields['field_b3']
		field_curlB_rho = IEq_fields['field_curlB_rho']
		field_u = IEq_fields['field_u']
		field_rho = IEq_fields['field_rho']
		current_time = IEq_fields['current_time']
		field_divB = IEq_fields['field_divB']

		if False: #other versions
			pass
		else:
			num_steps_run = I_TI['num_steps'];
			while num_steps_run > RK_block_size:
				kernel_list = RK_Step * RK_block_size
				t = cl.run_kernel(kernel_list, I_IEq['g_range'], I_IEq['l_range'],\
				field_b, field_b2, field_b3, field_curlB_rho,\
				field_u, field_rho, current_time)
				kernel_runtime =  kernel_runtime + t;
				num_steps_run = num_steps_run - RK_block_size;

			if num_steps_run > 0:
				kernel_list = RK_Step *num_steps_run
				t = cl.run_kernel(kernel_list, I_IEq['g_range'], I_IEq['l_range'],\
				field_b, field_b2, field_b3, field_curlB_rho,\
				field_u, field_rho, current_time);
				kernel_runtime = kernel_runtime + t;

	else:
		print('Wrong value of time_integrator_num_fields = ' + str(time_integrator_num_fields) + ' \n')
	#runtime = toc;
	runtime = 0
	print('Elapsed time is' + str(runtime) + 'seconds.\n')

	# Save total runtime and kernel runtime
	I_Results['runtime'] = runtime;
	I_Results['kernel_runtime'] = kernel_runtime;

	# save fields
	if I_RunOps['save_fields']:
		I_Results['field_b'] = field_b;

	field_b_ana = IEq_fields['field_b2']

	#Calculate analytical solution and error
	current_time[1] = I_TI['final_time']
	cl.run_kernel('analytical_b', I_IEq['g_range'], I_IEq['l_range'], field_b_ana, current_time);

	norm2_output[:] = 0;

	cl.run_kernel('norm2_diff', I_Tech['g_range'], I_Tech['l_range'], field_b, field_b_ana, norm2_output);
	I_Results['abs_err'] = math.sqrt(sum(sum(norm2_output)));


	norm2_output[:] = 0;
	cl.run_kernel( 'norm2', I_Tech['g_range'], I_Tech['l_range'], field_b_ana, norm2_output);
	I_Results['rel_err'] = I_Results['abs_err'] / math.sqrt(norm2_output.sum());


	#Calculate divergence
	cl.run_kernel('calc_div', I_IEq['g_range'], I_IEq['l_range'], field_b, field_divB);
	norm2_output[:] = 0;
	cl.run_kernel('norm2_S', I_Tech['g_range'], I_Tech['l_range'], field_divB, norm2_output);
	I_Results['divergence_norm'] = math.sqrt(norm2_output.sum());


	#Calculate energy
	norm2_output[:] = 0;
	cl.run_kernel('norm2', I_Tech['g_range'], I_Tech['l_range'], field_b, norm2_output);
	I_Results['energy'] = norm2_output.sum();


