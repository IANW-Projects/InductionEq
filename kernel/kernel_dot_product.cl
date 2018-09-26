//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

// Contains kernel to compute the dot product of two vectors. 
// By setting the define REDUCE the reduction will be performed by the kernel with the help of atomics.

#ifdef fp64_atomics
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

inline void atomic_add_global(volatile global double *source, const double operand) {
  union {
    unsigned long intVal;
    double floatVal;
  } newVal;
  union {
    unsigned long intVal;
    double floatVal;
  } prevVal;

  do {
    prevVal.floatVal = *source;
    newVal.floatVal = prevVal.floatVal + operand;
  } while (atomic_cmpxchg((volatile global unsigned long *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}

#else

inline void atomic_add_global(volatile global float *source, const float operand) {
  union {
    unsigned int intVal;
    float floatVal;
  } newVal;
  union {
    unsigned int intVal;
    float floatVal;
  } prevVal;

  do {
    prevVal.floatVal = *source;
    newVal.floatVal = prevVal.floatVal + operand;
  } while (atomic_cmpxchg((volatile global unsigned int *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}

#endif // fp64_atomics

// Computes the dot product of vector a_vec and b_vec and stores the result in output.
// Optionally the reduction can be perfomred by the kernel.
kernel void dot_product(global REAL *a_vec, global REAL *b_vec, global REAL *output) {

  const int gid = get_global_id(0);
  const int lid = get_local_id(0);
  const int group_size = get_local_size(0);
  const int wgid = get_group_id(0);

  local REAL partial_dot[(uint)W_SIZE];

  partial_dot[lid] = a_vec[gid] * b_vec[gid];

  barrier(CLK_LOCAL_MEM_FENCE);

  for(uint i = group_size/2; i > 0; i = i /2) {
    if(lid < i) {
      partial_dot[lid] += partial_dot[lid + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  if (lid == 0) {
    #ifdef REDUCE
	    atomic_add_global(&(output[0]), partial_dot[0]);
    #else
    	output[wgid]= partial_dot[0];
    #endif
  }
}
