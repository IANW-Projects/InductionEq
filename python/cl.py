import pyopencl as cl

def compile_kernels(src_path, cmd_line, *args):
    global program, ctx
    src = ' '
    for entry in src_path:
        with open(entry, 'r') as io:
            src = src + io.read()
    prg = cl.Program(ctx, src).build(options=cmd_line)
    program = prg
    return prg

def create_ctx():
    global ctx, queue, PYOPENCL_COMPILER_OUTPUT, devices, group_size
    PYOPENCL_COMPILER_OUTPUT=1
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    devices = ctx.devices
    group_size = devices[0].max_work_group_size
    return ctx, queue

def run_kernel(names, g_range, l_range, *arguments):
    global queue, program
    inputbuffers = []
    for arg in arguments:
        inputbuffers.append(cl.Buffer(ctx, cl.mem_flags.COPY_HOST_PTR, hostbuf = arg))
    callstrings = []
    if l_range[0] == 0 and len(l_range) == 1:
        l_range = None
    else:
        l_range = tuple(l_range)

    if type(names) is list:
        for name in names:
            callstrings.append('program.' + name + '(queue, tuple(g_range) , l_range, *inputbuffers)')
    else:
        callstrings.append('program.' + names + '(queue, tuple(g_range), l_range, *inputbuffers)')

    for callstring in callstrings:
        event = eval(callstring)
        event.wait()
    i = 0
    for arg in arguments:
        cl.enqueue_copy(queue, arg, inputbuffers[i]).wait()
        inputbuffers[i].release()
        i = i +1
    return 0

def run_kernel_ref(name, g_range, l_range, *arg_refs):
    global queue, program
    callstring = 'program.' + name + '(queue, None, None, *inputbuffers)'
    eval(callstring)
    return

def allocate_buf(array):
    return cl.Buffer(ctx, cl.mem_flags.COPY_HOST_PTR, hostbuf = array)

def malloc(size):
    return cl.Buffer(ctx, size=size)




