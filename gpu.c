#include "bh.h"
#include "err_code.h"

char *load_cl_file(char *filename)
{
    char *source;
    int fd;

    source = (char *)calloc(1,8192);
    fd = open(filename, O_RDONLY);
    read(fd, source, 8192);
    return (source);
}

t_context *setup_context(void)
{
     cl_uint numPlatforms;
     t_context *c = (t_context *)calloc(1, sizeof(t_context));
     int err;

    // Find number of platforms
    err = clGetPlatformIDs(0, NULL, &numPlatforms);

    // Get all platforms
    cl_platform_id Platform[numPlatforms];
    err = clGetPlatformIDs(numPlatforms, Platform, NULL);
    checkError(err, "Getting platforms");

    // Secure a GPU
    for (int i = 0; i < numPlatforms; i++)
    {
        err = clGetDeviceIDs(Platform[i], DEVICE, 1, &(c->device_id), NULL);
        if (err == CL_SUCCESS)
        {
            break;
        }
    }

    if (c->device_id == NULL)
        checkError(err, "Finding a device");

    // Create a compute context
    c->context = clCreateContext(0, 1, &(c->device_id), NULL, NULL, &err);
    checkError(err, "Creating context");

    // Create a command queue
    c->commands = clCreateCommandQueue(c->context, (c->device_id), 0, &err);
    checkError(err, "Creating command queue");
    return (c);
}

cl_kernel   make_kernel(t_context *c, char *sourcefile, char *name)
{
    cl_kernel k;
    cl_program p;
    int err;
    char *source;

    source = load_cl_file(sourcefile);
    p = clCreateProgramWithSource(c->context, 1, (const char **) & source, NULL, &err);
    checkError(err, "Creating program");

    // Build the program
    err = clBuildProgram(p, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048];

        printf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(p, c->device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        exit(0);
    }

    // Create the compute kernel from the program
    k = clCreateKernel(p, name, &err);
    checkError(err, "Creating kernel");
    return (k);
}

size_t nearest_pow_2(size_t n)
{
    return(pow(2, floor(log2(n)) + 1));
}

size_t nearest_mult_256(size_t n)
{
    return (((n / 256) + 1) * 256);
}

// vvv N CROSS M vvv

t_ret crunch_NxM(cl_float4 *N, cl_float4 *V, cl_float4 *M, size_t ncount, size_t mcount, cl_float4 force_bias)
{
    srand ( time(NULL) ); //before we do anything, seed rand() with the current time
    t_context   *context;
    cl_kernel   k_nbody;
    int err;

    context = setup_context();
    k_nbody = make_kernel(context, "nxm.cl", "nbody");

    cl_float4 *output_p = (cl_float4 *)calloc(ncount, sizeof(cl_float4));
    cl_float4 *output_v = (cl_float4 *)calloc(ncount, sizeof(cl_float4));

    cl_float4 *FB = (cl_float4 *)calloc(ncount, sizeof(cl_float4));
    for (int i = 0; i < ncount; i++)
        FB[i] = force_bias;
    // printf("\n\n\nBEFORE\n\n\n");
    // for (int i = 0; i < ncount; i++)
    //     printf("nx %f ny %f nz %f nw %f vx %f vy %f vz %f\n", N[i].x, N[i].y, N[i].z, N[i].w, V[i].x, V[i].y, V[i].z);
    // for (int i = 0; i < mcount; i++)
    //     printf("mx %f my %f mz %f mw %f\n", M[i].x, M[i].y, M[i].z, M[i].w);

    //device-side data
    cl_mem      d_N_start;
    cl_mem      d_M;
    cl_mem      d_A;
    cl_mem      d_V_start;
    cl_mem      d_V_end;
    cl_mem      d_N_end;
    //inputs
    d_N_start = clCreateBuffer(context->context, CL_MEM_READ_ONLY, sizeof(cl_float4) * ncount, NULL, NULL);
    d_M = clCreateBuffer(context->context, CL_MEM_READ_ONLY, sizeof(cl_float4) * mcount, NULL, NULL);
    d_V_start = clCreateBuffer(context->context, CL_MEM_READ_ONLY, sizeof(cl_float4) * ncount, NULL, NULL);
    //outputs
    d_A = clCreateBuffer(context->context, CL_MEM_READ_WRITE, sizeof(cl_float4) * ncount, NULL, NULL);
    d_V_end = clCreateBuffer(context->context, CL_MEM_WRITE_ONLY, sizeof(cl_float4) * ncount, NULL, NULL);
    d_N_end = clCreateBuffer(context->context, CL_MEM_WRITE_ONLY, sizeof(cl_float4) * ncount, NULL, NULL);
    //copy over initial data to device locations
    cl_event load;
    clEnqueueWriteBuffer(context->commands, d_N_start, CL_TRUE, 0, sizeof(cl_float4) * ncount, N, 0, NULL, NULL);
    clEnqueueWriteBuffer(context->commands, d_M, CL_TRUE, 0, sizeof(cl_float4) * mcount, M, 0, NULL, NULL);
    clEnqueueWriteBuffer(context->commands, d_A, CL_TRUE, 0, sizeof(cl_float4) * ncount, FB, 0, NULL, NULL);
    clEnqueueWriteBuffer(context->commands, d_V_start, CL_TRUE, 0, sizeof(cl_float4) * ncount, V, 0, NULL, &load);
    //NB need an event handling here (can I just trust they'll finish in order?)
    clFinish(context->commands);

    size_t global = ncount;
    size_t mscale = mcount;
    size_t local = GROUPSIZE < global ? GROUPSIZE : global;
    float soften = SOFTENING;
    float timestep = TIME_STEP;
    float grav = G;
    size_t count = ncount;

    clSetKernelArg(k_nbody, 0, sizeof(cl_mem), &d_N_start);
    clSetKernelArg(k_nbody, 1, sizeof(cl_mem), &d_N_end);
    clSetKernelArg(k_nbody, 2, sizeof(cl_mem), &d_M);
    clSetKernelArg(k_nbody, 3, sizeof(cl_mem), &d_V_start);
    clSetKernelArg(k_nbody, 4, sizeof(cl_mem), &d_V_end);
    clSetKernelArg(k_nbody, 5, sizeof(cl_mem), &d_A);
    clSetKernelArg(k_nbody, 6, sizeof(cl_float4) * GROUPSIZE, NULL);
    clSetKernelArg(k_nbody, 7, sizeof(float), &soften);
    clSetKernelArg(k_nbody, 8, sizeof(float), &timestep);
    clSetKernelArg(k_nbody, 9, sizeof(float), &grav);
    clSetKernelArg(k_nbody, 10, sizeof(unsigned int), &global);
    clSetKernelArg(k_nbody, 11, sizeof(unsigned int), &mscale);
    
    //printf("global is %zu, local is %zu\n", global, local);
    //printf("going onto the GPU\n");
    
    cl_event tick;
    cl_event tock;
    err = clEnqueueNDRangeKernel(context->commands, k_nbody, 1, NULL, &global, &local, 0, NULL, &tick);
    checkError(err, "Enqueueing kernel");
    clEnqueueReadBuffer(context->commands, d_N_end, CL_TRUE, 0, sizeof(cl_float4) * count, output_p, 1, &tick, &tock);
    clEnqueueReadBuffer(context->commands, d_V_end, CL_TRUE, 0, sizeof(cl_float4) * count, output_v, 1, &tick, &tock);
    clFinish(context->commands); //eventually, don't do this, do some cool async thing where I return an event (ie tock)
    //printf("gpu finished\n");
    clReleaseMemObject(d_N_start);
    clReleaseMemObject(d_N_end);
    clReleaseMemObject(d_M);
    clReleaseMemObject(d_V_start);
    clReleaseMemObject(d_V_end);
    clReleaseMemObject(d_A);

    // printf("\n\n\nAFTER\n\n\n");
    // for (int i = 0; i < ncount; i++)
    // {
    //     printf("nx %f ny %f nz %f nw %f vx %f vy %f vz %f\n", output_p[i].x, output_p[i].y, output_p[i].z, output_p[i].w, output_v[i].x, output_v[i].y, output_v[i].z);
    //     printf("message %f\n", output_v[i].w);
    // }

    return ((t_ret){output_p, output_v}); //whoops it needs to return positions and vels. maybe this shouldn't integrate?
}

cl_float4 vec_to_f4(t_vector v)
{
    return ((cl_float4){v.x, v.y, v.z, v.w});
}

t_ret gpu_magic(t_body **N0, t_body **M0, t_vector force_bias)
{
    static long crosstotal;
    cl_float4 fb = {force_bias.x, force_bias.y, force_bias.z, force_bias.w};
    size_t ncount = count_bodies(N0);
    size_t mcount = count_bodies(M0);
    size_t npadding = nearest_mult_256(ncount) - ncount;
    size_t mpadding = nearest_mult_256(mcount) - mcount; //padding should be to nearest multiple of 256, not nearest power of 2.
    crosstotal += ncount * mcount;
    cl_float4 *N = (cl_float4 *)calloc(ncount + npadding, sizeof(cl_float4));
    cl_float4 *M = (cl_float4 *)calloc(mcount + mpadding, sizeof(cl_float4));
    cl_float4 *V = (cl_float4 *)calloc(ncount + npadding, sizeof(cl_float4));
    for (int i = 0; N0[i]; i++)
    {
        N[i] = vec_to_f4(N0[i]->position);
        V[i] = vec_to_f4(N0[i]->velocity);
    }
    for (int i = 0; M0[i]; i++)
        M[i] = vec_to_f4(M0[i]->position);
    printf("crunching on GPU. ncount is %zu (+%zu), mcount %zu (+%zu)\n", ncount,npadding, mcount, mpadding);
    return(crunch_NxM(N, V, M, ncount + npadding, mcount + mpadding, fb));
}