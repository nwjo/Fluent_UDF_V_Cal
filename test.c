#include "udf.h"

DEFINE_ON_DEMAND(on_demnad_calc)
{
    #if !RP_HOST
    Domain *d;
    Thread *t;

    /* Integrate Velocity Magnitude */
    real sum_v=0.;
    cell_t c;
    d = Get_Domain(1);

    thread_loop_c(t, d)
    {
        if (FLUID_THREAD_P(t))
        {
            begin_c_loop(c,t)
            sum_v += sqrt(C_U(c,t)*C_U(c,t)+C_V(c,t)*C_V(c,t)+C_W(c,t)*C_W(c,t));
            end_c_loop(c,t)
        }

        #if RP_NODE
        sum_v = PRF_GRSUM1(sum_v);
        #endif
    }

    #if PARALLEL
    if(I_AM_NODE_ZERO_P)
    #endif
    printf("Sum of Velocity: %g\n", sum_v);
    fflush(stdout);

    #endif
}