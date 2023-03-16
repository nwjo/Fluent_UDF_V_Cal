#include "udf.h"
#include "metric.h"

DEFINE_ON_DEMAND(on_demnad_calc)
{
    #if !RP_HOST
    Domain *d;
    Thread *t;

    /* Integrate Velocity Magnitude */
    real sum_v;
    real x[ND_ND];
    real radius;
    cell_t c;
    d = Get_Domain(1);

    thread_loop_c(t, d)
    {
        if (FLUID_THREAD_P(t))
        {
            begin_c_loop_int(c,t)
            {
                radius = 0.0075;
                C_CENTROID(x,c,t);
                if (sqrt(pow(x[0],2) + pow(x[2],2))>=radius && sqrt(pow(x[0],2) + pow(x[2],2))< radius+0.0001 )
                {
                sum_v += sqrt(C_U(c,t)*C_U(c,t)+C_V(c,t)*C_V(c,t)+C_W(c,t)*C_W(c,t))*C_VOLUME[c,t];

                }
                else
                sum_v += 0;

            end_c_loop_int(c,t)
            }
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