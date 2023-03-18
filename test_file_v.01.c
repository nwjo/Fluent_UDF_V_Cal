#include "udf.h"
#include "metric.h"

DEFINE_ON_DEMAND(on_demnad_calc)
{
    #if !RP_HOST
    Domain *d;
    Thread *t;

    /* Integrate Velocity Magnitude */
    real sum_v;
    real sum_vol;
    real sum_ave_v[ND_ND];
    real x[ND_ND];
    real radius;
    int i;
    //real theta=atan(x[1]/x[0]);
    //real radial_v;
    //real tangetial_v;
    cell_t c;
    d = Get_Domain(1);

    for(i = 0; i<4; i++) 
    {

        thread_loop_c(t, d)
        {

            if (FLUID_THREAD_P(t))
            {
                radius = 0.0085;
                radius += (real) i*0.001; 
                begin_c_loop_int(c,t)
                {
                    
                    C_CENTROID(x,c,t);
                    if (sqrt(pow(x[0],2) + pow(x[1],2))>=radius && sqrt(pow(x[0],2) + pow(x[1],2))< radius+2.0e-6 )
                    {
                        //radial_v = sqrt(C_U(c,t)*C_U(c,t)+C_V(c,t)*C_V(c,t));
                        //tangential_v = atan2(C_W(c,t)/C_U(c,t));

                        sum_v += sqrt(C_U(c,t)*C_U(c,t)+C_V(c,t)*C_V(c,t))*C_VOLUME(c,t);
                        sum_vol += C_VOLUME(c,t);
                    }
                    else
                    sum_v += 0;
                    sum_vol += 0;
                end_c_loop_int(c,t)
                }
            sum_ave_v[i] = sum_v/sum_vol;
            }
        #if RP_NODE
        sum_ave_v[i] = PRF_GRSUM1(sum_ave_v[i]);
        #endif

        }

    #if PARALLEL
    if(I_AM_NODE_ZERO_P)
    #endif
    printf("Volume Averaged Velocity: %g\n", sum_ave_v[i]);
    fflush(stdout);

    #endif
    }
}
