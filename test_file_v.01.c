#include "udf.h"
#include "metric.h"

DEFINE_ON_DEMAND(on_demnad_calc)
{
    #if !RP_HOST
    Domain *d;
    Thread *t;

    /* Integrate Tangential Velocity */
    real sum_v;
    real sum_vol;
    real sum_ave_v[ND_ND];
    real x[ND_ND];
    real radius;
    real radial_coordinate;
    real theta_speed;
    real tangential_v;
    int i;

    cell_t c;
    d = Get_Domain(1);

        thread_loop_c(t, d)
        {

            if (FLUID_THREAD_P(t))
            {

                for(i = 1; i<5; i++) 
                {

                radius = 0.0085;
                radius += i*0.001+radius;

                begin_c_loop_int(c,t)
                {
                    C_CENTROID(x,c,t);
                    if (sqrt(pow(x[0],2) + pow(x[1],2))>=radius && sqrt(pow(x[0],2) + pow(x[1],2))< radius+2.0e-6 )
                    {
                        radial_coordinate = sqrt(pow(x[0],2) + pow(x[1],2));
                        theta_speed = (x[0]*C_V(c,t)-x[1]*C_U(c,t))/(pow(x[0],2)+pow(x[1],2));
                        tangential_v = radial_coordinate*theta_speed;

                        sum_v += tangential_v*C_VOLUME(c,t);
                        sum_vol += C_VOLUME(c,t);
                    }
                    else
                    sum_v += 0;
                    sum_vol += 0;

                end_c_loop_int(c,t)
                }

                sum_ave_v[i] = sum_v/sum_vol;

                #if RP_NODE
                sum_ave_v[i] = PRF_GRSUM1(sum_ave_v[i]);
                #endif /*RP_NODE*/

                #if PARALLEL
                if(I_AM_NODE_ZERO_P)
                #endif /*PARALLEL*/



                printf("Volume Averaged Velocity: %g\n", sum_ave_v[i]);
                fflush(stdout);
                }



            }
        } 

    #endif /*RP_HOST*/
}
