#include "udf.h"
#include "metric.h"

DEFINE_ON_DEMAND(on_demand_Vcalc)
{
    #if !RP_HOST
    Domain *d;
    Thread *tf;

    /* Integrate Velocity Magnitude */
    real sum_surf_vel;
    real sum_surf;
    real avg_vel;
    real x[ND_ND];
    real A[3];                 // Area vector
    real radius, delta_r;
    face_t f;
    int fcount;
    
    d = Get_Domain(1);
    
    fcount = 0;
    sum_surf_vel = 0.0;
    sum_surf = 0.0;
    radius = 0.07525; // half iso-surface between inner and outer cylinders
    delta_r = 4.0e-5;
            
    thread_loop_f(tf, d)
    {
        if (FLUID_THREAD_P(tf))
        {
            begin_f_loop(f,tf)
            {
                if PRINCIPAL_FACE_P(f, tf)
                {
                    F_CENTROID(x,f,tf);
                    F_AREA(A,f,tf);
                    if (sqrt(pow(x[0],2) + pow(x[1],2))>=radius && sqrt(pow(x[0],2) + pow(x[1],2))< radius + delta_r)
                    {
                        sum_surf_vel += sqrt(F_U(f,tf)*F_U(f,tf)+F_V(f,tf)*F_V(f,tf)+F_W(f,tf)*F_W(f,tf))*NV_MAG(A);
                        sum_surf += NV_MAG(A);
                        fcount++;
                    }
                    else
                    {
                        sum_surf_vel += 0;
                        sum_surf += 0;
                    }
                }
            }
            end_f_loop(f,tf)
            avg_vel = sum_surf_vel/sum_surf;
        }

        #if RP_NODE
        avg_vel = PRF_GRSUM1(avg_vel);
        #endif
    }

    #if PARALLEL
    if(I_AM_NODE_ZERO_P)
    #endif
    printf("== F Counts == Radius == Surface Area-Weighted Average Velocity\n");
    printf("%d\t%f\t%g\n", fcount, radius, avg_vel);
    fflush(stdout);

    #endif
}