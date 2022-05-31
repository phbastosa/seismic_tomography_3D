# include <cmath>
# include <ctime>
# include <cstdio>
# include "noble2014.hpp"

int main(int argc, char **argv)
{
    system("clear");
    float total_time;
    time_t t_0, t_f;
    t_0 = time(NULL);

    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int nz = atoi(argv[3]);

    float dx = atof(argv[4]);
    float dy = atof(argv[5]);
    float dz = atof(argv[6]);

    float xsrc = atof(argv[7]);
    float ysrc = atof(argv[8]);
    float zsrc = atof(argv[9]);

    float * velocity = importBinaryFloat(nx*ny*nz,argv[10]);

    float * travelTimes = eikonal3D(velocity,xsrc,ysrc,zsrc,nx,ny,nz,dx,dy,dz);    

    exportBinaryFloat(travelTimes,nx*ny*nz,argv[11]);

    t_f = time(NULL);
    total_time = difftime(t_f, t_0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);    

    return 0;
}
