#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>

#include "integrator.hh"

//#define phi_min 0
//#define phi_max 1.6
//#define d_phi  0.00000001
#define dd_phi 0.000000001
#define grid_size 10
#define d_start 0.0000001
#define N_steps 5
#define Times 100000

using namespace std;

struct angles {
  double phi;
  vector <double> proximity;
  vector <double> distance;
  vector <double> time;
  double d_min, t_min;
  bool coll;
};

vector <int> candidates;
vector <angles> angl;
angles temp;

double T = 100, tt = 20;
int candidates_size = 0;

double phi, d_phi, d_min, t_min;
int ngd, ns;

int new_ = 0;

void integrate(int n)
{
  double phi = angl[n].phi, _p;
  Pos p, p_init, ps1;

  p.init1(0, 0);
  _p = sqrt(-1 * p.penergy() / 3);

  double et=1e-12;
  double dt0=0.0001;
  double t=0;
  double se=0;
  double od=0;
  double md=1e10;

  bool da=false;

  p.init1(_p * cos(phi), _p * sin(phi));
  p_init.init1(_p * cos(phi), _p * sin(phi));

  for(int i=0;(t<t_min)&&(i<1000000);i++) {
      double dt=dt0;
      double e;

      p=p.integrate(dt,dt0,e,et);

      if (p.coll(0.001)){
        angl[n].coll = true;
        return;
      }

      if(dt>0) {
        t+=dt;
        se+=e;

        double a=(p_init-ps1)*(p-ps1)/((p-ps1)*(p-ps1));
        double d=(p-p_init).abs();
        if((a>0)&&(a<1.)) d=(ps1+(p-ps1)*a-p_init).abs();

        if(!da) {
          if(od>d) da=true;
          od=d;
        }
        else {
          if(d<angl[n].d_min) {
            angl[n].d_min = d;
            angl[n].t_min = t;
          }

          angl[n].time.push_back(t);
          angl[n].distance.push_back(d);
          angl[n].proximity.push_back(-log(d));
        }
      }

      ps1 = p;
    }
}

void gradient_descent()
{
  double phi_min = phi - d_phi*ns/2, phi_max = phi + d_phi*ns/2;
  int n = 0;
  for (double phi = phi_min; phi <= phi_max; phi += d_phi)
  {
    angl.push_back(temp);
    int m = angl.size()-1;
    angl[m].phi = phi;
    angl[m].d_min = 1e10;
    integrate(m);
    if (angl[m].d_min <= d_min && !angl[m].coll)
    {
      new_ = m;
      d_min = angl[m].d_min;
      //printf("new = %.10lf min = %.10lf\n\n\n", angl[m].phi, angl[m].d_min);
    }
  }
}

void print_all()
{
  for(int i = 0; i<angl.size(); i++)
  {
    printf("%.15lf %.10lf %lf %i\n", angl[i].phi, angl[i].d_min, angl[i].t_min, angl[i].coll);
  }
}

void print_max()
{
  printf("%lf %lf %lf %i", phi, d_min, t_min, 0);
}

int main(){
  scanf("%lf %lf %lf %lf %i %i", &phi, &d_phi, &d_min, &t_min, &ngd, &ns);

  for (int i = 0; i<ngd; i++){
    gradient_descent();
    d_phi /= 10;
    phi = angl[new_].phi;
    d_min = angl[new_].d_min;
    t_min = angl[new_].t_min;
  }

  print_all();

  printf("\n\n%.15lf %.15lf %.10lf\n", phi, d_min, t_min);

  return 0;
}
