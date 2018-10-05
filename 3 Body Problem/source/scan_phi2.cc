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

double phi_min, phi_max, d_phi;

void init_search();
void integrate(int n, bool g_d);
void search();
void gradient_descent(int n, double g_step, double g_size);
double distance(Pos p1, Pos p2);

void init_search()
{
  int n = 0;
  for(double phi = phi_min; phi<=phi_max; phi += d_phi, n++)
  {
    angl.push_back(temp);
    n = angl.size()-1;
    angl[n].phi = phi;
    angl[n].d_min = 99999999999;
    angl[n].coll = false;
    integrate(n, false);
  }
}

void integrate(int n, bool g_d)
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

  for(int i=0;(t<tt)&&(i<1000000);i++) {
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



        if (d < d_start && candidates[candidates.size() - 1] != n && !g_d) {
          candidates.push_back(n);
        }
      }

      ps1 = p;
    }
}
/*
void search()
{
  candidates_size = candidates.size();
  int curr = 0;
  while(curr < N_steps * candidates_size && curr < candidates.size())
  {
    int ex = int(curr/candidates_size);
    gradient_descent(candidates[curr], dd_phi/(pow(10, ex)), grid_size);
  }
}
*/


void print()
{
  int n = angl.size();
  for (int i = 0; i<n; i++)
  {
    printf("%lf %lf %lf %i\n", angl[i].phi, angl[i].d_min, angl[i].t_min, angl[i].coll);
  }
}

void print_min()
{
  double d_min = 0;
  //printf("Unesi minimalno rastojanje:\n");
  scanf("%lf", &d_min);
  int n = angl.size();
  for (int i = 0; i<n; i++)
  {
    if (i == 0){
      if (angl[i].d_min <= d_min && angl[i].d_min < angl[i+1].d_min)
        printf("%lf %lf %lf %i\n", angl[i].phi, angl[i].d_min, angl[i].t_min, angl[i].coll);
      continue;
    }
    if (i == n-1){
      if (angl[i].d_min <= d_min && angl[i].d_min < angl[i-1].d_min)
        printf("%lf %lf %lf %i\n", angl[i].phi, angl[i].d_min, angl[i].t_min, angl[i].coll);
      continue;
    }
    if (angl[i].d_min <= d_min  && angl[i].d_min < angl[i-1].d_min && angl[i].d_min < angl[i+1].d_min)
      printf("%lf %lf %lf %i\n", angl[i].phi, angl[i].d_min, angl[i].t_min, angl[i].coll);
  }
}

void print_all(int m)
{
  for (int i = 0; i < angl[m].distance.size(); i++)
    printf("%lf %lf %lf %lf %i\n", angl[m].phi, angl[m].distance[i], angl[m].proximity[i], angl[m].time[i], angl[m].coll);
}

int main() {

  //printf("Unesi pocetni i krajnji ugao i korak\n");
  int mod = 0;
  scanf("%lf %lf %lf %i", &phi_min, &phi_max, &d_phi, &mod);

  init_search(); //Finding candidates for periodical solutions
  //search();
  //print();
  /*if (mod)
    print_min();
  else
    print();*/

  //printf("%li %li %li %lf %lf\n", angl[1].distance.size(), angl[1].time.size(), angl[1].proximity.size(), angl[1].d_min, angl[1].phi);

  print_all(0);
  return 0;
}
