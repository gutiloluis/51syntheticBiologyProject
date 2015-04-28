#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define alpha_S 0.5
#define beta_S 30
#define K_S 15
#define n_S 2
#define gamma_rE 0.2
#define kE 20
#define gamma_E 1.0/30.0

#define h 1.0E-3
#define max_t 400
#define min_t 0

float func_rE_prime(float t, float rE, float E, float S);
float func_E_prime(float t, float rE, float E, float S);
float rk4_step(float t_old, float rE_old, float E_old, float S, float *t_new, float *rE_new, float *E_new);
float rk4(float t0, float rE0, float E0, char *filename);


/*
float func_B_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);
float func_rE_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);
float func_E_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);
float func_rF_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);
float func_F_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);
float func_rI_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);
float func_I_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);

void rk4_step(float t_old, float rA_old, float A_old, float rE_old, float E_old, float rF_old, float F_old, float rI_old, float I_old, float *t_new, float *rA_new, float *A_new, float *r
*/

float func_rE_prime(float t, float rE, float E, float S){
  return alpha_S + beta_S/(1.0 + pow((K_S/S),n_S)) - gamma_rE*rE; 
}

float func_E_prime(float t, float rE, float E, float S){
  return kE*rE - gamma_E*E;
}

float rk4_step(float t_old, float rE_old, float E_old, float S, float *t_new, float *rE_new, float *E_new){
  
  float krE_1 = func_rE_prime(t_old, rE_old, E_old, S);
  float kE_1 = func_E_prime(t_old, rE_old, E_old, S);

  float t_1 = t_old + h/2.0;
  float rE_1 = rE_old + (h/2.0)*krE_1;
  float E_1 = E_old + (h/2.0)*kE_1;


  float krE_2 = func_rE_prime(t_1, rE_1, E_1, S);
  float kE_2 = func_E_prime(t_1, rE_1, E_1, S);

  float t_2 = t_old + h/2.0;
  float rE_2 = rE_old + (h/2.0)*krE_2;
  float E_2 = E_old + (h/2.0)*kE_2;


  float krE_3 = func_rE_prime(t_2, rE_2, E_2, S);
  float kE_3 = func_E_prime(t_2, rE_2, E_2, S);

  float t_3 = t_old + h;
  float rE_3 = rE_old + h*krE_3;
  float E_3 = E_old + h*kE_3;

  
  float krE_4 = func_rE_prime(t_3, rE_3, E_3, S);
  float kE_4 = func_E_prime(t_3, rE_3, E_3, S);

  float krE_T = (1.0/6.0)*(krE_1 + 2.0*krE_2 + 2.0*krE_3 + krE_4);
  float kE_T = (1.0/6.0)*(kE_1 + 2.0*kE_2 + 2.0*kE_3 + kE_4);

  *t_new = t_old + h;
  *rE_new = rE_old + h*krE_T;
  *E_new = E_old + h*kE_T;
}

float rk4(float t0, float rE0, float E0, char *filename){

  float S;

  int n = (int) ((max_t - min_t)/h);
  
  float *t = malloc(sizeof(float)*n);
  float *rE = malloc(sizeof(float)*n);
  float *E = malloc(sizeof(float)*n);

  float t_new, rE_new, E_new;

  t[0] = t0;
  rE[0] = rE0;
  E[0] = E0;
 
  FILE *out = fopen(filename, "w");
  fprintf(out, "%f %f %f\n", t0, rE0, E0);

  int i;

  for(i=1;i<n;i++){

    rk4_step(t[i-1], rE[i-1], E[i-1], S, &t_new, &rE_new, &E_new);
    t[i] = t_new;

    if(t_new < 50.0){
      S = 1;
    } else if(t_new >= 50  && t_new <= 350){
      S = 5;
    } else {
      S = 1;
    }

    rE[i] = rE_new;
    E[i] = E_new;
    
    fprintf(out, "%f %f %f\n", t_new, rE_new, E_new);   
  }
  fclose(out);
}  

int main (void) {
  
  float t0 = 0.0;
  float rE0 = 0.0;
  float E0 = 0.0;
  
  char filename[50] = "data.dat";
  
  rk4(t0, rE0, E0, filename);

  return 0;

}
