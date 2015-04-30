#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Parametros del promotor de A */
#define alphaA 0.5
#define betaA 20
#define KA 20000
#define nA 2

/*Parametros del promotor de S */
#define alphaS 0.5
#define betaS 30
#define KS 1000
#define nS 2

/*Parametros del promotor de I y S juntos */
#define alphaIS 0
#define betaIS_I 2
#define betaIS_S 10
#define betaIS_IS 10
#define KI 20000
#define nI 2

/*Parametros de A*/
#define gammarA 0.2
#define kA 25
#define gammaA 1.0/30.0

/*Parametros de E*/
#define gammarE 0.2
#define kE 20
#define gammaE 1.0/30.0

/*Parametros de I*/
#define gammarI 0.2
#define kI 40
#define gammaI 1.0/30.0

#define h 1.0E-2
#define max_t 400
#define min_t 0

float func_rA_prime(float t, float S, float rA, float A, float rE, float E, float rI, float I);
float func_A_prime(float t, float S, float rA, float A, float rE, float E, float rI, float I);

float func_rE_prime(float t, float S, float rA, float A, float rE, float E, float rI, float I);
float func_E_prime(float t, float S, float rA, float A, float rE, float E, float rI, float I);

float func_rI_prime(float t, float S, float rA, float A, float rE, float E, float rI, float I);
float func_I_prime(float t, float S, float rA, float A, float rE, float E, float rI, float I);

void rk4_step(float t_old, float S, float rA_old, float A_old, float rE_old, float E_old, float rI_old, float I_old, float *t_new, float *rA_new, float *A_new, float *rE_new, float *E_new, float *rI_new, float *I_new);

void rk4_step_1(float t_old, float S, float rA_old, float A_old, float rE_old, float E_old, float rI_old, float I_old, float I_oldA, float *t_new, float *rA_new, float *A_new, float *rE_new, float *E_new, float *rI_new, float *I_new);

void rk4(float t0, float rA0, float A0, float rE0, float E0, float rI0, float I0, char *filename);


float func_rA_prime(float t, float S, float rA, float A, float rE, float E, float rI, float I){
  return alphaIS + betaIS_I/( pow((KI/I),nI) + 1 + pow((S/KS),nS)*pow((KI/I),nI) + pow((S/KS),nS) ) + betaIS_S/( pow((KS/S),nS) + pow((I/KI),nI)*pow((KS/S),nS) + 1 + pow((I/KI),nI) ) + betaIS_IS/( pow((KI/I),nI)*pow((KS/S),nS) + pow((KS/S),nS) + pow((KI/I),nI) + 1 ) - gammarA*rA;
}

float func_A_prime(float t, float S, float rA, float A, float rE, float E, float rI, float I){
  return kA*rA - gammaA*A;
}

float func_rE_prime(float t, float S, float rA, float A, float rE, float E, float rI, float I){
  return alphaS + betaS/(1.0 + pow((KS/S),nS)) - gammarE*rE; 
}

float func_E_prime(float t, float S, float rA, float A, float rE, float E, float rI, float I){
  return kE*rE - gammaE*E;
}

float func_rI_prime(float t, float S, float rA, float A, float rE, float E, float rI, float I){
  return alphaA + betaA/(1.0 + pow((KA/A),nA)) - gammarI*rI;
}

float func_I_prime(float t, float S, float rA, float A, float rE, float E, float rI, float I){
  return kI*rI - gammaI*I;
}

void rk4_step(float t_old, float S, float rA_old, float A_old, float rE_old, float E_old, float rI_old, float I_old, float *t_new, float *rA_new, float *A_new, float *rE_new, float *E_new, float *rI_new, float *I_new){
  

  float krA_1 = func_rA_prime(t_old, S, rA_old, A_old, rE_old, E_old, rI_old, I_old);
  float kA_1 = func_A_prime(t_old, S, rA_old, A_old, rE_old, E_old, rI_old, I_old);
  float krE_1 = func_rE_prime(t_old, S, rA_old, A_old, rE_old, E_old, rI_old, I_old);
  float kE_1 = func_E_prime(t_old, S, rA_old, A_old, rE_old, E_old, rI_old, I_old);
  float krI_1 = func_rI_prime(t_old, S, rA_old, A_old, rE_old, E_old, rI_old, I_old);
  float kI_1 = func_I_prime(t_old, S, rA_old, A_old, rE_old, E_old, rI_old, I_old);

  float t_1 = t_old + h/2.0;
  float rA_1 = rA_old + (h/2.0)*krA_1;
  float A_1 = A_old + (h/2.0)*kA_1;
  float rE_1 = rE_old + (h/2.0)*krE_1;
  float E_1 = E_old + (h/2.0)*kE_1;
  float rI_1 = rI_old + (h/2.0)*krI_1;
  float I_1 = I_old + (h/2.0)*kI_1;

  
  float krA_2 = func_rA_prime(t_1, S, rA_1, A_1, rE_1, E_1, rI_1, I_1);
  float kA_2 = func_A_prime(t_1, S, rA_1, A_1, rE_1, E_1, rI_1, I_1);
  float krE_2 = func_rE_prime(t_1, S, rA_1, A_1, rE_1, E_1, rI_1, I_1);
  float kE_2 = func_E_prime(t_1, S, rA_1, A_1, rE_1, E_1, rI_1, I_1);
  float krI_2 = func_rI_prime(t_1, S, rA_1, A_1, rE_1, E_1, rI_1, I_1);
  float kI_2 = func_I_prime(t_1, S, rA_1, A_1, rE_1, E_1, rI_1, I_1);

  float t_2 = t_old + h/2.0;
  float rA_2 = rA_old + (h/2.0)*krA_2;
  float A_2 = A_old + (h/2.0)*kA_2;
  float rE_2 = rE_old + (h/2.0)*krE_2;
  float E_2 = E_old + (h/2.0)*kE_2;
  float rI_2 = rI_old + (h/2.0)*krI_2;
  float I_2 = I_old + (h/2.0)*kI_2;

  
  float krA_3 = func_rA_prime(t_2, S, rA_2, A_2, rE_2, E_2, rI_2, I_2);
  float kA_3 = func_A_prime(t_2, S, rA_2, A_2, rE_2, E_2, rI_2, I_2);
  float krE_3 = func_rE_prime(t_2, S, rA_2, A_2, rE_2, E_2, rI_2, I_2);
  float kE_3 = func_E_prime(t_2, S, rA_2, A_2, rE_2, E_2, rI_2, I_2);
  float krI_3 = func_rI_prime(t_2, S, rA_2, A_2, rE_2, E_2, rI_2, I_2);
  float kI_3 = func_I_prime(t_2, S, rA_2, A_2, rE_2, E_2, rI_2, I_2);

  float t_3 = t_old + h;
  float rA_3 = rA_old + h*krA_3;
  float A_3 = A_old + h*kA_3;
  float rE_3 = rE_old + h*krE_3;
  float E_3 = E_old + h*kE_3;
  float rI_3 = rI_old + h*krI_3;
  float I_3 = I_old + h*kI_3;


  float krA_4 = func_rA_prime(t_3, S, rA_3, A_3, rE_3, E_3, rI_3, I_3);
  float kA_4 = func_A_prime(t_3, S, rA_3, A_3, rE_3, E_3, rI_3, I_3);
  float krE_4 = func_rE_prime(t_3, S, rA_3, A_3, rE_3, E_3, rI_3, I_3);
  float kE_4 = func_E_prime(t_3, S, rA_3, A_3, rE_3, E_3, rI_3, I_3);
  float krI_4 = func_rI_prime(t_3, S, rA_3, A_3, rE_3, E_3, rI_3, I_3);
  float kI_4 = func_I_prime(t_3, S, rA_3, A_3, rE_3, E_3, rI_3, I_3);

  float krA_T = (1.0/6.0)*(krA_1 + 2.0*krA_2 + 2.0*krA_3 + krA_4);
  float kA_T = (1.0/6.0)*(kA_1 + 2.0*kA_2 + 2.0*kA_3 + kA_4);
  float krE_T = (1.0/6.0)*(krE_1 + 2.0*krE_2 + 2.0*krE_3 + krE_4);
  float kE_T = (1.0/6.0)*(kE_1 + 2.0*kE_2 + 2.0*kE_3 + kE_4);
  float krI_T = (1.0/6.0)*(krI_1 + 2.0*krI_2 + 2.0*krI_3 + krI_4);
  float kI_T = (1.0/6.0)*(kI_1 + 2.0*kI_2 + 2.0*kI_3 + kI_4);

  *t_new = t_old + h;
  *rA_new = rA_old + h*krA_T;
  *A_new = A_old + h*kA_T;
  *rE_new = rE_old + h*krE_T;
  *E_new = E_old + h*kE_T;
  *rI_new = rI_old + h*krI_T;
  *I_new = I_old + h*kI_T;
}

void rk4_step_1(float t_old, float S, float rA_old, float A_old, float rE_old, float E_old, float rI_old, float I_old, float I_oldA, float *t_new, float *rA_new, float *A_new, float *rE_new, float *E_new, float *rI_new, float *I_new){
  
  float krA_1 = func_rA_prime(t_old, S, rA_old, A_old, rE_old, E_old, rI_old, I_oldA);
  float kA_1 = func_A_prime(t_old, S, rA_old, A_old, rE_old, E_old, rI_old, I_oldA);
  float krE_1 = func_rE_prime(t_old, S, rA_old, A_old, rE_old, E_old, rI_old, I_old);
  float kE_1 = func_E_prime(t_old, S, rA_old, A_old, rE_old, E_old, rI_old, I_old);
  float krI_1 = func_rI_prime(t_old, S, rA_old, A_old, rE_old, E_old, rI_old, I_old);
  float kI_1 = func_I_prime(t_old, S, rA_old, A_old, rE_old, E_old, rI_old, I_old);

  float t_1 = t_old + h/2.0;
  float rA_1 = rA_old + (h/2.0)*krA_1;
  float A_1 = A_old + (h/2.0)*kA_1;
  float rE_1 = rE_old + (h/2.0)*krE_1;
  float E_1 = E_old + (h/2.0)*kE_1;
  float rI_1 = rI_old + (h/2.0)*krI_1;
  float I_1 = I_old + (h/2.0)*kI_1;

  
  float krA_2 = func_rA_prime(t_1, S, rA_1, A_1, rE_1, E_1, rI_1, I_1);
  float kA_2 = func_A_prime(t_1, S, rA_1, A_1, rE_1, E_1, rI_1, I_1);
  float krE_2 = func_rE_prime(t_1, S, rA_1, A_1, rE_1, E_1, rI_1, I_1);
  float kE_2 = func_E_prime(t_1, S, rA_1, A_1, rE_1, E_1, rI_1, I_1);
  float krI_2 = func_rI_prime(t_1, S, rA_1, A_1, rE_1, E_1, rI_1, I_1);
  float kI_2 = func_I_prime(t_1, S, rA_1, A_1, rE_1, E_1, rI_1, I_1);

  float t_2 = t_old + h/2.0;
  float rA_2 = rA_old + (h/2.0)*krA_2;
  float A_2 = A_old + (h/2.0)*kA_2;
  float rE_2 = rE_old + (h/2.0)*krE_2;
  float E_2 = E_old + (h/2.0)*kE_2;
  float rI_2 = rI_old + (h/2.0)*krI_2;
  float I_2 = I_old + (h/2.0)*kI_2;

  
  float krA_3 = func_rA_prime(t_2, S, rA_2, A_2, rE_2, E_2, rI_2, I_2);
  float kA_3 = func_A_prime(t_2, S, rA_2, A_2, rE_2, E_2, rI_2, I_2);
  float krE_3 = func_rE_prime(t_2, S, rA_2, A_2, rE_2, E_2, rI_2, I_2);
  float kE_3 = func_E_prime(t_2, S, rA_2, A_2, rE_2, E_2, rI_2, I_2);
  float krI_3 = func_rI_prime(t_2, S, rA_2, A_2, rE_2, E_2, rI_2, I_2);
  float kI_3 = func_I_prime(t_2, S, rA_2, A_2, rE_2, E_2, rI_2, I_2);

  float t_3 = t_old + h;
  float rA_3 = rA_old + h*krA_3;
  float A_3 = A_old + h*kA_3;
  float rE_3 = rE_old + h*krE_3;
  float E_3 = E_old + h*kE_3;
  float rI_3 = rI_old + h*krI_3;
  float I_3 = I_old + h*kI_3;


  float krA_4 = func_rA_prime(t_3, S, rA_3, A_3, rE_3, E_3, rI_3, I_3);
  float kA_4 = func_A_prime(t_3, S, rA_3, A_3, rE_3, E_3, rI_3, I_3);
  float krE_4 = func_rE_prime(t_3, S, rA_3, A_3, rE_3, E_3, rI_3, I_3);
  float kE_4 = func_E_prime(t_3, S, rA_3, A_3, rE_3, E_3, rI_3, I_3);
  float krI_4 = func_rI_prime(t_3, S, rA_3, A_3, rE_3, E_3, rI_3, I_3);
  float kI_4 = func_I_prime(t_3, S, rA_3, A_3, rE_3, E_3, rI_3, I_3);

  float krA_T = (1.0/6.0)*(krA_1 + 2.0*krA_2 + 2.0*krA_3 + krA_4);
  float kA_T = (1.0/6.0)*(kA_1 + 2.0*kA_2 + 2.0*kA_3 + kA_4);
  float krE_T = (1.0/6.0)*(krE_1 + 2.0*krE_2 + 2.0*krE_3 + krE_4);
  float kE_T = (1.0/6.0)*(kE_1 + 2.0*kE_2 + 2.0*kE_3 + kE_4);
  float krI_T = (1.0/6.0)*(krI_1 + 2.0*krI_2 + 2.0*krI_3 + krI_4);
  float kI_T = (1.0/6.0)*(kI_1 + 2.0*kI_2 + 2.0*kI_3 + kI_4);

  *t_new = t_old + h;
  *rA_new = rA_old + h*krA_T;
  *A_new = A_old + h*kA_T;
  *rE_new = rE_old + h*krE_T;
  *E_new = E_old + h*kE_T;
  *rI_new = rI_old + h*krI_T;
  *I_new = I_old + h*kI_T;

}
void rk4(float t0, float rA0, float A0, float rE0, float E0, float rI0, float I0, char *filename){

  float S = 1;

  int n = (int) ((max_t - min_t)/h);

  float *t = malloc(sizeof(float)*n);
  float *rA = malloc(sizeof(float)*n);
  float *A = malloc(sizeof(float)*n);
  float *rE = malloc(sizeof(float)*n);
  float *E = malloc(sizeof(float)*n);
  float *rI = malloc(sizeof(float)*n);
  float *I = malloc(sizeof(float)*n);
  
  float t_new, rA_new, A_new, rE_new, E_new, rI_new, I_new;

  t[0] = t0;
  rA[0] = rA0;
  A[0] = A0;
  rE[0] = rE0;
  E[0] = E0;
  rI[0] = rI0;
  I[0] = I0;
 
  FILE *out = fopen(filename, "w");
  fprintf(out, "%f %f %f %f %f %f %f\n", t0, rA0, A0, rE0, E0, rI0, I0);

  int i;
  int Delta_i = 10000;
  for(i=1;i<n;i++){
    /*
    rk4_step(t[i-1], S, rA[i-1], A[i-1], rE[i-1], E[i-1], rI[i-1], I[i-1], &t_new, &rA_new, &A_new, &rE_new, &E_new, &rI_new, &I_new);
    */
    if(i <= Delta_i){
      rk4_step_1(t[i-1], S, rA[i-1], A[i-1], rE[i-1], E[i-1], rI[i-1], I[i-1],I0, &t_new, &rA_new, &A_new, &rE_new, &E_new, &rI_new, &I_new);
    } else {
      rk4_step_1(t[i-1], S, rA[i-1], A[i-1], rE[i-1], E[i-1], rI[i-1], I[i-1], I[i-Delta_i], &t_new, &rA_new, &A_new, &rE_new, &E_new, &rI_new, &I_new);
    }
    /*
    rk4_step_1(t[i-1], S, rA[i-1], A[i-1], rE[i-1], E[i-1], rI[i-1], I[i-1], I[i-10], &t_new, &rA_new, &A_new, &rE_new, &E_new, &rI_new, &I_new);
    */
    
    t[i] = t_new;
    rA[i] = rA_new;
    A[i] = A_new;
    rE[i] = rE_new;
    E[i] = E_new;
    rI[i] = rI_new;
    I[i] = I_new;

    /*
    if(t_new < 50.0){
      S = 1;
    } else if(t_new >= 50  && t_new <= 350){
      S = 5;
    } else {
      S = 1;
    }
    */
    fprintf(out, "%f %f %f %f %f %f %f\n", t_new, rA_new, A_new, rE_new, E_new, rI_new, I_new);
  }
  fclose(out);
}  

int main (void) {
  
  float t0 = 0.0;
  float rA0 = 0.0;
  float A0 = 0.0;
  float rE0 = 0.0;
  float E0 = 0.0;
  float rI0 = 0.0;
  float I0 = 0.0;
  
  char filename[50] = "data.dat";
  
  rk4(t0, rA0, A0, rE0, E0, rI0, I0, filename);

  return 0;

}
