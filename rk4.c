#include <stdio.h>

float func_rA_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);
float func_A_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);

float func_B_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);

float func_rE_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);
float func_E_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);
float func_rF_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);
float func_F_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);
float func_rI_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);
float func_I_prime(float t, float rA, float A, float rE, float E, float rF, float F, float rI, float I, float S);

void rk4_step(float t_old, float rA_old, float A_old, float rE_old, float E_old, float rF_old, float F_old, float rI_old, float I_old, float *t_new, float *rA_new, float *A_new, float *r

int main () {
  
  runge_kutta_4(t0, rE0, 
