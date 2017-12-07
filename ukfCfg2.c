 /******************************************************************************************************************************************************************************************************\
 *** 
 *** Description       : IMPLEMENTATION OF JOINT UKF for state and parameter estimation: UNDER CONSTRUCTION   
 ***                     This configuration describes vehicle kinematic equations augmented with effective tyre radius as additional unknown state.
 *** Codefile          : ukfCfg2.c
 *** Documentation(UNDER CONSTRUCTION) : https://github.com/ivo-georgiev/ukfLib/wiki/Vehicle-tracking-and-tyre-radius-estimation
 ***
 *** State vector
 *** state 0 = x[0] = x(k)
 *** state 1 = x[1] = y(k)
 *** state 2 = x[2] = theta(k)
 *** state 3 = x[3] = tyre_rad(k)
 ***
 *** tyres type 1: 245/40 R19
 *** tyres type 2: 245/45 R18
 *** Wheelbase: 2.8498 m
 *** Wheel track: 1.6002 m
 *** Steering ratio: ~16:1
 ***
 ***         #####               #####--
 ***           |                   |   ^
 ***           |                   |   |
 ***           |-------------------|   1.6002 m   
 ***           |                   |   |
 ***           |                   |   |
 ***         #####               #####--
 ***           |<-----2.8498 m---->|  
 ***
 *** MIT License
 ***
 *** Copyright (c) 2017 ivo-georgiev
 ***  
 *** Permission is hereby granted, free of charge, to any person obtaining a copy
 *** of this software and associated documentation files (the "Software"), to deal
 *** in the Software without restriction, including without limitation the rights
 *** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *** copies of the Software, and to permit persons to whom the Software is
 *** furnished to do so, subject to the following conditions:
 ***    
 *** The above copyright notice and this permission notice shall be included in all
 *** copies or substantial portions of the Software.
 ***      
 *** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *** OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *** SOFTWARE.
\******************************************************************************************************************************************************************************************************/
#include "ukfCfg.h"

static void Fx1(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT);
static void Fx2(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT);
static void Fx3(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT);
static void Fx4(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT);

static void Hy1(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx);
static void Hy2(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx);

static tPredictFcn PredictFcn[4] = {&Fx1,&Fx2,&Fx3,&Fx4};
static tObservFcn  ObservFcn[2] = {&Hy1,&Hy2};

//define state const for easy matrix indexing
#define xPosIdx (uint8)0
#define yPosIdx (uint8)1
#define thetaIdx (uint8)2
#define tyreRadIdx (uint8)3

//define 
#define angularSpeedIdx (uint8)0
#define steeringAngleIdx (uint8)1
//-----------------------
//UKF Processing matrix
//-----------------------
static float64 Sc_vector[1][3] = {{0.1,2,0}};
static float64 Wm_weight_vector[1][9] = {{0,0,0,0,0,0,0,0,0}};
static float64 Wc_weight_vector[1][9] = {{0,0,0,0,0,0,0,0,0}};
static float64 u_system_input[4][1] = {{0},{0},{0},{0}}; 
static float64 u_prev_system_input[4][1] = {{0},{0},{0},{0}};
static float64 y_meas[2][1] = {{0},{0}};
static float64 y_predicted_mean[2][1] = {{0},{0}};
static float64 x_system_states[4][1] = {{0},{0},{0},{0.3}};
static float64 x_system_states_ic[4][1] = {{0},{0},{0},{0.3}};
static float64 x_system_states_correction[4][1] = {{0},{0},{0},{0}};
static float64 X_sigma_points[4][9]=
{/*  s1  s2  s3  s4  s5  s6  s7  s8  s9        */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x2 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x3 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x4 */
};

//Sigma points Y(k|k-1) = y_m
static float64 Y_sigma_points[2][9]=
{/*  s1  s2  s3  s4  s5  s6  s7  s8  s9        */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* y1 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* y2 */
};

//State covariance  P(k|k-1) = P_m, P(k)= P  
static float64 Pxx_error_covariance[4][4]=
{/*  x1, x2, x3, x4        */
    {0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0}, /* x2 */ 
    {0,  0,  0,  0}, /* x3 */  
    {0,  0,  0,  0}, /* x4 */
};

//State covariance initial values
static float64 Pxx0_init_error_covariance[4][4]=
{/*  x1, x2, x3, x4        */
    {1,  0,  0,  0}, /* x1 */
    {0,  1,  0,  0}, /* x2 */ 
    {0,  0,  1,  0}, /* x3 */  
    {0,  0,  0,  1}, /* x4 */
};

//Process noise covariance Q : initial noise assumptions
static float64 Qxx_process_noise_cov[4][4]=
{/*  x1, x2, x3, x4        */
    {0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0}, /* x2 */ 
    {0,  0,  4,  0}, /* x3 */  
    {0,  0,  0,  4}, /* x4 */
};

//Output noise covariance: initial noise assumptions
static float64 Ryy0_init_out_covariance[2][2]=
{/*  y1, y2         */
    {1,  0},  /* y1 */
    {0,  1},  /* y2 */
};

//Output covariance Pyy = R (initial assumption)
static float64 Pyy_out_covariance[2][2]=
{/*  y1, y2         */
    {0,  0},  /* y1 */
    {0,  0},  /* y2 */
};

static float64 Pyy_out_covariance_copy[2][2]=
{/*  y1, y2         */
    {0,  0},  /* y1 */
    {0,  0},  /* y2 */
};

//cross-covariance of state and output
static float64 Pxy_cross_covariance[4][2]=
{/*  y1, y2         */
    {0,  0},  /* x1 */
    {0,  0},  /* x2 */
    {0,  0},  /* x3 */
    {0,  0},  /* x4 */
};

//Kalman gain matrix
static float64 K_kalman_gain[4][2]=
{  
    {0, 0},
    {0, 0},
    {0, 0},
    {0, 0},
};

static float64 Pxx_covariance_correction[4][4]=
{/*  x1, x2, x3, x4        */
    {0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0}, /* x2 */ 
    {0,  0,  0,  0}, /* x3 */  
    {0,  0,  0,  0}, /* x4 */
};

static float64 I_identity_matrix[2][2]=
{
    {0,  0},
    {0,  0},
};

tUkfMatrix UkfMatrixCfg2 = 
{
    {COLXROW(Sc_vector),NROWS(Sc_vector),NCOL(Sc_vector),&Sc_vector[0][0]},
    {COLXROW(Wm_weight_vector),NROWS(Wm_weight_vector),NCOL(Wm_weight_vector),&Wm_weight_vector[0][0]},
    {COLXROW(Wc_weight_vector),NROWS(Wc_weight_vector),NCOL(Wc_weight_vector),&Wc_weight_vector[0][0]},
    {COLXROW(x_system_states),NROWS(x_system_states),NCOL(x_system_states),&x_system_states[0][0]},
    {COLXROW(x_system_states_ic),NROWS(x_system_states_ic),NCOL(x_system_states_ic),&x_system_states_ic[0][0]},
    {0,0,0,NULL},//{COLXROW(x_system_states_limits),NROWS(x_system_states_limits),NCOL(x_system_states_limits),&x_system_states_limits[0][0]},
    {0,0,0,NULL},//{COLXROW(x_system_states_limits_enable),NROWS(x_system_states_limits_enable),NCOL(x_system_states_limits_enable),&x_system_states_limits_enable[0][0]},
    {COLXROW(x_system_states_correction),NROWS(x_system_states_correction),NCOL(x_system_states_correction),&x_system_states_correction[0][0]},
    {COLXROW(u_system_input),NROWS(u_system_input),NCOL(u_system_input),&u_system_input[0][0]},
    {COLXROW(u_prev_system_input),NROWS(u_prev_system_input),NCOL(u_prev_system_input),&u_prev_system_input[0][0]},
    {COLXROW(X_sigma_points),NROWS(X_sigma_points),NCOL(X_sigma_points),&X_sigma_points[0][0]},
    {COLXROW(Y_sigma_points),NROWS(Y_sigma_points),NCOL(Y_sigma_points),&Y_sigma_points[0][0]},
    {COLXROW(y_predicted_mean),NROWS(y_predicted_mean),NCOL(y_predicted_mean),&y_predicted_mean[0][0]},   
    {COLXROW(y_meas),NROWS(y_meas),NCOL(y_meas),&y_meas[0][0]},
    {COLXROW(Pyy_out_covariance),NROWS(Pyy_out_covariance),NCOL(Pyy_out_covariance),&Pyy_out_covariance[0][0]},
    {COLXROW(Pyy_out_covariance_copy),NROWS(Pyy_out_covariance_copy),NCOL(Pyy_out_covariance_copy),&Pyy_out_covariance_copy[0][0]},
    {COLXROW(Ryy0_init_out_covariance),NROWS(Ryy0_init_out_covariance),NCOL(Ryy0_init_out_covariance),&Ryy0_init_out_covariance[0][0]},
    {COLXROW(Pxy_cross_covariance),NROWS(Pxy_cross_covariance),NCOL(Pxy_cross_covariance),&Pxy_cross_covariance[0][0]},
    {COLXROW(Pxx_error_covariance),NROWS(Pxx_error_covariance),NCOL(Pxx_error_covariance),&Pxx_error_covariance[0][0]},
    {COLXROW(Pxx0_init_error_covariance),NROWS(Pxx0_init_error_covariance),NCOL(Pxx0_init_error_covariance),&Pxx0_init_error_covariance[0][0]},
    {COLXROW(Qxx_process_noise_cov),NROWS(Qxx_process_noise_cov),NCOL(Qxx_process_noise_cov),&Qxx_process_noise_cov[0][0]},
    {COLXROW(K_kalman_gain),NROWS(K_kalman_gain),NCOL(K_kalman_gain),&K_kalman_gain[0][0]},
    {COLXROW(I_identity_matrix),NROWS(I_identity_matrix),NCOL(I_identity_matrix),&I_identity_matrix[0][0]},  
    {COLXROW(Pxx_covariance_correction),NROWS(Pxx_covariance_correction),NCOL(Pxx_covariance_correction),&Pxx_covariance_correction[0][0]},   
    &PredictFcn[0],
    &ObservFcn[0],
    0.02
};
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void Fx1(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx)
 *** 
 ***  DESCRIPTION:
 ***       Calculate predicted state 0 for each sigma point
 ***       X_m[0][sigmaIdx] = f(X_p, u_p) = 
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tMatrix *          pu_p                                 
 ***      tMatrix *          pX_p                                 Pointer to the sigma points array at (k-1) moment 
 ***      tMatrix *          pX_m                                 Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      uint8              sigmaIdx                             Sigma point index.
 ***  RETURNS:
 ***           
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void Fx1(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT)
{
    const uint8 nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
    const float64 xPosPrev = pX_p->val[nCol*xPosIdx + sigmaIdx];
    const float64 angularSpeedPrev = pu_p->val[angularSpeedIdx];
    const float64 tyreRadPrev = pX_p->val[nCol*tyreRadIdx + sigmaIdx];
    const float64 thethaPrev = pX_p->val[nCol*thetaIdx + sigmaIdx];

    pX_m->val[nCol*xPosIdx + sigmaIdx] = xPosPrev + dT * angularSpeedPrev * tyreRadPrev * cos(thethaPrev);
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void Fx2(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx)
 *** 
 ***  DESCRIPTION:
 ***       Calculate predicted state 1 for each sigma point.
 ***       X_m[1][sigmaIdx] = f(X_p, u_p) = 
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tMatrix *          pu_p                                 
 ***      tMatrix *          pX_p                                 Pointer to the sigma points array at (k-1) moment 
 ***      tMatrix *          pX_m                                 Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      uint8              sigmaIdx                             Sigma point index.
 ***  RETURNS:
 ***           
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void Fx2(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT)
{
    const uint8 nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
    const float64 yPosPrev = pX_p->val[nCol*xPosIdx + sigmaIdx];
    const float64 angularSpeedPrev = pu_p->val[angularSpeedIdx];
    const float64 tyreRadPrev = pX_p->val[nCol*tyreRadIdx + sigmaIdx];
    const float64 thethaPrev = pX_p->val[nCol*thetaIdx + sigmaIdx];

    pX_m->val[nCol*yPosIdx + sigmaIdx] = yPosPrev + dT * angularSpeedPrev * tyreRadPrev * sin(thethaPrev);
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void Fx3(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx)
 *** 
 ***  DESCRIPTION:
 ***       Calculate predicted state 2 for each sigma point. 
 ***       X_m[2][sigmaIdx] = f(X_p, u_p) = 
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tMatrix *          pu_p                                 
 ***      tMatrix *          pX_p                                 Pointer to the sigma points array at (k-1) moment 
 ***      tMatrix *          pX_m                                 Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      uint8              sigmaIdx                             Sigma point index.
 ***  RETURNS:
 ***           
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void Fx3(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT)
{
    const uint8 nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
    const float64 xPosPrev = pX_p->val[nCol*xPosIdx + sigmaIdx];
    const float64 angularSpeedPrev = pu_p->val[angularSpeedIdx];
    const float64 steeringAnglePrev = pu_p->val[steeringAngleIdx];
    const float64 tyreRadPrev = pX_p->val[nCol*tyreRadIdx + sigmaIdx];
    const float64 thehtaPrev = pX_p->val[nCol*thetaIdx + sigmaIdx];
    const float64 WheelBase = 2.8498; //{m}
    
    pX_m->val[nCol*thetaIdx + sigmaIdx] = thehtaPrev + dT * angularSpeedPrev * tyreRadPrev * (1/WheelBase) * tan(steeringAnglePrev);

}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void Fx4(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx)
 *** 
 ***  DESCRIPTION:
 ***       Calculate predicted state 3 for each sigma point. 
 ***       X_m[3][sigmaIdx] = f(X_p, u_p) =
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tMatrix *          pu_p                                 
 ***      tMatrix *          pX_p                                 Pointer to the sigma points array at (k-1) moment 
 ***      tMatrix *          pX_m                                 Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      uint8              sigmaIdx                             Sigma point index.
 ***  RETURNS:
 ***           
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void Fx4(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT)
{
    const uint8 nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
    
    pX_m->val[nCol*tyreRadIdx + sigmaIdx] = pX_p->val[nCol*tyreRadIdx + sigmaIdx];

    pu_p = pu_p;
    dT = dT;
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void Hy1(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx)
 *** 
 ***  DESCRIPTION:
 ***       Calculate predicted state 3 for each sigma point. 
 ***       Y_m[0][sigmaIdx] = h1(X_m[0][sigmaIdx], u) = 
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tMatrix *          pu                                   
 ***      tMatrix *          pY_m                                 Pointer to the predicted output at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      tMatrix *          pX_m                                 Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      uint8              sigmaIdx                             Sigma point index.
 ***  RETURNS:
 ***           
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void Hy1(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx)
{
    const uint8 nCol = pX_m->ncol;
    pY_m->val[nCol*xPosIdx + sigmaIdx] = pX_m->val[nCol*xPosIdx + sigmaIdx];

    pu = pu;
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void Hy2(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx)
 *** 
 ***  DESCRIPTION:
 ***       Calculate predicted state 3 for each sigma point.
 ***       Y_m[1][sigmaIdx] = h2(X_m[1][sigmaIdx], u) = 
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tMatrix *          pu                                  
 ***      tMatrix *          pY_m                                 Pointer to the predicted output at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      tMatrix *          pX_m                                 Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      uint8              sigmaIdx                             Sigma point index.
 ***  RETURNS:
 ***           
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void Hy2(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx)
{
    const uint8 nCol = pX_m->ncol;
    pY_m->val[nCol*yPosIdx + sigmaIdx] = pX_m->val[nCol*yPosIdx + sigmaIdx];

    pu = pu;
}
