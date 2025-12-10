%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEVELOPMENT OF A DYNAMIC PHYSIOLOGICALLY-BASED MECHANISTIC KIDNEY MODEL TO PREDICT RENAL CLEARANCE
%% Authors: Weize Huang, Nina Isoherranen
%% Department of Pharmaceutics, School of Pharmacy, University of Washington, Seattle, WA, 98195
%% April 4th 2018


%%
clear all; close all; clc; 
%% 

Tfinish=240;
ivbolusDose1=	0;
IR = 1;
Tstart_infusion = 0;
Infusion_length = 240;
ivT1=	0;
ivT2=	6;
ivT3=	12;
poDose1= 10;
poDose2= 0;
poDose3= 0;
poT1=0;
poT2=0;
poT3=0;
Weight =74; 
BP0=	1;


%% Kidney Mechanistic Model
Q_GFR =7.2;

V_PT= 0.0305; % L
V_HT= 0.0027;
V_DT= 0.0194;
V_CDT= 0.237;   

V_PC= 0.0305; % L
V_HC= 0.0027;
V_DC= 0.0194;
V_CDC= 0.237;  

V_PB= 0.0305; % L
V_HB= 0.0027;
V_DB= 0.0194;
V_CDB= 0.237;  

Q_connectingtubule=11;
Q_initialCD=9;
Q_corticalCD=7;
Q_medullaryCD=5;
Q_papillary=3;
Q_urine=1;

CL_CT = 0; 
CL_BC = 0; 
CL_api_reabs = 0;
CL_bsl_reabs = 0;
CL_kidney_intrinsic = 0;

S_PT = 6107; % dm^2
S_HT = 61;
S_DT = 156;
S_CD = 6.7;

Qc0=	330; %L/hr
Q_kidney0=  60; 
Q_kidney1=  60;
V_ven0=	42; % L
CLh = 0; %L/hr

pH_1 = 7.2; % this is urine pH
pH_2 = 7.1; % this is urine pH
pH_3 = 7.0; % this is urine pH
pH_4 = 7.0; % this is urine pH
pH_5 = 7.0; % this is urine pH
pH_6 = 6.9; % this is urine pH
pH_7 = 6.8; % this is urine pH
pH_8 = 6.7; % this is urine pH
pH_9 = 6.6; % this is urine pH
pH_10 = 6.5; % this is urine pH
pH_11 = 6.5; % this is urine pH

%% For test compound: Neutral drug
alpha_1 = 1; 
alpha_2 = 1;
alpha_3 = 1;
alpha_4 = 1;
alpha_5 = 1;
alpha_6 = 1;
alpha_7 = 1;
alpha_8 = 1;
alpha_9 = 1;
alpha_10 = 1;
alpha_11 = 1;

beta = 1; 
gamma = 1;    
    
fu0=1; 
fu_k=fu0; 

P_drug = 10; 
Papp = 0.00036*P_drug/1.5; 


CL_PD_PX_apical = S_PT * Papp;
CL_PD_PX_basolateral = S_PT/30 * Papp;

CL_PD_LH =  S_HT * Papp;
CL_PD_DS =  S_DT * Papp;
CL_PD_CD =  S_CD * Papp;

sim('Huang_Isoherranen_Mechanistic_Kidney_Model');
    
CL_urine= IR/plasma(end)*1000/60
