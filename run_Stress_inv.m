clc
clear
close all

%% -------- Inversion Parameters ---------------------------------------
nB = 500;
dPAM = 0.1; % Controls the seperation between the stress vectors
MaxPAM = 30; % Maximum PAM angle
frics = 0.1:0.1:0.8; % Friction level range

%% ---------------------- Load FMS data ---------------------------
load('data/FMSDSfinal.mat')

%% ----------- Calc initial regional stress --------------------------------
[tblFinal,FMSDSA,fric_fms,m1,StressDir_init0S,sig123_err,plung123,strike123,Sratio,Sratio2,S1S2,S123] = friction_stress(frics,FMSDS,MaxPAM,dPAM,nB,600,1400);