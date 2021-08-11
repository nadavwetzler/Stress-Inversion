
function[StressDir,sig_Inv_vec,misfit1,SM,sig_Inv, sig_Inv_eig_d] = Stress_inv(FMS, frc, quiet, zcoord, rowm, Mag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fault Stress Inversion - Code for General Fault Stress Inversion Solution
% Code is based on the analysis by Reches 1987
% Author: Seth Busetti, 2013
% Modifues by Nadav Wetzler
%
% INPUTS:
% FMS:  (matrix) n x 3, n is number of focal mechanisms listed: strike, dip,
% and rake
% frc: (single value) the friction coefficients that you like to solve (e.g. 0.5)
% quiet: if set to '1' , no messeges are displayed
% zcoord: (vector) size of n, depth (positive) in km
% rowm: is the mean density in km/m^3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear('FMS')
%FMS = [5,   89,  1];
if frc == 0
    frc = 0.01;
end
%disp(['Calculating stress for mue ',num2str(frc)])
numevents = size(FMS,1);
%% Mechanical Input Data
plot_stereo = 0;
C = 0;


%%
if nargin == 2
    quiet = 1;
    Sv(1:numevents) = -100;
end

if nargin == 3
    Sv(1:numevents) = -100;
end

if nargin == 4
    if zcoord == -999
        Sv(1:numevents) = -100;
    else
        
        rowm = 2000; % Density
        Sv = zcoord.*9.8*rowm; % Pressure gradient in Pa/m
        Sv = -Sv/1000;
    end
end

if nargin >= 5
    if rowm == 0
        rowm = 2000;
    end
    if zcoord == -999
        Sv(1:numevents) = -100;
    else
        Sv = zcoord.*9.8*rowm; % Pressure gradient in Pa/m
        Sv = -Sv/1000;
    end
end
%% Solve for B vector
strike = make360(FMS(:,1))';
dip = FMS(:,2);
rake = FMS(:,3);

[Svec,Nvec] = FMS2SLN3(strike,dip,rake);
RT(1:numevents) = 0;
for ii = 1:numevents
    RT(ii) = abs(Svec(ii,1) * Nvec(ii,1) + Svec(ii,2) * Nvec(ii,2) + Svec(ii,3) * Nvec(ii,3));
    if RT(ii) >= 1
        RT(ii) = .99999;
    end
    RT(ii) = pi / 2 - atan(RT(ii) / sqrt(1 - RT(ii) ^ 2));
    RT(ii) = digit1(abs(RT(ii) * 180 / pi),1);
end

%Bvec = Nvec.*Svec;
Bvec(1:numevents,1:3) = 0;
for ii = 1:numevents
    Bvec(ii,1) = Nvec(ii,2)*Svec(ii,3) - Nvec(ii,3)*Svec(ii,2);
    Bvec(ii,2) = Nvec(ii,3)*Svec(ii,1) - Nvec(ii,1)*Svec(ii,3);
    Bvec(ii,3) = Nvec(ii,1)*Svec(ii,2) - Nvec(ii,2)*Svec(ii,1);
end


%% Stress Inversion Matrices for Loaded Faults
% (Using a linearization method after Reches, 1987)
%  Assemble A and F matrices for all faults
A = zeros(2*numevents,5);
Fv = zeros(2*numevents,1);

for i = 1:numevents
    
    
    A(i,1) = -frc*Nvec(i,1)^2 - Nvec(i,1)*Svec(i,1);
    A(i,2) = -frc*Nvec(i,2)^2 - Nvec(i,2)*Svec(i,2);
    A(i,3) = -2*frc*Nvec(i,2)*Nvec(i,3) - (Nvec(i,2)*Svec(i,3) + Svec(i,2)*Nvec(i,3));
    A(i,4) = -2*frc*Nvec(i,1)*Nvec(i,3) - (Nvec(i,3)*Svec(i,1) + Svec(i,3)*Nvec(i,1));
    A(i,5) = -2*frc*Nvec(i,1)*Nvec(i,2) - (Nvec(i,1)*Svec(i,2) + Svec(i,1)*Nvec(i,2));
    
    
    A(i+numevents,1) =  Nvec(i,1)*Bvec(i,1);
    A(i+numevents,2) =  Nvec(i,2)*Bvec(i,2);
    A(i+numevents,3) =  Nvec(i,2)*Bvec(i,3) + Nvec(i,3)*Bvec(i,2);
    A(i+numevents,4) =  Nvec(i,3)*Bvec(i,1) + Nvec(i,1)*Bvec(i,3);
    A(i+numevents,5) =  Nvec(i,1)*Bvec(i,2) + Nvec(i,2)*Bvec(i,1);
    
    
    Fv(i,1) = C+frc*Sv(i);
    Fv(i+numevents,1) = 0;
end


%% Solve for Global Inversion Matrix D (A X D = F) on Faults
if nargin == 6
    W = zeros(2*numevents,1);
    W(1:numevents) = Mag;
    W(numevents+1:end) = Mag;
    X = lscov(A,Fv,W);
else
    X = A\Fv;
end

% X = invZR(A,Fv);
%%

Sv_mean = mean(Sv);

sig_Inv = [X(1)+Sv_mean, X(5), X(4);...
    X(5), X(2)+Sv_mean, X(3);...
    X(4), X(3), Sv_mean];

[sig_Inv_vec,sig_Inv_eig]=eig(sig_Inv);

[~,ind] = sort(diag(sig_Inv_eig));
% ind = flipud(ind);
sig_Inv_vec = sig_Inv_vec(:,ind);
sig_Inv_eig_d = diag(sig_Inv_eig);


%% Calc stress directions
ALPHA = zeros(3,1);
DIP = zeros(3,1);
for ii=1:3
    [ALPHA(ii), DIP(ii)] = Eig2ang(sig_Inv_vec(:,ii));
end

StressDir(:,1) = ALPHA;
StressDir(:,2) = DIP;

%% Calc misfit
misfit1 = CalcMS_FMS_Stress(sig_Inv, FMS,frc);
SM = mean(misfit1);


%SM2 = mean(misfit1);
if quiet == 0
    disp(['SM  ',num2str(round(SM))])
end
%%

SHaz = atand(sig_Inv_vec(2,2)/sig_Inv_vec(1,2));
if SHaz < 0 SHaz = SHaz+180; end

pstress = 1/3*(sig_Inv_eig(1,1)+sig_Inv_eig(2,2)+sig_Inv_eig(3,3));
devstress = 1/2*((sig_Inv_eig(1,1)-sig_Inv_eig(3,3))^2+...
    (sig_Inv_eig(1,1)-sig_Inv_eig(2,2))^2+...
    (sig_Inv_eig(2,2)-sig_Inv_eig(3,3))^2)^.5;


Sv_mean = mean(Sv);


%% Stress Results + plot Planes
if plot_stereo == 1
    fz =14; %font size
    szs = 200; %size of stress dots
    Rr = 0.5; %size of stress circle
    ang = 0:1:360;
    
    dRr = 0.1*Rr;
    Xo = Rr*cosd(ang);
    Yo = Rr*sind(ang);
    % Plot stresses
    figure(2)
    % plot large circle
    plot(Xo,Yo)
    hold on
    plot([0 0],[Rr -Rr]);
    hold on
    plot([Rr -Rr],[0 0]);
    hold on
    ticks = [0 30 60 90];
    for ii=1:length(ticks)
        text(0,ticks(ii)/90*Rr,num2str(ticks(ii)),'fontsize',fz);
    end
    
    % plot internal circles
    for ii=2:(length(ticks)-1)
        Xo = ticks(ii)/90*Rr*cosd(ang);
        Yo = ticks(ii)/90*Rr*sind(ang);
        plot(Xo,Yo)
        hold on
    end
    
    text(Rr+dRr,0,num2str(90),'fontsize',fz);
    text(-(Rr+2*dRr),0,num2str(270),'fontsize',fz);
    text(0,-(Rr+dRr),num2str(180),'fontsize',fz);
    
    for ii=1:numevents
        
        azi = FMS(ii,1); %from dip dir to strike
        dip = FMS(ii,2);
        if (dip == 90)
            dip = 89.9;
        else
        end
        tpd = tand(180*0.5 - dip);
        for yy=1:90
            ang = (yy-1);
            arg = ((cosd(dip))^2*(sind(ang))^2)^0.5/cosd(ang);
            saz(yy) = atand(arg);
            taz = tand(saz(yy))^2;
            arg = (tpd + tpd*taz +taz)^0.5;
            ainp(yy) = acosd(tand(saz(yy))/arg);
        end
        saz(91) = 90;
        ainp(91) = 180*0.5 - dip;
        qq=1;
        for yy=1:2:180
            if (yy < 91)
                mi = yy;
                azz = saz(yy) + azi;
            else
                mi = 181 - yy;
                azz = 180 - saz(mi) + azi;
            end
            radius = Rr*sqrt(2.0)*sind(ainp(mi)*0.5);
            xSol(ii,qq) = radius*sind(azz);
            ySol(ii,qq) = radius*cosd(azz);
            qq=qq+1;
        end
        plot(xSol(ii,:),ySol(ii,:),'k','LineWidth',0.5)
        hold on
    end
    
    Sig_c = 'rgb';
    for ii = 1:3
        % plot Sigma
        scatter(Rr*sind(ALPHA(ii))*(1-DIP(ii)/90),Rr*cosd(ALPHA(ii))*(1-DIP(ii)/90),szs,Sig_c(ii),'filled')
    end
    set(gca,'fontsize',fz,'DataAspectRatio',[1 1 1])
    axis off
end
if quiet == 0
    disp('Fault Stress Parameters Successfully Calculated.')
    
    disp('---------FINISHED------------')
    disp('')
end



end



function [ALPHA, DIP] = Eig2ang(V)

if V(3) < 0
    V = -V;
end

ALPHA = atan(V(2) / V(1));

if V(1) < 0
    ALPHA = ALPHA + pi;
end
if V(2) <0 && V(1) >0
    ALPHA = ALPHA + 2 * pi;
end


if V(3) == 1
    V(3) = .9999;
end

ALPHA = rad2deg(ALPHA);
DIP = rad2deg(atan(V(3) / sqrt(1 - V(3) ^ 2)));


end
