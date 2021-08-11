function PAM = CalcPAM(FMS,fric, sig_Inv0)




strike = FMS(:,1);
dip = FMS(:,2);
rake = FMS(:,3);
[Svec,Nvec] = FMS2SLN3(strike,dip,rake);

numevents = size(FMS,1);

[sig_Inv_vec,sig_Inv_eig]=eig(sig_Inv0);

[~,ind] = sort(diag(sig_Inv_eig));
% ind = flipud(ind);
sig_Inv_vec = sig_Inv_vec(:,ind);
sig_Inv_eig_d = diag(sig_Inv_eig);
X1 = sig_Inv_eig_d(1);
X2 = sig_Inv_eig_d(2);
X3 = sig_Inv_eig_d(3);
Bvec(1:numevents,1:3) = 0;
for ii = 1:numevents
    Bvec(ii,1) = Nvec(ii,2)*Svec(ii,3) - Nvec(ii,3)*Svec(ii,2);
    Bvec(ii,2) = Nvec(ii,3)*Svec(ii,1) - Nvec(ii,1)*Svec(ii,3);
    Bvec(ii,3) = Nvec(ii,1)*Svec(ii,2) - Nvec(ii,2)*Svec(ii,1);
end
%% Arange sigmas X1<X2<X3 (l. 499)
%X1 = digit1(sig_Inv_eig_d(1),2); X2 = digit1(sig_Inv_eig_d(2),2); X3 = digit1(sig_Inv_eig_d(3),2);
% X1 = sig_Inv_eig_d(1); X2 = sig_Inv_eig_d(2); X3 = sig_Inv_eig_d(3);
% T1 = X1; T2 = X2; T3 = X3;
% 
% if X1 <= X2 && X2 <= X3
% 
% elseif X1 <= X3 && X3 <= X2
%     X1 = T1; X2 = T3; X3 = T2;
% 
% elseif X2 <= X3 && X3 <= X1
%     X1 = T2; X2 = T3; X3 = T1;
% 
% elseif X2 <= X1 && X1 <= X3
%     X1 = T2; X2 = T1; X3 = T3;
% 
% elseif X3 <= X1 && X1 <= X2
%     X1 = T3; X2 = T1; X3 = T2;
% 
% elseif X3 <= X2 && X2 <= X1
%     X1 = T3; X2 = T2; X3 = T1;
% else
% end
% if quiet == 0
%     disp('')
%     disp(['     Friction: ',num2str(frc)])
%     disp('Sigma1 Sigma2 Sigma3')
%     disp([round(X1),round(X2),round(X3)])
% end
% sig_Diag(1:3) = [X1,X2,X3];

%%
SX = sig_Inv0(1,1);
SY = sig_Inv0(2,2);
SZ = sig_Inv0(3,3);
SXY = sig_Inv0(1,2);
SXZ = sig_Inv0(1,3);
SYZ = sig_Inv0(2,3);
%% Directional cosines sig-1
A11 = SX - X1; A12 = SXY; A13 = SXZ;
A21 = SXY; A22 = SY - X1; A23 = SYZ;
A31 = SXZ; A32 = SYZ; A33 = SZ - X1;
Z1 = (A11 * A23 - A21 * A13) / (A12 * A23 - A22 * A13);
Z2 = (A11 - Z1 * A12) / A13;
K1 = 1 / (sqrt(1 + Z1 ^ 2 + Z2 ^ 2));
K2 = -K1 * Z1;
K3 = -K1 * Z2;
K1 = -K1; K2 = -K2; K3 = -K3;
RU = A11 * K1 + A12 * K2 + A13 * K3;
RV = A21 * K1 + A22 * K2 + A23 * K3;
RX = A31 * K1 + A32 * K2 + A33 * K3;
BW(1:3) = [K1, K2, K3];
sig1_dc = [K1, K2, K3];

if K3 < 0
    K1 = -K1; K2 = -K2; K3 = -K3;
end

ALPHA(1) = atan(K2 / K1);

if K1 < 0
    ALPHA(1) = ALPHA(1) + pi;
end
if K2 <0 && K1 >0
    ALPHA(1) = ALPHA(1) + 2 * pi;
end


if K3 == 1
    K3 = .9999;
end

ALPHA(1) = rad2deg(ALPHA(1));
DIP(1) = rad2deg(atan(K3 / sqrt(1 - K3 ^ 2)));
%disp('Plunge/Trend   Magnitude')
%disp([num2str(round(DIP(1))),'   / ',num2str(round(ALPHA(1))),'     ',num2str(round(X1))])
%% Directional cosines sig-2
A11 = SX - X2; A12 = SXY; A13 = SXZ;
A21 = SXY; A22 = SY - X2; A23 = SYZ;
A31 = SXZ; A32 = SYZ; A33 = SZ - X2;
Z1 = (A11 * A33 - A31 * A13) / (A12 * A33 - A32 * A13);
Z2 = (A31 - Z1 * A32) / A33;
if isnan(Z2)
    Z2 = 0.001;
end
L1 = 1 / (sqrt(1 + Z1 ^ 2 + Z2 ^ 2));
L2 = -L1 * Z1;
L3 = -L1 * Z2;
L1 = -L1; L2 = -L2; L3 = -L3;
RU = A11 * L1 + A12 * L2 + A13 * L3;
RV = A21 * L1 + A22 * L2 + A23 * L3;
RX = A31 * L1 + A32 * L2 + A33 * L3;
BW(4:6) = [L1, L2, L3];
sig2_dc = [L1, L2, L3];

if L3 < 0
    L1 = -L1; L2 = -L2; L3 = -L3;
end

ALPHA(2) = atan(L2 / L1);

if L1 < 0
    ALPHA(2) = ALPHA(2) + pi;
end
if L2 <0 && L1 >0
    ALPHA(2) = ALPHA(2) + 2 * pi;
end


if L3 == 1
    L3 = .9999;
end

ALPHA(2) = rad2deg(ALPHA(2));
DIP(2) = rad2deg(atan(L3 / sqrt(1 - L3 ^ 2)));
% disp([num2str(round(DIP(2))),'   / ',num2str(round(ALPHA(2))),'     ',num2str(round(X2))])



%% Directional cosines sig-3
A11 = SX - X3; A12 = SXY; A13 = SXZ;
A21 = SXY; A22 = SY - X3; A23 = SYZ;
A31 = SXZ; A32 = SYZ; A33 = SZ - X3;
Z1 = (A21 * A33 - A31 * A23) / (A22 * A33 - A32 * A23);
Z2 = (A21 - Z1 * A22) / A23;
M1 = 1 / (sqrt(1 + Z1 ^ 2 + Z2 ^ 2));
M2 = -M1 * Z1;
M3 = -M1 * Z2;
% M1 = -M1; M2 = -M2; M3 = -M3;
RU = A11 * M1 + A12 * M2 + A13 * M3;
RV = A21 * M1 + A22 * M2 + A23 * M3;
RX = A31 * M1 + A32 * M2 + A33 * M3;
BW(7:9) = [M1, M2, M3];
sig3_dc = [M1, M2, M3];

% this section might be problematic...

if M3 < 0
    M1 = -M1; M2 = -M2; M3 = -M3;
end

ALPHA(3) = atan(M2 / M1);

if M1 < 0
    ALPHA(3) = ALPHA(3) + pi;
end
if M2 <0 && M1 >0
    ALPHA(3) = ALPHA(3) + 2 * pi;
end


if M3 == 1
    M3 = .9999;
end

ALPHA(3) = rad2deg(ALPHA(3));
DIP(3) = rad2deg(atan(M3 / sqrt(1 - M3 ^ 2)));
% disp([num2str(round(DIP(3))),'   / ',num2str(round(ALPHA(3))),'     ',num2str(round(X3))])

StressDir(:,1) = ALPHA;
StressDir(:,2) = DIP;
%% stress ratio

PHI = (X2 - X3) / (X1 - X3);
% disp(['Stress-ratio  ',num2str(PHI)])

%% friction deviaion
NORM = 0;
COTOT =0;
FRIC(1: numevents) = 0;
B1(1: numevents) = 0;
B2(1: numevents) = 0;
B3(1: numevents) = 0;
for ii = 1: numevents
    A1 = Nvec(ii,1) ^ 2;
    B1(ii) = Nvec(ii,1) * Svec(ii,1);
    A2 = Nvec(ii,2) ^ 2;
    B2(ii) = Nvec(ii,2) * Svec(ii,2);
    A3 = Nvec(ii,3) ^ 2;
    B3(ii) = Nvec(ii,3) * Svec(ii,3);
    A4 = 2 * Nvec(ii,2) * Nvec(ii,3);
    B4 = (Nvec(ii,2) * Svec(ii,3) + Svec(ii,2) * Nvec(ii,3));
    A5 = 2 * Nvec(ii,1) * Nvec(ii,3);
    B5 = (Nvec(ii,3) * Svec(ii,1) + Svec(ii,3) * Nvec(ii,1));
    A6 = 2 * Nvec(ii,1) * Nvec(ii,2);
    B6 = Nvec(ii,1) * Svec(ii,2) + Svec(ii,1) * Nvec(ii,2);
    SN = A1 * SX + A2 * SY + A3 * SZ + A4 * SYZ + A5 * SXZ + A6 * SXY;
    SSH = -(B1(ii) * SX + B2(ii) * SY + B3(ii) * SZ + B4 * SYZ + B5 * SXZ + B6 * SXY);
    FRIC(ii) = SSH / SN; % Shear stress/normal stress in slip direction
    COH = SSH - SN * fric;
    COTOT = COTOT + COH;
    NORM = NORM + SN * fric;
end
MEANCO = COTOT / numevents;
NORM = NORM / numevents;
dfric = 0;
sumfric = 0;
for ii=1:numevents
    dfric = dfric + (fric - FRIC(ii));
    sumfric = sumfric + FRIC(ii);
end
fdfric = sqrt(dfric)/numevents;

%% Calculate mean angle between GENERAL and IDEAL principle stress axes: PAM


MIS(1:18) = 0;
MISS(1:18) = 0;
ROTS0(1:numevents) = 0;
MISANG = 90; MISSLP = 180;
NUMHF = 0; MISFIT = 0;
PAM(1:numevents) = 0;
for i = 1:numevents
    RX = (Nvec(i,1) * SX  + Nvec(i,2) * SXY + Nvec(i,3) * SXZ);
    RY = (Nvec(i,1) * SXY + Nvec(i,2) * SY  + Nvec(i,3) * SYZ);
    RZ = (Nvec(i,1) * SXZ + Nvec(i,2) * SYZ + Nvec(i,3) * SZ);
    SIGN = Nvec(i,1) * RX + Nvec(i,2) * RY  + Nvec(i,3) * RZ;
    RP = sqrt(RX ^ 2 + RY ^ 2 + RZ ^ 2);
    if abs(SIGN) > RP; SIGN = -.9999 * RP; end
    SIGS = sqrt(RP ^ 2 - SIGN ^ 2);
    FCMAX = abs(SIGS / SIGN); % Maximum shear stress/normal stress on this fault
    RX = RX / RP; RY = RY / RP; RZ = RZ / RP;
    RT = abs(Svec(i,1) * Nvec(i,1) + Svec(i,2) * Nvec(i,2) + Svec(i,3) * Nvec(i,3));
    if RT >= 1; RT = .99999; end
    RT = pi / 2 - atan(RT / sqrt(1 - RT ^ 2));
    RT = abs(RT * 180 / pi);
    ALP = atan(fric);
    ALP = pi / 4 - ALP / 2;
    ALP1 = ALP * RT / 90;
    ALP2 = RT * pi / 180 - ALP1;
    G(1) = -Nvec(i,1) * cos(ALP2) + Svec(i,1) * cos(ALP1);
    G(2) = -Nvec(i,2) * cos(ALP2) + Svec(i,2) * cos(ALP1);
    G(3) = -Nvec(i,3) * cos(ALP2) + Svec(i,3) * cos(ALP1);
    GTT = sqrt(G(1) ^ 2 + G(2) ^ 2 + G(3) ^ 2);
    G(1) = G(1) / GTT; G(2) = G(2) / GTT; G(3) = G(3) / GTT;
    if abs(G(1)) > abs(G(2)) && abs(G(1)) > abs(G(3)); GMARK = 1; end
    if abs(G(2)) > abs(G(1)) && abs(G(2)) > abs(G(3)); GMARK = 2; end
    if abs(G(3)) > abs(G(1)) && abs(G(3)) > abs(G(2)); GMARK = 3; end

    for k = 1:3
        % MARK = sign(BW(k));
        if abs(BW(k)) > 0.579
            if GMARK == 1
                if abs(BW(1)) > abs(BW(2)) || abs(BW(1)) > abs(BW(3))
                    if sign(G(GMARK)) == sign(BW(GMARK))
                        ROT = abs(BW(1) * G(1) + BW(2) * G(2) + BW(3) * G(3));
                    end
                end
            elseif GMARK == 2
                if abs(BW(2)) > abs(BW(1)) || abs(BW(2)) > abs(BW(3))
                    if sign(G(GMARK)) == sign(BW(GMARK))
                        ROT = abs(BW(1) * G(1) + BW(2) * G(2) + BW(3) * G(3));
                    end
                end
            else
                if abs(BW(3)) > abs(BW(2)) || abs(BW(3)) > abs(BW(1))
                    if sign(G(GMARK)) == sign(BW(GMARK))
                        ROT = abs(BW(1) * G(1) + BW(2) * G(2) + BW(3) * G(3));
                    end
                end

            end


        else
            G = -G;
            ROT = abs(BW(1) * G(1) + BW(2) * G(2) + BW(3) * G(3));
        end
    end



    ROT = rad2deg(pi / 2 - atan(ROT / sqrt(1 - ROT ^ 2)));
    G31 = G(2) * Bvec(i,3) - G(3) * Bvec(i,2);
    G32 = G(3) * Bvec(i,1) - G(1) * Bvec(i,3);
    G33 = G(1) * Bvec(i,2) - G(2) * Bvec(i,1);
    G21 = Bvec(i,1);
    G22 = Bvec(i,2);
    G23 = Bvec(i,3);

    ROT2 = abs(BW(4) * G21 + BW(5) * G22 + BW(6) * G23);
    if ROT2 >= 1; ROT2 = .99999; end
    ROT2 = pi / 2 - atan(ROT2 / sqrt(1 - ROT2 ^ 2));
    ROT2 = abs(ROT2 * 180 / pi);

    ROT3 = abs(BW(7) * G31 + BW(8) * G32 + BW(9) * G33);
    if ROT3 >= 1; ROT3 = .99999; end
    ROT3 = pi / 2 - atan(ROT3 / sqrt(1 - ROT3 ^ 2));
    ROT3 = abs(ROT3 * 180 / pi);
    ROTPA = (ROT * (1 - PHI) + ROT3 * PHI);
    % IMIS = int64((5 + ROTPA) / 5);
    % MIS(IMIS) = MIS(IMIS) + 1;
    if ROTPA > MISANG; NUMHF = NUMHF + 1; end
    MISFIT = MISFIT + ROTPA;
    PAM(i) = ROTPA;
    %     B1(i) = Nvec(i,2) * RZ - Nvec(i,3) * RY;
    %     B2(i) = Nvec(i,3) * RX - Nvec(i,1) * RZ;
    %     B3(i) = Nvec(i,1) * RY - Nvec(i,2) * RX;
    %     RB = sqrt(B1(i) ^ 2 + B2(i) ^ 2 + B3(i) ^ 2);
    %     B1(i) = B1(i) / RB; B2(i) = B2(i) / RB; B3(i) = B3(i) / RB;
    %     Q1 = Nvec(i,2) * B3(i) - Nvec(i,3) * B2(i);
    %     Q2 = Nvec(i,3) * B1(i) - Nvec(i,1) * B3(i);
    %     Q3 = Nvec(i,1) * B2(i) - Nvec(i,2) * B1(i);
    %     RQ = sqrt(Q1 ^ 2 + Q2 ^ 2 + Q3 ^ 2);
    %     Q1 = -Q1 / RQ; Q2 = -Q2 / RQ; Q3 = -Q3 / RQ;
    %     ROTS = (Q1 * Svec(i,1) + Q2 * Svec(i,2) + Q3 * Svec(i,3));
    %     if ROTS >= 1; ROTS = .999; end
    %     ROTS = pi / 2 - atan(ROTS / sqrt(1 - ROTS ^ 2));
    %     ROTS = abs(ROTS * 180 / pi);
    %     IMISS = int64((5 + abs(ROTS)) / 5);
    %     MISS(IMISS) = MISS(IMISS) + 1;
    %     if ROTS > MISSLP; NUMSLIP = NUMSLIP + 1; end
    %     ROTS0(i) = ROTS;
end
PAM = PAM';
end