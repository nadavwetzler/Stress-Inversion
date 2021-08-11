function [tblFinal,FMSDSA,fric_fms,m1,StressDir_init0S,sig123_err,plung123,strike123,SratioM,SratioM2,S1S2M,S123] = friction_stress(frics,FMSDSA,MaxPAM,dPAM,nB2,Fid1,Fid2)
frics = round(frics,1);

disp('Calc stress all FMS')
disp('----------------------------')
disp(['Calculating initial stress for mue ',num2str(frics(1))])


Rr = 0.5;
SB=numSubplots(length(frics));
[n1,m1] = size(FMSDSA);
m1 = m1 + 1;

FMSDSA(:,m1)= 0;
M1 = FMSDSA(:,8) - min(FMSDSA(:,8)) + 1; M1 = M1/max(M1);
M2 = M1;
FMS1 = FMSDSA(:,12:14);
FMS2 = FMSDSA(:,15:17);
[StressDir_init0,vectors_init0,~,~,sig_Inv0, ~] = Stress_inv([FMS1;FMS2],frics(1),1,-999,0,[M1;M2]);

if Fid2 > 0
    [MSIG_init,Pxb1_init,Pyb1_init] = BSstress([FMS1;FMS2], frics(1), nB2, Rr, [M1,M2]);
    figure(Fid2)
    plotStressInvResults([FMS1;FMS2], StressDir_init0,Pxb1_init,Pyb1_init,MSIG_init,'ALL PLANES - initial','k')
end
m2 = m1;
StressDir_init0S = zeros(3,2,length(frics));
selected_id = zeros(n1,1);
tblFinal = [];
SM = zeros(length(frics),1);
S123 = zeros(length(frics),3);
for ii=1:length(frics)
    [StressDir_init0,vectors_init0,~,~,sig_Inv0, ~] = Stress_inv([FMS1;FMS2],frics(1),1,-999,0,[M1;M2]);
    
    disp(['Calculating iterativly stress for mue ',num2str(frics(ii))])
    m2 = m2 + 1;
    
    % Calc PAMs
    PAM1 = CalcPAM(FMS1,frics(ii),sig_Inv0);
    PAM2 = CalcPAM(FMS2,frics(ii),sig_Inv0);
    
    % Select planes by misfit angle (SM (1) or PAM (2)
    [I1, I2] = selectfms(PAM1, PAM2,MaxPAM,dPAM);
    
    selected_id(I1) = 1;
    selected_id(I2) = 2;
    
    FMSDSA(I1,m1) = FMSDSA(I1,m1) + 1;
    FMSDSA(I2,m1) = FMSDSA(I2,m1) + 1;
    FMS3 = [FMS1(I1,:);FMS2(I2,:)];
    M3 = [FMSDSA(I1,8);FMSDSA(I2,8)];
    FMSDSA(:,m2) = selected_id;
    
    if size(FMS3,1) >0
        title0 = [num2str(size(FMS3,1)),'/',num2str(size(FMSDSA,1)),'        \mu = ',num2str(frics(ii))] ;
        
        [StressDir_init0,vectors_init0,~,SM(ii),sig_Inv0,sig_Inv_eig_d] = Stress_inv(FMS3,frics(ii),1,-999,0,M3); % Z3
        S123(ii,:) = sig_Inv_eig_d;
        StressDir_init0S(:,:,ii) = StressDir_init0;
        
        [MSIG_init,Pxb1_init,Pyb1_init] = BSstress(FMS3, frics(ii), nB2, Rr, M3);
        figure(Fid1)
        subplot(SB(1),SB(2),ii)
        plotStressInvResults(FMS3, StressDir_init0,Pxb1_init,Pyb1_init,MSIG_init,title0,'k')
        
        [tblFinal,Sratio,Sratio2,S1S2] = results2tbl(ii,tblFinal,StressDir_init0,sig_Inv_eig_d,sig_Inv0,frics(ii),MSIG_init,SM(ii));
        SratioM(ii) = Sratio;
        SratioM2(ii) = Sratio2;
        S1S2M(ii) = S1S2;
        for jj=1:3
            sig123_err(ii,jj) =  MSIG_init(jj);
            plung123(ii,jj) = StressDir_init0(jj,2);
            strike123(ii,jj) = StressDir_init0(jj,1);
        end
    else
        sig123_err(ii) = 0;
        plung123(ii) = 0;
        strike123(ii) = 0;
        SratioM(ii) = 0;
        SratioM2(ii) = 0;
        S1S2M(ii) = 0;
    end
end


fric_fms = zeros(n1,1);

end

function [tblFinal,Sratio,Sratio2,S1S2] = results2tbl(tt,tblFinal,StressDir_init0,sig_Inv_eig_d,sig_Inv0,frics,MSIG_init,SM)
% calc Shmax
theta = atand(2 * sig_Inv0(1,2) / (sig_Inv0(1,1) - sig_Inv0(2,2))) / 2;
Shmax =  (sig_Inv0(1,1) + sig_Inv0(2,2)) / 2 - sqrt(0.25 * (sig_Inv0(1,1) - sig_Inv0(2,2)) ^ 2 + sig_Inv0(1,2) ^ 2);
Shmin =  (sig_Inv0(1,1) + sig_Inv0(2,2)) / 2 + sqrt(0.25 * (sig_Inv0(1,1) - sig_Inv0(2,2)) ^ 2 + sig_Inv0(1,2) ^ 2);
Sratio = (sig_Inv_eig_d(2) - sig_Inv_eig_d(3)) / (sig_Inv_eig_d(1) - sig_Inv_eig_d(3));
S1S2 = sig_Inv_eig_d(1) / sig_Inv_eig_d(2);
Sratio2 = (sig_Inv_eig_d(2) - sig_Inv_eig_d(1)) / (sig_Inv_eig_d(3) - sig_Inv_eig_d(1));
tblFinal(tt).id = tt;
tblFinal(tt).mue = frics;
tblFinal(tt).s1a = round(StressDir_init0(1,1));
tblFinal(tt).s1b = round(StressDir_init0(1,2));
tblFinal(tt).s1r = round(MSIG_init(1));
tblFinal(tt).s2a = round(StressDir_init0(2,1));
tblFinal(tt).s2b = round(StressDir_init0(2,2));
tblFinal(tt).s2r = round(MSIG_init(2));
tblFinal(tt).s3a = round(StressDir_init0(3,1));
tblFinal(tt).s3b = round(StressDir_init0(3,2));
tblFinal(tt).s3r = round(MSIG_init(3));
sig1_TXT = [num2str(round(StressDir_init0(1,1))),' / ',num2str(round(StressDir_init0(1,2))),'  ± ',num2str(round(MSIG_init(1)))];
sig2_TXT = [num2str(round(StressDir_init0(2,1))),' / ',num2str(round(StressDir_init0(2,2))),'  ± ',num2str(round(MSIG_init(2)))];
sig3_TXT = [num2str(round(StressDir_init0(3,1))),' / ',num2str(round(StressDir_init0(3,2))),'  ± ',num2str(round(MSIG_init(3)))];
tblFinal(tt).s1s = sig1_TXT;
tblFinal(tt).s2s = sig2_TXT;
tblFinal(tt).s3s = sig3_TXT;
tblFinal(tt).Sratio = Sratio;
tblFinal(tt).Shmax = Shmax;
tblFinal(tt).Shmin = Shmin;
tblFinal(tt).theta = theta;
tblFinal(tt).sig_Inv0 = sig_Inv0;
tblFinal(tt).SM = SM;
end
