
function [I1, I2] = selectfms(PAM1,PAM2,MaxPAM,dPAM)

% Select planes from set-1
Iok1 = PAM1 < MaxPAM;
Iok2 = PAM2 < MaxPAM;

if dPAM > 1
    I1 = PAM1 < (PAM2 - dPAM);
    I2 = PAM2 < (PAM1 - dPAM);
else
    I1 = (PAM2 ./ PAM1 - 1) > dPAM & (PAM2 ./ PAM1 - 1) > 0;
    I2 = (PAM1 ./ PAM2 - 1) > dPAM & (PAM1 ./ PAM2 - 1) > 0;
end   

I1 = I1 & Iok1;
I2 = I2 & Iok2;


end