function [slick,Nvec] = FMS2SLN3(Strike,Dip,Rake)
% 
if Rake == 90
    Rake = 89.9999;
end

if Rake == -90
    Rake = -89.9999;
end

[DF,DD,SP,ST] = SDR2DDRPT(Strike,Dip,Rake);


ll = length(DF);
Nvec(1:ll,3) = 0;
slick(1:ll,3) = 0;
for ii=1:ll
    
% /* normal vector to fault plane */
    Nvec(ii,1) = cosd(DD(ii)) * sind(DF(ii)); 
    Nvec(ii,2) = sind(DD(ii)) * sind(DF(ii));
    Nvec(ii,3) = cosd(180-DF(ii));

    

%/* slickenside vector calculation */
    slick(ii,1) = cosd(ST(ii)) * sind(SP(ii)+90);
    slick(ii,2) = sind(ST(ii)) * sind(SP(ii)+90);
    slick(ii,3) = -cosd(SP(ii)+90);
end
end