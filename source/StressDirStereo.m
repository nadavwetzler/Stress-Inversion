function [Px,Py] = StressDirStereo(StressDir, Rr)

ALPHA = StressDir(:,1);
DIP = StressDir(:,2);

Px = Rr*sind(ALPHA).*(1-DIP/90);
Py = Rr*cosd(ALPHA).*(1-DIP/90);
