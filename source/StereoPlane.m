function [xSol,ySol] = StereoPlane(strike,dip,Rr) 
% Rr = 0.5;
n1 = length(strike);

for ii=1:n1

    azi = strike(ii); %from dip dir to strike
    dip1 = dip(ii);
    if (dip1 == 90)
        dip1 = 89.9;
    else
    end
    tpd = tand(180*0.5 - dip1);
    for yy=1:90
        ang = (yy-1);
        arg = ((cosd(dip1))^2*(sind(ang))^2)^0.5/cosd(ang);
        saz(yy) = atand(arg);
        taz = tand(saz(yy))^2;
        arg = (tpd + tpd*taz +taz)^0.5;
        ainp(yy) = acosd(tand(saz(yy))/arg);
    end
    saz(91) = 90;
    ainp(91) = 180*0.5 - dip1;
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
%     plot(xSol(ii,:),ySol(ii,:),'-k','LineWidth',1.5)
%     axis off
%     set(gca,'fontsize',fz,'DataAspectRatio',[1 1 1])

end