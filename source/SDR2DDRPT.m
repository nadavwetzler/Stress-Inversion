function [ Dip1, DipDir , Plunge, Trend] = SDR2DDRPT( Strike,Dip,Rake )
% The function converts convensional focal mechanism solution
% in the format of : Strike Dip Rake to 
% Dip, Dip direction, Plunge, and Trend
% Nadav Wetzler
 
l(1) = length(Strike);
l(2) = length(Dip);
l(3) = length(Rake);
checklengths = unique(l);
 
 
if (length(checklengths) > 1)
    error('Check vector lengths')
else
    RakeR(1:l(1)) = 0;
    Plunge(1:l(1)) = 0;
    Trend(1:l(1)) = 0;
    DipDir(1:l(1)) = 0;
    Dip1(1:l(1)) = 0;
    for ii=1:l(1)
        RakeR(ii) = 360-Rake(ii);
         
        Plunge(ii) = digit1(asind(sind(Dip(ii))*sind(RakeR(ii))),0);
         
        Trend(ii) = Strike(ii) + atand(tand(RakeR(ii))*cosd(Dip(ii)));
        %Trend(ii)  = digit1(make360(Strike(ii) + acosd(cosd(Rake(ii))/cosd(Plunge(ii)))),0);
         
        if (abs(Rake(ii)) > 90)
            Trend(ii) = Strike(ii) +180+ atand(tand(RakeR(ii))*cosd(Dip(ii)));
        else
            Trend(ii) = Strike(ii) + atand(tand(RakeR(ii))*cosd(Dip(ii)));
        end
         
        Trend(ii) = digit1(make360(Trend(ii)),0);
        DipDir(ii) = digit1(make360(Strike(ii) + 90),0);
        Dip1(ii) = Dip(ii);
    end
end
end
