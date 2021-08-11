
function [] = plotStressInvResults(FMS, StressDirF,Pxb1,Pyb1,MSIG,title0,clr1)
%--------------
ang = 0:1:360;
Rr = 0.5;
dRr = 0.1*Rr;
ticks = [90 60 30 0];
Xo1 = Rr*cosd(ang);
Yo1 = Rr*sind(ang);
% Sig_c = [1 0 0; 0 1 0; 0 0 1];%'rgb';
Sig_c = parula(3);
fz = 12;
%-----------------------
strike = FMS(:,1);
dip = FMS(:,2);
[xSol1,ySol1] = StereoPlane(strike,dip,Rr) ;


if nargin > 4
    okbt = 1;
else
    okbt = 0;
end

if nargin == 11
    clr = clr1;
else
    clr = 'k';
end

hold on
plot([0 0],[Rr -Rr],'-k');
plot([Rr -Rr],[0 0],'-k');
for ii=1:length(ticks)
    text(0,Rr-ticks(ii)/90*Rr,num2str(ticks(ii)),'fontsize',fz);
end
% plot internal circles
for ii=2:(length(ticks)-1)
    Xo = ticks(ii)/90*Rr*cosd(ang);
    Yo = ticks(ii)/90*Rr*sind(ang);
    plot(Xo,Yo,':k','LineWidth',0.5)
end
plot(Xo1,Yo1,'-k','LineWidth',3)
patch(Xo1,Yo1,ones(1,3)*0.95)

n2=size(xSol1,1);
for ii=1:n2
    plot(xSol1(ii,:),ySol1(ii,:),clr,'LineWidth',0.5)
end


if Pxb1 ~= -999
    
    if okbt == 1
        for ii = 1:3 % plot Bootstrap Sigma
            %scatter(Pxb1(:,ii),Pyb1(:,ii),30,Sig_c(ii),'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.0)
            scatter(Pxb1(:,ii),Pyb1(:,ii),'MarkerFaceColor',Sig_c(ii,:),'MarkerEdgeColor',Sig_c(ii,:),'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.0)
        end
    end
    
    
    
    if okbt == 1
        sig1_TXT = ['\sigma_1 ',num2str(round(StressDirF(1,1))),' / ',num2str(round(StressDirF(1,2))),'  +/- ',num2str(round(MSIG(1)))];
        sig2_TXT = ['\sigma_2 ',num2str(round(StressDirF(2,1))),' / ',num2str(round(StressDirF(2,2))),'  +/- ',num2str(round(MSIG(2)))];
        sig3_TXT = ['\sigma_3 ',num2str(round(StressDirF(3,1))),' / ',num2str(round(StressDirF(3,2))),'  +/- ',num2str(round(MSIG(3)))];
    else
        sig1_TXT = ['\sigma_1 ',num2str(round(StressDirF(1,1))),' / ',num2str(round(StressDirF(1,2)))];
        sig2_TXT = ['\sigma_2 ',num2str(round(StressDirF(2,1))),' / ',num2str(round(StressDirF(2,2)))];
        sig3_TXT = ['\sigma_3 ',num2str(round(StressDirF(3,1))),' / ',num2str(round(StressDirF(3,2)))];
        
    end
    
    
else
    sig1_TXT = ['\sigma_1 ',num2str(round(StressDirF(1,1))),' / ',num2str(round(StressDirF(1,2)))];
    sig2_TXT = ['\sigma_2 ',num2str(round(StressDirF(2,1))),' / ',num2str(round(StressDirF(2,2)))];
    sig3_TXT = ['\sigma_3 ',num2str(round(StressDirF(3,1))),' / ',num2str(round(StressDirF(3,2)))];
    
end
% plot stress Sigma
[Px,Py] = StressDirStereo(StressDirF, Rr);
stress_symbols = 'os^';
for ii = 1:3
    hp(ii) = scatter(Px(ii),Py(ii),250,'MarkerFaceColor',Sig_c(ii,:),'Marker',stress_symbols(ii),'MarkerEdgeColor','k','linewidth',2);
end

% simple
% legend(hp,'\sigma_1','\sigma_2','\sigma_3')
% title(title0, 'fontsize',20)
% text(0, -0.7, title0, 'fontsize',12, 'HorizontalAlignment','center')


% with error
%PositionL = [0.7 0.7 0.3 0.2];
%newPosition = [c(1) + PositionL(1)*c(3), c(2)+ PositionL(2)*c(4), PositionL(3)*c(3), PositionL(4)*c(4)];
%newUnits = 'normalized';
hL = legend(hp,sig1_TXT,sig2_TXT,sig3_TXT,'fontsize',12,'location','south');
%set(hL,'Position', newPosition,'Units', newUnits);
title(hL,title0)

axis off
axis equal

end