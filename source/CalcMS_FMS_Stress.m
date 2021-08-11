
function [misfit1n,instb_CFF,instb_pp,instb_critfric, misfit2, fault_tracR] = CalcMS_FMS_Stress(sig_Inv1, FMS,frc)
% Calculate angle between OBSERVED and CALCULATED Slip axis: SM
if ~isempty(FMS)
    strike = FMS(:,1);
    dip = FMS(:,2);
    rake = FMS(:,3);
    C = 0;
    [Svec,Nvec] = FMS2SLN3(strike,dip,rake);
    
    numevents = size(FMS,1);
    
    
    fault_trac(1:numevents,1:3)     = 0;
    fault_tracunit(1:numevents,1:3) = 0;
    fault_tracR(1:numevents)        = 0;
    normstress(1:numevents)         = 0;
    shearstress(1:numevents)        = 0;
    fault_shearvec(1:numevents,1:3) = 0;
    misfit1(1:numevents)            = 0;
    rakeS1(1:numevents)             = 0;
    misfit2(1:numevents)            = 0;
    instb_CFF(1:numevents)          = 0;
    instb_pp(1:numevents)           = 0;
    instb_critfric(1:numevents)     = 0;
    
    for i = 1:numevents
        
        % Solve for Tractions and Stress Maxima
        fault_trac(i,:) = (sig_Inv1*Nvec(i,:)')';
        fault_tracunit(i,:) = fault_trac(i,:)/norm(fault_trac(i,:));
        fault_tracR(i) = norm(fault_trac(i,:));
        
        normstress(i) = dot(Nvec(i,:)',fault_trac(i,:));
        shearstress(i) = (fault_tracR(i)^2-normstress(i)^2)^.5;
        
        % Solve for MaxiMum Shear Stress Direction and Rake
        fault_shearvec(i,:) = (fault_trac(i,:)-normstress(i)*Nvec(i,:))/shearstress(i);
        misfit1(i) = acosd(dot(fault_shearvec(i,:),Svec(i,:)'));
        rakeS1(i) = asind(-fault_shearvec(i,3)/sind(acosd(-Nvec(i,3))));        %Rake calculated from the solved shear vector
        misfit2(i) = abs(rake(i)-rakeS1(i));
        
        % determination of a quadrant
        cos_rakeS1 = fault_shearvec(i)*cos(strike(i)*pi/180)+fault_shearvec(i,2)*sin(strike(i)*pi/180);
        if (cos_rakeS1<0); rakeS1(i) = 180-rakeS1(i); end
        if (rakeS1(i) <-180); rakeS1(i) = rakeS1(i)+360; end
        if (rakeS1(i) > 180); rakeS1(i) = rakeS1(i)-360; end
        
        % Calculate Fault Instability in the Inversion Stress State
        % instb_CFF = Tau - fric*sigmanormal = Vertical distance to failure on Mohr plot (or critical cohesion)
        instb_CFF(i) = shearstress(i) - frc*(-1)*normstress(i) - C;
        
        % instb_pp = horizontal distance to failure, or pore pressure to failure
        instb_pp(i) = (-1)*normstress(i) - shearstress(i)/frc;
        
        % instb_critfric = friction coefficient that would cause failure
        instb_critfric(i) = abs(shearstress(i)/normstress(i));
        misfit1n(:,1) = misfit1;
    end
else
    misfit1n = [];
    
end

end
