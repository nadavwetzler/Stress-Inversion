function [x1] = make360(x)
n1= length(x);
x1(1:n1) = 0;
for ii=1:n1
    if (x(ii) < 0)
        x1(ii) = 360 + x(ii);
    elseif (x(ii) > 360)
        x1(ii) = x(ii) - 360;
    else
        x1(ii) = x(ii);
    end
    
end

