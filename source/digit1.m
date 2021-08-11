function [num] = digit1(x,n)

 dub = 10^n;
 temp1 = x*dub;
 num = round(temp1)/dub;
