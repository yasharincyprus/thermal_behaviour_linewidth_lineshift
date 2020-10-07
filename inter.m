%umax =upper de integration limit; -umax=lower de integration limit
% h=step size for DE; x0=demap(u0); x0=where you need the cpv integral
% in the old interval
%Between (u0+ h/2) and umax find out the integer(j) multiples of "h"
%The remainder is the fractional step size "hu"
%Integrate in u from (u0+h/2) to umax1;
%Add the contribution for the fractional width hu
% This is for the Left side of u0; For the Right side of u0
% Find the integer multiples of h between -umax and (u0+h/2);
%The remainder is the fraction step hl;
% Integrate from umin to umax1 & just add
% the contributions from the fractional widths on either side
function [umax1,hu,umin,hl]=inter(umax,u0,h)
format long;
z=umax-u0-(h/2);
hu=rem(z,h);
umax1=umax-hu;
z=u0+(h/2)+umax;
hl=rem(z,h);
umin=-umax+hl;