%Given the interval (a,b) and u, the new variable
% we get the old variable "x" & the transformation derivative
% x= (a+b)/2 +(b-a)/2 tanh(pi/2 sinh(u)); der= dx/du= (dx/dv)*(dv/du)
% Note that the transformation derivative is arrived at as follows.
% x= (b+a)/2 + (b-a)/2 tanh(v); v= (pi/2) sinh(u)
% dx/du= dx/dv * dv/du; dx/dv= 2(b-x)(x-a)/(b-a); dv/du = (pi/2) cosh(u)
% The term dvdu that is transferred here is just dv/du. It is multiplied by 
% dx/dv= 2(b-x)(x-a)/(b-a) in the CPV routine 
% 6-Feb-07
function [x,dvdu]=de(a,b,u);
format long;     
y=exp(u);
yi=exp(-u);
      x1=0.5*(y + yi);
      x2=0.5*(y - yi);
      x3=exp(pi*0.5*x2);
      x4=exp(-pi*0.5*x2);
      x=(b*x3+a*x4)/(x3+x4);
      dvdu=0.5*pi*x1;
      
      
