%%
function tdata = swintegral_effdos1(kelvin,cutoff)
    boltz = 0.695;
    
       a1 =      0.6765;
       b1 =        1026;
       c1 =       26.36;
       a2 =           0;
       b2 =         929;
       c2 =       26.55;
       a3 =      1.0000;
       b3 =       997.7;
       c3 =       88.54;
       a4 =       0.4112;
       b4 =       789.7;
       c4 =       19.18;
       a5 =      0.5986;
       b5 =         543;
       c5 =       15.79;
       a6 =     0.3451;
       b6 =       673.4;
       c6 =       7.465;
       a7 =      0.5574;
       b7 =       694.4;
       c7 =         333;
       a8 =      0.6014;
       b8 =       772.4;
       c8 =       41.33;
    
    func = @(x) (a1*exp(-((x-b1)./c1).^2) + a2*exp(-((x-b2)./c2).^2) + ...
              a3*exp(-((x-b3)./c3).^2) + a4*exp(-((x-b4)./c4).^2) + ...
              a5*exp(-((x-b5)./c5).^2) + a6*exp(-((x-b6)./c6).^2) + ...
              a7*exp(-((x-b7)./c7).^2) + a8*exp(-((x-b8)./c8).^2)) ./ (-2*pi.*(exp(x ./ (kelvin*boltz))-1));
    
    upb = cutoff;
%     sw_int = 0.*temp;
    
%     for i = 1:length(temp)
%         sw_int(i,1) = integral(func1,0,upb); % evaluates the integral
%     end
    sw_int = integral(func,1,upb);
    tdata = sw_int;

end