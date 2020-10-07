%%
function tdata = swintegral_effdos2(kelvin,cutoff)
    boltz = 0.695;
    
       a1 =      0.1156;
       b1 =        1022;
       c1 =       15.62;
       a2 =       0.897;
       b2 =        1011;
       c2 =       70.92;
       a3 =     -0.1243;
       b3 =       990.6;
       c3 =       9.912;
       a4 =      0.2389;
       b4 =       872.7;
       c4 =       35.25;
       a5 =      0.4437;
       b5 =       783.8;
       c5 =       46.23;
       a6 =       0.215;
       b6 =       705.1;
       c6 =       75.28;
    
    
    func = @(x) (a1*exp(-((x-b1)./c1).^2) + a2*exp(-((x-b2)./c2).^2) + ...
              a3*exp(-((x-b3)./c3).^2) + a4*exp(-((x-b4)./c4).^2) + ...
              a5*exp(-((x-b5)./c5).^2) + a6*exp(-((x-b6)./c6).^2)) ./ (-2*pi*(exp(x ./ (kelvin*boltz))-1));
    
    upb = cutoff;
    
    sw_int = integral(func,1,upb);
    tdata = sw_int;

end