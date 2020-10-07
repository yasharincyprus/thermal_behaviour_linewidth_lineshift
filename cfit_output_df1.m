%%
function lfit_val = cfit_output_df1(temp,tdata,debye,E1,epco,lw0,b1)
% calculate the fit output values given the parameters and temperatures
    k = 0.695;
    lfit_val = (epco/debye^7).*tdata + lw0 + b1./(exp(E1./(k.*temp))-1);
    
end
