%%
function lfit_val = lw_output_effdos(temp,tdata,E1,E2,epco,lw0,b1,b2)
% calculate the fit output values given the parameters and temperatures
    k = 0.695;
    lfit_val = epco.*tdata + lw0 + b1./(exp(E1./(k.*temp))-1); ...
        + b2./(exp(E2./(k.*temp))-1);
    
end
