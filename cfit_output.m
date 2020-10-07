%%
function lfit_val = cfit_output(temp,tdata,Er,Ez,Ez1,debye,epco,lw0,brr,bzz,bzz1)
% calculate the fit output values given the parameters and temperatures
    k = 0.695;
    lfit_val = (epco/debye^7).*tdata + lw0 + brr./(exp(Er./(k.*temp))-1) + bzz./(exp(Ez./(k.*temp))-1) + ...
        bzz1./(exp(Ez1./(k.*temp))-1);
    
end
