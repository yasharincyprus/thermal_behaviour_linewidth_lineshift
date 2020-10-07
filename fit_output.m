%%
function lfit_val = fit_output(tdata,epco,debye,lw0)
% calculate the fit output values given the parameters and temperatures

    lfit_val = (epco/debye^7)*tdata + lw0;
    
end
