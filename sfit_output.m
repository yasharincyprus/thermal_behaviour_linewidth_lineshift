%%
function lfit_val = sfit_output(tdata,epco,debye,lw0)
% calculate the fit output values given the parameters and temperatures

    lfit_val = (epco/debye^4)*tdata + lw0;
    
end
