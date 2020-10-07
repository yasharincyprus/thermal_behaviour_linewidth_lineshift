%%
function lfit_val = fit_output_df(tdata,epco,lw0,debye)
% calculate the fit output values given the parameters and temperatures

    lfit_val = (epco/debye^7)*tdata + lw0;
    
end
