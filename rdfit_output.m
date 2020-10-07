%%
function lfit_val = rdfit_output(tdata,debye,epco,lw0,Eij)
% calculate the fit output values given the parameters and temperatures

    lfit_val = (epco/debye^7)*tdata +  + lw0;
    
end