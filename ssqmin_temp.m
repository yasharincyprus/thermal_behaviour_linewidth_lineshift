%%
function ssq = peakfitmin(par,fwhm,temp)
% controls how ssq is minimised. constrains the values of Debye temperature and e-p coupling strength

    lfit_val = fit_output(temp,par(1,1),par(2,1),par(3,1));
    ssq = sum((lfit_val - fwhm).^2);
    
    % constrains Debye temperature to be positive
    pdw = 1e20*(par(1,1)<0);
    ssq = ssq*(1 + (par(1,1)^2)*pdw);
    
    paw = 1e20*(par(2,1)<0);
    ssq = ssq*(1 + (par(2,1)^2)*paw);
    
    %constrains baseline to be positive and below 10K FWHM
    pgf = 1e20*(par(2,1)<0);
    glf = 1e20*(par(2,1)>min(fwhm));
    ssq = ssq*(1 + (par(2,1)^2)*(pgf+glf));
    
end
