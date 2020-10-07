%%
function ssq = peakfitmin(par,fwhm,temp,debye)
% controls how ssq is minimised. constrains the values of Debye temperature and e-p coupling strength

    lfit_val = fit_output_df(temp,par(1,1),par(2,1),debye);
    ssq = sum((lfit_val - fwhm).^2);
    
    
    
    % constrains e-p coupling strength to be positive
    pepw = 1e20*(par(1,1)<0);
    ssq = ssq*(1 + (par(1,1)^2)*pepw);
    
    
    %constrains baseline to be positive and below 10K FWHM
    pgf = 1e20*(par(2,1)<0);
    glf = 1e20*(par(2,1)>min(fwhm));
    ssq = ssq*(1 + (par(2,1)^2)*(pgf+glf));
    
end
