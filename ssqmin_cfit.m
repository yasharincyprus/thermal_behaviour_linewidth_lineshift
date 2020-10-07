%%
function ssq = peakfitmin(par,fwhm,tdata,temp,Er,Ez,Ez1)
% controls how ssq is minimised. constrains the values of Debye temperature and e-p coupling strength

    lfit_val = cfit_output(temp,tdata,Er,Ez,Ez1,par(1,1),par(2,1),par(3,1),par(4,1),par(5,1),par(6,1));
    ssq = sum((lfit_val - fwhm).^2);
    
    % constrains Debye temperature to be positive
    pdw = 1e20*(par(1,1)<0);
    ssq = ssq*(1 + (par(1,1)^2)*pdw);
    
    
    % constrains e-p coupling strength to be positive
    pepw = 1e20*(par(2,1)<0);
    ssq = ssq*(1 + (par(2,1)^2)*pepw);
    
    %constrains baseline to be positive and below 10K FWHM
    pgf = 1e20*(par(3,1)<0);
    glf = 1e20*(par(3,1)>min(fwhm));
    ssq = ssq*(1 + (par(3,1)^2)*(pgf+glf));
    
    %constrains beta parameters to be positive
    pbrw = 1e20*(par(4,1)<0);
    ssq = ssq*(1 + (par(4,1)^2)*pbrw);
    
    pbzw = 1e20*(par(5,1)<0);
    ssq = ssq*(1 + (par(5,1)^2)*pbzw);
    
    pbz1w = 1e20*(par(6,1)<0);
    ssq = ssq*(1 + (par(6,1)^2)*pbz1w);
end
