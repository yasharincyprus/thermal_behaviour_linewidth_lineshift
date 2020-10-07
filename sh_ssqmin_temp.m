%%
function ssq = peakfitmin(par,pos,temp)
% controls how ssq is minimised. constrains the values of Debye temperature and e-p coupling strength

    lfit_val = sfit_output(temp,par(1,1),par(2,1),par(3,1));
    ssq = sum((lfit_val - pos).^2);
    
    % constrains Debye temperature to be positive
    pepw = 1e20*(par(1,1)>0);
    ssq = ssq*(1 + (par(1,1)^2)*pepw);
    
    
%     constrains e-p coupling strength to be positive
    pdw = 1e20*(par(2,1)<0);
    ssq = ssq*(1 + (par(2,1)^2)*pdw);
    
    %constrains baseline to be positive and below 10K FWHM
    pgf = 1e20*(par(3,1)<0);
    glf = 1e20*(par(3,1)<max(pos));
    ssq = ssq*(1 + (par(3,1)^2)*(pgf+glf));
    
end
