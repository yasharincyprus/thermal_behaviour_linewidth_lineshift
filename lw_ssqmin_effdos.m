%%
function ssq = peakfitmin(par,fwhm,tdata,temp,E1,E2)
% controls how ssq is minimised. constrains the values of Debye temperature and e-p coupling strength

    lfit_val = lw_output_effdos(temp,tdata,E1,E2,par(1,1),par(2,1),par(3,1),par(4,1));
    ssq = sum((lfit_val - fwhm).^2);
    
    % constrains e-p coupling parameter to be positive
    pepw = 1e20*(par(1,1)<0);
    pdpw = 1e20*(par(1,1)>500);
    ssq = ssq*(1 + (par(1,1)^2)*(pepw+pdpw));
    
    %constrains baseline to be positive and below 10K FWHM
    pgf = 1e20*(par(2,1)<0);
    glf = 1e20*(par(2,1)>min(fwhm));
    ssq = ssq*(1 + (par(2,1)^2)*(pgf+glf));
    
    %constrains beta parameters and phonon energies to be positive
    pb1w = 1e20*(par(3,1)<0);
    ssq = ssq*(1 + (par(3,1)^2)*pb1w);
    
    pb2w = 1e20*(par(4,1)<0);
    ssq = ssq*(1 + (par(4,1)^2)*pb2w);
%     
%     pb3w = 1e20*(par(5,1)<0);
%     ssq = ssq*(1 + (par(5,1)^2)*pb3w);
%     
%     pb4w = 1e20*(par(6,1)<0);
%     ssq = ssq*(1 + (par(6,1)^2)*pb4w);
%     
%     pb5w = 1e20*(par(7,1)<0);
%     ssq = ssq*(1 + (par(7,1)^2)*pb5w);
%     
%     pb6w = 1e20*(par(8,1)<0);
%     ssq = ssq*(1 + (par(8,1)^2)*pb6w);
    
end
