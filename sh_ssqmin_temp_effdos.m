%%
function ssq = peakfitmin(par,shift,tdata,temp,E1)
% controls how ssq is minimised. constrains the value of e-p coupling strength

    lfit_val = sfit_output_effdos(tdata,temp,E1,par(1,1),par(2,1),par(3,1));
    ssq = sum((lfit_val - shift).^2);
    
    
%     constrains e-p coupling strength to be positive
    pdw = 1e20*(par(1,1)<-500);
    pgw = 1e20*(par(1,1)>0);
    ssq = ssq*(1 + (par(1,1)^2)*(pdw+pgw));
    
    %constrains baseline to be less than zero and below 10K shift
    pgf = 1e20*(par(2,1)<0);
%     glf = 1e20*(par(3,1)<max(shift));
    ssq = ssq*(1 + (par(2,1)^2)*pgf);
    
    pbrw = 1e20*(par(3,1)>0);
    ssq = ssq*(1 + (par(3,1)^2)*pbrw);
% % % %     
%     psrw = 1e20*(par(4,1)>0);
%     ssq = ssq*(1 + (par(4,1)^2)*psrw);
% %     
%     pbrw = 1e20*(par(5,1)>0);
%     ssq = ssq*(1 + (par(5,1)^2)*pbrw);
end
