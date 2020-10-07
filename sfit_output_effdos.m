%%
function lfit_val = sfit_output_effdos(tdata,temp,E1,epco,sh0,b1)
% calculate the fit output values given the parameters and integral values
    k = 0.695;
%     gamma = 1084;
    lfit_val = epco*tdata + sh0 + b1.*exp(E1./(k.*temp)) ./ ((exp(E1./(k.*temp))-1).^2); %+ ...
%         b2.*exp(E2./(k.*temp)) ./ ((exp(E2./(k.*temp))-1).^2);%...
        %+ b2 ./ (exp(E2./(k.*temp))-1);
        
    
end
