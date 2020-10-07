%%
function lfit_val = nsfit_output(temp,lw0,gamma,brr,Er)
% calculate the fit output values given the parameters and temperatures
    k = 0.695;
    lfit_val =lw0+((2*brr/gamma).*exp(Er/k.*temp))./((exp(Er/k.*temp)-1).^2); %+ ...
%         ((2*bzz^2/gamma).*exp(Ez./k.*temp))./(exp(Ez./k.*temp)-1).^2 ...
%         + ((2*bz1^2/gamma).*exp(Ez1./k.*temp))./(exp(Ez1./k.*temp)-1).^2;
    
end
