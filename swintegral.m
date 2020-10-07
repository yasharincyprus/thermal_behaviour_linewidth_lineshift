%%
function tdata = swintegral(temp,debye)
    fun = @(x) (x.^3) ./ (exp(x) - 1) ;
    upb = debye ./temp;
    lw_int = 0.*temp;
    for i = 1:length(upb)
        lw_int(i,1) = integral(fun,0,upb(i,1)); % evaluates the integral
    end
    tdata = (temp.^4).*lw_int;

end