[a,b] = find(debfit>=575 & debfit<=585);
for i = 1:length(a)
    exc_debfit(i,1) = debfito(a(i,1),b(i,1));
    exc_alphit(i,1) = alphito(a(i,1),b(i,1));
    exc_ssqout(i,1) = ssq_outo(a(i,1),b(i,1));
    if (exc_ssqout(i,1) <= 0.0528) && (exc_alphit(i,1)>= -200)
        alpha_final = alepho(a(i,1),b(i,1));
        T_D = debo(a(i,1),b(i,1));
        temp_shift
        numerator = length(a) - i;
        disp(['number of remaining fits = ',num2str(numerator)])
        pause
        close all
    end
    exc_aleph(i,1) = alepho(a(i,1),b(i,1));
    exc_deb(i,1) = debo(a(i,1),b(i,1));
end

excavated = [a b exc_debfit exc_alphit exc_ssqout exc_aleph exc_deb];

