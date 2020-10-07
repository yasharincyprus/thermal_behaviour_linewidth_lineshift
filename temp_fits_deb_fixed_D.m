clearvars; clc; close all;

aleph(1,1) = 20;
alphit(1,1) = 20;
ssq_out(1,1) = 0;
incr = 1;

line_feats = csvread('temp_dep_fit.csv');
temp = line_feats(:,1);
fwhm = line_feats(:,4);
    
debye = 351;

for j = 2:800

    aleph(j,1) = aleph(j-1,1) + incr; 
    par(1,1) = aleph(j,1);
    par(2,1) = interp1(temp,fwhm,0,'linear','extrap');


    xData = lwintegral(temp,debye);
    yData = fwhm;

    tic;
    options=optimset('Display','final','TolFun',1e-4,'TolX',1e-5,...
                    'MaxFunEvals',5e4,'MaxIter',1e4);

    Niter=1;maxiter=50;ssq=100*ones(length(maxiter),1);ep1=1e-3;

    while Niter<50 && ssq(Niter)>0.001
        Niter=Niter+1;
       [parmin,fval,exitflag]=fminsearch(@ssqmin_temp,par,options,yData,xData,debye);
        npar=parmin;
        fit=fit_output(debye,xData,npar(1,1),npar(2,1));
        ssq(Niter)=sum((fit-yData).^2);
        if (ssq(Niter)-ssq(Niter-1))/ssq(Niter-1)>-ep1
            break
        end
        t1=toc;
        disp(['iteration ',num2str(Niter-1),',time ',num2str(t1),' secs,ssq ',num2str(ssq(Niter))])
    end

    alphit(j,1) = npar(1,1);
    ssq_out(j,1) = ssq(Niter);
end

figure
scatter(aleph,ssq_out,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75])

