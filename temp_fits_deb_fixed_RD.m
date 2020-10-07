clearvars; clc; close all;

aleph(1,1) = 20;
alphit(1,1) = 20;
ssq_out(1,1) = 0;
incr = 1;

line_feats = csvread('temp_dep_fit.csv');
temp = line_feats(:,1);
fwhm = line_feats(:,4);


Er = 173;
% Ez = 472;
% Ez1 = 472;
% Ez2 = 357;
% Ez3 = 472;
debye = 351;

for j = 2:800


    aleph(j,1) = aleph(j-1,1) + incr;
    par_n(1,1) = aleph(j,1);
    par_n(2,1) = interp1(temp,fwhm,0,'linear','extrap');
    par_n(3,1) = 0;
%     par_n(4,1) = 0;
%     par_n(5,1) = 0;
%     par_n(6,1) = 0;
%     par_n(7,1) = 0;

    xData = lwintegral(temp,debye);
    yData = fwhm;

 
    tic;
    options_n=optimset('Display','final','TolFun',1e-4,'TolX',1e-5,...
                    'MaxFunEvals',5e4,'MaxIter',1e3);

    Niter=1;maxiter=50;ssq_n=100*ones(length(maxiter),1);ep1=1e-3;

    while Niter<50 && ssq_n(Niter)>0.001
        Niter=Niter+1;
       [parmin_n,fval_n,exitflag_n]=fminsearch(@ssqmin_cfit_df,par_n,options_n,yData,xData,temp,debye,Er);
        npar_n=parmin_n;
        fit_n=cfit_output_df1(temp,xData,debye,Er,npar_n(1,1),npar_n(2,1),npar_n(3,1));
        ssq_n(Niter)=sum((fit_n-yData).^2);
        if (ssq_n(Niter)-ssq_n(Niter-1))/ssq_n(Niter-1)>-ep1
            break
        end
        t1=toc;
        disp(['iteration ',num2str(Niter-1),',time ',num2str(t1),' secs,ssq ',num2str(ssq_n(Niter))])
    end


    alphit(j,1) = npar_n(1,1);
    ssq_out(j,1) = ssq_n(Niter);
end

figure
scatter(aleph,ssq_out,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75])

