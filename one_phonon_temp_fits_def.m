clearvars; clc; close all;


line_feats = csvread('temp_dep_fit.csv');
temp = line_feats(:,1);
fwhm = line_feats(:,4);

debye = 580.3;

E1 = 179.3;



par(1,1) = 420;
par(2,1) = interp1(temp,fwhm,0,'linear','extrap');



par_n(1,1) = 539;
par_n(2,1) = interp1(temp,fwhm,0,'linear','extrap');
par_n(3,1) = 0;



xData = lwintegral(temp,debye);
yData = fwhm;

tic;
options=optimset('Display','final','TolFun',1e-4,'TolX',1e-5,...
                'MaxFunEvals',5e4,'MaxIter',1e4);
            
Niter=1;maxiter=50;ssq=100*ones(length(maxiter),1);ep1=1e-3;

while Niter<50 && ssq(Niter)>0.001
    Niter=Niter+1;
   [parmin,fval,exitflag]=fminsearch(@ssqmin_temp_df,par,options,yData,xData,debye);
    npar=parmin;
    fit=fit_output_df(xData,npar(1,1),npar(2,1),debye);
    ssq(Niter)=sum((fit-yData).^2);
    if (ssq(Niter)-ssq(Niter-1))/ssq(Niter-1)>-ep1
        break
    end
    t1=toc;
    disp(['iteration ',num2str(Niter-1),',time ',num2str(t1),' secs,ssq ',num2str(ssq(Niter))])
end
rmsd = sqrt(ssq(Niter)/(length(fwhm) - length(par)));

 fig5 = figure;
 sz2 = 10;
 scatter(temp,fwhm,sz2,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
 ylimLower = 0;
 ylim_curr = get(gca, 'ylim');
 set(gca, 'ylim', [ylimLower ylim_curr(2)]) 
 grid on
 xlim([temp(1,1), temp(end,1)])
 hold on;
%  errorbar(temp,fwhm,err,'LineStyle','none')
 p1 = plot(temp,fit,'r-','LineWidth',1);
 
tic;
options_n=optimset('Display','final','TolFun',1e-4,'TolX',1e-5,...
                'MaxFunEvals',5e4,'MaxIter',1e3);
            
Niter=1;maxiter=50;ssq_n=100*ones(length(maxiter),1);ep1=1e-3;

while Niter<50 && ssq_n(Niter)>0.001
    Niter=Niter+1;
   [parmin_n,fval_n,exitflag_n]=fminsearch(@one_phonon_ssqmin_cfit_def,par_n,options_n,yData,xData,temp,debye,E1);
    npar_n=parmin_n;
    fit_n=one_phonon_cfit_output_def(temp,xData,debye,E1,npar_n(1,1),npar_n(2,1),npar_n(3,1));
    ssq_n(Niter)=sum((fit_n-yData).^2);
    if (ssq_n(Niter)-ssq_n(Niter-1))/ssq_n(Niter-1)>-ep1
        break
    end
    t1=toc;
    disp(['iteration ',num2str(Niter-1),',time ',num2str(t1),' secs,ssq ',num2str(ssq_n(Niter))])
end
rmsd_n = sqrt(ssq_n(Niter)/(length(fwhm) - length(par_n)));

err = rmsd_n .* fwhm;
% errorbar(temp,fwhm,err,'LineStyle','none')
p2 = plot(temp,fit_n,'g-.','LineWidth',1);
legend([p1 p2],'Raman Scattering','Raman scattering + Direct Processes')