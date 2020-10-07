clearvars; clc; close all;


line_feats = csvread('temp_dep_fit.csv');
temp = line_feats(:,1);
fwhm = line_feats(:,4);

cutoff = 1084;

E1 = 88;
E2 = 111.5;


par_n(1,1) = 100;
par_n(2,1) = interp1(temp,fwhm,0,'linear','extrap');
par_n(3,1) = 1;
par_n(4,1) = 1;



xData = 0.*temp;
for i = 1:numel(temp)
    xData(i,1) = lwintegral_effdos(temp(i,1),cutoff);
end

yData = fwhm;

tic;
options_n=optimset('Display','final','TolFun',1e-4,'TolX',1e-5,...
                'MaxFunEvals',5e4,'MaxIter',1e3);
            
Niter=1;maxiter=50;ssq_n=100*ones(length(maxiter),1);ep1=1e-3;

while Niter<50 && ssq_n(Niter)>0.001
    Niter=Niter+1;
   [parmin_n,fval_n,exitflag_n]=fminsearch(@lw_ssqmin_effdos,par_n,options_n,yData,xData,temp,E1,E2);
    npar_n=parmin_n;
    fit_n=lw_output_effdos(temp,xData,E1,E2,npar_n(1,1),npar_n(2,1),npar_n(3,1),npar_n(4,1));
    ssq_n(Niter)=sum((fit_n-yData).^2);
    if (ssq_n(Niter)-ssq_n(Niter-1))/ssq_n(Niter-1)>-ep1
        break
    end
    t1=toc;
    disp(['iteration ',num2str(Niter-1),',time ',num2str(t1),' secs,ssq ',num2str(ssq_n(Niter))])
end

fig1 = figure;
 sz2 = 10;
 scatter(temp,fwhm,sz2,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
 ylimLower = 0;
 ylim_curr = get(gca, 'ylim');
 set(gca, 'ylim', [ylimLower ylim_curr(2)]) 
 grid on
 xlim([temp(1,1), temp(end,1)])
 hold on;
plot(temp,fit_n,'r-','LineWidth',1)
 



