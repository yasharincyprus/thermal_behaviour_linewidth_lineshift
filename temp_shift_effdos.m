clearvars; clc; close all;

cutoff = 1084;
line_feats = csvread('temp_dep_fit.csv');
temp = line_feats(:,1);
shift = line_feats(:,6);

E1 = 88;
% E2 = 173;
% E3 = 120.3;

par(1,1) = -10;
par(2,1) = interp1(temp,shift,0,'linear','extrap');
par(3,1) = -1;
% par(4,1) = -1;
% par(5,1) = -0.72;

xData = 0.*temp;
for i = 1:numel(temp)
    xData(i,1) = swintegral_effdos2(temp(i,1),cutoff);
end
yData = shift;

tic;
options=optimset('Display','final','TolFun',1e-4,'TolX',1e-5,...
                'MaxFunEvals',5e4,'MaxIter',1e4);
            
Niter=1;maxiter=50;ssq=100*ones(length(maxiter),1);ep1=1e-3;

while Niter<50 && ssq(Niter)>0.001
    Niter=Niter+1;
   [parmin,fval,exitflag]=fminsearch(@sh_ssqmin_temp_effdos,par,options,yData,xData,temp,E1);
    npar=parmin;
    fit=sfit_output_effdos(xData,temp,E1,npar(1,1),npar(2,1),npar(3,1));
    ssq(Niter)=sum((fit-yData).^2);
    if (ssq(Niter)-ssq(Niter-1))/ssq(Niter-1)>-ep1
        break
    end
    t1=toc;
    disp(['iteration ',num2str(Niter-1),',time ',num2str(t1),' secs,ssq ',num2str(ssq(Niter))])
end
rmsd = sqrt(ssq(Niter)/(length(shift) - length(par)));

 fig5 = figure;
 sz2 = 10;
 scatter(temp,shift,sz2,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
 ylimLower = min(shift)-1;
 ylim_curr = max(shift) + 1;
 set(gca, 'ylim', [ylimLower ylim_curr]) 
 grid on
 xlim([temp(1,1), temp(end,1)+10])
 hold on;

p1=plot(temp,fit,'r-','LineWidth',1);
% title(['(a,b) = (',num2str(a(i,1)),',',num2str(b(i,1)),'), T_D = ',num2str(npar(2,1)),', \alpha = ',num2str(npar(1,1))]) 

% err = rmsd .* shift;
% errorbar(temp,pos,err,'LineStyle','none')

% legend(p1,'Raman Scattering')