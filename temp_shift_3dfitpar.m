clearvars; clc; close all;
aleph(1,1) = 495;
alphit(1,1) = 0;
ssq_out(1,1) = 0;
deb(1,1) = 0;

for i = 2:100
line_feats = csvread('temp_dep_fit.csv');
temp = line_feats(:,1);
pos = line_feats(:,2);

aleph(i,1) = aleph(i-1,1) + 1;
par(1,1) = -30;
par(2,1) = aleph(i,1); 
par(3,1) = interp1(temp,pos,0,'linear','extrap');


xData = swintegral(temp,par(2,1));
yData = pos;

tic;
options=optimset('Display','final','TolFun',1e-4,'TolX',1e-5,...
                'MaxFunEvals',5e4,'MaxIter',1e3);
            
Niter=1;maxiter=50;ssq=100*ones(length(maxiter),1);ep1=1e-3;

while Niter<50 && ssq(Niter)>0.001
    Niter=Niter+1;
   [parmin,fval,exitflag]=fminsearch(@sh_ssqmin_temp,par,options,yData,xData);
    npar=parmin;
    fit=sfit_output(xData,npar(1,1),npar(2,1),npar(3,1));
    ssq(Niter)=sum((fit-yData).^2);
    if (ssq(Niter)-ssq(Niter-1))/ssq(Niter-1)>-ep1
        break
    end
    t1=toc;
    disp(['iteration ',num2str(Niter-1),',time ',num2str(t1),' secs,ssq ',num2str(ssq(Niter))])
end

alphit(i,1) = npar(1,1);
deb(i,1) = npar(2,1);
ssq_out(i,1) = ssq(Niter);


end



figure
scatter3(alphit,ssq_out,deb,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75])
set(gca,'XLim',[0 1000],'YLim',[0 4],'ZLim',[0 1000])
view(40,35)
% rmsd = sqrt(ssq(Niter)/(length(pos) - length(par)));
% 
%  fig5 = figure;
%  sz2 = 10;
%  scatter(temp,pos,sz2,'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
%  ylimLower = min(pos)-1;
%  ylim_curr = max(pos) + 1;
%  set(gca, 'ylim', [ylimLower ylim_curr]) 
%  grid on
%  xlim([temp(1,1), temp(end,1)+10])
%  hold on;
% 
%  p1 = plot(temp,fit,'r-','LineWidth',1);
%  
% 
% err = rmsd .* pos;
% errorbar(temp,pos,err,'LineStyle','none')

% legend(p1,'Raman Scattering')