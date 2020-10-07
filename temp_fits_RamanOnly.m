clearvars; clc; close all;


line_feats = csvread('temp_dep_fit.csv');
temp = line_feats(:,1);
fwhm = line_feats(:,4);
yData = fwhm;

par(3,1) = interp1(temp,fwhm,0,'linear','extrap');

sz = [200 430];
deb = 1.*ones(sz(1),sz(2));
aleph = 1.*ones(sz(1),sz(2));
alphit = 0.*ones(sz(1),sz(2));
debfit = 0.*ones(sz(1),sz(2));
ssq_out = 0.*ones(sz(1),sz(2));


for i = 1:sz(1)
    
    deb(i,:) = deb(i,:) + i;    
    par(2,1) = deb(i,1); 
    
    xData = lwintegral(temp,par(2,1));

    
    for j = 1:sz(2)

        aleph(i,j) = aleph(i,j) + j;
        par(1,1) = aleph(i,j);

        

        tic;
        options=optimset('Display','final','TolFun',1e-4,'TolX',1e-5,...
                        'MaxFunEvals',5e4,'MaxIter',1e4);

        Niter=1;maxiter=50;ssq=100*ones(length(maxiter),1);ep1=1e-3;

        while Niter<50 && ssq(Niter)>0.001
            Niter=Niter+1;
           [parmin,fval,exitflag]=fminsearch(@ssqmin_temp_df,par,options,yData,xData);
            npar=parmin;
            fit=fit_output_df(xData,npar(1,1),npar(2,1),npar(3,1));
            ssq(Niter)=sum((fit-yData).^2);
            if (ssq(Niter)-ssq(Niter-1))/ssq(Niter-1)>-ep1
                break
            end
            t1=toc;
            disp(['iteration ',num2str(Niter-1),',time ',num2str(t1),' secs,ssq ',num2str(ssq(Niter))])
        end

        alphit(i,j) = npar(1,1);
        debfit(i,j) = npar(2,1);
        ssq_out(i,j) = ssq(Niter);
    end
end

figure(1)
surf(alphit,debfit,ssq_out,'EdgeColor','none','facecolor','interp')
alpha 0.7
% set(gca,'XLim',[0 1000],'YLim',[0 4],'ZLim',[0 1000])
view(40,35)

figure(2)
surf(aleph,deb,ssq_out,'EdgeColor','none','facecolor','interp')
alpha 0.7
view(40,35)

t2=toc;
disp(['total time = ',num2str(t2)])