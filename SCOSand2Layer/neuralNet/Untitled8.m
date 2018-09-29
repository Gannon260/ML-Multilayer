clc
tic
addpath('H:\SCOSand2Layer\functions');
constants

%addpath('/Volumes/PETER1216/SCOSand2Layer/functions');
%change this to add noise
noise = true;

%ell = 1.0;
tau = DelayTime(1:1:100);
T = diff(DelayTime(1:1:101));
g= 0;
repnumber = 0;

Db1 = 5.00e-9:.01e-9:5.00e-9;
Ratio = .30:.01:1;
ell = 0.95:.002:1.05;
Beta = .5;

Rho = [10 12.5 15 17.5 20 22.5 25 27.5]./10

Db1s = Db1(:, randperm(size(Db1,2)));
ells = ell(:, randperm(size(ell,2)));
Ratio1 = Ratio(:, randperm(size(Ratio,2)));
j = 0;
Rep = 200;
siz = size(Db1s,2)*size(Ratio1,2)*size(ells,2)*Rep;
input = single(zeros(siz,size(tau,2)*2*4));
target = single(zeros(siz,3));

g2_25a = [];
inttime = 10;

load gauss_lag_5000.mat
for db1 = Db1s
    repnumber = repnumber+1
    tic
    for ratio = Ratio1
        for l = ells
            db2 = db1*10^ratio;

            g1_100 = diffusionforwardsolvergl(n,Reff,mua1,mus1,db1,tau,lambda,100e-2,w,l,mua2,mus2,db2,gl);
            normg1_10 = g1_100./g1_100(1);
            [b, index10] = min(abs(normg1_10-1/e)); %find where g1 = 1/e
            gamma = 1/tau(index10);
            sigma100 = getDCSNoise(500e3,T,inttime,beta,gamma,tau);
            
            g1_125 = diffusionforwardsolvergl(n,Reff,mua1,mus1,db1,tau,lambda,125e-2,w,l,mua2,mus2,db2,gl);
            normg1_125 = g1_125./g1_125(1);
            [b, index125] = min(abs(normg1_125-1/e)); %find where g1 = 1/e
            gamma = 1/tau(index125);
            sigma125 = getDCSNoise(500e3,T,inttime,beta,gamma,tau);
            
            sep150 = diffusionforwardsolvergl(n,Reff,mua1,mus1,db1,tau,lambda,150e-2,w,l,mua2,mus2,db2,gl);
            normsep150 = sep150/sep150(1);
            [b, index150] = min(abs(normsep150-1/exp(1))); %find where g1 = 1/e
            gamma = 1/tau(index150);
            nsep150 = getDCSNoise(200e3,T,inttime,beta,gamma,tau); %50 hz.
            
            g1_175 = diffusionforwardsolvergl(n,Reff,mua1,mus1,db1,tau,lambda,175e-2,w,l,mua2,mus2,db2,gl);
            normg1_175 = g1_175./g1_175(1);
            [b, index175] = min(abs(normg1_175-1/e)); %find where g1 = 1/e
            gamma = 1/tau(index175);
            sigma175 = getDCSNoise(500e3,T,inttime,beta,gamma,tau);
            
            sep200 = diffusionforwardsolvergl(n,Reff,mua1,mus1,db1,tau,lambda,200e-2,w,l,mua2,mus2,db2,gl);
            normsep200 = sep200/sep200(1);
            [b, index200] = min(abs(normsep200-1/exp(1))); %find where g1 = 1/e
            gamma200 = 1/tau(index200);
            nsep200 = getDCSNoise(90e3,T,inttime,beta,gamma200,tau); %50 hz.
            
            sep225 = diffusionforwardsolvergl(n,Reff,mua1,mus1,db1,tau,lambda,225e-2,w,l,mua2,mus2,db2,gl);
            normsep225 = sep225/sep225(1);
            [b, index225] = min(abs(normsep225-1/exp(1))); %find where g1 = 1/e
            gamma225 = 1/tau(index225);
            nsep225 = getDCSNoise(90e3,T,inttime,beta,gamma225,tau); %50 hz.
            
            g1_250 = diffusionforwardsolvergl(n,Reff,mua1,mus1,db1,tau,lambda,250e-2,w,l,mua2,mus2,db2,gl);
            normg1_250 = g1_250./g1_250(1);
            [b, index250] = min(abs(normg1_250-1/e));
            gamma = 1/tau(index250);
            sigma250 = getDCSNoise(30e3,T,inttime,beta,gamma,tau);
            
            g1_275 = diffusionforwardsolvergl(n,Reff,mua1,mus1,db1,tau,lambda,275e-2,w,l,mua2,mus2,db2,gl);
            normg1_275 = g1_275./g1_275(1);
            [b, index275] = min(abs(normg1_275-1/e));
            gamma = 1/tau(index275);
            sigma275 = getDCSNoise(30e3,T,inttime,beta,gamma,tau);

            for rep = 1:Rep
                j = j+1;
                noise100 = sigma100.*randn(length(tau),1)';
                noise125 = sigma125.*randn(length(tau),1)';
                noise150 = sigma100.*randn(length(tau),1)';
                noise175 = sigma175.*randn(length(tau),1)';
                noise100 = sigma100.*randn(length(tau),1)';
                noise100 = sigma100.*randn(length(tau),1)';
                noise100 = sigma100.*randn(length(tau),1)';
                noise100 = sigma100.*randn(length(tau),1)';
                
                
                input(j,:) = single([g2_30 g2_25a(:)' g2_20 g2_15 g2_10]);
                target(j,:) = single([db1*1e10 db2*1e9 l]);
            end
        end
    end
    toc
end
inputbeta = single(zeros(siz*size(Beta,2),size(input,2)+size(target,2)));
i = 0;
for bet = Beta
    indices = [1:siz] + siz*i;
    inputbet = (input - 1).*bet/0.5 + 1;
    inputbeta(indices, :) = [inputbet target];
    i = i + 1;
end
inputbeta = single(inputbeta(randperm(size(inputbeta,1)),:));
inputshuffle = single(inputbeta(:, 1:size(input,2)));

g2 = 0;
for detector = 1:ds
    indices = detector:ds:size(tau,2)*ds;
    g2 = g2 + inputbeta(:,indices);
end
avglong = g2./ds;
inputmean = single([avglong inputbeta(:, ds*size(tau,2)+1:size(input,2))]);
inputmean2 = single([avglong inputbeta(:, size(input,2)-size(tau,2):size(input,2))]);
inputshuffle2 = single(inputbeta(:, [1:size(tau,2)*ds size(input,2)-size(tau,2):size(input,2)]));

%targetshuffledb1 = inputtargetshuffle(:, 121);
targetshuffledb2 = inputbeta(:, size(input,2)+2);
%targetshuffleell = inputtargetshuffle(:, 123);
%targetshuffleall = inputtargetshuffle(:, 121:123);
%nnstart
clearvars -except inputshuffle targetshuffledb2