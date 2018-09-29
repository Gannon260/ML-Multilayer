clc
tic
addpath('H:\SCOSand2Layer\functions');
constants

%addpath('/Volumes/PETER128/SCOSand2Layer/functions');
%change this to add noise
noise = true;

%ell = 1.0;
tau = DelayTime(1:2:80);
T = diff(DelayTime(1:2:81));
g= 0
repnumber = 0;

Db1 = 4.5e-9:.01e-9:5e-9;
Ratio = .30:.01:1;
ell = 0.9:.01:1.1;
Beta = 0.45:.001:0.50;

Db1s = Db1(:, randperm(size(Db1,2)));
ells = ell(:, randperm(size(ell,2)));
Ratio1 = Ratio(:, randperm(size(Ratio,2)));
j = 0;
Rep = 10;
siz = size(Db1s,2)*size(Ratio1,2)*size(ells,2)*size(Rep,2);
input = zeros(siz,3*size(tau,2));
target = zeros(siz,1);

load gauss_lag_5000.mat
for db1 = Db1s
    repnumber = repnumber+1
    for ratio = Ratio1
        for l = ells
            db2 = db1*10^ratio;
            
            g1_10 = diffusionforwardsolvergl(n,Reff,mua1,mus1,db1,tau,lambda,rho(1),w,l,mua2,mus2,db2,gl);
            normg1_10 = g1_10./g1_10(1);
            [b, index10] = min(abs(normg1_10-1/e)); %find where g1 = 1/e
            gamma = 1/tau(index10);
            sigma10 = getDCSNoise(500e3,T,3,beta,gamma,tau);
            noise10 = sigma10.*randn(length(tau),1)';
            g2_10nn = 1 + beta.*normg1_10.^2;
            
            g2_10 = noise10 + 1 + beta.*normg1_10.^2;
            g2_10sf = strokefilter(g2_10,5);
            
            sep15 = diffusionforwardsolvergl(n,Reff,mua1,mus1,db1,tau,lambda,1.5,w,l,mua2,mus2,db2,gl);
            normsep15 = sep15/sep15(1);
            [b, index15] = min(abs(normsep15-1/exp(1))); %find where g1 = 1/e
            gamma = 1/tau(index15);
            nsep15 = getDCSNoise(300e3,T,3,beta,gamma,tau); %50 hz.
            noise15 = nsep15.*randn(length(tau),1)';
            g2_15 = beta.*normsep15.*normsep15 + noise15  + 1;
%             
%             sep20 = diffusionforwardsolvergl(n,Reff,mua1,mus1,db1,tau,lambda,1.5,w,l,mua2,mus2,db2,gl);
%             normsep20 = sep20/sep20(1);
%             [b, index20] = min(abs(normsep15-1/exp(1))); %find where g1 = 1/e
%             gamma20 = 1/tau(index15);
%             nsep20 = getDCSNoise(80e3,T,3,beta,gamma20,tau); %50 hz.
%             noise20 = nsep20.*randn(length(tau),1)';
%             g2_20 = beta.*normsep20.*normsep20 + noise20  + 1;
%             
            g1_25 = diffusionforwardsolvergl(n,Reff,mua1,mus1,db1,tau,lambda,rho(2),w,l,mua2,mus2,db2,gl);
            normg1_25 = g1_25./g1_25(1);
            if noise == true
                [b, index25] = min(abs(normg1_25-1/e));
                gamma = 1/tau(index25);
                sigma25 = getDCSNoise(20e3,T,3,beta,gamma,tau);
            end
            for rep = Rep
                j = j+1;
                if noise == true
                    for ds = 1:7
                        noise25 = sigma25.*randn(length(tau),1)';
                        g2_25a(ds,:) = noise25 + 1 + beta.*normg1_25.^2;
                    end
                    g2_25mean = mean(g2_25a, 1);
                else
                    g2_25mean = 1 + beta.*normg1_25.^2;
                end
                input(j,:) = [g2_25mean g2_15 g2_10];
                target(j,:) = [db2*1e8];
            end
        end
    end
end
inputbeta = zeros(siz*size(Beta,2),size(input,2));
targetbeta = zeros(siz*size(Beta,2),1);
i = 0;
for bet = Beta
    indices = [1:siz] + siz*i;
    inputbeta(indices, :) = (input - 1).*bet/0.5 + 1;
    targetbeta(indices,:) = target;
    i = i + 1;
end

inputb = inputbeta;
targetb = targetbeta;
inputtarget = [inputb targetb];
inputtargetshuffle = inputtarget(randperm(size(inputtarget,1)),:);
inputshuffle = inputtargetshuffle(:, 1:size(input,2));
targetshuffle = inputtargetshuffle(:, size(input,2) + 1);
nnstart