subplot(2,1,1);
noisy = inputshuffle2(1,:)
netfitted = denoisenet_10_500_80(inputshuffle2(1,:));
nonoise = dtargetshuffle2(1,:)
semilogx(inputshuffle2(1,:)); 
hold on; plot(netfitted); 
hold on; plot(dtargetshuffle2(1,:)); 

legend("noise","net corrected", "no-noise");

subplot(2,1,2);
noiseerror = (noisy-nonoise)./nonoise.*100;
neterror = (netfitted - nonoise)./nonoise.*100;
plot(noiseerror); hold on;
plot(neterror);