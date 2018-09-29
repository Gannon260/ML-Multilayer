 
net1 = feedforwardnet(10);
count = gpuDeviceCount;
%net1.trainFcn = 'trainbr';
a = 1:1e5;
x = double(inputshuffle(a,:))';
t =  double(targetshuffledb2(a,:))';
gpu1 = gpuDevice(1);

net2 = train(net1,x,t,'useGPU','yes');
y = net2(x,'useGPU','yes');