X = [a; a]
Y = a
[N, n] = size(X);

layers = [imageInputLayer([1 n]) fullyConnectedLayer(n) regressionLayer()];

opts = trainingOptions('sgdm');

net = trainNetwork(X, Y, layers, opts);
R = predict(net, X);