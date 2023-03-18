clear all;
load('UCI_3view.mat');
X=data;
gt = truth;
for i=1:3
    temp = X{i};
    X{i} = temp';
end

maxIters = 100;
lambda = 0.4;
beta = 0.0007;
C = size(unique(gt),1);

result = MvSCGE_solver(X, C, lambda, beta, maxIters, gt)
