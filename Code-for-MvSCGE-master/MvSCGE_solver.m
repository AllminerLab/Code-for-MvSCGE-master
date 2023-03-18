function [result] = MvSCGE_solver(X, C, lambda, beta, maxIters, y)
%X: 1 by n cell, each cell corresponds to each view
%C: number of classes

numOfViews = length(X);
numOfSamples = size(X{1}, 2);
MAXiter = 1000;
REPlic = 10; 

%Initialization Z
for i=1:numOfViews
    [Z{i}] = smooth_rep(X{i},lambda);
end

for it=1:maxIters
 
    %------update F--------
    F = solveF(Z, C);
    P = pdist2(F,F,'squaredeuclidean');
    
    %-------update Z---------
    for j=1:numOfViews
        Z{j} = SMR_mtv(X{j},P,lambda,beta);
    end
    
    %------print OBJ-------
    obj = 0;
    for j=1:numOfViews
        S = X{j}'*X{j};
        Ls = diag(sum(S)) - S;
        
        W = Z{j}'*Z{j};
        Dz=diag( 1./sqrt(sum(W, 2)+eps) ); 
        Lz = speye(numOfSamples) - Dz * W * Dz;
        
        obj = obj + norm((X{j} - X{j}*Z{j}),'fro')^2+lambda*trace(Z{j}*Ls*Z{j}')+beta*trace(F'*Lz*F);
    end
    Obj(it) = obj;
    if (it>1 && abs(Obj(it-1)-Obj(it)) < 10^-6)
        break;
    end
end

l = kmeans(F,C,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');

[ACC,NMI,PUR] = ClusteringMeasure(y,l);
[Fscore,Precision,R] = compute_f(y,l);
[AR,~,~,~]=RandIndex(y,l);
result = [ACC NMI PUR Fscore Precision R AR];






