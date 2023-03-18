function F = solveF(Z, numOfClasses)

numOfViews = length(Z);
numOfSamples = size(Z{1}, 1);
M = zeros(numOfSamples, numOfSamples);

%Compute M
for i=1:numOfViews
    %W = abs(Z{i}')*abs(Z{i});   %Compute similarity matrix,(abs(Z{i})+abs(Z{i}'))/2
    W = Z{i}'*Z{i};
    DN=diag( 1./sqrt(sum(W, 2)+eps) );   %Diagonal matrix
    LapN = speye(numOfSamples) - DN * W * DN;   %Normalized Laplacian matrix,LapN = diag(sum(W)) - W
   
    M = M + LapN;
end

%[V,D] = eig(M);
%[D_sort, ind] = sort(diag(D));
%ind2 = find(D_sort>1e-6);
%F = V(:, ind2(1:numOfClasses));  %the first c eigenvector

[~,~,vN] = svd(M);
FN = vN(:,numOfSamples-numOfClasses+1:numOfSamples);
for i = 1:numOfSamples
    F(i,:) = FN(i,:) ./ norm(FN(i,:)+eps);
end
