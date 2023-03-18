% By-Changqing Zhang
function [Z] =  SMR_mtv(X,P,lambda,beta)
%W = constructW_PKN(X,10);
W = X'*X;
L = diag(sum(W)) - W;
A_syl = X'*X;
B_syl = lambda*L + beta*P;
C_syl = -X'*X;

%X = lyap(A,B,C) solves the Sylvester equation,AX+XB+C=0 so here C is negative
Z = lyap(A_syl,B_syl,C_syl);

% CKSym = abs(Z)+abs(Z');
% C = SpectralClustering(CKSym,10);
% [A nmi avgent] = compute_nmi(gt,C);
% [f,p,r] = compute_f(gt,C);
end