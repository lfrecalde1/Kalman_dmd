function [sysmodel_DMDc] = Ident(X,U,dt)
X1 = X(:,1:end-1);
X2 = X(:,2:end);
% X3 = Xp(:,2:end);
numVar = size(X1,1); numOutputs = numVar; numInputs = size(U,1);
Gamma = U(:,1:end-1);
Omega = [X1; Gamma];

[U,S,V] = svd(Omega,'econ'); 
G = X2*V*S^(-1)*U';
% G2 = X3*V*S^(-1)*U';

% A = G2(:,1:numVar)
% B = G2(:,numVar+1:end)
% [Up,Sp,Vp] = svd(X2,'econ'); 
Ar = G(:,1:numVar);
Br = G(:,numVar+1:end);
Cr = eye(numVar); 
sysmodel_DMDc = ss(Ar,Br,Cr,zeros(numOutputs,numInputs),dt);