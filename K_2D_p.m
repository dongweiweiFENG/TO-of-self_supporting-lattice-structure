% Calculate compliance and corresponding gradient information
% ===output==========
% c is the compliance
% dcx is the gradient information about the density variable of the compliance

% ====input==========
% ssk -- the information about matrix of the struts
% edofMat -- the information about the relationship about the nodes and struts
% xPhys is the density of the struts
% num is the number of the struts
% F is the external force vector
% freedofs is the degree of freedom vector of nodes
% p is the penalization parameter
% U is the global displacement vector

function [c,dcx]=K_2D_p(ssk,num,F,freedofs,xPhys,edofMat,p,U)
iK = reshape(kron(edofMat,ones(4,1))',16*num,1);
jK = reshape(kron(edofMat,ones(1,4))',16*num,1);
[x,dx]=SIMP(xPhys,p);
x0=kron(x,ones(4,1))';
sk=ssk.*(x0+1e-12);
sK = full(reshape(sk(:),16*(num),1));
K = sparse(iK,jK,sK);
K = (K+K')/2;
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
clear K;
c=F(freedofs)'*U(freedofs);
u=U(edofMat)';
k1=ssk(:,1:4:end); k2=ssk(:,2:4:end); k3=ssk(:,3:4:end); k4=ssk(:,4:4:end);
u1=sum(u.*k1);u2=sum(u.*k2);u3=sum(u.*k3);u4=sum(u.*k4);
ku=dx'.*[u1;u2;u3;u4];
dcx=-sum(ku.*u);
clear ku;
