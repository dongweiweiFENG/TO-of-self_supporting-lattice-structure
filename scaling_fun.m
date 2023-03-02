function [x,dx]=scaling_fun(xPhys,aa,bb)
indx1=xPhys-aa>=0;
indx2=xPhys-aa<0;
x(indx1)=xPhys(indx1);
x(indx2)=xPhys(indx2).*(1-tanh(bb.*(xPhys(indx2)-aa(indx2))).*tanh(bb.*(xPhys(indx2)-aa(indx2))));
dx(indx1)=1;
dx(indx2)=1-tanh(bb*(xPhys(indx2)-aa(indx2))).*tanh(bb*(xPhys(indx2)-aa(indx2)))+...
    xPhys(indx2).*(-2*bb*tanh(bb.*(xPhys(indx2)-aa(indx2)))).*(1-tanh(bb.*(xPhys(indx2)-aa(indx2))).*tanh(bb.*(xPhys(indx2)-aa(indx2))));
x=x';
dx=dx';
% x2(indx1)=xPhys(indx1);
% eps=1e-6;
% x2(indx2)=(xPhys(indx2)+eps).*(1-tanh(bb.*((xPhys(indx2)+eps)-aa(indx2))).*tanh(bb.*((xPhys(indx2)+eps)-aa(indx2))));
% dx2=(x2'-x)/eps;
end