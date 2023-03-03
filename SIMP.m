%  the penalization of the struts density
%======input===========
% xPhys is the density of the struts
% p is the penalization parameter

%======output===========
% x and dx is the dnesity and gradient information after penalization
function [x,dx]=SIMP(xPhys,p)
indx1=xPhys-1>=0;
indx2=xPhys-1<0;
x(indx1)=xPhys(indx1);
x(indx2)=xPhys(indx2).^p;
dx(indx1)=1;
dx(indx2)=p*xPhys(indx2).^(p-1);
x=x';
dx=dx';
end
