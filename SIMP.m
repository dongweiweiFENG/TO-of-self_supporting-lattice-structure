%  the penalization of the struts density
% p is the penalization parameter
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
