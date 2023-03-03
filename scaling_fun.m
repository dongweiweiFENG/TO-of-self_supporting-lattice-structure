% size filter for simplification
% xPhys is the density of the struts
% aa is the lower bound of the virtual density value
% bb is the parameters controls the sharpness of the projection function
% x and dx is the density and its gradient information after the projection
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
end
