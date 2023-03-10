% the stiffness matrix of the struts and the information about the nodes and struts(edofMat)
% ssk -- the stiffness matrix of the struts
% edofMat -- the information about the nodes and struts

% x1,x3,x4 is the angles about the struts and printing direction
% nelx x nely is the resolution of the ground lattice structure
% E0 is the Young's modulus of the struts
% A0 is the initial cross-section of the struts
% l is the length of the minimum struts in the virtical direction
function [ssk,edofMat]=Bar2D_Stiffness(x1,x3,x4,nelx,nely,E0,A0,l) 
%x0=0
C1=cos(x1); 
S1=sin(x1); 
k01=[C1*C1 C1*S1 -C1*C1 -C1*S1; C1*S1 S1*S1 -C1*S1 -S1*S1; 
 -C1*C1 -C1*S1 C1*C1 C1*S1; -C1*S1 -S1*S1 C1*S1 S1*S1];
%x3=-pi/4
C3=cos(x3); 
S3=sin(x3); 
k03=[C3*C3 C3*S3 -C3*C3 -C3*S3; C3*S3 S3*S3 -C3*S3 -S3*S3; 
 -C3*C3 -C3*S3 C3*C3 C3*S3; -C3*S3 -S3*S3 C3*S3 S3*S3];
% x4=pi/4; 
C4=cos(x4); 
S4=sin(x4); 
k04=[C4*C4 C4*S4 -C4*C4 -C4*S4; C4*S4 S4*S4 -C4*S4 -S4*S4; 
 -C4*C4 -C4*S4 C4*C4 C4*S4; -C4*S4 -S4*S4 C4*S4 S4*S4];
K1=repmat(E0*A0/l*k01,1,nely);
K2=repmat(E0*A0/(sqrt(2)*l/2)*k03,1,nely);
K3=repmat(E0*A0/(sqrt(2)*l/2)*k04,1,nely);
ssk=[repmat([K1,K2,K2,K3,K3],1,nelx),K1];
e_0=0:2*nely+1:(2*nely+1)*(nelx-1);
e1=[2*(1:nely)'-1,2*(1:nely)',2*(1:nely)'+1,2*(1:nely)'+2];
ee2_1=[2*(1:nely)'-1,2*(1:nely)',2*(1:nely)'+2*nely+1,2*(1:nely)'+2*nely+2];
ee2_2=[2*(1:nely)'+2*nely+1,2*(1:nely)'+2*nely+2,2*(1:nely)'+4*nely+3,2*(1:nely)'+4*nely+4];
ee3_1=[2*(1:nely)'+1,2*(1:nely)'+2,2*(1:nely)'+2*nely+1,2*(1:nely)'+2*nely+2];
ee3_2=[2*(1:nely)'+2*nely+1,2*(1:nely)'+2*nely+2,2*(1:nely)'+4*nely+1,2*(1:nely)'+4*nely+2];
eee_0=[e1;ee2_1;ee2_2;ee3_1;ee3_2];
edofMat=[repmat(eee_0,nelx,1)+2*reshape(repmat(e_0,5*nely,1),5*nely*nelx,1);e1+repmat(2*(2*nely+1)*nelx,nely,1)];
