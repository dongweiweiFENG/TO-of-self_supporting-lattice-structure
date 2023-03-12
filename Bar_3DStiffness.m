% the stiffness matrix of the struts and the information about the nodes and struts(edofMat)
% ssk -- the stiffness matrix of the struts
% edofMat -- the information about the nodes and struts

% nelx x nely x nelz is the resolution of the ground lattice structure
% E0 is the Young's modulus of the struts
% A0 is the initial cross-section of the struts
% l is the length of the minimum struts in the virtical direction
function [ssk,edofMat]=Bar_3DStiffness(nelx,nely,nelz,E0,A0,l) 
k0=E0*A0/l*[1,-1;-1,1]; 
t1=[cos(pi/2),cos(pi/2),-1,0,0,0;0,0,0,cos(pi/2),cos(pi/2),-1];
t2=[0,-1/sqrt(3),sqrt(6)/3,0,0,0;0,0,0,0,-1/sqrt(3),sqrt(6)/3];
t3=[0,-1/sqrt(3),-sqrt(6)/3,0,0,0;0,0,0,0,-1/sqrt(3),-sqrt(6)/3];
t4=[1/sqrt(3),0,sqrt(6)/3,0,0,0;0,0,0,1/sqrt(3),0,sqrt(6)/3];
t5=[1/sqrt(3),0,-sqrt(6)/3,0,0,0;0,0,0,1/sqrt(3),0,-sqrt(6)/3];
t9=[1,1,-sqrt(2),0,0,0;0,0,0,1,1,-sqrt(2)]/2;
t8=[1,1,sqrt(2),0,0,0;0,0,0,1,1,sqrt(2)]/2;
t7=[1,-1,-sqrt(2),0,0,0;0,0,0,1,-1,-sqrt(2)]/2;
t6=[1,-1,sqrt(2),0,0,0;0,0,0,1,-1,sqrt(2)]/2;
k1=t1'*k0*t1;
k2=t2'*k0*t2/sqrt(6)*2;
k3=t3'*k0*t3/sqrt(6)*2;
k4=t4'*k0*t4/sqrt(6)*2;
k5=t5'*k0*t5/sqrt(6)*2;
k6=t6'*k0*t6/sqrt(2)*2;
k7=t7'*k0*t7/sqrt(2)*2;
k8=t8'*k0*t8/sqrt(2)*2;
k9=t9'*k0*t9/sqrt(2)*2;

K1=repmat(k1,1,nely+1);
K2=repmat(k2,1,nely);
K3=repmat(k3,1,nely);
K4=repmat(k4,1,nely+1);
K5=repmat(k5,1,nely+1);
K6=repmat(k6,1,2*nely);
K7=repmat(k7,1,2*nely);
K8=repmat(k8,1,2*nely);
K9=repmat(k9,1,2*nely);
sk=[repmat([K1,K2,K3,K4,K5,K6,K7,K8,K9],1,nelx),K1,K2,K3];
ssk=repmat(sk,1,nelz);
num0=13*nelx*nely+3*nelx;
num1=3*(nelx+1)*(nely+1)+3*nelx*nely;
num1_1=3*(nelx+1)*(nely+1);
num2=13*nelx*nely+3*nelx+3*nely+1;
e=(0:num1:num1*(nelz-1))';
e1=[3*(1:nely+1)'-2,3*(1:nely+1)'-1,3*(1:nely+1)',3*(1:nely+1)'+num1-2,3*(1:nely+1)'+num1-1,3*(1:nely+1)'+num1];
e2=[3*(1:nely)'-2,3*(1:nely)'-1,3*(1:nely)',3*(1:nely)'+num1+1,3*(1:nely)'+num1+2,3*(1:nely)'+num1+3];

e3=[3*(1:nely)'+num1-2,3*(1:nely)'+num1-1,3*(1:nely)'+num1,3*(1:nely)'+1,3*(1:nely)'+2,3*(1:nely)'+3];

e4=[3*(1:nely+1)'-2,3*(1:nely+1)'-1,3*(1:nely+1)',3*(1:nely+1)'+num1+3*nely+1,3*(1:nely+1)'+num1+3*nely+2,3*(1:nely+1)'+num1+3*nely+3];

e5=[3*(1:nely+1)'+num1-2,3*(1:nely+1)'+num1-1,3*(1:nely+1)'+num1,3*(1:nely+1)'+3*nely+1,3*(1:nely+1)'+3*nely+2,3*(1:nely+1)'+3*nely+3];


e6=[3*(1:nely)'-2,3*(1:nely)'-1,3*(1:nely)',3*(1:nely)'+num1_1-2,3*(1:nely)'+num1_1-1,3*(1:nely)'+num1_1];
e66=[3*(1:nely)'+num1_1-2,3*(1:nely)'+num1_1-1,3*(1:nely)'+num1_1,3*(1:nely)'+num1+3*nely+4,3*(1:nely)'+num1+3*nely+5,3*(1:nely)'+num1+3*nely+6];

e7=[3*(1:nely)'+num1-2,3*(1:nely)'+num1-1,3*(1:nely)'+num1,3*(1:nely)'+num1_1-2,3*(1:nely)'+num1_1-1,3*(1:nely)'+num1_1];
e77=[3*(1:nely)'+num1_1-2,3*(1:nely)'+num1_1-1,3*(1:nely)'+num1_1,3*(1:nely)'+3*nely+4,3*(1:nely)'+3*nely+5,3*(1:nely)'+3*nely+6];

e8=[3*(1:nely)'+1,3*(1:nely)'+2,3*(1:nely)'+3,3*(1:nely)'+num1_1-2,3*(1:nely)'+num1_1-1,3*(1:nely)'+num1_1];
e88=[3*(1:nely)'+num1_1-2,3*(1:nely)'+num1_1-1,3*(1:nely)'+num1_1,3*(1:nely)'+num1+3*nely+1,3*(1:nely)'+num1+3*nely+2,3*(1:nely)'+num1+3*nely+3];
e9=[3*(1:nely)'+num1+1,3*(1:nely)'+num1+2,3*(1:nely)'+num1+3,3*(1:nely)'+num1_1-2,3*(1:nely)'+num1_1-1,3*(1:nely)'+num1_1];
e99=[3*(1:nely)'+num1_1-2,3*(1:nely)'+num1_1-1,3*(1:nely)'+num1_1,3*(1:nely)'+3*nely+1,3*(1:nely)'+3*nely+2,3*(1:nely)'+3*nely+3];

ee0=[e1;e2;e3;e4;e5;e6;e66;e7;e77;e8;e88;e9;e99];
eee1=[e1+3*(nely+1)*nelx;e2+3*(nely+1)*nelx;...
    e3+3*(nely+1)*nelx];
e0=[repmat(nely+1,5*nely+3,1);...
    repmat([repmat(nely+1,nely,1);repmat(nely,nely,1)],4,1)];
e00=[repmat(nely+1,5*nely+3,1);...
    repmat([repmat(nely,nely,1);repmat(nely+1,nely,1)],4,1)];
e_0=[repmat(e0,nelx,3),repmat(e00,nelx,3)].*3.*reshape(repmat(0:nelx-1,13*nely+3,1),num0,1);
edofMat1=[repmat(ee0,nelx,1)+e_0;eee1];
edofMat=repmat(edofMat1,nelz,1)+reshape(repmat(reshape(repmat(e,1,6),1,6*nelz),num2,1),num2*nelz,6);