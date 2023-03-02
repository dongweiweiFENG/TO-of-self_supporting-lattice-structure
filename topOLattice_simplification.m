%[A,xPhys,L,V,E,cList,vList,sList]=topOLattice_simplification(128,64,1200,0.4);
function [A,xPhys,L,V,E,cList,vList,sList]=topOLattice_simplification(nelx,nely,ItMax,volfrac)
bBlackWhite=1;
E0=1;
folder = sprintf('images-%8.4f',volfrac);
mkdir(folder);
cList = zeros(ItMax,1);
vList = zeros(ItMax,1);
sList = zeros(ItMax,1);

load 2d_subdivision;% input the subdivision data
A0=pi;l=3*ceil(sqrt(A0/pi));
num=length(xPhys);
[ssk,edofMat]=Bar2D_Stiffness1_2_new(-pi/2,-pi/4,pi/4,nelx,nely,E0,A0,l);
iK = reshape(kron(edofMat,ones(4,1))',16*(5*nely*nelx+nely),1);
jK = reshape(kron(edofMat,ones(1,4))',16*(5*nely*nelx+nely),1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2*((2*nely+1)*nelx+nely+1),1,-1,2*((2*nely+1)*nelx+nely+1),1); %-nely
U = zeros(2*((2*nely+1)*nelx+nely+1),1);
fixeddofs = union([1:1:2*(nely+1)],[1]);  
alldofs = [1:2*((2*nely+1)*nelx+nely+1)];  
freedofs = setdiff(alldofs,fixeddofs);

fixed_E=find(xPhys<=0.1); % avoid the deleted struts in subdivision to participate in simplification
A(fixed_E)=0;
v0=sum(L*A0);

Nx_down=A_nodes_down(num,nelx,nely); % the relation matrix about the nodes and lower struts
Nx_up=A_nodes_up(num,nelx,nely); % the relation matrix about the nodes and upper struts
indx1=1:2*nely+1:(2*nely+1)*nelx+1;
indx2=nely+1:2*nely+1:(2*nely+1)*nelx+nely+1;
indxx=union(indx1,indx2);
Nx_down(indxx,:)=[];
Nx_up(indxx,:)=[];
L=L';
loop = 0;
loopbeta = 0;
p=3; 
xPhys=0.3*ones(num,1);
xPhys(fixed_E)=0;
alpha=600;
x_min=L.^2/alpha;
x_min=(max(0.16,x_min))';
x_min(fixed_E)=[];
x_min=x_min';
%% size filter
bb=100000;aa=x_min;
if bBlackWhite
    xTilde=xPhys;
    xTilde(fixed_E)=[];
    ind=setdiff(1:num,fixed_E);
    [xTilde,~]=scaling_fun(xTilde,aa,bb);
    xPhys(ind)=xTilde;
end

xx = xPhys(:);
xx(fixed_E) = [];
xold1 = xx;
xold2 = xx;
low = 0;
upp = 0;
while loop < ItMax
    loop = loop + 1;
    loopbeta = loopbeta + 1;
    A=A0*xPhys;
    
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    v=sum(A.*L);
    [Shx,dSx]=support_constriant(Nx_up,Nx_down,xPhys);
    if loop==1
        s_0=floor(10000*Shx)/10000; % can slightly change the value of s_0 to keep self-supporting
        s_1=Shx;
    end
    [c,dcx]=K_2D_p(ssk,iK,jK,num,F,freedofs,xPhys,edofMat,p,U);
    dvx=L*A0/(v0*volfrac);
    
    cList(loop) = c;
    vList(loop) = v/v0;
    sList(loop) = Shx;
    dcx(fixed_E) = [];
    dvx(fixed_E) = [];
    dSx(fixed_E) = [];
    %% heaviside
    if bBlackWhite
        [~,dx]=scaling_fun(xTilde,aa,bb);
        dcx(:) = dcx(:).*dx(:);
        dvx(:) = dvx(:).*dx(:);
        dSx(:)=dSx(:).*dx(:);
    end
    
    %% UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    % METHOD OF MOVING ASYMPTOTES (MMA)
    n = num-length(fixed_E);
    xval = xx;
    move = 0.01;
    xmin = max(xx-move, 0);
    xmax = min(xx+move, 1);
    %% constraints
    WScale = 1e-2;
    f0val = c*WScale;
    if (f0val > 100 || f0val < 1)% && maxIter > 1
        WScale = 10/c;
        f0val = c*WScale;
    end
    df0dx =dcx'*WScale;
    fval = [];
    dfdx = [];
    fval = [fval; v/(v0*volfrac)-1;Shx-s_0];
    dfdx = [dfdx,dvx,dSx]';
    m = length(fval);
    mdof = 1:m;
    
    %% solving
    a0 = 1;
    a = zeros(m,1);
    c_ = ones(m,1)*1000;
    d = zeros(m,1);
    iter = loop;
    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp] = ...
        mmasub(m, n, iter, xval, xmin, xmax, xold1, xold2,...
        f0val, df0dx, fval(mdof), dfdx(mdof,:),low, upp, a0, a, c_, d);

    xold2 = xold1;
    xold1 = xval;
    xx = xmma;
    xTilde(:) = xx(:);
    if bBlackWhite
        xT=scaling_fun(xTilde,aa,bb);
    else
        xT=xTilde;
    end
    xPhys = zeros(num, 1);
    index = setdiff(1:num, fixed_E);
    xPhys(1:nely)=1;
    xPhys(index) = xT;
    
    fprintf(' It.:%5i Obj.:%11.4f  vol.:%11.4f, Shx.:%11.4f\n',loop,c, ...
        v/v0,max(Shx));
        
    if Shx<=s_1&&loop>=1000&&min(xT(xT>0.01)-x_min(xT>0.01))>=0%&&change<=2e-3
        save mata_128_600
        clear VV FF;
        figure;
        xPhys(xPhys<1e-2)=0;
        D=2*sqrt(A0/pi*xPhys);
        drawframe(V,E,D,'k');
        break
    end
end
end

function [Shx,dSx]=support_constriant(Nx_up,Nx_down,A)
    p1=100;
    x_up=Nx_up.*exp(p1*A');
    x_down=Nx_down.*exp(p1*A');
    sx_upp=sum(x_up,2);
    sx_downn=sum(x_down,2);
    Sx_up=1/p1*log(sx_upp); 
    Sx_down=1/p1*log(sx_downn);
    p=100;
    SHx=Sx_up-Sx_down;
    SH=exp(1).^SHx;
    Shx=(sum(SH.^p))^(1/p);
    dSx=x_up'*(1./sx_upp.*((sum(SH.^p))^(1/p-1)*SH.^(p)))-x_down'*(1./sx_downn.*((sum(SH.^p))^(1/p-1)*SH.^(p)));
end

function Nx_down=A_nodes_down(num,nelx,nely)
num1=(2*nely+1)*nelx+nely+1;
iter = 1;
dim=num-3*(nelx-1)-4; %考虑悬空的杆数
iM = zeros(dim,1);
jM = zeros(dim,1);
vM = ones(dim,1);
i1=(2*nely+1)*nelx;
i2=5*nely*nelx-nely;
for i=2:nely %左右两边界
    it = i; % index of the element
    t=i;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it = i+nely; % index of the element
    t=i;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it = i+i2; % index of the element
    t=i+i1;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it = i+nely+i2; % index of the element
    t=i+i1;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
end
for j=1:nelx-1
    j1=j*(2*nely+1);
    for i=2:nely
        t=i+j1;
        
        it = (j-1)*5*nely+4*nely+i; % index of the element
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it = (j-1)*5*nely+5*nely+i; % index of the element
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it = (j-1)*5*nely+6*nely+i; % index of the element
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
    end
end
for j=1:nelx
    j1=(j-1)*(2*nely+1);
    for i=1:nely
        t=i+j1+nely+1;
        
        it = (j-1)*5*nely+3*nely+i; % index of the element
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it = (j-1)*5*nely+2*nely+i; % index of the element
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
    end
end
    Nx_down = sparse(jM',iM',vM',num1,num);  
end

function Nx_up=A_nodes_up(num,nelx,nely)
num1=(2*nely+1)*nelx+nely+1;
iter = 1;
dim=num-3*(nelx-1)-4;
iM = zeros(dim,1);
jM = zeros(dim,1);
vM = ones(dim,1);
i1=(2*nely+1)*nelx;
i2=5*nely*nelx-3*nely;
for i=2:nely
    it = i-1; 
    t=i;% index of the element
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it = i-1+3*nely; 
    t=i;% index of the element
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it = i-1+i2; 
    t=i+i1;% index of the element
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it = i-1+3*nely+i2; 
    t=i+i1;% index of the element
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
end
for j=1:nelx-1
    j1=j*(2*nely+1);
    for i=2:nely
        t=i+j1;
        
        it = (j-1)*5*nely+2*nely+i-1; % index of the element
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it = (j-1)*5*nely+5*nely+i-1; % index of the element
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it = (j-1)*5*nely+8*nely+i-1; % index of the element
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
    end
end
for j=1:nelx
    j1=(j-1)*(2*nely+1);
    for i=1:nely
        t=i+j1+nely+1;
        
        it = (j-1)*5*nely+nely+i; % index of the element
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it = (j-1)*5*nely+4*nely+i; % index of the element
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
    end
end
    Nx_up = sparse(jM',iM',vM',num1,num);  
end

