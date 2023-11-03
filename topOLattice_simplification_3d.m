%[V,E,A,L,xPhys]=topOLattice_simplification_3d(32,16,16,0.06,1,1000);
%====output=======
% V -- the coordinate of the nodes
% E -- the information about the nodes connected to the struts and the length of the struts
% A -- the vector about cross-section of struts
% L -- the vector about the length of struts
% xPhys -- the vector about the density of struts

%====input=======
% nelx x nely x nelz -- the resolution about the ground lattice structure
% volfrac -- the volume fraction
% ItMax -- the total iteration steps
% bBlackWhite -- value of 1 means using size filter
function [V,E,A,L,xPhys]=topOLattice_simplification_3d(nelx,nely,nelz,volfrac,bBlackWhite,ItMax)
E0=1;
load 3d_subdivision
num=nelz*((13*nely+3)*nelx+3*nely+1);
v0=A0*sum(L);

num_dofs=3*(nely+1)*(nelx+1)*(nelz+1)+3*nely*nelx*nelz;
U = zeros(num_dofs,1);

fix_A_0=find(xPhys<=1e-4);
fixed_E=fix_A_0; %do not consider the struts that has been optimized

Nx_down=A_nodes_down(num,nelx,nely,nelz);
Nx_up=A_nodes_up(num,nelx,nely,nelz);
indx1=1:(nely+1)*(nelx+1);
indx2=[(nely+1)*(nelx+1)+nely*nelx]*nelz+1:(nely+1)*(nelx+1)*(nelz+1)+nely*nelx*nelz;
indxx=union(indx1,indx2);
Nx_up(indxx,:)=[];
Nx_down(indxx,:)=[];
L=L';
loop = 0;
loopbeta = 0;
p=4;  
xPhys(:)=0.3; % initialize the density of every struts
xPhys(fix_A_0)=0;

alpha=600;
x_min=L.^2/alpha;
x_min=max(0.09,x_min);
x_min(fixed_E)=[];
%% size_filter
bb=1;aa=x_min;
if bBlackWhite    
    xTilde=xPhys;
    xTilde(fixed_E)=[];
    ind=setdiff(1:num,fixed_E);
    [xTilde,~]=size_filter(xTilde,aa,bb);
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
    if loop==1
        ssk0=ssk;
    end
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS 
    [c,dcx]=K_3D_p(ssk,num,F,freedofs,xPhys,edofMat,p,U);
    v=sum(A.*L); 
    [Shx,dSx]=support_constriant(Nx_up,Nx_down,xPhys);
    if loop==1
        s_0=floor(10000*Shx-2)/10000;
        s_1=roundn(Shx,-5)+1e-5;
    end
    dvx=L*A0/(v0*volfrac);
    dcx(fixed_E) = [];   
    dvx(fixed_E) = [];
    dSx(fixed_E) = [];
    %% Size filter
    if bBlackWhite
        [~,dx]=size_filter(xTilde,aa,bb);
        dcx(:) = dcx(:).*dx(:);
        dvx(:) = dvx(:).*dx(:);
        dSx(:)=dSx(:).*dx(:);
    end
    %%
    %% UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    % METHOD OF MOVING ASYMPTOTES (MMA)
    n = num-length(unique(fixed_E));
    xval = xx;
    move = 0.01;
    xmin = max(xx-move, 0);
    xmax = min(xx+move, 4.41);
    WScale = 1e-2;
    f0val = c*WScale;
    if (f0val > 100 || f0val < 1)
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
    change = max(abs(xmma(:)-xx(:)));
    xx = xmma;
    xTilde(:) = xx(:);
    if bBlackWhite
        xT=size_filter(xTilde,aa,bb);
    else
        xT=xx;
    end
    
    xPhys = zeros(num, 1);
    index = setdiff(1:num, fixed_E);
    xPhys(index) = xT;
    
    fprintf(' It.:%5i Obj.:%11.4f  V.:%11.4f,  Shx.:%11.4f\n',loop,c, ...
        v/v0,max(Shx));
    
    if loop>=800&&Shx<=s_1
        bb=100000;
    end
    if Shx<=s_1&&loop>=800&&min(xT(xT>0.09)-x_min(xT>0.09))>=0
        clear VV FF;
        xPhys(xPhys<1e-2)=0;
        r=2/3*sqrt(A0/pi*xPhys);
        [VV,FF]=draw_frame(V,E,r);
        saveOff(VV,FF,'3d_simplification.off');
        break
    end
end
end

function [Shx,dSx]=support_constriant(Nx_up,Nx_down,A)
    p1=120;
    x_up=Nx_up.*exp(p1*A');
    x_down=Nx_down.*exp(p1*A');
    sx_upp=sum(x_up,2);
    sx_downn=sum(x_down,2);
    Sx_up=1/p1*log(sx_upp); 
    Sx_down=1/p1*log(sx_downn);
    p=120;
    SHx=Sx_up-Sx_down;
    SH=exp(1).^SHx;
    Shx=(sum(SH.^p))^(1/p);
    dSx=x_up'*(1./sx_upp.*((sum(SH.^p))^(1/p-1)*SH.^(p)))-x_down'*(1./sx_downn.*((sum(SH.^p))^(1/p-1)*SH.^(p)));
end

% the matrix that considering the nodes and the lower struts that connected
function Nx_down=A_nodes_down(num,nelx,nely,nelz)
num0=(nely+1)*(nelx+1)*(nelz+1)+nelx*nely*nelz;
num1=13*nely+3;
num2=num1*nelx+3*nely+1;
iter = 1;
dim=num-num2+4*nely*nelx;
iM = zeros(dim,1);
jM = zeros(dim,1);
vM = ones(dim,1);
i1=(nely+1)*(nelx+1)+nelx*nely;

for k=2:nelz
    for j=2:nely
        for i=2:nelx
            t=(i-1)*(nely+1)+j+(k-1)*i1;
            it=3*nely+j+1+(i-2)*num1+(k-2)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;

            it=6*nely+j+2+(i-2)*num1+(k-2)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=10*nely+j+3+(i-2)*num1+(k-2)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=j+(i-1)*num1+(k-2)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=nely+j+(i-1)*num1+(k-2)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=2*nely+j+1+(i-1)*num1+(k-2)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=4*nely+2+j+(i-1)*num1+(k-2)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=7*nely+3+j+(i-1)*num1+(k-2)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=11*nely+2+j+(i-1)*num1+(k-2)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
        end
    end
    
    j=1;
    i=1;
    t=(i-1)*(nely+1)+j+(k-1)*i1;
    it=j+(i-1)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=2*nely+1+j+(i-1)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=4*nely+2+j+(i-1)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=7*nely+3+j+(i-1)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    for i=2:nelx
        t=(i-1)*(nely+1)+j+(k-1)*i1;
        it=3*nely+j+1+(i-2)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=10*nely+j+3+(i-2)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=2*nely+j+1+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=4*nely+2+j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=7*nely+3+j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
    end
    i=nelx+1;
    t=(i-1)*(nely+1)+j+(k-1)*i1;
    it=3*nely+j+1+(i-2)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=10*nely+j+3+(i-2)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=j+(i-1)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=2*nely+j+1+(i-1)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    j=nely+1;
    i=1;
    t=(i-1)*(nely+1)+j+(k-1)*i1;
    it=j+(i-1)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=nely+j+(i-1)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=4*nely+2+j+(i-1)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=11*nely+2+j+(i-1)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    for i=2:nelx
        t=(i-1)*(nely+1)+j+(k-1)*i1;
        it=3*nely+j+1+(i-2)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=6*nely+j+2+(i-2)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=nely+j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
                
        it=4*nely+2+j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=11*nely+2+j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
    end
    i=nelx+1;
    t=(i-1)*(nely+1)+j+(k-1)*i1;
    it=3*nely+j+1+(i-2)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=6*nely+j+2+(i-2)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=j+(i-1)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=nely+j+(i-1)*num1+(k-2)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    i=1;
    for j=2:nely
        t=(i-1)*(nely+1)+j+(k-1)*i1;
        it=j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=nely+j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=2*nely+1+j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=4*nely+2+j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=7*nely+3+j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=11*nely+2+j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
    end
    
    i=nelx+1;
    for j=2:nely
        t=(i-1)*(nely+1)+j+(k-1)*i1;
        it=3*nely+j+1+(i-2)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=6*nely+j+2+(i-2)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=10*nely+j+3+(i-2)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=nely+j+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=2*nely+j+1+(i-1)*num1+(k-2)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
    end
    
    
end
for k=1:nelz
    for i=1:nelx
        for j=1:nely
            t=(i-1)*nely+j+(nely+1)*(nelx+1)+i1*(k-1);
            it=5*nely+j+3+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=9*nely+j+3+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=8*nely+j+3+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=12*nely+3+j+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
        end
    end
end
    Nx_down = sparse(jM',iM',vM',num0,num);  
end

% the matrix that considering the nodes and the upper struts that connected
function Nx_up=A_nodes_up(num,nelx,nely,nelz)
num0=(nely+1)*(nelx+1)*(nelz+1)+nelx*nely*nelz;
num1=13*nely+3;
num2=num1*nelx+3*nely+1;
iter = 1;
dim=num-num2+4*nely*nelx;
iM = zeros(dim,1);
jM = zeros(dim,1);
vM = ones(dim,1);
i1=(nely+1)*(nelx+1)+nelx*nely;
for k=2:nelz  
    for j=2:nely
        for i=2:nelx
            t=(i-1)*(nely+1)+j+(k-1)*i1;
            it=4*nely+j+2+(i-2)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=8*nely+j+2+(i-2)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=12*nely+j+3+(i-2)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=j+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=nely+1+j+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=2*nely+j+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=3*nely+1+j+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=5*nely+3+j+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=9*nely+2+j+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;   
        end
    end
    j=1;
    i=1;
    t=(i-1)*(nely+1)+j+(k-1)*i1;
    it=j+(i-1)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=nely+1+j+(i-1)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=3*nely+1+j+(i-1)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=5*nely+3+j+(i-1)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    for i=2:nelx
        t=(i-1)*(nely+1)+j+(k-1)*i1;
        it=4*nely+j+2+(i-2)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=12*nely+j+3+(i-2)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=nely+1+j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=3*nely+1+j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=5*nely+3+j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
    end
    i=nelx+1;
    t=(i-1)*(nely+1)+j+(k-1)*i1;
    it=4*nely+j+2+(i-2)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=12*nely+j+3+(i-2)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=j+(i-1)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=nely+1+j+(i-1)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;%j=1;
    
    j=nely+1;
    i=1;
    t=(i-1)*(nely+1)+j+(k-1)*i1;
    it=j+(i-1)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=2*nely+j+(i-1)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=3*nely+1+j+(i-1)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=9*nely+2+j+(i-1)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    for i=2:nelx
        t=(i-1)*(nely+1)+j+(k-1)*i1;
        it=4*nely+j+2+(i-2)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=8*nely+j+2+(i-2)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=2*nely+j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=3*nely+1+j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=9*nely+2+j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
    end
    i=nelx+1;
    t=(i-1)*(nely+1)+j+(k-1)*i1;
    it=4*nely+j+2+(i-2)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=8*nely+j+2+(i-2)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=j+(i-1)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    it=2*nely+j+(i-1)*num1+(k-1)*num2;
    iM(iter) = it;
    jM(iter) = t;
    iter = iter + 1;
    
    i=1;
    for j=2:nely
        t=(i-1)*(nely+1)+j+(k-1)*i1;
        it=j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=nely+1+j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=2*nely+j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=3*nely+1+j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=5*nely+3+j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
        
        it=9*nely+2+j+(i-1)*num1+(k-1)*num2;
        iM(iter) = it;
        jM(iter) = t;
        iter = iter + 1;
    end
    i=nelx+1;
    for j=2:nely
       t=(i-1)*(nely+1)+j+(k-1)*i1;
       it=4*nely+j+2+(i-2)*num1+(k-1)*num2;
       iM(iter) = it;
       jM(iter) = t;
       iter = iter + 1;
       
       it=8*nely+j+2+(i-2)*num1+(k-1)*num2;
       iM(iter) = it;
       jM(iter) = t;
       iter = iter + 1;
       
       it=12*nely+j+3+(i-2)*num1+(k-1)*num2;
       iM(iter) = it;
       jM(iter) = t;
       iter = iter + 1;
       
       it=j+(i-1)*num1+(k-1)*num2;
       iM(iter) = it;
       jM(iter) = t;
       iter = iter + 1;
       
       it=nely+1+j+(i-1)*num1+(k-1)*num2;
       iM(iter) = it;
       jM(iter) = t;
       iter = iter + 1;
       
       it=2*nely+j+(i-1)*num1+(k-1)*num2;
       iM(iter) = it;
       jM(iter) = t;
       iter = iter + 1;
    end
   
end
for k=1:nelz
    for i=1:nelx
        for j=1:nely
            t=(i-1)*nely+j+(k-1)*i1+(nely+1)*(nelx+1);
            it=6*nely+j+3+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=7*nely+j+3+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=10*nely+j+3+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
            
            it=11*nely+3+j+(i-1)*num1+(k-1)*num2;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
        end
    end
end
    Nx_up = sparse(jM',iM',vM',num0,num);  
end

%size filter
function [x,dx]=size_filter(xPhys,aa,bb)
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
