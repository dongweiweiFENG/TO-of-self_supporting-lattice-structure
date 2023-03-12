% [V,E,xPhys,L]=topOLattice_subdivision_3d(32,16,16,0.5,1,400)

% V -- the coordinate of the nodes
% E -- the information about the nodes connected to the struts and the length of the struts
% A -- the vector about cross-section of struts
% L -- the vector about the length of struts
% xPhys -- the vector about the density of struts
% cList -- the iterative compliance information
% vList -- the iterative volume information
% nelx x nely x nelz -- the resolution about the ground lattice structure
% volfraction -- the volume fraction
% bBlackWite -- 0 for no heaviside projection and 1 for heaviside projection
% ItMax -- the total iteration steps
function [V,E,xPhys,L]=topOLattice_subdivision_3d(nelx,nely,nelz,volfrac,bBlackWhite,ItMax)
clc;
% close all;
bSubdivision=2;
MMA = 1;
self_s=1; % self_supporting filter
folder = sprintf('images-%8.4f',volfrac);
mkdir(folder);
cList = zeros(ItMax,1);
vList = zeros(ItMax,1);
sList = zeros(ItMax,1);

penal = 3;
pNorm = 16;

A0=pi;
l=6;
%storage the length of the struts
num0=13*nelx*nely+3*nelx+3*nely+1;
L0=ones(1,num0)*l;
for i=1:4*nely+2
    L0(1,nely+1+i:13*nely+3:end)=l*sqrt(6)/2;
end
for i=1:8*nely
    L0(1,5*nely+3+i:13*nely+3:end)=l*sqrt(2)/2;
end
L=repmat(L0,1,nelz);
v0=sum(L.*A0);
clear L0;
%coordinate V
X = repmat(repmat(0:l:nelx*l, (nely+1), 1),1,nelz+1)/sqrt(2);
X1 = repmat(repmat(l/2:l:nelx*l, nely, 1),1,nelz)/sqrt(2);
Y = repmat([nely*l:-l:0]', 1, (nelx+1)*(nelz+1))/sqrt(2);
Y1 = repmat([nely*l-l/2:-l:0]', 1, nelx*nelz)/sqrt(2);
Z=repmat(0:l:nelz*l,(nelx+1)*(nely+1),1);
Z1=repmat(l/2:l:nelz*l,nelx*nely,1);
V0 = [X(:),Y(:),Z(:)];
V1 = [X1(:),Y1(:),Z1(:)];
V=zeros(length(V0)+length(V1),3);
num_V_=(nely+1)*(nelx+1)+nely*nelx;
for i=1:length(V)
    j=floor((i-1)/(num_V_));
    if i-j*num_V_<=(nely+1)*(nelx+1)
        V(i,:)=V0(j*(nely+1)*(nelx+1)+i-j*num_V_,:);
    else
        V(i,:)=V1(j*(nely*nelx)+i-j*num_V_-(nely+1)*(nelx+1),:);
    end
end
E0 = 1;
[ssk,edofMat]=Bar_3DStiffness(nelx,nely,nelz,E0,A0,l);
e=edofMat(:,[3,6])/3;
E=[e,L']; %the struts information about its connected nodes and  length

nLevels = min(log2(nelx),min(log2(nely),log2(nelz)));  %k_bar
% Cell dimension
nc0 = 2^nLevels;   % Coarset cell, e.g., nc0*nc0 = 32*32  --2^m
nc = zeros(nLevels,1); 
for level = 1:nLevels   %k
    nc(level) = nc0/(2^(level-1));
end

% Elements within cell
nx0 = nelx/nc0;
ny0 = nely/nc0;
nz0=nelz/nc0;
nx = zeros(nLevels,1);
ny = zeros(nLevels,1); 
nz = zeros(nLevels,1); 
for level = 1:nLevels
    nx(level) = nx0 * 2^(level-1);
    ny(level) = ny0 * 2^(level-1); 
    nz(level) = nz0 * 2^(level-1); 
end

h_min = 0;

% Design variables
h0 = ones(nx0*ny0*nz0,1);
h = {zeros(nx(1)*ny(1)*nz(1),1)};
for level = 2:nLevels
    h{level} = zeros(nx(level)*ny(level)*nz(level),1);
end

% A copy of design variables, for caculating the initial volume
h1 = {ones(nx(1)*ny(1)*nz(1),1)};
for level = 2:nLevels
    h1{level} = ones(nx(level)*ny(level)*nz(level),1);
end

% Design variables projected by the subdivision rule
hp = h; 
% Design variables converted to black-white design
hbw = h;

%%
num_x = 13*nelx*nely*nelz+3*nelx*nelz+3*nely*nelz+nelz;
[M0x] = buildSparseM0x(nx,ny,nz,nc0,nelx,nely,num_x);
[Mx] = buildSparseMx(nLevels,nc,nx,ny,nz,nelx,nely,num_x);  

[xFix] = designToX(M0x,h0,Mx,nLevels,h);

xFixx=A0*L'.*xFix;

initFrac = (volfrac*v0 - sum(xFixx(:))) / (v0 - sum(xFixx(:)));
for level = 1:nLevels
    h{level} = initFrac * h1{level};
end

[hlinear] = linearize(h,nLevels);
[hplinear] = linearize(hp,nLevels);

numDesign = 0; % Number of total design variables
for level = 1:nLevels
    numDesign = numDesign + size(h{level},1);
end

[hde,hdeCo] = initializeDependency(nx,ny,nz,nLevels);

[hdeEx,hdeExiLevels,hdeExNum,hdeExCo] = initializeExtendedDependency(nx,ny,nz,nLevels,hde);

%% subdivisionProjection
[hp,hdeExCo] = restrictedSubdivisionProjectionPNorm(h,hp,h0,hdeEx,hdeExCo,hdeExiLevels,hdeExNum,pNorm,nLevels,nx,ny,nz);

beta = 1;       % beta continuation
eta = 0.5;      % projection threshold, fixed at 0.5

%% Adapted from 88 lines
% %% PREPARE FINITE ELEMENT ANALYSIS
num_node=(nely+1)*(nelx+1)+nely*nelx;
num_nodes=(nely+1)*(nelx+1)*(nelz+1)+nelz*nely*nelx;
iK = reshape(kron(edofMat,ones(6,1))',36*(num_x),1);
jK = reshape(kron(edofMat,ones(1,6))',36*(num_x),1);
F = sparse(3*(nely+1)*nelx+3*[1:nely+1],1,-0.1,3*num_nodes,1);
U = zeros(3*num_nodes,1);
fixednid=(1:nely+1)+[0:num_node:num_nodes]';
fixeddof=[3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2];
fixeddofs=fixeddof(:);
alldofs = [1:3*num_nodes];
freedofs = setdiff(alldofs,fixeddofs);
%% INITIALIZE ITERATION
[x] = designToX(M0x,h0,Mx,nLevels,h);
h01=h;
for i=1:num_x
    if abs(x(i)-2*initFrac)<1e-3
        x(i)=initFrac;
    end
end
xT=x;

if self_s
    xT=self_supporting(xT,nely,nelx,nx,ny,nz,nLevels,nc,pNorm);
end

xPhys = xT;
loop = 0;
loopbeta = 0;
change = 1;
hlinearold1 = hlinear;
hlinearold2 = hlinear;
low = 0;
upp = 0;
p=4;
%% START ITERATION     
while loop < ItMax && change>0.00001
    loop = loop + 1;
    loopbeta = loopbeta + 1;
    %% FE-ANALYSIS
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    [c,dcx]=K_3D_p(ssk,num_x,F,freedofs,xPhys,edofMat,p,U);
    v=A0.*L*xPhys;   
    dvx=L*A0/(v0*volfrac);
    if self_s
        dx_new=self_supporting_dx(xT,nely,nelx,nx,ny,nz,nLevels,nc,pNorm);
        dcx=dcx.*dx_new;
        dvx=dvx.*dx_new;
    end
%     dvx = ones(nely,nelx)/(nelx*nely*volfrac);

    cList(loop) = c;
    vList(loop) = v;
    
	dch = {zeros(1,size(h{1},1))};
    dvh = {zeros(1,size(h{1},1))};
    for level = 2:nLevels
        dch{level} = zeros(1,size(h{level},1));
        dvh{level} = zeros(1,size(h{level},1));
    end
    dcxlinear = reshape(dcx',1,num_x);
    dvxlinear = reshape(dvx',1,num_x);
    for level = 1:nLevels
        dch{level} = dcxlinear*Mx{level};
        dvh{level} = dvxlinear*Mx{level};
    end
    dc = dch{1};
    dv = dvh{1};
    for level = 2:nLevels
        dc = [dc dch{level}];
        dv = [dv dvh{level}];
    end

    if bSubdivision > 0
        dhpdh = zeros(numDesign,numDesign);
        offset = zeros(nLevels,1);
        offset(1,1) = 0;
        for level = 2:nLevels
            offset(level,1) = offset(level-1,1) + size(h{level-1},1);
        end
        
        if bSubdivision == 1
            for level = 1:nLevels
                for i=1:nx(level)
                    for j=1:ny(level)
                        for s=1:nz(level)
                            it = (j-1)*nx(level)+(i-1)+1+(s-1)*nx(level)*ny(level);
                            for ii = 1:(level-1)
                                ic = hde{level}(it,ii);
                                dhpdh(offset(level,1)+it,offset(ii,1)+ic) = hdeCo{level}(it,ii);
                            end
                            dhpdh(offset(level,1)+it,offset(level,1)+it) = hdeCo{level}(it,level);
                        end
                    end
                end
            end
        elseif bSubdivision == 2        
            for level = 1:nLevels
                for i=1:nx(level)
                    for j=1:ny(level)
                        for s=1:nz(level)
                            it = (j-1)*nx(level)+(i-1)+1+(s-1)*nx(level)*ny(level);
                            numEx = hdeExNum{level}(it,1);
                            for ii = 1:numEx
                                ic = hdeEx{level}(it,ii);
                                iLevel = hdeExiLevels{level}(it,ii);
                                dhpdh(offset(level,1)+it,offset(iLevel,1)+ic) = hdeExCo{level}(it,ii);
                            end
                            dhpdh(offset(level,1)+it,offset(level,1)+it) = hdeExCo{level}(it,numEx+1);
                        end
                    end
                end
            end
        end

        dc = dc*dhpdh;
        dv = dv*dhpdh;
    end
    
    if bBlackWhite
        dx = beta * (1-tanh(beta*(hlinear-eta)).*tanh(beta*(hlinear-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
        dc = dc.*dx';
        dv = dv.*dx';
    end
    
    if MMA == 1
        %% MMA
        cScale = 1e-2;
        f0val = c*cScale;
        if (f0val > 100 || f0val < 1)% && maxIter > 1
            cScale = 10/c;
            f0val = c*cScale;
        end

        df0dx = reshape(dc,numDesign,1)*cScale;
        dfdx = reshape(dv,1,numDesign);
        iter = loopbeta;
        xval = hlinear;
        move = 0.01;
        xmin=max(h_min,xval-move);
        xmax=min(1,xval+move);

        fval = v/ (v0*volfrac) - 1;
        
        m = 1;
        n = numDesign; 
        [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
            mmasub(m,n,iter,xval,xmin,xmax,hlinearold1,hlinearold2,...
            f0val,df0dx,fval,dfdx,low,upp,1,0,1000,0);
        hlinearold2 = hlinearold1;
        hlinearold1 = xval;
        
        hnewlinear = xmma;
        offset = 0;
        for level = 1:nLevels
            if level == 1
                hnew = {hnewlinear((offset+1):(offset+size(h{level},1)), 1)};
            else
                hnew{level} = hnewlinear((offset+1):(offset+size(h{level},1)), 1);
            end
            offset = offset + size(h{level},1);
        end
        
        if bBlackWhite
            for level = 1:nLevels
                hbw{level} = (tanh(beta*eta) + tanh(beta*(hnew{level}-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
            end
        else
            hbw = hnew;
        end
        
        [hp,hdeExCo] = restrictedSubdivisionProjectionPNorm(hbw,hp,h0,hdeEx,hdeExCo,hdeExiLevels,hdeExNum,pNorm,nLevels,nx,ny,nz);
       
        [xnew] = designToX(M0x,h0,Mx,nLevels,hp);
        [x_test] = designToX(M0x,h0,Mx,nLevels,h01);
        for i=1:num_x
            if abs(x_test(i)-2*initFrac)<1e-3
                xnew(i)=xnew(i)/2;
            end
        end
        if self_s
            xnew=self_supporting(xnew,nely,nelx,nx,ny,nz,nLevels,nc,pNorm);
        end
    end
    
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    xT=x;    
    xPhys = xT;
    h = hnew;
    [hlinear] = linearize(h,nLevels);
    [hplinear] = linearize(hp,nLevels);

    %% Beta-continuation for Heaviside projection 
    if bBlackWhite
        if beta < 64 && loopbeta == 60
            beta = 2*beta;
            loopbeta = 0;
            change = 1;
            fprintf('Parameter beta increased to %g.\n',beta);
        end
    end
    
    tmp = (x.*(1-x)*4);
    sharp = sum(tmp(:))/num_x;
    sList(loop) = sharp;

    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%11.4f ch.:%7.3f, Sharpness: %7.3f\n',loop,c, ...
        v/ (sum(A0*L)),change,sharp);
    %% PLOT DENSITIES
    if rem(loop,200)==0
        clearvars -except V E L xPhys A A0 l ssk edofMat freedofs alldofs F
        save 3d_subdivision;
        r=sqrt(A0/pi*xPhys);
        [VV,FF]=draw_frame(V,E,r);
        saveOff(VV,FF,'3d_SUBDIVISION.off');
    end
end


%% 
function [hde,hdeCo] = initializeDependency(nx,ny,nz,nLevels)
hde = {zeros(nx(1)*ny(1)*nz(1),1)}; % Design variables' dependence
for level = 2:nLevels % From second level upwards, the dependence to zeroth is neglected.
    hde{level} = zeros(nx(level)*ny(level)*nz(level),level-1);
end
% The first level doesn't depend on the zeroth. The following is added for
% consistency
for i=1:nx(1)
    for j=1:ny(1)
        for s=1:nz(1)
            it = (j-1)*nx(1)+(i-1)+1+(s-1)*nx(1)*ny(1);
            ic = it;
            hde{1}(it) = ic;
        end
    end
end

hdeCo = {zeros(nx(1)*ny(1)*nz(1),1)}; % Design variables' dependence, coefficient of derivative
for level = 2:nLevels % From second level upwards, the dependence to zeroth is neglected.
    hdeCo{level} = zeros(nx(level)*ny(level)*nz(level),level);
end

for level = 2:nLevels
   for i=1:nx(level)
        for j=1:ny(level)
            for s=1:nz(level)
                it = (j-1)*nx(level)+(i-1)+1+(s-1)*nx(level)*ny(level);
                ic = (floor((j+1)/2)-1)*nx(level-1) + (floor((i+1)/2)-1) + 1+(floor((s+1)/2)-1)*nx(level-1)*ny(level-1);
                hde{level}(it,level-1) = ic;
                if level > 2
                    for ii = 1:(level-2)
                        hde{level}(it, ii) = hde{level-1}(ic,ii);
                    end
                end
            end
        end
    end
end

function [hdeEx,hdeExiLevels,hdeExNum,hdeExCo] = initializeExtendedDependency(nx,ny,nz,nLevels,hde)
hdeEx = {zeros(nx(1)*ny(1)*nz(1),1)};
for level = 2:nLevels
    hdeEx{level} = zeros(nx(level)*ny(level)*nz(level),3^(level-1));
end

hdeExiLevels = {zeros(nx(1)*ny(1)*nz(1),1)};
for level = 2:nLevels
    hdeExiLevels{level} = zeros(nx(level)*ny(level)*nz(level),3^(level-1));
end

hdeExNum = {zeros(nx(1)*ny(1)*nz(1),1)};
for level = 2:nLevels
    hdeExNum{level} = zeros(nx(level)*ny(level)*nz(level),1);
end
for level = 2:nLevels
    for i=1:nx(level)
        for j=1:ny(level)
            for s=1:nz(level)
                it = (j-1)*nx(level)+(i-1)+1+(s-1)*nx(level)*ny(level);
                ic = (floor((j+1)/2)-1)*nx(level-1) + (floor((i+1)/2)-1) + 1+(floor((s+1)/2)-1)*nx(level-1)*ny(level-1);
                if rem(i,2) == 0
                    ix = i+1;
                else
                    ix = i-1;
                end
                if rem(j,2) == 0
                    jy = j+1;
                else
                    jy = j-1;
                end
                if rem(s,2)==0
                    sz=s+1;
                else
                    sz=s-1;
                end
                
                icx = 0;
                icy = 0;
                icz=0;
                
                if ix >= 1 && ix <= nx(level)
                    itx = (j-1)*nx(level)+(ix-1)+1+(s-1)*nx(level)*ny(level);
                    icx = hde{level}(itx,level-1);
                end
                if jy >= 1 && jy <= ny(level)
                    ity = (jy-1)*nx(level)+(i-1)+1+(s-1)*nx(level)*ny(level);
                    icy = hde{level}(ity,level-1);
                end
                if sz >= 1 && sz <= nz(level)
                    itz = (j-1)*nx(level)+(i-1)+1+(sz-1)*nx(level)*ny(level);
                    icz = hde{level}(itz,level-1);
                end
                hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
                hdeEx{level}(it,hdeExNum{level}(it,1)) = ic;
                hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = level-1;
                if icx ~= 0
                    hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
                    hdeEx{level}(it,hdeExNum{level}(it,1)) = icx;
                    hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = level-1;
                end
                if icy ~= 0
                    hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
                    hdeEx{level}(it,hdeExNum{level}(it,1)) = icy;
                    hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = level-1;
                end
                if icz ~= 0
                    hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
                    hdeEx{level}(it,hdeExNum{level}(it,1)) = icz;
                    hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = level-1;
                end
                if level > 2
                    for ii = 1:(level-2)
                        for iii = 1:hdeExNum{level-1}(ic,1)
                            bExist = 0;
                            for jj = 1:hdeExNum{level}(it,1)
                                if hdeEx{level}(it,jj) == hdeEx{level-1}(ic,iii)
                                    if hdeExiLevels{level}(it,jj) == hdeExiLevels{level-1}(ic,iii)
                                        bExist = 1;
                                    end
                                end
                            end
                            if bExist == 0
                                hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
                                hdeEx{level}(it,hdeExNum{level}(it,1)) = hdeEx{level-1}(ic,iii);
                                hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = hdeExiLevels{level-1}(ic,iii);
                            end
                        end
                        
                        if icx ~= 0
                            for iii = 1:hdeExNum{level-1}(icx,1)
                                bExist = 0;
                                for jj = 1:hdeExNum{level}(it,1)
                                    if hdeEx{level}(it,jj) == hdeEx{level-1}(icx,iii)
                                        if hdeExiLevels{level}(it,jj) == hdeExiLevels{level-1}(icx,iii)
                                            bExist = 1;
                                        end
                                    end
                                end
                                if bExist == 0
                                    hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
                                    hdeEx{level}(it,hdeExNum{level}(it,1)) = hdeEx{level-1}(icx,iii);
                                    hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = hdeExiLevels{level-1}(icx,iii);
                                end
                            end
                        end
                        
                        if icy ~= 0
                            for iii = 1:hdeExNum{level-1}(icy,1)
                                bExist = 0;
                                for jj = 1:hdeExNum{level}(it,1)
                                    if hdeEx{level}(it,jj) == hdeEx{level-1}(icy,iii)
                                        if hdeExiLevels{level}(it,jj) == hdeExiLevels{level-1}(icy,iii)
                                            bExist = 1;
                                        end
                                    end
                                end
                                if bExist == 0
                                    hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
                                    hdeEx{level}(it,hdeExNum{level}(it,1)) = hdeEx{level-1}(icy,iii);
                                    hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = hdeExiLevels{level-1}(icy,iii);
                                end
                            end
                        end
                        if icz ~= 0
                            for iii = 1:hdeExNum{level-1}(icz,1)
                                bExist = 0;
                                for jj = 1:hdeExNum{level}(it,1)
                                    if hdeEx{level}(it,jj) == hdeEx{level-1}(icz,iii)
                                        if hdeExiLevels{level}(it,jj) == hdeExiLevels{level-1}(icz,iii)
                                            bExist = 1;
                                        end
                                    end
                                end
                                if bExist == 0
                                    hdeExNum{level}(it,1) = hdeExNum{level}(it,1) + 1;
                                    hdeEx{level}(it,hdeExNum{level}(it,1)) = hdeEx{level-1}(icz,iii);
                                    hdeExiLevels{level}(it,hdeExNum{level}(it,1)) = hdeExiLevels{level-1}(icz,iii);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

hdeExCo = {zeros(nx(1)*ny(1)*nz(1),1)}; % Design variables' dependence, coefficient of derivative
for level = 2:nLevels % From second level upwards, the dependence to zeroth is neglected.
    hdeExCo{level} = zeros(nx(level)*ny(level)*nz(level),3^(level-1));
end

function [x] = designToX(M0x,h0,Mx,nLevels,h)
x = M0x*h0;
for level = 1:nLevels
    x = x + Mx{level}*h{level};
end

function [hp,hdeExCo] = restrictedSubdivisionProjectionPNorm(h,hp,h0,hdeEx,hdeExCo,hdeExiLevels,hdeExNum,pNorm,nLevels,nx,ny,nz)
for i=1:nx(1)
    for j=1:ny(1)
        for s=1:nz(1)
            it = (j-1)*nx(1)+(i-1)+1+(s-1)*nx(1)*ny(1);
            ic = it;
            hp{1}(it) = min(h{1}(it), h0(ic));
            hdeExCo{1}(it,1) = 1;
        end
    end
end

for level = 2:nLevels
    de = zeros(3^(level-1),1);
    for i=1:nx(level)
        for j=1:ny(level)
            for s=1:nz(level)
                it = (j-1)*nx(level)+(i-1)+1+(s-1)*nx(level)*ny(level);
                numEx = hdeExNum{level}(it,1);
                de(numEx+1) = h{level}(it);
                for ii = 1:numEx
                    ic = hdeEx{level}(it,ii);
                    iLevel = hdeExiLevels{level}(it,ii);
                    de(ii) = h{iLevel}(ic);
                end
                hp{level}(it) = 0;
                for ii = 1:(numEx+1)
                    hp{level}(it) = hp{level}(it) + (de(ii))^(-pNorm);
                end
                
                A = hp{level}(it) / (numEx+1);
                hp{level}(it) = (hp{level}(it) / (numEx+1))^(-1/pNorm);
                
                for ii = 1:(numEx+1)
                    hdeExCo{level}(it,ii) = A^((-1/pNorm)-1)*(de(ii)^(-pNorm-1))/(numEx+1);
                end
            end
        end
    end
end

function [hlinear] = linearize(h,nLevels)
hlinear = h{1};
for level = 2:nLevels
    hlinear = [hlinear; h{level}];
end

function [M0x] = buildSparseM0x(nx,ny,nz,nc0,nelx,nely,num)
M0x = sparse(num,nx(1)*ny(1)*nz(1));
dim = (13*nx(1)*ny(1)*nz(1)+3*ny(1)*nz(1)+3*nx(1)*nz(1)+nz(1))*nc0;
num2=13*nely+3;
num3=13*nelx*nely+3*nelx+3*nely+1;
iter = 1;
iM = zeros(dim,1);
jM = zeros(dim,1);
vM = ones(dim,1);
for i=1:nx(1)
    for j=1:ny(1)
        for s=1:nz(1)
            ii=(j-1)*nc0+1+(i-1)*nc0*(13*ny(1)*nc0+3)+(s-1)*nc0*num3;
            t=i+(j-1)*nx(1)+(s-1)*nx(1)*ny(1);
            for k=1:nc0
                it=ii+(k-1)*num3;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
            end
            if j==ny(1)
                for k=1:nc0
                    it=ii+(k-1)*num3+nc0;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
                if i==nx(1)
                    for k=1:nc0
                        it=ii+(k-1)*num3+num2*nc0+nc0;
                        iM(iter) = it;
                        jM(iter) = t;
                        iter = iter + 1;
                    end
                end
            end
            if i==nx(1)
                for k=1:nc0
                    it=ii+(k-1)*num3+num2*nc0;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
            end
            k11=0;k12=0;
            for k=1:nc0
                it=ii+nely+1+(k-1)*num3+k11;
                k11=k11+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
                
                it=ii+2*nely+nc0+(k-1)*num3-k12;
                k12=k12+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
            end
            k13=0;k14=0;
            if i==nx(1)
                for k=1:nc0
                    it=ii+nely+1+nc0*num2+(k-1)*num3+k13;
                    k13=k13+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                   
                    it=ii+2*nely+nc0+nc0*num2+(k-1)*num3-k14;
                    k14=k14+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
            end
          
            for k=1:nc0
                it=ii+3*nely+1+(k-1)*(num3+num2);
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
                it=ii+4*nely+2+(nc0-1)*num2+(k-1)*(num3-num2);
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
            end

            if j==ny(1)
                for k=1:nc0
                    it=ii+3*nely+1+(k-1)*(num3+num2)+nc0;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                  
                    it=ii+4*nely+2+(nc0-1)*num2+(k-1)*(num3-num2)+nc0;
                    k14=k14+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
            end
            
            iii=0;iiii=0;
            for k=1:nc0
                it=ii+5*nely+3+(k-1)*(num2+num3)+iii;
                iii=iii+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
                it=ii+6*nely+3+(k-1)*(num2+num3)+iiii;
                iiii=iiii+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
            end
            jjj=0;jjjj=0;
            for k=1:nc0
                it=(nc0-k)*num3+(k-1)*num2+ii+7*nely+3+jjj;
                jjj=jjj+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
                it=(nc0-k)*num3+(k-1)*num2+ii+8*nely+3+jjjj;
                jjjj=jjjj+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
            end
            jj1=0;jjj1=0;
            for k=1:nc0
                it=(k-1)*(num2+num3)+ii+9*nely+nc0+2-jj1;
                jj1=jj1+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
                it=(k-1)*(num2+num3)+ii+10*nely+nc0+2-jjj1;
                jjj1=jjj1+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
            end
            jj2=0;jjj2=0;
            for k=1:nc0
                it=(nc0-k)*num3+(k-1)*num2+ii+11*nely+nc0+2-jj2;
                jj2=jj2+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
                it=(nc0-k)*num3+(k-1)*num2+ii+12*nely+2+nc0-jjj2;
                jjj2=jjj2+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
            end
        end
    end
end
M0x = sparse(iM',jM',vM',num,nx(1)*ny(1)*nz(1));

function [Mx] = buildSparseMx(nLevels,nc,nx,ny,nz,nelx,nely,num)
Mx = {sparse(num,nx(1)*ny(1)*nz(1))};
for level = 2:nLevels
    Mx{level} = sparse(num,nx(level)*ny(level)*nz(level));
end
num2=13*nely+3;
num3=13*nelx*nely+3*nelx+3*nely+1;
for level = 1:nLevels
    nc_l = nc(level);
    dim = 45*nc_l*nx(level)*ny(level)*nz(level);
    iter = 1;
    iM = zeros(dim,1);
    jM = zeros(dim,1);
    vM = ones(dim,1);
    for i=1:nx(level)
        for j=1:ny(level)
            for s=1:nz(level)
                ii=(j-1)*nc_l+1+(i-1)*nc_l*num2+(s-1)*nc_l*num3;
                t=(i+(j-1)*nx(level))+(s-1)*nx(level)*ny(level);
                for k=1:nc_l
                    it = ii+(k-1)*num3+nc_l/2*num2;    % index of the element
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it = ii+(k-1)*num3+nc_l*num2/2+nc_l/2;    % index of the element
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it = ii+(k-1)*num3+nc_l*num2/2+nc_l;    % index of the element
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
                for k=1:nc_l
                    it=ii+nc_l/2+(k-1)*num3;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+nc_l/2+(k-1)*num3+nc_l*num2;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
                k111=0;k112=0;
                for k=1:nc_l
                    it=ii+nely+1+(k-1)*num3+k111+nc_l/2*num2;
                    k111=k111+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+2*nely+nc_l+(k-1)*num3-k112+nc_l/2*num2;
                    k112=k112+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
                for k=1:nc_l
                    it=ii+3*nely+1+(k-1)*(num3+num2)+nc_l/2;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+4*nely+2+(nc_l-1)*num2+(k-1)*(num3-num2)+nc_l/2;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
                kk1=0;kk2=0;kk3=0;kk4=0;kk5=0;kk6=0;
                for k=1:nc_l/2
                    it=ii+nely+1+nc_l/2+(k-1)*num3+kk1;
                    kk1=kk1+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+nely+1+nc_l/2*num3+(k-1)*num3+kk2;
                    kk2=kk2+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+nely+1+nc_l/2*(num2+1)+(k-1)*num3+kk3;
                    kk3=kk3+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+nely+1+nc_l/2*(num3+num2)+(k-1)*num3+kk4;
                    kk4=kk4+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+nely+1+nc_l/2*(2*num2+1)+(k-1)*num3+kk5;
                    kk5=kk5+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+nely+1+nc_l/2*(num3+2*num2)+(k-1)*num3+kk6;
                    kk6=kk6+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
                kk11=0;kk12=0;kk13=0;kk14=0;kk15=0;kk16=0;
                for k=1:nc_l/2
                    it=ii+2*nely+nc_l/2+(k-1)*num3-kk11;
                    kk11=kk11+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+2*nely+(k-1+nc_l/2)*num3+nc_l-kk12;
                    kk12=kk12+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+2*nely+nc_l/2+nc_l/2*num2+(k-1)*num3-kk13;
                    kk13=kk13+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+2*nely+nc_l/2*(num3+num2)+(k-1)*num3+nc_l-kk14;
                    kk14=kk14+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+2*nely+nc_l/2+nc_l*num2+(k-1)*num3-kk15;
                    kk15=kk15+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+2*nely+nc_l/2*(num3+2*num2)+(k-1)*num3+nc_l-kk16;
                    kk16=kk16+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end  
                 
                for k=1:nc_l/2
                    it=ii+4*nely+2+(k-1)*(num3-num2)+(nc_l/2-1)*num2;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+4*nely+2+(k-1)*(num3-num2)+nc_l/2*num3+(nc_l-1)*num2;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+4*nely+2+(k-1)*(num3-num2)+nc_l/2+(nc_l/2-1)*num2;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+4*nely+2+(k-1)*(num3-num2)+nc_l/2*num3+(nc_l-1)*num2+nc_l/2;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+4*nely+2+(k-1)*(num3-num2)+nc_l+(nc_l/2-1)*num2;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+4*nely+2+(k-1)*(num3-num2)+nc_l/2*num3+(nc_l-1)*num2+nc_l;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end 
                
                for k=1:nc_l/2
                    it=ii+3*nely+1+(k-1)*(num3+num2)+nc_l/2*num2;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+3*nely+1+(k-1)*(num3+num2)+nc_l/2*num3;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+3*nely+1+(k-1)*(num3+num2)+nc_l/2*num2+nc_l/2;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+3*nely+1+(k-1)*(num3+num2)+nc_l/2*num3+nc_l/2;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+3*nely+1+(k-1)*(num3+num2)+nc_l/2*num2+nc_l;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+3*nely+1+(k-1)*(num3+num2)+nc_l/2*num3+nc_l;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end 
                
                k1=0;k2=0;k11=0;k22=0;
                k3=0;k4=0;k33=0;k44=0;
                for k=1:nc_l/2
                    it=ii+7*nely+2+nc_l/2+(k-1)*(num3-num2)+(nc_l/2-1)*num2-k1;
                    k1=k1+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+8*nely+2+nc_l/2+(k-1)*(num3-num2)+(nc_l/2-1)*num2-k11;
                    k11=k11+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+7*nely+2+nc_l+(k-1)*(num3-num2)+(nc_l-1)*num2+num3*nc_l/2-k2;
                    k2=k2+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+8*nely+2+nc_l+(k-1)*(num3-num2)+(nc_l-1)*num2+num3*nc_l/2-k22;
                    k22=k22+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+(num2+1)*nc_l/2+5*nely+3+(k-1)*(num2+num3)+k3;
                    k3=k3+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+(num2+1)*nc_l/2+6*nely+3+(k-1)*(num2+num3)+k33;
                    k33=k33+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+5*nely+3+nc_l/2*num3+(k-1)*(num2+num3)+k4;
                    k4=k4+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+6*nely+3+nc_l/2*num3+(k-1)*(num2+num3)+k44;
                    k44=k44+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
            
                k11=0;k12=0;k1=0;k2=0;
                k13=0;k14=0;k3=0;k4=0;
                for k=1:nc_l/2
                    it=ii+11*nely+nc_l/2+3+(k-1)*(num3-num2)+k11+(nc_l/2-1)*num2;
                    k11=k11+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+12*nely+nc_l/2+3+(k-1)*(num3-num2)+k1+(nc_l/2-1)*num2;
                    k1=k1+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+11*nely+3+(k-1)*(num3-num2)+(nc_l-1)*num2+nc_l/2*num3+k12;
                    k12=k12+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+12*nely+3+(k-1)*(num3-num2)+(nc_l-1)*num2+nc_l/2*num3+k2;
                    k2=k2+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+(num2+1)*nc_l/2+9*nely+2+(k-1)*(num2+num3)-k13;
                    k13=k13+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+(num2+1)*nc_l/2+10*nely+2+(k-1)*(num2+num3)-k3;
                    k3=k3+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+9*nely+2+nc_l+nc_l/2*num3+(k-1)*(num2+num3)-k14;
                    k14=k14+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+10*nely+2+nc_l+nc_l/2*num3+(k-1)*(num2+num3)-k4;
                    k4=k4+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
                k21=0;k22=0;k23=0;k24=0;
                k1=0;k2=0;k3=0;k4=0;
                for k=1:nc_l/2
                    it=ii+5*nely+3+nc_l/2+(k-1)*(num2+num3)+k21;
                    k21=k21+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+6*nely+3+nc_l/2+(k-1)*(num2+num3)+k1;
                    k1=k1+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+5*nely+3+nc_l/2+nc_l/2*num3+(k-1)*(num2+num3)+k22;
                    k22=k22+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+6*nely+3+nc_l/2+nc_l/2*num3+(k-1)*(num2+num3)+k2;
                    k2=k2+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+7*nely+2+nc_l+(k-1)*(num3-num2)+(nc_l/2-1)*num2-k23;
                    k23=k23+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+8*nely+2+nc_l+(k-1)*(num3-num2)+(nc_l/2-1)*num2-k3;
                    k3=k3+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+7*nely+2+nc_l+nc_l/2*num3+(k-1)*(num3-num2)+(nc_l/2-1)*num2-k24;
                    k24=k24+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+8*nely+2+nc_l+nc_l/2*num3+(k-1)*(num3-num2)+(nc_l/2-1)*num2-k4;
                    k4=k4+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
                
                k31=0;k32=0;k33=0;k34=0;
                k1=0;k2=0;k3=0;k4=0;
                for k=1:nc_l/2
                    it=ii+9*nely+2+nc_l+num2*nc_l/2+(k-1)*(num2+num3)-k31;
                    k31=k31+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+10*nely+2+nc_l+num2*nc_l/2+(k-1)*(num2+num3)-k1;
                    k1=k1+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+9*nely+2+nc_l+num2*nc_l/2+num3*nc_l/2+(k-1)*(num2+num3)-k32;
                    k32=k32+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+10*nely+2+nc_l+num2*nc_l/2+num3*nc_l/2+(k-1)*(num2+num3)-k2;
                    k2=k2+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+11*nely+3+nc_l/2+num2*(nc_l-1)+(k-1)*(num3-num2)+k33;
                    k33=k33+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+12*nely+3+nc_l/2+num2*(nc_l-1)+(k-1)*(num3-num2)+k3;
                    k3=k3+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+11*nely+3+nc_l/2+num2*(nc_l-1)+num3*nc_l/2+(k-1)*(num3-num2)+k34;
                    k34=k34+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+12*nely+3+nc_l/2+num2*(nc_l-1)+num3*nc_l/2+(k-1)*(num3-num2)+k4;
                    k4=k4+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
                
                k41=0;k42=0;k43=0;k44=0;
                k1=0;k2=0;k3=0;k4=0;
                for k=1:nc_l/2
                    it=ii+5*nely+3+nc_l/2*num2+(k-1)*(num2+num3)+k41;
                    k41=k41+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+6*nely+3+nc_l/2*num2+(k-1)*(num2+num3)+k1;
                    k1=k1+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+5*nely+3+nc_l/2*(num2+num3)+(k-1)*(num2+num3)+k42;
                    k42=k42+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+6*nely+3+nc_l/2*(num2+num3)+(k-1)*(num2+num3)+k2;
                    k2=k2+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+7*nely+2+nc_l/2+(nc_l-1)*num2+(k-1)*(num3-num2)-k43;
                    k43=k43+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+8*nely+2+nc_l/2+(nc_l-1)*num2+(k-1)*(num3-num2)-k3;
                    k3=k3+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+7*nely+2+nc_l/2+(nc_l-1)*num2+nc_l/2*num3+(k-1)*(num3-num2)-k44;
                    k44=k44+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+8*nely+2+nc_l/2+(nc_l-1)*num2+nc_l/2*num3+(k-1)*(num3-num2)-k4;
                    k4=k4+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end
                
                k51=0;k52=0;k53=0;k54=0;
                k1=0;k2=0;k3=0;k4=0;
                for k=1:nc_l/2
                    it=ii+9*nely+2+nc_l/2+(k-1)*(num2+num3)-k51;
                    k51=k51+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+10*nely+2+nc_l/2+(k-1)*(num2+num3)-k1;
                    k1=k1+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+9*nely+2+nc_l/2+num3*nc_l/2+(k-1)*(num2+num3)-k52;
                    k52=k52+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+10*nely+2+nc_l/2+num3*nc_l/2+(k-1)*(num2+num3)-k2;
                    k2=k2+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+11*nely+3+(nc_l/2-1)*num2+(k-1)*(num3-num2)+k53;
                    k53=k53+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+12*nely+3+(nc_l/2-1)*num2+(k-1)*(num3-num2)+k3;
                    k3=k3+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    
                    it=ii+11*nely+3+(nc_l/2-1)*num2+num3*nc_l/2+(k-1)*(num3-num2)+k54;
                    k54=k54+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                    it=ii+12*nely+3+(nc_l/2-1)*num2+num3*nc_l/2+(k-1)*(num3-num2)+k4;
                    k4=k4+1;
                    iM(iter) = it;
                    jM(iter) = t;
                    iter = iter + 1;
                end    
            end
        end
    end
    Mx{level} = sparse(iM',jM',vM',num,nx(level)*ny(level)*nz(level));
end
function x_new=self_supporting(xnew,nely,nelx,nx,ny,nz,nlevels,nc,pNorm)
ii0=13*nely+3;
ii1=13*nely*nelx+3*nelx+3*nely+1;
x_new=xnew;
for level=1:nlevels
    for s=2:nz(level)
        for j=1:ny(level)
            for i=1:nx(level)
                nc_l = nc(level);
                kk=[(nc_l-1)*ii1:-ii1:0]+1;
                k=kk+(i-1)*nc_l*ii0+(j-1)*nc_l+(ii0+1-ii1)*nc_l/2+(s-1)*nc_l*ii1;
                x_kncl=((1/nc_l)*sum(xnew(k).^pNorm))^(1/pNorm);
                x_new(k(nc_l/2+1:nc_l))=x_kncl;
            end
        end
         for j=1:ny(level)
            for i=1:nx(level)+1
                nc_l = nc(level);
                kk=[(nc_l-1)*ii1:-ii1:0]+1;
                k =kk+(i-1)*nc_l*ii0+(j-1)*nc_l+(1-ii1)*nc_l/2+(s-1)*nc_l*ii1;
                x_kncl=((1/nc_l)*sum(xnew(k).^pNorm))^(1/pNorm);
                x_new(k(nc_l/2+1:nc_l))=x_kncl;
            end
         end
        for j=1:ny(level)+1
            for i=1:nx(level)
                nc_l = nc(level);
                kk=[(nc_l-1)*ii1:-ii1:0]+1;
                k =kk+(i-1)*nc_l*ii0+(j-1)*nc_l+(ii0-ii1)*nc_l/2+(s-1)*nc_l*ii1;
                x_kncl=((1/nc_l)*sum(xnew(k).^pNorm))^(1/pNorm);
                x_new(k(nc_l/2+1:nc_l))=x_kncl;
            end
        end
    end
end
function dx_new=self_supporting_dx(xnew,nely,nelx,nx,ny,nz,nlevels,nc,pNorm)
ii0=13*nely+3;
ii1=13*nely*nelx+3*nelx+3*nely+1;
x_new=xnew;
dx_new=ones(1,length(xnew));
dx1=dx_new;x1=x_new;
eps=1e-6;
for level=2:nlevels
    for s=2:nz(level)
        for j=1:ny(level)
            for i=1:nx(level)
                nc_l = nc(level);
                kk=[(nc_l-1)*ii1:-ii1:0]+1;
                k=kk+(i-1)*nc_l*ii0+(j-1)*nc_l+(ii0+1-ii1)*nc_l/2+(s-1)*nc_l*ii1;
                x_kncl=((1/nc_l)*sum(xnew(k).^pNorm))^(1/pNorm);
                x_new(k(nc_l/2+1:nc_l))=x_kncl;
                dx_new(k(nc_l/2+1:nc_l))=((1/nc_l)*sum(xnew(k).^pNorm))^(1/pNorm-1).*(xnew(k(nc_l/2+1:nc_l))).^(pNorm-1)/nc_l;
                for ii=nc_l/2+1:nc_l
                    xnew1=xnew;
                    xnew1(k(ii))=xnew(k(ii))+eps;
                    x_1=((1/nc_l)*sum((xnew1(k)).^pNorm))^(1/pNorm);
                    x1(k(nc_l/2+1:nc_l))=x_1;
                end
                dx1(k(nc_l/2+1:nc_l))=(x1(k(nc_l/2+1:nc_l))-x_new(k(nc_l/2+1:nc_l)))./eps;
            end
        end
         for j=1:ny(level)
            for i=1:nx(level)+1
                nc_l = nc(level);
                kk=[(nc_l-1)*ii1:-ii1:0]+1;
                k =kk+(i-1)*nc_l*ii0+(j-1)*nc_l+(1-ii1)*nc_l/2+(s-1)*nc_l*ii1;
                x_kncl=((1/nc_l)*sum(xnew(k).^pNorm))^(1/pNorm);
                x_new(k(nc_l/2+1:nc_l))=x_kncl;
                dx_new(k(nc_l/2+1:nc_l))=((1/nc_l)*sum(xnew(k).^pNorm))^(1/pNorm-1).*(xnew(k(nc_l/2+1:nc_l))).^(pNorm-1)/nc_l;
                for ii=nc_l/2+1:nc_l
                    xnew1=xnew;
                    xnew1(k(ii))=xnew(k(ii))+eps;
                    x_1=((1/nc_l)*sum((xnew1(k)).^pNorm))^(1/pNorm);
                    x1(k(nc_l/2+1:nc_l))=x_1;
                end
                dx1(k(nc_l/2+1:nc_l))=(x1(k(nc_l/2+1:nc_l))-x_new(k(nc_l/2+1:nc_l)))./eps;
            end
         end
        for j=1:ny(level)+1
            for i=1:nx(level)
                nc_l = nc(level);
                kk=[(nc_l-1)*ii1:-ii1:0]+1;
                k =kk+(i-1)*nc_l*ii0+(j-1)*nc_l+(ii0-ii1)*nc_l/2+(s-1)*nc_l*ii1;
                x_kncl=((1/nc_l)*sum(xnew(k).^pNorm))^(1/pNorm);
                x_new(k(nc_l/2+1:nc_l))=x_kncl;
                dx_new(k(nc_l/2+1:nc_l))=((1/nc_l)*sum(xnew(k).^pNorm))^(1/pNorm-1).*(xnew(k(nc_l/2+1:nc_l))).^(pNorm-1)/nc_l;
                for ii=nc_l/2+1:nc_l
                    xnew1=xnew;
                    xnew1(k(ii))=xnew(k(ii))+eps;
                    x_1=((1/nc_l)*sum((xnew1(k)).^pNorm))^(1/pNorm);
                    x1(k(nc_l/2+1:nc_l))=x_1;
                end
                dx1(k(nc_l/2+1:nc_l))=(x1(k(nc_l/2+1:nc_l))-x_new(k(nc_l/2+1:nc_l)))./eps;
            end
        end
    end
end


