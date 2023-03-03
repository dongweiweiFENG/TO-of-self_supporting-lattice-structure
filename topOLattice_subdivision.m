% topOLattice_subdivision(128,64,0.5,1,400)

% V -- the coordinate of the nodes
% E -- the information about the nodes connected to the struts and the length of the struts
% A -- the vector about cross-section of struts
% L -- the vector about the length of struts
% xPhys -- the vector about the density of struts
% cList -- the iterative compliance information
% vList -- the iterative volume information
% nelx x nely -- the resolution about the ground lattice structure
% volfraction -- the volume fraction
% bBlackWite -- 0 for no heaviside projection and 1 for heaviside projection
% ItMax -- the total iteration steps

function [V,E,A,L,xPhys,cList,vList]=topOLattice_subdivision(nelx,nely,volfrac,bBlackWhite,ItMax)

clc;
close all;
MMA = 1;

folder = sprintf('images-%8.4f',volfrac);
mkdir(folder);
cList = zeros(ItMax,1);
vList = zeros(ItMax,1);
sList = zeros(ItMax,1);

penal = 3;
pNorm = 16;

A0=pi; % the strut's initial cross-section
l=3*ceil(sqrt(A0/pi)); % the length of the minimum vertical struts
%the vector of struts length L
L=ones(1,5*nely*nelx+nely)*l;
for i=1:4*nely
    L(1,nely+i:5*nely:end)=l*sqrt(2)/2;
end
v0=sum(L*A0);

%Coordinates corresponding to nodes V
X = repmat(0:l:nelx*l, nely+1, 1);
X1=repmat(l/2:l:nelx*l, nely, 1);
X = flipud(X);
Y = repmat([nely*l:-l:0]', 1, nelx+1);
Y1 = repmat([nely*l-l/2:-l:0]', 1, nelx);

V0 = [X(:), Y(:), zeros(length(X(:)), 1)];
V1=[X1(:), Y1(:), zeros(length(X1(:)), 1)];
V=zeros(length(V0)+length(V1),3);
for i=1:length(V)
    if rem(i,nelx+1)==0 || rem(i,nelx+1)>nelx/2+1
        if rem(i,nelx+1)==0
            V(i,:)=V1(floor(i/(nelx+1))*nelx/2+rem(i,(nelx+1)),:);
        else  
            V(i,:)=V1(floor(i/(nelx+1))*nelx/2+rem(i,(nelx+1))-nelx/2-1,:);
        end
    else
        V(i,:)=V0(floor(i/(nelx+1))*(nelx/2+1)+rem(i,(nelx+1)),:);
    end
end

nLevels = log2(nely);   %The maximal subdivision level
% Cell dimension
nc0 = 2^nLevels;   % Coarset cell, e.g., nc0*nc0 = 32*32  --2^m
nc = zeros(nLevels,1); 
for level = 1:nLevels   %k
    nc(level) = nc0/(2^(level-1));
end

% Elements within cell
nx0 = nelx/nc0;
ny0 = nely/nc0;
nx = zeros(nLevels,1);
ny = zeros(nLevels,1); 
for level = 1:nLevels
    nx(level) = nx0 * 2^(level-1);
    ny(level) = ny0 * 2^(level-1); 
end

h_min = 0;%0.00001;

% Design variables
h0 = ones(nx0*ny0,1);
h = {zeros(nx(1)*ny(1),1)};
for level = 2:nLevels
    h{level} = zeros(nx(level)*ny(level),1);
end

% A copy of design variables, for caculating the initial volume
h1 = {ones(nx(1)*ny(1),1)};
for level = 2:nLevels
    h1{level} = ones(nx(level)*ny(level),1);
end

% Design variables projected by the subdivision rule
hp = h; 
% Design variables converted to black-white design
hbw = h;

%%
num_x = 5*nely*nelx+nely;% the total number of struts

% the matrix that connect the state variables and struts density
[M0x] = buildSparseM0x(nx,ny,nc0,nelx,nely,num_x);
[Mx] = buildSparseMx(nLevels,nc,nx,ny,nelx,nely,num_x);  

[xFix] = designToX(M0x,h0,Mx,nLevels,h);
[xFix1] = designToX(M0x,h0,Mx,nLevels,h1);
xFixx=A0*L'.*xFix;
xFixx1=A0*L'.*xFix1;
initFrac = (volfrac*v0 - sum(xFixx(:))) / (sum(xFixx1(:)) - sum(xFixx(:)));
for level = 1:nLevels
    h{level} = initFrac * h1{level};
end

[hlinear] = linearize(h,nLevels);
[hplinear] = linearize(hp,nLevels);

numDesign = 0; % Number of total design variables
for level = 1:nLevels
    numDesign = numDesign + size(h{level},1);
end

[hde,hdeCo] = initializeDependency(nx,ny,nLevels);

[hdeEx,hdeExiLevels,hdeExNum,hdeExCo] = initializeExtendedDependency(nx,ny,nLevels,hde);

%% subdivisionProjection
[hp,hdeExCo] = restrictedSubdivisionProjectionPNorm(h,hp,h0,hdeEx,hdeExCo,hdeExiLevels,hdeExNum,pNorm,nLevels,nx,ny);

beta = 1;       % beta continuation
eta = 0.5;      % projection threshold, fixed at 0.5

%% Adapted from 88 lines
E0 = 1;
% nu = 0.3;

%% PREPARE FINITE ELEMENT ANALYSIS

[ssk,edofMat]=Bar2D_Stiffness1_2_new(-pi/2,-pi/4,pi/4,nelx,nely,E0,A0,l);

e=edofMat(:,[2,4])/2;
E=[e,L']; %the struts information about its connected nodes and  length
num=5*nely*nelx+nely;
iK = reshape(kron(edofMat,ones(4,1))',16*num,1);
jK = reshape(kron(edofMat,ones(1,4))',16*num,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2*((2*nely+1)*nelx+nely+1),1,-1,2*((2*nely+1)*nelx+nely+1),1); %-nely
U = zeros(2*((2*nely+1)*nelx+nely+1),1);
fixeddofs = union([1:1:2*(nely+1)],[1]);  
alldofs = [1:2*((2*nely+1)*nelx+nely+1)];  
freedofs = setdiff(alldofs,fixeddofs);

%% INITIALIZE ITERATION
[x] = designToX(M0x,h0,Mx,nLevels,h);
xPhys = x;
loop = 0;
loopbeta = 0;
change = 1;
hlinearold1 = hlinear;
hlinearold2 = hlinear;
low = 0;
upp = 0;
p=4; %  the penalization parameter
kkk=1;
A=A0*xPhys;
%% START ITERATION     
while loop < ItMax && change>=0.0001
    loop = loop + 1;
    loopbeta = loopbeta + 1;

    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    [c,dcx]=K_2D_p(ssk,iK,jK,num,F,freedofs,xPhys,edofMat,p,U);
    v=A0*L*xPhys;   
    dvx=L*A0/(sum(A0*L)*volfrac);
    
    %calculate the gradient information after self-supporting filter
    [dcx,dvx]=self_dc_dv(dcx,dvx,xPhys,nely,nx,ny,nLevels,nc,pNorm); 

    cList(loop) = c;
    vList(loop) = v/v0;
    
	dch = {zeros(1,size(h{1},1))};
    dvh = {zeros(1,size(h{1},1))};
    for level = 2:nLevels
        dch{level} = zeros(1,size(h{level},1));
        dvh{level} = zeros(1,size(h{level},1));
    end
    dcxlinear = reshape(dcx',1,num);
    dvxlinear = reshape(dvx',1,num);
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
    
    % subdivision
    dhpdh = zeros(numDesign,numDesign);
    offset = zeros(nLevels,1);
    offset(1,1) = 0;
    for level = 2:nLevels
        offset(level,1) = offset(level-1,1) + size(h{level-1},1);
    end
    for level = 1:nLevels
        for i=1:nx(level)
            for j=1:ny(level)
                it = (j-1)*nx(level)+(i-1)+1;
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
    dc = dc*dhpdh;
    dv = dv*dhpdh;

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

        fval = v/ (sum(A0*L)*volfrac) - 1;
        
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
        % subdivision
        [hp,hdeExCo] = restrictedSubdivisionProjectionPNorm(hbw,hp,h0,hdeEx,hdeExCo,hdeExiLevels,hdeExNum,pNorm,nLevels,nx,ny);
        % covert state variables to struts denstity
        [xnew] = designToX(M0x,h0,Mx,nLevels,hp);
        
        % self-supporting filter
        xnew=self_ff(xnew,nely,nx,ny,nLevels,nc,pNorm);
        
    end
    
    change = max(abs(xnew(:)-x(:)));
    x = xnew; 
    xPhys = x;
    h = hnew;
    [hlinear] = linearize(h,nLevels);

    %% Beta-continuation for Heaviside projection
    if bBlackWhite
        if beta < 64 && loopbeta >= 60
            beta = 2*beta;
            loopbeta = 0;
            change = 1/beta;
            fprintf('Parameter beta increased to %g.\n',beta);
        end
    end
    
    tmp = (x.*(1-x)*4);
    sharp = sum(tmp(:))/(4*nelx*nely+nelx+nely);
    sList(loop) = sharp;

    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%11.4f ch.:%7.3f, Sharpness: %7.3f\n',loop,c, ...
        v/ (sum(A0*L)),change,sharp);
end
%% PLOT DENSITIES
clearvars -except V E A L xPhys A0
save 2d_subdivision;
A=A0*xPhys;
figure;
clf;
D=2*sqrt(A0/pi*xPhys);
drawframe(V,E,D,[0.5,0.5,0.5]);
set(gcf,'color','w');

%% 
function [hde,hdeCo] = initializeDependency(nx,ny,nLevels)
hde = {zeros(nx(1)*ny(1),1)}; % Design variables' dependence
for level = 2:nLevels % From second level upwards, the dependence to zeroth is neglected.
    hde{level} = zeros(nx(level)*ny(level),level-1);
end
% The first level doesn't depend on the zeroth. The following is added for
% consistency
for i=1:nx(1)
    for j=1:ny(1)
        it = (j-1)*nx(1)+(i-1)+1;
        ic = it;
        hde{1}(it) = ic;
    end
end

hdeCo = {zeros(nx(1)*ny(1),1)}; % Design variables' dependence, coefficient of derivative
for level = 2:nLevels % From second level upwards, the dependence to zeroth is neglected.
    hdeCo{level} = zeros(nx(level)*ny(level),level);
end

for level = 2:nLevels
   for i=1:nx(level)
        for j=1:ny(level)
            it = (j-1)*nx(level)+(i-1)+1;
            ic = (floor((j+1)/2)-1)*nx(level-1) + (floor((i+1)/2)-1) + 1;
            hde{level}(it,level-1) = ic;
            if level > 2
                for ii = 1:(level-2)
                    hde{level}(it, ii) = hde{level-1}(ic,ii);
                end
            end
        end
    end
end

function [hdeEx,hdeExiLevels,hdeExNum,hdeExCo] = initializeExtendedDependency(nx,ny,nLevels,hde)
hdeEx = {zeros(nx(1)*ny(1),1)};
for level = 2:nLevels
    hdeEx{level} = zeros(nx(level)*ny(level),3^(level-1));
end

hdeExiLevels = {zeros(nx(1)*ny(1),1)};
for level = 2:nLevels
    hdeExiLevels{level} = zeros(nx(level)*ny(level),3^(level-1));
end

hdeExNum = {zeros(nx(1)*ny(1),1)};
for level = 2:nLevels
    hdeExNum{level} = zeros(nx(level)*ny(level),1);
end
for level = 2:nLevels
   for i=1:nx(level)
        for j=1:ny(level)
            it = (j-1)*nx(level)+(i-1)+1;
            ic = (floor((j+1)/2)-1)*nx(level-1) + (floor((i+1)/2)-1) + 1;
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

            icx = 0;
            icy = 0;
            
            if ix >= 1 && ix <= nx(level)
                itx = (j-1)*nx(level)+(ix-1)+1;
                icx = hde{level}(itx,level-1);
            end
            if jy >= 1 && jy <= ny(level)
                ity = (jy-1)*nx(level)+(i-1)+1;
                icy = hde{level}(ity,level-1);
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
                end
            end
        end
   end
end 

hdeExCo = {zeros(nx(1)*ny(1),1)}; % Design variables' dependence, coefficient of derivative
for level = 2:nLevels % From second level upwards, the dependence to zeroth is neglected.
    hdeExCo{level} = zeros(nx(level)*ny(level),3^(level-1));
end

function [x] = designToX(M0x,h0,Mx,nLevels,h)
x = M0x*h0;
for level = 1:nLevels
    x = x + Mx{level}*h{level};
end

function [hp,hdeCo] = subdivisionProjectionPNorm(h,hp,h0,hde,hdeCo,pNorm,nLevels,nx,ny)
for i=1:nx(1)
    for j=1:ny(1)
        it = (j-1)*nx(1)+(i-1)+1;
        ic = it;
        hp{1}(it) = min(h{1}(it), h0(ic));
        hdeCo{1}(it,1) = 1;
    end
end

for level = 2:nLevels
    de = zeros(level,1);
    for i=1:nx(level)
        for j=1:ny(level)
            it = (j-1)*nx(level)+(i-1)+1;
            de(level,1) = h{level}(it);
            for ii = 1:(level-1)
                ic = hde{level}(it, ii);
                de(ii,1) = h{ii}(ic);
            end
            hp{level}(it) = 0;
            for ii = 1:level
                hp{level}(it) = hp{level}(it) + (de(ii,1))^(-pNorm);
            end
            A = hp{level}(it) / level;
            hp{level}(it) = A^(-1/pNorm);

            for ii = 1:level
                hdeCo{level}(it,ii) = A^((-1/pNorm)-1)*(de(ii,1)^(-pNorm-1))/level;
            end
        end
    end
end

function [hp,hdeExCo] = restrictedSubdivisionProjectionPNorm(h,hp,h0,hdeEx,hdeExCo,hdeExiLevels,hdeExNum,pNorm,nLevels,nx,ny)
for i=1:nx(1)
    for j=1:ny(1)
        it = (j-1)*nx(1)+(i-1)+1;
        ic = it;
        hp{1}(it) = min(h{1}(it), h0(ic));
        hdeExCo{1}(it,1) = 1;
    end
end

for level = 2:nLevels
    de = zeros(3^(level-1),1);
    for i=1:nx(level)
        for j=1:ny(level)
            it = (j-1)*nx(level)+(i-1)+1;
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

function [hlinear] = linearize(h,nLevels)
hlinear = h{1};
for level = 2:nLevels
    hlinear = [hlinear; h{level}];
end

function [M0x] = buildSparseM0x(nx,ny,nc0,nelx,nely,num)
M0x = sparse(num,nx(1)*ny(1));
dim = (5*nx(1)*ny(1)+ny(1))*nc0;
iter = 1;
iM = zeros(dim,1);
jM = zeros(dim,1);
vM = ones(dim,1);
for i=1:nx(1)
    for j=1:ny(1)
        ii = (i-1)*nc0*5*nely+(j-1)*nc0;
        jj=(nc0-1)*5*nely;
        t=j+(i-1)*ny(1);
        for k=1:nc0
            it=ii+k;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
        end
        if i==nx(1)
            for k=1:nc0
                it=ii+k+nc0*5*nely;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
            end
        end
        iii=0;
        for k=nely+1:5*nely:jj+nely+1
            it=ii+k+iii;
            iii=iii+1;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
        end
        iii=0;
        for k=2*nely+1:5*nely:jj+2*nely+1
            it=ii+k+iii;
            iii=iii+1;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
        end
        jjj=0;
        for k=3*nely+nc0:5*nely:jj+3*nely+nc0
            it=ii+k-jjj;
            jjj=jjj+1;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
        end
        jjj=0;
        for k=4*nely+nc0:5*nely:jj+4*nely+nc0
            it=ii+k-jjj;
            jjj=jjj+1;
            iM(iter) = it;
            jM(iter) = t;
            iter = iter + 1;
        end
    end
end
M0x = sparse(iM',jM',vM',num,nx(1)*ny(1));

function [Mx] = buildSparseMx(nLevels,nc,nx,ny,nelx,nely,num)
Mx = {sparse(num,nx(1)*ny(1))};
for level = 2:nLevels
    Mx{level} = sparse(num,nx(level)*ny(level));
end
for level = 1:nLevels
    nc_l = nc(level);
    dim = nx(level)*ny(level)*(5*nc_l);
    iter = 1;
    iM = zeros(dim,1);
    jM = zeros(dim,1);
    vM = ones(dim,1);
    for i=1:nx(level)
        for j=1:ny(level)
            ii0=5*nely;
            ii =(i-1)*nc_l*ii0+(j-1)*nc_l;
            ii1 = nely+nc_l/2+1;
            ii2=ii0*nc_l/2;
            ii3=4*nely+nc_l/2;
            ii4=2*nely+1;
            t=i+(j-1)*nx(level);
            
            for k=1:nc_l
                it = ii+ii2+k;    % index of the element
                %Mx{level}(it,ic) = 1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
            end
            k0=0;
            k1=0;
            k00=0;
            k11=0;
            for k=1:nc_l/2
                it=ii+ii3+ii0*(k-1)-k0;
                k0=k0+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
                
                it=ii+ii3+ii0*(k-1)-k00-nely;
                k00=k00+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
                
                it=ii+ii3++nc_l/2+ii0*(nc_l/2+k-1)-k11-nely;
                k11=k11+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
                
                it=ii+ii3++nc_l/2+ii0*(nc_l/2+k-1)-k1;
                k1=k1+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
            end
            k2=0;
            k22=0;
            k3=0;
            k33=0;
            for k=1:nc_l/2
                it=ii+ii4+nc_l/2+ii0*(k-1)+k22-nely;
                k22=k22+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
                
                it=ii+ii4+nc_l/2+ii0*(k-1)+k2;
                k2=k2+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
                
                it=ii+ii4+ii0*(nc_l/2+k-1)+k33-nely;
                k33=k33+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
                
                it=ii+ii4+ii0*(nc_l/2+k-1)+k3;
                k3=k3+1;
                iM(iter) = it;
                jM(iter) = t;
                iter = iter + 1;
            end
        end
    end
    Mx{level} = sparse(iM',jM',vM',num,nx(level)*ny(level));
end

