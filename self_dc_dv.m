%calculate the gradient information after self-supporting filter in subdivision
function [df0,dfv]=self_dc_dv(df0dx,dfdx,xPhys,nely,nx,ny,nlevels,nc,pNorm)
ii0=5*nely;
df0=df0dx;
dfv=dfdx;
for level=2:nlevels
%     hc=hp{level};
    for j=1:ny(level)-1
        for i=1:nx(level)
            nc_l = nc(level);
            kk=1:nc_l;
            k=kk+(i-1)*nc_l*ii0+(j-1)*nc_l+(ii0+1)*nc_l/2;
            df0(k(nc_l/2+1:nc_l))=df0dx(k(nc_l/2+1:nc_l)).*((1/nc_l)*sum(xPhys(k).^pNorm))^(1/pNorm-1).*(xPhys(k(nc_l/2+1:nc_l)).^(pNorm-1))'/nc_l;
            dfv(k(nc_l/2+1:nc_l))=dfdx(k(nc_l/2+1:nc_l)).*((1/nc_l)*sum(xPhys(k).^pNorm))^(1/pNorm-1).*(xPhys(k(nc_l/2+1:nc_l)).^(pNorm-1))'/nc_l;
        end
    end
end
