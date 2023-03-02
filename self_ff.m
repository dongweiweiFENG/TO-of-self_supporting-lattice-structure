function x_new=self_ff(xnew,nely,nx,ny,nlevels,nc,pNorm)
ii0=5*nely;
x_new=xnew;
for level=2:nlevels
%     hc=hp{level};
    for j=1:ny(level)-1
        for i=1:nx(level)
            nc_l = nc(level);
            kk=1:nc_l;
            k=kk+(i-1)*nc_l*ii0+(j-1)*nc_l+(ii0+1)*nc_l/2;
            x_kncl=(1/nc_l*sum(xnew(k).^pNorm))^(1/pNorm);
            x_new(k(nc_l/2+1:nc_l))=x_kncl;
        end
    end
end