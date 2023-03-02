function drawframe6__1(VV,EE,xPhys,color)
aa=length(EE);
k=0;
for i=1:aa 
    if xPhys(i-k)<=0.09*pi
        xPhys(i-k)=[];
        EE(i-k,:)=[];
        k=k+1;
        aa=aa-1;
    end
end
linesize=xPhys;
for i=1:length(EE)
    xx1=VV(EE(i,1),:);
    xx2=VV(EE(i,2),:);
    a=[xx1(1,1);xx2(1,1)];
    b=[xx1(1,2);xx2(1,2)];
    if linesize(i)>=1
        plot(a,b,'linewidth', linesize(i), 'color', color);
    else
        plot(a,b,'linewidth', linesize(i), 'color', color);
    end
    hold on
end
axis off
end