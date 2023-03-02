function drawframe(VV,EE,xPhys,color)
aa=length(EE);
k=0;
for i=1:aa 
    if xPhys(i-k)<2e-1
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
    plot(a,b,'linewidth', linesize(i), 'color', color);
    hold on
end
axis off