% IO008
dat=bnetout.coord;
dat=horzcat(dat,bnetout.resected);
dat=horzcat(dat,bnetout.rate);
writematrix(dat,'IO008e','Delimiter','\t');
for i=1:numel(bnetout.mi(:,1))
    for j=1:numel(bnetout.mi(:,1))
        if ~isnan(bnetout.mi(i,j))
            if bnetout.mi(i,j)~=0
            bnetout.mi(i,j)=1/bnetout.mi(i,j);
            end;
        end;
    end;
end;
dat=bnetout.mi;
writematrix(dat,'IO008mi','Delimiter','\t');
dat=bnetout.coord;
dmy0=zeros(numel(dat(:,1)),1);
dmy1=dmy0
dmy1(bnetout.rchan)=45;
dat=horzcat(dat,dmy0,dmy1)
writematrix(dat,'IO008r','Delimiter','\t');

