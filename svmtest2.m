load('fourvarSVM.mat')
if numel(simoutput(1,:))>9
    simoutput=simoutput(:,1:9);
end;
for i=1:numel(simoutput(:,1))
if (isnan(table2array(simoutput(i,5))))
simoutput(i,5)=table(0);
end;
end;

for i=1:numel(simoutput(:,1))
if ~(isreal(table2array(simoutput(i,5))))
simoutput(i,5)=table(0);
end;
end;

for i=1:numel(simoutput(:,1))
if (isnan(table2array(simoutput(i,7))))
simoutput(i,7)=table(1);
end;
end;

for i=1:numel(simoutput(:,1))
if ~(isreal(table2array(simoutput(i,7))))
simoutput(i,7)=table(1);
end;
end;
szfree=[];
responder=[];
for i=1:numel(simoutput(:,1))
 szfree(i)=predict(mdlSVM,table2array(simoutput(i,4:7)));
 responder(i)=predict(RmdlSVM,table2array(simoutput(i,4:7)));
end;
szfree=szfree';
responder=responder';
szfree=table(szfree);
responder=table(responder);
simoutput=horzcat(simoutput,szfree);
simoutput=horzcat(simoutput,responder);
