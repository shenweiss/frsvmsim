%R01 SVM
load('workingdata0817.mat')
load('frmi0825.mat')
load('ratedistance.mat')
load('FRtables.mat')
A=[];
szfree=vertcat(ones(10,1),zeros(13,1));
resp=vertcat(ones(13,1),zeros(10,1));
resp(21)=1;
nulls=[4 13 15 17 19];

rr=[];
rr350=[];
for i=1:23
    [idx,~]=find(frono_table.patients_frono==i);
    [idx2,~]=find(frons_table.patients_frons==i);
    [idx3,~]=find(frono_table.resected_frono==1);
    [idx4,~]=find(frons_table.resected_frons==1);
    [idx5,~]=find(frono_table.freq_frono>=350);
    resectedfrono=intersect(idx,idx3);
    resectedfrons=intersect(idx2,idx4);
    rr(i)=((numel(resectedfrono)+numel(resectedfrons))/(numel(idx)+numel(idx2)));
    total350frono=intersect(idx,idx5);
    resectedfrono350=intersect(resectedfrono,idx5);
    rr350(i)=((numel(resectedfrono350)+numel(resectedfrons))/(numel(total350frono)+numel(idx2)));
end;
rr350=rr350';

[idx]=isnan(ratedistance.fripple_rate_resected);
ratedistance.fripple_rate_resected(idx)=0;
radius_rate_sub=ratedistance.fripple_rate-ratedistance.fripple_rate_resected;
radius_rate_sub_sqrt=sqrt(radius_rate_sub);

[frmi]=calc_localmeasures_0819(frmi);
temp=frmi.lambda_r;
temp(isnan(temp))=0;
lambda_nr=temp./frmi.lambda;
[idx]=isnan(frmi.lambda);
lambda_nr(nulls)=NaN;

ur_mLE=frmi.nonresected_ur_leff_mean';

rate=frmi.rate_array;
patients=frmi.patient_array;
e_l=frmi.e_l;
e_l_nr=frmi.e_l_nr;
resected=frmi.resected_array
responder=frmi.responder_array';
szfree=vertcat(zeros(1194,1),ones(1375,1));
glmm_table=table(rate,e_l,e_l_nr,patients,szfree,responder,resected);
glmm_table2=glmm_table;
idx=find(rate==0);
glmm_table2(idx,:)=[];
idx2=find(e_l==0);
idx=vertcat(idx,idx2);
glmm_table(idx,:)=[];
X=[glmm_table.e_l log10(glmm_table.rate) ];
[kmeans_idx,centroid]=kmeans(X,3,'distance','cityblock');
load('kmeans_backup.mat')

patientlist=unique(glmm_table2.patients);
totalnode_r=[];
totalnode_nr=[];
zeronode_r=[];
zeronode_nr=[];
szfree=[];
clust1_r=[];
clust2_r=[];
clust3_r=[];
clust1_nr=[];
clust2_nr=[];
clust3_nr=[];
for i=1:numel(patientlist)
    [idx]=find(glmm_table2.patients==patientlist(i));
    [idx2]=find(glmm_table2.resected==1);
    totalnode_r(i)=numel(intersect(idx,idx2));
    [idx2]=find(glmm_table2.resected==0);
    totalnode_nr(i)=numel(intersect(idx,idx2));
    [idx3]=find(glmm_table2.e_l==0);
    [idx2]=find(glmm_table2.resected==1);
    zeronode_r(i)=numel(intersect(intersect(idx,idx2),idx3));
    [idx2]=find(glmm_table2.resected==0);
    zeronode_nr(i)=numel(intersect(intersect(idx,idx2),idx3));
    if i<=10
        szfree(i)=1;
    else
        szfree(i)=0;
    end;
    [idx]=find(glmm_table.patients==patientlist(i));
    [idx2]=find(glmm_table.resected==1);
    [idx3]=find(kmeans_idx==1);
    temp=intersect(intersect(idx,idx2),idx3);
    if ~isempty(temp)
    clust1_r(i)=numel(temp);
    else
    clust1_r(i)=NaN;
    end;
    [idx3]=find(kmeans_idx==2);
    temp=intersect(intersect(idx,idx2),idx3);
    if ~isempty(temp)
    clust2_r(i)=numel(temp);
    else
    clust2_r(i)=NaN;
    end;
    [idx3]=find(kmeans_idx==3);
    temp=intersect(intersect(idx,idx2),idx3);
    if ~isempty(temp)
    clust3_r(i)=numel(temp);
    else
    clust3_r(i)=NaN;
    end;

    [idx2]=find(glmm_table.resected==0);
    [idx3]=find(kmeans_idx==1);
    temp=intersect(intersect(idx,idx2),idx3);
    if ~isempty(temp)
    clust1_nr(i)=numel(temp);
    else
    clust1_nr(i)=NaN;
    end;
    [idx3]=find(kmeans_idx==2);
    temp=intersect(intersect(idx,idx2),idx3);
    if ~isempty(temp)
    clust2_nr(i)=numel(temp);
    else
    clust2_nr(i)=NaN;
    end;
    [idx3]=find(kmeans_idx==3);
    temp=intersect(intersect(idx,idx2),idx3);
    if ~isempty(temp)
    clust3_nr(i)=numel(temp);
    else
    clust3_nr(i)=NaN;
    end;
end;
p_vals=[];
p_vals2=[];
zeronode_r(isnan(zeronode_r))=0;
clust1_r(isnan(clust1_r))=0;
clust2_r(isnan(clust2_r))=0;
clust3_r(isnan(clust3_r))=0;
zeronode_nr(isnan(zeronode_nr))=0;
clust1_nr(isnan(clust1_nr))=0;
clust2_nr(isnan(clust2_nr))=0;
clust3_nr(isnan(clust3_nr))=0;

zero_ratio_r=(zeronode_r./totalnode_r)';
zero_ratio_nr=(zeronode_nr./totalnode_nr)';
clust1_ratio_r=(clust1_r./totalnode_r)';
clust2_ratio_r=(clust2_r./totalnode_r)';
clust3_ratio_r=(clust3_r./totalnode_r)';
clust1_ratio_nr=(clust1_nr./totalnode_nr)';
clust2_ratio_nr=(clust2_nr./totalnode_nr)';
clust3_ratio_nr=(clust3_nr./totalnode_nr)';
[idx]=find(lambda_nr==0);
lambda_nr(idx)=1;
X=[rr350 radius_rate_sub_sqrt lambda_nr];
ssamp=[4 13 15 17 19];
X(ssamp,:)=[];
szfree(ssamp)=[];
resp(ssamp)=[];
mdlSVM = fitcsvm(X,szfree,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto')
CVMdl = crossval(mdlSVM,'kfold',100)
A(1)=kfoldLoss(CVMdl)
RmdlSVM = fitcsvm(X,resp,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto')
CVMdl = crossval(mdlSVM,'kfold',100)
A(2)=kfoldLoss(CVMdl)
X=[radius_rate_sub_sqrt lambda_nr ur_mLE rr350];
mdlSVM = fitcsvm(X,szfree,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto')
CVMdl = crossval(mdlSVM,'kfold',10)
A(3)=kfoldLoss(CVMdl)
mdlSVM = fitcsvm(X,resp,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto')
CVMdl = crossval(mdlSVM,'kfold',10)
A(4)=kfoldLoss(CVMdl)