function [simoutput] = frresectmodels(RESP,szfreeMDL,respMDL) 

server='localhost';
username='admin';
password='';
dbname='deckard_new';
port=27017;
conn = mongo(server,port,dbname,'UserName',username,'Password',password);
collection = "HFOs";
%Input: Patient
%Output : 
%1) Array of rank ordered LE electrodes
%2) Array of LE values
%3) Array of resection radius
%5) Array of SOZ RR
%6) Array of FR RR
%7) Array of RDRRD
%8) Array of gammaRR
%9) Array of urMLE
%10) contras

%Initialize variables
rate_array=[];
rate_array2=[];
%build Euclidian distance map
x=[];
y=[];
z=[];
    deleted_electrodes=[];
    test_query=['{"patient_id":"' RESP{1} '"}'];
    electrodes = distinct(conn,collection,"electrode",'Query',test_query);
    collection = "Electrodes";
    electrodes_2 = distinct(conn,collection,"electrode",'Query',test_query);
    total_electrodes=[electrodes electrodes_2];
    unique_electrodes=unique(total_electrodes);
  for j=1:numel(unique_electrodes)
      collection = "HFOs";
      test_query=['{"patient_id":"' RESP{1} '","electrode":"' unique_electrodes{j} '"}'];
      x_temp = distinct(conn,collection,"x",'Query',test_query);
      x_temp=cell2mat(x_temp);
      if ~isempty(x_temp)
          x(j)=x_temp(1);          
          y_temp = distinct(conn,collection,"y",'Query',test_query);
          y_temp = cell2mat(y_temp);
          y(j)=y_temp(1);
          z_temp = distinct(conn,collection,"z",'Query',test_query);
          z_temp = cell2mat(z_temp);
          z(j)=z_temp(1); 
      else
         collection = "Electrodes";
         x_temp = distinct(conn,collection,"x",'Query',test_query);
         x_temp(cellfun(@ischar,x_temp)) = {nan};
         x_temp = cell2mat(x_temp);         
         if isempty(x_temp)
             edited_unique_electrodes=unique_electrodes;
             deleted_electrodes=[deleted_electrodes j];
             x(j)=NaN;
             y(j)=NaN;
             z(j)=NaN;
         else
             x(j)=x_temp(1);
             y_temp = distinct(conn,collection,"y",'Query',test_query);
             y_temp(cellfun(@ischar,y_temp)) = {nan};
             y_temp = cell2mat(y_temp);
             y(j)=y_temp(1);
             z_temp = distinct(conn,collection,"z",'Query',test_query);
             z_temp(cellfun(@ischar,z_temp)) = {nan};
             z_temp = cell2mat(z_temp);
             z(j)=z_temp(1);
         end;
      end;  
       if x(j) == -1
         collection = "Electrodes";
         x_temp = distinct(conn,collection,"x",'Query',test_query);
         x_temp = cell2mat(x_temp);
         x(j)=x_temp(1);
         y_temp = distinct(conn,collection,"y",'Query',test_query);
         y_temp = cell2mat(y_temp);
         y(j)=y_temp(1);
         z_temp = distinct(conn,collection,"z",'Query',test_query);
         z_temp = cell2mat(z_temp);
         z(j)=z_temp(1);
       end;
       if x(j) == -1
           deleted_electrodes=[deleted_electrodes j];
           x(j)=NaN;
           y(j)=NaN;
           z(j)=NaN;
       end;
  end; 
  unique_electrodes(deleted_electrodes)=[];
  x(deleted_electrodes)=[];
  y(deleted_electrodes)=[];
  z(deleted_electrodes)=[];
  distance_matrix=inf(numel(unique_electrodes), numel(unique_electrodes));
  for j=1:numel(unique_electrodes)
    for k=1:numel(unique_electrodes)
        distance_matrix(j,k)=sqrt(((x(j)-x(k))^2)+((y(j)-y(k))^2)+((z(j)-z(k))^2));
    end;
  end;
  z2=z;

% Calculate FR MI local efficiency of each node
 fripple_mi_matrix=zeros(numel(unique_electrodes), numel(unique_electrodes));
 fripple_dmi_matrix=zeros(numel(unique_electrodes), numel(unique_electrodes));
 
% Find number of file blocks
 collection = "HFOs";
 test_query=['{"patient_id":"' RESP{1} '"}'];
 blocks=distinct(conn, collection, 'file_block','query',test_query);
 blocks=cellfun(@str2num, blocks);
 
 for j=1:numel(unique_electrodes) % j iteration
 in=[];    
 in_frono=[];  
 for z=1:numel(blocks)
     test_query=['{"patient_id":"' RESP{1} '","electrode":"' unique_electrodes{j} '","file_block":"' num2str(blocks(z)) '", "type":"' num2str(4) '","freq_pk": {$gt:350} }'];
     in_t=distinct(conn,collection,'start_t','query',test_query);
     in_t=cell2mat(in_t);
     test_query=['{"patient_id":"' RESP{1} '","electrode":"' unique_electrodes{j} '","file_block":"' num2str(blocks(z)) '", "type":"' num2str(2) '" }'];
     in_0=distinct(conn,collection,'start_t','query',test_query);
     in_0=cell2mat(in_0);
     if z == 1
         in_frono = in_t;
         in_rons = in_0;
     else 
         in_t=(in_t+(600*(z-1)));
         in_frono = [in_frono in_t];
         in_rons = [in_rons in_0];
     end;
 end;    
 
 in_frons=[];  
 for z=1:numel(blocks)
     test_query=['{"patient_id":"' RESP{1} '","electrode":"' unique_electrodes{j} '","file_block":"' num2str(blocks(z)) '", "type":"' num2str(5) '" }'];
     in_t=distinct(conn,collection,'start_t','query',test_query);
     in_t=cell2mat(in_t);
     if z == 1
         in_frons = in_t;
     else 
         in_t=(in_t+(600*(z-1)));
         in_frons = [in_frons in_t];
     end;
 end;    

 in = [in_frono in_frons];
 in0 = [in_rons in_frons];
 in = sort(in,'ascend');
 rate=numel(in)/(numel(blocks)*10);
 rate0=numel(in0);
 rate_array = [rate_array rate];
 rate_array2 = [rate_array2 rate0];


 if numel(in) == 0
   for k=1:numel(unique_electrodes)
      fripple_mi_matrix(j,k)=0;     
  end;
 else  
 
 for k=1:numel(unique_electrodes)  % k iteration
    collection = "HFOs";
    out=[];  
    out_frono=[];  
       for z=1:numel(blocks)
         test_query=['{"patient_id":"' RESP{1} '","electrode":"' unique_electrodes{k} '","file_block":"' num2str(blocks(z)) '", "type":"' num2str(4) '","freq_pk": {$gt:350} }'];
         out_t=distinct(conn,collection,'start_t','query',test_query);
         out_t=cell2mat(out_t);
       if z == 1
         out_frono = out_t;
       else 
         out_t=(out_t+(600*(z-1)));
         out_frono = [out_frono out_t];
     end;
 end;    
 
 out_frons=[];  
 for z=1:numel(blocks)
     test_query=['{"patient_id":"' RESP{1} '","electrode":"' unique_electrodes{k} '","file_block":"' num2str(blocks(z)) '", "type":"' num2str(5) '" }'];
     out_t=distinct(conn,collection,'start_t','query',test_query);
     out_t=cell2mat(out_t);
     if z == 1
         out_frons = out_t;
     else 
         out_t=(out_t+(600*(z-1)));
         out_frons = [out_frons out_t];
     end;
 end;    
 
 out = [out_frono out_frons];
 out = sort(out,'ascend');
 
 if numel(out) == 0
      fripple_mi_matrix(j,k)=0;
      fripple_dmi_matrix(j,k)=inf;
 else
     mi = AIMIE(in,out);
     if isnan(mi)
        fripple_mi_matrix(j,k)=0;
        fripple_dmi_matrix(j,k)=inf;
     else
     if mi==0
        fripple_mi_matrix(j,k)=0;
        fripple_dmi_matrix(j,k)=inf;
     else   
     fripple_mi_matrix(j,k)=mi;
     fripple_dmi_matrix(j,k)=(1/mi);
    end;
    end;
  end;
 if j==k
     fripple_mi_matrix(j,k)=0;
     fripple_dmi_matrix(j,k)=0;
 end;    
 end;
 end;
 end;
 % calculate efficiency
 Eloc = efficiency_wei(fripple_mi_matrix,2);
 
 % RD measures
             %find events    
             fripple_rate_matrix=distance_matrix;
             for j=1:numel(unique_electrodes)
             for k=1:numel(unique_electrodes) 
                 j
                 k
               test_query=['{"patient_id":"' RESP{1} '","electrode":"' unique_electrodes{j} '","type":"' num2str(4) '","freq_pk": {$gt:350} }'];
               frono_events = count(conn,collection,'Query',test_query);  
               test_query=['{"patient_id":"' RESP{1} '","electrode":"' unique_electrodes{j} '","type":"' num2str(5) '" }'];
               frons_events = count(conn,collection,'Query',test_query); 
               test_query=['{"patient_id":"' RESP{1} '"}'];
               blocks = distinct(conn,collection,"file_block",'Query',test_query);
               electrode_1_fr_rate=(frono_events+frons_events)/(numel(blocks)*10);
               test_query=['{"patient_id":"' RESP{1} '","electrode":"' unique_electrodes{k} '","type":"' num2str(4) '","freq_pk": {$gt:350} }'];
               frono_events = count(conn,collection,'Query',test_query);  
               test_query=['{"patient_id":"' RESP{1} '","electrode":"' unique_electrodes{k} '","type":"' num2str(5) '" }'];
               frons_events = count(conn,collection,'Query',test_query); 
               electrode_2_fr_rate=(frono_events+frons_events)/(numel(blocks)*10);
               if electrode_1_fr_rate == 0
                  fripple_rate_matrix(j,k)=inf;
               else
               if electrode_2_fr_rate == 0
                  fripple_rate_matrix(j,k)=inf;
               else
                  fripple_rate_matrix(j,k)=fripple_rate_matrix(j,k)*((electrode_1_fr_rate+electrode_2_fr_rate)/2);
               end;
               end;
             end;
             end;
 
 h_rate=[];
 l_rate=[];

 [hrate]=find(rate_array>=.9);
 h_rate=rate_array(hrate);
 [lrate]=find(rate_array<.9);
 l_rate=rate_array(lrate);
 channel_hrate=unique_electrodes(hrate);
 channel_lrate=unique_electrodes(lrate);
 x_hrate=x(hrate);
 y_hrate=y(hrate);
 z_hrate=z2(hrate);
 x_lrate=x(lrate);
 y_lrate=y(lrate);
 z_lrate=z2(lrate);
 le_hrate=Eloc(hrate);
 le_lrate=Eloc(lrate);

 % high rate channels sort by ascending LE or rates if LE>0.1
 le_hrate_t=le_hrate;
 [zidx]=find(le_hrate_t==0);
 le_hrate_t(zidx)=1;
 [testidx]=find(le_hrate_t<=0.1);
 if isempty(testidx)
     [~,idx]=sort(h_rate,'descend')
 else
     [~,idx]=sort(le_hrate_t,'ascend');
 end;
 le_hrate=le_hrate(idx);
 channel_hrate=channel_hrate(idx);
 channel_id_hrate=hrate(idx);
 x_hrate=x_hrate(idx);
 y_hrate=y_hrate(idx);
 z_hrate=z_hrate(idx);

 % low rate channels sort by descending rates and ascending distance.
 rad_t=[];
 [~,idx_max]=max(l_rate);
 for i=1:numel(lrate)
        rad_t(i)=sqrt((x(lrate(idx_max))-x(lrate(i)))^2+(y(lrate(idx_max))-y(lrate(i)))^2+(z2(lrate(idx_max))-z2(lrate(i)))^2); 
 end;
 temp=[lrate',l_rate',rad_t'];
 [~,idx]=sortrows(temp,[-2,3]);
 %   For example, SORTROWS(A,[2 -3]) first sorts the rows in ascending order
 %   according to column 2; then, rows with equal entries in column 2 get
 %   sorted in descending order according to column 3.
 channel_id_lrate=lrate(idx);
 l_rate=l_rate(idx);
 channel_lrate=channel_lrate(idx);
 channel_lrate=channel_lrate(idx);
 channel_id_lrate=lrate(idx);
 x_lrate=x_lrate(idx);
 y_lrate=y_lrate(idx);
 z_lrate=z_lrate(idx);

sz_free=0;
counter=1;
counter2=1;
contra=0;
rnode={};
rnode_id=[];
% 
hratenum=numel(channel_id_hrate);
lratenum=numel(channel_id_lrate);
tratenum=hratenum+lratenum;
% % Derive and rank order the simulated resected nodes.
% % Don't select contra nodes 
if ~isempty(channel_id_hrate)
for i=1:tratenum
     if i<=numel(channel_id_hrate)  
         if i==1
             rlat=(x_hrate(1)>0);
             rnode{counter}=channel_hrate(1);
             rnode_id(counter)=horzcat(rnode_id,channel_id_hrate(1));
         else
           if (x_hrate(i)>0)==rlat
             counter=counter+1;
             rnode{counter}=channel_hrate(i);
             rnode_id=horzcat(rnode_id,channel_id_hrate(i));  
           end;
         end; 
     else
       if (x_lrate(counter2)>0)==rlat
          counter=counter+1;
          rnode{counter}=channel_lrate(counter2);
          rnode_id=horzcat(rnode_id,channel_id_lrate(counter2));
          counter2=counter2+1;
       end;
     end;
end;
else
for i=1:lratenum
         if i<=numel(channel_id_lrate)  
         if i==1
             rlat=(x_lrate(1)>0);
             rnode{counter}=channel_lrate(1);
             rnode_id=horzcat(rnode_id,channel_id_lrate(counter));
         else
           if (x_lrate(i)>0)==rlat
             counter=counter+1;
             rnode{counter}=channel_lrate(i);
             rnode_id=horzcat(rnode_id,channel_id_lrate(i));  
           end;
         end; 
 
     end
 end;
end;

% find radius from FR desync node using rank order
x_matrix=zeros(numel(unique_electrodes),numel(rnode_id));
y_matrix=zeros(numel(unique_electrodes),numel(rnode_id));
z_matrix=zeros(numel(unique_electrodes),numel(rnode_id));

for i=1:numel(rnode_id)
    if i~=1
        for j=1:i
            x_matrix(rnode_id(j),i)=abs(x(rnode_id(1))-x(rnode_id(j)));
            y_matrix(rnode_id(j),i)=abs(y(rnode_id(1))-y(rnode_id(j)));
            z_matrix(rnode_id(j),i)=abs(z2(rnode_id(1))-z2(rnode_id(j)));
        end;
    end;
end;

max_r=[];
tempvar=[];
for i=1:numel(rnode)
    tempvar(1)=max(x_matrix(:,i));
    tempvar(2)=max(y_matrix(:,i));
    tempvar(3)=max(z_matrix(:,i));
    max_r(i)=max(tempvar);
end;

max_r=max_r+10;

% Find channels within max R sphere centered on FR desync node
inspere_num=[];
insphere_id={''};
for i=1:numel(max_r)
    counter=0;
for j=1:numel(unique_electrodes)
    radius=sqrt((x(rnode_id(1))-x(j))^2+(y(rnode_id(1))-y(j))^2+(z2(rnode_id(1))-z2(j))^2);
    if radius <= max_r(i)
       counter=counter+1;
       insphere_num(i,counter)=j;
       insphere_id{i,counter}=unique_electrodes{j};
    end;
end;
end;

delidx=[];
% Remove contra channels from FR desync nodes
if x(rnode_id(1))>0
    for i=1:numel(inspere_num)
        if x(inspere_num(i))<0
            delidx=vertcat(delidx,i);
        end;
    end;
else
    for i=1:numel(inspere_num)
        if x(inspere_num(i))>0
            delidx=vertcat(delidx,i);
        end;
    end;
end;

if ~isempty(delidx)
inspere_num(delidx,:)=[];
inspere_id(delidx,:)=[];
end;

% Calculate RDRRD
hRDRRD=zeros(numel(unique_electrodes),numel(unique_electrodes),numel(max_r));
for i=1:numel(max_r)
hRDRRD(:,:,i)=fripple_rate_matrix;
for j=1:numel(unique_electrodes)
    for k=1:numel(unique_electrodes)
         if isempty(intersect(insphere_num(i,:),j))
             hRDRRD(j,k,i)=inf;
         end;    
         if isempty(intersect(insphere_num(i,:),k))
             hRDRRD(j,k,i)=inf;
         end;
         if j==k
             hRDRRD(j,k,i)=0;
         end;
    end;
end;
end;

mRDRRD=[];
for i=1:numel(max_r)
    [~,~,~,r_radius_all,~] = charpath(fripple_rate_matrix,0,0);
    [~,~,~,r_radius_resected,~] = charpath(hRDRRD(:,:,i),0,0);
    if isnan(r_radius_resected)
        r_radius_resected=0;
    end;
    mRDRRD(i)=sqrt(r_radius_all-r_radius_resected);
end;

% Calculate gRR
hgRR=zeros(numel(unique_electrodes),numel(unique_electrodes),numel(max_r));
for i=1:numel(max_r)
hgRR(:,:,i)=fripple_dmi_matrix;
for j=1:numel(unique_electrodes)
    for k=1:numel(unique_electrodes)
         if isempty(intersect(insphere_num(i,:),j))
                if isempty(intersect(insphere_num(i,:),k))
                    hgRR(j,k,i)=inf;
                end;
         end;
         if j==k
             hgRR(j,k,i)=0;
         end;
    end;
end;
end;

mgRR=[]
for i=1:numel(max_r)
    [lambda_all,~,~,~,~] = charpath(fripple_dmi_matrix,0,0);
    [lambda_r,~,~,~,~] = charpath(hgRR(:,:,i),0,0);
    mgRR(i)=lambda_r./lambda_all;
end;

% uRLE
% Calculate uRmatrix
uR=zeros(numel(unique_electrodes),numel(unique_electrodes),numel(max_r));
for i=1:numel(max_r)
uR(:,:,i)=fripple_mi_matrix;
for j=1:numel(unique_electrodes)
    for k=1:numel(unique_electrodes)
         if ~isempty(intersect(insphere_num(i,:),j))
         uR(j,k,i)=0;
         end;
         if ~isempty(intersect(insphere_num(i,:),k))
         uR(j,k,i)=0;
         end;
         if j==k
             uR(j,k,i)=1;
         end;
    end;
end;
end;

urEloc=[];
for i=1:numel(max_r)
     temp = efficiency_wei(uR(:,:,i),2);
     urEloc=horzcat(urEloc, temp);
end;

urmLE=[];
for i=1:numel(max_r)
    temp=urEloc(:,i);
    [idx]=find(isnan(temp)==1);
    [idx2]=find(temp==0);
    idx3=vertcat(idx,idx2);
    temp(idx3)=[];
    urmLE(i)=mean(temp);
end;

% SOZ
soz_array=[];
for k=1:numel(unique_electrodes)
              test_query=['{"patient_id":"' RESP{1} '","electrode":"' unique_electrodes{k} '"}'];
              collection = "HFOs";
              soz = distinct(conn,collection,"soz",'Query',test_query);
              soz=cell2mat(soz);
              if ~isempty(soz)       soz=str2double(soz);                     
              soz=soz(1);
              else
              collection = "Electrodes";
              soz = distinct(conn,collection,"soz",'Query',test_query);
              soz=cell2mat(soz);
              if ~isempty(soz)       soz=str2double(soz);       end;
              soz=soz(1);
              end;
              if soz == 0
                 soz_array(k)=0;
              else
                 soz_array(k)=1;
              end;
end;

% Resected
resected=[];
for k=1:numel(unique_electrodes)
              test_query=['{"patient_id":"' RESP{1} '","electrode":"' unique_electrodes{k} '"}'];
              collection = "HFOs";
              resected = distinct(conn,collection,"resected",'Query',test_query);
              resected=cell2mat(resected);
              if ~isempty(resected)       resected=str2double(resected);                     
              resected=resected(1);
              else
              collection = "Electrodes";
              resected = distinct(conn,collection,"resected",'Query',test_query);
              resected=cell2mat(resected);
              if ~isempty(resected)       resected=str2double(resected);       end;
              resected=resected(1);
              end;
              if resected == 0
                 resected_array(k)=0;
              else
                 resected_array(k)=1;
              end;
end;

%of resected electrodes predicted by model
percent_r=[];
novel_r=[];
for i=1:numel(max_r)
[idx1]=find(resected_array==1);
temp=insphere_num(i,:);
[idxx]=find(temp==0);
temp(idxx)=[];
overlapchans=numel(intersect(idx1,temp));
[int,ia,ib]=intersect(idx1,temp);
temp2=temp;
temp2(ib)=[];
percent_r(i)=overlapchans./numel(idx1);
novel_r(i)=numel(temp2)./numel(temp);
end;

% uRLE
vSOZRR=[];
for i=1:numel(max_r)
  denom=sum(soz_array);
  temp=insphere_num(i,:);
  [idx]=find(temp==0);
  temp(idx)=[];
  sozchan=find(soz_array==1);
  numerator=numel(intersect(temp,sozchan));
  vSOZRR(i)=numerator./denom
end
fprintf('to here')

% FR RR
vFRRR=[];
for i=1:numel(max_r)
  temp=insphere_num(i,:);
  [idx]=find(temp==0);
  temp(idx)=[];
  temp2=rate_array(temp);
  numerator=sum(temp2);
  denom=sum(rate_array);
  vFRRR(i)=numerator./denom;
end;

%HFOonS resection ratio

vRonSRR=[];
for i=1:numel(max_r)
  temp=insphere_num(i,:);
  [idx]=find(temp==0);
  temp(idx)=[];
  temp2=rate_array2(temp);
  numerator=sum(temp2);
  denom=sum(rate_array2);
  vRonSRR(i)=numerator./denom;
end;

max_r=max_r';
vSOZRR=vSOZRR';
vFRRR=vFRRR';
vRonSRR=vRonSRR';
mRDRRD=mRDRRD';
mgRR=mgRR';
urmLE=urmLE';
percent_r=percent_r';
novel_r=novel_r';
simoutput=table(max_r,vSOZRR,vRonSRR,vFRRR,mRDRRD,mgRR,urmLE,percent_r,novel_r);

% szfree=[];
% responder=[];
% for i=1:numel(simoutput(:,1))
%  szfree(i)=predict(szfreeMDL,table2array(simoutput(i,4:6)))
%  responder(i)=predict(respMDL,table2array(simoutput(i,4:6)))
% end;
% szfree=szfree';
% responder=responder';
% szfree=table(szfree);
% responder=table(responder);
% simoutput=horzcat(simoutput,szfree);
% simoutput=horzcat(simoutput,responder);
