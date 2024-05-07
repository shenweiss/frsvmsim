%% 
function [bnetout] = preparebrainnet(RESP) 
bnetout=[];
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
  coord=horzcat(x',y',z2');
  bnetout.coord=coord;

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
bnetout.mi=fripple_mi_matrix;
bnetout.rate=rate_array';

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

 bnetout.resected=resected_array';
 % calculate efficiency
 Eloc = efficiency_wei(fripple_mi_matrix,2);

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
 h_rate=h_rate(idx);
 channel_hrate=channel_hrate(idx);
 channel_id_hrate=hrate(idx);
 x_hrate=x_hrate(idx);
 y_hrate=y_hrate(idx);
 z_hrate=z_hrate(idx);


 le_lrate_t=le_lrate;
 [zidx]=find(le_lrate_t==0);
 le_lrate_t(zidx)=1;
 [testidx]=find(le_lrate_t<=0.1);
 if isempty(testidx)
     [~,idx]=sort(l_rate,'descend')
 else
     [~,idx]=sort(le_lrate_t,'ascend');
 end;
 le_lrate=le_lrate(idx);
 l_rate=l_rate(idx);
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

hratenum=numel(channel_id_hrate);
lratenum=numel(channel_id_lrate)
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
bnetout.rchan=rnode_id(1);


