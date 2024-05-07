function dataout = main_nresected_types_0810(RESP) 

server='localhost';
username='admin';
password='';
dbname='deckard_new';
port=27017;
conn = mongo(server,port,dbname,'UserName',username,'Password',password);
collection = "HFOs";

% RESP = {'2061', '4122', '4124', '4145', '4163', '4166', '449', '456', '458', '463', '466', '468', '473', '477', '478', '481','IO004', 'IO006', 'IO013', 'IO014'};
% 
% NONRESP = {'4100', '4110', '479', '480', 'IO010', 'IO015', 'IO027'};
% 
% NOSURG = {'451','474','475','IO009','IO022'};
% For each subject
% 
% Identify all the unique electrodes in the HFO database, identify all the unique electrodes in the Electrode database. Merge unique electrodes
% Query MNI coordinates for each unique electrode in HFO database, and Electrode database. If not found in one, use the other. If MNI coordinates do not exist, remove electrode.
% Create two adjacency matrix for all electrodes by calculating geometric distance between the electrodes
% Query SOZ electrodes in the HFO database, Query SOZ electrodes in the Electrode database. Merge unique SOZ electrodes.
% For adjacency matrix one, delete all non-SOZ electrodes.
% Query HFO database for electrodes with total fRonO (>350 Hz) > 0, query HFO database for electrodes with total fRonS (>350 Hz) > 0. Merge these electrodes
% For adjacency matrix two, deletif ~isempty(resected)       resected=str2double(resected);       end;                             resected=resected(1);e all non-fast ripple electrodes
% Query HFO database for number of blocks, use block number to calculate the fast ripple (>350 Hz) rate
% Duplicate adjacency matrix two creating adjacency matrix three and multiply geometric distance by the respective inverse fast ripple rates in each electrode pair.
% Calculate the distance for the three adjacency matrices

r_radius_resected=[];
r_radius_fripple=[];
r_radius_fripple_rate=[];
r_radius_fripple_nresected=[];
r_radius_fripple_rate_nresected=[];
r_radius_fripple_resected=[];
r_radius_fripple_rate_resected=[];

for i=1:numel(RESP)
    x=[];
    y=[];
    z=[];
    resected_matrix=[];
    fripple_matrix=[];
    fripple_rate_matrix=[];
    deleted_electrodes=[];
    test_query=['{"patient_id":"' RESP{i} '"}'];
    electrodes = distinct(conn,collection,"electrode",'Query',test_query);
    collection = "Electrodes";
    electrodes_2 = distinct(conn,collection,"electrode",'Query',test_query);
    total_electrodes=[electrodes electrodes_2];
    unique_electrodes=unique(total_electrodes);
  for j=1:numel(unique_electrodes)
      collection = "HFOs";
      test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{j} '"}'];
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
         if numel(x_temp) > 1
         x_temp(cellfun(@ischar,x_temp)) = {nan};
         end;
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
             if numel(y_temp) > 1
             y_temp(cellfun(@ischar,y_temp)) = {nan};
             end;
             y_temp = cell2mat(y_temp);
             y(j)=y_temp(1);
             z_temp = distinct(conn,collection,"z",'Query',test_query);
             if numel(z_temp) > 1
             z_temp(cellfun(@ischar,z_temp)) = {nan};
             end;
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
  
  resected_matrix=distance_matrix;
  for j=1:numel(unique_electrodes)
      test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{j} '"}'];
      collection = "HFOs";
      resected = distinct(conn,collection,"resected",'Query',test_query);
      resected=cell2mat(resected);
      if ~isempty(resected)       
          resected=str2double(resected);       
          resected=resected(1);
      else
      collection = "Electrodes";
      resected = distinct(conn,collection,"resected",'Query',test_query);
      resected=cell2mat(resected);
      if ~isempty(resected)       
          resected=str2double(resected);      
          resected=resected(1);
      end;
      end;
      if resected == 0
          collection = "Electrodes";
          resected = distinct(conn,collection,"resected",'Query',test_query);
          resected=cell2mat(resected);
          if ~isempty(resected)      resected=str2double(resected);    resected=resected(1);     end;   
      end;    
      if resected == 0
          for k=1:numel(unique_electrodes)
              resected_matrix(j,k)=inf;
          end;
      else
          for k=1:numel(unique_electrodes)
              test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{k} '"}'];
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
                  collection = "Electrodes";
                  resected = distinct(conn,collection,"resected",'Query',test_query);
                  resected=cell2mat(resected);
                  if ~isempty(resected)       resected=str2double(resected);       end;
                  resected=resected(1);
                end;
             if resected == 0
                 resected_matrix(j,k)=inf;
             end;
          end; 
      end;
  end;
    
 
 fripple_matrix=distance_matrix;
 fripple_rate_matrix=distance_matrix;
 fripple_nresected_matrix=distance_matrix;
 fripple_rate_nresected_matrix=distance_matrix;
 fripple_resected_matrix=distance_matrix;
 fripple_rate_resected_matrix=distance_matrix;

 fprintf('fripple ');
 collection = "HFOs";
 for j=1:numel(unique_electrodes)
   electrode1_resected=0;
   electrode2_resected=0;
     % find if resected electrode 1 is resected
              test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{j} '"}'];
              collection = "HFOs";
              resected = distinct(conn,collection,"resected",'Query',test_query);
              resected=cell2mat(resected);
              if ~isempty(resected)       resected=str2double(resected);       end;
              if ~isempty(resected)
              resected=resected(1);
              else
                  collection = "Electrodes";
                  resected = distinct(conn,collection,"resected",'Query',test_query);
                  resected=cell2mat(resected);
                  if ~isempty(resected)       resected=str2double(resected);       end;
                  resected=resected(1);
              end;
                if resected == 0
                  collection = "Electrodes";
                  resected = distinct(conn,collection,"resected",'Query',test_query);
                  resected=cell2mat(resected);
                  if ~isempty(resected)       resected=str2double(resected);       end;
                  resected=resected(1);
                end;
             if resected == 0
                 electrode1_resected=0;
             else
                 electrode1_resected=1;
             end;
     % Find Events   
     collection = "HFOs";
     test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{j} '","type":"' num2str(4) '","freq_pk": {$gt:350} }'];
     num_records = count(conn,collection,'Query',test_query);
     if num_records == 0
         test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{j} '","type":"' num2str(5) '" }'];
         num_records = count(conn,collection,'Query',test_query);
     end;
     if num_records == 0
         for k=1:numel(unique_electrodes)
             fripple_matrix(j,k)=inf;
             fripple_rate_matrix(j,k)=inf;
             fripple_nresected_matrix(j,k)=inf;
             fripple_rate_nresected_matrix(j,k)=inf;
             fripple_resected_matrix(j,k)=inf;
             fripple_rate_resected_matrix(j,k)=inf;
         end;
     else
         for k=1:numel(unique_electrodes)
             electrode2_resected=0;
             
             % find if resected electrode 2 is resected
              test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{k} '"}'];
              collection = "HFOs";
              resected = distinct(conn,collection,"resected",'Query',test_query);
              resected=cell2mat(resected);
              if ~isempty(resected)       resected=str2double(resected);       end;
              if ~isempty(resected)
              resected=resected(1);
              else
                  collection = "Electrodes";
                  resected = distinct(conn,collection,"resected",'Query',test_query);
                  resected=cell2mat(resected);
                  if ~isempty(resected)       resected=str2double(resected);       end;
                  resected=resected(1);
              end;
                if resected == 0
                  collection = "Electrodes";
                  resected = distinct(conn,collection,"resected",'Query',test_query);
                  resected=cell2mat(resected);
                  if ~isempty(resected)       resected=str2double(resected);       end;
                  resected=resected(1);
                end;
             if resected == 0
                 electrode2_resected=0;
             else
                 electrode2_resected=1;
             end;
             
             %find events
             collection = "HFOs";
             test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{k} '","type":"' num2str(4) '","freq_pk": {$gt:350} }'];
              num_records = count(conn,collection,'Query',test_query);
                 if num_records == 0
                    test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{k} '","type":"' num2str(5) '" }'];
                    num_records = count(conn,collection,'Query',test_query);
                 end;
             if num_records == 0
               fripple_matrix(j,k)=inf;
               fripple_rate_matrix(j,k)=inf;
               fripple_nresected_matrix(j,k)=inf;
               fripple_rate_nresected_matrix(j,k)=inf;
               fripple_resected_matrix(j,k)=inf;
               fripple_rate_resected_matrix(j,k)=inf;
             else                       
               test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{j} '","type":"' num2str(4) '","freq_pk": {$gt:350} }'];
               frono_events = count(conn,collection,'Query',test_query);  
               test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{j} '","type":"' num2str(5) '" }'];
               frons_events = count(conn,collection,'Query',test_query); 
               test_query=['{"patient_id":"' RESP{i} '"}'];
               blocks = distinct(conn,collection,"file_block",'Query',test_query);
               electrode_1_fr_rate=(frono_events+frons_events)/(numel(blocks)*10);
               test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{k} '","type":"' num2str(4) '","freq_pk": {$gt:350} }'];
               frono_events = count(conn,collection,'Query',test_query);  
               test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{k} '","type":"' num2str(5) '" }'];
               frons_events = count(conn,collection,'Query',test_query); 
               electrode_2_fr_rate=(frono_events+frons_events)/(numel(blocks)*10);
               fripple_rate_matrix(j,k)=fripple_rate_matrix(j,k)*((electrode_1_fr_rate+electrode_2_fr_rate)/2);
               if electrode1_resected == 1
                   fripple_nresected_matrix(j,k)=inf;
                   fripple_rate_nresected_matrix(j,k)=inf;
                   fripple_rate_resected_matrix(j,k)=fripple_rate_resected_matrix(j,k)*((electrode_1_fr_rate+electrode_2_fr_rate)/2);
               else
                   if electrode2_resected == 1
                      fripple_nresected_matrix(j,k)=inf;
                      fripple_rate_nresected_matrix(j,k)=inf;
                      fripple_rate_resected_matrix(j,k)=fripple_rate_resected_matrix(j,k)*((electrode_1_fr_rate+electrode_2_fr_rate)/2);
                   else    
                      fripple_rate_nresected_matrix(j,k)=fripple_rate_nresected_matrix(j,k)*((electrode_1_fr_rate+electrode_2_fr_rate)/2);
                      fripple_resected_matrix(j,k)=inf;
                      fripple_rate_resected_matrix(j,k)=inf;
                   end;
               end;    
             end;
         end;
     end;
 end; 

 i
 [~,~,~,r_radius_resected(i),~] = charpath(resected_matrix,0,0);
 [~,~,~,r_radius_fripple(i),~] = charpath(fripple_matrix,0,0);
 [~,~,~,r_radius_fripple_rate(i),~] = charpath(fripple_rate_matrix,0,0);
 [~,~,~,r_radius_fripple_nresected(i),~] = charpath(fripple_nresected_matrix,0,0);
 [~,~,~,r_radius_fripple_rate_nresected(i),~] = charpath(fripple_rate_nresected_matrix,0,0);
 [~,~,~,r_radius_fripple_resected(i),~] = charpath(fripple_resected_matrix,0,0);
 [~,~,~,r_radius_fripple_rate_resected(i),~] = charpath(fripple_rate_resected_matrix,0,0);

end;

dataout.resected=r_radius_resected'
dataout.fripple=r_radius_fripple'
dataout.fripple_rate=r_radius_fripple_rate';
dataout.fripple_nresected=r_radius_fripple_nresected';
dataout.fripple_rate_nresected=r_radius_fripple_rate_nresected';
dataout.fripple_resected=r_radius_fripple_resected';
dataout.fripple_rate_resected=r_radius_fripple_rate_resected';