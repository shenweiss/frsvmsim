function dataout = main_functional_fr_resect_0825(RESP, responder) 
% 
% I think I can use a measure like path length for the global network to compare across patients in the functional connectivity analysis. 
% Another analysis I should do is to compare centrality or hubs.
% I can do this analysis at the electrode level, and ask if their is a difference in hubness for all electrodes, 
% SOZ electrodes, and non-SOZ electrodes in classifying non-reponders.
% 
server='localhost';
username='admin';
password='';
dbname='deckard_new';
port=27017;
conn = mongo(server,port,dbname,'UserName',username,'Password',password);
collection = "HFOs";

% For each subject
% 
% Identify all the unique electrodes in the HFO database, identify all the unique electrodes in the Electrode database. Merge unique electrodes
% Query MNI coordinates for each unique electrode in HFO database, and Electrode database. If not found in one, use the other. If MNI coordinates do not exist, remove electrode.
% Create two adjacency matrix for all electrodes by calculating geometric distance between the electrodes
% Query SOZ electrodes in the HFO database, Query SOZ electrodes in the Electrode database. Merge unique SOZ electrodes.
% For adjacency matrix one, delete all non-SOZ electrodes.
% Query HFO database for electrodes with total fRonO (>350 Hz) > 0, query HFO database for electrodes with total fRonS (>350 Hz) > 0. Merge these electrodes
% For adjacency matrix two, delete all non-fast ripple electrodes
% Query HFO database for number of blocks, use block number to calculate the fast ripple (>350 Hz) rate
% Duplicate adjacency matrix two creating adjacency matrix three and multiply geometric distance by the respective inverse fast ripple rates in each electrode pair.
% Calculate the distance for the three adjacency matrices

e_l=[];
e_l_uw=[];
cC=[];
cC_uw=[];

e_l_nr=[];
e_l_uw_nr=[];
cC_nr=[];
cC_uw_nr=[];

Cnt=[];
Cnt_uw=[];
Cnt_nr=[];
Cnt_nr_uw=[];

lambda=[];
lambda_nr=[];
efficiency=[];
efficiency_nr=[];

radius=[]; 
radius_r=[];
diameter=[];
diameter_r=[];

is=[];
is_nr=[];
os=[];
os_nr=[];
strength=[];
strength_nr=[];
strength_uw=[];
strength_nr_uw=[];

patient_array=[];
responder_array=[];
rate_array=[];

mi_patient_array=[];
mi_resected_array=[];
mi_array=[];
mi_responder=[];
resected_all=[];

for i=1:numel(RESP)
    resected_array=[];
    x=[];
    y=[];
    z=[];
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
 
  for j=1:numel(unique_electrodes)
      test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{j} '"}'];
      collection = "HFOs";
      resected = distinct(conn,collection,"resected",'Query',test_query);
      resected=cell2mat(resected);
      if ~isempty(resected)
      resected=str2double(resected);
      end;
      if ~isempty(resected)
      resected=resected(1);
      else
      collection = "Electrodes";
      resected = distinct(conn,collection,"resected",'Query',test_query);      
      resected=cell2mat(resected);
      if ~isempty(resected)
      resected=str2double(resected);
      end;
      resected=resected(1);
      end;
      if resected == 0
          collection = "Electrodes";
          resected = distinct(conn,collection,"resected",'Query',test_query);
          resected=cell2mat(resected);
          if ~isempty(resected)
          resected=str2double(resected);
          end;
          resected=resected(1);
      end;
      if resected == 0
          resected_array=[resected_array 0];
          patient_array=[patient_array i];
          responder_array=[responder_array responder(i)];
      else
          resected_array=[resected_array 1];
          patient_array=[patient_array i];
          responder_array=[responder_array responder(i)];
      end;
  end;      
 
 fripple_dmi_matrix=inf(numel(unique_electrodes), numel(unique_electrodes));
 fripple_wmi_matrix=inf(numel(unique_electrodes), numel(unique_electrodes));
 fripple_mi_matrix=zeros(numel(unique_electrodes), numel(unique_electrodes));
 
 fripple_uw_dmi_matrix=inf(numel(unique_electrodes), numel(unique_electrodes));
 fripple_uw_wmi_matrix=inf(numel(unique_electrodes), numel(unique_electrodes));
 fripple_uw_mi_matrix=zeros(numel(unique_electrodes), numel(unique_electrodes));

 fripple_resected_dmi_matrix=inf(numel(unique_electrodes), numel(unique_electrodes));

 fripple_nresected_dmi_matrix=inf(numel(unique_electrodes), numel(unique_electrodes));
 fripple_nresected_wmi_matrix=inf(numel(unique_electrodes), numel(unique_electrodes));
 fripple_nresected_mi_matrix=zeros(numel(unique_electrodes), numel(unique_electrodes));
 
 fripple_nresected_uw_dmi_matrix=inf(numel(unique_electrodes), numel(unique_electrodes));
 fripple_nresected_uw_wmi_matrix=inf(numel(unique_electrodes), numel(unique_electrodes));
 fripple_nresected_uw_mi_matrix=zeros(numel(unique_electrodes), numel(unique_electrodes));

 % Find number of file blocks
 collection = "HFOs";
 test_query=['{"patient_id":"' RESP{i} '"}'];
 blocks=distinct(conn, collection, 'file_block','query',test_query);
 blocks=cellfun(@str2num, blocks);
 
 for j=1:numel(unique_electrodes)
 in=[];    
 in_frono=[];  
 for z=1:numel(blocks)
     test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{j} '","file_block":"' num2str(blocks(z)) '", "type":"' num2str(4) '","freq_pk": {$gt:350} }'];
     in_t=distinct(conn,collection,'start_t','query',test_query);
     in_t=cell2mat(in_t);
     if z == 1
         in_frono = in_t;
     else 
         in_t=(in_t+(600*(z-1)));
         in_frono = [in_frono in_t];
     end;
 end;    
 
 in_frons=[];  
 for z=1:numel(blocks)
     test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{j} '","file_block":"' num2str(blocks(z)) '", "type":"' num2str(5) '" }'];
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
 in = sort(in,'ascend');
 rate=numel(in)/(numel(blocks)*10);
 rate_array = [rate_array rate];
 
 if numel(in) == 0
   for k=1:numel(unique_electrodes)
      fripple_dmi_matrix(j,k)=inf;
      fripple_uw_dmi_matrix(j,k)=inf;
      fripple_wmi_matrix(j,k)=inf;
      fripple_uw_wmi_matrix(j,k)=inf;
      fripple_mi_matrix(j,k)=0;
      fripple_uw_mi_matrix(j,k)=0;
      
  end;
 else     
  for k=1:numel(unique_electrodes)      
   % check for resected
      test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{j} '"}'];
      collection = "HFOs";
      resected2 = distinct(conn,collection,"resected",'Query',test_query);
      resected2=cell2mat(resected2);
      if ~isempty(resected2)
      resected2=str2double(resected2);
      end;
      if ~isempty(resected2)
      resected2=resected2(1);
      else
      collection = "Electrodes";
      resected2 = distinct(conn,collection,"resected",'Query',test_query);      
      resected2=cell2mat(resected2);
      if ~isempty(resected2)
      resected2=str2double(resected2);
      end;
      resected2=resected2(1);
      end;
      if resected2 == 0
          collection = "Electrodes";
          resected2 = distinct(conn,collection,"resected",'Query',test_query);
          resected2=cell2mat(resected2);
          if ~isempty(resected2)
          resected2=str2double(resected2);
          end;
          resected2=resected2(1);
      end;
    collection = "HFOs";
    out=[];  
    out_frono=[];  
       for z=1:numel(blocks)
         test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{k} '","file_block":"' num2str(blocks(z)) '", "type":"' num2str(4) '","freq_pk": {$gt:350} }'];
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
     test_query=['{"patient_id":"' RESP{i} '","electrode":"' unique_electrodes{k} '","file_block":"' num2str(blocks(z)) '", "type":"' num2str(5) '" }'];
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
     fripple_dmi_matrix(j,k)=inf;
     fripple_uw_dmi_matrix(j,k)=inf;
     fripple_wmi_matrix(j,k)=inf;
     fripple_uw_wmi_matrix(j,k)=inf;
     fripple_mi_matrix(j,k)=0;
     fripple_uw_mi_matrix(j,k)=0;
 else
     mi = AIMIE(in,out);
     if isnan(mi)
        fripple_dmi_matrix(j,k)=inf;
        fripple_uw_dmi_matrix(j,k)=inf;
        fripple_wmi_matrix(j,k)=inf;
        fripple_uw_wmi_matrix(j,k)=inf;
        fripple_mi_matrix(j,k)=0;
        fripple_uw_mi_matrix(j,k)=0;
     else
     if mi==0
        fripple_dmi_matrix(j,k)=inf;
        fripple_uw_dmi_matrix(j,k)=inf;
        fripple_wmi_matrix(j,k)=inf;
        fripple_uw_wmi_matrix(j,k)=inf;
        fripple_mi_matrix(j,k)=0;
        fripple_uw_mi_matrix(j,k)=0;
     else   
     fripple_dmi_matrix(j,k)=(1/mi);
     fripple_uw_dmi_matrix(j,k)=1;
     fripple_wmi_matrix(j,k)=distance_matrix(j,k)*(mi);
     fripple_uw_wmi_matrix(j,k)=1;
     fripple_mi_matrix(j,k)=mi;
     fripple_uw_mi_matrix(j,k)=1;
     if ((resected_array(j)==1)||(resected_array(k)==1))  
             fripple_resected_dmi_matrix(j,k)=(1/mi);
             fripple_nresected_mi_matrix(j,k)=0;
             fripple_nresected_uw_mi_matrix(j,k)=0;
             fripple_nresected_dmi_matrix(j,k)=inf;
             fripple_nresected_uw_dmi_matrix(j,k)=inf;
             mi_patient_array=[mi_patient_array i];
             mi_resected_array=[mi_resected_array 1];
             mi_array=[mi_array mi];
             mi_responder=[mi_responder responder(i)];
     else
             fripple_nresected_mi_matrix(j,k)=mi;
             fripple_nresected_uw_mi_matrix(j,k)=1;
             fripple_nresected_dmi_matrix(j,k)=(1/mi);
             fripple_nresected_uw_dmi_matrix(j,k)=1;
             mi_patient_array=[mi_patient_array i];
             mi_resected_array=[mi_resected_array 0];
             mi_array=[mi_array mi];
             mi_responder=[mi_responder responder(i)];
     end;
     end;
     end;
  end;
 if j==k
     fripple_dmi_matrix(j,k)=inf;
     fripple_wmi_matrix(j,k)=inf;
     fripple_mi_matrix(j,k)=0;
 end;    
 end;
 end;
 end;

resected_all=[resected_all resected_array];
[lambda(i),efficiency(i),~,radius(i),diameter(i)] = charpath(fripple_dmi_matrix,0,0);
[uw_lambda(i),uw_efficiency(i),~,uw_radius(i),uw_diameter(i)] = charpath(fripple_uw_dmi_matrix,0,0);

[w_lambda(i),w_efficiency(i),~,w_radius(i),w_diameter(i)] = charpath(fripple_wmi_matrix,0,0); 
[lambda_r(i),efficiency_r(i),~,radius_r(i),diameter_r(i)] = charpath(fripple_resected_dmi_matrix,0,0);
[uw_lambda_r(i),uw_efficiency_r(i),~,uw_radius_r(i),uw_diameter_r(i)] = charpath(fripple_nresected_uw_dmi_matrix,0,0);

Eloc = efficiency_wei(fripple_mi_matrix,2);
Eloc_uw = efficiency_bin(fripple_uw_mi_matrix,2);
Eloc_nr=efficiency_wei(fripple_nresected_mi_matrix,2);
Eloc_uw_nr=efficiency_bin(fripple_nresected_uw_mi_matrix,2);

[is_t, os_t, strength_t]=strengths_dir(fripple_mi_matrix);
[strength_uw_t]=strengths_und(fripple_uw_mi_matrix);
[is_nr_t, os_nr_t, strength_nr_t]=strengths_dir(fripple_nresected_mi_matrix);
[strength_nr_uw_t]=strengths_und(fripple_nresected_uw_mi_matrix);

fripple_mi_matrix = weight_conversion(fripple_mi_matrix,'normalize');
fripple_dmi_matrix = weight_conversion(fripple_mi_matrix,'lengths');
fripple_nresected_mi_matrix = weight_conversion(fripple_nresected_mi_matrix,'normalize');
fripple_nresected_dmi_matrix = weight_conversion(fripple_nresected_mi_matrix,'lengths');

C = clustering_coef_wu(fripple_mi_matrix);
C_uw = clustering_coef_bu(fripple_uw_mi_matrix);
C_nr=clustering_coef_wu(fripple_nresected_mi_matrix);
C_uw_nr = clustering_coef_bu(fripple_nresected_uw_mi_matrix);

[~,Cent]=edge_betweenness_wei(fripple_dmi_matrix);    
[~,Cent_uw]=edge_betweenness_bin(fripple_uw_mi_matrix);
[~,Cent_nr]=edge_betweenness_wei(fripple_nresected_dmi_matrix);
[~,Cent_nr_uw]=edge_betweenness_bin(fripple_nresected_uw_mi_matrix);

e_l = vertcat(e_l,Eloc);
e_l_uw = vertcat(e_l_uw,Eloc_uw);
e_l_nr = vertcat(e_l_nr, Eloc_nr);
e_l_uw_nr = vertcat(e_l_uw_nr, Eloc_uw_nr);

cC = vertcat(cC, C);
cC_uw = vertcat(cC_uw, C_uw);
cC_nr = vertcat(cC_nr, C_nr);
cC_uw_nr = vertcat(cC_uw_nr, C_uw_nr);

is = [is, is_t];
is_nr = [is_nr, is_nr_t];
os = [os, os_t];
os_nr = [os_nr, os_nr_t];
strength = [strength, strength_t];
strengh_uw = [strength_uw, strength_uw_t];
strength_nr = [strength_nr, strength_nr_t];
strength_nr_uw = [strength_nr_uw, strength_nr_uw_t];


Cnt=vertcat(Cnt,Cent);
Cnt_uw=vertcat(Cnt_uw, Cent_uw);
Cnt_nr=vertcat(Cnt_nr,Cent_nr);
Cnt_nr_uw=vertcat(Cnt_nr_uw,Cent_nr_uw);
i

end;

% For braingraph demos
dataout.x=x';
dataout.y=y';
dataout.z=z2';
dataout.fripple_mi_matrix=fripple_mi_matrix;
dataout.fripple_nresected_mi_matrix=fripple_nresected_mi_matrix; 

% for data out
dataout.e_l=e_l;
dataout.e_l_uw=e_l_uw;
dataout.e_l_nr=e_l_nr;
dataout.e_l_uw_nr=e_l_uw_nr;

dataout.is=is';
dataout.is_nr=is_nr';
dataout.os=os';
dataout.os_nr=os_nr';
dataout.strength=strength';
dataout.strength_nr=strength_nr';
dataout.strength_uw=strength_uw';
dataout.strength_nr_uw=strength_nr_uw';

dataout.cC=cC;
dataout.cC_uw=cC_uw;
dataout.cC_nr=cC_nr;
dataout.cC_uw_nr=cC_uw_nr;

dataout.Cnt=Cnt;
dataout.Cnt_uw=Cnt_uw;
dataout.Cnt_nr=Cnt_nr;
dataout.Cnt_nr_uw=Cnt_nr_uw;

dataout.lambda=lambda';         
dataout.efficiency=efficiency';     
dataout.radius=radius';         
dataout.diameter=diameter';

dataout.uw_lambda=uw_lambda';         
dataout.uw_efficiency=uw_efficiency';     
dataout.uw_radius=uw_radius';         
dataout.uw_diameter=uw_diameter';

dataout.lambda_r=lambda_r';
dataout.efficiency_r=efficiency_r';
dataout.radius_r=radius_r';
dataout.diameter_r=diameter_r';

dataout.uw_lambda_r=uw_lambda_r';
dataout.uw_efficiency_r=uw_efficiency_r';
dataout.uw_radius_r=uw_radius_r';
dataout.uw_diameter_r=uw_diameter_r';

dataout.w_lambda=w_lambda';         
dataout.w_efficiency=w_efficiency';     
dataout.w_radius=w_radius';         
dataout.w_diameter=w_diameter';

dataout.resected_array=resected_all';
dataout.patient_array=patient_array';
dataout.responder_array=responder_array;
dataout.rate_array=rate_array';

dataout.mi_patient_array=mi_patient_array;
dataout.mi_resected_array=mi_resected_array
dataout.mi_array=mi_array;
dataout.mi_responder=mi_responder;