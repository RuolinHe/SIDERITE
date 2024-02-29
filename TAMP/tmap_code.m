%% load data
Sid_raw=readcell('D:\课题组\zhiyuan_Lab\10-Database_resource\Program\Sid_structure_output3.xlsx','Sheet','Unique1','Range','A:BV'); %改
%% SID
SID=cell(size(Sid_raw,1)-1,1);
for i = 1:size(Sid_raw,1)-1
    SID{i}=['SID',num2str(cell2mat(Sid_raw(i+1,1)),'%05d')];
end
%% Biosynthetic Type
Bio_type_index = find(ismember(Sid_raw(1,:),'Biosynthetic Type'));
sider_bio_type=Sid_raw(2:end,Bio_type_index);
for i = 1:length(sider_bio_type)
    sider_bio_type{i}=strrep(sider_bio_type{i},'Putative ','');
end
%% Kingdom
Kingdom=cell(size(Sid_raw,1)-1,1);
Bacteria_index = find(ismember(Sid_raw(1,:),'Bacteria'));
Animal_index = find(ismember(Sid_raw(1,:),'Animal'));
Kingdom_type=Sid_raw(1,Bacteria_index:Animal_index);
for i = 1:length(Kingdom)
    Kingdom(i)=join(Kingdom_type(cell2mat(Sid_raw(i+1,Bacteria_index:Animal_index))>0),'/');
end
%% Precursor
Precursor=cell(size(Sid_raw,1)-1,1);
Source1_index = find(ismember(Sid_raw(1,:),'Hydroxamate source'));
Source2_index = find(ismember(Sid_raw(1,:),'2-Nitrosophenol source'));
for i = 1:length(Precursor)
    loc_index=[];
    for j = Source1_index:Source2_index
        if ~ismissing(Sid_raw{i+1,j})
            loc_index=[loc_index;j];
        end
    end
    Precursor(i)=join(Sid_raw(i+1,loc_index),'/');
    Precursor{i}=strrep(Precursor{i},'; ','/');
    if isempty(Precursor{i})
        Precursor{i}='Unknown';
    end
end
%% Siderophore Class
Large_index=find(ismember(Sid_raw(1,:),'Siderophore Large Class'));
Sid_Class=cell(size(Sid_raw,1)-1,1);
for i = 1:length(Sid_Class)
    Sid_Class{i}=[num2str(Sid_raw{i+1,Large_index}),'.',num2str(Sid_raw{i+1,Large_index+1}),'.',num2str(Sid_raw{i+1,Large_index+2})];
end
%% Siderophore Small Class
Sid_Small_Class=cell(size(Sid_raw,1)-1,1);
for i = 1:length(Sid_Small_Class)
    Sid_Small_Class{i}=[num2str(Sid_raw{i+1,Large_index}),'.',num2str(Sid_raw{i+1,Large_index+1})];
end
%% Ligand Type
Ligand_index=find(ismember(Sid_raw(1,:),'Ligand Type'));
Ligand_Type=cell(size(Sid_raw,1)-1,1);
for i = 1:length(Ligand_Type)
    Ligand_Type{i}=strrep(Sid_raw{i+1,Ligand_index},'/','/__');
end
%% in chemistry space
SID_URL='https://www.test.org/';
SID_Con_SMILES=readcell('D:\课题组\zhiyuan_Lab\10-Database_resource\Program\Sid_structure_unique_Canonical.xlsx','Range','D:D');
SID_Con_SMILES(1)=[];
file_name={'NPAtlas','COCONUT','LOTUS','SuperNatural','NPASS','chemdiv'};
for i = 1:length(file_name)
    if ~exist(['./chem_space/',file_name{i}],'dir')
        mkdir(['./chem_space/',file_name{i}])
    end
end
SID_overlap_list=[];
%%
SMILES_index=find(ismember(Sid_raw(1,:),'Canonical SMILES'));
Denticity_index=find(ismember(Sid_raw(1,:),'Theoretical denticity'));
Count_index=find(ismember(Sid_raw(1,:),'Count'));
Monomer_index=find(ismember(Sid_raw(1,:),'monomers Number'));
MW_index=find(ismember(Sid_raw(1,:),'Molecular Weight'));
Phylum_index=find(ismember(Sid_raw(1,:),'Phylum'));
Hydroxamate_index=find(ismember(Sid_raw(1,:),'Hydroxamate'));
Nitrosophenol_index=find(ismember(Sid_raw(1,:),'2-Nitrosophenol'));
%% Predicted diffusion coefficient 
Pre_D = cell(size(Sid_raw,1)-1,1);
for i = 2:size(Sid_raw,1)
    Pre_D{i-1}=segwe_d(cell2mat(Sid_raw(i,MW_index)));
end
%% 
result_tab3=[Sid_raw(1,[1:3,SMILES_index,Bio_type_index]),{'Kingdom'},Sid_raw(1,[Denticity_index,Count_index,Monomer_index,Monomer_index+1,MW_index]),{'Ligand Type','Precursor'},Sid_raw(1,[Large_index,Phylum_index,Hydroxamate_index:Nitrosophenol_index]),{'Siderophore Class','Siderophore Small Class','Predicted logS','Predicted diffusion coefficient'}];
result_tab3=[result_tab3;SID,Sid_raw(2:end,2:3),SID_Con_SMILES,sider_bio_type,Kingdom,Sid_raw(2:end,[Denticity_index,Count_index,Monomer_index,Monomer_index+1,MW_index]),Ligand_Type,Precursor,Sid_raw(2:end,[Large_index,Phylum_index,Hydroxamate_index:Nitrosophenol_index]),Sid_Class,Sid_Small_Class,Sid_raw(2:end,end-1),Pre_D];
for i = 1:size(result_tab3,1)
    if ismissing(result_tab3{i,3}) % Siderophore other name
        result_tab3{i,3}='None';
    end
    if isempty(result_tab3{i,15}) % Phylum
        result_tab3{i,15}='Other';
    end
end
%%
writecell(result_tab3,'D:\课题组\zhiyuan_Lab\10-Database_resource\Program\sider_tmap\Sid_tmap20230306.csv','Delimiter','tab');
%% COCONUT
NPA_raw=readcell('D:\课题组\zhiyuan_Lab\10-Database_resource\Program\COCONUT4MetFrag_Canonical.xlsx','Range','A:D');
NPA_raw(1,:)=[];
NPA_ID=NPA_raw(:,2);
NPA_SMILES=NPA_raw(:,1);
NPA_Valid=cell2mat(NPA_raw(:,3));
NPA_Con_SMILES=NPA_raw(:,4);
for i = 1:length(NPA_Con_SMILES)
    if ismissing(NPA_Con_SMILES{i})
        NPA_Con_SMILES{i}='';
    end
end
% remove the same SMILES
SID_overlap_list=[SID_overlap_list,ismember(SID_Con_SMILES,NPA_Con_SMILES)];
remove_list=ismember(NPA_Con_SMILES,SID_Con_SMILES);
fprintf('%d siderophores in %s are known in our database\n',sum(remove_list),file_name{2})% 326
fprintf('%d siderophores in %s are invalid\n',length(NPA_Valid)-sum(NPA_Valid),file_name{2})% 351
remove_list=remove_list|NPA_Valid==0;
NPA_SMILES(remove_list)=[];
NPA_ID(remove_list)=[];
NPA_Con_SMILES(remove_list)=[];
NPA_URL='https://coconut.naturalproducts.net/compound/coconut_id/';
Header=[{'ID','SMILES','URL','In_database','Siderophore Large Class'},result_tab3(1,[29,2,3,11,6,5,7,9,13,12])];
Chem_raw=[Header;[[NPA_ID;SID],[NPA_Con_SMILES;SID_Con_SMILES],[repmat({NPA_URL},length(NPA_ID),1);repmat({SID_URL},length(SID),1)],num2cell([zeros(length(NPA_ID),1);ones(length(SID),1)])],[repmat({'0'},length(NPA_ID),1);Sid_raw(2:end,Large_index)],[[repmat({' '},length(NPA_ID),3),repmat({'0'},length(NPA_ID),1),repmat({' '},length(NPA_ID),6)];result_tab3(2:end,[29,2,3,11,6,5,7,9,13,12])]];
% writecell(Chem_raw,'D:\课题组\zhiyuan_Lab\10-Database_resource\Program\sider_tmap\chem_space\COCONUT\COCONUT.csv','Delimiter','tab');
%% COCONUT for new siderophore discovery
% Specify the file name
fileName = 'COCONUT_r.txt';

% Open the file for writing
fileID = fopen(fileName, 'w');

% Check if the file was opened successfully
if fileID == -1
    error('Could not open the file for writing.');
end

% Write each cell string to the file
for i = 1:length(NPA_SMILES)
    fprintf(fileID, '%s\n', NPA_SMILES{i});
end

% Close the file
fclose(fileID);

P_NPA_SIMLES=readcell('D:\课题组\zhiyuan_Lab\10-Database_resource\Program\output\summary.xlsx','Sheet','COCONUT-01','Range','B:B');
P_NPA_SIMLES(1)=[];
P_NPA_ID=cell(length(P_NPA_SIMLES),1);
P_NPA_Link=cell(length(P_NPA_SIMLES),1);
for i = 1:length(P_NPA_SIMLES)
    P_NPA_ID(i)=NPA_ID(ismember(NPA_SMILES,P_NPA_SIMLES{i}));
    P_NPA_Link{i}=[NPA_URL,P_NPA_ID{i}];
end
%% function
function answer = segwe_d(argument)
    % modified from segwe_dmw_v2
    % https://www.nmr.chemistry.manchester.ac.uk/?q=node/432
    % http://dx.doi.org/10.1021/acs.analchem.7b05032
    % argument is MW in g/mol
    % answer is predicted diffusion coefficient in 10^(-10) m^2 s-1
    % Constant
    kB = 1.380644e-23;
    NA = 6.02214e23;
    T = 298.15; % 室温
    
    % Water
    MWs = 18.01;
    a = -14.0236;
    b = 2092.95;
    peff=627;

    dosydata.temperature=T;

    n=exp(a+(b/dosydata.temperature));

    dosydata.argument=argument;
    
    dosydata.viscositydisplay=n*1000;
    dosydata.MWs=MWs;
    dosydata.peff=peff;
    dosydata.packing=1;
    dosydata.density=dosydata.peff/1000;
    
%     radius=(((dosydata.argument/(10^(3)))/((4*pi/3)*NA*dosydata.peff))^(1/3))/(10^-12); % Hydrodynamic radius / (10^-12 m)
    alpha = ((dosydata.MWs/dosydata.argument)^(1/3));
    dosydata.answer = ((kB*dosydata.temperature*(3*alpha/2 + 1/(1+alpha)))/(6*pi*n*(3*dosydata.argument*10^(-3)/(4*pi*dosydata.peff*NA))^(1/3)))/(10^-10);
    answer = dosydata.answer;
end