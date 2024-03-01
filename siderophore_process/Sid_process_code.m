%%
Last_row_num = 873;
sider_Molecular_Formula = deblank(readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range',['AB2:AB',num2str(Last_row_num)]));
Element = 'CHNOS';
Ele_matrix = zeros(length(sider_Molecular_Formula),5);
for i = 1:length(sider_Molecular_Formula)
    loc_formular = strrep(sider_Molecular_Formula{i},'+','');
    loc_formular = strrep(loc_formular,'Cl','');
    Ele_list = [];
    for j = 1:length(Element)
       if contains(loc_formular,Element(j))
           Ele_list = [Ele_list;Element(j)];
       end
    end
    for j = 1:length(Ele_list)-1
        loc_index = strfind(loc_formular,Ele_list(j));
        next_index = strfind(loc_formular,Ele_list(j+1));
        if next_index - loc_index==1
           Ele_matrix(i,strfind(Element, Ele_list(j)))=1;
        else
            Ele_matrix(i,strfind(Element, Ele_list(j)))=str2double(loc_formular(loc_index+1:next_index-1));
        end
    end
    last_index = strfind(loc_formular,Ele_list(end));
    if last_index == length(loc_formular)
        Ele_matrix(i,strfind(Element, Ele_list(end)))=1;
    else
        Ele_matrix(i,strfind(Element, Ele_list(end)))=str2double(loc_formular(last_index+1:end));
    end
end
%% siderophore ligand type
sider_ligand_type_matrix = readmatrix('Sid_structure_output3.xlsx','Sheet','1-2','Range',['F2:R',num2str(Last_row_num)]);
sider_ligand_type_matrix_all = readmatrix('Sid_structure_output3.xlsx','Sheet','1-2','Range',['F2:T',num2str(Last_row_num)]);
sider_ref = readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range',['X2:X',num2str(Last_row_num)]);
sider_ref_bio_type = readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range',['AH2:AH',num2str(Last_row_num)]);
sider_ligand_type_header_all = deblank(readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range','F1:T1'));
sider_ligand_type_header = deblank(readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range','F1:R1'));
sider_ligand_source_header = readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range','AX1:BE1');
sider_record_source = readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range',['AG2:AG',num2str(Last_row_num)]);
sider_record_source_tab = tabulate(sider_record_source);
% writecell([{'Source','Count','Frequency (%)'};sider_record_source_tab],'Sid_structure_output3.xlsx','Sheet','C_source','Range','A:C');
sider_ligand_type_header{5}=strrep(sider_ligand_type_header{5},' in Citrate','');
sider_ligand_type_header{11}=strrep(sider_ligand_type_header{11},' in Citrate','');
sider_ligand_type=cell(size(sider_ligand_type_matrix,1),1);
for i = 1:size(sider_ligand_type_matrix,1)
    if sum(sider_ligand_type_matrix(i,:))==0
        sider_ligand_type{i} = 'Other';
    else
        loc_sider_ligand_type = unique(sider_ligand_type_header(sider_ligand_type_matrix(i,:)~=0));
        sider_ligand_type(i) = join(loc_sider_ligand_type,'/');
    end
end
for i = 1:length(sider_ref_bio_type)
    if ismissing(sider_ref_bio_type{i})
        sider_ref_bio_type{i}=[];
    end
end
%% siderophore SMILES
sider_name = deblank(readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range',['D2:D',num2str(Last_row_num)]));
sider_other_name = readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range',['E2:E',num2str(Last_row_num)]);
sider_con_SMILES = readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range',['Z2:Z',num2str(Last_row_num)]);
sider_bio_type = readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range',['AA2:AA',num2str(Last_row_num)]);
sider_microbe = deblank(readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range',['V2:V',num2str(Last_row_num)]));
sider_Rich_infor_matrix = readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range',['AI2:BP',num2str(Last_row_num)]);
sider_Molecular_Weight = readcell('Sid_structure_output3.xlsx','Sheet','1-2','Range',['AD2:AD',num2str(Last_row_num)]);
%%
for i = 1:size(sider_Rich_infor_matrix,1)
   for j = [1:4,16:23,33,34] % 搜索改
      if ismissing(sider_Rich_infor_matrix{i,j})
          sider_Rich_infor_matrix{i,j} = '';
      end
   end
end
sider_microbe_phylum_genus = join(sider_Rich_infor_matrix(:,[9,10,14]),'/'); % 改
index = 0;
index_list = [];
last_str = '';
for i = 1:length(sider_name)
    loc_microbe = split(sider_microbe{i});
    if ~strcmp(sider_name{i},last_str)
       last_str = sider_name{i};
       index = index + 1;
    end
    index_list = [index_list;index];
end
sider_microbe_genus_phylum_tab = tabulate(sider_microbe_phylum_genus);
[uni_sider_name,sider_name_ia,sider_name_ic] = unique(sider_name);
uni_sider_other_name=sider_other_name(sider_name_ia);
assert(length(unique(sider_con_SMILES(sider_name_ia)))==length(sider_name_ia)) % check smiles are unique
for i = 1:length(uni_sider_other_name)
    if ismissing(uni_sider_other_name{i})
        uni_sider_other_name{i}=[];
    end
end
uni_sider_name_count = zeros(length(uni_sider_name),1);
uni_sider_name_microbe = cell(length(uni_sider_name),1);
uni_sider_phylum_str = cell(length(uni_sider_name),1);
uni_sider_phylum_header={'Actinobacteria','Cyanobacteria','Firmicutes','Proteobacteria','Bacteroidetes','Ascomycota','Basidiomycota','Mucoromycota','Streptophyta','Chordata'};
uni_sider_phylum_freq = zeros(length(uni_sider_name),length(uni_sider_phylum_header));
uni_sider_kingdom_header = {'Bacteria','Fungi','Plant','Animal'};
uni_sider_kingdom_freq = zeros(length(uni_sider_name),length(uni_sider_kingdom_header));
for i = 1:length(uni_sider_name)
    loc_index=find(sider_name_ic == i);
    uni_sider_name_count(i) = sum(sider_name_ic == i);
    uni_sider_name_microbe(i) = join(sider_microbe(sider_name_ic == i),'; ');
    loc_Phylum_tab=tabulate(sider_Rich_infor_matrix(sider_name_ic == i,10)); % 改
    uni_sider_phylum_str(i)=join(loc_Phylum_tab(:,1),'/');
    for j = 1:length(uni_sider_phylum_header)
        uni_sider_phylum_freq(i,j)=sum(ismember(sider_Rich_infor_matrix(sider_name_ic == i,10),uni_sider_phylum_header{j}))/uni_sider_name_count(i); % 改
    end
    uni_sider_kingdom_freq(i,4) = uni_sider_phylum_freq(i,ismember(uni_sider_phylum_header,'Chordata'));
    uni_sider_kingdom_freq(i,3) = uni_sider_phylum_freq(i,ismember(uni_sider_phylum_header,'Streptophyta'));
    uni_sider_kingdom_freq(i,2) = sum(uni_sider_phylum_freq(i,ismember(uni_sider_phylum_header,'Ascomycota')|ismember(uni_sider_phylum_header,'Basidiomycota')|ismember(uni_sider_phylum_header,'Mucoromycota')));
    uni_sider_kingdom_freq(i,1) = 1-sum(uni_sider_kingdom_freq(i,2:4));
end
uni_sider_name_ligand_type = sider_ligand_type(sider_name_ia);
uni_sider_name_ligand_type_tab = tabulate(uni_sider_name_ligand_type);
uni_sider_name_bio_type = sider_bio_type(sider_name_ia);
uni_sider_name_bio_type_tab = tabulate(uni_sider_name_bio_type);
result_tab = [{'Siderophore ID','Siderophore name','Siderophore other name','Count'},sider_ligand_type_header_all,{'Ligand Type','Biosynthetic Type','Canonical SMILES','monomers Name', 'monomers Number', 'coverage', 'missing Monomers Number','Molecular Formula','Molecular Weight','Microorganisms'},uni_sider_kingdom_header,uni_sider_phylum_header,'Phylum',sider_ligand_source_header,'References','Biosynthetic Type Reference','Fe3+ Affinity Reference','logKf Fe3+','pFe3+','Redox potential','C','H','N','O','S','C/H','C/N','C/O','C/S','Biosynthetic Gene','Receptor Gene'];
% result_tab=[result_tab;num2cell([1:length(uni_sider_name)]'),uni_sider_name,uni_sider_other_name,num2cell(uni_sider_name_count),num2cell(sider_ligand_type_matrix_all(sider_name_ia,:)),uni_sider_name_ligand_type,uni_sider_name_bio_type,sider_con_SMILES(sider_name_ia),sider_Rich_infor_matrix(sider_name_ia,1:4),sider_Molecular_Formula(sider_name_ia),sider_Molecular_Weight(sider_name_ia),uni_sider_name_microbe,num2cell(uni_sider_kingdom_freq),num2cell(uni_sider_phylum_freq),uni_sider_phylum_str,sider_Rich_infor_matrix(sider_name_ia,12:17),uni_sider_ref,uni_sider_bio_ref,sider_Rich_infor_matrix(sider_name_ia,18:26)];
result_tab=[result_tab;num2cell([1:length(uni_sider_name)]'),uni_sider_name,uni_sider_other_name,num2cell(uni_sider_name_count),num2cell(sider_ligand_type_matrix_all(sider_name_ia,:)),uni_sider_name_ligand_type,uni_sider_name_bio_type,sider_con_SMILES(sider_name_ia),sider_Rich_infor_matrix(sider_name_ia,5:8),sider_Molecular_Formula(sider_name_ia),sider_Molecular_Weight(sider_name_ia),uni_sider_name_microbe,num2cell(uni_sider_kingdom_freq),num2cell(uni_sider_phylum_freq),uni_sider_phylum_str,sider_Rich_infor_matrix(sider_name_ia,16:23),sider_ref(sider_name_ia),sider_ref_bio_type(sider_name_ia),sider_Rich_infor_matrix(sider_name_ia,[4,1:3,end-10:end])];
% writecell(result_tab,'Sid_structure_output3.xlsx','Sheet','Unique');
uni_sider_ligand_type_matrix = sider_ligand_type_matrix(sider_name_ia,:);
uni_sider_ligand_type_num = zeros(size(uni_sider_ligand_type_matrix,2),1);
uni_sider_ligand_type_num_count = zeros(size(uni_sider_ligand_type_matrix,2),1);
uni_sider_ligand_type_num_count1 = zeros(size(uni_sider_ligand_type_matrix,2),1);
Bio_type_index = ismember(result_tab(1,:),'Biosynthetic Type');
for i = 1:size(uni_sider_ligand_type_matrix,2)
    uni_sider_ligand_type_num(i) = sum(uni_sider_ligand_type_matrix(:,i)~=0); % Ligand type (including)
    uni_sider_ligand_type_num_count(i) = sum(uni_sider_ligand_type_matrix(uni_sider_ligand_type_matrix(:,i)~=0,i)); % Ligand type (sum)
    uni_sider_ligand_type_num_count1(i) = sum(uni_sider_ligand_type_matrix((uni_sider_ligand_type_matrix(:,i)~=0)&(ismember(result_tab(2:end,Bio_type_index),'NRPS')|ismember(result_tab(2:end,Bio_type_index),'Putative NRPS')),i)); % Ligand type (sum) in NRPS
end
%% Distance matrix
% Tc=readcell('Tanimoto.xlsx');
Dc=readcell('Dice.xlsx');
% Tc(1,:)=[];Tc(:,1)=[];Tc=cell2mat(Tc);
Dc(1,:)=[];Dc(:,1)=[];Dc=cell2mat(Dc);
%%
Dist = pdist(Dc);
Dc_tree = linkage(Dc,'average');
leafOrder = optimalleaforder(Dc_tree,Dist);
Dc_sorted=Dc(leafOrder,leafOrder);
figure;
h = heatmap(Dc_sorted,'GridVisible','off');
h.XDisplayLabels = repmat({''},size(Dc,1),1);
h.YDisplayLabels = repmat({''},size(Dc,2),1);
colormap("parula")
title('Dice (average)')
Dc_cluster = cluster(Dc_tree,'Cutoff',0.6);
cgo=clustergram(Dc_sorted,'colormap','parula','linkage','average','Symmetric','false');
colorbar;
%%
% result_tab=[uni_sider_name,num2cell(uni_sider_name_count),num2cell(sider_ligand_type_matrix_all(sider_name_ia,:)),uni_sider_name_ligand_type,uni_sider_name_bio_type,sider_con_SMILES(sider_name_ia),sider_Rich_infor_matrix(sider_name_ia,[10:13,1:9]),sider_Molecular_Formula(sider_name_ia),sider_Molecular_Weight(sider_name_ia),uni_sider_name_microbe,num2cell(uni_sider_kingdom_freq),num2cell(uni_sider_phylum_freq),sider_Rich_infor_matrix(sider_name_ia,[21:26])];
result_tab_header=result_tab(1,:);
result_tab1=result_tab;
result_tab1(1,:)=[];
result_tab1=result_tab1(leafOrder,:);
result_tab1(:,1)=[];
result_tab1=[num2cell([1:length(uni_sider_name)]'),result_tab1];
result_tab1=[result_tab_header;result_tab1];
% result_tab1 = [{'Siderophore ID','Siderophore name','Count'},sider_ligand_type_header_all,{'Ligand Type','Biosynthetic Type','Canonical SMILES','monomers Name', 'monomers Number', 'coverage', 'missing Monomers Number','Molecular Formula','Molecular Weight','Microorganisms'},uni_sider_kingdom_header,uni_sider_phylum_header,sider_ligand_type_header_all([1,2,6,7,8,10]),'C','H','N','O','S','C/H','C/N','C/O','C/S';result_tab];
% writecell(result_tab1,'Sid_structure_output3.xlsx','Sheet','Unique');
%% siderophore cluster
% result_tab=result_tab1;
result_tab=readcell('Sid_structure_output3.xlsx','Sheet','Unique1','Range','A:BQ'); % 改
for i = 1:size(result_tab,1)
    for j = 1:size(result_tab,2)
        if ismissing(result_tab{i,j})
            result_tab{i,j}=[];
        end
    end
end
Sid_type=readcell('Sid_structure_output3.xlsx','Sheet','Unique1','Range','BR:BS'); % 改
Sid_type(1,:)=[];
Sid_type=cell2mat(Sid_type);
%%
MW_index = ismember(result_tab(1,:),'Molecular Weight');
Ligand_index = ismember(result_tab(1,:),'Ligand Type');
Bio_type_index = ismember(result_tab(1,:),'Biosynthetic Type');
Source=cell2mat(result_tab(2:end,find(ismember(result_tab(1,:),'Bacteria')):find(ismember(result_tab(1,:),'Animal')))); 
Source_type=result_tab(1,find(ismember(result_tab(1,:),'Bacteria')):find(ismember(result_tab(1,:),'Animal')));
uni_sider_name_count_sorted=cell2mat(result_tab(2:end,ismember(result_tab(1,:),'Count')));
Sid_type_uni=unique(Sid_type,'rows','stable');
Sid_type_uni_n=zeros(length(Sid_type_uni),1);
Sid_type_uni_total_count=zeros(length(Sid_type_uni),1);
Source_type_uni=cell(length(Sid_type_uni),1);
Sid_type_result_tab=result_tab(1,:);
for i = 1:length(Sid_type_uni)
% for i = 12
    loc_index=find(Sid_type(:,1)==Sid_type_uni(i,1)&Sid_type(:,2)==Sid_type_uni(i,2));
    loc_Source=sum(Source(loc_index,:),1);
    loc_Source_type=Source_type(loc_Source>0);
    loc_Bio_type=result_tab(1+loc_index,Bio_type_index);
    for j = 1:length(loc_Bio_type)
        loc_Bio_type{j}=strrep(loc_Bio_type{j},'Putative ','');
    end
    loc_Bio_type_uni=unique(loc_Bio_type);
    loc_Bio_type_uni=join(loc_Bio_type_uni,'; ');
    Source_type_uni(i)=join(loc_Source_type,'/');
    Sid_type_uni_n(i)=length(loc_index);
    Sid_type_uni_total_count(i)=sum(uni_sider_name_count_sorted(loc_index));
    loc_Dc_sorted=Dc_sorted(loc_index,loc_index);
    [~,I]=max(sum(loc_Dc_sorted));
    Sid_type_result_tab=[Sid_type_result_tab;result_tab(1+loc_index(I),:)];
    Sid_type_result_tab(end,21)=loc_Bio_type_uni;
end
Sid_type_result_tab=[Sid_type_result_tab,[{'Siderophore Large Class','Siderophore Class'};num2cell(Sid_type_uni)],[{'Source in Siderophore Class'};Source_type_uni],[{'Count in Siderophore Class'};num2cell(Sid_type_uni_n)],[{'All Counts in Siderophore Class'};num2cell(Sid_type_uni_total_count)]];
Sid_type_uni_c=zeros(size(Sid_type_uni,1),1);
for i = 1:length(Sid_type_uni_c)
    Sid_type_uni_c(i)=Sid_type_uni(i,1)+0.01*Sid_type_uni(i,2);
end
[~,I]=sort(Sid_type_uni_c);
Sid_type_result_tab1=Sid_type_result_tab(2:end,:);
Sid_type_result_tab=[Sid_type_result_tab(1,:);Sid_type_result_tab1(I,:)];
% writecell(Sid_type_result_tab,'Sid_structure_output3.xlsx','Sheet','S_Class');
Sid_large_type_uni=unique(Sid_type_uni(:,1));
Sid_large_type_uni_n=zeros(length(Sid_large_type_uni),1);
Sid_large_type_uni_total_count=zeros(length(Sid_large_type_uni),1);
for i = 1:length(Sid_large_type_uni)
    Sid_large_type_uni_n(i)=sum(Sid_type_uni_n(Sid_type_uni(:,1)==Sid_large_type_uni(i)));
    Sid_large_type_uni_total_count(i)=sum(Sid_type_uni_total_count(Sid_type_uni(:,1)==Sid_large_type_uni(i)));
end
Sid_large_type_uni_n_sorted=sort(Sid_large_type_uni_n,'descend');
Sid_large_type_uni_total_count_sorted=sort(Sid_large_type_uni_total_count,'descend');