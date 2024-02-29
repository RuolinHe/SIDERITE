%% sider ligand
Ligand_type_raw=readcell('D:\课题组\zhiyuan_Lab\10-Database_resource\Program\Sid_structure_output3.xlsx','Sheet','L_type','Range','A:B');
Ligand_type_raw(1,:)=[];
%%
Ligand_type_count=cell2mat(Ligand_type_raw(:,2));
[Ligand_type_count,I]=sort(Ligand_type_count,'descend');
Ligand_type_label=Ligand_type_raw(I,1);
%     Ligand_type_label=Ligand_type_label(1:find(ismember(Ligand_type_label,'Other')));
Ligand_type_label=Ligand_type_label(1:11); % TOP10
Ligand_type_count(length(Ligand_type_label))=sum(Ligand_type_count(length(Ligand_type_label):end));
Ligand_type_count(length(Ligand_type_label)+1:end)=[];
Ligand_type_label{end}='Others';

for i = 1:length(Ligand_type_label)
    Ligand_type_label{i}=strrep(Ligand_type_label{i},'/','+');
    Ligand_type_label{i}=[Ligand_type_label{i},' (',num2str(100*(round(Ligand_type_count(i)/sum(Ligand_type_count),4))),'%)'];
end
%%
% color_matrix=[255,235,205;255,255,0;152,251,152;221,160,221;0,255,0;0,206,209;184,134,11;255,0,255;255,20,147;178,34,34;47,79,79;0,100,0;0,0,255;138,138,138]./255;
figure;
% colormap(color_matrix);
colormap('parula')
pie(Ligand_type_count,Ligand_type_label); %Sider_ligand.svg
%% sider ligand contain
Ligand_type_raw=readcell('D:\课题组\zhiyuan_Lab\10-Database_resource\Program\Sid_structure_output3.xlsx','Sheet','L_type2','Range','B:C');
Ligand_type_raw(1,:)=[];
%%
Ligand_type_count=cell2mat(Ligand_type_raw(:,2));
[Ligand_type_count,I]=sort(Ligand_type_count,'descend');
Ligand_type_label=Ligand_type_raw(I,1);
% Ligand_type_count=Ligand_type_count./(size(result_tab,1)-1);
Ligand_type_label1 = categorical(Ligand_type_label);
Ligand_type_label1 = reordercats(Ligand_type_label1,Ligand_type_label);
figure;
hold on
% b=bar(Ligand_type_label1,Ligand_type_count*100);
b=bar(Ligand_type_label1,Ligand_type_count);
b.FaceColor='#0072BD';
xtips1 = b.XEndPoints;
ytips1 = b.YEndPoints;
labels1 = string(b.YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
% loc_xlim=xlim;
% ylim([0 size(result_tab,1)-1])
% ytickformat('percentage')
xlabel('Ligand type')
ylabel('Count in siderophores')
% line(loc_xlim,[size(result_tab,1)-1 size(result_tab,1)-1],'Color','red','LineStyle','--')
% xlim(loc_xlim)
% xtickangle(45)
hold off
%%
saveas(gcf,'D:\课题组\zhiyuan_Lab\10-Database_resource\Progress\figure\Sider_ligand_contain.svg','svg');
%% sider ligand sum
Ligand_type_raw=readcell('D:\课题组\zhiyuan_Lab\10-Database_resource\Program\Sid_structure_output3.xlsx','Sheet','L_type3','Range','B:C');
Ligand_type_raw(1,:)=[];
%%
Ligand_type_count=cell2mat(Ligand_type_raw(:,2));
[Ligand_type_count,I]=sort(Ligand_type_count,'descend');
Ligand_type_label=Ligand_type_raw(I,1);
% Ligand_type_count=Ligand_type_count./(size(result_tab,1)-1);
Ligand_type_label1 = categorical(Ligand_type_label);
Ligand_type_label1 = reordercats(Ligand_type_label1,Ligand_type_label);
figure;
% b=bar(Ligand_type_label1,Ligand_type_count*100);
b=bar(Ligand_type_label1,Ligand_type_count);
b.FaceColor='#0072BD';
xtips1 = b.XEndPoints;
ytips1 = b.YEndPoints;
labels1 = string(b.YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylim([0 1050])
% ytickformat('percentage')
xlabel('Ligand type')
ylabel('Sum of ligand')
saveas(gcf,'D:\课题组\zhiyuan_Lab\10-Database_resource\Progress\figure\Sider_ligand_sum.svg','svg');
%% sider source
Ligand_type_raw=readcell('D:\课题组\zhiyuan_Lab\10-Database_resource\Program\Sid_structure_output3.xlsx','Sheet','Source','Range','A:B');
Ligand_type_raw(1,:)=[];
%%
Ligand_type_count=cell2mat(Ligand_type_raw(:,2));
[Ligand_type_count,I]=sort(Ligand_type_count,'descend');
Ligand_type_label=Ligand_type_raw(I,1);
for i = 1:length(Ligand_type_label)
    Ligand_type_label{i}=[Ligand_type_label{i},' (',num2str(100*(round(Ligand_type_count(i)/sum(Ligand_type_count),4))),'%)'];
end
% color_matrix=[0.8500 0.3250 0.0980;0.3010 0.7450 0.9330;0 100/255 0;1 1 0];
color_matrix=[255,165,0;255,20,147;0,255,0;0,128,128]./255;
figure;
colormap(color_matrix);
pie(Ligand_type_count,Ligand_type_label);
%%
saveas(gcf,'D:\课题组\zhiyuan_Lab\10-Database_resource\Progress\figure\Sider_source.svg','svg');
%% Biosynthetic Type
Ligand_type_raw=readcell('D:\课题组\zhiyuan_Lab\10-Database_resource\Program\Sid_structure_output3.xlsx','Sheet','B_Type','Range','E11:F15');
% Ligand_type_raw(1,:)=[];
%%
Ligand_type_count=cell2mat(Ligand_type_raw(:,2));
[Ligand_type_count,I]=sort(Ligand_type_count,'descend');
Ligand_type_label=Ligand_type_raw(I,1);
for i = 1:length(Ligand_type_label)
    Ligand_type_label{i}=[Ligand_type_label{i},' (',num2str(100*(round(Ligand_type_count(i)/sum(Ligand_type_count),4))),'%)'];
end
% color_matrix=[61,139,93;213,0,57;240,162,104;255,254,138;177,197,220]./255;
color_matrix=[255,20,147;0,0,255;255,165,0;0,128,128;0,255,0]./255;
figure;
colormap(color_matrix);
pie(Ligand_type_count,Ligand_type_label);
colormap(color_matrix(1:length(Ligand_type_count),:));
%%
saveas(gcf,'D:\课题组\zhiyuan_Lab\10-Database_resource\Progress\figure\Sider_Biosynthetic_Type.svg','svg');
%% Theoretical denticity
Denticity_index = ismember(result_tab(1,:),'Theoretical denticity');
figure;
histogram(cell2mat(result_tab(2:end,Denticity_index)));
xlabel('Denticity number')
ylabel('Count in siderophores')
saveas(gcf,'D:\课题组\zhiyuan_Lab\10-Database_resource\Progress\figure\Sider_denticity.svg','svg');
%% MW
MW_index = ismember(result_tab(1,:),'Molecular Weight');
figure;
histogram(cell2mat(result_tab(2:end,MW_index)));
xlabel('Molecular weight (Da)')
ylabel('Count in siderophores')
saveas(gcf,'D:\课题组\zhiyuan_Lab\10-Database_resource\Progress\figure\Sider_MW.svg','svg');
%% Theoretical denticity vs MW.
Denticity = cell2mat(loc_result_tab(:,Denticity_index));
MW = cell2mat(loc_result_tab(:,MW_index));
%%
figure;
boxplot(MW,Denticity);
xlabel('Denticity number')
ylabel('Molecular weight (Da)')
saveas(gcf,'D:\课题组\zhiyuan_Lab\10-Database_resource\Progress\figure\Denticity_vs_MW.svg','svg');
%% logS
logS=readmatrix('D:\课题组\zhiyuan_Lab\10-Database_resource\Program\Sid_structure_output3.xlsx','Sheet','Unique1','Range','BU:BU'); % 改
figure;
histogram(logS(2:end,:));
xlabel('Predicted aqueous solubility (logS)')
ylabel('Count in siderophores')
saveas(gcf,'D:\课题组\zhiyuan_Lab\10-Database_resource\Progress\figure\Sider_logS.svg','svg');
%% Predicted diffusion coefficient (10^(-10) m2 s-1) in H2O, 298.15K
diff_coef=readmatrix('D:\课题组\zhiyuan_Lab\10-Database_resource\Program\Sid_structure_output3.xlsx','Sheet','Unique1','Range','BV:BV'); % 改
%% 
figure;
histogram(diff_coef(2:end,:));
xlabel('Predicted diffusion coefficient (10^{-10} m^2 s^{-1})')
ylabel('Count in siderophores')
saveas(gcf,'D:\课题组\zhiyuan_Lab\10-Database_resource\Progress\figure\Sider_diff_coef.svg','svg');