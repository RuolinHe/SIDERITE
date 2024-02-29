%% COCONUT new siderophore clustering
Dist_matrix=readcell('Tanimoto_COCONUT_SIDERITE.xlsx');
%% 
Header_index=Dist_matrix(2:end,1);
SMILES=Dist_matrix(2:end,2);
Dist_matrix(1,:)=[];Dist_matrix(:,[1,2])=[];Dist_matrix=cell2mat(Dist_matrix);
Dist = pdist(Dist_matrix);
Dist_tree = linkage(Dist_matrix,'average');
leafOrder = optimalleaforder(Dist_tree,Dist);
Dist_sorted=Dist_matrix(leafOrder,leafOrder);
Header_index_sorted=Header_index(leafOrder);
SMILES_sorted=SMILES(leafOrder);
% writematrix(Dist_sorted,'D:\课题组\zhiyuan_Lab\10-Database_resource\Program\Tanimoto_COCONUT_SIDERITE_sorted.xlsx','Sheet','Sheet1');
%% 
figure;
h = heatmap(Dist_sorted,'GridVisible','off');
h.XDisplayLabels = repmat({''},size(Dist_matrix,1),1);
h.YDisplayLabels = repmat({''},size(Dist_matrix,2),1);
colormap("jet")
title('Tanimoto distance clustering')
% saveas(gcf,'D:\课题组\zhiyuan_Lab\10-Database_resource\Progress\figure\New_siderophore_clustering.svg','svg');
%% 
figure;
hold on
for i = 1:length(Header_index_sorted)
    if strcmp(Header_index_sorted{i},'SIDERITE')
        patch([0,0.2,0.2,0],[-(i-1),-(i-1),-(i),-(i)],[0,1,0],'EdgeColor','none');
    else
        patch([0,0.2,0.2,0],[-(i-1),-(i-1),-(i),-(i)],[0,0,1],'EdgeColor','none');
    end
end
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
axis off
hold off
% saveas(gcf,'D:\课题组\zhiyuan_Lab\10-Database_resource\Progress\figure\New_siderophore_clustering_legend.svg','svg');