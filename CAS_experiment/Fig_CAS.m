% Figure CAS of siderophore experiments

% H20
h2o_raw1=readmatrix('CAS-data-2023.07.30.xlsx','Sheet','figure','Range','D62:S67');
h2o_raw2=readmatrix('CAS-data-2023.07.30.xlsx','Sheet','figure','Range','D68:G68');
h2o_control_raw=readmatrix('CAS-data-2023.07.30.xlsx','Sheet','figure','Range','I68:K68');

h2o_matrix=zeros(25,3);
h2o_id=zeros(25,1);
for i = 1:size(h2o_raw1,1)
    for j = 1:size(h2o_raw1,2)/4
        h2o_id((i-1)*size(h2o_raw1,2)/4+j)=h2o_raw1(i,4*j-3);
        h2o_matrix((i-1)*size(h2o_raw1,2)/4+j,:)=h2o_raw1(i,4*j-2:4*j);
    end
end
h2o_id(end)=h2o_raw2(1);
h2o_matrix(end,:)=h2o_raw2(2:end);
h2o_mean=mean(h2o_matrix,2);
h2o_std=std(h2o_matrix,0,2);
my_xticklabels_h2o=cell(length(h2o_id),1);
for i = 1:length(h2o_id)
    my_xticklabels_h2o{i}=['\bf ',num2str(h2o_id(i))];
end

f=figure;
f.Position(3)=800;
bar(1:length(h2o_id),h2o_mean)
hold on
errorbar(1:length(h2o_id),h2o_mean, h2o_std, 'k.', 'LineWidth', 1);
xticks(1:length(h2o_id));
xticklabels(my_xticklabels_h2o);
xtickangle(0)
xlabel('Compound')
ylabel('CAS')
my_xlim=xlim;
line(my_xlim,[mean(h2o_control_raw),mean(h2o_control_raw)],'Color','red','LineStyle','-')
line(my_xlim,[mean(h2o_control_raw)+std(h2o_control_raw),mean(h2o_control_raw)+std(h2o_control_raw)],'Color','red','LineStyle','--')
line(my_xlim,[mean(h2o_control_raw)-std(h2o_control_raw),mean(h2o_control_raw)-std(h2o_control_raw)],'Color','red','LineStyle','--')
title('Water-soluble siderophore candidates(n=25)')
hold off
saveas(gcf,'Experiment_H2O.svg','svg');
DMSO
DMSO_raw1=readmatrix('CAS-data-2023.07.30.xlsx','Sheet','figure','Range','D72:S78');
DMSO_raw2=readmatrix('CAS-data-2023.07.30.xlsx','Sheet','figure','Range','D79:k79');
DMSO_control_raw=readmatrix('CAS-data-2023.07.30.xlsx','Sheet','figure','Range','M79:O79');

DMSO_matrix=zeros(30,3);
DMSO_id=zeros(30,1);
for i = 1:size(DMSO_raw1,1)
    for j = 1:size(DMSO_raw1,2)/4
        DMSO_id((i-1)*size(DMSO_raw1,2)/4+j)=DMSO_raw1(i,4*j-3);
        DMSO_matrix((i-1)*size(DMSO_raw1,2)/4+j,:)=DMSO_raw1(i,4*j-2:4*j);
    end
end
DMSO_id(end-1)=DMSO_raw2(1);
DMSO_matrix(end-1,:)=DMSO_raw2(2:4);
DMSO_id(end)=DMSO_raw2(5);
DMSO_matrix(end,:)=DMSO_raw2(6:end);
DMSO_mean=mean(DMSO_matrix,2);
DMSO_std=std(DMSO_matrix,0,2);
my_xticklabels_DMSO=cell(length(DMSO_id),1);
for i = 1:length(DMSO_id)
    my_xticklabels_DMSO{i}=['\bf ',num2str(DMSO_id(i))];
end

f=figure;
f.Position(3)=1000;
bar(1:length(DMSO_id),DMSO_mean)
hold on
errorbar(1:length(DMSO_id),DMSO_mean, DMSO_std, 'k.', 'LineWidth', 1);
xticks(1:length(DMSO_id));
xticklabels(my_xticklabels_DMSO);
xtickangle(0)
xlabel('Compound')
ylabel('CAS')
my_xlim=xlim;
line(my_xlim,[mean(DMSO_control_raw),mean(DMSO_control_raw)],'Color','red','LineStyle','-')
line(my_xlim,[mean(DMSO_control_raw)+std(DMSO_control_raw),mean(DMSO_control_raw)+std(DMSO_control_raw)],'Color','red','LineStyle','--')
line(my_xlim,[mean(DMSO_control_raw)-std(DMSO_control_raw),mean(DMSO_control_raw)-std(DMSO_control_raw)],'Color','red','LineStyle','--')
title('Water-insoluble siderophore candidates(n=30, in DMSO)')
hold off
saveas(gcf,'Experiment_DMSO.svg','svg');
