clc
clear
close all;
%%
% reg_folder = '/run/user/1001/gvfs/smb-share:server=172.16.19.80,share=proj4/TracerX/SqCCs_fi_20x/serial_3pairs/data/Reg_manual';
reg_folder = '/Users/hzhang/Documents/project/sum/CellPos_all/HE_r/alignment_acc';
reg_files = dir(fullfile(reg_folder, '*cell.ndpi'));
%%
npoint = [];
sum_d = [];
rTRE = [];
mtre = [];
sdtre = [];
for i = 1:length(reg_files)
    reg_files(i).name
    Vars = {'fixedPoints','movingPoints','tform','ref'};
    S = load(fullfile(reg_folder, reg_files(i).name, 'reg.mat'), Vars{:});
    %S = load(fullfile(reg_folder, reg_files(i).name, 'tform.mat'), Vars{:});
    
%     pre_points = transformPointsForward(S.tform, S.movingPoints);
    pre_points = S.movingPoints;
    diagonal = (size(S.ref.Ss1,1)^2 + size(S.ref.Ss1,2)^2)^0.5;
    dis = [];

    for r = 1:size(pre_points,1)  
        d = ((pre_points(r,1) - S.fixedPoints(r,1))^2 + ...
            (pre_points(r,2) - S.fixedPoints(r,2))^2)^0.5;
        tre = d/diagonal;
        dis = cat(1,dis,d);
    end
    npoint = cat(1,npoint,length(S.movingPoints));
    sum_d = cat(1,sum_d, sum(dis));
    mtre = cat(1,mtre,mean(dis/diagonal));
    sdtre = cat(1,sdtre,std(dis/diagonal));
    rTRE = cat(1, rTRE, dis/diagonal);

end
%%
sumdf = [string(cat(1,reg_files.name)),npoint,sum_d, mtre,sdtre];
df = array2table(sumdf, 'VariableNames', {'Slide','npoint','sum_d','mean_rTRE','SD_rTRE'});
writetable(df,fullfile('/Users/hzhang/Documents/project/sum/CellPos_all/HE_r/alignment_acc/','reg_accuracy_T_new_1.csv'),...
'Delimiter',',','QuoteStrings',true)
histogram(rTRE)