function [liverMask] = GenerateliverMask1107(S0,savepath,Thmask_datapath)
[nx,ny] = size(S0);
try 
    load(fullfile(Thmask_datapath, 'liverMask.mat'))
    liverMask = imresize(liverMask,size(S0));
catch
    figure; imagesc(S0); axis off; 
    'please plot the liver region'
    liverMask = roipoly;
    ROI = find(liverMask);  
    show_ROI(nx,ny,{ROI},fullfile(savepath,'liverMask.fig'))
    save(fullfile(savepath, 'liverMask.mat'), 'liverMask')
end