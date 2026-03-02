function show_ROI(nx, ny, ROI_total,savename)
% p: handle with nROI
% nx, ny: matrix size
% ROI_total: stored ROI matrix
% savename: the path & name to save the ROI image
%

nROI = length(ROI_total);
colorRange = {'b-','-r','k-','k--','g','b','w','k'};
for k = 1 : nROI
    hold all
    temp = zeros(nx,ny);
    ind = ROI_total{k};
    temp(ind) = 1; 
%     contour(gca, temp, [0.5 0.5],colorRange(k),'LineWidth',2);
    contour(gca, temp, [0.5 0.5],'b-','LineWidth',3);
    clear ind temp
%     pause
end

savefig(savename)

