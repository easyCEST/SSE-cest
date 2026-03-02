function ROIname= ROI_name(nROI)
%%选择ROI的数量
ROIname = [];
for ii = 1:nROI
    ROIname = [ROIname, string(['ROI', num2str(ii)])];
end
