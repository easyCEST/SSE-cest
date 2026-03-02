function [MTRasymM0_LiverRange] = MTRcontrast_LiverRange_RHQ(S0B0path,V_as,x_as,S0,Thmask,displayrange,offset_range,liverMask,alpha)
close all;
[bx, by]=find(liverMask);
Nx = size(V_as,1)
Ny = size(V_as,2)
MTRasymM0_LiverRange = zeros(Nx,Ny);
MTR_mRange = zeros(Nx,Ny);
for i = 1:length(offset_range)
    offset = offset_range(i);
    if offset > 0
        [~,index]= min(abs(x_as-offset));
        [~,index_m]= min(abs(x_as-(-offset))); %寻找与-offset最近的index，以便计算aysm
        MTRasymM0 = (V_as(:,:,index_m)-V_as(:,:,index))./(S0+1e-5).*Thmask;
        MTRasymM0_LiverRange = MTRasymM0_LiverRange + MTRasymM0;
    else
        [~,index]= min(abs(x_as-offset));
        MTR_m = (V_as(:,:,index))./(S0+1e-5).*Thmask;
        MTR_mRange = MTR_mRange + MTR_m;
%         strOff2 = ['MTRm',num2str(offset),'ppm.mat'];
%         save(fullfile(S0B0path, strOff2),'MTR_m')
%         savefigRegion(MTR_m,displayrange.MTRpm, bx, by, strOff2, S0B0path, Thmask)
    end
end
        MTRasymM0_LiverRange = MTRasymM0_LiverRange./length(offset_range);
        save(fullfile(S0B0path,'MTRasymM0_LiverRange.mat'),'MTRasymM0_LiverRange')
        savefigRegion(MTRasymM0_LiverRange,displayrange.MTRasymM0, bx, by, S0B0path, Thmask)
        MTR_mRange = MTR_mRange./length(offset_range);
        save(fullfile(S0B0path, 'MTR_mRange.mat'),'MTR_mRange')
%         savefigRegion(MTR_mRange,displayrange.MTRpm, bx, by, S0B0path, Thmask)
end
function savefigRegion(data,dispRange, bx, by, path,mask)
    tmp = medfilt2(data.*mask,[2, 2]);
    tmp = data.*mask;
    tmp(mask==0) = dispRange(1);
    try
        figure; imagesc(tmp(min(bx)-3:max(bx)+3,min(by)-3:max(by)+3),dispRange); 
    catch
        figure; imagesc(tmp(min(bx):max(bx),min(by):max(by)),dispRange); 
    end
    axis off; 
    %colormap(slanCM('gist_rainbow'));
    colorbar; altered_colormap();
    title('MTRasymM0 Range of 0.5-1.5ppm')
    savefig(fullfile(path,'MTRasymM0_LiverRange.fig'))
end