% function ROIplotSong(nROI,legends,Thmask,ROIpath,Img)
% function [ROI_total]=plotROI(nROI,legends,Thmask,ROIpath,V_exp,Img)
function [ROI_total]= ROIplotSong(nROI,legends,Thmask,ROIpath,Img)
close all
% this function plot ROI and store it as a struct
str=[num2str(nROI),'ROI']
    try load(fullfile(ROIpath,'ROI_total'))
    catch
        figure
%        imagesc(Img,[-0.1 0.1]); axis off; colormap('jet'); colorbar
% imshow(Img,[]);
%         imagesc(Img,[-0.1 0.1]); axis off; colormap('jet'); colorbar
        imagesc(Img); axis off; colormap('gray'); caxis([0, 1]); colorbar
%        imagesc(Img); axis off; colormap('gray'); caxis([min(Img(:)), max(Img(:))]); colorbar
       [Nx Ny]=size(Img)
%        imagesc(S0_use); axis off; colormap('jet'); colorbar
       ROI_total = cell(nROI,1);
       for k=1:nROI
           ['Please select ', legends(k), ' now>>>']
           ROI = roipoly;
           save(fullfile(ROIpath, ['ROI',num2str(k)]),'ROI')
           index = find(ROI.*Thmask);
            if isempty(index)
               ['ROI' legends(k) 'too small, please plot again']
               ROI = roipoly;	%matrix of 0's outside the roi and 1's inside
               index = find(ROI.*Thmask);
            end
           ROI_total{k} = index;
           ROIname = ['ROI' num2str(k)]
           show_ROI(Nx, Ny, {index}, fullfile(ROIpath, ROIname));         
       end
%        ROI1=zeros(Nx,Ny);ROI2=zeros(Nx,Ny);ROI3=zeros(Nx,Ny);ROI4=zeros(Nx,Ny);ROI5=zeros(Nx,Ny);
%        ROI1(44:45, 28:29) = 1; ROI2(44:45, 42:43) = 1;ROI3(46:47, 53:54) = 1;ROI4(33:34, 30:31) = 1;ROI5(33:34, 24:25) = 1;
%        ROI_total{1} = find(ROI1.*Thmask);
%        ROI_total{2} = find(ROI2.*Thmask);
%        ROI_total{3} = find(ROI3.*Thmask);
%        ROI_total{4} = find(ROI4.*Thmask);
%        ROI_total{5} = find(ROI5.*Thmask);
    end
    save(fullfile(ROIpath, 'ROI_total'),'ROI_total')
