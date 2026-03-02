function [subtract_Zspc,subtract_Zspc_uncorrect] = MTRcontrast(slice,S0B0path,single_offset,V_as,x_as,V_as_uncorrect,S0,Thmask,displayrange,brainMask,alpha,x_lw)
close all;
[bx, by]=find(brainMask);
[Nx,Ny,Nz]=size(V_as);
Zspc = zeros(Nx,Ny,Nz);
Zspc_uncorrect = zeros(Nx,Ny,Nz);
for i=1:size(V_as,3)
        Zspc(:,:,i) = (V_as(:,:,i))./(S0+1e-5).*Thmask;
        Zspc_uncorrect(:,:,i) = (V_as_uncorrect(:,:,i))./(S0+1e-5).*Thmask;
end
subtract_Zspc = 1-Zspc;
subtract_Zspc_uncorrect = 1-Zspc_uncorrect;
save(fullfile(S0B0path, 'subtract_Zspc'),'subtract_Zspc')
save(fullfile(S0B0path, 'subtract_Zspc_uncorrect'),'subtract_Zspc_uncorrect')

for i = 1:length(single_offset)
    offset = single_offset(i);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%调整ST的ColorBar
    if(abs(offset)==3.5 || abs(offset)==2)
        displayrange_ST = [0.55 0.85];
    else
        displayrange_ST = [0.6 0.9];
    end


    if offset > 0
        [~,index]= min(abs(x_as-offset));
        if(offset==2 || offset==3.5)
            subtract_Zspc_single = mean(subtract_Zspc(:,:,(index-(0.1/x_lw)):(index+(0.1/x_lw))),3);
            subtract_Uncor_Zspc_single = mean(subtract_Zspc_uncorrect(:,:,(index-(0.1/x_lw)):(index+(0.1/x_lw))),3);
            strOff1 = ['subtract_Zspc',num2str(offset),'ppm_Mean.mat'];
            strOff5 = ['subtract_Uncor_Zspc',num2str(offset),'ppm_Mean.mat'];
        else
            subtract_Zspc_single = subtract_Zspc(:,:,index);
            subtract_Uncor_Zspc_single = subtract_Zspc_uncorrect(:,:,index);
            strOff1 = ['subtract_Zspc',num2str(offset),'ppm.mat'];
            strOff5 = ['subtract_Uncor_Zspc',num2str(offset),'ppm.mat'];
        end
        save(fullfile(S0B0path, strOff1),'subtract_Zspc_single')
        save(fullfile(S0B0path, strOff5),'subtract_Uncor_Zspc_single')
        savefigRegion(subtract_Zspc_single,displayrange_ST, bx, by, strOff1, S0B0path, Thmask)
        savefigRegion(subtract_Uncor_Zspc_single,displayrange_ST, bx, by, strOff5, S0B0path, Thmask)
        

%         MTRasym = 1-(V_as(:,:,index)+1e-5)./(V_as(:,:,-index)+1e-5);
%         strOff3 = ['MTRasym',num2str(offset),'ppm.mat'];
%         save(fullfile(S0B0path, strOff3),'MTRasym')
%         savefigRegion(MTRasym,displayrange.MTRasym, bx, by, strOff3, S0B0path, Thmask)
        [~,index_m]= min(abs(x_as-(-offset))); %寻找与-offset最近的index，以便计算aysm
        MTRasymM0 = (V_as(:,:,index_m)-V_as(:,:,index))./(S0+1e-5).*Thmask;
        strOff4 = ['MTRasymM0',num2str(offset),'ppm.mat'];
        save(fullfile(S0B0path, strOff4),'MTRasymM0')
        savefigRegion(MTRasymM0,displayrange.MTRasymM0, bx, by, strOff4, S0B0path, Thmask)

%         MTRasym_uncorrect = 1-(V_as_uncorrect(:,:,in_p)+1e-5)./(V_as(:,:,in_m)+1e-5);
%         strOff7 = ['MTRasymUncorrect',num2str(offset),'ppm.mat'];
%         save(fullfile(S0B0path, strOff7),'MTRasym')
%         savefigRegion(MTRasym_uncorrect,displayrange.MTRasym, bx, by, strOff7, S0B0path, Thmask)
        MTRasymM0_uncorrect = (V_as_uncorrect(:,:,index_m)-V_as_uncorrect(:,:,index))./(S0+1e-5).*Thmask;
        strOff8 = ['MTRasymM0uncorrect',num2str(offset),'ppm.mat'];
        save(fullfile(S0B0path, strOff8),'MTRasymM0_uncorrect')
        savefigRegion(MTRasymM0_uncorrect,displayrange.MTRasymM0, bx, by, strOff8, S0B0path, Thmask)
        
%         alpha = 1;
        MTRasymM0FC = alpha*MTRasymM0.*(V_as(:,:,index_m)+V_as(:,:,index))./(S0+1e-5).*Thmask;
        strOff9 = ['MTRasymM0FC',num2str(offset),'ppm.mat'];
        save(fullfile(S0B0path, strOff9),'MTRasymM0FC')
        savefigRegion(MTRasymM0FC,displayrange.MTRasymM0, bx, by, strOff9, S0B0path, Thmask)
        
        MTRasymM0FC_uncorrect = alpha*MTRasymM0_uncorrect.*(V_as_uncorrect(:,:,index_m)+V_as_uncorrect(:,:,index))./(S0+1e-5).*Thmask;
        strOff10 = ['MTRasymM0FC_uncorrect',num2str(offset),'ppm.mat'];
        save(fullfile(S0B0path, strOff10),'MTRasymM0FC_uncorrect')
        savefigRegion(MTRasymM0FC_uncorrect,displayrange.MTRasymM0, bx, by, strOff10, S0B0path, Thmask)
        
    else

        [~,index]= min(abs(x_as-offset));
        if(offset==-3.5)
            subtract_Zspc_single = mean(subtract_Zspc(:,:,(index-(0.5/x_lw)):(index+(0.5/x_lw))),3);
            subtract_Uncor_Zspc_single = mean(subtract_Zspc_uncorrect(:,:,(index-(0.5/x_lw)):(index+(0.5/x_lw))),3);
            strOff2 = ['subtract_Zspc',num2str(offset),'ppm_Mean.mat'];
            strOff6 = ['subtract_Uncor_Zspc',num2str(offset),'ppm_Mean.mat'];
        else
            subtract_Zspc_single = subtract_Zspc(:,:,index);
            subtract_Uncor_Zspc_single = subtract_Zspc_uncorrect(:,:,index);
            strOff2 = ['subtract_Zspc',num2str(offset),'ppm.mat'];
            strOff6 = ['subtract_Uncor_Zspc',num2str(offset),'ppm.mat'];
        end
        save(fullfile(S0B0path, strOff2),'subtract_Zspc_single')
        savefigRegion(subtract_Zspc_single,displayrange_ST, bx, by, strOff2, S0B0path, Thmask)
        save(fullfile(S0B0path, strOff6),'subtract_Uncor_Zspc_single')
        savefigRegion(subtract_Uncor_Zspc_single,displayrange_ST, bx, by, strOff6, S0B0path, Thmask)
        %%mean(-0.8~-2)
        [~,indexm]= min(abs(x_as-(-2)));
        [~,indexn]= min(abs(x_as-(-0.8)));
        subtract_Zspc_single = mean(subtract_Zspc(:,:,indexm:indexn),3);
        subtract_Uncor_Zspc_single = mean(subtract_Zspc_uncorrect(:,:,indexm:indexn),3);
        strOff2 = ['subtract_Zspc(-2~-0.8)ppm_Mean.mat'];
        strOff6 = ['subtract_Uncor_Zspc(-2~-0.8)ppm_Mean.mat'];
        save(fullfile(S0B0path, strOff2),'subtract_Zspc_single')
        savefigRegion(subtract_Zspc_single,displayrange_ST, bx, by, strOff2, S0B0path, Thmask)
        save(fullfile(S0B0path, strOff6),'subtract_Uncor_Zspc_single')
        savefigRegion(subtract_Uncor_Zspc_single,displayrange_ST, bx, by, strOff6, S0B0path, Thmask)
    end
    
end

originalPos = [3,559,1912,221];
originalUnits = 'pixels';
fig = figure;
fig.Units = originalUnits;
fig.Position = originalPos;
t = tiledlayout(1, length(single_offset)+2,'TileSpacing', 'tight', 'Padding', 'compact');
len = length(single_offset);
for i=1:(len+2)
    if i == (len+1)
        load(fullfile(S0B0path,'MTRasymM03.5ppm.mat'))
        subtract_Zspc_single = MTRasymM0;
        displayrange_ST = [-0.05 0.05];
    elseif i == (len+2)
        load(fullfile(S0B0path,'MTRasymM01.2ppm.mat'))
        subtract_Zspc_single = MTRasymM0;
        displayrange_ST = [-0.05 0.05];
    else
        try
           load(fullfile(S0B0path,['subtract_Zspc',num2str(single_offset(i)),'ppm.mat']));
        catch
           load(fullfile(S0B0path,['subtract_Zspc',num2str(single_offset(i)),'ppm_Mean.mat']));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %调整ST的ColorBar
        if(abs(single_offset(i))==3.5 || abs(single_offset(i))==2)
            displayrange_ST = [0.55 0.85];
        else
            displayrange_ST = [0.6 0.9];
        end
    end


%     subplot(1,length(single_offset),i);
    temp = subtract_Zspc_single.*Thmask;
    temp(Thmask==0) = nan;
    nexttile;
    try
        imagesc(temp(min(bx)-3:max(bx)+3,min(by)-3:max(by)+3),displayrange_ST); 
    catch
        imagesc(temp(min(bx):max(bx),min(by):max(by)),displayrange_ST); 
    end
    axis off; 
    %colormap(slanCM('gist_rainbow'));   % rainbow colorbar 需要导入文件夹
%     colorbar; 
    altered_colormap();
    if i==len+1
        title(['MTRasymM03.5ppm']);
    elseif i==len+2
        title(['MTRasymM01.2ppm']);
    else
        title(['subtract_Zspc',num2str(single_offset(i)),'ppm']);
    end
end
layer_number = slice; % 替换为你需要的层编号
annotation(fig, 'textbox', [0.005, 0.5, 0.05, 0.1], ...
    'String', ['slice ', num2str(layer_number)], ...
    'FitBoxToText', 'on', ...
    'EdgeColor', 'none', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left');
% colorbar; altered_colormap();
saveas(gcf, fullfile(S0B0path,'ST_All.fig'));
saveas(gcf, fullfile(S0B0path,'ST_All.jpg'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%逐像素对Z谱拉伸到[0，1]
Zspc = zeros(Nx,Ny,Nz);
Zspc_uncorrect = zeros(Nx,Ny,Nz);
for i=1:size(V_as,3)
        Zspc(:,:,i) = (V_as(:,:,i))./(S0+1e-5).*Thmask;
        Zspc_uncorrect(:,:,i) = (V_as_uncorrect(:,:,i))./(S0+1e-5).*Thmask;
end
for i=1:size(Zspc,1)
    for j=1:size(Zspc,2)
        spectra = squeeze(Zspc(i,j,:));
        spectra = spectra - min(spectra);
        spectra = spectra./max(spectra);
        Zspc(i,j,:) = spectra;

        spectra = squeeze(Zspc_uncorrect(i,j,:));
        spectra = spectra - min(spectra);
        spectra = spectra./max(spectra);
        Zspc_uncorrect(i,j,:) = spectra;

        spectra = squeeze(V_as(i,j,:));
        spectra = spectra - min(spectra);
        spectra = spectra./max(spectra);
        V_as(i,j,:) = spectra;

        spectra = squeeze(V_as_uncorrect(i,j,:));
        spectra = spectra - min(spectra);
        spectra = spectra./max(spectra);
        V_as_uncorrect(i,j,:) = spectra;
    end
end
subtract_Zspc = 1-Zspc;
subtract_Zspc_uncorrect = 1-Zspc_uncorrect;
save(fullfile(S0B0path, 'subtract_Zspc_0to1'),'subtract_Zspc')
save(fullfile(S0B0path, 'subtract_Zspc_uncorrect_0to1'),'subtract_Zspc_uncorrect')

for i = 1:length(single_offset)
    offset = single_offset(i);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%调整ST的ColorBar
    if(abs(offset)==3.5 || abs(offset)==2)
        displayrange_ST = [0.55 0.85];
    else
        displayrange_ST = [0.6 0.9];
    end

    if offset > 0
        [~,index]= min(abs(x_as-offset));
        if(offset==2 || offset==3.5)
            subtract_Zspc_single = mean(subtract_Zspc(:,:,(index-(0.1/x_lw)):(index+(0.1/x_lw))),3);
            subtract_Uncor_Zspc_single = mean(subtract_Zspc_uncorrect(:,:,(index-(0.1/x_lw)):(index+(0.1/x_lw))),3);
            strOff1 = ['subtract_Zspc',num2str(offset),'ppm_Mean_0to1.mat'];
            strOff5 = ['subtract_Uncor_Zspc',num2str(offset),'ppm_Mean_0to1.mat'];
        else
            subtract_Zspc_single = subtract_Zspc(:,:,index);
            subtract_Uncor_Zspc_single = subtract_Zspc_uncorrect(:,:,index);
            strOff1 = ['subtract_Zspc',num2str(offset),'ppm_0to1.mat'];
            strOff5 = ['subtract_Uncor_Zspc',num2str(offset),'ppm_0to1.mat'];
        end
        save(fullfile(S0B0path, strOff1),'subtract_Zspc_single')
        save(fullfile(S0B0path, strOff5),'subtract_Uncor_Zspc_single')
        savefigRegion(subtract_Zspc_single,displayrange_ST, bx, by, strOff1, S0B0path, Thmask)
        savefigRegion(subtract_Uncor_Zspc_single,displayrange_ST, bx, by, strOff5, S0B0path, Thmask)
        

%         MTRasym = 1-(V_as(:,:,index)+1e-5)./(V_as(:,:,-index)+1e-5);
%         strOff3 = ['MTRasym',num2str(offset),'ppm.mat'];
%         save(fullfile(S0B0path, strOff3),'MTRasym')
%         savefigRegion(MTRasym,displayrange.MTRasym, bx, by, strOff3, S0B0path, Thmask)
        [~,index_m]= min(abs(x_as-(-offset))); %寻找与-offset最近的index，以便计算aysm
        MTRasymM0 = (V_as(:,:,index_m)-V_as(:,:,index))./(S0+1e-5).*Thmask;
        strOff4 = ['MTRasymM0',num2str(offset),'ppm_0to1.mat'];
        save(fullfile(S0B0path, strOff4),'MTRasymM0')
        savefigRegion(MTRasymM0,displayrange.MTRasymM0, bx, by, strOff4, S0B0path, Thmask)

%         MTRasym_uncorrect = 1-(V_as_uncorrect(:,:,in_p)+1e-5)./(V_as(:,:,in_m)+1e-5);
%         strOff7 = ['MTRasymUncorrect',num2str(offset),'ppm.mat'];
%         save(fullfile(S0B0path, strOff7),'MTRasym')
%         savefigRegion(MTRasym_uncorrect,displayrange.MTRasym, bx, by, strOff7, S0B0path, Thmask)
        MTRasymM0_uncorrect = (V_as_uncorrect(:,:,index_m)-V_as_uncorrect(:,:,index))./(S0+1e-5).*Thmask;
        strOff8 = ['MTRasymM0uncorrect',num2str(offset),'ppm_0to1.mat'];
        save(fullfile(S0B0path, strOff8),'MTRasymM0_uncorrect')
        savefigRegion(MTRasymM0_uncorrect,displayrange.MTRasymM0, bx, by, strOff8, S0B0path, Thmask)
        
%         alpha = 1;
        MTRasymM0FC = alpha*MTRasymM0.*(V_as(:,:,index_m)+V_as(:,:,index))./(S0+1e-5).*Thmask;
        strOff9 = ['MTRasymM0FC',num2str(offset),'ppm_0to1.mat'];
        save(fullfile(S0B0path, strOff9),'MTRasymM0FC')
        savefigRegion(MTRasymM0FC,displayrange.MTRasymM0, bx, by, strOff9, S0B0path, Thmask)
        
        MTRasymM0FC_uncorrect = alpha*MTRasymM0_uncorrect.*(V_as_uncorrect(:,:,index_m)+V_as_uncorrect(:,:,index))./(S0+1e-5).*Thmask;
        strOff10 = ['MTRasymM0FC_uncorrect',num2str(offset),'ppm_0to1.mat'];
        save(fullfile(S0B0path, strOff10),'MTRasymM0FC_uncorrect')
        savefigRegion(MTRasymM0FC_uncorrect,displayrange.MTRasymM0, bx, by, strOff10, S0B0path, Thmask)
        
    else

        [~,index]= min(abs(x_as-offset));
        if(offset==-3.5)
            subtract_Zspc_single = mean(subtract_Zspc(:,:,(index-(0.5/x_lw)):(index+(0.5/x_lw))),3);
            subtract_Uncor_Zspc_single = mean(subtract_Zspc_uncorrect(:,:,(index-(0.5/x_lw)):(index+(0.5/x_lw))),3);
            strOff2 = ['subtract_Zspc',num2str(offset),'ppm_Mean_0to1.mat'];
            strOff6 = ['subtract_Uncor_Zspc',num2str(offset),'ppm_Mean_0to1.mat'];
        else
            subtract_Zspc_single = subtract_Zspc(:,:,index);
            subtract_Uncor_Zspc_single = subtract_Zspc_uncorrect(:,:,index);
            strOff2 = ['subtract_Zspc',num2str(offset),'ppm_0to1.mat'];
            strOff6 = ['subtract_Uncor_Zspc',num2str(offset),'ppm_0to1.mat'];
        end
        save(fullfile(S0B0path, strOff2),'subtract_Zspc_single')
        savefigRegion(subtract_Zspc_single,displayrange_ST, bx, by, strOff2, S0B0path, Thmask)
        save(fullfile(S0B0path, strOff6),'subtract_Uncor_Zspc_single')
        savefigRegion(subtract_Uncor_Zspc_single,displayrange_ST, bx, by, strOff6, S0B0path, Thmask)
        %%mean(-0.8~-2)
        [~,indexm]= min(abs(x_as-(-2)));
        [~,indexn]= min(abs(x_as-(-0.8)));
        subtract_Zspc_single = mean(subtract_Zspc(:,:,indexm:indexn),3);
        subtract_Uncor_Zspc_single = mean(subtract_Zspc_uncorrect(:,:,indexm:indexn),3);
        strOff2 = ['subtract_Zspc(-2~-0.8)ppm_Mean_0to1.mat'];
        strOff6 = ['subtract_Uncor_Zspc(-2~-0.8)ppm_Mean_0to1.mat'];
        save(fullfile(S0B0path, strOff2),'subtract_Zspc_single')
        savefigRegion(subtract_Zspc_single,displayrange_ST, bx, by, strOff2, S0B0path, Thmask)
        save(fullfile(S0B0path, strOff6),'subtract_Uncor_Zspc_single')
        savefigRegion(subtract_Uncor_Zspc_single,displayrange_ST, bx, by, strOff6, S0B0path, Thmask)
    end
    
end

originalPos = [3,559,1912,221];
originalUnits = 'pixels';
fig = figure;
fig.Units = originalUnits;
fig.Position = originalPos;
t = tiledlayout(1, length(single_offset),'TileSpacing', 'tight', 'Padding', 'compact');

for i=1:length(single_offset)
    try
       load(fullfile(S0B0path,['subtract_Zspc',num2str(single_offset(i)),'ppm_0to1.mat']));
    catch
       load(fullfile(S0B0path,['subtract_Zspc',num2str(single_offset(i)),'ppm_Mean_0to1.mat']));
    end
    
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%调整ST的ColorBar
    if(abs(single_offset(i))==3.5 || abs(single_offset(i))==2)
        displayrange_ST = [0.55 0.85];
    else
        displayrange_ST = [0.6 0.9];
    end

%     subplot(1,length(single_offset),i);
    temp = subtract_Zspc_single.*Thmask;
    temp(Thmask==0) = nan;
    nexttile;
    try
        imagesc(temp(min(bx)-3:max(bx)+3,min(by)-3:max(by)+3),displayrange_ST); 
    catch
        imagesc(temp(min(bx):max(bx),min(by):max(by)),displayrange_ST); 
    end
    axis off; 
%     colormap(slanCM('gist_rainbow'));   % rainbow colorbar 需要导入文件夹
    colorbar; altered_colormap();
    title(['subtract_Zspc',num2str(single_offset(i)),'ppm_0to1.mat']);
end
%     colorbar; altered_colormap();
    saveas(gcf, fullfile(S0B0path,'ST_0to1_All.fig'));
    saveas(gcf, fullfile(S0B0path,'ST_0to1_All.jpg'));
end





function savefigRegion(data,dispRange, bx, by, NameStr, path,mask)
%     tmp = medfilt2(data.*mask,[2, 2]);
    tmp = data.*mask;
%     tmp(mask==0) = dispRange(1);
    tmp(mask==0) = nan;
    try
        figure; imagesc(tmp(min(bx)-3:max(bx)+3,min(by)-3:max(by)+3),dispRange); 
    catch
        figure; imagesc(tmp(min(bx):max(bx),min(by):max(by)),dispRange); 
    end
    axis off; 
    %colormap(slanCM('gist_rainbow'));   % rainbow colorbar 需要导入文件夹
    colorbar; altered_colormap();
    title(NameStr)
    savefig(fullfile(path,[NameStr 'liver.fig']))
end
%     


