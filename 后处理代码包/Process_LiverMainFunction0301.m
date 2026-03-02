%% 针对人体肝脏CEST磁共振图像的后处理代码
%% By Ren Haiqi   --2024.6.25
clc
clear all
close all
disp('************************* Notice **************************');
NegOrPos = input('Please enter the sequence type that needs to be processed:\n(0 for POS/1 for NEG/2 for APT/3 for Glyco/4 for phantom3T)');
mode = input('Please enter the sequence type that needs to be processed:\n(0 for no lb/1 for with lb/)');
%%%%%%%%%%%%%% root folder with Scanned data and also put analysis results
datapath='C:\Users\likaixiang\Desktop\2026\03\后处理代码包\数据\fasting';
resultpath='C:\Users\likaixiang\Desktop\2026\03\后处理代码包\结果\fasting';
DataFormat = 1;   % 0:bruker;1:dicom
register = 0;   % register ：  0:no register;  1:rigid ;  2:RPCA
B0correct = 2;   % 0: no B0 correction; 1: WASSR correction(with B0map); 2: self correction 3:2-Pools LD Corr
%%%%%%%%%%%%% the slice selected for analysis，如果是0则分析所有的slice
N_slice=0;   % 选择要分析的slice，如果是0则分析所有的slice
total_slice_num = 41;   % 总层数
multi_slice = 1;   % data dimension: 0: single-slice; 1: multi-slice
snr=15;
nROI=2;   % ROI数量
use_exist_ROI = 1;   % 是否用已有ROI路径 1: Yes; 0: No
Thmask_datapath ='C:\Users\likaixiang\Desktop\2026\03\后处理代码包\结果\fasting\no_lb_neg_interp_0227\slice21';
ROI_datapath=Thmask_datapath;
% 可根据需要修改勾画ROI的背景图

if NegOrPos == 0
    background_frequency = 155;
else 
    background_frequency = 85;
end

date_suffix = '0227';  % 根据日期和其他因素修改这里即可

if NegOrPos == 0 && mode == 0
    filename = ['no_lb_pos_interp_' date_suffix];
elseif NegOrPos == 1 && mode == 0
    filename = ['no_lb_neg_interp_' date_suffix];
elseif NegOrPos == 0 && mode == 1
    filename = ['with_lb_pos_interp_' date_suffix];
else  
    filename = ['with_lb_neg_interp_' date_suffix];
end

SavepathRoot=fullfile(resultpath,filename);
if ~exist (SavepathRoot,'dir')
    mkdir(SavepathRoot)
end
if ~multi_slice
    slice_ind_list = 1;
elseif N_slice == 0
    slice_ind_list = 1:total_slice_num;
    tempIndex = min(slice_ind_list)-1;
elseif N_slice <= total_slice_num
    slice_ind_list = N_slice;
    tempIndex = N_slice-1; % for select ROI path
else
    error(['Only N_slice <= ',num2str(total_slice_num),' is valid.'])
end
b0corr_type_strings = ["No_corr","B0map_corr","Self_corr"];
if use_exist_ROI
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %每层使用不同的ROI
    % ROIpath_list = [];
    % for Ns = slice_ind_list
    %     Savepath = fullfile(SavepathRoot,['slice',num2str(Ns)]);
    %     ROIpath = uigetdir(ROI_datapath, ['Choose the ROIpath for slice ',num2str(Ns)]);
    %     ROIpath_list = [ROIpath_list, string(ROIpath)];
    % end

    %每层使用相同的ROI
        ROIpath_list = [];
        for Ns = slice_ind_list
            Savepath = fullfile(SavepathRoot,['slice',num2str(Ns)]);
        end
        ROIpath = uigetdir(ROI_datapath, ['Choose the ROIpath for slice ',num2str(Ns)]);
        ROIpath_list = [ROIpath_list, string(ROIpath)];
end
% offset choice MTRasym定量时使用
offset_range = 0.5:0.1:1.5;
if NegOrPos == 2 
    Single_off = [3.5 2 1.2 -1.2 -3.5];
%     Single_off = [3.5 -3.5 -1.2 1.2 2 -0.8 0.8];
elseif NegOrPos == 3 
    Single_off=[-1 1];

% add by liuhao,用于处理phantom3T数据
elseif NegOrPos == 4
    Single_off = [3.5 -3.5];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif NegOrPos == 1 
    Single_off = [-3.5 -0.8 0.8 -1.2 1.2 2 3.5];
elseif NegOrPos == 0 
    Single_off = [-0.8 0.8 -1.2 1.2 2 3.5];
end

orientation=0;    % readout direction
% optimized display range 设置定量图Colorbar大小
displayrange.MTRpm = [0 1];
displayrange.MTRpmS0 = [0.5 1];
displayrange.MTRasym = [-0.15, 0.05];
displayrange.MTRasymM0 = [-0.05, 0.05];
displayrange.b0map = [-0.5, 0.5];
FS=128;

%%%%%%%%%%%%%% load Raw-data：img_data
if NegOrPos == 0
    load(fullfile(datapath,'0p7P_interp.mat'));
else
    load(fullfile(datapath,'0p7N_interp.mat'));
end
img_interp = 1.37*img_interp;  % compensate for S0 readout timepoint
[Nx,Ny,Nz] = size(img_interp);   % 四维数组：n_slice*n_offset
%%确保w_offsetOri 的矩阵维度为1*x
if size(img_interp,4) == 201
   w_original = -10:0.1:10; 
else
    w_original = -12:0.1:12; 
end
w = w_original;
w_offsetOri = w;
Mz_ori = zeros(Nx,Ny,total_slice_num,length(w_offsetOri));   

for in=1:total_slice_num
    Mz_ori(:,:,in,:) = img_interp(:,:,in,:);
end
%%设置S0
S0_temp = ones(80,80);   % 重建时已做归一化处理，后续均使用全1矩阵S0_temp
M0_stack = zeros(Nx,Ny,total_slice_num);
for in=1:total_slice_num
    M0_stack(:,:,in) = S0(:,:,in);
end
%% 逐层后处理
for Ns = slice_ind_list
    % if Ns == 1 || Ns == 2 || Ns == 40 || Ns == 41
    %     continue
    % end

    Savepath = fullfile(SavepathRoot,['slice',num2str(Ns)]);    
    S0B0path=fullfile(Savepath,b0corr_type_strings{(B0correct+1)});
    if exist (S0B0path)~=7
        mkdir(S0B0path);
    end    
    [Nx,Ny,Nz] = size(Mz_ori);
    clear img_interp
    % load B0map(2D) or B0map stack(3D) 单次屏气序列没有B0map
    if B0correct ~= 1
        B0map = zeros(size(Mz_ori,1),size(Mz_ori,2));
    else
        B0map = b0mapWASSR(:,:,Ns);   % Ns第几层
    end 
    % 后处理部分
    if use_exist_ROI==0
        Thmask_datapath=Savepath;
    end
        
    x_lw=0.025;   % 插值步长
    Zdata = squeeze(Mz_ori(:,:,Ns,:));   % 即未归一化的除S0外的data
    S0=squeeze(M0_stack(:,:,Ns));   % S0
    S0m=squeeze(M0_stack(:,:,Ns));    
    [Zdata] = registration(Zdata,S0,register,S0B0path);   % 运动校正
    [Thmask] = GenerateThmask1107(Savepath,Thmask_datapath,S0,snr);   % 画小三角滤掉噪点
    liverMask = GenerateliverMask1107(S0,Savepath,Thmask_datapath);
    Thmask = liverMask.*Thmask;   % 去掉所画脑区mask的噪点
    save(fullfile(Savepath, 'Thmask'), 'Thmask'); 
    % 插值|B0校正 获得经过B0校正后未插值和插值的Z谱数据
    [x_as,V_as,V_as_uncorrect,w_offset,V_exp]=Zspec_InterpandB0Corr_RHQ(NegOrPos,x_lw,w_offsetOri,Zdata,Thmask,B0map,B0correct,Ns);

    save(fullfile(Savepath, 'w_offset'), 'w_offset');
    save(fullfile(Savepath, 'V_as'), 'V_as');
    save(fullfile(Savepath, 'V_exp'), 'V_exp');
    save(fullfile(Savepath, 'S0'), 'S0');
    save(fullfile(Savepath, 'S0m'), 'S0m');
    save(fullfile(Savepath, 'S0_temp'), 'S0_temp');
    [V_exp_mask] = NormZ(S0_temp,V_exp,liverMask);   % 重建时已做归一化操作，固使用全1矩阵的S0_temp
    [V_as_mask] = NormZ(S0_temp,V_as,liverMask);   % 同上
    save(fullfile(Savepath, 'w_offset'), 'w_offset');
    save(fullfile(Savepath, 'V_exp_mask'), 'V_exp_mask');
    save(fullfile(Savepath, 'V_as_mask'), 'V_as_mask');
    save(fullfile(Savepath, 'x_as'), 'x_as');
%% 量化部分
if NegOrPos ~= 1 && NegOrPos ~= 0
    %%MTRasym量化
    alpha = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MTRcontrast(Ns,S0B0path,Single_off,V_as,x_as,V_as_uncorrect,S0_temp,Thmask,displayrange,liverMask, alpha,x_lw);
    [MTRasymM0_LiverRange] = MTRcontrast_LiverRange_RHQ(S0B0path,V_as,x_as,S0_temp,Thmask,displayrange,offset_range,liverMask, alpha);
    %%新三频偏法量化结果
    [NOE_3pt_new] = ThreeOffset_NEW_RHQ(V_as,x_as,x_lw,-1.75,-0.8,0.8,S0_temp,Thmask);
    display_image_save(Savepath,Thmask,NOE_3pt_new,'Glycogen NEW3pT');
    save(fullfile(Savepath, 'Glycogen_NEW3pT'), 'NOE_3pt_new');
    [Glyco_3pt_new] = ThreeOffset_NEW_RHQ(V_as,x_as,x_lw,2,0.8,-0.8,S0_temp,Thmask);
    display_image_save(Savepath,Thmask,Glyco_3pt_new,'Glucose NEW3pT');
    save(fullfile(Savepath, 'Glucose_NEW3pT'), 'Glyco_3pt_new');
end
%% ROI 分析
    ROIname=ROI_name(nROI);
    if ~use_exist_ROI
        ROIpath=fullfile(S0B0path,['ROI ',strrep(datestr(datetime('now')),':','-')]);
        mkdir(ROIpath) ;
        ROI_total = ROIplotSong(nROI,ROIname,Thmask,ROIpath,squeeze(Zdata(:,:,background_frequency)));
        ROIprocessSXL1(nROI,ROIpath,w_offset,x_lw,Savepath,ROI_total);
        clear img
    else
        ROIpath = fullfile(S0B0path,['ROI ',strrep(datestr(datetime('now')),':','-')]);
        mkdir(ROIpath)
        sourcePath = ROIpath_list{1};  
        targetPath = ROIpath;  
        fileList = dir(fullfile(sourcePath, '*'));  % 获取所有文件
        validFiles = regexp({fileList.name}, '^(ROI|FocusROI)');  % 匹配以 ROI 或 FocusROI 开头的文件
        
        % 遍历并复制符合条件的文件
        for k = 1:length(fileList)
            if ~isempty(validFiles{k}) && ~fileList(k).isdir
                copyfile(fullfile(sourcePath, fileList(k).name), targetPath);
            end
        end
        ROIprocessSXL1_readROI(nROI,ROIpath,w_offset,x_lw,Savepath);         
    end
%% 获得z谱|洛伦兹拟合曲线|LD谱线
    if NegOrPos == 1 || NegOrPos == 0
        alpha = 1;
        % MTRcontrast_0p7N(S0B0path,Single_off,V_as,x_as,V_as_uncorrect,S0_temp,Thmask,displayrange,liverMask, alpha,x_lw);
        % [MTRasymM0_LiverRange] = MTRcontrast_LiverRange_RHQ(S0B0path,V_as,x_as,S0_temp,Thmask,displayrange,offset_range,liverMask, alpha);
        w_offset_LD = -12:0.1:12;
        V_exp_mask_LD = zeros(size(V_exp_mask,1),size(V_exp_mask,2),length(w_offset_LD));
        [~,midneg]=min(abs(w_offset_LD+6));
        [~,midpos]=min(abs(w_offset_LD-6));
        if NegOrPos == 1
            for ii = 1:size(V_exp_mask,1)
                for jj = 1:size(V_exp_mask,2)
                    V_exp_mask_LD(ii,jj,Ns+10:130+Ns) = V_exp_mask(ii,jj,:);
                    temp = 121-Ns;
                    V_exp_mask_LD(ii,jj,midpos:111+temp) = flip(V_exp_mask_LD(ii,jj,Ns+10:midneg));
                end
            end
        else
            for ii = 1:size(V_exp_mask,1)
                for jj = 1:size(V_exp_mask,2)
                    V_exp_mask_LD(ii,jj,70+Ns:190+Ns) = V_exp_mask(ii,jj,:);
                    temp = 200+Ns-121;
                    V_exp_mask_LD(ii,jj,131-temp:midneg) = flip(V_exp_mask_LD(ii,jj,midpos:190+Ns));
                end
            end
        end

        load(fullfile(ROIpath,['Zspc_',num2str(nROI),'ROI_Raw.mat']))

        %% 洛伦兹拟合参数设置
        if mode==0
            lb = [];
        else
            lb = [0, -Inf,  0,   1];
        end
        ub = [0,  Inf,  1.8, Inf];
        par0 = [0,  -50, 20, 1];  % [x0, a, b, c]

        if NegOrPos == 1  % 负扫描参数
            l1ppm = -0.4;  % 左侧第一个峰位置(ppm)
            l2ppm = -6;     % 左侧第二个峰位置(ppm)
            l3ppm = -10;    % 左侧第三个峰位置(ppm)
            h1ppm = 0.6;    % 右侧第一个峰位置(ppm)
            h2ppm = 0.7;    % 右侧第二个峰位置(ppm)
            h3ppm = 1.5;    % 右侧第三个峰位置(ppm)
        else  % 正扫描参数
            l1ppm = -0.3;   % 左侧第一个峰位置(ppm)
            l2ppm = -0.4;   % 左侧第二个峰位置(ppm)
            l3ppm = -0.6;   % 左侧第三个峰位置(ppm)
            h1ppm = 1.1;    % 右侧第一个峰位置(ppm)
            h2ppm = 6;      % 右侧第二个峰位置(ppm)
            h3ppm = 10;     % 右侧第三个峰位置(ppm)
        end


        [Zexp_trimmed,w_offset_trimmed,Lorentzian_fit,LFori,LD,w_offset_1, LD_MTRrex] = ...
            LD_Zspec_SXL1105_no_normalize(NegOrPos,Ns,ROIpath,mSigE_Raw,w_offset,nROI,...
            l1ppm,l2ppm,l3ppm,h1ppm,h2ppm,h3ppm, lb, ub, par0);

        save(fullfile(ROIpath, 'LFori'), 'LFori');
    %% For Voxel-by-voxel LD fitting /逐体素LD
    [voxel_LD,w_offset_1,Lorentzian_fits,MTRrex_map] = Voxel_LD_filterV1105(...
        NegOrPos,Ns,V_exp_mask,w_offset,...
        l1ppm,l2ppm,l3ppm,h1ppm,h2ppm,h3ppm, lb, ub, par0);


                         save(fullfile(Savepath, 'Lorentzian_fits.mat'), 'Lorentzian_fits');
                         save(fullfile(Savepath, 'MTRrex_map.mat'), 'MTRrex_map');
                         save(fullfile(Savepath, 'voxel_LD'), 'voxel_LD');
                         save(fullfile(Savepath, 'w_offsetLD'), 'w_offset_1');

             if NegOrPos == 1
                 n = get_closest_index(w_offset_1 +1.2, voxel_LD, 3);%这块是与-1.2做差，实际上不是1.2
                 NOEmean = mean(voxel_LD(:,:,(n-3):(n+3)),3);
                 display_image_save1103(Savepath, Thmask, NOEmean, sprintf('slice%d LD NOE-1.2',Ns), [0 0.2]);
                 save(fullfile(Savepath, sprintf('slice%d_NOEmean',Ns)), 'NOEmean');

                 n = get_closest_index(w_offset_1 +3.5, voxel_LD, 3);
                 mean_neg3p5 = mean(voxel_LD(:,:,(n-3):(n+3)),3);
                 display_image_save1103(Savepath, Thmask, mean_neg3p5, sprintf('slice%d LD NOE-3.5', Ns), [0 0.15]);
                 save(fullfile(Savepath, sprintf('slice%d_mean_neg3p5',Ns)), 'mean_neg3p5');
             else
                 n = get_closest_index(w_offset_1 -3.5, voxel_LD, 2);
                 APTmean = mean(voxel_LD(:,:,(n-1):(n+2)), 3);%n+2的话上面最后一个参数填2
                 display_image_save1103(Savepath, Thmask, APTmean, sprintf('slice%d LD 3.5ppm', Ns), [0 0.15]);
                 save(fullfile(Savepath, sprintf('slice%d_APTmean', Ns)), 'APTmean');
                 
                 n = get_closest_index(w_offset_1 -1.2, voxel_LD, 2);
                 mean_1p2 = mean(voxel_LD(:,:,(n-1):(n+2)), 3);
                 display_image_save1103(Savepath, Thmask, mean_1p2, sprintf('slice%d LD 1.2ppm', Ns), [0 0.15]);
                 save(fullfile(Savepath, sprintf('slice%d_mean_1p2', Ns)), 'mean_1p2');
                 
                 n = get_closest_index(w_offset_1 -2, voxel_LD, 3);
                 mean_2 = mean(voxel_LD(:,:,(n-1):(n+2)), 3);
                 display_image_save1103(Savepath, Thmask, mean_2, sprintf('slice%d LD 2ppm', Ns), [0 0.15]);
                 save(fullfile(Savepath, sprintf('slice%d_mean_2', Ns)), 'mean_2');
                 close all
             end
    else
    end
end

