function [mSigE_trimmed,w_offset_trimmed,Lorentzian_fit,LFori,LD,w_offset_1, MTRrex]=...
    LD_Zspec_SXL1105_no_normalize(NegOrPos,N_slice,ROIpath,mSigE,w_offset,nROI,l1ppm,l2ppm,l3ppm,h1ppm, h2ppm, h3ppm,lb,ub,par0)
%%for +Zpositive CEST:LD_Zspec_SXL1027_no_normalize(0,mSigE_Raw,w_offset,4,-0.3,-0.4,-0.6,1.1,6,10);
%%for -Znegtive NOE:LD_Zspec_SXL1027_no_normalize(1,Z(:,2:5),Z(:,1),4,-0.4,-6,-10,0.6,0.7,1.5);
%%if load ROI
% str=[num2str(nROI),'ROI_Raw'];
% load(fullfile(ROIpath,['Zspc_',str,]))
% mSigE = mSigE_Raw;
Nw=length(w_offset);
w_offset_1=(min(w_offset):0.1:max(w_offset))'; 
[~ , ~]=size(w_offset_1);
[~,~]=size(w_offset);

LD = cell(1,nROI); 
Lorentzian_fit = cell(1,nROI);
Z = cell(1,nROI);
ZL = cell(1,nROI);
MTRrex = cell(1,nROI);

for k = 1 : nROI %%%%%%%%%%需要画的ROI的个数
  
    %差分判断，剔除两端递减的数据
    %z_spectrum = 1.37*mSigE(:, k);  % scaling of S0 according to steadystate value
    z_spectrum = mSigE(:, k);  % scaling of S0 according to steadystate value
    offset = w_offset(:);  % 确保是列向量
    assert(length(offset) == length(z_spectrum), 'w_offset 尺寸与 Z 谱不一致');

    if(NegOrPos == 1)
        left_indices = find(offset < -5);
        right_indices = find(offset > 1.5);
    elseif(NegOrPos == 0)
        left_indices = find(offset < -1.5);
        right_indices = find(offset > 5);
    end

    left_violation_count = 0;
    for i = 1:length(left_indices)-1
        cur = left_indices(end - i + 1);     
        prev = left_indices(end - i);        
        if z_spectrum(cur) > z_spectrum(prev)
            left_violation_count = left_violation_count + 1;
        else
            left_violation_count = 0;  
        end
        if left_violation_count >= 3
            z_spectrum(1:(cur+3)) = [];
            offset(1:(cur+3)) = [];
            break;
        end
    end

    if(NegOrPos == 1)
        right_indices = find(offset > 1.5);
    elseif(NegOrPos == 0)
        right_indices = find(offset > 5);
    end

    right_violation_count = 0;
    for i = 1:length(right_indices)-1
        cur = right_indices(i);       
        prev = right_indices(i+1);    
        if z_spectrum(cur) > z_spectrum(prev)
             right_violation_count = right_violation_count + 1;
        else
            right_violation_count = 0;
        end
        if right_violation_count >= 3
            z_spectrum((cur-3):end) = [];
            offset((cur-3):end) = [];
            break;
        end
    end

    mSigE_trimmed{k} = z_spectrum;     
    w_offset_trimmed{k} = offset;

    z_raw = mSigE_trimmed{k};
    w = w_offset_trimmed{k};
    %% 6ppm 至最边缘offset均值作为9.5,10ppm 的Z谱先验
  %% mean(Z(6ppm：max(abs(w))) as MTC
     %%% compensate for MTC，consider as a constant
        [~, max_idx] = max(abs(w));  % 找到w_offset最大值的索引
        MTC6_idx = find (abs(w) == 6);
  %        MTC6_idx = []       %for porcine liver
          if(NegOrPos == 1)  %%neg
              Z_edge = mean(z_raw(max_idx:MTC6_idx)) ;
              % 如果找不到 w = 6，则使用 z_raw 最大值对应的索引
              if isempty(MTC6_idx)
                  Z_edge  = z_raw(max_idx);
              end
          elseif(NegOrPos == 0) %%pos
              Z_edge  = mean(z_raw(MTC6_idx:max_idx)) ;  % 使用6-最远频偏的average作为MTC
              % 如果找不到 w = 6，则使用 z_raw 最大值对应的索引
              if isempty(MTC6_idx)
                  Z_edge  = z_raw(max_idx);
              end
          end
%     z_ds=z_raw+MTC*ones(size(z_raw));
    %% -归一化处理 [0,maxZ]-----------------------------
%    [~, max_idx] = max(abs(w)); 
   % z_norm = z_raw(max_idx)*(z_raw - min(z_raw)) / (z_raw(max_idx)-min(z_raw));
   z_norm = z_raw;
    %% 边缘加入Lorentzian先验约束
    if(NegOrPos == 1)  %%neg
        extra_w = [-9.5, -9.8, -10];
     
    elseif(NegOrPos == 0)
        extra_w = abs([-9.5,-9.8, -10]);
    
    end
   % [~, index] = max(abs(w))
    extra_z = ones(size(extra_w)).*Z_edge ; % [使用最大值填充模拟点
    w_aug = [w(:); extra_w(:)];
    z_aug = [z_norm(:); extra_z(:)];

    [w_sorted, sort_idx] = sort(w_aug);
    z_sorted = z_aug(sort_idx);

    [w_unique, ia] = unique(w_sorted, 'stable');
    z_unique = z_sorted(ia);
    mSigE_processed{k} = z_unique;
    w_offset_processed{k} = w_unique;

    Z{k} = mSigE_processed{k};
    Z{k}(isnan(Z{k})) = 0;
    ZL{k} = spline(w_offset_processed{k}, Z{k}, w_offset_processed{k});

    options = optimset('MaxFunEvals',1000000,'TolFun',1e-10,'TolX',1e-10, 'Display',  'off' );
    options.Algorithm = 'levenberg-marquardt';

    [~,minus10]=min(abs(w_unique-l3ppm));  
    [~,minusP2]=min(abs(w_unique-l2ppm));  
    [~,minusP3]=min(abs(w_unique-l1ppm));  
    [~,plus10]=min(abs(w_unique-h3ppm));  
    [~,plusP2]=min(abs(w_unique-h2ppm));  
    [~,plusP3]=min(abs(w_unique-h1ppm));  

    toFit = ZL{k}([minus10:minusP2 minusP3:plusP3 plusP2:plus10]);
    xdata = w_offset_processed{k}([minus10:minusP2 minusP3:plusP3 plusP2:plus10]);
    ydata = toFit;

    LFori(:,1)= (-11:0.1:11)';
    xall_da = w_offset_processed{k};
    x0 = 0;

    par = lsqcurvefit(@lorentz_N,par0, xdata, ydata, lb, ub, options);
    disp(par)
    Lorentzian_fit{k} = lorentz_N(par, xall_da);
    LFori(:,k+1) = lorentz_N(par, LFori(:,1));
    
    % z=z_unique+(1-z_raw(N7ppm))*ones(size(z_unique));
    %LD{k}= Lorentzian_fit{k} - z_raw - (1-z_raw(N7ppm))*ones(size(z_raw));
    %LD{k}= Lorentzian_fit{k} - z_unique;
    LD{k}= Lorentzian_fit{k} - z_unique(1:length(Lorentzian_fit{k}));
%    MTRrex{k} = 1./(z_raw+1e-5) - 1./(Lorentzian_fit{k}+1e-5);

    % 创建一个不可见的图形窗口
    fig = figure('Visible', 'off');
    str1 = [num2str(k), ' ROIZ谱图'];
    plot(w_offset_processed{k}, LD{k}, 'b--');
    title(str1, 'FontWeight', 'bold', 'FontSize', 14);
    set(gca, 'Xdir', 'reverse', 'FontWeight', 'bold', 'FontSize', 14);
    savefig(fullfile(ROIpath,['LDspectra_ROI', num2str(k), '.fig']));
    close(fig);

    clear right_indices;
    clear left_indices;

    figure
%    plot(w_offset',1.37*mSigE(:,k),'r--')
%    hold on
     plot(w_offset_trimmed{k},z_raw,'k-',w_offset_trimmed{k},z_norm,'b-');
    %plot(w_offset_trimmed{k},z_raw,'k-',w_offset_trimmed{k},z_ds,'b-',w_offset_trimmed{k},z_norm,'g-');
    legend('Raw Z_l_a_b','Normalized Z_l_a_b')
    hold on
    plot(w_offset_processed{k}, [z_unique(1:length(Lorentzian_fit{k})),Lorentzian_fit{k},LD{k}]);
    legend('normalized Z_l_a_b','Lorentzian fitting', 'LD')
    % 根据NegOrPos值确定标题前缀
    if NegOrPos == 1
        prefix = 'neg-slice';
    else
        prefix = 'pos-slice';
    end

    title(sprintf('ROI%d, LD谱 - %s%d', k, prefix, N_slice), ...
        'FontWeight','bold','FontSize',14);

    set(gca, 'Xdir', 'reverse', 'FontWeight','bold','FontSize',14);
    savefig(fullfile(ROIpath, ['Lorentzian_Difference_ROI', num2str(k), '.fig']));
end
%% ====== 补全对齐并写入 Excel ======
Nw_full = length(w_offset);
LD_full = zeros(Nw_full, nROI);
Lorentzian_fit_full = zeros(Nw_full, nROI);

for k = 1:nROI
    w_trim = w_offset_trimmed{k};
    LD_trim = LD{k};
    Lorentzian_trim = Lorentzian_fit{k};

    [tf, idx_in_full] = ismember(round(w_trim,3), round(w_offset,3));
    valid_idx = idx_in_full(tf);

    LD_full(valid_idx, k) = LD_trim(tf);
    Lorentzian_fit_full(valid_idx, k) = Lorentzian_trim(tf);
end
xlswrite(fullfile(ROIpath, 'Lorentzian_Difference_raw.xls'), LD_full, 1, 'A1');
xlswrite(fullfile(ROIpath, 'Lorentzian_fit_raw.xls'), Lorentzian_fit_full, 1, 'A1');


