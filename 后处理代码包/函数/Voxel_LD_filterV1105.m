function [voxel_LD,w_offset_1,Lorentzian_fits,MTRrex_map]=...
    Voxel_LD_filterV1105(NegOrPos,Ns,V_exp_ROI,w_offset,l1ppm,l2ppm,l3ppm,h1ppm, h2ppm, h3ppm,lb,ub,par0)

% 初始化变量
Row=size(V_exp_ROI,1);
Column=size(V_exp_ROI,2);
Nw = size(V_exp_ROI,3);
w_offset_1=[min(w_offset):0.1:max(w_offset)]';

voxel_LD=zeros(Row,Column,length(w_offset_1));
Lorentzian_fits = zeros(Row,Column,length(w_offset_1));
MTRrex_map = zeros(Row,Column,length(w_offset_1));
Z_as=zeros(Row,Column,length(w_offset_1));

% 处理NaN值
V_exp_ROI(isnan(V_exp_ROI)) = 0;

% 主循环处理每个体素
for i=1:Row
    for j=1:Column
        % 获取单个体素的Z谱数据
%         z_spectrum = 1.37 * squeeze(V_exp_ROI(i,j,:));  % scaling of S0 according to steadystate value
        z_spectrum = squeeze(V_exp_ROI(i,j,:)); 
        offset = w_offset(:);  % 确保是列向量
        
        % 检查数据有效性
        if all(z_spectrum == 0) || length(offset) ~= length(z_spectrum)
            continue;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 差分判断，剔除两端递减的数据
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if(NegOrPos == 1)
            left_indices = find(offset < -5);
            right_indices = find(offset > 1.5);
        elseif(NegOrPos == 0)
            left_indices = find(offset < -1.5);
            right_indices = find(offset > 5);
        end

        % 处理左侧（从中心往更负的方向）
        left_violation_count = 0;
        for ii = 1:length(left_indices)-1
            cur = left_indices(end - ii + 1);     
            prev = left_indices(end - ii);        
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

        % 重新计算右侧索引
        if(NegOrPos == 1)
            right_indices = find(offset > 1.5);
        elseif(NegOrPos == 0)
            right_indices = find(offset > 5);
        end

        % 处理右侧（从中心往更正的方向）
        right_violation_count = 0;
        for ii = 1:length(right_indices)-1
            cur = right_indices(ii);       
            prev = right_indices(ii+1);    
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

        % 保存修剪后的数据
        mSigE_trimmed = z_spectrum;     
        w_offset_trimmed = offset;

        z_raw = mSigE_trimmed;
        w = w_offset_trimmed;

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % 数据增强：添加额外的偏移点
        % 6ppm 至最边缘offset均值作为9.5,10ppm 的Z谱先验
      % mean(Z(6ppm：max(abs(w))) as MTC
      %%% compensate for MTC，consider as a constant
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
       
        mSigE_processed = z_unique;
        w_offset_processed = w_unique;
    
        Z = mSigE_processed;
        Z(isnan(Z)) = 0;
        ZL = spline(w_offset_processed, Z, w_offset_processed);





        
        options = optimset('MaxFunEvals',1000000,'TolFun',1e-10,'TolX',1e-10, 'Display',  'off');
        options.Algorithm = 'levenberg-marquardt';

        % 找到拟合区间索引
        [~,minus10]=min(abs(w_unique-l3ppm));  
        [~,minusP2]=min(abs(w_unique-l2ppm));  
        [~,minusP3]=min(abs(w_unique-l1ppm));  
        [~,plus10]=min(abs(w_unique-h3ppm));  
        [~,plusP2]=min(abs(w_unique-h2ppm));  
        [~,plusP3]=min(abs(w_unique-h1ppm));  

        toFit = ZL([minus10:minusP2 minusP3:plusP3 plusP2:plus10]);
        xdata = w_offset_processed([minus10:minusP2 minusP3:plusP3 plusP2:plus10]);
        ydata = toFit;
        xall_da = w_offset_processed;
  
        % 执行拟合 - 直接调用外部的 lorentz_N 函数
        par = lsqcurvefit(@lorentz_N, par0, xdata, ydata, lb, ub, options);
        
        % 生成拟合结果 - 直接调用外部的 lorentz_N 函数
        xall_da = w_offset_trimmed;
        Lorentzian_fit = lorentz_N(par, xall_da);
        LD_result = Lorentzian_fit - z_unique(1:(length(Lorentzian_fit)));
        MTRrex_result = 1./(z_raw+1e-5) - 1./(Lorentzian_fit+1e-5);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 对齐到统一频偏网格并存储结果
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        w_offset_trimmed_rounded = round(w_offset_trimmed, 4);
        w_offset_1_rounded = round(w_offset_1, 4);
        [tf, idx_in_w_offset_1] = ismember(w_offset_trimmed_rounded, w_offset_1_rounded);
        
        % 存储结果
        valid_indices = idx_in_w_offset_1(tf);
        voxel_LD(i,j,valid_indices) = LD_result(tf)';
        Lorentzian_fits(i,j,valid_indices) = Lorentzian_fit(tf)';
        MTRrex_map(i,j,valid_indices) = MTRrex_result(tf)';

    end
end

end

% 注意：这里不再包含 lorentz_N 函数的定义
% MATLAB会自动调用同目录下的 lorentz_N.m 文件
