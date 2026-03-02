function   [x_as,V_as,V_as_uncorrect,w_offset,V_exp]=Zspec_InterpandB0Corr_RHQ(NegOrPos,x_lw,w1,Vo_ori_use,Thmask,B0map,B0correct,Ns)
    %% 取有效值
    if NegOrPos == 2
        offsetnum = 121;
    elseif NegOrPos == 3
        offsetnum = 61;
    elseif NegOrPos == 0 || NegOrPos == 1
        offsetnum = 121;

    %add by liuhao,用于处理phantom3T数据
    elseif NegOrPos == 4
%         offsetnum = 241;
        offsetnum = 201;
%         offsetnum = 191;
    end

    if length(w1) == 201
        i = 1;
        for temp = 1:length(w1)
            if mean(Vo_ori_use(:,:,temp),[1,2]) ~= 0 
                Vo_tmp(:,:,i) = Vo_ori_use(:,:,temp);
                w_tmp(1,i) = temp;
                i = i+1;
            end
        end
        i = 1;
        for temp = 1:length(w1)
            if i > offsetnum
                break
            end
            if temp == w_tmp(1,i)
                w_offsetTemp(1,i) = w1(1,temp);
                i = i+1;
            end
        end
    end
    if length(w1) == 241
        i = 1;
%         for temp = 1:length(w1)
% %                 if mean(Vo_ori_use(:,:,temp),[1,2]) ~= 0 && temp>=26 && temp<=216
%                 if mean(Vo_ori_use(:,:,temp),[1,2]) ~= 0
%                     Vo_tmp(:,:,i) = Vo_ori_use(:,:,temp);
%                     w_tmp(1,i) = temp;
%                     i = i+1;
%                 end
%          end

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %phantom浓度梯度，频偏没有规律
        if NegOrPos==4
            for temp = 1:length(w1)
                if mean(Vo_ori_use(:,:,temp),[1,2]) ~= 0 && temp>=21 && temp<=221
%                 if mean(Vo_ori_use(:,:,temp),[1,2]) ~= 0
                    Vo_tmp(:,:,i) = Vo_ori_use(:,:,temp);
                    w_tmp(1,i) = temp;
                    i = i+1;
                end
            end
        else
            for temp = 1:length(w1)
                if mean(Vo_ori_use(:,:,temp),[1,2]) ~= 0
                    Vo_tmp(:,:,i) = Vo_ori_use(:,:,temp);
                    w_tmp(1,i) = temp;
                    i = i+1;
                end
            end
        end
        
        i = 1;
        for temp = 1:length(w1)
            if i > offsetnum 
                break
            end
            if temp == w_tmp(1,i)
                w_offsetTemp(1,i) = w1(1,temp);
                i = i+1;
            end
        end
    end
    w_offsetOri = w_offsetTemp;
    [Ny,Nx,Nz1]=size(Vo_ori_use);
    Nw = length(w_offsetOri); 
    [w_offset, I] = sort(w_offsetOri);
%     for in = 1 : Nw
%         Vo_cest(:,:,in) = medfilt2(Vo_tmp(:,:,in),[3 3]);
%     end
    Vo_cest = Vo_tmp; %去掉平滑
    [px,py]=find(Thmask);
    [Ny,Nx,Nz1]=size(Vo_cest);
    %% 插值和B0校正  
    %%%%%V_as:the shifted Z-spectra//_as:插值后的 _exp插值前的
    if length(w1) == 201
        if NegOrPos == 2
            x_as = [w1(20+Ns):x_lw:w1(140+Ns)];%APT
        elseif NegOrPos == 3
            x_as = [w1(50+Ns):x_lw:w1(110+Ns)];%Glyco
        end
    else 
        if NegOrPos == 1
            x_as = [w1(Ns+10):x_lw:w1(130+Ns)];% NEG
        elseif NegOrPos == 0
            x_as = [w1(70+Ns):x_lw:w1(190+Ns)];% POS

        %add by liuhao,用于处理phantom3T数据
        elseif NegOrPos == 4
%             x_as = [w1(26):x_lw:w1(216)];
%             x_as = [w1(10+Ns):x_lw:w1(210+Ns)];
%             x_as = [w1(16):x_lw:w1(224)];
%             x_as = [min(w1):x_lw:max(w1)];
            x_as = [w1(21):x_lw:w1(221)];

        end
    end
    V_as=zeros(Nx,Ny,length(x_as));
    V_exp=zeros(Nx,Ny,length(w_offsetOri));
    V_as_uncorrect=zeros(Nx,Ny,length(x_as));
    V_exp_uncorrect=zeros(Nx,Ny,length(w_offsetOri));

    if B0correct ~= 0
    end
    for in=1:length(px)
        spectra=squeeze(Vo_cest(px(in),py(in),:));
        if B0correct == 1
        %%%% load existed B0 map，using pchip B0 correction
            sft=B0map(px(in),py(in));
            V_as(px(in),py(in),:)= pchip(w_offset-sft,spectra,x_as); 
            V_exp(px(in),py(in),:)= pchip(w_offset-sft,spectra,w_offset);
        elseif B0correct == 0
        %%%% NO B0 correction, only interpolation of Zspectra
            if length(w_offset) > 2
                V_as(px(in),py(in),:)= pchip(w_offset,spectra,x_as); 
                V_exp(px(in),py(in),:)= pchip(w_offset,spectra,w_offset); 
            else 
                V_as(px(in),py(in),:) = nonzeros(spectra);
                V_exp(px(in),py(in),:) = nonzeros(spectra);
            end   
        elseif B0correct == 2
       %%%% with Self B0 correction
       %%%%1.find the lowest point using spline interpolation of original Z-spectra
            yi = spline(w_offset,spectra,x_as);
            %%%calculate spectrum-shift according to Z-spectra itself
            if (NegOrPos == 2)
                [value_min,ind_min] = min(yi);
                sft=x_as(ind_min); %%shifted frequency  
            elseif NegOrPos == 3
                [value_min,ind_min] = min(yi(45:165));
                sft=x_as(ind_min+44); %%shifted frequency 
            elseif (NegOrPos == 0 || NegOrPos == 1)
                % 设定搜索频偏范围，例如 [-3, 3] ppm
                search_min = -1.5;
                search_max = 1.5;
                
                % 找到落在这个范围内的索引
                valid_idx = find(x_as >= search_min & x_as <= search_max);
                
                % 从限制范围内找最小值点作为主峰
                [~, local_min_idx] = min(yi(valid_idx));
                sft = x_as(valid_idx(local_min_idx));  % 获取该点对应的频偏作为 shift
            else 
                %%%%%%%%% add by liuhao,用于处理phantom3T数据
                [value_min,ind_min] = min(yi);
                sft=x_as(ind_min); %%shifted frequency 
            end
            V_as(px(in),py(in),:)= pchip(w_offset-sft,spectra,x_as); 
            V_exp(px(in),py(in),:)= pchip(w_offset-sft,spectra,w_offset); 
        end
        V_as_uncorrect(px(in),py(in),:)= pchip(w_offset,spectra,x_as); 
        V_exp_uncorrect(px(in),py(in),:)= pchip(w_offset,spectra,w_offset);
    end