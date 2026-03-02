function ROIprocessSXL1(nROI,ROIpath,w_offset,x_lw,Savepath,ROI_total)
% load(fullfile(S0B0path, 'ROI_total'));
load(fullfile(Savepath,['V_as.mat']));
load(fullfile(Savepath,['V_exp.mat']));
% load(fullfile(Savepath,['V_as_mask.mat']));
% load(fullfile(Savepath,['V_exp_mask.mat']));
% V_as = V_as_mask;
% V_exp = V_exp_mask;
N_echo=1;%%暂取echo=1.
x_as=[min(w_offset):x_lw:max(w_offset)];
[Nx,Ny,Nas]= size(V_as);
Nw = length(w_offset);
str=[num2str(nROI),'ROI'];
%%MTRasym no S0
% % 记得修改S0 这里是全1矩阵
%  try load(fullfile(Savepath,['S0.mat']))
%         S0_use = S0;
%     catch
%         load(fullfile(Savepath,['S0m.mat']))
%         S0_use = S0m;
%  end
  try load(fullfile(Savepath,['S0_temp.mat']))
        S0_use = S0_temp;
    catch
        load(fullfile(Savepath,['S0_temp.mat']))
        S0_use = S0_temp;
 end
 for k=1:nROI
     for n_offset=1:Nas
         temp = squeeze(V_as(:,:,n_offset));
         mSig(n_offset,k) = mean2(temp(ROI_total{k}));
         eSig(n_offset,k) = std(temp(ROI_total{k}));
     end
     MTRasym_as(:,k)= 1- mSig(Nas:-1:((Nas+3)/2),k)./mSig(1:(Nas-1)/2,k);
     for n_offset=1:Nw   
         temp = squeeze(V_exp(:,:,n_offset));
         mSigE(n_offset,k) = mean2(temp(ROI_total{k}));
         eSigE(n_offset,k) = std(temp(ROI_total{k}));
     end
     if min(w_offset) < 0
         MTRasym_exp(:,k)=1-mSigE(Nw:-1:((Nw+3)/2),k)./mSigE(1:(Nw-1)/2,k);
     end
         MTR_exp(:,k) = mSigE(Nw:-1:1,k); 
 end
 temp = flip(w_offset);
if min(w_offset) < 0
      MTRe = temp(1:length(w_offset)/2);
else
      MTRe = temp;
end
 if min(w_offset) < 0
      MTRw=[max(x_as):-x_lw:x_lw];
else
      MTRw=[max(x_as):-x_lw:0];
end
%%MTRasymM0
for k=1:nROI
    for n_offset=1:Nas
        temp = squeeze(V_as(:,:,n_offset)./S0_use);
        mSig(n_offset,k) = mean2(temp(ROI_total{k}));%存每个ppm的平均z谱信号
        eSig(n_offset,k) = std(temp(ROI_total{k}));
    end
    if min(w_offset) < 0
         MTRasymM0_as(:,k)= mSig(1:(Nas-1)/2,k)- mSig(Nas:-1:((Nas+3)/2),k);%求MTRasym
         m_offset = w_offset(w_offset<0);
    else
         MTRasymM0_as(:,k)= mSig(Nas:-1:1,k);
    end
%       temp_as = flip([-3:0.025:3]);
      temp_as = flip(x_as);
       if min(x_as) < 0
            MTRe_as= temp_as(1:length(x_as)/2);
        else
            MTRe_as = temp_as;
        end
    for n_offset=1:Nw
        temp = squeeze(V_exp(:,:,n_offset))./S0_use;%%?
        mSigE(n_offset,k) = mean2(temp(ROI_total{k}));
        eSigE(n_offset,k) = std(temp(ROI_total{k}));
    end
  w_offset=sort(w_offset);
  positive=sum(w_offset>0);
  negative=sum(w_offset<0);
  index_po=find(w_offset>0);
  po_w_offset=w_offset(index_po);
  index_ne=find(w_offset<0);
  ne_w_offset=w_offset(index_ne);
  if (positive>negative)%正的ppm大于负的ppm个数时执行
       ne_w_offset=-fliplr(po_w_offset);
       try
            w_offset_1 = [ne_w_offset 0 po_w_offset];%让负的ppm个数变多
       catch
            w_offset_1 = [ne_w_offset' 0 po_w_offset'];
       end
       mSigE_1(:,k) = spline(w_offset,mSigE(:,k),w_offset_1)';   
       [Row_mSigE_1 column_mSigE_1]=size(w_offset_1); 
       MTRasymM0_exp(:,k)=mSigE_1(1:(column_mSigE_1-1)/2,k)-mSigE_1(column_mSigE_1:-1:((column_mSigE_1+3)/2),k); 
    elseif(negative>positive)
       po_w_offset=-fliplr(ne_w_offset);
       try
            w_offset_1 = [ne_w_offset 0 po_w_offset]
       catch
            w_offset_1 = [ne_w_offset' 0 po_w_offset']
       end
       mSigE_1(:,k) = spline(w_offset,mSigE(:,k),w_offset_1)';  
       [Row_mSigE_1 column_mSigE_1]=size(w_offset_1); 
       MTRasymM0_exp(:,k)=mSigE_1(1:(column_mSigE_1-1)/2,k)-mSigE_1(column_mSigE_1:-1:((column_mSigE_1+3)/2),k); 
     else 
       MTRasymM0_exp(:,k)=mSigE(1:(Nw-1)/2,k)-mSigE(Nw:-1:((Nw+3)/2),k); 
  end
   MTRM0_exp(:,k) = mSigE(Nw:-1:1,k);
end
if min(w_offset) < 0
     mSigM= mSigE(w_offset<0,:);
end
po_w_offset=flip(po_w_offset);

% %%%%% add by liuhao 
% %对8个ROI的Z谱数据归一化
% for k=1:nROI
%      mSigE(:,k) = mSigE(:,k)./max(mSigE(:,k));
%      mSig(:,k) = mSig(:,k)./max(mSig(:,k));
%      mSigM(:,k) = mSigM(:,k)./max(mSigM(:,k));
% end

if strfind(computer, 'WIN')
    if min(w_offset) < 0
        xlswrite( fullfile(ROIpath, ['MTRasym_',str,'.xls']),MTRasym_exp,N_echo,['C2']);
        xlswrite( fullfile(ROIpath, ['MTRasym_',str,'.xls']),MTRasym_as,N_echo,['C2']);
        xlswrite( fullfile(ROIpath, ['MTRasym_',str,'.xls']),MTRw',N_echo,['B2']);
        xlswrite( fullfile(ROIpath, ['MTRasymM0_',str,'.xls']),MTRasymM0_exp,N_echo,['C2']);
        xlswrite( fullfile(ROIpath, ['MTRasymM0_',str,'.xls']), po_w_offset',N_echo,['B2']);
        xlswrite( fullfile(ROIpath, ['MTRasymM0_as_',str,'.xls']),MTRasymM0_as,N_echo,['C2']);
        xlswrite( fullfile(ROIpath, ['MTRasymM0_as_',str,'.xls']),MTRe_as',N_echo,['B2']);
        xlswrite( fullfile(ROIpath, ['MTRm_',str,'.xls']),mSigM,N_echo,['C2']);
        xlswrite( fullfile(ROIpath, ['MTRm_',str,'.xls']),m_offset',N_echo,['B2']);
    end
        if size(w_offset,1)==1
            w_offset = w_offset';
        end


        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mSigE_0to1 = mSigE;
        mSig_0to1 = mSig;

        mSigE_Raw = mSigE;
        mSig_Raw = mSig;
        for k=1:nROI

        %对Z谱数据归一化
        mSigE(:,k) = mSigE(:,k)./max(mSigE(:,k));
        mSig(:,k) = mSig(:,k)./max(mSig(:,k));

            mSigE_0to1(:,k) = mSigE_0to1(:,k)-min(mSigE_0to1(:,k));
            mSigE_0to1(:,k) = mSigE_0to1(:,k)./max(mSigE_0to1(:,k))

            mSig_0to1(:,k) = mSig_0to1(:,k)-min(mSig_0to1(:,k));
            mSig_0to1(:,k) = mSig_0to1(:,k)./max(mSig_0to1(:,k))
        end
        
        %%归一化的Z谱
        save(fullfile(ROIpath,['Zspc_as_',str]),'mSig')
        save(fullfile(ROIpath,['Zspc_',str]),'mSigE')
        xlswrite( fullfile(ROIpath, ['Zspc_',str,'.xls']),w_offset,N_echo,['B1']);
        xlswrite( fullfile(ROIpath, ['Zspc_',str,'.xls']),mSigE,N_echo,['C1']);    
        xlswrite( fullfile(ROIpath, ['Zspc_as_',str,'.xls']),x_as',N_echo,['B1']);
        xlswrite( fullfile(ROIpath, ['Zspc_as_',str,'.xls']),mSig,N_echo,['C1']);
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%0to1的Z谱
        save(fullfile(ROIpath,['Zspc_as_',str,'_0to1']),'mSig_0to1')
        save(fullfile(ROIpath,['Zspc_',str,'_0to1']),'mSigE_0to1')
        xlswrite( fullfile(ROIpath, ['Zspc_',str,'_0to1.xls']),w_offset,N_echo,['B2']);
        xlswrite( fullfile(ROIpath, ['Zspc_',str,'_0to1.xls']),mSigE_0to1,N_echo,['C2']);    
        xlswrite( fullfile(ROIpath, ['Zspc_as_',str,'_0to1.xls']),x_as',N_echo,['B2']);
        xlswrite( fullfile(ROIpath, ['Zspc_as_',str,'_0to1.xls']),mSig_0to1,N_echo,['C2']);

        %%原始未做任何处理的Z谱
        save(fullfile(ROIpath,['Zspc_as_',str,'_Raw']),'mSig_Raw')
        save(fullfile(ROIpath,['Zspc_',str,'_Raw']),'mSigE_Raw')
        xlswrite( fullfile(ROIpath, ['Zspc_',str,'_Raw.xls']),w_offset,N_echo,['B1']);
        xlswrite( fullfile(ROIpath, ['Zspc_',str,'_Raw.xls']),mSigE_Raw,N_echo,['C1']);    
        xlswrite( fullfile(ROIpath, ['Zspc_as_',str,'_Raw.xls']),x_as',N_echo,['B1']);
        xlswrite( fullfile(ROIpath, ['Zspc_as_',str,'_Raw.xls']),mSig_Raw,N_echo,['C1']);

        
 elseif strfind(computer, 'MAC')    
     if min(w_offset) < 0
         csvwrite( fullfile(ROIpath, ['MTRasym_',str,'.xls']),MTRasym_exp,0,0);
         csvwrite( fullfile(ROIpath, ['MTRasym_',str,'.xls']),MTRe',0,0);
         csvwrite( fullfile(ROIpath, ['MTRasymM0_',str,'.xls']),MTRasymM0_exp,0,0);
         csvwrite( fullfile(ROIpath, ['MTRasymM0_',str,'.xls']),MTRe',0,0);
     end
     csvwrite( fullfile(ROIpath, ['Zspc_',str,'.csv']),mSigE,0,0);
end
str1=[num2str(k),' ROIZ谱图'];
for k=1:nROI
    figure
    str1=[num2str(k),' ROIZ谱图'];
%     plot(w_offset,mSigE(:,k),'b--o');
    plot(w_offset,mSigE_Raw(:,k),'r:.');
%     axis([-5 5,-inf,inf]);
    title(str1,'FontWeight','bold','FontSize',14);
    set(gca, 'Xdir', 'normal', 'FontWeight','bold','FontSize',14);
    savefig(fullfile(ROIpath,['Zspc',num2str(k),'ROI']))
%    if min(w_offset)<0   %显示MATRasym的，需要看的就是加上
%        figure(k+2)
%        str1=[num2str(k),'MTRasym_exp'];
%        plot(MTRe,MTRasym_exp(:,k),'b--o');
%        title(str1,'FontWeight','bold','FontSize',18);
%        set(gca, 'Xdir', 'reverse', 'FontWeight','bold','FontSize',18);
%        savefig(fullfile(ROIpath,['MTRasym_exp',num2str(k),'ROI']))
%    end
end
