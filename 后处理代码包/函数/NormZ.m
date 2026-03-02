function [V_exp_mask]=NormZ(S0,V_exp,brainMask) 
    S0_use=S0;
    Nw = size(V_exp,3);
    Row=size(V_exp,1);        
    Column=size(V_exp,2); 
    V_exp_mask=zeros(Row,Column,Nw);
    for n_offset=1:Nw
              m1=V_exp(:,:,n_offset)./S0_use;%对每个ppm的zdata的每个像素求z谱信号
              m1_ROI= m1.*brainMask;
              V_exp_mask(:,:,n_offset)=m1_ROI;%存每个像素的z谱信号
    end
% % %   save(fullfile(ROIpath,'V_exp_ROI'),'V_exp_ROI')    
% % %   save(fullfile(ROIpath,'w_offset'),'w_offset')      