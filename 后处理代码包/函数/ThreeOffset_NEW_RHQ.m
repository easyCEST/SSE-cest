function [NOE_3pt] = ThreeOffset_NEW_RHQ(V_as,x_as,x_lw,Lp,Mp,Rp,S0_temp,Thmask)
%ThreeOffset_NEW_RHQ(V_as,x_as,x_lw,-2,-0.8,0.3,S0_temp,Thmask);    
    [~,x1] = min(abs(x_as-(Lp)));
    [~,x2] = min(abs(x_as-(Rp)));
    [~,x3] = min(abs(x_as-(Mp)));
    V_as_temp = V_as./(S0_temp+1e-5).*Thmask;
    %% 求Sref
    D = V_as_temp(:,:,x1).*abs(Lp-Rp)./(V_as_temp(:,:,x1)-V_as_temp(:,:,x2));%所有点的D：大三角x
    % S1 = D.*V_as(:,:,x1)./2;%所有点的大三角面积
    d = D - abs(Lp-Mp);%小三角的x
    H0p5 = d.*(V_as_temp(:,:,x1)-V_as_temp(:,:,x2))./abs(Lp-Rp);%每个像素点-0.5ppm的Zref
    H2 = V_as_temp(:,:,x1);%每个像素点-2ppm的Zref
    Sref = (H0p5+H2).*abs(Lp-Mp)./2;%每个像素点Zref的面积

    %% 求Srel
    for ii = 1:size(V_as,1)
        for jj = 1:size(V_as,2)
            if x1 < x3
                Srel(ii,jj,:) = trapz(Lp:x_lw:Mp,squeeze(V_as_temp(ii,jj,x1:x3)));
            else
                Srel(ii,jj,:) = trapz(Mp:x_lw:Lp,squeeze(V_as_temp(ii,jj,x3:x1)));
            end
        end
    end

    %% 求NOE_3pt
    NOE_3pt = Sref - Srel;