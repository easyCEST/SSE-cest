function [Vo_ori_reg] = registration(Vo_ori_use,S0,register,zspecDir)
[Nx,Ny,Nz1] = size(Vo_ori_use);
if register == 0 
    Vo_ori_reg = Vo_ori_use;
elseif register == 1 % 自带的rigid配准
        Vo_ori_reg = zeros(size(Vo_ori_use));
    transType = 'rigid';
    for ind = 1 : Nz1
        [tmp, ~, ~] = coRegister(S0, Vo_ori_use(:,:,ind), transType);
        Vo_ori_reg(:,:,ind) = tmp;
    end
elseif register == 2 % RPCA
    % RPCA
    mtype = 'ssd';  %mtype 和 ttype：这两个变量是设置后续图像配准方法的类型。'ssd'代表最小平方误差（Sum of Squared Differences），而'rigid'表示刚性变换（没有旋转或缩放，仅平移）
    ttype = 'rigid';
    depth = 4;  %设定了图像配准的深度，通常是用于图像金字塔（multi-resolution）策略，逐渐细化配准。
    %save I0 in a 2D matrix M(Nx*Ny Nz1)
    M = zeros(Nx*Ny,Nz1);
    for iz = 1:Nz1
         Vo = Vo_ori_use(:,:,iz);   %每一列是每个信号的值
         M(:,iz) = Vo(:);
    end
    [m,n] = size(M);
    if exist(fullfile(zspecDir,['lambda.mat']),'file')~=2
        i = 1;  %记录当前的迭代次数
        lambda = zeros(Nz1,1);  %存储每次迭代中计算得到的正则化参数
        rankL = zeros(Nz1,1);   %存储每次迭代中低秩矩阵的秩（rank）
        non0S = 1;   %记录稀疏矩阵（sparse）中非零元素的比例。初始化为 1（表示最大比例），以确保进入 while 循环。
         while non0S > 0.1  %迭代计算正则化参数 λ
                if i == 1
                    lambda(i) =1 / sqrt(max(m,n))/5;   %这里的初始值计算方式通常是为了避免过大或过小的正则化值，保证初步的正则化适中。
                else
                    lambda(i) = (1 + 0.05) *lambda(i - 1);  %对于后续的迭代，lambda(i) 会逐步增加，每次增加 5%（(1 + 0.05) * lambda(i - 1)）。这种递增方式是为了使正则化越来越强，从而使稀疏部分的非零元素逐渐减少。
                end

                % RPCA and claculate the no 0 of sparse
                [lowrank,sparse] = RPCA_iALM(M,lambda(i));
                rankL(i) = rank(lowrank);   %记录低秩矩阵的秩
                non0S = sum(sum(sparse ~= 0))/(numel(sparse));   %计算稀疏矩阵的非零元素比例
                i = i + 1;
         end
         lambda(find(lambda==0)) = [];
         rankL(find(rankL==0)) = [];

         save(fullfile(zspecDir,['rankL.mat']),'rankL');
         save(fullfile(zspecDir,['lambda.mat']),'lambda');
    else
        load(fullfile(zspecDir,['lambda.mat']));
        load(fullfile(zspecDir,['rankL.mat']));
    end  
%     figure(1)
%     rawlambda = 1 / sqrt(max(m,n));
%     plot(lambda / rawlambda,rankL,'-o')
%     savefig(fullfile(zspecDir,'lambda-rankL.fig'))

    % registration among the range of lambda_use
    [~,rank3]=min(abs(rankL-3));
    lam0 = lambda(rank3);
    lambda_use =lam0;   %选择最优的正则化参数lambda

    [lowrank, sparse] = RPCA_iALM(M,lambda_use);    
    L = zeros(size(Vo_ori_use));
    for iy = 1:n
        tmp = lowrank(:,iy);
        L(:,:,iy) = reshape(tmp,Nx,Ny);
    end

    % 3. compute Reference Image: mean of lowrank
    L_Ref = mean(lowrank(:,:),2);
    L_Ref = reshape(L_Ref,Nx,Ny);

    %4. Registration
    Vo_ori_reg = zeros(size(Vo_ori_use));
    for ind = 1:Nz1
        tmp = image_registration_pairwise(Vo_ori_use(:,:,ind),L_Ref,mtype,ttype,depth);
        Vo_ori_reg(:,:,ind) = tmp;
    end
    save(fullfile(zspecDir,['Vo_ori_reg.mat']),'Vo_ori_reg');
    save(fullfile(zspecDir,['lowrank.mat']),'lowrank');
    save(fullfile(zspecDir,['sparse.mat']),'sparse');
end
