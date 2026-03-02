function [Thmask] = GenerateThmask1107(Savepath,Thmask_datapath,S0,snr)
%[Thmask] = GenerateThmask(Savepath,S0,snr);
[Nx,Ny] = size(S0);
try 
    load(fullfile(Thmask_datapath, 'Thmask.mat'))
    Thmask = imresize(Thmask,size(S0));
catch
    Thmask = ones(Nx,Ny);
    temp= ThreMask_SNR(S0, snr);
    Thmask = Thmask.*temp;
    save(fullfile(Savepath, 'Thmask'),'Thmask');
end
end

function [image_mask]=ThreMask_SNR(ref_MI, snr)
[m,n]= size(ref_MI);
fh=figure;
imagesc(ref_MI)
colormap(gray)
disp('please draw the ROI of Noise...')
Noise_mask= roipoly;
index=find( Noise_mask >0);
Noise=ref_MI(index);
N = std(Noise) ;

image_mask= zeros(m,n);
[indX, indY]= find( ref_MI >snr*N) ;

for ind= 1: length(indX)
    image_mask(indX(ind), indY(ind)) = 1;
end

close(fh) 
end