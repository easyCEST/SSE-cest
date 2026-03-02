function altered_colormap()
cmap = colormap(jet);
nzeros = 2;
cmap(1:nzeros,:) = zeros(nzeros,3);
colormap(cmap); 