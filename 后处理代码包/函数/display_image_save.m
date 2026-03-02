function display_image_save(Savepath,mask,MatName,FigTitle)
    index_ROI = find(mask < 1);
%     MatName(index_ROI) = -0.5;
    MatName(index_ROI) = NaN;
    [px,py]=find(mask);
    figure
    imagesc(MatName(min(px):max(px),min(py):max(py)),[0 0.35]);
    colormap(jet(256))
    colorbar
    altered_colormap()
    axis off
    set(gca, 'FontWeight','bold','FontSize',14)
    title(FigTitle,'FontWeight','bold','FontSize',14)
    savefig(strcat(Savepath,'\',FigTitle,'.fig'))
end
