function display_image_save1103(Savepath, mask, MatName, FigTitle, colorRange)
    % 如果没有传入colorRange参数，则使用默认值[0 0.35]
    if nargin < 5
        colorRange = [0 0.35]; 
    end
    
    index_ROI = find(mask < 1);
    MatName(index_ROI) = NaN;
    [px,py] = find(mask);
    
    figure
    imagesc(MatName(min(px):max(px), min(py):max(py)), colorRange); % 使用传入的颜色范围
    colormap(jet(256))
    colorbar
    altered_colormap()
    axis off
    set(gca, 'FontWeight', 'bold', 'FontSize', 14)
    title(FigTitle, 'FontWeight', 'bold', 'FontSize', 14)
    savefig(strcat(Savepath, '\', FigTitle, '.fig'))
end
