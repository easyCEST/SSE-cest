function n = get_closest_index(diff_array, voxel_LD, add_num)
    % 功能：根据输入的差值数组（如w_offset_1 - 3.5）的绝对值排序，返回满足n+add_num不超出voxel_LD第三维度范围的索引n
    % 输入：
    %   diff_array - 差值数组（如w_offset_1 - 目标值，由外部计算后传入）
    %   voxel_LD - 用于判断范围的三维数组
    %   add_num - 指定需判断的"n+add_num"中的偏移量
    % 输出：
    %   n - 符合条件的索引
    
    % 按差值的绝对值排序，获取索引（接近程度排序）
    [~, idx] = sort(abs(diff_array));
    
    % 初始化n为最接近的索引
    n = idx(1);
    
    % 检查是否越界，若越界则依次选取下一个接近的索引
    if n + add_num > size(voxel_LD, 3)
        n = idx(2);
        if n + add_num > size(voxel_LD, 3)
            n = idx(3);
        end
    end
end