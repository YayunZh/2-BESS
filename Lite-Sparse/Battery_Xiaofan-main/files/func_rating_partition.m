function [conv_rating_array, conv_rating_matrix] = func_rating_partition(conv_power_matrix,partition_num)
    % initializations
    conv_rating_matrix = zeros(size(conv_power_matrix,1),size(conv_power_matrix,2));% 10*10
    conv_power_array = [];
    for i = 1:size(conv_power_matrix,1)
        for j = i:size(conv_power_matrix,2)
            if (conv_power_matrix(i,j)~=0)
                conv_power_array(end+1) = conv_power_matrix(i,j);
            end    
        end
    end
    % sortings
    sorted_power = sort(conv_power_array);
    conv_num = length(conv_power_array); % 3
    devide_index_array = linspace(0,conv_num,partition_num+1);% 0 3
    devide_index_array_ceil = ceil(devide_index_array(2:end));% 3
    conv_rating_array = zeros();
    
    for i = 1:conv_num
        for j = 1:partition_num
            if(i <= devide_index_array_ceil(j))
                conv_rating_array(i) = sorted_power(devide_index_array_ceil(j));
                break;
            end
        end
    end
    for i = 1:size(conv_power_matrix,1) % input, 10
        for k = 1:size(conv_power_matrix,2)
            if (conv_power_matrix(i,k) == 0)
                conv_rating_matrix(i,k) = 0;
            else
                for j = 1:partition_num
                    if(abs(conv_power_matrix(i,k)) <= (sorted_power(devide_index_array_ceil(j))+1e-6))
                        conv_rating_matrix(i,k) = sorted_power(devide_index_array_ceil(j));
                        break;
                    end
                end
            end
        end
    end
end

