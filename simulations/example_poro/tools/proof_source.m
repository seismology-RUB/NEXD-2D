% version='$Rev: 20 $ ($Date: 2017-06-07 18:22:15 +0200 (Mi, 07 Jun 2017) $, $Author: Janis Heuel, Marc S. Boxberg $)'
% Proofs if the source was calculated before. If cmp_elements is false,
% then the source exists twice or more, else the answer is true, then
% the source is used before.

sourcematrix = [frequency_vec, stf];
sort_extw = find(stf==4);

if src == 1
    for i = 1:length(stf)
        if stf(i) ~= 4
            extwavelet{i} = 0;
        end
    end
end

if src >= 2
    source_proof = [frequency_vec(1:src-1), stf(1:src-1)];
    for j = 1:src-1
        cmp_elements = find(sourcematrix(src,:) == source_proof(j,:));
        if wavelet_type ~= 4
            if length(cmp_elements) == 2
                cmp_elements = false;
                break
            elseif isempty(cmp_elements)
                cmp_elements = true;
            end
        elseif wavelet_type == 4
            cmp_str = zeros(1,j);
            for k = 1:j
                cmp_str(k) = double(strcmp(extwavelet{k},extwavelet{src}));
            end
            if norm(cmp_str) == 0
                cmp_elements = true;
            elseif norm(cmp_str) ~= 0
                cmp_elements = false;
            end
        end
    end
end

if wavelet_type == 4 && cmp_elements == false
    fprintf('The external wavelet is used before. \n');
end

cd (workpath)
cd ../data
