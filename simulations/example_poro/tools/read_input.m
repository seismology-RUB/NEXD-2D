% version='$Rev: 20 $ ($Date: 2017-06-07 18:22:15 +0200 (Mi, 07 Jun 2017) $, $Author: Janis Heuel, Marc S. Boxberg $)'
%
% Function to read the parameters from input files.
%   The following parameters are needed to use this function in the
%       main program:
%       filename:   contains relative path and name of the input file.
%       paramter:   contains the name of the input parameter
%       true_false: if the parameter is element of logical type true or
%                   false, this parameter has to be set as 1, else 0.
%       del:        if a special entry in filename should be deleted, this
%                   parameter has to be set as 1, else 0.
%       strdel:     if del=1 then strdel contains the name of the string
%                   that should deleted, else 0.
%
%   -> filename, parameter and strdel are both strings. This function
%      inputs has to be called with str2func('name') when calling the
%      function.
%
%       out: the output, which can further be used in main program.

function [out] = read_input(filename,parameter,true_false,del,strdel)
    filename = func2str(filename);
    parameter = func2str(parameter);
    workpath = pwd;
    read = table2cell(readtable(filename,'Delimiter','\t'));

    if del == 1
        strdel = func2str(strdel);
        find   = strfind(read,strdel);                        % finds the keyword strdel
        b_loop = 0;
        l_loop = 1;
        while b_loop ~= 1                                     % loop finds the line with the keyword
              a_loop  = cell2mat(find(l_loop));
              b_loop  = length(a_loop);
              l_loop  = l_loop+1;
        end
        read = strrep(read,read(l_loop-1),'');                % Deletes the entry
    end
% The following lines find the keyword and convert the input to
% logical type, number or string.
    find = strfind(read,parameter); % finds the keyword parameter
    b_loop = 0;
    l_loop = 1;
    while b_loop == 0                                                 % loop finds the line with the keyword
        a_loop  = cell2mat(find(l_loop));
        b_loop  = length(a_loop);
        l_loop  = l_loop+1;
    end
% Extract string and convert it to logical type
    convert = cell2mat(read(l_loop-1));
    find_first_sign   = strfind(convert,'=');
    del_first_sign    = convert(1+find_first_sign:end);
    find_second_sign  = strfind(del_first_sign,'#');
    del_second_sign   = del_first_sign(2:find_second_sign-2);
    if true_false == 1
        out_true    = strfind(del_second_sign,'.true.');
        out_false   = strfind(del_second_sign,'.false.');
        if out_true == 1
            out = 1;
        elseif out_false == 1
            out = 0;
        elseif out_true ~= 1 && out_false ~= 1
            out = 2;
        end
    elseif true_false == 0
        if strcmp(parameter,'extwavelet') == false
            out = str2double(del_second_sign);
        elseif strcmp(parameter,'extwavelet') == true
            out = cellstr(del_second_sign);

        end
    end
    cd (workpath)
    cd ../tools
end

