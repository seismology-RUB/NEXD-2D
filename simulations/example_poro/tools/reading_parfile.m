% Script to read the parameters from parfile, to read the number of
% sources and the wavelet type from 'source' and to read the number of
% materials from 'matpropporo'.
% The function read_input.m is used to read from different files the
% parameter. For further information see description in read_input.m
%
% version='$Rev: 20 $ ($Date: 2017-06-07 18:22:15 +0200 (Mi, 07 Jun 2017) $, $Author: Janis Heuel, Marc S. Boxberg $)'

cd (workpath)

%% Checking if calculation is elastic or poroelastic.
poroelastic = read_input(str2func('../data/parfile'),str2func('poroelastic'),1,1,str2func('Poroelasticity'));
fprintf('\n')
if poroelastic == 1
    fprintf('Calculation is poroelastic. \n')
elseif poroelastic == 0
    fprintf('Calcualtion is elastic. \n')
else
    fprintf('Error in parfile at point poroelastic!\n');
    cd (workpath)
    cd ../tools
    return
end


%% Reading the number of fluids
number_of_fluids = read_input(str2func('../data/parfile'),str2func('fluidn'),0,0,0);
if ~ ismember(number_of_fluids,[0,1,2])
    fprintf('Number of fluids is greater than 2. The script stops here. \n')
    cd (workpath)
    cd ../tools
    return
end
number_of_fluids_save = number_of_fluids;
fprintf('\n');
fprintf('Number of fluids: %g \n',number_of_fluids);


%% Checking if tortuosity is given in 'matpropporo' or if tortuosity
%% will be calculated by a further parameter r.
tortuosity = read_input(str2func('../data/parfile'),str2func('calculate_tortuosity'),1,0,0);
if tortuosity == 1
    calculate_tortuosity = true;
elseif tortuosity == 0
    calculate_tortuosity = false;
else
    fprintf('calculate_tortuosity has to be set as true or false in parfile. \n');
    cd (workpath)
    cd ../tools
    return
end


%% Checking if fractional step is used for the calculation
%fractional_step = read_input(str2func('../data/parfile'),str2func('fractional_step'),1,0,0);
fractional_step = false;
if ~ ismember(fractional_step,[0,1])
    fprintf('fractional_step has to be set as true or false in parfile. \n');
    cd (workpath)
    cd ../tools
    return
end


%% Reading number of sources from 'source'
nsrc = read_input(str2func('../data/source'),str2func('nsrc'),0,0,0);
fprintf('Number of sources: %g \n', nsrc);

%% Indetify center frequency (frequency_vec) and type of wavelet (stf) in 'source'.
%% Both frequency_vec and stf are vectors with length nsrc.
%% extwavelet contains the string with the path of each external wavelet.
frequency_vec = zeros(nsrc,1);
stf           = zeros(nsrc,1);

copyfile('../data/source','../data/source_tmp')
for f = 1:nsrc
    stf(f)           = read_input(str2func('../data/source'),str2func('stf'),0,0,0);
    frequency_vec(f) = read_input(str2func('../data/source'),str2func('f0'),0,0,0);
    extwavelet{f}    = read_input(str2func('../data/source'),str2func('extwavelet'),0,0,0);

    %% Deletes the lines in 'source' which are unimportant to read the next source
    if f <= nsrc-1
        source_str       = sprintf('Source %d', f+1);
        read_source      = table2cell(readtable('../data/source'));
        find_line_source = strfind(read_source,source_str);
        b_loop = 0;
        l_loop = 1;
        while b_loop ~= 1 % finds the line with the keyword 'source'
              a_loop  = cell2mat(find_line_source(l_loop));
              b_loop  = length(a_loop);
              l_loop  = l_loop+1;
        end
        read_source = read_source((l_loop-1):end);
        fileID = '../data/source';
        writetable(cell2table(read_source),fileID)
        movefile('../data/source.txt','../data/source')
    end
end

copyfile('../data/source_tmp','../data/source')


%% Reading number of materials from 'matpropporo'
cd (workpath)
cd ../cubit
read_matpropporo    = table2cell(readtable('matpropporo','Delimiter','\t'));
find_begin          = find(strcmp(read_matpropporo, 'BEGIN'));
find_end            = find(strcmp(read_matpropporo, 'END'));
matpropporo         = read_matpropporo(find_begin+1:find_end-1);
number_of_materials = str2double(cell2mat(matpropporo(1)));
if number_of_materials <= 0
    fprintf('Number of materials is zero. The script stops here! \n');
    cd (workpath)
    cd ../tools
    return
elseif number_of_materials > 0
    fprintf('Number of materials: %g \n', number_of_materials);
end
