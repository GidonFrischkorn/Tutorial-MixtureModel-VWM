%% Setup
% Add folders with functions
addpath('MATLAB_data','data');

% specfiy experiments
experiments = {
    'E1' ; 'E2'; 'E3'; 'E4'; 'E5';
    'E6';'E7';'E8';'E9';'E10'
    };

% get all expfiles
expfiles = dir('MATLAB_data\');
expfiles_names   = {expfiles.name};

% Loop thourgh experiments
for e = 1:size(experiments,1)
    exp = experiments{e};

    this_expFiles_IDs = regexp(expfiles_names, exp, 'start');
    this_expFiles_IDs  = ~cellfun('isempty', this_expFiles_IDs);

    print_header = 0;

    for d = 1:size(expfiles,1)
        if this_expFiles_IDs(d) == 0
            continue
        else
            data = load([expfiles(d).folder,'\',expfiles(d).name]);

            subID = erase(expfiles(d).name, exp);
            subID = erase(subID, '_subject_');
            subID = erase(subID, '.mat');
            subID = str2double(subID);

            expName = data.experiment_name{1};
            expName = erase(expName,' et al. ');
            expName = erase(expName,' ');
            subData = data.data;

            dataFile = ['data\', expName, '.txt'];

            if print_header == 0
                fid = fopen(dataFile, 'a');
                % Print header into data file
                fprintf (fid, '%s %s %s %s ', 'subID', 'trial', 'setsize', 'RespErr');
                for i = 1:10
                    % Print retrieval specific inforamtion into data file
                    fprintf (fid, '%s ', ['Pos_Lure',num2str(i)]);
                end
                % Print new Line after each trial
                fprintf (fid, '\n');
                fclose (fid);

                % Set print_header to 1
                print_header = 1;
            end



            for t = 1:size(subData.error_vec,2)
                error = subData.error_vec(t);
                setsize = subData.N(t);
                lures = subData.dist_error_vec{t};

                fid = fopen(dataFile, 'a');
                % Print general information (constant over Trial(expTrial)) into data file
                fprintf (fid, '%d %i %i %d ', subID, t, setsize, error);
                for i = 1:10
                    if i > size(lures,1)
                        % Print NA for lure position into data file
                        fprintf (fid, '%s ', 'NA');
                    else
                        % Print lure positions into data file
                        fprintf (fid, '%d ', lures(i));
                    end
                end
                % Print new Line after each trial
                fprintf (fid, '\n');
                fclose (fid);
            end
        end
    end
end