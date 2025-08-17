function chosenFile = load_SSmatching_file(folderPath,targetSuffix)
    files         = dir(fullfile(folderPath,'*.mat'));
    matchingFiles = [];
    for i = 1:length(files)
        fname = files(i).name;
        % Check if the file matches the target suffix
        if endsWith(fname,targetSuffix)
            % Extract prefix before suffix
            prefix = erase(fname,targetSuffix);
            % Check if prefix is a valid date
            if length(prefix) == 9 && ~isempty(datetime(prefix,'InputFormat','MMMddyyyy','Locale','en_US'))
                matchingFiles = [matchingFiles;files(i)];
            end
        end
    end
    if isempty(matchingFiles)
        error('No matching file found.');
    elseif length(matchingFiles) > 1
        % Sort by date and pick the most recent
        [~, idx]   = max([matchingFiles.datenum]);
        chosenFile = matchingFiles(idx).name;
    else
        chosenFile = matchingFiles.name;
    end
end
