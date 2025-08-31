function SSfile = load_SSmatching_file(folderPath,targetSuffix)
    files         = dir(fullfile(folderPath,'*.mat'));
    matchingFiles = [];
    for i = 1:length(files)
        fname = files(i).name;
        % Check if the file matches the target suffix
        if endsWith(fname,targetSuffix)
            prefix = erase(fname,targetSuffix);
            if length(prefix) == 9 && ~isempty(datetime(prefix,'InputFormat','MMMddyyyy','Locale','en_US'))
                matchingFiles = [matchingFiles;files(i)];
            end
        end
    end
    if isempty(matchingFiles)
        error('No matching file found.');
    elseif length(matchingFiles) > 1
        [~, idx]   = max([matchingFiles.datenum]);            % Sort by date and pick the most recent
        SSfile = matchingFiles(idx).name;
    else
        SSfile = matchingFiles.name;
    end
end
