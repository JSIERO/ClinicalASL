function ChangeElastixParameterFileEntry(filename,oldTextLine,newTextLine, newfilename )

% Define file path and text to search for
%filename =  % Replace with your actual file name
%oldTextLine = % The text you want to replace
%ewTextLine = ' % The text you want to replace it with
%newfilename = new name for modified file

% Open the file for reading
fid = fopen(filename, 'r');

if fid == -1
    error('File not found or permission denied');
end

% Read the file line by line into a cell array
fileContents = {};
lineIndex = 1;
while ~feof(fid)
    fileContents{lineIndex} = fgetl(fid);
    lineIndex = lineIndex + 1;
end

% Close the file after reading
fclose(fid);

% Replace the old text line with the new one
for i = 1:length(fileContents)
    if strcmp(fileContents{i}, oldTextLine)
        fileContents{i} = newTextLine;
    end
end

% Open the new file for writing
fid = fopen(newfilename, 'w');

if fid == -1
    error('File cannot be opened for writing');
end

% Write the modified content back to the file
for i = 1:length(fileContents)
    fprintf(fid, '%s\n', fileContents{i});
end

% Close the file after writing
fclose(fid);

disp('Text replacement completed successfully.');