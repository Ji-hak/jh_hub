%% RUN ALL SCRIPTS
profile on;
generateSection3results;
generateSection4results;
plotFigures2and3B1B3B4;
plotFigures4C4C5;
plotFigureB2;
plotFigureC1;
plotFigureC2;
plotFigureC3;
computeTBFGlosses;
profile off;

%% HOUSEKEEPING
close all;
restoredefaultpath;
clear variables;
delete('*.asv');
clc; 

%% CHOOSE TO DELETE UNUSED FUNCTIONS (OR JUST REPORT)
deleteUnusedFuncs = false;

%% ADD DIRECTORIES
fullPathNameForThisFile = mfilename('fullpath');
fullPathFileNameSplitByFinalBackSlash = regexp(...
    fullPathNameForThisFile,'[^\\]+$','split');
thisDir = fullPathFileNameSplitByFinalBackSlash{1};
addpath(genpath(thisDir));

%% DIRECTORIES TO CHECK THE FUNCTION CALLS
dirsToCheck = {[thisDir,'MAPSlite'],[thisDir,'Functions']};

%% EXTRACT INFO ABOUT THE FULL PATH FILE NAMES OF THE FUNCTIONS CALLED
% Get profile info.
P = profile('INFO');
fullPathFileNames = {P.FunctionTable.FileName}';
functionNames = {P.FunctionTable.FunctionName}';
% Find function calls for files living in the directories specified above
patMatchCellStr = strrep(strcat('^(',dirsToCheck,')'),'\','\\');
nDirsToCheck = size(dirsToCheck,2);
patMatchStr = patMatchCellStr{1};
for iDir = 2:nDirsToCheck
    patMatchStr = [patMatchStr,'|',patMatchCellStr{iDir}];                  %#ok<AGROW>
end
nFuncCalls = size(fullPathFileNames,1);
isFuncCallInSpecifiedDir = false(nFuncCalls,1);
for iCall = 1:nFuncCalls
    if ~isempty(regexp(fullPathFileNames{iCall},patMatchStr,'match'))
        isFuncCallInSpecifiedDir(iCall) = true;
    end
end
namesOfCalledFuncsIncSubFuncs = functionNames(isFuncCallInSpecifiedDir);
% Remove the sub function calls to leave just the function names
nFuncsIncSubFuncs = size(namesOfCalledFuncsIncSubFuncs,1);
isCallAsubFunc = false(nFuncsIncSubFuncs,1);
for iCall = 1:nFuncsIncSubFuncs
    if ~isempty(regexp(namesOfCalledFuncsIncSubFuncs{iCall},'>','match'))
        isCallAsubFunc(iCall) = true;
    end
end
namesOfCalledFuncs = namesOfCalledFuncsIncSubFuncs(~isCallAsubFunc);

%% FIND FUNCTIONS IN THE ABOVE DIRECTORIES THAT NOT CALLED IN THE CODE
% Create a list of all subdirectories
allFoldersToCheck = {};
for iDir = 1:nDirsToCheck
    iDirAllFoldersStr = genpath(dirsToCheck{iDir});
    iDirAllFoldersCellStr = regexp(iDirAllFoldersStr,';','split');
    iDirAllFoldersCellStr = iDirAllFoldersCellStr(1:end-1);
    allFoldersToCheck = [allFoldersToCheck,iDirAllFoldersCellStr];           %#ok<AGROW>
end
allFoldersToCheck = allFoldersToCheck';
% Loop through each folder
nFoldersToCheck = size(allFoldersToCheck,1);
for iFolder = 1:nFoldersToCheck
    % Extract complete list of .m functions in that folder
    iFolderToCheck = allFoldersToCheck{iFolder};
    fprintf(['Checking: ',strrep(iFolderToCheck,'\','\\'),'\n']);
    iFolderFuncInfo = dir([iFolderToCheck,'\*.m']);
    % Loop through each function found
    nFuncsToCheck = size(iFolderFuncInfo,1);    
    for iFunc = 1:nFuncsToCheck
        iFuncToCheck = strrep(iFolderFuncInfo(iFunc).name,'.m','');
        % Check function has been called in the code
        if ~ismember(iFuncToCheck,namesOfCalledFuncs)
            % If not report
            fprintf(['\t',iFuncToCheck,' is unused\n']);
            % And delete as appropriate
            if deleteUnusedFuncs
                fprintf('\t\t deleting ...\n');
                fullPathNameOfFileToDelete = ...
                    [iFolderToCheck,'\',iFuncToCheck,'.m'];
                delete(fullPathNameOfFileToDelete);
            end     
        end
    end
end