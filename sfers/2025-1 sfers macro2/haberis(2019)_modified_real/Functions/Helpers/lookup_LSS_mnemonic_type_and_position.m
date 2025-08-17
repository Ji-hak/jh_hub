function variableInfo = ...
    lookup_LSS_mnemonic_type_and_position(Model,mnemonicList)
% This function analyses the structure of an LSS model and extracts
% information about the position of mnemonics in a forecast run data
% structure.
%
% INPUTS
% Model: LSS model
% mnemonicList: nV x 1 cell vector of strings containing valid mnemonics 
% for variables in the model
%
% OUTPUTS
% variableInfo: nV x 3 cell array of information about the variables. First
% column is the mnemonic, second column is the type (eg "rawObservables"),
% third column is the index number within that type.
%
% Author: Richard Harrison
% Last modified: 18 April 2013

%% INITIAL CHECKS AND VALIDATION


%% GATHER INFORMATION ABOUT LSS MODEL
[xMnems,zMnems,modelHasMeasurementEqs,modelHasDataTransformationEqs] = ...
    unpack_model(Model,{'xMnems';'zMnems';'modelHasMeasurementEqs';...
    'modelHasDataTransformationEqs'});
if modelHasMeasurementEqs
    Ymnems = unpack_model(Model,'Ymnems');
end
if modelHasDataTransformationEqs
    YtildeMnems = unpack_model(Model,'YtildeMnems');
end
%% INITIALISE OUTPUTS
nVars = size(mnemonicList,1);
variableInfo = cell(nVars,3);

%% LOOP THROUGH VARIABLES AND POPULATE CELL ARRAY
for iVar = 1:nVars
    iMnem = mnemonicList{iVar};
    rowToSelect = nan;
    if ~isempty(intersect(iMnem,xMnems))
        rowToSelect = lookup_model_index_numbers(xMnems,iMnem);
        varType = 'modelVariables';
    elseif ~isempty(intersect(iMnem,zMnems))
        varType = 'Shocks';
        rowToSelect = lookup_model_index_numbers(zMnems,iMnem);
    elseif modelHasMeasurementEqs
        if ~isempty(intersect(iMnem,Ymnems))
            varType = 'modelObservables';
            rowToSelect = lookup_model_index_numbers(Ymnems,iMnem);
        end
        if modelHasDataTransformationEqs && ~isempty(intersect(iMnem,YtildeMnems))
            varType = 'rawObservables';
            rowToSelect = lookup_model_index_numbers(YtildeMnems,iMnem);
        end
    end
    if isnan(rowToSelect)
        error('Variable mnemonic not found as an endogenous variable in model structure.');
    end
    variableInfo{iVar,1} = iMnem;
    variableInfo{iVar,2} = varType;
    variableInfo{iVar,3} = rowToSelect;
end

