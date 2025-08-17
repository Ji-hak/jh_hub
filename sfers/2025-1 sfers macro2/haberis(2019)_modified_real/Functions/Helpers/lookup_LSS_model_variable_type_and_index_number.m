function [varTypeStr,varIndexNumber] = ...
    lookup_LSS_model_variable_type_and_index_number(Model,varMnem)
% This function looks-up an LSS model variable type and index number.
%
% INPUTS:
%   -> Model: LSS model structure
%   -> varMnem: model variable, model observable or raw observable
%
% OUTPUTS:
%   -> varTypeStr: identifier for the type of variable consistent with
%      wider MAPS use: 'modelVariables' or 'modelObservables' or
%      'timeVaryingTrends' or 'Shocks'
%   -> varIndexNumber: integer for position of variable within type 
%
% DETAILS:
%   -> This model utility function looks-up an LSS model variable type and
%      index number. It is useful when used in conjunction with code that
%      allows users to identify variables by mnemonics alone.
%
% NOTES:
%   -> See the MAPS user guide for more information on model analysis 
%      functionality in MAPS.
%
% This version: 05/07/2013
% Author(s): Rich Harrison Matt Waldron

%% CHECK INPUTS
if nargin < 2
    
elseif ~isstruct(Model)
    
elseif ~ischar(varMnem)
    
end

%% GATHER RELEVANT INFORMATION
[xMnems,zMnems,modelHasMeasurementEqs,modelHasDataTransformationEqs,...
    modelHasTimeVaryingTrends] = unpack_model(Model,{...
    'xMnems','zMnems','modelHasMeasurementEqs',...
    'modelHasDataTransformationEqs','modelHasTimeVaryingTrends'});

%% WORK OUT VARIABLE TYPE
% Note the hierachical logic of this cell. Models without measurement
% equations cannot have data transformation equations. And models with data
% transformation equations cannot have time-varying trends.
if any(ismember(varMnem,xMnems))
    varTypeStr = 'modelVariables';
    varTypeMnems = xMnems;
elseif any(ismember(varMnem,zMnems))
    varTypeStr = 'Shocks';
    varTypeMnems = zMnems;
elseif modelHasMeasurementEqs 
    Ymnems = unpack_model(Model,'Ymnems');
    if any(ismember(varMnem,Ymnems))
        varTypeStr = 'modelObservables';
        varTypeMnems = Ymnems;
    elseif modelHasDataTransformationEqs
        YtildeMnems = unpack_model(Model,'YtildeMnems');
        if any(ismember(varMnem,YtildeMnems))
            varTypeStr = 'rawObservables';
            varTypeMnems = YtildeMnems;
        elseif modelHasTimeVaryingTrends
            ttMnems = unpack_model(Model,'etatMnems');
            if any(ismember(varMnem,ttMnems))
                varTypeStr = 'timeVaryingTrends';
                varTypeMnems = ttMnems;
            else
                error(['The variable mnemonics passed in as input to '...
                    'this function does not exist in the model input'])
            end
        end
    end
end

%% LOOKUP INDEX NUMBER
varIndexNumber = lookup_model_index_numbers(varTypeMnems,varMnem);

end