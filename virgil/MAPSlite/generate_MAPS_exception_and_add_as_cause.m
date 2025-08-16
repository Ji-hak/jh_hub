function updatedExc = generate_MAPS_exception_and_add_cause(masterExc, errId, varargin)
% Wrapper to add a new MAPS exception as a cause to an existing exception.
% Usage:
%   updatedExc = generate_MAPS_exception_and_add_cause(masterExc, errId)
%   updatedExc = generate_MAPS_exception_and_add_cause(masterExc, errId, errArgsCell)

    % Validate inputs
    if nargin < 2
        error('generate_MAPS_exception_and_add_cause:BadNargin', ...
              'Requires at least two inputs: master exception and errId.');
    end
    if ~isa(masterExc, 'MException')
        error('generate_MAPS_exception_and_add_cause:BadInput1', ...
              'First input must be an MException object.');
    end

    % Handle optional errArgs
    if nargin < 3
        errArgs = {};
    else
        errArgs = varargin{1};
        if ~iscell(errArgs)
            errArgs = {errArgs};
        end
    end

    % Generate the new exception and add as cause
    newExc = generate_MAPS_exception(errId, errArgs);
    updatedExc = addCause(masterExc, newExc);
end
