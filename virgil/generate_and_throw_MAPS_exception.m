function generate_and_throw_MAPS_exception(errId, varargin)
% Wrapper to create and immediately throw a MAPS exception.
% Accepts either:
%   1) generate_and_throw_MAPS_exception(errId)
%   2) generate_and_throw_MAPS_exception(errId, errArgsCell)

    % Handle optional errArgs
    if nargin < 2
        errArgs = {};
    else
        errArgs = varargin{1};
        if ~iscell(errArgs)
            errArgs = {errArgs};
        end
    end

    % Create the exception and throw it
    exc = generate_MAPS_exception(errId, errArgs);
    throw(exc);
end
