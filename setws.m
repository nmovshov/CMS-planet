function setws(dbg)
% Set workspace (paths and/or common variables).

if nargin == 0, dbg = false; end

addpath(pwd)
addpath(fullfile(pwd, 'extern'))
addpath(fullfile(pwd, 'extern', 'physunits'))

if dbg % use dimensionally aware preal class to help formula debugging
    evalin('base', 'si = setUnits;')
else   % use purely numerical funits for speed
    evalin('base', 'si = setFUnits;')
end
