


%% Subfunction: get field if exist, return default value otherwise
function val = tryGetField(s, field, dftVal)
if isfield(s, field), val = s.(field); 
elseif nargin>2, val = dftVal;
else, val = [];
end
