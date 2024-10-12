function assignUserInputs(funcParamsNames, varargin)
%% assignUserInputs
% Load user inputs into calling environment variables/parameters easily.
%
%% Syntax
% assignUserInputs({'param1', 'param1', 'param2'}, StructureOfInputs, 'param1',  val,...
%   'param2',  val)
%
%% Description
% The goal of this function is allow easy and intuitive parsing of user
%   inputs. The function supports three common way to pass variables to
%   functions:
%     (1) "Name value" pairs- the "variable name" string is followed by the
%           "variable value". The inputs orders is arbitrary, and the user
%           can specify the relevant inputs.
%     (2) Structure input- the structure field names will specify the variable name, and
%           field content is the variable value.
%     (3) Regular case, where variables name is defined by their position in function call.  
%   Both methods "name-value" pairs and structure can be combined. User are advised to
%   load default values, before calling this function, to guarantee all variable are
%   initialized. 
%   User inputs will be interpreted by the default Matlab way, only if none
%   of the inputs was interpreted by the other two methods.
%

%% Example I- transfering variables via structure and "name value" pairs
% clear all;
% S=struct( 'param1', true, 'param2', 'bla bla bla', 'param3', {10} );
% funcParamsNames=fieldnames(S);
% funcParamsNames=cat(1, funcParamsNames, 'param4');
% who
% assignUserInputs(funcParamsNames, S, 'param4', rand(3,2) );
% fprintf('Now see the newly added variables: %s, %s, %s, %s.\n', funcParamsNames{:});

%% Load uses params, overifding default ones

% The list of all legal input parameter names. Others will be ignored
if ~(iscell(funcParamsNames) && all( cellfun(@ischar, funcParamsNames) ))
    % if no cell array of string was specified as funcParamsNames input, consider it to
    % be empty, and append varargin
    varargin=cat(2, funcParamsNames, varargin); % varargin is a 1xN cell array
    funcParamsNames=[];
end

if numel(varargin)==1 && isempty(varargin{1}) % functioned called wihout arguments
    return;
end

% verify if no funcParamsNames input was specified
if isempty(funcParamsNames)
    isNoFuncParamsNames=true;
else
    if verLessThan('matlab', '7.14') % again, old version do not support 'stable'. 
        funcParamsNames=unique(funcParamsNames); % This can lead to bugs :(
    else
        funcParamsNames=unique(funcParamsNames , 'stable');
    end%
    isNoFuncParamsNames=false;
end

unUsedVarargin=varargin;
isUsedParamsName=false( size(funcParamsNames) );
%% Load "param name"- "param value" pairs.
nVarargin=length(varargin);
if nVarargin > 1
    isSpecificParam=false(1, nVarargin);
    iArg=1;
    while iArg <= (nVarargin-1)
        % automatically get all RELEVANT input pairs and store in local vairables
        isFuncParamsNames=strcmpi(varargin{iArg}, funcParamsNames);
        if isNoFuncParamsNames || any( isFuncParamsNames  )
            assignin('caller', varargin{iArg}, varargin{iArg+1});
            
            isSpecificParam( [iArg, iArg+1] )=true; % parameters come in Name-Value pairs
            iArg=iArg+1;
            isUsedParamsName(isFuncParamsNames)=true;
        end
        iArg=iArg+1;
    end % while iArg < (nVarargin-1)
    unUsedVarargin=varargin(~isSpecificParam); % Save varargin elements that were not used
    funcParamsNames=funcParamsNames(~isUsedParamsName); % do not allow repeating 
                                                        % initialisation of variables
    isUsedParamsName=false( size(funcParamsNames) );                                                 
end % if nargin>1

%% Attempt loading users parameters from input structures.
% In this case structure field name will be parameter name, and filed content will be
% parameter value.

iStructures=find( cellfun(@isstruct, unUsedVarargin) );
if ~isempty(iStructures)
    isSpecificParam=false(iStructures);
end
for iStruct=iStructures
    % analyze each structure unattained by previous "Load param name- param value pairs"
    % process
    CurrStruct=unUsedVarargin{iStruct};
    if numel(CurrStruct)>1
        CurrStruct=CurrStruct(1);
        warning('Structure arrays are not supported, 1''st element will be used.');
    elseif isempty(CurrStruct) % Ignore empty structures
        continue;
    end
    currFieldNames=fieldnames(CurrStruct);
    if isNoFuncParamsNames
        funcParamsNames=currFieldNames;
    end

    nFields=length(currFieldNames);
    for iFieldStr=1:nFields
        % Find relevant structure field names supported by function legal input names
        isFuncParamsNames=strcmpi(currFieldNames{iFieldStr}, funcParamsNames);
        if sum(isFuncParamsNames) > 1 % if several names were found try case sensitive match
            isFuncParamsNames=strcmp(currFieldNames{iFieldStr}, funcParamsNames);
        end
        
        % in case of ambiguty, use the first appearing name, as they are identical.
        iFirstFittingName=find(isFuncParamsNames, 1, 'first'); 
        
        if ~isempty(iFirstFittingName) % Load parameters into current environment
            assignin('caller',  funcParamsNames{iFirstFittingName},...
                CurrStruct.(currFieldNames{iFieldStr}) );
            isSpecificParam(iStruct)=true; % mark used input
            isUsedParamsName(iFirstFittingName)=true;
        end
    end % for iFieldStr=1:nFields
    if isNoFuncParamsNames
        funcParamsNames=[];
    else
        funcParamsNames=funcParamsNames(~isUsedParamsName); % do not allow repeating
                                                        % initialisation of variables
    end
end % for iStruct=find( cellfun(@isstruct, unUsedVarargin) )


if ~isempty(iStructures)
    % remove used input elements
    unUsedVarargin=unUsedVarargin( iStructures(~isSpecificParam) ); 
end
if isequal(unUsedVarargin, varargin) % neither inputs were used to extract user inputs
    % Preserve custom Matlab input parameters transfer scheme. Here inpus order defines
    % the variable destination.
    nInputs=min( numel(varargin), numel(funcParamsNames) );
    for iVargin=1:nInputs
        assignin( 'caller',  funcParamsNames{iVargin}, varargin{iVargin} );
    end
end