function outArr=uniqueUnsortedVec(inArr)
%% uniqueUnsortedVec
% The function returns the same data as in input, but with no repetitions and in orignal
%   order.
%
%% Syntax
%  outArr=uniqueUnsortedVec(inArr);
%
%% Description
% The function is just a layer over the excellent and well known Matlab 'unique' function.
%  The unique function sotrs the returned results- whihc is not always desirable. So this
%  function does the same, but the elemenst order doe not change. This function is highly
%  similar to 'unique', but it operates on arrays only (numerical, character and cell
%  arays). It returns only the resulting array.
%
%% Input arguments (defaults exist):
% inArr- a 1D vector of: numericals, logicals, characters, categoricals, or a cell array
%   of strings.
%
%% Output arguments
% outArr- a 1D vector composed of unique elements of inArr, without cvhnaging their
%   interal order.
%
%% Example
% inArr={'a', 'c', 'a', 'a', 'b', 'c'}
% unique(inArr)
% uniqueUnsortedVec(inArr)
%
[~, iOrigA, ~]=unique(inArr);
if length(iOrigA) == length(inArr)  % unique found no repetitions- no need to do anything
    outArr=inArr;
else
    outArr=inArr( sort(iOrigA) );   % remove repeating elements, keeping same ordering
end