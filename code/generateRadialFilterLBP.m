function [radInterpFilt]=generateRadialFilterLBP(p, r)
%% generateRadialFilterLBP
% The function returns a filter with indexes arranged clock-wise in a circular shape.

%% Description
% The function is aimed to return a radial filter for LMP generation. 
%% Input arguments (defaults exist):
% inArr- a 1D vector of: numericals, logicals, characters, categoricals, or a cell array
%   of strings.
%
%% Output arguments
% p- an integer number specifying the number of neigbours- number of enabled filter
%   elemnts.
% r- a positive number specifying the raduis of the filter. Can be a non integer

%% Default params
if nargin<2
    r=1;
    if nargin<1
        p=8;
    end
end

%% verify params leget values
r=max(1, r);    % radius below 1 is illegal
p=round(p);     % non integer number of neighbours sound oucward
p=max(1, p);    % number of neighbours below 1 is illegal


theta=linspace(0, 2*pi, p+1)+pi/2;   
theta=theta(1:end-1);           % remove obsolite last element (0=2*pi)

%% Find relevant coordinates
[rowsFilt, colsFilt] = pol2cart(theta, repmat(r, size(theta) )); % convert to cartesian
nEps=-3;
rowsFilt=roundnS(rowsFilt, nEps);
colsFilt=roundnS(colsFilt, nEps);

% Matrix indexes should be integers
rowsFloor=floor(rowsFilt);
rowsCeil=ceil(rowsFilt);

colsFloor=floor(colsFilt);
colsCeil=ceil(colsFilt);

rowsDistFloor=1-abs( rowsFloor-rowsFilt );
rowsDistCeil=1-abs( rowsCeil-rowsFilt );
colsDistFloor=1-abs( colsFloor-colsFilt );
colsDistCeil=1-abs( colsCeil-colsFilt );

% Find minimal filter dimentions, based on indexes
filtDims=[ceil( max(rowsFilt) )-floor( min(rowsFilt) ),...
    ceil( max(colsFilt) )-floor( min(colsFilt) ) ];
filtDims=filtDims+mod(filtDims+1, 2); % verify filter dimentions are odd

filtCenter=(filtDims+1)/2;

%% Convert cartesian coordinates to matrix elements coordinates via simple shift
rowsFloor=rowsFloor+filtCenter(1);
rowsCeil=rowsCeil+filtCenter(1);
colsFloor=colsFloor+filtCenter(2);
colsCeil=colsCeil+filtCenter(2);


%% Generate the filter- each 2D slice for filter element  
radInterpFilt=zeros( [filtDims,  p], 'single'); % initate filter with zeros
for iP=1:p
    radInterpFilt( rowsFloor(iP), colsFloor(iP), iP )=...
        radInterpFilt( rowsFloor(iP), colsFloor(iP), iP )+rowsDistFloor(iP)+colsDistFloor(iP);
    
    radInterpFilt( rowsFloor(iP), colsCeil(iP), iP )=...
        radInterpFilt( rowsFloor(iP), colsCeil(iP), iP )+rowsDistFloor(iP)+colsDistCeil(iP);
    
    radInterpFilt( rowsCeil(iP), colsFloor(iP), iP )=...
        radInterpFilt( rowsCeil(iP), colsFloor(iP), iP )+rowsDistCeil(iP)+colsDistFloor(iP);
   
    radInterpFilt( rowsCeil(iP), colsCeil(iP), iP )=...
        radInterpFilt( rowsCeil(iP), colsCeil(iP), iP )+rowsDistCeil(iP)+colsDistCeil(iP);
    
    radInterpFilt( :, :, iP )=radInterpFilt( :, :, iP )/sum(sum(radInterpFilt( :, :, iP )));
end
% imshow(sum(radInterpFilt,3), []);

% Substract 1 at central element to get difference between central element and relevant
% neighbours: (5) T=p{s(g1-g0), s(g2-g0),...,s(gn-g0)}
radInterpFilt( filtCenter(1), filtCenter(2), : )=...
    radInterpFilt( filtCenter(1), filtCenter(2), : )-1; 


