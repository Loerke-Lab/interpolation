function [ result ] = interpolatePointVals2Vec_bin( valuematrix, bincenters, binwidth, displayVariable )
% construct a value vector using values at unevenly spaced locations; the 
% parameter vector is interpolated using gaussian weighting of distances; 
% in addition it is possible to weight contributions using an external
% weight vector
%
% INPUT:    valuematrix     n x 2 matrix containing:
%                           first column:   features' x-coordinates
%                           second column:  features' y-coordinates
%                           third column: (optional) weight
%                           
%
%           bincenters      vector of the desired value bin centers of the
%                           x-vector to which the value will be
%                           interpolated
%
%           binwidth        total width of the bins (width/2 on the left,
%                           and width/2 on the right of the bin center)
%                           if no value is specified, the bin width will be
%                           equivalent to the distance between bin centers;
%                           if width is =0, then only the values with an
%                           exact match to the bin centers will be
%                           considered
%
% OUPUT:   result          n x 2 matrix containing:
%                           first column:   bin x-coordinates
%                           second column:  interpolated y-coordinates
%
% written 01/09/2017 Dinah Loerke

display = 0;
if nargin>3
    if displayVariable==1
        display = 1;
    end
end

bindiff = nanmean(diff(bincenters));
if nargin>2
    if ~isempty(binwidth)
        bindiff = binwidth;
        binedges_left = bincenters-(bindiff/2);
        binedges_right = bincenters+(bindiff/2);
    end
end

%   extract x,y coordinate and parameter vectors from input
xvec = valuematrix(:,1);
yvec = valuematrix(:,2);


% loop over bins
for i=1:length(bincenters)
    
    % display progress
    % display(['processing mask section ',num2str(i),' of ', num2str(length(stepvector))]);
      
    % distance vector (from current bin) for each data point
    fpos = find( (xvec==bincenters(i)) | ( (xvec>=binedges_left(i)) & (xvec<binedges_right(i) ) ) ) ;
    
    % weighted average of the data parameter (yvec) for this bin
    valuvec = nanmean(yvec(fpos));
    
    % write results (averaged parameter values for the current pixels) into
    % the results at the specified positions 
    result(i,1:2) = [bincenters(i) valuvec];
    
end % of for i-loop

% display results
if display == 1
    figure;
    hold off; plot(xvec,yvec,'c.');
    hold on; plot(result(:,1),result(:,2),'r.-');
end

end

