function [ result ] = interpolatePointVals2Vec_gauss( valuematrix, bins, sigma, displayVariable )
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
%           bins            vector of the desired value bins of the
%                           x-vector to which the value will be
%                           interpolated
%
%           sigma           sigma of the Gaussian weighting function (pix)
%                           smaller sigmas will of course allow you to 
%                           resolve smaller spatial variations, but it
%                           makes no sense to choose a sigma that is much
%                           below the typical distance between individual
%                           features
%
% OUPUT:   result          n x 2 matrix containing:
%                           first column:   bin x-coordinates
%                           second column:  interpolated y-coordinates
%
% written 02/06/2014 Dinah Loerke

display = 0;
if nargin>3
    if displayVariable==1
        display = 1;
    end
end

%   extract x,y coordinate and parameter vectors from input
xvec = valuematrix(:,1);
yvec = valuematrix(:,2);

weightvec_sig = 1+0*xvec;
if size(valuematrix,2)>2
    weightvec_sig = valuematrix(:,3);
end


% loop over bins
for i=1:length(bins)
    
    % display progress
    % display(['processing mask section ',num2str(i),' of ', num2str(length(stepvector))]);
      
    % distance vector (from current bin) for each data point
    distvec = abs(xvec-bins(i));
    
    % weighting vector for each data point (Gaussian function with distance
    % from bin)
    weightvec_dist = exp(-(distvec.^2)/(2*sigma^2));
    
    nanpos = isnan(yvec);
    weightvec_dist(nanpos) = nan;
    
    % multiply this weight vector with the 'external' weights, which may
    % represent e.g. statistical significance of the point
    weightvec = weightvec_dist.*weightvec_sig;
    
    % weighted average of the data parameter (yvec) for this bin
    % valuvec = nansum((yvec.*weightvec),1)./nansum(weightvec,1); % OLD CODE -- nansum no longer supported
    valuvec = sum((yvec.*weightvec),1,'omitnan')./sum(weightvec,1,'omitnan'); % added 6/13/22 KB


    % write results (averaged parameter values for the current pixels) into
    % the results at the specified positions 
    result(i,1:2) = [bins(i) valuvec];
    
end % of for i-loop

% display results
if display == 1
    figure;
    hold off; plot(xvec,yvec,'c.');
    hold on; plot(result(:,1),result(:,2),'r.-');
end

end

