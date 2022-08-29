function [ result ] = interpolatePointVals2Mat_gauss2018( valuematrix, binsX, binsY, sigmaX, sigmaY, threshold, displayVariable )
% construct an interpolated value matrix at grid positions using original
% values at unevenly spaced locations; the parameter vector is
% interpolated using gaussian weighting of distances (of points from the
% corresponding bin); in addition it is possible to weight contributions
% using an external weight vector; in order to avoid over-interpretation of
% sparsely distributed data, it is also possible to constrain the
% interpolation to the bins points that are within the specified threshold 
% distance of existing points
%
% INPUT:    
% valuematrix   n x 3 matrix containing:
%               first column:       features' x-coordinates (e.g. distance)
%               second column:      features' y-coordinates (e.g. angle)
%               third column:       features' value (e.g. correlation magnitude)
%               fourth column:      (optional) weight            
%
% binsX         vector of the desired value bins of the x-vector to which
%               the value will be interpolated
% binsY         vector of the desired value bins of the y-vector to which
%               the value will be interpolated
%
% sigmaX        sigma of the Gaussian weighting function in x
% sigmaY        sigma of the Gaussian weighting function in y
%
% threshold (OPTIONAL) 1x2 vector of threshold distances; the threshold
%               distances respresent what distance ((in units of multiples
%               of sigma in x and y-direction)) from real data points is
%               still allowable for interpolation (so e.g. threshold = [3 3] 
%               means that only those bins within 3*sigma distance of real
%               data are calculated), so that bin points too far from real
%               data points could be assigned to nan
%
% displayVariable (OPTINAL) set to 1 if plotting of results is desired
%                           
% OUPUT:   
% result        n x m  x 3 matrix containing:
%               rows:       interpolated values for x bins (for given y)
%               columns:    interpolated values for y bins (for given x)
%               layer 1:    interpolated values at these positions (from
%                           third row)
%               layers 2:   standard deviation associated with average
%               layer 3:    true number of points within rectangular area
%               layer 4:    local point density for interpolation
%               layer 5:    x-bin values
%               layer 6:    y-bin values
%
% last modified 10/07/2014 Dinah Loerke
% last modified 06/01/2018 Dinah Loerke

display = 0;
if nargin>6
    if displayVariable==1
        display = 1;
    end
end

% determine if threhsolding is desired and set thresholds
thresholdVal = 0;
if nargin>5
    if ~isempty(threshold)
        thresholdVal = 1;
        if length(thresholdVal)>1
            thrX = threshold(1);
            thrY = threshold(2);
        else
            thrX = threshold(1);
            thrY = threshold(1);
        end        
    end
end
%   extract x,y coordinate and parameter vectors from input
xvec = valuematrix(:,1);
yvec = valuematrix(:,2);
zvec = valuematrix(:,3);

% default weight vector is equal to 1; if an external weight vector exists
% (fourth column), then that weight vector is used
weightvec_sig = 1+0*xvec;
if size(valuematrix,2)>3
    weightvec_sig = valuematrix(:,4);
end

% make binWidth vectors
binWidthX = diff(binsX); binWidthX(length(binWidthX)+1)=binWidthX(length(binWidthX));
binWidthY = diff(binsY); binWidthY(length(binWidthY)+1)=binWidthY(length(binWidthY));

% loop over x-bins
for i=1:length(binsX)
    
    % display progress
    % display(['processing mask section ',num2str(i),' of ', num2str(length(stepvector))]);
    
    % distance vector 1 (from current x bin) for each data point
    distvec1 = abs(xvec-binsX(i));
    
    % internal weighting vector is Gaussian function with distance of data
    % points from bin center, meaning that data points farther away from
    % the bin center contribute progressively less
    weightvec_dist1 = exp(-(distvec1.^2)/(2*sigmaX^2));
    
    % true number of points within bin (for output): distance is within bin/2
    fd1 = find(distvec1<binWidthX(i)/2);
    
    % threshold if desired
    if thresholdVal==1
        findThrX = find( distvec1 > (thrX*sigmaX) );
        weightvec_dist1(findThrX) = nan;
    end
    
   % loop over y-bins
    for k=1:length(binsY)
        
        % distance vector (from current y-bin) for each data point
        distvec2            = abs(yvec-binsY(k));

        % weighting vector for each data point is Gaussian function with 
        % distance from bin (multiplying factors for both bins)
        weightvec_dist2     = exp(-(distvec2.^2)/(2*sigmaY^2));
        
        % threshold if desired
        if thresholdVal==1
            findThrY        = find( distvec2 > (thrY*sigmaY) );
            weightvec_dist2(findThrY) = nan;
        end
        
        % combined internal weight vecyor
        weightvec_distTotal = weightvec_dist1 .* weightvec_dist2;
        
        % total weight vector (including external if desired)
        nanpos              = isnan(zvec);
        weightvec_distTotal(nanpos) = nan;

        % multiply this weight vector with the 'external' weights, which may
        % represent e.g. statistical significance of the point
        weightvec           = weightvec_distTotal.*weightvec_sig;

        % weighted average of the data parameter (yvec) for this bin
        valuvec             = nansum((zvec.*weightvec),1)./nansum(weightvec,1);
        
        % calculate standard deviation as the (weighted) root mean variance
        variancevec         = ((zvec - valuvec).^2);
        varianceave         = nansum((variancevec.*weightvec),1)./nansum(weightvec,1);
        standard_dev        = sqrt(varianceave);

        % true number of points within bin: distance is within bin/s
        fd2                 = find(distvec2<binWidthY(k)/2);
        % combine fd1 (from i-loop) with fd2 to find number of points in
        % current bin
        fd                  = intersect(fd1, fd2);
        true_number         = length(fd);
        
        % write results (averaged parameter values for the current pixels)
        % into the results at the specified positions 
        % first layer: interpolated weighted intensities
        result(i,k,1) = valuvec;
        % second layer: standard deviation
        result(i,k,2) = standard_dev;
        % third layer: true number of points fully within each bin
        result(i,k,3) = true_number;
        % fourth layer: local point density from weight vector
        result(i,k,4) = nansum(weightvec,1);
        % fifth layer: x bin values
        result(i,k,5) = binsX(i);
        % sixth layer: y bin values
        result(i,k,6) = binsY(k);
              
    end % of for k-loop
    
end % of for i-loop

% display results
if display == 1
    figure; 
    subplot(1,2,1);
    plot(xvec,zvec,'b.');
    subplot(1,2,2);
    plot(yvec,zvec,'r.');
    
    figure;
    subplot(2,2,1); imshow(result(:,:,1),[]); 
                    colormap jet; title('value');
                    xlabel('channel 2'); ylabel('channel 1');
    subplot(2,2,2); imshow(result(:,:,2),[]); 
                    colormap jet; title('std');
    subplot(2,2,3); imshow(result(:,:,3),[]); 
                    colormap jet; title('true #');
    subplot(2,2,4); imshow(result(:,:,4),[]); 
                    colormap jet; title('density');
    % set(gca,'Position',[0 0 1 1]);
end

end

