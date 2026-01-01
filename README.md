# interpolation
Gaussian interpolation functions for 2D and 3D datasets.

constructs an interpolated value matrix at grid positions using original values at unevenly spaced locations; the parameter vector is interpolated using gaussian weighting of distances (of points from the corresponding bin); in addition it is possible to weight contributions using an external weight vector; in order to avoid over-interpretation of sparsely distributed data, it is also possible to constrain the interpolation to the bins points that are within the specified threshold distance of existing points
