% My idea:
% 1. just give lower weights/confidence to data from neighboring vehicles 
% outside of a "better commutable" range
% 2. give lower weights to data from neighboring vehicles within a certain
% range of road degree.
% 3. set a block probability: randomly assign [0,1] to each vehicle, block
% those with label < p, i.e., p = 0.2
