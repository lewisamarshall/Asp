% Buffer Capacity Example

TRIS=ion('Tris', 1, 8.076, 29.5e-9);
CL=ion('Chloride', -1, -2, -79.1e-9);

BUFFER=solution({TRIS, CL}, [0.1, .05])
conductivity=BUFFER.conductivity
buffering_capacity=BUFFER.buffering_capacity
transference_numbers=BUFFER.transference

%From Peakmaster, w/o Ionic Strength Correction
% Ionic strength = 50.001 mM
% Conductivity = 0.524 S/m
%Buffer Capacity = 57.567 mM