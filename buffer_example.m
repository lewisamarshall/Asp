% Buffer Capacity Example

TRIS=ion('Tris', 1, 8.076, 29.5e-9);
CL=ion('Chloride', -1, -2, -79.1e-9);

BUFFER=solution({TRIS, CL}, [.1, .05]);
BUFFER.buffering_capacity
BUFFER.ionic_strength
BUFFER.cond

%From Peakmaster, w/o Ionic Strength Correction
% Ionic strength = 50.001 mM
% Conductivity = 0.524 S/m
%Buffer Capacity = 57.567 mM