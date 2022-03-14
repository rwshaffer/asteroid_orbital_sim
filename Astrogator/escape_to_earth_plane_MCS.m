function [] = escape_to_earth_plane_MCS(root,sat)

ASTG = sat.Propagator;

MCS = ASTG.MainSequence;
MCS.RemoveAll;



end