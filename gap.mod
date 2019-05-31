NEURON {
	POINT_PROCESS Gap
	NONSPECIFIC_CURRENT i_gap
	RANGE r, i_gap
	POINTER vgap
}

PARAMETER {
	v (millivolt)
	vgap (millivolt)
      r =100000 (megohm) 
}

ASSIGNED {
	i_gap (nanoamp)
}

BREAKPOINT {
	i_gap = (v - vgap)/r
}

