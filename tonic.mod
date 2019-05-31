NEURON{
SUFFIX tonic
NONSPECIFIC_CURRENT igaba
RANGE  igaba, e, g
}
PARAMETER {
g=0.001 (siemens/cm2) <0, 1e9>
e = -70 (millivolt)
}
ASSIGNED{
igaba (milliamp/cm2)
v (millivolt)
}
BREAKPOINT {igaba=g*(v-e)}
