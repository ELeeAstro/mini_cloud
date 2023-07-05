# mini_cloud

A tidy up and more in-depth guide is in progress. Additional species of interest are also being added.


NOTE: McCormack solver will be changed to pressure vertical coordinates

## using the code

src_mini_cloud is the main version using dvode as the ODE solver, *_Hairer, *_L and *_limex are test codes and should not be used.

nk_tables contains the optical constant database used in DIHRT and mini-cloud.

main.f90 shows an example interface to mini-cloud and how it should appear in the GCM

You can compile the source code in src_mini_cloud and entering make.
Then run the example code in the directory above ./mini_cloud

Outputs are in tracers.txt and opac.txt, results of the test can be plotted using the plot_examples.py and plot_examples_opac.py codes.

The old old version tried to use the budaj tables to interpolate the opacity, kept here for posterity, however, it was found using the Mie theory was more accurate and faster in the end.

NOTE: The test is deliberately designed to push the integrator to the limits with very strong evaporation and condensation rates with rapidly changing temperatures. This typically does not occur in the GCM models (that are reasonably thermally stabilized).

