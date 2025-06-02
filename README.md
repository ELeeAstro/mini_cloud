# mini_cloud

`mini-cloud' is a selection of time-dependent cloud microphysical models primarily using the mass moment method for exoplanet atmosphere cloud formation.



## using the code

Each source code directory contains a makefile to build the model. Each model is slightly different, but generally src_mini_cloud_2_mono_mix is the current recommended version for most exoplanet mineral cloud schemes.

This model is being developed and is dynamically changing, expect frequent bug fixes, improvements and updates.

### src_mini_cloud_sat_adj

Single moment tracer saturation adjustment scheme

### src_mini_cloud_2_mono_mix

-This is the main `stable' code used for coupling to dynamics-

Two-moment monodisperse size distribution code 

### src_mini_cloud_2_mono

Research/development code for two-moment monodisperse size distribution.
This code is generally used for developing theory/modelling using the monodisperse moment method and not much else.


### src_mini_cloud_2_exp

Research/development code for two-moment exponential size distribution.
This code is generally used for developing theory/modelling using the exponential moment method and not much else.

### src_mini_cloud_3_gamma

Research/development code for three-moment gamma size distribution.
This code is generally used for developing theory/modelling using the gamma moment method and not much else.

### src_mini_cloud_3_lognormal

Research/development code for three-moment log-normal size distribution.
This code is generally used for developing theory/modelling using the log-normal moment method and not much else.

### src_mini_cloud_ori

Legacy original mini-cloud source code.


