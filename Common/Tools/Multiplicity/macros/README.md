# Example multiplicity calibration macros

You will find some example macros in this directory that will allow for the calculation of a glauber fit and the corresponding calibration histograms.

A simplified description of the procedure necessary to generate a
Glauber + NBD fit to a certain histogram is described below.

## Performing a Glauber + NBD fit

### First step: calculation of Glauber MC sample

The Glauber + NBD model assumes that the multiplicity / amplitude
distribution that one intends to fit can be described as coming
from a certain number N_{ancestor} of particle-emitting sources called 'ancestors'. Each ancestor is assumed to produce particles
according to a negative binominal distribution. Further,
the number of ancestors is assumed to be related to the basic
glauber quantities N_{part} and N_{coll} via the formula:

N_{ancestor} = f * N_{part}  + (1-f) N_{coll}

Usually, the value f is fixed at 0.8, and thus the range of
N_{ancestors} is typically 0-700 in Pb-Pb collisions.

In order to allow for Glauber + NBD fitting, the correlation of
(N_{part}, N_{coll}) needs to be known. For that purpose,
a tree of Glauber MC needs to be generated using TGlauberMC using the relevant nuclei, which also involves choosing an appropriate nuclear profile.

Once TGlauberMC has been used to produce a tree with N_{part} and
N_{coll} values, their correlation needs to be saved to a 2D histogram that serves as input to the Glauber + NBD fitter used in
O2/O2Physics. This is done using the macro called `saveCorrelation.C` contained in this directory, which produces a file named `baseHistos.root`. The file `saveCorrelation.C` serves as
an example for the Pb-Pb case and minor adaptation may be needed
in case other nuclei are to be used.

### Second step: execution of Glauber + NBD fitter

The fitting procedure is done via the macro ``. Because
the numerical fitting utilizes a convolution of a N_{ancestor}
distribution and NBD distributions, it will not be as fast
as typical one-dimensional fits: typical times of 10-100 seconds
are not unusual. The macro will produce an output file
that contains:

* the original fitted histogram
* the actual glauber fit distribution, plotted all the way
to zero multiplicity

This output can then be used to extract percentiles.

### Third step: calculation of percentile boundaries

Once both the data and glauber fit distributions are known,
the next step is to create a 'stitched' data/glauber distribution in
which the distribution at low multiplicities follows
the glauber fit and the distribution at higher multiplicities
follows the data. The multiplicity value in which the switch from
glauber to data takes place is called 'anchor point'.

Because of the fact that this 'stitching' procedure may need to be tuned and the actual glauber fit is slow, the stitching is done
in a separate macro called `runCalibration.C`. It is provided
in this directory as well and it is the third and last step in calculating percentile boundaries. The macro will printout some
key boundaries as well as save an output file with the calibration.

*Bonus*: at the stage in which the glauber fit has been done
and all information is available, a de-convolution process
can be used to calculate the average N_{part} and N_{coll}
in centrality slices. This functionality is also provided
in O2Physics as part of the `multGlauberNBDFitter` and the
`runCalibration.C` macro can optionally also perform that
deconvolution. *Warning*: this procedure might take a mooment.