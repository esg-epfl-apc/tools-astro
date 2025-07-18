Axis priors have different format according to whether the model parameter is a numerical or categorical variable.

- SED and Reddening Curve Axes Priors: here the first column is a string, representing the name of the SED template or the Reddening Curve accordingly. The second column is a double precision decimal number, representing the prior weight to multiply the likelihood with.

- E(B-V) and Redshift Axes Priors: the first column contains E(B-V) or redshift values accordingly, and the second column contains the prior weight to multiply the likelihood with.
