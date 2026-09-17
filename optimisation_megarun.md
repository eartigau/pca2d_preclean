NIRPS data per objet is here :

/cosmos99/nirps/apero-data/nirps_he_online/objects

SPIRou objects are here:

/cosmos99/spirou/apero-data/spirou_offline/objects

We want to make per-instrument/per-objects folder in /home/artigau/pca2d_optimisation and do symlinks to the relevant t.fits files

We want to understand and optimise the code. The ultimate optimisation metric is the reduction of the robust-sigma of the RV sequences after the LBL analysis. Here are things that need to be played-with and understood how the affect the ultimate performance of the method:

correction fit metric : F or dF/dv^2 ??
how many star components : 0 (mean stellar flux), 1 (single stellar variablity accepted), more?
how many observer components : 1, 3, 7, more?
what is the best high-pass scale length : 100km/s? That's the current best guess, how about 50 or 200 km/s?

fit a velocity per exposure : yes or no? still to be understood

sets you will want to try : 

GL406 alone, with NIRPS and with SPIRou, one instrument at a time (never
mixed in one fit: two instruments are two runs)
TOI4552 TOIM4508 TOI782 with NIRPS together vs alone
TOI2120 TOI6091 with SPIRou together vs alone
GL725B, GL251, GL48 with SPIRou together vs alone (in particular GL725B that has little BERV coverage)
GJ1 and GJ707 with NIRPS, alone/together

All tested scenario should generate a full report

At the end, present a full scientific paper, with references. Use the consensus app to fully get info that others worked on, such as the yarara and references therein.

Prepare figures and tables to uphold your findings and conclusions.