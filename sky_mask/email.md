Hello Etienne,

A super quick one, I just got the results of the pca2d calculated on the t.fits masked with 3sigma sky and 50% recon transmission.
I used --n-star 0 --n-earth 10 

We improved, the 6d signal now is clear. We have no evidence for the 3 nor 18d signals. 
However, the pca2d lbl is not corrected by the fp drift nor zero point

When i compare the pca2d with NO drift corrected lbl2,
an mcmc on the data show the same planet with the same amplitude for both methods 

NO drift masked tcorr: P=6.02s K=5.9m/s  rms series: 8.5 m/s
pca2d masked tfits: P=6.02s K=5.9 m/s rms series: 12.4 m/s

for reference:
drift masked tcorr: P=6.02d K=5.7 m/s   rms series: 7.5 m/s

I dont why the RV rms gets much higher on the pca2d rdb data. Ideas?

I will run tomorrow  GL251 and I will send you the code to mask the t.fits. Here the full pca2d report.
https://filesender.renater.fr/?s=download&token=e2ff86cc-dd8d-4d90-b1f9-87c0ebb0baac


Cheers!

Andres


PS. Periodograms pca2d from t. fits masked 

PastedGraphic-6.png


PS2. periodograms NO drift lbl of the tcorr masked data.

PastedGraphic-5.png




***************************************
Andres CARMONA GONZALEZ

Institut de Recherche en Astrophysique et Planétologie (IRAP)
Observatoire de Midi-Pyrénées
14, avenue Edouard Belin
31400 Toulouse
FRANCE

andres.carmona-gonzalez@univ-tlse3.fr
andres.carmona.astronomy@gmail.com

Project SPIP: https://spip.omp.eu/fr/exoplanetes
Project SPIRou: https://spirou.omp.eu/science/

***************************************