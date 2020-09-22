from astropy.cosmology import FlatLambdaCDM

def timeEstimate(delayTime): 
	cosmo = FlatLambdaCDM(H0=67.90, Om0=0.3065)
	return cosmo.lookback_time(delayTime)
	
print(timeEstimate(0.41))
