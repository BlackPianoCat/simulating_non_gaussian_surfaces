
##################################################################
#
# Coded in Python by SEBASTIAN KORSAK © 2021 (sebykorsak@gmail.com)
# Original file in MATLAB by Dr. V. Costantoudis
#
##################################################################

import numpy as np
import statistics as st
import cv2
import PIL
import matplotlib.pyplot as plt
from scipy.stats import johnsonsu

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from skimage import filters
from PIL import Image
from skimage.filters import threshold_otsu
from matplotlib.pyplot import figure

# Special thanks to Max Pierini who translated the MATLAB code to python for me
from j_johnson_M import f_johnson_M

def SAimage_fft_2(N_total=500,rms=3,skewness=1,kurtosis=3,corlength_x=10,corlength_y=10,alpha=0.9,th=1000.0,non_Gauss=True,corr=False):
	'''
	Input:
	
	Generation of non-Gaussian surfaces with predetermined 
	1. N_total: The dimension of our plane
	2. Standard deviation (rms)
	3. Skewness (skew)
	4. Kurtosis (kurt)
	5. Correlation lengths (ksix, ksiy)
	6. Roughness exponent (aplha)
	7. The threshold we choose to binarize our image
	8. non_Gaus: is True if we want the non-Gaussian simulation, and False when we want the Gaussian one
	9. corr: True if you want to illustrate correlation functions 

	The input values of kurtosis and skewness should satisfy the inequality
	kurtosis > 1 + skewness^2
	The method is described in Yang et al. CMES, vol.103, no.4, pp.251-279,2014
	We use pearson translation to transform random noise to non-gaussian series and 
	inverse Fourier transform to generate Gaussian surfaces.
	
	The method is implemented in three steps
	1. Generation of a Gaussian self-affine surface z_gs with rms, ksix,ksiy,alpha
	2. Generation of non-Gaussian noise z_ngn with mu,rms,skew,kurt
	3. Arranging z_ngn according to the spatial arrangement of z_g to get a non-Gaussian self-affine surface z_ngs
	'''
	
	##################################################################
	#
	# Coded in Python by SEBASTIAN KORSAK © 2021 (sebykorsak@gmail.com)
	# Original file in MATLAB by Dr. V. Costantoudis
	#
	##################################################################
	
	# Input
	N = N_total
	numpoints = N
	hhcf_y = np.zeros(int(N/2))
	Pyy = 0.0
	
	# input moments and spatial parameters
	skew = skewness;
	kurt = kurtosis
	ksix = corlength_x
	ksiy = corlength_y
	alpha = alpha
	mu = 0

	# 1st step: Generation of a Gaussian surface

	# Determine the autocorrelation function R(tx,ty) 
	R = np.zeros([numpoints+1,numpoints+1])
	txmin, txmax, tymin, tymax = -numpoints/2,numpoints/2,-numpoints/2,numpoints/2
	dtx, dty = (txmax-txmin)/(numpoints), (tymax-tymin)/(numpoints)
	
	tx = np.arange(txmin,txmax,dtx)
	ty = np.arange(tymin,tymax,dty)

	for txx in tx:
		for tyy in ty:
			R[int(txx+txmax+1),int(tyy+tymax+1)]=((rms**2)*np.exp(-(abs(np.sqrt((txx/ksix)**2+(tyy/ksiy)**2)))**(2*alpha)))
			
	# According to the Wiener-Khinchine theorem FR is the power spectrum of the desired profile
	FR = np.fft.fft2(R, s=[numpoints,numpoints])
	AMPR = np.sqrt(dtx**2+dty**2)*abs(FR)
	
	# 2nd step: Generate a white noise, normalize it and take its Fourier transform
	X = np.random.rand(numpoints,numpoints)
	aveX = X.mean(axis=0).mean(axis=0)
	dif2X = (X-aveX)**2;
	stdX = np.sqrt(dif2X.mean(axis=0).mean(axis=0));
	X = X/stdX;
	XF = np.fft.fft2(X, s=[numpoints,numpoints])
	
	# 3nd step: Multiply the two Fourier transforms
	YF = XF*np.sqrt(AMPR)
	
	# 4th step: Perform the inverse Fourier transform of YF and get the desired surface
	zaf = np.fft.ifft2(YF,s=[numpoints,numpoints])
	z = np.real(zaf)
	avez = z.mean(axis=0).mean(axis=0)
	dif2z = (z-avez)**2
	stdz = np.sqrt(dif2z.mean(axis=0).mean(axis=0))
	z = ((z-avez)*rms)/stdz

	# Define the fraction of the surface to be analysed
	xmin, xmax, ymin, ymax= 0, N, 0, N
	z_gs=z[xmin:xmax,ymin:ymax]
	print(z_gs.shape)
	Nh=xmax-xmin+1
	
	# For the Gaussian Surface Simulation
	if not non_Gauss:
		X = np.arange(0,N)
		Y = np.arange(0,N)
		X, Y = np.meshgrid(X, Y)

		fig = plt.figure()
		fig.set_figwidth(10)
		fig.set_figheight(10)
		ax = plt.axes(projection='3d')
		ax.plot_surface(X,Y,z_gs,cmap='viridis', edgecolor='none')
		ax.view_init(30, 150)
		plt.title('Gaussian Surface',fontsize=16)
		plt.show()

		plt.imshow(z_gs<th, cmap='gray', origin='lower', interpolation='none')
		plt.title('Binarized Image from Gaussian Surface')
		plt.show()

	
	del R, AMPR, z

	# For the Non-Gaussian Surface Simulator
	if non_Gauss:
		#2nd step: Generation of a non-Gaussian noise NxN
		## Finding the parameters of JohnsonSU distribution
		## Special thank to Max Pierini for this part of the code
		coef, j_type, err = f_johnson_M(mu, rms, skew, kurt)
		
		gamma, delta, xi, lam = coef
		
		## Simulating from JohnsonSU distribution
		z_ngn = johnsonsu.rvs(a=gamma, b=delta, loc=xi, scale=lam,size=[N,N])
		
		#3rd step: Combination of z_gs with z_ngn to output a z_ngs
		v_gs = z_gs.reshape(-1)
		v_ngn = z_ngn.reshape(-1)
		Igs = np.argsort(v_gs)
		Ingn = np.argsort(v_ngn)
		vs_gs = np.sort(v_gs)
		vs_ngn = np.sort(v_ngn)
		
		v_ngs = np.zeros(np.max(Igs)+1)
		
		for iv in range(N*N):
			ivs = Igs[iv]
			v_ngs[ivs] = vs_ngn[iv]
		
		X = np.arange(0,N)
		Y = np.arange(0,N)
		X, Y = np.meshgrid(X, Y)
		z_ngs = -v_ngs.reshape(N,N)
		# Creating color map
		my_cmap = plt.get_cmap('hot')

		fig = plt.figure()
		fig.set_figwidth(10)
		fig.set_figheight(10)
		ax = plt.axes(projection='3d')
		ax.set_xlim(0,N)
		ax.set_ylim(0,N)
		ax.set_zlim(-1000,5000)
		ax.plot_surface(X,Y,z_ngs,cmap=my_cmap, edgecolor='black')
		ax.view_init(30, 150)
		plt.title('Non-Gaussian Surface',fontsize=16)
		# Add a color bar which maps values to colors.
		plt.show()

		plt.imshow(z_ngs>th, cmap='gray', origin='lower', interpolation='none')
		plt.title('Binarized Image from Non-Gaussian Surface',fontsize=16)
		plt.show()


		if corr:
			# Correlation Functions
			# a. 1-D height-height correlation function
			inpsur=z_ngs
			hhcf1d=np.zeros(N//2)
			rdif=np.zeros(N//2)
			for ndif in range(N//2):
				surf1=inpsur[0:N,0:(N-ndif)]
				surf2=inpsur[0:N,ndif:N]
				difsur2=(surf1-surf2)**2;
				hhcf1d[ndif]=np.sqrt(difsur2.mean(axis=0).mean(axis=0));
				rdif[ndif]=ndif
				

			plt.loglog(rdif,hhcf1d)
			plt.grid()
			plt.xlabel('log(r(nm))')
			plt.ylabel('log(G(r) (nm))')
			plt.title('1-D height-height correlation function (non-Gaussian surface)')
			plt.show()