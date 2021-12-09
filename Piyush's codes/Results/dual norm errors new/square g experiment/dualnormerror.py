import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("/scratch/users/ppanchal/mpltools/mpltools")
from mpltools import annotation

base_fname1 = "mp5polyxymn1_"
base_fname2 = "mp5polyxymn2_"

gramat = np.loadtxt('sqsq_poly_gramat.txt')
Ginv = np.linalg.inv(gramat)

# length of the error matrices
size = 8

cvg_val = 8

# Matrix containing both errors
errors = np.empty([size,2])
h = np.empty([])
for k in range(size):
	# boundary formula
	# creating the f vector
	fbem = np.empty([50])
	fbnd = np.empty([50])
	for i in range(5):
		for j in range(5):
			fname1 = base_fname1 + str(i+1) + "_" + str(j+1)
			converged_fname1 = "converged/"+fname1
			temp1 = np.loadtxt(fname1)
			temp1_converged = np.loadtxt(converged_fname1)
			idx_ij = 5*i + j
			fbnd[idx_ij] = temp1[k,4]-temp1[cvg_val,2]
			fbem[idx_ij] = temp1[k,2]-temp1[cvg_val,2]
			
			
			fname2 = base_fname2 + str(i+1) + "_" + str(j+1)
			converged_fname2 = "converged/"+fname2
			temp2 = np.loadtxt(fname2)
			temp2_converged = np.loadtxt(converged_fname2)
			fbnd[25+idx_ij] = temp2[k,4]-temp2[cvg_val,2]
			fbem[25+idx_ij] = temp2[k,2]-temp2[cvg_val,2]
			if k == 0:
				h = 1/temp1[:size,0]
	#print("fbem", fbem)
	#print("fbnd", fbnd)
	errors[k,0] = np.sqrt(np.dot(fbem, np.matmul(Ginv, fbem))) 
	errors[k,1] = np.sqrt(np.dot(fbnd, np.matmul(Ginv, fbnd))) 
	#print(errors[k,:])

fig = plt.figure()
ax = fig.add_subplot(111)
plt.loglog(h,abs(errors[:,0]),'.-',color="red",zorder=2)
plt.loglog(h,abs(errors[:,1]),'.-',color="blue",zorder=2)
#plt.legend(["BEM shape derivative","Boundary shape derivative"])
plt.xlabel("meshwidth h")
plt.ylabel("Dual norm error")

# Adding best fit lines
fit_bnd = np.polyfit(np.log(h)[:],np.log(errors[:,1]),1)
fit_bem = np.polyfit(np.log(h)[:],np.log(errors[:,0]),1)
y_bnd = fit_bnd[0] * np.log(h) + fit_bnd[1]
y_bem = fit_bem[0] * np.log(h) + fit_bem[1]

h1 = h
y_bnd1 = y_bnd
y_bem1 = y_bem

bestfitcolor = "silver"

print("poly Bnd rate:",fit_bnd[0])
print("poly BEM rate:",fit_bem[0])

# Doing the same thing for sin fields
base_fname1 = "scg1_"
base_fname2 = "scg2_"

#gramat = np.loadtxt('sqsq_POLY_gramat.txt')
Ginv = np.linalg.inv(gramat)
# Matrix containing both errors
errors = np.empty([size,2])
h = np.empty([])
for k in range(size):
	# boundary formula
	# creating the f vector
	fbem = np.empty([50])
	fbnd = np.empty([50])
	for i in range(5):
		for j in range(5):
			fname1 = base_fname1 + str(i+1) + "_" + str(j+1)
			converged_fname1 = "converged/"+fname1
			temp1 = np.loadtxt(fname1)
			temp1_converged = np.loadtxt(converged_fname1)
			idx_ij = 5*i + j
			fbnd[idx_ij] = temp1[k,4]-temp1[cvg_val,2]
			fbem[idx_ij] = temp1[k,2]-temp1[cvg_val,2]
			#fbnd[idx_ij] = temp1[k,4]-temp1[-1,2]
			#fbem[idx_ij] = temp1[k,2]-temp1[-1,2]
			
			
			fname2 = base_fname2 + str(i+1) + "_" + str(j+1)
			converged_fname2 = "converged/"+fname2
			temp2 = np.loadtxt(fname2)
			temp2_converged = np.loadtxt(converged_fname2)
			fbnd[25+idx_ij] = temp2[k,4]-temp2[cvg_val,2]
			fbem[25+idx_ij] = temp2[k,2]-temp2[cvg_val,2]
			#fbnd[25+idx_ij] = temp2[k,4]-temp2[-1,2]
			#fbem[25+idx_ij] = temp2[k,2]-temp2[-1,2]
			if k == 0:
				h = 1/temp1[:size,0]
	#print("fbem", fbem)
	#print("fbnd", fbnd)
	errors[k,0] = np.sqrt(np.dot(fbem, np.matmul(Ginv, fbem))) 
	errors[k,1] = np.sqrt(np.dot(fbnd, np.matmul(Ginv, fbnd))) 
	#print(errors[k,:])

plt.loglog(h,abs(errors[:,0]),'.-',color="green",zorder=2)
plt.loglog(h,abs(errors[:,1]),'.-',color="orange",zorder=2)
plt.legend(["Pullback approach (smooth g)","Stress tensor (smooth g)","Pullback approach (non smooth g)","Stress tensor (non smooth g)"],loc="lower right")

plt.loglog(h1[:],np.exp(y_bnd1)[:],'--',color=bestfitcolor,zorder = 1)
plt.loglog(h1[:],np.exp(y_bem1)[:],'--',color=bestfitcolor,zorder = 1)

# Adding best fit lines
fit_bnd = np.polyfit(np.log(h),np.log(errors[:,1]),1)
fit_bem = np.polyfit(np.log(h),np.log(errors[:,0]),1)
y_bnd = fit_bnd[0] * np.log(h) + fit_bnd[1]
y_bem = fit_bem[0] * np.log(h) + fit_bem[1]

bestfitcolor = "silver"
plt.loglog(h,np.exp(y_bnd),'--',color=bestfitcolor,zorder = 1)
plt.loglog(h,np.exp(y_bem),'--',color=bestfitcolor,zorder = 1)

# Adding the slope triangles
annotation.slope_marker((3e-3, 1e-2), (0.42, 1),ax = ax)
annotation.slope_marker((3e-3, 8e-4), (1.23, 1),ax = ax)
#annotation.slope_marker((3e-3, 2.5e-2), (0.56, 1),ax = ax)
annotation.slope_marker((5.2e-3, 1.1e-1), (0.50, 1),invert=True,ax = ax)


print("g Bnd rate:",fit_bnd[0])
print("g BEM rate:",fit_bem[0])
plt.savefig('sqsq_dualnorm.eps')
plt.show()

