#******************************************************************************************************************************************
#******************************************************************************************************************************************
# Generates N realizations of a random field with correlation structure defined 
#******************************************************************************************************************************************
#******************************************************************************************************************************************
setwd("D:/berenice/Documents/GitHub/1D_land_subsidence_with_ensemble_Kalman_filter/src/SSG")

datos<-read.table("D:/berenice/Documents/GitHub/1D_land_subsidence_with_ensemble_Kalman_filter/src/SSG/IN/diff_ecs.dat")
kRF<-read.table("D:/berenice/Documents/GitHub/1D_land_subsidence_with_ensemble_Kalman_filter/src/SSG/IN/K.txt")


L=datos[5,]							
nnd=datos[1,]						
dz=L/(nnd-1)						
N=2010							

#******************************************************************************************************************************************
# Mexico Plain:
media=-8.964 						# Average random field (from Vargas y Ortega,2003)
sd=(0.871) 	# define it in logarithm base 10	# Standard deviation (from Vargas y Ortega,2003)
sd=sd*2.303	# lnK=2.303*log10K
var=sd*sd							# Variance random field

		library(RandomFields) 			
		library(abind)

		RFoptions(seed=NA)

								# RMspheric(var, scale, Aniso, proj) 
								# RMexp(var, scale, Aniso, proj) 
								# RMnugget(tol, vdim, var, scale, Aniso, proj) 
								# RFvariogram(model, x, y = NULL, z = NULL, T = NULL, grid, distances, dim, ...)

		model<-RMexp(var=var,scale=2.1)	
		z<-seq(0,L,dz)
		Re<-RFsimulate(model=RPsequential(phi=model),z,n=N) 
		#Re<-RFsimulate(model=model,z,n=N)	# Re are n realizations of ln(K)
		lnK=as.matrix(Re)				# ln(K)
		mean(lnK)			
		sd(lnK)*sd(lnK)
		Kin<-matrix(kRF[,1],nnd,N)		# kRF of RyF (1991)
		lnK=(lnK)+(log(Kin))			# lnK + ln average RyF(1991)
		k=exp(lnK)					# k
		

# To graph the re-scaled realizations:

Mk<-matrix(0,1,nnd)
for(i in 1:(nnd)){			
	Mk[i]=sum(log(k[i,]))/N
}

	
plot(log(k[,1]),z,type="l",lwd=2,col="NA",xlim=c(-30,-10),ylim=c(0,15),yaxp=c(0,15,5),cex.axis=1,xlab="Y = ln K' (m/s)",ylab="Posici�n en z (m)")
for(i in 1:N){
	points(log(k[,i]),z,type="l",lwd=0.5,col="gray40")#,log="x")
}
points(log(kRF[,1]),z,type="l",lwd=2,col="red")#,log="x")
points(Mk[1,],z,type="l",lwd=3,col="blue")#,log="x")
	

# To export the realizations:

		rea_k=abind(Re(k[,1]),Re(k[,2]),along=1)
		for(i in 3:N){
			rea_k=abind(rea_k,Re(k[,i]),along=1)
		}

write.table(rea_k,file="realizaciones_k.txt",row.names=FALSE,col.names=FALSE)


#************************************************************************************************************
#************************************************************************************************************
# An alternative way to produce the Y realizations: gstat 
#************************************************************************************************************
#************************************************************************************************************
library(sp)
library(gstat)
library(abind)

# Mexico Plain:
media=-8.964 						# Average random field (from Vargas y Ortega,2003)
sd=(0.871) 	# define it in logarithm base 10	# Standard deviation (from Vargas y Ortega,2003)
sd=sd*2.303	# lnK=2.303*log10K
var=sd*sd							# Variance random field

z<-seq(0,L,dz)


# Sequential Gaussian Simulation
xy <- expand.grid(z, seq(1, 1, by = 1))
names(xy) <- c("x","y")
gridded(xy) = ~x+y
x<-xy@coords[,1]
y<-xy@coords[,2]
variogram<-vgm(var, "Sph", range=2.1, nugget=0)
g.dummy <- gstat(formula = z~1, dummy = TRUE, beta = 0, model = variogram, nmax = 30, maxdist = 2.1) # for speed -- 10 is too small!!
yy <- predict(g.dummy, xy, nsim = N)

sim<-yy@data
sim

# To graph the Y realizations with mean cero:

Msim<-Vsim<-matrix(0,1,nnd)
for(i in 1:(nnd)){						
	Msim[i]=sum(sim[i,])/N
}
Msim

for(i in 1:nnd){
	Vsim[i]=(sum((sim[i,1:N]-Msim[i])^2))/N
}
Vsim

SDsim<-sqrt(Vsim)


plot(sim[,1],z,col=NA,xlim=c(-10,10))
for(i in 1:N){
	points(sim[,i],z,type="l",col="grey40")
}
points(Msim,z,type="l",col="red")
points(Msim+2*SDsim,z,type="l",col="orange")
points(Msim-2*SDsim,z,type="l",col="orange")


# To re-scale the realizations:

Kin<-matrix(kRF[,1],nnd,N)		# kRF of RyF (1991)
lnK=(sim)+(log(Kin))			# lnK + ln average RyF(1991)
k=exp(lnK)					# k


# To graph the rescaled realizations:

Mk<-matrix(0,1,nnd)
for(i in 1:(nnd)){				
	Mk[i]=sum(log(k[i,]))/N
}

	
plot(log(k[,1]),z,type="l",lwd=2,col="NA",xlim=c(-30,-10),ylim=c(0,15),yaxp=c(0,15,5),cex.axis=1,xlab="Y = ln K' (m/s)",ylab="Posici�n en z (m)")
for(i in 1:N){
	points(log(k[,i]),z,type="l",lwd=0.5,col="gray40")#,log="x")
}
points(log(kRF[,1]),z,type="l",lwd=2,col="red")#,log="x")
points(Mk[1,],z,type="l",lwd=3,col="blue")#,log="x")


# To export the realizations:

		rea_k=abind(Re(k[,1]),Re(k[,2]),along=1)
		for(i in 3:N){
			rea_k=abind(rea_k,Re(k[,i]),along=1)
		}

write.table(rea_k,file="realizaciones_k.txt",row.names=FALSE,col.names=FALSE)



