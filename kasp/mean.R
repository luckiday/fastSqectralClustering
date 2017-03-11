library(kernlab);
library(MASS);
library(gtools);
library(cluster);
library(wskm);
# check the result
cRate=function(sp0, sp1, nc, N)
{
	tr = 0;
	seqs=seq(1,nc);
	spx=matrix(0,N,1);
	spy=matrix(0,N,1);

	if(nc <8) 
	{ perms=permutations(n=nc,r=nc);
		np=dim(perms)[1];
		for(i in 1:np)
		{
			for(j in 1:nc) { spx[sp1==j]=perms[i,j]; }
			tmp=sum(sp0==spx)/N; 
			if(tr<tmp) {tr=tmp;spy=spx;}
		}
	} else {
			for(i in 1:10000)
			{ 
				permx=sample(seqs,nc, replace=FALSE);
				for(j in 1:nc) { spx[sp1==j]=permx[j]; }
				tmp=sum(as.integer(sp0)==spx)/N; if(tr<tmp) {tr=tmp;spy=spx;}
			}
	}
	cat("The rate is ",tr,"\n");
}



espectral=function(x,sp,ncluster,dim)
{


	N=nrow(x);	#The # of observations
	m=dim;		#The # of features
	sp0=sp;
	y=x;
	alpha = 200

	n=floor(length(y[,1])/alpha);

	n1=N*0.05;
	idx= sample((1:N),n1,replace = FALSE);
	xx=y[idx,];

	ptm = proc.time()
	
	# # use coarse k-means
	# cat("begin coarse K-means\n");
	# xxkms=kmeans(xx[,1:m],centers = n, iter.max = 200, nstart = 20,algorithm = c("Hartigan-Wong")); 
	# cat("begin K-means\n")
	# xkms= kmeans(y[,1:m],centers = xxkms$centers, iter.max = 200, nstart = 1, algorithm = c("Hartigan-Wong"));
	# tmp = xkms$cluster;
	# x = xkms$centers;
	# sp = specc(x, centers=3);
	# sp = sp@.Data;
	# spp=xkms$cluster;
	# for(i in 1:n)
	# {
	# 	spp[spp==i]=sp[i];
	# }

	# # use ewkm only
	# cat("begin ewkm \n");
	# xkms = ewkm(y[,1:m],3, lambda = 0.55, maxiter = 100);
	# spp = xkms$cluster;

	# # ewkm + spec
	# xkms = ewkm(y[,1:m],n, lambda = 0.55, maxiter = 200);
	# tmp = xkms$cluster;
	# x = xkms$centers;

	# sp = specc(x, centers=3);
	# sp = sp@.Data;

	# spp=xkms$cluster;
	# for(i in 1:n)
	# {
	# 	spp[spp==i]=sp[i];
	# }


	# # use Clustering Large Applications
	# cat("begin clara \n");
	# xkms = clara(y[,1:m],3, sample = 200);
	# spp = xkms$cluster;


	# # Clustering Large Applications + spec
	cat("begin clara\n");
	xkms = clara(y[,1:m],n,pamLike = TRUE);
	tmp = xkms$clustering;
	x = xkms$medoids;
	sp = specc(x, centers=3);
	sp = sp@.Data;
	spp=xkms$cluster;
	for(i in 1:n)
	{
		spp[spp==i]=sp[i];
	}

	# tmp = clara(y[,1:m],3);
	# spp = tmp$clustering;

	cRate(sp0,spp,ncluster,N);
	

	cat("run time:", proc.time() - ptm,"\n")
}



# loader
main=function(){
	x=read.table(file="connect4.Rdata",header=FALSE,sep=",",strip.white=FALSE);
	x$label=x$V43;
	sp=matrix(0,nrow(x),1);
	sp[x$label=="win"]="1";
	sp[x$label=="loss"]="2";
	sp[x$label=="draw"]="3";
	tmp=c("V43");
	rmcls=match(tmp,names(x));
	x=x[,-rmcls];
	for(i in 1:42) {x[,i]=as.integer(x[,i])};
	nc=3;
	m=42;

	# Normalization
	for(i in 1:m)
	{
		if(sd(x[,i])==0) {x[,i]=0;}
		else {x[,i]=(x[,i]-mean(x[,i]))/sd(x[,i]); }
	}
	# x: data+ label, sp: result, ncluster: 3, 42 feature
	z=espectral(x,sp,ncluster=nc,42);
}

main()