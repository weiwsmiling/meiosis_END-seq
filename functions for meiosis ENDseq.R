BinMat<-function(mat,binsize){
nbins<-floor(ncol(mat)/binsize)
bmat<-matrix(0,nrow(mat),nbins)
for (i in 1:nbins){
	bmat[,i]=rowSums(mat[,((i-1)*binsize+1):(i*binsize)])
	}
	return(bmat)
}

MeanResection<-function(ba,bs){
	res_mean<-matrix(0,nrow(ba),2)
	for(i in 1:nrow(ba)){
		tmp<-ba[,1000:1]
		res_mean[i,1]=(((tmp[i,521:900]/sum(tmp[i,521:900]))%*%matrix(1:380,ncol=1))+20)*10
		res_mean[i,2]=(((bs[i,521:900]/sum(bs[i,521:900]))%*%matrix(1:380,ncol=1))+20)*10
}
	return(res_mean)
}
Resection<-function(ba,bs){
	res<-matrix(0,nrow(bs),2)
	for (i in 1:nrow(res)){
		ba[i,]<-smooth.spline((ba)[i,],spar=0.2)$y
		bs[i,]<-smooth.spline((bs)[i,],spar=0.2)$y
		ind<-order(ba[i,],decreasing=T)[1:5]
		k<-cutree(hclust(dist(ind)),h=100)
		ba.max=ind[k==which.max(table(k))][1]
		ba.bg=max(ba[i,1:50])
		ind<-order(bs[i,],decreasing=T)[1:5]
		k<-cutree(hclust(dist(ind)),h=100)
		bs.max=ind[k==which.max(table(k))][1]
		bs.bg=max(bs[i,551:600])
		if ((bs.max>551)|(ba.max<50)|(ba.max>bs.max)|(ba.max>301)|(bs.max<300)){res[i,]=(-1)}else{
			j=bs.max
			a=0
			bs.tmp=bs[i,(j):(j+19)]
			while(length(which(bs.tmp>bs.bg))>=10){
				j<-(j+1)
				bs.tmp=bs[i,(j):(j+19)]
				a=max(which(bs.tmp>bs.bg))
			}
			res[i,2]=(j+a)*20-6000
			m=ba.max
			ba.tmp=ba[i,(m-19):(m)]
			b=20
			while(length(which(ba.tmp>ba.bg))>=10){
				m<-(m-1)
				ba.tmp=ba[i,(m-19):(m)]
				b=min(which(ba.tmp>ba.bg))
			}
			res[i,1]=6000-(b-20+m)*20
		}
	}
	return(res)
}

GetGap2<-function(ba,bs,co=0.05){
	gp<-matrix(0,nrow(bs),2)
	for (i in 1:nrow(gp)){
		ind<-order(ba[i,],decreasing=T)[1:5]
		k<-cutree(hclust(dist(ind)),h=100)
		ba.max=ind[k==which.max(table(k))][1]
		ind<-order(bs[i,],decreasing=T)[1:5]
		k<-cutree(hclust(dist(ind)),h=100)
		bs.max=ind[k==which.max(table(k))][1]
		if ((bs.max>900)|(ba.max<100)|(ba.max>bs.max)){gp[i,]=(-1)}else{
		ba.mu=mean(log(ba[i,(ba.max-50):(ba.max+50)]+0.01))
		bs.mu=mean(log(bs[i,(bs.max-50):(bs.max+50)]+0.01))
		p=1
		j=ba.max
		while(p>co){
			p=t.test(log(ba[i,(j):(j+10)]+0.01)+rnorm(11,0,0.0001),mu=ba.mu+rnorm(1,0,0.0001),alternative="less")$p.value
			j=j+1
		}
		gp[i,1]=j+5
			j=bs.max
			p=1
		while(p>co){
			p=t.test(log(bs[i,(j-10):(j)]+0.01)+rnorm(11,0,0.0001),mu=bs.mu+rnorm(1,0,0.0001),alternative="less")$p.value
			j=j-1
		}
		gp[i,2]=j-5
			
	}
	}
	return(gp)
}
