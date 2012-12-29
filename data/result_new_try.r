.libPaths("/mnt/Storage/home/wangxj/R/x86_64-unknown-linux-gnu-library/2.13")
library(Hotelling)

read.table("/mnt/Storage/home/wangxj/GSE6477_rma.txt",row.names=1,header=T)->rma
unique(rma)->rma




c(1,5,7,8,9,10,11,12,15,18,20,21,23,24,25,28,29,32,35,39,40,45,47,48,49,51,53,56,58,61,63,64,69,70,71,72,73,76,77,79,80,81,82,83,84,86,87,91,93,94,95,96,97,104,105,106,107,113,119,122,123,127,128,129,133,136,137,140,143,150)->Hyper1
c(2,3,4,6,13,14,16,17,22,27,30,31,33,34,36,37,39,41,42,43,44,46,50,52,54,55,57,60,62,65,67,68,74,75,78,88,89,90,92,98,99,100,101,103,108,109,111,112,114,115,116,117,120,121,124,126,130,131,132,135,138,139,141,142,144,145,146,147,148,149)->Nonhyper1
read.table("/mnt/Storage/home/wangxj/ref_tf.txt")->TF
TF<-as.vector(t(TF))
intersect(rownames(rma),TF)->TF
read.table("/mnt/Storage/home/wangxj/pathway/arrest.txt")->arrest
arrest<-as.vector(t(arrest))
intersect(rownames(rma),arrest)->arrest


rma[,Hyper1]->Hyper
rma[,Nonhyper1]->Nonhyper
matrix=c()
list(arrest)->ttt
p=c()
for(ref in TF)
{
#ref="NM_002467_at"
	print (ref)
	for (pathway in ttt)
	{
		cor(t(Hyper[ref,]),t(Hyper[,]),method="spearman")->cor_hyper
		cor(t(Nonhyper[ref,]),t(Nonhyper[,]),method="spearman")->cor_nonhyper
		rbind(cor_hyper,cor_nonhyper)->a
		rownames(a)<-c("cor_hyper","cor_nonhyper")
#		plot(t(a) , pch=18 , col='#00001010')
#		abline(a=0, b=1, col = 2)
#		abline(a=0,b=-1,col=2)
		cor(t(Hyper[ref,]),t(Hyper[pathway,]),method="spearman")->cor_hyper2
		cor(t(Nonhyper[ref,]),t(Nonhyper[pathway,]),method="spearman")->cor_nonhyper2
		rbind(cor_hyper2,cor_nonhyper2)->b
		rownames(b)<-c("cor_hyper2","cor_nonhyper2")
#		lines(t(b),pch=18,col='blue',type="p")
#		kmeans(t(b),centers=centers)->cc
#		points(cc$centers, col = 1:2, pch = 8, cex=2)
		
#po=cc$centers
		new1=c()
		new=0
		for (i in c(1:length(p)))
		{
			distance=abs((b[1,i]-b[2,i])/sqrt(2))*sign(abs(b[1,i])-abs(b[2,i]))
			tb=t(b)
			for(i in c(1:length(tb[,1]))){tb[i,1]=(t(b)[i,1]+t(b)[i,2])/sqrt(2);tb[i,2]=(-t(b)[i,1]+t(b)[i,2])/sqrt(2)}
			summary(lm(tb[,1]~tb[,2]))
			
			
			new=distance
#			len=length(which(cc$cluster ==i))
#po[i,]=len/length(cc$cluster)*po[i,]
#new=new*abs(p[i])
			new1=c(new1,new)
			
		}
		new1=new1*min((-log(hotel.test(t(a),t(b))$pval,10)),15) 
#new1=sum(new1)
		
		
				

	}
new1=sum(new1)
	matrix=c(matrix,new1)
}
write.table(matrix,"vect.txt",quote=F,sep="\t")


########Bayerson model

.libPaths("/mnt/Storage/home/wangxj/R/x86_64-unknown-linux-gnu-library/2.13")
library(Hotelling)

read.table("/mnt/Storage/home/wangxj/GSE6477_rma.txt",row.names=1,header=T)->rma
unique(rma)->rma
read.table("/mnt/Storage/home/wangxj/rma2.txt",row.names=1,header=T)->rma2
read.table("/mnt/Storage/home/wangxj/GSE6365ExprData.txt",row.names=1,header=T)->rma3
read.table("/mnt/Storage/home/wangxj/Hyper2.txt",header=F)->Hyper2.names
read.table("/mnt/Storage/home/wangxj/Hyper5.txt")->Hyper3.names
read.table("/mnt/Storage/home/wangxj/Nonhyper5.txt")->Nonhyper3.names
read.table("/mnt/Storage/home/wangxj/Nonhyper2.txt",header=F)->Nonhyper2.names




c(1,5,7,8,9,10,11,12,15,18,20,21,23,24,25,28,29,32,35,39,40,45,47,48,49,51,53,56,58,61,63,64,69,70,71,72,73,76,77,79,80,81,82,83,84,86,87,91,93,94,95,96,97,104,105,106,107,113,119,122,123,127,128,129,133,136,137,140,143,150)->Hyper1
c(2,3,4,6,13,14,16,17,22,27,30,31,33,34,36,37,39,41,42,43,44,46,50,52,54,55,57,60,62,65,67,68,74,75,78,88,89,90,92,98,99,100,101,103,108,109,111,112,114,115,116,117,120,121,124,126,130,131,132,135,138,139,141,142,144,145,146,147,148,149)->Nonhyper1
read.table("/mnt/Storage/home/wangxj/ref_tf.txt")->TF
TF<-as.vector(t(TF))
intersect(rownames(rma),TF)->TF
read.table("/mnt/Storage/home/wangxj/pathway/cellcyclearrest.txt")->arrest
arrest<-as.vector(t(arrest))
intersect(rownames(rma),arrest)->arrest


rma[,Hyper1]->Hyper
rma[,Nonhyper1]->Nonhyper
as.vector(t(Hyper2.names))->Hyper2.names

rma2[,Hyper2.names]->Hyper2
as.vector(t(Nonhyper2.names))->Nonhyper2.names

rma2[,Nonhyper2.names]->Nonhyper2
as.vector(t(Hyper3.names))->Hyper3.names

rma3[,Hyper3.names]->Hyper3
as.vector(t(Nonhyper3.names))->Nonhyper3.names


rma3[,Nonhyper3.names]->Nonhyper3
intersect(rownames(rma),rownames(rma3))->all
intersect(rownames(rma2),all)->all
matrix=c()
matrix2=c()
list(arrest)->ttt
p=c()
arrest<-intersect(all,arrest)
arrest<-intersect(arrest,rownames(rma))
TF<-intersect(TF,rownames(rma3))
for(ref in TF)
{
#ref="NM_002467_at"
	print (ref)
	
	cor(t(Hyper[ref,]),t(Hyper[all,]),method="spearman")->cor_hyper1
	cor(t(Nonhyper[ref,]),t(Nonhyper[all,]),method="spearman")->cor_nonhyper1
	rbind(cor_hyper1,cor_nonhyper1)->b1
	cor(t(Hyper2[ref,]),t(Hyper2[all,]),method="spearman")->cor_hyper2
	cor(t(Nonhyper2[ref,]),t(Nonhyper2[all,]),method="spearman")->cor_nonhyper2
	rbind(cor_hyper2,cor_nonhyper2)->b2
	
	cor(t(Hyper3[ref,]),t(Hyper3[all,]),method="spearman")->cor_hyper3
	cor(t(Nonhyper3[ref,]),t(Nonhyper3[all,]),method="spearman")->cor_nonhyper3
	rbind(cor_hyper3,cor_nonhyper3)->b3
	
	
	
	z=0.5*log((1+b1)/(1-b1))
	z[z>100]=100
	z[1,]->Zr1.h
	z[2,]->Zr1.n
	
	z=0.5*log((1+b2)/(1-b2))
	z[z>100]=100
	
	z[1,]->Zr2.h
	z[2,]->Zr2.n
	
	
	z=0.5*log((1+b3)/(1-b3))
	z[z>100]=100
	
	z[1,]->Zr3.h
	z[2,]->Zr3.n
	
	
	s1=(length(Hyper[1,])-3)^(-1/2)
	s2=(length(Hyper2[1,])-3)^(-1/2)
	s3=(length(Hyper3[1,])-3)^(-1/2)
	cbind(s1,s2,s3)->S.h
	s1.2=(length(Nonhyper[1,])-3)^(-1/2)
	s2.2=(length(Nonhyper2[1,])-3)^(-1/2)
	s3.2=(length(Nonhyper3[1,])-3)^(-1/2)
	
	
	cbind(s1.2,s2.2,s3.2)->S.n
	wi.h=S.h^-2
	wi.n=S.n^-2
	result.h=c()
	result.n=c()
	for(i in c(1:length(all)))
	{
#U.h=sum(wi.h[1]*Zr1.h[i]+wi.h[2]*Zr2.h[i])/sum(wi.h)
		U.h=sum(c(Zr1.h[i],Zr2.h[i],Zr3.h[i])*wi.h)/sum(wi.h)
#U.n=sum(wi.n[1]*Zr1.n[i]+wi.n[2]*Zr2.n[i])/sum(wi.n)
		U.n=sum(c(Zr1.n[i],Zr2.n[i],Zr3.n[i])*wi.n)/sum(wi.n)
#Q.h=wi.h[1]*(Zr1.h[i]-U.h)^2+wi.h[2]*(Zr2.h[i]-U.h)^2
		Q.h=sum(wi.h*(c((Zr1.h[i]-U.h)^2,(Zr2.h[i]-U.h)^2,(Zr3.h[i]-U.h)^2)))
#Q.n=wi.n[1]*(Zr1.n[i]-U.n)^2+wi.n[2]*(Zr2.n[i]-U.n)^2
		Q.n=sum(wi.n*(c((Zr1.n[i]-U.n)^2,(Zr2.n[i]-U.n)^2,(Zr3.n[i]-U.n)^2)))
		tao.sqrh=max((Q.h-1)/(sum(wi.h)-sum(wi.h^2)/sum(wi.h)),0)
		tao.sqrn=max((Q.n-1)/(sum(wi.n)-sum(wi.n^2)/sum(wi.n)),0)
#esti.uh=sum(((S.h[1]^2+tao.sqrh)^-1)*Zr1.h[i]+((S.h[2]^2+tao.sqrh)^-1)*Zr2.h[i])/sum(((S.h[1]^2+tao.sqrh)^-1)+((S.h[2]^2+tao.sqrh)^-1))
		esti.uh=sum(((S.h^2+tao.sqrh)^-1)*c(Zr1.h[i],Zr2.h[i],Zr3.h[i]))/sum((S.h^2+tao.sqrh)^-1)
#esti.un=sum(((S.n[1]^2+tao.sqrn)^-1)*Zr1.n[i]+((S.n[2]^2+tao.sqrn)^-1)*Zr2.n[i])/sum(((S.n[1]^2+tao.sqrn)^-1)+((S.n[2]^2+tao.sqrn)^-1))
		esti.un=sum(((S.n^2+tao.sqrh)^-1)*c(Zr1.n[i],Zr2.n[i],Zr3.n[i]))/sum((S.n^2+tao.sqrh)^-1)
		
		result.h=c(result.h, esti.uh)
		result.n=c(result.n,esti.un)
		
		
		
	}
	
	
	result.hall=(exp(2*result.h)-1)/(exp(2*result.h)+1)
	names(result.hall)<-all
	result.nall=(exp(2*result.n)-1)/(exp(2*result.n)+1)
	names(result.nall)<-all
	
	
	result.h1<-result.hall[arrest]
	result.n1<-result.nall[arrest]
	b=rbind(result.h1,result.n1)
	
	a=rbind(result.hall,result.nall)
	
	new1=c()
	new2=c()
	new=0
	new0=0
	
	distance=abs((b[1,]-b[2,])/sqrt(2))*sign(abs(b[1,])-abs(b[2,]))
	tb=t(b)
	
	for(i in c(1:length(tb[,1]))){tb[i,1]=(t(b)[i,1]+t(b)[i,2])/sqrt(2);tb[i,2]=(-t(b)[i,1]+t(b)[i,2])/sqrt(2)}
	new=summary(lm(tb[,1]~tb[,2]))$coefficients[2,4]
	
	new0=sum(distance)
	new0=new0*min((-log(hotel.test(t(a),t(b))$pval,10)),15)
	matrix=c(matrix,new0)
	

	
	

	matrix2=c(matrix2,new)
}
write.table(matrix,"vect3_new.txt",quote=F,sep="\t")
write.table(matrix2,"vect_P_new.txt",quote=F,sep="\t")


#####Permuation

.libPaths("/mnt/Storage/home/wangxj/R/x86_64-unknown-linux-gnu-library/2.13")
library(Hotelling)

read.table("/mnt/Storage/home/wangxj/GSE6477_rma.txt",row.names=1,header=T)->rma
unique(rma)->rma
read.table("/mnt/Storage/home/wangxj/rma2.txt",row.names=1,header=T)->rma2
read.table("/mnt/Storage/home/wangxj/GSE6365ExprData.txt",row.names=1,header=T)->rma3
read.table("/mnt/Storage/home/wangxj/Hyper2.txt",header=F)->Hyper2.names
read.table("/mnt/Storage/home/wangxj/Hyper5.txt")->Hyper3.names
read.table("/mnt/Storage/home/wangxj/Nonhyper5.txt")->Nonhyper3.names
read.table("/mnt/Storage/home/wangxj/Nonhyper2.txt",header=F)->Nonhyper2.names




#rma[,c(1,5,7,8,9,10,11,12,15,18,20,21,23,24,25,28,29,32,35,39,40,45,47,48,49,51,53,56,58,61,63,64,69,70,71,72,73,76,77,79,80,81,82,83,84,86,87,91,93,94,95,96,97,104,105,106,107,113,119,122,123,127,128,129,133,136,137,140,143,150)]->Hyper
#rma[,c(2,3,4,6,13,14,16,17,22,27,30,31,33,34,36,37,39,41,42,43,44,46,50,52,54,55,57,60,62,65,67,68,74,75,78,88,89,90,92,98,99,100,101,103,108,109,111,112,114,115,116,117,120,121,124,126,130,131,132,135,138,139,141,142,144,145,146,147,148,149)]->Nonhyper

c(1,5,7,8,9,10,11,12,15,18,20,21,23,24,25,28,29,32,35,39,40,45,47,48,49,51,53,56,58,61,63,64,69,70,71,72,73,76,77,79,80,81,82,83,84,86,87,91,93,94,95,96,97,104,105,106,107,113,119,122,123,127,128,129,133,136,137,140,143,150)->Hyper1
c(2,3,4,6,13,14,16,17,22,27,30,31,33,34,36,37,39,41,42,43,44,46,50,52,54,55,57,60,62,65,67,68,74,75,78,88,89,90,92,98,99,100,101,103,108,109,111,112,114,115,116,117,120,121,124,126,130,131,132,135,138,139,141,142,144,145,146,147,148,149)->Nonhyper1
#c(1,5,7,8,9,10,11,12,15,18,19,20,22,23,24,26,27,30,33,37,38,43,45,46,47,49,51,54,56,58,60,61,65,66,67,68,69,72,73,75,76,77,78,79,80,81,82,86,88,89,90,91,92,98,99,100,101,106,111,114,115,118,119,120,124,126,127,130,133,140)->Hyper1
#c(2,3,4,6,13,14,16,17,21,25,28,29,31,32,34,35,36,39,40,41,42,44,48,50,52,53,55,57,59,62,63,64,70,71,74,83,84,85,87,93,94,95,96,97,102,103,104,105,107,108,109,110,112,113,116,117,121,122,123,125,128,129,131,132,134,135,136,137,138,139)->Nonhyper1
read.table("/mnt/Storage/home/wangxj/ref_tf.txt")->TF
TF<-as.vector(t(TF))
intersect(rownames(rma),TF)->TF
read.table("/mnt/Storage/home/wangxj/pathway/telomere maintenance.txt")->telomere_maintenance
telomere_maintenance<-as.vector(t(telomere_maintenance))
intersect(rownames(rma),telomere_maintenance)->telomere_maintenance
read.table("/mnt/Storage/home/wangxj/pathway/growth factor activity.txt")->growth_factor_activity

growth_factor_activity<-as.vector(t(growth_factor_activity))
intersect(rownames(rma),growth_factor_activity)->growth_factor_activity
read.table("/mnt/Storage/home/wangxj/pathway/cell cycle.txt")->cell_cycle
cell_cycle<-as.vector(t(cell_cycle))
intersect(rownames(rma),cell_cycle)->cell_cycle
read.table("/mnt/Storage/home/wangxj/pathway/response to DNA damage stimulus.txt")->response_to_DNA_damage_stimulus
response_to_DNA_damage_stimulus<-as.vector(t(response_to_DNA_damage_stimulus))
intersect(rownames(rma),response_to_DNA_damage_stimulus)->response_to_DNA_damage_stimulus
read.table("/mnt/Storage/home/wangxj/pathway/cell division.txt")->cell_division
cell_division<-as.vector(t(cell_division))
intersect(rownames(rma),cell_division)->cell_division
read.table("/mnt/Storage/home/wangxj/pathway/Wnt signaling.txt")->Wnt_signaling

Wnt_signaling<-as.vector(t(Wnt_signaling))
intersect(rownames(rma),Wnt_signaling)->Wnt_signaling
read.table("/mnt/Storage/home/wangxj/pathway/cell differentiation.txt")->cell_differentiation

cell_differentiation<-as.vector(t(cell_differentiation))
intersect(rownames(rma),cell_differentiation)->cell_differentiation

read.table("/mnt/Storage/home/wangxj/pathway/Apoptosis.txt")->Apoptosis
Apoptosis<-as.vector(t(Apoptosis))
intersect(rownames(rma),Apoptosis)->Apoptosis
read.table("/mnt/Storage/home/wangxj/pathway/cellcyclearrest.txt")->arrest
arrest<-as.vector(t(arrest))
intersect(rownames(rma),arrest)->arrest
read.table("/mnt/Storage/home/wangxj/pathway/mit.txt")->mit
mit<-as.vector(t(mit))
intersect(rownames(rma),mit)->mit


#rma[,Hyper1]->Hyper
#rma[,Nonhyper1]->Nonhyper
as.vector(t(Hyper2.names))->Hyper2.names

#rma2[,Hyper2.names]->Hyper2
as.vector(t(Nonhyper2.names))->Nonhyper2.names

#rma2[,Nonhyper2.names]->Nonhyper2
as.vector(t(Hyper3.names))->Hyper3.names

#rma3[,Hyper3.names]->Hyper3
as.vector(t(Nonhyper3.names))->Nonhyper3.names


#rma3[,Nonhyper3.names]->Nonhyper3
matrix3=c()


for(count in c(1:100))

{print (count)
	
	sample2<-c(Hyper2.names,Nonhyper2.names)
	sample3<-c(Hyper3.names,Nonhyper3.names)
	
	
	sample1<-c(Hyper1,Nonhyper1)
	sample(sample1,size=70)->random1
	rma[,random1]->Hyper
	rma[,setdiff(sample1,random1)]->Nonhyper
	random2<-sample(sample2,length(Hyper2.names))
	rma2[,random2]->Hyper2
	rma2[,setdiff(sample2,random2)]->Nonhyper2
	random3<-sample(sample3,length(Hyper3.names))
	rma3[,random3]->Hyper3
	rma3[,setdiff(sample3,random3)]->Nonhyper3
	
	
	intersect(rownames(rma),rownames(rma3))->all
	intersect(rownames(rma2),all)->all
	matrix=c()
	matrix2=c()
#list(mit,telomere_maintenance,growth_factor_activity,cell_cycle,response_to_DNA_damage_stimulus,cell_division,Wnt_signaling,cell_differentiation,Apoptosis,arrest)->ttt
	list(arrest)->ttt
	p=c()
	arrest<-intersect(all,arrest)
	arrest<-intersect(arrest,rownames(rma))
#for(m in arrest){t=t.test(Hyper[m,],Nonhyper[m,])$p.value;t=-log(t,10);t=t*sign(rowMeans(Hyper[m,])-rowMeans(Nonhyper[m,]));p=c(p,t)}
	TF<-intersect(TF,rownames(rma3))
	for(ref in TF)
	{
		print (ref)
		
		cor(t(Hyper[ref,]),t(Hyper[all,]),method="spearman")->cor_hyper1
		cor(t(Nonhyper[ref,]),t(Nonhyper[all,]),method="spearman")->cor_nonhyper1
		rbind(cor_hyper1,cor_nonhyper1)->b1
		cor(t(Hyper2[ref,]),t(Hyper2[all,]),method="spearman")->cor_hyper2
		cor(t(Nonhyper2[ref,]),t(Nonhyper2[all,]),method="spearman")->cor_nonhyper2
		rbind(cor_hyper2,cor_nonhyper2)->b2
		
		cor(t(Hyper3[ref,]),t(Hyper3[all,]),method="spearman")->cor_hyper3
		cor(t(Nonhyper3[ref,]),t(Nonhyper3[all,]),method="spearman")->cor_nonhyper3
		rbind(cor_hyper3,cor_nonhyper3)->b3
		
		
		
		z=0.5*log((1+b1)/(1-b1))
		z[z>100]=100
		z[1,]->Zr1.h
		z[2,]->Zr1.n
		
		z=0.5*log((1+b2)/(1-b2))
		z[z>100]=100
		
		z[1,]->Zr2.h
		z[2,]->Zr2.n
		
		
		z=0.5*log((1+b3)/(1-b3))
		z[z>100]=100
		
		z[1,]->Zr3.h
		z[2,]->Zr3.n
		
		
		s1=(length(Hyper[1,])-3)^(-1/2)
		s2=(length(Hyper2[1,])-3)^(-1/2)
		s3=(length(Hyper3[1,])-3)^(-1/2)
		cbind(s1,s2,s3)->S.h
		s1.2=(length(Nonhyper[1,])-3)^(-1/2)
		s2.2=(length(Nonhyper2[1,])-3)^(-1/2)
		s3.2=(length(Nonhyper3[1,])-3)^(-1/2)
		
		
		cbind(s1.2,s2.2,s3.2)->S.n
		wi.h=S.h^-2
		wi.n=S.n^-2
		result.h=c()
		result.n=c()
		for(i in c(1:length(all)))
		{
#U.h=sum(wi.h[1]*Zr1.h[i]+wi.h[2]*Zr2.h[i])/sum(wi.h)
			U.h=sum(c(Zr1.h[i],Zr2.h[i],Zr3.h[i])*wi.h)/sum(wi.h)
#U.n=sum(wi.n[1]*Zr1.n[i]+wi.n[2]*Zr2.n[i])/sum(wi.n)
			U.n=sum(c(Zr1.n[i],Zr2.n[i],Zr3.n[i])*wi.n)/sum(wi.n)
#Q.h=wi.h[1]*(Zr1.h[i]-U.h)^2+wi.h[2]*(Zr2.h[i]-U.h)^2
			Q.h=sum(wi.h*(c((Zr1.h[i]-U.h)^2,(Zr2.h[i]-U.h)^2,(Zr3.h[i]-U.h)^2)))
#Q.n=wi.n[1]*(Zr1.n[i]-U.n)^2+wi.n[2]*(Zr2.n[i]-U.n)^2
			Q.n=sum(wi.n*(c((Zr1.n[i]-U.n)^2,(Zr2.n[i]-U.n)^2,(Zr3.n[i]-U.n)^2)))
			tao.sqrh=max((Q.h-1)/(sum(wi.h)-sum(wi.h^2)/sum(wi.h)),0)
			tao.sqrn=max((Q.n-1)/(sum(wi.n)-sum(wi.n^2)/sum(wi.n)),0)
#esti.uh=sum(((S.h[1]^2+tao.sqrh)^-1)*Zr1.h[i]+((S.h[2]^2+tao.sqrh)^-1)*Zr2.h[i])/sum(((S.h[1]^2+tao.sqrh)^-1)+((S.h[2]^2+tao.sqrh)^-1))
			esti.uh=sum(((S.h^2+tao.sqrh)^-1)*c(Zr1.h[i],Zr2.h[i],Zr3.h[i]))/sum((S.h^2+tao.sqrh)^-1)
#esti.un=sum(((S.n[1]^2+tao.sqrn)^-1)*Zr1.n[i]+((S.n[2]^2+tao.sqrn)^-1)*Zr2.n[i])/sum(((S.n[1]^2+tao.sqrn)^-1)+((S.n[2]^2+tao.sqrn)^-1))
			esti.un=sum(((S.n^2+tao.sqrh)^-1)*c(Zr1.n[i],Zr2.n[i],Zr3.n[i]))/sum((S.n^2+tao.sqrh)^-1)
			
			result.h=c(result.h, esti.uh)
			result.n=c(result.n,esti.un)
			
			
			
		}
		
		
		result.hall=(exp(2*result.h)-1)/(exp(2*result.h)+1)
		names(result.hall)<-all
		result.nall=(exp(2*result.n)-1)/(exp(2*result.n)+1)
		names(result.nall)<-all
		
		
		result.h1<-result.hall[arrest]
		result.n1<-result.nall[arrest]
		b=rbind(result.h1,result.n1)
		
		a=rbind(result.hall,result.nall)
		
#po=cc$centers
		new1=c()
		new2=c()
		new=0
		new0=0
		
		distance=abs((b[1,]-b[2,])/sqrt(2))*sign(abs(b[1,])-abs(b[2,]))
		tb=t(b)
		
		for(i in c(1:length(tb[,1]))){tb[i,1]=(t(b)[i,1]+t(b)[i,2])/sqrt(2);tb[i,2]=(-t(b)[i,1]+t(b)[i,2])/sqrt(2)}
		new=summary(lm(tb[,1]~tb[,2]))$coefficients[2,4]
		
		new0=sum(distance)
		new0=new0*min((-log(hotel.test(t(a),t(b))$pval,10)),15)
		matrix=c(matrix,new0)
		
		matrix2=c(matrix2,new)
	}
	matrix3=cbind(matrix3,matrix)	
}
write.table(matrix3,"vectper.txt",quote=F,sep="\t")
