standardized_mean<-function(m.1,sd.1,n.1,m.2,sd.2,n.2){
  sd_pooled=sqrt(((n.1-1)*sd.1^2+(n.2-1)*sd.2^2)/(n.1+n.2-2))
  (m.1-m.2)/sd_pooled
}
library(MASS)

complex<-read.table("~/Desktop/Complexity_FIle.txt", header=1, sep='\t')
head(complex)

sub1<- complex[1:24,]
sub2<- complex[25:48,]
sub3<-complex[49:120,]
case1<-apply(sub1[1:4,3:22],2,sum)
case2<-apply(sub1[5:8,3:22],2,sum)
case3<-apply(sub1[9:12,3:22],2,sum)
case4<-apply(sub1[13:16,3:22],2,sum)
case5<-apply(sub1[17:20,3:22],2,sum)
case6<-apply(sub1[21:24,3:22],2,sum)


s1<-rbind(case1,case2,case3,case4,case5,case6)


case7<-apply(sub2[1:4,3:22],2,sum)
case8<-apply(sub2[5:8,3:22],2,sum)
case9<-apply(sub2[9:12,3:22],2,sum)
case10<-apply(sub2[13:16,3:22],2,sum)
case11<-apply(sub2[17:20,3:22],2,sum)
case12<-apply(sub2[21:24,3:22],2,sum)

s2<-rbind(case7,case8,case9,case10,case11,case12)


case13<-apply(sub3[1:4,3:22],2,sum)
case14<-apply(sub3[5:8,3:22],2,sum)
case15<-apply(sub3[9:12,3:22],2,sum)
case16<-apply(sub3[13:16,3:22],2,sum)
case17<-apply(sub3[17:20,3:22],2,sum)
case18<-apply(sub3[21:24,3:22],2,sum)
case19<-apply(sub3[25:28,3:22],2,sum)
case20<-apply(sub3[29:32,3:22],2,sum)
case21<-apply(sub3[33:36,3:22],2,sum)
case22<-apply(sub3[37:40,3:22],2,sum)
case23<-apply(sub3[41:44,3:22],2,sum)
case24<-apply(sub3[45:48,3:22],2,sum)
case25<-apply(sub3[48:52,3:22],2,sum)
case26<-apply(sub3[53:56,3:22],2,sum)
case27<-apply(sub3[57:60,3:22],2,sum)
case28<-apply(sub3[61:64,3:22],2,sum)
case29<-apply(sub3[65:68,3:22],2,sum)
case30<-apply(sub3[69:72,3:22],2,sum)


s3<-rbind(case13,case14,case15,case16,case17,case18,case19,case20,case21,case22,case23,case24,case25,case26,case27,case28,case29,case30)


chisq.test(margin.table(s1, 2))
chisq.test(margin.table(s2, 2))
chisq.test(margin.table(s3, 2))

s_all<-rbind(s1,s2,s3)
summary(s_all)
phi_coeff<-function(mat){
a=sum(mat[,1]==0&mat[,2]==0)
b=sum(mat[,1]==0&mat[,2]!=0)
c=sum(mat[,1]!=0&mat[,2]==0)
d=sum(mat[,1]!=0&mat[,2]!=0)
phi(c(a,b,c,d))
}
colnames(s_all)
res<-matrix(0,ncol(s_all),ncol(s_all))
rownames(res)=colnames(res)=colnames(s_all)
for(i in 1:ncol(s_all)){
  for(j in 1:ncol(s_all)){
   res[i,j]=phi_coeff(cbind(s_all[,i],s_all[,j]))
    }
  }
heatmap(res)
write.table(res,file="~/Desktop/Phi.txt",sep='\t',quote=F,col.names=NA)



