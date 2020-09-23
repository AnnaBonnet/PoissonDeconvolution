

h_hat_func=function(Hseq,Y,a,n,nb_cores,Tmax=1, etaseq,gamma=0.01)
{
  hh=matrix(rep(Hseq,each=length(Hseq)),length(Hseq),length(Hseq))
  
  Bth=matrix(mcmapply(FUN=function(tt,h)
  {
    norm2_t=calc_norm2_t_sym_Rcpp(tt,Y,a,n,Tmax)
    norm2_t2h2=calc_norm2_t2h2_sym_Rcpp(tt,h,Y,a,n,Tmax)
    norm2_crois=calc_norm2_crois_sym_Rcpp(tt,h,Y,a,n,Tmax)
    g_mat_t=3*a/n*sqrt(norm2_t/tt^6-2*norm2_crois/((h^2+tt^2)^1.5*tt^3)+norm2_t2h2/(h^2+tt^2)^3)
    return(g_mat_t)
  },hh,t(hh),mc.cores=nb_cores),length(Hseq),length(Hseq))
  
  pen = matrix(NA,length(etaseq),length(Hseq))
  for (hh in 1:length(etaseq))
  {
    for (kk in 1:length(Hseq)){
      pen[hh,kk] = ( (1+etaseq[hh])*2*sqrt(3/5)*sqrt(a*Tmax/2) ) * sqrt(N) / (Hseq[kk]^(3/2)*n)		
    }
  }
  w = mapply(1:length(etaseq), FUN = function(ee){
    Ah = apply(sweep(Bth,2,pen[ee,],"-"),1,max)
    Ah + gamma*pen[ee,]          	
  })
  h=Hseq[which.min(apply(w,1,min))]
  return(h)
  
}
