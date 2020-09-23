
fh_estim=function(h_hat,Y,a,n,seq_x)
{
    F_h=function(x,h_hat,Y,a,n)
    {
      N=length(Y)
      sum_f=0
      sum_f_sym=0
      s_max=matrix(0,N)
      s_min=matrix(0,N)
      s_min_ceil=matrix(0,N)
      s_max_floor=matrix(0,N)
      s_max_sym=matrix(0,N)
      s_min_sym=matrix(0,N)
      s_min_ceil_sym=matrix(0,N)
      s_max_floor_sym=matrix(0,N)
      for(i in 1:N)
      {
        s_min_sym[i]=(-x+Y[i]-h_hat)/(2*a)-1/2
        s_min_ceil_sym[i]=max(ceiling(s_min_sym[i]),0)
        s_max_sym[i]=(-x+Y[i]+h_hat)/(2*a)-1/2
        s_max_floor_sym[i]=floor(s_max_sym[i])
        
        
        if(s_max_floor_sym[i]>=0 && s_max_floor_sym[i]>=s_min_ceil_sym[i]){
          for (k in s_min_ceil_sym[i]:s_max_floor_sym[i])
          { 
            sum_f_sym=sum_f_sym- 3*a/(n*h_hat^3)*(Y[i]-(2*k+1)*a-x)
          }
          
          
          
        }
        
        s_min[i]=(x-Y[i]-h_hat)/(2*a)-1/2
        s_min_ceil[i]=max(ceiling(s_min[i]),0)
        s_max[i]=(x-Y[i]+h_hat)/(2*a)-1/2
        s_max_floor[i]=floor(s_max[i])
        if(s_max_floor[i]>=0 && s_max_floor[i]>=s_min_ceil[i]){
          
          for (k in s_min_ceil[i]:s_max_floor[i])
          { 
            sum_f=sum_f+ 3*a/(n*h_hat^3)*(Y[i]+(2*k+1)*a-x)
            #cat(sum_f)
            # cat("nb ",count)
          }
        }
      }
      fh=0.5*sum_f+ 0.5*sum_f_sym
      return(sum_f)
    }
  
  
  
  f_func_r=function(x){F_h(x,h_hat,Y,a,n)}
  lambda_hat=sapply(seq_x,f_func_r)
  return(lambda_hat)
}