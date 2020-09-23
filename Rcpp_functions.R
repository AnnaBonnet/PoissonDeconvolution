
cppFunction('
            double calc_norm2_crois_sym_Rcpp(double t, double h, NumericVector Y, double a,int n, double Tmax){
            int N=Y.size();
            double sum_tot=0;
            double test=0;
            for(int i=0; i<N; i++) 
            for(int j=0; j<N; j++)
            {  
            double s_max1=(1-Y[i]+t)/(2*a)-0.5 ;
            double s_max2=(1-Y[j]+sqrt(h*h+t*t))/(2*a)-0.5 ;
            for(int k1=-s_max1; k1<=s_max1; k1++)
            for(int k2=-s_max2; k2<=s_max2; k2++){
            int s_k1=-1 ;
            int s_k2=-1 ;
            if(k1 >=0){s_k1=1;} 
            if(k2 >=0){ s_k2=1;} 
            double int_sup=min(NumericVector::create((Y[i]+(2*k1+1)*a+t), (Y[j]+(2*k2+1)*a+sqrt(h*h+t*t)),Tmax)) ;
            double int_inf=max(NumericVector::create((Y[i]+(2*k1+1)*a-t), (Y[j]+(2*k2+1)*a-sqrt(h*h+t*t)),0));
            if(int_sup>int_inf){
            sum_tot=sum_tot+ 0.25*s_k1*s_k2*(((Y[i]+(2*k1+1)*a)*(Y[j]+(2*k2+1)*a)*(int_sup-int_inf)-0.5*(int_sup*int_sup-int_inf*int_inf)*(Y[i]+Y[j]+2*(k1+k2+1)*a)+(int_sup*int_sup*int_sup-int_inf*int_inf*int_inf)/3));
            }
            }
            }
            return sum_tot;
            }
            ')




cppFunction('
            double calc_norm2_t_sym_Rcpp(double t, NumericVector Y, double a,int n, double Tmax){
            int N=Y.size();
            double sum_tot=0;
            double test=0;
            for(int i=0; i<N; i++) 
            for(int j=0; j<N; j++)
            {  
            double s_max1=(1-Y[i]+t)/(2*a)-0.5 ;
            double s_max2=(1-Y[j]+t)/(2*a)-0.5 ;
            for(int k1=-s_max1; k1<=s_max1; k1++)
            for(int k2=-s_max2; k2<=s_max2; k2++){
            int s_k1=-1 ;
            int s_k2=-1 ;
            if(k1 >=0){s_k1=1;} 
            if(k2 >=0){ s_k2=1;} 
            double int_sup=min(NumericVector::create((Y[i]+(2*k1+1)*a+t), (Y[j]+(2*k2+1)*a+t),Tmax)) ;
            double int_inf=max(NumericVector::create((Y[i]+(2*k1+1)*a-t), (Y[j]+(2*k2+1)*a-t),0));
            if(int_sup>int_inf){
            sum_tot=sum_tot+ 0.25*s_k1*s_k2*(((Y[i]+(2*k1+1)*a)*(Y[j]+(2*k2+1)*a)*(int_sup-int_inf)-0.5*(int_sup*int_sup-int_inf*int_inf)*(Y[i]+Y[j]+2*(k1+k2+1)*a)+(int_sup*int_sup*int_sup-int_inf*int_inf*int_inf)/3));
            }
            }
            }
            return sum_tot;
            }
            ')


cppFunction('
            double calc_norm2_t2h2_sym_Rcpp(double t, double h, NumericVector Y, double a,int n, double Tmax){
            int N=Y.size();
            double sum_tot=0;
            double test=0;
            for(int i=0; i<N; i++) 
            for(int j=0; j<N; j++)
            {  
            double s_max1=(1-Y[i]+sqrt(h*h+t*t))/(2*a)-0.5 ;
            double s_max2=(1-Y[j]+sqrt(h*h+t*t))/(2*a)-0.5 ;
            for(int k1=-s_max1; k1<=s_max1; k1++)
            for(int k2=-s_max2; k2<=s_max2; k2++){
            int s_k1=-1 ;
            int s_k2=-1 ;
            if(k1 >=0){s_k1=1;} 
            if(k2 >=0){ s_k2=1;} 
            double int_sup=min(NumericVector::create((Y[i]+(2*k1+1)*a+sqrt(h*h+t*t)), (Y[j]+(2*k2+1)*a+sqrt(h*h+t*t)),Tmax)) ;
            double int_inf=max(NumericVector::create((Y[i]+(2*k1+1)*a-sqrt(h*h+t*t)), (Y[j]+(2*k2+1)*a-sqrt(h*h+t*t)),0));
            if(int_sup>int_inf){
            sum_tot=sum_tot+0.25*s_k1*s_k2*(((Y[i]+(2*k1+1)*a)*(Y[j]+(2*k2+1)*a)*(int_sup-int_inf)-0.5*(int_sup*int_sup-int_inf*int_inf)*(Y[i]+Y[j]+2*(k1+k2+1)*a)+(int_sup*int_sup*int_sup-int_inf*int_inf*int_inf)/3));
            }
            }
            }
            return sum_tot;
            }
            ')


