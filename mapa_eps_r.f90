program chimera
integer::i,N,t,tmax,seed,trans,j,v1,v2,Nl,aa,bb
parameter(N=301)
real::xold(N+5),xnew(N+5),soma
real::sigma,mi,C,delta,gama,nn
real:: gamamax,gamamin,epsmax,pi
external f,ran0
pi=4.0d0*atan(1.0d0)
!open(10,file='chimera.dat')
!open(12,file='plato.dat')
open(1,file='mapa_eps_gama_0_5_r.dat')
sigma=N/1.0
tmax=2000
mi=N/2.0

Nl=(N-1)/2
 
gama=0.5
 
 
 Nl=(N-1)/2
 
 
 rmax = 4.0
 epsmax = 1.5
 
 eps = 0.05
 do while(eps.LE.epsmax)
 r= 2.5
 do while(r.LE.rmax)
 nn = 0.0
 do i=1,Nl
 	nn = nn + 1/(i**gama)
 end do
 nn=2.0*nn
 do i=1,N
	xnew(i)=0.1+0.01*(1/sqrt(2*pi)*exp(-((i-mi)*(i-mi))/(2*(sigma*sigma))))
 end do


!calculo do transiente
  do t=1,tmax
   	do i=1,N
 	   	xold(i) =xnew(i) 
   	end do
    
	do i=1,N
 	  	 soma =0
 	   
 	  	 do j=1,Nl
 	   
 	    		aa = i+j
	  		bb = i-j 
	   		if(aa>N) then
	   			aa = aa - N
	  		end if
 	   		if(bb < 1) then
		 		bb = bb + N
	   		end if
	   	
	   		soma=soma + (1/j**gama)*(f(r,xold(aa)) - 2*f(r,xold(i))+f(r,xold(bb)))
	   		
	   		
 	   	 end do
 	   	xnew(i)=f(r,xold(i))+ eps/nn *soma  
 	end do
 end do
 
	do i=1,N
	write(1,*) xnew(i)
	end do
	r = r+0.005
	end do
	write(1,*)
 eps = eps+0.0025
   end do
 
 close(1)
 !close(12)
  !close(13)
end program 
function f(r,x)
real:: f,x,r
f = r*x*(1-x)
return 
end

