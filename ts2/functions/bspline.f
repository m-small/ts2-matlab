C-------------------------------------------------------------------------------
C  	bspline.f
C-------------------------------------------------------------------------------
C
C	This is a mex file for MATLAB. It calculates the values of the cardinal 
C	B-spline scaling function for a given (vector) set of parameters (for a 
C 	vector of input values), used to specify the dilation, translation and 
C	order parameter of the function. The function is normalised to have 
C	constant L^2 norm for different dilation and translation paramters.
C	   See "Fundamentals of Wavelets", pp.106 by Goswami abd Chan.
C
C	MATLAB call: y = bspline(x,pms);
C
C	 INPUT:
C	   'x' vector of input arguments to function
C	   'pms' vector of function parameters
C	
C	 OUTPUT:
C	   'y' function values
C
C	Notes:
C	   'x' is a 1-by-n or n-by-1 vector
C	   'y' is 1-by-n
C	   y(i) is the evaluation of the function at x(i) with parameters values
C	   pms(:,i) (or pms(:) if pms is 2-by-1).
C
C	The format of the vector 'pms' is as follows:
C		pms = [a; b];
C	with
C	   'a' dilation parameters (size = [1 n])
C	   'b' translation parameters (size = [1 n])
C
C	The order of the spline, parameter 'ord', is taken in from a global 
C	variable 'SPLINE_ORDER'	in the global MATLAB workspace.
C	   'ord' order of spline (polynomial of degree ord-1) 
C
C	NOTES: If a vector 'x' is given with only one paramter 'ord'. We 
C	calculate the output of the B-spline for values 'x' and order 'ord'. If
C	'pms' (2-by-1) is given as only two parameters, those parameters are 
C	taken for all values of 'x'. The functions are currently centred at zero
C	too.
C
C	Brian Fleming, 11/6/00.
C-------------------------------------------------------------------------------

C-------------------------------------------------------------------------------
C	Set mex requirments
C-------------------------------------------------------------------------------

	subroutine mexFunction(nlhs,plhs,nrhs,prhs)
	
	integer plhs(*), prhs(*)
	integer mxCreateFull, mxGetPr
	integer y_pr, x_pr, pms_pr
	
	integer nlhs, nrhs
	integer mxGetM, mxGetN
	integer size_x, size_pms, m, n
	integer NMAX, NO_PARAMS, ORD_MAX
	
C-------------------------------------------------------------------------------
	PARAMETER (NMAX = 20000, NO_PARAMS = 2, ORD_MAX = 99)
C-------------------------------------------------------------------------------

	real*8  y(NMAX), x(NMAX), pms(NO_PARAMS,NMAX), c(ORD_MAX-1)
	integer ord, tot
	
	real*8  x1, temp
	integer i,k,q,p,so_pr
	
C-------------------------------------------------------------------------------
C	Error checking
C-------------------------------------------------------------------------------

	if (nrhs.LT.2) then
	   call mexErrMsgTxt('2 inputs required!')
	else if (nlhs.NE.1) then
	   call mexErrMsgTxt('1 output required!')
	endif

C-------------------------------------------------------------------------------
C	Sort out input arguments
C-------------------------------------------------------------------------------
C	Values of x
C-------------------------------------------------------------------------------

	m = mxGetM(prhs(1))
	n = mxGetN(prhs(1))
	size_x = m*n
	
	if ((m.ne.1) .AND. (n.ne.1)) then
	   call mexErrMsgTxt('Input is not a vector!')
	end if
	
	if (size_x.GT.NMAX) then
	   call mexErrMsgTxt('Vector length exceeds maximum!')
	end if
	
C-------------------------------------------------------------------------------
C	Length of vector x
C-------------------------------------------------------------------------------

	if ((m.EQ.1) .AND. (n.EQ.1)) then
	   tot = 1
	else
	   tot = max(m,n)
	endif
	
C-------------------------------------------------------------------------------
	
	x_pr = mxGetPr(prhs(1))
	call mxCopyPtrToReal8(x_pr,x,size_x)
	
C-------------------------------------------------------------------------------
C	Get pms
C-------------------------------------------------------------------------------

	m = mxGetM(prhs(2))
	n = mxGetN(prhs(2))
	size_pms = m*n
	
	if (size_pms.NE.1) then 
	   if (m.NE.NO_PARAMS) then
	      call mexErrMsgTxt('Parameter matrix orientation incorrect!')
	   endif
	
	   if ((size_pms.NE.(NO_PARAMS*size_x)) .AND. (n.NE.1))then
	      call mexErrMsgTxt('Incorrect number of parameters to function!')
	   end if
	
	   pms_pr = mxGetPr(prhs(2))
	   call mxCopyPtrToReal8(pms_pr,pms,size_pms)
	   
C-------------------------------------------------------------------------------
C	Get spline-order from MATLAB global variable, 'SPLINE_ORDER'.
C-------------------------------------------------------------------------------
	
	   so_pr = mexGetGlobal('SPLINE_ORDER')
	   if (so_pr.EQ.0) then
	      call mexErrMsgTxt('global variable spline-order not defined!')
	   else
	      call mxCopyPtrToReal8(mxGetPr(so_pr),temp,1)
	      ord =  int(temp)
	   endif
	else
	   pms_pr = mxGetPr(prhs(2))
	   call mxCopyPtrToReal8(pms_pr,temp,size_pms)
	   ord = int(temp)
	endif
	
	if ((ord.LT.1) .OR. (ord.GT.ORD_MAX)) then
	   call mexErrMsgTxt('Order exceeds maximum (= 99)!')
	endif
	
C-------------------------------------------------------------------------------

	do 10 i = 1,tot
	
	   y(i) = 0.0
	   
	   if ((size_pms.NE.1) .AND. (n.NE.1)) then
C-------------------------------------------------------------------------------
C	   Centre Functions
	      x(i) = x(i) + (pms(1,i)*ord)/2.0
C-------------------------------------------------------------------------------
	      x(i) = (x(i)-pms(2,i))/pms(1,i)
	   elseif ((size_pms.NE.1) .AND. (n.EQ.1)) then
C-------------------------------------------------------------------------------
C	   Centre Functions
	      x(i) = x(i) + (pms(1,1)*ord)/2.0
C-------------------------------------------------------------------------------
	      x(i) = (x(i)-pms(2,1))/pms(1,1)
	   else
C-------------------------------------------------------------------------------
C	   Centre Functions
	      x(i) = x(i) + ord/2.0
C-------------------------------------------------------------------------------
	   endif	
	         
C-------------------------------------------------------------------------------
C	   Characteristic function
C-------------------------------------------------------------------------------

	   if (ord.EQ.1) then
	      if ((x(i).GE.0.0) .AND. (x(i).LT.1.0)) then
	         y(i) = 1.0
	      else
	         y(i) = 0.0
	      endif
	   endif
	   
C-------------------------------------------------------------------------------
C	   Higher Order
C-------------------------------------------------------------------------------

	   if ((ord.GE.2) .AND. (ord.LE.ORD_MAX)) then
	   
	      do 20 k = 1,ord-1
	      
	         c(k) = 0.0
	         x1 = x(i) - k + 1
	         
	         if ((x1.GE.0.0) .AND. (x1.LT.1.0)) then
	            c(k) = x1
	         endif
	         
	         if ((x1.GE.1.0) .AND. (x1.LT.2.0)) then
	            c(k) = 2 - x1
	         endif
	         
20	      continue

	      do 30 p = 1,ord-2
	         do 40 q = 1,ord-1-p
	            c(q) = ((x(i)-q+1) * c(q)+(p+q+1-x(i))*c(q+1)) / (p+1)
40	         continue
30	      continue

	      if ((size_pms.NE.1) .AND. (n.NE.1)) then
	         y(i) = (1/sqrt(pms(1,i)))*c(1)
	      elseif ((size_pms.NE.1) .AND. (n.EQ.1)) then
	         y(i) = (1/sqrt(pms(1,1)))*c(1)
	      else
	         y(i) = c(1)
	      endif
	      
	   endif
10	continue

C-------------------------------------------------------------------------------
C	Arrange Output
C-------------------------------------------------------------------------------
	         
	plhs(1) = mxCreateFull(1,tot,0)
	y_pr = mxGetPr(plhs(1))
	call mxCopyReal8ToPtr(y,y_pr,tot)
	
	return
	end
	
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
