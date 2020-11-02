# ZPL for Nelde-Mead (Downhill-Simplex) Minimization algorithm
# Works with radius, thickness, and conic variables only
# Dr. Daniel P. Sablowski, dsablowski@aip.de
# October 2020
# 
# Define number of iterations:
NMIter = 1000
# Define minimal simplex sigma
sigma = 1e-13
# The algorithm will stop if one of the criteria is reached

# Load variables
nthic=0
nrad=0
nconi=0

GETVARDATA 1
nVar = VEC1(0)
PRINT "Number of variables: ", nVar

for i, 0, nVar-1, 1
	PRINT VEC1(5*i+1), "; ", VEC1(5*i+2), "; ", VEC1(5*i+3), "; ",  VEC1(5*i+4), "; ",  VEC1(5*i+5)
	IF(VEC1(5*i+1)==1)
		nrad=nrad+1
	ENDIF
	IF(VEC1(5*i+1)==2)
		nthic=nthic+1
	ENDIF
	IF(VEC1(5*i+1)==3)
		nconi=nconi+1
	ENDIF
	#PRINT 
NEXT
#PRINT nrad, "; ", nconi

# If there is no member of one of the allowed variable classes (radii, thicknesses and conic)
# we need to make sure that the dedicated vectors have at least a length of unity - otherwise
# we'll get an error message.
# After the declaration we need to change the counter of the accosiated variable class back to zero
#
# check if variable class is not present
rind=0
tind=0
cind=0
IF(nrad==0)		# if yes, no radius variables set
	rind=1
	nrad=1
ENDIF
IF(nthic==0)	# if yes, no thickness variables set
	tind=1
	nthic=1
ENDIF
IF(nconi==0)	# if yes, no conic variables set
	cind=1
	nconi=1
ENDIF

# declare vectors for each variable class with values of variables and the surface number
DECLARE RVar, DOUBLE, 1, nrad	# Vector with radii variables
DECLARE RSur, INT, 1, nrad		# Vector with radii variables surface number
DECLARE TVar, DOUBLE, 1, nthic	# Vector with thicknes variables
DECLARE TSur, INT, 1, nthic		# Vector with thicknes variables surface number
DECLARE CVar, DOUBLE, 1, nconi	# Vector with conic variables
DECLARE CSur, INT, 1, nconi		# Vector with conic variables surface number
rcount=1
tcount=1
ccount=1

# Change back the counter for the variable classes
IF(rind==1)
	nrad=0
ENDIF
IF(tind==1)
	nthic=0
ENDIF
IF(cind==1)
	nconi=0
ENDIF

# load variable values and surface numbers
FOR i, 0, nVar-1, 1
	IF(VEC1(5*i+1)==1)
		RVar(rcount)=VEC1(5*i+5)
		RSur(rcount)=VEC1(5*i+2)
		#PRINT RVar(rcount), "; ", RSur(rcount)
		rcount=rcount+1
	ENDIF
	IF(VEC1(5*i+1)==2)
		TVar(tcount)=VEC1(5*i+5)
		TSur(tcount)=VEC1(5*i+2)
		#PRINT TVar(tcount), "; ", TSur(tcount)
		tcount=tcount+1
	ENDIF
	IF(VEC1(5*i+1)==3)
		CVar(ccount)=VEC1(5*i+5)
		CSur(ccount)=VEC1(5*i+2)
		#PRINT CVar(ccount), "; ", CSur(ccount)
		ccount=ccount+1
	ENDIF
NEXT

# vector containing all variable values
DECLARE optval, DOUBLE, 1, nVar

FOR i, 1, nVar, 1
	IF(i<=nrad)
		optval(i)=RVar(i)
	else
		if(i<=nrad+nthic)
			optval(i)=TVar(i-nrad)
		else
			optval(i)=CVar(i-nrad-nthic)
		endif
	ENDIF
	PRINT optval(i)
NEXT
PRINT

# Simplex tranformation coefficients:
gamma = 2.0
alpha = 1.0
beta = 0.5
btot = 0.5

DECLARE y, DOUBLE, 1, nVar+1
DECLARE Z, Double, 1, nVar
DECLARE C, Double, 1, nVar
DECLARE S, Double, 1, nVar
DECLARE Em, Double, 1, nVar
DECLARE step, Double, 1, nVar
DECLARE Pm, Double, 2, nVar+1, nVar
yh=0
yl=0
ysh=0
ym=0
ys=0
yi=0
yt=0
Ph=0
Pl=0
Psh=0
Iters=0


# Create initial simplex points
# The step size defines the "size" of simplex within the initial space
FOR i, 1, nVar, 1
	Pm(1,i) = optval(i)
	if(optval(i)==0)
		step(i)=0.01
		# or use randomized steps
		#step(i)=RAND(0.01)*SIGN(RAND(1)-0.5)
	else
		step(i)=optval(i)*0.1
		# or use randomized steps
		#step(i) = optval(i)*RAND(0.1)*SIGN(RAND(1)-0.5)
	endif
	#PRINT "Pm(1,",i,"): ", Pm(1,i), "; step: ", step(i)
NEXT
PRINT
FOR i, 1, nVar+1, 1
	FOR j, 1, nVar, 1
		if((i>1) & (i==j+1))
			aa = 1.0
		else
			aa = 0.0
		endif
		if(i==1)
			optval(j)=Pm(1,j)
		else
			Pm(i,j)=Pm(1,j)+step(j)*aa
		endif
		optval(j)=Pm(i,j)
		if(j<=nrad)
			SURP RSur(j), 2, optval(j)
		else
			if(j<=nrad+nthic)
				SURP TSur(j-nrad), 3, optval(j)
			else
				SURP CSur(j-nrad-nthic), 5, optval(j)
			endif
		endif
	NEXT
	y(i)=MFCN()
	#PRINT y(i)
NEXT

# Start main loop
FOR tc, 1, NMIter, 1
	#initialize next interation
	ym=0
	ys=0
	FOR i, 1, nVar, 1
		Z(i)=0.0
	NEXT
	# looking for worst point
	yh=y(1)
	Ph=1
	FOR j, 1, nVar+1, 1
		if(y(j)>=yh)
			yh=y(j)
			Ph = j
		ENDIF
	NEXT
	# looking for best point
	yl=yh
	Pl=Ph
	FOR j, 1, nVar+1, 1
		if(y(j)<yl)
			yl=y(j)
			Pl=j
		endif
	NEXT
	# looking for second highest value
	ysh=yl
	Psh=Pl
	FOR j, 1, nVar+1, 1
		if((y(j)>ysh) & (y(j)<yh) & (y(j)>yl) & (j!=Pl))
			ysh=y(j)
			Psh=j
		endif
	NEXT
	#PRINT Pl, "; ", Psh, "; ", Ph
	# compute mean and sigma
	FOR i, 1, nVar+1, 1
		ym=ym+y(i)/(nVar+1)
	NEXT
	FOR i, 1, nVar+1, 1
		ys=ys+SQRT((y(i)-ym)*(y(i)-ym))
	NEXT
	ys = ys/nVar
	#compute centroid point
	FOR j, 1, nVar, 1
		Z(j)=0
		FOR i, 1, nVar+1, 1
			if(i!=Ph)
				Z(j)=Z(j)+Pm(i,j)/nVar
			endif
		NEXT
	NEXT
	# reflect highest value at centroid
	FOR i, 1, nVar, 1
		C(i)=Z(i)+alpha*(Z(i)-Pm(Ph,i))
		optval(i)=C(i)
		if(i<=nrad)
			SURP RSur(i), 2, optval(i)
		else
			if(i<=nrad+nthic)
				SURP TSur(i-nrad), 3, optval(i)
			else
				SURP CSur(i-nrad-nthic), 5, optval(i)
			endif
		endif
	NEXT
	yi = MFCN()
	# Alpha Branch
	if(yi<yl)
	#PRINT "exp"
		FOR i, 1, nVar, 1
			Em(i)=Z(i)+gamma*(C(i)-Z(i))
			optval(i)=Em(i)
			if(i<=nrad)
				SURP RSur(i), 2, optval(i)
			else
				if(i<=nrad+nthic)
					SURP TSur(i-nrad), 3, optval(i)
				else
					SURP CSur(i-nrad-nthic), 5, optval(i)
				endif
			endif
		NEXT
		yt=MFCN()
		if(yt<yl)
			FOR i, 1, nVar, 1
				Pm(Ph,i)=Em(i)
			NEXT
			y(Ph)=yt
		else
			FOR i, 1, nVar, 1
				Pm(Ph,i)=C(i)
			NEXT
			y(Ph)=yi
		endif
	else	# end alpha / start beta branch
		if(yi<=ysh)
			FOR i, 1, nVar, 1
				Pm(Ph,i) = C(i)
			NEXT
			y(Ph)=yi
		endif
		if(yi>ysh)	# beta a
			if(yi<=yh)
				FOR i, 1, nVar, 1
					Pm(Ph,i) = C(i)
				NEXT
				y(Ph)=yi
				yh=y(Ph)
			endif
			#PRINT "single contr."
			FOR i, 1, nVar, 1,
				S(i)=Z(i)+beta*(Pm(Ph,i)-Z(i))
				optval(i)=S(i)
				if(i<=nrad)
					SURP RSur(i), 2, optval(i)
				else 
					if(i<=nrad+nthic)
						SURP TSur(i-nrad), 3, optval(i)
					else
						SURP CSur(i-nrad-nthic), 5, optval(i)
					endif
				endif
			NEXT
			yt=MFCN()
			if(yt>yh)	# beta b
				#PRINT "total contr."
				FOR j, 1, nVar+1, 1
					FOR i, 1, nVar, 1
						Pm(j,i)=Pm(Pl,i)+btot*(Pm(j,i)-Pm(Pl,i))
						optval(i)=Pm(j,i)
						if(i<=nrad)
							SURP RSur(i), 2, optval(i)
						else
							if(i<=nrad+nthic)
								SURP TSur(i-nrad), 3, optval(i)
							else
								SURP CSur(i-nrad-nthic), 5, optval(i)
							endif
						endif
					NEXT
					y(j)=MFCN()
				NEXT
			endif	# end beta b
			if(yt<=yh)
				FOR i, 1, nVar, 1
					Pm(Ph,i)=S(i)
				NEXT
				y(Ph)=yt
			endif
		endif	# end beta a
	endif	# end beta Branch
	#PRINT yl
	# compute mean and sigma
	ym=0
	ys=0
	FOR i, 1, nVar+1, 1
		ym=ym+y(i)/(nVar+1)
	NEXT
	FOR i, 1, nVar+1, 1
		ys=ys+SQRT((y(i)-ym)*(y(i)-ym))
	NEXT
	ys = ys/nVar
	Iters=tc
	if(ys<=sigma)
		tc=NMiter+1
	endif
NEXT	# end main loop

# looking for smallest value
yl=y(1)
FOR j, 1, nVar+1, 1
	if(y(j)<yl)
		yl=y(j)
		Pl = j
	endif
NEXT
FOR i, 1, nVar, 1
	optval(i)=Pm(Pl,i)
	if(i<=nrad)
		SURP RSur(i), 2, optval(i)
	else
		if(i<=nrad+nthic)
			SURP TSur(i-nrad), 3, optval(i)
		else
			SURP CSur(i-nrad-nthic), 5, optval(i)
		endif
	endif
NEXT
PRINT "Converged after ", Iters, " number of iterations."

UPDATE ALL
