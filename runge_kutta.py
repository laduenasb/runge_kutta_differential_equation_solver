def Runge_kutta_sm(U,w,xa,xb,de,h): #### Codigo del metodo Runge Kutta Felhberg para solucionar sistemas de ecuaciones diferenciales ordinarias numericamente con U el vector de funciones igualadas du2(t)/dt,du3(t)/dt,...,dun(t)/dt y el vector w con las condiciones iniciales w2(xa),w3(xa),...,wn(xa) en el intervalo xa,xb con "de" la cantidad de decimales exactos, donde la matriz es mostrada verticalmente donde de izquierda a derecha en forma vertical se muestra las soluciones u1,u2(t),...,un(t), con u1 el tiempo, siempre ubicar u1 como la variable temporal y h el deltax, advertencia: poner las derivadas en orden en el vector U, es decir que las derivadas queden en orden creciente, al ubicarlas en el vector U.  
	m=len(U) ###numero de ecuaciones
	def definirn1(U):
		a="u1,"
		for j in range(1,len(U)):
			a+="u"+str(j+1)+","
		a+="u"+str(len(U)+1)
		for i in range(0,len(U)):
			print("        def funcion"+str(i+1)+"("+a+"):")
			print("                f="+str(U[i]))
			print("                return f")
		return ""
	print("from definiciones import crear_vector,crear_matriz,vertical2,minimo,Runge_kutta_sm,equalv")
	print("from math import *")
	print("xa="+str(xa))
	print( "xb="+str(xb))
	print("de="+str(de))
	print("h="+str(h))
	print("w="+str(w))
	print("def Runge_kutta_s2(w,xa,xb,de,h):")
	print("        m="+str(m))
	print("        xa=float(xa)")
	print("        xb=float(xb)")
	print("        npuntos=int(ceil((xb-xa)/h+1.0))")
	print("        vector_t=crear_vector(xa,xb,npuntos)")
	print("        h=vector_t[1]-vector_t[0]")
	print("        W1=crear_matriz(len(vector_t),m)")
	print("        w1=[]")
	print("        w2=[]")
	print("        w9=[]")
	print("        h_new=[]")
	print("        error=10**(-de-1)")
	print("        for i in range(0,len(w)):")
	print("                w1.append(float(w[i]))")
	print("                w2.append(float(w[i]))")
	print("                w9.append(float(w[i]))")
	print("                W1[0][i]=float(w[i])")
	print(definirn1(U))
	print('        for i in range(1,len(vector_t)):')
	######### k1
	a="		"
	b="vector_t[i-1],"
	for j in range(1,int(m)):
		b+="w1["+str(j-1)+"]"+","
	b+="w1["+str(m-1)+"]"
	for i in range(0,int(m)):
		y1=a+"k1"+str(i+1)+"="+"h*funcion"+str(i+1)+"("+b+")"
		print(y1)
	######### k2
	c="		"
	d="vector_t[i-1]+(h/4.0),"
	for j in range(1,int(m)):
		d+="w1["+str(j-1)+"]"+"+(1.0/4.0)*k1"+str(j)+","
	d+="w1["+str(m-1)+"]"+"+(1.0/4.0)*k1"+str(m)
	for i in range(0,int(m)):
		y2=c+"k2"+str(i+1)+"="+"h*funcion"+str(i+1)+"("+d+")"
		print(y2)
	######### k3
	e="		"
	f="vector_t[i-1]+(3.0*h/8.0),"
	for j in range(1,int(m)):
		f+="w1["+str(j-1)+"]"+"+(3.0/32.0)*k1"+str(j)+"+(9.0/32.0)*k2"+str(j)+","
	f+="w1["+str(m-1)+"]"+"+(3.0/32.0)*k1"+str(m)+"+(9.0/32.0)*k2"+str(m)
	for i in range(0,int(m)):
		y3=e+"k3"+str(i+1)+"="+"h*funcion"+str(i+1)+"("+f+")"
		print(y3)
	######## k4
	g="		"
	h="vector_t[i-1]+(12.0*h/13.0),"
	for j in range(1,int(m)):
		h+="w1["+str(j-1)+"]"+"+(1932.0/2197.0)*k1"+str(j)+"-(7200.0/2197.0)*k2"+str(j)+"+(7296.0/2197.0)*k3"+str(j)+","
	h+="w1["+str(m-1)+"]"+"+(1932.0/2197.0)*k1"+str(m)+"-(7200.0/2197.0)*k2"+str(m)+"+(7296.0/2197.0)*k3"+str(m)
	for i in range(0,int(m)):
		y4=g+"k4"+str(i+1)+"="+"h*funcion"+str(i+1)+"("+h+")"
		print(y4)
	####### k5
	k="		"
	l="vector_t[i-1]+h,"
	for j in range(1,int(m)):
		l+="w1["+str(j-1)+"]"+"+(439.0/216.0)*k1"+str(j)+"-8.0*k2"+str(j)+"+(3680.0/513.0)*k3"+str(j)+"-(845.0/4104.0)*k4"+str(j)+","
	l+="w1["+str(m-1)+"]"+"+(439.0/216.0)*k1"+str(m)+"-8.0*k2"+str(m)+"+(3680.0/513.0)*k3"+str(m)+"-(845.0/4104.0)*k4"+str(m)
	for i in range(0,int(m)):
		y5=k+"k5"+str(i+1)+"="+"h*funcion"+str(i+1)+"("+l+")"
		print(y5)
	####### k6
	q="		"
	w="vector_t[i-1]+(h/2.0),"
	for j in range(1,int(m)):
		w+="w1["+str(j-1)+"]"+"-(8.0/27.0)*k1"+str(j)+"+2.0*k2"+str(j)+"-(3544.0/2565.0)*k3"+str(j)+"+(1859.0/4104.0)*k4"+str(j)+"-(11.0/40.0)*k5"+str(j)+","
	w+="w1["+str(m-1)+"]"+"-(8.0/27.0)*k1"+str(m)+"+2.0*k2"+str(m)+"-(3544.0/2565.0)*k3"+str(m)+"+(1859.0/4104.0)*k4"+str(m)+"-(11.0/40.0)*k5"+str(m)
	for i in range(0,int(m)):
		y6=k+"k6"+str(i+1)+"="+"h*funcion"+str(i+1)+"("+w+")"
		print(y6)
	####### w
	for j in range(0,int(m)):
		y7="		w1["+str(j)+"]"+"="+"w1["+str(j)+"]"+"+(25.0/216.0)*k1"+str(j+1)+"+(1408.0/2565.0)*k3"+str(j+1)+"+(2197.0/4104.0)*k4"+str(j+1)+"-(1.0/5.0)*k5"+str(j+1)	
		print(y7)
		print("		W1[i]["+str(j)+"]"+"="+"w1["+str(j)+"]")
	for j in range(0,int(m)):
		y8="		w2["+str(j)+"]"+"="+"w2["+str(j)+"]"+"+(16.0/135.0)*k1"+str(j+1)+"+(6656.0/12825.0)*k3"+str(j+1)+"+(28561.0/56430.0)*k4"+str(j+1)+"-(9.0/50.0)*k5"+str(j+1)+"+(2.0/55.0)*k6"+str(j+1)
		print(y8)
		print("                if abs("+"w1["+str(j)+"]-"+"w2["+str(j)+"])!=0.0:")
		print("                        lo=(0.5*(h**5)*error/abs("+"w1["+str(j)+"]-"+"w2["+str(j)+"]"+"))**(1.0/4.0)")
		print("                        h_new.append(lo)")
	print("        equalv(w9,w)")
	print('        hkl=""')
	print("        for i in range(0,m+1):")
	print("                d=str(i+1)")
	print('                g="u"')
	print('                ghj=" "')
	print("                hkl+=g+d+ghj")
	print("        print hkl")
	print('        h=vector_t[1]-vector_t[0]')
	print('        if len(h_new)==0:')
	print("                h_new.append(h)")
	print("        return vector_t,W1,h,minimo(h_new)")
	print("F=Runge_kutta_s2(w,xa,xb,de,h)")
	print("vt=F[0]")
	print("H=F[1]")
	print("hn=F[2]")
	print("hmin=F[3]")
	print("print vertical2(vt,H)")
	print("print 'h=',hn")
	print("print 'hmin=',hmin")
	return ""
