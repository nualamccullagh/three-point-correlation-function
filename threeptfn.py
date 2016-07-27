import numpy as N
import scipy.special as S
import scipy.interpolate as I


global xi00i
global xi0n2i
global xi02i
global xi11i
global xi31i
global xi51i
global xi1n1i
global xi3n1i
global xi5n1i
global xi20i
global xi2n2i
global xi22i
global xi40i


def setup(pklinfile):
    
    
    pklin=N.loadtxt(pklinfile)
    k=pklin[:,0]
    pk=pklin[:,1]
    pklinI=I.interp1d(k, pk)
    
    global xi00i
    global xi0n2i
    global xi02i
    global xi11i
    global xi31i
    global xi51i
    global xi1n1i
    global xi3n1i
    global xi5n1i
    global xi20i
    global xi2n2i
    global xi22i
    global xi40i
    rxi=N.arange(0.0, 300, 0.5)
    xx00=xinm(0, 0, rxi, pklinI)
    xx0n2=xinm(0, -2, rxi, pklinI)
    xx02=xinm(0, 2, rxi, pklinI)
    xx11=xinm(1, 1, rxi, pklinI)
    xx31=xinm(3, 1, rxi, pklinI)
    xx51=xinm(5, 1, rxi, pklinI)
    xx1n1=xinm(1, -1, rxi, pklinI)
    xx3n1=xinm(3, -1, rxi, pklinI)
    xx5n1=xinm(5, -1, rxi, pklinI)
    xx20=xinm(2, 0, rxi, pklinI)
    xx2n2=xinm(2, -2, rxi, pklinI)
    xx22=xinm(2, 2, rxi, pklinI)
    xx40=xinm(4, 0, rxi, pklinI)
    xi00i=I.interp1d(rxi, xx00, kind='cubic')
    xi02i=I.interp1d(rxi, xx02, kind='cubic')
    xi0n2i=I.interp1d(rxi, xx0n2, kind='cubic')
    xi20i=I.interp1d(rxi, xx20, kind='cubic')
    xi2n2i=I.interp1d(rxi, xx2n2, kind='cubic')
    xi22i=I.interp1d(rxi, xx22, kind='cubic')
    xi40i=I.interp1d(rxi, xx40, kind='cubic')
    xi11i=I.interp1d(rxi, xx11, kind='cubic')
    xi31i=I.interp1d(rxi, xx31, kind='cubic')
    xi51i=I.interp1d(rxi, xx51, kind='cubic')
    xi1n1i=I.interp1d(rxi, xx1n1, kind='cubic')
    xi3n1i=I.interp1d(rxi, xx3n1, kind='cubic')
    xi5n1i=I.interp1d(rxi, xx5n1, kind='cubic')


        

def zeta_eq_rs(rr, theta, phi, b1, b2, bs, f):
    global xi00i
    global xi11i
    global xi31i
    global xi51i
    global xi1n1i
    global xi3n1i
    global xi5n1i
    global xi20i
    global xi40i
    global rxi
    
    rr=N.asarray(rr)
    theta=N.asarray(theta)
    phi=N.asarray(phi)
    
    t1=phi
    
    xx=N.cos(N.pi/3)
    
    xi00=xi00i(rr)
    xi11=xi11i(rr)
    xi31=xi31i(rr)
    xi51=xi51i(rr)
    xi1n1=xi1n1i(rr)
    xi3n1=xi3n1i(rr)
    xi5n1=xi5n1i(rr)
    xi20=xi20i(rr)
    xi40=xi40i(rr)
    
    if (t1.shape[0]>1):
        xi00=N.tile(xi00, (t1.shape[0], 1))
        xi20=N.tile(xi20, (t1.shape[0], 1))
        xi40=N.tile(xi40, (t1.shape[0], 1))
        xi11=N.tile(xi11, (t1.shape[0], 1))
        xi31=N.tile(xi31, (t1.shape[0], 1))
        xi51=N.tile(xi51, (t1.shape[0], 1))
        xi1n1=N.tile(xi1n1, (t1.shape[0], 1))
        xi3n1=N.tile(xi3n1, (t1.shape[0], 1))
        xi5n1=N.tile(xi5n1, (t1.shape[0], 1))
    
    
    def do_calc(i, j):
        t00=34./21.+b2/b1+(8.*f)/3.-(44.*f)/(63.*b1)+(2.*b2*f)/(3.*b1)-(4.*bs*f)/(9.*b1**2)+(4.*bs*f)/(9.*b1)+(22.*f**2)/15.-(88.*f**2)/(225.*b1**2)+(2.*f**2)/(9.*b1)+(b2*f**2)/(9.*b1)-(14.*bs*f**2)/(225.*b1**3)+(2.*bs*f**2)/(27.*b1)+(4.*f**3)/15.-(4.*f**3)/(75.*b1**3)+(4.*f**3)/(15.*b1)+(2.*f**4)/(25.*b1)
        t02=(1/(105840.*b1**3))*f*(-21.*f*(16.*bs+15.*f)-6.*b1*(196.*bs+235.*f)+126.*b1**3*(70.+65.*f+17.*f**2)+7.*b1**2*(105.*b2*(3.+f)+70.*bs*(3.+f)+3.*(-67.+70.*f+102.*f**2+36.*f**3)))*(8.-6.*N.cos(2.*t1)-3.*N.cos(2.*(t1-theta))-24.*N.cos(2.*theta)-3.*N.cos(2*(t1+theta))+6.*N.sqrt(3.)*N.sin(2.*t1)+3.*N.sqrt(3) *N.sin(2.*(t1-theta))+3.*N.sqrt(3)*N.sin(2.*(t1+theta)))
        t04=(1./(940800.*b1**3.))*f**2.*(-27.*b1+70.*b1**3.*(3.+f)+7.*(2.*bs+f)+14.*b1**2.*f*(5.+3.*f))*(216.-120.*N.cos(2.*t1)+210.*N.cos(4.*t1)+140.*N.cos(2.*(t1-2.*theta))+140.*N.cos(4.*t1-2.*theta)+80.*N.cos(2.*(t1-theta))+35.*N.cos(4.*(t1-theta))-480.*N.cos(2.*theta)+840.*N.cos(4.*theta)+80.*N.cos(2.*(t1+theta))+35.*N.cos(4.*(t1+theta))+140.*N.cos(2.*(2.*t1+theta))+140.*N.cos(2.*(t1+2.*theta))+120.*N.sqrt(3.)*N.sin(2.*t1)+210.*N.sqrt(3.)*N.sin(4.*t1)-140.*N.sqrt(3.)*N.sin(2.*(t1-2.*theta))+140.*N.sqrt(3.)*N.sin(4.*t1-2.*theta)-80.*N.sqrt(3.)*N.sin(2.*(t1-theta))+35.*N.sqrt(3.)*N.sin(4.*(t1-theta))-80.*N.sqrt(3.)*N.sin(2.*(t1+theta))+35.*N.sqrt(3.)*N.sin(4.*(t1+theta))+140.*N.sqrt(3.)*N.sin(2.*(2.*t1+theta))-140.*N.sqrt(3.)*N.sin(2.*(t1+2.*theta)))
        t22=(1./(37044.*b1**3.))*(-21.*f**2.*(59.*bs+48.*f)+6.*b1*f*(637.*bs+862.*f)+1764.*b1**3.*(-1.+14.*f**2.+8.*f**3.)+7.*b1**2.*(49.*bs*(-9.+8.*f**2.)+12.*f*(124.+49.*(2.+b2)*f+168.*f**2.+72.*f**3.))-9.*f**2.*(-654.*b1-630.*(bs+f)+882.*b1**3.*(7.+4.*f)+7.*b1**2.*(147.*b2+98.*bs+6.*(7.+6.*f)**2.))*N.cos(t1)**2.*N.cos(theta)**4.*(-2.+N.cos(2.*t1)+N.sqrt(3.)*N.sin(2.*t1))-3.*f*N.cos(theta)**2.*(6825.*b1**2.+2940.*b1*bs+2274.*b1*f+8232.*b1**2.*f+24696.*b1**3.*f+4116.*b1**2.*b2*f-1680.*bs*f+2744.*b1**2.*bs*f-1575.*f**2.+14112.*b1**2.*f**2.+14112.*b1**3.*f**2.+6048.*b1**2.*f**3.+(882.*b1**3.*f*(7.+4.*f)-21.*f*(38.*bs+39.*f)-12.*b1*(49.*bs+136.*f)+7.*b1**2.*(-195.+49.*(6.+3.*b2+2.*bs)*f+504.*f**2.+216.*f**3.))*N.cos(2.*t1)-N.sqrt(3.)*(882.*b1**3.*f*(7.+4.*f)-21.*f*(38.*bs+39.*f)-12.*b1*(49.*bs+136.*f)+7.*b1**2.*(-195.+49.*(6.+3.*b2+2.*bs)*f+504.*f**2.+216.*f**3.))*N.sin(2.*t1)))
        t24=(1./(82320.*b1**3.))*f*(24.*(392.*b1**3.*f**2.+7.*f*(17.*bs+20.*f)+b1*(49.*bs+460.*f)+28.*b1**2.*(1.+14.*f**2.+12.*f**3.))+735.*f**2.*(-1.+14.*b1**3.+2.*b1**2.*(7.+6.*f))*N.cos(t1)**2.*N.cos(theta)**6.*(-9.+2.*N.cos(2.*t1)+2.*N.cos(4.*t1)+6.*N.sqrt(3.)*N.sin(2.*t1))-3.*N.cos(theta)**2.*(4.*(5096.*b1**3.*f**2.+7.*f*(151.*bs+170.*f)+3.*b1*(49.*bs+1410.*f)+28.*b1**2.*(3.+182.*f**2.+156.*f**3.))+(5096.*b1**3.*f**2.-57.*b1*(49.*bs+20.*f)-7.*f*(89.*bs+100.*f)+28.*b1**2.*(-57.+182.*f**2.+156.*f**3.))*N.cos(2.*t1)-N.sqrt(3.)*(5096.*b1**3.*f**2.-57.*b1*(49.*bs+20.*f)-7.*f*(89.*bs+100.*f)+28.*b1**2.*(-57.+182.*f**2.+156.*f**3.))*N.sin(2.*t1))-5.*f*N.cos(theta)**4.*(-8550.*b1-2100.*bs-1806.*f-18816.*b1**2.*f-18816.*b1**3.*f-16128.*b1**2.*f**2.-(1311.*b1+7.*(46.*bs-31.*f)+9800.*b1**3.*f+1400.*b1**2.*f*(7.+6.*f))*N.cos(2.*t1)+(57.*b1+14.*(bs-13.*f)+2842.*b1**3.*f+406.*b1**2.*f*(7.+6.*f))*N.cos(4.*t1)+1311.*N.sqrt(3.)*b1*N.sin(2.*t1)+322.*N.sqrt(3.)*bs*N.sin(2.*t1)-217.*N.sqrt(3.)*f*N.sin(2.*t1)+9800.*N.sqrt(3.)*b1**2.*f*N.sin(2.*t1)+9800.*N.sqrt(3.)*b1**3.*f*N.sin(2.*t1)+8400.*N.sqrt(3.)*b1**2.*f**2.*N.sin(2.*t1)+57.*N.sqrt(3.)*b1*N.sin(4.*t1)+14.*N.sqrt(3.)*bs*N.sin(4.*t1)-182.*N.sqrt(3.)*f*N.sin(4.*t1)+2842.*N.sqrt(3.)*b1**2.*f*N.sin(4.*t1)+2842.*N.sqrt(3.)*b1**3.*f*N.sin(4.*t1)+2436.*N.sqrt(3.)*b1**2.*f**2.*N.sin(4.*t1)))
        
        t44=(1./(137200.*b1**3.))*f**2.*(4.*(116.*b1+203.*bs+504.*f+504.*b1**2.*f**2.)+5.*N.cos(theta)**2.*(-4.*(92.*b1+161.*bs+693.*f+1008.*b1**2.*f**2.)+(388.*b1+679.*bs+252.*f-1008.*b1**2.*f**2.)*N.cos(2.*t1)+N.sqrt(3.)*(-388.*b1-679.*bs-252.*f+1008.*b1**2.*f**2.)*N.sin(2.*t1))+1225.*f*N.cos(t1)**2.*N.cos(theta)**6.*((19.+24.*b1**2.*f)*N.cos(2.*t1)+(-2.+24.*b1**2.*f)*N.cos(4.*t1)+3.*(-11.-36.*b1**2.*f+N.sqrt(3.)*(5.+24.*b1**2.*f)*N.sin(2.*t1)))-17150.*b1**2.*f**2.*N.cos(t1)**4.*N.cos(theta)**8.*(-6.+4.*N.cos(2.*t1)+N.cos(4.*t1)+4.*N.sqrt(3.)*N.sin(2.*t1)-N.sqrt(3.)*N.sin(4.*t1))-15./4.*N.cos(theta)**4.*(-36.*b1-63.*bs-6384.*f-14784.*b1**2.*f**2.-2.*(-324.*b1-567.*bs+854.*f+4144.*b1**2.*f**2.)*N.cos(2.*t1)+(12.*b1+21.*bs+1148.*f+2968.*b1**2.*f**2.)*N.cos(4.*t1)-648.*N.sqrt(3.)*b1*N.sin(2.*t1)-1134.*N.sqrt(3.)*bs*N.sin(2.*t1)+1708.*N.sqrt(3.)*f*N.sin(2.*t1)+8288.*N.sqrt(3.)*b1**2.*f**2.*N.sin(2.*t1)+12.*N.sqrt(3.)*b1*N.sin(4.*t1)+21.*N.sqrt(3.)*bs*N.sin(4.*t1)+1148.*N.sqrt(3.)*f*N.sin(4.*t1)+2968.*N.sqrt(3.)*b1**2.*f**2.*N.sin(4.*t1)))
        t11=-(1./(175.*b1**4.))*(b1*(175.*b1**3.+105.*b1**2.*f+29.*b1*f**2.+3.*f**3.)+f*(385.*b1**3.+487.*b1**2.*f+243.*b1*f**2.+45.*f**3.)*N.cos(t1)**2.*N.cos(theta)**2.-N.sqrt(3.)*f*(385.*b1**3.+487.*b1**2.*f+243.*b1*f**2.+45.*f**3.)*N.cos(t1)*N.cos(theta)**2.*N.sin(t1))
        t31=-(1./(1400.*b1**4.))*f*(8.*b1*(35.*b1**2.+14.*b1*f+3.*f**2.)+3.*N.cos(theta)**2.*(-5.*b1*(35.*b1**2.+14.*b1*f+3.*f**2.)+2.*(35.*b1**3.+266.*b1**2.*f+219.*b1*f**2.+60.*f**3.)*N.cos(t1)**2.-70.*N.sqrt(3.)*b1**3.*N.cos(t1)*N.sin(t1)-N.sqrt(3.)*f*(266.*b1**2.+219.*b1*f+60.*f**2.)*N.sin(2.*t1))+15.*f*(21.*b1**2.+18.*b1*f+5.*f**2.)*N.cos(t1)*N.cos(theta)**4.*(-6.*N.cos(t1)+N.cos(3.*t1)+N.sqrt(3.)*(4.*N.sin(t1)+N.sin(3.*t1))))
        
        t13=-(1./(900.*b1**4.))*f*(24.*b1*(15.*b1**2.+8.*b1*f+f**2.)+3.*N.cos(theta)**2.*(-15.*b1*(15.*b1**2.+8.*b1*f+f**2.)+2.*(45.*b1**3.+316.*b1**2.*f+259.*b1*f**2.+60.*f**3.)*N.cos(t1)**2.-90.*N.sqrt(3.)*b1**3.*N.cos(t1)*N.sin(t1)-N.sqrt(3.)*f*(316.*b1**2.+259.*b1*f+60.*f**2.)*N.sin(2.*t1))+5.*f*(73.*b1**2.+64.*b1*f+15.*f**2.)*N.cos(t1)*N.cos(theta)**4.*(-6.*N.cos(t1)+N.cos(3.*t1)+N.sqrt(3.)*(4.*N.sin(t1)+N.sin(3.*t1))))
        t33=(1./(900.*b1**4.))*f**2.*(-24.*b1*(3.*b1+f)+50.*f*(11.*b1+5.*f)*N.cos(t1)**3.*N.cos(theta)**6.*(-3.*N.cos(t1)+2.*N.cos(3.*t1)+3.*N.sqrt(3.)*N.sin(t1))+15.*N.cos(t1)*N.cos(theta)**4.*(3.*(9.*b1**2.+47.*b1*f+20.*f**2.)*N.cos(t1)+(3.*b1**2.-21.*b1*f-10.*f**2.)*N.cos(3.*t1)+2.*N.sqrt(3.)*(-5.*(3.*b1**2.+12.*b1*f+5.*f**2.)+(3.*b1**2.-21.*b1*f-10.*f**2.)*N.cos(2.*t1))*N.sin(t1))-6.*N.cos(theta)**2.*(-15.*b1*(3.*b1+f)+2.*(33.*b1**2.+77.*b1*f+30.*f**2.)*N.cos(t1)**2.-N.sqrt(3.)*(33.*b1**2.+77.*b1*f+30.*f**2.)*N.sin(2.*t1)))
        t15=-(1./(10080.*b1**4.))*f**2.*(96.*b1*(5.*b1+f)+24.*N.cos(theta)**2.*(-21.*b1*(5.*b1+f)+2.*(5.*b1**2.+71.*b1*f+30.*f**2.)*N.cos(t1)**2.-N.sqrt(3.)*(5.*b1**2.+71.*b1*f+30.*f**2.)*N.sin(2.*t1))+63.*f*(7.*b1+3.*f)*N.cos(t1)*N.cos(theta)**6.*(20.*N.cos(t1)-5.*N.cos(3.*t1)+2.*N.cos(5.*t1)-12.*N.sqrt(3.)*N.sin(t1)-3.*N.sqrt(3.)*N.sin(3.*t1))+7.*N.cos(theta)**4.*(300.*b1**2.-780.*b1*f-360.*f**2.-4.*(5.*b1**2.+176.*b1*f+75.*f**2.)*N.cos(2.*t1)+(85.*b1**2.+157.*b1*f+60.*f**2.)*N.cos(4.*t1)+20.*N.sqrt(3.)*b1**2.*N.sin(2.*t1)+704.*N.sqrt(3.)*b1*f*N.sin(2.*t1)+300.*N.sqrt(3.)*f**2.*N.sin(2.*t1)+85.*N.sqrt(3.)*b1**2.*N.sin(4.*t1)+157.*N.sqrt(3.)*b1*f*N.sin(4.*t1)+60.*N.sqrt(3.)*f**2.*N.sin(4.*t1)))
        t35=-(1./(20160.*b1**4.))*f**3.*(192.*b1+24.*N.cos(theta)**2.*(-57.*b1+2.*(37.*b1+60.*f)*N.cos(t1)**2.-N.sqrt(3.)*(37.*b1+60.*f)*N.sin(2.*t1))+14.*N.cos(t1)*N.cos(theta)**6.*(40.*(7.*b1+31.*f)*N.cos(t1)-5.*(8.*b1+47.*f)*N.cos(3.*t1)+b1*N.cos(5.*t1)-146.*f*N.cos(5.*t1)-201.*N.sqrt(3.)*b1*N.sin(t1)-624.*N.sqrt(3.)*f*N.sin(t1)-39.*N.sqrt(3.)*b1*N.sin(3.*t1)-381.*N.sqrt(3.)*f*N.sin(3.*t1))+6.*N.cos(theta)**4.*(-240.*b1-1140.*f-2.*(198.*b1+475.*f)*N.cos(2.*t1)+(33.*b1+190.*f)*N.cos(4.*t1)+396.*N.sqrt(3.)*b1*N.sin(2.*t1)+950.*N.sqrt(3.)*f*N.sin(2.*t1)+33.*N.sqrt(3.)*b1*N.sin(4.*t1)+190.*N.sqrt(3.)*f*N.sin(4.*t1))-630.*f*N.cos(t1)**3.*N.cos(theta)**8.*(17.*N.cos(t1)-11.*N.cos(3.*t1)-N.cos(5.*t1)-13.*N.sqrt(3.)*N.sin(t1)-3.*N.sqrt(3.)*N.sin(3.*t1)+N.sqrt(3.)*N.sin(5.*t1)))
        return t00*xi00**2+t02*xi00*xi20+t04*xi00*xi40+t22*xi20**2+t24*xi20*xi40+t44*xi40**2+t11*xi11*xi1n1+t31*xi3n1*xi11+t13*xi31*xi1n1+t33*xi3n1*xi31+t15*xi1n1*xi51+t35*xi3n1*xi51
    return (do_calc(0, 1)+do_calc(1, 2)+do_calc(0, 2))*b1**3
        
        





def besph(n, z):
	ans=S.jn(n+0.5, z)*N.sqrt(N.pi/(2*z))
	if n==0:
		ans[N.where(z==0)]=1
	else:
	    ans[N.where(z==0)]=0
	return ans


def xinm(n, m, r, pklin):
    L=20000
    maxk=4*N.pi
    deltak=maxk/L
    #kk=N.arange(deltak, deltak*L, deltak)
    upper=maxk
    deltak=upper/L
    kk=N.arange(10**-3, deltak*L, deltak)
    sk=r.shape
    jrk=besph(n, N.outer(r,kk))
    rsm=0.5
    
	
	#gaussian smoothing at 0.1*r
    #w=N.exp(-N.outer(r,kk)**2*0.1**2/2)
    w=N.exp(-(kk*rsm)**2/2.0)
    #w=N.sin(kk*N.pi/(2*kny))**2/(kk*N.pi/(2*kny))**2
    #w=1.0
    #w=N.exp(-rsm*rsm*kk*kk/2)
    ans=N.sum((pklin(kk)*w**2*kk**(2+m)*jrk*deltak)/(2*N.pi**2), 1)
    ans=N.reshape(ans, sk)
    #ans[N.where(r==0)]=1 if n==0 else 0
    return ans
    
    
def mu1mu2_to_thetaphi(r1, r2, r3, mu1, mu2):
    r1=N.asarray(r1)
    r2=N.asarray(r2)
    r3=N.asarray(r3)
    mu1=N.asarray(mu1)
    mu2=N.asarray(mu2)
    
    
    mu12=(r1**2+r2**2-r3**2)/(2*r1*r2)
    #th12=N.arccos(mu12)
    sth12=N.sqrt(1-mu12**2)
    
    if (r1.size==1):
        phi=N.arctan((mu2-mu1*mu12)/(sth12*mu1))
        theta=N.arccos(mu1/N.cos(phi))
    else:
        phi=N.zeros([mu1.size, r1.size], dtype=N.float)
        theta=N.zeros([mu1.size, r1.size], dtype=N.float)
        mu1mat=N.transpose(N.tile(mu1, (r1.size, 1)))
        mu2mat=N.transpose(N.tile(mu2, (r1.size, 1)))
        phi=N.arctan((mu2mat-N.outer(mu1, mu12))/(N.outer(mu1, sth12)))
        theta=N.arccos(mu1mat/N.cos(phi))
    
    
    return theta, phi
