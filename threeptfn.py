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
    rxi=N.arange(0.0, 400, 0.5)
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
    xi00i=I.interp1d(rxi, xx00)
    xi02i=I.interp1d(rxi, xx02)
    xi0n2i=I.interp1d(rxi, xx0n2)
    xi20i=I.interp1d(rxi, xx20)
    xi2n2i=I.interp1d(rxi, xx2n2)
    xi22i=I.interp1d(rxi, xx22)
    xi40i=I.interp1d(rxi, xx40)
    xi11i=I.interp1d(rxi, xx11)
    xi31i=I.interp1d(rxi, xx31)
    xi51i=I.interp1d(rxi, xx51)
    xi1n1i=I.interp1d(rxi, xx1n1)
    xi3n1i=I.interp1d(rxi, xx3n1)
    xi5n1i=I.interp1d(rxi, xx5n1)




def zeta_eq_real(rr, b1, b2, bs):
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
    
    
    xi00=xi00i(rr)
    xi11=xi11i(rr)
    xi1n1=xi1n1i(rr)
    xi20=xi20i(rr)
    
    
    def do_calc(i, j):
        t00=34./21.+b2/b1
        t22=-1./21.-bs/(12.*b1)
        t11=-1.
       
        return t00*xi00**2+t22*xi20**2+t11*xi11*xi1n1
    return (do_calc(0, 1)+do_calc(1, 2)+do_calc(0, 2))*b1**3

        
def zeta_eq_mon(rr, b1, b2, bs, f):
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
    
    
    def do_calc(i, j):
        t00=34./21.+b2/b1+(2.*f)/3.+(82.*f)/(63.*b1)+(2.*b2*f)/(3.*b1**2.)+(34.*f**2.)/(75.*b1**2.)+(38.*f**2.)/(45.*b1)+(b2*f**2.)/(9.*b1**3.)+(8.*bs*f**2.)/(675.*b1**3.)+(2.*f**3.)/(25.*b1**3.)+(2.*f**3.)/(5.*b1**2.)+(2.*f**4.)/(25.*b1**3.)
        t02=f/3.+(248.*f)/(315.*b1)+(b2*f)/(3.*b1**2.)+(2.*bs*f)/(45.*b1**2.)+(76.*f**2.)/(147.*b1**2.)+(46.*f**2.)/(63.*b1)+(b2*f**2.)/(9.*b1**3.)+(22.*bs*f**2.)/(945.*b1**3.)+(4.*f**3.)/(35.*b1**3.)+(17.*f**3.)/(35.*b1**2.)+(4.*f**4.)/(35.*b1**3.)
        t04=(13.*f**2.)/(490.*b1**2.)+f**2./(28.*b1)+(bs*f**2.)/(210.*b1**3.)+f**3./(70.*b1**3.)+f**3./(28.*b1**2.)+f**4./(70.*b1**3.)
        t22=-(1./21.)-bs/(12.*b1)+(19.*f)/(196.*b1)+(bs*f)/(42.*b1**2.)+(8413.*f**2.)/(41160.*b1**2.)+(3.*f**2.)/(20.*b1)+(3.*b2*f**2.)/(80.*b1**3.)+(41.*bs*f**2.)/(2940.*b1**3.)+(111.*f**3.)/(1960.*b1**3.)+(27.*f**3.)/(140.*b1**2.)+(27.*f**4.)/(490.*b1**3.)
        t24=f/(245.*b1)+(bs*f)/(140.*b1**2.)+(89.*f**2.)/(2744.*b1**2.)+(13.*bs*f**2.)/(1470.*b1**3.)+(43.*f**3.)/(1960.*b1**3.)+f**3./(28.*b1**2.)+f**4./(49.*b1**3.)
        
        t44=-((367.*f**2.)/(411600.*b1**2.))-(367.*bs*f**2.)/(235200.*b1**3.)+(73.*f**3.)/(78400.*b1**3.)+(2539.*f**4.)/(1411200.*b1**3.)
        t11=-1.-f/6.-(29.*f)/(30.*b1)-(661.*f**2.)/(1050.*b1**2.)-(7.*f**2.)/(15.*b1)-(87.*f**3.)/(350.*b1**3.)-(86.*f**3.)/(175.*b1**2.)-(3.*f**4.)/(70.*b1**4.)-(41.*f**4.)/(175.*b1**3.)-(3.*f**5.)/(70.*b1**4.)
        t31=-(f/(10.*b1))-(17.*f**2.)/(200.*b1**2.)-f**2./(40.*b1)-(33.*f**3.)/(700.*b1**3.)-(11.*f**3.)/(200.*b1**2.)-(3.*f**4.)/(280.*b1**4.)-(57.*f**4.)/(1400.*b1**3.)-(3.*f**5.)/(280.*b1**4.)
        
        t13=-(f/(5.*b1))-(169.*f**2.)/(900.*b1**2.)-f**2./(20.*b1)-(19.*f**3.)/(225.*b1**3.)-(97.*f**3.)/(900.*b1**2.)-f**4./(60.*b1**4.)-(67.*f**4.)/(900.*b1**3.)-f**5./(60.*b1**4.)
        t33=-(f**2./(120.*b1**2.))-(7.*f**3.)/(1200.*b1**3.)-f**3./(400.*b1**2.)-f**4./(720.*b1**4.)-(7.*f**4.)/(1800.*b1**3.)-f**5./(720.*b1**4.)
        t15=-((5.*f**2.)/(504.*b1**2.))-(17.*f**3.)/(5040.*b1**3.)-f**3./(1008.*b1**2.)-f**4./(1680.*b1**4.)-f**4./(630.*b1**3.)-f**5./(1680.*b1**4.)
        t35=-(f**3./(840.*b1**3.))-f**4./(10080.*b1**4.)-f**4./(10080.*b1**3.)-f**5./(10080.*b1**4.)
        return t00*xi00**2+t02*xi00*xi20+t04*xi00*xi40+t22*xi20**2+t24*xi20*xi40+t44*xi40**2+t11*xi11*xi1n1+t31*xi3n1*xi11+t13*xi31*xi1n1+t33*xi3n1*xi31+t15*xi1n1*xi51+t35*xi3n1*xi51
    return (do_calc(0, 1)+do_calc(1, 2)+do_calc(0, 2))*b1**3
    
    
    


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
        t00=34./21.+b2/b1+(2.*f)/3.+(82.*f)/(63.*b1)+(2.*b2*f)/(3.*b1**2.)+(34.*f**2.)/(75.*b1**2.)+(38.*f**2.)/(45.*b1)+(b2*f**2.)/(9.*b1**3.)+(8.*bs*f**2.)/(675.*b1**3.)+(2.*f**3.)/(25.*b1**3.)+(2.*f**3.)/(5.*b1**2.)+(2.*f**4.)/(25.*b1**3.)
        t02=(1./(105840.*b1**3.))*f*(2205.*b1**3.+42.*b1**2.*(124.+115.*f)+3.*b1*(735.*b2+98.*bs+3.*f*(380.+357.*f))+7.*f*(105.*b2+2.*(11.*bs+54.*f*(1.+f))))*(8.-6.*N.cos(2.*t1)-3.*N.cos(2.*(t1-theta))-24.*N.cos(2.*theta)-3.*N.cos(2.*(t1+theta))+6.*N.sqrt(3.)*N.sin(2.*t1)+3.*N.sqrt(3.)*N.sin(2.*(t1-theta))+3.*N.sqrt(3.)*N.sin(2.*(t1+theta)))
        t04=(1./(940800.*b1**3.))*f**2.*(105.*b1**2.+3.*b1*(26.+35.*f)+14.*(bs+3.*f*(1.+f)))*(216.-120.*N.cos(2.*t1)+210.*N.cos(4.*t1)+140.*N.cos(2.*(t1-2.*theta))+140.*N.cos(4.*t1-2.*theta)+80.*N.cos(2.*(t1-theta))+35.*N.cos(4.*(t1-theta))-480.*N.cos(2.*theta)+840.*N.cos(4.*theta)+80.*N.cos(2.*(t1+theta))+35.*N.cos(4.*(t1+theta))+140.*N.cos(2.*(2.*t1+theta))+140.*N.cos(2.*(t1+2.*theta))+120.*N.sqrt(3.)*N.sin(2.*t1)+210.*N.sqrt(3.)*N.sin(4.*t1)-140.*N.sqrt(3.)*N.sin(2.*(t1-2.*theta))+140.*N.sqrt(3.)*N.sin(4.*t1-2.*theta)-80.*N.sqrt(3.)*N.sin(2.*(t1-theta))+35.*N.sqrt(3.)*N.sin(4.*(t1-theta))-80.*N.sqrt(3.)*N.sin(2.*(t1+theta))+35.*N.sqrt(3.)*N.sin(4.*(t1+theta))+140.*N.sqrt(3.)*N.sin(2.*(2.*t1+theta))-140.*N.sqrt(3.)*N.sin(2.*(t1+2.*theta)))
        t22=-(1./(37044.*b1**3.))*(1764.*b1**3.-7.*f**2.*(588.*b2+215.*bs+864.*f*(1.+f))+21.*b1**2.*(147.*bs-16.*f*(31.+49.*f))-6.*b1*f*(637.*bs+6.*f*(601.+588.*f))+9.*f**2.*(4116.*b1**2.+6.*b1*(577.+882.*f)+7.*(147.*b2+8.*bs+54.*f*(3.+4.*f)))*N.cos(t1)**2.*N.cos(theta)**4.*(-2.+N.cos(2.*t1)+N.sqrt(3.)*N.sin(2.*t1))+3.*f*N.cos(theta)**2.*(6825.*b1**2.+2940.*b1*bs+18738.*b1*f+16464.*b1**2.*f+4116.*b2*f+1064.*bs*f+5481.*f**2.+21168.*b1*f**2.+6048.*f**3.+(21.*b1**2.*(-65.+196.*f)+7.*f*(147.*b2-16.*bs+27.*f*(5.+8.*f))+b1*(-588.*bs+108.*f*(23.+49.*f)))*N.cos(2.*t1)-N.sqrt(3.)*(21.*b1**2.*(-65.+196.*f)+7.*f*(147.*b2-16.*bs+27.*f*(5.+8.*f))+b1*(-588.*bs+108.*f*(23.+49.*f)))*N.sin(2.*t1)))
        t24=(1./(82320.*b1**3.))*f*(24.*(28.*b1**2.+7.*f*(17.*bs+48.*f*(1.+f))+b1*(49.*bs+4.*f*(115.+147.*f)))+2205.*f**2.*(2.+7.*b1+4.*f)*N.cos(t1)**2.*N.cos(theta)**6.*(-9.+2.*N.cos(2.*t1)+2.*N.cos(4.*t1)+6.*N.sqrt(3.)*N.sin(2.*t1))-3.*N.cos(theta)**2.*(4.*(84.*b1**2.+7.*f*(151.*bs+6.*f*(89.+104.*f))+3.*b1*(49.*bs+2.*f*(705.+1274.*f)))+(-1596.*b1**2.-3.*b1*(931.*bs+4.*(95.-637.*f)*f)+7.*f*(-89.*bs+24.*f*(11.+26.*f)))*N.cos(2.*t1)+N.sqrt(3.)*(1596.*b1**2.+3.*b1*(931.*bs+4.*(95.-637.*f)*f)+7.*f*(89.*bs-24.*f*(11.+26.*f)))*N.sin(2.*t1))-5.*f*N.cos(theta)**4.*(-8550.*b1-2100.*bs-11214.*f-28224.*b1*f-16128.*f**2.-(3.*b1*(437.+4900.*f)+7.*(46.*bs+3.*f*(223.+400.*f)))*N.cos(2.*t1)+(b1*(57.+4263.*f)+7.*(2.*bs+3.*f*(59.+116.*f)))*N.cos(4.*t1)+1311.*N.sqrt(3.)*b1*N.sin(2.*t1)+322.*N.sqrt(3.)*bs*N.sin(2.*t1)+4683.*N.sqrt(3.)*f*N.sin(2.*t1)+14700.*N.sqrt(3.)*b1*f*N.sin(2.*t1)+8400.*N.sqrt(3.)*f**2.*N.sin(2.*t1)+57.*N.sqrt(3.)*b1*N.sin(4.*t1)+14.*N.sqrt(3.)*bs*N.sin(4.*t1)+1239.*N.sqrt(3.)*f*N.sin(4.*t1)+4263.*N.sqrt(3.)*b1*f*N.sin(4.*t1)+2436.*N.sqrt(3.)*f**2.*N.sin(4.*t1)))
        
        t44=(1./(137200.*b1**3.))*f**2.*(4.*(116.*b1+7.*(29.*bs+72.*f*(1.+f)))-5.*N.cos(theta)**2.*(4.*(92.*b1+7.*(23.*bs+9.*f*(11.+16.*f)))-(388.*b1+7.*(97.*bs+36.*(1.-4.*f)*f))*N.cos(2.*t1)+N.sqrt(3.)*(388.*b1+7.*(97.*bs+36.*(1.-4.*f)*f))*N.sin(2.*t1))+1225.*f*N.cos(t1)**2.*N.cos(theta)**6.*((19.+24.*f)*N.cos(2.*t1)+(-2.+24.*f)*N.cos(4.*t1)+3.*(-11.-36.*f+N.sqrt(3.)*(5.+24.*f)*N.sin(2.*t1)))-17150.*f**2.*N.cos(t1)**4.*N.cos(theta)**8.*(-6.+4.*N.cos(2.*t1)+N.cos(4.*t1)+4.*N.sqrt(3.)*N.sin(2.*t1)-N.sqrt(3.)*N.sin(4.*t1))-15./4.*N.cos(theta)**4.*(-36.*b1-63.*bs-6384.*f-14784.*f**2.+2.*(324.*b1+7.*(81.*bs-2.*f*(61.+296.*f)))*N.cos(2.*t1)+(12.*b1+7.*(3.*bs+4.*f*(41.+106.*f)))*N.cos(4.*t1)-648.*N.sqrt(3.)*b1*N.sin(2.*t1)-1134.*N.sqrt(3.)*bs*N.sin(2.*t1)+1708.*N.sqrt(3.)*f*N.sin(2.*t1)+8288.*N.sqrt(3.)*f**2.*N.sin(2.*t1)+12.*N.sqrt(3.)*b1*N.sin(4.*t1)+21.*N.sqrt(3.)*bs*N.sin(4.*t1)+1148.*N.sqrt(3.)*f*N.sin(4.*t1)+2968.*N.sqrt(3.)*f**2.*N.sin(4.*t1)))
        t11=-(1./(175.*b1**4.))*(b1*(175.*b1**3.+105.*b1**2.*f+29.*b1*f**2.+3.*f**3.)+f*(175.*b1**4.+45.*f**3.*(1.+f)+35.*b1**3.*(11.+14.*f)+3.*b1*f**2.*(81.+82.*f)+b1**2.*f*(487.+516.*f))*N.cos(t1)**2.*N.cos(theta)**2.-N.sqrt(3.)*f*(175.*b1**4.+45.*f**3.*(1.+f)+35.*b1**3.*(11.+14.*f)+3.*b1*f**2.*(81.+82.*f)+b1**2.*f*(487.+516.*f))*N.cos(t1)*N.cos(theta)**2.*N.sin(t1))
        t31=-(1./(1400.*b1**4.))*f*(8.*b1*(35.*b1**2.+14.*b1*f+3.*f**2.)+3.*N.cos(theta)**2.*(-5.*b1*(35.*b1**2.+14.*b1*f+3.*f**2.)+2.*(60.*f**3.*(1.+f)+35.*b1**3.*(1.+4.*f)+14.*b1**2.*f*(19.+22.*f)+3.*b1*f**2.*(73.+76.*f))*N.cos(t1)**2.-N.sqrt(3.)*(60.*f**3.*(1.+f)+35.*b1**3.*(1.+4.*f)+14.*b1**2.*f*(19.+22.*f)+3.*b1*f**2.*(73.+76.*f))*N.sin(2.*t1))+5.*f*(35.*b1**3.+15.*f**2.*(1.+f)+7.*b1**2.*(9.+11.*f)+3.*b1*f*(18.+19.*f))*N.cos(t1)*N.cos(theta)**4.*(-6.*N.cos(t1)+N.cos(3.*t1)+N.sqrt(3.)*(4.*N.sin(t1)+N.sin(3.*t1))))
        
        t13=-(1./(900.*b1**4.))*f*(24.*b1*(15.*b1**2.+8.*b1*f+f**2.)+3.*N.cos(theta)**2.*(-15.*b1*(15.*b1**2.+8.*b1*f+f**2.)+2.*(60.*f**3.*(1.+f)+45.*b1**3.*(1.+4.*f)+4.*b1**2.*f*(79.+97.*f)+b1*f**2.*(259.+268.*f))*N.cos(t1)**2.-N.sqrt(3.)*(60.*f**3.*(1.+f)+45.*b1**3.*(1.+4.*f)+4.*b1**2.*f*(79.+97.*f)+b1*f**2.*(259.+268.*f))*N.sin(2.*t1))+5.*f*(45.*b1**3.+15.*f**2.*(1.+f)+b1*f*(64.+67.*f)+b1**2.*(73.+97.*f))*N.cos(t1)*N.cos(theta)**4.*(-6.*N.cos(t1)+N.cos(3.*t1)+N.sqrt(3.)*(4.*N.sin(t1)+N.sin(3.*t1))))
        t33=(1./(900.*b1**4.))*f**2.*(-24.*b1*(3.*b1+f)+50.*f*(9.*b1**2.+5.*f*(1.+f)+b1*(11.+14.*f))*N.cos(t1)**3.*N.cos(theta)**6.*(-3.*N.cos(t1)+2.*N.cos(3.*t1)+3.*N.sqrt(3.)*N.sin(t1))-15.*N.cos(t1)*N.cos(theta)**4.*(-3.*(20.*f**2.*(1.+f)+9.*b1**2.*(1.+4.*f)+b1*f*(47.+56.*f))*N.cos(t1)+(10.*f**2.*(1.+f)+7.*b1*f*(3.+4.*f)+3.*b1**2.*(-1.+6.*f))*N.cos(3.*t1)+2.*N.sqrt(3.)*(5.*(5.*f**2.*(1.+f)+2.*b1*f*(6.+7.*f)+b1**2.*(3.+9.*f))+(10.*f**2.*(1.+f)+7.*b1*f*(3.+4.*f)+3.*b1**2.*(-1.+6.*f))*N.cos(2.*t1))*N.sin(t1))-6.*N.cos(theta)**2.*(-15.*b1*(3.*b1+f)+2.*(30.*f**2.*(1.+f)+7.*b1*f*(11.+12.*f)+b1**2.*(33.+54.*f))*N.cos(t1)**2.-N.sqrt(3.)*(30.*f**2.*(1.+f)+7.*b1*f*(11.+12.*f)+b1**2.*(33.+54.*f))*N.sin(2.*t1)))
        t15=-(1./(10080.*b1**4.))*f**2.*(96.*b1*(5.*b1+f)+24.*N.cos(theta)**2.*(-21.*b1*(5.*b1+f)+2.*(30.*f**2.*(1.+f)+b1**2.*(5.+50.*f)+b1*f*(71.+80.*f))*N.cos(t1)**2.-N.sqrt(3.)*(30.*f**2.*(1.+f)+b1**2.*(5.+50.*f)+b1*f*(71.+80.*f))*N.sin(2.*t1))+63.*f*(5.*b1**2.+3.*f*(1.+f)+b1*(7.+8.*f))*N.cos(t1)*N.cos(theta)**6.*(20.*N.cos(t1)-5.*N.cos(3.*t1)+2.*N.cos(5.*t1)-12.*N.sqrt(3.)*N.sin(t1)-3.*N.sqrt(3.)*N.sin(3.*t1))+7.*N.cos(theta)**4.*(300.*b1**2.-780.*b1*f-600.*b1**2.*f-360.*f**2.-960.*b1*f**2.-360.*f**3.-4.*(75.*f**2.*(1.+f)+5.*b1**2.*(1.+25.*f)+8.*b1*f*(22.+25.*f))*N.cos(2.*t1)+(60.*f**2.*(1.+f)+5.*b1**2.*(17.+20.*f)+b1*f*(157.+160.*f))*N.cos(4.*t1)+20.*N.sqrt(3.)*b1**2.*N.sin(2.*t1)+704.*N.sqrt(3.)*b1*f*N.sin(2.*t1)+500.*N.sqrt(3.)*b1**2.*f*N.sin(2.*t1)+300.*N.sqrt(3.)*f**2.*N.sin(2.*t1)+800.*N.sqrt(3.)*b1*f**2.*N.sin(2.*t1)+300.*N.sqrt(3.)*f**3.*N.sin(2.*t1)+85.*N.sqrt(3.)*b1**2.*N.sin(4.*t1)+157.*N.sqrt(3.)*b1*f*N.sin(4.*t1)+100.*N.sqrt(3.)*b1**2.*f*N.sin(4.*t1)+60.*N.sqrt(3.)*f**2.*N.sin(4.*t1)+160.*N.sqrt(3.)*b1*f**2.*N.sin(4.*t1)+60.*N.sqrt(3.)*f**3.*N.sin(4.*t1)))
        t35=-(1./(20160.*b1**4.))*f**3.*(192.*b1+24.*N.cos(theta)**2.*(-57.*b1+2.*(60.*f*(1.+f)+b1*(37.+60.*f))*N.cos(t1)**2.-N.sqrt(3.)*(60.*f*(1.+f)+b1*(37.+60.*f))*N.sin(2.*t1))-14.*N.cos(t1)*N.cos(theta)**6.*(-40.*(31.*f*(1.+f)+b1*(7.+31.*f))*N.cos(t1)+5.*(47.*f*(1.+f)+b1*(8.+47.*f))*N.cos(3.*t1)-b1*N.cos(5.*t1)+146.*f*N.cos(5.*t1)+146.*b1*f*N.cos(5.*t1)+146.*f**2.*N.cos(5.*t1)+201.*N.sqrt(3.)*b1*N.sin(t1)+624.*N.sqrt(3.)*f*N.sin(t1)+624.*N.sqrt(3.)*b1*f*N.sin(t1)+624.*N.sqrt(3.)*f**2.*N.sin(t1)+39.*N.sqrt(3.)*b1*N.sin(3.*t1)+381.*N.sqrt(3.)*f*N.sin(3.*t1)+381.*N.sqrt(3.)*b1*f*N.sin(3.*t1)+381.*N.sqrt(3.)*f**2.*N.sin(3.*t1))+6.*N.cos(theta)**4.*(-240.*b1-1140.*f-1140.*b1*f-1140.*f**2.-2.*(475.*f*(1.+f)+b1*(198.+475.*f))*N.cos(2.*t1)+(190.*f*(1.+f)+b1*(33.+190.*f))*N.cos(4.*t1)+396.*N.sqrt(3.)*b1*N.sin(2.*t1)+950.*N.sqrt(3.)*f*N.sin(2.*t1)+950.*N.sqrt(3.)*b1*f*N.sin(2.*t1)+950.*N.sqrt(3.)*f**2.*N.sin(2.*t1)+33.*N.sqrt(3.)*b1*N.sin(4.*t1)+190.*N.sqrt(3.)*f*N.sin(4.*t1)+190.*N.sqrt(3.)*b1*f*N.sin(4.*t1)+190.*N.sqrt(3.)*f**2.*N.sin(4.*t1))-630.*f*(1.+b1+f)*N.cos(t1)**3.*N.cos(theta)**8.*(17.*N.cos(t1)-11.*N.cos(3.*t1)-N.cos(5.*t1)-13.*N.sqrt(3.)*N.sin(t1)-3.*N.sqrt(3.)*N.sin(3.*t1)+N.sqrt(3.)*N.sin(5.*t1)))
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
    kk=N.arange(10**-4, deltak*L, deltak)
    sk=r.shape
    jrk=besph(n, N.outer(r,kk))
    rsm=1.0
    
	
	#gaussian smoothing at 0.1*r
    #w=N.exp(-N.outer(r,kk)**2*0.1**2/2)
    w=N.exp(-(kk*rsm)**2/2.0)
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
