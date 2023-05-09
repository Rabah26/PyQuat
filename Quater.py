# Quaternion Numbers class ..

from math import sqrt, log, acos, cos, sin, exp

class Quat:
    """
    This file creates a class quaternion type in python, 
    with methods manipulate usual algebraic operations 
    like: +,-,*,**,norm, unit,abs,exp,log .. 
    """
    def __init__(self,a,b,c,d):
        self.a=a; self.b=b; self.c=c; self.d=d

    def __repr__(self):
        aff='Quat('+str(self.a)+','+str(self.b)+','+str(self.c)+','+str(self.d)+')'
        return aff
    
    def __neg__(self):
        return Quat(-self.a,-self.b,-self.c,-self.d)
    
    def __add__(self,other):
        if not isinstance(other,Quat): 
            other=Quat(other,0,0,0)
        return Quat(self.a+other.a,self.b+other.b,self.c+other.c,self.d+other.d)
    
    def __radd__(self,other):
        if not isinstance(self,Quat): 
            self=Quat(self,0,0,0)
        return self+other
    
    def __sub__(self,other):
        other=-other
        return self+other
    
    def __rsub__(self,other):
        self=-self
        return self+other
    
    def __mul__(self,other):
        if not isinstance(other,Quat): 
            other=Quat(other,0,0,0)
        a1=self.a*other.a-self.b*other.b-self.c*other.c-self.d*other.d
        b1=self.a*other.b+self.b*other.a+self.c*other.d-self.d*other.c
        c1=self.a*other.c+self.c*other.a+self.d*other.b-self.b*other.d
        d1=self.a*other.d+self.d*other.a+self.b*other.c-self.c*other.b
        return Quat(a1,b1,c1,d1)
    
    def __rmul__(self,other):
        if not isinstance(self,Quat): 
            self=Quat(self,0,0,0)
        return self*other
    
    def __truediv__(self,other):
        if not isinstance(other,Quat): 
            return self*(1/other)
        return self*(other.qinv())
    
    def qinv(self):
        sqr=self*self.conj()
        return self.conj()*(1/sqr.a)
    
    def unit(self):
        return self*(1/self.qabs())
    
    def conj(self):
        return Quat(self.a,-self.b,-self.c,-self.d)
    
    def qabs(self):
        sqr=self*self.conj()
        return (sqr.a)**0.5
    
    def qrac(self):
        s1=self+self.qabs()
        if s1.a==0: 
            return None
        s2=((s1.a)*2)**(-0.5)
        return s1*s2
    
    def qr(self):
        return Quat(self.a,0,0,0)
    
    def qi(self):
        return Quat(0,self.b,self.c,self.d)
    
    def __pow__(self,n):
        if n % 0.5: 
            return 'Error.. the arg not a half int!'
        a, b = abs(int(n)), n % 1
        r = Quat(1,0,0,0)
        for i in range(a): r = self * r
        if b: r = r * self.qrac()
        if n < 0: r = r.qinv()
        return r

    def log(self):
        # log(q) = (ln|q|, (acos(s/|q|)(v/|v|)
        if self.qabs() == 0:
            return 
        if self.qi().qabs() == 0:
            return Quat(log(self.a),0,0,0)        
        ss, aq = self.a, self.qabs()
        uv, ac = self.qi().unit(), acos(ss/aq)
        return Quat(log(aq),0,0,0)+uv*ac
        
    def exp(self):
        if self.qabs() == 0:
            return Quat(1,0,0,0)
        if self.qi().qabs() == 0:
            return Quat(exp(self.a),0,0,0)
        ss, vv = self.a, self.qi()
        av, uv = vv.qabs(), vv.unit()
        return exp(ss)*(cos(av)+uv*sin(av))
      
