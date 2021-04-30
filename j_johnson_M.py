import numpy as np


################################################################
# Ported from original MATLAB ToolBox "Johnson Curve Toolbox"
# Dave (2021). Johnson Curve Toolbox 
# (https://www.mathworks.com/matlabcentral/fileexchange/46123-johnson-curve-toolbox), 
# MATLAB Central File Exchange. Retrieved April 29, 2021.
#
# Coded in Python by MAX PIERINI © 2021 EpiData.it (info@epidata.it)
#
################################################################


def sub_sbfit(xbar,sigma,rtb1,b2):
    
    def n_fault():
        # - assigns something to output arguments when fault=1
        gamma = 0
        delta = 0
        xlam  = 0
        xi    = 0
        return gamma, delta, xlam, xi
    
    # Preallocate:
    deriv = np.empty(shape=(4))
    dd = np.empty(shape=(4))

    # Define constants:
    tt = 1.0e-4; tol = 0.01; limit = 50; zero = 0.0; one = 1.0; two = 2.0
    three = 3.0; four = 4.0; six = 6.0; half = 0.5; quart = 0.25; one5 = 1.5
    a1 = 0.0124; a2 = 0.0623; a3 = 0.4043; a4 = 0.408; a5 = 0.479; a6 = 0.485
    a7 = 0.5291; a8 = 0.5955; a9 = 0.626; a10 = 0.64; a11 = 0.7077; a12 = 0.7466
    a13 = 0.8; a14 = 0.9281; a15 = 1.0614; a16 = 1.25; a17 = 1.7973; a18 = 1.8
    a19 = 2.163; a20 = 2.5; a21 = 8.5245; a22 = 11.346; rb1 = abs(rtb1)
    b1 = rb1 * rb1; neg = (rtb1 < zero)
    
    # Get d as first estimate of delta:
    e = b1 + one
    x = half*b1 + one
    y = abs(rb1)*np.sqrt(quart*b1+one)
    u = (x+y)**(one/three)
    w = u + one/u - one
    f = w*w*(three+w*(two+w)) - three
    e = (b2-e)/(f-e)
    if (abs(rb1)>tol):
        d = one/np.sqrt(np.log(w))
        if (d < a10):
            f = a16*d
        else:
            f = two - a21/(d*(d*(d-a19)+a22))
    else:
        f = two

    f = e*f + one
    if (f < a18):
        d = a13*(f-one)
    else:
        d = (a9*f-a4)*(three-f)**(-a5)
        
    # Get g as first estimate of gamma:
    g = zero
    if (b1 >= tt):
        if (d > one):
            if (d <= a20):
                u = a2
                y = a3
            else:
                u = a1
                y = a7
            g = b1**(u*d+y)*(a14+d*(a15*d-a11))
        else:
            g = (a12*d**a17+a8)*b1**a6
    
    # -----Main iteration starts here:-----
    stopWhile = 0 # initialize
    m    = 0
    while (stopWhile==0):
        m     = m + 1
        fault = (m > limit)
        if (fault):
            # TODO: ==================
            gamma, delta, xlam, xi = n_fault()
            return gamma,delta,xlam,xi,fault

        # Get first six moments for latest g and d values:

        hmu, fault = sub_mom(g,d)

        if (fault):
            # TODO: ==================
            gamma, delta, xlam, xi = n_fault()
            return gamma,delta,xlam,xi,fault
        s     = hmu[0]*hmu[0]
        h2    = hmu[1] - s
        fault = (h2 <= zero)
        if (fault):
            # TODO: ==================
            gamma, delta, xlam, xi = n_fault()
            return gamma,delta,xlam,xi,fault
        t    = np.sqrt(h2)
        h2a  = t*h2
        h2b  = h2*h2
        h3   = hmu[2] - hmu[0]*(three*hmu[1]-two*s)
        rbet = h3/h2a
        h4   = hmu[3] - hmu[0]*(four*hmu[2]-hmu[0]*(six*hmu[1]-three*s))
        bet2 = h4/h2b
        w    = g*d
        u    = d*d
    
        # Get derivatives:
        for j in range(2):
            for k in range(4):
                t = k + 1
                if (j==0):
                    s = hmu[k+1] - hmu[k]
                else:
                    s = ((w-t)*(hmu[k]-hmu[k+1])+(t+one)*(hmu[k+1]-hmu[k+2]))/u
                dd[k] = t*s/d

            t          = two*hmu[0]*dd[0]
            s          = hmu[0]*dd[1]
            y          = dd[1] - t

            deriv[j]   = (dd[2]-three*(s+hmu[1]*dd[0]-t*hmu[0])-one5*h3*y/h2)/h2a
            deriv[j+2] = (dd[3]-four*(dd[2]*hmu[0]+dd[0]*hmu[2])+six*(hmu[1]*t\
                         +hmu[0]*(s-t*hmu[0]))-two*h4*y/h2)/h2b

        t = one/(deriv[0]*deriv[3]-deriv[1]*deriv[2])
        u = (deriv[3]*(rbet-rb1)-deriv[1]*(bet2-b2))*t
        y = (deriv[0]*(bet2-b2)-deriv[2]*(rbet-rb1))*t
        
        # Form new estimates of g and d:
        g = g - u
        if ((b1 == zero) or (g < zero)):
            g = zero
        d = d - y

        # Assess WHILE loop:
        if ((abs(u) <= tt) and (abs(y) <= tt)):
            delta = d
            xlam  = sigma/np.sqrt(h2)
            if (neg):
                gamma  = -g
                hmu[0] = one - hmu[0]
            else:
                gamma  = g
            xi = xbar - xlam*hmu[0]
            break
    
    return gamma,delta,xlam,xi,fault


def sub_mom(g,d):
    # - evaluates 1st six moments of a johnson SB distribution, using Goodwin method

    # -----Notes:-----
    # rttwo : sqrt(2.0)
    # rrtpi : reciprocal of sqrt(pi)
    # expa  : a value such that exp(expa) does not quite cause overflow
    # expb  : a value such that 1.0 + exp(-expb) may be taken to be 1.0

    # Define constants:
    zz = 1.0e-5; vv = 1.0e-8; limit = 500; rttwo = 1.414213562; rrtpi = 0.5641895835
    expa = 80.0; expb = 23.7; zero = 0.0; quart = 0.25; half = 0.5; p75 = 0.75
    one = 1.0; two = 2.0; three = 3.0; w = g/d

    # Preallocate or initialize:
    a     = np.empty(shape=(6))
    b     = np.empty(shape=(6))
    fault = 0
    c     = np.zeros(shape=(6))

    # Trial value of h:
    if (w > expa):
        fault = 1
        return a,fault
    e = np.exp(w) + one
    r = rttwo/d
    h = p75
    if (d < three):
        h = quart*d
    k = 1
    
    # -----OUTER WHILE loop:-----
    skip     = 1 # skip 1st block of outer WHILE loop on 1st run
    stop_out = 0 # initiialize
    while (stop_out==0):
        # -----Skip this block on 1st run:-----
        if (skip==1):
            skip = 0 # don't skip anymore:
        else:
            k = k + 1
            if ( k > limit ):
                fault = 1
                return a,fault
            for i in range(6):
                c[i] = a[i]
            #  No convergence yet - try smaller h:
            h = half*h
        # -------------------------------------
        t = w
        u = t
        y = h*h
        x = two*y
        a[0] = one/e
        for i in np.arange(1,6):
            a[i] = a[i-1]/e
        v = y
        f = r*h
        m = 0
        # -----INNER WHILE loop evaluates infinite series:-----
        stop_inn = 0
        
        while (stop_inn==0):
            m = m + 1
            if (m > limit):
                break # terminate INNER WHILE loop
            for i in range(6):
                b[i] = a[i]
            u = u - f
            z = one
            if (u > -expb):
                z = np.exp(u) + z
            t = t + f
            l = (t > expb)
            if (l==0):
                s = np.exp(t) + one
            p = np.exp(-v)
            q = p
            for i in range(6):
                aa = a[i]
                p  = p/z
                ab = aa
                aa = aa + p
                if (aa == ab):
                    break # terminate this FOR loop
                if (l==0):
                    q  = q/s
                    ab = aa
                    aa = aa + q
                    l = (aa==ab)
                a[i] = aa
            y = y + x
            v = v + y
            for i in range(6):
                if (a[i]==zero):
                    fault = 1
                    return a,fault
                if (abs((a[i]-b[i])/a[i]) > vv):
                    cont_inn = 1 # continue next iteration of INNER WHILE loop
                    break        # terminate this FOR loop
                else:
                    cont_inn = 0 # don't continue
            if (cont_inn==1):
                continue # skip to start of INNER WHILE loop
            v = rrtpi*h
            for i in range(6):
                a[i] = v*a[i]
            
            for i in range(6):
                if (a[i] == zero):
                    fault = 1
                    return a,fault
                if (abs((a[i]-c[i])/a[i]) > zz):
                    break_out = 1 # signal OUTER WHILE loop should be terminated
                    break         # terminate this FOR loop
                else:
                    break_out = 0 # don't terminate OUTER WHILE loop
            if (break_out==1):
                break # terminate INNER WHILE loop
                # -----------------------------------------------------
        if (break_out==1):
            break # terminate OUTER WHILE loop

    return a,fault


def sub_sufit(xbar,sd,rb1,b2):
    # - finds parameters of Johnson SU curve with given first four moments

    # Define constants:
    tol = 0.01; zero = 0.0; one = 1.0; two = 2.0; three = 3.0; four = 4.0; six = 6.0
    seven = 7.0; eight = 8.0; nine = 9.0; ten = 10.0; sixten = 16.0; half = 0.5
    one5 = 1.5; two8 = 2.8

    b1 = rb1 * rb1; b3 = b2 - three

    # w is first estimate of exp(delta^(-2)):
    w = np.sqrt(two*b2-two8*b1-two)
    w = np.sqrt(w-one)

    # Initialize:
    stopWhile = 0
    
    if (abs(rb1)>tol):
        while (stopWhile==0): # Johnson iteration:
            w1  = w + one
            wm1 = w - one
            z   = w1*b3
            v   = w*(six+w*(three+w))
            a   = eight*(wm1*(three+w*(seven+v))-z)
            b   = sixten*(wm1*(six+v)-b3)
            y   = (np.sqrt(a*a-two*b*(wm1*(three+w*(nine+w*(ten+v)))-two*w1*z))-a)/b
            z   = y*wm1*(four*(w+two)*y+three*w1*w1)**2/(two*(two*y+w1)**3)
            v   = w*w
            w   = np.sqrt(one-two*(one5-b2+(b1*(b2-one5-v*(one+half*v)))/z))
            w   = np.sqrt(w-one)
            if (abs(b1-z) <= tol):
                y = y/w
                y = np.log(np.sqrt(y)+np.sqrt(y+one))
                if (rb1 > zero):
                    y = -y
                break # terminate WHILE loop
            else:
                continue #next iteration of WHILE loop
    else:
        # Symmetrical case - results are known
        y = zero
        
    x     = np.sqrt(one/np.log(w))
    delta = x
    gamma = y*x
    y     = np.exp(y)
    z     = y*y
    x     = sd/np.sqrt(half*(w-one)*(half*w*(z+one/z)+one))
    xlam  = x
    xi    = (half*np.sqrt(w)*(y-one/y))*x + xbar
    
    return gamma,delta,xlam,xi


def sub_sign(A,B):
    # - port of SIGN statement from FORTRAN
    #
    # If B\ge 0 then the result is ABS(A), else it is -ABS(A).
    A      = abs(A)
    if B<0:
        A = A * -1
    return A


def sub_jnsn(xbar,sd,rb1,bb2):
    
    def n_GOTO(where, gamma,delta,xlam,xi,itype,ifault):
        # - nested function of sub_jnsn to emulate 'GOTO' statements from FORTRAN
        if where == 'SN':
            # SN (Normal) distribution:
            itype = 4
            delta = one/sd
            gamma = -xbar/sd
            xlam  = one # after Simonato (2011)
        elif where == 'ST':
            # ST distribution:
            itype = 5
            y     = half + half*np.sqrt(one-four/(b1+four))
            if (rb1 > zero):
                y = one - y
            x     = sd/np.sqrt(y*(one-y))
            xi    = xbar - y*x
            xlam  = xi + x
            delta = y
        else:
            raise Exception('Unknown parameter for WHERE!')
        return gamma,delta,xlam,xi,itype,ifault
    
    # constants
    tol = 0.01; zero = 0.0; quart = 0.25
    half = 0.5; one = 1.0; two = 2.0
    three = 3.0; four  = 4.0
    itype = np.nan
    
    # Check for negative SD:
    if (sd < zero):
        itype  = np.nan; gamma = np.nan; delta = np.nan
        xlam = np.nan; xi = np.nan; ifault = 1
        return gamma,delta,xlam,xi,itype,ifault
    else:
        ifault = 0
        xi     = zero
        xlam   = zero
        gamma  = zero
        delta  = zero
        
    if (sd > zero):
        b1    = rb1*rb1
        b2    = bb2
        fault = 0
        # Test whether Lognormal (or Normal) requested:
        if (b2 >= zero):
            # Test for position relative to boundary line:
            if (b2 > b1+tol+one):
                if ((abs(rb1) <= tol) and (abs(b2-three) <= tol)):
                    gamma,delta,xlam,xi,itype,ifault = n_GOTO('SN', gamma,delta,xlam,xi,itype,ifault)
                    return gamma,delta,xlam,xi,itype,ifault
                else:
                    stopWhile = 0 # proceed to WHILE loop
                    skip      = 1 # skip first line of WHILE loop:
            else:
                if (b2 >= b1+one):
                    gamma,delta,xlam,xi,itype,ifault = n_GOTO('ST', gamma,delta,xlam,xi,itype,ifault)
                    return gamma,delta,xlam,xi,itype,ifault
                itype  = 5 # no 'itype' was included here in original FORTRAN 
                ifault = 2
                return gamma,delta,xlam,xi,itype,ifault
        else:
            stopWhile = 0 # proceed to WHILE loop
            skip      = 0 # don't skip 1st line of WHILE loop
    else:
        itype = 5
        xi    = xbar
        return gamma,delta,xlam,xi,itype,ifault
    
    while (stopWhile==0):
        # -----Skip this block on 1st run:-----
        if (skip==1):
            skip = 0 # don't skip anymore
        else:
            if (abs(rb1) < tol):
                gamma,delta,xlam,xi,itype,ifault = n_GOTO('SN', gamma,delta,xlam,xi,itype,ifault)
                return gamma,delta,xlam,xi,itype,ifault
        # Test for position relative to Lognormal line:
        x = half*b1 + one
        y = abs(rb1) * np.sqrt(quart*b1+one)
        u = (x+y)**(one/three)
        w = u + one/u - one
        u = w*w*(three+w*(two+w)) - three
        if ((b2 < zero) or (fault)):
            b2 = u
        x = u - b2
        
        if (abs(x) <= tol):
            # Lognormal (SL) distribution:
            itype = 1
            xlam  = sub_sign(one,rb1)
            u     = xlam*xbar
            x     = one/np.sqrt(np.log(w))
            delta = x
            y     = half*x*np.log(w*(w-one)/(sd*sd))
            gamma = y
            xi    = xlam*(u-np.exp((half/x-y)/x))
            return gamma,delta,xlam,xi,itype,ifault
            # SB or SU distribution:
        elif (x > zero):
            itype = 3
            gamma,delta,xlam,xi,fault = sub_sbfit(xbar,sd,rb1,b2)
            if (fault==0):
                return gamma,delta,xlam,xi,itype,ifault
            # Failure - try to fit approximate result:
            ifault = 3
            if (b2 <= b1+two):
                gamma,delta,xlam,xi,itype,ifault = n_GOTO('ST', gamma,delta,xlam,xi,itype,ifault)
                return gamma,delta,xlam,xi,itype,ifault
            else:
                skip = 0 # don't skip 1st line of WHILE loop
                continue # next iteration of WHILE loop
        else:
            itype = 2
            gamma,delta,xlam,xi = sub_sufit(xbar,sd,rb1,b2)
            stopWhile=1
            
    return gamma,delta,xlam,xi,itype,ifault


def f_johnson_M(mu,sd,skew,kurt):
    """
    Use moments to estimate parameters of a Johnson distribution.

    ARGS:
        - mu   [float] : mean
        - sd   [float] : standard deviation
        - skew [float] : skewness
        - kurt [float] : kurtosis

    RETURNS:
        - coef [tuple] : coefficients of Johnson distributions
            - gamma  [float]
            - delta  [float]
            - xi     [float]
            - lambda [float]
        - type   [str] : Johnson distribution type
                         (SL, SU, SB, SN, SY)
        - error  [str] : any error occured

    Ported from original MATLAB ToolBox "Johnson Curve Toolbox"
    Dave (2021). Johnson Curve Toolbox 
    (https://www.mathworks.com/matlabcentral/fileexchange/46123-johnson-curve-toolbox), 
    MATLAB Central File Exchange. Retrieved April 29, 2021.

    ######################################################################
    # Coded in Python by MAX PIERINI © 2021 EpiData.it (info@epidata.it) #
    ######################################################################
    """
    if not np.all((np.isscalar(mu), 
                   np.isscalar(sd),
                   np.isscalar(skew),
                   np.isscalar(kurt))):
        raise Exception('All inputs must be scalars!')
    if sd<0:
        raise Exception('Cannot have a negative SD!')
    
    gamma,delta,lambd,xi,itype,ifault = sub_jnsn(mu,sd,skew,kurt)
    coef = (gamma, delta, xi, lambd)

    if itype == 1:
        _type = 'SL'
    elif itype == 2:
        _type = 'SU'
    elif itype == 3:
        _type = 'SB'
    elif itype == 4:
        _type = 'SN'
    elif itype == 5:
        _type = 'ST'
    else:
        raise Exception('Unknown ITYPE!')
        
    if ifault == 0:
        msg = ''
    elif ifault == 1:
        msg = 'Negative SD'
    elif ifault == 2:
        msg = '(b2 < b1+two)'
    else:
        msg = 'SB failure, SL or ST used instead'
    
    return coef, _type, msg
