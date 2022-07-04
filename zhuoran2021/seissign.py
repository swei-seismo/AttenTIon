## Functions to analyse seismic signal
## Written by S. Wei, July 2014
import scipy.signal as sg
import numpy as np
from scipy.odr import odrpack as odr
from scipy.odr import models

def polyodr(x,y,n,verbose=False,itmax=200):
    ''' Performs a polynomial least squares fit to the data,
    with errors! Uses scipy odrpack, but for least squares.
    IN:
       x,y (arrays) - data to fit
       n (int)      - polinomial order
       verbose      - can be 0,1,2 for different levels of output
                      (False or True are the same as 0 or 1)
       itmax (int)  - optional maximum number of iterations
    OUT:
       coeff -  polynomial coefficients, lowest order first
       err   - standard error (1-sigma) on the coefficients
    --Tiago, 20071114
    '''
    # http://www.scipy.org/doc/api_docs/SciPy.odr.odrpack.html
    # see models.py and use ready made models!!!!
    func   = models.polynomial(n)
    mydata = odr.Data(x, y)
    myodr  = odr.ODR(mydata, func, maxit=itmax)
    # Set type of fit to least-squares:
    myodr.set_job(fit_type=2)
    if verbose == 2: myodr.set_iprint(final=2)
    fit = myodr.run()
    # Display results:
    if verbose: fit.pprint()
    if fit.stopreason[0] == 'Iteration limit reached':
        print('(WWW) poly_lsq: Iteration limit reached, result not reliable!')
    # Results and errors
    coeff = fit.beta[::-1]
    err   = fit.res_var
    
    return coeff,err

def lincorrcoef(x,y):
## Calculate Linear Correlation Coefficient
## Input: x - 1D numpy array
##        y - 1D numpy array
## Output: r - scalar
    denom = np.sqrt(np.sum((x-np.mean(x))**2)*np.sum((y-np.mean(y))**2))
    numer = np.sum((x-np.mean(x))*(y-np.mean(y)))
    r = numer/denom
    
    return r


def seisfilter(data0,dt,btype,cutfreq):
## Filter seismic signal with a Butterworth filter
## Input:   data0   - original data (1-D numpy array)
##          dt      - sample interval (float)
##          btype   - 'lowpass', 'highpass', 'bandpass', 'bandstop'
##          cutfreq - critical frequency (list of two floats for bandpass, one float otherwise)
## Output:  data    - filtered data (1-D numpy array)
    Wn = np.array(cutfreq)/(1.0/2.0/dt)
    b, a = sg.butter(4,Wn,btype)
    data = sg.lfilter(b,a,data0)
    
    return data


def taper(data0,width=0.05,ttype='hanning',side='both'):
## Apply a symmetric taper to each end of data
## Input:   data0   - original data (1-D numpy array)
##          width   - width on each end to taper (0.0 to 0.5)
##          ttype   - 'hanning', 'hamming', 'blackman'
##          side    - sample interval (float)
## Output:  data    - taperd data (1-D numpy array)
    npts = data0.size
    if width >= 0 and width < 0.5:
        wlen = int(npts*width)
    else:
        raise ValueError("taper width should between 0.0 and 0.5")
    if ttype == 'hanning':
        window = np.hanning(2*wlen)
    elif ttype == 'hamming':
        window = np.hamming(2*wlen)
    elif ttype == 'blackman':
        window = np.blackman(2*wlen)
    else:
        raise ValueError("taper type needs to be hanning/hamming/blackman")
    if side == 'left':
        taper = np.hstack((window[:wlen], np.ones(npts - wlen)))
    elif side == 'right':
        taper = np.hstack((np.ones(npts-wlen),window[len(window)-wlen:]))
    else:
        taper = np.hstack((window[:wlen],np.ones(npts-2*wlen),window[len(window)-wlen:]))
    data = data0 * taper
    
    return data


def spec(data,dt,maxfreq=0,stype='amp'):
## Calculate spectrum of signal
## Input:   data    - time domain signal (1-D numpy array)
##          dt      - time interval (float)
##          maxfreq - maximum frequency
##          stype   - 'amp','power',or 'phase' (string)
## Output:  spec    - spectrum (1-D numpy array)
##          freq    - frequency (1-D numpy array)
    if maxfreq == 0:
        maxfreq = 1.0/2.0/dt
    freq=np.fft.fftfreq(data.size,d=dt)
    if stype == 'amp':
        spec=(np.absolute(np.fft.fft(data)))*dt
    elif stype == 'power':
        spec=((np.absolute(np.fft.fft(data)))*dt)**2
    elif stype == 'phase':
        spec = (np.angle(np.fft.fft(data)))*dt
    else:
        raise ValueError("Spectrum type: amp/power/phase")
    ind = np.argsort(freq)
    freq = freq[ind]
    spec = spec[ind]
    pind = (freq>0)&(freq<maxfreq)
    freq = freq[pind]
    spec = spec[pind]
    
    return spec,freq


def smooth(x,window_len,window='hanning'):
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    if window_len % 2 == 0:
        window_len = window_len + 1
    s = np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = eval('np.'+window+'(window_len)')
    y = np.convolve(w/w.sum(),s,mode='valid')
    
    return y[int(window_len/2):-int(window_len/2)]