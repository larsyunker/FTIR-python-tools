"""
This function processes an FTIR interferogram, converting it into an IR spectrum
"""

def FTIR(lst,pad_factor=16,apdzr=None,apodize_method='Blackman-Harris-3',pcpoints=512,**kwargs):
    """
    applies a FFT to the supplied interferogram
    
    factor is the power of 2 used to zero pad the interferogram
    apodize_method is the apodization window method to be used
    apdzr (if provided) is an instance of the apodizer class (this can avoid redefinition of the same method for multiple calls of FTIR)
    apodize_method defines the apodization method to be used if no apodizer instance is supplied
    pcpoints defines the number of points around the zpd that will be used for phase correction
    kwargs allows for the function to be called from another function, allowing that function to supply the keywords
    """
    def extract_centerburst(interferogram,points=512):
        """
        extracts the n number of points surrounding the centerburst
        if the script cannot find n points centered around the zpd, the number of points will be narrowed to keep the zpd centered
        """
        zpd = np.where(interferogram == max(interferogram))[0][0] # find the zpd
        if zpd < points/2:
            return interferogram[:zpd*2]
        else:
            return interferogram[zpd-points/2:zpd+points/2]
    
    def Mertz_phase_correct(interferogram,ftresult):
        """
        phase corrects the interferogram using the Mertz method
        the interferogram is expected to be rotated
        the ftresult is the output of the fft of the interferogram
        
        series of steps:
            - FFT the apodized, rotated interferogram
            - calculate the ratio of imaginary/real along the axis
            - take the arctan of that ratio (this is the power spectrum)
            - take the cosine and sine of the power spectrum
            - multiply the real part of the FT result by the cosine
            - multiply the imaginary part of the FT result by the sine
            - return the sum of the these as the phase corrected spectrum
        """
        narrowed = extract_centerburst(interferogram,pcpoints) # extract centerburst
        narrowed = apdzr.apodize(narrowed,about='center') # apodize about the center
        #narrowed.resize(2**pad_factor) # zero pad
        rotated = rotate(narrowed) # rotate
        fftpc = np.fft.fft(rotated) # apply discrete FFT
        powspec = np.arctan(fftpc.imag/fftpc.real) # calculate the power spectrum
        powspec = np.interp(np.arange(len(ftresult)),np.linspace(0,len(ftresult),len(powspec)),powspec) # apply a linear interpolation to make the array the appropriate size
        return (ftresult*np.exp(-np.complex(0,1)*powspec)).real # phase correct the supplied ftresult
    
    def pipramp(interferogram):
        """
        multiplies the double-sided portion of the interferogram by a ramp to avoid wiggles
        the ramp goes from 0 at the beginning of the interferogram to 1 at 2 times the zpd location, and remain at one for the remainder of the pulse
        """
        zpd = np.where(interferogram == max(interferogram))[0][0] # find the zpd
        y = np.linspace(-0.5,0.5,zpd*2) # generate the ramp
        y -= 0.5 # shift so that all values are negative
        y.resize(len(interferogram)) # pad with zeros to the length of the interferogram
        y += 1. # shift y
        return interferogram * y
    
    try:
        np
    except NameError: # import numpy if not already defined
        import numpy as np
    try:
        rotate
    except NameError: # import rotate if not already defined
        from rotate import rotate
    
    if apdzr is None: # if the function has not been handed an apodizer instance, create one
        from _apodizer import apodizer
        apdzr = apodizer(apodize_method)
    ifg = lst - sum(lst)/len(lst) # shift the interferogram to be centered at zero
    #ifg = lst # do not apply the zeroing
    apdzd = apdzr.apodize(ifg,about='zpd') # apodize about the zpd
    apdzd = pipramp(apdzd) # ramp to remove wiggles
    apdzd.resize(2**pad_factor) # zero pad
    rtd = rotate(apdzd) # rotate
    ft = np.fft.fft(rtd) # fft
    pc = Mertz_phase_correct(ifg,ft) # phase correct
    return pc[:len(pc)/2]