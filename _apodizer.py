"""
This is a class that generates and stores an apodization function based on the supplied method
The instance.apodize() method can then apodize a supplied list across the supplied width
"""

class apodizer(object):
    def __init__(self,method='blackman-harris-3',sigma=0.25):
        """
        builds and stores an apodizing function based on the specified method
        instance.apodize can then be called as many times as required
        sigma is only used when the Gauss apodization is used
        """
        self.sym = __import__('sympy')
        self.np = __import__('numpy')
        from sympy.abc import t
        self.t = t
        try: # check to see if the rotate function is already loaded
            rotate
        except NameError:
            from rotate import rotate
            self.rotate = rotate
        available = [
        'none', # for no apodization
        'blackman-harris-3', # Blackman-Harris 3 term
        'blackman-harris-4', # Blackman-Harris 4 term
        'blackman-harris-7', # Blackman-Harris 7 term
        'flattop', # flat top
        'gauss', # gaussian
        'hamming', # Hamming
        'hann', # Hann
        'rectangular', # apply the same intensity across the bin width
        ]
        self.method = method
        if method is None:
            self.fn = self.sym.lambdify(self.t,1.*self.t,'numpy')
            self.method = 'None'
        if self.method.lower() not in available:
            raise KeyError('The method %s is not defined in the apodizer class' %method)
        if self.method.lower().startswith('blackman-harris'):
            self.fn = self.bh(self.method)
        elif self.method.lower() == 'flattop':
            self.fn = self.flattop()
        elif self.method.lower() == 'gauss':
            self.fn = self.gauss(sigma)
        elif self.method.lower() == 'hamming':
            self.fn = self.hamming()
        elif self.method.lower() == 'hann':
            self.fn = self.hann()
        elif self.method.lower() == 'rectangular':
            self.fn = self.rectangular()
        
    def __str__(self):
        detailed_names = {
        'none':'No',
        'blackman-harris-3': 'Blackman-Harris 3-Term',
        'blackman-harris-4': 'Blackman-Harris 4-Term',
        'blackman-harris-7': 'Blackman-Harris 7-Term',
        'flattop': 'Flattop',
        'gauss': 'Gaussian',
        'hamming': 'Hamming',
        'hann': 'Hann',
        'rectangular': 'Rectangular', # apply the same intensity across the bin width
        }
        return 'Apodizer instance using %s window method' %detailed_names[self.method.lower()]
    def __repr__(self):
        return "%s('%s')" %(self.__class__.__name__,self.method) 
    
    def apodize(self,lst,binwidth=None,about='edges'):
        """
        apodizes a list of intensities based on the specified method
        the function assumes that the burst is rotated
        the about keyword specifies whether the function should be applied to the center of the list of the edges of the list (if the list is rotated)
        """
        if binwidth is None: # if binwidth is none, set to the width of the list
            binwidth = len(lst)
        x = self.np.linspace(-0.5,0.5,binwidth)
        rampy = self.fn(x)
        if len(lst) > len(rampy): # if the list sizes are not equal, pad edges with zeros to get the lists to be the same size
            rampy = self.np.pad(rampy,int((len(lst)-len(rampy))/2),'constant',constant_values=0.)
            if len(lst) != len(rampy): # if the lists aren't equal in size, add an extra zero to the end of rampy
                rampy.resize(len(lst))
        if about == 'center': # if the apodization is to be applied around the center of the spectrum
            return lst*rampy
        elif about == 'edges': # if the apodization is to be applied to the edges
            return lst * self.rotate(rampy)
        elif about == 'zpd': # if the apodization is to be applied around the zpd
            """applies equivalent ramps to either side of the zpd"""
            zpd = self.np.where(lst == max(lst))[0][0] # find the zpd
            lx = self.np.linspace(-0.5,0,len(lst[zpd:]))[-zpd:]
            rx = self.np.linspace(0,0.5,len(lst[zpd:]))
            ly = self.fn(lx)
            ry = self.fn(rx)
            rampy = self.np.append(ly,ry)
            return lst * rampy
        else:
            raise ValueError('The provided about parameter is not supported by this function')
        
    def bh(self,method):
        if method.endswith('3') or method == 'Blackman-Harris':
            a = [0.21557895,0.4973406,0.0782793,0.]
        elif method.endswith('4'):
            a = [0.35875,0.48829,0.14128,0.01168]
        elif method.endswith('7'):
            a = [0.27105140069342,0.43329793923448,0.21812299954311,0.06592544638803,0.01081174209837,0.00077658482522,0.00001388721735]
        BH = 0
        for ind,val in enumerate(a): # build sympy function
            BH += val*self.sym.cos(2*ind*self.sym.pi*self.t)
        return self.sym.lambdify(self.t,BH,'numpy') # convert into a lambda function

    def flattop(self):
        """generates the flat top function"""
        b = [0.21557895,0.41663158,0.277263158,0.083578947,0.006947368]
        fn = 0
        for ind,val in enumerate(b): # build sympy function
            fn += val*self.sym.cos(2*ind*self.sym.pi*self.t)
        return self.sym.lambdify(self.t,fn,'numpy')
    
    def gauss(self,sigma=0.25):
        """generates the gauss function based on the specified sigma"""
        fn = self.sym.exp(-self.t**2/(2.*sigma**2))
        return self.sym.lambdify(self.t,fn,'numpy')
    
    def hamming(self):
        """generates the Hamming function"""
        fn = (27./50.)+(23./50.)*self.sym.cos(2.*self.sym.pi*self.t)
        return self.sym.lambdify(self.t,fn,'numpy')
    
    def hann(self):
        """generates the Hann function"""
        fn = (self.sym.cos(self.sym.pi*self.t))**2
        return self.sym.lambdify(self.t,fn,'numpy')
    
    def rectangular(self):
        """applies a rectangular function"""
        return self.sym.lambdify(self.t,1.*self.t,'numpy')
    
    def triangular(self):
        """applies a triangular function"""
        def tri(val):
            if val < 0:
                return 2*val + 1
            elif val > 0:
                return -2*val + 1
        return self.np.vectorize(tri) # return vectorized function
        
    def plotwindow(self):
        """plots the window function for quick visualization"""
        import pylab as pl
        fig = pl.figure()
        ax = fig.add_subplot(111)
        x = self.np.linspace(-0.5,0.5,1000)
        y = self.fn(x)
        ax.axhline(0,color='k',linewidth=0.5,linestyle='dashed')
        ax.plot(x,y)
        ax.fill_between(x,0,y,alpha=0.25)
        ax.set_xlim(-0.5,0.5)
        ax.set_ylim(-0.1,1)
        pl.title(self.method)
        pl.show()

if __name__ == '__main__':
    apdz = apodizer(
    method='blackman-harris-3'
    )
    print apdz
    apdz.plotwindow()
