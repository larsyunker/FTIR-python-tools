"""
features to add:
- break into usable subfunctions?(one to fit kinetics, one to generate figures, one to save xl)
- orrrrrrr create a specific plot function for the fit plots
- comment all lines to indicate what's going on
- set up the githug repository
- change output to calculate the values in excel
- add an x-shifter to not consider t < 0 values
"""
book = 'Book1 - Copy'
#book = 'Book1 - padded - Copy'
#book = 'Suzuki ex - pad removed - Copy'
sheet = 'Sheet1'


def fit_kinetics(xin,yin=None,**kwargs):
    """
    tests the provided x and y lists for kinetic behaviour
    """
    ks = { # default keyword arguments
    'verbose': True, # toggle verbose
    'order':1, # what order to check the data for
    'output': False, #whether to save the output
    'outname': 'kinetic fit', # name for the output excel file
    'genplot': False, # whether to generate a plot
    'saveplot': False, # whether to save a plot
    'sheet': None, # if excel book being used, specify sheet
    'Y0guess': None, # guess for Y0
    'ybounds': None, # bounds for Y0 (2 item tuple or list)
    'kguess': None, # guess for kobserved
    'kbounds': None, # bounds for kobs (2 item tuple or list)
    'yadjust': None, # guess for y axis offset
    'yadjbounds': None, # bounds for kobs
    'dx': None, # guess for delta x
    'eqnoutput': False, # output equation format for cells in excel
    }
    if set(kwargs.keys()) - set(ks.keys()): # check for invalid keyword arguments
        string = ''
        for i in set(kwargs.keys()) - set(ks.keys()):
            string += ` i`
        raise KeyError('Unsupported keyword argument(s): %s' %string)
    ks.update(kwargs) # update defaules with provided keyword arguments
    
    from _XLSX import XLSX
    import numpy as np
    from scipy.optimize import least_squares
    #from math import exp,log
    if ks['verbose'] is True:
        import sys
    
    def zero_order(guess,x):
        """
        simulates zero order behaviour based on guess values
        """
        Y0,kobs = guess # unpack guess tuple
        ysim = np.asarray([Y0 + -kobs*i for i in x]) # calculate simulated array
        return ysim
        
    def first_order(guess,x):
        """simulates first order behaviour based on guess values
        guess should be a tuple of the form (y value at time zero, observed rate constant)"""
        Y0,kobs = guess # unpack guess tuple
        ysim = np.asarray([Y0*np.exp(-kobs*i) for i in x]) # calculate simulated array
        return ysim
    
    def second_order(guess,x):
        """
        simulates second order behaviour based on guess values"""
        Y0,kobs = guess # unpack guess tuple
        ysim = np.asarray([1/(1/Y0-kobs*i) for i in x]) # calculate simulated array
        return ysim
    
    def error(guess,x,yact,ysim,order=None):
        """calculates the error between the actual y values and the simulated y values"""
        Y0,kobs,yadj = guess # unpack guess tuple
        if ysim is None: # if function has not been handed a set of simulated y values
            if order == 0:
                ysim = zero_order((Y0,kobs),x)
            if order == 1:
                ysim = first_order((Y0,kobs),x)
            if order == 2:
                ysim = second_order((Y0,kobs),x)
        yadjusted = yact + yadj # shift the ys in order to check fit
        resarray = yadjusted - ysim # calculate residuals
        return resarray
    
    def rsquared(ysimulated,yadjusted):
        """calculates the rsquared for the simulation"""
        avgy = sum(yadjusted)/len(yadjusted) # average y value
        SStot = sum([(i - avgy)**2 for i in ysimulated]) # sum of squares total
        SSres = sum((yadjusted - ysimulated)**2) # sum of squares residual
        #SSreg = sum((ysimulated - avgy)**2) # explained sum of squares
        Rsq = 1-(SSres/SStot)
        return Rsq,SStot,SSres
    
    if type(xin) is str and yin is None: # if the function has been handed a 
        xl = XLSX(book) # load xlsx workbook
        if ks['sheet'] is None: # if sheet is not specified, use first sheet
            ks['sheet'] = xl.wb.get_sheet_names()[0]
        (x,y),xunit,yunit = xl.pullspectrum(sheet) # import data
        ks['outname'] = xl.bookname # set outname to the input excel workbook name
    else:
        x = xin
        y = yin
        xl = None # otherwise set xl object to None
        xunit = None
        yunit = None
    
    if isinstance(xin,np.ndarray) is False: # convert to numpy arrays for ease of manipulation
        x = np.asarray(x)
    if isinstance(yin,np.ndarray) is False:
        y = np.asarray(y)
    
    #if min(y) < 0: # if any y values are negative
    #    yadj = abs(min(y)) # initial guess for y adjustment
    #else:
    #    yadj = 0.
    #Y0guess = max(y)+yadj
    
    
    bounds = [[None,None,None],[None,None,None]]
    # parameters for Y0 guess
    if ks['Y0guess'] is None: # if there is no guess for Y0
        if sum(y) < 0: # if the majority of the values are negative
            Y0guess = abs(min(y)) # initial guess for initial y value
        else: #otherwise use the maximum y values
            Y0guess = max(y)
    else:
        Y0guess = ks['Y0guess']
    if ks['ybounds'] is None: # if no bounds were specified
        bounds[0][0] = 0.
        bounds[1][0] = np.inf
    else:
        bounds[0][0] = ks['ybounds'][0]
        bounds[1][0] = ks['ybounds'][1]
    
    # parameters for kobserved guess
    if ks['kguess'] is None: # if there is no guess for kobs
        kguess = abs(max(y)-min(y))/(max(x)-min(x))
        #kguess = 1E-03 # initial guess for observed rate constant
    else:
        kguess = ks['kguess']
    if ks['kbounds'] is None:
        bounds[0][1] = kguess
        bounds[1][1] = np.inf
    else:
        bounds[0][1] = ks['kbounds'][0]
        bounds[1][1] = ks['kbounds'][1]
    
    # parameters for yadjust guess
    if ks['yadjust'] is None:
        yadjust = 0.
    else:
        yadjust = ks['yadjust']
    if ks['yadjbounds'] is None:
        bounds[0][2] = 0.
        bounds[1][2] = max(abs(min(y)),abs(max(y)))
    else:
        bounds[0][2] = ks['yadjbounds'][0]
        bounds[1][2] = ks['yadjbounds'][1]
    
    if type(ks['dx']) == int: # if a list narrower is specified
        x = x[:-ks['dx']]
        y = y[ks['dx']:]
    elif ks['dx'] == 'Auto': # if automatic narrower is specified
        rsqlst = []
    
    if ks['verbose'] is True:
        sys.stdout.write('Using guesses\nY0: %f\nk: %f\nyadjust: %f\nWith bounds: %s\n' %(Y0guess,kguess,yadjust,str(bounds)))
    
    #kmin = (abs(min(y)-max(y)))/max(x) # minimum rate to achieve the overall decrease in absorbance
    #bounds = ((max(y),kmin,0),(15,1,3))
    #res = least_squares(error,(Y0guess,kguess,yadj),bounds=bounds,args=(x,y,None),verbose=1)
    if ks['order'] == 0:
        res = least_squares(error,(Y0guess,kguess,yadjust),args=(x,y,None,0),bounds=bounds)
        ysimulated = zero_order(res.x[0:2],x)
        if ks['verbose'] is True:
            sys.stdout.write('Zero order fit:\nY0: %f\nkobs: %f\nyadjust: %f\n' %tuple(res.x))
    if ks['order'] == 1:
        if ks['dx'] is None:
            res = least_squares(error,(Y0guess,kguess,yadjust),args=(x,y,None,1),bounds=bounds)
        elif type(ks['dx']) is int:
            res = least_squares(error,(Y0guess,kguess,yadjust),args=(x,y,None,1),bounds=bounds)
        elif ks['dx'] == 'Auto':
            for ind,val in enumerate(x):
                sys.stdout.write('\r%d %d %d' %(ind,len(x)-1,len(y[ind:])))
                if ind == 0:
                    res = least_squares(error,(Y0guess,kguess,yadjust),args=(x,y,None,1),bounds=bounds)
                else:
                    res = least_squares(error,(Y0guess,kguess,yadjust),args=(x[:-ind],y[ind:],None,1),bounds=bounds)
                ysimulated = first_order(res.x[0:2],x[ind:])
                yadjusted = y[ind:] + res.x[2]
                Rsq,SStot,SSres = rsquared(ysimulated,yadjusted)
                if ind > 0 and Rsq < rsqlst[-1]:
                    print 'breaking'
                    ind = ind-1
                    break
                rsqlst.append(Rsq)
            res = least_squares(error,(Y0guess,kguess,yadjust),args=(x[:-(ind+1)],y[ind:],None,1),bounds=bounds)
        if ks['verbose'] is True:
            sys.stdout.write('First order fit:\nY0: %f\nkobs: %f\nyadjust: %f\n' %tuple(res.x))
        ysimulated = first_order(res.x[0:2],x) # calculate y values based on least squares output
    yadjusted = y + res.x[2] # adjust y values
    Rsq,SStot,SSres = rsquared(ysimulated,yadjusted) # calculate R-squared statistics
    print 'R-squared: %.4f' %Rsq
    
    if ks['genplot'] is True or ks['saveplot'] is True: # import pyplot if necessary
        import pylab as pl
    #reslin = least_squares(error,(Y0guess,kguess),args=(x,y,None,0))
    #liny = zero_order(reslin.x,x)
    #print reslin.x
    #
    #res2 = least_squares(error,(Y0guess,kguess),args=(x,y,None,2))
    #y2 = second_order(res2.x,x)
    #print rsquared(y2,y+res2.x[0])
    #print rsquared(first_order((14293.,0.172557),x),y)
    if ks['genplot'] is True or ks['saveplot'] is True: # if a shown plot is called for
        axwidth = 1.5
        linewidth = 2.5
        fig = pl.figure(figsize=(7.87,4.87),dpi=300)
        ax = fig.add_subplot(111)
        font = {'fontname':'Arial','fontsize':12} #font parameters for axis/text labels
        tickfont = pl.matplotlib.font_manager.FontProperties(family='Arial',size=12) # font parameters for axis ticks
        ax.spines["top"].set_visible(False) # hide bottom axis
        ax.spines["right"].set_visible(False) # hide left axis
        for axis in ["top","bottom","left","right"]:
            ax.spines[axis].set_linewidth(axwidth)
        ax.tick_params(axis='x', length=axwidth*3, width=axwidth, direction='out',top='off')
        ax.tick_params(axis='y', length=axwidth*3, width=axwidth, direction='out',right='off')
        for label in ax.get_xticklabels():
            label.set_fontproperties(tickfont) 
        for label in ax.get_yticklabels():
            label.set_fontproperties(tickfont)
        ax.plot(x,yadjusted, # plot actual values
            label='values',
            marker='o',
            markersize=3.,
            linestyle='None',
            color='#2c7fb8',
            markeredgewidth=0.25
            )
        ax.plot(x,ysimulated, # plot fitted line
            'r-',
            #label=('least squares simulated\n$R^2 = %.4f$' %(Rsq)), # label with Rsq incorporated
            label='least squares simulated', # plain label
            linewidth=linewidth,
            alpha=0.80
            )
        if ks['order'] == 0: # put the rate equation in
            ax.text(
                max(x),
                (max(y)-min(y))/5.,
                ('$f_{x} = -%.4fx+%.4f$\n$R^2 = %.4f$' %(res.x[1],res.x[0],Rsq)), # with equation
                horizontalalignment='right',
                fontsize=12,
                ) # insert rate equation in LaTeX format
        
        if ks['order'] == 1: # put the rate equation in
            ax.text(
                max(x),
                (max(y)-min(y))/5.,
                ('$f_{x} = %.4fe^{-%.4fx}$\n$R^2 = %.4f$' %(res.x[0],res.x[1],Rsq)), # with equation
                horizontalalignment='right',
                fontsize=12,
                ) # insert rate equation in LaTeX format
        #ax.set_xlim((0,3000))
        #ax.set_ylim((0,0.04))
        ax.set_xlabel(xunit,**font) # x label
        ax.set_ylabel(yunit,**font) # y label
        pl.legend( # legend properties
            loc='best',
            fontsize=12,
            fancybox=True,
            framealpha=0.
            )
        pl.tight_layout()
        if ks['saveplot'] is True:
            if ks['outname'].endswith('.xlsx'):
                out = ks['outname'].split('.')[0]
            else:
                out = ks['outname']
            pl.savefig(
                out,
                #format='svg', # scalable vector graphics
                dpi=300,
                transparent=True,
                )
        if ks['genplot'] is True:
            pl.show()
    
    if ks['output'] is True:
        if xl is None:
            xl = XLSX(ks['outname'],create=True)
            if ks['sheet'] is None: # if sheet is not specified, use first sheet
                ks['sheet'] = xl.wb.get_sheet_names()[0]
        
        cs = xl.wb.get_sheet_by_name(sheet)
        cs.cell(row=1,column=3).value = 'Adjusted Y'
        cs.cell(row=1,column=4).value = 'Simulated Y'
        for ind,val in enumerate(ysimulated):
            if ks['eqnoutput'] is False:
                cs.cell(row=ind+2,column=3).value = yadjusted[ind]
                cs.cell(row=ind+2,column=4).value = val
            if ks['eqnoutput'] is True:
                cs.cell(row=ind+2,column=3).value = ('=$B%d+$G$3' %(ind+2))
                cs.cell(row=ind+2,column=4).value = ('=$G$1*EXP(-$G$2*$A%d)' %(ind+2))
        
        cs.cell(row=1,column=6).value = 'Y0'
        cs.cell(row=1,column=7).value = res.x[0]
        cs.cell(row=2,column=6).value = 'kobs'
        cs.cell(row=2,column=7).value = res.x[1]
        #cs.cell(row=3,column=6).value = 'Yadjust'
        #cs.cell(row=3,column=7).value = res.x[2]
        cs.cell(row=4,column=6).value = 'T1/2'
        cs.cell(row=4,column=7).value = np.log(2)/res.x[1]
        cs.cell(row=5,column=6).value = 'SStot'
        cs.cell(row=5,column=7).value = SStot
        cs.cell(row=6,column=6).value = 'SSres'
        cs.cell(row=6,column=7).value = SSres
        cs.cell(row=7,column=6).value = 'R-squared'
        cs.cell(row=7,column=7).value = Rsq
        #cs.cell(row=7,column=7).value = '=1-(%s/%s)' %(xl.ind_to_rowandcolumn(5,6),xl.ind_to_rowandcolumn(4,6))
        
        xl.save()
    return res.x,Rsq


if __name__ == '__main__':
    fit_kinetics(book,
    sheet=sheet,
    output=True,
    #dx=299,
    )