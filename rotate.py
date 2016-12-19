"""
this function rotates the supplied list
(rotation in the case of an interferogram refers to moving the left part of the burst to the end of the list and the 
right part of the burst, including the maximum of the burst, to the beginning of the list
"""


def rotate(burst,lenlist=None):
    """
    rotates the provided burst list into a list of length lenlist
    if lenlist is not specified, the rotated list will be the same size as the burst list
    """
    try:
        np
    except NameError:
        import numpy as np
    zpd = np.where(burst == max(burst))[0][0] # find the zpd
    if lenlist is None: # if there is no list length specified, cheat and return the reordered list
        return np.append(burst[zpd:],burst[:zpd])
        #lenlist = len(burst)
    zer = np.zeros(lenlist) # create a zeros array
    for ind,val in enumerate(burst): # extract the n points around the zpd into the zeros array
        if ind < zpd:
            zer[-(zpd-ind)] = val # put the zpd and lower data at the end of the zeros array
        elif ind >= zpd:
            zer[ind-zpd] = val # put the zpd and higher data at the beginning of zeros array
    if np.where(zer == max(zer))[0][0] != 0: # error check just in case the zpd is not located at index 0 using this method
        print np.where(zer == max(zer))[0][0]
        raise IndexError('The ZPD is not located at index zero, fix the script')
    return zer

