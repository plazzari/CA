#!/usr/bin/python

import os,sys, getopt
import glob
import numpy  as np
import random
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def dist(x,y):   
    return np.sqrt(np.sum((x-y)**2))    
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print 'test.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'test.py -i <inputfile> -o <outputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    print 'Input file is "', inputfile
    print 'Output file is "', outputfile

    face_np_list =[]
    face_point_list =[]
    face_normal_list=[]
    filename =inputfile + '_b.out' + '*.pcd'

    listFileName = glob.glob( filename )

    for i,filein in enumerate(listFileName):
         a=np.loadtxt(filein,skiprows=11,usecols=(0, 1, 2))
         b=np.loadtxt(filein,skiprows=11,usecols=(3, 4, 5))
         face_np_list.append(a.shape[0])
         face_point_list.append(a)
         face_normal_list.append(b)

    seq = sorted(range(len(face_np_list)), key=lambda k: face_np_list[k])
    seq.reverse()
    
    aux = [ face_np_list[i] for i in seq];    face_np_list=list(aux)
    aux = [ face_point_list[i] for i in seq]; face_point_list=list(aux)
    aux = [ face_normal_list[i] for i in seq];face_normal_list=list(aux)
    
    nn=min(30,len(face_np_list)) 

# Constructing incidence matrix
    IM=np.zeros((nn, nn))   
    IM[:,:]=np.nan
   

    for i,f1 in enumerate(face_point_list[:nn]):
        for p1 in f1:
            for j,f2 in enumerate(face_point_list[:nn]):
                for p2 in f2:
                    if dist(p1,p2) < 3.:
                        IM[i,j]=1

    masked_array = np.ma.masked_invalid(IM)
    fig=plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 1,
                       width_ratios=[1,1],
                       height_ratios=[1,5]
                       )
    gs.update(wspace = 0, hspace = 0,bottom=0.42)
    ax1=plt.subplot(gs[0])
    ax1.bar(np.arange(0,30),face_np_list[0:30], alpha=0.4, color='black')
    ax1.set_xlim((0, 30))
    ax1.axis('off')
    ax1.set_xticklabels(())
    ax1.set_yticklabels(())
 

# Constructing incidence matrix extended
    IME=np.zeros((nn, nn))   
    IME[:,:]=np.nan

    for i,n1 in enumerate(face_normal_list[:nn]):
        for j,n2 in enumerate(face_normal_list[:nn]):
            if IM[i,j] == 1:

               a =face_normal_list[i]
               b =face_normal_list[j]
               v1= np.mean(a,axis=0)
               v2= np.mean(b,axis=0)
 
               IME[i,j]=angle_between(v1, v2)/3.141592653589793*180.


    masked_array = np.ma.masked_invalid(IME)
    plt.subplot(gs[1])
    plt.pcolormesh(masked_array.T, cmap = 'RdBu', edgecolors = 'None')
#   plt.colorbar()
    plt.gca().set_aspect('equal')
    theOutputFile= "PICTURES/" + inputfile + ".png"
    fig.savefig(theOutputFile)

# Generate report
    rpt_name=inputfile + '.report2'
    f =open(rpt_name, 'w')

    string1 = "Number of faces = " + str(nn) +  "\n"  ; f.write(string1);
#   string1 = "Surface=   " + str(SS)    +  "\n"  ; f.write(string1);
#   string1 = "S/V    =   " + str(SV)    +  "\n"  ; f.write(string1);

    f.close()
  
    
if __name__ == "__main__":
    main(sys.argv[1:])


