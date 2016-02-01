#!/usr/bin/python

import os,sys, getopt
import numpy  as np
import random

def APPLY_RULE(xyz,xyz_n):
    SMALL_CUBE=np.zeros((3,3,3),np.int)
    for jk in np.arange(1,100):
        for jj in np.arange(1,100):
            for ji in np.arange(1,100):
                SMALL_CUBE      = xyz[jk-1:jk+2,jj-1:jj+2,ji-1:ji+2];
                xyz_n[jk,jj,ji] = RULE_CA(SMALL_CUBE);

def RULE_CA(SMALL_CUBE):
   SM=SMALL_CUBE.copy() 
   b=SM[1,1,1].copy()

   c = IS_value2(SM,-1) 
   d = IS_value2(SM,1) 

   if (b  == -1)           : return -1

   if (b  ==  0) & ~c & ~d : return  0
   if (b  ==  0) &  c & ~d : return -1
   if (b  ==  0) & ~c &  d : return  0
   if (b  ==  0) &  c &  d : return -1 

   if (b  ==  1) & ~c & ~d : return  1
   if (b  ==  1) &  c & ~d : return  2
   if (b  ==  1) & ~c &  d : return  1
   if (b  ==  1) &  c &  d : return  2

   if (b  ==  2)           : return  2

def IS_value(SM,v):
    for jk in range(3):
        for jj in range(3):
            for ji in range(3):
                if ( jk != 1) & ( jj != 1) & ( ji != 1):
                    if SM[jk,jj,ji] == v : return 1
    return 0

def IS_value2(SM,v):
    if SM[0,1,1] == v : return 1
    if SM[2,1,1] == v : return 1
    if SM[1,0,1] == v : return 1
    if SM[1,2,1] == v : return 1
    if SM[1,1,0] == v : return 1
    if SM[1,1,2] == v : return 1
    return 0


def fill_hole3D(L,xyz):
   aux=xyz.copy() 
   for k in np.arange(1,L):
       for j in np.arange(1,L):
           for i in np.arange(1,L):
               if aux[k,j,i-1] == 1 and aux[k,j,i+1] == 1: xyz[k,j,i] =1
               if aux[k,j-1,i] == 1 and aux[k,j+1,i] == 1: xyz[k,j,i] =1
               if aux[k-1,j,i] == 1 and aux[k+1,j,i] == 1: xyz[k,j,i] =1

               if aux[k,j-1,i-1] == 1 and aux[k,j+1,i+1] == 1: xyz[k,j,i] =1
               if aux[k,j-1,i+1] == 1 and aux[k,j+1,i-1] == 1: xyz[k,j,i] =1

               if aux[k-1,j,i-1] == 1 and aux[k+1,j,i+1] == 1: xyz[k,j,i] =1
               if aux[k-1,j,i+1] == 1 and aux[k+1,j,i-1] == 1: xyz[k,j,i] =1

               if aux[k-1,j-1,i] == 1 and aux[k+1,j+1,i] == 1: xyz[k,j,i] =1
               if aux[k-1,j+1,i] == 1 and aux[k+1,j-1,i] == 1: xyz[k,j,i] =1

               if aux[k-1,j-1,i-1] == 1 and aux[k+1,j+1,i+1] == 1: xyz[k,j,i] =1
               if aux[k-1,j-1,i+1] == 1 and aux[k+1,j+1,i-1] == 1: xyz[k,j,i] =1
               if aux[k-1,j+1,i+1] == 1 and aux[k+1,j-1,i-1] == 1: xyz[k,j,i] =1
               if aux[k-1,j+1,i-1] == 1 and aux[k+1,j-1,i+1] == 1: xyz[k,j,i] =1


def find_nearest(mylist,p_in):

    p=mylist[0]
    myd=np.sqrt( (p[0]-p_in[0])**2 + (p[1]-p_in[1])**2 + (p[2]-p_in[2])**2 )

    d= myd
    res=p

    for p in mylist:
      myd=np.sqrt( (p[0]-p_in[0])**2 + (p[1]-p_in[1])**2 + (p[2]-p_in[2])**2)
      if (myd < d):
          d=myd
          res=p

    return res    

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

    a=np.loadtxt(inputfile,skiprows=1,usecols=(1, 2, 3))
    
    aa=a[:,0];# x
    bb=a[:,1];# y
    cc=a[:,2];# z
    
    nn=len(aa)
    
    px=list(set(aa)); px.sort() 
    py=list(set(bb)); py.sort()
    pz=list(set(cc)); pz.sort()
    
    x= np.zeros((nn,),dtype=np.int)
    y= np.zeros((nn,),dtype=np.int)
    z= np.zeros((nn,),dtype=np.int)
    
    for i,ia in enumerate(aa):
        for j,ja in enumerate (px):   
            if ia == ja:
                x[i] = j  
    
    for i,ib in enumerate(bb):
        for j,jb in enumerate (py):   
            if ib == jb:
                y[i] = j  
    
    for i,ic in enumerate(cc):
        for j,jc in enumerate (pz):   
            if ic == jc:
                z[i] = j  
    
    xyz  =np.zeros((101,101,101),np.int)
    xyz_n=np.zeros((101,101,101),np.int)
    
    for idx in range(nn):
        i = 50 + (x[idx] - len(px)/2)
        j = 50 + (y[idx] - len(py)/2)
        k = 50 + (z[idx] - len(pz)/2)
        xyz[i,j,k] = 1

# GAP filling
    fill_hole3D(100,xyz)
    fill_hole3D(100,xyz)

# Compute number of occupied points (Volume) 
    VV = np.sum(xyz)

# Initialize border values of the cube to -1   
    res  = xyz.copy()
    for idx in range(101):
        if np.allclose(xyz[idx,:,:],0) : res[idx,:,:]=-1.
        if np.allclose(xyz[:,idx,:],0) : res[:,idx,:]=-1.
        if np.allclose(xyz[:,:,idx],0) : res[:,:,idx]=-1.

    xyz=res.copy()

    for i in range(50):
        print i
        APPLY_RULE(xyz,xyz_n)
        if np.allclose(xyz-xyz_n,0) : break # Steady state reached, stop looping
        xyz = xyz_n.copy()

    out_b=[]
    for jk in np.arange(1,100):
        for jj in np.arange(1,100):
            for ji in np.arange(1,100):
                if xyz[jk,jj,ji] == 2:
                   point=[ji,jj,jk]
                   out_b.append(point)

# Compute number of occupied points at Surface  
    SS = len(out_b)


    headPC= "# .PCD v0.7 - Point Cloud Data file format\n" + \
            "VERSION 0.7\n"  + \
            "FIELDS x y z\n" + \
            "SIZE 4 4 4\n"   + \
            "TYPE F F F\n"   + \
            "COUNT 1 1 1\n"  + \
            "WIDTH " +  str(len(out_b)) + "\n"+ \
            "HEIGHT 1\n" + \
            "VIEWPOINT 0 0 0 1 0 0 0\n" + \
            "POINTS "  +  str(len(out_b)) + "\n" \
            "DATA ascii"

    a=np.savetxt(outputfile, out_b, delimiter=' ', header=headPC, newline='\n', comments='')
    out_b1=list(out_b)

    for jk in np.arange(1,100):
        for jj in np.arange(1,100):
            for ji in np.arange(1,100):
                if xyz[jk,jj,ji] == 2: xyz[jk,jj,ji] = -1

    for i in range(10):
        print i
        APPLY_RULE(xyz,xyz_n)
        if np.allclose(xyz-xyz_n,0) : break
        xyz = xyz_n.copy()

    out_b2=[]
    for jk in np.arange(1,100):
        for jj in np.arange(1,100):
            for ji in np.arange(1,100):
                if xyz[jk,jj,ji] == 2:
                   point=[ji,jj,jk]
                   out_b2.append(point)

    out_b_inner=[]
    for p in out_b1:
        res=find_nearest(out_b2,p)
        out_b_inner.append(res)
        
    headPC= "# .PCD v0.7 - Point Cloud Data file format\n" + \
            "VERSION 0.7\n"  + \
            "FIELDS x y z\n" + \
            "SIZE 4 4 4\n"   + \
            "TYPE F F F\n"   + \
            "COUNT 1 1 1\n"  + \
            "WIDTH " +  str(len(out_b_inner)) + "\n"+ \
            "HEIGHT 1\n" + \
            "VIEWPOINT 0 0 0 1 0 0 0\n" + \
            "POINTS "  +  str(len(out_b_inner)) + "\n" \
            "DATA ascii"
    outputfile=outputfile + '.inner'
    a=np.savetxt(outputfile, out_b_inner, delimiter=' ', header=headPC, newline='\n', comments='')
    
# Generate report
    basename=os.path.splitext(inputfile)[0]
    rpt_name=os.path.splitext(inputfile)[0]+'.report1'
    SV = float(SS)/float(VV)
    f =open(rpt_name, 'w')

    string1 = "Report for " + basename + " CA\n"  ; f.write(string1);
    string1 = "Volume =   " + str(VV)    +  "\n"  ; f.write(string1);
    string1 = "Surface=   " + str(SS)    +  "\n"  ; f.write(string1);
    string1 = "S/V    =   " + str(SV)    +  "\n"  ; f.write(string1);

    f.close()
  
    
if __name__ == "__main__":
    main(sys.argv[1:])


