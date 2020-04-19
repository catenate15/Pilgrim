#!/usr/bin/python3.6
'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2020.2
License     : MIT/x11

Copyright (c) 2020, David Ferro Costas (david.ferro@usc.es) and
Antonio Fernandez Ramos (qf.ramos@usc.es)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software
is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
---------------------------

*----------------------------------*
| Module     :  common             |
| Sub-module :  gaussian           |
| Last Update:  2020/02/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

Interface for the Electronic Structure calculation
using GAUSSIAN
'''

#==========================================================#
import os                                                  #
import time
import common.Exceptions as Exc
import subprocess
from   common.fncs      import xyz
from   common.fncs      import clean_lines
from   common.fncs      import flatten_llist
from   common.fncs      import symbols2atonums
from   common.fncs      import atonums2masses
from   common.physcons  import ANGSTROM
from   common.files     import read_file
from   common.files     import write_file
from   common.pgs       import get_pgs
#==========================================================#


#==========================================================#
#                  TO BE MODIFIED BY USER                  #
#==========================================================#
#EXE  = "/home/programs/G09_64david/g09/g09"               #
#FCHK = "/home/programs/G09_64david/g09/formchk"           #
#----------------------------------------------------------#
# in .bashrc:                                              #
#  export GauExe="/home/programs/G09_64david/g09/g09"      #
#  export GauFchk="/home/programs/G09_64david/g09/formchk" #
#==========================================================#




#=======================================================#
# FUNCTIONS FOR READING THE LOG FILE OF GAUSSIAN PROG.  #
#=======================================================#
def split_gaulog_into_gaublocks(filename):
    # Key to split the system into blocks
    str_end='Normal termination'
    # Read log/out file and join lines as string
    with open(filename,'r') as asdf: lines = asdf.readlines()
    text  = "".join(lines)
    # Divide by str_end (last element has to be excluded)
    blocks = text.split(str_end)[:-1]
    # For some reason, sometimes Gaussian prints a block
    # without information. In these cases, the block consists
    # of a few lines. Here, we exclude that cases
    #print [len(block.split("\n")) for block in blocks]
    blocks = [block for block in blocks if len(block.split("\n")) > 300]
    # Remove the lines list and the whole text
    del lines, text
    # No normal termination?
    return blocks
#-------------------------------------------------------#
def get_data_from_maintext(mtext):
    key_geom1 = 'Z-Matrix orientation'
    key_geom2 = 'Input orientation'
    key_geom3 = 'Standard orientation'
    key_force = "     Forces ("
    key_zmat  = 'Final structure in terms of initial Z-matrix:'
    key_end   = "------------------"
    key_1bar  = "\\"
    key_oniom = "ONIOM: extrapolated energy"

    # (a) Find cartesian coordinates
    if   key_geom1 in mtext: geom = mtext.split(key_geom1)[-1]
    elif key_geom2 in mtext: geom = mtext.split(key_geom2)[-1]
    elif key_geom3 in mtext: geom = mtext.split(key_geom3)[-1]
    else                   : geom,xcc = None, None
    if geom is not None:
       # convert to list of lines and get the lines associated to geometry
       geom = "\n".join(geom.split("\n")[5:])
       idx  = geom.find(key_end)
       geom = geom[:idx].strip()
       # convert to list of floats
       geom = [line.split() for line in geom.split("\n") if line.strip() != ""]
       xcc  = [[float(x),float(y),float(z)] for (_,atnum,_,x,y,z) in geom]
       xcc  = flatten_llist(xcc)
    # (b) Find forces --> gradient
    if key_force in mtext:
       force = mtext.split(key_force)[-1]
       # convert to list of lines and get the lines associated to forces
       force = "\n".join(force.split("\n")[3:])
       idx = force.find(key_end)
       force = force[:idx]
       # convert to list of floats
       force = [line.split() for line in force.split("\n") if line.strip() != ""]
       gcc  = [[-float(gx),-float(gy),-float(gz)] for (_,atnum,gx,gy,gz) in force]
       gcc  = flatten_llist(gcc)
    else: gcc = None
    # (c) Find z-matrix
    if  key_zmat in mtext:
        lines_zmat = mtext.split(key_zmat)[-1].strip().split("\n")
        for idx,line in enumerate(lines_zmat):
            line = line.strip()
            if line == "" or key_1bar in line:
               idx1 = idx
               break
        zmat = [line.strip() for line in lines_zmat[:idx1] if "Variables:" not in line]
    else: zmat = None
    # (d) ONIOM energy?
    E_ONIOM = None
    for line in mtext.split("\n")[::-1]:
        if key_oniom in line:
           E_ONIOM = float(line.split()[-1])
           break
    # Convert xcc to bohr
    if xcc is not None: xcc = [xi/ANGSTROM for xi in xcc]
    # Return data
    return xcc, gcc, zmat, E_ONIOM
#-------------------------------------------------------#
def get_data_from_archive(summary):
    # Keywords to look for
    key_hag  = "#"
    key_1bar ='\\'
    key_2bar ='\\\\'
    key_ver  ='Version'
    key_en1  ='State='
    key_en2  ='RMSD='
    key_imag ='NImag='
    # logical variables
    Lzmat = False
    Lxyz  = False
    Lhess = False
    # (a) the command line
    idx1 = summary.find(key_hag,0)+len(key_hag)
    idx2 = summary.find(key_2bar,idx1)
    commands = summary[idx1:idx2]
    # (b) the comment line
    idx3 = summary.find(key_2bar,idx2)+len(key_2bar)
    idx4 = summary.find(key_2bar,idx3)
    comment = summary[idx3:idx4]
    # (c) charge and multiplicity
    idx5 = summary.find(key_2bar,idx4)+len(key_2bar)
    idx6 = summary.find(key_1bar,idx5)
    ch,mtp = [int(value) for value in summary[idx5:idx6].split(",")]
    # (d) z-matrix or Cartesian coordinates
    idx7 = summary.find(key_ver,idx6)
    geom = [string for string in summary[idx6+len(key_1bar):idx7].split(key_1bar) if string != ""]
    if len(geom[0]) <= 4:
       Lzmat   = True
       zmat    = list(geom)
       symbols = [line.split(",")[0] for line in zmat if "=" not in line]
       xcc     = None
    else:
       Lxyz    = True
       zmat    = None
       symbols = [line.split(",")[0]   for line in geom]
       # sometimes, this line has 5 elements instead of 4
       # for this reason, coordinates are extracted with [-3:]
       # instead of [1:]
       xyz = [line.split(",")[-3:] for line in geom]
       xyz = [[float(x),float(y),float(z)] for (x,y,z) in xyz]
       xcc = flatten_llist(xyz)
    # (e) Energy and other info
    idx8a = summary.find(key_ver,idx7)
    idx8b = summary.find(key_en1,idx7)
    idx8  = max(idx8a,idx8b)
    idx9  = summary.find(key_1bar,idx8)+len(key_1bar)
    idx10 = summary.find(key_en2,idx9)
    str_energies = summary[idx9:idx10].replace(key_1bar," ")
    energies = str_energies.split()
    energies = [line.split("=") for line in energies]
    # remove S**2 (for open-shell)
    energies = [(float(energy),level.strip()) for level,energy in energies if not level.strip().startswith("S2")]
    # (f) Hessian matrix
    Lhess = key_imag in summary
    if Lhess:
       idx11=summary.find(key_imag,0)+len(key_imag)
       idx12=summary.find(key_2bar,idx11)
       num_imag = int(summary[idx11:idx12])
       idx12 += len(key_2bar)
       idx13=summary.find(key_2bar,idx12)
       # low-triangle hessian
       Fcc = [float(value) for value in summary[idx12:idx13].split(",")]
    else:
       num_imag = -1
       Fcc  = None
    # (g) Gradient may appera after hessian matrix
    if Lhess:
       idx13 += len(key_2bar)
       idx14=summary.find(key_2bar,idx13)
       gcc = [float(value) for value in summary[idx13:idx14].split(",")]
    else:
       gcc = None
    # remove dummies from symbols
    symbols = [symbol for symbol in symbols if not symbol.lower().startswith("x")]
    # Convert xcc to bohr
    if xcc is not None: xcc = [xi/ANGSTROM for xi in xcc]
    return commands,comment,ch,mtp,symbols,xcc,gcc,Fcc,energies,num_imag,zmat
#-------------------------------------------------------#
def get_data_from_gaublock(gaublock):
    '''
    gaublock is a string 
    gaublock contains the info from begining till Normal termination
    '''
    # Divide the block into the summary part and the rest
    key_start='GINC'
    key_end  ='@'
    mtext    = gaublock.split(key_start)[0]
    summary  = gaublock.split(key_start)[1].split(key_end)[0]
    # Remove the initial blank space of each line in summary
    # Also remove the line breaks
    summary  = "".join([line.strip() for line in summary.split("\n")])
    # Data in summary (aka archive)
    commands,comment,ch,mtp,symbols,xcc,gcc,Fcc,energies,num_imag,zmat = get_data_from_archive(summary)
    # Data in main text (excluding archive)
    xcc_mt, gcc_mt, zmat_mt, E_oniom = get_data_from_maintext(mtext)
    # If data not in archive --> get from main text
    if xcc  is None and  xcc_mt is not None: xcc  =  xcc_mt
    if gcc  is None and  gcc_mt is not None: gcc  =  gcc_mt
    if zmat is None and zmat_mt is not None: zmat = zmat_mt
    # Return data
    return commands,comment,ch,mtp,symbols,xcc,gcc,Fcc,energies,E_oniom,num_imag,zmat
#-------------------------------------------------------#
def read_gaussian_log(filename,target_level=None):
    if not os.path.exists(filename): return
    # split lines into blocks (in case of Link1)
    blocks = split_gaulog_into_gaublocks(filename)
    # Get info of each block
    data = [get_data_from_gaublock(block) for block in blocks]
    # There is nothing to return
    if data == []: return [None]*12
    # Localize data with hessian matrix
    IDX = -1
    for idx,data_i in enumerate(data):
        Fcc =  data_i[7]
        if Fcc is not None:
           IDX = idx
           break
    # Return the best set of data (the last with the hessian or the last block)
    commands,comment,ch,mtp,symbols,xcc,gcc,Fcc,energies,E_oniom,num_imag,zmat = data[IDX]
    # If user does not ask for level, send one of lowest energy
    if target_level is None:
       energies.sort()
       energy,level = energies[0]
    else:
       IDX = None
       exception = Exc.LevelNotFound()
       exception._var = target_level
       for idx,(energy,level) in enumerate(energies):
           if level.lower() == target_level.lower():
               IDX = idx
               break
       if IDX is None: raise exception
       energy, level = energies[IDX]
    # oniom?
    if E_oniom is not None:
       energy = E_oniom
       level  = "ONIOM"
    # Return data
    return commands,comment,ch,mtp,symbols,xcc,gcc,Fcc,energy,num_imag,zmat,level
#-------------------------------------------------------#
# reading method for Pilgrim                            #
#-------------------------------------------------------#
def read_gauout(filename):
    # read gaussian file
    data_gaulog = read_gaussian_log(filename)
    # split data
    ch      = data_gaulog[2]
    mtp     = data_gaulog[3]
    symbols = data_gaulog[4]
    xcc     = data_gaulog[5]
    gcc     = data_gaulog[6]
    Fcc     = data_gaulog[7]
    V0      = data_gaulog[8]
    level   = data_gaulog[11]
    # symbols to atomic numbers
    atonums = symbols2atonums(symbols)
    # atomic mass
    atomasses = atonums2masses(atonums)
    # return data
    return xcc, atonums, ch, mtp, V0, gcc, Fcc, atomasses, level
#=======================================================#


   

#==========================================================#
def set_EXE():
    global EXE
    # Defined in this file
    if 'EXE' in globals():
        return
    # Try to export it from bashrc
    elif "GauExe" in os.environ:
        EXE = os.environ["GauExe"]
        return
    # Not found
    else: raise Exc.ExeNotDef(Exception)
#----------------------------------------------------------#
def set_FCHK():
    global FCHK
    # Defined in this file
    if 'FCHK' in globals():
        return
    # Try to export it from bashrc
    elif "GauFchk" in os.environ:
        FCHK = os.environ["GauFchk"]
        return
    # Not found
    else: raise Exc.ExeNotDef(Exception)
#----------------------------------------------------------#
def check_EXE():
    if not os.path.exists(EXE): raise Exc.ExeNotFound(Exception)
#----------------------------------------------------------#
def check_FCHK():
    if not os.path.exists(FCHK): raise Exc.ExeNotFound(Exception)
#----------------------------------------------------------#
def execute(ifile,ofile,err,folder=None):
    # Add folder to names
    if folder is not None:
       if not folder.endswith("/"): folder = folder + "/"
       ifile = folder + ifile
       ofile = folder + ofile
       err   = folder + err
    # Execution file?
    set_EXE()
    check_EXE()
    # Exception
    exception = Exc.CalcFails(Exception)
    exception._var = ofile
    # Execute Gaussian
    command = "%s <%s 1>%s 2>%s"%(EXE,ifile,ofile,err)
    # Try up to NN times
    NN = 5
    for ii in range(NN):
        try   : status = os.system(command)
        except: status = 0
        # interrupted by ctrl+c (status = 2?)
        if status == 2: raise KeyboardInterrupt
        # Get log status
        logstatus = log_status(ofile)
        # Act according log status
        # (a) Normal termination
        if   logstatus == 1       : break
        # (b) empty file, no file or open-new-file error
        elif logstatus in (-1,0,2):
             time.sleep(1)
             continue
        # (c) other type of error
        else: raise exception
    return status
#----------------------------------------------------------#
def log_status(ofile):
    '''
    -1 --> file does not exists
     0 --> files exists but it is empty
     1 --> normal termination
     2 --> open-new-file
     3 --> other...
    '''
    if not os.path.exists(ofile)      : return -1
    with open(ofile,'r') as asdf: lines = asdf.readlines()
    if lines == []                    : return  0
    lline = lines[-1].lower()
    if   "normal termination" in lline: return  1
    elif "open-new-file"      in lline: return  2
    else                              : return  3
#----------------------------------------------------------#
def normal_termination(ofile):
    lines = read_file(ofile)
    if len(lines) == 0: return False
    lastline = lines[-1].lower()
    if "normal termination" in lastline: return True
    else                               : return False
#----------------------------------------------------------#
def genfchk(chk,fchk,err,folder=None):
    # Add folder to names
    if folder is not None:
       if not folder.endswith("/"): folder = folder + "/"
       chk   = folder + chk
       fchk  = folder + fchk
       err   = folder + err
    # fchk tool?
    set_FCHK()
    check_FCHK()
    # Exception
    exception = Exc.CalcFails(Exception)
    exception._var = chk
    # Execute fchk tool
    command = "%s %s %s 1>%s 2>&1"%(FCHK,chk,fchk,err)
    try   : status = os.system(command)
    except: raise exception
    return status
#----------------------------------------------------------#
def iofiles(name,folder=None):
    '''
    For a given name, it returns the name of all files
    for the Gaussian calculation
    '''

    if   folder is None          : folder = ""
    elif not folder.endswith("/"): folder = folder+"/"
    else                         : folder = folder
    wname  = folder + name # whole name
    ifile  = wname + ".gjf"
    ofile  = wname + ".log"
    chk    = wname + ".chk"
    fchk   = wname + ".fchk"
    err    = wname + ".err"
    return wname, ifile, ofile, chk, fchk, err
#==========================================================#


