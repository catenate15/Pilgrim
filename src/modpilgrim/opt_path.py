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
| Module     :  modpilgrim         |
| Sub-module :  opt_path           |
| Last Update:  2020/04/18 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#--------------------------------------------------#
import copy
import gc
import os
import readline
import sys
import time
import datetime
import numpy             as np
import shutil
#--------------------------------------------------#
import modpilgrim.names      as PN
import modpilgrim.strings    as PS
import modpilgrim.rdwr       as RW
#--------------------------------------------------#
from   modpilgrim.itf        import get_spc_fnc
from   modpilgrim.diverse    import get_input_data
from   modpilgrim.diverse    import status_check
from   modpilgrim.diverse    import ffchecking
from   modpilgrim.diverse    import find_label_in_rst
#--------------------------------------------------#
from   modpilgrim.plotting   import manage_data_for_plot_mep
from   modpilgrim.plotting   import manage_data_for_plot_cvt
from   modpilgrim.plotting   import manage_data_for_plot_sct
from   modpilgrim.plotting   import write_plotfile
from   modpilgrim.exceptions import deal_with_exception
#--------------------------------------------------#
import common.fncs       as fncs
import common.adipot     as ap
import common.cvtst      as cv
import common.files      as ff
import common.internal   as intl
import common.physcons   as pc
import common.cag        as cag
import common.sct        as sct
import common.steepdesc  as sd
import common.Exceptions as Exc
import common.ispe       as ispe
#--------------------------------------------------#
from   common.Logger     import Logger
from   common.Molecule   import Molecule
from   common.steepdesc  import TSLABEL
from   common.criteria   import EPS_MEPS
#--------------------------------------------------#

WARNINGS = []

#===============================================================#
def get_itargets(targets,dpath,dctc):
    '''
    get individual targets
    '''
    # Targets?
    if (len(targets) == 0) or ("*" in targets): targets = dpath.keys()

    # generate list of individual targets
    itargets = []
    for target in targets:
        ctc, itc = PN.name2data(target)
        # check ctc
        if ctc not in dctc.keys() :
           fncs.print_string("* '%s' not in '%s'"%(target,PN.IFILE1),3)
           print("")
           continue
        if ctc not in dpath.keys():
           fncs.print_string("* '%s' not in '%s'"%(target,PN.IFILE3),3)
           print("")
           continue
        # list of itcs for the ctc
        itclist  = [itc_i for itc_i,weight_i in dctc[ctc]._itcs]
        # add itc
        if   itc is None    : itargets += [PN.struckey(ctc,itc) for itc,weight in dctc[ctc]._itcs]
        elif itc in itclist : itargets += [PN.struckey(ctc,itc)]
    return sorted(itargets)
#---------------------------------------------------------------#
def in_interval(s,si,sf,eps=1e-8):
    return si-eps <= s <= sf+eps
#---------------------------------------------------------------#
def get_masses(target,dctc,dimasses):
    ctc, itc = PN.name2data(target)
    diso     = dctc[ctc]._diso
    if   itc in diso.keys(): imod = diso[itc]
    elif "*" in diso.keys(): imod = diso["*"]
    else                   : imod = None
    if imod is None: return None
    gtsTS  = dctc[ctc].gtsfile(itc)
    TS     = Molecule()
    TS.set_from_gts(gtsTS)
    TS.apply_imods(imod,dimasses)
    masses = list(TS._masses)
    return masses
#===============================================================#


#===============================================================#
def compare_tpath(tpath,tpath2,rstfile,eps=1e-8):
    exception = Exc.RstDiffVar(Exception)
    exception._rst = rstfile
    if tpath2 is None: return
    the_vars = ("method","mu","ds","hsteps","cubic")
    for idx,var in enumerate(the_vars):
        exception._var = var
        var1, var2 = tpath[idx],tpath2[idx]
        if type(var1) != type(var2): raise exception
        if type(var1) == type(""):
            if var1 != var2: raise exception
        if type(var1) == type(1.0):
            if abs(var1-var2) > eps: raise exception
#---------------------------------------------------------------#
def compare_tcommon(tcommon,tcommon2,rstfile,eps=0.01):
    exception = Exc.RstDiffVar(Exception)
    exception._rst = rstfile
    if tcommon2 is None: return
    the_vars = ("ch","mtp","atonums","masses","mu")
    for idx,var in enumerate(the_vars):
        exception._var = var
        var1, var2 = tcommon[idx],tcommon2[idx]
        if type(var1) == type(1) or type(var1) == type(1.0):
            var1 = [var1]
            var2 = [var2]
        for ii,jj in zip(var1,var2):
            if abs(ii-jj) > eps: raise exception
#===============================================================#

#===============================================================#
def calc_mep(itarget,gtsTS,pathvars,tsoftware,TMP):
    # data in name
    ctc, itc = PN.name2data(itarget)
    # Low-Level
    if pathvars._dlevel is None:
       tcommon,drst,pathvars = obtain_mep(itarget,gtsTS,pathvars,tsoftware,TMP)
       fncs.print_string(PS.smep_table(drst,pathvars._eref),4)
    # Dual-Level
    else:
       rstfile = PN.get_rst(ctc,itc)
       xyzfile = PN.get_rstxyz(ctc,itc)
       software, tes = tsoftware
       fncs.print_string(PS.smep_init(itarget,software,pathvars._paral,pathvars.tuple_first(),\
                    pathvars.tuple_sdbw(),pathvars.tuple_sdfw()),4)
       fncs.print_string(PS.smep_ff(TMP,PN.DIR4,PN.DIR5,rstfile,xyzfile),4)

       # read rst
       tpath, tcommon, drst = ff.read_rst(rstfile)
       fncs.print_string(PS.smep_rst(rstfile,drst),4)

       fncs.print_string("Applying Dual-Level...",4)
       print("")
       dlevel_xy = [(find_label_in_rst(x,drst)[0],y) for x,y in pathvars._dlevel.items()]
       dlevel_xy = [(x,y) for x,y in dlevel_xy if x is not None]
       # Print points (sorted by s value)
       dummy = [(drst[xx][0],idx) for idx,(xx,yy) in enumerate(dlevel_xy)]
       dummy.sort()
       for s_i,idx_i in dummy:
           xx_i, yy_i = dlevel_xy[idx_i]
           fncs.print_string("%+8.4f bohr (%-6s) --> %.7f hartree"%(s_i,xx_i,yy_i),13)
       print("")
       # interpolation
       drst, points, xx, yyll, yyhl = ispe.ispe(tcommon,drst,dlevel_xy,tension=0.0)
       tdleveldata = (points,xx,yyll,yyhl)
       # table new values
       fncs.print_string(PS.smep_tableDLEVEL(drst,tdleveldata,pathvars._eref),4)
    return tcommon, drst, pathvars
#---------------------------------------------------------------#
def onesidemep(ivars,rstfile,drst):
    points = []
    for geom in sd.steepest(*ivars):
        (point,si), E, xms, gms, Fms, t = geom
        points.append(point)
        if point not in drst.keys():
           drst[point] = (si,E,xms,gms,Fms,None,None,t)
           ff.write_rst_add(rstfile,point,drst[point])
    # check convergence of last point
    if len(points) > 2:
       masses = ivars[3]
       mu     = ivars[4][1]
       epse   = ivars[4][5]
       epsg   = ivars[4][6]
       sA, EA, xA, gA = drst[points[-2]][0:4]
       sB, EB, xB, gB = drst[points[-1]][0:4]
       # check convergence of MEP
       dE   = abs(EB-EA)
       ngB  = fncs.norm(fncs.ms2cc_g(gB,masses,mu))
       conv = (dE  < epse) and (ngB < epsg)
       s_y = "%+6.3f"%drst[points[-2]][0]
       s_z = "%+6.3f"%drst[points[-1]][0]
       tt1 = "|E(s=%s)-E(s=%s)|"%(s_z,s_y)
       tt2 = "|grad(s=%s)|"%(s_z)
       ml = max(max([len(tt1),len(tt2)]),20)
       line1 = ("%%-%is = %%.2E hartree     "%ml)%(tt1,dE)
       line2 = ("%%-%is = %%.2E hartree/bohr"%ml)%(tt2,ngB)
       if (dE  < epse): line1 += " < %.2E  YES"%epse
       else           : line1 += " > %.2E  NO"%epse
       if (ngB < epsg): line2 += " < %.2E  YES"%epsg
       else           : line2 += " > %.2E  NO"%epsg
    else:
       conv  = False
       line1 = ""
       line2 = ""
    tconv = (conv,line1,line2)
    return (drst,points,tconv)
#---------------------------------------------------------------#
def obtain_mep(target,gtsTS,pathvars,tsoftware,TMP):
    ctc, itc = PN.name2data(target)

    # create folder now and not in the parallel process
    if not os.path.exists(TMP):
       try   : os.mkdir(TMP)
       except: pass

    # Files
    rstfile = PN.get_rst(ctc,itc)
    xyzfile = PN.get_rstxyz(ctc,itc)

    # print
    software, tes = tsoftware
    fncs.print_string(PS.smep_init(target,software,pathvars._paral,pathvars.tuple_first(),\
                 pathvars.tuple_sdbw(),pathvars.tuple_sdfw()),4)
    fncs.print_string(PS.smep_ff(TMP,PN.DIR4,PN.DIR5,rstfile,xyzfile),4)

    # read rst
    try: tpath2, tcommon2, drst = ff.read_rst(rstfile)
    except:
        exception = Exc.RstReadProblem(Exception)
        exception._var = rstfile
        raise exception
    fncs.print_string(PS.smep_rst(rstfile,drst),4)
    
    # correct MEP direction?
    if TSLABEL in drst.keys():
       ii_s, ii_V, ii_x, ii_g, ii_F, ii_v0, ii_v1, ii_t = drst[TSLABEL]
       ii_ic, ii_sign = pathvars._fwdir
       if not intl.ics_correctdir(ii_x,ii_v0,ii_ic,ii_sign,tcommon2[3],tcommon2[4]):
           fncs.print_string("'fwdir' variable differs from MEP direction in rst file!",4)
           fncs.print_string("* modifying rst internal dictionary...",8)
           new_drst = {}
           for key in drst.keys():
               ii_s, ii_V, ii_x, ii_g, ii_F, ii_v0, ii_v1, ii_t = drst[key]
               ii_s = -ii_s
               if ii_v0 is not None: ii_v0 = [-ii for ii in ii_v0]
               if ii_v1 is not None: ii_v1 = [-ii for ii in ii_v1]
               if "bw" in key: newkey = key.replace("bw","fw")
               else          : newkey = key.replace("fw","bw")
               new_drst[newkey] = (ii_s, ii_V, ii_x, ii_g, ii_F, ii_v0, ii_v1, ii_t)
           drst = new_drst
           del new_drst
           fncs.print_string("* rewriting rst file...",8)
           ff.write_rst(rstfile,tpath2,tcommon2,drst)

    # Extension of MEP in rst is bigger
    if drst != {}:
       lbw,lfw,sbw,sfw,V0bw,V0fw = sd.rstlimits(drst)
       pathvars._sbw = min(pathvars._sbw,sbw)
       pathvars._sfw = max(pathvars._sfw,sfw)

    # Read gts
    ts = Molecule()
    ts.set_from_gts(gtsTS)
    # scaling of frequencies
    ts.setvar(fscal=pathvars._freqscal)
    # apply iso mod
    if pathvars._masses is not None: ts.mod_masses(pathvars._masses)
    # setup
    ts.setup(mu=pathvars._mu)
    ts.ana_freqs(case="cc")

    fncs.print_string(PS.smep_ts(ts),4)

    tcommon = (ts._ch,ts._mtp,ts._atnums,ts._masses,ts._mu)
    compare_tpath(pathvars.tuple_rst(),tpath2,rstfile)
    compare_tcommon(tcommon,tcommon2,rstfile)

    # data for single-point calculation
    frozen  = RW.read_frozen(gtsTS+".frozen")
    oniom_layers = (list(pathvars._oniomh),list(pathvars._oniomm),list(pathvars._onioml))
    spc_args  = (tes,TMP,False,frozen,oniom_layers)
    spc_fnc   = get_spc_fnc(software)

    #------------#
    # First step #
    #------------#
    print("")
    fncs.print_string("Performing first step of MEP...",4)
    print("")
    inputvars = (ts._xcc,ts._gcc,ts._Fcc,ts._symbols,ts._masses,pathvars.tuple_first(),\
                 spc_fnc,spc_args,drst,pathvars._paral)
    (xms,gms,Fms),(v0,v1),(xms_bw,xms_fw) = sd.mep_first(*inputvars)
    s1bw = -float(pathvars._ds)
    s1fw = +float(pathvars._ds)

    # oniom layers?
    oniom_ok  = pathvars.isONIOMok(ts._natoms,software)
    layers    = pathvars.get_layers()
    fncs.print_string(PS.smep_oniom(layers,ts._natoms,software),8)
    if not oniom_ok: raise Exc.WrongONIOMlayers(Exception)

    # Print MEP info
    fncs.print_string(PS.smep_first(ts._symbols,ts._xms,v0,v1,layers),8)

    # write rst file
    if TSLABEL not in drst.keys():
        drst[TSLABEL] = (0.0,ts._V0,xms,gms,Fms,v0,v1,None)
        ff.write_rst_head(rstfile,pathvars.tuple_rst(),tcommon)
        ff.write_rst_add(rstfile,TSLABEL,drst[TSLABEL])

    #------------#
    # The MEP    #
    #------------#
    print("")
    fncs.print_string("Calculating MEP...",4)
    print("")
    fncs.print_string("* REMEMBER: data of each step will be saved at %s"%rstfile,7)
    fncs.print_string("            a summary will be printed when finished",7)
    # preparation
    xcc_bw = fncs.ms2cc_x(xms_bw,ts._masses,pathvars._mu)
    xcc_fw = fncs.ms2cc_x(xms_fw,ts._masses,pathvars._mu)
    args_bw = ((xcc_bw,s1bw,ts._symbols,ts._masses,pathvars.tuple_sdbw(),\
                spc_fnc,spc_args,drst,ts._Fms,"bw%i") , rstfile,drst)
    args_fw = ((xcc_fw,s1fw,ts._symbols,ts._masses,pathvars.tuple_sdfw(),\
                spc_fnc,spc_args,drst,ts._Fms,"fw%i") , rstfile,drst)
    # execution
    if pathvars._paral:
        import multiprocessing
        pool = multiprocessing.Pool()
        out  = [pool.apply_async(onesidemep,args=args) for args in [args_bw,args_fw]]
        drstbw,pointsbw,convbw = out[0].get()
        drstfw,pointsfw,convfw = out[1].get()
        del out
        # clean up pool
        pool.close()
        pool.join()
    else:
        drstbw,pointsbw,convbw = onesidemep(*args_bw)
        drstfw,pointsfw,convfw = onesidemep(*args_fw)
    # update drst
    fncs.print_string("* FINISHED!",7)
    print("")
    drst.update( drstbw )
    drst.update( drstfw )
    points = [TSLABEL]+pointsbw+pointsfw
    # empty variables
    del drstbw,pointsbw
    del drstfw,pointsfw
    # Rewrite rst
    fncs.print_string("* (re)writing file: %s (sorted version)"%rstfile,7)
    ff.write_rst(rstfile,pathvars.tuple_rst(),tcommon,drst)
    print("")
   ## restrict drst to points
   #if restrict: drst = {point:drst[point] for point in points}
    # Get limits of rst
    lbw,lfw,sbw,sfw,V0bw,V0fw = sd.rstlimits(drst)
    convbw,l1bw,l2bw = convbw
    convfw,l1fw,l2fw = convfw
    if l1bw+l1fw != "":
        fncs.print_string("* MEP convergence criteria (epse and epsg):",7)
        print("")
        if l1bw != "": fncs.print_string("%s"%l1bw,9)
        if l2bw != "": fncs.print_string("%s"%l2bw,9)
        if l1fw != "": fncs.print_string("%s"%l1fw,9)
        if l2fw != "": fncs.print_string("%s"%l2fw,9)
        print("")
        if convbw:
           pathvars.converged_in_bw(sbw)
           fncs.print_string("CRITERIA FULFILLED in backward dir.!",9)
           fncs.print_string("path stopped at sbw = %+8.4f bohr"%sbw,9)
           print("")
        if convfw:
           pathvars.converged_in_fw(sfw)
           fncs.print_string("CRITERIA FULFILLED in forward dir.!",9)
           fncs.print_string("path stopped at sfw = %+8.4f bohr"%sfw,9)
           print("")
    # write molden file
    fncs.print_string("* writing file: %s"%xyzfile,7)
    ff.rst2xyz(rstfile,xyzfile,onlyhess=True)
    print("")
    # reference energy
    if pathvars._eref is None: pathvars._eref = V0bw
    # return data
    return tcommon, drst, pathvars
#---------------------------------------------------------------#
def obtain_drp(): pass
#===============================================================#


#===============================================================#
def calc_coefs(itarget,tcommon,drst,pathvars,ltemp,symmetry=None,plotfile=None):
    # initialize dcfs
    dcfs = {}
    # sorted points
    points = sd.sorted_points(drst,hess=True)
    # Calculate adiabatic potential
    dMols, Vadi, pathvars = obtain_adipot(tcommon,drst,pathvars,symmetry)
    # Calculate CVT correction factor
    if pathvars._cvt == "yes":
       dcfs, lscvt, gibbs,lnew = obtain_cvt(dMols,points,Vadi,ltemp,pathvars,dcfs=dcfs)
    else:
       lscvt = None
    # Calculate SCT correction factor
    if pathvars._sct == "yes":
       # v1 vector along MEP
       if   pathvars._v1mode == "grad": dv1 = sct.get_numv1(drst)
       elif pathvars._v1mode == "hess": dv1 = {}
       # calculate SCT coefficient
       dcfs, tplot_sct, E0, VAG = obtain_sct(dMols,points,Vadi,ltemp,dv1,pathvars,dcfs=dcfs)
       # Calculate CAG correction factor
       dcfs = obtain_cag(Vadi,ltemp,E0,VAG,lscvt,dcfs=dcfs)
    # save data for plotting
    if plotfile is not None:
       plotdata = {}
       plotdata.update(manage_data_for_plot_mep(itarget,drst,pathvars._eref,Vadi))
       if "cvt" in dcfs.keys():
         plotdata.update(manage_data_for_plot_cvt(itarget,ltemp,Vadi.xx(),(gibbs,lnew),lscvt,dcfs["cvt"]))
       if "sct" in dcfs.keys() and tplot_sct is not None:
         plotdata.update(manage_data_for_plot_sct(itarget,dcfs["zct"],dcfs["sct"],*tplot_sct))
       write_plotfile(plotfile,plotdata)
    # Delete variables
    palpha = Vadi.get_alpha()
    pomega = Vadi.get_omega()
    del dMols
    del Vadi
    # Return data
    return dcfs, pathvars, palpha, pomega
#---------------------------------------------------------------#
def obtain_adipot(tcommon,drst,pathvars,symmetry):
    print("")
    fncs.print_string("Calculating adiabatic potential...",4)
    print ()
    ics   = pathvars.get_ics()
    icsbw = pathvars._icsbw
    icsfw = pathvars._icsfw
    if pathvars._useics == "yes" and (ics is None or len(ics) == 0): raise Exc.NoICS(Exception)
    # calculate adiabatic potential
    idata = (tcommon,drst,pathvars._eref,ics,pathvars._useics,pathvars._lowfq,pathvars._freqscal,icsbw,icsfw,symmetry)
    dMols, ccVadi, icVadi, tuple_cc, tuple_ic, listV0 = ap.path2vadi(*idata)
    # expand tuples
    data_x, lcc_frqs , lcc_tzpe = tuple_cc
    data_x, lic_frqs , lic_tzpe = tuple_ic
    # check data
    fncs.print_string(PS.sadipot_ics(ics,pathvars._useics,icsbw,icsfw),8)
    if pathvars._useics == "yes":
       ok1 = ap.ccvsic_checknfreqs(lcc_frqs,lic_frqs)
       ok2 = ap.ccvsic_checkts(data_x,lcc_frqs,lic_frqs)
       ok3 = ap.ccvsic_checkzpe(lcc_tzpe,lic_tzpe)
       checks = (ok1,ok2,ok3)
       # select spline
       if not ok1 or not ok2          : Vadi, pathvars._useics = ccVadi, "no"
       elif pathvars._useics == "yes" : Vadi, pathvars._useics = icVadi, "yes"
       else                           : Vadi, pathvars._useics = ccVadi, "no"
    else:
       ok1,ok2,ok3 = True, True, True
       Vadi, pathvars._useics = ccVadi, "no"
    # Any imag freq along the MEP?
    if pathvars._useics == "yes": lfqs = list(lic_frqs)
    else                        : lfqs = list(lcc_frqs)
    nfreqs  = min([len(xx) for xx in lcc_frqs])
    allfqreal = True
    for lfq in lfqs: 
        if len(lfq) > nfreqs: lfq = lfq[1:]
        if True in [fq<0.0 for fq in lfq]: allfqreal = False
    del lfqs
    # print checks
    fncs.print_string(PS.sadipot_checks(ok1,ok2,ok3,pathvars._useics,allfqreal),8)
    # setup spline
    Vadi.setup()
    data4print = (Vadi.xx(), Vadi.yy(), Vadi._sAG, Vadi._VAG, listV0, lcc_tzpe, lic_tzpe, pathvars._eref)
    fncs.print_string(PS.sadipot_table(*data4print),8)
    # print freqs
    fncs.print_string("- Vibrational frequencies summary: cc (ic) [in cm^-1]:",8)
    print("")
    fncs.print_string(PS.sadipot_freqs(Vadi.xx(),lcc_frqs,lic_frqs),8)
    return dMols, Vadi, pathvars
#---------------------------------------------------------------#
def obtain_cvt(dMols,points,VadiSpl,temps,pathvars,si=-float("inf"),sj=+float("inf"),dcfs={}):
    print("")
    fncs.print_string("Calculating CVT variational coefficient...",4)
    print("")
    useics = pathvars._useics
    if len(temps) == 0: raise Exc.NoTemps(Exception)
    # Only points between si and sj
    points = [pp for pp in points if si <= dMols[pp][0] <= sj]
    lcvt_s, lcvt_gamma, gibbs_matrix, gibbsTS, lnew = cv.get_cvt(dMols,points,VadiSpl,temps,useics)
    # print gibbs
    svals = [dMols[point][0] for point in points]
    fncs.print_string(PS.scvt_gibbs(svals,temps,gibbs_matrix.copy(),pathvars,gibbsTS),8)
    # print cvt coefs
    fncs.print_string(PS.scvt_coefs(lcvt_s, lcvt_gamma, temps),8)
    # save data
    dcfs["cvt"] = lcvt_gamma
    return dcfs, lcvt_s, gibbs_matrix, lnew
#---------------------------------------------------------------#
def obtain_cag(Vadi,ltemp,E0,VAG,lscvt=None,dcfs={}):
    print("")
    fncs.print_string("Calculating CAG coefficient...",4)
    print("")
    # calculate CAGTST
    dE_cagtst, cagtst = cag.calc_cag(ltemp,Vadi)
    # calculate CAGCVT
    if lscvt is not None: dE_cagcvt, cagcvt = cag.calc_cag(ltemp,Vadi,lscvt)
    else                : dE_cagcvt, cagcvt = None, None
    # In case E0 > VAG
    if E0 >= VAG:
       cagtst = list([1.0 for T in ltemp])
       if cagcvt is not None: cagcvt = list([1.0 for T in ltemp])
    # print data
    fncs.print_string(PS.scag_table(ltemp,dE_cagtst,cagtst,dE_cagcvt,cagcvt),8)
    # save in dict and return
    if cagtst is not None: dcfs["cagtst"] = cagtst
    if cagcvt is not None: dcfs["cagcvt"] = cagcvt
    return dcfs
#---------------------------------------------------------------#
def obtain_sct(dMols,points,VadiSpl,temps,dv1,pathvars,dcfs={}):
    print("")
    fncs.print_string("Calculating SCT transmission coefficient...",4)
    print("")
    # data from  pathvars
    useics = pathvars._useics
    v1mode = pathvars._v1mode
    # E0 value
    if pathvars._e0 is None:
       V1bw = pathvars._eref + VadiSpl.get_alpha()[1]
       V1fw = pathvars._eref + VadiSpl.get_omega()[1]
       if pathvars._V1R is not None: E0bw = pathvars._V1R
       else                        : E0bw = V1bw
       if pathvars._V1P is not None: E0fw = pathvars._V1P
       else                        : E0fw = V1fw
       E0 = max(E0bw,E0fw) - pathvars._eref
    else:
       E0 = pathvars._e0 - pathvars._eref
    # some checks
    if len(temps) == 0: raise Exc.NoTemps(Exception)
    if useics in ["yes",True]: case = "ic"
    else                     : case = "cc"
    # MEP LIMITS
    sbw,sfw = VadiSpl.get_alpha()[0],VadiSpl.get_omega()[0]
    # Part I - Get E0 and VAG
    E0      = sct.get_sct_part1(points,VadiSpl,E0)
    sAG,VAG = VadiSpl.get_max()
    fncs.print_string(PS.ssct_init(E0,VadiSpl,pathvars,v1mode),8)
    # Part II - Calculate tbar, bmfs and mueff
    tuple_part2 = (dMols,points,dv1,case,pathvars._muintrpl)
    svals, lkappa, ltbar, ldtbar, mu, lmueff, toignore = sct.get_sct_part2(*tuple_part2)
    fncs.print_string(PS.ssct_mueff(svals,VadiSpl,lkappa,ltbar,lmueff,mu,toignore),8)
    #----------#
    # E0 < VAG #
    #----------#
    if E0 < VAG:
       # Part III - Quantum reaction coordinate
       fncs.print_string(PS.ssct_E0VAG(E0,VAG),8)
       if pathvars._qrc is not None:
          afreq   = pathvars._qrcafreq
          lEquant = [E0+E_i for E_i in pathvars._qrclE]
          fncs.print_string(PS.ssct_qrc(pathvars),8)
          if pathvars._qrccase != 0: return dcfs, None, E0, VAG
          qrc_ZCT = sct.get_sct_part3(svals, mu   ,VadiSpl,afreq,lEquant,E0,VAG,temps)
          qrc_SCT = sct.get_sct_part3(svals,lmueff,VadiSpl,afreq,lEquant,E0,VAG,temps)
          fncs.print_string(PS.ssct_probs(qrc_SCT[1],qrc_ZCT[2],qrc_SCT[2],qrc_SCT[3],sbw,sfw),12)
          kappaI1_zct = qrc_ZCT[0]
          kappaI1_sct = qrc_SCT[0]
          qrc_Elim = lEquant[1]
       else:
          kappaI1_zct = None
          kappaI1_sct = None
          qrc_Elim    = None
       # apply QRC always?
       if not pathvars._qrcauto: qrc_Elim = None
       # Part IV - calculate thetas and probs
       fncs.print_string("Transmission probabilities for Kappa^SAG calculation:",8)
       print("")
       outZCT = sct.get_sct_part4(svals,mu    ,VadiSpl,E0)
       outSCT = sct.get_sct_part4(svals,lmueff,VadiSpl,E0)
       weights_ZCT,lE_ZCT,probs_ZCT,rpoints_ZCT,diffs_ZCT,(pZCT0,rpZCT0) = outZCT
       weights_SCT,lE_SCT,probs_SCT,rpoints_SCT,diffs_SCT,(pSCT0,rpSCT0) = outSCT
       # include also prob at E=E0 (pZCT0 and pSCT0)
       fncs.print_string(PS.ssct_probs(         [E0]+lE_SCT   ,\
                                       [pZCT0 ]+probs_ZCT,\
                                       [pSCT0 ]+probs_SCT,\
                                       [rpSCT0]+rpoints_SCT,sbw,sfw),8)
       fncs.print_string(PS.ssct_diffs(lE_SCT,diffs_SCT),8)
       # Part V - calculate coefficients
       ZCTdata = sct.get_sct_part5(lE_ZCT,probs_ZCT,weights_ZCT,E0,VAG,temps,kappaI1_zct,qrc_Elim)
       SCTdata = sct.get_sct_part5(lE_SCT,probs_SCT,weights_SCT,E0,VAG,temps,kappaI1_sct,qrc_Elim)
       ZCT,lIi_ZCT, RTE_ZCT, INTG_ZCT, bqrcZCT = ZCTdata
       SCT,lIi_SCT, RTE_SCT, INTG_SCT, bqrcSCT = SCTdata
       fncs.print_string(PS.ssct_kappa(temps,ZCT,lIi_ZCT,RTE_ZCT,E0,bqrcZCT,case="zct"),8)
       fncs.print_string(PS.ssct_kappa(temps,SCT,lIi_SCT,RTE_SCT,E0,bqrcSCT,case="sct"),8)
    #----------#
    # E0 > VAG #
    #----------#
    else:
       ZCT = [1.0 for T in temps]
       SCT = [1.0 for T in temps]
       INTG_ZCT = None
       INTG_SCT = None
       RTE_ZCT  = None
       RTE_SCT  = None
       fncs.print_string(PS.ssct_E0_above_VAG(E0,VAG),8)
       fncs.print_string(PS.ssct_onlykappa(temps,ZCT,SCT),8)
    # data for the plot
    forplot = (svals,lmueff,temps,INTG_ZCT,INTG_SCT,RTE_ZCT,RTE_SCT,E0,VAG)
    # save data
    dcfs["zct"] = ZCT
    dcfs["sct"] = SCT
    return dcfs, forplot, E0, VAG
#===============================================================#

def get_path_sctconv(itarget,gtsTS,pathvars,tsoftware,ltemp,TMP,symmetry,plotfile):
    # assert no DLEVEL will be done (just in case)
    pathvars._dlevel = None
    # The step to be increased
    hessds = abs(pathvars._hsteps*pathvars._ds)
    # some initializations
    convlist_sct = []
    convlist_lim = []
    # do we have products?
    bool_prods = pathvars._V1P is not None
    ntimesdecrease = 0
    # (a) Define the limits of MEP to calculate SCT
    sbw0, sfw0 = 0.0, 0.0
    while sbw0 > -0.5: sbw0 -= hessds
    while sfw0 < +0.5: sfw0 += hessds
    sbwN, sfwN = pathvars._sbw,pathvars._sfw
    if fncs.is_greatereq(sbwN,sbw0,EPS_MEPS): sbw0 = sbwN +   hessds
    if fncs.is_smallereq(sfwN,sfw0,EPS_MEPS): sfw0 = sfwN -   hessds
    # (b) Update pathvars and calculate MEP
    pathvars._sbw = sbw0
    pathvars._sfw = sfw0
    fncs.print_string("--> scterr keyword was detected! <--",4)
    fncs.print_string("    MEP limited to: [%.4f,%.4f] (bohr)\n"%(sbw0,sfw0),4)
    tcommon,drst,pathvars = calc_mep(itarget,gtsTS,pathvars,tsoftware,TMP)
    # (c) if MEP was calculated previously, it may be longer
    lbw,lfw,sbw,sfw,V0bw,V0fw = sd.rstlimits(drst)
    if fncs.is_smaller(sbw,sbw0,EPS_MEPS): sbw0 = sbw + 2*hessds
    if fncs.is_greater(sfw,sfw0,EPS_MEPS): sfw0 = sfw - 2*hessds
    # (c) Go ahead
    pathvars._sbw,pathvars._sfw = sbw0, sfw0
    while True:
        # Reduce MEP to current limits
        sbwi = pathvars._sbw
        sfwi = pathvars._sfw
        fncs.print_string("--> CALCULATING COEFFICIENTS IN [%.4f,%.4f] (in bohr)<--"%(sbwi,sfwi),4)
        drst = {label:data for label,data in drst.items() if in_interval(data[0],sbwi,sfwi)}
        # Calculate corresponding coefficients
        dcfs,pathvars,palpha,pomega = calc_coefs(itarget,tcommon,drst,pathvars,ltemp,symmetry,plotfile)
        # expand data
        try   : SCT = dcfs["sct"]
        except: return dcfs
        # Save data only for lower temperature
        convlist_sct.append(SCT[0])
        convlist_lim.append((sbwi,sfwi))
        # Check convergence
        if len(convlist_sct) > 1:
           SCT_a  = convlist_sct[-1]
           SCT_b  = convlist_sct[-2]
           dif100 = 100*abs(SCT_a-SCT_b)/SCT_a
           # reduce value only to two decimal places
           dif100 = float("%9.2f"%dif100)
           # compare
           converged = (dif100 <= pathvars._scterr)
           if float("%.4E"%SCT_a) <= float("%.4E"%SCT_b): ntimesdecrease += 1
        else: converged = False
        # Be more strict with convergence if no products
        if not bool_prods and ntimesdecrease < 2: converged = False
        # Print convergence table
        tuple_table = (convlist_sct,convlist_lim,pathvars._scterr,converged)
        fncs.print_string(PS.ssct_convergence(*tuple_table),4)
        if converged: break
        # If current limits are outside required, break
        if fncs.is_smallereq(sbwi,sbwN,EPS_MEPS) or fncs.is_greatereq(sfwi,sfwN,EPS_MEPS):
           fncs.print_string("WARNING: SCT transmission coefficient is NOT converged YET!",4)
           fncs.print_string("         Unfortunately, the MEP limit defined by the user",4)
           fncs.print_string("         was reached...",4)
           fncs.print_string("         [sbw,sfw] = [%.4f,%.4f]"%(sbwN,sfwN),4)
           print("")
           break
        # Increase MEP
        sbwi = pathvars._sbw
        sfwi = pathvars._sfw
        V1bw = palpha[1]+pathvars._eref
        V1fw = pomega[1]+pathvars._eref
        lbw,lfw,sbw,sfw,V0bw,V0fw = sd.rstlimits(drst)
        pathvars.increase_svals(V0bw,V0fw,V1bw,V1fw)
        sbwj = pathvars._sbw
        sfwj = pathvars._sfw
        fncs.print_string("--> MEP will be calculated between [%.4f,%.4f] (in bohr) <--"%(sbwj,sfwj),4)
        print("")
        fncs.print_string("sbw: %+8.4f --> %+8.4f bohr"%(sbwi,sbwj),8)
        fncs.print_string("sfw: %+8.4f --> %+8.4f bohr"%(sfwi,sfwj),8)
        print("")
        tcommon,drst,pathvars = calc_mep(itarget,gtsTS,pathvars,tsoftware,TMP)
        pathvars._sbw = sbwj
        pathvars._sfw = sfwj
    del drst
    return dcfs, converged
#---------------------------------------------------------------#
def deal_with_path(target,dlevel,software,ltemp,dctc,pathvars,dtes,dchem,dhighlvl,dimasses):
    dof      = PN.get_dof(dlevel)
    plotfile = PN.get_plf(dlevel)
    # gts file for this TS
    ctc, itc = PN.name2data(target)
    gtsTS    = dctc[ctc].gtsfile(itc)
    # rotational symmetry
    moleculeTS = Molecule()
    moleculeTS.set_from_gts(gtsTS)
    symmetry = str(moleculeTS._pgroup) , int(moleculeTS._rotsigma)
    # temporal folder
    TMP = PN.TMPi%(target)
    # if exists,remove content (unless keeptmp is activated)
    # (to avoid reading old files from different levels of calculation)
    if os.path.exists(TMP) and not pathvars._keeptmp: shutil.rmtree(TMP,ignore_errors=True)
    # split target
    ctc, itc = PN.name2data(target)
    # name of rst file
    rstfile = PN.get_rst(ctc,itc)
    # tuple software
    tes    = dtes.get(software,{}).get(ctc,None)
    tsoftw = (software,tes)
    # internal coordinates
    if   itc in dctc[ctc]._dics.keys(): ics = dctc[ctc]._dics[itc]
    elif "*" in dctc[ctc]._dics.keys(): ics = dctc[ctc]._dics["*"]
    else                              : ics = None
    if   itc in dctc[ctc]._dicsbw.keys(): icsbw = dctc[ctc]._dicsbw[itc]
    elif "*" in dctc[ctc]._dicsbw.keys(): icsbw = dctc[ctc]._dicsbw["*"]
    else                                : icsbw = None
    if   itc in dctc[ctc]._dicsfw.keys(): icsfw = dctc[ctc]._dicsfw[itc]
    elif "*" in dctc[ctc]._dicsfw.keys(): icsfw = dctc[ctc]._dicsfw["*"]
    else                                : icsfw = None
    # path variables
    pathvars.set_ics(ics,icsbw,icsfw) # before setup3!!
    pathvars.apply_specific(itc)
    pathvars.setup1()
    pathvars.setup2()
    pathvars.setup3()
    # Set Eref (from reaction)
    pathvars.set_eref_from_reaction(target,dchem,dof)
    # Quantum reaction coordinate qrc
    pathvars.prepare_qrc(dchem,dctc,dimasses)
    # frequency scaling factor
    pathvars._freqscal = float(dctc[ctc]._fscal)
    # if dlevel --> no convergence and dlevel data
    if dlevel:
       exception = Exc.NoDLEVELdata(Exception)
       pathvars._scterr = None
       keydhl = "%s.%s.path"%(ctc,itc)
       if keydhl not in dhighlvl.keys():
          # maybe only TS
          keydhl = "%s.%s"%(dctc[ctc]._root,itc)
          if keydhl in dhighlvl.keys():
             dictE  = {0.0:dhighlvl[keydhl]}
          else:
            global WARNINGS
            WARNINGS.append("No high-level data for %s"%target)
            raise exception
       else: dictE = dhighlvl[keydhl]
       pathvars._dlevel = dictE
    # LOGGER
    pof = PN.get_pof(dlevel,"path",target)
    sys.stdout = Logger(pof,"w",True)
    sys.stdout.writeinfile(PS.init_txt())
    #string
    fncs.print_string(PS.smep_title(target,pathvars,pof),2)
    # Was previously converged???
    ffloat = "%.3f"
    drstconv = RW.read_rstconv()
    if target in drstconv.keys() and not dlevel and os.path.exists(rstfile):
       lowtemp, useics, scterr, converged = drstconv[target]
       b0 = (converged == "yes")
       b1 = (pathvars._useics == useics)
       b2 = (ffloat%pathvars._scterr == ffloat%scterr)
       b3 = (ffloat%min(ltemp) == ffloat%lowtemp)
       if b0 and b1 and b2 and b3:
          pathvars._scterr = None
          fncs.print_string("THIS PATH IS STORED AS CONVERGED!\n",4)
          tpath, tcommon, drst = ff.read_rst(rstfile)
          lbw,lfw,sbw,sfw,V0bw,V0fw = sd.rstlimits(drst)
          pathvars._sbw = sbw
          pathvars._sfw = sfw
          del drst
    #----------------#
    # calculate path #
    #----------------#

    # 1. Only MEP
    if not pathvars._beyondmep:
       common,drst,pathvars = calc_mep(target,gtsTS,pathvars,tsoftw,TMP)
       dcoefs = {}
       del drst
       # raise error
       raise Exc.OnlyMEP(Exception)

    # 2. MEP expanded till SCT convergence
    elif pathvars.sct_convergence():
         dcoefs, converged = get_path_sctconv(target,gtsTS,pathvars,tsoftw,ltemp,TMP,symmetry,plotfile)
         # save convergence in file!
         drstconv = RW.read_rstconv()
         if converged: drstconv[target] = (min(ltemp),pathvars._useics,pathvars._scterr,"yes")
         else        : drstconv[target] = (min(ltemp),pathvars._useics,pathvars._scterr,"no")
         RW.write_rstconv(drstconv)

    # 3. Coefs with the current MEP extension
    else:
         tcommon,drst,pathvars = calc_mep(target,gtsTS,pathvars,tsoftw,TMP)
         dcoefs,pathvars,palpha,pomega = calc_coefs(target,tcommon,drst,pathvars,ltemp,symmetry,plotfile)
         del drst

    # print summary with the coefficients
    fncs.print_string(PS.spath_allcoefs(ltemp,dcoefs),3)
    # return data
    return dcoefs, pathvars
#===============================================================#

#===============================================================#
def execute_pfn(itargets,dchem,idata,status,case):
    reactions = set([])
    for target in itargets:
        ctc1, itc1 = PN.name2data(target)
        the_reaction = None
        for reaction,(Rs,TS,Ps) in dchem.items():
            try   : ctc2, itc2 = PN.name2data(TS)
            except: continue
            if ctc1 == ctc2:
               the_reaction = reaction
               if itc1 == itc2: break
        if the_reaction is not None: reactions.add(the_reaction)
    if len(reactions) == 0: return
    reactions = sorted(list(reactions))

    fncs.print_string("The selected transitions states are involved in defined reactions!",3)
    pfn_targets = set([])
    for reaction in reactions:
        fncs.print_string("* %s"%reaction,6)
        Rs,TS,Ps = dchem[reaction]
        for xx in Rs+[TS]+Ps: pfn_targets.add(PN.name2data(xx)[0])
    if len(pfn_targets) == 0: return
    pfn_targets = sorted(list(pfn_targets))
    print("")

    import opt_pfn as pfn
    fncs.print_string("Calculating partition functions for the next targets:",3)
    for target in pfn_targets:
        fncs.print_string("* %s"%target,6)
    print("")
    pfn.main(idata,status,case,targets=pfn_targets)
#===============================================================#

#===============================================================#
def main(idata,status,case,targets="*"):

    stat2check = [1,2,3,4]
    mustexist  = [PN.DIR1]
    tocreate   = [PN.DIR4,PN.DIR2,PN.DIR3,PN.DIR5,PN.DIR6,PN.TMP]

    #-------------------------------------------------------#
    # Read Pilgrim input files, check file/folder status    #
    # and expand tuple 'case'                               #
    #-------------------------------------------------------#
    # expand data
    (dctc,dimasses), ltemp, dpath, (dtesLL,dtesHL), dchem, tkmc, ddlevel = idata
    # status ok?
    fstatus = status_check(status,stat2check)
    if fstatus == -1: exit()
    # existency of folders
    fstatus = ffchecking(mustexist,tocreate)
    if fstatus == -1: exit()
    # expand case
    (dof,hlf,plotfile),dlevel,software = case
    #-------------------------------------------------------#


    # print selected software
    fncs.print_string("Selected software: %s"%software,3)
    print("")

    # read high level file
    if dlevel: dhighlvl = RW.read_highlevelfile(hlf)
    else     : dhighlvl = {}

    # update targets
    if targets == [] or targets == "*": targets = dpath.keys()
    itargets = get_itargets(targets,dpath,dctc)

    # Print targets
    if len(itargets) != 0:
       fncs.print_string("The paths will be calculated for:",3)
       for idx in range(0,len(itargets),4):
           fncs.print_string("%s"%(", ".join(itargets[idx:idx+4])),7)
    else: fncs.print_string("No target fits for the selected options...",3)
    print("")

    # check if TSs are involved in any reaction
#   try   : execute_pfn(itargets,dchem,idata,status,case)
#   except: pass

    # loop over each target
    dall = RW.read_alldata(dof,ltemp)[0]
    for target in itargets:
        # collect garbage
        gc.collect()
        t1 = time.time()
        # Get ctc and itc from name
        ctc, itc = PN.name2data(target)
        # initialize pathvars
        pathvars = copy.deepcopy(dpath[ctc])
        # isotopic modification?
        masses = get_masses(target,dctc,dimasses)
        pathvars.set_masses(masses)
        # calculate path + coefficients
        idata = (target,dlevel,software,ltemp,dctc,pathvars,dtesLL,dchem,dhighlvl,dimasses)
        try: dcoefs, pathvars = deal_with_path(*idata)
        except Exception as exception:
           deal_with_exception(exception)
           sys.stdout = Logger(None)
           continue
        # print end of file
        t2    = time.time()
        etime = fncs.time2human(t2-t1,"secs")
        sys.stdout.writeinfile(PS.end_txt(*etime))
        # update file (read, check ltemp, update, write)
        sys.stdout = Logger(None)
        fncs.print_string("Updating data file: %s"%dof,3)
        dall = RW.read_alldata(dof,ltemp)[0]
        for coef in dcoefs.keys(): dall[coef][target] = dcoefs[coef]
        RW.write_alldata(dof,ltemp,dall)
        print("")
        # collect garbage
        gc.collect()

    # print WARNINGS
    if len(WARNINGS) != 0:
       fncs.print_string("WARNINGS:",3)
       for warning in WARNINGS: print_string("      *",warning)
       print("")
#===============================================================#
