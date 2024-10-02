#!/usr/bin/env python
# vim: set fileencoding=UTF-8 :

"""Various flux derivative stuff."""

import numpy as np
import math
import scipy.interpolate
import sys
import matplotlib.pyplot as plt
import matplotlib.backends
import re

class DataError(Exception):
    """Do some thing."""
    
    def __init__(self, value):
        """Define some value."""
        self.value = value
    def __str__(self):
        """Define some value."""
        return repr(self.value)

def wheref(array, value):
    """Find where a floating point value is equal to another."""
    #Floating point inaccuracy.
    eps=1e-7
    return np.where((array > value-eps)*(array < value+eps))

def rebin(data, xaxis,newx):
    """Just rebins the data."""
    if newx[0] < xaxis[0] or newx[-1]> xaxis[-1]:
            raise ValueError("A value in newx is beyond the interpolation range")
    intp=scipy.interpolate.InterpolatedUnivariateSpline(np.log(xaxis),data)
    newdata=intp(np.log(newx))
    return newdata

def save_figure(path):
    """Save the figure, automatically determining file extension."""
    bk=matplotlib.backends.backend
    if path == "":
            return
    elif bk == 'TkAgg' or bk == 'Agg' or bk == 'GTKAgg':
            path = path+".png"
    elif bk == 'PDF' or bk == 'pdf':
            path = path+".pdf"
    elif bk == 'PS' or bk == 'ps':
            path = path+".ps"
    return plt.savefig(path)

def corr_table(table, dvecs,table_name):
    """Little function to adjust a table so it has a different central value."""
    new=np.array(table)
    new[12:,:] = table[12:,:]+2*table[0:12,:]*dvecs
    pkd="/home/spb41/cosmomc-src/cosmomc/data/lya-interp/"
    np.savetxt(pkd+table_name,new,("%1.3g","%1.3g","%1.3g","%1.3g","%1.3g","%1.3g","%1.3g","%1.3g","%1.3g","%1.3g","%1.3g"))
    return new


class PowerSpec:
    """A class to be derived from by flux and matter PowerSpec classes. Stores various helper methods."""
    
    # snapshots
    snaps=()
    #SDSS redshift bins.
    zz=np.array([1.0])
    #SDSS kbins, in s/km units.
    sdsskbins=np.array([1.0])
    # Omega_matter
    om=0.267
    #hubble constant
    h0=0.71
    #Boxsize in Mpc/h
    box=60.0
    #Size of the best-fit box, for testing varying boxsize
    bfbox=60.0
    #Some paths
    knotpos=np.array([1.0])
    base=""
    pre=""
    suf=""
    ext=""
    #For plotting
    ymin=0.8
    ymax=1.2
    figprefix="/figure"
    def __init__(self, snaps=("snapshot_000", "snapshot_001","snapshot_002","snapshot_003","snapshot_004","snapshot_005","snapshot_006","snapshot_007","snapshot_008","snapshot_009","snapshot_010","snapshot_011"),
             zz=np.array([4.2,4.0,3.8,3.6,3.4,3.2,3.0,2.8,2.6,2.4,2.2,2.0]),
             sdsskbins=np.array([0.00141,0.00178,0.00224,0.00282,0.00355,0.00447,0.00562,0.00708,0.00891,0.01122,0.01413,0.01778]),
             knotpos=np.array([0.07,0.15,0.475, 0.75, 1.19, 1.89,4,25]), om=0.267, h0=0.71,box=60.0,
             base="/home/spb41/Lyman-alpha/MinParametricRecon/runs/",suf="/", ext=".txt"):
        """Plot the snapshots."""
        if len(snaps) != np.size(zz):
            raise DataError("There are "+str(len(snaps))+" snapshots, but "+str(np.size(zz))+"redshifts given.")
        self.snaps=snaps
        self.zz=zz
        self.sdsskbins=sdsskbins
        self.om=om
        self.h0=h0
        self.box=box
        self.knotpos=knotpos*h0
        self.base=base
        self.suf=suf
        self.ext=ext
        return

    def getz(self, snap):
        """Get redshift associated with a snapshot."""
        ind=np.where(np.array(self.snaps) == snap)
        if np.size(ind):
            return self.zz[ind]
        else:
            raise DataError(str(snap)+" does not exist!")

    def getsnap(self, redshift):
        """Get snapshot associated with a redshift."""
        ind=wheref(self.zz, redshift)
        if np.size(ind):
            return str(np.asarray(self.snaps)[ind][0])
        else:
            raise DataError("No snapshot at redshift "+str(redshift))

    def getsdsskbins(self, redshift):
        """Get the k bins at a given redshift in h/Mpc units."""
        return self.sdsskbins*self.hubble(redshift)/(1.0+redshift)
        #Corr is /sqrt(self.h0)...

    def hubble(self, zz):
        """Hubble parameter."""
        return 100*self.h0*math.sqrt(self.om*(1+zz)**3+(1-self.om))
        #Conversion factor between s/km and h/Mpc is (1+z)/H(z)

    def loaddata(self, file, box):
        """Do correct units conversion to return k and one-d power."""
        #Adjust Fourier convention.
        fluxpower=np.loadtxt(file)
        scale=self.h0/box
        k=fluxpower[1:,0]*scale*2.0*math.pi
        pf=fluxpower[1:,1]/scale
        return (k, pf)

    def plot_z(self,knot,redshift,title="",ylabel="", legend=True):
        """Plot comparisons between a bunch of sims on one graph."""
        """plot_z(Redshift, sims to use ( eg, A1.14)."""
        """Note this will clear current figures."""
        (simk,bfpk) = self.loadpk(knot.bstft + self.suf + self.pre + self.getsnap(redshift) + self.ext,self.bfbox)
        #Setup figure plot.
        ind=wheref(self.zz, redshift)
        plt.figure(ind[0][0])
        plt.clf()
        if title != '':
            plt.title(title+" at z="+str(redshift))
        plt.ylabel(ylabel)
        plt.xlabel(r"$k\; (\mathrm{Mpc}^{-1})$")
        line=np.array([])
        legname=np.array([])
        for sim in knot.names:
            (k,pk)=self.loadpk(sim+self.suf+self.pre+self.getsnap(redshift)+self.ext,self.box)
            oi = np.where(simk <= k[-1])
            ti = np.where(simk[oi] >= k[0])
            relp=rebin(pk, k, simk[oi][ti])
            relp=relp/rebin(bfpk, simk, simk[oi][ti])
            line=np.append(line, plt.semilogx(simk[oi][ti]/self.h0,relp,linestyle="-"))
            legname=np.append(legname,sim)
        if legend:
            plt.legend(line,legname)
            plt.semilogx(self.knotpos,np.ones(len(self.knotpos)),"ro")
        plt.ylim(self.ymin,self.ymax)
        plt.xlim(simk[0]*0.8, 10)
        return

    def plot_all(self, knot,zzz=np.array([]), out=""):
        """Plot a whole suite of snapshots."""
        if np.size(zzz) == 0:
            zzz=self.zz    #lolz
        for z in zzz:
            self.plot_z(knot,z)
            if out != "":
                save_figure(out+self.figprefix+str(z))
        return

    def plot_power(self,path, redshift,colour="black"):
        """Plot absolute power spectrum, not relative."""
        (k_g,pk_g) = self.loadpk(path + self.suf + self.pre + self.getsnap(redshift) + self.ext,self.box)
        plt.loglog(k_g,pk_g, color=colour)
        plt.xlim(0.01,k_g[-1]*1.1)
        plt.ylabel("P(k) /(h-3 Mpc3)")
        plt.xlabel("k /(h MPc-1)")
        plt.title("Power spectrum at z="+str(redshift))
        return(k_g, pk_g)

    def plot_power_all(self, knot,zzz=np.array([]), out=""):
        """Plot absolute power for all redshifts."""
        if np.size(zzz) == 0:
            zzz=self.zz    #lolz
        for z in zzz:
            ind=wheref(self.zz, z)
            plt.figure(ind[0][0])
            for sim in knot.names:
                self.plot_power(sim,z)
            if out != "":
                save_figure(out+self.figprefix+str(z))
        return

    def plot_compare_two(self, one, onebox, two,twobox,colour=""):
        """Compare two power spectra directly. Smooths result."""
        """plot_compare_two(first P(k), second P(k))."""
        (onek,onep)=self.loadpk(one,onebox)
        (twok,twop)=self.loadpk(two,twobox)
        onei = np.where(onek <= twok[-1])
        twoi= np.where (onek[onei] >= twok[0])
        relp=rebin(twop, twok, onek[onei][twoi])
        relp=relp/rebin(onep, onek, onek[onei][twoi])
        onek=onek[onei][twoi]
        plt.title("Relative Power spectra "+one+" and "+two)
        plt.ylabel(r"$P_2(k)/P_1(k)$")
        plt.xlabel(r"$k\; (h\,\mathrm{Mpc}^{-1})$")
        if colour == "":
            line=plt.semilogx(onek,relp)
        else:
            line=plt.semilogx(onek,relp,color=colour)
        plt.semilogx(self.knotpos,np.ones(len(self.knotpos)),"ro")
        ind=np.where(onek < 10)
        plt.ylim(min(relp[ind])*0.98,max(relp[ind])*1.01)
        plt.xlim(onek[0]*0.8, 10)
        return line

    def compare_two(self, one, two,redshift):
        """Get the difference between two simulations on scales probed by the SDSS power spectrum."""
        (onek,onep)=self.loadpk(one,self.bfbox)
        (twok,twop)=self.loadpk(two,self.box)
        onei = np.where(onek <= twok[-1])
        twoi= np.where (onek[onei] >= twok[0])
        relp=rebin(twop, twok, onek[onei][twoi])
        relp=relp/rebin(onep, onek, onek[onei][twoi])
        onek=onek[onei][twoi]
        sdss=self.getsdsskbins(redshift)
        relp_r=np.ones(np.size(sdss))
        ind = np.where(sdss > onek[0])
        relp_r[ind]=rebin(relp,onek,sdss[ind])
        return relp_r

    def compare_two_table(self,onedir, twodir):
        """Do the above 12 times to get a correction table."""
        nk=np.size(self.sdsskbins)
        nz=np.size(self.zz)-1
        table=np.empty([nz,nk])
        for i in np.arange(0,nz):
            sim=self.pre+self.getsnap(self.zz[i])+self.ext
            table[-1-i,:]=self.compare_two(onedir+sim,twodir+sim,self.zz[i])
        return table

    def plot_compare_two_sdss(self, onedir, twodir, zzz=np.array([]), out="", title="", ylabel="", ymax=0,ymin=0, colour="",legend=False):
        """Plot a whole redshift range of relative power spectra on the same figure."""
        """plot_all(onedir, twodir)."""
        """Pass onedir and twodir as relative to basedir. ie, for default settings something like best-fit/flux-power/."""
        if np.size(zzz) == 0:
            zzz=self.zz    #lolz
        line=np.array([])
        legname=np.array([])
        sdss=self.sdsskbins
        plt.xlabel(r"$k_v\; (s\,\mathrm{km}^{-1})$")
        for z in zzz:
            sim=self.pre+self.getsnap(z)+self.ext
            pk=self.compare_two(onedir+sim,twodir+sim,z)
            if colour == "":
                line=np.append(line, plt.semilogx(sdss,pk))
            else:
                line=np.append(line, plt.semilogx(sdss,pk,color=colour))
            legname=np.append(legname,"z="+str(z))
        plt.title(title)
        if ylabel != "":
            plt.ylabel(ylabel)
        if legend:
            plt.legend(line, legname,bbox_to_anchor=(0., 0, 1., .25), loc=3,ncol=3, mode="expand", borderaxespad=0.)
        if ymax != 0 and ymin !=0:
            plt.ylim(ymin,ymax)
        plt.xlim(sdss[0],sdss[-1])
        plt.xticks(np.array([sdss[0],3e-3,5e-3,0.01,sdss[-1]]),("0.0014","0.003","0.005","0.01","0.0178"))
        if out != "":
            save_figure(out)
        return plt.gcf()

    def plot_compare_two_all(self, onedir, twodir, zzz=np.array([]), out="", title="", ylabel="", ymax=0,ymin=0, colour="",legend=False):
        """Plot a whole redshift range of relative power spectra on the same figure."""
        """plot_all(onedir, twodir)"""
        """Pass onedir and twodir as relative to basedir. ie, for default settings something like best-fit/flux-power/ onedir uses bfbox, twodir uses box."""
        if np.size(zzz) == 0:
            zzz=self.zz    #lolz
        line=np.array([])
        legname=np.array([])
        for z in zzz:
            line=np.append(line, self.plot_compare_two(onedir+self.pre+self.getsnap(z)+self.ext,self.bfbox,twodir+self.pre+self.getsnap(z)+self.ext,self.box,colour))
            legname=np.append(legname,"z="+str(z))
        if title == "":
            plt.title("Relative Power spectra "+onedir+" and "+twodir)
        else:
            plt.title(title)
        if ylabel != "":
            plt.ylabel(ylabel)
        if legend:
            plt.legend(line, legname,bbox_to_anchor=(0., 0, 1., .25), loc=3,ncol=3, mode="expand", borderaxespad=0.)
        if ymax != 0 and ymin !=0:
            plt.ylim(ymin,ymax)
        if out != "":
            save_figure(out)
        return plt.gcf()

    def getflat(self,dir, si=0.0):
        """Get a power spectrum in the flat format we use."""
        """for inputting some cosmomc tables."""
        pk_sdss=np.empty([11, 12])
        #Note this omits the z=2.0 bin
        #SiIII corr now done on the fly in lya_sdss_viel.f90
        for i in np.arange(0,np.size(self.snaps)-1):
            scale=self.hubble(self.zz[i])/(1.0+self.zz[i])
            (k,pk)=self.loadpk(dir+self.snaps[i]+self.ext, self.box)
            fbar=math.exp(-0.0023*(1+self.zz[i])**3.65)
            a=si/(1-fbar)
            #The SiIII correction is kind of oscillatory, so we want
            #to average over the whole interval being probed.
            sdss=self.sdsskbins
            kmids=np.zeros(np.size(sdss)+1)
            sicorr=np.empty(np.size(sdss))
            for j in np.arange(0,np.size(sdss)-1):
                kmids[j+1]=math.exp((math.log(sdss[j+1])+math.log(sdss[j]))/2.)
            #Final segment should make no difference
            kmids[-1]=2*math.pi/2271+kmids[-2]
            #Near zero bad things will happen
            sicorr[0]=1+a**2+2*a*math.cos(2271*sdss[0])
            for j in np.arange(1,np.size(sdss)):
                sicorr[j]=1+a**2+2*a*(math.sin(2271*kmids[j+1])-math.sin(2271*kmids[j]))/(kmids[j+1]-kmids[j])/2271
            sdss=self.getsdsskbins(self.zz[i])
            pk_sdss[-i-1,:]=rebin(pk, k, sdss)*scale*sicorr
        return pk_sdss

    def calc_z(self, redshift,s_knot, kbins):
        """Calculate the flux derivatives for a single redshift. Output: (kbins d2P...kbins dP (flat vector of length 2xkbins))."""
        #Array to store answers.
        #Format is: k x (dP, d²P, χ²)
        kbins=np.array(kbins)
        nk=np.size(kbins)
        snap=self.getsnap(redshift)
        results=np.zeros(2*nk)
        if np.size(s_knot.qvals) > 0:
            results = np.zeros(4*nk)
        pdifs=s_knot.pvals-s_knot.p0
        qdifs=np.array([])
        #If we have somethign where the parameter is redshift-dependent, eg, gamma.
        if np.size(s_knot.p0) > 1:
            i=wheref(redshift, self.zz)
            pdifs = s_knot.pvals[:,i]-s_knot.p0[i]
            if np.size(s_knot.qvals) > 0:
                qdifs=s_knot.qvals[:,i] - s_knot.q0[i]
        npvals=np.size(pdifs)
        #Load the data
        (k,pfp0)=self.loadpk(s_knot.bstft+self.suf+snap+self.ext,s_knot.bfbox)
        #This is to rescale by the mean flux, for generating mean flux tables.
        ###
        #tau_eff=0.0023*(1+redshift)**3.65
        #tmin=0.2*((1+redshift)/4.)**2
        #tmax=0.5*((1+redshift)/4.)**4
        #teffs=tmin+s_knot.pvals*(tmax-tmin)/30.
        #pdifs=teffs/tau_eff-1.
        ###
        powerfluxes=np.zeros((npvals,np.size(k)))
        for i in np.arange(0,np.size(s_knot.names)):
            (k,powerfluxes[i,:])=self.loadpk(s_knot.names[i]+self.suf+snap+self.ext, s_knot.bfbox)
        #So now we have an array of data values, which we want to rebin.
        ind = np.where(kbins >= k[0])
        difpf_rebin=np.ones((npvals,np.size(kbins)))
        for i in np.arange(0, npvals):
            difpf_rebin[i,ind]=rebin(powerfluxes[i,:]/pfp0,k,kbins[ind])
            #Set the things beyond the range of the interpolator
            #equal to the final value.
            if ind[0][0] > 0:
                difpf_rebin[i,0:ind[0][0]]=difpf_rebin[i,ind[0][0]]
        #So now we have an array of data values.
        #Pass each k value to flux_deriv in turn.
        for k in np.arange(0,np.size(kbins)):
            #Format of returned data is:
            # y = ax**2 + bx + cz**2 +dz +e xz
            # derives = (a,b,c,d,e)
            derivs=self.flux_deriv(difpf_rebin[:,k], pdifs,qdifs)
            results[k]=derivs[0]
            results[nk+k]=derivs[1]
            if np.size(derivs) > 2:
                results[2*nk+k]=derivs[2]
                results[3*nk+k]=derivs[3]
        return results

    def calc_all(self, s_knot,kbins):
        """Calculate the flux derivatives for all redshifts. Input: sims to load, parameter values, mean parameter value. Output: (2*kbins) x (zbins)."""
        flux_derivatives=np.zeros((2*np.size(kbins),np.size(self.zz)))
        if np.size(s_knot.qvals) > 1:
            flux_derivatives=np.zeros((4*np.size(kbins),np.size(self.zz)))
        #Call flux_deriv_const_z for each redshift.
        for i in np.arange(0,np.size(self.zz)):
            flux_derivatives[:,i]=self.calc_z(self.zz[i], s_knot,kbins)
        return flux_derivatives

    def flux_deriv(self, pfdif, pdif, qdif=np.array([])):
        """Calculate the flux-derivative for a single redshift and k bin."""
        pdif=np.ravel(pdif)
        if np.size(pdif) != np.size(pfdif):
            raise DataError(str(np.size(pdif))+" parameter values, but "+str(np.size(pfdif))+" P_F values")
        if np.size(pdif) < 2:
            raise DataError(str(np.size(pdif))+" pvals given. Need at least 2.")
        pfdif=pfdif-1.0
        if np.size(qdif) > 2:
            qdif=np.ravel(qdif)
            mat=np.vstack([pdif**2, pdif, qdif**2, qdif] ).T
        else:
            mat=np.vstack([pdif**2, pdif] ).T
        (derivs, residues,rank, sing)=np.linalg.lstsq(mat, pfdif)
        return derivs

    def get_error_z(self, sim, bstft,box, derivs, params, redshift,qarams=np.empty([])):
        """Get the error on one test simulation at single redshift."""
        #Need to load and rebin the sim.
        (k, test)=self.loadpk(sim+self.suf+self.getsnap(redshift)+self.ext, box)
        (k,bf)=self.loadpk(bstft+self.suf+self.getsnap(redshift)+self.ext,box)
        kbins=self.getkbins()
        ind = np.where(kbins >= 1.0*2.0*math.pi*self.h0/box)
        test2=np.ones(np.size(kbins))
        test2[ind]=rebin(test/bf,k,kbins[ind])#
        if ind[0][0] > 0:
            test2[0:ind[0][0]]=test2[ind[0][0]]
        if np.size(qarams) > 0:
            guess=derivs.Getpf(params, redshift,qarams)+1.0
        else:
            guess=derivs.Getpf(params, redshift)+1.0
        return np.array((test2/guess))[0][0]

class FluxPow(PowerSpec):
    """A class written to store the various methods related to calculating of the flux derivatives and plotting of the flux power spectra."""
    
    figprefix="/flux-figure"
    kbins=np.array([])
    def __init__(self, snaps=("snapshot_000", "snapshot_001","snapshot_002","snapshot_003","snapshot_004","snapshot_005","snapshot_006","snapshot_007","snapshot_008","snapshot_009","snapshot_010","snapshot_011"),
         zz=np.array([4.2,4.0,3.8,3.6,3.4,3.2,3.0,2.8,2.6,2.4,2.2,2.0]),
         sdsskbins=np.array([0.00141,0.00178,0.00224,0.00282,0.00355,0.00447,0.00562,0.00708,0.00891,0.01122,0.01413,0.01778]),
         knotpos=np.array([0.07,0.15,0.475, 0.75, 1.19, 1.89,4,25]), om=0.266, h0=0.71,box=60.0,kmax=4.0,
         base="/home/spb41/Lyman-alpha/MinParametricRecon/runs/",bf="best-fit/",suf="flux-power/", ext="_fluxpower.txt"):
        """Do something with the snapshots."""
        PowerSpec.__init__(self, snaps,zz,sdsskbins,knotpos, om, h0,box,base,suf, ext)
        (k_bf,pk_bf)= self.loadpk(bf+suf+"snapshot_000"+self.ext,self.bfbox)
        ind=np.where(k_bf <= kmax)
        self.kbins=k_bf[ind]

    def plot_z(self,sims,redshift,title="Relative Flux Power",ylabel=r"$\mathrm{P}_\mathrm{F}(k,p)\,/\,\mathrm{P}_\mathrm{F}(k,p_0)$", legend=True):
        """Plot comparisons between a bunch of sims on one graph."""
        PowerSpec.plot_z(self,sims,redshift,title,ylabel,legend)
        if legend:
                kbins=self.getsdsskbins(redshift)
                plt.axvspan(kbins[0], kbins[-1], color="#B0B0B0")
        plt.ylim(self.ymin,self.ymax)
        plt.xlim(self.kbins[0]*0.8, 10)

    def getkbins(self):
        """Get the kbins to interpolate onto."""
        return self.kbins

    def macdonaldpf(self,sdss, zz):
        """Load the SDSS power spectrum."""
        psdss=sdss[np.where(sdss[:,0] == zz)][:,1:3]
        fbar=math.exp(-0.0023*(1+zz)**3.65)
        #multiply by the hubble parameter to be in 1/(km/s)
        scale=self.hubble(zz)/(1.0+zz)
        pf=psdss[:,1]*fbar**2/scale
        k=psdss[:,0]*scale
        return (k, pf)

    def loadpk(self, path,box):
        """Load a pk. Different function due to needing to be different for each class."""
        #Adjust Fourier convention.
        fluxpower=np.loadtxt(self.base+path)
        scale=self.h0/box
        k=(fluxpower[1:,0]-0.5)*scale*2.0*math.pi
        pf=fluxpower[1:,1]/scale
        return (k, pf)

class MatterPow(PowerSpec):
    """A class to plot matter power spectra."""
    
    ob=0.0
    #For plotting
    ymin=0.4
    ymax=1.6
    figprefix="/matter-figure"
    def __init__(self, snaps=("snapshot_000", "snapshot_001","snapshot_002","snapshot_003","snapshot_004","snapshot_005","snapshot_006","snapshot_007","snapshot_008","snapshot_009","snapshot_010","snapshot_011"),
         zz=np.array([4.2,4.0,3.8,3.6,3.4,3.2,3.0,2.8,2.6,2.4,2.2,2.0]),
         sdsskbins=np.array([0.00141,0.00178,0.00224,0.00282,0.00355,0.00447,0.00562,0.00708,0.00891,0.01122,0.01413,0.01778]),
         knotpos=np.array([0.07,0.15,0.475, 0.75, 1.19, 1.89,4,25]), om=0.266,ob=0.0449, h0=0.71,box=60.0,
         base="/home/spb41/Lyman-alpha/MinParametricRecon/runs/",suf="matter-power/", ext=".0", matpre="pk-by-"):
        """Do something with the snapshots."""
        PowerSpec.__init__(self, snaps,zz,sdsskbins,knotpos, om, h0,box,base,suf,ext)
        self.ob=ob
        self.pre=matpre

    def plot_z(self,sims,redshift,title="Relative Matter Power",ylabel=r"$\mathrm{P}(k,p)\,/\,\mathrm{P}(k,p_0)$"):
        """Plot comparisons between a bunch of sims on one graph."""
        PowerSpec.plot_z(self,sims,redshift,title,ylabel)

    def loadpk(self, path,box):
        """Load a pk. Different function due to needing to be different for each class."""
        #Load baryon P(k)
        matterpower=np.loadtxt(self.base+path)
        scale=self.h0/box
        #Adjust Fourier convention.
        simk=matterpower[1:,0]*scale*2.0*math.pi
        pkbar=matterpower[1:,1]/scale**3
        #Load DM P(k)
        matterpower=np.loadtxt(self.base+re.sub("by","DM",path))
        pkdm=matterpower[1:,1]/scale**3
        pk=(pkbar*self.ob+pkdm*(self.om-self.ob))/self.om
        return (simk,pk)

    def plot_power(self,path, redshift,camb_filename=""):
        """Plot absolute power spectrum, not relative."""
        (k_g, pk_g)=PowerSpec.plot_power(self,path,redshift)
        sigma=2.0
        pkg=np.loadtxt(self.base+path+self.suf+self.pre+self.getsnap(redshift)+self.ext)
        samp_err=pkg[1:,2]
        sqrt_err=np.array(np.sqrt(samp_err))
        plt.loglog(k_g,pk_g*(1+sigma*(2.0/sqrt_err+1.0/samp_err)),linestyle="-.",color="black")
        plt.loglog(k_g,pk_g*(1-sigma*(2.0/sqrt_err+1.0/samp_err)),linestyle="-.",color="black")
        if camb_filename != "":
                camb=np.loadtxt(camb_filename)
                #Adjust Fourier convention.
                k=camb[:,0]*self.h0
                #NOW THERE IS NO h in the T anywhere.
                pk=camb[:,1]
                plt.loglog(k/self.h0, pk, linestyle="--")
        plt.xlim(0.01,k_g[-1]*1.1)
        return(k_g, pk_g)

class FluxPdf(PowerSpec):
    """The PDF is an instance of the PowerSpec class. Perhaps poor naming."""
    
    def __init__(self, snaps=("snapshot_006","snapshot_007","snapshot_008","snapshot_009","snapshot_010","snapshot_011"),
        zz=np.array([3.0,2.8,2.6,2.4,2.2,2.0]),
        sdsskbins=np.arange(0,20),
        knotpos=np.array([]), om=0.266,ob=0.0449, h0=0.71,box=48.0,
        base="/home/spb41/Lyman-alpha/MinParametricRecon/runs/",suf="flux-pdf/", ext="_fluxpdf.txt",):
        """Do something with the snapshots."""
        PowerSpec.__init__(self, snaps,zz,sdsskbins,knotpos, om, h0,box,base,suf,ext)

    def loadpk(self, path, box):
        """Load a pk. Different function due to needing to be different for each class."""
        fluxpdf = np.loadtxt(self.base+path)
        return(fluxpdf[:,0], fluxpdf[:,1])

    def plot_compare_two(self, one, onebox, two,twobox,colour=""):
        """Compare two power spectra directly. Smooths result."""
        """plot_compare_two(first P(k), second P(k))."""
        (onek,onep)=self.loadpk(one,onebox)
        (twok,twop)=self.loadpk(two,twobox)
        relp=onep/twop
        plt.title("Relative flux PDF "+one+" and "+two)
        plt.ylabel(r"$F_2(k)/F_1(k)$")
        plt.xlabel(r"$Flux$")
        line=plt.semilogy(onek,relp)
        return line

    def plot_power(self,path, redshift):
        """Plot absolute power spectrum, not relative."""
        (k,pdf)=self.loadpk(path+self.suf+self.pre+self.getsnap(redshift)+self.ext,self.box)
        plt.semilogy(k,pdf, color="black", linewidth="1.5")
        plt.ylabel("P(k) /(h-3 Mpc3)")
        plt.xlabel("k /(h MPc-1)")
        plt.title("PDF at z="+str(redshift))
        return(k, pdf)

    def calc_z(self, redshift,s_knot):
        """Calculate the flux derivatives for a single redshift. Output: (kbins d2P...kbins dP (flat vector of length 2x21))."""
        #Array to store answers.
        #Format is: k x (dP, d²P, χ²)
        npvals=np.size(s_knot.pvals)
        nk=21
        results=np.zeros(2*nk)
        pdifs=s_knot.pvals-s_knot.p0
        #This is to rescale by the mean flux, for generating mean flux tables.
        ###
        #tau_eff=0.0023*(1+redshift)**3.65
        #tmin=0.2*((1+redshift)/4.)**2
        #tmax=0.5*((1+redshift)/4.)**4
        #teffs=tmin+s_knot.pvals*(tmax-tmin)/30.
        #pdifs=teffs/tau_eff-1.
        ###
        ured=np.ceil(redshift*5)/5.
        lred=np.floor(redshift*5)/5.
        usnap=self.getsnap(ured)
        lsnap=self.getsnap(lred)
        #Load the data
        (k,upfp0)=self.loadpk(s_knot.bstft+self.suf+usnap+self.ext,s_knot.bfbox)
        upower=np.zeros((npvals,np.size(k)))
        for i in np.arange(0,np.size(s_knot.names)):
            (k,upower[i,:])=self.loadpk(s_knot.names[i]+self.suf+usnap+self.ext, s_knot.bfbox)
        (k,lpfp0)=self.loadpk(s_knot.bstft+self.suf+lsnap+self.ext,s_knot.bfbox)
        lpower=np.zeros((npvals,np.size(k)))
        for i in np.arange(0,np.size(s_knot.names)):
            (k,lpower[i,:])=self.loadpk(s_knot.names[i]+self.suf+lsnap+self.ext, s_knot.bfbox)
        powerfluxes=5*((redshift-lred)*upower+(ured-redshift)*lpower)
        pfp0=5*((redshift-lred)*upfp0+(ured-redshift)*lpfp0)
        #So now we have an array of data values.
        #Pass each k value to flux_deriv in turn.
        for k in np.arange(0,nk):
            (dpf, d2pf,chi2)=self.flux_deriv(powerfluxes[:,k]/pfp0[k], pdifs)
            results[k]=d2pf
            results[nk+k]=dpf
        return results

        def getkbins(self):
            """Get the kbins to interpolate onto."""
            return np.arange(0,20,1)+0.5
        
    def plot_z(self,knot,redshift,title="",ylabel=""):
        """Plot comparisons between a bunch of sims on one graph."""
        """plot_z(Redshift, sims to use ( eg, A1.14)."""
        """Note this will clear current figures."""
        #Load best-fit
        (simk,bfpk)=self.loadpk(knot.bstft+self.suf+self.pre+self.getsnap(redshift)+self.ext,self.bfbox)
        #Setup figure plot.
        ind=wheref(self.zz, redshift)
        plt.figure(ind[0][0])
        plt.clf()
        if title != '':
            plt.title(title+" at z="+str(redshift),)
        plt.ylabel(ylabel)
        plt.xlabel(r"$\mathcal{F}$")
        line=np.array([])
        legname=np.array([])
        for sim in knot.names:
            (k,pk)=self.loadpk(sim+self.suf+self.pre+self.getsnap(redshift)+self.ext,self.box)
            line=np.append(line, plt.semilogy(simk,pk/bfpk,linestyle="-", linewidth=1.5))
            legname=np.append(legname,sim)
        plt.legend(line,legname)
        return

    def getflat(self,dir):
        """Get a power spectrum in the flat format we use."""
        """for inputting some cosmomc tables."""
        # pk_sdss=np.empty([11, 12])
        #For z=2.07 we need to average snap_011 and snap_010
        z=2.07
        (k,pf_a)=self.loadpk(dir+self.suf+"snapshot_011"+self.ext, self.box)
        (k,pf_b)=self.loadpk(dir+self.suf+"snapshot_010"+self.ext, self.box)
        pf1=(z-2.0)*5*(pf_b-pf_a)+pf_a
        z=2.52
        (k,pf_a)=self.loadpk(dir+self.suf+"snapshot_009"+self.ext, self.box)
        (k,pf_b)=self.loadpk(dir+self.suf+"snapshot_008"+self.ext, self.box)
        pf2=(z-2.4)*5*(pf_b-pf_a)+pf_a
        z=2.94
        (k,pf_a)=self.loadpk(dir+self.suf+"snapshot_007"+self.ext, self.box)
        (k,pf_b)=self.loadpk(dir+self.suf+"snapshot_006"+self.ext, self.box)
        pf3=(z-2.8)*5*(pf_b-pf_a)+pf_a
        pdf = np.array([pf1,pf2,pf3])
        np.savetxt(sys.stdout,pdf.T("%1.8f","%1.8f","%1.8f"))
        return (pf1, pf2, pf3)

if __name__=='__main__':
        # flux=FluxPow()
        # matter=MatterPow()
        # fpdf=fluxpdf()
        # A_knot=knot(("A0.54/","A0.74/","A0.84/","A1.04/", "A1.14/","A1.34/"), (0.54,0.74,0.84,1.04,1.14,1.34),0.94,"best-fit/", 60)
        # AA_knot=knot(("AA0.54/","AA0.74/","AA1.14/","AA1.34/"), (0.54,0.74,1.14,1.34),0.94,"boxcorr400/", 120)
        # B_knot=knot(("B0.33/","B0.53/","B0.73/", "B0.83/", "B1.03/","B1.13/", "B1.33/"), (0.33,0.53,0.73,0.83,1.03, 1.13,1.33),0.93,"best-fit/", 60)
        # C_knot=knot(("C0.11/", "C0.31/","C0.51/","C0.71/","C1.11/","C1.31/","C1.51/"),(0.11, 0.31,0.51,0.71, 1.11,1.31,1.51),0.91,"bf2/", 60)
        # D_knot=knot(("D0.50/","D0.70/","D1.10/","D1.20/", "D1.30/", "D1.50/", "D1.70/"),(0.50, 0.70,1.10,1.20,1.30, 1.50, 1.70),0.90,"bfD/", 48)
        # interp=flux_interp(flux, (AA_knot, B_knot, C_knot, D_knot))
#Thermal parameters for models

#Models where gamma is varied
        bf2zg=np.array([1.5704182,1.5754168,1.5797292,1.5836439,1.5870027,1.5915027,1.5943816,1.5986097,1.6027860,1.6073484,1.6116327,1.6158375])
        bf2zt=np.array([21221.521,21915.172,22303.908,22432.774,22889.881,22631.245,22610.916,22458.702,22020.437,21752.465,21278.358,20720.927])/1e3
        bf2g=np.array([1.4260277,1.4332204,1.4354097,1.4422902,1.4455976,1.4502961,1.4547729,1.4593188,1.4637958,1.4682989,1.4728728,1.4769461 ])
        bf2t=np.array([21453.313,21917.668,22655.974, 22600.650, 23061.612, 22798.62,22694.608,22721.762,22208.044, 21826.779, 21552.183,20908.159])/1e3
        bf2bg=np.array([1.1853212,1.1950833,1.2041183,1.2128127,1.2183342,1.2275667,1.2328254,1.2391863,1.2463307,1.2518918,1.2576843,1.263106])
        bf2bt=np.array([21464.978,22181.427,22594.685,22749.661,23218.541,22985.179,22974.422,22822.704,22387.685,22101.844,21612.527,21031.671])/1e3
        bf2cg=np.array([1.0216708,1.0334782,1.0445199,1.0552473,1.0616736,1.0729254,1.0792092,1.0865544,1.0950485,1.1010989,1.1077290,1.1138842])
        bf2ct=np.array([21561.621,22289.611,22714.780,22882.540,23358.805,23138.204,23134.987,22988.727,22560.401,22273.657,21786.965,21207.151])/1e3
        bf2dg=np.array([0.69947395,0.71491958,0.72959149,0.74404393,0.75198024,0.76702529,0.77516276,0.78435663,0.79540249,0.80256333,0.81069041,0.81829536])
        bf2dt=np.array([21769.984,22520.422,22969.096,23161.952,23655.146,23460.672,23475.417,23345.052,22934.061,22657.571,22182.125,21612.365])/1e3

#Models where T0 is varied; gamma changes a bit also
        T15g=np.array([1.3485881,1.3577712,1.3659720,1.3736938,1.3795726,1.3878026,1.3931488,1.3998762,1.4070048,1.4137031,1.4203951,1.4270162])
        T15t=np.array([13985.271,14489.076,14788.379,14914.136,15257.258,15123.863,15145.278,15078.820,14821.798,14677.305,14395.929,14059.430])/1e3
        T35g=np.array([1.3398174,1.3480600,1.3555022,1.3624323,1.3671665,1.3741617,1.3781567,1.3829295,1.3879474,1.3919720,1.3960214,1.3998402])
        T35t=np.array([32144.682,33160.696,33728.095,33910.458,34564.844,34159.548,34094.253,33809.077,33095.110,32603.121,31807.908,30881.237])/1e3
        T45g=np.array([1.3248931,1.3350548,1.3440339,1.3521382,1.3579792,1.3654530,1.3699611,1.3748911,1.3798712,1.3838409,1.3877698,1.3915019])
        T45t=np.array([40273.406,41586.734,42333.225,42588.319,43436.152,42930.669,42854.484,42486.193,41569.709,40929.328,39905.296,38717.452])/1e3

#Check knot
        bf2ag=np.array([1.0184905 ,1.0332849 ,1.0404092 ,1.0538235 ,1.0598905 ,1.0695964 ,1.0771414 ,1.0827464 ,1.0904784 ,1.0964669 ,1.1009428,1.1068357])
        bf2at=np.array([26914.856, 27497.555, 28407.030, 28357.785, 28919.376, 28610.038, 28478.749, 28484.856, 27841.160, 27335.328, 26933.453,26099.643])/1e3

        # G_knot=knot(("bf2z/","bf2b/","bf2c/","bf2d/", "bf2T15/","bf2T35/","bf2T45/","bf2a/"),(bf2zg,bf2bg,bf2cg,bf2dg,T15g,T35g,T45g,bf2ag),bf2g,"bf2/",60, (bf2zt,bf2bt,bf2ct,bf2dt,T15t,T35t,T45t,bf2at),bf2t)
        # g_int=flux_interp(flux,(G_knot))

