#!/usr/bin/env python
# vim: set fileencoding=UTF-8 :
#pylint:disable=R,C
"""
Various flux derivative stuff
"""
import math
import sys
import re
import matplotlib.pyplot as plt
import matplotlib.backends
import scipy.interpolate
import numpy as np

class DataError(Exception):
    """
    A class to check for data errors
    """
    def __init__(self, value):
        super().__init__()
        self.value = value
    def __str__(self):
        return repr(self.value)


def wheref(array, value):
    """
    A where function to find where a floating point value is equal to another
    """
    #Floating point inaccuracy.
    eps=1e-7
    return np.where((array > value-eps)*(array < value+eps))


def rebin(data, xaxis,newx):
    """Just rebins the data"""
    if newx[0] < xaxis[0] or newx[-1]> xaxis[-1]:
        raise ValueError("A value in newx is beyond the interpolation range")
    intp=scipy.interpolate.InterpolatedUnivariateSpline(np.log(xaxis),data)
    newdata=intp(np.log(newx))
    return newdata

def save_figure(path):
    """Saves the figure, automatically determining file extension"""
    bka=matplotlib.backends.backend
    if path == "":
        return()
    if bka in ('TkAgg', 'Agg', 'GTKAgg'):
        path = path+".png"
    if bka in ('PDF', 'pdf'):
        path = path+".pdf"
    if bka in ('PS', 'ps'):
        path = path+".ps"
    return plt.savefig(path)


def corr_table(table, dvecs,table_name):
    """Little function to adjust a table so it has a different central value"""
    new=np.array(table)
    new[12:,:] = table[12:,:]+2*table[0:12,:]*dvecs
    pkd="/home/spb41/cosmomc-src/cosmomc/data/lya-interp/"
    np.savetxt(pkd+table_name,new,("%1.3g","%1.3g","%1.3g","%1.3g","%1.3g"
    ,"%1.3g","%1.3g","%1.3g","%1.3g","%1.3g","%1.3g"))
    return new


class PowerSpec():
    """ A class to be derived from by flux and matter PowerSpec classes.
 Stores various helper methods."""
    # pylint: disable=too-many-instance-attributes
    # Snapshots
    snaps=()
    #SDSS redshift bins.
    zza=np.array([1.0])
    #SDSS kbins, in s/km units.
    sdsskbins=np.array([1.0])
    # Omega_matter
    omm=0.267
    #Hubble constant
    hco=0.71
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

    def __init__(self, snaps=("snapshot_000", "snapshot_001","snapshot_002","snapshot_003"
    ,"snapshot_004", "snapshot_005","snapshot_006","snapshot_007",
    "snapshot_008","snapshot_009","snapshot_010","snapshot_011"),
         zza=np.array([4.2,4.0,3.8,3.6,3.4,3.2,3.0,2.8,2.6,2.4,2.2,2.0]),
         sdsskbins=np.array([0.00141,0.00178,0.00224,0.00282,0.00355,
         0.00447,0.00562,0.00708,0.00891,0.01122,0.01413,0.01778]),
         knotpos=np.array([0.07,0.15,0.475, 0.75, 1.19, 1.89,4,25]), omm=0.267, hco=0.71,box=60.0,
         base="/home/spb41/Lyman-alpha/MinParametricRecon/runs/",suf="/", ext=".txt"):
        #pylint: disable=too-many-arguments
        if len(snaps) != np.size(zza):
            raise DataError("There are "+str(len(snaps))+" snapshots, but "+
            str(np.size(zza))+"redshifts given.")
        self.snaps=snaps
        self.zza=zza
        self.sdsskbins=sdsskbins
        self.omm=omm
        self.hco=hco
        self.box=box
        self.knotpos=knotpos*hco
        self.base=base
        self.suf=suf
        self.ext=ext

    def get(self, snap):
        """ Get redshift associated with a snapshot """
        ind=np.where(np.array(self.snaps) == snap)
        if np.size(ind):
            return self.zza[ind]
        raise DataError(str(snap)+" does not exist!")

    def get_snap(self, redshift):
        """ Get snapshot associated with a redshift """
        ind=wheref(self.zza, redshift)
        if np.size(ind):
            return str(np.asarray(self.snaps)[ind][0])
        raise DataError("No snapshot at redshift "+str(redshift))

    def get_sdss_kbins(self, redshift):
        """Get the k bins at a given redshift in h/Mpc units"""
        return self.sdsskbins*self.hubble(redshift)/(1.0+redshift)
#Corr is /sqrt(self.hco)...

    def hubble(self, zzb):
        """ Hubble parameter. Hubble(Redshift) """
        return 100*self.hco*math.sqrt(self.omm*(1+zzb)**3+(1-self.omm))
    #Conversion factor between s/km and h/Mpc is (1+z)/H(z)

    def loadpk(self, file, box):
        """ Do correct units conversion to return k and one-d power """
        #Adjust Fourier convention.
        flux_power=np.loadtxt(file)
        scale=self.hco/box
        k=flux_power[1:,0]*scale*2.0*math.pi
        pff=flux_power[1:,1]/scale
        return (k, pff)

    def plot_z(self,knot,redshift,title="",ylabel="", legend=True):
        """ Plot comparisons between a bunch of sims on one graph
    plot_z(Redshift, Sims to use ( eg, A1.14).
    Note this will clear current figures."""
        #pylint: disable=too-many-arguments
        #pylint: disable=too-many-locals
        #Load best-fit
        (simk,bfpk)=self.loadpk(knot.bstft+self.suf+self.pre+self.get_snap(redshift)+
        self.ext,self.bfbox)
        #Setup figure plot.
        ind=wheref(self.zza, redshift)
        plt.figure(ind[0][0])
        plt.clf()
        if title != '':
            plt.title(title+" at z="+str(redshift))
        plt.ylabel(ylabel)
        plt.xlabel(r"$k\; (\mathrm{Mpc}^{-1})$")
        line=np.array([])
        legname=np.array([])
        for sim in knot.names:
            (k,pka)=self.loadpk(sim+self.suf+self.pre+self.get_snap(redshift)+
            self.ext,self.box)
            oin = np.where(simk <= k[-1])
            tin = np.where(simk[oin] >= k[0])
            relp=rebin(pka, k, simk[oin][tin])
            relp=relp/rebin(bfpk, simk, simk[oin][tin])
            line=np.append(line, plt.semilogx(simk[oin][tin]/self.hco,relp,linestyle="-"))
            legname=np.append(legname,sim)
        if legend:
            plt.legend(line,legname)
            plt.semilogx(self.knotpos,np.ones(len(self.knotpos)),"ro")
        plt.ylim(self.ymin,self.ymax)
        plt.xlim(simk[0]*0.8, 10)

    def plot_all(self, knot,zzz=np.array([]), out=""):
        """ Plot a whole suite of snapshots: plot_all(knot, outdir) """
        if np.size(zzz) == 0:
            zzz=self.zza    #lolz
        for zaa in zzz:
            self.plot_z(knot,zaa)
            if out != "":
                save_figure(out+self.figprefix+str(zaa))

    def plot_power(self,path, redshift,colour="black"):
        """ Plot absolute power spectrum, not relative"""
        (k_g,pk_g)=self.loadpk(path+self.suf+self.pre+self.get_snap(redshift)+
        self.ext,self.box)
        plt.loglog(k_g,pk_g, color=colour)
        plt.xlim(0.01,k_g[-1]*1.1)
        plt.ylabel("P(k) /(h-3 Mpc3)")
        plt.xlabel("k /(h MPc-1)")
        plt.title("Power spectrum at z="+str(redshift))
        return(k_g, k_g)

    def plot_power_all(self, knot,zzz=np.array([]), out=""):
        """ Plot absolute power for all redshifts """
        if np.size(zzz) == 0:
            zzz=self.zza    #lolz
        for zaa in zzz:
            ind=wheref(self.zza, zaa)
            plt.figure(ind[0][0])
            for sim in knot.names:
                self.plot_power(sim,zaa)
            if out != "":
                save_figure(out+self.figprefix+str(zaa))

    def plot_compare_two(self, one, onebox, two,twobox,colour=""):
        """ Compare two power spectra directly. Smooths result.
       plot_compare_two(first P(k), second P(k))"""
        #pylint: disable=too-many-arguments
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
        """Get the difference between two simulations on scales probed by
       the SDSS power spectrum"""
        (onek,onep)=self.loadpk(one,self.bfbox)
        (twok,twop)=self.loadpk(two,self.box)
        onei = np.where(onek <= twok[-1])
        twoi= np.where (onek[onei] >= twok[0])
        relp=rebin(twop, twok, onek[onei][twoi])
        relp=relp/rebin(onep, onek, onek[onei][twoi])
        onek=onek[onei][twoi]
        sdss=self.get_sdss_kbins(redshift)
        relpr=np.ones(np.size(sdss))
        ind = np.where(sdss > onek[0])
        relpr[ind]=rebin(relp,onek,sdss[ind])
        return relpr

    def compare_two_table(self,onedir, twodir):
        """Do the above 12 times to get a correction table"""
        nka=np.size(self.sdsskbins)
        nza=np.size(self.zza)-1
        table=np.empty([nza,nka])
        for i in np.arange(0,nza):
            sim=self.pre+self.get_snap(self.zza[i])+self.ext
            table[-1-i,:]=self.compare_two(onedir+sim,twodir+sim,self.zza[i])
        return table

    def plot_compare_two_sdss(self, onedir, twodir, zzz=np.array([]), out="",
    title="", ylabel="", ymax=0,ymin=0, colour="",legend=False):
        """ Plot a whole redshift range of relative power spectra on the same figure.
    plot_all(onedir, twodir)Pass onedir and twodir as relative to basedir.
    ie, for default settings something like best-fit/flux-power/"""
        #pylint: disable=too-many-locals
        #pylint: disable=too-many-arguments
        if np.size(zzz) == 0:
            zzz=self.zza    #lolz
        line=np.array([])
        legname=np.array([])
        sdss=self.sdsskbins
        plt.xlabel(r"$k_v\; (s\,\mathrm{km}^{-1})$")
        for zaa in zzz:
            sim=self.pre+self.get_snap(zaa)+self.ext
            pka=self.compare_two(onedir+sim,twodir+sim,zaa)
            if colour == "":
                line=np.append(line, plt.semilogx(sdss,pka))
            else:
                line=np.append(line, plt.semilogx(sdss,pka,color=colour))
            legname=np.append(legname,"z="+str(zaa))
        plt.title(title)
        if ylabel != "":
            plt.ylabel(ylabel)
        if legend:
            plt.legend(line, legname,bbox_to_anchor=(0., 0, 1., .25), loc=3,
            ncol=3,mode="expand", borderaxespad=0.)
        if ymax != 0 and ymin !=0:
            plt.ylim(ymin,ymax)
        plt.xlim(sdss[0],sdss[-1])
        plt.xticks(np.array([sdss[0],3e-3,5e-3,0.01,sdss[-1]]),("0.0014","0.003",
        "0.005","0.01","0.0178"))
        if out != "":
            save_figure(out)
        return plt.gcf()

    def plot_compare_two_all(self, onedir, twodir, zzz=np.array([]), out="", title="",
    ylabel="", ymax=0,ymin=0, colour="",legend=False):
        """ Plot a whole redshift range of relative power spectra on the same figure.
            plot_all(onedir, twodir)
            Pass onedir and twodir as relative to basedir.
            ie, for default settings something like
            best-fit/flux-power/
            onedir uses bfbox, twodir uses box"""
        #pylint: disable=too-many-arguments
        if np.size(zzz) == 0:
            zzz=self.zza    #lolz
        line=np.array([])
        legname=np.array([])
        for zaa in zzz:
            line=np.append(line, self.plot_compare_two(onedir+self.pre+self.get_snap(zaa)+
            self.ext,self.bfbox,twodir+self.pre+self.get_snap(zaa)+self.ext,self.box,colour))
            legname=np.append(legname,"z="+str(zaa))
        if title == "":
            plt.title("Relative Power spectra "+onedir+" and "+twodir)
        else:
            plt.title(title)
        if ylabel != "":
            plt.ylabel(ylabel)
        if legend:
            plt.legend(line, legname,bbox_to_anchor=(0., 0, 1., .25), loc=3,
            ncol=3,mode="expand", borderaxespad=0.)
        if ymax != 0 and ymin !=0:
            plt.ylim(ymin,ymax)
        if out != "":
            save_figure(out)
        return plt.gcf()

    def get_flat(self,dire, sil=0.0):
        """Get a power spectrum in the flat format we use
        for inputting some cosmomc tables"""
        pk_sdss=np.empty([11, 12])
        #Note this omits the z=2.0 bin
        #SiIII corr now done on the fly in lya_sdss_viel.f90
        for i in np.arange(0,np.size(self.snaps)-1):
            scale=self.hubble(self.zza[i])/(1.0+self.zza[i])
            (kaa,pka)=self.loadpk(dire+self.snaps[i]+self.ext, self.box)
            fbar=math.exp(-0.0023*(1+self.zza[i])**3.65)
            avg=sil/(1-fbar)
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
            sicorr[0]=1+avg**2+2*avg*math.cos(2271*sdss[0])
            for j in np.arange(1,np.size(sdss)):
                sicorr[j]=1+avg**2+2*avg*(math.sin(2271*kmids[j+1])-
                math.sin(2271*kmids[j]))/(kmids[j+1]-kmids[j])/2271
            sdss=self.get_sdss_kbins(self.zza[i])
            pk_sdss[-i-1,:]=rebin(pka, kaa, sdss)*scale*sicorr
        return pk_sdss

    def calc_z(self, redshift,s_knot, kbins):
        """ Calculate the flux derivatives for a single redshift
        Output: (kbins d2P...kbins dP (flat vector of length 2xkbins))"""
        #Array to store answers.
        #Format is: k x (dP, d²P, χ²)
        #pylint: disable=too-many-locals
        kbins=np.array(kbins)
        nka=np.size(kbins)
        snap=self.get_snap(redshift)
        results=np.zeros(2*nka)
        if np.size(s_knot.qvals) > 0:
            results = np.zeros(4*nka)
        pdifs=s_knot.pvals-s_knot.p0
        qdifs=np.array([])
        #If we have somethign where the parameter is redshift-dependent, eg, gamma.
        if np.size(s_knot.p0) > 1:
            i=wheref(redshift, self.zza)
            pdifs = s_knot.pvals[:,i]-s_knot.p0[i]
            if np.size(s_knot.qvals) > 0:
                qdifs=s_knot.qvals[:,i] - s_knot.q0[i]
        npvals=np.size(pdifs)
        #Load the data
        (k,pfpo)=self.loadpk(s_knot.bstft+self.suf+snap+self.ext,s_knot.bfbox)
        #This is to rescale by the mean flux, for generating mean flux tables.
        ###
        #tau_eff=0.0023*(1+redshift)**3.65
        #tmin=0.2*((1+redshift)/4.)**2
        #tmax=0.5*((1+redshift)/4.)**4
        #teffs=tmin+s_knot.pvals*(tmax-tmin)/30.
        #pdifs=teffs/tau_eff-1.
        ###
        power_fluxes=np.zeros((npvals,np.size(k)))
        for i in np.arange(0,np.size(s_knot.names)):
            (k,power_fluxes[i,:])=self.loadpk(s_knot.names[i]+self.suf+snap+self.ext, s_knot.bfbox)
        #So now we have an array of data values, which we want to rebin.
        ind = np.where(kbins >= k[0])
        difpf_rebin=np.ones((npvals,np.size(kbins)))
        for i in np.arange(0, npvals):
            difpf_rebin[i,ind]=rebin(power_fluxes[i,:]/pfpo,k,kbins[ind])
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
            results[nka+k]=derivs[1]
            if np.size(derivs) > 2:
                results[2*nka+k]=derivs[2]
                results[3*nka+k]=derivs[3]
        return results

    def calc_all(self, s_knot,kbins):
        """ Calculate the flux derivatives for all redshifts
    Input: Sims to load, parameter values, mean parameter value
    Output: (2*kbins) x (zbins)"""
        flux_derivatives=np.zeros((2*np.size(kbins),np.size(self.zza)))
        if np.size(s_knot.qvals) > 1:
            flux_derivatives=np.zeros((4*np.size(kbins),np.size(self.zza)))
        #Call flux_deriv_const_z for each redshift.
        for i in np.arange(0,np.size(self.zza)):
            flux_derivatives[:,i]=self.calc_z(self.zza[i], s_knot,kbins)
        return flux_derivatives

    def flux_deriv(self, pfdif, pdif, qdif=np.array([])):
        """Calculate the flux-derivative for a single redshift and k bin"""
        pdif=np.ravel(pdif)
        if np.size(pdif) != np.size(pfdif):
            raise DataError(str(np.size(pdif))+" parameter values, but "+
            str(np.size(pfdif))+" P_F values")
        if np.size(pdif) < 2:
            raise DataError(str(np.size(pdif))+" pvals given. Need at least 2.")
        pfdif=pfdif-1.0
        if np.size(qdif) > 2:
            qdif=np.ravel(qdif)
            mat=np.vstack([pdif**2, pdif, qdif**2, qdif] ).T
        else:
            mat=np.vstack([pdif**2, pdif] ).T
        derivs=np.linalg.lstsq(mat, pfdif)
        return derivs

    def get_error_z(self, sim, bstft,box, derivs, params, redshift,
    qarams=np.empty([])):
        """ Get the error on one test simulation at single redshift """
        #Need to load and rebin the sim.
        #pylint: disable=too-many-arguments
        (k, test)=self.loadpk(sim+self.suf+self.get_snap(redshift)+self.ext, box)
        (k,bfa)=self.loadpk(bstft+self.suf+self.get_snap(redshift)+self.ext,box)
        kbins=np.array(redshift)
        ind = np.where(kbins >= 1.0*2.0*math.pi*self.hco/box)
        test2=np.ones(np.size(kbins))
        test2[ind]=rebin(test/bfa,k,kbins[ind])#
        if ind[0][0] > 0:
            test2[0:ind[0][0]]=test2[ind[0][0]]
        if np.size(qarams) > 0:
            guess=derivs.GetPF(params, redshift,qarams)+1.0
        else:
            guess=derivs.GetPF(params, redshift)+1.0
        return np.array((test2/guess))[0][0]

class FluxPow(PowerSpec):
    """ A class written to store the various methods related to calculating
of the flux derivatives and plotting of the flux power spectra"""
    figprefix="/flux-figure"
    kbins=np.array([])
    def __init__(self, snaps=("snapshot_000", "snapshot_001","snapshot_002","snapshot_003",
    "snapshot_004","snapshot_005","snapshot_006","snapshot_007","snapshot_008",
    "snapshot_009","snapshot_010","snapshot_011"),
         zza=np.array([4.2,4.0,3.8,3.6,3.4,3.2,3.0,2.8,2.6,2.4,2.2,2.0]),
         sdsskbins=np.array([0.00141,0.00178,0.00224,0.00282,0.00355,0.00447,
0.00562,0.00708,0.00891,0.01122,0.01413,0.01778]),
         knotpos=np.array([0.07,0.15,0.475, 0.75, 1.19, 1.89,4,25]), omm=0.266,
         hco=0.71,box=60.0,kmax=4.0,
         base="/home/spb41/Lyman-alpha/MinParametricRecon/runs/",bf="best-fit/",
         suf="flux-power/", ext="_flux_power.txt"):
        #pylint: disable=too-many-arguments
        PowerSpec.__init__(self, snaps,zza,sdsskbins,knotpos, omm, hco,
        box,base,suf, ext)
        k_bf= self.loadpk(bf+suf+"snapshot_000"+self.ext,self.bfbox)
        ind=np.where(k_bf <= kmax)
        self.kbins=ind

    def plot_za(self,sims,redshift,title="Relative Flux Power",
    ylabel=r"$\mathrm{P}_\mathrm{F}(k,p)\,/\,\mathrm{P}_\mathrm{F}(k,p_0)$", legend=True):
        """Plotting z"""
        #pylint: disable=too-many-arguments
        PowerSpec.plot_z(self,sims,redshift,title,ylabel,legend)
        if legend:
            kbins=self.get_sdss_kbins(redshift)
            plt.axvspan(kbins[0], kbins[-1], color="#B0B0B0")
        plt.ylim(self.ymin,self.ymax)
        plt.xlim(self.kbins[0]*0.8, 10)

    def get_kbins(self):
        """Get the kbins to interpolate onto"""
        return self.kbins

    def mac_donald_pf(self,sdss, zzb):
        """Load the SDSS power spectrum"""
        psdss=sdss[np.where(sdss[:,0] == zzb)][:,1:3]
        fbar=math.exp(-0.0023*(1+zzb)**3.65)
        #multiply by the hubble parameter to be in 1/(km/s)
        scale=self.hubble(zzb)/(1.0+zzb)
        pfb=psdss[:,1]*fbar**2/scale
        k=psdss[:,0]*scale
        return (k, pfb)

    def loadpka(self, path,box):
        """Load a Pk. Different function due to needing to be different for each class"""
        #Adjust Fourier convention.
        flux_power=np.loadtxt(self.base+path)
        scale=self.hco/box
        k=(flux_power[1:,0]-0.5)*scale*2.0*math.pi
        pfb=flux_power[1:,1]/scale
        return (k, pfb)

class MatterPow(PowerSpec):
    """ A class to plot matter power spectra """
    obs=0.0
    #For plotting
    ymin=0.4
    ymax=1.6
    figprefix="/matter-figure"
    def __init__(self, snaps=("snapshot_000", "snapshot_001", "snapshot_002",
    "snapshot_003","snapshot_004","snapshot_005","snapshot_006","snapshot_007",
    "snapshot_008","snapshot_009","snapshot_010","snapshot_011"),
         zza=np.array([4.2,4.0,3.8,3.6,3.4,3.2,3.0,2.8,2.6,2.4,2.2,2.0]),
         sdsskbins=np.array([0.00141,0.00178,0.00224,0.00282,0.00355,
         0.00447,0.00562,0.00708,0.00891,0.01122,0.01413,0.01778]),
         knotpos=np.array([0.07,0.15,0.475, 0.75, 1.19, 1.89,4,25]), omm=0.266,
obs=0.0449, hco=0.71,box=60.0,
         base="/home/spb41/Lyman-alpha/MinParametricRecon/runs/",suf="matter-power/",
         ext=".0", matpre="PK-by-"):
        #pylint: disable=too-many-arguments
        PowerSpec.__init__(self, snaps,zza,sdsskbins,knotpos, omm, hco,box,base,suf,ext)
        self.obs=obs
        self.pre=matpre

    def plot_zb(self,sims,redshift,title="Relative Matter Power",
ylabel=r"$\mathrm{P}(k,p)\,/\,\mathrm{P}(k,p_0)$"):
        """Plotting zb"""
        PowerSpec.plot_z(self,sims,redshift,title,ylabel)

    def loadpkb(self, path,box):
        """Load a Pk. Different function due to needing to be different for each class"""
        #Load baryon P(k)
        matter_power=np.loadtxt(self.base+path)
        scale=self.hco/box
        #Adjust Fourier convention.
        simk=matter_power[1:,0]*scale*2.0*math.pi
        pk_bar=matter_power[1:,1]/scale**3
        #Load DM P(k)
        matter_power=np.loadtxt(self.base+re.sub("by","DM",path))
        pkdm=matter_power[1:,1]/scale**3
        pkc=(pk_bar*self.obs+pkdm*(self.omm-self.obs))/self.omm
        return (simk,pkc)

    def plot_powera(self,path, redshift,camb_filename=""):
        """ Plot absolute power spectrum, not relative"""
        (k_g, pk_g)=PowerSpec.plot_power(self,path,redshift)
        sigma=2.0
        pkg=np.loadtxt(self.base+path+self.suf+self.pre+self.get_snap(redshift)+self.ext)
        samp_err=pkg[1:,2]
        sqrt_err=np.array(np.sqrt(samp_err))
        plt.loglog(k_g,pk_g*(1+sigma*(2.0/sqrt_err+1.0/samp_err)),
        linestyle="-.",color="black")
        plt.loglog(k_g,pk_g*(1-sigma*(2.0/sqrt_err+1.0/samp_err)),
        linestyle="-.",color="black")
        if camb_filename != "":
            camb=np.loadtxt(camb_filename)
            #Adjust Fourier convention.
            pkc=camb[:,0]*self.hco
            #NOW THERE IS NO h in the T anywhere.
            pkc=camb[:,1]
            plt.loglog(k_g/self.hco, pkc, linestyle="--")
        plt.xlim(0.01,k_g[-1]*1.1)
        return(k_g, pk_g)

class FluxPdf(PowerSpec):
    """The PDF is an instance of the PowerSpec class. Perhaps poor naming"""
    def __init__(self, snaps=("snapshot_006","snapshot_007",
    "snapshot_008","snapshot_009","snapshot_010","snapshot_011"),
         zza=np.array([3.0,2.8,2.6,2.4,2.2,2.0]),
         sdsskbins=np.arange(0,20),
         knotpos=np.array([]), omm=0.266, hco=0.71,box=48.0,
         base="/home/spb41/Lyman-alpha/MinParametricRecon/runs/",
         suf="flux-pdf/", ext="_FluxPdf.txt",):
        #pylint: disable=too-many-arguments
        PowerSpec.__init__(self, snaps,zza,sdsskbins,knotpos,omm,hco,box,base,suf,ext)

    def loadpkc(self, path):
        """Returning Fluxpdf"""
        flux_pdf = np.loadtxt(self.base+path)
        return(flux_pdf[:,0], flux_pdf[:,1])

    def plot_compare_twoa(self, one, onebox, two,twobox):
        """ Compare two power spectra directly. Smooths result.
        plot_compare_two(first P(k), second P(k))"""
        (onek,onep)=self.loadpk(one,onebox)
        (twok,twop)=self.loadpk(two,twobox)
        relp=onep/twop
        dummy=twok
        plt.title("Relative flux PDF "+one+" and "+two)
        plt.ylabel(r"$F_2(k)/F_1(k)$")
        plt.xlabel(r"$Flux$")
        line=plt.semilogy(onek,relp)
        return line

    def plot_powerb(self,path, redshift):
        """ Plot absolute power spectrum, not relative"""
        (k,pdf)=self.loadpk(path+self.suf+self.pre+self.get_snap(redshift)+
        self.ext,self.box)
        plt.semilogy(k,pdf, color="black", linewidth="1.5")
        plt.ylabel("P(k) /(h-3 Mpc3)")
        plt.xlabel("k /(h MPc-1)")
        plt.title("PDF at z="+str(redshift))
        return(k, pdf)

    def calc_za(self, redshift,s_knot):
        """ Calculate the flux derivatives for a single redshift
    Output: (kbins d2P...kbins dP (flat vector of length 2x21))"""
        #Array to store answers.
        #Format is: k x (dP, d²P, χ²)
        #pylint: disable=too-many-locals
        npvals=np.size(s_knot.pvals)
        nka=21
        results=np.zeros(2*nka)
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
        usnap=self.get_snap(ured)
        lsnap=self.get_snap(lred)
        #Load the data
        (k,upfpo)=self.loadpk(s_knot.bstft+self.suf+usnap+self.ext,s_knot.bfbox)
        u_power=np.zeros((npvals,np.size(k)))
        for i in np.arange(0,np.size(s_knot.names)):
            (k,u_power[i,:])=self.loadpk(s_knot.names[i]+self.suf+usnap+self.ext,s_knot.bfbox)
        (k,lpfpo)=self.loadpk(s_knot.bstft+self.suf+lsnap+self.ext,s_knot.bfbox)
        l_power=np.zeros((npvals,np.size(k)))
        for i in np.arange(0,np.size(s_knot.names)):
            (k,l_power[i,:])=self.loadpk(s_knot.names[i]+self.suf+lsnap+self.ext,s_knot.bfbox)
        power_fluxes=5*((redshift-lred)*u_power+(ured-redshift)*l_power)
        pfpo=5*((redshift-lred)*upfpo+(ured-redshift)*lpfpo)
        qdifs=np.array([])
        #So now we have an array of data values.
        #Pass each k value to flux_deriv in turn.
        for k in np.arange(0,nka):
            deriv_pf=self.flux_deriv(power_fluxes[:,k]/pfpo[k], pdifs,qdifs)
            results[k]=deriv_pf[0]
            results[nka+k]=deriv_pf[1]
        return results

    def get_kbins(self):
        """Get the kbins to interpolate onto"""
        return np.arange(0,20,1)+0.5

    def plot_zc(self,knot,redshift,title="",ylabel=""):
        """ Plot comparisons between a bunch of sims on one graph
    plot_z(Redshift, Sims to use ( eg, A1.14).
    Note this will clear current figures."""
        #Load best-fit
        (simk,bfpk)=self.loadpk(knot.bstft+self.suf+self.pre+self.get_snap(redshift)+
        self.ext,self.bfbox)
        #Setup figure plot.
        ind=wheref(self.zza, redshift)
        plt.figure(ind[0][0])
        plt.clf()
        if title != '':
            plt.title(title+" at z="+str(redshift),)
        plt.ylabel(ylabel)
        plt.xlabel(r"$\mathcal{F}$")
        line=np.array([])
        legname=np.array([])
        for sim in knot.names:
            pka=self.loadpk(sim+self.suf+self.pre+self.get_snap(redshift)+self.ext,self.box)
            line=np.append(line, plt.semilogy(simk,pka/bfpk,linestyle="-", linewidth=1.5))
            legname=np.append(legname,sim)
        plt.legend(line,legname)

    def get_flata(self,dire):
        """Get a power spectrum in the flat format we use
        for inputting some cosmomc tables"""
        #pk_sdss=np.empty([11, 12])
        #For z=2.07 we need to average snap_011 and snap_010
        zzd=2.07
        pf_a=self.loadpk(dire+self.suf+"snapshot_011"+self.ext,self.box)
        pf_b=self.loadpk(dire+self.suf+"snapshot_010"+self.ext,self.box)
        pf_one=(zzd-2.0)*5*(pf_b-pf_a)+pf_a
        zzd=2.52
        pf_a=self.loadpk(dire+self.suf+"snapshot_009"+self.ext,self.box)
        pf_b=self.loadpk(dire+self.suf+"snapshot_008"+self.ext,self.box)
        pf_two=(zzd-2.4)*5*(pf_b-pf_a)+pf_a
        zzd=2.94
        pf_a=self.loadpk(dire+self.suf+"snapshot_007"+self.ext,self.box)
        pf_b=self.loadpk(dire+self.suf+"snapshot_006"+self.ext,self.box)
        pf_three=(zzd-2.8)*5*(pf_b-pf_a)+pf_a
        pdf = np.array([pf_one,pf_two,pf_three])
        np.savetxt(sys.stdout,pdf,("%1.8f","%1.8f","%1.8f"))
        return (pf_one, pf_two, pf_three)

if __name__=='__main__':
    flux=FluxPow()
    matter=MatterPow()
    fpdf=FluxPdf()
    A_knot=(("A0.54/","A0.74/","A0.84/","A1.04/", "A1.14/","A1.34/"),
    (0.54,0.74,0.84,1.04,1.14,1.34),0.94,"best-fit/", 60)
    AA_knot=(("AA0.54/","AA0.74/","AA1.14/","AA1.34/"), (0.54,0.74,1.14,1.34),
    0.94,"boxcorr400/", 120)
    B_knot=(("B0.33/","B0.53/","B0.73/", "B0.83/", "B1.03/","B1.13/", "B1.33/"),
    (0.33,0.53,0.73,0.83,1.03, 1.13,1.33),0.93,"best-fit/", 60)
    C_knot=(("C0.11/", "C0.31/","C0.51/","C0.71/","C1.11/","C1.31/","C1.51/"),
    (0.11, 0.31,0.51,0.71, 1.11,1.31,1.51),0.91,"bf2/", 60)
    D_knot=(("D0.50/","D0.70/","D1.10/","D1.20/", "D1.30/", "D1.50/", "D1.70/"),
    (0.50, 0.70,1.10,1.20,1.30, 1.50, 1.70),0.90,"bfD/", 48)
    #interp=FluxPow.__init__(flux, (AA_knot, B_knot, C_knot, D_knot))
#Thermal parameters for models

#Models where gamma is varied
    bf2zg=np.array([1.5704182,1.5754168,1.5797292,1.5836439,1.5870027,1.5915027,1.5943816,
    1.5986097,1.6027860,1.6073484,1.6116327,1.6158375])
    bf2zt=np.array([21221.521,21915.172,22303.908,22432.774,22889.881,22631.245,
    22610.916,22458.702,22020.437,21752.465,21278.358,20720.927])/1e3
    bf2g=np.array([1.4260277,1.4332204,1.4354097,1.4422902,1.4455976,1.4502961,
    1.4547729,1.4593188,1.4637958,1.4682989,1.4728728,1.4769461 ])
    bf2t=np.array([21453.313,21917.668,22655.974, 22600.650, 23061.612, 22798.62,
    22694.608,22721.762,22208.044, 21826.779, 21552.183,20908.159])/1e3
    bf2bg=np.array([1.1853212,1.1950833,1.2041183,1.2128127,1.2183342,1.2275667,
    1.2328254,1.2391863,1.2463307,1.2518918,1.2576843,1.263106])
    bf2bt=np.array([21464.978,22181.427,22594.685,22749.661,23218.541,22985.179,
    22974.422,22822.704,22387.685,22101.844,21612.527,21031.671])/1e3
    bf2cg=np.array([1.0216708,1.0334782,1.0445199,1.0552473,1.0616736,1.0729254,
    1.0792092,1.0865544,1.0950485,1.1010989,1.1077290,1.1138842])
    bf2ct=np.array([21561.621,22289.611,22714.780,22882.540,23358.805,23138.204,
    23134.987,22988.727,22560.401,22273.657,21786.965,21207.151])/1e3
    bf2dg=np.array([0.69947395,0.71491958,0.72959149,0.74404393,0.75198024,
    0.76702529,0.77516276,0.78435663,0.79540249,0.80256333,0.81069041,0.81829536])
    bf2dt=np.array([21769.984,22520.422,22969.096,23161.952,23655.146,23460.672,
    23475.417,23345.052,22934.061,22657.571,22182.125,21612.365])/1e3

#Models where T0 is varied; gamma changes a bit also
    T15g=np.array([1.3485881,1.3577712,1.3659720,1.3736938,1.3795726,1.3878026,
    1.3931488,1.3998762,1.4070048,1.4137031,1.4203951,1.4270162])
    T15t=np.array([13985.271,14489.076,14788.379,14914.136,15257.258,15123.863,
    15145.278,15078.820,14821.798,14677.305,14395.929,14059.430])/1e3
    T35g=np.array([1.3398174,1.3480600,1.3555022,1.3624323,1.3671665,1.3741617,
    1.3781567,1.3829295,1.3879474,1.3919720,1.3960214,1.3998402])
    T35t=np.array([32144.682,33160.696,33728.095,33910.458,34564.844,34159.548,
    34094.253,33809.077,33095.110,32603.121,31807.908,30881.237])/1e3
    T45g=np.array([1.3248931,1.3350548,1.3440339,1.3521382,1.3579792,1.3654530,
    1.3699611,1.3748911,1.3798712,1.3838409,1.3877698,1.3915019])
    T45t=np.array([40273.406,41586.734,42333.225,42588.319,43436.152,42930.669,
    42854.484,42486.193,41569.709,40929.328,39905.296,38717.452])/1e3

#Check knot
    bf2ag=np.array([1.0184905 ,1.0332849 ,1.0404092 ,1.0538235 ,1.0598905
    ,1.0695964 ,1.0771414 ,1.0827464 ,1.0904784 ,1.0964669 ,1.1009428,1.1068357])
    bf2at=np.array([26914.856, 27497.555, 28407.030, 28357.785, 28919.376,
    28610.038, 28478.749, 28484.856, 27841.160, 27335.328, 26933.453,26099.643])/1e3

    G_knot=(("bf2z/","bf2b/","bf2c/","bf2d/", "bf2T15/","bf2T35/",
    "bf2T45/","bf2a/"),(bf2zg,bf2bg,bf2cg,bf2dg,T15g,T35g,T45g,bf2ag),bf2g,"bf2/",60, (bf2zt,bf2bt,
    bf2ct,bf2dt,T15t,T35t,T45t,bf2at),bf2t)
   # g_int=FluxPow__init__(flux,(G_knot))
