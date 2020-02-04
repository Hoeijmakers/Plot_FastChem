



#This is a script to plot the FastChem output into a spaghettiplot.
class FastChem_output(object):
    def __init__(self,p,plot_species=['Fe','C1O1','H2O1','H'],outpath='',noshow=False):
        """
        Initialise the FastChem output class used to plot the output of FastChem.
        Parameters
        ----------
        p : file path
            File to read the output from.
        plot_species: list of strs
            A list containing the FastChem labels of the species to be plotted.
            If a label is set wrongly, the code will abort, printing the available
            labels as recognised from the map-file.
        outpath: file path
            Optional outpath of the figure that is going to be plotted. If not set,
            or set to '', no figure will be written, and the plot is shown on screen instead (default).
        Returns
        -------
        FastChem_output : `FastChem_output` object
            Upon init, this object will plot the default species, either to screen if outpath is
            not set (default), or to file if outpath is set.

        Attributes
        -------
        self.plot_species: list of str
            Same as plot_species above. Used to change the list of species to plot
            after init.

        self.ymin: int,float
            Maximum x-axis range.
        self.ymax: int,float
            Maximum x-axis range.
        self.xmin: int,float
            Maximum x-axis range.
        self.xmax: int,float
            Maximum x-axis range.


        self.outpath: file path
            See above.
        self.fontsize: int
            Font size of the axes and labels.
        self.dpi: int
            DPI of the output figure.
        self.figsize: tuple
            Figsize, as in the plt.subplots figsize keyword.


        self.colors: list of str
            Matplotlib colour names in a list of arbitrary length. The species lines to plot
            will cycle through this list rather than the default. Can be used to specify the colour
            of each single element.
        self.styles: list of str
            Matplotlib short-hand line styles (':','-.','--','-') in a list of arbitrary length.
            The species lines to plot will cycle through this list. Set to '-' by default, meaning
            that all lines will be solid.
        self.labels: list of str
            Manually set the species labels. You probably need this for making paper-ready figures. If you set a label to '', the line will be ignored in the legend.
            If you set fewer labels than species, the trailing species will not make it into the legend.
        self.alpha: list of float
            Transparency values of lines. Same philosophy as colours and styles.
        self.linewidth: list of float
            Line thicknesses.
        self.Tcolor: str
            Matplotlib colour of the temperature profile on the second axis.
        self.Tstyle: str
            Matplotlib short-hand line style for the TP profile.
        self.plot_TP: Bool
            Turn plotting of the TP profile on or off. On by default.


        Methods
        -------
        self.print_species()
            Print the FastChem labels of all available species that can be included
            in plot_species.
        self.plot(show==False)
            Re-plot the Figure (after changing the plot parameters for instance).
            set show==True to force plotting the Fig. on-screen, even if outpath is set.


        Example
        -------
        You can make use of this class in e.g. the following way.
        my_environment$ python3
        >>>from plotting import FastChem_output as FCO
        >>>a = FCO('output/chem_output_Kepler-7b.dat')
        >>>a.plot_Tcolor = 'green'
        >>>a.Tstyle = '-'
        >>>a.styles = ['--']
        >>>a.outpath = 'Kepler_7_chemistry.png'
        >>>a.plot()
        """



        import numpy as np
        import astropy.io.ascii as ascii
        from astropy.table import Table
        import pdb
        f = open(p,'r').read().splitlines()
        header = f[0].split('\t')


        self.species = [i.replace(' ','') for i in header][5:]
        # self.species = ascii.read(map_path,delimiter=' = ',data_start=0,names=['Name','Col_Index'])#,headerstart=None)

        self.data = ascii.read(p,delimiter='\t',data_start=1,header_start=0,names=['P (bar)','T (K)','n_H (cm-3)','n_g (cm-3)','m (u)']+self.species)
        self.plot_species = plot_species

        #The following encodes defaults for setting up the plotting window. These
        #can be changed by the user.
        P = self.data['P (bar)']
        self.ymin = max(P)
        self.ymax = min(P)
        self.xmin = 0
        self.xmax = 0
        self.xlabel = 'Volume mixing ratio'
        self.ylabel = 'P (bar)'
        self.fontsize = 9
        self.outpath=outpath#The out path of the figure.
        self.dpi = 400
        self.figsize = (6,4.5)
        #All linestyle stuff:
        self.colors = []#If this has nonzero length, lines will cycle through this list.
        self.styles = ['-']
        self.labels = []
        self.alpha = [1.0]
        self.linewidth = [1.0]

        self.Tcolor = 'red'
        self.Tstyle = '--'
        self.plot_TP = True


        if noshow == False:
            self.plot()


    def print_species(self):
        for i in self.species:
            print(i)


    def resolve(self,label):
        """ This resolves a species label generated by FastChem."""
        import copy
        import math

        Natoms = sum(int(x) for x in label if x.isdigit())#Sum the numbers in the string. If greater than 1, we have a molecule.
        Nion = int(label.count('+'))
        Nanion = int(label.count('-'))
        if Natoms >= 2:
            molecule = True
        else:
            molecule = False
        outlabel = copy.deepcopy(label).replace('1','')

        if molecule == True:
            for i in range (99,0,-1):
                outlabel=outlabel.replace(str(i),'$_{'+str(i)+'}$')#Replace each character that remains with  LaTeX.
            #We go from 99 downwards in order to take into account hypothetical molecules with double-digit numbers.
            outlabel = outlabel.replace('+','$^{+}$')
            outlabel = outlabel.replace('-','$^{-}$')
        else:
            if Nanion >= 1: #We are an anion.
                outlabel=outlabel.replace('-','$^-$')
            else:
                outlabel=outlabel.replace('+','')
                outlabel+=r'{\fontsize{%spt}{3em}\fontfamily{cursif}\selectfont{}{ %s}'%(int(math.floor(self.fontsize*0.85)),(Nion+1)*'I')
        return(outlabel)


    def plot(self,show=False):
        import matplotlib.pyplot as plt
        from matplotlib import rc
        import sys

        for s in self.plot_species:
            if s not in self.species:
                self.print_species()
                print('')
                print("Error: %s is not available. Above are the available species."%s)
                return

        plt.rcParams.update({'font.size': self.fontsize})
        rc('text', usetex=True)

        fig, ax1 = plt.subplots(figsize=self.figsize)
        P = self.data['P (bar)']
        lines=[]
        for i,name in enumerate(self.plot_species):

            #Deal with customly set labels.
            if len(self.labels) == 0:
                    label = self.resolve(name)
            else:
                label = 'empty'
                if i < len(self.labels):#If not, the label is left empty at ''.
                    label = self.labels[i]
                    if label == '':
                        label = 'empty'#This is a placeholder word really.

            #Deal with line colours, styles, thicknesses and transparencies:
            if len(self.colors) == 0:
                l,=ax1.plot(self.data[name],P,self.styles[i%(len(self.styles))],label=label,linewidth=self.linewidth[i%(len(self.linewidth))],alpha=self.alpha[i%(len(self.alpha))])
            else:
                l,=ax1.plot(self.data[name],P,self.styles[i%(len(self.styles))],label=label,color=self.colors[i%(len(self.colors))],linewidth=self.linewidth[i%(len(self.linewidth))],alpha=self.alpha[i%(len(self.alpha))])
            lines.append(l)

        ax1.set_yscale('log')
        ax1.set_xscale('log')

        if self.xmin > 0:
            xmin = self.xmin
        else:
            xmin = ax1.get_xlim()[0]
        if self.xmax > self.xmin:
            xmax = self.xmax
        else:
            xmax = ax1.get_xlim()[1]

        ax1.set_xlabel(self.xlabel)
        ax1.set_ylabel(self.ylabel)
        ax1.set_ylim(self.ymin,self.ymax)
        ax1.set_xlim(xmin,xmax)
        plt.gca().invert_yaxis()

        if self.plot_TP == True:
            ax2=ax1.twiny()
            ax2.set_xlabel('T (K)')  # we already handled the x-label with ax1
            l,=ax2.plot(self.data['T (K)'],P,self.Tstyle,color=self.Tcolor,label='T (K)')
            lines.append(l)
            plt.gca().invert_yaxis()
        fig.tight_layout()  # otherwise the right y-label is slightly clipped

        #This is for dealing with the empty labels that may have remained:
        labels_accepted = []
        lines_accepted = []
        for l in lines:
            la = l.get_label()
            if la != 'empty':
                labels_accepted.append(la)
                lines_accepted.append(l)
        leg = plt.legend(lines_accepted, labels_accepted, loc=0)
        if show == True:#If show == True, always show.
            self.interactive_legend(fig,ax1,lines,leg)
            plt.show()
        elif len(self.outpath) > 0:#Otherwise, check that the outpath is set, if so, we save.
            plt.savefig(self.outpath,dpi=self.dpi)
        else:#If show == False but the outpath is not set (default), we still show:
            self.interactive_legend(fig,ax1,lines,leg)
            plt.show()




    def interactive_legend(self,fig,ax,lines,leg):
        #leg = ax.legend(loc='upper left', fancybox=False, shadow=False)
        leg.get_frame().set_alpha(0.4)
# we will set up a dict mapping legend line to orig line, and enable
# picking on the legend line
        lined = dict()
        for legline, origline in zip(leg.get_lines(), lines):
            legline.set_picker(5)  # 5 pts tolerance
            lined[legline] = origline
        def onpick(event):
    # on the pick event, find the orig line corresponding to the
    # legend proxy line, and toggle the visibility
            legline = event.artist
            origline = lined[legline]
            vis = not origline.get_visible()
            origalpha = origline.get_alpha()
            origline.set_visible(vis)
    # Change the alpha on the line in the legend so we can see what lines
    # have been toggled
            if vis:
                legline.set_alpha(origalpha)
            else:
                legline.set_alpha(origalpha/3.0)
            fig.canvas.draw()
        fig.canvas.mpl_connect('pick_event', onpick)

# a=FastChem_output('output/chem_output.dat')
