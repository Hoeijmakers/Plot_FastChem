#This is a script to plot the FastChem output into a spaghettiplot.
class FastChem_output(object):
    def __init__(self,p,plot_species=['Fe','C1O1','H2O1','H'],outpath=''):
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
        #All colour stuff:
        self.colors = []#If this has nonzero length, lines will cycle through this list.
        self.styles = ['-']
        self.Tcolor = 'red'
        #All linestyle stuff:
        self.Tstyle = '--'
        self.plot_TP = True

        self.plot()

    def print_species(self):
        for i in self.species:
            print(i)

    def resolve_label(self,label):
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

        #
        #
        #
        # elif any(char.isdigit() for char in outlabel) == True:
        #     print('%s is a molecule' % label)
        #
        #
        #
        # if len(label) <= 2 and label[-1].isdigit() == False:
        #     print('%s is a neutral atom' % label)
        # elif len(label) <= 4 and label.count('+') == 1:
        #     print('%s is a single ion' % label)
        #
        # if len(outlabel) <= 2 and outlabel[-1].isdigit() == False and outlabel[-1] not in ['-','+']:#We are a single neutral atom species.
        #     print('%s is a neutral atom' % outlabel)
        # elif any(char.isdigit() for char in outlabel) == True:#If not, are we a molecule?
        #     for i in range (99,0,-1):
        #         if i > 1:
        #             outlabel=outlabel.replace(str(i),'$_'+str(i)+'$')#Paste LaTeX into the label for each character found.
        #             #We go from 99 downwards in order to take into account hypothetical molecules with double-digit numbers.
        #         else:
        #             outlabel=outlabel.replace(str(i),'')#Remove 1's.
        #     #Now we still have potential plusses and minuses. We have also all the atomic ions, which are like Ab+ now because the 1's are removed.
        #     #But these are all short:
        #
        # else:
        #     if any((char == '+') for char in outlabel) == True:
        #         outlabel+=r'{\fontsize{%spt}{3em}\fontfamily{cursif}\selectfont{}{ II}'%int(math.floor(self.fontsize*0.85))
        #     elif any((char == '-') for char in outlabel) == True:
        #         outlabel=outlabel.replace('-','$^-$')
        #     else:
        #         outlabel+=r'{\fontsize{%spt}{3em}\fontfamily{cursif}\selectfont{}{ I}'%int(math.floor(self.fontsize*0.85))
        return(outlabel)



    def plot(self,show=False):
        import matplotlib.pyplot as plt
        from matplotlib import rc

        for s in self.plot_species:
            if s not in self.species:
                self.print_species()
                print('')
                print("Error: %s is not available. Above are the available species."%s)
                return

        plt.rcParams.update({'font.size': self.fontsize})
        rc('text', usetex=True)

        fig, ax1 = plt.subplots()
        P = self.data['P (bar)']
        lines=[]
        for i,name in enumerate(self.plot_species):
            label = self.resolve_label(name)
            if len(self.colors) == 0:
                l,=ax1.plot(self.data[name],P,self.styles[i%(len(self.styles))],label=label)
            else:
                l,=ax1.plot(self.data[name],P,self.styles[i%(len(self.styles))],label=label,color=self.colors[i%(len(self.colors))])
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
        labels = [l.get_label() for l in lines]
        plt.legend(lines, labels, loc=0)
        if show == True:#If show == True, always show.
            plt.show()
        elif len(self.outpath) > 0:#Otherwise, check that the outpath is set, if so, we save.
            plt.savefig(self.outpath,dpi=self.dpi)
        else:#If show == False but the outpath is not set (default), we still show:
            plt.show()


a=FastChem_output('output/chem_output.dat')
