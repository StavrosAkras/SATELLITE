# Introduction

The Spectroscopic Analysis Tool for intEgraL fieLd unIt daTacubEs
(<span>Satellite</span>) is a newly developed code for the spectroscopic
characterization of extended photo-ionised nebulae such as planetary
nebulae, H II regions or galaxies in the optical regime observed with
any integral field unit (IFU). SATELLITE has been written in
<span>python</span> and was developed to be an automatic and
user-friendly code. The user does not need any programming or coding
knowledge.

The capabilities and performance of the <span>satellite</span> code
(v0.1) were first presented in  using the IFU data of the Abell 14 PN
obtained with VIMOS@ESO. <span>satellite</span> v1.2 has been applied to
three more PNe: Hen 2-108 (VIMOS) (Miranda Marques et al. 2021,
submitted), NGC 7009 (MUSE), and NGC 6778 (MUSE) (Akras et al. 2021,
submitted).

<span>satellite</span> carries out a spectroscopic analysis of extended
ionized nebulae through 1D and 2D approach on a list of 35 emission line
(the brightest and more frequently detected in ionized nebulae) via a
number of pseudo-slits that simulate slit spectrometry and 2D emission
line imaging. The analysis is performed in four different modules:

  - (I) rotation analysis,

  - (II) radial analysis,

  - (III) specific slits analysis and

  - (IV) 2D analysis

For each module, <span>satellite</span> computes all the typically used
nebular parameters and their uncertainties using as input information
the available emission lines maps. For all four modules, the
uncertainties of the line intensities as well as those of the nebular
parameters are computed following a Monte Carlo approach considering a
number of spectra.

  - extinction coefficient (c(H\(\beta\))),

  - electron temperatures and densities for different diagnostic lines,

  - ionic and elemental abundances, abundances ratio relative to oxygen
    and the ionization correction factors,

  - emission line ratios from a pre-defined list

The input parameters necessary to run the code are provided by the user
in four ASCII files:

  - <span>*input.txt*</span>

  - <span>*numerical\_input.txt*</span>

  - <span>*output.txt (the requested outputs from the code\!)*</span>

  - <span>*diagnostic\_diagrams\_input.txt*</span>

  - <span>*plots\_parameters\_input.txt*</span>

and they are:

  - emission line and error maps (or additional error as percentage of
    its pixel flux. )

  - the pixel scale of the IFU

  - the interstellar extinction law

  - number of replicate spectra

  - the width, length, position angle (PA) and coordinates of the
    pseudo-slits for the 1D analysis

  - the coordinates of the central star or the centre of the nebula

  - atomic data

  - Te and Ne from different diagnostic lines for the proper calculation
    of the ionic abundances

  - the line ratios that the code will compute

  - the emission line diagnostic diagrams that the code will construct

  - the module or modules that the code will execute

The philosophy behind the development of the <span>satellite</span> code
is, besides the unique 2D imaging spectroscopy that IFU technology
provides, to carry out a detailed 1D spectroscopic analysis through a
number of pseudo-slits that simulate slit spectrometry and emission line
imaging in order to properly compare the results presented in previous
studies. The outputs from each of <span>satellite</span>’s modules are
saved in five different folders:

  - output\_angles\_plots \(\leftarrow\) rotation analysis module

  - output\_Diagnostic\_Diagrams \(\leftarrow\) 2D analysis module

  - output\_images \(\leftarrow\) 2D analysis module ( & slit position
    testing)

  - output\_plots \(\leftarrow\) 2D analysis module

  - output\_radial\_plots\(\leftarrow\) radial analysis module

The description of each module as well as the input parameters and
outcomes are described in more details in the following sections using
as a representative example the analysis of the planetary nebulae
NGC 7009 and its MUSE data from the Science Verification phase .

# Set up and Run the SATELLITE Code

<span>satellite</span> has been successfully run/tested in Fedora 23
operation systems and <span>Python</span> version 2.7.11.  
To use the <span>satellite</span> code, it is also necessary to install
a number of libraries such as <span>matplotlib, numpy, scipy,
astropy</span>, and <span>seaborn</span> libraries.  
The version of the libraries was: matplotlib v2.2.5, numpy v1.11.1 and
scipy v0.18.0, astropy v2.0.12, and seaborn v0.9.1.
<span>satellite</span> also make use of the <span>PyNeb</span> package
version 1.1.15 .  
The installation of a library can be done via <span>*pip*</span> and the
command: <span>*pip install library.*</span>

To get the version of the installed libraries add the following lines in
a <span>python</span> script and run it.

\(\\
import\ matplotlib\\
import\ numpy\\
import\ seaborn\\
import\ astropy\\
import\ scipy\)  
\(\\
print(matplotlib.\_\_version\_\_)\\
print(numpy.\_\_version\_\_)\\
print(seaborn.\_\_version\_\_)\\
print(astropy.\_\_version\_\_)\\
print(scipy.\_\_version\_\_)\\
\)

For the atomic data from the Chianti database, it is necessary to
download the data from the website <https://www.chiantidatabase.org/>.

## GitHub repository

1\) Download <span>satellite</span> from GitHub
(https://github.com/StavrosAkras/SATELLITE)  
command: ‘‘<span>**git clone
https://github.com/StavrosAkras/SATELLITE.git**</span> ’’  
A new folder will be created namely SATELLITE with all the necessary
files.  
2\) run setup.py in the satellite folder  
command: ‘‘<span>**python setup.py install**</span> ’’  
Setup.py also installs all the necessary python packages need to use
<span>satellite</span> included <span>pyneb</span>.  
3\) run the code command: ‘‘<span>**./satellite \(>\)
outputLog.txt.**</span> ’’  
A number of general comments provided from the different modules of the
code are written in the <span>*outputLog.txt*</span> ASCII file.

# READ INPUT DATA

The number of input data that the user has to provide the code.

## read\_input\_script

In the current version of <span>satellite</span>, there is a pre-defined
list of 35 emission lines from which the user can select the emission
line that will be used by the code for the analysis. The list contains
the most commonly detected and relatively bright emission line in
ionized nebula such as the H and He recombination lines (recombination
lines from O or N will be implemented to future versions) and
collisionaly excited lines from O+, O++, N+, N++ ,S+, S++, etc. (see
Figure [1](#fig1)).

![Example of the <span>*input.txt*</span> file from which the user can
select the emission line that will be used for 1D spectroscopic analysis
or radial analysis.](input_lines.pdf)

The user can access the list from the <span>*input.txt*</span> file and
just add ‘‘ yes ’’ or ‘‘ no ’’ in the second column for both flux and
error maps. The first column in the <span>*input.txt*</span> lists the
name of the emission lines. <span>**The input "FIT files" of the
emission line map must have the same name**</span>. Moreover, all the
input "FIT files" (maps) must be located in the
<span>*image\_data*</span> folder. <span>**Hint:**</span> In case error
maps are not available, the user can create some fake error maps by
multiplying the flux maps by e.g. 0.1 to replicate an uncertainty of 10
percent for each emission line. The option of an extra error is also
possible. In the forth column, the user can add an extra error for each
line as a percentage of the flux (e.g. 1% of the total flux of
\[O III\] 5007Å line). The even rows for each line (error lines) can
take two numerical values "0" and "any integer". (i) "O" means that the
final error for each pixel and each emission line is equal to the
percentage of the flux, and (ii) "any number" means that the final error
for each pixel and each emission line is equal to the error from the
error map + the percentage of the flux.

The third column in the <span>*input.txt*</span> deals with the
<span>*radial\_analysis*</span> module and the user adds,
‘‘ radial\_yes ’’ or ‘‘ radial\_no ’’, for the lines(with their
error maps) are available and will be used for this module.
<span>**Note:**</span> The H\(\alpha\) and H\(\beta\) lines must always
be selected for the determination of the interstellar extinction
coefficient and corrected line intensities.

The <span>read\_input\_script</span> reads the flux and error maps of
each line defined in the <span>*input.txt*</span> file pixel-by-pixel
and save the values in 2D arrays.

For the example of NGC 7009, the selected emission lines are presented
in Figure. [1](#fig1)

## read\_input\_lines\_parameters\_script

Besides the selection of the emission line and error maps, there is also
a list of numerical parameters that the code needs to run properly and
they are provided by the user in different files. All these parameters
are read by the code via the
<span>read\_input\_lines\_parameters\_script</span>. Below, the most
important parameters are listed together with the ascii files:

  - the pixel scale of the IFU in arcsec (x10\(^{2}\)) \(\rightarrow\)
    <span>*numerical\_input.txt*</span>

  - the interstellar extinction law (Rx10\(^{1}\)) \(\rightarrow\)
    <span>*numerical\_input.txt*</span>

  - atomic data \(\rightarrow\) <span>*numerical\_input.txt*</span>

  - the width, length, position angle (PA) and coordinates of the
    pseudo-slits for the 1D analysis in the
    <span>*rotation\_analysis*</span>, <span>*radial\_analysis*</span>,
    and <span>*specific\_slits\_analysis*</span> modules \(\rightarrow\)
    <span>*numerical\_input.txt*</span>

  - the coordinates of the central star or the centre of the nebula
    \(\rightarrow\) <span>*numerical\_input.txt*</span> file

  - number of replicate spectra for the determination of the
    uncertainties \(\rightarrow\) <span>*numerical\_input.txt*</span>
    file

  - the minimum radius from which the maximum line flux will be
    determined and the profiles will be normalized to 1. \(\rightarrow\)
    <span>*numerical\_input.txt*</span> file

  - the number of column and row pixels that have to be added to the raw
    maps in order to put the centre of the nebula at the center of the
    map \(\rightarrow\) <span>*numerical\_input.txt*</span>  

  - the module or modules that the code will execute \(\rightarrow\)
    <span>*outputs.txt*</span> file

  - the emission line ratios that the code will compute based on the
    available line maps \(\rightarrow\) <span>*outputs.txt*</span>

  - the physical parameters (c, T\(_{\rm e}\), N\(_{\rm e}\)) that the
    code will compute based on the available line maps \(\rightarrow\)
    <span>*outputs.txt*</span>

  - the T\(_{\rm e}\) and N\(_{\rm e}\) from specific diagnostic lines
    that will be used to compute ionic abundances \(\rightarrow\)
    <span>*outputs.txt*</span>

  - the elemental abundances and ICFs that the code will compute based
    on the available line maps \(\rightarrow\)
    <span>*outputs.txt*</span>  

  - the diagnostic diagrams that the code will construct \(\rightarrow\)
    <span>*diagnostic\_diagrams\_input.txt*</span>

  - the xmin/xmax and ymin/ymax for the diagnostic diagrams
    \(\rightarrow\) <span>*diagnostic\_diagrams\_input.txt*</span>

  - if the user wants to over-plot the selection criteria from
    Kewley2001 and Kauffmann2003 \(\rightarrow\)
    <span>*diagnostic\_diagrams\_input.txt*</span>

  - ymin/ymax for the scatter plots \(\rightarrow\)
    <span>*plots\_parameters\_input.txt*</span> file

## Reorder line flux and error maps

The first task the user must perform before run <span>satellite</span>
is to compute the number of columns and rows of pixels that need to be
added to the raw line flux and error maps in order to put the central
star or centre of the nebula at the centre of the new map. This task is
very important in order to properly rotate the maps and measure the line
fluxes in each pseudo-slit. <span>**Note:**</span> <span>**The final
maps have to have the same number of columns and rows\!**</span>

Figure [2](#fig2) presents a cartoon that illustrates the numbers of
rows of pixels need to be added to the top and bottom of the raw map as
well as the columns of pixels to the left and right parts of the raw
maps. The value these extra pixels have is <span>**zero**</span> and
they do not have any impact on the spectroscopic analysis. The extra row
and columns added in the flux and error maps are provided by the user in
the <span>*numerical\_input.txt*</span> file.

CM (15,15) is the centre of the raw map and CS (19,12) is the central
star of an ellipsoidal planetary nebula (orange colour). The offset
between the CS and CM is 4 pixels in row and 3 pixels is column. It is
necessary to move the CS to the left and up by adding some extra columns
and rows. So, 8 columns are added to the right side and the CS is now in
the middle (19,19) and 6 rows at the bottom so the CS is now at the
position (18,18). But the map still does not have the same number of
columns and rows (38,36). For this reason, we have to add 1 extra row to
the top and 1 extra to the bottom. Finally, the CS is the centre of the
new map at position (19,19). These numbers are added in the
<span>*numerical\_input.txt*</span> file as it is show below:

  - add\_pixels\_above=1

  - add\_pixels\_below=7

  - add\_pixels\_left=0

  - add\_pixels\_right=8

  - total\_num\_pixels\_verti=30

  - total\_num\_pixels\_horiz=38

As for NGC 7009, the number of rows to the top and bottom of the raw
maps are 35 and 35, respectively, and the extra columns to the left and
right side are 4 and 6, respectively.

![Illustrative example of a line flux map and the rows/columns that have
to be added in order to coincide the CS with the CM of the new
map.](reorder_map_example.pdf)

## Extinction Coefficient

<span>satellite</span> determines the extinction coefficient
(c(H\(\beta\))) using the <span>PyNeb</span> package  (see: 
<https://github.com/Morisset/PyNeb_devel/blob/master/docs/Notebooks/PyNeb_manual_5.ipynb>)
and constructs the 2D map and its error map. All the extinction law
available in the the <span>PyNeb</span> can be selected (’No
correction’, ’CCM89’, ’CCM89\_Bal07’, ’CCM89\_oD94’,
’S79\_H83\_CCM89’, ’K76’, ’SM79\_Gal’, ’G03\_LMC’,
’MCC99\_FM90\_LMC’, ’F99’, ’F88\_F99\_LMC’\]). The R parameter is
also given by the user (e.g. R=3.1) in the
<span>*numerical\_input.txt*</span> file as an integer number
(R\*10=031). The <span>save\_FITSimages\_script</span> saves the 2D
arrays as FITS image.

The code does not take into account all the pixels in the pseudo-slits
or estimates the extinction coefficient for all the pixels
<span>**BUT**</span> only for those that satisfy the following criteria:

  - F(\(Ha\))\(>\)0,

  - F(\(Hb\))\(>\)0

  - F(\(Ha\))\(>\)F(\(Hb\))\*2.86

  - <span>**otherwise a value ’’zero‘‘ is applied.**</span>

There is also an option to find the outliers and exclude them from the
calculation of the mean/median values. The outliers are defined as those
pixels with values that do not satisfy the criteria:
value\(>\)per25-1.5\*iqr and value\(<\) per75+1.5\*iqr), where per25 and
per75 are the 25% and 75% percentiles, respectively, and iqr = per75 -
per25 the inter-quartile-range. These calculations are made in the
<span>calculations\_excluding\_outliers\_script</span>. In the current
version, the outliers are not excluded from the physical parameters that
the code computes.

## Atomic Data

The atomic data can also be changed by the user. All the options
available in the <span>PyNeb</span> package can be used by adding in the
<span>*numerical\_input.txt*</span> the words: IRAF\_09\_orig, IRAF\_09,
PYNEB\_13\_01, PYNEB\_14\_01, PYNEB\_14\_02, PYNEB\_14\_03’,
’PYNEB\_16\_01, , PYNEB\_17\_01, PYNEB\_17\_02, PYNEB\_18\_01,
PYNEB\_20\_01, and PYNEB\_21\_01. The atomic data from Chianti group can
also be used (adding the word ’’Chianti‘‘in the
<span>*numerical\_input.txt*</span>) but they have to be downloaded from
the website ([
http://www.chiantidatabase.org/chianti\_download.html.](%20http://www.chiantidatabase.org/chianti_download.html.)).
For more information the user must refer to the website of PyNeb package
(<https://github.com/Morisset/PyNeb_devel/blob/master/docs/Notebooks/PyNeb_manual_3.ipynb>)

# SATELLITE’s Modules

The modules that the code executes are defined in the second column
(‘‘ yes ’’ or ‘‘ no ’’) of the <span>*outputs.txt*</span> file.
Each module is described below. <span>**Note:**</span> It is recommended
to run the <span>radial\_analysis</span> module separately from the
other three modules.

## rotation\_analysis module

The <span>*rotation\_analysis*</span> module deals with the
spectroscopic analysis of a number of radial placed pseudo-slits from
the centre to outer parts with position angles (PA) between 0 and 360.
Figure [3](#fig4) illustrates as example the position of these
pseudo-slits on the \[N II\] 6584 emission line map of NGC 7009. The
minimum and maximum values of the PA, the step in PA , the width and the
length of the pseudo-slits are provided by the user in the file
<span>*numerical\_input.txt*</span>.

![An illustrative image of the radial positioned pseudo-slits with PA
from 0 to 360 every 20 degrees overlaid the \[N II\] 6584 image of
NGC 7009.](NGC709rotate_slits.png)

The code first computes the integrated H\(\beta\) flux and line fluxes
in the <span>rotate\_line\_fluxes\_script</span>. Then, line intensities
(relative to H\(\beta\) and corrected for the interstellar excitation)
as well as all the nebular parameters defined by the user in the
<span>*output.txt*</span> file for all the pseudo-slits are computed by
the <span>PyNeb</span> package. This analysis allows to explore the
variation of line intensities, line ratios and physical parameters (Te,
Ne), chemical abundances as functions of PA.

The outcomes from this module are multiple and are saved in different
files:

  - an ASCII file with the c(H\(\beta\)) and the intensity of each
    emission line and for each PA \(\rightarrow\)  
    <span>*output\_linesintensities\_per\_angles.txt*</span>

  - an ASCII file with various line ratios defined by the user for each
    PA \(\rightarrow\)<span>*output\_lineratios\_per\_angles.txt*</span>

  - an ASCII file with T\(_{\rm e}\) and N\(_{\rm e}\) defined by the
    user for each PA
    \(\rightarrow\)<span>*PyNeb\_output\_Te\_and\_Ne\_per\_angles.txt*</span>

  - an ASCII file with the ionic abundances for each ion defined by the
    user for each PA \(\rightarrow\)  
    <span>*PyNeb\_output\_total\_abund\_ICFs\_per\_angles.txt*</span>

  - an ASCII file with the elemental abundances and ICFs for each
    element defined by the user for each PA
    \(\rightarrow\) <span>*PyNeb\_output\_ionic\_abund\_per\_angles.txt*</span>

  - plots of c(H\(\beta\)), T\(_{\rm e}\), N\(_{\rm e}\), ionic,
    elemental abundances, ICFs and abundances ratios as function of the
    PA                                                   \(\rightarrow\)
      <span>*output\_angles\_plots*</span> folder

The user can use the ASCII files to carry out any further analysis
and/or construct his/her own proper plots.

Figure [6](#fig6a) and [9](#fig6b) present the plots of c(H\(\beta\)),
T\(_{\rm e}\), N\(_{\rm e}\), ionic, elemental abundances and ICFs of N
as functions of the position angle of the pseudo-slits for the analysis
of NGC 7009. T\(_{\rm e}\) and N\(_{\rm e}\) are shown in the same plot
as in Figures [6](#fig6a) and [9](#fig6b) or in two separate plots.

![Representative plot of c(\(H\beta\)) (upper panel), T\(_{\rm e}\) and
N\(_{\rm e}\) for two different diagnostic lines (middle/lower panels)
as function of the PA of the pseudo-slits.](fig_cHbeta_angles.png)
![Representative plot of c(\(H\beta\)) (upper panel), T\(_{\rm e}\) and
N\(_{\rm e}\) for two different diagnostic lines (middle/lower panels)
as function of the PA of the
pseudo-slits.](fig_Te\(NII6548_84\)_Ne\(SII6716_31\)_both_angles.png)
![Representative plot of c(\(H\beta\)) (upper panel), T\(_{\rm e}\) and
N\(_{\rm e}\) for two different diagnostic lines (middle/lower panels)
as function of the PA of the
pseudo-slits.](fig_Te\(SIII6312_9069\)_Ne\(ClIII5517_38\)_both_angles.png)

![Representative plot of the ionic/total abundance of N (upper panel),
the ICF(N) (middle panel) and the N/O ratio as functions of the
pseudo-slits’ PA.](N_abundances_angles.png) ![Representative plot of the
ionic/total abundance of N (upper panel), the ICF(N) (middle panel) and
the N/O ratio as functions of the pseudo-slits’ PA.](N_ICFs_angle.png)
![Representative plot of the ionic/total abundance of N (upper panel),
the ICF(N) (middle panel) and the N/O ratio as functions of the
pseudo-slits’ PA.](NO_abundances_ratio_angles.png)

### slit\_line\_flux\_script

This script calculates the flux for each emission line along a
pseudo-slit with a specific position angle, width and length given by
the user in the <span>*numerical\_input.txt*</span> file.
<span>**Note1:**</span> The width should always be an integer number.
The pseudo-slit starts from the centre of the image or of the nebula,
and it covers only one half of the nebula. The code sums up the values
of all the pixels within the area defined by the width and length of the
pseudo-slit, except those pixels which have F(\(Ha\))\(<\)0,
F(\(Hb\))\(<\)0 and/or F(\(Ha\))\(<\)F(\(Hb\))\*2.86 (negative or
unrealistic c(Hb)).

For the calculations, the script first rotates the entire image/table by
a given angle, then calculates the new size (x,y) of the rotated imaged
as well as the center of the new image. The pseudo-slit is always
oriented along the up-down direction of the image. The script computes
the total flux for each emission line and the total number of pixels.

<span>**Note2:**</span> It is necessary the image be large enough to be
sure that after the rotation, the entire nebula or galaxy remains inside
the image. Sometimes an elongated nebula (or even a galaxy) is observed
in a specific PA, so the entire nebula fits the field of view of the
instrument (<span>**rotation angle=0 in the <span>satellite</span> code
corresponds to observed PA=0**</span>).  
<span>**Note3:**</span> The orientation of the maps/images should always
be north up and east to the left. Otherwise the user has to take into
account the offset between the sky (North) and image orientation.

In case, the slit width and length are equal to the number of the pixels
in the raw map (parameter total\_num\_pixels\_horiz in the
<span>*numerical\_input.txt*</span> file), the code calculates the
integrated fluxes of the entire nebula for all the position angles. This
specific task was used to verify if the rotation of the images affects
the integrated line fluxes.

If the slit width and length are larger than the maximum number the
software return the following message: ‘‘Sorry, your slit width or/and
length are larger than the true size of the image’’.

### TeNe\_angles\_script, ionicabundances\_angles\_script and element\_abundances\_ICFs\_angles\_script

T\(_{\rm e}\), N\(_{\rm e}\), ionic/elemental abundances, ICFs and
abundance ratios are also computed for each pseudo-slit. Various
diagnostic lines can be used for T\(_{\rm e}\)/N\(_{\rm e}\). The user
can also choose which T\(_{\rm e}\)/N\(_{\rm e}\) combination will be
applied for the abundances of each ion (see Figure [11](#fig7)). All
these parameters are defined by the user in the
<span>*outputs.txt*</span> file. T\(_{\rm e}\) and N\(_{\rm e}\)
parameters are computed in the <span>*TeNe\_angles\_script*</span>,
while the ionic,elemental abundances and ICFs are computed in the
<span>*ionicabundances\_angles\_script*</span> and
<span>*element\_abundances\_ICFs\_angles\_script*</span>. All the
scripts make use of the <span>*PyNeb*</span> package.

![An example of the <span>*outputs.txt*</span> file for NGC 7009 and the
parameters that the user has to define for the calculations of
T\(_{\rm e}\), N\(_{\rm e}\), ionic/elemental abundances, ICFs and
abundance ratios.](outputsfile.pdf) ![An example of the
<span>*outputs.txt*</span> file for NGC 7009 and the parameters that the
user has to define for the calculations of T\(_{\rm e}\), N\(_{\rm e}\),
ionic/elemental abundances, ICFs and abundance
ratios.](output_abund.pdf)

## specific\_slit\_analysis module

In this module, the user can define 10 pseudo-slits for a spectroscopic
analysis of specific regions/structures in PN (e.g. knots, blobs, inset
or outer regions) or in any extended nebula. All the input information
from the 10 specific pseudo-slits are given by the user in the
<span>*numerical\_input.txt*</span> file:

![An example of the <span>*numerical\_input.txt*</span> file for
NGC 7009 and the parameters that the user has to define for the 10
pseudo-slits in the <span>*specific\_slit\_analysis*</span>
module.](spec_sl_output.pdf)

  - PA\_for\_specific\_slit\_n (in angle)

  - width\_for\_specific\_slit\_n (in pixels)

  - length\_for\_specific\_slit\_n (in pixels)

  - x\_coor\_of\_spec\_slit\_n (in pixels)

  - y\_coor\_of\_spec\_slit\_n (in pixels)

where n is the number of the pseudo-slit from 1 to 10 (see
Figure [12](#fig8)). The x\_coor\_of\_spec\_slit\_n, and
y\_coor\_of\_spec\_slit\_n parameters refer to the centre of each slit.

![Ten selected regions in NGC 7009 overlaid on the \[N II\] 6584Å image.
The position of the centre (x, y coordinates), position angle, width and
length of the slits are free parameters provided by the user. Slits 1
and 10 represent the slits position from  and , respectively. Numbered
regions from 2 to 9 correspond to the sub-structures of knots and
jet-like, or sub-regions of rims defined in . ](NGC7009regions.pdf)

<span>satellite</span> calculates the H\(\beta\) flux, line intensities
(normalized to Hb=100 and corrected for interstellar extinction),
emission line ratios (from the <span>*output.txt*</span> file), nebular
parameters (T\(_{\rm e}\), N\(_{\rm e}\)), ionic/total elemental
abundances, ICFs and abundance ratios for all 10 pseudo-slits. This
module is executed from the
<span>*specificPA\_line\_fluxes\_script*</span>. The scripts
<span>*slit\_line\_flux\_script*</span> is employed in this module for
each pseudo-slit.

c(H\(\beta\)), emission line intensities and line ratios are saved in
the <span>*output\_linesintensities\_per\_angles.txt*</span> and
<span>*output\_lineratios\_per\_angles.txt*</span> files, respectively.
So, the user can also perform any extra analysis he/she wants.
Similarly, T\(_{\rm e}\) and N\(_{\rm e}\) parameters are computed in
the <span>*TeNe\_specific\_slits\_script*</span> and are saved in the
<span>*PyNeb\_output\_Te\_and\_Ne\_specific\_slits*</span> file, the
ionic abundances are computed in the
<span>*ionicabundances\_specific\_slits\_script*</span> and are saved in
the <span>*PyNeb\_output\_ionic\_abund\_specific\_slits*</span> file,
finally and the elemental abundances, ICFs and abundance ratios are
computed in the
<span>*element\_abundances\_ICFs\_specific\_slits\_script*</span> and
are saved in the
<span>PyNeb\_output\_total\_abund\_ICFs\_specific\_slits</span> file.  
Figure [13](#fig9) shows the position of the specific areas/regions
selected for the study of NGC 7009 overlaid on the \[N II\] flux map.
The selected regions are the same as those defined by  for a direct
comparison of the results from the <span>*specific\_slit\_analysis
module*</span> with 1D long-slit spectroscopic data. Possible
differences between the two studies can be associated with the position
of the pseudo-slits.  
At this point, it is worth mentioning the
<span>*slit\_position\_testing*</span> module of the
<span>satellite</span> code. This module is used to verify the position
of the pseudo-slits before use the software. <span>**Hint:**</span> When
the <span>*slit\_position\_testing*</span> module is used first deselect
all other modules. Moreover, at least an emission line has to be used
and defined in the <span>*input.txt*</span> file in order to properly
use this module (e.g., H\(\alpha\) and H\(\beta\) to avoid multiple
maps). The output of this module is 10 figures (in png and pdf formats)
with the position of each pseudo-slit overlaid on the emission line map
(see Figure [15](#fig10)).

![Representative examples of the output figures produced by the
<span>*slit\_position\_testing*</span> module. The slit from  (left
panel) and R1 slit from  (right panel) are shown overlaid on the
H\(\alpha\) emission line maps. Scale is in "python counting; usual
pixel counting starts in 1.](fig_slit0.pdf) ![Representative examples of
the output figures produced by the
<span>*slit\_position\_testing*</span> module. The slit from  (left
panel) and R1 slit from  (right panel) are shown overlaid on the
H\(\alpha\) emission line maps. Scale is in "python counting; usual
pixel counting starts in 1.](fig_slit4.pdf)

# 2D analysis module

Besides the 1D spectroscopic analysis, <span>satellite</span> also
performs a spectroscopic analysis in both spatial dimensions
simultaneously using the entire maps. For this module, the
<span>analysis2D\_script</span>, <span>TeNe\_2D\_script</span>,
<span>generate\_2D\_lineratio\_maps\_script</span>,
<span>ionicabundances\_2D\_script</span>,
<span>element\_abundances\_ICFs\_2D\_script</span> and
<span>diagnotic\_diagrams\_script</span> are employed.

c(H\(\beta\)), line intensities, line ratios, T\(_{\rm e}\),
N\(_{\rm e}\), ionic, elemental abundances, ICFs and abundances ratios
are computed for each individual pixels, if the criteria
F(\(Ha\))\(>\)0, F(\(Hb\))\(>\)0 and F(\(Ha\))\(>\)F(\(Hb\))\*2.86 are
satisfied, otherwise a value equals to <span>**’’zero‘‘**</span> is
applied.

The main outcomes from this module are 2D maps for all the
aforementioned nebular parameters saved in the
<span>*output\_images*</span> folder. In
Figures [16](#fig14), [18](#fig15) and [20](#fig16), the maps of
c(H\(\beta\)), T\(_{\rm e}\) and N\(_{\rm e}\) using the \[S III\] and
\[S II\] diagnostic lines, the line ratios log(\[N II\]/\[O III\]) and
log(\[S II\]/S III\]) are presented as representative examples of the
outcomes from this module.

The <span>*2D\_spectroscopic\_analysis*</span> module also calculates
and returns the distribution of each maps (histogram plots), e.g.,
c(H\(\beta\)), N\(_{\rm e}\) and T\(_{\rm e}\) maps (see
Figures [21](#fig17) and [23](#fig18)) as well as emission line
diagnostic diagrams (see Figures [24](#fig19a) and [25](#fig19b)) using
the <span>diagnotic\_diagrams\_script</span> which are selected by the
user in the <span>*diagnostic\_diagrams\_input*</span> file (see
Figure [26](#fig13)).  
<span>**Note1:**</span> At this point, it has to be clarified that when
the <span>*specific\_slits\_analysis*</span> and/or
<span>*rotation\_analysis*</span> modules are used together with the
<span>*2D\_analysis*</span> module, the emission lines ratios for the
three modules are plotted on the same diagnostic diagrams for a direct
comparison between an 1D and 2D analysis.  
Last but not least, <span>satellite</span> computes and returns an ASCII
file with the mean value, standard deviation and the percentiles of 5%,
25% (Q1), 50% (median), 75% (Q3), 95% for all the nebular parameters and
emission line ratios for a thorough statistical analysis of the observed
nebula.

![c(H\(\beta\)) map of NGC 7009.](chb_2Dimage.png)

![N\(_{\rm e}\) and T\(_{\rm e}\) maps obtained from the \[S II\] and
\[S III\] diagnostic lines of
NGC 7009.](Ne\(SII6716_31\)_Te\(NII6548_84\)_2Dimage.png)
![N\(_{\rm e}\) and T\(_{\rm e}\) maps obtained from the \[S II\] and
\[S III\] diagnostic lines of
NGC 7009.](Te\(SIII6312_9069\)_Ne\(ClIII5517_38\)_2Dimage.png)

![Log(\[N II\]/\[O III\]) and log(\[S II\]/S III\]) line ratio maps of
NGC 7009.](log\(\(N2_6548s+N2_6583s\)_\(O3_4959s+O3_5007s\)\)_2Dimage.png)
![Log(\[N II\]/\[O III\]) and log(\[S II\]/S III\]) line ratio maps of
NGC 7009.](log\(\(S2_6716s+S2_6731s\)_\(S3_6312s+S3_9069s\)\)_2Dimage.png)

![The histogram of c(H\(\beta\) map.](histc_Hb.png)

![The histograms of N\(_{\rm e}\)\[S II\] and T\(_{\rm e}\)\[S III\]
maps.](histNe\(SII6716_31\)_Te\(NII6548_84\).png) ![The histograms of
N\(_{\rm e}\)\[S II\] and T\(_{\rm e}\)\[S III\]
maps.](histTe\(NII6548_84\)_Ne\(SII6716_31\).png)

![A representative example of emission line diagnostic diagram: \[S II\]
6716/6731 versus H\(\alpha\)/\[N II\] 6548+6584. Cyan dots correspond to
the values of individual pixels, pink circles and yellow diamonds show
the values obtained from the simulated long-slits of the rotational
analysis task with position angles from 0 to 360 degrees with 10 degrees
step and the values from the 10 simulated slits in the specific slits
task, respectively. The inset plot illustrate the variation of the line
ratios with the position angle of the simulated
slits.](fig_SII_HaNII.png)

![A representative example of emission line diagnostic diagram:
\[O III\] 5007/H\(\beta\) versus \[S II\] 6716+6731/H\(\alpha\). Cyan
dots correspond to the values of individual pixels, pink circles and
yellow diamonds show the values obtained from the simulated long-slits
of the rotational analysis task with position angles from 0 to 360
degrees with 10 degrees step and the values from the 10 simulated slits
in the specific slits task, respectively. The inset plot illustrate the
variation of the line ratios with the position angle of the simulated
slits. The regimes of the PNe, H <span>ii</span> regions and supernova
remnants are also drawn.](fig_OIIIHb_SIIHa.png)

![An example of the <span>*diagnostic\_diagrams\_input.txt*</span> file
for NGC 7009 and the parameters that the user can select for the
diagnostic diagrams in the <span>*2D\_analysis\_module*</span>
module.](DD_input.pdf)

## radial analysis module

The final module in the current version of <span>satellite</span> (v1.2)
conducts a radial spectroscopic analysis considering a pseudo-slit with
specific width, length and position angle
(parameter=angle\_for\_radial\_flux) provided by the user in the
<span>*numerical\_input.txt*</span> file. The user must also select the
emission lines that will be used for this analysis. This can be made in
the third column of the <span>*input.txt*</span> file: radial\_yes, or
radial\_no.  
<span>**Note1:**</span> It is recommended to disable all other modules
when the <span>*radial\_analysis*</span> module is executed.  
The main outcomes of this module are:

  - (I) the radial profiles of all the selected emission lines in the
    <span>*input.txt*</span> normalized by the peak flux.

  - (II) the calculation of all the nebular parameters (c(H\(\beta\)),
    line intensities, line ratios, T\(_{\rm e}\), N\(_{\rm e}\), ionic,
    elemental abundances, ICFs and abundances ratios) as functions of
    the distance from the central star or the central point of the
    nebula or galaxy in general.

The normalization of the radial profiles is made using the peak of the
flux of each emission line. However, the user can also select the range
from where this peak can be obtained by providing the code with the
minimum radius (<span>*limit\_radial\_in\_arcsec*</span> parameter).
This option permit to investigate the radial distribution of emission
lines for regions/substructures with specific interest.

The radial profile of various emission lines for the example of NGC 7009
are shown in Figure [28](#fig11) (<span>**Hint:**</span> It is
recommended to use maximum 4-6 lines for the construction of more
illustrative plots.). All radial profiles are normalized to a peak flux
found for distances r\(>\)20 arcsec (\_radial\_in\_arcsec\(>\)20)
focused to the low-ionization structures/knot of NGC 7009. The
calculation are made in the <span>*find\_maxvlaue\_script*</span>.
Hence, <span>satellite</span> returns the distance between the peak of
each selected line and the central star in arcsec. Table 1 lists the
distances for the example of NGC 7009 and it can be seen that there is a
spatial offset of 1 arcsec between the high/moderate- and low-ionization
lines. The values of each radial step (pixel scale of the IFU) are also
saved in an ASCII file, so the user can build his/her own radial
profiles.

The radial variation of c(H\(\beta\)), T\(_{\rm e}\), and N\(_{\rm e}\)
parameters of NGC 7009 are shown in Figure [31](#fig12)

![Radial profiles for several emission lines of NGC 7009 at
PA=79 degrees. Upper panel shows all the radial profiles, while the
lower panel zooms-in to the much fainter emission
lines.](radial_distribution_1.png) ![Radial profiles for several
emission lines of NGC 7009 at PA=79 degrees. Upper panel shows all the
radial profiles, while the lower panel zooms-in to the much fainter
emission lines.](radial_distribution_2.png)

![Representative examples of the radial distribution of c(H\(\beta\))
upper panel and T\(_{\rm e}\), N\(_{\rm e}\) (lower
panel).](c\(Hb\)_radial.png) ![Representative examples of the radial
distribution of c(H\(\beta\)) upper panel and T\(_{\rm e}\),
N\(_{\rm e}\) (lower
panel).](Te\(NII6548_84\)_Ne\(SII6716_31\)_both_radial.png)
![Representative examples of the radial distribution of c(H\(\beta\))
upper panel and T\(_{\rm e}\), N\(_{\rm e}\) (lower
panel).](Te\(SIII6312_9069\)_Ne\(ClIII5517_38\)_both_radial.png)

One again, it has to be pointed out that the code sums up the values of
the pixels which have F(H\(\alpha\))\(>\)0, F(H\(\beta\))\(>\)0 and/or
F(H\(\alpha\))\(>\)F(H\(\beta\))\*2.86.

<div id="distancepeak">

|                          |                     |                         |          |
| :----------------------: | :-----------------: | :---------------------: | :------: |
|           Line           | distance\(^{\dag}\) |          Line           | distance |
|                          |      (arcsec)       |                         | (arcsec) |
|  H <span>i</span> 4861Å  |        23.6         |     \[N II\] 6548Å      |   24.8   |
|  \([\)O III\(]\) 4959Å   |        23.8         | H <span>i</span> 6563Å  |   23.6   |
|   \([\)N I\(]\) 5199Å    |        24.8         |     \[N II\] 6584Å      |   24.8   |
| He <span>ii</span> 5412Å |        20.2         | He <span>i</span> 6678Å |   23.8   |
|  \([\)Cl III\(]\) 5517Å  |        24.2         |     \[S II\] 6716Å      |   24.8   |
|  \([\)Cl III\(]\) 5538Å  |        24.4         |     \[S II\] 6731Å      |   24.8   |
|     \([\)N I\] 5755Å     |        24.2         |    \[Ar III\] 7136Å     |   24.8   |
| He <span>i</span> 5876Å  |        23.6         |     \[O II\] 7320Å      |   24.8   |
|   \([\)O II\(]\) 6300Å   |        24.8         |     \[O II\] 7330Å      |   24.8   |
|  \([\)S III\(]\) 6312Å   |        23.8         |     \[S III\] 9069Å     |   23.8   |

Distances from the central stars of emission line’s peak for a
pseudo-slit at 79 degree position angle

</div>

\(^{\dag}\) The spacial resolution of MUSE maps is 0.2 arcsec.

### radial distance calculations

At this point, it is necessary to further explain how
<span>satellite</span> calculates the fluxes of the emission lines as
function of the distance from the central star. Figure [32](#fig20)
shown an example of a pseudo-slit at PA=90 degrees.

The width of the pseudo-slit defines how many pixels will be taken into
consideration for the flux at each distance. For instance, the fluxes
(and errors) of 5 pixels are summed up for the first column (or radial
distance r=0.2 arcsec). Then, the code moves to the second column and
computes the flux and the corresponding error from the next 5 pixels at
the radial distance r=0.4 arcsec and so on (see Figure [32](#fig20)).

After finishing the computation of the fluxes and errors for all the
lines, the code computes the extinction coefficient (c(H\(\beta\))) and
corrected line intensities (relative to H\(\beta\)=100) as function of
the radial distance from the central star or geometric centre as well as
all nebular parameters (T\(_{\rm e}\), N\(_{\rm e}\), ionic, elemental
abundances, ICFs and abundance ratios) (Figures [28](#fig11) and
[31](#fig12)).

![An illustrative example of how <span>satellite</span> computes the
fluxes and radial stances from the central star of the nebula or the
geometric centre of the nebula.](NGC7009radial.png)

# Uncertainty calculations

The uncertainties of emission lines and all nebular parameters for all
four modules are computed following the same Monte Carlo approach. In
particular, <span>satellite</span>, first, computes the total error of
the flux for each pseudo-slit or pixel, using the provided error maps +
an extra error as the percentage of the flux.

\[\Delta F_{tot}=\Sigma_{i=1}^{N}~(\sigma_{F_i}+\lambda*F_i)^{1/2}\]
where i corresponds to the pseudo-slit or pixel and range from 1 to the
total number N, \(\sigma_{F_i}\) is the uncertainties of the flux in the
pseudo-slit or pixel i based on the provided error maps, F\(_i\) is the
flux in the pseudo-slit or pixel i and \(\lambda\) corresponds to the
percentage of the flux (from 0 to 1.0). If \(\lambda\)=0, then the code
takes into account only the errors from the maps. The \(\lambda\)=0
parameter is given to the code by the user in the
<span>*input.txt*</span> file (forth column). There is also the option
not to use the error maps. This is defined in the
<span>*input.txt*</span> file (forth column) even rows (the rows of
errors). If a non-zero value is provided, the code uses the formula (1)
while for a "0" value , the code uses the formula (2).

\[\Delta F_{tot}=\Sigma_{i=1}^{N}~(\lambda*F_i)^{1/2}\]

These resultant uncertainties of the line fluxes are then used to
replicate the spectrum of a pseudo-slit or pixel, and a number of
additional spectra are generated using a Monte Carlo algorithm. The
number of the replicate spectra is given by the user in the
<span>*numerical\_input.txt*</span>. <span>satellite</span> computes all
the nebular parameters, e.g. T\(_{\rm e}\), N\(_{\rm e}\), ionic,
elemental abundances and ICFs for all the replicate spectra and the
standard deviation of each parameters is the uncertainty of the
parameter that the <span>satellite</span> code provides.

# General Notes:

  - It has to be clarified that the <span>satellite</span> code takes as
    input a list of emission line fluxes and error maps extracted from
    the datacubes obtained from any IFU. It does not extract the maps
    from the datacudes. Therefore, this is a step that has to be done
    before the use of <span>satellite</span>.

  - Moreover, <span>satellite</span> can also be applied to individual
    emission lines images obtained with narrow band filters (if there
    are available) or the emission line images obtained from 3D
    photo-ionization models.

# Possible error messages

In this section, a number of possible errors that may come out are
described.

  - In case an emission line is missing, the user has to define that in
    the <span>*input.txt*</span> file by writing "no" in the second and
    third columns of the corresponding line. If the user has forgotten
    to properly change the <span>*output.txt*</span> file or the
    <span>*diagnostic\_diagram\_input.txt*</span> file, an error will be
    return, see Figure [33](#figmissingline).

  - In case an emission line is missing, but the user has forgotten to
    properly change the <span>*output.txt*</span> file or the
    <span>*diagnostic\_diagram\_input.txt*</span> file and a physical
    parameters has to be computed such as Te, Ne, abundances, an error
    will be return like in Figure [34](#figmissingline2).

  - In case, the arrays of the emission lines have different sizes an
    error will be return like in Figure [35](#figerrorimages). Moreover,
    the error line may also be related to the parameters
    <span>*total\_num\_pixels\_verti*</span> and
    <span>*total\_num\_pixels\_horiz*</span> in the
    <span>*numerical\_input.txt*</span> file.

![An example of error in case the He II line is missing (in the
<span>*input.txt*</span> file, it has been set as "no") but the
He I/He II ratio is required to be computed (in the
<span>*output.txt*</span>, the He I/He II ratio is still
"yes").](linemissing.png)

![An example of error in case a physical parameter is required to be
computed and the corresponding line is missing (in the
<span>*input.txt*</span> file, it has been set as "no". In this example,
the T\(_{e}\) and N\(_{e}\) from the \[N II\] and \[Cl III\] diagnostic
lines have to be computed but a diagnostic line is
missing.](linemissing2.png)

![An example of error in case there is a problem with the dimensions of
the arrays that correspond to the emission lines. In this case, the
error is because the size of the flux maps is not consistent with the
<span>*total\_num\_pixels\_verti*</span> and
<span>*total\_num\_pixels\_horiz*</span> parameters in the
<span>*numerical\_input.txt*</span> file. ](error_with_tables.png)

# How to Cite <span>satellite</span> in a publication? 

There are two papers appropriate as references for
<span>satellite</span> in a publication. They are:  
1\) Akras, Stavros; Monteiro, Hektor; Aleman, Isabel; Farias, Marcos A.
F. ; May, Daniel ; Pereira, Claudio B., 2020, MNRAS, 493, 2238A  
  
bibtex code:  
@ARTICLE<span>2020MNRAS.493.2238A,  
author = <span><span>Akras</span>, Stavros and <span>Monteiro</span>,
Hektor and <span>Aleman</span>, Isabel and <span>Farias</span>, Marcos
A. F. and <span>May</span>, Daniel and <span>Pereira</span>, Claudio
B.</span>,  
title = "<span>Exploring the differences of integrated and spatially
resolved analysis using integral field unit data: the case of Abell
14</span>",  
journal = <span>\(\setminus\)mnras</span>,  
keywords = <span>techniques: imaging spectroscopy; techniques:
spectroscopic; (stars:) binaries: general, ISM: abundances, (ISM:)
planetary nebulae: individual: Abell 14, Astrophysics - Astrophysics of
Galaxies, Astrophysics - Solar and Stellar Astrophysics</span>,  
year = 2020,  
month = apr,  
volume = <span>493</span>,  
number = <span>2</span>,  
pages = <span>2238-2252</span>,  
doi = <span>10.1093/mnras/staa383</span>,  
archivePrefix = <span>arXiv</span>,  
eprint = <span>2002.12380</span>,  
primaryClass = <span>astro-ph.GA</span>,  
adsurl =
<span>https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.2238A</span>,  
adsnote = <span>Provided by the SAO/NASA Astrophysics Data
System</span>  
</span>  
2\) S. Akras; H. Monteiro; J. R. Walsh; J. García-Rojas; I. Aleman; H.
Boffin; P. Boumis; A. Chiotelis; R. M. L. Corradi; D. R. Gonçalves; L.
A. Gutiérrez-Soto; D. Jones; C. Morisset, 2021, MNRAS, submitted  

# Licence and Copyright Information

<span>satelite</span> is freely available under the General Public
License (GPL).

# Questions of problems

For questions please write an email to Dr. Stavros Akras
(stavrosakras@gmail.com)
