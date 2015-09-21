# Looking at local environment of Type Ia supernova hosts via Hubble Space Telescope

Distance is the hardest thing to measure in Astronomy. Finding a standard candle, an object of known brightness, is a simple way to solve this. Type Ia supernova are our standard candles. This project is an attempt to look at the relationship of local galactic properties, available through HST, and the known scatter of Ia absolute magnitude. 

This file is primarily for me to document what I have done, need to do, and my other thoughts. Secondary it is for others to look at this project and understand what is happening. For them, the most useful sections are [Procedure](#procedure), [Analysis Method](#analysis-method), and [Future Work](#future-work). 

## Table of Contents

* [Background](#background)
* [Project Goals](#project-goal)
* [Data](#data)
* [Procedure](#procedure)
* [Questions & Answers](#questions--answers)
* [Analysis Method](#analysis-method)
* [Future Work](#future-work)

## Background

* Mass step (Childress 2013 is recent)
* Mentality (Hayden 2013)
* Local H$_{\alpha}$ (Rigault 2013)
* Stuff in emails between Bryan & Peter (with me in there but not commenting.) In early spring 2015

## Project Goal

SN are used in cosmology and need to be a standard candle, we need to remove all variance in their absolute magnitude. Host mass & metallicity have been seen to effect the absolute magnitude of SN, but intuition says it should be a local property not a global galactic property.

## Data

### Data can be found at 

* HST data from proposals 
    * [11670](http://archive.stsci.edu/proposal_search.php?id=11670&mission=hst) for Cycle 17 and 
    * [12969](http://archive.stsci.edu/proposal_search.php?id=12969&mission=hst) for Cycle 20 
* SDSS [transient info](http://sdssdp62.fnal.gov/sdsssn/DataRelease/index.html)

### Data description

* .drz
    - final image after HST processing (astrodrizzel?)
* .flt
    - This comes with two copies. These can be looked at to see each actual observation. But these will have cosmic ray hits and other imperfections. 

### Data issues

| SN | Issue |
|----|----|
|SN13038| the host does not seem to exist. I am thinking about removing it form the analysis or performing its analyses by hand |
|SN19023| host seems to be far away |
|SN4019| does not seem to have a star in its field. I used a larger elliptical galaxy as apposed to a small (& faint) almost point sized unknown in type galaxy to calculate the WCS shift. |
|SN6491| has a star on top of the host |
|SN6614| might have a "star" that is really a "galaxy", check SDSS object ID 1237663785282764968. Also SN6614 uses a galaxy for its WCS shift calibration. |
|SN8297| has a galaxy very far away and I do not think I found it with the code. Code says that I found something at `r < 300 pixels` |

### Variables

This section is to keep me consistent

| name | variable | Units or range | description |
|------|:------:|-------|------|
|`Raw' pixel value | $\dfrac{e^{-1}}{s}$ | some flux or something | The value in the fits file for the processed (drz).
| fractional pixel rank | $r_f$ | [0-1] | The fractional rank of the pixel value of the SN[^fractional pixel rank] compared to the galaxy.|
| WCS shift | $\Delta_{WCS}$  | arcsec | The difference between the position a star in the .fits file and SDSS DR12. [More details below](#calculating-delta)

## Procedure

### Processing Data

#### HST Data rename & combining - `renameFITS.py`

+ has two functions `rename(doc)` & `merge(img1, img2)`
+ takes all HST data and renames in formate `data/HST - renamed/<<SN>>_<<filter>>_<<img>>_<<count>>.fits`
    * `<<SN>>` is the SDSS tranient number
    * `<<filter>>` is the HST filter (by number)
    * `<<img>>` is wether its drz or flt
    * `<<count>>` optional and only for flt, since there are two images per filter.
+ combines the two drz images
    * `<<filter>>` becomes `combined`
    * `<<img>>` and `<<count>>` are dropped
    * saves as `data/HST - combined/*.fits`

### Repeated function

#### General propose  - `ancillary.py`

* `import_fits`
    - imports either the science data or the full fits object (hdu) from fits files given. Only does one file at a time. Call via map()
    - `biteswap` what in the world is this and why do I use it! Does this affect the data[y,x] paradigm I am use to??
* `get_sn_names`
    - Reads SDSS transient number from files that were renamed & combined. Now we know the numbers associated with the data we currently have. 

#### Plotting - `ploting.py`

* `cdf`

#### Testing - `tests.py`

* I need to test to visually inspect if I found the right thing using `sep`.

### Generating fractional rank variable

#### Defining the galaxy - `defGalaxy.py`

The main goal is to get an elliptical shape that defines a galaxy.
    
* output
    - This should output a `.csv` file whose columns are `SN name, x, y, a, b, theta`
    - file is save at `resources/galaxies_#.csv` where `#` is the sigma used to define galaxy.
        + note larger sigma equals smaller galaxy.

This currently does not do all galaxies just a hard code of SN1415!

* can I define for a galaxy edge by surface brightness (there is a standard galaxy edge definition in this units!)
    - if not, can I find objects by say 5-sigma, but define their edges by 1-sigma over background?

`sep.extract()` settings and why

* thresh
    - why is this based off of background rms? It should be based off of surface brightness. 

#### Calculating $\Delta$

Oh too much work went into doing this by hand 'cause code just would not do.

`resources/shift.csv` contains the results. SN13038 has an ultra faint (or low surface birghtness?) host so it has no data and an SDSS object ID of 0. Just an FYI

More in section [What are we using to do this $\Delta$](#what-are-we-using-to-do-this-delta)

#### Calculating fractional pixel rank - `fractionalRank.py`

* Get SN's pixel position
    - read from SDSS file SN's RA & Dec
    - then convert it to HST's RA & Dec
* Rank Galaxie's pixels ()
    - read in galaxies's properties
        + seems to have an issue with the csv created. (cant seem to have 3 levels of header? Or that the first line has no 'comma')
        + can't have extra newline at the end
        + using `format='ascii.commented_header', header_start=1` seem to be ok, even with new line at end and now needs three header comments!
    - get pixels for this shape
    - import .fits
        + need to change this, get this outside of `fractionalRank.py`. See question and answer below.
    - rank pixles
    - find SN's pixel (should be another function but called here)
    - calculate fractional pixel rank
* save data
    - save galaxie's pixel's values
    - save fractional pixel rank

Questions: I want `import_image` from `defGalaxy.py` in `fractionalRank.py`. Also I think I want `getSNNames` in both (its currently in `findSN.py` that will be deleted.)
Possible Answer: Create a file-IO file and put those two functions in them. AKA make a shared header file for these two mostly separate application files. 

**Question: What to do if SN is outside the ellipse?**

Question: How do we know we have the correct object?

Question: How do we get the whole galaxy? Looking at SN11860, with ds9 in histogram, we see way more galaxy then the resulting shape from sep.

## Questions & Answers

Does it matter if a galaxy is universally dimmer. Rigault will preferentially put that in passive (cause they use an absolute H-alpha measurement) we will not. Our passive/active is relative to the Galaxy. Also we have all the same biases talked about in appendix B of Rigault 2014. 

#### $\Delta$

We need to know the pixel in the HST image where the SN Ia was located. We are using HST for it small pixel size (?). yada yada

* Do we need spherical trig when calculating shift?
    - These are all located in Stripe 82 (RA range and Dec range). Therefore since they are along the equator of this coordinate system the shifts are of order $1~\text{arcsec}$. A simple $\Delta$ is just fine.
* Do we need to worry about the effects of proper motion 
* How accurate is SDSS's position of the SN?
    - Since files are from after everything (have years of data). then the SN's position can be known by looking at it in all available images. Accuracy should be ~0.1 arcsec (unlikely ~0.05 arcsec) says Peter who know these things and the pixel size of SDSS!

#### What are we using to do this $\Delta$?
DS9 7.3.2 centroid's on an object in the fits image (preferably a star) compared to SDSS DR12's location. $\Delta = \text{.fits} - \text{DR12}$, using ICRS frame, saved in a csv file `?.csv`. 

* Using ds9 7.3.2 mac version for the centroid in the fits image. I make the region about the same size as the object (preferably star, but galaxy if there are no stars in the field) but where it does include some sky. Record the center's centroid-ed position. 
    * As [previously published](http://arxiv.org/abs/1307.3539) centroid-ing is good, but different then a Gaussian fit.
* SDSS - should we use DR12 (most recent) or DR7 (the last of SDSS-II where SN survey was from). For a star near SN1415 they appear to be basically identical (different ID's). DR12 has sexagesimal location (that is helpful).
    - http://skyserver.sdss.org/dr7/en/tools/explore/obj.asp?id=588015509808873611
    - http://skyserver.sdss.org/dr12/en/tools/explore/Summary.aspx?id=1237663716016980092
    - SDSS said there was an error in [astonomy in DR8](http://www.sdss3.org/dr8/algorithms/astrometry.php#dr9)
    - http://sdssdp62.fnal.gov/sdsssn/DataRelease/updates.txt says that they posted it in 2014 & used BOSS redshifts (final update in 2015) the DR around then was DR12, July 2014 (dr7 was '08')!
    - to go back: DR12 data can be found at http://skyserver.sdss.org/dr12/en/tools/explore/Summary.aspx?id=1237678617938821531 where 1237678617938821531 is the ID of the object.
* what about proper motion?
    - don't worry about the host they are too far away 
    - what about local stars used for $\Delta$? 
        + I think that it is ok (or the best we can do), because they are recent HST & recent DR12 data. 

## Analysis Method

### Ideas
* Compare fractional pixel rank (~ how many stars are near a SN) to SN parameters ($\Delta M_B$)[^1]
    - Analyses, statistically, if each group (high and low brightness)produce the same SN properties or not.
    - See if we can analyses the continuum of the brightness rather then discreet groups. 
* Repeat with color
* Repeat with SDSS pixel color (found when SN was first analyzed)
* Use what we find to make a better Hubble diagram or apply it somehow to cosmology.
    - sections 5 and 6 of Rigault. 
    - talking about how this would look like a "Mass Step" will likely be enough. 

### Issues
* Low luminosity hosts in SDSS - like SN13038's host
    - What to do with IAU 2006gn (SN13038)? the host is very faint and very low surface brightness and it seems ~0.5 arcmin away
* hot pixels?
    - in SN17886, with linear scale and scale parameters: -1.15171, & 123.66848. 
* How to check if there is a cosmic ray hit on the galaxy? or star?

## Future Work
* include KMOS as an option for IFU, 
* look up hex pack (a pi instrument) on the WIN telescope. 
    * (I can’t spell the instrument or the telescope name).
    * http://www.wiyn.org/Instruments/wiynhexgrad.html


[^fractional pixel rank]: This might be the 3x3 pixel square around the SN rather then just the SN.
[^1]: like figure 6 in Rigault 2013