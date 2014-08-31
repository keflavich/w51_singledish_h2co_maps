from astroquery.simbad import Simbad
from astropy import coordinates
from astropy import units as u

S = Simbad()
S.add_votable_fields('id(Gal|RAFGL|IRAS|1)')
S.add_votable_fields('otype')

tbl = S.query_region(coordinates.Galactic(49*u.deg,-0.3*u.deg), radius=1.5*u.deg)

with open('/Users/adam/work/w51/simbad_HII_regions.reg','w') as outf:
    print >>outf,'global color=white\nfk5'
    for line in tbl[(tbl['OTYPE'] == 'HII')]:
        ra,dec = line['RA'].replace(" ",":"), line['DEC'].replace(" ",":")
        if ra.count(":") == 2 and dec.count(":") == 2:
            print >>outf,"point(%s,%s) # text={%s} point=cross" % (ra,dec,line['ID_Gal_RAFGL_IRAS_1'] if '[' in line['MAIN_ID'] else line['MAIN_ID'])
