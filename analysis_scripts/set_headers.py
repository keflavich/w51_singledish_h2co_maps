def set_header_keywords(header, telescope='Arecibo'):
    # Remove starlink-related headers: STARLINK headers were copied, but they
    # are no longer appropriate for the python-modified data
    for starlink_kw in ['HDUCALS1','HDUCLAS2','HDSTYPE']:
        if starlink_kw in header:
            del header[starlink_kw]

    header['OBSERVER'] = 'Adam Ginsburg'
    header['CONTACT'] = 'adam.g.ginsburg@gmail.com'
    header['TELESCOP'] = telescope
    header['ORIGIN'] = 'sdpy: W51 H2CO pipeline'

    return header

