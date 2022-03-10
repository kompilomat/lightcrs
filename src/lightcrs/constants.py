import math

# datastructures, that only need to be loaded once
mgrs_lat_bands = "CDEFGHJKLMNPQRSTUVWXX"
easting_100k_letters = [ 'ABCDEFGH', 'JKLMNPQR', 'STUVWXYZ' ]
northing_100k_Letters = [ 'ABCDEFGHJKLMNPQRSTUV', 'FGHJKLMNPQRSTUVABCDE' ]

# WGS84 constants
semimajor_axis = 6378137.0
flattening = 1 / 298.257223563
scale = 0.9996  



# easting, northing: Karney 2011 Eq 7-14, 29, 35
eccentricity = math.sqrt(flattening * (2 - flattening))
n = flattening / (2 - flattening)
n2 = n**2
n3 = n**3
n4 = n**4
n5 = n**5
n6 = n**6
A = semimajor_axis/(1+n) * (1 + 1/4*n2 + 1/64*n4 + 1/256*n6)

alpha = [ None,
        1/2*n - 2/3*n2 + 5/16*n3 +   41/180*n4 -     127/288*n5 +      7891/37800*n6,
                13/48*n2 -  3/5*n3 + 557/1440*n4 +     281/630*n5 - 1983433/1935360*n6,
                        61/240*n3 -  103/140*n4 + 15061/26880*n5 +   167603/181440*n6,
                                49561/161280*n4 -     179/168*n5 + 6601661/7257600*n6,
                                                    34729/80640*n5 - 3418889/1995840*n6,
                                                                212378941/319334400*n6 ]
# from Karney 2011 Eq 15-22, 36
beta = [ None, 
        1/2*n - 2/3*n2 + 37/96*n3 -    1/360*n4 -   81/512*n5 +    96199/604800*n6,
                1/48*n2 +  1/15*n3 - 437/1440*n4 +   46/105*n5 - 1118711/3870720*n6,
                        17/480*n3 -   37/840*n4 - 209/4480*n5 +      5569/90720*n6,
                                    4397/161280*n4 -   11/504*n5 -  830251/7257600*n6,
                                                4583/161280*n5 -  108847/3991680*n6,
                                                                20648693/638668800*n6 ]


