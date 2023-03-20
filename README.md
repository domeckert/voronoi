# voronoi
2D Voronoi binning for X-ray images

Usage:
$ voronoi img=image.fits minc=20 outimg=voronoi.fits [expimg=exp.fits backimg=bkg.fits srcreg=src.reg wvt=Y]
- img: input image
- minc: number of counts for Voronoi tessellation
- outimg: output image
- expimg: exposure map
- backimg: background map
- srcreg: region file for source filtering
- wvt: use weighted Voronoi tessellation (default), Diehl & Stattle 06; if set to N, use Cappellari & Copin 03
