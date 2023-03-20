#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "math.h"
#include "fitsio.h"
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <vector>
#include "functions.h"

using namespace std;
using namespace boost;

void region(char *fnam,double *exposure,long *axes){
    int nlin=line_num(fnam);
    if (nlin==0) {
        printf("    Unable to read file %s\n",fnam);
    }
    else {
        FILE *ff=fopen(fnam,"r");
        double *xpos=new double[nlin-3];
        double *ypos=new double[nlin-3];
        double *rad=new double[nlin-3];
        double *r2=new double[nlin-3];
        double *angles=new double[nlin-3];
        char c;
        // Go to the 5th line
        int nn=0;
        while (nn<2){
            c = fgetc(ff);
            if(c == '\n') nn++;
        }
        char temp[200];
        fscanf(ff,"%s\n",temp);
        if (!strcmp(temp,"physical")||!strcmp(temp,"wcs")){
            printf("    Region file must be in image coordinates\n");
        }
        else if (strcmp(temp,"image")){
            printf("    Invalid region file %s\n",temp);
        }
        else {
            printf("    %d regions will be ignored\n",nlin-3);
            char type[100];
            for (int i=0;i<nlin-3;i++){
                fscanf(ff,"%[^\n]\n",type);
                char *pch=strtok(type,"-(,)");
                //printf("pch: %s\n",pch);
                if (!strcmp(pch,"circle")) {
                    pch = strtok (NULL, "(,)");
                    xpos[i]=atof(pch)-1;
                    pch = strtok (NULL, "(,)");
                    ypos[i]=atof(pch)-1;
                    pch = strtok (NULL, "(,)");
                    rad[i]=atof(pch);
                    r2[i]=rad[i];
                    angles[i]=0.0;
                }
                else if (!strcmp(pch,"ellipse")) {
                    pch = strtok (NULL, "(,)");
                    xpos[i]=atof(pch)-1;
                    pch = strtok (NULL, "(,)");
                    ypos[i]=atof(pch)-1;
                    pch = strtok (NULL, "(,)");
                    rad[i]=atof(pch);
                    pch = strtok (NULL, "(,)");
                    r2[i]=atof(pch);
                    pch = strtok (NULL, "(,)");
                    angles[i]=atof(pch)*3.14159/180.+3.14159/2.;
                }
            }
            fclose(ff);
            // Modified 01-24-2017 ; improved performance
            for (int ns=0;ns<nlin-3;ns++){
                int boxsize=(int)round(rad[ns]+0.5);
                if (r2[ns]>rad[ns]) {
                    boxsize=(int)round(r2[ns]+0.5);
                }
                int cx=(int)round(xpos[ns]);
                int cy=(int)round(ypos[ns]);
                for (int i=cx-boxsize; i<cx+boxsize+1; i++) {
                    for (int j=cy-boxsize; j<cy+boxsize+1; j++) {
                        if (i>=0 && j>=0 && i<axes[0] && j<axes[1]) {
                            double posx=(i-xpos[ns]);
                            double posy=(j-ypos[ns]);
                            double xtil=cos(angles[ns])*posx+sin(angles[ns])*posy;
                            double ytil=-sin(angles[ns])*posx+cos(angles[ns])*posy;
                            double aoverb=rad[ns]/r2[ns];
                            double dist=aoverb*sqrt(xtil*xtil+ytil*ytil/aoverb/aoverb);
                            if (dist<rad[ns]){
                                exposure[j*axes[0]+i]=0.0;
                            }
                        }
                    }
                }
            }
            delete [] xpos;
            delete [] ypos;
            delete [] rad;
            delete [] r2;
            delete [] angles;
        }
    }
}


int main(int argc, char **argv){
do {
	int status=0;
	if (argc<4 || argc>8){
        printf("Tool to create a Voronoi tessellation from a counts image\n");
		printf("Usage:\n");
		printf("$ voronoi img=image.fits minc=20 outimg=voronoi.fits [expimg=exp.fits backimg=bkg.fits srcreg=src.reg wvt=Y]\n");
        printf("img: input image\n");
        printf("minc: number of counts for Voronoi tessellation\n");
        printf("outimg: output image\n");
        printf("expimg: exposure map\n");
        printf("backimg: background map\n");
        printf("srcreg: region file for source filtering\n");
        printf("wvt: use weighted Voronoi tessellation (default), Diehl & Stattle 06; if set to N, use Cappellari & Copin 03\n");
	}
	else {
        int mincounts;
        char *arg1=new char[200];
        char *arg2=new char[200];
        char *expfile=new char[200];
        char *backfile=new char[200];
        char *srcreg=new char[200];
        bool isimg=false;
        bool isoutimg=false;
        bool isminc=false;
        bool isexp=false;
        bool isback=false;
        bool isreg=false;
        bool isbacksub=false;
        bool wvt=true;
        int narg=1;
		*argv++;
		char comm[200];
        while (narg<argc) {
            string s=*argv;
            vector<string> fields;
            split(fields,s,is_any_of("="));
            char *f1=new char[200];
            char *f2=new char[200];
            strcpy(f1, fields[0].c_str());
            strcpy(f2, fields[1].c_str());
            if (!strcmp(f1,"img")) {
                strcpy(arg1, f2);
                isimg=true;
                printf("Input image is %s\n",arg1);
            }
            else if (!strcmp(f1,"minc")) {
                mincounts=atof(f2);
                printf("We will use %d number of counts for the Voronoi tessellation\n",mincounts);
                isminc=true;
            }
            else if (!strcmp(f1,"outimg")) {
                strcpy(arg2, f2);
                snprintf(comm,100,"if [ -e %s ]; then rm %s; fi",arg2,arg2);
                system(comm);
                isoutimg=true;
                printf("Output image is %s\n",arg2);
            }
            else if (!strcmp(f1,"expimg")) {
                strcpy(expfile, f2);
                isexp=true;
                printf("Exposure map is %s\n",expfile);
            }
            else if (!strcmp(f1,"backimg")) {
                strcpy(backfile, f2);
                isback=true;
                printf("Background map is %s\n",backfile);
            }
            else if (!strcmp(f1,"srcreg")) {
                strcpy(srcreg, f2);
                isreg=true;
                printf("Sources in file %s will be filtered\n",srcreg);
            }
            else if (!strcmp(f1,"wvt")) {
                if (!strcmp(f2,"n")||!strcmp(f2,"N")) {
                    wvt=false;
                    printf("We will use the Cappellari & Copin technique\n");
                }
                else {
                    printf("We will use the Diehl & Stattler technique\n");
                }
            }
            else {
                printf("Unknown parameter %s\n",f1);
            }
            delete [] f1;
            delete [] f2;
            *argv++;
            narg++;
        }
        if (isimg==false) {
            printf("Missing mandatory parameter img\n");
            break;
        }
        if (isminc==false) {
            printf("Missing mandatory parameter minc\n");
            break;
        }
        if (isoutimg==false) {
            printf("Missing mandatory parameter outimg\n");
            break;
        }
        
		int status=0;
		fitsfile *x,*x2;
		fits_open_file(&x,arg1,READONLY,&status);
		if (status!=0){
			printf("Program exiting with status %d\n",status);
			break;
		}
		long axes[2];
		int bitpix,naxis;
		fits_get_img_param(x,2,&bitpix,&naxis,axes,&status);
		long start[2]={1,1};
		int anynul;
        double *img=new double[axes[0]*axes[1]];
        fits_read_pix(x,TDOUBLE,start,axes[0]*axes[1],NULL,img,&anynul,&status);
        if (status!=0) {
            printf("Error %d\n",status);
            return status;
        }
        double pixsize,cdelt1,crval1,crval2,crpix1,crpix2;
        
        fits_read_key(x,TDOUBLE,(char *)"CDELT2",&pixsize,NULL,&status);
        if (status!=0) {
            pixsize=60.;
            status=0;
        }
        else {
            pixsize*=3600.;//arcsec
            fits_read_key(x,TDOUBLE,(char *)"CDELT1",&cdelt1,NULL,&status);
            fits_read_key(x,TDOUBLE,(char *)"CRVAL1",&crval1,NULL,&status);
            fits_read_key(x,TDOUBLE,(char *)"CRVAL2",&crval2,NULL,&status);
            fits_read_key(x,TDOUBLE,(char *)"CRPIX1",&crpix1,NULL,&status);
            fits_read_key(x,TDOUBLE,(char *)"CRPIX2",&crpix2,NULL,&status);
        }
 		fits_close_file(x,&status);
		if (status!=0) {
			printf("Error %d\n",status);
			return status;
		}
        
        double *expo,*backimg;
        if (isexp) {
            fitsfile *y;
            fits_open_file(&y,expfile,READONLY,&status);
            if (status!=0){
                printf("Program exiting with status %d\n",status);
                break;
            }
            expo=new double[axes[0]*axes[1]];
            fits_read_pix(y,TDOUBLE,start,axes[0]*axes[1],NULL,expo,&anynul,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_close_file(y,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
        }
        else {
            expo=new double[axes[0]*axes[1]];
            for (int i=0; i<axes[0]*axes[1]; i++) {
                expo[i]=1.0;
            }
        }
        if (isback) {
            fitsfile *y;
            fits_open_file(&y,backfile,READONLY,&status);
            if (status!=0){
                printf("Program exiting with status %d\n",status);
                break;
            }
            backimg=new double[axes[0]*axes[1]];
            fits_read_pix(y,TDOUBLE,start,axes[0]*axes[1],NULL,backimg,&anynul,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_close_file(y,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
        }
        double back=0.0;
        double eback=0.0;
        if (isreg) {
            region(srcreg,expo,axes);
        }
        
        long npix=axes[0]*axes[1];
        
        double cx,cy;//image pixels of pixel with largest number of counts
        double maximg=maxel(npix,img);
        for (int i=0;i<axes[0];i++){
            for (int j=0;j<axes[1];j++){
                //if (img[j*axes[0]+i]==maximg){
                if (img[j*axes[0]+i]==maximg){
                    if (isexp) {
                        if (expo[j*axes[0]+i]>0.0) {
                            cx=i;
                            cy=j;                            
                        }
                    }
                    else {
                        cx=i;
                        cy=j;                                                    
                    }
                }
            }
        }
        int *binning=new int[npix];//mapping of pixels to bins
        bool *isbinned=new bool[npix];
        long npt=0;
        if (!isexp) {
            npt=npix;
        }
        for (int i=0; i<npix; i++) {
            isbinned[i]=false;
            if (isexp) {
                if (expo[i]>0.0) {
                    binning[i]=0;
                    npt++;
                }
                else {
                    binning[i]=-1;
                }
            }
            else {
                binning[i]=0;
            }
        }
        printf("Number of pixels: %ld\n",npt);
        int cpix=cy*axes[0]+cx;
        //printf("Step 0: %g  %g\n",cx,cy);
        int **tpix=new int*[npix];//array to store pixels present in each bin
        int *npixbin=new int[npix];//number of pixels per bin
        int *binnedpix=new int[npix];//array to store all pixels already binned
        double *cra=new double[npix];//centroids of the bins
        double *cdec=new double[npix];
        bool *isgood=new bool[npix];
        long nok=0;
        long nbinned=0;
        int bin=0;
        int *tarray=new int[npix];
        printf("Creating initial condition...\n");
        while (nok<npt-2) {
            //printf("Working in bin %d...\n",bin);
            tarray[0]=cpix;
            npixbin[bin]=1;
            int nbin=bin+1;
            binning[cpix]=nbin;
            double ncb=totcounts(npixbin[bin],tarray,img);//number of counts in the current pixel
            getmean(npixbin[bin],axes,tarray,cra[bin],cdec[bin]);//Centroid of one pixel
            bool cont=true;
            while (cont) {
                binnedpix[nok]=cpix;
                nok++;
                if (nok%10000==0) {
                    printf("Working with pixel %ld...\n",nok);
                }
                cpix=getclosest(axes,binning,cra[bin],cdec[bin]);
                if (cpix==-1) break;
                bool crit1=topologicalcrit(cpix,npixbin[bin],tarray,axes);
                bool crit2=morphologicalcrit(cpix,npixbin[bin],tarray,axes);//Estimate the roundness of the possible new bin
                bool crit3=uniformitycrit(cpix,npixbin[bin],tarray,img,mincounts);
                if (!crit1 || !crit2 || !crit3){
                    /*if (!crit1) {
                        printf("Crit 1 failed\n");
                    }
                    if (!crit2) {
                        printf("Crit 2 failed\n");
                    }
                    if (!crit3) {
                        printf("Crit 3 failed\n");
                    }*/
                    ncb=totcounts(npixbin[bin],tarray,img);
                    if (ncb>0.8*mincounts) {//The bin has reached the target number of counts
                        //printf("Bin %d has reached the target number of counts\n",nbin);
                        for (int i=0; i<npixbin[bin]; i++) {
                            int pix=tarray[i];
                            isbinned[pix]=true;
                            nbinned++;
                        }
                        isgood[bin]=true;
                    }
                    else {
                        isgood[bin]=false;
                    }
                    tpix[bin]=new int[npixbin[bin]];
                    for (int i=0; i<npixbin[bin]; i++) {
                        tpix[bin][i]=tarray[i];
                    }
                    //printf("npixbin, ncounts: %d  %g\n",npixbin[bin],ncb);
                    bin++;
                    cont=false;
                }
                else {//we accrete
                    int tnb=npixbin[bin];
                    tarray[tnb]=cpix;                                           
                    npixbin[bin]++;
                    binning[cpix]=nbin;
                    //Update centroid of the current bin
                    getmean(npixbin[bin],axes,tarray,cra[bin],cdec[bin]);
                }
            }
            //Calculate new centroid of binned pixels and determine first pixel of next bin
            double tcx,tcy;
            getmean(nok,axes,binnedpix,tcx,tcy);
            /*if (tcx==cx && tcy==cy) {
                nok++;
            }*/
            cx=tcx;
            cy=tcy;
            cpix=getclosest(axes,binning,cx,cy);
            if (cpix==-1) break;
        }

        printf("Adding unsuccessfully binned pixels...\n");
        for (int i=0; i<npix; i++) {
            if (!isbinned[i] && binning[i]>-1) {
                int tbin=findbin(i,axes,bin,isgood,cra,cdec);
                binning[i]=tbin+1;
                int tnb=npixbin[tbin];
                int *transarr=new int[tnb+1];
                for (int j=0; j<tnb; j++) {
                    transarr[j]=tpix[tbin][j];
                }
                transarr[tnb]=i;
                delete [] tpix[tbin];
                tpix[tbin]=transarr;
                npixbin[tbin]++;
            }
        }
        printf("Recomputing centroids...\n");
        int nbgood=0;
        for (int i=0; i<bin; i++) {
            if (isgood[i]) {
                getmean(npixbin[i],axes,tpix[i],cra[i],cdec[i]);
                nbgood++;
            }
        }
        printf("Number of bins for Voronoi tessellation: %d\n",nbgood);
        
        //Sort good bins
        double *cxn=new double[nbgood];
        double *cyn=new double[nbgood];
        int ng=0;
        for (int i=0; i<bin; i++) {
            if (isgood[i]) {
                cxn[ng]=cra[i];
                cyn[ng]=cdec[i];
                ng++;
            }            
        }
        delete [] cra;
        delete [] cdec;
        delete [] isbinned;
        delete [] tarray;
        for (int i=0; i<npix; i++) {
            delete [] tpix[i];
        }
        delete [] tpix;
        delete [] npixbin;
        delete [] isgood;
        
        printf("Compute Voronoi tessellation...\n");
        //Scale factor for the use of the Diehl & Statler 2006 method
        //If wvt=false these are fixed to 1.0
        double *scale=new double[nbgood];
        for (int i=0; i<nbgood; i++) {
            scale[i]=1.0;
        }
        double *cpb=new double[nbgood];
        double *sb=new double[nbgood];
        double *esb=new double[nbgood];
        double totdiff=1e10;
        int maxiter=20;
        int iter=0;
        while (totdiff>0.2 && iter<maxiter) {//threshold for fixed point
            totdiff=compute_voronoi(axes,img,binning,nbgood,scale,wvt,cxn,cyn,cpb)/nbgood;
            printf("Iteration %d, difference %g\n",iter+1,totdiff);
            iter++;
        }
        
        if (ng==nbgood) {
            printf("No empty bins\n");
        }
        else {
            printf("Some bins are empty\n");
        }
        
        printf("Binning completed, now we compute the surface-brightness map...\n");
        double *outimg=new double[npix];
        double *error=new double[npix];
        int *npb=new int[nbgood];
        compute_outimg(binning,nbgood,axes,cxn,cyn,pixsize,img,expo,isexp,backimg,isback,
                       back,eback,isbacksub,outimg,error,npix,scale,sb,esb,npb);

        
        FILE *fc=fopen("centroids.txt","w");
        ng=0;
        for (int i=0; i<nbgood; i++) {
            fprintf(fc,"%d  %g  %g  %g  %g  %g  %d  %g\n",i,cxn[i],cyn[i],cpb[i],sb[i],esb[i],npb[i],scale[i]);
            if (cpb[i]>0.0) {
                ng++;                
            }
        }
        fclose(fc);
        
        

        fits_create_file(&x2,arg2,&status);
		if (status!=0) {
			printf("Error %d\n",status);
			return status;
		}
        fits_create_img(x2,-32,2,axes,&status);
		fits_write_pix(x2,TDOUBLE,start,axes[0]*axes[1],outimg,&status);
		if (status!=0){
			printf("Error %d\n",status);
			return status;
		}	
        double tps=pixsize/3600.;
        fits_write_key(x2,TSTRING,(char *)"CREATOR",(void *)"proffit 1.1",NULL,&status);
        fits_write_key(x2,TSTRING,(char *)"CTYPE1",(void *)"RA---TAN",(char *)"LONGPROJ where LONG can be RA, GLON, ELON and PROJ can be CAR, TAN or AIT",&status);
        fits_write_key(x2,TDOUBLE,(char *)"CRPIX1",&crpix1,(char *)"Pixel at reference point",&status);
        fits_write_key(x2,TDOUBLE,(char *)"CRVAL1",&crval1,(char *)"LONG at the reference value",&status);
        fits_write_key(x2,TSTRING,(char *)"CUNIT1",(void *)"deg",(char *)"Physical units of axis 1",&status);
        fits_write_key(x2,TSTRING,(char *)"CTYPE2",(void *)"DEC--TAN",(char *)"LAT-PROJ where LAT can be DEC, GLAT, ELAT and PROJ can be CAR, TAN or AIT",&status);
        fits_write_key(x2,TDOUBLE,(char *)"CRPIX2",&crpix2,(char *)"Pixel at reference point",&status);
        fits_write_key(x2,TDOUBLE,(char *)"CRVAL2",&crval2,(char *)"LAT at the reference value",&status);
        fits_write_key(x2,TSTRING,(char *)"CUNIT2",(void *)"deg",(char *)"Physical units of axis 2",&status);
        fits_write_key(x2,TDOUBLE,(char *)"CDELT1",&cdelt1,(char *)"Element (1,1) of coordinate transf. matrix (default 1)",&status);
        fits_write_key(x2,TDOUBLE,(char *)"CDELT2",&tps,(char *)"Element (2,2) of coordinate transf. matrix (default 1)",&status);
        fits_write_key(x2,TSTRING,(char *)"RADECSYS",(void *)"FK5",(char *)"Stellar reference frame",&status);
        fits_write_key(x2,TSTRING,(char *)"EQUINOX",(void *)"2000.0",(char *)"Coordinate system equinox",&status);
        fits_create_img(x2,-32,2,axes,&status);
        if (status!=0) {
            printf("Error %d\n",status);
            return status;
        }
        fits_write_pix(x2,TDOUBLE,start,npix,error,&status);
        if (status!=0) {
            printf("    Error %d\n",status);
        }
        fits_write_key(x2,TSTRING,(char *)"CREATOR",(void *)"proffit 1.1",NULL,&status);
        fits_write_key(x2,TSTRING,(char *)"CTYPE1",(void *)"RA---TAN",(char *)"LONGPROJ where LONG can be RA, GLON, ELON and PROJ can be CAR, TAN or AIT",&status);
        fits_write_key(x2,TDOUBLE,(char *)"CRPIX1",&crpix1,(char *)"Pixel at reference point",&status);
        fits_write_key(x2,TDOUBLE,(char *)"CRVAL1",&crval1,(char *)"LONG at the reference value",&status);
        fits_write_key(x2,TSTRING,(char *)"CUNIT1",(void *)"deg",(char *)"Physical units of axis 1",&status);
        fits_write_key(x2,TSTRING,(char *)"CTYPE2",(void *)"DEC--TAN",(char *)"LAT-PROJ where LAT can be DEC, GLAT, ELAT and PROJ can be CAR, TAN or AIT",&status);
        fits_write_key(x2,TDOUBLE,(char *)"CRPIX2",&crpix2,(char *)"Pixel at reference point",&status);
        fits_write_key(x2,TDOUBLE,(char *)"CRVAL2",&crval2,(char *)"LAT at the reference value",&status);
        fits_write_key(x2,TSTRING,(char *)"CUNIT2",(void *)"deg",(char *)"Physical units of axis 2",&status);
        fits_write_key(x2,TDOUBLE,(char *)"CDELT1",&cdelt1,(char *)"Element (1,1) of coordinate transf. matrix (default 1)",&status);
        fits_write_key(x2,TDOUBLE,(char *)"CDELT2",&pixsize,(char *)"Element (2,2) of coordinate transf. matrix (default 1)",&status);
        fits_write_key(x2,TSTRING,(char *)"RADECSYS",(void *)"FK5",(char *)"Stellar reference frame",&status);
        fits_write_key(x2,TSTRING,(char *)"EQUINOX",(void *)"2000.0",(char *)"Coordinate system equinox",&status);
		fits_close_file(x2,&status);
		if (status!=0) {
			printf("Program exiting with status %d\n",status);
			break;
		}

		printf("Program exited successfully\n");
        
        delete [] binning;
        delete [] img;
        delete [] cpb;
        delete [] cxn;
        delete [] cyn;
        delete [] outimg;
        delete [] error;
        delete [] scale;
        delete [] sb;
        delete [] esb;
        if (isexp) delete [] expo;
        if (isback) delete [] backimg;
	}
}
while (0);
}	
