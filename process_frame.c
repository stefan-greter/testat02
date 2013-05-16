/* Copying and distribution of this file, with or without modification,
 * are permitted in any medium without royalty. This file is offered as-is,
 * without any warranty.
 */

/*! @file process_frame.c
 * @brief Contains the actual algorithm and calculations.
 */

/* Definitions specific to this application. Also includes the Oscar main header file. */
#include "template.h"
#include <string.h>
#include <stdlib.h>
#define DEBUG 1				//1 = Debug mode enabled (Threshold shown in first row of GRAYSCALE picture)

OSC_ERR OscVisDrawBoundingBoxBW(struct OSC_PICTURE *picIn, struct OSC_VIS_REGIONS *regions, uint8 Color);

void ProcessFrame(uint8 *pInputImg)
{
	int c, r;
	int nc = OSC_CAM_MAX_IMAGE_WIDTH/2;
	int siz = sizeof(data.u8TempImage[GRAYSCALE]);

	struct OSC_PICTURE Pic1, Pic2;//we require these structures to use Oscar functions
	struct OSC_VIS_REGIONS ImgRegions;//these contain the foreground objects

	//Variables used for OSTU algorithm
	int k,i;
	uint32 omega0, omega1, my0, my1;
	uint32 hist[256];
	float sigma, max_sigma;
	uint8 otsu_thres;
	memset(hist, 0 , sizeof(hist));


	for(r = 0; r < siz; r+= nc)
	{
		for(c = 0; c < nc; c++)
		{
			//calc histogram
			hist[(int) data.u8TempImage[GRAYSCALE][r+c]]++;
		}
	}


	max_sigma = 0;

	//calc otsu
	for(k = 0; k < 256; k++){
		omega0 = 0;
		omega1 = 0;
		my0 = 0;
		my1 = 0;
		for(i = 0; i <= k ; i++){
			omega0 += hist[i];
			my0 += hist[i] * i;
		}
		for(i = k+1 ; i < 256 ; i++){
			omega1 += hist[i];
			my1 += hist[i] * i;
		}

		sigma = (my0 / omega0 - my1 / omega1 ) * (my0 / omega0 - my1 / omega1 );
		sigma *= omega0 * omega1;

		//set maximum to otsu_thres
		if(sigma > max_sigma){
			max_sigma = sigma;
			otsu_thres = k;
		}
	}



	if(data.ipc.state.nStepCounter == 1)
	{
		/* this is the first time we call this function */
		/* first time we call this; index 1 always has the background image */
		memcpy(data.u8TempImage[BACKGROUND], data.u8TempImage[GRAYSCALE], sizeof(data.u8TempImage[GRAYSCALE]));
		/* set foreground counter to zero */
		memset(data.u8TempImage[FGRCOUNTER], 0, sizeof(data.u8TempImage[FGRCOUNTER]));
	}
	else
	{
		for(r = 0; r < siz; r+= nc)/* we strongly rely on the fact that them images have the same size */
		{
			for(c = 0; c < nc; c++)
			{
				data.u8TempImage[THRESHOLD][r +c ] = data.u8TempImage[GRAYSCALE][r + c] > otsu_thres ? 0 : 0xFF;


#if DEBUG		//plot otsu_thres in first row
				data.u8TempImage[GRAYSCALE][otsu_thres] = 0;
				data.u8TempImage[GRAYSCALE][otsu_thres+1] = 0xFF;
				data.u8TempImage[GRAYSCALE][otsu_thres+2] = 0;
#endif
			}
		}
		/*
		{
			//for debugging purposes we log the background values to console out
			//we chose the center pixel of the image (adaption to other pixel is straight forward)
			int offs = nc*(OSC_CAM_MAX_IMAGE_HEIGHT/2)/2+nc/2;

			OscLog(INFO, "%d %d %d %d %d\n", (int) data.u8TempImage[GRAYSCALE][offs], (int) data.u8TempImage[BACKGROUND][offs], (int) data.u8TempImage[BACKGROUND][offs]-data.ipc.state.nThreshold,
											 (int) data.u8TempImage[BACKGROUND][offs]+data.ipc.state.nThreshold, (int) data.u8TempImage[FGRCOUNTER][offs]);
		}
		*/

		for(r = nc; r < siz-nc; r+= nc)/* we skip the first and last line */
		{
			for(c = 1; c < nc-1; c++)/* we skip the first and last column */
			{
				unsigned char* p = &data.u8TempImage[THRESHOLD][r+c];
				data.u8TempImage[EROSION][r+c] = *(p-nc-1) & *(p-nc) & *(p-nc+1) &
												 *(p-1)    & *p      & *(p+1)    &
												 *(p+nc-1) & *(p+nc) & *(p+nc+1);
			}
		}

		for(r = nc; r < siz-nc; r+= nc)/* we skip the first and last line */
		{
			for(c = 1; c < nc-1; c++)/* we skip the first and last column */
			{
				unsigned char* p = &data.u8TempImage[EROSION][r+c];
				data.u8TempImage[DILATION][r+c] = *(p-nc-1) | *(p-nc) | *(p-nc+1) |
												  *(p-1)    | *p      | *(p+1)    |
												  *(p+nc-1) | *(p+nc) | *(p+nc+1);
			}
		}

		//wrap image DILATION in picture struct
		Pic1.data = data.u8TempImage[DILATION];
		Pic1.width = nc;
		Pic1.height = OSC_CAM_MAX_IMAGE_HEIGHT/2;
		Pic1.type = OSC_PICTURE_GREYSCALE;
		//as well as EROSION (will be used as output)
		Pic2.data = data.u8TempImage[EROSION];
		Pic2.width = nc;
		Pic2.height = OSC_CAM_MAX_IMAGE_HEIGHT/2;
		Pic2.type = OSC_PICTURE_BINARY;//probably has no consequences
		//have to convert to OSC_PICTURE_BINARY which has values 0x01 (and not 0xff)
		OscVisGrey2BW(&Pic1, &Pic2, 0x80, false);

		//now do region labeling and feature extraction
		OscVisLabelBinary( &Pic2, &ImgRegions);
		OscVisGetRegionProperties( &ImgRegions);

		//OscLog(INFO, "number of objects %d\n", ImgRegions.noOfObjects);
		//plot bounding boxes both in gray and dilation image
		Pic2.data = data.u8TempImage[GRAYSCALE];
		OscVisDrawBoundingBoxBW( &Pic2, &ImgRegions, 255);
		OscVisDrawBoundingBoxBW( &Pic1, &ImgRegions, 128);
	}



}


/* Drawing Function for Bounding Boxes; own implementation because Oscar only allows colored boxes; here in Gray value "Color"  */
/* should only be used for debugging purposes because we should not drawn into a gray scale image */
OSC_ERR OscVisDrawBoundingBoxBW(struct OSC_PICTURE *picIn, struct OSC_VIS_REGIONS *regions, uint8 Color)
{
	 uint16 i, o;
	 uint8 *pImg = (uint8*)picIn->data;
	 const uint16 width = picIn->width;
	 for(o = 0; o < regions->noOfObjects; o++)//loop over regions
	 {
		 /* Draw the horizontal lines. */
		 for (i = regions->objects[o].bboxLeft; i < regions->objects[o].bboxRight; i += 1)
		 {
				 pImg[width * regions->objects[o].bboxTop + i] = Color;
				 pImg[width * (regions->objects[o].bboxBottom - 1) + i] = Color;
		 }

		 /* Draw the vertical lines. */
		 for (i = regions->objects[o].bboxTop; i < regions->objects[o].bboxBottom-1; i += 1)
		 {
				 pImg[width * i + regions->objects[o].bboxLeft] = Color;
				 pImg[width * i + regions->objects[o].bboxRight] = Color;
		 }
	 }
	 return SUCCESS;
}


