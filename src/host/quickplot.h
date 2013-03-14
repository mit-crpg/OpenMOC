/*
 * quickplot.h
 *
 *  Created on: Apr 5, 2012
 *      Author: samuelshaner
 */


#ifndef QUICKPLOT_H_
#define QUICKPLOT_H_


#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <utility>
#include <sstream>
#include "Magick++.h"
#include "silo.h"


typedef enum colortypes {
	SCALED,
	RANDOM,
	BLACKWHITE
}colortype;


/* define BitMap struct */
template <typename U>
struct BitMap {
	U* pixels;
	colortype color_type;
	int pixel_x;
	int pixel_y;
	double geom_x;
	double geom_y;
	double center_x;
	double center_y;
	std::list<Magick::Drawable> drawList;
};


/*
 *  function definitions
 */
template <typename U>
void plot(BitMap<U>* bitMap, std::string name, std::string extension);

template <typename U>
void plotSilo(BitMap<U>* bitMap, float* pixMap, std::string type, std::string extension);

template <typename U>
void plotMagick(BitMap<U>* bitMap, float* pixMap, std::string type, std::string extension);

template <typename U>
void copyBitMap(BitMap<U>* bitMap, float* pixMap);

template <typename U>
void normalize(BitMap<U>* bitMap, float* pixMap);

template <typename U>
void getBounds(BitMap<U>* bitMap, float* pixMap, float* bounds);

template <typename U>
void getColor(BitMap<U>* bitMap, float value, float* color);

template <typename U>
void randomize(BitMap<U>* bitMap, float* pixMap);

template <typename U>
void initialize(BitMap<U>* bitMap);

template <typename U, typename V>
void drawLine(BitMap<U>* bitMap, V xIn, V yIn, V xOut, V yOut, U color);

template <typename U, typename V>
int convertToBitmapY(BitMap<U>* bitMap, V y);

template <typename U, typename V>
int convertToBitmapX(BitMap<U>* bitMap, V x);

template <typename U>
void deleteBitMap(BitMap<U>* bitMap);

template <typename U>
void drawText(BitMap<U>* bitMap, std::string text, int x, int y);

template <typename U>
void addScalebar(BitMap<U>* bitMap, float* pixMap, std::list<Magick::Drawable>* drawList);



/*
 * function declarations
 */

/* templated general plot function */
template <typename U>
void plot(BitMap<U>* bitMap, std::string name, std::string extension){

	/* create array to store color values */
	float* pixMap = new float[bitMap->pixel_x * bitMap->pixel_y];
	copyBitMap(bitMap, pixMap);

	/* decide which plot function to call */
	if (extension == "png" || extension == "tiff" || extension == "jpg"){
		plotMagick(bitMap, pixMap, name, extension);
	}
	else if (extension == "pdb" || extension == "h5"){
		plotSilo(bitMap, pixMap, name, extension);
	}

	delete [] pixMap;
}


template <typename U>
void plotSilo(BitMap<U>* bitMap, float* pixMap, std::string name, std::string extension){
	printf("plotting silo mesh...\n");

	/* Create file pointer */
    DBfile *file;

    /* create filename with correct extension */
	std::stringstream string;
	string << name << "." << extension;
	std::string title_str = string.str();
	const char* title = title_str.c_str();

	/* Create file */
	if (extension == "h5"){
		file = DBCreate(title, DB_CLOBBER, DB_LOCAL, "structured mesh bitmap", DB_HDF5);
	}
	else{
		file = DBCreate(title, DB_CLOBBER, DB_LOCAL, "structured mesh bitmap", DB_PDB);
	}

	/* if color_type is RANDOM, randomize bitMapRGB */
	if (bitMap->color_type == RANDOM){
		normalize(bitMap, pixMap);
		randomize(bitMap, pixMap);
	}

    /* create mesh point arrays */
	double mesh_x[bitMap->pixel_x + 1];
	double mesh_y[bitMap->pixel_y + 1];

	/* generate structured mesh */
	for (int i = 0; i < (bitMap->pixel_x + 1); i++){
		mesh_x[i] = (double(i) - double(bitMap->pixel_x)/2.0 + 1.0) * (bitMap->geom_x/double(bitMap->pixel_x));
	}
	for (int i = 0; i < (bitMap->pixel_y + 1); i++){
		mesh_y[i] = (double(bitMap->pixel_y)/2.0 - double(i)) * (bitMap->geom_y/double(bitMap->pixel_y));
	}

	/* descriptions of mesh */
	double *coords[] = {mesh_x, mesh_y};
	int dims[] = {bitMap->pixel_x + 1, bitMap->pixel_y + 1};
	int ndims = 2;

	/* Write structured mesh to file */
	DBPutQuadmesh(file, "quadmesh", NULL, coords, dims, ndims, DB_DOUBLE, DB_COLLINEAR, NULL);

	/* dimensions of mesh */
	int dimsvar[] = {bitMap->pixel_x, bitMap->pixel_y};

	/* description of what is being plotted */
	const char* type_char = name.c_str();

	/* write pixMap data to file */
	DBPutQuadvar1(file, type_char, "quadmesh", pixMap, dimsvar, ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL);

	/* close file */
    DBClose(file);
	printf("done plotting silo mesh...\n");
}


/**
 * Generic function for plotting pixMap in png, tiff, or jpg file
 * using Magick++
 */
template <typename U>
void plotMagick(BitMap<U>* bitMap, float* pixMap, std::string name, std::string extension){
	printf("Writing Magick bitmap...\n");

	/* declare variables */
	float* color = new float[3];
	std::list<Magick::Drawable> drawList;

	/* add scale bar to SCALED bitMaps */
	if (bitMap->color_type == SCALED){
		addScalebar(bitMap, pixMap, &drawList);
	}

	/* normalize pixMap*/
	normalize(bitMap, pixMap);

	/* if color_type is RANDOM, randomize numbers */
	if (bitMap->color_type == RANDOM){
		randomize(bitMap, pixMap);
	}

	/* create image and open for modification */
	Magick::Image image(Magick::Geometry(bitMap->pixel_x,bitMap->pixel_y), "white");
	image.modifyImage();

	/* Make pixel cache */
	Magick::Pixels pixel_cache(image);
	Magick::PixelPacket* pixels;
	pixels = pixel_cache.get(0,0,bitMap->pixel_x,bitMap->pixel_y);

	/* Write pixMapRGB array to Magick pixel_cache */
	for (int y=0;y<bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			/* if pixel is not blank, color pixel */
			if (bitMap->pixels[y * bitMap->pixel_x + x] != -1){

				if (bitMap->color_type == BLACKWHITE){
					*(pixels+(y * bitMap->pixel_x + x)) = Magick::ColorRGB(0, 0, 0);
				}
				else{
					getColor(bitMap, pixMap[y * bitMap->pixel_x + x], color);
					*(pixels+(y * bitMap->pixel_x + x)) = Magick::ColorRGB(color[0], color[1], color[2]);
				}
			}
		}
	}

	/* Sync pixel cache with Magick image */
	pixel_cache.sync();

	/* Draw items in drawList */
	image.draw(bitMap->drawList);

	if (bitMap->color_type == SCALED){
		image.draw(drawList);
	}


	/* create filename with correct extension */
	std::stringstream string;
	string << name << "." << extension;
	std::string title = string.str();

	/* write Magick image to file */
	image.write(title);

	delete [] color;
}


/* copy elements in bitMap to bitMapRGB */
template <typename U>
void copyBitMap(BitMap<U>* bitMap, float* pixMap){

	/* copy bitMap to bitMapRGB */
	for (int y=0;y<bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			pixMap[y * bitMap->pixel_x + x] = (float)bitMap->pixels[y * bitMap->pixel_x + x];
		}
	}

}


/* normalize bitMapRGB to numbers between 0 and 1 */
template <typename U>
void normalize(BitMap<U>* bitMap, float* pixMap){

	float* bounds = new float[2];
	getBounds(bitMap, pixMap, bounds);

	/* copy bitMap to bitMapRGB and normalize */
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x=0;x< bitMap->pixel_x; x++){
			if (pixMap[y * bitMap->pixel_x + x] == -1){
				pixMap[y * bitMap->pixel_x + x] = bounds[0];
			}
			pixMap[y * bitMap->pixel_x + x] = (pixMap[y * bitMap->pixel_x + x] - bounds[0]) /  (bounds[1] - bounds[0]);
		}
	}

	delete [] bounds;
}


/* get min and max bounds of bitMapRGB */
template <typename U>
void getBounds(BitMap<U>* bitMap, float* pixMap, float* bounds){

	bounds[0] = pixMap[0];
	bounds[1] = pixMap[0];

	/* find max */
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			bounds[1] = std::max(bounds[1], pixMap[y * bitMap->pixel_x + x]);
		}
	}

	bounds[0] = bounds[1] - 1e-10;

	/* find min */
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			if (pixMap[y * bitMap->pixel_x + x] != -1){
				bounds[0] = std::min(bounds[0], pixMap[y * bitMap->pixel_x + x]);
			}
		}
	}
}


/* write RGB triplet to color using HOT color scheme */
template <typename U>
void getColor(BitMap<U>* bitMap, float value, float* color){

	if (bitMap->color_type == SCALED){
		if (value <= 1.0/3.0){
			color[0] = 0.0;
			color[1] = 3.0 * value;
			color[2] = 1.0;
		}
		else if (value <= 2.0/3.0){
			color[0] = 3.0 * value - 1.0;
			color[1] = 1.0;
			color[2] = -3.0 * value + 2.0;
		}
		else {
			color[0] = 1.0;
			color[1] = -3.0 * value + 3.0;
			color[2] = 0.0;
		}
	}
	else{
		color[0] = int(value * 100)     % 100 / 100.0;
		color[1] = int(value * 10000)   % 100 / 100.0;
		color[2] = int(value * 1000000) % 100 / 100.0;
	}
}


/* pseudorandomize bitMapRGB with number between 0 and 1 */
template <typename U>
void randomize(BitMap<U>* bitMap, float* pixMap){

	/* make array to store random numbers */
	float* myRandoms = new float[131];

	/* make random numbers */
	srand(1);
	for (int i=0;i< 131; i++){
		myRandoms[i] = rand() / float(RAND_MAX);
	}

	/* randomize bitMapRGB */
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			pixMap[y * bitMap->pixel_x + x] = myRandoms[abs(int(pixMap[y * bitMap->pixel_x + x] / 1e-6)) % 131];
		}
	}

	delete [] myRandoms;
}


/* initialize values to -1 */
template <typename U>
void initialize(BitMap<U>* bitMap){

	try{
		bitMap->pixels = new U[bitMap->pixel_x * bitMap->pixel_y];
	}
	catch (std::exception &e){
		printf("Could not allocate memory for BitMap pixels. "
				"Backtrace:\n%s", e.what());
	}

	/* initialize pixMap to -1 */
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			bitMap->pixels[y * bitMap->pixel_x + x] = -1;
		}
	}

	/* initialize parameters to default values */
	bitMap->center_x = 0;
	bitMap->center_y = 0;
	bitMap->color_type = RANDOM;
}


/**
 * Bresenham's line drawing algorithm. Takes in the start and end coordinates
 * of line (in geometry coordinates), pointer to pixMap array, and line color.
 * "Draws" the line on pixMap array.
 * Taken from "Simplificaiton" code at link below
 * http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
*/
template <typename U, typename V>
void drawLine(BitMap<U>* bitMap, V xIn, V yIn, V xOut, V yOut, U color){

	/* initialize variables */
	int x0, y0, x1,y1;

	/* convert geometry coordinates to bitmap coordinates */
	x0 = convertToBitmapX(bitMap, xIn);
	y0 = convertToBitmapY(bitMap, yIn);
	x1 = convertToBitmapX(bitMap, xOut);
	y1 = convertToBitmapY(bitMap, yOut);

	/* "draw" line on pixMap array */
	int dx = abs(x1-x0);
	int dy = abs(y1-y0);
	int sx, sy;
	if (x0 < x1){
		sx = 1;
	}
	else{
		sx = -1;
	}
	if (y0 < y1){
		sy = 1;
	}
	else{
		sy = -1;
	}
	int error = dx - dy;
	bitMap->pixels[y0 * bitMap->pixel_x + x0] = color;
	bitMap->pixels[y1 * bitMap->pixel_x + x1] = color;
	while (x0 != x1 && y0 != y1){
		bitMap->pixels[y0 * bitMap->pixel_x + x0] = color;
		int e2 = 2 * error;
		if (e2 > -dy){
			error = error - dy;
			x0 = x0 + sx;
		}
		if (e2 < dx){
			error = error + dx;
			y0 = y0 + sy;
		}
	}

}


/**
 * Convert a x value our from geometry coordinates to Bitmap coordinates.
 */
template <typename U, typename V>
int convertToBitmapX(BitMap<U>* bitMap, V x){
	return int((x - bitMap->center_x) * (bitMap->pixel_x - 1) / bitMap->geom_x + (bitMap->pixel_x - 1) / 2.0);
}


/**
 * Convert a y value our from geometry coordinates to Bitmap coordinates.
 */
template <typename U, typename V>
int convertToBitmapY(BitMap<U>* bitMap, V y){
	return int(-(y - bitMap->center_y) * (bitMap->pixel_y - 1) / bitMap->geom_x + (bitMap->pixel_y - 1) / 2.0);
}


/**
 * delete BitMap
 */
template <typename U>
void deleteBitMap(BitMap<U>* bitMap){
	delete [] bitMap->pixels;
	delete bitMap;
}


/**
 * Write text on BitMap
 */
template <typename U>
void drawText(BitMap<U>* bitMap, std::string text, int x, int y){

	/* add item to drawlist	 */
	bitMap->drawList.push_back(Magick::DrawableText(x, y, text));
}


/**
 * add Scalebar to bitMap
 */
template <typename U>
void addScalebar(BitMap<U>* bitMap, float* pixMap, std::list<Magick::Drawable>* drawList){

	drawList->push_back(Magick::DrawableStrokeWidth(0));
	drawList->push_back(Magick::DrawablePointSize(20));
	float* color = new float[3];
	std::stringstream text_stream;
	text_stream.width(9);
	std::string text;

	/* make box */
	drawList->push_back(Magick::DrawableStrokeColor("black"));
	drawList->push_back(Magick::DrawableRectangle(bitMap->pixel_x - 35, 18, bitMap->pixel_x - 10, 122));

	/* make scaled bar */
	for (int y = 20; y < 121; y++){
		getColor(bitMap, 1.0 - (y - 20) / 100.0, color);
		drawList->push_back(Magick::DrawableStrokeColor(Magick::ColorRGB(color[0], color[1], color[2])));
		drawList->push_back(Magick::DrawableLine(bitMap->pixel_x - 33, y, bitMap->pixel_x - 12, y));
	}

	float* bounds = new float[2];
	getBounds(bitMap, pixMap, bounds);
	drawList->push_back(Magick::DrawableStrokeColor("black"));

	/* draw text for max and min */

	text_stream << bounds[1];
	text = text_stream.str();
	text_stream.str("");
	drawList->push_back(Magick::DrawableText(bitMap->pixel_x - 145, 26, text));
	text.clear();

	if (bounds[0] > 1E-15)
		text_stream << bounds[0];
	else
		text_stream << 0;
	text = text_stream.str();
	drawList->push_back(Magick::DrawableText(bitMap->pixel_x - 145, 126, text));

}



#endif /* QUICKPLOT_H_ */
