#include <vector>
#include <iostream>
#include "Image.h"

namespace Image
{
	ImageAoS::ImageAoS(const Projection::ProjectionParameters &params, ImageOrder order_)
		: rows(params.num_pixels)
		, cols(params.num_pixels)
		, half_width(0.5*params.phantom_radius/(double)params.num_pixels)
		, order (order_)
	{
		pixels.reserve(rows*cols);
		for(int i=0; i<rows; i++)
			for(int j=0; j<cols; j++)
				setPixelCoords(this, j, i, order);
	}
	
	void makeUnitDisk(ImageAoS &image, int radius)
	{
		int r2 = radius*radius;
		for(int i=0; i<image.rows; i++)
			for(int j=0; j<image.cols; j++)
			{
				int d2 = (i-image.rows/2)*(i-image.rows/2) + (j-image.cols/2)*(j-image.cols/2);
				setPixelValue(image, i, j, (int)(d2<=r2));
			}
	}
	
	void printImage(ImageAoS &image)
	{
		for(int i=0; i<image.rows; i++)
		{
			for(int j=0; j<image.cols; j++)
				std::cout << Image::getPixel(image, i, j).value;
			std::cout << "\n";
		}
	}
	
	void printImageCoords(ImageAoS &image_row_order, ImageAoS &image_col_order)
	{
		for(int i = 0; i<image_row_order.rows; i++)
		{
			Projection::PointValue pixel0 = getPixel(image_row_order, i, 0);
			Projection::PointValue pixel1 = getPixel(image_col_order, i, 0);
			std::cout << "(" << pixel0.x  << "," << pixel0.y << "," << pixel0.value << "),(" << pixel1.x  << "," << pixel1.y << "," << pixel1.value << ")\n";
		}
		std::cout << std::endl;
	}
	
	void setPixelValue(ImageAoS &image, int row, int col, int value)
	{
		image.pixels[row*image.rows + col].value = value;
	}
	
	void setPixelCoords(ImageAoS *image, int row, int col, ImageOrder order_)
	{
		if(order_ == ImageOrder::ROW_MAJOR)
		{
			image->pixels[row*image->rows + col].x = ((double)(row - image->rows/2) + 0.5)*image->half_width*2;
			image->pixels[row*image->rows + col].y = ((double)(col - image->cols/2) + 0.5)*image->half_width*2;
		}
		else
		{
			image->pixels[col*image->rows + row].x = ((double)(row - image->rows/2) + 0.5)*image->half_width*2;
			image->pixels[col*image->rows + row].y = ((double)(col - image->cols/2) + 0.5)*image->half_width*2;
		}
	}
	
	Projection::PointValue getPixel(const ImageAoS &image, int row, int col)
	{
		return image.pixels[row*image.rows + col];
	}
}
