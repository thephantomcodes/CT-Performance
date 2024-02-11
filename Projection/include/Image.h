#ifndef IMAGE_H
#define IMAGE_H

#include <vector>
#include "ProjectionParameters.h"

namespace Image
{
	enum class ImageOrder {ROW_MAJOR, COL_MAJOR};
	
	class ImageAoS
	{
		public:
			ImageAoS(const Projection::ProjectionParameters &params, ImageOrder order_);
			int rows;
			int cols;
			double half_width;
			ImageOrder order;
			std::vector<Projection::PointValue> pixels;
	};
	
	void makeUnitDisk(ImageAoS &image, int radius);
	void setPixelValue(ImageAoS &image, int row, int col, int value);
	void printImage(ImageAoS &image);
	void printImageCoords(ImageAoS &image_row_order, ImageAoS &image_col_order);
	void setPixelCoords(ImageAoS *image, int row, int col, ImageOrder order_);
	Projection::PointValue getPixel(const ImageAoS &image, int row, int col);
}

#endif //IMAGE_H
