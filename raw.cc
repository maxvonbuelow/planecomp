#include "image.h"
#include "libraw/libraw.h"
#include <cstring>
#include "raw.h"

template <typename Tin>
void copy_mosaiced2(Tin *in, int ow, int oh, int stride, image_s &out, int w, int h, int nc, int shift, int bitdepth)
{
	uint16_t mask = (1 << bitdepth) - 1;
	for (int j = 0; j < std::min(nc, 4); ++j) {
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				int ox = j & 1, oy = j >> 1;
				int xx = std::min(x, ow / 2 - 1) * 2 + ox, yy = std::min(y, oh / 2 - 1) * 2 + oy;
				out.at2d(x, y, j) = (in[xx + yy * stride] >> shift) & mask;
				if (shift != 0 && (in[xx + yy * stride] & ((1 << shift) - 1)) != 0) {
					std::cout << "Not lossless..." << std::endl;
					std::exit(1);
				}
			}
		}
	}
}

image_s load_raw(const char *filename, int &bitdepth)
{
	/* NOTE: Sorry for all of these ugly hacks. LibRaw seems to have serious issues with compatibility between versions... */
	/* ... we validate that we don't destroy everything on the raw image */
	bool iiq_support = LIBRAW_CHECK_VERSION(0, 18, 0);
	std::size_t l = std::strlen(filename);
	bool is_iiq = filename[l - 4] == '.' && std::tolower(filename[l - 3]) == 'i' && std::tolower(filename[l - 2]) == 'i' && std::tolower(filename[l - 1]) == 'q';
	bool is_cr2 = filename[l - 4] == '.' && std::tolower(filename[l - 3]) == 'c' && std::tolower(filename[l - 2]) == 'r' && std::tolower(filename[l - 1]) == '2';
	int shift = is_iiq && !iiq_support ? 2 : 0;

	LibRaw img;
	img.open_file(filename);
	if (is_cr2) {
		img.imgdata.sizes.left_margin = 104;

		img.imgdata.sizes.width = img.imgdata.sizes.raw_width - img.imgdata.sizes.left_margin * 2;
		img.imgdata.sizes.height = img.imgdata.sizes.raw_height - img.imgdata.sizes.top_margin * 2;
	}
	int ow = img.imgdata.sizes.width;
	int oh = img.imgdata.sizes.height;
	int w = ow / 2;
	int h = oh / 2;
	img.unpack();

	int stride = img.imgdata.sizes.raw_pitch / 2;

	int max = 0;
	uint32_t ored = 0;
	for (int y = img.imgdata.sizes.top_margin; y < img.imgdata.sizes.height; ++y) {
		for (int x = img.imgdata.sizes.left_margin; x < img.imgdata.sizes.width; ++x) {
			int cur = img.imgdata.rawdata.raw_image[x + y * stride] >> shift;
			max = std::max(max, cur);
			ored |= cur;
		}
	}
	bitdepth = 0;
	while (max > 0) {
		max >>= 1;
		++bitdepth;
	}
	if (bitdepth & 1) ++bitdepth;

	image_s image(w, h, 4);
	copy_mosaiced2<uint16_t>(img.imgdata.rawdata.raw_image + img.imgdata.sizes.top_margin * stride + img.imgdata.sizes.left_margin, ow, oh, stride, image, w, h, 4, shift, bitdepth);

	return image;
}
