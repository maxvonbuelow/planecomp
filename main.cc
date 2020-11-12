#include <iostream>
#include <vector>
#include <random>
#include <limits>
#include <cmath>
#include <core/core.hpp>
#include "slic.h"
#include <unordered_set>
#include <fstream>
#include "image.h"
#include "raw.h"
#include "args.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <omp.h>
#include <atomic>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <chrono>

struct MiddleburryCalib {
	float cam0[9];
	float cam1[9];
	float doffs;
	uint32_t w, h;
	uint32_t ndisp, vmin, vmax;
	float dyavg, dymax;
	float baseline;
};

MiddleburryCalib read_calib(const char *fn)
{
	MiddleburryCalib calib;
	std::ifstream is(fn);
	while (is) {
		std::string key, val;
		std::getline(is, key, '=');
		if (key == "cam0") {
			for (int i = 0; i < 9; ++i) {
				while (is.peek() < '0' || is.peek() > '9') is.get();
				is >> calib.cam0[i];
			}
		} else if (key == "cam1") {
			for (int i = 0; i < 9; ++i) {
				while (is.peek() < '0' || is.peek() > '9') is.get();
				is >> calib.cam1[i];
			}
		} else if (key == "doffs") {
			is >> calib.doffs;
		} else if (key == "width") {
			is >> calib.w;
		} else if (key == "height") {
			is >> calib.h;
		} else if (key == "ndisp") {
			is >> calib.ndisp;
		} else if (key == "vmin") {
			is >> calib.vmin;
		} else if (key == "vmax") {
			is >> calib.vmax;
		} else if (key == "dyavg") {
			is >> calib.dyavg;
		} else if (key == "dymax") {
			is >> calib.dymax;
		} else if (key == "baseline") {
			is >> calib.baseline;
		}
		is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	return calib;
}

struct Plane {
	float a, b,/* c,*/ d;

	bool fit(const float *p0, const float *p1, const float *p2)
	{
		float a0[3], a1[3];
		for (int i = 0; i < 3; ++i) {
			a0[i] = p2[i] - p0[i];
			a1[i] = p1[i] - p0[i];
		}

		a = a0[1] * a1[2] - a0[2] * a1[1];
		b = a0[2] * a1[0] - a0[0] * a1[2];
		float c = a0[0] * a1[1] - a0[1] * a1[0];
		d = -(a * p2[0] + b * p2[1] + c * p2[2]);

		a /= -c;
		b /= -c;
		d /= -c;

		if (std::abs(a) < std::numeric_limits<float>::epsilon() && std::abs(b) < std::numeric_limits<float>::epsilon() && std::abs(c) < std::numeric_limits<float>::epsilon() && std::abs(d) < std::numeric_limits<float>::epsilon()) return false;
		return true;
	}

	float dist(const float *p) const
	{
		return std::abs(get_val_at(p[0], p[1]) - p[2]);
	}
	float get_val_at(float x, float y) const
	{
		return /*-(*/a * x + b * y + d/*) / c*/;
	}
};

void stereo(const MiddleburryCalib &calib, const image_f &im0g, const image_f &im1g, std::vector<float> &errors)
{
	image_f im1p(calib.w, calib.h);
	errors.clear();
	errors.resize(calib.w * calib.h * calib.ndisp, std::numeric_limits<float>::infinity());

	for (uint32_t i = 0; i < calib.ndisp; ++i) {
		for (uint32_t y = 0; y < calib.h; ++y) {
			for (uint32_t x = 0; x < calib.w; ++x) {
				float xx = x - i;
				float yy = y;
				im1p.at2d(x, y) = std::numeric_limits<float>::infinity();
				if (xx < 0 || yy < 0 || xx >= calib.w || yy >= calib.h) continue;

				im1p.at2d(x, y) = im1g.at_lin(xx, yy);
			}
		}
		for (uint32_t y = 0; y < calib.h; ++y) {
			for (uint32_t x = 0; x < calib.w; ++x) {
				float a[9], b[9];
				int n = 0;
				float avg_a = 0, avg_b = 0;
				for (int yy = y == 0; yy < 3 - (y == calib.h - 1); ++yy) {
					for (int xx = x == 0; xx < 3 - (x == calib.w - 1); ++xx) {
						avg_a += a[n] = im0g.at2d(x + xx - 1, y + yy - 1);
						avg_b += b[n] = im1p.at2d(x + xx - 1, y + yy - 1);
						++n;
					}
				}
				avg_a /= n;
				avg_b /= n;

				float var_a = 0, var_b = 0;
				for (int i = 0; i < n; ++i) {
					float aa = a[i] - avg_a;
					float bb = b[i] - avg_b;
					var_a += aa * aa;
					var_b += bb * bb;
				}
				var_a = std::sqrt(var_a / n);
				var_b = std::sqrt(var_b / n);

				float ncc = 0;
				for (int i = 0; i < n; ++i) {
					float aa = (a[i] - avg_a) / var_a;
					float bb = (b[i] - avg_b) / var_b;
					ncc += aa * bb;
				}
				ncc = std::max(0.f, ncc / n);
				float err = 1 - ncc;

				errors[y * calib.w * calib.ndisp + x * calib.ndisp + i] = err;
			}
		}
	}
}

void lsolve(int w, int h, int steps, const float *scores, int *result)
{
	// compute result depth map (get the best depth for each pixel)
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			float best = std::numeric_limits<float>::infinity();
			int beststep = -1;
			for (int s = 0; s < steps; ++s) {
				float score = scores[y * w * steps + x * steps + s];
				if (score <= best) {
					best = score;
					beststep = s;
				}
			}
			result[y * w + x] = beststep;
		}
	}
}

template <typename T>
struct CompImg {
	image<T> img;
	image<T> lossless;
	image_s mask;
	uint32_t ncomp = 0, nuncomp = 0;
	std::vector<Plane> planes;
	std::vector<float> stddevs;
    std::random_device rd;
    std::mt19937 g;
	CompImg() : g(rd())
	{}

	void alloc(uint32_t w, uint32_t h, uint32_t nc)
	{
		img = image<T>(w, h, nc);
		lossless = image<T>(w, h, nc);
		mask = image_s(w, h);
	}

	void recon(int bd)
	{
		img = lossless;
		for (uint32_t y = 0; y < img.height(); ++y) {
			for (uint32_t x = 0; x < img.width(); ++x) {
				uint16_t pid = mask.at2d(x, y);
				if (pid == 0) continue; // already copied

				float stddev = stddevs[pid];
				std::normal_distribution<float> dist(0, stddev);
				for (uint32_t c = 0; c < img.channels(); ++c) {
					float v = planes[pid * img.channels() + c].get_val_at((float)x / (img.width() - 1), (float)y / (img.height() - 1));
					img.at2d(x, y, c) = std::min(std::max(v + dist(g), 0.f), 1.f) * (1 << bd);
				}
			}
		}
	}

	void copy(uint32_t x, uint32_t y, const T *val)
	{
		for (int c = 0; c < img.channels(); ++c) {
			lossless.at2d(x, y, c) = img.at2d(x, y, c) = val[c];
		}
		++nuncomp;
		mask.at2d(x, y) = 0;
	}
	uint32_t sp = 1;
	uint32_t plane(Plane *pl, float stddev)
	{
		if (sp >= planes.size() / img.channels()) {
			planes.resize((sp + 1) * img.channels());
			stddevs.resize(sp + 1);
		}
		for (int c = 0; c < img.channels(); ++c) {
			planes[sp * img.channels() + c] = pl[c];
		}
		stddevs[sp] = stddev;
		return sp++;
	}
	void superpixel(uint32_t x, uint32_t y, uint32_t sp)
	{
		++ncomp;
		const double mean = 0.0;
		const double stddev = stddevs[sp];
		std::normal_distribution<float> dist(mean, stddev);
		for (int c = 0; c < img.channels(); ++c) {
			float v = planes[sp * img.channels() + c].get_val_at((float)x / (img.width() - 1), (float)y / (img.height() - 1));
			img.at2d(x, y, c) = std::min(std::max(v + dist(g), 0.f), 1.f) * (1 << (sizeof(T) * 8));
			lossless.at2d(x, y, c) = (1 << (sizeof(T) * 8 - 1));
		}
		mask.at2d(x, y) = sp;
	}

};

void write_contours(const char *fn, const cv::Mat &in)
{
	cv::Mat img;
	cv::erode(in, img, cv::Mat(), cv::Point(-1, -1), 5, 1, 1);

	image_b out(img.cols, img.rows, 4);
	for (uint32_t y = 0; y < out.height(); ++y) {
		for (uint32_t x = 0; x < out.width(); ++x) {
			cv::Vec3b v = img.at<cv::Vec3b>(y, x);
			if (v[0] == 255 && v[1] == 255 && v[2] == 255) {
				out.at2d(x, y, 0) = 0;
				out.at2d(x, y, 1) = 0;
				out.at2d(x, y, 2) = 0;
				out.at2d(x, y, 3) = 0;
			} else {
				out.at2d(x, y, 0) = 0;
				out.at2d(x, y, 1) = 0;
				out.at2d(x, y, 2) = 0;
				out.at2d(x, y, 3) = 255;
			}
		}
	}
	image_io::save(out, fn);
}

template <typename T>
float planecomp(double thressigma, const image<T> &im0, CompImg<T> &out, int bitdepth, int n, const image_b *mask = nullptr)
{
	uint32_t w = im0.width(), h = im0.height();

	out.alloc(w, h, im0.channels());
	
	cv::Mat inputMatrix(h, w, CV_8UC3, cv::Scalar(0, 0, 0));
	for (uint32_t y = 0; y < h; ++y) {
		for (uint32_t x = 0; x < w; ++x) {
			int c0, c1, c2;
			if (im0.channels() == 3) {
				c0 = im0.at2d(x, y, 0) / (1 << (bitdepth - 8));
				c1 = im0.at2d(x, y, 1) / (1 << (bitdepth - 8));
				c2 = im0.at2d(x, y, 2) / (1 << (bitdepth - 8));
			} else {
				c0 = im0.at2d(x, y, 0) / (1 << (bitdepth - 8));
				c1 = ((int)im0.at2d(x, y, 1) + (int)im0.at2d(x, y, 2)) / 2 / (1 << (bitdepth - 8));
				c2 = im0.at2d(x, y, 3) / (1 << (bitdepth - 8));
			}
			inputMatrix.at<cv::Vec3b>(y, x) = cv::Vec3b(c2, c1, c0);
		}
	}

    SLIC slic;
	slic.GenerateSuperpixels(inputMatrix, n);
	int *result = slic.GetLabel();
	n = 0;
	for (int i = 0; i < w * h; ++i) n = std::max(result[i] + 1, n);
	if (mask) ++n; // all masked pixel will end in patch n-1
	uint32_t tobedone = 0;

	for (uint32_t y = 0; y < h; ++y) {
		for (uint32_t x = 0; x < w; ++x) {
			if (mask && mask->at2d(x, y) != 0) result[y * w + x] = n - 1;
			else ++tobedone;
		}
	}

	std::vector<uint32_t> cluster_cnt(n, 0);
	for (uint32_t y = 0; y < h; ++y) {
		for (uint32_t x = 0; x < w; ++x) {
			int c = result[y * w + x];
			++cluster_cnt[c];
		}
	}
	std::vector<uint32_t> cluster_off(n, 0);
	for (int i = 1; i < n; ++i) {
		cluster_off[i] = cluster_off[i - 1] + cluster_cnt[i - 1];
	}
	cluster_off.push_back(cluster_off.back() + cluster_cnt.back());
	std::vector<std::pair<uint32_t, uint32_t>> cluster_pixels(cluster_off.back());
	std::fill(cluster_cnt.begin(), cluster_cnt.end(), 0);
	for (uint32_t y = 0; y < h; ++y) {
		for (uint32_t x = 0; x < w; ++x) {
			int c = result[y * w + x];
			cluster_pixels[cluster_off[c] + cluster_cnt[c]++] = std::make_pair(x, y);
		}
	}

	float sc = 1 << bitdepth;

	uint32_t x0, y0, x1, y1, x2, y2;
    std::random_device rd;
    std::mt19937 g(rd());
	std::vector<float> mses(n, 0);

	uint32_t done = 0;
	uint32_t nnn = n - (mask ? 1 : 0);
	uint32_t cntsucc = 0;
	for (uint32_t i = 0; i < nnn; ++i) {
		bool last = i == nnn - 1;
		uint32_t o = cluster_off[i];
		uint32_t nn = cluster_off[i + 1] - o;

		float min_peakerr[128]; for (int j = 0; j < 128; ++j) min_peakerr[j] = std::numeric_limits<float>::infinity();
		float min_avgerr[128];
		Plane min_plane[128];
		std::uniform_int_distribution<> dist(0, nn - 1);
		if (nn >= 3) {
			for (int c = 0; c < im0.channels(); ++c) {
				std::vector<T> cache(nn);
				for (int i = 0; i < nn; ++i) {
					std::tie(x0, y0) = cluster_pixels[o + i];
					cache[i] = im0.at2d(x0, y0, c);
				}

				std::vector<Eigen::Vector3d> coords;
				for (uint32_t j = 0; j < nn; ++j) {
					std::tie(x0, y0) = cluster_pixels[o + j];
					coords.emplace_back((float)x0 / (w - 1), (float)y0 / (h - 1), cache[j] / sc);
				}
				Eigen::Matrix<Eigen::Vector3d::Scalar, Eigen::Dynamic, Eigen::Dynamic> coord(3, nn);
				for (size_t i = 0; i < nn; ++i) coord.col(i) = coords[i];
				// calculate centroid
				Eigen::Vector3d centroid(coord.row(0).mean(), coord.row(1).mean(), coord.row(2).mean());
				// subtract centroid
				coord.row(0).array() -= centroid(0); coord.row(1).array() -= centroid(1); coord.row(2).array() -= centroid(2);
				auto svd = coord.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeThinV);
				Eigen::Vector3d plane_normal = svd.matrixU().rightCols<1>();
				double err = svd.singularValues()(2);

				Plane p = { -plane_normal(0) / plane_normal(2), -plane_normal(1) / plane_normal(2), plane_normal.dot(centroid) / plane_normal(2) };
				float peakerr = 0;
				double avgerr = 0;
				for (uint32_t j = 0; j < nn; ++j) {
					std::tie(x0, y0) = cluster_pixels[o + j];
					float p0[] = { (float)x0 / (w - 1), (float)y0 / (h - 1), cache[j] / sc };
					float e = std::abs(p.dist(p0));
					peakerr = std::max(peakerr, e);
					avgerr += e;
				}
				avgerr /= nn;
				if (last) std::cout << "Channel " << c << " error: " << avgerr << std::endl;
				min_plane[c] = p;
				min_avgerr[c] = avgerr;
				min_peakerr[c] = peakerr;
			}
		}
		mses[i] = 0;
		for (int j = 0; j < im0.channels(); ++j) mses[i] += min_avgerr[j];
		mses[i] /= im0.channels();
		float thres = 6.f / 255.f;

		thres = std::numeric_limits<float>::infinity();
		if (min_peakerr[0] > thres || min_peakerr[1] > thres || min_peakerr[2] > thres || mses[i] > thressigma) {
			for (uint32_t j = 0; j < nn; ++j) {
				std::tie(x0, y0) = cluster_pixels[o + j];
				out.copy(x0, y0, im0.ptr2d(x0, y0));
			}
		} else {
			uint32_t sp;
			sp = out.plane(min_plane, mses[i]);
			for (uint32_t j = 0; j < nn; ++j) {
				std::tie(x0, y0) = cluster_pixels[o + j];
				out.superpixel(x0, y0, sp);
			}
			++cntsucc;
		}
		done += nn;
	}

	// coppy masked pixels
	if (mask) {
		uint32_t o = cluster_off[n - 1];
		uint32_t nn = cluster_off[n] - o;
		for (uint32_t j = 0; j < nn; ++j) {
			std::tie(x0, y0) = cluster_pixels[o + j];
			out.copy(x0, y0, im0.ptr2d(x0, y0));
		}
	}

	return (float)cntsucc / n;
}


struct Args {
	std::vector<std::string> in;
	bool oc = false;
	bool decomp = false;
	int n = 1000;
	std::string validate;
	std::string out;
	int bitdepth = 14;
	double threshold = 1.5 / 255.f / 1.f;
	Args(int argc, const char **argv)
	{
		args::parser args(argc, argv, "Blubcompression");
		const int ARG_IN = args.add_nonopt("INPUT"); args.range(1);
		const int ARG_OC = args.add_opt('c', "only-comp", "Do only compression, no stereo");
		const int ARG_DE = args.add_opt('d', "decompress", "Decompress planepng file");
		const int ARG_VA = args.add_opt('v', "validate", "Compare lossless parts with");
		const int ARG_N  = args.add_opt('n', "superpixels", "Number of superpixels");
		const int ARG_OU = args.add_opt('o', "out", "Write output");
		const int ARG_BD = args.add_opt('b', "bitdepth", "Bitdepth for output (default 14)");
		const int ARG_TH = args.add_opt('t', "threshold", "Threshold (default 1.5 / 255.f / 1.f)");

		for (int arg = args.next(); arg != args::parser::end; arg = args.next()) {
			if (arg == ARG_IN)      in.push_back(args.val<std::string>());
			else if (arg == ARG_OC) oc = true;
			else if (arg == ARG_DE) decomp = true;
			else if (arg == ARG_N)  n = args.val<int>();
			else if (arg == ARG_OU) out = args.val<std::string>();
			else if (arg == ARG_VA) validate = args.val<std::string>();
			else if (arg == ARG_BD) bitdepth = args.val<int>();
			else if (arg == ARG_TH) threshold = args.val<double>();
		}
	}
};
long GetFileSize(const char *filename)
{
    struct stat stat_buf;
    int rc = stat(filename, &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}

#define BUFSIZE ((uint32_t)1024 * 1024)
uint32_t write_final(std::ostream &os, std::istream &ismask, uint32_t mask_s, std::istream &islossless, uint32_t lossless_s, const std::vector<Plane> &planes, const std::vector<float> &stddevs)
{
	uint32_t s = stddevs.size(), nc = stddevs.empty() ? 0 : planes.size() / stddevs.size();
	os.write((const char*)&s, 4);
	os.write((const char*)&nc, 4);
	os.write((const char*)planes.data(), planes.size() * sizeof(Plane));
	os.write((const char*)stddevs.data(), stddevs.size() * sizeof(float));

	os.write((const char*)&mask_s, 4);
	uint32_t c = 0;
	char buf[BUFSIZE];
	while (c < mask_s) {
		uint32_t wr = std::min(BUFSIZE, mask_s - c);
		ismask.read(buf, wr);
		os.write(buf, wr);
		c += wr;
	}
	os.write((const char*)&lossless_s, 4);
	c = 0;
	while (c < lossless_s) {
		uint32_t wr = std::min(BUFSIZE, lossless_s - c);
		islossless.read(buf, wr);
		os.write(buf, wr);
		c += wr;
	}
	os.flush();
	return os.tellp();
}
uint32_t write_final(const char *out_file, const char *mask, const char *lossless, const std::vector<Plane> &planes, const std::vector<float> &stddevs)
{
	uint32_t mask_s = GetFileSize(mask), lossless_s = GetFileSize(lossless);
	std::ofstream os(out_file);
	std::ifstream masks(mask), losslesss(lossless);
	write_final(os, masks, mask_s, losslesss, lossless_s, planes, stddevs);
}

void read_final(const char *in_file, const char *mask, const char *lossless, std::vector<Plane> &planes, std::vector<float> &stddevs)
{
	std::ifstream is(in_file);
	uint32_t s, nc;
	is.read((char*)&s, 4);
	is.read((char*)&nc, 4);
	planes.resize(s * nc);
	stddevs.resize(s);
	is.read((char*)planes.data(), planes.size() * sizeof(Plane));
	is.read((char*)stddevs.data(), stddevs.size() * sizeof(float));

	is.read((char*)&s, 4);
	uint32_t c = 0;
	char buf[BUFSIZE];
	std::ofstream osmask(mask);
	while (c < s) {
		uint32_t wr = std::min(BUFSIZE, s - c);
		is.read(buf, wr);
		osmask.write(buf, wr);
		c += wr;
	}

	is.read((char*)&s, 4);
	c = 0;
	std::ofstream oslossless(lossless);
	while (c < s) {
		uint32_t wr = std::min(BUFSIZE, s - c);
		is.read(buf, wr);
		oslossless.write(buf, wr);
		c += wr;
	}
}

void calcavgerr(const image_f &gt, const image_b &noocc, int *result) {
	float avgerr = 0;
	float avgerr_noocc = 0;
	int pxcnt = 0, noocccnt = 0;
	for (uint32_t y = 0; y < gt.height(); ++y) {
		for (uint32_t x = 0; x < gt.width(); ++x) {
			float vgt = gt.at2d(x, y);
			float v = result[x + y * gt.width()];

			float err = std::abs(v - vgt);
			if (std::isinf(err) || std::isnan(err)) continue;
			avgerr += err;
			++pxcnt;
			if (noocc.at2d(x, y) == 255) {
				avgerr_noocc += err;
				++noocccnt;
			}
		}
	}
	avgerr /= pxcnt;
	avgerr_noocc /= noocccnt;

	std::cout << "avgerr: " << avgerr << std::endl << "avgerr noocc: " << avgerr_noocc << std::endl;
}
bool isdir(const char *fn)
{
	struct stat st;
	if(stat(fn,&st) == 0)
		return st.st_mode & S_IFDIR != 0;
	return false;
}
int main(int argc, const char **argv)
{
	Args args(argc, argv);
	std::cout << "Threshold: " << args.threshold << std::endl;
	int n = args.n;

	if (args.oc) {
		std::vector<uint32_t> bytes(512, 0);
		std::vector<uint32_t> bytes_copy(512, 0);
		std::vector<float> timer(512, 0);
		std::vector<float> timerpng(512, 0);
#pragma omp parallel for
		for (int i = 0; i < args.in.size(); ++i) {
			const std::string &file = args.in[i];
			std::string prefix = std::to_string(omp_get_thread_num()) + "_";
			int bitdepth;
			image_s im0;
			im0 = load_raw(file.c_str(), bitdepth);
			std::string maskname = file.substr(0, file.find_last_of('.')) + "_mask.png";
			std::cout << "File: " << file << " mask: " << maskname << std::endl;
			image_b *mask0 = nullptr;
			image_b mask0_dta(im0.width(), im0.height(), 1);
			if (std::ifstream(maskname.c_str()).good()) {
				std::cout << "mask exists" << std::endl;
				image_b mask_big = image_io::load(maskname.c_str());
				for (uint32_t y = 0; y < mask0_dta.height(); ++y) {
					for (uint32_t x = 0; x < mask0_dta.width(); ++x) {
						mask0_dta.at2d(x, y) = (mask_big.at2d(x * 2, y * 2) == 0 || mask_big.at2d(x * 2 + 1, y * 2) == 0 || mask_big.at2d(x * 2, y * 2 + 1) == 0 || mask_big.at2d(x * 2 + 1, y * 2 + 1) == 0) ? 0 : 255;
					}
				}
				mask0 = &mask0_dta;
			}
			std::cout << "Use mask: " << !!mask0 << std::endl;
// #define DEBUG
			std::ofstream final0_png("final0.ppng");
#ifdef DEBUG
			std::ofstream mask0_png("mask0.png"), lossless0_png("lossless0.png"), copy0_png("copy0.png");
#else
			std::stringstream mask0_png, lossless0_png, copy0_png;
#endif
			auto startpng = std::chrono::high_resolution_clock::now();
			image_io::save_png(im0, copy0_png);
			auto endpng = std::chrono::high_resolution_clock::now();
			auto start = std::chrono::high_resolution_clock::now();
			CompImg<uint16_t> comp0;
			planecomp(args.threshold, im0, comp0, bitdepth, n, mask0);
			std::cout << "Finished pc" << std::endl;
			image_io::save_png(comp0.mask, mask0_png);
			image_io::save_png(comp0.lossless, lossless0_png);
			uint32_t mask_s = mask0_png.tellp(), lossless_s = lossless0_png.tellp();
#ifdef DEBUG
			mask0_png.close();
			lossless0_png.close();
			std::ifstream imask0_png("mask0.png"), ilossless0_png("lossless0.png");
			uint32_t cur = write_final(final0_png, imask0_png, mask_s, ilossless0_png, lossless_s, comp0.planes, comp0.stddevs);
#else
			uint32_t cur = write_final(final0_png, mask0_png, mask_s, lossless0_png, lossless_s, comp0.planes, comp0.stddevs);
#endif
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<float> elapsed = end - start;
			float cur_sec = elapsed.count();
			std::chrono::duration<float> elapsedpng = endpng - startpng;
			float cur_secpng = elapsedpng.count();
			timer[omp_get_thread_num()] += cur_sec;
			timerpng[omp_get_thread_num()] += cur_secpng;
			uint32_t cur_copy = copy0_png.tellp();
			bytes[omp_get_thread_num()] += cur;
			bytes_copy[omp_get_thread_num()] += cur_copy;
#pragma omp critical
			{
				std::cout << "### " << file << " ### Cur rate: 1:" << (float)cur_copy / cur << " [" << cur_copy << "/" << cur << "] took " << cur_sec << " sec" << std::endl;
			}
		}
		uint32_t bytes_sum = 0;
		uint32_t bytes_copy_sum = 0;
		float timer_sum = 0, timer_sumpng = 0;
		for (int i = 0; i < omp_get_num_threads(); ++i) {
			bytes_sum += bytes[i];
			bytes_copy_sum += bytes_copy[i];
			timer_sum += timer[i];
			timer_sumpng += timerpng[i];
		}
		std::cout << bytes_sum << " bytes = " << bytes_sum / 1024.f / 1024.f << " MiB" << std::endl;
		std::cout << bytes_copy_sum << " bytes = " << bytes_copy_sum / 1024.f / 1024.f << " MiB" << std::endl;
		std::cout << "Rate: 1:" << (float)bytes_copy_sum / bytes_sum << std::endl;
		std::cout << "Time [sec]: " << timer_sum << std::endl;
		std::cout << "Time PNG [sec]: " << timer_sumpng << std::endl;
		return 0;
	}

	if (args.decomp) {
		CompImg<uint16_t> comp0;
		read_final(args.in[0].c_str(), "__tmp_mask.png", "__tmp_lossless.png", comp0.planes, comp0.stddevs);
		comp0.mask = image_io::load("__tmp_mask.png");
		comp0.lossless = image_io::load("__tmp_lossless.png");
		if (!args.validate.empty()) {
			int bitdepth;
			std::cout << "Starting validation" << std::endl;
			std::cout << "Open " << args.validate.c_str() << std::endl;
			image_s im0 = load_raw(args.validate.c_str(), bitdepth);
			std::cout << im0.width() << " " << im0.height() << std::endl;
			std::cout << comp0.lossless.width() << " " << comp0.lossless.height() << std::endl;
			for (uint32_t y = 0; y < comp0.mask.height(); ++y) {
				for (uint32_t x = 0; x < comp0.mask.width(); ++x) {
					if (comp0.mask.at2d(x, y) != 0) continue;
					for (int j = 0; j < 4; ++j) {
						if (comp0.lossless.at2d(x, y, j) != im0.at2d(x, y, j)) {
							std::cerr << "Found mismatch @" << x << " " << y << " " << j << ": " << comp0.lossless.at2d(x, y, j) << " != " << im0.at2d(x, y, j) << std::endl;
						}
					}
				}
			}
		}
		comp0.recon(args.bitdepth);
		image_io::save(comp0.img, "oc_recon.png");
		image_s bayered(comp0.img.width() * 2, comp0.img.height() * 2);
		for (uint32_t y = 0; y < comp0.img.height(); ++y) {
			for (uint32_t x = 0; x < comp0.img.width(); ++x) {
				bayered.at2d(x * 2, y * 2) = comp0.img.at2d(x, y, 0);
				bayered.at2d(x * 2 + 1, y * 2) = comp0.img.at2d(x, y, 1);
				bayered.at2d(x * 2, y * 2 + 1) = comp0.img.at2d(x, y, 2);
				bayered.at2d(x * 2 + 1, y * 2 + 1) = comp0.img.at2d(x, y, 3);
			}
		}
		image_io::save(bayered, "oc_bayer.tiff");
		return 0;
	}

	for (int i = 0; i < args.in.size(); ++i) {
		std::cout << "###########" << std::endl << args.in[i] << std::endl;
		std::string prefix = std::to_string(i) + "_";

		MiddleburryCalib calib = read_calib((args.in[i] + "/calib.txt").c_str());
		image_b im0 = image_io::load((args.in[i] + "/im0.png").c_str());
		image_b im1 = image_io::load((args.in[i] + "/im1.png").c_str());
		image_b noocc = image_io::load((args.in[i] + "/mask0nocc.png").c_str());
		image_f gt = image_io::load((args.in[i] + "/disp0GT.pfm").c_str());
		std::stringstream comp0_png, comp1_png, mask0_png, mask1_png, lossless0_png, lossless1_png, copy0_png, copy1_png;
		image_io::save_png(im0, copy0_png);
		image_io::save_png(im1, copy1_png);

		calib.vmin = 0;
		calib.vmax = calib.ndisp - 1;

		CompImg<uint8_t> comp0, comp1;
		float p0 = planecomp(args.threshold, im0, comp0, 8, n);
		float p1 = planecomp(args.threshold, im1, comp1, 8, n);
		std::cout << "Ratio success: " << (int)((p0 + p1) / 2.f * 10000 + 0.5) / 100.f << std::endl;

		image_f im0g = image_manip::grayscale(im0);
		image_f im1g = image_manip::grayscale(im1);

		image_f comp0g = image_manip::grayscale(comp0.img);
		image_f comp1g = image_manip::grayscale(comp1.img);
		image_io::save_png(comp0.img, comp0_png);
		image_io::save_png(comp1.img, comp1_png);
		image_io::save_png(comp0.mask, mask0_png);
		image_io::save_png(comp1.mask, mask1_png);

		image_io::save_png(comp0.lossless, lossless0_png);
		image_io::save_png(comp1.lossless, lossless1_png);

		std::vector<float> errors, errorscomp;
		stereo(calib, im0g, im1g, errors);
		stereo(calib, comp0g, comp1g, errorscomp);
		std::vector<int> result(im0g.width() * im0g.height()), resultcomp(im0g.width() * im0g.height());
		lsolve(im0g.width(), im0g.height(), calib.ndisp, errors.data(), result.data());
		lsolve(im0g.width(), im0g.height(), calib.ndisp, errorscomp.data(), resultcomp.data());

		std::cout << "Original:" << std::endl;
		calcavgerr(gt, noocc, result.data());
		std::cout << "Compressed:" << std::endl;
		calcavgerr(gt, noocc, resultcomp.data());

		std::stringstream final0_ss, final1_ss;
		std::ofstream final0_f, final1_f;
		std::ostream *final0_png = &final0_ss, *final1_png = &final1_ss;
		if (!args.out.empty()) {
			final0_f.open((prefix + "_0_" + args.out).c_str());
			final1_f.open((prefix + "_1_" + args.out).c_str());
			final0_png = &final0_f;
			final1_png = &final1_f;
		}
		std::cout << "save 0" << std::endl;
		uint32_t cur0 = write_final(*final0_png, mask0_png, mask0_png.tellp(), lossless0_png, lossless0_png.tellp(), comp0.planes, comp0.stddevs);
		uint32_t cur_copy0 = copy0_png.tellp();

		std::cout << "save 1" << std::endl;
		uint32_t cur1 = write_final(*final1_png, mask1_png, mask1_png.tellp(), lossless1_png, lossless1_png.tellp(), comp1.planes, comp1.stddevs);
		uint32_t cur_copy1 = copy1_png.tellp();


		std::cout << "Quality:" << std::endl;
		std::cout << "Original:" << std::endl;
		calcavgerr(gt, noocc, result.data());
		std::cout << "Compressed:" << std::endl;
		calcavgerr(gt, noocc, resultcomp.data());
		std::cout << "Size LHS:" << std::endl;
		std::cout << cur0 << " bytes = " << cur0 / 1024.f / 1024.f << " MiB" << std::endl;
		std::cout << cur_copy0 << " bytes = " << cur_copy0 / 1024.f / 1024.f << " MiB" << std::endl;
		std::cout << "Rate: 1:" << (float)cur_copy0 / cur0 << std::endl;
		std::cout << "Size RHS:" << std::endl;
		std::cout << cur1 << " bytes = " << cur1 / 1024.f / 1024.f << " MiB" << std::endl;
		std::cout << cur_copy1 << " bytes = " << cur_copy1 / 1024.f / 1024.f << " MiB" << std::endl;
		std::cout << "Rate: 1:" << (float)cur_copy1 / cur1 << std::endl << std::endl;
	}

	return 0;
}
