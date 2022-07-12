Segmentation-Based Near-Lossless Compression of Multi-View Cultural Heritage Image Data
======

This is the official implementation of the near-lossless image compression algorithm from the paper "Segmentation-Based Near-Lossless Compression of Multi-View Cultural Heritage Image Data" that has been presented and accepted at *18 Eurographics Workshop on Graphics and Cultural Heritage* in Granada, Spain (actually online...).

Build instructions
------
```
cd planecomp
mkdir build && cd $_
cmake ..
make
```

License & Reference
------
Our program (except of slic.cc) is licensed under the liberal BSD 3-Clause license included as LICENSE.md file.

If you decide to use our code or code based on this project in your application, please make sure to cite our GCH 2020 paper:

```
@inproceedings{10.2312/gch.20201294,
	title = {{Segmentation-Based Near-Lossless Compression of Multi-View Cultural Heritage Image Data}},
	author = {von Buelow, Max and Tausch, Reimar and Knauthe, Volker and Wirth, Tristan and Guthe, Stefan and Santos, Pedro and Fellner, Dieter W.},
	booktitle = {Eurographics Workshop on Graphics and Cultural Heritage},
	editor = {Spagnuolo, Michela and Melero, Francisco Javier},
	year = {2020},
	month = {nov},
	publisher = {The Eurographics Association},
	doi = {10.2312/gch.20201294},
	issn = {2312-6124},
	isbn = {978-3-03868-110-6},
	location = {Online},
	orcid = {0000-0002-0036-319X, 0000-0002-4474-1578, 0000-0001-6993-5099, 0000-0002-2445-9081, 0000-0001-5539-9096, 0000-0003-1813-1714, 0000-0001-7756-0901}
}
```
