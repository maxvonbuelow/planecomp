Segmentation-Based Near-Lossless Compression of Multi-View Cultural Heritage Image Data
======

This is the official implementation of the near-lossless image compression algorithm from the paper "Segmentation-Based Near-Lossless Compression of Multi-View Cultural Heritage Image Data" that has been presented and accepted at *18 Eurographics Workshop on Graphics and Cultural Heritage* in Granada, Spain. 

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
@inproceedings{Buelow2020Imagecomp,
	title = {Segmentation-Based Near-Lossless Compression of Multi-View Cultural Heritage Image Data},
	author = {{von Buelow}, Max and Tausch, Reimar and Knauthe, Volker and Wirth, Tristan and Guthe, Stefan and Santos, Pedro and Fellner, Dieter W.},
	booktitle = {Proceedings of the Eurographics Workshop on Graphics and Cultural Heritage},
	series = {GCH},
	year = {2020},
	month = {November}
}  
```
